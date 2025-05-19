#include "DataFormats/PortableVertex/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/PortableVertex/interface/VertexHostCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include "BlockAlgo.h"
#include "ClusterizerAlgo.h"
#include "FitterAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  /**
   * This class does vertexing by
   * - consuming set of portablevertex::Track
   * - clusterizing them into track clusters
   * - fitting cluster properties to vertex coordinates
   * - produces a device vertex product (portablevertex::Vertex)
   */
  class PrimaryVertexProducer_Alpaka : public stream::EDProducer<> {
  public:
    PrimaryVertexProducer_Alpaka(edm::ParameterSet const& config) : stream::EDProducer<>(config) {
      trackToken_ = consumes(config.getParameter<edm::InputTag>("TrackLabel"));
      devicePutToken_ = produces();
      blockSize = config.getParameter<int32_t>("blockSize");
      blockOverlap = config.getParameter<double>("blockOverlap");
      clusterParams = {
          .Tmin = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("Tmin"),
          .Tpurge = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("Tpurge"),
          .Tstop = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("Tstop"),
          .vertexSize = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("vertexSize"),
          .coolingFactor =
              config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("coolingFactor"),
          .d0CutOff = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("d0CutOff"),
          .dzCutOff = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("dzCutOff"),
          .uniquetrkweight =
              config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("uniquetrkweight"),
          .uniquetrkminp =
              config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("uniquetrkminp"),
          .zmerge = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("zmerge"),
          .sel_zrange = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("zrange"),
          .convergence_mode =
              config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<int>("convergence_mode"),
          .delta_lowT = config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("delta_lowT"),
          .delta_highT =
              config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<double>("delta_highT")};
      cParams = std::make_shared<portablevertex::ClusterParamsHostCollection>(1, cms::alpakatools::host());
      auto cpview = cParams->view();
      cpview.TMin() = clusterParams.Tmin;
      cpview.Tpurge() = clusterParams.Tpurge;
      cpview.Tstop() = clusterParams.Tstop;
      cpview.vertexSize() = clusterParams.vertexSize;
      cpview.coolingFactor() = clusterParams.coolingFactor;
      cpview.d0CutOff() = clusterParams.d0CutOff;
      cpview.dzCutOff() = clusterParams.dzCutOff;
      cpview.uniquetrkweight() = clusterParams.uniquetrkweight;
      cpview.uniquetrkminp() = clusterParams.uniquetrkminp;
      cpview.zmerge() = clusterParams.zmerge;
      cpview.zrange() = clusterParams.sel_zrange;
      cpview.convergence_mode() = clusterParams.convergence_mode;
      cpview.delta_lowT() = clusterParams.delta_lowT;
      cpview.delta_highT() = clusterParams.delta_highT;
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) {
      printf("Start produce\n");
      const portablevertex::TrackDeviceCollection& inputtracks = iEvent.get(trackToken_);
      printf("Get tracks\n");
      int32_t nT = inputtracks.view().metadata().size();
      printf("Compute block number\n");
      int32_t nBlocks = nT > blockSize ? nT/blockSize : 1;
	      //((nT - 1) / (blockOverlap * blockSize)) 
              //                         : 1;  // If the block size is big enough we process everything at once
      // Now the device collections we still need
      printf("Reserve block memory \n");
      portablevertex::TrackDeviceCollection tracksInBlocks{nBlocks * blockSize, iEvent.queue()};  // As high as needed
      printf("nBlocks %i, nT %i\n", nBlocks, nBlocks * blockSize);
      portablevertex::VertexDeviceCollection deviceVertex{
          1024, iEvent.queue()};  // Hard capped to 1024, though we might want to restrict it for low PU cases

      // run the algorithm
      //// First create the individual blocks
      BlockAlgo blockKernel_{};
      blockKernel_.createBlocks(iEvent.queue(), inputtracks, tracksInBlocks, blockSize, blockOverlap);
      // Need to have the blocks created before launching the next step
      alpaka::wait(iEvent.queue());
      //// Then run the clusterizer per blocks
      ClusterizerAlgo clusterizerKernel_{iEvent.queue(), blockSize};
      clusterizerKernel_.clusterize(iEvent.queue(), tracksInBlocks, deviceVertex, cParams, nBlocks, blockSize);
      // clusterizerKernel_.resplit_tracks(iEvent.queue(), tracksInBlocks, deviceVertex, cParams, nBlocks, blockSize);
      // clusterizerKernel_.reject_outliers(iEvent.queue(), tracksInBlocks, deviceVertex, cParams, nBlocks, blockSize);
      // Need to have all vertex before arbitrating and deciding what we keep
      alpaka::wait(iEvent.queue());
      clusterizerKernel_.arbitrate(iEvent.queue(), tracksInBlocks, deviceVertex, cParams, nBlocks, blockSize);
      alpaka::wait(iEvent.queue());
      //// And then fit
      // Put the vertices in the event as a portable collection
      iEvent.emplace(devicePutToken_, std::move(deviceVertex));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("TrackLabel");
      desc.add<double>("blockOverlap");
      desc.add<int32_t>("blockSize");
      edm::ParameterSetDescription parc0;
      parc0.add<double>("d0CutOff", 3.0);
      parc0.add<double>("Tmin", 2.0);
      parc0.add<double>("delta_lowT", 0.001);
      parc0.add<double>("zmerge", 0.01);
      parc0.add<double>("dzCutOff", 3.0);
      parc0.add<double>("Tpurge", 2.0);
      parc0.add<int32_t>("convergence_mode", 0);
      parc0.add<double>("delta_highT", 0.01);
      parc0.add<double>("Tstop", 0.5);
      parc0.add<double>("coolingFactor", 0.6);
      parc0.add<double>("vertexSize", 0.006);
      parc0.add<double>("uniquetrkweight", 0.8);
      parc0.add<double>("uniquetrkminp", 0.0);
      parc0.add<double>("zrange", 4.0);
      desc.add<edm::ParameterSetDescription>("TkClusParameters", parc0);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<portablevertex::TrackDeviceCollection> trackToken_;
    device::EDPutToken<portablevertex::VertexDeviceCollection> devicePutToken_;
    int32_t blockSize;
    double blockOverlap;
    fitterParameters fitterParams;
    clusterParameters clusterParams;
    std::shared_ptr<portablevertex::ClusterParamsHostCollection> cParams;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PrimaryVertexProducer_Alpaka);
