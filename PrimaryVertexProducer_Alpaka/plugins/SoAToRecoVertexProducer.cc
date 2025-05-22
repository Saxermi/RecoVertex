#include "DataFormats/PortableVertex/interface/VertexHostCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"


#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/SequentialPrimaryVertexFitterAdapter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/AdaptiveChisquarePrimaryVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/WeightedMeanFitter.h"

#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include <algorithm>
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include <cmath>  // for std::isnan

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

#define cputime
#ifdef cputime
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
#endif




/**
   * This plugin takes the SoA portableVertex and converts them to reco::Vertex, for usage within other workflows
   * - consuming set of reco::Tracks and portablevertex SoA
   * - produces a host reco::vertexCollection
 */
class SoAToRecoVertexProducer : public edm::stream::EDProducer<> {

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  const edm::EDGetTokenT<portablevertex::VertexHostCollection> portableVertexToken_;
  const edm::EDGetTokenT<reco::TrackCollection> recoTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> trkToken_; // can only be const when assigen i
  const edm::EDPutTokenT<reco::VertexCollection> recoVertexToken_;
  // added WE, copied from  PrimaryVertexProducer
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
  TrackFilterForPVFindingBase* theTrackFilter;
  TrackClusterizerInZ* theTrackClusterizer;
  // vtx fitting algorithms
  struct algo {
    PrimaryVertexFitterBase* pv_fitter;
    VertexCompatibleWithBeam* vertexSelector;
    std::string label;
    bool useBeamConstraint;
    double minNdof;
  };

  std::vector<algo> algorithms;
  edm::ParameterSet theConfig;
  bool fVerbose;



public:
  SoAToRecoVertexProducer(edm::ParameterSet const& config)
      : portableVertexToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        recoTrackToken_(consumes(config.getParameter<edm::InputTag>("srcTrack"))),
        recoVertexToken_(produces<reco::VertexCollection>()),
	theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
	theConfig(config)
  {
    trkToken_ = consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("TrackLabel"));
    bsToken_ = consumes<reco::BeamSpot>(config.getParameter<edm::InputTag>("beamSpotLabel"));
    theTrackFilter = new TrackFilterForPVFinding(config.getParameter<edm::ParameterSet>("TkFilterParameters"));
    // select and configure the track selection
    std::string trackSelectionAlgorithm =
      config.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
    if (trackSelectionAlgorithm == "filter") {
      theTrackFilter = new TrackFilterForPVFinding(config.getParameter<edm::ParameterSet>("TkFilterParameters"));
    } else {
      throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
    }
    
    // select and configure the track clusterizer
    std::string clusteringAlgorithm =
      config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
    if (clusteringAlgorithm ==  "DA_vect") {
      theTrackClusterizer = new DAClusterizerInZ_vect(
	config.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
    } else {
      throw VertexException("SoAToPrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
    }

    // select and configure the vertex fitters
    std::vector<edm::ParameterSet> vertexCollections =
      config.getParameter<std::vector<edm::ParameterSet> >("vertexCollections");

    for (std::vector<edm::ParameterSet>::const_iterator algoconf = vertexCollections.begin();
	 algoconf != vertexCollections.end();
	 algoconf++) {
      algo algorithm;

      algorithm.label = algoconf->getParameter<std::string>("label");

    // configure the fitter and selector
    std::string fitterAlgorithm = algoconf->getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.pv_fitter = new SequentialPrimaryVertexFitterAdapter(new KalmanVertexFitter());
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      auto fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
      algorithm.pv_fitter = new SequentialPrimaryVertexFitterAdapter(fitter);
    } else if (fitterAlgorithm.empty()) {
      algorithm.pv_fitter = nullptr;
    } else if (fitterAlgorithm == "AdaptiveChisquareVertexFitter") {
      algorithm.pv_fitter = new AdaptiveChisquarePrimaryVertexFitter(algoconf->getParameter<double>("chi2cutoff"),
                                                                     algoconf->getParameter<double>("zcutoff"),
                                                                     algoconf->getParameter<double>("mintrkweight"),
                                                                     false);
    } else if (fitterAlgorithm == "MultiPrimaryVertexFitter") {
      algorithm.pv_fitter = new AdaptiveChisquarePrimaryVertexFitter(algoconf->getParameter<double>("chi2cutoff"),
                                                                     algoconf->getParameter<double>("zcutoff"),
                                                                     algoconf->getParameter<double>("mintrkweight"),
                                                                     true);
    } else if (fitterAlgorithm == "WeightedMeanFitter") {
      algorithm.pv_fitter = new WeightedMeanPrimaryVertexEstimator();
    } else {
      throw VertexException("PrimaryVertexProducer: unknown algorithm: " + fitterAlgorithm);
    }
    algorithm.minNdof = algoconf->getParameter<double>("minNdof");
    algorithm.useBeamConstraint = algoconf->getParameter<bool>("useBeamConstraint");
    algorithm.vertexSelector =
      new VertexCompatibleWithBeam(VertexDistanceXY(), algoconf->getParameter<double>("maxDistanceToBeam"));
    
    algorithms.push_back(algorithm);
    // commmented this in again in attempt to solve the ievent not registered failure
    //produces<reco::VertexCollection>(algorithm.label);
    // implemented  this in again in attempt to solve the ievent not registered failure

    if (!algorithm.label.empty()) {
      produces<reco::VertexCollection>(algorithm.label);
    }
    
    }
#ifdef cputime
    produces<float>("clusteringCPUtime");
#endif
  } // constructor
  
  

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("soaVertex");
    desc.add<edm::InputTag>("srcTrack");
    // >> added from PrimaryVertexProducer
    // vertex collections
    {
      edm::ParameterSetDescription vpsd1;
      vpsd1.add<double>("maxDistanceToBeam", 1.0);
      vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
      vpsd1.add<bool>("useBeamConstraint", false);
      vpsd1.add<std::string>("label", "");
      vpsd1.add<double>("chi2cutoff", 2.5);
      vpsd1.add<double>("zcutoff", 1.0);
      vpsd1.add<double>("mintrkweight", 0.0);
      vpsd1.add<double>("minNdof", 0.0);
      //vpsd1.add<edm::ParameterSetDescription>("vertexTimeParameters", psd_pv_time);
      
      // two default values : with- and without beam constraint
      std::vector<edm::ParameterSet> temp1;
      temp1.reserve(2);
      {
	edm::ParameterSet temp2;
	temp2.addParameter<double>("maxDistanceToBeam", 1.0);
	temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
	temp2.addParameter<bool>("useBeamConstraint", false);
	temp2.addParameter<std::string>("label", "");
	temp2.addParameter<double>("chi2cutoff", 2.5);
	temp2.addParameter<double>("zcutoff", 1.0);
	temp2.addParameter<double>("mintrkweight", 0.);
	temp2.addParameter<double>("minNdof", 0.0);
	//edm::ParameterSet temp_vertexTime;
	//temp_vertexTime.addParameter<std::string>("algorithm", "");
	//temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
	temp1.push_back(temp2);
      }
      {
	edm::ParameterSet temp2;
	temp2.addParameter<double>("maxDistanceToBeam", 1.0);
	temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
	temp2.addParameter<bool>("useBeamConstraint", true);
	temp2.addParameter<std::string>("label", "WithBS");
	temp2.addParameter<double>("chi2cutoff", 2.5);
	temp2.addParameter<double>("zcutoff", 1.0);
	temp2.addParameter<double>("mintrkweight", 0.);
	temp2.addParameter<double>("minNdof", 2.0);
	//edm::ParameterSet temp_vertexTime;
	//temp_vertexTime.addParameter<std::string>("algorithm", "");
	//temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
	temp1.push_back(temp2);
      }
      desc.addVPSet("vertexCollections", vpsd1, temp1);
    }
    desc.addUntracked<bool>("verbose", false);
    {
      edm::ParameterSetDescription psd0;
      TrackFilterForPVFinding::fillPSetDescription(psd0);
      desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
    }
    desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
    desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));

    {
      edm::ParameterSetDescription psd0;
      {
	edm::ParameterSetDescription psd1;
	DAClusterizerInZ_vect::fillPSetDescription(psd1);
	
	//edm::ParameterSetDescription psd2;
	//DAClusterizerInZT_vect::fillPSetDescription(psd2);
	
	edm::ParameterSetDescription psd3;
	GapClusterizerInZ::fillPSetDescription(psd3);
	
	psd0.ifValue(
		     edm::ParameterDescription<std::string>("algorithm", "DA_vect", true),
		     "DA_vect" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkDAClusParameters", psd1, true) or
		     "gap" >> edm::ParameterDescription<edm::ParameterSetDescription>("TkGapClusParameters", psd3, true));
      }
      desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
    }
    // <<< 
    descriptions.addWithDefaultLabel(desc);
  }
};
  
void SoAToRecoVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Book inputs and space for outputs
  const portablevertex::VertexHostCollection& hostVertex = iEvent.get(portableVertexToken_);
  const portablevertex::VertexHostCollection::ConstView& hostVertexView = hostVertex.const_view();
  std::cout << "im in soatorecovertec"<<std::endl;
  std::cout << "hostVertexView[0].nV() "<< hostVertexView[0].nV() << std::endl;

  // get the BeamSpot, it will always be needed, even when not used as a constraint
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken_, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }
  bool validBS = true;
  VertexState beamVertexState(beamSpot);
  if ((beamVertexState.error().cxx() <= 0.) || (beamVertexState.error().cyy() <= 0.) ||
      (beamVertexState.error().czz() <= 0.)) {
    edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors " << beamVertexState.error().matrix();
    validBS = false;
  }

  //auto tks =  iEvent.getHandle(recoTrackToken_).product();
   edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(trkToken_, tks);
  // interface RECO tracks to vertex reconstruction
  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;
  // ... skipped track timing for now
  t_tks = (*theB).build(tks, beamSpot);


  // repackage the clustering result from the gpu
  //struct VertexProtoType
  //{
  //  std::vector<double> Vtx_proto_z_vec;
  //  std::vector<double> Vtx_proto_rho_vec;
   // unsigned int getSize() const { return Vtx_proto_z_vec.size(); }
  
  //};
  DAClusterizerInZ_vect::VertexProtoType  prototypes;
  /*for (int iV = 0; iV < hostVertexView[0].nV(); iV++) {
    if (not(hostVertexView[iV].isGood()))
      continue;
    std::cout << "SoAToRecoVertexProducer: " << iV <<" z=" << hostVertexView[iV].z() << "  rho=" <<  hostVertexView[iV].rho()  <<  " nV=" <<  hostVertexView[iV].nV() << std::endl;
    prototypes.Vtx_proto_z_vec.push_back(hostVertexView[iV].z());
    prototypes.Vtx_proto_rho_vec.push_back(hostVertexView[iV].rho());
  }
  */



  for (int iV = 0; iV < hostVertexView[0].nV(); iV++) {
    int ivv = hostVertexView[iV].order();

    double z_val = hostVertexView[ivv].z();
    double rho_val = hostVertexView[ivv].rho();

    // Always print values, even if they are NaN
    std::cout << "== " << iV << " " << ivv << " of " << hostVertexView[0].nV() << std::endl;
    std::cout << "iVertex " << ivv << ": z=" << z_val << ", rho=" << rho_val << std::endl;

    // Only push values into the prototype if neither z nor rho is NaN
    if (!std::isnan(z_val) && !std::isnan(rho_val)) {
        prototypes.Vtx_proto_z_vec.push_back(z_val);
        prototypes.Vtx_proto_rho_vec.push_back(rho_val);
    }
}

for (size_t i = 0; i < prototypes.Vtx_proto_z_vec.size(); ++i) {
  std::cout << "Prototype " << i << ": z=" << prototypes.Vtx_proto_z_vec[i]
            << ", rho=" << prototypes.Vtx_proto_rho_vec[i] << std::endl;
}
// Zuerst alle mit isGood() == true
std::cout << "Vertices with isGood() == true:" << std::endl;
for (int iV = 0; iV < hostVertexView[0].nV(); ++iV) {
    int ivv = hostVertexView[iV].order();
    if (hostVertexView[ivv].isGood()) {
        std::cout << "iVertex " << ivv << ": z=" << hostVertexView[ivv].z()
                  << ", rho=" << hostVertexView[ivv].rho() << std::endl;
    }
}

// Danach alle mit isGood() == false
std::cout << "Vertices with isGood() == false:" << std::endl;
for (int iV = 0; iV < hostVertexView[0].nV(); ++iV) {
    int ivv = hostVertexView[iV].order();
    if (!hostVertexView[ivv].isGood()) {
        std::cout << "iVertex " << ivv << ": z=" << hostVertexView[ivv].z()
                  << ", rho=" << hostVertexView[ivv].rho() << std::endl;
    }
}


double rhoSum = 0.0;
for (size_t i = 0; i < prototypes.Vtx_proto_rho_vec.size(); ++i) {
    double val = prototypes.Vtx_proto_rho_vec[i];
    if (!std::isnan(val)) {
        rhoSum += val;
    }
}

std::cout << "Sum of all non-NaN rho values: " << rhoSum << std::endl;
prototypes.normalizeRho();
prototypes.sortByZ();

rhoSum = 0.0;
for (size_t i = 0; i < prototypes.Vtx_proto_rho_vec.size(); ++i) {
    double val = prototypes.Vtx_proto_rho_vec[i];
    if (!std::isnan(val)) {
        rhoSum += val;
    }
}
  // carlos: sortByZ
//for (int iV = 0; iV < vertexView[0].nV() ; iV++){
 //   printf("Vertex %i: x=%1.5f, y=%1.5f, z=%1.5f, rho=%1.5f\n",vertexView[iV].x(), vertexView[iV].y(), vertexView[iV].z(), vertexView[iV].rho());
//}

  // select tracks for post clustering
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);

//calculate blockborders
auto block_boundaries = theTrackClusterizer->get_block_boundaries(seltks);

#ifdef cputime
  auto start_clustering = std::chrono::high_resolution_clock::now();
#endif

  // now do the global post-clustering
 // std::vector<TransientVertex>&& clusters = theTrackClusterizer->vertices(seltks, prototypes);
 // std::vector<TransientVertex>&& clusters = theTrackClusterizer->vertices(seltks);
  //std::vector<TransientVertex> clusters;  // to make it compile
  auto* daClusterizer = dynamic_cast<DAClusterizerInZ_vect*>(theTrackClusterizer);
  //if (!daClusterizer) {
  //    throw cms::Exception("LogicError") << "theTrackClusterizer is not a DAClusterizerInZ_vect!";
  
  std::cout << "made it to the global da after DAB "  << std::endl;

  std::vector<TransientVertex> clusters = daClusterizer->Global_DA_after_DAB(prototypes, seltks);
  std::cout << "made it trough the global da after DAB "  << std::endl;


  for (size_t i = 0; i < clusters.size(); ++i) {
    const auto& vertex = clusters[i];
    if (vertex.isValid()) {
        auto pos = vertex.position();
        std::cout << "Cluster " << i << ": "
                  << "x = " << pos.x() << ", "
                  << "y = " << pos.y() << ", "
                  << "z = " << pos.z() << ", ";
    } else {
      std::cout << "Cluster " << i << ": invalid vertex" << std::endl;
    }
}

  
#ifdef cputime
  auto stop_clustering = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int, std::micro> tcpu_clustering = std::chrono::duration_cast<std::chrono::microseconds>(stop_clustering - start_clustering);
#endif

//creating extra info
/*
auto extraInfo = std::make_unique<std::vector<float>>();
for(auto & z : block_boundaries){
  extraInfo->push_back(z);
}
 iEvent.put(std::move(extraInfo), "extraInfo");
*/


   // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);
    std::vector<TransientVertex> pvs;

    if (algorithm->pv_fitter == nullptr) {
      pvs = clusters;
    } else {
      pvs = algorithm->pv_fitter->fit(seltks, clusters, beamSpot, algorithm->useBeamConstraint);
    }

    /*
      if (algorithm->pv_time_estimator != nullptr) {
      algorithm->pv_time_estimator->fill_vertex_times(pvs);
    }
    */
    
    // sort vertices by pt**2  vertex
    if (pvs.size() > 1) {
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }

    // select and convert transient vertices to (reco) vertices
    for (std::vector<TransientVertex>::const_iterator iv = pvs.begin(); iv != pvs.end(); iv++) {
      if (iv->isValid() && (iv->degreesOfFreedom() >= algorithm->minNdof)) {
        reco::Vertex v = *iv;
        if (!validBS || ((*(algorithm->vertexSelector))(v, beamVertexState))) {
          vColl.push_back(v);
        }
      }
    }
    if (fVerbose) {
      edm::LogPrint("PrimaryVertexProducer") << "PrimaryVertexProducer \"" << algorithm->label << "\" contains "
                                             << pvs.size() << " reco::Vertex candidates";
    }

    if (clusters.size() > 2 && clusters.size() > 2 * pvs.size()) {
      edm::LogWarning("PrimaryVertexProducer")
          << "More than 50% of candidate vertices lost (" << pvs.size() << " out of " << clusters.size() << ")";
    }

    if (pvs.empty() && seltks.size() > 5) {
      edm::LogWarning("PrimaryVertexProducer")
          << "No vertex found with " << seltks.size() << " tracks and " << clusters.size() << " vertex candidates";
    }

    if (vColl.empty()) {
      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        edm::LogWarning("PrimaryVertexProducer") << "Zero recostructed vertices, will put reco::Vertex derived from "
                                                    "dummy/fake BeamSpot into Event, BeamSpot has invalid errors: "
                                                 << bse.matrix();
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          edm::LogWarning("PrimaryVertexProducer")
              << "Zero recostructed vertices, will put reco::Vertex derived from BeamSpot into Event.";
        }
      }
    }

    if (fVerbose) {
      int ivtx = 0;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        edm::LogPrint("PrimaryVertexProducer")
            << "recvtx " << std::setw(3) << std::fixed << ivtx++ << " #trk " << std::setw(3) << v->tracksSize()
            << " chi2 " << std::setw(5) << std::setprecision(1) << v->chi2() << " ndof " << std::setw(5)
            << std::setprecision(1) << v->ndof() << " x " << std::setw(7) << std::setprecision(4) << v->position().x()
            << " dx " << std::setw(6) << std::setprecision(4) << v->xError() << " y " << std::setw(7)
            << std::setprecision(4) << v->position().y() << " dy " << std::setw(6) << std::setprecision(4)
            << v->yError() << " z " << std::setw(8) << std::setprecision(4) << v->position().z() << " dz "
            << std::setw(6) << std::setprecision(4) << v->zError();
        if (v->tError() > 0) {
          edm::LogPrint("PrimaryVertexProducer") << " t " << std::setw(6) << std::setprecision(3) << v->t() << " dt "
                                                 << std::setw(6) << std::setprecision(3) << v->tError();
        }
      }
    }
    iEvent.put(std::move(result), algorithm->label);
  }
 
  
#ifdef cputime
  auto clusteringCPUtime = std::make_unique<float>(tcpu_clustering.count() * 1.e-3);
  iEvent.put(std::move(clusteringCPUtime), "clusteringCPUtime");
#endif

}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SoAToRecoVertexProducer);
