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

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

/**
   * This plugin takes the SoA portableVertex and converts them to reco::Vertex, for usage within other workflows
   * - consuming set of reco::Tracks and portablevertex SoA
   * - produces a host reco::vertexCollection
 */
class SoAToVertex_tProducer : public edm::stream::EDProducer<> {
public:
  SoAToVertex_tProducer(edm::ParameterSet const& config)
      : portableVertexToken_(consumes(config.getParameter<edm::InputTag>("soaVertex"))),
        recoTrackToken_(consumes(config.getParameter<edm::InputTag>("srcTrack"))),
        recoVertexToken_(produces<reco::VertexCollection>()) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("soaVertex");
    desc.add<edm::InputTag>("srcTrack");

    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  const edm::EDGetTokenT<portablevertex::VertexHostCollection> portableVertexToken_;
  const edm::EDGetTokenT<reco::TrackCollection> recoTrackToken_;
  const edm::EDPutTokenT<reco::VertexCollection> recoVertexToken_;
};

void SoAToVertex_tProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Book inputs and space for outputs
  const portablevertex::VertexHostCollection& hostVertex = iEvent.get(portableVertexToken_);
  const portablevertex::VertexHostCollection::ConstView& hostVertexView = hostVertex.const_view();
  auto tracks =
      iEvent.getHandle(recoTrackToken_)
          .product();  // Note that we need reco::Tracks for building the track Reference vector inside the reco::Vertex

  // This is an annoying conversion as the vertex expects a transient track here, which is a dataformat which we otherwise bypass
  auto result = std::make_unique<reco::VertexCollection>();

  // Do the conversion back to reco::Vertex
  reco::VertexCollection& vColl = (*result);
  for (int iV = 0; iV < hostVertexView[0].nV(); iV++) {
    int ivv = hostVertexView[iV].order();
    printf("== %i %i of %i\n", iV, ivv, hostVertexView[0].nV());
    if (not(hostVertexView[ivv].isGood()))
      continue;
    printf("iVertex %i: z=%1.5f, rho=%1.5f\n", ivv, hostVertexView[ivv].z(), hostVertexView[ivv].rho());
  }
  // And finally put the collection in the event
  iEvent.put(std::move(result));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SoAToVertex_tProducer);
