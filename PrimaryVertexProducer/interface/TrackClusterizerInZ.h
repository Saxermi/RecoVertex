#ifndef TrackClusterizerInZ_h
#define TrackClusterizerInZ_h

/**\class TrackClusterizerInZ 
 
  Description: interface/base class for track clusterizers that separate event tracks into clusters along the beam line

*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

class TrackClusterizerInZ {
public:
  TrackClusterizerInZ() = default;
  TrackClusterizerInZ(const edm::ParameterSet& conf){};
  virtual std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks) const = 0;
  virtual std::vector<std::vector<reco::TransientTrack> > clusterize(
      const std::vector<reco::TransientTrack>& tracks) const = 0;
  // for testing the clustering in blocks only
  virtual std::vector<float> get_block_boundaries(const std::vector<reco::TransientTrack>& tracks) const{
    std::vector<float> values;
    values.push_back(42.);
    return values;
  }

  virtual ~TrackClusterizerInZ() = default;
};

#endif
