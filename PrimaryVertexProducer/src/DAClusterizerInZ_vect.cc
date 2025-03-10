#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"


#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include "FWCore/Utilities/interface/isFinite.h"
#include "vdt/vdtMath.h"
// using this for debugging
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
// using f stream and sstream for debugging 

#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
using namespace std;

//#define DEBUG
#ifdef DEBUG
#define DEBUGLEVEL 0
#endif

DAClusterizerInZ_vect::DAClusterizerInZ_vect(const edm::ParameterSet& conf) {
  // hardcoded parameters
  maxIterations_ = 1000;
  mintrkweight_ = 0.5;

  // configurable debug output
#ifdef DEBUG
  zdumpcenter_ = conf.getUntrackedParameter<double>("zdumpcenter", 0.);
  zdumpwidth_ = conf.getUntrackedParameter<double>("zdumpwidth", 20.);
#endif

  // configurable parameters
  double Tmin = conf.getParameter<double>("Tmin");
  double Tpurge = conf.getParameter<double>("Tpurge");
  double Tstop = conf.getParameter<double>("Tstop");
  vertexSize_ = conf.getParameter<double>("vertexSize");
  coolingFactor_ = conf.getParameter<double>("coolingFactor");
  d0CutOff_ = conf.getParameter<double>("d0CutOff");
  dzCutOff_ = conf.getParameter<double>("dzCutOff");
  uniquetrkweight_ = conf.getParameter<double>("uniquetrkweight");
  uniquetrkminp_ = conf.getParameter<double>("uniquetrkminp");
  zmerge_ = conf.getParameter<double>("zmerge");
  sel_zrange_ = conf.getParameter<double>("zrange");
  convergence_mode_ = conf.getParameter<int>("convergence_mode");
  delta_lowT_ = conf.getParameter<double>("delta_lowT");
  delta_highT_ = conf.getParameter<double>("delta_highT");
  runInBlocks_ = conf.getParameter<bool>("runInBlocks");
  block_size_ = conf.getParameter<unsigned int>("block_size");
  overlap_frac_ = conf.getParameter<double>("overlap_frac");
  // temporary kludge: negative overlap_frac values turn on the unblocking test
  if (overlap_frac_ < 0){
    overlap_frac_ = -overlap_frac_;
    unblock_ = true;
    std::cout << "DAClusterizerinZ_vect unblocking test enabled " << std::endl;
  }else{
    unblock_ = false;
  }
  
#ifdef DEBUG
  std::cout << "DAClusterizerinZ_vect: mintrkweight = " << mintrkweight_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: uniquetrkweight = " << uniquetrkweight_ << std::endl;
  std::cout << "DAClusterizerInZ_vect: uniquetrkminp = " << uniquetrkminp_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: zmerge = " << zmerge_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: Tmin = " << Tmin << std::endl;
  std::cout << "DAClusterizerinZ_vect: Tpurge = " << Tpurge << std::endl;
  std::cout << "DAClusterizerinZ_vect: Tstop = " << Tstop << std::endl;
  std::cout << "DAClusterizerinZ_vect: vertexSize = " << vertexSize_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: coolingFactor = " << coolingFactor_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: d0CutOff = " << d0CutOff_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: dzCutOff = " << dzCutOff_ << std::endl;
  std::cout << "DAClusterizerInZ_vect: zrange = " << sel_zrange_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: convergence mode = " << convergence_mode_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: delta_highT = " << delta_highT_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: delta_lowT = " << delta_lowT_ << std::endl;

  std::cout << "DAClusterizerinZ_vect: run in blocks = " << runInBlocks_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: block_size = " << block_size_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: overlap_fraction = " << overlap_frac_ << std::endl;
  std::cout << "DAClusterizerinZ_vect: DEBUGLEVEL " << DEBUGLEVEL << std::endl;
#endif

  if (convergence_mode_ > 1) {
    edm::LogWarning("DAClusterizerinZ_vect")
        << "DAClusterizerInZ_vect: invalid convergence_mode " << convergence_mode_ << "  reset to default " << 0;
    convergence_mode_ = 0;
  }

  if (Tmin == 0) {
    betamax_ = 1.0;
    edm::LogWarning("DAClusterizerinZ_vect")
        << "DAClusterizerInZ_vect: invalid Tmin " << Tmin << "  reset to default " << 1. / betamax_;
  } else {
    betamax_ = 1. / Tmin;
  }

  if ((Tpurge > Tmin) || (Tpurge == 0)) {
    edm::LogWarning("DAClusterizerinZ_vect")
        << "DAClusterizerInZ_vect: invalid Tpurge " << Tpurge << "  set to " << Tmin;
    Tpurge = Tmin;
  }
  betapurge_ = 1. / Tpurge;

  if ((Tstop > Tpurge) || (Tstop == 0)) {
    edm::LogWarning("DAClusterizerinZ_vect")
        << "DAClusterizerInZ_vect: invalid Tstop " << Tstop << "  set to  " << max(1., Tpurge);
    Tstop = max(1., Tpurge);
  }
  betastop_ = 1. / Tstop;
}

namespace {
  inline double local_exp(double const& inp) { return vdt::fast_exp(inp); }

  inline void local_exp_list_range(double const* __restrict__ arg_inp,
                                   double* __restrict__ arg_out,
                                   const int kmin,
                                   const int kmax) {
    for (auto i = kmin; i != kmax; ++i)
      arg_out[i] = vdt::fast_exp(arg_inp[i]);
  }

}  // namespace

void DAClusterizerInZ_vect::verify(const vertex_t& v, const track_t& tks, unsigned int nv, unsigned int nt) const {
  if (!(nv == 999999)) {
    assert(nv == v.getSize());
  } else {
    nv = v.getSize();
  }

  if (!(nt == 999999)) {
    assert(nt == tks.getSize());
  } else {
    nt = tks.getSize();
  }

  assert(v.zvtx_vec.size() == nv);
  assert(v.rho_vec.size() == nv);
  assert(v.swz_vec.size() == nv);
  assert(v.exp_arg_vec.size() == nv);
  assert(v.exp_vec.size() == nv);
  assert(v.se_vec.size() == nv);
  assert(v.swz_vec.size() == nv);
  assert(v.swE_vec.size() == nv);

  assert(v.zvtx == &v.zvtx_vec.front());
  assert(v.rho == &v.rho_vec.front());
  assert(v.exp_arg == &v.exp_arg_vec.front());
  assert(v.sw == &v.sw_vec.front());
  assert(v.swz == &v.swz_vec.front());
  assert(v.se == &v.se_vec.front());
  assert(v.swE == &v.swE_vec.front());

  for (unsigned int k = 0; k < nv - 1; k++) {
    if (v.zvtx_vec[k] <= v.zvtx_vec[k + 1])
      continue;
    cout << " Z, cluster z-ordering assertion failure   z[" << k << "] =" << v.zvtx_vec[k] << "    z[" << k + 1
         << "] =" << v.zvtx_vec[k + 1] << endl;
  }

  assert(nt == tks.zpca_vec.size());
  assert(nt == tks.dz2_vec.size());
  assert(nt == tks.tt.size());
  assert(nt == tks.tkwt_vec.size());
  assert(nt == tks.sum_Z_vec.size());
  assert(nt == tks.kmin.size());
  assert(nt == tks.kmax.size());

  assert(tks.zpca == &tks.zpca_vec.front());
  assert(tks.dz2 == &tks.dz2_vec.front());
  assert(tks.tkwt == &tks.tkwt_vec.front());
  assert(tks.sum_Z == &tks.sum_Z_vec.front());

  for (unsigned int i = 0; i < nt; i++) {
    if ((tks.kmin[i] < tks.kmax[i]) && (tks.kmax[i] <= nv))
      continue;
    cout << "track vertex range assertion failure" << i << "/" << nt << "   kmin,kmax=" << tks.kmin[i] << ", "
         << tks.kmax[i] << "  nv=" << nv << endl;
  }

  for (unsigned int i = 0; i < nt; i++) {
    assert((tks.kmin[i] < tks.kmax[i]) && (tks.kmax[i] <= nv));
  }
}

DAClusterizerInZ_vect::track_t DAClusterizerInZ_vect::fill(const vector<reco::TransientTrack>& tracks) const {
  // prepare track data for clustering
  track_t tks;
  double sumtkwt = 0.;
  for (auto it = tracks.begin(); it != tracks.end(); it++) {
    if (!(*it).isValid())
      continue;
    double t_tkwt = 1.;
    double t_z = ((*it).stateAtBeamLine().trackStateAtPCA()).position().z();
    if (std::fabs(t_z) > 1000.)
      continue;
    auto const& t_mom = (*it).stateAtBeamLine().trackStateAtPCA().momentum();
    //  get the beam-spot
    reco::BeamSpot beamspot = (it->stateAtBeamLine()).beamSpot();
    double t_dz2 = std::pow((*it).track().dzError(), 2)  // track errror
                   + (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
                         std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
                   + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
    t_dz2 = 1. / t_dz2;
    if (edm::isNotFinite(t_dz2) || t_dz2 < std::numeric_limits<double>::min())
      continue;
    if (d0CutOff_ > 0) {
      Measurement1D atIP = (*it).stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
      t_tkwt = 1. / (1. + local_exp(std::pow(atIP.value() / atIP.error(), 2) -
                                    std::pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
      if (edm::isNotFinite(t_tkwt) || t_tkwt < std::numeric_limits<double>::epsilon()) {
        edm::LogWarning("DAClusterizerinZ_vect") << "rejected track t_tkwt " << t_tkwt;
        continue;  // usually is > 0.99
      }
    }
    tks.addItemSorted(t_z, t_dz2, &(*it), t_tkwt);
    sumtkwt += t_tkwt;
  }

  if (sumtkwt > 0) {
    tks.extractRaw();
    tks.osumtkwt = 1. / sumtkwt;
  } else {
    tks.osumtkwt = 0.;
  }

#ifdef DEBUG
  if (DEBUGLEVEL > 0) {
    std::cout << "Track count (Z) " << tks.getSize() << std::endl;
  }
#endif
  return tks;
}

namespace {
  inline double Eik(double t_z, double k_z, double t_dz2) { return std::pow(t_z - k_z, 2) * t_dz2; }
}  // namespace

void DAClusterizerInZ_vect::set_vtx_range(double beta, track_t& gtracks, vertex_t& gvertices) const {
  const unsigned int nv = gvertices.getSize();
  const unsigned int nt = gtracks.getSize();

  if (nv == 0) {
    edm::LogWarning("DAClusterizerinZ_vect") << "empty cluster list in set_vtx_range";
    return;
  }

  for (auto itrack = 0U; itrack < nt; ++itrack) {
    double zrange = max(sel_zrange_ / sqrt(beta * gtracks.dz2[itrack]), zrange_min_);

    double zmin = gtracks.zpca[itrack] - zrange;
    unsigned int kmin = min(nv - 1, gtracks.kmin[itrack]);
    // find the smallest vertex_z that is larger than zmin
    if (gvertices.zvtx[kmin] > zmin) {
      while ((kmin > 0) && (gvertices.zvtx[kmin - 1] > zmin)) {
        kmin--;
      }
    } else {
      while ((kmin < (nv - 1)) && (gvertices.zvtx[kmin] < zmin)) {
        kmin++;
      }
    }

    double zmax = gtracks.zpca[itrack] + zrange;
    unsigned int kmax = min(nv - 1, gtracks.kmax[itrack] - 1);
    // note: kmax points to the last vertex in the range, while gtracks.kmax points to the entry BEHIND the last vertex
    // find the largest vertex_z that is smaller than zmax
    if (gvertices.zvtx[kmax] < zmax) {
      while ((kmax < (nv - 1)) && (gvertices.zvtx[kmax + 1] < zmax)) {
        kmax++;
      }
    } else {
      while ((kmax > 0) && (gvertices.zvtx[kmax] > zmax)) {
        kmax--;
      }
    }

    if (kmin <= kmax) {
      gtracks.kmin[itrack] = kmin;
      gtracks.kmax[itrack] = kmax + 1;
    } else {
      gtracks.kmin[itrack] = max(0U, min(kmin, kmax));
      gtracks.kmax[itrack] = min(nv, max(kmin, kmax) + 1);
    }
  }
}

void DAClusterizerInZ_vect::clear_vtx_range(track_t& gtracks, vertex_t& gvertices) const {
  const unsigned int nt = gtracks.getSize();
  const unsigned int nv = gvertices.getSize();
  for (auto itrack = 0U; itrack < nt; ++itrack) {
    gtracks.kmin[itrack] = 0;
    gtracks.kmax[itrack] = nv;
  }
}

double DAClusterizerInZ_vect::update(
    double beta, track_t& gtracks, vertex_t& gvertices, const double rho0, const bool updateTc) const {
  // update weights and vertex positions
  // returns the maximum of changes of vertex positions
  // sums needed for Tc are only updated if updateTC == true

  const unsigned int nt = gtracks.getSize();
  const unsigned int nv = gvertices.getSize();
  auto osumtkwt = gtracks.osumtkwt;

  double Z_init = 0;
  if (rho0 > 0) {
    Z_init = rho0 * local_exp(-beta * dzCutOff_ * dzCutOff_);  // cut-off
  }

  // define kernels
  auto kernel_calc_exp_arg_range = [beta](const unsigned int itrack,
                                          track_t const& tracks,
                                          vertex_t const& vertices,
                                          const unsigned int kmin,
                                          const unsigned int kmax) {
    const double track_z = tracks.zpca[itrack];
    const double botrack_dz2 = -beta * tracks.dz2[itrack];

    // auto-vectorized
    for (unsigned int ivertex = kmin; ivertex < kmax; ++ivertex) {
      auto mult_res = track_z - vertices.zvtx[ivertex];
      vertices.exp_arg[ivertex] = botrack_dz2 * (mult_res * mult_res);
    }
  };

  auto kernel_add_Z_range = [Z_init](
                                vertex_t const& vertices, const unsigned int kmin, const unsigned int kmax) -> double {
    double ZTemp = Z_init;
    for (unsigned int ivertex = kmin; ivertex < kmax; ++ivertex) {
      ZTemp += vertices.rho[ivertex] * vertices.exp[ivertex];
    }
    return ZTemp;
  };

  auto kernel_calc_normalization_range = [updateTc](const unsigned int track_num,
                                                    track_t& tracks,
                                                    vertex_t& vertices,
                                                    const unsigned int kmin,
                                                    const unsigned int kmax) {
    auto o_trk_sum_Z = tracks.tkwt[track_num] / tracks.sum_Z[track_num];
    auto o_trk_dz2 = tracks.dz2[track_num];
    auto tmp_trk_z = tracks.zpca[track_num];

    // auto-vectorized
    if (updateTc) {
#pragma GCC ivdep
      for (unsigned int k = kmin; k < kmax; ++k) {
        vertices.se[k] += vertices.exp[k] * o_trk_sum_Z;
        auto w = vertices.rho[k] * vertices.exp[k] * (o_trk_sum_Z * o_trk_dz2);
        vertices.sw[k] += w;
        vertices.swz[k] += w * tmp_trk_z;
        vertices.swE[k] += w * vertices.exp_arg[k];
      }
    } else {
      // same loop but without updating sWE
#pragma GCC ivdep
      for (unsigned int k = kmin; k < kmax; ++k) {
        vertices.se[k] += vertices.exp[k] * o_trk_sum_Z;
        auto w = vertices.rho[k] * vertices.exp[k] * (o_trk_sum_Z * o_trk_dz2);
        vertices.sw[k] += w;
        vertices.swz[k] += w * tmp_trk_z;
      }
    }
  };

  if (updateTc) {
    for (auto ivertex = 0U; ivertex < nv; ++ivertex) {
      gvertices.se[ivertex] = 0.0;
      gvertices.sw[ivertex] = 0.0;
      gvertices.swz[ivertex] = 0.0;
      gvertices.swE[ivertex] = 0.0;
    }
  } else {
    for (auto ivertex = 0U; ivertex < nv; ++ivertex) {
      gvertices.se[ivertex] = 0.0;
      gvertices.sw[ivertex] = 0.0;
      gvertices.swz[ivertex] = 0.0;
    }
  }

  // loop over tracks
  for (auto itrack = 0U; itrack < nt; ++itrack) {
    const unsigned int kmin = gtracks.kmin[itrack];
    const unsigned int kmax = gtracks.kmax[itrack];

    kernel_calc_exp_arg_range(itrack, gtracks, gvertices, kmin, kmax);
    local_exp_list_range(gvertices.exp_arg, gvertices.exp, kmin, kmax);
    gtracks.sum_Z[itrack] = kernel_add_Z_range(gvertices, kmin, kmax);

    if (edm::isNotFinite(gtracks.sum_Z[itrack]))
      gtracks.sum_Z[itrack] = 0.0;

    if (gtracks.sum_Z[itrack] > 1.e-100) {
      kernel_calc_normalization_range(itrack, gtracks, gvertices, kmin, kmax);
    }
  }

  // (un-)apply the factor -beta which is needed in exp_arg, but not in swE
  if (updateTc) {
    auto obeta = -1. / beta;
    for (auto ivertex = 0U; ivertex < nv; ++ivertex) {
      gvertices.swE[ivertex] *= obeta;
    }
  }

  // now update z and rho
  auto kernel_calc_z = [osumtkwt, nv](vertex_t& vertices) -> double {
    double delta = 0;
    // does not vectorize
    for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
      if (vertices.sw[ivertex] > 0) {
        auto znew = vertices.swz[ivertex] / vertices.sw[ivertex];
        // prevents from vectorizing if
        delta = max(std::abs(vertices.zvtx[ivertex] - znew), delta);
        vertices.zvtx[ivertex] = znew;
      }
    }

    for (unsigned int ivertex = 0; ivertex < nv; ++ivertex)
      vertices.rho[ivertex] = vertices.rho[ivertex] * vertices.se[ivertex] * osumtkwt;

    return delta;
  };

  double delta = kernel_calc_z(gvertices);

  // return how much the prototypes moved
  return delta;
}

unsigned int DAClusterizerInZ_vect::thermalize(
    double beta, track_t& tks, vertex_t& v, const double delta_max0, const double rho0) const {
  unsigned int niter = 0;
  double delta = 0;
  double delta_max = delta_lowT_;

  if (convergence_mode_ == 0) {
    delta_max = delta_max0;
  } else if (convergence_mode_ == 1) {
    delta_max = delta_lowT_ / sqrt(std::max(beta, 1.0));
  }

  set_vtx_range(beta, tks, v);
  double delta_sum_range = 0;  // accumulate max(|delta-z|) as a lower bound
  std::vector<double> z0 = v.zvtx_vec;

  while (niter++ < maxIterations_) {
    delta = update(beta, tks, v, rho0);
    delta_sum_range += delta;

    if (delta_sum_range > zrange_min_) {
      for (unsigned int k = 0; k < v.getSize(); k++) {
        if (std::abs(v.zvtx_vec[k] - z0[k]) > zrange_min_) {
          set_vtx_range(beta, tks, v);
          delta_sum_range = 0;
          z0 = v.zvtx_vec;
          break;
        }
      }
    }

    if (delta < delta_max) {
      break;
    }
  }

#ifdef DEBUG
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect.thermalize niter = " << niter << " at T = " << 1 / beta
              << "  nv = " << v.getSize() << std::endl;
    if (DEBUGLEVEL > 2)
      dump(beta, v, tks, 0, rho0);
  }
#endif

  return niter;
}

bool DAClusterizerInZ_vect::merge(vertex_t& y, track_t& tks, double& beta) const {
  // merge clusters that collapsed or never separated,
  // only merge if the estimated critical temperature of the merged vertex is below the current temperature
  // return true if vertices were merged, false otherwise
  const unsigned int nv = y.getSize();

  if (nv < 2)
    return false;

  // merge the smallest distance clusters first
  std::vector<std::pair<double, unsigned int>> critical;
  for (unsigned int k = 0; (k + 1) < nv; k++) {
    if (std::fabs(y.zvtx[k + 1] - y.zvtx[k]) < zmerge_) {
      critical.push_back(make_pair(std::fabs(y.zvtx[k + 1] - y.zvtx[k]), k));
    }
  }
  if (critical.empty())
    return false;

  std::stable_sort(critical.begin(), critical.end(), std::less<std::pair<double, unsigned int>>());

  for (unsigned int ik = 0; ik < critical.size(); ik++) {
    unsigned int k = critical[ik].second;
    double rho = y.rho[k] + y.rho[k + 1];

#ifdef DEBUG
    assert((k + 1) < nv);
    if (DEBUGLEVEL > 1) {
      std::cout << "merging " << fixed << setprecision(4) << y.zvtx[k + 1] << " and " << y.zvtx[k]
                << "  sw = " << y.sw[k] + y.sw[k + 1] << std::endl;
    }
#endif

    if (rho > 0) {
      y.zvtx[k] = (y.rho[k] * y.zvtx[k] + y.rho[k + 1] * y.zvtx[k + 1]) / rho;
      if (not edm::isFinite(y.zvtx[k]))
        y.zvtx[k] = 0.5 * (y.zvtx[k] + y.zvtx[k + 1]);
    } else {
      y.zvtx[k] = 0.5 * (y.zvtx[k] + y.zvtx[k + 1]);
    }
    y.rho[k] = rho;
    y.sw[k] += y.sw[k + 1];
    y.removeItem(k + 1, tks);
    set_vtx_range(beta, tks, y);
    y.extractRaw();
    return true;
  }

  return false;
}

bool DAClusterizerInZ_vect::purge(vertex_t& y, track_t& tks, double& rho0, const double beta) const {
  constexpr double eps = 1.e-100;
  // eliminate clusters with only one significant/unique track
  const unsigned int nv = y.getSize();
  const unsigned int nt = tks.getSize();

  if (nv < 2)
    return false;

  std::vector<double> sump_v(nv), arg_cache_v(nv), exp_cache_v(nv), pcut_cache_v(nv);
  std::vector<int> nUnique_v(nv);
  double* __restrict__ parg_cache;
  double* __restrict__ pexp_cache;
  double* __restrict__ ppcut_cache;
  double* __restrict__ psump;
  int* __restrict__ pnUnique;
  int constexpr nunique_min_ = 2;

  set_vtx_range(beta, tks, y);

  parg_cache = arg_cache_v.data();
  pexp_cache = exp_cache_v.data();
  ppcut_cache = pcut_cache_v.data();
  psump = sump_v.data();
  pnUnique = nUnique_v.data();

  const auto rhoconst = rho0 * local_exp(-beta * dzCutOff_ * dzCutOff_);
  for (unsigned int k = 0; k < nv; k++) {
    const double pmax = y.rho[k] / (y.rho[k] + rhoconst);
    ppcut_cache[k] = uniquetrkweight_ * pmax;
  }

  for (unsigned int i = 0; i < nt; i++) {
    const auto invZ = ((tks.sum_Z[i] > eps) && (tks.tkwt[i] > uniquetrkminp_)) ? 1. / tks.sum_Z[i] : 0.;
    const auto track_z = tks.zpca[i];
    const auto botrack_dz2 = -beta * tks.dz2[i];
    const auto kmin = tks.kmin[i];
    const auto kmax = tks.kmax[i];

    for (unsigned int k = kmin; k < kmax; k++) {
      const auto mult_resz = track_z - y.zvtx[k];
      parg_cache[k] = botrack_dz2 * (mult_resz * mult_resz);
    }

    local_exp_list_range(parg_cache, pexp_cache, kmin, kmax);

    for (unsigned int k = kmin; k < kmax; k++) {
      const double p = y.rho[k] * pexp_cache[k] * invZ;
      psump[k] += p;
      pnUnique[k] += (p > ppcut_cache[k]) ? 1 : 0;
    }
  }

  double sumpmin = nt;
  unsigned int k0 = nv;
  for (unsigned k = 0; k < nv; k++) {
    if ((pnUnique[k] < nunique_min_) && (psump[k] < sumpmin)) {
      sumpmin = psump[k];
      k0 = k;
    }
  }

  if (k0 != nv) {
#ifdef DEBUG
    assert(k0 < y.getSize());
    if (DEBUGLEVEL > 1) {
      std::cout << "eliminating prototype at " << std::setw(10) << std::setprecision(4) << y.zvtx[k0]
                << " with sump=" << sumpmin << "  rho*nt =" << y.rho[k0] * nt << " pnUnique=" << pnUnique[k0] << endl;
    }
#endif

    y.removeItem(k0, tks);
    set_vtx_range(beta, tks, y);
    return true;
  } else {
    return false;
  }
}

double DAClusterizerInZ_vect::beta0(double betamax, track_t const& tks, vertex_t const& y) const {
  double T0 = 0;  // max Tc for beta=0
  // estimate critical temperature from beta=0 (T=inf)
  const unsigned int nt = tks.getSize();
  const unsigned int nv = y.getSize();

  for (unsigned int k = 0; k < nv; k++) {
    // vertex fit at T=inf
    double sumwz = 0;
    double sumw = 0;
    for (unsigned int i = 0; i < nt; i++) {
      double w = tks.tkwt[i] * tks.dz2[i];
      sumwz += w * tks.zpca[i];
      sumw += w;
    }

    y.zvtx[k] = sumwz / sumw;

    // estimate Tcrit
    double a = 0, b = 0;
    for (unsigned int i = 0; i < nt; i++) {
      double dx = tks.zpca[i] - y.zvtx[k];
      double w = tks.tkwt[i] * tks.dz2[i];
      a += w * std::pow(dx, 2) * tks.dz2[i];
      b += w;
    }
    double Tc = 2. * a / b;  // the critical temperature of this vertex

    if (Tc > T0)
      T0 = Tc;

  }  // vertex loop (normally there should be only one vertex at beta=0)

#ifdef DEBUG
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect.beta0:   Tc = " << T0 << std::endl;
    int coolingsteps = 1 - int(std::log(T0 * betamax) / std::log(coolingFactor_));
    std::cout << "DAClusterizerInZ_vect.beta0:   nstep = " << coolingsteps << std::endl;
  }
#endif

  if (T0 > 1. / betamax) {
    int coolingsteps = 1 - int(std::log(T0 * betamax) / std::log(coolingFactor_));

    return betamax * std::pow(coolingFactor_, coolingsteps);
  } else {
    // ensure at least one annealing step
    return betamax * coolingFactor_;
  }
}

bool DAClusterizerInZ_vect::split(const double beta, track_t& tks, vertex_t& y, double threshold) const {
  // split only critical vertices (Tc >~ T=1/beta   <==>   beta*Tc>~1)
  // returns true if at least one cluster was split

  // update the sums needed for Tc, rho0 is never non-zero while splitting is still active
  update(beta, tks, y, 0., true);  // the "true" flag enables Tc evaluation

  constexpr double epsilon = 1e-3;  // minimum split size
  unsigned int nv = y.getSize();

  // avoid left-right biases by splitting highest Tc first; Handwavey: Chi-squared to bad then split vertex

  std::vector<std::pair<double, unsigned int>> critical;
  for (unsigned int k = 0; k < nv; k++) {
    double Tc = 2 * y.swE[k] / y.sw[k];
    if (edm::isFinite(Tc) and beta * Tc > threshold) {
      critical.push_back(make_pair(Tc, k));
    }
  }
  if (critical.empty())
    return false;

  std::stable_sort(critical.begin(), critical.end(), std::greater<std::pair<double, unsigned int>>());

  bool split = false;
  const unsigned int nt = tks.getSize();

  for (unsigned int ic = 0; ic < critical.size(); ic++) {
    unsigned int k = critical[ic].second;

    // estimate subcluster positions and weight
    double p1 = 0, z1 = 0, w1 = 0;
    double p2 = 0, z2 = 0, w2 = 0;
    for (unsigned int i = 0; i < nt; i++) {
      if (tks.sum_Z[i] > 1.e-100) {
        // winner-takes-all, usually overestimates splitting; What does this mean??
        double tl = tks.zpca[i] < y.zvtx[k] ? 1. : 0.;
        double tr = 1. - tl;

        // soften it, especially at low T
        double arg = (tks.zpca[i] - y.zvtx[k]) * sqrt(beta * tks.dz2[i]);
        if (std::fabs(arg) < 20) {
          double t = local_exp(-arg);
          tl = t / (t + 1.);
          tr = 1 / (t + 1.);
        }
        // calculate distance between the 2 new vertices, this is not constant; handwavey not clear if it helps
        double p = y.rho[k] * tks.tkwt[i] * local_exp(-beta * Eik(tks.zpca[i], y.zvtx[k], tks.dz2[i])) / tks.sum_Z[i];
        double w = p * tks.dz2[i];
        p1 += p * tl;
        z1 += w * tl * tks.zpca[i];
        w1 += w * tl;
        p2 += p * tr;
        z2 += w * tr * tks.zpca[i];
        w2 += w * tr;
      }
    }

    if (w1 > 0) {
      z1 = z1 / w1;
    } else {
      z1 = y.zvtx[k] - epsilon;
    }
    if (w2 > 0) {
      z2 = z2 / w2;
    } else {
      z2 = y.zvtx[k] + epsilon;
    }

    // reduce split size if there is not enough room
    if ((k > 0) && (z1 < (0.6 * y.zvtx[k] + 0.4 * y.zvtx[k - 1]))) {
      z1 = 0.6 * y.zvtx[k] + 0.4 * y.zvtx[k - 1];
    }
    if ((k + 1 < nv) && (z2 > (0.6 * y.zvtx[k] + 0.4 * y.zvtx[k + 1]))) {
      z2 = 0.6 * y.zvtx[k] + 0.4 * y.zvtx[k + 1];
    }

#ifdef DEBUG
    assert(k < nv);
    if (DEBUGLEVEL > 1) {
      if (std::fabs(y.zvtx[k] - zdumpcenter_) < zdumpwidth_) {
        std::cout << " T= " << std::setw(8) << 1. / beta << " Tc= " << critical[ic].first << "    splitting "
                  << std::fixed << std::setprecision(4) << y.zvtx[k] << " --> " << z1 << "," << z2 << "     [" << p1
                  << "," << p2 << "]";
        if (std::fabs(z2 - z1) > epsilon) {
          std::cout << std::endl;
        } else {
          std::cout << "  rejected " << std::endl;
        }
      }
    }
#endif

    // split if the new subclusters are significantly separated
    if ((z2 - z1) > epsilon) {
      split = true;
      double pk1 = p1 * y.rho[k] / (p1 + p2);
      double pk2 = p2 * y.rho[k] / (p1 + p2);

      if (not(edm::isFinite(pk1) and edm::isFinite(pk2)))
        continue;

      y.zvtx[k] = z2;
      y.rho[k] = pk2;
      y.insertItem(k, z1, pk1, tks);
      if (k == 0)
        y.extractRaw();

      nv++;

      // adjust remaining pointers
      for (unsigned int jc = ic; jc < critical.size(); jc++) {
        if (critical[jc].second >= k) {
          critical[jc].second++;
        }
      }
    } else {
#ifdef DEBUG
      std::cout << "warning ! split rejected, too small." << endl;
#endif
    }
  }

  return split;
}

// Standard method DA pure
vector<TransientVertex> DAClusterizerInZ_vect::vertices_no_blocks(const vector<reco::TransientTrack>& tracks) const {

  //initalize csv file
//const std::string directory = "/work/msaxer/ba/data/";
const std::string directory = "./"; // Current directory for slurm jobs

const std::string base_filename = "daten_";
auto now = std::chrono::system_clock::now();
std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
std::tm* now_tm = std::localtime(&now_time_t);
char buffer[80];
std::strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", now_tm);
std::string filename = directory + base_filename + buffer + ".csv";
std::ofstream daten_csv(filename, std::ios::app);
if (!daten_csv.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
}
std::ostringstream oss;
oss << "Current date and time: " << std::asctime(now_tm) << "running the vertices in  no blocks method" << std::endl;
oss << "Checkpoint;Time it took (microseconds); Number of clusters after checkpoint; comment" << std::endl;







  track_t&& tks = fill(tracks);
  vector<TransientVertex> clusters;
  if (tks.getSize() == 0)
    return clusters;
  tks.extractRaw();

  double rho0 = 0.0;  // start with no outlier rejection

  vertex_t y;  // the vertex prototypes

  // initialize:single vertex at infinite temperature
  y.addItem(0, 1.0); // First vertex is set here at z=0
  clear_vtx_range(tks, y);

  // estimate first critical temperature; betamax = Tmin, beta = 1/T
  double beta = beta0(betamax_, tks, y);
#ifdef DEBUG
  if (DEBUGLEVEL > 0)
    std::cout << "Beta0 is " << beta << std::endl;
#endif

  thermalize(beta, tks, y, delta_highT_); // Iteration at T=const, takes for-f-ever for high, till stable for clusters

  // annealing loop, stop when T<Tmin  (i.e. beta>1/Tmin)
  auto start_clustering_second_loop = std::chrono::high_resolution_clock::now();

  double betafreeze = betamax_ * sqrt(coolingFactor_);

  // main loop which takes a long time for high T; this runs until stable
  while (beta < betafreeze) {
    while (merge(y, tks, beta)) {
      update(beta, tks, y, rho0, false);
    }
    split(beta, tks, y);

    beta = beta / coolingFactor_;
    thermalize(beta, tks, y, delta_highT_);
          cout << "beta is" << beta << std::endl;
      cout << "betafreeze is" << betafreeze << std::endl;

  }
      auto stop_clustering_second_loop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> second_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(stop_clustering_second_loop - start_clustering_second_loop);
std::cout<<"the the second loop clustering took ms :"<< second_loop_clustering.count() << std::endl;
cout << "made it trough the second loop" << std::endl;
cout << "size after second loop" << y.getSize() << std::endl;
oss << "global Annealing;" << second_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;

#ifdef DEBUG
  verify(y, tks);

  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices_no_blocks :"
              << "last round of splitting" << std::endl;
  }
#endif

  set_vtx_range(beta, tks, y);
  update(beta, tks, y, rho0, false);

  while (merge(y, tks, beta)) {
    set_vtx_range(beta, tks, y);
    update(beta, tks, y, rho0, false);
  }
    auto further_cooling_timer_start = std::chrono::high_resolution_clock::now();


  // Right place?; Further cooling post spiltting to get closer to T=1, closer to Gauss-dist
  unsigned int ntry = 0;
  double threshold = 1.0;
  while (split(beta, tks, y, threshold) && (ntry++ < 10)) {
    thermalize(beta, tks, y, delta_highT_, rho0);  // rho0 = 0. here
    while (merge(y, tks, beta)) {
      update(beta, tks, y, rho0, false);
    }

    // relax splitting a bit to reduce multiple split-merge cycles of the same cluster
    threshold *= 1.1;
  }
        auto further_cooling_timer_stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> further_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(further_cooling_timer_stop - further_cooling_timer_start);
std::cout<<"further cooling took ms :"<< further_loop_clustering.count() << std::endl;
cout << "size after further cooling loop" << y.getSize() << std::endl;
oss << "further cooling;" << further_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;

#ifdef DEBUG
  verify(y, tks);
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices_no_blocks :"
              << "turning on outlier rejection at T=" << 1 / beta << std::endl;
  }
#endif

  // switch on outlier rejection at T=Tmin
  if (dzCutOff_ > 0) {
    rho0 = y.getSize() > 1 ? 1. / y.getSize() : 1.;
    for (unsigned int a = 0; a < 5; a++) {
      update(beta, tks, y, a * rho0 / 5.);  // adiabatic turn-on
    }
  }

  thermalize(beta, tks, y, delta_lowT_, rho0);
    auto collapsing_outliers_timer_start = std::chrono::high_resolution_clock::now();

#ifdef DEBUG
  verify(y, tks);
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices_no_blocks :"
              << "merging with outlier rejection at T=" << 1 / beta << std::endl;
  }
  if (DEBUGLEVEL > 2)
    dump(beta, y, tks, 2, rho0);
#endif

  // merge again (some cluster split by outliers collapse here)
  while (merge(y, tks, beta)) {
    set_vtx_range(beta, tks, y);
    update(beta, tks, y, rho0, false);
  }
    auto collapsing_outliers_timer_stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> collapsing_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(collapsing_outliers_timer_stop - collapsing_outliers_timer_start);
std::cout<<"collapsing outliers  took ms :"<< collapsing_loop_clustering.count() << std::endl;
cout << "size after collapsing outliers loop" << y.getSize() << std::endl;
oss << "collapsing outliers;" << collapsing_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;
    auto eliminate_insignificant_time_start = std::chrono::high_resolution_clock::now();
#ifdef DEBUG
  verify(y, tks);
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices_no_blocks :"
              << "after merging with outlier rejection at T=" << 1 / beta << std::endl;
  }
  if (DEBUGLEVEL > 2)
    dump(beta, y, tks, 2, rho0);
#endif

  // go down to the purging temperature (if it is lower than Tmin)
  while (beta < betapurge_) {
    beta = min(beta / coolingFactor_, betapurge_);
    thermalize(beta, tks, y, delta_lowT_, rho0);
  }

#ifdef DEBUG
  verify(y, tks);
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices :"
              << "purging at T=" << 1 / beta << std::endl;
  }
#endif

  // eliminate insignificant vertices, this is more restrictive at higher T
  while (purge(y, tks, rho0, beta)) {
    thermalize(beta, tks, y, delta_lowT_, rho0);
  }

#ifdef DEBUG
  verify(y, tks);
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices_no_blocks :"
              << "last cooling T=" << 1 / beta << std::endl;
  }
#endif
      auto eliminate_insignificant_time_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int, std::micro> eliminate_insignificant_duration = std::chrono::duration_cast<std::chrono::microseconds>(eliminate_insignificant_time_stop - eliminate_insignificant_time_start);
std::cout<<"eliminating insignificant (purging) took ms :"<< eliminate_insignificant_duration.count() << std::endl;
cout << "size after eliminating insignificant" << y.getSize() << std::endl;
oss << "eliminating insignificant (purge);" << eliminate_insignificant_duration.count() << ";" << y.getSize() << ";none" << std::endl;

    auto cool_some_more_time_start = std::chrono::high_resolution_clock::now();

  // optionally cool some more without doing anything, to make the track assignment harder (harder = sharper more clear)
  while (beta < betastop_) {
    beta = min(beta / coolingFactor_, betastop_);
    thermalize(beta, tks, y, delta_lowT_, rho0);
  }
    auto cool_some_more_time_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int, std::micro> cool_some_more_duration = std::chrono::duration_cast<std::chrono::microseconds>(cool_some_more_time_stop - cool_some_more_time_start);
std::cout<<"some more cooling  took ms :"<< cool_some_more_duration.count() << std::endl;
cout << "size after some more cooling" << y.getSize() << std::endl;
oss << "some more cooling;" << cool_some_more_duration.count() << ";" << y.getSize() << ";none" << std::endl;

#ifdef DEBUG
  verify(y, tks);
  if (DEBUGLEVEL > 0) {
    std::cout << "DAClusterizerInZ_vect::vertices_no_blocks :"
              << "stop cooling at T=" << 1 / beta << std::endl;
  }
  if (DEBUGLEVEL > 2)
    dump(beta, y, tks, 2, rho0);
#endif
daten_csv<< oss.str();

daten_csv.close();

  // assign tracks and fill into transient vertices
  return fill_vertices(beta, rho0, tks, y);
}

// Wolfram implementation of getting blocks, use this in vertices_in_blocks
std::vector<float> DAClusterizerInZ_vect::get_block_boundaries(const std::vector<reco::TransientTrack>& tracks) const {
  /* get the block boundaries as used by vertices_in_blocks
     the code is a direct copy from that method, but does not produce tracklists,
     it only stores the position of the first and the last track of a block
  */
  std::vector<float> values; // two numbers for each block: z of the first and the last track
  if (tracks.size() == 0) return values;  // don't want to handle this case
  
  std::vector<reco::TransientTrack> sorted_tracks; // Can we change direction here?
  for (unsigned int i = 0; i < tracks.size(); i++) {
    sorted_tracks.push_back(tracks[i]);
  }
  std::sort(sorted_tracks.begin(),
            sorted_tracks.end(),
            [](const reco::TransientTrack& a, const reco::TransientTrack& b) -> bool {
              return (a.stateAtBeamLine().trackStateAtPCA()).position().z() <
                     (b.stateAtBeamLine().trackStateAtPCA()).position().z();
            });

  unsigned int nBlocks = (unsigned int)std::floor(sorted_tracks.size() / (block_size_ * (1 - overlap_frac_)));
  if (nBlocks < 1) {
    nBlocks = 1;
  }
  for (unsigned int block = 0; block < nBlocks; block++) {
    unsigned int begin = (unsigned int)(block * block_size_ * (1 - overlap_frac_));
    unsigned int end = (unsigned int)std::min(begin + block_size_, (unsigned int)sorted_tracks.size());
    // instead of copying the tracks we just note the coordinate of the first and the last one here
    values.push_back((sorted_tracks[begin].stateAtBeamLine().trackStateAtPCA()).position().z());
    if (end > 0){
      values.push_back((sorted_tracks[end-1].stateAtBeamLine().trackStateAtPCA()).position().z());
    }else{
      values.push_back(values[0]);
    }
  }
  return values;
}


// DA in blocks
vector<TransientVertex> DAClusterizerInZ_vect::vertices_in_blocks(const vector<reco::TransientTrack>& tracks) const {
  //initalize csv file
//const std::string directory = "/work/msaxer/ba/data/";
const std::string directory = "./"; // Current directory for slurm jobs

const std::string base_filename = "daten_";
auto now = std::chrono::system_clock::now();
std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
std::tm* now_tm = std::localtime(&now_time_t);
char buffer[80];
std::strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", now_tm);
std::string filename = directory + base_filename + buffer + ".csv";
std::ofstream daten_csv(filename, std::ios::app);
if (!daten_csv.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
}
std::ostringstream oss;
oss << "Current date and time: " << std::asctime(now_tm) << "running the vertices in blocks method" << std::endl;
oss << "Checkpoint;Time it took (microseconds); Number of clusters after checkpoint; comment" << std::endl;


//csv end


  cout<<"running in blocks"<<std::endl;


  vector<reco::TransientTrack> sorted_tracks; // initalizes empty vectors and coppies all tracks into it
  vector<pair<float, float>> vertices_tot;  // z, rho for each vertex

// using this vector we collect all vertices protoypes form 
vertex_t combined_vertex_prototypes;
double betasave = 0.0;
for (unsigned int i = 0; i < tracks.size(); i++)
{
  sorted_tracks.push_back(tracks[i]);
  }
  double rho0, beta; // get blocborders
  auto blockBoundaries = get_block_boundaries(sorted_tracks);  
  // starting timer for clustering in blocks
    auto start_clustering_first_loop = std::chrono::high_resolution_clock::now();
  /* Reneval of code here. THis now gets block boundaries and then works for DA in blocks*/
  // iterates over each block defined in blocboundaries for each block it finds the range of tracks that fall into it and collects those tracks
  for (unsigned int b = 0; b < blockBoundaries.size(); b += 2) {
    float zBegin = blockBoundaries[b];
    float zEnd   = blockBoundaries[b+1];
    
    // find iBegin as the first track with z >= zBegin
    // find iEnd   as the first track with z > zEnd
    auto itBegin = std::lower_bound(
        sorted_tracks.begin(), sorted_tracks.end(), zBegin,
        [](auto const& tk, float zVal) {
            return tk.stateAtBeamLine().trackStateAtPCA().position().z() < zVal;
        }
    );
    auto itEnd = std::upper_bound(
        sorted_tracks.begin(), sorted_tracks.end(), zEnd,
        [](float zVal, auto const& tk) {
            return zVal < tk.stateAtBeamLine().trackStateAtPCA().position().z();
        }
    );
    unsigned int beginIdx = itBegin - sorted_tracks.begin();
    unsigned int endIdx   = itEnd   - sorted_tracks.begin();

    // gather block tracks
    std::vector<reco::TransientTrack> block_tracks;
    block_tracks.reserve(endIdx - beginIdx);
  
    for (unsigned int i = beginIdx; i < endIdx; i++) {
      block_tracks.push_back(sorted_tracks[i]); // 
      block_tracks.push_back(sorted_tracks[i]); // vielleicht mal weglassen? // apparently this is used to increase weight eg its a weighting method and should be left in place
    }
// maybe further discuss this above conclusion wether or not this makes sense???
    if (block_tracks.empty()) {
      continue;
    }


#ifdef DEBUG
    std::cout << "Running vertices_in_blocks on" << std::endl;
    std::cout << "- block no." << block << " on " << nBlocks << " blocks " << std::endl;
    std::cout << "- block track size: " << sorted_tracks.size() << " - block size: " << block_size_ << std::endl;
#endif
    track_t&& tks = fill(block_tracks);
    tks.extractRaw();

    rho0 = 0.0;  // start with no outlier rejection

    vertex_t y;  // the vertex prototypes

    // initialize:single vertex at infinite temperature
    y.addItem(0, 1.0);
    clear_vtx_range(tks, y);

    // estimate first critical temperature
    beta = beta0(betamax_, tks, y);
      cout << "beta before 1 loop" << beta<< std::endl;

#ifdef DEBUG
    if (DEBUGLEVEL > 0)
      std::cout << "Beta0 is " << beta << std::endl;
#endif

    thermalize(beta, tks, y, delta_highT_);

    // annealing loop, stop when T<Tmin  (i.e. beta>1/Tmin)

    double betafreeze = 1.e-5; // 0.5; // seting betafreeze to T=20 betamax_ * sqrt(coolingFactor_);
    int iterations = 0;


    while (beta < betafreeze)
    {
      iterations++;
      while (merge(y, tks, beta))
      {
        update(beta, tks, y, rho0, false);
      }
      split(beta, tks, y);
     // cout << "iteration is " << iterations << std::endl;
    //  cout << "beta is" << beta << std::endl;
    //  cout << "betafreeze is" << betafreeze << std::endl;

      beta = beta / coolingFactor_;
      thermalize(beta, tks, y, delta_highT_);
    }
    // store vertex prototypes of the processed block
        // Add the refined vertex prototypes for this block to the combined vertex prototype
    for (unsigned int i = 0; i < y.getSize(); ++i) {
        combined_vertex_prototypes.addItem(y.zvtx_vec[i], y.rho_vec[i]);
    }

    betasave = beta;
    // what if we thermalize ombined_vertex_prototypes now?

    std::cout << "and the following niter" << iterations << std::endl;
    // closes  loop starting 1005
  }
  auto stop_clustering_first_loop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> first_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(stop_clustering_first_loop - start_clustering_first_loop);
std::cout<<"the first loop clustering took ms:"<< first_loop_clustering.count() << std::endl;


// Output the combined vertex prototype's cluster positions
//std::cout << "Combined Vertex Prototype Cluster Positions:" << std::endl;
//for (unsigned int i = 0; i < combined_vertex_prototypes.getSize(); ++i) {
 //   std::cout << "Cluster " << i << ": z = " << combined_vertex_prototypes.zvtx_vec[i]
 //             << ", rho = " << combined_vertex_prototypes.rho_vec[i] << std::endl;
//}

// (re)defining variables to fit to classic da
vertex_t y;  // the vertex prototypes

y = combined_vertex_prototypes;
cout << "size before" << y.getSize() << std::endl;
oss << "Annealing_in_blocks;" << first_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;


vector<TransientVertex> clusters;
 track_t tks = fill(sorted_tracks); // track_t&& doesent work it then says out of scope for while merge loop

 if (tks.getSize() == 0){
    return clusters;
 }





#ifdef cputime
  auto stop_clustering = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int, std::micro> tcpu_clustering = std::chrono::duration_cast<std::chrono::microseconds>(stop_clustering - start_clustering);
#endif


/* thermalization not necessary at such high temperatures
//trying to thermalize
  auto thermalizing_inbetween_loop_start = std::chrono::high_resolution_clock::now();

    thermalize(beta, tks, y, delta_highT_);
     while (merge(y, tks, beta))
    {
      update(beta, tks, y, rho0, false); // check update method again wether merging done correctly
    } // merging after thermalizing
  auto thermalizing_inbetween_loop_stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> thermalize_inbetween_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(thermalizing_inbetween_loop_stop - thermalizing_inbetween_loop_start);
std::cout<<"thermalizing innbetween took ms:"<< thermalize_inbetween_loop_clustering.count() << std::endl;
oss << "thermalizing_between_loops;" << thermalize_inbetween_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;

*/ 
// insert da global code here

// global annealing loop, stop when T<Tmin  (i.e. beta>1/Tmin)
//hardcoding beta to 0.005 not sure if this is necessary but maybe?
  //

  beta = betasave; //*0.5;
  cout << "beta before 2 loop" << beta << std::endl;

  // beta = 0.005;
  double betafreeze =  betamax_ * sqrt(coolingFactor_);/// betamax_ * sqrt(coolingFactor_);
  cout << "betafreeze second loop" << betafreeze << std::endl;
  cout << "made it to the second loop" << std::endl;
  auto start_clustering_second_loop = std::chrono::high_resolution_clock::now();

  // main loop which takes a long time for high T; this runs until stable
  while (beta < betafreeze)
  {
    while (merge(y, tks, beta))
    {
      update(beta, tks, y, rho0, false);
    }
    split(beta, tks, y);
    cout << "beta is" << beta << std::endl;
    cout << "betafreeze is" << betafreeze << std::endl;

    beta = beta / coolingFactor_;
    thermalize(beta, tks, y, delta_highT_);
  }
    auto stop_clustering_second_loop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> second_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(stop_clustering_second_loop - start_clustering_second_loop);
std::cout<<"the the second loop clustering took ms :"<< second_loop_clustering.count() << std::endl;
cout << "made it trough the second loop" << std::endl;
cout << "size after second loop" << y.getSize() << std::endl;
oss << "global Annealing;" << second_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;
#ifdef DEBUG
    verify(y, tks);

    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "last round of splitting" << std::endl;
    }
#endif



//trying to thermalize
  auto thermalizing_after_loop_start = std::chrono::high_resolution_clock::now();

    thermalize(beta, tks, y, delta_highT_);
     while (merge(y, tks, beta))
    {
      update(beta, tks, y, rho0, false); // check update method again wether merging done correctly
    } // merging after thermalizing
  auto thermalizing_after_loop_stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> thermalize_after_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(thermalizing_after_loop_stop - thermalizing_after_loop_start);
std::cout<<"thermalizing after loops took ms:"<< thermalize_after_loop_clustering.count() << std::endl;
oss << "thermalizing_after_loops;" << thermalize_after_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;

  set_vtx_range(beta, tks, y);
  update(beta, tks, y, rho0, false);

  while (merge(y, tks, beta)) {
    set_vtx_range(beta, tks, y);
    update(beta, tks, y, rho0, false);
  }
    auto further_cooling_timer_start = std::chrono::high_resolution_clock::now();


  // Right place?; Further cooling post spiltting to get closer to T=1, closer to Gauss-dist
  unsigned int ntry = 0;
  double threshold = 1.0;
  while (split(beta, tks, y, threshold) && (ntry++ < 10)) {
    thermalize(beta, tks, y, delta_highT_, rho0);  // rho0 = 0. here
    while (merge(y, tks, beta)) {
      update(beta, tks, y, rho0, false);
    }

    // relax splitting a bit to reduce multiple split-merge cycles of the same cluster
    threshold *= 1.1;
  }
      auto further_cooling_timer_stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> further_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(further_cooling_timer_stop - further_cooling_timer_start);
std::cout<<"further cooling took ms :"<< further_loop_clustering.count() << std::endl;
cout << "size after further cooling loop" << y.getSize() << std::endl;
oss << "further cooling;" << further_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;
    auto final_cooling_start = std::chrono::high_resolution_clock::now();

  // switch on outlier rejection at T=Tmin
  if (dzCutOff_ > 0) {
    rho0 = y.getSize() > 1 ? 1. / y.getSize() : 1.;
    for (unsigned int a = 0; a < 5; a++) {
      update(beta, tks, y, a * rho0 / 5.);  // adiabatic turn-on
    }
  }

  thermalize(beta, tks, y, delta_lowT_, rho0);

    auto collapsing_outliers_timer_start = std::chrono::high_resolution_clock::now();

  // merge again (some cluster split by outliers collapse here)
  while (merge(y, tks, beta)) {
    set_vtx_range(beta, tks, y);
    update(beta, tks, y, rho0, false);
  }
    auto collapsing_outliers_timer_stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<int, std::micro> collapsing_loop_clustering = std::chrono::duration_cast<std::chrono::microseconds>(collapsing_outliers_timer_stop - collapsing_outliers_timer_start);
std::cout<<"collapsing outliers  took ms :"<< collapsing_loop_clustering.count() << std::endl;
cout << "size after collapsing outliers loop" << y.getSize() << std::endl;
oss << "collapsing outliers;" << collapsing_loop_clustering.count() << ";" << y.getSize() << ";none" << std::endl;
    auto eliminate_insignificant_time_start = std::chrono::high_resolution_clock::now();
  // go down to the purging temperature (if it is lower than Tmin)
  while (beta < betapurge_) {
    beta = min(beta / coolingFactor_, betapurge_);
    thermalize(beta, tks, y, delta_lowT_, rho0);
  }






  // eliminate insignificant vertices, this is more restrictive at higher T
  while (purge(y, tks, rho0, beta)) {
    thermalize(beta, tks, y, delta_lowT_, rho0);
  }
      auto eliminate_insignificant_time_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int, std::micro> eliminate_insignificant_duration = std::chrono::duration_cast<std::chrono::microseconds>(eliminate_insignificant_time_stop - eliminate_insignificant_time_start);
std::cout<<"eliminating insignificant (purging) took ms :"<< eliminate_insignificant_duration.count() << std::endl;
cout << "size after eliminating insignificant" << y.getSize() << std::endl;
oss << "eliminating insignificant (purge);" << eliminate_insignificant_duration.count() << ";" << y.getSize() << ";none" << std::endl;

    auto cool_some_more_time_start = std::chrono::high_resolution_clock::now();

  // optionally cool some more without doing anything, to make the track assignment harder (harder = sharper more clear)
  while (beta < betastop_) {
    beta = min(beta / coolingFactor_, betastop_);
    thermalize(beta, tks, y, delta_lowT_, rho0);
  }
    auto cool_some_more_time_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<int, std::micro> cool_some_more_duration = std::chrono::duration_cast<std::chrono::microseconds>(cool_some_more_time_stop - cool_some_more_time_start);
std::cout<<"some more cooling  took ms :"<< cool_some_more_duration.count() << std::endl;
cout << "size after some more cooling" << y.getSize() << std::endl;
oss << "some more cooling;" << cool_some_more_duration.count() << ";" << y.getSize() << ";none" << std::endl;


    auto final_cooling_stop = std::chrono::high_resolution_clock::now();


  std::chrono::duration<int, std::micro> final_cooling_clustering = std::chrono::duration_cast<std::chrono::microseconds>(final_cooling_stop - final_cooling_start);
std::cout<<"total time cooling took ms :"<< final_cooling_clustering.count() << std::endl;
cout << "size at the end" << y.getSize() << std::endl;

daten_csv<< oss.str();

daten_csv.close();
cout << "wrote csv" << y.getSize() << std::endl;

  return fill_vertices(beta, rho0, tks, y);

/*
    set_vtx_range(beta, tks, y);
    update(beta, tks, y, rho0, false);

    while (merge(y, tks, beta)) {
      set_vtx_range(beta, tks, y);
      update(beta, tks, y, rho0, false);
    }

    unsigned int ntry = 0;
    double threshold = 1.0;
    while (split(beta, tks, y, threshold) && (ntry++ < 10)) {
      thermalize(beta, tks, y, delta_highT_, rho0);  // rho0 = 0. here
      while (merge(y, tks, beta)) {
        update(beta, tks, y, rho0, false);
      }

      // relax splitting a bit to reduce multiple split-merge cycles of the same cluster
      threshold *= 1.1;
    }

#ifdef DEBUG
    verify(y, tks);
    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "turning on outlier rejection at T=" << 1 / beta << std::endl;
    }
#endif

    // switch on outlier rejection at T=Tmin
    if (dzCutOff_ > 0) {
      rho0 = y.getSize() > 1 ? 1. / y.getSize() : 1.;
      for (unsigned int a = 0; a < 5; a++) {
        update(beta, tks, y, a * rho0 / 5.);  // adiabatic turn-on
      }
    }

    thermalize(beta, tks, y, delta_lowT_, rho0);

#ifdef DEBUG
    verify(y, tks);
    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "merging with outlier rejection at T=" << 1 / beta << std::endl;
    }
    if (DEBUGLEVEL > 2)
      dump(beta, y, tks, 2, rho0);
#endif

    // merge again  (some cluster split by outliers collapse here)
    while (merge(y, tks, beta)) {
      set_vtx_range(beta, tks, y);
      update(beta, tks, y, rho0, false);
    }

#ifdef DEBUG
    verify(y, tks);
    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "after merging with outlier rejection at T=" << 1 / beta << std::endl;
    }
    if (DEBUGLEVEL > 2)
      dump(beta, y, tks, 2, rho0);
#endif

    // go down to the purging temperature (if it is lower than tmin)
    while (beta < betapurge_) {
      beta = min(beta / coolingFactor_, betapurge_);
      thermalize(beta, tks, y, delta_lowT_, rho0);
    }

#ifdef DEBUG
    verify(y, tks);
    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "purging at T=" << 1 / beta << std::endl;
    }
#endif

    // eliminate insignificant vertices, this is more restrictive at higher T
    while (purge(y, tks, rho0, beta)) {
      thermalize(beta, tks, y, delta_lowT_, rho0);
    }

#ifdef DEBUG
    verify(y, tks);
    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "last cooling T=" << 1 / beta << std::endl;
    }
#endif

    // optionally cool some more without doing anything, to make the assignment harder
    while (beta < betastop_) {
      beta = min(beta / coolingFactor_, betastop_);
      thermalize(beta, tks, y, delta_lowT_, rho0);
    }

#ifdef DEBUG
    verify(y, tks);
    if (DEBUGLEVEL > 0) {
      std::cout << "DAClusterizerInZSubCluster_vect::vertices :"
                << "stop cooling at T=" << 1 / beta << std::endl;
    }
    if (DEBUGLEVEL > 2)
      dump(beta, y, tks, 2, rho0);
#endif

    // simple attempt to merge some duplicates in the overlap regions
    // a better or additional approach might be to do this during the track assignment step based
    if(unblock_){
      for (unsigned int ivertex = 0; ivertex < y.getSize(); ivertex++) {
	if (y.zvtx_vec[ivertex] != 0 && y.rho_vec[ivertex] != 0) {
	  bool merged = false;
	  for(auto & vtx : vertices_tot){
	    if (std::abs(vtx.first - y.zvtx_vec[ivertex]) < 20e-4){ // 20 microns 
	      double rho_new =  vtx.second  + y.rho_vec[ivertex];
	      double znew = (vtx.first * vtx.second + y.zvtx_vec[ivertex] * y.rho_vec[ivertex] ) / rho_new;
	      vtx.first = znew;
	      vtx.second = rho_new;
	      merged = true;
	    }
	  }
	  if (! merged){
	    vertices_tot.push_back(pair(y.zvtx_vec[ivertex], y.rho_vec[ivertex]));
	  }
	}
      }
    }else{
      // default code
      for (unsigned int ivertex = 0; ivertex < y.getSize(); ivertex++) {
	if (y.zvtx_vec[ivertex] != 0 && y.rho_vec[ivertex] != 0) {
	  vertices_tot.push_back(pair(y.zvtx_vec[ivertex], y.rho_vec[ivertex]));
#ifdef DEBUG
	  std::cout << "Found new vertex " << y.zvtx_vec[ivertex] << " , " << y.rho_vec[ivertex] << std::endl;
#endif
	}
      }
    }
  }

  std::sort(vertices_tot.begin(),
            vertices_tot.end(),
            [](const pair<float, float>& a, const pair<float, float>& b) -> bool { return a.first < b.first; });
  

  
  // reassign tracks to vertices
  track_t&& tracks_tot = fill(tracks);
  const unsigned int nv = vertices_tot.size();
  const unsigned int nt = tracks_tot.getSize();

  for (auto itrack = 0U; itrack < nt; ++itrack) {
    double zrange = max(sel_zrange_ / sqrt(beta * tracks_tot.dz2[itrack]), zrange_min_);

    double zmin = tracks_tot.zpca[itrack] - zrange;
    unsigned int kmin = min(nv - 1, tracks_tot.kmin[itrack]);
    // find the smallest vertex_z that is larger than zmin
    if (vertices_tot[kmin].first > zmin) {
      while ((kmin > 0) && (vertices_tot[kmin - 1].first > zmin)) {
        kmin--;
      }
    } else {
      while ((kmin < (nv - 1)) && (vertices_tot[kmin].first < zmin)) {
        kmin++;
      }
    }

    double zmax = tracks_tot.zpca[itrack] + zrange;
    unsigned int kmax = min(nv - 1, tracks_tot.kmax[itrack] - 1);
    // note: kmax points to the last vertex in the range, while gtracks.kmax points to the entry BEHIND the last vertex
    // find the largest vertex_z that is smaller than zmax
    if (vertices_tot[kmax].first < zmax) {
      while ((kmax < (nv - 1)) && (vertices_tot[kmax + 1].first < zmax)) {
        kmax++;
      }
    } else {
      while ((kmax > 0) && (vertices_tot[kmax].first > zmax)) {
        kmax--;
      }
    }

    if (kmin <= kmax) {
      tracks_tot.kmin[itrack] = kmin;
      tracks_tot.kmax[itrack] = kmax + 1;
    } else {
      tracks_tot.kmin[itrack] = max(0U, min(kmin, kmax));
      tracks_tot.kmax[itrack] = min(nv, max(kmin, kmax) + 1);
    }
  }

  rho0 = nv > 1 ? 1. / nv : 1.;
  const auto z_sum_init = rho0 * local_exp(-beta * dzCutOff_ * dzCutOff_);

  std::vector<std::vector<unsigned int>> vtx_track_indices(nv);
  for (unsigned int i = 0; i < nt; i++) {
    const auto kmin = tracks_tot.kmin[i];
    const auto kmax = tracks_tot.kmax[i];
    double p_max = -1;
    unsigned int iMax = 10000;  //just a "big" number w.r.t. number of vertices
    float sum_Z = z_sum_init;
    for (auto k = kmin; k < kmax; k++) {
      float v_exp = local_exp(-beta * Eik(tracks_tot.zpca[i], vertices_tot[k].first, tracks_tot.dz2[i]));
      sum_Z += vertices_tot[k].second * v_exp;
    }
    double invZ = sum_Z > 1e-100 ? 1. / sum_Z : 0.0;
    for (auto k = kmin; k < kmax && invZ != 0.0; k++) {
      float v_exp = local_exp(-beta * Eik(tracks_tot.zpca[i], vertices_tot[k].first, tracks_tot.dz2[i]));
      double p = vertices_tot[k].second * v_exp * invZ;
      if (p > p_max && p > mintrkweight_) {
        p_max = p;
        iMax = k;
      }
    }
    if (iMax < vtx_track_indices.size()) {
      vtx_track_indices[iMax].push_back(i);
    }
  }
#ifdef DEBUG
  for (auto itrack = 0U; itrack < nt; ++itrack) {
    std::cout << "itrack " << itrack << " , " << tracks_tot.kmin[itrack] << " , " << tracks_tot.kmax[itrack]
              << std::endl;
  }
#endif

  vector<TransientVertex> clusters;
  if (nv == 0) {
    return clusters;
  }

  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
  vector<reco::TransientTrack> vertexTracks;

  for (unsigned int k = 0; k < nv; k++) {
    if (!vtx_track_indices[k].empty()) {
      for (auto i : vtx_track_indices[k]) {
        vertexTracks.push_back(*(tracks_tot.tt[i]));
#ifdef DEBUG
        std::cout << vertices_tot[k].first << ","
                  << (*(tracks_tot.tt[i])).stateAtBeamLine().trackStateAtPCA().position().z() << std::endl;
#endif
      }
    }

    // implement what clusterize() did before : merge left-to-right if distance < 2 * vertexSize_
    if ((k + 1 == nv) || (abs(vertices_tot[k + 1].first - vertices_tot[k].first) > (2 * vertexSize_))) {
      // close a cluster
      if (vertexTracks.size() > 1) {
        GlobalPoint pos(0, 0, vertices_tot[k].first);  // only usable with subsequent fit
        TransientVertex v(pos, dummyError, vertexTracks, 0);
        clusters.push_back(v);
      }
      vertexTracks.clear();
    }
  }
  */



 // return clusters;

}  // end of vertices_in_blocks





vector<TransientVertex> DAClusterizerInZ_vect::fill_vertices(double beta, double rho0, track_t& tks, vertex_t& y) const {
  // select significant tracks and use a TransientVertex as a container

  set_vtx_range(beta, tks, y);
  const unsigned int nv = y.getSize();
  for (unsigned int k = 0; k < nv; k++) {
    if (edm::isNotFinite(y.rho[k]) || edm::isNotFinite(y.zvtx[k])) {
      y.rho[k] = 0;
      y.zvtx[k] = 0;
    }
  }

  // ensure consistent assignment probabillities and make a hard assignment
  const unsigned int nt = tks.getSize();
  const auto z_sum_init = rho0 * local_exp(-beta * dzCutOff_ * dzCutOff_);
  std::vector<std::vector<unsigned int>> vtx_track_indices(nv);
  std::vector<std::vector<float>> vtx_track_weights(nv);
  for (unsigned int i = 0; i < nt; i++) {
    const auto kmin = tks.kmin[i];
    const auto kmax = tks.kmax[i];
    for (auto k = kmin; k < kmax; k++) {
      y.exp_arg[k] = -beta * Eik(tks.zpca[i], y.zvtx[k], tks.dz2[i]);
    }

    local_exp_list_range(y.exp_arg, y.exp, kmin, kmax);

    tks.sum_Z[i] = z_sum_init;
    for (auto k = kmin; k < kmax; k++) {
      tks.sum_Z[i] += y.rho[k] * y.exp[k];
    }
    const double invZ = tks.sum_Z[i] > 1e-100 ? 1. / tks.sum_Z[i] : 0.0;

    double pmax = -1;
    unsigned int k_pmax = 0;
    for (auto k = kmin; k < kmax; k++) {
      double p = y.rho[k] * y.exp[k] * invZ;
      if (p > pmax) {
        pmax = p;
        k_pmax = k;
      }
    }

    if (pmax > mintrkweight_) {
      // assign to the cluster with the highest assignment weight, if it is at least mintrkweight_
      vtx_track_indices[k_pmax].push_back(i);
      vtx_track_weights[k_pmax].push_back(pmax);
    }
  }

  // fill transient vertices
  // the position is normally not used, probably not optimal when Tstop <> 2, anyway
  vector<TransientVertex> clusters;
  for (unsigned int k = 0; k < nv; k++) {
    double sump = 0;
    double sumw = 0;
    double sumwp = 0, sumwz = 0;
    if (!vtx_track_indices[k].empty()) {
      vector<reco::TransientTrack> vertexTracks;
      TransientVertex::TransientTrackToFloatMap trkWeightMap;
      unsigned int j = 0;
      for (auto i : vtx_track_indices[k]) {
        auto p = vtx_track_weights[k][j];
        vertexTracks.push_back(*(tks.tt[i]));
        trkWeightMap[vertexTracks[j]] = p;
        auto w = p * tks.dz2[i];
        sump += p;
        sumw += w;
        sumwp += w * p;
        sumwz += w * tks.zpca[i];
        j++;
      }
      float zerror_squared = 1.;  //
      if ((sumw > 0) && (sumwp > 0)) {
        zerror_squared = sumwp / (sumw * sumw);
        y.zvtx[k] = sumwz / sumw;
      }

      reco::BeamSpot bs = vertexTracks[0].stateAtBeamLine().beamSpot();
      GlobalPoint pos(bs.x(y.zvtx[k]), bs.y(y.zvtx[k]), y.zvtx[k]);
      const float xerror_squared = pow(bs.BeamWidthX(), 2);
      const float yerror_squared = pow(bs.BeamWidthY(), 2);
      GlobalError err(xerror_squared, 0, yerror_squared, 0., 0., zerror_squared);
      TransientVertex v(pos, err, vertexTracks, 0, 2 * sump - 3.);
      v.weightMap(trkWeightMap);
      clusters.push_back(v);
    }
  }

  return clusters;
}

vector<TransientVertex> DAClusterizerInZ_vect::vertices(const vector<reco::TransientTrack>& tracks) const {
  if (runInBlocks_ and (block_size_ < tracks.size()))  //doesn't bother if low number of tracks
    return vertices_in_blocks(tracks);
  else
    return vertices_no_blocks(tracks);
}

vector<vector<reco::TransientTrack>> DAClusterizerInZ_vect::clusterize(  // OBSOLETE
    const vector<reco::TransientTrack>& tracks) const {
  vector<vector<reco::TransientTrack>> clusters;
  vector<TransientVertex>&& pv = vertices(tracks);

#ifdef DEBUG
  if (DEBUGLEVEL > 0) {
    std::cout << "###################################################" << endl;
    std::cout << "# vectorized DAClusterizerInZ_vect::clusterize   nt=" << tracks.size() << endl;
    std::cout << "# DAClusterizerInZ_vect::clusterize   pv.size=" << pv.size() << endl;
    std::cout << "###################################################" << endl;
  }
#endif

  if (pv.empty()) {
    return clusters;
  }

  // fill into clusters and merge
  vector<reco::TransientTrack> aCluster = pv.begin()->originalTracks();

  for (auto k = pv.begin() + 1; k != pv.end(); k++) {
    if (std::abs(k->position().z() - (k - 1)->position().z()) > (2 * vertexSize_)) {
      // close a cluster
      if (aCluster.size() > 1) {
        clusters.push_back(aCluster);
      }
#ifdef DEBUG
      else {
        std::cout << " one track cluster at " << k->position().z() << "  suppressed" << std::endl;
      }
#endif
      aCluster.clear();
    }
    for (unsigned int i = 0; i < k->originalTracks().size(); i++) {
      aCluster.push_back(k->originalTracks()[i]);
    }
  }
  clusters.emplace_back(std::move(aCluster));

  return clusters;
}

void DAClusterizerInZ_vect::dump(
    const double beta, const vertex_t& y, const track_t& tks, const int verbosity, const double rho0) const {
#ifdef DEBUG
  const unsigned int nv = y.getSize();
  const unsigned int nt = tks.getSize();
  // make a local copies of clusters and tracks to update Tc  [ = y_local.swE / y_local.sw ]
  vertex_t y_local = y;
  track_t tks_local = tks;
  update(beta, tks_local, y_local, rho0, true);

  std::vector<unsigned int> iz;
  for (unsigned int j = 0; j < nt; j++) {
    iz.push_back(j);
  }
  std::sort(iz.begin(), iz.end(), [tks](unsigned int a, unsigned int b) { return tks.zpca[a] < tks.zpca[b]; });
  std::cout << std::endl;
  std::cout << "-----DAClusterizerInZ::dump ----" << nv << "  clusters " << std::endl;
  std::cout << "                                                                   ";
  for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
    if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) < zdumpwidth_) {
      std::cout << "   " << setw(3) << ivertex << "  ";
    }
  }
  std::cout << endl;
  std::cout << "                                                                z= ";
  std::cout << setprecision(4);
  for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
    if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) < zdumpwidth_) {
      std::cout << setw(8) << fixed << y.zvtx[ivertex];
    }
  }
  std::cout << endl
            << "T=" << setw(15) << 1. / beta << " Tmin =" << setw(10) << 1. / betamax_
            << "                             Tc= ";
  for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
    if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) < zdumpwidth_) {
      double Tc = 2 * y_local.swE[ivertex] / y_local.sw[ivertex];
      std::cout << setw(8) << fixed << setprecision(1) << Tc;
    }
  }
  std::cout << endl;

  std::cout << "                                                               pk= ";
  double sumpk = 0;
  for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
    sumpk += y.rho[ivertex];
    if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) > zdumpwidth_)
      continue;
    std::cout << setw(8) << setprecision(4) << fixed << y.rho[ivertex];
  }
  std::cout << endl;

  std::cout << "                                                               nt= ";
  for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
    if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) > zdumpwidth_)
      continue;
    std::cout << setw(8) << setprecision(1) << fixed << y.rho[ivertex] * nt;
  }
  std::cout << endl;

  if (verbosity > 0) {
    double E = 0, F = 0;
    std::cout << endl;
    std::cout << "----        z +/- dz                ip +/-dip       pt    phi  eta    weights  ----" << endl;
    std::cout << setprecision(4);
    for (unsigned int i0 = 0; i0 < nt; i0++) {
      unsigned int i = iz[i0];
      if (tks.sum_Z[i] > 0) {
        F -= std::log(tks.sum_Z[i]) / beta;
      }
      double tz = tks.zpca[i];

      if (std::fabs(tz - zdumpcenter_) > zdumpwidth_)
        continue;
      std::cout << setw(4) << i << ")" << setw(8) << fixed << setprecision(4) << tz << " +/-" << setw(6)
                << sqrt(1. / tks.dz2[i]);
      if ((tks.tt[i] == nullptr)) {
        std::cout << "          effective track                             ";
      } else {
        if (tks.tt[i]->track().quality(reco::TrackBase::highPurity)) {
          std::cout << " *";
        } else {
          std::cout << "  ";
        }
        if (tks.tt[i]->track().hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1)) {
          std::cout << "+";
        } else {
          std::cout << "-";
        }
        std::cout << setw(1)
                  << tks.tt[i]
                         ->track()
                         .hitPattern()
                         .pixelBarrelLayersWithMeasurement();  // see DataFormats/TrackReco/interface/HitPattern.h
        std::cout << setw(1) << tks.tt[i]->track().hitPattern().pixelEndcapLayersWithMeasurement();
        std::cout << setw(1) << hex
                  << tks.tt[i]->track().hitPattern().trackerLayersWithMeasurement() -
                         tks.tt[i]->track().hitPattern().pixelLayersWithMeasurement()
                  << dec;
        std::cout << "=" << setw(1) << hex
                  << tks.tt[i]->track().hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) << dec;

        Measurement1D IP = tks.tt[i]->stateAtBeamLine().transverseImpactParameter();
        std::cout << setw(8) << IP.value() << "+/-" << setw(6) << IP.error();
        std::cout << " " << setw(6) << setprecision(2) << tks.tt[i]->track().pt() * tks.tt[i]->track().charge();
        std::cout << " " << setw(5) << setprecision(2) << tks.tt[i]->track().phi() << " " << setw(5) << setprecision(2)
                  << tks.tt[i]->track().eta();
      }  // not a dummy track

      double sump = 0.;
      for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
        if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) > zdumpwidth_)
          continue;

        if ((tks.tkwt[i] > 0) && (tks.sum_Z[i] > 0)) {
          //double p=pik(beta,tks[i],*k);
          double p = y.rho[ivertex] * local_exp(-beta * Eik(tks.zpca[i], y.zvtx[ivertex], tks.dz2[i])) / tks.sum_Z[i];
          if (p > 0.0001) {
            std::cout << setw(8) << setprecision(3) << p;
          } else {
            std::cout << "    .   ";
          }
          E += p * Eik(tks.zpca[i], y.zvtx[ivertex], tks.dz2[i]);
          sump += p;
        } else {
          std::cout << "        ";
        }
      }
      std::cout << "  ( " << std::setw(3) << tks.kmin[i] << "," << std::setw(3) << tks.kmax[i] - 1 << " ) ";
      std::cout << endl;
    }
    std::cout << "                                                                   ";
    for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
      if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) < zdumpwidth_) {
        std::cout << "   " << setw(3) << ivertex << "  ";
      }
    }
    std::cout << endl;
    std::cout << "                                                                z= ";
    std::cout << setprecision(4);
    for (unsigned int ivertex = 0; ivertex < nv; ++ivertex) {
      if (std::fabs(y.zvtx[ivertex] - zdumpcenter_) < zdumpwidth_) {
        std::cout << setw(8) << fixed << y.zvtx[ivertex];
      }
    }
    std::cout << endl;
    std::cout << endl
              << "T=" << 1 / beta << " E=" << E << " n=" << y.getSize() << "  F= " << F << endl
              << "----------" << endl;
  }
#endif
}

void DAClusterizerInZ_vect::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.addUntracked<double>("zdumpcenter", 0.);
  desc.addUntracked<double>("zdumpwidth", 20.);
  desc.add<double>("d0CutOff", 3.0);
  desc.add<double>("Tmin", 2.0);
  desc.add<double>("delta_lowT", 0.001);
  desc.add<double>("zmerge", 0.01);
  desc.add<double>("dzCutOff", 3.0);
  desc.add<double>("Tpurge", 2.0);
  desc.add<int>("convergence_mode", 0);
  desc.add<double>("delta_highT", 0.01);
  desc.add<double>("Tstop", 0.5);
  desc.add<double>("coolingFactor", 0.6);
  desc.add<double>("vertexSize", 0.006);
  desc.add<double>("uniquetrkweight", 0.8);
  desc.add<double>("uniquetrkminp", 0.0);
  desc.add<double>("zrange", 4.0);
  desc.add<bool>("runInBlocks", false);
  desc.add<unsigned int>("block_size", 10000);
  desc.add<double>("overlap_frac", 0.0);
}
