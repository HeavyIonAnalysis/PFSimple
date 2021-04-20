#ifndef KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_

#include <limits>

#include "Constants.hpp"

class Mother {

 public:
  Mother() = default;
  Mother(const Mother&) = default;
  Mother(Mother&&) = default;
  Mother& operator=(Mother&&) = default;
  Mother& operator=(const Mother&) = default;
  ~Mother() = default;

  explicit Mother(Pdg_t pdg) : pdg_(pdg) {}

//  Mother(float distance, float chi_2_geo, float ldl, float chi_2_topo)
//      : distance_(distance), chi2_geo_(chi_2_geo), ldl_(ldl), chi2_topo_(chi_2_topo) {}

  void SetDistance(float distance) { distance_ = distance; }
  void SetChi2Geo(float chi_2_geo) { chi2_geo_ = chi_2_geo; }
  void SetLdL(float ldl) { ldl_ = ldl; }
  void SetChi2Topo(float chi_2_topo) { chi2_topo_ = chi_2_topo; }

  Pdg_t GetPdg() const { return pdg_; }
  float GetDistance() const { return distance_; }
  float GetChi2Geo() const { return chi2_geo_; }
  float GetLdL() const { return ldl_; }
  float GetChi2Topo() const { return chi2_topo_; }

 protected:
  Pdg_t pdg_{-1};

  float distance_{std::numeric_limits<float>::max()}; ///< lower value
  float chi2_geo_{std::numeric_limits<float>::min()}; ///< upper value
  float ldl_{std::numeric_limits<float>::min()};      ///< upper value
  float chi2_topo_{std::numeric_limits<float>::max()};///< lower value
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
