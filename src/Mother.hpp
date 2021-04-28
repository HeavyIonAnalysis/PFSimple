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

  void SetCutDistance(float distance) { distance_ = distance; }
  void SetCutChi2Geo(float chi_2_geo) { chi2_geo_ = chi_2_geo; }
  void SetCutLdL(float ldl) { ldl_ = ldl; }
  void SetCutChi2Topo(float chi_2_topo) { chi2_topo_ = chi_2_topo; }
  void CancelCutDistance() { this->SetCutDistance(huge_value); }
  void CancelCutChi2Geo() { this->SetCutChi2Geo(huge_value); }
  void CancelCutLdL() { this->SetCutLdL(-huge_value); }
  void CancelCutChi2Topo() { this->SetCutChi2Topo(huge_value); }
  void CancelCuts();

  Pdg_t GetPdg() const { return pdg_; }
  float GetCutDistance() const { return distance_; }
  float GetCutChi2Geo() const { return chi2_geo_; }
  float GetCutLdL() const { return ldl_; }
  float GetCutChi2Topo() const { return chi2_topo_; }

 protected:
  Pdg_t pdg_{-1};

  float distance_{1.f};         ///< lower value
  float chi2_geo_{3.f};         ///< upper value
  float ldl_{5.f};              ///< upper value
  float chi2_topo_{huge_value}; ///< lower value
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
