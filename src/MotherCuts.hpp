#ifndef KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_

#include <limits>

class MotherCuts {

 public:
  MotherCuts() = default;
  MotherCuts(const MotherCuts&) = default;
  MotherCuts(MotherCuts&&) = default;
  MotherCuts& operator=(MotherCuts&&) = default;
  MotherCuts& operator=(const MotherCuts&) = default;
  ~MotherCuts() = default;

  MotherCuts(float distance, float chi_2_geo, float ldl, float chi_2_topo)
    : distance_(distance), chi2_geo_(chi_2_geo), ldl_(ldl), chi2_topo_(chi_2_topo) {}

  void SetDistance(float distance) { distance_ = distance; }
  void SetChi2Geo(float chi_2_geo) { chi2_geo_ = chi_2_geo; }
  void SetLdL(float ldl) { ldl_ = ldl; }
  void SetChi2Topo(float chi_2_topo) { chi2_topo_ = chi_2_topo; }

  float GetDistance() const { return distance_; }
  float GetChi2Geo() const { return chi2_geo_; }
  float GetLdL() const { return ldl_; }
  float GetChi2Topo() const { return chi2_topo_; }

 protected:
  float distance_{std::numeric_limits<float>::max()};
  float chi2_geo_{0.f};
  float ldl_{0.f};
  float chi2_topo_{std::numeric_limits<float>::max()};
};

#endif //KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
