#ifndef KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_

#include <limits>
#include <string>
#include <vector>

#include "Constants.hpp"

class Mother {

 public:
  Mother() = default;
  Mother(const Mother&) = default;
  Mother(Mother&&) = default;
  Mother& operator=(Mother&&) = default;
  Mother& operator=(const Mother&) = default;
  ~Mother() = default;

  explicit Mother(Pdg_t pdg) : pdg_(pdg) {
    mass_pdg_ = -1.0;
    mass_pdg_sigma_ = -1.0;
  }

  void SetMassPdg(float mass_pdg) { mass_pdg_ = mass_pdg; }
  void SetMassPdgSigma(float mass_pdg_sigma) { mass_pdg_sigma_ = mass_pdg_sigma; }

  void SetCutDistance(float distance) { distance_ = distance; }
  void SetCutDistanceToSV(float distance_sv) { distance_sv_ = distance_sv; }
  void SetCutChi2Geo(float chi_2_geo) { chi2_geo_.at(0) = chi_2_geo; }
  void SetCutChi2GeoSM(std::vector<float> chi_2_geo_sm) {
    for (std::size_t i = 1; i < chi_2_geo_sm.size() + 1; ++i) chi2_geo_.at(i) = chi_2_geo_sm.at(i - 1);
  }
  void SetCutLdL(float ldl) { ldl_ = ldl; }
  void SetCutDecayLength(float l) { l_ = l; }
  void SetCutDistancePVLine(float distance_pv) { distance_pv_ = distance_pv; }
  void SetCutChi2Topo(float chi_2_topo) { chi2_topo_.at(0) = chi_2_topo; }
  void SetCutChi2TopoSM(std::vector<float> chi_2_topo_sm) {
    for (std::size_t i = 1; i < chi_2_topo_sm.size() + 1; ++i) chi2_topo_.at(i) = chi_2_topo_sm.at(i - 1);
  }
  void SetCutCosTopo(float cos_topo) { cos_topo_.at(0) = cos_topo; }
  void SetCutCosTopoSM(std::vector<float> cos_topo_sm) {
    for (std::size_t i = 1; i < cos_topo_sm.size() + 1; ++i) cos_topo_.at(i) = cos_topo_sm.at(i - 1);
  }
  void SetCutCosOpen(float cos_open) { cos_topo_.at(0) = cos_open; }
  void SetCutCosOpenSM(std::vector<float> cos_open_sm) {
    for (std::size_t i = 1; i < cos_open_sm.size() + 1; ++i) cos_open_.at(i) = cos_open_sm.at(i - 1);
  }
  void SetCutInvMass(float value) { invmass_discrepancy_ = value; }

  void CancelCutDistance() { this->SetCutDistance(huge_value); }
  void CancelCutDistanceToSV() { this->SetCutDistanceToSV(huge_value); }
  void CancelCutChi2Geo() { this->SetCutChi2Geo(huge_value); }
  void CancelCutChi2GeoSM() { this->SetCutChi2GeoSM({huge_value, huge_value, huge_value}); }
  void CancelCutLdL() { this->SetCutLdL(-huge_value); }
  void CancelCutDecayLength() { this->SetCutDecayLength(-huge_value); }
  void CancelCutDistancePVLine() { this->SetCutDistancePVLine(-huge_value); }
  void CancelCutChi2Topo() { this->SetCutChi2Topo(huge_value); }
  void CancelCutChi2TopoSM() { this->SetCutChi2TopoSM({-huge_value, -huge_value, -huge_value}); }
  void CancelCutCosTopo() { this->SetCutCosTopo(-huge_value); }
  void CancelCutCosTopoSM() { this->SetCutCosTopoSM({huge_value, huge_value, huge_value}); }
  void CancelCutCosOpen() { this->SetCutCosOpen(-huge_value); }
  void CancelCutCosOpenSM() { this->SetCutCosOpenSM({-huge_value, -huge_value, -huge_value}); }
  void CancelCutInvMass() { this->SetCutInvMass(huge_value); }
  void CancelCuts();

  Pdg_t GetPdg() const { return pdg_; }
  float GetMassPdg() {
    if (mass_pdg_ < 0) mass_pdg_ = GetMassKFDatabase();
    return mass_pdg_;
  };
  float GetMassPdgSigma() {
    if (mass_pdg_sigma_ < 0) mass_pdg_sigma_ = GetMassSigmaKFDatabase();
    return mass_pdg_sigma_;
  }
  float GetMassKFDatabase() const;
  float GetMassSigmaKFDatabase() const;
  float GetCutDistance() const { return distance_; }
  float GetCutDistanceToSV() const { return distance_sv_; }
  std::array<float, 4> GetCutChi2Geo() const { return chi2_geo_; }
  std::array<float, 4> GetCutCosOpen() const { return cos_open_; }
  float GetCutLdL() const { return ldl_; }
  float GetCutDecayLength() const { return l_; }
  float GetCutDistancePVLine() const { return distance_pv_; }
  std::array<float, 4> GetCutChi2Topo() const { return chi2_topo_; }
  std::array<float, 4> GetCutCosTopo() const { return cos_topo_; }
  float GetCutInvMass() const { return invmass_discrepancy_; }

 protected:
  Pdg_t pdg_{-1};
  float mass_pdg_{-1.0};
  float mass_pdg_sigma_{-1.0};
  float mass_kaon_{0.497614};// default mass value in KFParticleDatabase

  float distance_{1.f};///< lower value
  float distance_sv_{1.f};
  std::array<float, 4> chi2_geo_{{3.f, huge_value, huge_value, huge_value}};///< upper value
  float ldl_{5.f};                                                          ///< upper value
  float l_{-huge_value};
  float distance_pv_{-huge_value};
  std::array<float, 4> chi2_topo_{{huge_value, -huge_value, -huge_value, -huge_value}};///< lower value
  std::array<float, 4> cos_topo_{{-huge_value, huge_value, huge_value, huge_value}};
  std::array<float, 4> cos_open_{{-huge_value, -huge_value, -huge_value, -huge_value}};
  float invmass_discrepancy_{huge_value};
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
