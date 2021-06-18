#ifndef KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_

#include <vector>

#include "Constants.hpp"

class Daughter {
 public:
  Daughter() = default;
  Daughter(const Daughter&) = default;
  Daughter(Daughter&&) = default;
  Daughter& operator=(Daughter&&) = default;
  Daughter& operator=(const Daughter&) = default;
  ~Daughter() = default;

  Daughter(Pdg_t pdg_hypo, std::vector<Pdg_t> pids={}) : pdg_hypo_(pdg_hypo),
                                                         pids_(pids) {
  if(pids.empty())
    pids_ = {pdg_hypo};                                                           
  }

  Pdg_t GetPdgHypo() const { return pdg_hypo_; }
  const std::vector<Pdg_t>& GetPids() const { return pids_; }
  float GetCutChi2Prim() const { return chi2_prim_; }
  float GetCutCos() const { return cos_; }
  float GetCutInvMass() const { return invmass_discrepancy_; }
  int GetId() const { return id_; }

  void SetCutChi2Prim(float value) { chi2_prim_ = value; }
  void SetCutCos(float value) { cos_ = value; }
  void SetCutInvMass(float value) { invmass_discrepancy_ = value; }
  void CancelCutChi2Prim() { this->SetCutChi2Prim(-huge_value); }
  void CancelCutCos() { this->SetCutCos(-huge_value); }
  void CancelCutInvMass()  { this->SetCutInvMass(huge_value); }
  void SetId(int id) { id_ = id; }
  void CancelCuts();

 protected:
  Pdg_t pdg_hypo_{-1};                    ///< PDG code hypothesis
  std::vector<Pdg_t> pids_{};             ///< vector of PDG codes to use
  float chi2_prim_{18.4207};              ///< \f$\chi^2\f$ lower value
  float cos_{0.f};                        ///< cosine lower value
  float invmass_discrepancy_{huge_value};
  int id_{-1};                            ///< daughther number (0, 1, 2)
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
