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

  Daughter(Pdg_t pdg_hypo, std::vector<Pdg_t> pids, float chi_2_prim, float cos) : pdg_hypo_(pdg_hypo),
                                                                                   pids_(std::move(pids)),
                                                                                   chi2_prim_(chi_2_prim),
                                                                                   cos_(cos) {}

  Pdg_t GetPdgHypo() const { return pdg_hypo_; }
  const std::vector<Pdg_t>& GetPids() const { return pids_; }
  float GetChi2Prim() const { return chi2_prim_; }
  float GetCos() const { return cos_; }
  int GetId() const { return id_; }

  void SetId(int id) { id_ = id; }

 protected:
  Pdg_t pdg_hypo_{-1};       ///< PDG code hypothesis
  std::vector<Pdg_t> pids_{};///< vector of PDG codes to use
  float chi2_prim_{0.f};     ///< \f$\chi^2\f$ lower value
  float cos_{0.f};           ///< cosine lower value
  int id_{-1};               ///< daughther number (0, 1, 2)
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
