#ifndef KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_

#include <vector>

#include "Constants.hpp"

class DaughterCuts {
 public:
  DaughterCuts() = default;
  DaughterCuts(const DaughterCuts&) = default;
  DaughterCuts(DaughterCuts&&) = default;
  DaughterCuts& operator=(DaughterCuts&&) = default;
  DaughterCuts& operator=(const DaughterCuts&) = default;
  ~DaughterCuts() = default;

  DaughterCuts(Pdg_t pdg_hypo, std::vector<Pdg_t> pids, float chi_2_prim) : pdg_hypo_(pdg_hypo), pids_(std::move(pids)), chi2_prim_(chi_2_prim) {}

  Pdg_t GetPdgHypo() const { return pdg_hypo_; }
  const std::vector<Pdg_t>& GetPids() const { return pids_; }
  float GetChi2Prim() const { return chi2_prim_; }

 protected:
  Pdg_t pdg_hypo_{-1};
  std::vector<Pdg_t> pids_{};
  float chi2_prim_{0};
};

#endif //KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
