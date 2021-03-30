#ifndef KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_

#include <utility>
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

  DaughterCuts(std::vector<Pdg_t> pids, float chi_2_prim) : pids_(std::move(pids)), chi2_prim_(chi_2_prim) {}

 protected:
  std::vector<Pdg_t> pids_{};
  float chi2_prim_{-1};
};

#endif //KFPARTICLESIMPLE_KFSIMPLE_DAUGHTERCUTS_HPP_
