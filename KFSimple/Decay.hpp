#ifndef KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_

#include <vector>

#include "DaughterCuts.hpp"
#include "MotherCuts.hpp"

class Decay{

 public:
  Decay() = default;
  Decay(const Decay&) = default;
  Decay(Decay&&) = default;
  Decay& operator=(Decay&&) = default;
  Decay& operator=(const Decay&) = default;
  ~Decay() = default;

 protected:
  std::vector<DaughterCuts> daughters_{};
  MotherCuts mother_;
};



#endif //KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
