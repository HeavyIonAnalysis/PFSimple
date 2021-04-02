#ifndef KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_

#include <utility>
#include <vector>
#include <string>

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

  Decay(std::string name, const MotherCuts& mother, std::vector<DaughterCuts> daughters)
    : name_(std::move(name)), mother_(mother), daughters_(std::move(daughters)) {}

  const std::vector<DaughterCuts>& GetDaughters() const { return daughters_; }
  const MotherCuts& GetMother() const { return mother_; }

  void SetDaughters(const std::vector<DaughterCuts>& daughters) { daughters_ = daughters; }
  void SetMother(const MotherCuts& mother) { mother_ = mother; }

  int GetNDaughters() const { return daughters_.size(); }

 protected:
  std::string name_;

  MotherCuts mother_;
  std::vector<DaughterCuts> daughters_{};
};



#endif //KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
