#ifndef KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_

#include <utility>
#include <vector>
#include <string>
#include <stdexcept>

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

  Decay(std::string name, Pdg_t pdg, const MotherCuts& mother, std::vector<DaughterCuts> daughters) :
    name_(std::move(name)),
    pdg_(pdg),
    mother_(mother),
    daughters_(std::move(daughters)) {}

  const std::vector<DaughterCuts>& GetDaughters() const { return daughters_; }
  const MotherCuts& GetMother() const { return mother_; }

  void SetDaughters(const std::vector<DaughterCuts>& daughters) {
    if(!daughters_.empty()){
      throw std::runtime_error("Daughters are already set!");
    }
    daughters_ = daughters;
    int i{0};
    for(auto& daughter : daughters_){
      daughter.SetId(i++);
    }
  }
  void SetMother(const MotherCuts& mother) { mother_ = mother; }

  int GetNDaughters() const { return daughters_.size(); }
  Pdg_t GetPdg() const { return pdg_; }

 protected:
  std::string name_;  ///< decay name, to be used in branch name for example
  Pdg_t pdg_{-1};

  MotherCuts mother_;  ///< cuts for mother particle
  std::vector<DaughterCuts> daughters_{}; ///< cuts for daughter particles
};



#endif //KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
