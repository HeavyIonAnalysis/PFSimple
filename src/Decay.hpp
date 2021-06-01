#ifndef KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Daughter.hpp"
#include "Mother.hpp"

class Decay {

 public:
  Decay() = default;
  Decay(const Decay&) = default;
  Decay(Decay&&) = default;
  Decay& operator=(Decay&&) = default;
  Decay& operator=(const Decay&) = default;
  ~Decay() = default;

  Decay(std::string name, const Mother& mother, std::vector<Daughter> daughters) : name_(std::move(name)),
                                                                                   mother_(mother),
                                                                                   daughters_(std::move(daughters)) {
    int i{0};
    for(auto& daughter : daughters_){
      daughter.SetId(i++);
    }
  }

  const std::vector<Daughter>& GetDaughters() const { return daughters_; }
  const Mother& GetMother() const { return mother_; }

  void SetDaughters(const std::vector<Daughter>& daughters) {
    if (!daughters_.empty()) {
      throw std::runtime_error("Daughters are already set!");
    }
    daughters_ = daughters;
    int i{0};
    for (auto& daughter : daughters_) {
      daughter.SetId(i++);
    }
  }
  void SetMother(const Mother& mother) { mother_ = mother; }

  int GetNDaughters() const { return daughters_.size(); }

 protected:
  std::string name_;///< decay name, to be used in branch name for example

  Mother mother_;                    ///< cuts for mother particle
  std::vector<Daughter> daughters_{};///< cuts for daughter particles
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
