#ifndef KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_

#include <iostream>
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
    for (auto& daughter : daughters_) {
      daughter.SetId(i++);
    }
  }

  const std::vector<Daughter>& GetDaughters() const { return daughters_; }
  const Mother& GetMother() const { return mother_; }
  const std::string GetName() const { return name_; }

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

  void SetIsApplyMassConstraint(bool is = true) { is_apply_mass_constraint_ = is; }
  void SetIsTransportToPV(bool is = true) { is_transport_to_pv_ = is; }
  void SetIsDoNotWriteMother(bool is = true) {
    std::cout << "WARNING!! If the mother particle is not written in the output-tree, MC-matching is not possible for upper-level mothers from cascade decays. " << std::endl;
    is_do_not_write_mother_ = is;
  }

  bool GetIsApplyMassConstraint() const { return is_apply_mass_constraint_; }
  bool GetIsTransportToPV() const { return is_transport_to_pv_; }
  bool GetIsDoNotWriteMother() const { return is_do_not_write_mother_; }

 protected:
  std::string name_;///< decay name, to be used in branch name for example

  Mother mother_;                    ///< cuts for mother particle
  std::vector<Daughter> daughters_{};///< cuts for daughter particles

  bool is_apply_mass_constraint_{false};
  bool is_transport_to_pv_{false};
  bool is_do_not_write_mother_{false};
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_DECAY_HPP_
