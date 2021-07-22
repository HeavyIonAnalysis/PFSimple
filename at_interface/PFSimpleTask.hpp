#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_

#include "Decay.hpp"
#include "SimpleFinderNew.hpp"

#include "AnalysisTree/Task.hpp"

class ConverterIn;

class PFSimpleTask : public AnalysisTree::Task {

 public:
  PFSimpleTask() = default;
  ~PFSimpleTask() override = default;

  void Init() override{};
  void Exec() override;
  void Finish() override{};
  void SetInTask(ConverterIn* in_task) { in_task_ = in_task; }
  SimpleFinderNew* GetSimpleFinder() { return pf_simple_; }
  void SetDecays(const std::vector<Decay>& decays) { decays_ = decays; }
  void AddDecay(const Decay& decay) { decays_.emplace_back(decay); }

 protected:
  ConverterIn* in_task_{nullptr};
  SimpleFinderNew* pf_simple_{nullptr};
  std::vector<Decay> decays_{};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
