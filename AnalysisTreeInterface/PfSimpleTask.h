#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_

#include "SimpleFinder.h"

#include "AnalysisTree/Task.hpp"
#include "ConverterIn.h"
#include "ConverterOut.h"

class PFSimpleTask : public AnalysisTree::Task {

 public:
  PFSimpleTask() = default;
  ~PFSimpleTask() override = default;

  void Init() override {};
  void Exec() override {

    SimpleFinder pf_simple_;
    pf_simple_.Init(in_task_->CreateInputContainer());
    pf_simple_.SetDecay(in_task_->GetDecay());
    pf_simple_.SetCuts(in_task_->GetCuts());
    pf_simple_.SortTracks();
    pf_simple_.FindParticles();

    out_task_->SetDecay(in_task_->GetDecay());
    out_task_->SetCandidates(pf_simple_.GetMotherCandidates());

  };
  void Finish() override {};
  void SetInTask(ConverterIn* in_task) { in_task_ = in_task; }
  void SetOutTask(ConverterOut* out_task) { out_task_ = out_task; }

 protected:

  ConverterIn* in_task_{nullptr};
  ConverterOut* out_task_{nullptr};

};

#endif //KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
