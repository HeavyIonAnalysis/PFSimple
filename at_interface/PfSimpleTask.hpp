#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_

#include "SimpleFinderNew.hpp"
//#include "SimpleFinder.hpp"

#include "AnalysisTree/Task.hpp"
#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

class PFSimpleTask : public AnalysisTree::Task {

 public:
  PFSimpleTask() = default;
  ~PFSimpleTask() override = default;

  void Init() override{};
  void Exec() override {

    SimpleFinderNew pf_simple_;

    DaughterCuts proton(2212, {2212}, 18.6);
    DaughterCuts pion(-211, {-211}, 18.6);
    Decay lambda("lambda", MotherCuts(1, 3, 5, 3), {proton, pion});
    pf_simple_.AddDecay(lambda);

    pf_simple_.Init(in_task_->GetInputContainer());

    pf_simple_.FindParticles();

//    out_task_->SetDecay(in_task_->GetDecay());
    out_task_->SetCandidates(pf_simple_.GetCandidates());
  };
  void Finish() override{};
  void SetInTask(ConverterIn* in_task) { in_task_ = in_task; }
  void SetOutTask(ConverterOut* out_task) { out_task_ = out_task; }

 protected:
  ConverterIn* in_task_{nullptr};
  ConverterOut* out_task_{nullptr};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
