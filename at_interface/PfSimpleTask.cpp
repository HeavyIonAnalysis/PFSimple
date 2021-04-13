#include "PfSimpleTask.hpp"

#include "SimpleFinderNew.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

void PFSimpleTask::Exec() {
  SimpleFinderNew pf_simple_;

  DaughterCuts proton(2212, {2212}, 18.6, 0.99);
  DaughterCuts pion(-211, {-211}, 18.6, 0.f);
  Decay lambda("lambda", MotherCuts(1, 3, 5, 3), {proton, pion});
  pf_simple_.AddDecay(lambda);

  pf_simple_.Init(in_task_->GetInputContainer());
  pf_simple_.FindParticles();
  out_task_->SetCandidates(pf_simple_.GetCandidates());
}
