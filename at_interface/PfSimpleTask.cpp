#include "PfSimpleTask.hpp"

#include "SimpleFinderNew.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

void PFSimpleTask::Exec() {
  SimpleFinderNew pf_simple_;

  pf_simple_.SetDecays(decays_);

  pf_simple_.Init(in_task_->GetInputContainer());
  pf_simple_.FindParticles();
  out_task_->SetCandidates(pf_simple_.GetCandidates());
}
