#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

void PFSimpleTask::Exec() {
  pf_simple_ = new SimpleFinderNew();
  
  pf_simple_->SetDecays(decays_);

  pf_simple_->Init(in_task_->GetInputContainer());
  pf_simple_->FindParticles();
}
