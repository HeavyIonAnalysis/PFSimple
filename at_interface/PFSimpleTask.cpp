#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

void PFSimpleTask::Init() {
  pf_simple_ = new SimpleFinder();
  pf_simple_->SetDecays(decays_);
}

void PFSimpleTask::Exec() {
  pf_simple_->Init(in_task_->GetInputContainer());
  pf_simple_->FindParticles();
}

void PFSimpleTask::Finish() {
  delete pf_simple_;
}
