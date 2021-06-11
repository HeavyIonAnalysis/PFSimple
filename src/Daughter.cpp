#include "Daughter.hpp"


void Daughter::CancelCuts(){
  this->CancelCutChi2Prim();
  this->CancelCutCos();
  this->CancelCutInvMass();
}