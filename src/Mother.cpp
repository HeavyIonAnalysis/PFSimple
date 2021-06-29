#include "Mother.hpp"

void Mother::CancelCuts(){
  this->CancelCutDistance();
  this->CancelCutChi2Geo();
  this->CancelCutLdL();
  this->CancelCutChi2Topo();
  this->CancelCutChi2TopoLower();
  this->CancelCutInvMass();
}
