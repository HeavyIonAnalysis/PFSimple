#include "Mother.hpp"

void Mother::CancelCuts() {
  this->CancelCutDistance();
  this->CancelCutDistanceToSV();
  this->CancelCutChi2Geo();
  this->CancelCutChi2GeoSM();
  this->CancelCutLdL();
  this->CancelCutDecayLength();
  this->CancelCutDistancePVLine();
  this->CancelCutChi2Topo();
  this->CancelCutChi2TopoSM();
  this->CancelCutCosTopo();
  this->CancelCutCosTopoSM();
  this->CancelCutCosOpen();
  this->CancelCutCosOpenSM();
}
