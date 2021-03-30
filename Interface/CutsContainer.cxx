#include "CutsContainer.h"

void CutsContainer::CancelCuts() {
  CancelCutChi2PrimPos();
  CancelCutChi2PrimNeg();
  CancelCutDistance();
  CancelCutCosineDaughterPos();
  CancelCutCosineDaughterNeg();
  CancelCutChi2Geo();
  CancelCutLUp();
  CancelCutLDown();
  CancelCutLdL();
  CancelCutIsFromPV();
  CancelCutCosineTopoDown();
  CancelCutCosineTopoUp();
  CancelCutSigmaMassRatio();
  CancelCutChi2Topo();

  CancelCutChi2PrimThird();
  CancelCutDistanceThird();
  CancelCutCosineDaughterThird();
  CancelCutChi2GeoThree();
  CancelCutLThreeUp();
  CancelCutLThreeDown();
  CancelCutLdLThree();
  CancelCutIsFromPVThree();
  CancelCutCosineTopoThreeDown();
  CancelCutCosineTopoThreeUp();
  CancelCutChi2TopoThree();
}
