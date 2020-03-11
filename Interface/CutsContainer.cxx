#include "CutsContainer.h"

void CutsContainer::CancelCuts()
{
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
  CancelCutCosineTopo();
  CancelCutSigmaMassRatio();
  CancelCutChi2Topo();  
}
