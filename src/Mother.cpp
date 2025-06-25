#include "Mother.hpp"
#include <KFParticleDatabase.h>
using std::to_string;

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
  this->CancelCutInvMass();
}

float Mother::GetMassKFDatabase() const {
  if (mass_pdg_ > 0)
    return mass_pdg_;

  float massPDG, massPDGsigma;
  KFParticleDatabase::Instance()->GetMotherMass(pdg_, massPDG, massPDGsigma);

  if (massPDG == mass_kaon_ && pdg_ != 310)
    throw std::runtime_error("Mass for mother pdg " + to_string(pdg_) + " is not availabe in KFParticleDatabase. Set mass manually with Mother::SetMassPdg(mass).");

  return massPDG;
};

float Mother::GetMassSigmaKFDatabase() const {
  if (mass_pdg_sigma_ > 0)
    return mass_pdg_sigma_;

  float massPDG, massPDGsigma;
  KFParticleDatabase::Instance()->GetMotherMass(pdg_, massPDG, massPDGsigma);

  if (massPDG == mass_kaon_ && pdg_ != 310)
    throw std::runtime_error("Mass sigma for mother pdg " + to_string(pdg_) + "  is not availabe in KFParticleDatabase. Set mass sigma manually with Mother::SetMassPdgSigma(mass_sigma).");

  return massPDGsigma;
};
