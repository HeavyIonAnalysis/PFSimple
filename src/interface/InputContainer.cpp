#include "InputContainer.hpp"
#include <Constants.hpp>

#include <stdexcept>

#include "TMath.h"

void InputContainer::SetPV(float x, float y, float z) {
  KFPVertex primVtx_tmp;
  primVtx_tmp.SetXYZ(x, y, z);
  primVtx_tmp.SetCovarianceMatrix(0, 0, 0, 0, 0, 0);//NOTE
  primVtx_tmp.SetNContributors(0);
  primVtx_tmp.SetChi2(-100);

  vtx_ = KFVertex(primVtx_tmp);
}

void InputContainer::AddTrack(const std::vector<float>& par,
                              const std::vector<float>& cov,
                              const std::vector<float>& field,
                              int charge,
                              int pdg,
                              int id) {
  if (par.size() != kNumberOfTrackPars || cov.size() != NumberOfCovElements || field.size() != NumberOfFieldPars) {
    throw std::runtime_error("Wrong size of input vector!");
  }

  KFParticle particle;
  particle.X() = par[kX];
  particle.Y() = par[kY];
  particle.Z() = par[kZ];
  particle.Px() = par[kPx];
  particle.Py() = par[kPy];
  particle.Pz() = par[kPz];

  for (int i = 0; i < NumberOfCovElements; i++)
    particle.Covariance(i) = cov[i];

  for (int i = 0; i < NumberOfFieldPars; i++)
    particle.SetFieldCoeff(field[i], i);

  particle.Q() = char(charge);//NOTE: is not safe
  particle.SetPDG(pdg);
  particle.SetId(id);

  tracks_.push_back(particle);
}

void InputContainer::Reserve(size_t n) {
  tracks_.reserve(n);
}
