#include "InputContainer.hpp"
#include <Constants.hpp>

#include <iostream>
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
                              int id,
                              int nhits) {
  if (pdg == 0 || pdg == -2)
    return;

  if (par.size() != kNumberOfTrackPars || cov.size() != NumberOfCovElements || field.size() != kNumberOfFieldPars) {
    throw std::runtime_error("Wrong size of input vector!");
  }

  KFParticle particle;
  particle.X() = par[kX];
  particle.Y() = par[kY];
  particle.Z() = par[kZ];
  particle.Px() = par[kPx];
  particle.Py() = par[kPy];
  particle.Pz() = par[kPz];

  for (int i = 0; i < 21; i++)
    particle.Covariance(i) = cov[i];

  for (int i = 0; i < kNumberOfFieldPars; i++)
    particle.SetFieldCoeff(field[i], i);

  particle.Q() = char(charge);//NOTE: is not safe
  particle.SetPDG(pdg);
  particle.SetId(id);

  tracks_.push_back(particle);
}

double InputContainer::InversedChi2Prob(double p, int ndf) {
  const double epsilon = 1.e-14;
  double chi2Left = 0.f;
  double chi2Right = 10000.f;

  double probLeft = p - TMath::Prob(chi2Left, ndf);

  double chi2Centr = (chi2Left + chi2Right) / 2.f;
  double probCentr = p - TMath::Prob(chi2Centr, ndf);

  while (TMath::Abs(chi2Right - chi2Centr) / chi2Centr > epsilon) {
    if (probCentr * probLeft > 0.f) {
      chi2Left = chi2Centr;
      probLeft = probCentr;
    } else {
      chi2Right = chi2Centr;
    }

    chi2Centr = (chi2Left + chi2Right) / 2.f;
    probCentr = p - TMath::Prob(chi2Centr, ndf);
  }

  return chi2Centr;
}

void InputContainer::Reserve(size_t n) {
  tracks_.reserve(n);
}

// KFParticleTopoReconstructor* InputContainer::CreateTopoReconstructor()
// {
//   /*
//    * Creates the pointer on the KFParticleTopoReconstructor object
//    * with all necessary input information in order to perform particle selection using
//    * non-simplified "standard" KFParticle algorithm.
//    */
//   auto* TR = new KFParticleTopoReconstructor;
//
//   // cuts setting
//   TR -> GetKFParticleFinder() -> SetChiPrimaryCut2D(cuts_.GetCutChi2PrimPos());
//   TR -> GetKFParticleFinder() -> SetMaxDistanceBetweenParticlesCut(cuts_.GetCutDistance());
//   TR -> GetKFParticleFinder() -> SetChi2Cut2D(cuts_.GetCutChi2Geo());
//   TR -> GetKFParticleFinder() -> SetLCut(cuts_.GetCutLDown());
//   TR -> GetKFParticleFinder() -> SetLdLCut2D(cuts_.GetCutLdL());
//
//   KFPTrackVector track_tmp, track_empty;
//   track_tmp.Resize(tracks_.size());
//   for(int iTr=0; iTr<tracks_.size(); iTr++)
//   {
//     for(Int_t iP=0; iP<6; iP++)
//       track_tmp.SetParameter(tracks_[iTr].GetParameter(iP), iP, iTr);
//     for(Int_t iC=0; iC<21; iC++)
//       track_tmp.SetCovariance(tracks_[iTr].GetCovariance(iC), iC, iTr);
//     for(Int_t iF=0; iF<10; iF++)
//       track_tmp.SetFieldCoefficient(tracks_[iTr].GetFieldCoeff()[iF], iF, iTr);
//     track_tmp.SetPDG(tracks_[iTr].GetPDG(), iTr);
//     track_tmp.SetQ(tracks_[iTr].GetQ(), iTr);
//     track_tmp.SetPVIndex(-1, iTr);
//     track_tmp.SetId(tracks_[iTr].Id(), iTr);
//   }
//   TR->Init(track_tmp, track_empty);
//
//   TR->AddPV(vtx_);
//
//   return TR;
// }

// SimpleFinder InputContainer::CreateSimpleFinder()
// {
//   /*
//    * Creates the SimpleFinder object with all necessary input information in order to
//    * perform particle selection using src algorithm.
//    */
//
//   SimpleFinder FCF;
//
//   KFPTrackVector track_tmp;
//   track_tmp.Resize(tracks_.size());
//
//   for(int iTr=0; iTr<tracks_.size(); iTr++)
//   {
//     for(Int_t iP=0; iP<6; iP++)
//       track_tmp.SetParameter(tracks_[iTr].GetParameter(iP), iP, iTr);
//     for(Int_t iC=0; iC<21; iC++)
//       track_tmp.SetCovariance(tracks_[iTr].GetCovariance(iC), iC, iTr);
//     for(Int_t iF=0; iF<10; iF++)
//       track_tmp.SetFieldCoefficient(tracks_[iTr].GetFieldCoeff()[iF], iF, iTr);
//     track_tmp.SetPDG(tracks_[iTr].GetPDG(), iTr);
//     track_tmp.SetQ(tracks_[iTr].GetQ(), iTr);
//     track_tmp.SetPVIndex(-1, iTr);
//     track_tmp.SetId(tracks_[iTr].Id(), iTr);
//   }
//   FCF.Init(track_tmp, vtx_);
//   FCF.SetCuts(cuts_);
//
//   return FCF;
// }
