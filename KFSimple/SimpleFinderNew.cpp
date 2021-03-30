#include "SimpleFinderNew.hpp"

void SimpleFinderNew::Init(KFPTrackVector&& tracks, const KFVertex& pv) {
  tracks_ = tracks;
  prim_vx_ = pv;
}

void SimpleFinderNew::Init(const InputContainer& input) {
  const std::vector<KFParticle>& tracks = input.GetTracks();
  KFPTrackVector track_tmp;
  track_tmp.Resize(tracks.size());

  for (size_t iTr = 0; iTr < tracks.size(); iTr++) {
    for (Int_t iP = 0; iP < 6; iP++)
      track_tmp.SetParameter(tracks[iTr].GetParameter(iP), iP, iTr);
    for (Int_t iC = 0; iC < 21; iC++)
      track_tmp.SetCovariance(tracks[iTr].GetCovariance(iC), iC, iTr);
    for (Int_t iF = 0; iF < 10; iF++)
      track_tmp.SetFieldCoefficient(tracks[iTr].GetFieldCoeff()[iF], iF, iTr);
    track_tmp.SetPDG(tracks[iTr].GetPDG(), iTr);
    track_tmp.SetQ(tracks[iTr].GetQ(), iTr);
    track_tmp.SetPVIndex(-1, iTr);
    track_tmp.SetId(tracks[iTr].Id(), iTr);
  }

  Init(std::move(track_tmp), input.GetVertex());
//  SetDecay(input.GetDecay());
//  SetCuts(input.GetCuts());
}
