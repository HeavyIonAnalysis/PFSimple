#include <KFParticleSIMD.h>
#include "SimpleFinderNew.hpp"

void SimpleFinderNew::Init(KFPTrackVector&& tracks, const KFVertex& pv) {
  tracks_ = tracks;
  prim_vx_ = pv;
  InitIndexesMap();
}

void SimpleFinderNew::SetTrack(const KFParticle& particle, int id, KFPTrackVector& tracks){
  for (Int_t iP = 0; iP < 6; iP++){
    tracks.SetParameter(particle.GetParameter(iP), iP, id);
  }
  for (Int_t iC = 0; iC < 21; iC++){
    tracks.SetCovariance(particle.GetCovariance(iC), iC, id);
  }
  for (Int_t iF = 0; iF < 10; iF++){
    tracks.SetFieldCoefficient(particle.GetFieldCoeff()[iF], iF, id);
  }
  tracks.SetPDG(particle.GetPDG(), id);
  tracks.SetQ(particle.GetQ(), id);
  tracks.SetPVIndex(-1, id);
  tracks.SetId(particle.Id(), id);
}


void SimpleFinderNew::Init(const InputContainer& input) {
  const std::vector<KFParticle>& tracks = input.GetTracks();
  KFPTrackVector track_tmp;
  track_tmp.Resize(tracks.size());
  for (size_t iTr = 0; iTr < tracks.size(); iTr++) {
    SetTrack(tracks[iTr], iTr, track_tmp);
  }
  Init(std::move(track_tmp), input.GetVertex());
}

float SimpleFinderNew::CalculateDistanceBetweenParticles(const Parameters_t& parameters) {
  const float dx = parameters.at(0).at(kX) - parameters.at(1).at(kX);
  const float dy = parameters.at(0).at(kY) - parameters.at(1).at(kY);
  const float dz = parameters.at(0).at(kZ) - parameters.at(1).at(kZ);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void SimpleFinderNew::CalculateParamsInPCA(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2) {
  params_ = Parameters_t(2);

  KFParticle particle1(track1, pid1);
  KFParticle particle2(track2, pid2);
  KFParticleSIMD particleSIMD1(particle1);// the same particle is copied to each SIMD element
  KFParticleSIMD particleSIMD2(particle2);

  float_v dS[2];
  float_v params1[8], params2[8];
  particleSIMD1.GetDStoParticleFast(particleSIMD2, dS);
  particleSIMD1.TransportFast(dS[0], params1);
  particleSIMD2.TransportFast(dS[1], params2);

  for (int i = 0; i < 8; i++) {
    params_.at(0).at(i) = params1[i][0];
    params_.at(1).at(i) = params2[i][0];
  }
}

KFParticleSIMD SimpleFinderNew::ConstructMother(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs) {
  const auto n = tracks.size();

  KFParticleSIMD mother;
  std::vector<KFParticle> particles{};
  std::vector<KFParticleSIMD> particles_simd{};

  for(size_t i=0; i<n; ++i){
    particles.emplace_back( KFParticle(tracks.at(i), pdgs.at(i)) );
    particles.at(i).SetId(tracks.at(i).Id());
    particles_simd.emplace_back(particles.at(i));
  }

  auto sv = GetSecondaryVertex();
  float_v vertex[3] = {sv[0], sv[1], sv[2]};
  for(size_t i=0; i<n; ++i) {
    float_v ds, dsdr[6];
    ds = particles_simd.at(i).GetDStoPoint(vertex, dsdr);
    particles_simd.at(i).TransportToDS(ds, dsdr);
  }

  if(n == 2){
    const KFParticleSIMD* vDaughtersPointer[2] = {&particles_simd.at(0), &particles_simd.at(1)};
    mother.Construct(vDaughtersPointer, 2, nullptr);
  }
  else if(n == 3){
    const KFParticleSIMD* vDaughtersPointer[3] = {&particles_simd.at(0), &particles_simd.at(1),  &particles_simd.at(2)};
    mother.Construct(vDaughtersPointer, 3, nullptr);
  }

  return mother;
}

float SimpleFinderNew::CalculateChiToPrimaryVertex(const KFPTrack& track, Pdg_t pid) const {
  // SIMD'ized version
  KFParticle tmpPart(track, pid);
  KFParticleSIMD tmpPartSIMD(tmpPart);

  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);
  const float_v point[3] = {prim_vx_Simd.X(), prim_vx_Simd.Y(), prim_vx_Simd.Z()};
  tmpPartSIMD.TransportToPoint(point);

  float_v chi2vec = tmpPartSIMD.GetDeviationFromVertex(prim_vx_Simd);

  return chi2vec[0];
}

void SimpleFinderNew::CalculateMotherProperties(const KFParticleSIMD& mother) {
  float_v l_Simd;
  float_v dl_Simd;
  float_m isFromPV_Simd;
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);

  mother.GetDistanceToVertexLine(prim_vx_Simd, l_Simd, dl_Simd, &isFromPV_Simd);

  KFParticleSIMD motherTopo = mother;
  motherTopo.SetProductionVertex(prim_vx_Simd);
  motherTopo.KFParticleBaseSIMD::GetDecayLength(l_Simd, dl_Simd);

  values_.l = l_Simd[0];
  values_.l_over_dl = l_Simd[0] / dl_Simd[0];
  values_.is_from_PV = isFromPV_Simd[0];
}

std::vector<int> SimpleFinderNew::GetIndexes(const DaughterCuts& cuts) {
  std::vector<int> result{};
  for(auto pid : cuts.GetPids()){
    auto it = indexes_.find(pid);
    if(it != indexes_.end()){
      for(auto i_track : it->second){
        auto track = GetTrack(i_track);
        if(IsGoodDaughter(track, cuts)){
          result.emplace_back(i_track);
        }
      }
    }
  }
  return result;
}

void SimpleFinderNew::InitIndexesMap() {
  for (int i = 0; i < tracks_.Size(); i++) {
    auto pdg = tracks_.PDG()[i];
    auto it = indexes_.find(pdg);
    if(it != indexes_.end()) {
      it->second.emplace_back(i);
    }
    else{
      indexes_[pdg] = {i};
    }
  }
}

float SimpleFinderNew::CalculateChi2Geo(const KFParticleSIMD& mother) {
  float_v chi2 = mother.Chi2() / simd_cast<float_v>(mother.NDF());
  return chi2[0];
}

bool SimpleFinderNew::IsGoodDaughter(const KFPTrack& track, const DaughterCuts& cuts) {
  int id = cuts.GetId();
  values_.chi2_prim[id] = CalculateChiToPrimaryVertex(track, cuts.GetPdgHypo());
  if (values_.chi2_prim[id] < cuts.GetChi2Prim()){ return false; }
  return true;
}

bool SimpleFinderNew::IsGoodPair(const KFPTrack& track1,
                                 const KFPTrack& track2,
                                 const Decay& decay) {
  const auto& daughters = decay.GetDaughters();
  CalculateParamsInPCA(track1, daughters[0].GetPdgHypo(), track2, daughters[1].GetPdgHypo());
  values_.distance[0] = CalculateDistanceBetweenParticles(params_);

  if(values_.distance[0] > decay.GetMother().GetDistance()){ return false; }
  return true;
}

bool SimpleFinderNew::IsGoodMother(const KFParticleSIMD& mother, const MotherCuts& cuts) {
  values_.chi2_geo = CalculateChi2Geo(mother);
  if(values_.chi2_geo > cuts.GetChi2Geo()){ return false; }

  CalculateMotherProperties(mother);
  if(values_.l_over_dl < cuts.GetLdL()){ return false; }

  return true;
}

void SimpleFinderNew::SaveParticle(KFParticleSIMD& particle_simd) {
  KFParticle particle;

  particle_simd.GetKFParticle(particle, 0);
  particle.SetPDG(particle.GetPDG()); //TODO: not needed?

  OutputContainer mother(particle);
  mother.SetSelectionValues(values_);

  output_.emplace_back(mother);

  tracks_.Resize(tracks_.Size()+1);
  SetTrack(particle, tracks_.Size(), tracks_);
}
