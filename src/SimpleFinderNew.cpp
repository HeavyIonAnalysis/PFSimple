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
  const float dx = parameters.first.at(kX) - parameters.second.at(kX);
  const float dy = parameters.first.at(kY) - parameters.second.at(kY);
  const float dz = parameters.first.at(kZ) - parameters.second.at(kZ);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

SimpleFinderNew::Parameters_t SimpleFinderNew::CalculateParamsInPCA(const KFPTrack& track1, int pid1,
                                                                    const KFPTrack& track2, int pid2) {
  SimpleFinderNew::Parameters_t pars;

  KFParticle particle1(track1, pid1);
  KFParticle particle2(track2, pid2);
  KFParticleSIMD particleSIMD1(particle1);// the same particle is copied to each SIMD element
  KFParticleSIMD particleSIMD2(particle2);

  float_v dS[2];
  float_v params1[8], params2[8];
  particleSIMD1.GetDStoParticleFast(particleSIMD2, dS);
  particleSIMD1.TransportFast(dS[0], params1);
  particleSIMD2.TransportFast(dS[1], params2);

//  float_v parbuf;
  for (int i = 0; i < 8; i++) {
//    parbuf = params1[i];
    pars.first.at(i) = params1[i][0];
//    parbuf = params2[i];
    pars.second.at(i) = params2[i][0];
  }
  return pars;
}

KFParticleSIMD SimpleFinderNew::ConstructMother(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs) {
  KFParticleSIMD mother;

  std::vector<KFParticle> particles{};
  std::vector<KFParticleSIMD> particles_simd{};

  for(size_t i=0; i<tracks.size(); ++i){
    particles.emplace_back( KFParticle(tracks.at(i), pdgs.at(i)) );
    particles.at(i).SetId(tracks.at(i).Id());
    particles_simd.emplace_back(particles.at(i));
  }

  if(tracks.size() == 2){
    float_v ds[2] = {0.f, 0.f};
    float_v dsdr[4][6];
    particles_simd.at(0).GetDStoParticle(particles_simd.at(1), ds, dsdr);
    particles_simd.at(0).TransportToDS(ds[0], dsdr[0]);
    particles_simd.at(1).TransportToDS(ds[1], dsdr[3]);
    const KFParticleSIMD* vDaughtersPointer[2] = {&particles_simd.at(0), &particles_simd.at(1)};
    mother.Construct(vDaughtersPointer, 2, nullptr);
  }
  else if(tracks.size() == 3){

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

std::vector<size_t> SimpleFinderNew::GetIndexes(const DaughterCuts& cuts) {
  std::vector<size_t> result{};
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
      indexes_[pdg] = {static_cast<unsigned long>(i)};
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
                                 Pdg_t pdg1,
                                 const KFPTrack& track2,
                                 Pdg_t pdg2,
                                 const Decay& decay) {
  Parameters_t parameters = CalculateParamsInPCA(track1, pdg1, track2, pdg2);
  values_.distance = CalculateDistanceBetweenParticles(parameters);

  if(values_.distance > decay.GetMother().GetDistance()){ return false; }

  return true;
}

bool SimpleFinderNew::IsGoodMother(const KFParticleSIMD& mother, const MotherCuts& cuts) {
  values_.chi2_geo = CalculateChi2Geo(mother);

  if(values_.chi2_geo > cuts.GetChi2Geo()){ return false; }

  CalculateMotherProperties(mother);
  if(values_.l_over_dl < cuts.GetLdL()){ return false; }

  return true;
}

std::array<float_v, 3> SimpleFinderNew::CalculateCoordinatesSecondaryVertex(const SimpleFinderNew::Parameters_t& pars) {
  std::array<float_v, 3> sv;
  //** Calculate coordinates of the secondary vertex of the first two daughters
  sv.at(0) = (pars.first.at(kX) + pars.second.at(kX)) / 2;
  sv.at(1) = (pars.first.at(kY) + pars.second.at(kY)) / 2;
  sv.at(2) = (pars.first.at(kZ) + pars.second.at(kZ)) / 2;
  return sv;
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
