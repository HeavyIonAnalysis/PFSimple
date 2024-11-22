#include "SimpleFinder.hpp"
#include <KFParticleSIMD.h>

void SimpleFinder::Init(std::vector<KFParticle>&& tracks, const KFVertex& pv) {
  tracks_ = tracks;
  prim_vx_ = pv;
  InitTrackIndexesMap();
}

void SimpleFinder::Init(const InputContainer& input) {// TODO know how to de-const without time lose for copying by value
  std::vector<KFParticle> tracks = input.GetTracks();
  Init(std::move(tracks), input.GetVertex());
}

float SimpleFinder::CalculateDistanceBetweenParticles() {
  const float dx = params_.at(0).at(kX) - params_.at(1).at(kX);
  const float dy = params_.at(0).at(kY) - params_.at(1).at(kY);
  const float dz = params_.at(0).at(kZ) - params_.at(1).at(kZ);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void SimpleFinder::CalculateParamsInPCA(const KFParticle& track1, int pid1, const KFParticle& track2, int pid2) {
  params_ = Parameters_t(2);

  KFParticle particle1(track1);
  KFParticle particle2(track2);
  particle1.SetPDG(pid1);
  particle2.SetPDG(pid2);
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

void SimpleFinder::CalculateParamsInSV(const KFParticle& track, int pid) {

  CalculateSecondaryVertex();

  params_.resize(3);

  KFParticle particle(track);
  particle.SetPDG(pid);
  KFParticleSIMD particleSIMD(particle);

  float_v dS;
  float_v dsdr[6];
  float_v params[8];
  float_v sv[3];

  for (int i = 0; i < 3; i++)
    sv[i] = sec_vx_.at(i);

  dS = particleSIMD.GetDStoPoint(sv, dsdr);
  particleSIMD.TransportFast(dS, params);

  for (int i = 0; i < 8; i++)
    params_.at(2).at(i) = params[i][0];
}

float SimpleFinder::CalculateDistanceToSV() const {
  const float dx = params_.at(2).at(0) - sec_vx_.at(0);
  const float dy = params_.at(2).at(1) - sec_vx_.at(1);
  const float dz = params_.at(2).at(2) - sec_vx_.at(2);
  const float dr = sqrt(dx * dx + dy * dy + dz * dz);
  return dr;
}

KFParticleSIMD SimpleFinder::ConstructMother(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs) {
  const auto n = tracks.size();

  KFParticleSIMD mother;
  std::vector<KFParticle> particles = tracks;
  std::vector<KFParticleSIMD> particles_simd{};

  for (size_t i = 0; i < n; ++i) {
    SetKFParticleEnergy(particles.at(i), pdgs.at(i));
    particles.at(i).SetPDG(pdgs.at(i));
    particles.at(i).SetId(tracks.at(i).Id());// TODO rm obsolet ID copying?
    particles_simd.emplace_back(particles.at(i));
  }

  if (n == 2) {
    float_v ds[2] = {0.f, 0.f};
    float_v dsdr[4][6];
    particles_simd.at(0).GetDStoParticle(particles_simd.at(1), ds, dsdr);
    particles_simd.at(0).TransportToDS(ds[0], dsdr[0]);
    particles_simd.at(1).TransportToDS(ds[1], dsdr[3]);
  } else {
    float_v vertex[3] = {sec_vx_[0], sec_vx_[1], sec_vx_[2]};
    for (size_t i = 0; i < n; ++i) {
      float_v ds, dsdr[6];
      ds = particles_simd.at(i).GetDStoPoint(vertex, dsdr);
      particles_simd.at(i).TransportToDS(ds, dsdr);
    }
  }

  if (n == 2) {
    const KFParticleSIMD* vDaughtersPointer[2] = {&particles_simd.at(0), &particles_simd.at(1)};
    mother.Construct(vDaughtersPointer, 2, nullptr);
  } else if (n == 3) {
    const KFParticleSIMD* vDaughtersPointer[3] = {&particles_simd.at(0), &particles_simd.at(1), &particles_simd.at(2)};
    mother.Construct(vDaughtersPointer, 3, nullptr);
  }

  return mother;
}

float SimpleFinder::CalculateChiToPrimaryVertex(const KFParticle& track, Pdg_t pid) const {
  // SIMD'ized version
  KFParticle tmpPart = track;
  tmpPart.SetPDG(pid);
  KFParticleSIMD tmpPartSIMD(tmpPart);

  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);
  const float_v point[3] = {prim_vx_Simd.X(), prim_vx_Simd.Y(), prim_vx_Simd.Z()};
  tmpPartSIMD.TransportToPoint(point);

  float_v chi2vec = tmpPartSIMD.GetDeviationFromVertex(prim_vx_Simd);
  return chi2vec[0];
}

std::vector<int> SimpleFinder::GetGoodDaughtersIndexes(const Daughter& daughter) {
  std::vector<int> result{};
  for (auto pid : daughter.GetPids()) {
    auto it = indexes_.find(pid);
    if (it != indexes_.end()) {
      for (auto i_track : it->second) {
        auto track = GetTrack(i_track);
        if (IsGoodDaughter(track, daughter)) {
          result.emplace_back(i_track);
        }
      }
    }
  }
  return result;
}

void SimpleFinder::InitTrackIndexesMap() {
  for (size_t i = 0; i < tracks_.size(); i++) {
    auto pdg = tracks_.at(i).GetPDG();
    auto it = indexes_.find(pdg);
    if (it != indexes_.end()) {
      it->second.emplace_back(i);
    } else {
      indexes_[pdg] = {static_cast<int>(i)};
    }
  }
}

KFPTrack SimpleFinder::ToKFPTrack(const KFParticle& particle) {
  KFPTrack track;

  track.SetX(particle.GetX());
  track.SetY(particle.GetY());
  track.SetZ(particle.GetZ());
  track.SetPx(particle.GetPx());
  track.SetPy(particle.GetPy());
  track.SetPz(particle.GetPz());

  for (Int_t iC = 0; iC < NumberOfCovElements; iC++) {
    track.SetCovariance(iC, particle.GetCovariance(iC));
  }
  for (Int_t iF = 0; iF < 10; iF++) {
    track.SetFieldCoeff(particle.GetFieldCoeff()[iF], iF);
  }
  track.SetCharge(particle.GetQ());
  track.SetId(particle.Id());

  return track;
}

void SimpleFinder::SetKFParticleEnergy(KFParticle& particle, int pdg) const {
  const KFPTrack& track = ToKFPTrack(particle);
  particle = KFParticle(track, pdg);
}

bool SimpleFinder::IsGoodDaughter(const KFParticle& track, const Daughter& daughter) {
  int id = daughter.GetId();
  values_.chi2_prim.at(id) = CalculateChiToPrimaryVertex(track, daughter.GetPdgHypo());
  if (values_.chi2_prim.at(id) < daughter.GetCutChi2Prim() || std::isnan(values_.chi2_prim.at(id))) { return false; }
  return true;
}

bool SimpleFinder::IsGoodPair(const KFParticle& track1,
                              const KFParticle& track2,
                              const Decay& decay) {
  if (track1.Id() == track2.Id()) { return false; }

  values_.distance = CalculateDistanceBetweenParticles();

  if (values_.distance > decay.GetMother().GetCutDistance() || std::isnan(values_.distance)) { return false; }
  return true;
}

bool SimpleFinder::IsGoodThree(const KFParticle& track1,
                               const KFParticle& track2,
                               const KFParticle& track3,
                               const Decay& decay) {

  if (track1.Id() == track2.Id()) { return false; }
  if (track1.Id() == track3.Id()) { return false; }
  if (track2.Id() == track3.Id()) { return false; }

  values_.distance_sv = CalculateDistanceToSV();
  if (values_.distance_sv > decay.GetMother().GetCutDistanceToSV() || std::isnan(values_.distance_sv)) { return false; }
  return true;
}

bool SimpleFinder::IsGoodMother(const KFParticleSIMD& kf_mother, const Mother& mother, const int id_mother) {
  values_.chi2_geo[id_mother] = kf_mother.Chi2()[0] / simd_cast<float_v>(kf_mother.NDF())[0];
  if (values_.chi2_geo[id_mother] > mother.GetCutChi2Geo()[id_mother] || std::isnan(values_.chi2_geo[id_mother])) { return false; }
  return true;
}

bool SimpleFinder::IsMotherFromPV(const KFParticleSIMD& kf_mother, const Mother& mother, const int id_mother) {
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);
  KFParticleSIMD motherTopo = kf_mother;
  motherTopo.SetProductionVertex(prim_vx_Simd);

  values_.chi2_topo[id_mother] = motherTopo.GetChi2()[0] / simd_cast<float_v>(motherTopo.GetNDF())[0];
  values_.cos_topo[id_mother] = CalculateCosTopo(kf_mother);

  if (id_mother > 0) {
    std::array<bool, 2> is_from_pv = {true, true};
    if (values_.chi2_topo[id_mother] > mother.GetCutChi2Topo()[id_mother] || std::isnan(values_.chi2_topo[id_mother])) is_from_pv.at(0) = false;
    if (values_.cos_topo[id_mother] < mother.GetCutCosTopo()[id_mother] || std::isnan(values_.cos_topo[id_mother])) is_from_pv.at(1) = false;
    if ((is_from_pv.at(0) == false && is_from_pv.at(1) == false) || (is_from_pv.at(0) == false && mother.GetCutCosTopo()[id_mother] == huge_value) || (mother.GetCutChi2Topo()[id_mother] == -huge_value && is_from_pv.at(1) == false))
      return false;
  } else {
    if (values_.chi2_topo[id_mother] > mother.GetCutChi2Topo()[id_mother] || std::isnan(values_.chi2_topo[id_mother])) return false;
    if (values_.cos_topo[id_mother] < mother.GetCutCosTopo()[id_mother] || std::isnan(values_.cos_topo[id_mother])) return false;
  }

  return true;
}

bool SimpleFinder::IsGoodDecayLength(const KFParticleSIMD& kf_mother, const Mother& mother) {

  float_v d_pv_Simd, dd_pv_Simd, l_Simd, dl_Simd;
  float_m isFromPV_Simd;
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);

  kf_mother.GetDistanceToVertexLine(prim_vx_Simd, d_pv_Simd, dd_pv_Simd, &isFromPV_Simd);
  values_.distance_pv = d_pv_Simd[0];
  if (values_.distance_pv < mother.GetCutDistancePVLine() || std::isnan(values_.distance_pv)) { return false; }

  KFParticleSIMD motherTopo = kf_mother;
  motherTopo.SetProductionVertex(prim_vx_Simd);
  motherTopo.KFParticleBaseSIMD::GetDecayLength(l_Simd, dl_Simd);

  values_.l_over_dl = l_Simd[0] / dl_Simd[0];
  if (values_.l_over_dl < mother.GetCutLdL() || std::isnan(values_.l_over_dl)) { return false; }

  values_.l = l_Simd[0];
  if (values_.l < mother.GetCutDecayLength() || std::isnan(values_.l)) { return false; }

  values_.is_from_PV = isFromPV_Simd[0];

  return true;
}

void SimpleFinder::SaveParticle(KFParticleSIMD& particle_simd, const Decay& decay) {
  KFParticle particle;

  particle_simd.GetKFParticle(particle, 0);
  particle.SetPDG(decay.GetMother().GetPdg());

  OutputContainer mother(particle);
  mother.SetSelectionValues(values_);

  output_.emplace_back(mother);
  tracks_.emplace_back(particle);
}

float SimpleFinder::CalculateCosTopo(const KFParticleSIMD& mother) const {
  const float px = mother.GetPx()[0];
  const float py = mother.GetPy()[0];
  const float pz = mother.GetPz()[0];

  const float delta_x = mother.GetX()[0] - prim_vx_.GetX();
  const float delta_y = mother.GetY()[0] - prim_vx_.GetY();
  const float delta_z = mother.GetZ()[0] - prim_vx_.GetZ();

  const float sp = delta_x * px + delta_y * py + delta_z * pz;
  const float norm = std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z) * mother.GetP()[0];
  return sp / norm;
}

bool SimpleFinder::IsGoodCos(const KFParticleSIMD& mother, const Decay& decay) {
  for (int i = 0; i < decay.GetNDaughters(); ++i) {
    const auto cut = decay.GetDaughters()[i].GetCutCos();
    const auto& par = params_.at(i);
    const float norm = mother.GetP()[0] * std::sqrt(par[kPx] * par[kPx] + par[kPy] * par[kPy] + par[kPz] * par[kPz]);
    values_.cos[i] = (mother.GetPx()[0] * par[kPx] + mother.GetPy()[0] * par[kPy] + mother.GetPz()[0] * par[kPz]) / norm;
    if (values_.cos[i] < cut || std::isnan(values_.cos[i])) {
      return false;
    }
  }
  return true;
}

void SimpleFinder::CalculateSecondaryVertex() {
  if (params_.size() != 2) {
    throw std::runtime_error("Daughter parameters size is wrong");
  }
  for (int i = 0; i < 3; ++i)
    sec_vx_.at(i) = (params_[0].at(kX + i) + params_[1].at(kX + i)) / 2;
}

void SimpleFinder::ReconstructDecay(const Decay& decay) {

  std::vector<std::vector<int>> indexes{};
  std::vector<Pdg_t> pdgs{};
  for (const auto& daughter : decay.GetDaughters()) {
    indexes.emplace_back(GetGoodDaughtersIndexes(daughter));
    pdgs.emplace_back(daughter.GetPdgHypo());
  }

  for (auto index_1 : indexes.at(0)) {

    std::array<KFParticle, 3> track;

    track.at(0) = GetTrack(index_1);
    if (std::abs(decay.GetDaughters().at(0).GetPdgHypo()) == 1000020030 || std::abs(decay.GetDaughters().at(0).GetPdgHypo()) == 1000020040) {
      int charge = (int) track.at(0).Q();
      charge *= 2;
      track.at(0).Q() = charge;
      track.at(0).Px() *= std::abs(charge);
      track.at(0).Py() *= std::abs(charge);
      track.at(0).Pz() *= std::abs(charge);
    }

    for (auto index_2 : indexes.at(1)) {
      track.at(1) = GetTrack(index_2);
      if (std::abs(decay.GetDaughters().at(1).GetPdgHypo()) == 1000020030 || std::abs(decay.GetDaughters().at(1).GetPdgHypo()) == 1000020040) {
        int charge = (int) track.at(1).Q();
        charge *= 2;
        track.at(1).Q() = charge;
        track.at(1).Px() *= std::abs(charge);
        track.at(1).Py() *= std::abs(charge);
        track.at(1).Pz() *= std::abs(charge);
      }

      CalculateParamsInPCA(track.at(0), decay.GetDaughters().at(0).GetPdgHypo(), track.at(1), decay.GetDaughters().at(1).GetPdgHypo());
      if (!IsGoodPair(track.at(0), track.at(1), decay)) continue;

      KFParticleSIMD kf_mother = ConstructMother({track.at(0), track.at(1)}, pdgs);
      if (decay.GetNDaughters() == 2) {

        int id_mother = 0;
        if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
        if (!IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;
        if (!IsGoodDecayLength(kf_mother, decay.GetMother())) continue;
        if (!IsGoodCos(kf_mother, decay)) continue;

        FillDaughtersInfo({track.at(0), track.at(1)}, pdgs);
        SaveParticle(kf_mother, decay);

      } else if (decay.GetNDaughters() == 3) {

        int id_mother = 1;

        if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
        if (IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;

        for (auto index_3 : indexes.at(2)) { 
		
	  track.at(2) = GetTrack(index_3);
	  if (std::abs(decay.GetDaughters().at(2).GetPdgHypo()) == 1000020030 || std::abs(decay.GetDaughters().at(2).GetPdgHypo()) == 1000020040) {
	    int charge = (int) track.at(2).Q();
	    charge *= 2;
	    track.at(2).Q() = charge;
	    track.at(2).Px() *= std::abs(charge);
	    track.at(2).Py() *= std::abs(charge);
	    track.at(2).Pz() *= std::abs(charge);
	  }
	  
          CalculateParamsInSV(track.at(2), decay.GetDaughters().at(0).GetPdgHypo());
	  if (!IsGoodThree(track.at(0), track.at(1), track.at(2), decay)) continue;

          kf_mother = ConstructMother({track.at(0), track.at(2)}, {pdgs.at(0), pdgs.at(2)});
          id_mother = 2;
          if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
          if (IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;

          kf_mother = ConstructMother({track.at(1), track.at(2)}, {pdgs.at(1), pdgs.at(2)});
          id_mother = 3;
          if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
          if (IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;

          kf_mother = ConstructMother({track.at(0), track.at(1), track.at(2)}, pdgs);
          id_mother = 0;
          if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
          if (!IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;
          if (!IsGoodDecayLength(kf_mother, decay.GetMother())) continue;
          if (!IsGoodCos(kf_mother, decay)) continue;

          FillDaughtersInfo({track.at(0), track.at(1), track.at(2)}, pdgs);

          SaveParticle(kf_mother, decay);
        }
      } else {
        throw std::runtime_error("Number of daughters should be 2 or 3. Current number is " + std::to_string(decay.GetNDaughters()));
      }
    }
  }
}

void SimpleFinder::FillDaughtersInfo(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs) {
  for (size_t i = 0; i < tracks.size(); ++i) {
    values_.chi2_prim[i] = CalculateChiToPrimaryVertex(tracks.at(i), pdgs.at(i));
  }
}
