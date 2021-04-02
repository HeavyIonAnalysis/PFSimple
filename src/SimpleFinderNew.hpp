/**
 * @class SimpleFinder
 * @brief V0 particle (Lambda) reconstruction algorithm
 * @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov, Susanne Glaessel
 *
 * Simplified version of KFParticleFinder class. At the current moment it is developed to reconstruct 2- and 3-body-decays. \n
 * SimpleFinder is based on KFParticle package, it uses mathematical apparatus implemented in KFParticle.
 * The basic idea and reconstruction algorithm also reproduce the KFParticle, but SimpleFinder is free of overloading of too complicated procedure as in KFParticleFinder. \n
 * Also the advantage of SimpleFinder is that reconstruction procedure is clear and under full control of user, almost in the "hand mode".
 * It gives a possibility of detailed analysis of V0 reconstruction, in particular of decay parameters distributions in order to optomize selection criterias (cuts).\n
 */

#ifndef KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_

#include <vector>
#include <limits>

#include <KFPTrackVector.h>
#include <KFVertex.h>
#include <KFParticleSIMD.h>

#include "InputContainer.hpp"
#include "OutputContainer.hpp"

#include "Constants.hpp"
#include "Decay.hpp"

struct SelectionValues {
  SelectionValues() = default;

  float chi2_prim[3]{-1.f, -1.f, -1.f};
  float distance{std::numeric_limits<float>::max()};
  float l{-1.f};
  float l_over_dl{-1.f};
  float chi2_geo{-1.f};
  float chi2_topo{-1.f};

  bool is_from_PV{false};
};

class SimpleFinderNew{

  typedef std::pair<std::array<float, 8>, std::array<float, 8>> Parameters_t;

 public:

  SimpleFinderNew() = default;
  SimpleFinderNew(const SimpleFinderNew&) = default;
  SimpleFinderNew(SimpleFinderNew&&) = default;
  SimpleFinderNew& operator=(SimpleFinderNew&&) = default;
  SimpleFinderNew& operator=(const SimpleFinderNew&) = default;
  ~SimpleFinderNew() = default;

  void Init(KFPTrackVector&& tracks, const KFVertex& pv); ///< Initialize SimpleFinder object with PV and set of tracks of the current event
  void Init(const InputContainer& input);

  bool IsGoodDaughter(const KFPTrack& track, const DaughterCuts& cuts, int id=0) {
    values_.chi2_prim[id] = CalculateChiToPrimaryVertex(track, cuts.GetPdgHypo());
    if (values_.chi2_prim[id] < cuts.GetChi2Prim()){ return false; }
    return true;
  }

  bool IsGoodPair(const KFPTrack& track1, Pdg_t pdg1, const KFPTrack& track2, Pdg_t pdg2, const Decay& decay) {
    Parameters_t parameters = CalculateParamsInPCA(track1, pdg1, track2, pdg2);
    values_.distance = CalculateDistanceBetweenParticles(parameters);

    if(values_.distance > decay.GetMother().GetDistance()){ return false; }

    return true;
  }

  bool IsGoodMother(const KFParticleSIMD& mother, const MotherCuts& cuts) {
    values_.chi2_geo = CalculateChi2Geo(mother);

    if(values_.chi2_geo > cuts.GetChi2Geo()){ return false; }

    CalculateMotherProperties(mother);
    if(values_.l_over_dl < cuts.GetLdL()){ return false; }

    return true;
  }

  void InitIndexesMap();

  void FindParticles() {
    for(const auto& decay : decays_){
      ReconstructDecay(decay);
    }
  }


  std::array<float_v, 3> CalculateCoordinatesSecondaryVertex(const Parameters_t& pars) {
    std::array<float_v, 3> sv;
    //** Calculate coordinates of the secondary vertex of the first two daughters
    sv.at(0) = (pars.first.at(kX) + pars.second.at(kX)) / 2;
    sv.at(1) = (pars.first.at(kY) + pars.second.at(kY)) / 2;
    sv.at(2) = (pars.first.at(kZ) + pars.second.at(kZ)) / 2;
    return sv;
  }

  void SaveParticle(KFParticleSIMD& particle_simd){
    KFParticle particle;

    particle_simd.GetKFParticle(particle, 0);
    particle.SetPDG(particle.GetPDG()); //TODO: not needed?

    OutputContainer mother(particle);
    mother.SetL(values_.l);
    mother.SetLdL(values_.l_over_dl);
    mother.SetIsFromPV(values_.is_from_PV);
    mother.SetChi2Topo(values_.chi2_topo);
    mother.SetChi2Geo(values_.chi2_geo);

    output_.emplace_back(mother);
  }

  void ReconstructDecay(const Decay& decay) {

    const auto& mother_cuts = decay.GetMother();

    const auto& cuts = decay.GetDaughters();

    std::vector<std::vector<size_t>> indexes{};
    std::vector<Pdg_t> pdgs{};
    for(const auto& cut : cuts) {
      indexes.emplace_back(GetIndexes(cut.GetPids()));
      pdgs.emplace_back(cut.GetPdgHypo());
    }
    for (auto index_1 : indexes.at(0)){
      KFPTrack track_1;
      tracks_.GetTrack(track_1, index_1);
      if(!IsGoodDaughter(track_1, cuts.at(0), 0)) continue;

      for (auto index_2 : indexes.at(1)){
        KFPTrack track_2;
        tracks_.GetTrack(track_2, index_2);
        if(!IsGoodDaughter(track_2, cuts.at(1), 1)) continue;
        if(!IsGoodPair(track_1, pdgs.at(0), track_2, pdgs.at(1), decay)) continue;

        if(decay.GetNDaughters() == 2){
          KFParticleSIMD kf_mother = ConstructMother({track_1, track_2}, pdgs);
          if(!IsGoodMother(kf_mother, mother_cuts)) continue;

          SaveParticle(kf_mother);
        }
        else if (decay.GetNDaughters() == 3){
          for (auto index_3 : indexes.at(2)) {
            KFPTrack track_3;
            tracks_.GetTrack(track_3, index_3);
            if(!IsGoodDaughter(track_3, cuts.at(2), 2)) continue;
//            KFParticleSIMD kf_mother = ConstructMother(particleSIMDPos, particleSIMDNeg, particleSIMDThird,sec_vx);
//            SaveParticle(kf_mother);
          }
        }
        else{
          throw std::runtime_error("Number of daughters should be 2 or 3. Current number is " + std::to_string(decay.GetNDaughters()));
        }
      }
    }
  }

  void AddDecay(const Decay& decay){ decays_.emplace_back(decay); }


  std::vector<size_t> GetIndexes(const std::vector<Pdg_t>& pids);
  static float CalculateDistanceBetweenParticles(const Parameters_t& parameters);
  static Parameters_t CalculateParamsInPCA(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2);
  static KFParticleSIMD ConstructMother(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs);
  float CalculateChiToPrimaryVertex(const KFPTrack& track, Pdg_t pid) const;
  static float CalculateChi2Geo(const KFParticleSIMD& mother);
  void CalculateMotherProperties(const KFParticleSIMD& mother);

 private:
  KFPTrackVector tracks_;
  KFVertex prim_vx_;

  std::map<Pdg_t, std::vector<size_t>> indexes_{};

  std::vector <Decay> decays_{};

  SelectionValues values_{};

  std::vector<OutputContainer> output_{};

};

#endif //KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
