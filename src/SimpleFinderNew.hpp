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

#include <KFPTrackVector.h>
#include <KFVertex.h>
#include <KFParticleSIMD.h>

#include "InputContainer.hpp"
#include "OutputContainer.hpp"

#include "Constants.hpp"
#include "Decay.hpp"

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

  bool IsGoodDaughter(const KFPTrack& track, const DaughterCuts& cuts);

  bool IsGoodPair(const KFPTrack& track1, Pdg_t pdg1, const KFPTrack& track2, Pdg_t pdg2, const Decay& decay);

  bool IsGoodMother(const KFParticleSIMD& mother, const MotherCuts& cuts);

  void InitIndexesMap();

  void FindParticles() {
    for(const auto& decay : decays_){
      ReconstructDecay(decay);
    }
  }

  KFPTrack GetTrack(int i) {
    KFPTrack track;
    tracks_.GetTrack(track, i);
    return track;
  }

  std::array<float_v, 3> CalculateCoordinatesSecondaryVertex(const Parameters_t& pars);

  void SaveParticle(KFParticleSIMD& particle_simd);

  void ReconstructDecay(const Decay& decay) {
    const auto& mother_cuts = decay.GetMother();
    const auto& daughters_cuts = decay.GetDaughters();

    std::vector<std::vector<size_t>> indexes{};
    std::vector<Pdg_t> pdgs{};
    for(const auto& daughter : daughters_cuts) {
      indexes.emplace_back(GetIndexes(daughter));
      pdgs.emplace_back(daughter.GetPdgHypo());
    }

    for (auto index_1 : indexes.at(0)){
      auto track_1 = GetTrack(index_1);

      for (auto index_2 : indexes.at(1)){
        auto track_2 = GetTrack(index_2);
        if(!IsGoodPair(track_1, pdgs.at(0), track_2, pdgs.at(1), decay)) continue;

        if(decay.GetNDaughters() == 2){
          KFParticleSIMD kf_mother = ConstructMother({track_1, track_2}, pdgs);
          if(!IsGoodMother(kf_mother, mother_cuts)) continue;
          SaveParticle(kf_mother);
        }
        else if (decay.GetNDaughters() == 3){
          for (auto index_3 : indexes.at(2)) {
            auto track_3 = GetTrack(index_3);
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


  std::vector<size_t> GetIndexes(const DaughterCuts& cuts);
  static float CalculateDistanceBetweenParticles(const Parameters_t& parameters);
  static Parameters_t CalculateParamsInPCA(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2);
  static KFParticleSIMD ConstructMother(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs);
  float CalculateChiToPrimaryVertex(const KFPTrack& track, Pdg_t pid) const;
  static float CalculateChi2Geo(const KFParticleSIMD& mother);
  void CalculateMotherProperties(const KFParticleSIMD& mother);

  const std::vector<OutputContainer>& GetCandidates() const { return output_; }
 private:
  KFPTrackVector tracks_;
  KFVertex prim_vx_;

  std::map<Pdg_t, std::vector<size_t>> indexes_{};

  std::vector <Decay> decays_{};

  SelectionValues values_{};

  std::vector<OutputContainer> output_{};

};

#endif //KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
