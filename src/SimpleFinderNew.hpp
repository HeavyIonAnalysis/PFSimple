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
#include <TVector3.h>

#include "InputContainer.hpp"
#include "OutputContainer.hpp"

#include "Constants.hpp"
#include "Decay.hpp"
#include "NonLinearCutBase.hpp"

class SimpleFinderNew{

  typedef std::array<float, 8> Param_t;
  typedef std::vector<Param_t> Parameters_t;

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

  bool IsGoodPair(const KFPTrack& track1, const KFPTrack& track2, const Decay& decay);
  bool IsGoodThree(const KFPTrack& track1, const KFPTrack& track2, const KFPTrack& track3, const Decay& decay){
    return true;
  }

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

  bool IsGoodCos(const KFParticleSIMD& mother, const Parameters_t& daughter_pars, const Decay& decay){
    for(int i=0; i<decay.GetNDaughters(); ++i){
      const auto cut = decay.GetDaughters()[i].GetCos();
      const auto& par = daughter_pars.at(i);
      const float norm = mother.GetP()[0] * std::sqrt(par[kPx]*par[kPx] + par[kPy]*par[kPy] + par[kPz]*par[kPz]);
      values_.cos[i] = (mother.GetPx()[0]*par[kPx] + mother.GetPy()[0]*par[kPy] + mother.GetPz()[0]*par[kPz])/norm;
      if(values_.cos[i] < cut){
        return false;
      }
    }
    return true;
  }

  std::array<float, 3> GetSecondaryVertex(){
    if(params_.size() < 2 || params_.size() > 3){
      throw std::runtime_error("Daughter parameters size is wrong");
    }
    std::array<float, 3> sv{};
    for(int i=0; i<3; ++i){
      sv.at(i) = (params_[0].at(kX+i) + params_[1].at(kX+i)) / 2;
    }
    return sv;
  }

  void SaveParticle(KFParticleSIMD& particle_simd);

  void ReconstructDecay(const Decay& decay) {
    const auto& mother_cuts = decay.GetMother();
    const auto& daughters_cuts = decay.GetDaughters();

    std::vector<std::vector<int>> indexes{};
    std::vector<Pdg_t> pdgs{};
    for(const auto& daughter : daughters_cuts) {
      indexes.emplace_back(GetIndexes(daughter));
      pdgs.emplace_back(daughter.GetPdgHypo());
    }

    for (auto index_1 : indexes.at(0)){
      auto track_1 = GetTrack(index_1);

      for (auto index_2 : indexes.at(1)){
        auto track_2 = GetTrack(index_2);
        if(!IsGoodPair(track_1, track_2, decay)) continue;

        if(decay.GetNDaughters() == 2){
          KFParticleSIMD kf_mother = ConstructMother({track_1, track_2}, pdgs);
          if(!IsGoodMother(kf_mother, mother_cuts)) continue;
          if(!IsGoodCos(kf_mother, params_, decay)) continue;
          FillDaughtersInfo({track_1, track_2}, pdgs);
          SaveParticle(kf_mother);
        }
        else if (decay.GetNDaughters() == 3){
          for (auto index_3 : indexes.at(2)) {
            auto track_3 = GetTrack(index_3);
            if(!IsGoodThree(track_1, track_2, track_3, decay)) continue;

            KFParticleSIMD kf_mother = ConstructMother({track_1, track_2, track_3}, pdgs);
            if(!IsGoodMother(kf_mother, mother_cuts)) continue;
            if(!IsGoodCos(kf_mother, params_, decay)) continue;
            FillDaughtersInfo({track_1, track_2, track_3}, pdgs);
            SaveParticle(kf_mother);
          }
        }
        else{
          throw std::runtime_error("Number of daughters should be 2 or 3. Current number is " + std::to_string(decay.GetNDaughters()));
        }
      }
    }
  }

  void AddDecay(const Decay& decay){ decays_.emplace_back(decay); }

  const std::vector<OutputContainer>& GetCandidates() const { return output_; }

 private:

  KFPTrackVector tracks_;  ///< input information: vector of tracks
  KFVertex prim_vx_;  ///< input information: primiry vertex
  std::vector <Decay> decays_{};  ///< input information: list of decays to reconstruct
  NonLinearCutBase* ml_cuts_{nullptr};  ///< input information: non-linear cuts class (optional)

  std::map<Pdg_t, std::vector<int>> indexes_{};  ///< map of indexes for a given particle specie

  Parameters_t params_{};   ///< vector of daughter parameters at current SV estimation
  SelectionValues values_{};  ///< struct with mother and daughters properties used to apply cuts

  std::vector<OutputContainer> output_{};  ///< output information: vector of candidates

  /**
  * Find indexes of good daughters
  * @param cuts daughter particle cuts container
  * @return vector of indexes
  */
  std::vector<int> GetIndexes(const DaughterCuts& cuts);

  static float CalculateDistanceBetweenParticles(const Parameters_t& parameters);
  void CalculateParamsInPCA(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2);
  KFParticleSIMD ConstructMother(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs);
  float CalculateChiToPrimaryVertex(const KFPTrack& track, Pdg_t pid) const;
  static float CalculateChi2Geo(const KFParticleSIMD& mother);
  void CalculateMotherProperties(const KFParticleSIMD& mother);
  static void SetTrack(const KFParticle& particle, int id, KFPTrackVector& tracks);

  void FillDaughtersInfo(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs){
    for(int i=0; i<tracks.size(); ++i){
      values_.chi2_prim[i] = CalculateChiToPrimaryVertex(tracks.at(i), pdgs.at(i));
    }
  }

  float CalculateCosTopo(const KFParticleSIMD& mother) const;

  bool ApplyNonLinearCut() const {
    return ml_cuts_ == nullptr || ml_cuts_->ApplyCut(values_);
  }


};

#endif //KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
