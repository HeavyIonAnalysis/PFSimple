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

  void InitIndexesMap();

  void FindParticles() {
    for(const auto& decay : decays_){
      ReconstructDecay(decay);
    }
  }

  void ReconstructDecay(const Decay& decay);
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

  bool IsGoodDaughter(const KFPTrack& track, const DaughterCuts& cuts);
  bool IsGoodPair(const KFPTrack& track1, const KFPTrack& track2, const Decay& decay);
  bool IsGoodThree(const KFPTrack& track1, const KFPTrack& track2, const KFPTrack& track3, const Decay& decay){ return true; } //TODO
  bool IsGoodMother(const KFParticleSIMD& mother, const MotherCuts& cuts);
  bool IsGoodCos(const KFParticleSIMD& mother, const Parameters_t& daughter_pars, const Decay& decay);

  KFParticleSIMD ConstructMother(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs);
  std::array<float, 3> GetSecondaryVertex();

  void CalculateParamsInPCA(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2);
  float CalculateChiToPrimaryVertex(const KFPTrack& track, Pdg_t pid) const;
  float CalculateCosTopo(const KFParticleSIMD& mother) const;
  static float CalculateDistanceBetweenParticles(const Parameters_t& parameters);
  static void SetTrack(const KFParticle& particle, int id, KFPTrackVector& tracks);

  void FillDaughtersInfo(const std::vector<KFPTrack>& tracks, const std::vector<Pdg_t>& pdgs);
  void SaveParticle(KFParticleSIMD& particle_simd);

  bool ApplyNonLinearCut() const {
    return ml_cuts_ == nullptr || ml_cuts_->ApplyCut(values_);
  }

  KFPTrack GetTrack(int i) {
    KFPTrack track;
    tracks_.GetTrack(track, i);
    return track;
  }

};

#endif //KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
