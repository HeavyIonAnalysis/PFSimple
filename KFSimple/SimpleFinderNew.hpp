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

#include "InputContainer.h"
#include "OutputContainer.h"

#include "Constants.hpp"
#include "Decay.hpp"

class SimpleFinderNew{

 public:

  SimpleFinderNew() = default;
  SimpleFinderNew(const SimpleFinderNew&) = default;
  SimpleFinderNew(SimpleFinderNew&&) = default;
  SimpleFinderNew& operator=(SimpleFinderNew&&) = default;
  SimpleFinderNew& operator=(const SimpleFinderNew&) = default;
  ~SimpleFinderNew() = default;

  void Init(KFPTrackVector&& tracks, const KFVertex& pv);///< Initialize SimpleFinder object with PV and set of tracks of the current event
  void Init(const InputContainer& input);

  bool IsGoodDaughter(const KFPTrack& track, const DaughterCuts& cuts){
    return false;
  }

  void ReconstructDecay(const Decay& decay) {

    DaughterCuts cut_1;
    DaughterCuts cut_2;

    std::vector<size_t> indexes_1{};
    std::vector<size_t> indexes_2{};

    for (auto index_1 : indexes_1){

      KFPTrack track_1;
      tracks_.GetTrack(track_1, index_1);
      int pdg_1 = tracks_.PDG()[index_1];

      if(!IsGoodDaughter(track_1, cut_1)) continue;

      for (auto index_2 : indexes_2){

        KFPTrack track_2;
        tracks_.GetTrack(track_2, index_2);
        int pdg_2 = tracks_.PDG()[index_2];

        if(!IsGoodDaughter(track_2, cut_2)) continue;
      }
    }

  }

 private:
  KFPTrackVector tracks_;
  KFVertex prim_vx_;

  std::map<Pdg_t, std::vector<size_t>> indexes_{};

  std::vector <Decay> decays_{};

  std::vector<OutputContainer> output_{};




};

#endif //KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
