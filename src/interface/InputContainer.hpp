/**
 ** @class InputContainer
 ** @brief Container with input information about event (vertex and tracks) and cuts used in the certain event
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov
 **
 ** Each event is characterized with primary vertex, set of tracks and set of cuts
 ** (in general case cuts for different events can be different).\n
 ** Primary vertex is characterized with its coordinates {x, y, z} and optionally with
 ** corresponding covariation matrix.\n
 ** In order to store primary vertex the KFVertex object is used.
 ** Each track is characterized with:\n
 ** x, y, z - three coordinates of point where it is defined;\n
 ** px, py, pz - three components of its momentum in this point;\n
 ** cov - covariation matrix of parameters mentioned above;\n
 ** field[0] - field[8] - magnetic field approximation coefficients along the track's trajectory
 ** field[9] - reference point for MF coefficients;\n
 ** charge - its charge;\n
 ** pdg - PID hypothesis for the track;\n
 ** id - its unique number (conserves through all the algorithm);\n
 ** passcuts - flag variable which shows whether track satisfies pre-selection criteria. By default
 ** passcuts=1 should be used. \n
 ** In order to store tracks the vector of KFParticle objects is used.\n
 ** In order to store cut values used for current event the CutsContainer object is used.
 **/

#ifndef InputContainer_H
#define InputContainer_H

#include "KFParticle.h"
#include "KFVertex.h"

class InputContainer {
 public:
  InputContainer() = default;
  virtual ~InputContainer() = default;

  void SetPV(float x, float y, float z);
  void SetPV(KFVertex vertex);
  void SetPV(KFPVertex vertex);
  void AddTrack(const std::vector<float>& par, const std::vector<float>& cov, const std::vector<float>& field, int charge, int pdg, int id);

  const KFVertex& GetVertex() const { return vtx_; };
  const std::vector<KFParticle>& GetTracks() const { return tracks_; };

  void Clear() {
    if (!tracks_.empty()) {
      tracks_.clear();
    }
  }

  void Reserve(size_t n);

 protected:
  KFVertex vtx_;
  std::vector<KFParticle> tracks_{};
};

#endif//InputContainer_H
