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

#ifndef SimpleFinder_H
#define SimpleFinder_H

#include "Constants.hpp"
#include "InputContainer.hpp"
#include "OutputContainer.hpp"

#include "KFPTrack.h"
#include "KFPTrackVector.h"
#include "KFParticleSIMD.h"
#include "KFVertex.h"

class SimpleFinder {

 public:
  SimpleFinder() = default;
  virtual ~SimpleFinder() = default;

  void Init(KFPTrackVector&& tracks, const KFVertex& pv);///< Initialize SimpleFinder object with PV and set of tracks of the current event
  void Init(const InputContainer& input);

  void SortTracks();
  void FindParticles();

  const KFPTrackVector* GetTracks() const { return &tracks_; };

  const std::vector<float>& GetMass() const { return vec_mass_; };// TODO remove after debug procedure

  const std::vector<OutputContainer>& GetMotherCandidates() const { return vec_mother_; };

  //  void SetCuts(const CutsContainer& cuts) { cuts_ = cuts; }
  //  void SetDecay(const DecayContainer& decay) { decay_ = decay; }

 protected:
  float CalculateChiToPrimaryVertex(const KFPTrack& track, int pid) const;///< Calculates \f$\chi^2\f$ of the track to the primary vertex (PV)
  static void CalculateParamsInPCA(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2, std::array<float, 8>& pars1, std::array<float, 8>& pars2);
  ///< Recalculates daughters tracks' parameters in the point of their closest approach

  static float CalculateDistanceBetweenParticles(const std::array<float, 8>& pars1, const std::array<float, 8>& pars2);///< Calculates the distance between daughter tracks in their closest approach
  static float CalculateCosMomentumSum(const std::array<float, 8>& pars1, const std::array<float, 8>& pars2);          ///< Calculates the cosine of the angle between daughter's and mother's momenta
  static KFParticleSIMD ConstructMother(const KFPTrack& track1, int pid1, const KFPTrack& track2, int pid2);           ///< Creates mother particle as the KFParticleSIMD object
  static float CalculateChi2Geo(const KFParticleSIMD& mother);                                                         ///< Calculates \f$\chi^2\f$ of daughters' tracks in their closest approach
  void CalculateMotherProperties(const KFParticleSIMD& mother, float& l, float& ldl, int& isFromPV) const;
  ///< Calculates distance between primary and secondary vertices with error and determines whether mother comes from the PV

  float CalculateCosTopo(const KFParticleSIMD& mother) const; ///< Calculates cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV
  float CalculateChi2Topo(const KFParticleSIMD& mother) const;///< Calculates \f$\chi^2\f$ of the mother's track to the PV
  void SaveParticle(const OutputContainer& Lambda);           ///< Saves selected particle with set of geometrical decay parameters

  static std::array<float_v, 3> CalculateCoordinatesSecondaryVertex(const std::array<float, 8>& pars1, const std::array<float, 8>& pars2);   //Calculates the coordinates of the secondary vertex for the first two daugthers
  static void CalculateParamsInSecondaryVertex(const KFParticleSIMD& particleSIMD1, std::array<float_v, 3> xyz, std::array<float, 8>& pars1);//Recalculates 3rd daugthers track parameters in secondary vertex
  static float CalculateDistanceToSecondaryVertex(const std::array<float, 8>& pars1, std::array<float_v, 3>& sec_vx);                        //Calculates Distance of third daughter from secondary vertex
  static float CalculateCosMomentumSumThird(const std::array<float, 8>& pars1, const std::array<float, 8>& pars2, const std::array<float, 8>& pars3);
  static KFParticleSIMD ConstructMotherThree(KFParticleSIMD& particleSIMD1, KFParticleSIMD& particleSIMD2, KFParticleSIMD& particleSIMD3, std::array<float_v, 3> sec_vx);///< Creates mother particle as KFParticleSIMD object for three daughters

  bool IsGoodDaughter(const KFPTrack& track, int pid);

  KFPTrackVector tracks_;
  KFVertex prim_vx_;

  std::array<std::vector<int>, kNumberOfTrackTypes> trIndex_;

  float mass_{0.};               // TODO remove after debug procedure
  std::vector<float> vec_mass_{};// TODO remove after debug procedure

  std::vector<OutputContainer> vec_mother_{};
};

#endif//SimpleFinder_H
