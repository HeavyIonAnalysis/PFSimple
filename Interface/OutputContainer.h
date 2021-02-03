/**
 ** @class OutputContainer
 ** @brief Container with output information about reconstructed particles and geometrical decay parameters (quantities to be cut in order to select particles)
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov, Susanne Glaessel
 **
 ** Each particle candidate is characterized with set of geometrical decay parameters. Depending on the
 ** value of each parameter the candidate is saved or rejected.
 ** In order to save the reconstructed particle, the KFParticle object is used. It contains all
 ** information about the particle (mass, momentum etc), and access to this information is
 ** possible via KFParticle methods.
 **/ 

#ifndef OutputContainer_H
#define OutputContainer_H

#include "KFParticle.h"

class OutputContainer
{
  
 public:
  
  OutputContainer() = default;
  virtual ~OutputContainer() = default;  
  
  //  candidate parameters setters for two daugthers
  void SetChi2PrimPos(float value) {chi2_prim_pos_ = value;};
  void SetChi2PrimNeg(float value) {chi2_prim_neg_ = value;};
  void SetDistance(float value) {distance_ = value;};
  void SetCosineDaughterPos(float value) {cosine_daughter_pos_ = value;};
  void SetCosineDaughterNeg(float value) {cosine_daughter_neg_ = value;};
  void SetChi2Geo(float value) {chi2_geo_ = value;};
  void SetL(float value) {l_ = value;};
  void SetLdL(float value) {ldl_ = value;};
  void SetIsFromPV(int value) {is_from_pv_ = value;};
  void SetCosineTopo(float value) {cosine_topo_ = value;};
  void SetSigmaMassRatio(float value) {sigma_mass_ratio_ = value;};
  void SetChi2Topo(float value) {chi2_topo_ = value;};
  void SetNHitsPos(int value) {nhits_pos_ = value;};
  void SetNHitsNeg(int value) {nhits_neg_ = value;};

  //  candidate parameters setters for third daugther
  void SetChi2PrimThird(float value) {chi2_prim_third_ = value;};
  void SetDistanceThird(float value) {distance_third_ = value;};
  void SetCosineDaughterThird(float value) {cosine_daughter_third_ = value;};
  void SetChi2GeoThree(float value) {chi2_geo_three_ = value;};
  void SetCosineTopoThree(float value) {cosine_topo_three_ = value;};
  void SetChi2TopoThree(float value) {chi2_topo_three_ = value;};
  void SetNHitsThird(int value) {nhits_third_ = value;};
  
  void SetParticle(const KFParticle& particle) {particle_ = particle;};
  
  //  candidate parameters getters for two daugthers
  float GetChi2PrimPos() const {return chi2_prim_pos_;};
  float GetChi2PrimNeg() const {return chi2_prim_neg_;};
  float GetDistance() const {return distance_;};
  float GetCosineDaughterPos() const {return cosine_daughter_pos_;};
  float GetCosineDaughterNeg() const {return cosine_daughter_neg_;};
  float GetChi2Geo() const {return chi2_geo_;};
  float GetL() const {return l_;};
  float GetLdL() const {return ldl_;};
  int   GetIsFromPV() const {return is_from_pv_;};
  float GetCosineTopo() const {return cosine_topo_;};
  float GetSigmaMassRatio() const {return sigma_mass_ratio_;};
  float GetChi2Topo() const {return chi2_topo_;};
  int   GetNHitsPos() const {return nhits_pos_;};
  int   GetNHitsNeg() const {return nhits_neg_;};

  //  candidate parameters getters for third daugther
  float GetChi2PrimThird() const {return chi2_prim_third_;};
  float GetDistanceThird() const {return distance_third_;};
  float GetCosineDaughterThird() const {return cosine_daughter_third_;};
  float GetChi2GeoThree() const {return chi2_geo_three_;};
  float GetCosineTopoThree() const {return cosine_topo_three_;};
  float GetChi2TopoThree() const {return chi2_topo_three_;};
  int   GetNHitsThird() const {return nhits_third_;};

  const KFParticle& GetParticle() const {return particle_;};

 protected:
   
  // candidate selection parameters (to be cut) for two daughters
  float chi2_prim_pos_ {-1.};       ///< \f$\chi^2\f$ of the positive track to the primary vertex (PV)
  float chi2_prim_neg_ {-1.};       ///< \f$\chi^2\f$ of the negative track to the PV
  float distance_ {-1.};            ///< Distance between daughter tracks in their closest approach
  float cosine_daughter_pos_ {-1.}; ///< Cosine of the angle between positive daughter's and mother's momenta
  float cosine_daughter_neg_ {-1.}; ///< Cosine of the angle between negative daughter's and mother's momenta
  float chi2_geo_ {-1.};            ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
  float l_ {-1.};                   ///< Lenght of interpolated track from secondary to primary vertex
  float ldl_ {-1.};                 ///< Distance between primary and secondary vertices divided by error 
  int   is_from_pv_ {-1};           ///< Flag variable whether mother particle comes from the PV (1-yes, 0-no)
  float cosine_topo_{-1.};          ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV
  float sigma_mass_ratio_ {-1.};    ///< Difference between invariant and real mother's mass divided by the error (not used now)
  float chi2_topo_ {-1.};           ///< \f$\chi^2\f$ of the mother's track to the PV
  int   nhits_pos_{-1};
  int   nhits_neg_{-1};

  // candidate selection parameters (to be cut) for third daughter
  float chi2_prim_third_ {-1.};       ///< \f$\chi^2\f$ of the third track to the PV
  float distance_third_ {-1.};        ///< Distance between third daughter track and SV in their closest approach
  float cosine_daughter_third_ {-1.}; ///< Cosine of the angle between third daughter's and mother's momenta
  float chi2_geo_three_ {-1.};        ///< \f$\chi^2\f$ of all three daughters' tracks in their closest approach
  float cosine_topo_three_{-1.};      ///< Cosine of the angle between reconstructed mother's momentum of three particles and mother's radius vector beginning in the PV
  float chi2_topo_three_ {-1.};       ///< \f$\chi^2\f$ of the mother's track of three particles to the PV
  int   nhits_third_{-1};

  KFParticle particle_;
  
};

#endif // OutputContainer_H
