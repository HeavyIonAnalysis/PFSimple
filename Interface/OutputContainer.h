/**
 ** @class OutputContainer
 ** @brief Container with output information about reconstructed particles and geometrical decay parameters (quantities to be cut in order to select particles)
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov
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
  
  //  lambda candidate parameters setters
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
  
  void SetParticle(KFParticle particle) {particle_ = particle;};
  
  //  lambda candidate parameters getters
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
             
  const KFParticle& GetParticle() const {return particle_;};

 protected:
   
  //  lambda candidate selection parameters (to be cut)
  float chi2_prim_pos_ {-1.};       ///< \f$\chi^2\f$ of the positive track to the primary vertex (PV)
  float chi2_prim_neg_ {-1.};       ///< \f$\chi^2\f$ of the negative track to the PV
  float distance_ {-1.};            ///< Distance between daughter tracks in their closest approach
  float cosine_daughter_pos_ {-1.}; ///< Cosine of the angle between positive daughter's and mother's momenta
  float cosine_daughter_neg_ {-1.}; ///< Cosine of the angle between negative daughter's and mother's momenta
  float chi2_geo_ {-1.};            ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
  float l_ {-1.};                   ///< Distance between primary and secondary vertices
  float ldl_ {-1.};                 ///< Distance between primary and secondary vertices divided by error 
  int   is_from_pv_ {-1};           ///< Flag variable whether mother particle comes from the PV (1-yes, 0-no)
  float cosine_topo_{-1.};          ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV
  float sigma_mass_ratio_ {-1.};    ///< Difference between invariant and real mother's mass divided by the error (not used now)
  float chi2_topo_ {-1.};           ///< \f$\chi^2\f$ of the mother's track to the PV

  KFParticle particle_;
  
};

#endif // OutputContainer_H