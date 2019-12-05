/**
 ** @class CutsContainer
 ** @brief Container with values of cuts.
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov
 **
 ** The meaning of quantities to be cut is described in the OutputContainer.h
 **/


#ifndef CutsContainer_H
#define CutsContainer_H

class CutsContainer
{
  
 public:
   
  CutsContainer() = default;
  virtual ~CutsContainer() = default;  
  
  //  lambda candidate parameters setters
  void SetCutChi2PrimPos(float value){cut_chi2_prim_pos_ = value;};
  void SetCutChi2PrimNeg(float value){cut_chi2_prim_neg_ = value;};
  void SetCutDistance(float value){cut_distance_ = value;};
  void SetCutCosineDaughterPos(float value){cut_cosine_daughter_pos_ = value;};
  void SetCutCosineDaughterNeg(float value){cut_cosine_daughter_neg_ = value;};
  void SetCutChi2Geo(float value){cut_chi2_geo_ = value;};
  void SetCutLUp(float value){cut_l_up_ = value;};
  void SetCutLDown(float value){cut_l_down_ = value;};
  void SetCutLdL(float value){cut_ldl_ = value;};
  void SetCutIsFromPV(int value){cut_is_from_pv_ = value;};
  void SetCutCosineTopo(float value){cut_cosine_topo_ = value;};
  void SetCutSigmaMassRatio(float value){cut_sigma_mass_ratio_ = value;};
  void SetCutChi2Topo(float value){cut_chi2_topo_ = value;};
  
  void  CancelCuts();     ///< Sets cuts very large (very small) in order to select all the candidates (reject none of them)
  
  //  lambda candidate parameters getters
  float GetCutChi2PrimPos() const {return cut_chi2_prim_pos_;};
  float GetCutChi2PrimNeg() const {return cut_chi2_prim_neg_;};
  float GetCutDistance() const {return cut_distance_;};
  float GetCutCosineDaughterPos() const {return cut_cosine_daughter_pos_;};
  float GetCutCosineDaughterNeg() const {return cut_cosine_daughter_neg_;};
  float GetCutChi2Geo() const {return cut_chi2_geo_;};
  float GetCutLUp() const {return cut_l_up_;};
  float GetCutLDown() const {return cut_l_down_;};
  float GetCutLdL() const {return cut_ldl_;};
  int   GetCutIsFromPV() const {return cut_is_from_pv_;};
  float GetCutCosineTopo() const {return cut_cosine_topo_;};
  float GetCutSigmaMassRatio() const {return cut_sigma_mass_ratio_;};
  float GetCutChi2Topo() const {return cut_chi2_topo_;};

 protected:
   
  // Cuts with their default values
  float cut_chi2_prim_pos_{18.4207};
  float cut_chi2_prim_neg_{18.4207};
  float cut_distance_{1.};
  float cut_cosine_daughter_pos_{0.};
  float cut_cosine_daughter_neg_{0.};
  float cut_chi2_geo_{3.};
  float cut_l_up_{200.}; float cut_l_down_{-5.};
  int   cut_is_from_pv_{0};
  float cut_cosine_topo_{0.};
  float cut_ldl_{5.}; float cut_ldl_sec_{10.};
  float cut_sigma_mass_ratio_{3.};
  float cut_chi2_topo_{5.};   

};
#endif//CutsContainer_H