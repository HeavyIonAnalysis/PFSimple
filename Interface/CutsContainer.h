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
  void SetCutChi2PrimPos(float value){cut_chi2_prim_pos_ = value; is_apply_cut_chi2_prim_pos_ = true;};
  void SetCutChi2PrimNeg(float value){cut_chi2_prim_neg_ = value; is_apply_cut_chi2_prim_neg_ = true;};
  void SetCutDistance(float value){cut_distance_ = value; is_apply_cut_distance_ = true;};
  void SetCutCosineDaughterPos(float value){cut_cosine_daughter_pos_ = value; is_apply_cut_cosine_daughter_pos_ = true;};
  void SetCutCosineDaughterNeg(float value){cut_cosine_daughter_neg_ = value; is_apply_cut_cosine_daughter_neg_ = true;};
  void SetCutChi2Geo(float value){cut_chi2_geo_ = value; is_apply_cut_chi2_geo_ = true;};
  void SetCutLUp(float value){cut_l_up_ = value; is_apply_cut_l_up_ = true;};
  void SetCutLDown(float value){cut_l_down_ = value; is_apply_cut_l_down_ = true;};
  void SetCutLdL(float value){cut_ldl_ = value; is_apply_cut_ldl_ = true;};
  void SetCutIsFromPV(int value){cut_is_from_pv_ = value; is_apply_cut_is_from_pv_ = true;};
  void SetCutCosineTopo(float value){cut_cosine_topo_ = value; is_apply_cut_cosine_topo_ = true;};
  void SetCutSigmaMassRatio(float value){cut_sigma_mass_ratio_ = value; is_apply_cut_sigma_mass_ratio_ = true;};
  void SetCutChi2Topo(float value){cut_chi2_topo_ = value; is_apply_cut_chi2_topo_ = true;};
  
  void CancelCutChi2PrimPos(){is_apply_cut_chi2_prim_pos_ = false;};
  void CancelCutChi2PrimNeg(){is_apply_cut_chi2_prim_neg_ = false;};
  void CancelCutDistance(){is_apply_cut_distance_ = false;};
  void CancelCutCosineDaughterPos(){is_apply_cut_cosine_daughter_pos_ = false;};
  void CancelCutCosineDaughterNeg(){is_apply_cut_cosine_daughter_neg_ = false;};
  void CancelCutChi2Geo(){is_apply_cut_chi2_geo_ = false;};
  void CancelCutLUp(){is_apply_cut_l_up_ = false;};
  void CancelCutLDown(){is_apply_cut_l_down_ = false;};
  void CancelCutLdL(){is_apply_cut_ldl_ = false;};
  void CancelCutIsFromPV(){is_apply_cut_is_from_pv_ = false;};
  void CancelCutCosineTopo(){is_apply_cut_cosine_topo_ = false;};
  void CancelCutSigmaMassRatio(){is_apply_cut_sigma_mass_ratio_ = false;};
  void CancelCutChi2Topo(){is_apply_cut_chi2_topo_ = false;};
  
  void CancelCuts();
  
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
  
  bool GetIsApplyCutChi2PrimPos() const {return is_apply_cut_chi2_prim_pos_;};
  bool GetIsApplyCutChi2PrimNeg() const {return is_apply_cut_chi2_prim_neg_;};
  bool GetIsApplyCutDistance() const {return is_apply_cut_distance_;};
  bool GetIsApplyCutCosineDaughterPos() const {return is_apply_cut_cosine_daughter_pos_;};
  bool GetIsApplyCutCosineDaughterNeg() const {return is_apply_cut_cosine_daughter_neg_;};
  bool GetIsApplyCutChi2Geo() const {return is_apply_cut_chi2_geo_;};
  bool GetIsApplyCutLUp() const {return is_apply_cut_l_up_;};
  bool GetIsApplyCutLDown() const {return is_apply_cut_l_down_;};
  bool GetIsApplyCutLdL() const {return is_apply_cut_ldl_;};
  bool GetIsApplyCutIsFromPV() const {return is_apply_cut_is_from_pv_;};
  bool GetIsApplyCutCosineTopo() const {return is_apply_cut_cosine_topo_;};
  bool GetIsApplyCutSigmaMassRatio() const {return is_apply_cut_sigma_mass_ratio_;};
  bool GetIsApplyCutChi2Topo() const {return is_apply_cut_chi2_topo_;};
  
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
  float cut_ldl_{5.};
  float cut_sigma_mass_ratio_{3.};
  float cut_chi2_topo_{5.};   
  
  bool is_apply_cut_chi2_prim_pos_{true};
  bool is_apply_cut_chi2_prim_neg_{true};
  bool is_apply_cut_distance_{true};
  bool is_apply_cut_cosine_daughter_pos_{true};
  bool is_apply_cut_cosine_daughter_neg_{true};
  bool is_apply_cut_chi2_geo_{true};
  bool is_apply_cut_l_up_{true}; bool is_apply_cut_l_down_{true};
  bool is_apply_cut_is_from_pv_{true};
  bool is_apply_cut_cosine_topo_{true};
  bool is_apply_cut_ldl_{true};
  bool is_apply_cut_sigma_mass_ratio_{true};
  bool is_apply_cut_chi2_topo_{true}; 
  
};
#endif//CutsContainer_H