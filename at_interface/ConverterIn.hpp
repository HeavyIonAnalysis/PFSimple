#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_

#include <utility>

#include "AnalysisTree/Task.hpp"
#include <AnalysisTree/Detector.hpp>

#include "Decay.hpp"
#include "InputContainer.hpp"

namespace AnalysisTree {
class EventHeader;
class Matching;
class Cuts;
}// namespace AnalysisTree

class ConverterIn : public AnalysisTree::Task {

 public:
  explicit ConverterIn() {
    rec_event_header_name_ = "RecEventHeader";
    sim_event_header_name_ = "SimEventHeader";
    kf_tracks_name_ = "VtxTracks";
    sim_tracks_name_ = "SimParticles";
  }
  ~ConverterIn() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override{};

  const InputContainer& GetInputContainer() const { return container_; }
  void SetIsShine(bool is = true) { is_shine_ = is; }
  void SetTrackCuts(AnalysisTree::Cuts* const cuts) { track_cuts_ = cuts; };
  void SetMotherPdgsToBeConsidered(std::vector<int>&& pdgs) { mother_pdgs_to_be_considered_ = pdgs; };
  void SetPidMode(int value) { pid_mode_ = value; }/// selection mode for pid:
                                                   /// 0 = no pid
                                                   /// 1 = mc pid
                                                   /// 2 = rec pid with max. purity & purity > min. requested purity;
                                                   /// 3 = rec pid with purity > min. requested purity, pdg-specific purity is possible
  void UseNoPID() { pid_mode_ = 0; }
  void UseMcPID() { pid_mode_ = 1; }
  void UseRecPIDPurityMax() { pid_mode_ = 2; }
  void UseRecPIDPurityMin() { pid_mode_ = 3; }
  void SetPidPurity(const float pid_purity) {
    if (pid_mode_ != 2 && pid_mode_ != 3)
      throw std::runtime_error("Pdg purity can only be set in pidmode 2 or 3");
    for (int ipid = 0; ipid < pid_purity_.size(); ipid++)
      pid_purity_.at(ipid) = pid_purity;
  }
  void SetPidPurityProton(const float pid_purity) {
    if (pid_mode_ != 3) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 3");
    } else
      pid_purity_.at(0) = pid_purity;
  };
  void SetPidPurityPion(const float pid_purity) {
    if (pid_mode_ != 3) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 3");
    } else
      pid_purity_.at(1) = pid_purity;
  };
  void SetPidPurityKaon(const float pid_purity) {
    if (pid_mode_ != 3) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 3");
    } else
      pid_purity_.at(2) = pid_purity;
  };
  void SetPidPurityDeuteron(const float pid_purity) {
    if (pid_mode_ != 3) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 3");
    } else
      pid_purity_.at(3) = pid_purity;
  };
  void SetPidPurityBG(const float pid_purity) {
    if (pid_mode_ != 3) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 3");
    } else
      pid_purity_.at(4) = pid_purity;
  };

 protected:
  std::vector<float> GetCovMatrixCbm(const AnalysisTree::BranchChannel&) const;
  std::vector<float> GetCovMatrixShine(const AnalysisTree::BranchChannel&) const;
  void FillParticle(const AnalysisTree::BranchChannel&);
  bool IsGoodTrack(const AnalysisTree::BranchChannel& rec_track) const;
  bool CheckMotherPdgs(const AnalysisTree::BranchChannel& rec_track) const;

  AnalysisTree::Branch kf_tracks_;
  AnalysisTree::Branch sim_tracks_;
  AnalysisTree::Matching* kf2sim_tracks_{nullptr};
  AnalysisTree::Branch rec_event_header_;
  AnalysisTree::Branch sim_event_header_;
  AnalysisTree::Cuts* track_cuts_{nullptr};
  std::vector<int> mother_pdgs_to_be_considered_;

  InputContainer container_;

  std::string rec_event_header_name_;
  std::string sim_event_header_name_;
  std::string kf_tracks_name_;
  std::string sim_tracks_name_;

  // fields for input parameters
  
  std::vector<AnalysisTree::Field> mf_field_;
  
  AnalysisTree::Field x_field_;
  AnalysisTree::Field y_field_;
  AnalysisTree::Field z_field_;
  AnalysisTree::Field px_field_;
  AnalysisTree::Field py_field_;
  AnalysisTree::Field pz_field_;
  
  AnalysisTree::Field q_field_;
  AnalysisTree::Field mc_pdg_field_;
  AnalysisTree::Field nhits_field_;
  
  AnalysisTree::Field rec_pdg_field_;
  AnalysisTree::Field prob_p_field_;
  AnalysisTree::Field prob_pi_field_;
  AnalysisTree::Field prob_K_field_;
  AnalysisTree::Field prob_d_field_;
  AnalysisTree::Field prob_bg_field_;
  
  AnalysisTree::Field vtx_x_field_;
  AnalysisTree::Field vtx_y_field_;
  AnalysisTree::Field vtx_z_field_;
  
  AnalysisTree::Field tx_field_;
  AnalysisTree::Field ty_field_;
  AnalysisTree::Field qp_field_;
  
  std::vector<AnalysisTree::Field> cov_field_;
  
  AnalysisTree::Field mother_id_field_;
  AnalysisTree::Field sim_pdg_field_;
  
  int pid_mode_{1};
  std::array<float, NumberOfPids> pid_purity_{0.5, 0.5, 0.5, 0.5, 0.5};
  bool is_shine_{false};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
