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
  explicit ConverterIn() {}
  ~ConverterIn() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override{};

  const InputContainer& GetInputContainer() const { return container_; }
  void SetTrackCuts(AnalysisTree::Cuts* const cuts) { track_cuts_ = cuts; };
  void SetAncestorMotherPdgsToBeConsidered(std::vector<int>&& pdgs) { ancestor_pdgs_to_be_considered_ = pdgs; };
  void SetPidMode(int value) { pid_mode_ = value; }/// selection mode for pid:
                                                   /// 0 = no pid
                                                   /// 1 = mc pid
                                                   /// 2 = rec pid provided by Pid framework
                                                   /// 3 = rec pid with max. purity & purity > min. requested purity;
                                                   /// 4 = rec pid with purity > min. requested purity, pdg-specific purity is possible

  void SetRecEventHeaderName(const std::string& name) { rec_event_header_name_ = name; }
  void SetRecTracksName(const std::string& name) { kf_tracks_name_ = name; }
  void SetSimTracksName(const std::string& name) { sim_tracks_name_ = name; }

  void UseNoPID() { pid_mode_ = 0; }
  void UseMcPID() { pid_mode_ = 1; }
  void UseRecPID() { pid_mode_ = 2; }
  void UseRecPIDPurityMax() { pid_mode_ = 3; }
  void UseRecPIDPurityMin() { pid_mode_ = 4; }
  void SetPidPurity(const float pid_purity) {
    if (pid_mode_ != 3 && pid_mode_ != 4)
      throw std::runtime_error("Pdg purity can only be set in pidmode 3 or 4. Set purity to -1 or change pidmode.");
    for (int ipid = 0; ipid < pid_purity_.size(); ipid++)
      pid_purity_.at(ipid) = pid_purity;
  }
  void SetPidPurityProton(const float pid_purity) {
    if (pid_mode_ != 4) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 4. Set purity to -1 or change pidmode.");
    } else
      pid_purity_.at(0) = pid_purity;
  };
  void SetPidPurityPion(const float pid_purity) {
    if (pid_mode_ != 4) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 4. Set purity to -1 or change pidmode.");
    } else
      pid_purity_.at(1) = pid_purity;
  };
  void SetPidPurityKaon(const float pid_purity) {
    if (pid_mode_ != 4) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 4. Set purity to -1 or change pidmode.");
    } else
      pid_purity_.at(2) = pid_purity;
  };
  void SetPidPurityDeuteron(const float pid_purity) {
    if (pid_mode_ != 4) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 4. Set purity to -1 or change pidmode.");
    } else
      pid_purity_.at(3) = pid_purity;
  };
  void SetPidPurityBG(const float pid_purity) {
    if (pid_mode_ != 4) {
      throw std::runtime_error("Pdg-specific purity can only be set in pidmode 4. Set purity to -1 or change pidmode.");
    } else
      pid_purity_.at(4) = pid_purity;
  };

 protected:
  std::vector<float> GetCovMatrixCbm(const AnalysisTree::BranchChannel&) const;
  void FillParticle(const AnalysisTree::BranchChannel&);
  bool IsGoodTrack(const AnalysisTree::BranchChannel& rec_track) const;
  bool CheckAncestorPdgs(const AnalysisTree::BranchChannel& rec_track) const;

  AnalysisTree::Branch kf_tracks_;
  AnalysisTree::Branch sim_tracks_;
  AnalysisTree::Matching* kf2sim_tracks_{nullptr};
  AnalysisTree::Branch rec_event_header_;
  AnalysisTree::Cuts* track_cuts_{nullptr};
  std::vector<int> ancestor_pdgs_to_be_considered_;

  InputContainer container_;

  std::string rec_event_header_name_;
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
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
