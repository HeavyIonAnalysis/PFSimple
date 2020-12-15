#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Matching.hpp>
#include "AnalysisTree/FillTask.hpp"

#include "Interface/InputContainer.h"

class ConverterIn : public AnalysisTree::FillTask {

  enum eInputBranches {
    kRecEventHeader = 0,
    kSimEventHeader,
    kKfpfTracks,
    kSimTracks,
    kNumberOfInputBranches
  };

 public:
  explicit ConverterIn(const CutsContainer& cuts) : cuts_(cuts)
  {  // Use default names for the input branches
    in_branches_.resize(kNumberOfInputBranches);
    in_branches_[kRecEventHeader] = "RecEventHeader";
    in_branches_[kSimEventHeader] = "SimEventHeader";
    in_branches_[kKfpfTracks] = "VtxTracks";
    in_branches_[kSimTracks] = "SimParticles";
  }
  ~ConverterIn() override = default;

  void Init(std::map<std::string, void*>& branches) override;

  void Exec() override {};

  void Finish() override {};

  const CutsContainer& GetCuts() const { return cuts_; }
//  const InputContainer& GetInputContainer() const { return input_container_; }

  void SetIsShine(bool is=true) { is_shine_ = is; }

  void SetCuts(const CutsContainer& cuts) { cuts_ = cuts; };
  void SetTrackCuts(AnalysisTree::Cuts* const cuts) { track_cuts_ = cuts; };

  InputContainer CreateInputContainer() const;
 protected:

  std::vector<float> GetCovMatrixCbm(const AnalysisTree::Track&) const;
  std::vector<float> GetCovMatrixShine(const AnalysisTree::Track&) const;
  void FillParticle(const AnalysisTree::Track&, InputContainer&) const;
  bool IsGoodTrack(const AnalysisTree::Track& rec_track) const;

  AnalysisTree::TrackDetector* kf_tracks_{nullptr};
  AnalysisTree::Particles* sim_tracks_{nullptr};
  AnalysisTree::Matching* kf2sim_tracks_{nullptr};
  AnalysisTree::EventHeader* rec_event_header_{nullptr};
  AnalysisTree::EventHeader* sim_event_header_{nullptr};
  AnalysisTree::Cuts* track_cuts_{nullptr};

  CutsContainer cuts_;
//  InputContainer input_container_;

  // field ids for input parameters
  int q_field_id_{AnalysisTree::UndefValueInt};
  int pdg_field_id_{AnalysisTree::UndefValueInt};
  int par_field_id_ {AnalysisTree::UndefValueInt};
  int mf_field_id_ {AnalysisTree::UndefValueInt};
  int cov_field_id_{AnalysisTree::UndefValueInt};
  int mother_id_field_id_{AnalysisTree::UndefValueInt};
  int sim_pdg_field_id_{AnalysisTree::UndefValueInt};
  int passcuts_field_id_{AnalysisTree::UndefValueInt};
  int nhits_field_id_{AnalysisTree::UndefValueInt};

  bool is_shine_{false};

};

#endif //KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
