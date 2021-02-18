#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_

#include "AnalysisTree/Task.hpp"
#include <AnalysisTree/Detector.hpp>

#include "DecayContainer.h"
#include "InputContainer.h"

namespace AnalysisTree{
  class EventHeader;
  class Matching;
  class Cuts;
}

class ConverterIn : public AnalysisTree::Task {

 public:
  ConverterIn(const DecayContainer& decay, const CutsContainer& cuts) : decay_(decay), cuts_(cuts) {
    rec_event_header_name_ = "RecEventHeader";
    sim_event_header_name_ = "SimEventHeader";
    kf_tracks_name_ = "VtxTracks";
    sim_tracks_name_ = "SimParticles";
  }
  ~ConverterIn() override = default;

  void Init() override;
  void Exec() override{};
  void Finish() override{};

  const DecayContainer& GetDecay() const { return decay_; }
  const CutsContainer& GetCuts() const { return cuts_; }
  //  const InputContainer& GetInputContainer() const { return input_container_; }

  void SetIsShine(bool is = true) { is_shine_ = is; }

  void SetDecay(const DecayContainer& decay) { decay_ = decay; };
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

  DecayContainer decay_;
  CutsContainer cuts_;

  std::string rec_event_header_name_;
  std::string sim_event_header_name_;
  std::string kf_tracks_name_;
  std::string sim_tracks_name_;

  // field ids for input parameters
  int q_field_id_{AnalysisTree::UndefValueInt};
  int pdg_field_id_{AnalysisTree::UndefValueInt};
  int par_field_id_{AnalysisTree::UndefValueInt};
  int mf_field_id_{AnalysisTree::UndefValueInt};
  int cov_field_id_{AnalysisTree::UndefValueInt};
  int mother_id_field_id_{AnalysisTree::UndefValueInt};
  int sim_pdg_field_id_{AnalysisTree::UndefValueInt};
  int passcuts_field_id_{AnalysisTree::UndefValueInt};
  int nhits_field_id_{AnalysisTree::UndefValueInt};

  bool is_shine_{false};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTERIN_H_
