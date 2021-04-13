#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_

#include "OutputContainer.hpp"

#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/EventHeader.hpp"
#include "AnalysisTree/Task.hpp"

class ConverterOut : public AnalysisTree::Task {
 public:
  explicit ConverterOut() = default;
  ~ConverterOut() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override {}

  void SetCandidates(const std::vector<OutputContainer>& candidates) { candidates_ = candidates; }

  void CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) const;

    protected:
  void InitIndexes();
  void MatchWithMc();

  // output branches
  AnalysisTree::EventHeader* events_{nullptr};
  AnalysisTree::Particles* lambda_reco_{nullptr};
  AnalysisTree::Particles* lambda_sim_{nullptr};
  AnalysisTree::Matching* lambda_reco2sim_{nullptr};

  // input branches
  std::string mc_particles_name_{"SimParticles"};
  std::string rec_tracks_name_{"VtxTracks"};
  std::string sim_events_name_{"SimEventHeader"};

  AnalysisTree::Particles* mc_particles_{nullptr};
  AnalysisTree::TrackDetector* rec_tracks_{nullptr};
  AnalysisTree::Matching* rec_to_mc_{nullptr};
  AnalysisTree::EventHeader* sim_events_{nullptr};

  std::vector<OutputContainer> candidates_;

  // field ids of simulated events
  int b_field_id_{-1};

  // field ids of input simulated mother
  int mother_id_field_id_{-1};

  int x_field_id_{-1};
  int daughter_id_field_id_{-1};
  int is_signal_field_id_{-1};
  int pt_err_field_id_{-1};

  int chi2prim_field_id_{-1};
  int distance_field_id_{-1};
  int cosine_field_id_{-1};
  int chi2geo_field_id_{-1};

};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
