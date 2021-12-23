#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_

#include "Decay.hpp"
#include "OutputContainer.hpp"

#include "AnalysisTree/Cuts.hpp"
#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/EventHeader.hpp"
#include "AnalysisTree/Task.hpp"

class PFSimpleTask;

class ConverterOut : public AnalysisTree::Task {
 public:
  explicit ConverterOut() = default;
  ~ConverterOut() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override {}

  void SetPFSimpleTask(PFSimpleTask* pfsimple_task) { pfsimple_task_ = pfsimple_task; }

  void CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) const;
  void SetDecay(const Decay& decay) { decay_ = decay; }
  void SetOutputCuts(AnalysisTree::Cuts* output_cuts) { output_cuts_ = output_cuts; }
  void SetPidMode(const int pid_mode) { pid_mode_ = pid_mode; }

 protected:
  void InitIndexes();
  void MatchWithMc(AnalysisTree::Particle& particle);
  int GetMothersSimId(AnalysisTree::Particle& lambdarec);
  int DetermineGeneration(int mother_sim_id);

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
  AnalysisTree::Cuts* output_cuts_{nullptr};
  Decay decay_{};
  int pid_mode_{0};

  std::vector<OutputContainer> candidates_;

  PFSimpleTask* pfsimple_task_{nullptr};

  // field ids of simulated events
  int b_field_id_{-1};

  // field ids of input simulated mother
  int mother_id_field_id_{-1};

  int x_field_id_{-1};
  int daughter_id_field_id_{-1};
  int generation_field_id_{-1};
  int g4process_field_id_{-1};
  int g4process_field_id_w_{-1};
  int pt_err_field_id_{-1};

  int chi2prim_field_id_{-1};
  int distance_field_id_{-1};
  int cosine_field_id_{-1};
  int chi2geo_field_id_{-1};

  int chi2geo_sm_field_id_{-1};
  int chi2topo_sm_field_id_{-1};
  int cosine_topo_sm_field_id_{-1};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
