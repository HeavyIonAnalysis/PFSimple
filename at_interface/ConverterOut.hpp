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

  void CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle, int const Ndaughters) const;
  void SetDecays(const std::vector<Decay>& decays) { decays_ = decays; }
  void SetOutputCuts(AnalysisTree::Cuts* output_cuts) { output_cuts_ = output_cuts; }

  void SetSimEventHeaderName(const std::string& name) { sim_events_name_ = name; }
  void SetRecTracksName(const std::string& name) { rec_tracks_name_ = name; }
  void SetSimTracksName(const std::string& name) { mc_particles_name_ = name; }
  void SetIsWriteDetailedBG(bool is=true) { is_detailed_bg_ = is; }

 protected:
  void InitIndexes();
  void MatchWithMc(const OutputContainer& candidate, AnalysisTree::Particle& particle, int const Ndaughters);
  int GetMothersSimId(const OutputContainer& candidate, AnalysisTree::Particle& particlerec, const int Ndaughters) const;
  int DetermineGeneration(int mother_sim_id) const;

  // first returned value is:
  // 1 - reco daughter is unmatched to mc
  // 2 - reco daughter is matched, but primary
  // 3 - reco daughter is secondary, produced not in decay from mother with not expected pdg
  // 4 - reco daughter is secondary, produced not in decay from mother with expected pdg
  // 5 - reco daughter is secondary, produced in decay from mother with not expected pdg
  // 6 - reco daughter is secondary, produced in decay from mother with expected pdg
  // second returned value is mother's sim id (if the first value is >2) or -999 (if mother doesn't exist)
  std::pair<int, int> DetermineDaughtersMCStatus(int daughter_rec_id, Pdg_t mother_expected_pdg, int daughter_gen) const;

  // 1 - daughters have the same mother
  // 2 - daughters have different mothers
  // 0 - at least one daughter does not have mother (e.g. primary)
  static int DetermineMotherMCStatus(int mid1, int mid2) ;

  int DetermineBGType(const OutputContainer& candidate, AnalysisTree::Particle& particle, int const Ndaughters);

  // output branches
  AnalysisTree::EventHeader* events_{nullptr};
  AnalysisTree::Particles* particle_reco_{nullptr};
  AnalysisTree::Particles* particle_sim_{nullptr};
  AnalysisTree::Matching* particle_reco2sim_{nullptr};

  // input branches
  std::string mc_particles_name_{"SimParticles"};
  std::string rec_tracks_name_{"RecParticles"};
  std::string sim_events_name_{"SimEventHeader"};

  AnalysisTree::Particles* mc_particles_{nullptr};
  AnalysisTree::TrackDetector* rec_tracks_{nullptr};
  AnalysisTree::Matching* rec_to_mc_{nullptr};
  AnalysisTree::EventHeader* sim_events_{nullptr};
  AnalysisTree::Cuts* output_cuts_{nullptr};
  std::vector<Decay> decays_{};

  std::vector<OutputContainer> candidates_;

  std::map<int, long unsigned int> rec2sim_daughters_;

  PFSimpleTask* pfsimple_task_{nullptr};

  int Ndaughters_max_{0};
  
  // field ids of simulated events
  int b_field_id_{-1};

  // field ids of input simulated mother
  int mother_id_field_id_{-1};
  //int x_sim_field_id_{-1};

  int x_field_id_{-1};
  int particle_id_field_id_{-1};
  int daughter_id_field_id_{-1};
  int generation_field_id_{-1};
  int g4process_field_id_{-1};
  int g4process_field_id_w_{-1};
  int mc_id_field_id_{-1};

  int pt_err_field_id_{-1};

  int mother_id_field_id_w_{-1};

  int chi2prim_field_id_{-1};
  int distance_field_id_{-1};
  int cosine_field_id_{-1};
  int chi2geo_field_id_{-1};

  int chi2geo_sm_field_id_{-1};
  int cosopen_sm_field_id_{-1};
  int chi2topo_sm_field_id_{-1};
  int cosine_topo_sm_field_id_{-1};

  int chi2prim_mother_field_id_{-1};
  int invmass_discr_field_id_{-1};

  bool is_detailed_bg_{false};
  bool is_write_mc_{true};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
