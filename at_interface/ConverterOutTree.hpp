#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUTTREE_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUTTREE_H_

#include "Decay.hpp"
#include "OutputContainer.hpp"

#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/EventHeader.hpp"
#include "AnalysisTree/Task.hpp"

class PFSimpleTask;

class ConverterOutTree : public AnalysisTree::Task {
 public:
  explicit ConverterOutTree() = default;
  ~ConverterOutTree() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override;

  void SetPFSimpleTask(PFSimpleTask* pfsimple_task) { pfsimple_task_ = pfsimple_task; }
  void SetOutFilename(const std::string out_file_name) { out_file_name_ = out_file_name; }
  void CopyParticle(const OutputContainer& kf_particle);
  void SetDecay(const Decay& decay) { decay_ = decay; }

 protected:
  void MatchWithMc();
  int GetMothersSimId();
  int DetermineGeneration(int mother_sim_id);

  // output branches
  std::string out_file_name_{""};
  TFile* out_file_{nullptr};
  TTree* out_events_{nullptr};
  TTree* out_reco_{nullptr};
  TTree* out_sim_{nullptr};
  TTree* out_match_{nullptr};
  
  // input branches
  std::string mc_particles_name_{"SimParticles"};
  std::string rec_tracks_name_{"VtxTracks"};
  std::string sim_events_name_{"SimEventHeader"};

  AnalysisTree::Particles* mc_particles_{nullptr};
  AnalysisTree::TrackDetector* rec_tracks_{nullptr};
  AnalysisTree::Matching* rec_to_mc_{nullptr};
  AnalysisTree::EventHeader* sim_events_{nullptr};
  Decay decay_{};

  std::vector<OutputContainer> candidates_;

  PFSimpleTask* pfsimple_task_{nullptr};
  
  // field ids of simulated events
  int b_field_id_{-1};

  // field ids of input simulated mother
  int mother_id_field_id_{-1};
  int g4process_field_id_{-1};

  // output variables
  float b_;
  int id_rec_{-1};
  int pid_;
  float mass_, mass_err_;
  float px_, py_, pz_;
  float pt_err_;
  float x_, y_, z_;
  float x_err_, y_err_, z_err_;
  float phi_err_, eta_err_;

  int daughter_id_1_{-1}, daughter_id_2_{-1}, daughter_id_3_{-1};
  float chi2prim_1_{-1.f}, chi2prim_2_{-1.f}, chi2prim_3_{-1.f};
  float cos_1_{-1.f}, cos_2_{-1.f}, cos_3_{-1.f};  
  float distance_{-1.f}, distance_sv_{-1.f};
  float chi2geo_sm_1_{-1.f}, chi2geo_sm_2_{-1.f}, chi2geo_sm_3_{-1.f};
  float chi2topo_sm_1_{-1.f}, chi2topo_sm_2_{-1.f}, chi2topo_sm_3_{-1.f};
  float costopo_sm_1_{-1.f}, costopo_sm_2_{-1.f}, costopo_sm_3_{-1.f};  
  float chi2geo_{-1.f};
  float chi2topo_{-1.f};
  float costopo_{-1.f};
  float L_{-1.f}, LdL_{-1.f}, distance_pv_line_{-1.f};

  float mass_mc_, px_mc_, py_mc_, pz_mc_;
  int generation_{-1}; 
  int g4process_{-1};
  int is_mc_{0};  
  int id_mc_{-1};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUTTREE_H_
