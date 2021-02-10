#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_

#include <Interface/DecayContainer.h>
#include <Interface/OutputContainer.h>

#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/EventHeader.hpp"
#include "AnalysisTree/Task.hpp"

class ConverterOut : public AnalysisTree::Task {
 public:
  explicit ConverterOut(const DecayContainer& decay) : decay_(decay) {}
  //ConverterOut() = default;
  ~ConverterOut() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override {}

  void SetCandidates(const std::vector<OutputContainer>& canditates) { canditates_ = canditates; }
  void SetDecay(const DecayContainer& decay) { decay_ = decay; };

 protected:
  void InitIndexes();
  void MatchWithMc();

  // output branches
  AnalysisTree::EventHeader* events_{nullptr};
  AnalysisTree::Particles* lambda_reco_{nullptr};
  AnalysisTree::Particles* lambda_sim_{nullptr};
  AnalysisTree::Matching* lambda_reco2sim_{nullptr};

  // input branches
  AnalysisTree::Particles* mc_particles_{nullptr};
  AnalysisTree::TrackDetector* rec_tracks_{nullptr};
  AnalysisTree::Matching* rec_to_mc_{nullptr};
  AnalysisTree::EventHeader* sim_events_{nullptr};

  std::vector<OutputContainer> canditates_;
  DecayContainer decay_;

  // field ids of simulated events
  int b_field_id_{-1};

  // field ids of input simulated mother
  int mother_id_field_id_{-1};

  // field ids for selected mother candidates kinematic parameters
  int x_field_id_{-1};
  int y_field_id_{-1};
  int z_field_id_{-1};
  int daughter1_id_field_id_{-1};
  int daughter2_id_field_id_{-1};
  int daughter3_id_field_id_{-1};
  int is_signal_field_id_{-1};
  int px_err_field_id_{-1};
  int py_err_field_id_{-1};
  int pz_err_field_id_{-1};
  int mass_err_field_id_{-1};

  // field ids for mother candidate cutting parameters for two daughters
  int chi2primpos_field_id_{-1};
  int chi2primneg_field_id_{-1};
  int distance_field_id_{-1};
  int cosinepos_field_id_{-1};
  int cosineneg_field_id_{-1};
  int chi2geo_field_id_{-1};
  int l_field_id_{-1};
  int ldl_field_id_{-1};
  int isfrompv_field_id_{-1};
  int cosinetopo_field_id_{-1};
  int chi2topo_field_id_{-1};
  int nhits_pos_field_id_{-1};
  int nhits_neg_field_id_{-1};

  // field ids for mother candidate cutting parameters for third daughter
  int chi2primthird_field_id_{-1};
  int distancethird_field_id_{-1};
  int cosinethird_field_id_{-1};
  int chi2geothree_field_id_{-1};
  int cosinetopothree_field_id_{-1};
  int chi2topothree_field_id_{-1};
  int nhits_third_field_id_{-1};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
