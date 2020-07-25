#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_

#include <Interface/OutputContainer.h>

#include "AnalysisTree/FillTask.hpp"
#include "AnalysisTree/Detector.hpp"

class ConverterOut : public AnalysisTree::FillTask {
 public:
  ConverterOut() = default;
  ~ConverterOut() override = default;

  void Init(std::map<std::string, void*>& branches) override;
  void Exec() override;
  void Finish() override {}

  void SetCandidates(const std::vector<OutputContainer>& canditates) { canditates_ = canditates; }

 protected:

  void InitIndexes();
  void MatchWithMc();

  AnalysisTree::Particles* lambda_reco_{nullptr};

  AnalysisTree::Particles* mc_particles_{nullptr};
  AnalysisTree::TrackDetector* rec_tracks_{nullptr};

  std::vector<OutputContainer> canditates_;

  // field ids for selected lambda candidates kinematic parameters
  int x_field_id_{-1};
  int y_field_id_{-1};
  int z_field_id_{-1};
  int mass_field_id_{-1};
//  int rap_lab_field_id_{-1};
  int rap_cm_field_id_{-1};
//  int pdg_field_id_w_{-1};
  int daughter1_id_field_id_{-1};
  int daughter2_id_field_id_{-1};

  // field ids for lambda candidate cutting parameters
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
};

#endif //KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
