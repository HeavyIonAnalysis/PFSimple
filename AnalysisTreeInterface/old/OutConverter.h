#ifndef KFC_out_AT_HH
#define KFC_out_AT_HH

#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "AnalysisTree/Track.h"
#include "AnalysisTree/Hit.h"
#include "AnalysisTree/Detector.h"
#include "AnalysisTree/Configuration.h"
#include "AnalysisTree/EventHeader.h"
#include "AnalysisTree/Matching.h"

#include "OutputContainer.h"

class OutConverter
{
public:
  
  // Constructors/Destructors ---------
  OutConverter() = default;
  virtual ~OutConverter() = default;
  
  void InitAT();
  void SetOutFileName(std::string name){out_file_name_ = std::move(name);};
  void WriteCandidates(const std::vector<OutputContainer>& canditates);
  
  void Finish();
  
protected:
  
  float RapidityShift(float pbeam);
  
  std::string out_file_name_{"KFPS_output.root"};
  TFile* out_file_{nullptr};
  TTree* out_tree_{nullptr};
  
  AnalysisTree::Configuration out_config_;
  AnalysisTree::Particles* lambda_reco_{nullptr};
   
  // field ids for selected lambda candidates kinematic parameters
  int x_field_id_{-1};
  int y_field_id_{-1};
  int z_field_id_{-1};
  int mass_field_id_{-1};
  int rap_lab_field_id_{-1};
  int rap_cm_field_id_{-1};
  int pdg_field_id_w_{-1};
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
#endif