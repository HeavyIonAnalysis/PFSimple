#ifndef AT_KFC_in_HH
#define AT_KFC_in_HH

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "AnalysisTree/Track.h"
#include "AnalysisTree/Hit.h"
#include "AnalysisTree/Detector.h"
#include "AnalysisTree/Configuration.h"
#include "AnalysisTree/EventHeader.h"
#include "AnalysisTree/Matching.h"
#include "AnalysisTree/Cuts.h"

#include "InputContainer.h"

class InConverter
{
public:
  
  // Constructors/Destructors ---------
  InConverter() = default;
  virtual ~InConverter() = default;
  
  void InitAnalysisTree(const std::string& file_name, const std::string& tree_name);
  InputContainer CreateInputContainer(int iEvent);

  void SetIsShine(bool is=true) { is_shine_ = is; }
  int GetNEvents(){return in_chain_ -> GetEntries();};  

  void SetCuts(const CutsContainer& cuts) { cuts_ = cuts; };
  void SetTrackCuts(AnalysisTree::Cuts* const cuts) { track_cuts_ = cuts; };

protected:
  std::vector<float> GetCovMatrixCbm(const AnalysisTree::Track& track) const;
  std::vector<float> GetCovMatrixShine(const AnalysisTree::Track& track) const;
  bool IsGoodTrack(const AnalysisTree::Track& rec_track);
  void FillTrack(const AnalysisTree::Track& rec_track, InputContainer& input_info);
  int GetTrackPid(const AnalysisTree::Track& rec_track);

  TFile* in_file_{nullptr};
  TChain* in_chain_{nullptr};

  AnalysisTree::TrackDetector* kf_tracks_{nullptr};
  AnalysisTree::TrackDetector* sim_tracks_{nullptr};
  AnalysisTree::Matching* kf2sim_tracks_{nullptr};
  AnalysisTree::EventHeader* rec_event_header_{nullptr};
  AnalysisTree::EventHeader* sim_event_header_{nullptr};
  AnalysisTree::Configuration* config_{nullptr};
// <<<<<<< HEAD
  
// =======
// 
  AnalysisTree::Cuts* track_cuts_{nullptr};
// 
//   const int nParametersKF{6};
// >>>>>>> lubynets-order
  CutsContainer cuts_;

  // field ids for input parameters
  int q_field_id_{AnalysisTree::UndefValueInt};
  int pdg_field_id_{AnalysisTree::UndefValueInt};
  int par_field_id_ {AnalysisTree::UndefValueInt};
  int mf_field_id_ {AnalysisTree::UndefValueInt};
  int cov_field_id_{AnalysisTree::UndefValueInt};
  int mother_id_field_id_{AnalysisTree::UndefValueInt};
  int sim_pdg_field_id_{AnalysisTree::UndefValueInt};
  int passcuts_field_id_{AnalysisTree::UndefValueInt};

  bool is_shine_{false};
  
  
};

#endif
