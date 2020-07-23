#ifndef KFPS_manager_HH
#define KFPS_manager_HH

#include <string>
#include <utility>

#include "AnalysisTreeInterface/old/InConverter.h"
#include "AnalysisTreeInterface/old/OutConverter.h"
#include "SimpleFinder.h"

class Manager
{
public:

  void Init() {};
  void Run(int n_events=-1);
  void Finish(){};

  void SetInFileName(const std::string& name){in_file_name_ = name;};
  void SetInTreeName(const std::string& name){in_tree_name_ = name;};
  void SetIsShine(bool is=true) { is_shine_ = is; }

  void SetCuts(const CutsContainer& cuts) { cuts_ = cuts; };
  void SetTrackCuts(AnalysisTree::Cuts* const cuts) { track_cuts_ = cuts; };

 protected:
  std::string in_file_name_;
  std::string in_tree_name_;
  CutsContainer cuts_;
  AnalysisTree::Cuts* track_cuts_{nullptr};

  bool is_shine_{false};

};
#endif