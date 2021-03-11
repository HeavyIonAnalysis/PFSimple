#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFTASKMANAGER_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFTASKMANAGER_H_

#include "AnalysisTree/TaskManager.hpp"

#include <utility>

#include "ConverterIn.h"
#include "ConverterOut.h"

class PFTaskManager {

  enum eTasks {
    kInConverter = 0,
    kOutConverter
  };

 public:
  PFTaskManager() = default;
  PFTaskManager(PFTaskManager& other) = delete;
  void operator=(const PFTaskManager&) = delete;
  ~PFTaskManager() = default;

  void Init(const std::vector<std::string>& filelists, const std::vector<std::string>& in_trees) {
    man_->Init(filelists, in_trees);
  }

  void Run(long long nEvents);

  void Finish(){
    man_->Finish();
  }

  void AddTasks(ConverterIn* in_task, ConverterOut* out_task) {
    man_->AddTask(in_task);
    man_->AddTask(out_task);
  }

  void SetDecay(const DecayContainer& decay) { decay_ = decay; };
  void SetCuts(const CutsContainer& cuts) { cuts_ = cuts; };

  void SetOutputName(std::string file, std::string tree) {
    man_->SetOutputName(std::move(file), std::move(tree));
  }

 protected:
  AnalysisTree::TaskManager* man_{ AnalysisTree::TaskManager::GetInstance() };
  DecayContainer decay_;
  CutsContainer cuts_;
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFTASKMANAGER_H_
