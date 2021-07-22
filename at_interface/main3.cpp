//#include "AnalysisTree/PlainTreeFiller.hpp"
#include <PlainTreeFiller.hpp>
#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

#include "AnalysisTree/TaskManager.hpp"

int main(int argc, char** argv) {
   
  if (argc < 2) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  const bool make_plain_tree{false};

  const std::string& filename = argv[1];

  std::vector<Daughter> daughters = {2212,-211,1000010020};
  for (int idaughter = 0; idaughter < daughters.size() ; idaughter++) {
    daughters.at(idaughter).CancelCuts();
    daughters.at(idaughter).SetCutChi2Prim(18.42);
    //daughters.at(idaughter).SetCutCos(0.99825);
  }                                                                                                            
  Mother mother(3004);
  mother.CancelCuts();
  mother.SetCutDistance(1.0);
  mother.SetCutChi2GeoSM({3.0});
  mother.SetCutChi2TopoSM({29});
  mother.SetCutDistanceToSV(1.0);
  mother.SetCutChi2Geo(3.0);
  //mother.SetCutChi2Topo(29);
  //mother.SetCutLdL(3);
  Decay decay("H3L", mother, {daughters});
  
  auto* man = AnalysisTree::TaskManager::GetInstance();
  man->SetOutputName("PFSimpleOutput.root", "pTree");
  
  auto* in_converter = new ConverterIn();
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name
  
  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({decay});
  
  auto* out_converter = new ConverterOut();
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks", "SimEventHeader", "RecEventHeader"});
  out_converter->SetDecay(decay);

  //man.AddTasks(in_converter, out_converter);
  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  man->Init({filename}, {"rTree"});
  man->Run(-1);// -1 = all events
  man->Finish();
  //man->ClearTasks();

  if (make_plain_tree) {
    std::ofstream filelist;
    filelist.open("filelist.txt");
    filelist << "PFSimpleOutput.root\n";
    filelist.close();

    man->SetOutputName("PFSimplePlainTree.root", "plain_tree");

//    auto* tree_task_events = new AnalysisTree::PlainTreeFiller();
//    std::string branchname_events = "Events";
//    tree_task_events->SetInputBranchNames({branchname_events});
//    tree_task_events->AddBranch(branchname_events);

    auto* tree_task = new AnalysisTree::PlainTreeFiller();
    std::string branchname_rec = "Candidates";
    tree_task->SetInputBranchNames({branchname_rec});
    tree_task->AddBranch(branchname_rec);

    man->AddTask(tree_task);

    man->Init({"filelist.txt"}, {"pTree"});
    man->Run(-1);// -1 = all events
    man->Finish();
  }

  return EXIT_SUCCESS;
}
