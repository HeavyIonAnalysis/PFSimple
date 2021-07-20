//#include "AnalysisTree/PlainTreeFiller.hpp"
#include <PlainTreeFiller.hpp>
#include "PfSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

#include "AnalysisTree/TaskManager.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  const bool make_plain_tree{true};

  const std::string& filename = argv[1];

  Daughter proton(2212);
  Daughter pion(-211);
  
  proton.SetCutChi2Prim(26);
  proton.SetCutCos(0.99825);
  
  pion.SetCutChi2Prim(110);

  Mother lambda(3122);
  lambda.SetCutChi2Geo(11);
  lambda.SetCutChi2Topo(29);
  lambda.SetCutDistance(0.15);
  lambda.SetCutLdL(4);

  Decay lambda_pi_p("lambda", lambda, {pion, proton});

  auto* man = AnalysisTree::TaskManager::GetInstance();
  man->SetOutputName("PFSimpleOutput.root", "pTree");

  auto* in_converter = new ConverterIn();
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name

  auto* out_converter = new ConverterOut();
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks", "SimEventHeader", "RecEventHeader"});
  out_converter->SetDecay(lambda_pi_p);

  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetOutTask(out_converter);
  pf_task->SetDecay(lambda_pi_p);

  //  man.AddTasks(in_converter, out_converter);
  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  man->Init({filename}, {"rTree"});
  man->Run(-1);// -1 = all events
  man->Finish();
  man->ClearTasks();

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