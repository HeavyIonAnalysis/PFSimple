//#include "AnalysisTree/PlainTreeFiller.hpp"
#include "PFSimpleTask.hpp"
#include <PlainTreeFiller.hpp>

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  const bool make_plain_tree{false};

  const std::string& filename = argv[1];
  
  const int pid_mode = 1;       // MC-PID
  Daughter proton(2212);
  Daughter pion(-211);

  proton.SetCutChi2Prim(2);
  pion.SetCutChi2Prim(2);

  Mother lambda(3122);
  lambda.SetCutChi2Geo(3);
  lambda.SetCutDistance(1);
  lambda.SetCutLdL(5);

//  proton.CancelCuts();
//  pion.CancelCuts();
//  lambda.CancelCuts();

  Decay lambda_pi_p("lambda", lambda, {pion, proton});


  Daughter pion_plus(211);
  Daughter pion_minus(-211);
  pion_plus.SetCutChi2Prim(2);
  pion_minus.SetCutChi2Prim(2);

  Mother kshort(310);
  kshort.SetCutChi2Geo(3);
  kshort.SetCutDistance(1);
  kshort.SetCutLdL(5);

  Decay kshort_pi_pi("kshort", kshort, {pion_minus, pion_plus});

  auto* man = TaskManager::GetInstance();
  man->SetOutputName("PFSimpleOutput.root", "pTree");
  
  std::string tree_name;
  std::string rec_tracks_name;
  if(pid_mode < 2) {
    tree_name = "rTree";
    rec_tracks_name = "VtxTracks";
  }
  else {
    tree_name = "aTree";
    rec_tracks_name = "RecParticles";
  }

  auto* in_converter = new ConverterIn();
  in_converter->SetRecEventHeaderName("RecEventHeader");
  in_converter->SetRecTracksName(rec_tracks_name);
  in_converter->SetSimTracksName("SimParticles");
  
  in_converter->SetTrackCuts(new Cuts("Cut to reproduce KFPF", {EqualsCut((rec_tracks_name + ".pass_cuts").c_str(), 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name
  in_converter->SetPidMode(pid_mode);
  //   in_converter->SetPidPurity(min_pur);

  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({lambda_pi_p, kshort_pi_pi});

  auto* out_converter = new ConverterOut();
  out_converter->SetSimEventHeaderName("SimEventHeader");
  out_converter->SetRecTracksName(rec_tracks_name);
  out_converter->SetSimTracksName("SimParticles");
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetDecay(lambda_pi_p);

  //   Cuts* post_cuts = new Cuts("post_cuts", {RangeCut("Candidates.generation", 0.9, 100)});
  //   Cuts* post_cuts = new Cuts("post_cuts", {EqualsCut("Candidates.generation", 0)});
  //   Cuts* post_cuts = new Cuts("post_cuts", {RangeCut("Candidates.mass", 1.09, 1.14)});
  //   out_converter->SetOutputCuts(post_cuts);

  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  man->Init({filename}, {tree_name});
  man->Run(-1);// -1 = all events
               //   man->Run(9900);// -1 = all events
  man->Finish();
  man->ClearTasks();

  if (make_plain_tree) {
    std::ofstream filelist;
    filelist.open("filelist.txt");
    filelist << "PFSimpleOutput.root\n";
    filelist.close();

    man->SetOutputName("PFSimplePlainTree.root", "plain_tree");

    //    auto* tree_task_events = new PlainTreeFiller();
    //    std::string branchname_events = "Events";
    //    tree_task_events->SetInputBranchNames({branchname_events});
    //    tree_task_events->AddBranch(branchname_events);

    auto* tree_task = new PlainTreeFiller();
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