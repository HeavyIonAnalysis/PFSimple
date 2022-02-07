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

  //   // ******** optimized cuts ***************
  //     Daughter proton(2212);
  //     Daughter pion(-211);
  //
  //     proton.SetCutChi2Prim(26);
  //     proton.SetCutCos(0.99825);
  //
  //     pion.SetCutChi2Prim(110);
  //
  //     Mother lambda(3122);
  //     lambda.SetCutChi2Geo(11);
  //     lambda.SetCutChi2Topo(29);
  //     lambda.SetCutDistance(0.15);
  //     lambda.SetCutLdL(4);
  //   //****************************************

  // ******** default kfpf cuts *************
  //   const int pid_mode = 1;
  const int pid_mode = 2;
  Daughter proton(2212, {2212, 2});// for MC-PID
  Daughter pion(-211, {-211, -2});
  Daughter pion_plus(211, {211, 2});
  Daughter pion_minus(-211, {-211, -2});

  //   const int pid_mode = 0;
  //   Daughter proton(2212, {1});       // for no-PID
  //   Daughter pion(-211, {-1});
  //   Daughter pion_plus(211, {1});
  //   Daughter pion_minus(-211, {-1});

  proton.SetCutChi2Prim(18.42);
  pion_plus.SetCutChi2Prim(18.42);

  pion.SetCutChi2Prim(18.42);
  pion_minus.SetCutChi2Prim(18.42);

  Mother lambda(3122);
  lambda.SetCutChi2Geo(3);
  lambda.SetCutDistance(1);
  lambda.SetCutLdL(5);

  Mother kshort(310);
  kshort.SetCutChi2Geo(3);
  kshort.SetCutDistance(1);
  kshort.SetCutLdL(5);
  // ***************************************

  //   // ******** no cuts **********************
  //   Daughter proton(2212);
  //   Daughter pion(-211);
  //   Mother lambda(3122);
  //   proton.CancelCuts();
  //   pion.CancelCuts();
  //   lambda.CancelCuts();
  //   // ***************************************

  Decay lambda_pi_p("lambda", lambda, {pion, proton});
  Decay kshort_pi_pi("kshort", kshort, {pion_minus, pion_plus});

  auto* man = TaskManager::GetInstance();
  man->SetOutputName("PFSimpleOutput.root", "pTree");

  auto* in_converter = new ConverterIn();

  //   SimpleCut kfpf_cut = EqualsCut("VtxTracks.pass_cuts", 1);
  //   SimpleCut mother_cut = EqualsCut("VtxTracks.mother_pdg", 3122);
  //   SimpleCut several_mother_cut = SimpleCut({"VtxTracks.mother_pdg"}, []( std::vector<double>& var ) { return std::fabs(var.at(0)-3122)<0.1 || std::fabs(var.at(0)-3312)<0.1; });
  //   Cuts* cuts = new Cuts("cuts", {kfpf_cut, several_mother_cut});
  //   in_converter->SetTrackCuts(cuts);

  //   in_converter->SetMotherPdgsToBeConsidered({3122});

  in_converter->SetTrackCuts(new Cuts("Cut to reproduce KFPF", {EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name
  in_converter->SetPidMode(pid_mode);
  //   in_converter->SetPidPurity(min_pur);

  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({lambda_pi_p, kshort_pi_pi});

  auto* out_converter = new ConverterOut();
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks", "RecParticles", "SimEventHeader", "RecEventHeader"});
  out_converter->SetDecay(lambda_pi_p);

  //   Cuts* post_cuts = new Cuts("post_cuts", {RangeCut("Candidates.generation", 0.9, 100)});
  //   Cuts* post_cuts = new Cuts("post_cuts", {EqualsCut("Candidates.generation", 0)});
  //   Cuts* post_cuts = new Cuts("post_cuts", {RangeCut("Candidates.mass", 1.09, 1.14)});
  //   out_converter->SetOutputCuts(post_cuts);

  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  man->Init({filename}, {"aTree"});
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