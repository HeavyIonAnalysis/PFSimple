#include "PFTaskManager.h"

#include "AnalysisTree/PlainTreeFiller.hpp"

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt mother_pdg\n";
    return EXIT_FAILURE;
  }

  const bool make_plain_tree{true};

  DecayContainer decay;
  if (!decay.SetDecay(atoi(argv[2]))){
    return EXIT_FAILURE;
  }

  CutsContainer cuts;
  cuts.CancelCuts();
  cuts.SetCutChi2PrimPos(18.6);
  cuts.SetCutChi2PrimNeg(18.6);
  cuts.SetCutDistance(1.);
  cuts.SetCutChi2Geo(3.);
  cuts.SetCutLdL(5.);

  if (decay.GetNdaughters() == 3) {
    cuts.SetCutChi2PrimThird(18.42);
    cuts.SetCutDistanceThird(.6);
    //cuts.SetCutChi2GeoThree(3.);
    //cuts.SetCutCosineTopoDown(0.);
    //cuts.SetCutCosineTopoUp(0.99996);
    //cuts.SetCutLThreeDown(16.);
  }

  const std::string& filename = argv[1];

  PFTaskManager man{};

  man.SetOutputName("PFSimpleOutput.root", "sTree");

  auto* in_converter = new ConverterIn(decay, cuts);
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name

  auto* out_converter = new ConverterOut(decay);
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks", "SimEventHeader"});

  man.AddTasks(in_converter, out_converter);
  man.Init({filename}, {"aTree"});
  man.Run(-1);// -1 = all events
  man.Finish();

//  if (make_plain_tree) {
//    std::ofstream filelist;
//    filelist.open("filelist.txt");
//    filelist << "PFSimpleOutput.root\n";
//    filelist.close();
//
//    AnalysisTree::TaskManager pl_man({{"filelist.txt"}}, {{"sTree"}});
//    pl_man.SetOutFileName("PFSimplePlainTree.root");
//
//    auto* tree_task_events = new AnalysisTree::PlainTreeFiller();
//    std::string branchname_events = std::string("Events");
//    tree_task_events->SetInputBranchNames({branchname_events});
//    tree_task_events->SetOutputBranchName(branchname_events);
//
//    auto* tree_task = new AnalysisTree::PlainTreeFiller();
//    std::string branchname_rec = decay.GetNameMother() + std::string("Candidates");
//    tree_task->SetInputBranchNames({branchname_rec});
//    tree_task->SetOutputBranchName(branchname_rec);
//
//    //     pl_man.AddTask(tree_task_events);
//    pl_man.AddTask(tree_task);
//
//    pl_man.Init();
//    pl_man.Run(-1);// -1 = all events
//    pl_man.Finish();
//  }

  return EXIT_SUCCESS;
}