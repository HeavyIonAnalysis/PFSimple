#include "PFTaskManager.h"

#include "AnalysisTree/PlainTreeFiller.hpp"

int main(int argc, char** argv)
{
  if(argc < 2){
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  const bool make_plain_tree{true};

  CutsContainer cuts;
  cuts.CancelCuts();
  cuts.SetCutChi2PrimPos(18.6);
  cuts.SetCutChi2PrimNeg(18.6);
  cuts.SetCutDistance(1.);
  cuts.SetCutChi2Geo(3.);
  cuts.SetCutLdL(5.);

  const std::string& filename = argv[1];

  PFTaskManager man(filename, "aTree");
  man.SetOutFileName("PFSimpleOutput.root");
  man.SetOutTreeName("sTree");

  auto* in_converter = new ConverterIn(cuts);
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {{{"VtxTracks", "pass_cuts"}, 1}}));
  in_converter->SetIsShine(false); //TODO maybe change name

  auto* out_converter = new ConverterOut;
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks"});

  man.AddTasks(in_converter, out_converter);
  man.Init();
  man.Run(-1); // -1 = all events
  man.Finish();

  if(make_plain_tree){
    std::ofstream filelist;
    filelist.open ("filelist.txt");
    filelist << "PFSimpleOutput.root\n";
    filelist.close();

    AnalysisTree::TaskManager pl_man({{"filelist.txt"}}, {{"sTree"}});
    pl_man.SetOutFileName("PFSimplePlainTree.root");

    auto* tree_task = new AnalysisTree::PlainTreeFiller();
    tree_task->SetInputBranchNames({"LambdaCandidates"});
    pl_man.AddTask(tree_task);

    pl_man.Init();
    pl_man.Run(-1); // -1 = all events
    pl_man.Finish();
  }

  return EXIT_SUCCESS;
}