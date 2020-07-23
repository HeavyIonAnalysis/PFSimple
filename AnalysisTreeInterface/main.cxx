#include "PFTaskManager.h"

int main(int argc, char** argv)
{
  if(argc < 2){
    std::cout << "Wrong number of arguments! Please use:\n  ./main filename\n";
    return EXIT_FAILURE;
  }

  CutsContainer cuts;
  cuts.CancelCuts();
  cuts.SetCutChi2PrimPos(3.);
  cuts.SetCutChi2PrimNeg(3.);
  cuts.SetCutDistance(1.);
  cuts.SetCutChi2Geo(3.);
  cuts.SetCutLdL(5.);

  const std::string& filename = argv[1];

  PFTaskManager man(filename, "aTree");
  man.SetOutFileName("PFSimpleOutput.root");
  man.SetOutTreeName("sTree");

  auto* in_converter = new ConverterIn(cuts);
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {{{"KfpfTracks", "pass_cuts"}, 1}}));
  in_converter->SetIsShine(false); //TODO maybe change name

  man.AddTasks(in_converter, new ConverterOut);
  man.Init();
  man.Run(-1); // -1 = all events
  man.Finish();

  return EXIT_SUCCESS;
}