#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOutTree.hpp"

#include "AnalysisTree/TaskManager.hpp"

int main(int argc, char** argv) {

  if (argc < 2) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  const std::string& filename = argv[1];

  //const int pid_mode = 1;// mc-pid
  std::vector<Daughter> daughters = {2212, -211, 1000010020};

  const int pid_mode = 2;// rec-pid
  float purity = 0.7;
  //   std::array<float,3> purity_pdg = {0.7,0.7,0.1};
  //   const int pid_mode = 0; // no-pid
  //   std::vector<Daughter> daughters = {{2212, {1}}, {-211, {-1}}, {1000010020, {1}}}

  for (size_t idaughter = 0; idaughter < daughters.size(); ++idaughter) {
    daughters.at(idaughter).CancelCuts();
    daughters.at(idaughter).SetCutChi2Prim(18.42);
    //daughters.at(idaughter).SetCutCos(0.99825);
  }
  Mother mother(3004);
  mother.CancelCuts();
  mother.SetCutDistance(1.0);
  mother.SetCutChi2GeoSM({3.0});
  //mother.SetCutChi2TopoSM({29});
  mother.SetCutDistanceToSV(1.0);
  //mother.SetCutChi2Geo(3.0);
  //mother.SetCutChi2Topo(29);
  //mother.SetCutLdL(3);
  Decay decay("H3L", mother, {daughters});

  std::string outname;
  if (argc > 2) outname = argv[2];
  else
    outname = "PFSimpleOutput";
  std::string outfilename;
  if (pid_mode == 1) outfilename = outname + "_pid_mc.root";
  if (pid_mode > 1) outfilename = outname + "_pid_rec.root";

  std::string tracks_name;
  if (pid_mode < 2) tracks_name = "VtxTracks";
  else
    tracks_name = "RecTracks";

  auto* man = AnalysisTree::TaskManager::GetInstance();
  //man->SetOutputName("PFSimpleOutput.root", "pTree");

  auto* in_converter = new ConverterIn();
  in_converter->SetRecEventHeaderName("RecEventHeader");
  in_converter->SetRecTracksName(tracks_name);
  in_converter->SetSimTracksName("SimParticles");
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut(tracks_name + ".pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name
  in_converter->SetPidMode(pid_mode);
  if (pid_mode > 2) {
    in_converter->SetPidPurity(purity);
    //in_converter->SetPidPurityProton(purity_pdg.at(0)); // in pid-mode 3 pdg-specific purity possible
    //in_converter->SetPidPurityPion(purity_pdg.at(1));
    //in_converter->SetPidPurityDeuteron(purity_pdg.at(2));
  }

  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({decay});

  auto* out_converter = new ConverterOutTree();
  out_converter->SetSimEventHeaderName("SimEventHeader");
  out_converter->SetRecTracksName("RecParticles");
  out_converter->SetSimTracksName("SimParticles");
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetDecay(decay);
  //   out_converter->SetPidMode(pid_mode);
  out_converter->SetOutFilename(outfilename);

  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  if (pid_mode < 2)
    man->Init({filename}, {"rTree"});
  else
    man->Init({filename}, {"pTree"});
  man->Run(-1);// -1 = all events
  man->Finish();
  man->ClearTasks();

  return EXIT_SUCCESS;
}
