#include "PFTaskManager.h"
#include "AnalysisTree/PlainTreeFiller.hpp"

int main(int argc, char** argv)
{
  
  /**
   * argv[1]: list with inputfiles \n
   * optional arguments: \n
   * argv[2]: name of outputfile \n
   * argv[3]: mother pdg - if decay from list is selected (set_user_decay{false}) \n
   */
  
  if(argc < 2){
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  const bool set_user_decay{true};
  
  const bool make_plain_tree{false};  // only works with updated PlainTreeFiller

  DecayContainer decay;

  if(argc < 4 && set_user_decay == false)
    std::cout << "Default decay is set. If you want to change the decay, please use:\n  ./main filelist.txt outname pdgcode\n";
    
  //Set decay from list
  if (argc == 4){
    if (set_user_decay == true)
      {
	std::cout << "There is a conflict: pdgcode for list-decay is given and user-decay is set! Please reset user decay or use:\n  ./main filelist.txt\n";
	return EXIT_FAILURE;
      }
    if (decay.SetDecay(atoi(argv[3])) == false) return EXIT_FAILURE;
  }
  
  //Set user decay
  if (set_user_decay == true) {
    decay.SetNameMother("H3L");
    decay.SetPdgMother(3004);
    decay.SetNdaughters(3);
    decay.SetPdgDaughterPos(2212);
    decay.SetPdgDaughterNeg(-211);
    decay.SetPdgDaughterThird(1000010020);
  
    /* Optional: User can set additional particle pids for a daughter that will be considered as candidates for the reconstruction of the mother.
       eg.: decay.SetPdgDaughterPos({2212,1,211,-1});
       vector<int> with: 1st position: daughter pid, 2nd position: =1: nuclei pids (pid > 1000000000) will be considered as daughter candidates, =0: no nuclei pids considered, 3rd, 4th, etc. & positions: alternative pids
    */

    if (decay.GetNdaughters() > 2 && decay.GetPdgDaughterThird() == 0)
      {
	std::cout <<"Number of daughters was set to 3, but only 2 daughter pdgs were given!\n";
	return EXIT_FAILURE;
      }

    if (decay.GetNdaughters() == 2 && decay.GetPdgDaughterThird() != 0)
      {
	std::cout <<"Number of daughters was set to 2, but 3 daughter pdgs were given!\n";
	return EXIT_FAILURE;
      }

    if (decay.GetNdaughters() > 3)
      {
	std::cout <<"Number of daughters was set to > 3, but only 2 or 3 daughters are allowed!\n";
	return EXIT_FAILURE;
      }
    
    std::cout << "User decay is set."<<std::endl;
  }

  if (decay.GetPdgDaughterPosCandidates().size() > 1)
    std::cout << "Alternative candidates for positive daughter are set."<<std::endl;

  if (decay.GetPdgDaughterNegCandidates().size() > 1)
    std::cout << "Alternative candidates for negative daughter are set."<<std::endl;

  if (decay.GetPdgDaughterThirdCandidates().size() > 1)
    std::cout << "Alternative candidates for third daughter are set."<<std::endl;

  if (decay.GetNdaughters() == 2) std::cout<<"Reconstruction of "<<decay.GetNameMother()<<"-decay: "<<decay.GetPdgMother()<<" -> "<<decay.GetPdgDaughterPos()<<" + "<<decay.GetPdgDaughterNeg()<<std::endl;
  if (decay.GetNdaughters() == 3) std::cout<<"Reconstruction of "<<decay.GetNameMother()<<"-decay: "<<decay.GetPdgMother()<<" -> "<<decay.GetPdgDaughterPos()<<" + "<<decay.GetPdgDaughterNeg()<<" + "<<decay.GetPdgDaughterThird()<<std::endl;

  CutsContainer cuts;
  cuts.CancelCuts();
  cuts.SetCutChi2PrimPos(18.42);
  cuts.SetCutChi2PrimNeg(18.42);
  cuts.SetCutDistance(1.);
  cuts.SetCutChi2Geo(3.);
  //cuts.SetCutLdL(5.);

  if (decay.GetNdaughters() == 3)
    {
      cuts.SetCutChi2PrimThird(18.42);
      cuts.SetCutDistanceThird(.6);
      //cuts.SetCutChi2GeoThree(3.);
      //cuts.SetCutCosineTopoDown(0.);
      //cuts.SetCutCosineTopoUp(0.99996); 
      //cuts.SetCutLThreeDown(16.);
    }


  std::string outname;
  if (argc > 2) outname = argv[2];
  else outname = "PFSimpleOutput";
  
  const std::string& filename = argv[1];
  std::string outfilename = outname + ".root";

  std::cout<<"Output file:\n"<<outfilename<<std::endl;

  PFTaskManager man(filename, "aTree");
  man.SetOutTreeName("sTree");
  man.SetOutFileName(outfilename);

  auto* in_converter = new ConverterIn(decay,cuts);
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {{{"VtxTracks", "pass_cuts"}, 1}}));
  in_converter->SetIsShine(false); //TODO maybe change name

  auto* out_converter = new ConverterOut(decay);
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks","SimEventHeader"});

  man.AddTasks(in_converter, out_converter);
  man.Init();
  man.Run(-1); // -1 = all events
  man.Finish();
  

  if(make_plain_tree){
    std::ofstream filelist;
    std::string filelist_stree = outname + "_list.txt";
	   
    filelist.open (filelist_stree);
    filelist << outfilename;
    filelist << "\n";
    filelist.close();

    AnalysisTree::TaskManager pl_man({{filelist_stree}}, {{"sTree"}});
    std::string filename_pl = outname + "_PlainTree.root";
    pl_man.SetOutFileName(filename_pl);
    std::cout<<"PlainTree file:\n"<<outfilename<<std::endl;

    auto* tree_task_events = new AnalysisTree::PlainTreeFiller();
    std::string branchname_events = std::string("Events"); 
    tree_task_events->SetInputBranchNames({branchname_events});
    tree_task_events->SetOutputBranchName(branchname_events);

    auto* tree_task = new AnalysisTree::PlainTreeFiller();
    std::string branchname_rec = decay.GetNameMother() + std::string("Candidates"); 
    tree_task->SetInputBranchNames({branchname_rec}); 
    tree_task->SetOutputBranchName(branchname_rec);  

    pl_man.AddTask(tree_task_events);
    pl_man.AddTask(tree_task);

    pl_man.Init();
    pl_man.Run(-1); // -1 = all events
    pl_man.Finish();
  }

  return EXIT_SUCCESS;
}
