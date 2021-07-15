//#include "AnalysisTree/PlainTreeFiller.hpp"
#include "PFSimpleTask.hpp"

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

  // ******** Lambda *********************************
  Daughter proton(2212, {2212, 211, -1});
  Daughter pion(-211, {-211});
  Mother lambda(3122);
  
//   proton.CancelCuts();
//   pion.CancelCuts();
//   lambda.CancelCuts();
  
//   proton.SetCutChi2Prim(18.42);
//   pion.SetCutChi2Prim(18.42);
//   lambda.SetCutDistance(1.);
//   lambda.SetCutChi2Geo(3.);
//   lambda.SetCutLdL(10.);
//   lambda.SetCutChi2TopoLower(5.);
//   lambda.SetCutInvMass(3.);
  
// // ***** optimized ******************
//   proton.SetCutChi2Prim(26);
//   proton.SetCutCos(0.99825);
//   
//   pion.SetCutChi2Prim(110);
// 
//   lambda.SetCutChi2Geo(11);
// //   lambda.SetCutChi2Topo(29);
//   lambda.SetCutDistance(0.15);
//   lambda.SetCutLdL(4);
//   
//   lambda.SetCutInvMass(6.);
// // **********************************

  
  Decay lambda_pi_p("lambda", lambda, {pion, proton});
  //**************************************************
  
  // ******* Xi - ************************************
  Daughter pion_from_xi(-211, {-211});
  Daughter lambda_from_xi(3122, {3122});
  Mother xi(3312);
  
  pion_from_xi.CancelCuts();
  lambda_from_xi.CancelCuts();
  xi.CancelCuts();
  
  pion_from_xi.SetCutChi2Prim(18.42);
  xi.SetCutDistance(1.);
  xi.SetCutChi2Geo(6.);
  xi.SetCutLdL(5.);
  xi.SetCutChi2Topo(5.);
  
//   
  Decay xi_pi_lambda("xi", xi, {pion_from_xi, lambda_from_xi});
  //**************************************************

  auto* man = AnalysisTree::TaskManager::GetInstance();
  man->SetOutputName("PFSimpleOutput.root");
  
  auto* in_converter = new ConverterIn();
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name

  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({lambda_pi_p});  
//   pf_task->SetDecays({lambda_pi_p, xi_pi_lambda});  
  
  auto* out_converter = new ConverterOut();
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks", "SimEventHeader", "RecEventHeader"});

  //  man.AddTasks(in_converter, out_converter);
  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  man->Init({filename}, {"rTree"});
  man->Run(-1);// -1 = all events
  man->Finish();

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