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
  
//   const int pid_mode = 0; // no pid (topo)
//   const int pid_mode = 1; // mc pid
  const int pid_mode = 2; // rec pid

  // ******** Lambda *********************************
  Daughter proton(2212, {2212, 2});
  Daughter pion(-211, {-211, -2});
  Mother lambda(3122);
  
  proton.CancelCuts();
  pion.CancelCuts();
  lambda.CancelCuts();
  
  proton.SetCutChi2Prim(18.42);
  pion.SetCutChi2Prim(18.42);
  lambda.SetCutDistance(1.);
  lambda.SetCutChi2Geo(3.);
  lambda.SetCutLdL(10.);
  lambda.SetCutChi2TopoLower(5.);
  lambda.SetCutInvMass(3.);
  
  Decay lambda_pi_p("lambda", lambda, {pion, proton});
  lambda_pi_p.SetIsApplyMassConstraint();
  //**************************************************
  
  // ******* Xi - ************************************
  Daughter pion_from_xi(-211, {-211, -2});
  Daughter lambda_from_xi(3122);
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
  xi_pi_lambda.SetIsTransportToPV();
  //**************************************************

  // ******* Omega - ************************************
  Daughter kaon_from_omega(-321, {-321, -2});
  Daughter lambda_from_omega(3122);
  Mother omega(3334);
  
  kaon_from_omega.CancelCuts();
  lambda_from_omega.CancelCuts();
  omega.CancelCuts();
  
  kaon_from_omega.SetCutChi2Prim(18.42);
  omega.SetCutDistance(1.);
  omega.SetCutChi2Geo(6.);
  omega.SetCutLdL(5.);
  omega.SetCutChi2Topo(5.);
  
//   
  Decay omega_K_lambda("omega", omega, {kaon_from_omega, lambda_from_omega});
  omega_K_lambda.SetIsTransportToPV();
  //**************************************************

  auto* man = AnalysisTree::TaskManager::GetInstance();
  man->SetOutputName("PFSimpleOutput.root");
  
  auto* in_converter = new ConverterIn();
  in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name
  in_converter->SetPidMode(pid_mode);
  
//   in_converter->SetAncestorPdgsToBeConsidered({3312});

  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({lambda_pi_p, xi_pi_lambda, omega_K_lambda});  
  
  auto* out_converter = new ConverterOut();
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetInputBranchNames({"SimParticles", "VtxTracks", "RecParticles", "SimEventHeader", "RecEventHeader"});
  
//   AnalysisTree::Cuts* post_cuts = new AnalysisTree::Cuts("post_cuts", {AnalysisTree::RangeCut("Candidates.generation", 0.9, 100)});
//   AnalysisTree::Cuts* post_cuts = new AnalysisTree::Cuts("post_cuts", {AnalysisTree::EqualsCut("Candidates.generation", 0)});
//   AnalysisTree::Cuts* post_cuts = new AnalysisTree::Cuts("post_cuts", {AnalysisTree::RangeCut("Candidates.mass", 1.09, 1.14)});
//   out_converter->SetOutputCuts(post_cuts);

  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  man->Init({filename}, {"aTree"});
  man->Run(-1);// -1 = all events
  man->Finish();

  return EXIT_SUCCESS;
}