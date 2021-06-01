void lambda_pfmatch_check(const TString infile="/home/user/cbmdir/working/KFPS_match.1.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  
  auto* lambda_reco = new AnalysisTree::TrackDetector();
  auto* lambda_sim = new AnalysisTree::TrackDetector();
  auto* matching = new AnalysisTree::Matching();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) file->Get("Configuration");

  tree->SetBranchAddress("LambdaSimulated", &lambda_sim);
  tree->SetBranchAddress("LambdaCandidates", &lambda_reco);
  tree->SetBranchAddress(config->GetMatchName("LambdaCandidates", "LambdaSimulated").c_str(), &matching);
  
  const int nEntry = 104;
  tree->GetEntry(nEntry);
  
  for(int i=0; i<lambda_reco->GetNumberOfChannels(); i++)
  {
    auto& rec_lam = lambda_reco -> GetChannel(i);
    if(rec_lam.GetField<int>(4) == 1)
    {
      std::cout << "Dids: " << rec_lam.GetField<int>(2) << "\t" << rec_lam.GetField<int>(3) << std::endl;
      const int sim_lam_id = matching -> GetMatch(i);
      const auto& sim_lambda = lambda_sim->GetChannel(sim_lam_id);
      std::cout << sim_lambda.GetField<int>(1) << std::endl;      
    }
  }
  
}
