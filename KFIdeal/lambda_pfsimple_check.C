void lambda_pfsimple_check(const TString infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/na61.aTree.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  
  auto* sim_tracks = new AnalysisTree::TrackDetector();
  auto* rec_tracks = new AnalysisTree::TrackDetector();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* matching_vtx_sim = new AnalysisTree::Matching();
  auto* matching_rec_sim = new AnalysisTree::Matching();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) file->Get("Configuration");
  
  tree->SetBranchAddress("SimTracks", &sim_tracks);
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);
  tree->SetBranchAddress("VtxTracks", &vtx_tracks);
  tree->SetBranchAddress(config->GetMatchName("VtxTracks", "SimTracks").c_str(), &matching_vtx_sim);
  tree->SetBranchAddress(config->GetMatchName("KfpfTracks", "SimTracks").c_str(), &matching_rec_sim);
  
  const int nEntry = 103;
  const int nReco = 141;
//   const int nMother = 1603;
  
  tree -> GetEntry(nEntry);
  auto& rec_track = rec_tracks -> GetChannel(nReco);
  const int simtrackid1 = matching_rec_sim->GetMatch(nReco);
  
  
  const auto& sim_track = sim_tracks->GetChannel(simtrackid1);
  
  std::cout << simtrackid1 << "\t" << sim_track.GetField<int>(2) << "\t" << sim_track.GetField<int>(3) << std::endl;  
  
//   for(int i=0; i<sim_tracks->GetNumberOfChannels(); i++)
//   {
//     auto& sim_track = sim_tracks -> GetChannel(i);
//     if(sim_track.GetField<int>(2) == nMother)
//     {
//       std::cout << i << "\t" << matching_rec_sim -> GetMatchInverted(i) << "\n";
//       const auto& daughter = rec_tracks->GetChannel(matching_rec_sim -> GetMatchInverted(i));
//       const int nhits_vtpc1 = daughter.GetField<int>(2);
//       const int nhits_vtpc2 = daughter.GetField<int>(3);
//       const int nhits_mtpc = daughter.GetField<int>(4);
//       const int nhits_pot_vtpc1 = daughter.GetField<int>(6);
//       const int nhits_pot_vtpc2 = daughter.GetField<int>(7);
//       const int nhits_pot_mtpc = daughter.GetField<int>(8);
//       const float chi2ndf = daughter.GetField<float>(1);
//       
//       std::cout << nhits_vtpc1 + nhits_vtpc2 << " [>=15]\n";
//       std::cout << nhits_vtpc1 + nhits_vtpc2 + nhits_mtpc << " [ >=30]\n";
//       std::cout << 1.*(nhits_vtpc1 + nhits_vtpc2 + nhits_mtpc)/(nhits_pot_vtpc1 + nhits_pot_vtpc2 + nhits_pot_mtpc) << " [0.5-1.1]\n";
//       std::cout << chi2ndf << " [<=3]\n\n";
//     }
//     
//   }
  
  
  
}