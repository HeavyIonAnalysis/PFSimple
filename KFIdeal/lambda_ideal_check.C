std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks);

void lambda_ideal_check(const TString infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/na61.aTree.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  
  auto* sim_tracks = new AnalysisTree::TrackDetector();
  auto* rec_tracks = new AnalysisTree::TrackDetector();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* matching_vtx_sim = new AnalysisTree::Matching();
  auto* matching_rec_sim = new AnalysisTree::Matching();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) file->Get("Configuration");
  
  auto* vtx_matching = new AnalysisTree::Matching();
  tree->SetBranchAddress("SimTracks", &sim_tracks);
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);
  tree->SetBranchAddress("VtxTracks", &vtx_tracks);
  tree->SetBranchAddress(config->GetMatchName("VtxTracks", "SimTracks").c_str(), &matching_vtx_sim);
  tree->SetBranchAddress(config->GetMatchName("KfpfTracks", "SimTracks").c_str(), &matching_rec_sim);
  
  tree -> GetEntry(104);  
  auto& sim_track = sim_tracks -> GetChannel(1603);
  
  auto daughters_ids = FindDaughters(sim_track, sim_tracks);
  std::cout << daughters_ids.size() << " daugherts\n";
  
  const int rec_ids[2] = {matching_rec_sim->GetMatchInverted(daughters_ids[0]), matching_rec_sim->GetMatchInverted(daughters_ids[1])};
  
  for(int i=0; i<daughters_ids.size(); i++)
  {
    std::cout << daughters_ids[i] << "\t" << rec_ids[i] << std::endl;
    
    const auto& daughter = rec_tracks->GetChannel(rec_ids[i]);
    
    const int nhits_vtpc1 = daughter.GetField<int>(2);
    const int nhits_vtpc2 = daughter.GetField<int>(3);
    const int nhits_mtpc = daughter.GetField<int>(4);
    const int nhits_pot_vtpc1 = daughter.GetField<int>(6);
    const int nhits_pot_vtpc2 = daughter.GetField<int>(7);
    const int nhits_pot_mtpc = daughter.GetField<int>(8);
    const float chi2ndf = daughter.GetField<float>(1);
    
    std::cout << nhits_vtpc1 + nhits_vtpc2 << " [>=15]\n";
    std::cout << nhits_vtpc1 + nhits_vtpc2 + nhits_mtpc << " [ >=30]\n";
    std::cout << 1.*(nhits_vtpc1 + nhits_vtpc2 + nhits_mtpc)/(nhits_pot_vtpc1 + nhits_pot_vtpc2 + nhits_pot_mtpc) << " [0.5-1.1]\n";
    std::cout << chi2ndf << " [<=3]\n\n";
  }
}

std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks)
{
  std::vector<int> daughters_ids{};
  const int lambda_id = mc_track.GetId();
  for(int id=lambda_id; id<sim_tracks->GetNumberOfChannels(); ++id)
  {
    const auto& daughter = sim_tracks->GetChannel(id);
    if( daughter.GetField<int>(2) == lambda_id )
      daughters_ids.emplace_back(daughter.GetId());
  }
  return daughters_ids;
}