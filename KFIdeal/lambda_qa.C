std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks);

void lambda_qa(const TString infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/na61.aTree.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  
  TH1F h_mc_mother("mc_mother", "", 1000, 1., 2.);
  TH2F h_pion_px("pion_px", "x=sim", 1000, -3, 3, 1000, -3, 3);
  TH2F h_pion_py("pion_py", "x=sim", 1000, -3, 3, 1000, -3, 3);
  TH2F h_pion_pz("pion_pz", "x=sim", 1000, -2, 10, 1000, -2, 10);
  TH2F h_pion_e("pion_e", "x=sim", 1000, 0, 10, 1000, 0, 10);
  TH2F h_proton_px("proton_px", "x=sim", 1000, -3, 3, 1000, -3, 3);
  TH2F h_proton_py("proton_py", "x=sim", 1000, -3, 3, 1000, -3, 3);
  TH2F h_proton_pz("proton_pz", "x=sim", 1000, -2, 10, 1000, -2, 10);
  TH2F h_proton_e("proton_e", "x=sim", 1000, 0, 10, 1000, 0, 10);
  
  auto* sim_tracks = new AnalysisTree::TrackDetector();
  auto* rec_tracks = new AnalysisTree::TrackDetector();
  auto* matching = new AnalysisTree::Matching();
  
  tree->SetBranchAddress("SimTracks", &sim_tracks);
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);
  tree->SetBranchAddress("VKfpfTracks2SimTracks", &matching);
  
  const int n_entries = tree->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    tree->GetEntry(i_event);

    for(const auto& mc_track : *(sim_tracks->GetChannels()) )
    {
      if( mc_track.GetField<int>(1) == 3122 )
      {
        auto daughters_ids = FindDaughters(mc_track, sim_tracks);
        if(daughters_ids.size() != 2) continue;

        const auto& mc_daughter0 = sim_tracks->GetChannel(daughters_ids[0]);
        const auto& mc_daughter1 = sim_tracks->GetChannel(daughters_ids[1]);
        
        const auto mc_daughter0_mom = mc_daughter0.GetMomentum(mc_daughter0.GetField<int>(1));
        const auto mc_daughter1_mom = mc_daughter1.GetMomentum(mc_daughter1.GetField<int>(1));
                
        const auto mc_mother_mom = mc_daughter0_mom + mc_daughter1_mom;
        h_mc_mother.Fill(mc_mother_mom.M());
        
        const int rec_ids[2] = {matching->GetMatchInverted(daughters_ids[0]), matching->GetMatchInverted(daughters_ids[1])};
        if (rec_ids[0] < 0 || rec_ids[1] < 0) continue;
        if (rec_ids[0] >= rec_tracks->GetNumberOfChannels() || rec_ids[1] >= rec_tracks->GetNumberOfChannels()) continue;
        
        const auto& reco_daughter0 = rec_tracks->GetChannel(rec_ids[0]);
        const auto& reco_daughter1 = rec_tracks->GetChannel(rec_ids[1]);
        
        const auto reco_daughter0_mom = reco_daughter0.GetField<int>(1) < 0 ? reco_daughter0.GetMomentum(0.140f) : reco_daughter0.GetMomentum(0.938f);
        const auto reco_daughter1_mom = reco_daughter1.GetField<int>(1) < 0 ? reco_daughter1.GetMomentum(0.140f) : reco_daughter1.GetMomentum(0.938f);
        
        if(reco_daughter0.GetField<int>(1) < 0)
        {
          h_pion_px.Fill(mc_daughter0_mom.Px(), reco_daughter0_mom.Px());
          h_pion_py.Fill(mc_daughter0_mom.Py(), reco_daughter0_mom.Py());
          h_pion_pz.Fill(mc_daughter0_mom.Pz(), reco_daughter0_mom.Pz());
          h_pion_e.Fill(mc_daughter0_mom.E(), reco_daughter0_mom.E());
          h_proton_px.Fill(mc_daughter1_mom.Px(), reco_daughter1_mom.Px());
          h_proton_py.Fill(mc_daughter1_mom.Py(), reco_daughter1_mom.Py());
          h_proton_pz.Fill(mc_daughter1_mom.Pz(), reco_daughter1_mom.Pz());
          h_proton_e.Fill(mc_daughter1_mom.E(), reco_daughter1_mom.E());
        }
        else
        {
          h_proton_px.Fill(mc_daughter0_mom.Px(), reco_daughter0_mom.Px());
          h_proton_py.Fill(mc_daughter0_mom.Py(), reco_daughter0_mom.Py());
          h_proton_pz.Fill(mc_daughter0_mom.Pz(), reco_daughter0_mom.Pz());
          h_proton_e.Fill(mc_daughter0_mom.E(), reco_daughter0_mom.E());
          h_pion_px.Fill(mc_daughter1_mom.Px(), reco_daughter1_mom.Px());
          h_pion_py.Fill(mc_daughter1_mom.Py(), reco_daughter1_mom.Py());
          h_pion_pz.Fill(mc_daughter1_mom.Pz(), reco_daughter1_mom.Pz());
          h_pion_e.Fill(mc_daughter1_mom.E(), reco_daughter1_mom.E());
        }
        
      }
    }
  }
  
  auto *out_file = TFile::Open("out.lambdaQA.root", "recreate");
  h_mc_mother.Write();
  h_pion_px.Write();
  h_pion_py.Write();
  h_pion_pz.Write();
  h_pion_e.Write();
  h_proton_px.Write();
  h_proton_py.Write();
  h_proton_pz.Write();
  h_proton_e.Write();
  out_file -> Close();  
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

//   if(daughters_ids.size() != 2){
//     std::cout << "number of daughters = " << daughters_ids.size() << std::endl;
//   }
//  else
//    std::cout << "sim: " << daughters_ids[0] << "  " << daughters_ids[1] << std::endl;

  return daughters_ids;
}