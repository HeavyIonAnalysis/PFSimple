std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks);
AnalysisTree::Cuts* GetShineTrackCuts();

void lambda_qa(const TString infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/na61.aTree.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  
  TH1F h_mc_mother("mc_mother", "", 1000, 1., 2.);
  TH2F h_pion_px("pion_px", "x=sim", 1000, -1, 1, 1000, -1, 1);
  TH2F h_pion_py("pion_py", "x=sim", 1000, -0.5, 0.5, 1000, -0.5, 0.5);
  TH2F h_pion_pz("pion_pz", "x=sim", 1000, -1, 4, 1000, -1, 4);
  TH2F h_pion_e("pion_e", "x=sim", 1000, 0, 10, 1000, 0, 10);
  TH2F h_proton_px("proton_px", "x=sim", 1000, -2, 2, 1000, -2, 2);
  TH2F h_proton_py("proton_py", "x=sim", 1000, -1, 1, 1000, -1, 1);
  TH2F h_proton_pz("proton_pz", "x=sim", 1000, -1, 10, 1000, -1, 10);
  TH2F h_proton_e("proton_e", "x=sim", 1000, 0, 10, 1000, 0, 10);
  TH2F h_reco_charge("reco_charge", "x=reco0", 3, -1.5, 1.5, 3, -1.5, 1.5);
  
  auto* sim_tracks = new AnalysisTree::TrackDetector();
  auto* rec_tracks = new AnalysisTree::TrackDetector();
  auto* matching = new AnalysisTree::Matching();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) file->Get("Configuration");
  config->Print();
  
  tree->SetBranchAddress("SimTracks", &sim_tracks);
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);
  tree->SetBranchAddress(config->GetMatchName("KfpfTracks", "SimTracks").c_str(), &matching);        // SHINE
  
  AnalysisTree::Cuts* track_cuts = GetShineTrackCuts();
  track_cuts->Init(*config);
    
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
        
        if(mc_daughter0.GetField<int>(1)!=2212 && mc_daughter0.GetField<int>(1)!=-211) continue;
        if(mc_daughter1.GetField<int>(1)!=2212 && mc_daughter1.GetField<int>(1)!=-211) continue;
                
        const auto mc_daughter0_mom = mc_daughter0.GetMomentum(mc_daughter0.GetField<int>(1));
        const auto mc_daughter1_mom = mc_daughter1.GetMomentum(mc_daughter1.GetField<int>(1));
        
        if(mc_daughter0.GetField<int>(1)==2212 && mc_daughter0_mom.E() < 0.938)
          std::cout << mc_daughter0_mom.E() << "\n";
        if(mc_daughter1.GetField<int>(1)==2212 && mc_daughter1_mom.E() < 0.938)
          std::cout << mc_daughter1_mom.E() << "\n";          
        
        const auto mc_mother_mom = mc_daughter0_mom + mc_daughter1_mom;
        h_mc_mother.Fill(mc_mother_mom.M());
        
        const int rec_ids[2] = {matching->GetMatchInverted(daughters_ids[0]), matching->GetMatchInverted(daughters_ids[1])};
        if (rec_ids[0] < 0 || rec_ids[1] < 0) continue;
        if (rec_ids[0] >= rec_tracks->GetNumberOfChannels() || rec_ids[1] >= rec_tracks->GetNumberOfChannels()) continue;

        const auto& reco_daughter0 = rec_tracks->GetChannel(rec_ids[0]);
        const auto& reco_daughter1 = rec_tracks->GetChannel(rec_ids[1]);
                
//         if ( !(track_cuts->Apply(reco_daughter0) && track_cuts->Apply(reco_daughter0)) ) continue;
        
        const auto reco_daughter0_mom = mc_daughter0.GetField<int>(0) < 0 ? reco_daughter0.GetMomentum(0.140f) : reco_daughter0.GetMomentum(0.938f);
        const auto reco_daughter1_mom = mc_daughter1.GetField<int>(0) < 0 ? reco_daughter1.GetMomentum(0.140f) : reco_daughter1.GetMomentum(0.938f);
        
//         if(!(  (reco_daughter0.GetField<int>(1)==1 && reco_daughter1.GetField<int>(1)==-1 ) || (reco_daughter0.GetField<int>(1)==-1 && reco_daughter1.GetField<int>(1)==1 ) ))
//         {
//           std::cout << reco_daughter0.GetField<int>(1) << "\t" << reco_daughter1.GetField<int>(1) << std::endl;
//           std::cout << mc_daughter0.GetField<int>(0) << "\t" << mc_daughter1.GetField<int>(0) << std::endl;
//           std::cout << mc_daughter0.GetField<int>(1) << "\t" << mc_daughter1.GetField<int>(1) << std::endl << std::endl;
//         }
        
//         std::cout << reco_daughter0.GetField<int>(1) << "\t" << reco_daughter1.GetField<int>(1) << std::endl;
        h_reco_charge.Fill(reco_daughter0.GetField<int>(1), reco_daughter1.GetField<int>(1));

        if(mc_daughter0.GetField<int>(0) < 0)
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
  h_reco_charge.Write();
  out_file -> Close();  
}

std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks)
{
  std::vector<int> daughters_ids{};
  const int lambda_id = mc_track.GetId();
  for(int id=lambda_id; id<sim_tracks->GetNumberOfChannels(); ++id)
  {
    const auto& daughter = sim_tracks->GetChannel(id);
    if( daughter.GetField<int>(2) == lambda_id )          // SHINE
//     if( daughter.GetField<int>(0) == lambda_id )          // CBM
      daughters_ids.emplace_back(daughter.GetId());
  }

//   if(daughters_ids.size() != 2){
//     std::cout << "number of daughters = " << daughters_ids.size() << std::endl;
//   }
//  else
//    std::cout << "sim: " << daughters_ids[0] << "  " << daughters_ids[1] << std::endl;

  return daughters_ids;
}

AnalysisTree::Cuts* GetShineTrackCuts() {
//   AnalysisTree::SimpleCut dcax("dcax", -2-0.083, 2-0.083);
//   AnalysisTree::SimpleCut dcay("dcay", -1+0.006, 1+0.006);
  AnalysisTree::SimpleCut hits({"nhits_vtpc1", "nhits_vtpc2", "nhits_mtpc"},
                               [](std::vector<double>& hits) { return hits[0]+hits[1]>=15 && hits[0]+hits[1]+hits[2]>=30; } );

  AnalysisTree::SimpleCut hits_pot({"nhits_vtpc1", "nhits_vtpc2", "nhits_mtpc", "nhits_pot_vtpc1", "nhits_pot_vtpc2", "nhits_pot_mtpc"},
                                   [](std::vector<double>& hits) {
                                     const double total = hits[0]+hits[1]+hits[2];
                                     const double total_pot = hits[3]+hits[4]+hits[5];
                                     return total/total_pot>0.55 && total/total_pot<1.1;  } );

  AnalysisTree::Cuts* rec_track_cuts{ new AnalysisTree::Cuts("KfpfTracks")};
  rec_track_cuts->AddCuts({hits, hits_pot});
  return rec_track_cuts;
};
