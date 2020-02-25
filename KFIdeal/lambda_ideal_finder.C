std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks);
AnalysisTree::Cuts* GetShineTrackCuts();

void lambda_ideal_finder(const TString infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/na61.aTree.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  TH1F inv_mass("inv_mass", "", 1000, 1.0, 2.);

  TH3F h3_rec_tracks("rec_tracks", "", 50, -3.2, 3.2, 50, 0, 3, 50, 0.5, 4);
  TH3F h3_not_rec_tracks("not_rec_tracks", "", 50, -3.2, 3.2, 50, 0, 3, 50, -0.5, 4);

  auto* sim_tracks = new AnalysisTree::TrackDetector();
  auto* rec_tracks = new AnalysisTree::TrackDetector();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* matching = new AnalysisTree::Matching();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) file->Get("Configuration");

  auto* vtx_matching = new AnalysisTree::Matching();
  tree->SetBranchAddress("SimTracks", &sim_tracks);
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);
  tree->SetBranchAddress("VtxTracks", &vtx_tracks);
  tree->SetBranchAddress(config->GetMatchName("KfpfTracks", "SimTracks").c_str(), &matching);
  
//   AnalysisTree::Cuts* track_cuts = GetShineTrackCuts();
//   track_cuts->Init(*config);

  const int n_entries = tree->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    tree->GetEntry(i_event);
    int number_of_lambdas = 0;

    for(const auto& mc_track : *(sim_tracks->GetChannels()) )
    {
      if( mc_track.GetField<int>(1) == 3122 )
      {
        auto daughters_ids = FindDaughters(mc_track, sim_tracks);
        if(daughters_ids.size() != 2) continue;

        const auto& mc_daughter0 = sim_tracks->GetChannel(daughters_ids[0]);
        const auto& mc_daughter1 = sim_tracks->GetChannel(daughters_ids[1]);
        
//         if(mc_daughter0.GetField<int>(1)!=2212 && mc_daughter0.GetField<int>(1)!=-211) continue;
//         if(mc_daughter1.GetField<int>(1)!=2212 && mc_daughter1.GetField<int>(1)!=-211) continue;

        const int rec_ids[2] = {matching->GetMatchInverted(daughters_ids[0]), matching->GetMatchInverted(daughters_ids[1])};

        if(rec_ids[0] < 0 || rec_ids[0] >= rec_tracks->GetNumberOfChannels())
          h3_not_rec_tracks.Fill( mc_daughter0.GetPhi(), mc_daughter0.GetPt(), mc_daughter0.GetRapidity(mc_daughter0.GetField<int>(1)) );
        else
          h3_rec_tracks.Fill( mc_daughter0.GetPhi(), mc_daughter0.GetPt(), mc_daughter0.GetRapidity(mc_daughter0.GetField<int>(1)) );

        if(rec_ids[1] < 0 || rec_ids[1] >= rec_tracks->GetNumberOfChannels())
          h3_not_rec_tracks.Fill( mc_daughter1.GetPhi(), mc_daughter1.GetPt(), mc_daughter1.GetRapidity(mc_daughter1.GetField<int>(1)) );
        else
          h3_rec_tracks.Fill( mc_daughter1.GetPhi(), mc_daughter1.GetPt(), mc_daughter1.GetRapidity(mc_daughter1.GetField<int>(1)) );
        
//         std::cout << rec_ids[0] << "\t" << rec_ids[1] << "\n";
        if (rec_ids[0] < 0 || rec_ids[1] < 0) continue;
        if (rec_ids[0] >= rec_tracks->GetNumberOfChannels() || rec_ids[1] >= rec_tracks->GetNumberOfChannels()) continue;

//         std::cout << "rec: " << rec_ids[0] << "  " << rec_ids[1] << std::endl;

        const auto& daughter0 = rec_tracks->GetChannel(rec_ids[0]);
        const auto& daughter1 = rec_tracks->GetChannel(rec_ids[1]);

//         const auto daughter0_mom = mc_daughter0.GetField<int>(0) < 0 ? daughter0.GetMomentum(0.14f) : daughter0.GetMomentum(0.938f);  // SHINE
//         const auto daughter1_mom = mc_daughter1.GetField<int>(0) < 0 ? daughter1.GetMomentum(0.14f) : daughter1.GetMomentum(0.938f);  // SHINE
        const auto daughter0_mom = mc_daughter0.GetField<int>(1) < 0 ? daughter0.GetMomentum(0.14f) : daughter0.GetMomentum(0.938f);  // CBM
        const auto daughter1_mom = mc_daughter1.GetField<int>(1) < 0 ? daughter1.GetMomentum(0.14f) : daughter1.GetMomentum(0.938f);  // CBM
        
//         if ( !(track_cuts->Apply(daughter0) && track_cuts->Apply(daughter1)) ) continue;
        
        if(!(  (daughter0.GetField<int>(1)==1 && daughter1.GetField<int>(1)==-1 ) || (daughter0.GetField<int>(1)==-1 && daughter1.GetField<int>(1)==1 ) )) continue;

        const auto mc_daughter0_mom = mc_daughter0.GetMomentum(mc_daughter0.GetField<int>(1));
        const auto mc_daughter1_mom = mc_daughter1.GetMomentum(mc_daughter1.GetField<int>(1));

        const auto mc_mother_mom = mc_daughter0_mom + mc_daughter1_mom;
        const auto mother_mom = daughter0_mom + daughter1_mom;

//         std::cout << "mass: " << mother_mom.M() << "  " << mc_mother_mom.M() << std::endl;
        inv_mass.Fill(mother_mom.M());
        number_of_lambdas++;
      }
    }
//     std::cout << i_event << "\t" << number_of_lambdas << "\n";
  }
  auto *out_file = TFile::Open("out.root", "recreate");
  inv_mass.Write("mass");
  h3_rec_tracks.Write("h3_rec_tracks");
  h3_not_rec_tracks.Write("h3_not_rec_tracks");
  out_file->Close();
}

std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks)
{
  std::vector<int> daughters_ids{};
  const int lambda_id = mc_track.GetId();
  for(int id=lambda_id; id<sim_tracks->GetNumberOfChannels(); ++id)
  {
    const auto& daughter = sim_tracks->GetChannel(id);
//     if( daughter.GetField<int>(2) == lambda_id )          // SHINE
    if( daughter.GetField<int>(0) == lambda_id )          // CBM
      daughters_ids.emplace_back(daughter.GetId());
  }

  if(daughters_ids.size() != 2 && daughters_ids.size() != 0){
    std::cout << "number of daughters = " << daughters_ids.size() << std::endl;
    for(int iD=0; iD<daughters_ids.size(); iD++)
      std::cout << sim_tracks->GetChannel(daughters_ids[iD]).GetField<int>(1) << std::endl;
  }
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
  
  AnalysisTree::SimpleCut chi2_ndf("chi2ndf", 0., 3.);

  AnalysisTree::Cuts* rec_track_cuts{ new AnalysisTree::Cuts("KfpfTracks")};
  rec_track_cuts->AddCuts({hits, hits_pot, chi2_ndf});
  return rec_track_cuts;
}
