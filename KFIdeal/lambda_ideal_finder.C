std::vector<int> FindDaughters(const AnalysisTree::Track& mc_track, const AnalysisTree::TrackDetector* sim_tracks);

void lambda_ideal_finder(const TString infile="/home/vklochkov/Data/na61/pbpb/13gev/mc/dcmqgsm_150.analysistree.root")
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("aTree");
  TH1F inv_mass("inv_mass", "", 100, 1.05, 1.15);

  TH3F h3_rec_tracks("rec_tracks", "", 50, -3.2, 3.2, 50, 0, 3, 50, 0.5, 4);
  TH3F h3_not_rec_tracks("not_rec_tracks", "", 50, -3.2, 3.2, 50, 0, 3, 50, -0.5, 4);

  auto* sim_tracks = new AnalysisTree::TrackDetector();
  auto* rec_tracks = new AnalysisTree::TrackDetector();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* matching = new AnalysisTree::Matching();
  auto* vtx_matching = new AnalysisTree::Matching();
  tree->SetBranchAddress("SimTracks", &sim_tracks);
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);
  tree->SetBranchAddress("VtxTracks", &vtx_tracks);
  tree->SetBranchAddress("KfpfracksToSimTracks", &matching);
  tree->SetBranchAddress("VtxTracksToSimTracks", &vtx_matching);

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

        const int rec_ids[2] = {matching->GetMatchInverted(daughters_ids[0]), matching->GetMatchInverted(daughters_ids[1])};

        if(rec_ids[0] < 0 || rec_ids[0] >= rec_tracks->GetNumberOfChannels())
          h3_not_rec_tracks.Fill( mc_daughter0.GetPhi(), mc_daughter0.GetPt(), mc_daughter0.GetRapidity(mc_daughter0.GetField<int>(1)) );
        else
          h3_rec_tracks.Fill( mc_daughter0.GetPhi(), mc_daughter0.GetPt(), mc_daughter0.GetRapidity(mc_daughter0.GetField<int>(1)) );

        if(rec_ids[1] < 0 || rec_ids[1] >= rec_tracks->GetNumberOfChannels())
          h3_not_rec_tracks.Fill( mc_daughter1.GetPhi(), mc_daughter1.GetPt(), mc_daughter1.GetRapidity(mc_daughter1.GetField<int>(1)) );
        else
          h3_rec_tracks.Fill( mc_daughter1.GetPhi(), mc_daughter1.GetPt(), mc_daughter1.GetRapidity(mc_daughter1.GetField<int>(1)) );

        if (rec_ids[0] < 0 || rec_ids[1] < 0) continue;
        if (rec_ids[0] >= rec_tracks->GetNumberOfChannels() || rec_ids[1] >= rec_tracks->GetNumberOfChannels()) continue;

        std::cout << "rec: " << rec_ids[0] << "  " << rec_ids[1] << std::endl;

        const auto& daughter0 = rec_tracks->GetChannel(rec_ids[0]);
        const auto& daughter1 = rec_tracks->GetChannel(rec_ids[1]);

        const auto daughter0_mom = daughter0.GetField<int>(1) < 0 ? daughter0.GetMomentum(0.14f) : daughter0.GetMomentum(0.938f);
        const auto daughter1_mom = daughter1.GetField<int>(1) < 0 ? daughter1.GetMomentum(0.14f) : daughter1.GetMomentum(0.938f);

        const auto mc_daughter0_mom = mc_daughter0.GetMomentum(mc_daughter0.GetField<int>(1));
        const auto mc_daughter1_mom = mc_daughter1.GetMomentum(mc_daughter1.GetField<int>(1));

        const auto mc_mother_mom = mc_daughter0_mom + mc_daughter1_mom;
        const auto mother_mom = daughter0_mom + daughter1_mom;

        std::cout << "mass: " << mother_mom.M() << "  " << mc_mother_mom.M() << std::endl;
        inv_mass.Fill(mother_mom.M());
      }
    }
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
    if( daughter.GetField<int>(2) == lambda_id )
      daughters_ids.emplace_back(daughter.GetId());
  }

  if(daughters_ids.size() != 2){
    std::cout << "number of daughters = " << daughters_ids.size() << std::endl;
  }
//  else
//    std::cout << "sim: " << daughters_ids[0] << "  " << daughters_ids[1] << std::endl;

  return daughters_ids;
}
