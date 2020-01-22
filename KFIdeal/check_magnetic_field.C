#include "../KFSimple/Constants.h"

void check_magnetic_field(const TString& infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/reference.aTree.root",
                          const TString& field_file_name="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/field_mapF.root",
                          const TString& field_name="histoBy")
{
  std::unique_ptr<TFile> field_file {TFile::Open(field_file_name)};
  auto* field_map = field_file->Get<TH3F>(field_name);

  std::unique_ptr<TFile> file {TFile::Open(infile)};
  std::unique_ptr<TTree> tree {file->Get<TTree>("aTree")};

  AnalysisTree::Configuration *config = file->Get<AnalysisTree::Configuration>("Configuration");
  const int i_field_coff = config->GetBranchConfig("KfpfTracks").GetFieldId("cx0");
  const int i_par = config->GetBranchConfig("KfpfTracks").GetFieldId("x");

  auto* rec_tracks = new AnalysisTree::TrackDetector();
  tree->SetBranchAddress("KfpfTracks", &rec_tracks);

  TProfile2D out("diff", "", 50, -50, 50, 50, -50, 50);

  const int n_entries = tree->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    tree->GetEntry(i_event);
    const int n_tracks = rec_tracks->GetNumberOfChannels();

    for(int i_track=0; i_track<n_tracks; ++i_track)
    {
      const auto& track = rec_tracks->GetChannel(i_track);

      const float z0 = track.GetField<float>(i_field_coff + kZ0);
      const float cy0 = track.GetField<float>(i_field_coff + kCy0);
      const float cy1 = track.GetField<float>(i_field_coff + kCy1);
      const float cy2 = track.GetField<float>(i_field_coff + kCy2);
      const float x = track.GetField<float>(i_par);
      const float y = track.GetField<float>(i_par+1);
      const float z = track.GetField<float>(i_par+2);

//       const float mf_map = field_map->GetBinContent(field_map->FindBin(x, y, z));
      const float mf_map = field_map->GetBinContent(field_map->FindBin(fabs(x), fabs(y), z));
      const float mf_y0 = 0.1*(cy0 + cy1*(z-z0) + cy2*(z-z0)*(z-z0));

      out.Fill(x, y, (mf_map-mf_y0) / (mf_map+mf_y0) );
//      std::cout << mf_map << " " << mf_y0 << std::endl;
    }
  }

  auto *out_file = TFile::Open("out.root", "recreate");
  out.Write();
  out_file->Close();
}
