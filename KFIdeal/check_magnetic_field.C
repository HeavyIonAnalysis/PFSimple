#include "../KFSimple/Constants.h"

float signum(float x)
{
  if(x>0)
    return 1.;
  else
    return -1.;
}

void check_magnetic_field(const TString& infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/dcmqgsm_150.analysistree_more.root",
//                           const TString& field_file_name="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/field_x.root",
//                           const TString& field_file_name="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/field_y.root",
                          const TString& field_file_name="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/field_z.root",
//                           const TString& field_name="field_x"
//                           const TString& field_name="field_y"
                          const TString& field_name="field_z"
                         )
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

  TProfile2D out("diff", "", 300, -150, 150, 120, -60, 60);

  float Range = 0.004;
  TH1F hMF_bin("MFbin", "", 1000, -Range, Range);
  TH1F hMF_Inter("MFInter", "", 1000, -Range, Range);
  TH1F hMF_Param("MFParam", "", 1000, -Range, Range);
  
  TH2F hInt_Par("IntPar", "", 1000, -Range, Range, 1000, -Range, Range);
  hInt_Par.GetXaxis()->SetTitle("Interpolation");
  hInt_Par.GetYaxis()->SetTitle("Parametrization");
  
  const int n_entries = tree->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    tree->GetEntry(i_event);
    const int n_tracks = rec_tracks->GetNumberOfChannels();

    for(int i_track=0; i_track<n_tracks; ++i_track)
    {
      const auto& track = rec_tracks->GetChannel(i_track);

      const float z0 = track.GetField<float>(i_field_coff + kZ0);
      
//       const float cy0 = track.GetField<float>(i_field_coff + kCx0);
//       const float cy1 = track.GetField<float>(i_field_coff + kCx1);
//       const float cy2 = track.GetField<float>(i_field_coff + kCx2);
      
//       const float cy0 = track.GetField<float>(i_field_coff + kCy0);
//       const float cy1 = track.GetField<float>(i_field_coff + kCy1);
//       const float cy2 = track.GetField<float>(i_field_coff + kCy2);
      
      const float cy0 = track.GetField<float>(i_field_coff + kCz0);
      const float cy1 = track.GetField<float>(i_field_coff + kCz1);
      const float cy2 = track.GetField<float>(i_field_coff + kCz2);
      
      const float x = track.GetField<float>(i_par);
      const float y = track.GetField<float>(i_par+1);
      const float z = track.GetField<float>(i_par+2);
      
//       // NO interpolation
// 
//       const float mf_map = (field_map->GetBinContent(field_map->FindBin(fabs(x), fabs(y), fabs(z-40.))))*signum(x)*signum(y);  //X
//       const float mf_map = field_map->GetBinContent(field_map->FindBin(fabs(x), fabs(y), fabs(z-40.)));  //Y
//       const float mf_map = (field_map->GetBinContent(field_map->FindBin(fabs(x), fabs(y), fabs(z-40.))))*signum(y)*signum(z-40.);  //Z
      
//       hMF_bin.Fill((field_map->GetBinContent(field_map->FindBin(x, y, z)));    //X
      hMF_bin.Fill(field_map->GetBinContent(field_map->FindBin(x, y, z)));        //Y
//       hMF_bin.Fill((field_map->GetBinContent(field_map->FindBin(x, y, z)));    //Z
      
      // BEGIN root interpolation
      
//       const float x_l = fabs(x);
//       const float y_l = fabs(y);
//       const float z_l = fabs(z-40.);
//       if(x_l<=0 || x_l>=300. || y_l<=0 || y_l>=300. || z_l<=0. || z_l>=500.) continue;
//       if(x_l!=x_l || y_l!=y_l || z_l!=z_l) continue;
      
//       const float mf_map = field_map -> Interpolate(x, y, z);       //X
      const float mf_map = field_map -> Interpolate(x, y, z);          //Y
//       const float mf_map = field_map -> Interpolate(x, y, z);       //Z
           
      // END root interpolation
      
      hMF_Inter.Fill(mf_map);
      
      const float mf_y0 = 0.1*(cy0 + cy1*(z-z0) + cy2*(z-z0)*(z-z0));
      hMF_Param.Fill(mf_y0);
      hInt_Par.Fill(mf_map, mf_y0);

      out.Fill(x, y, fabs(mf_map-mf_y0) / fabs(mf_map+mf_y0) );
//      std::cout << mf_map << " " << mf_y0 << std::endl;
      
    }
  }

  auto *out_file = TFile::Open("out.root", "recreate");
  out.Write();
  hMF_bin.Write();
  hMF_Inter.Write();
  hMF_Param.Write();
  hInt_Par.Write();
  out_file->Close();
}
