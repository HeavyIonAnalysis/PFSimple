#include "../KFSimple/Constants.h"

float signum(float x)
{
  if(x>0)
    return 1.;
  else
    return -1.;
}

void check_magnetic_field(const TString& infile="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/reference.aTree.root",
                          const TString& field_file_name="/home/user/cbmdir/kfpf/kfpf_analysis_tree_converter/input/field_mapF.root",
//                           const TString& field_name="histoBx"
//                           const TString& field_name="histoBy"
                          const TString& field_name="histoBz"
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

  TProfile2D out("diff", "", 50, -50, 50, 50, -50, 50);
  TH1F histodeltaX("deltaX", "", 1000, -0.1, 1.1);
  TH1F histodeltaY("deltaY", "", 1000, -0.1, 1.1);
  TH1F histodeltaZ("deltaZ", "", 1000, -0.1, 1.1);
  
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
      
      // BEGIN interpolation ##################################################################################
      const float stepX = 2., stepY = 2., stepZ = 2.;
      const float x_l = fabs(x);
      const float y_l = fabs(y);
      const float z_l = fabs(z-40.);
      
      if(x_l<0 || x_l>300. || y_l<0 || y_l>300. || z_l<0. || z_l>500.) continue;
      
      float delta_X = (x_l - field_map->GetXaxis()->GetBinCenter(field_map->GetXaxis()->FindBin(x_l)))/stepX;
      float delta_Y = (y_l - field_map->GetYaxis()->GetBinCenter(field_map->GetYaxis()->FindBin(y_l)))/stepY;
      float delta_Z = (z_l - field_map->GetZaxis()->GetBinCenter(field_map->GetZaxis()->FindBin(z_l)))/stepZ;
      
      if(delta_X < 0.)
        delta_X = 1. + delta_X;
      if(delta_Y < 0.)
        delta_Y = 1. + delta_Y;
      if(delta_Z < 0.)
        delta_Z = 1. + delta_Z;
      
//       if(z>1.)
      {
        histodeltaX.Fill(delta_X);
        histodeltaY.Fill(delta_Y);
        histodeltaZ.Fill(delta_Z);
      }
      
      float Ha[2][2][2];
      
      
      Ha[0][0][0] = field_map -> GetBinContent(field_map->FindBin(x_l-stepX/2., y_l-stepY/2., z_l-stepZ/2.));
      Ha[0][0][1] = field_map -> GetBinContent(field_map->FindBin(x_l-stepX/2., y_l-stepY/2., z_l+stepZ/2.));
      Ha[0][1][0] = field_map -> GetBinContent(field_map->FindBin(x_l-stepX/2., y_l+stepY/2., z_l-stepZ/2.));
      Ha[0][1][1] = field_map -> GetBinContent(field_map->FindBin(x_l-stepX/2., y_l+stepY/2., z_l+stepZ/2.));
      Ha[1][0][0] = field_map -> GetBinContent(field_map->FindBin(x_l+stepX/2., y_l-stepY/2., z_l-stepZ/2.));
      Ha[1][0][1] = field_map -> GetBinContent(field_map->FindBin(x_l+stepX/2., y_l-stepY/2., z_l+stepZ/2.));
      Ha[1][1][0] = field_map -> GetBinContent(field_map->FindBin(x_l+stepX/2., y_l+stepY/2., z_l-stepZ/2.));
      Ha[1][1][1] = field_map -> GetBinContent(field_map->FindBin(x_l+stepX/2., y_l+stepY/2., z_l+stepZ/2.));      
        
      
            // x-axis interpolation
      float Hb[2][2];
      Hb[0][0] = Ha[0][0][0]*(1.-delta_X) + Ha[1][0][0]*delta_X;
      Hb[0][1] = Ha[0][0][1]*(1.-delta_X) + Ha[1][0][1]*delta_X;
      Hb[1][0] = Ha[0][1][0]*(1.-delta_X) + Ha[1][1][0]*delta_X;
      Hb[1][1] = Ha[0][1][1]*(1.-delta_X) + Ha[1][1][1]*delta_X;
      
      // y-axis interpolation
      float Hc[2];
      Hc[0] = Hb[0][0]*(1.-delta_Y) + Hb[1][0]*delta_Y;
      Hc[1] = Hb[0][1]*(1.-delta_Y) + Hb[1][1]*delta_Y;
      
//       const float mf_map = (Hc[0]*(1.-delta_Z) + Hc[1]*delta_Z)*signum(x)*signum(y);         //X
//       const float mf_map = Hc[0]*(1.-delta_Z) + Hc[1]*delta_Z;                               //Y
      const float mf_map = (Hc[0]*(1.-delta_Z) + Hc[1]*delta_Z)*signum(y)*signum(z-40.);     //Z

      // END interpolation ####################################################################################
      
      
      const float mf_y0 = 0.1*(cy0 + cy1*(z-z0) + cy2*(z-z0)*(z-z0));

      out.Fill(x, y, fabs(mf_map-mf_y0) / fabs(mf_map+mf_y0) );
//      std::cout << mf_map << " " << mf_y0 << std::endl;
      
    }
  }

  auto *out_file = TFile::Open("out.root", "recreate");
  out.Write();
  histodeltaX.Write();
  histodeltaY.Write();
  histodeltaZ.Write();
  out_file->Close();
}
