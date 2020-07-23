#include "InConverter.h"

//KF Particle headers
#include "KFParticle.h"
#include "KFPTrackVector.h"

//KF Simple headers
#include "SimpleFinder.h"
#include "Constants.h"

//c++ and std headers
#include <iostream>
#include <cmath>
#include <vector>

void InConverter::InitAnalysisTree(const std::string& file_name, const std::string& tree_name)
{
  in_file_ = TFile::Open(file_name.c_str(), "read");
  config_ = (AnalysisTree::Configuration*) in_file_->Get("Configuration");

  in_chain_ =  new TChain(tree_name.c_str());
  in_chain_ -> Add(file_name.c_str());
  in_chain_ -> SetBranchAddress("KfpfTracks", &kf_tracks_);
  in_chain_ -> SetBranchAddress("SimTracks", &sim_tracks_);
  in_chain_ -> SetBranchAddress("RecEventHeader", &rec_event_header_);
  in_chain_ -> SetBranchAddress("SimEventHeader", &sim_event_header_);
  in_chain_ -> SetBranchAddress(config_->GetMatchName("KfpfTracks", "SimTracks").c_str(), &kf2sim_tracks_);

  auto branch_conf_kftr = config_->GetBranchConfig( "KfpfTracks" );
  q_field_id_ = branch_conf_kftr.GetFieldId("q");

  par_field_id_ = branch_conf_kftr.GetFieldId("x");   // par0
  mf_field_id_ = branch_conf_kftr.GetFieldId("cx0");  // magnetic field par0
  cov_field_id_ = branch_conf_kftr.GetFieldId("cov1"); // cov matrix 0
  
  passcuts_field_id_ = branch_conf_kftr.GetFieldId("pass_cuts");
  pdg_field_id_ = branch_conf_kftr.GetFieldId("mc_pdg");

  auto branch_conf_simtr = config_->GetBranchConfig( "SimTracks" );
  mother_id_field_id_ = branch_conf_simtr.GetFieldId("mother_id");
  sim_pdg_field_id_ = branch_conf_simtr.GetFieldId("pdg");

  if (track_cuts_ != nullptr)
    track_cuts_->Init(*config_);
}

int InConverter::GetTrackPid(const AnalysisTree::Track& rec_track)
{
  int pdg = -1;
  if(is_shine_) {
    pdg = rec_track.GetField<int>(q_field_id_) > 0 ? 2212 : -211;
  }
  else {
    pdg = rec_track.GetField<int>(pdg_field_id_);
  }
  return pdg;
}

void InConverter::FillTrack(const AnalysisTree::Track& rec_track, InputContainer& input_info)
{
  std::vector <float> mf(kNumberOfFieldPars, 0.f);
  for(int iF=0; iF<kNumberOfFieldPars; iF++)
    mf[iF] = rec_track.GetField<float>(mf_field_id_+iF);

  auto cov_matrix = is_shine_ ? GetCovMatrixShine(rec_track) : GetCovMatrixCbm(rec_track);
  std::vector <float> par(kNumberOfTrackPars,0.f);

  par[kX] = rec_track.GetField<float>(par_field_id_);
  par[kY] = rec_track.GetField<float>(par_field_id_ + 1);
  par[kZ] = rec_track.GetField<float>(par_field_id_ + 2);
  par[kPx] = rec_track.GetPx();
  par[kPy] = rec_track.GetPy();
  par[kPz] = rec_track.GetPz();

  const int pdg = GetTrackPid(rec_track);

  input_info.AddTrack(par, cov_matrix, mf, rec_track.GetField<int>(q_field_id_), pdg, rec_track.GetId(), rec_track.GetField<int>(passcuts_field_id_));
}


InputContainer InConverter::CreateInputContainer(int iEvent)
{
  in_chain_->GetEntry(iEvent);
  const int n_tracks = kf_tracks_->GetNumberOfChannels();
  std::cout << "Event # " << iEvent << " with Ntracks = " << n_tracks << std::endl;

  InputContainer input_info;
  input_info.SetCuts(cuts_);
  input_info.SetPV(rec_event_header_->GetVertexX(), rec_event_header_->GetVertexY(), rec_event_header_->GetVertexZ());

  int n_good_tracks{0};
  for(int i_track=0; i_track<n_tracks; ++i_track)
  {
// <<<<<<< HEAD
    const auto& rec_track = kf_tracks_->GetChannel(i_track);
    if(!IsGoodTrack(rec_track)) continue;
    FillTrack(rec_track, input_info);
    n_good_tracks++;
// =======
//     const AnalysisTree::Track& rec_track = kf_tracks_->GetChannel(i_track);
//     
//     if (track_cuts_ != nullptr)
//       if(!track_cuts_->Apply(rec_track))
//         continue;
// 
//     std::vector<float> par;
//     for(int iP=0; iP<3; iP++)
//       par.push_back(rec_track.GetField<float>(par_field_id_ + iP));
//     par.push_back(rec_track.GetPx());
//     par.push_back(rec_track.GetPy());
//     par.push_back(rec_track.GetPz());      
//           
//     std::vector<float> mf;
//     for(int iF=0; iF<10; iF++)
//       mf.push_back(rec_track.GetField<float>(mf_field_id_+iF));
//     
//     int pdg = -1;
//     if(is_shine_) {
//       pdg = rec_track.GetField<int>(q_field_id_) > 0 ? 2212 : -211;
//     }
//     else {
//       pdg = rec_track.GetField<int>(pdg_field_id_);
//     }
//     auto cov_matrix = is_shine_ ? GetCovMatrixShine(rec_track) : GetCovMatrixCbm(rec_track);
//     
//   //       if(ststrack.GetField<int>(m_pdg_field_id_)!=3122) continue;
//     inputInfo.AddTrack( par,
//                         cov_matrix,
//                         mf,
//                         rec_track.GetField<int>(q_field_id_),
//                         pdg,
//                         rec_track.GetId(),
//                         rec_track.GetField<int>(passcuts_field_id_));
// //                        ststrack.GetField<int>(nhits_field_id_),
// //                        ststrack.GetField<int>(passcuts_field_id_));
// >>>>>>> lubynets-order
  }

  std::cout << "Good tracks = " << n_good_tracks << std::endl;
  return input_info;
}

bool InConverter::IsGoodTrack(const AnalysisTree::Track& rec_track)
{
  return true;
  bool is_good{false};

  const int sim_id = kf2sim_tracks_->GetMatch(rec_track.GetId());
//    std::cout<< "sim  "  << sim_id << std::endl;
  if (sim_id >= 0 && sim_id < sim_tracks_->GetNumberOfChannels()) {
    const AnalysisTree::Track& sim_track = sim_tracks_->GetChannel(sim_id);
    const int mother_id = sim_track.GetField<int>(mother_id_field_id_);
//      std::cout << "mother " << mother_id << std::endl;

    if (mother_id >= 0 && mother_id < sim_tracks_->GetNumberOfChannels() )
    {
      const AnalysisTree::Track& mother_track = sim_tracks_->GetChannel(mother_id);
      const int mother_pdg = mother_track.GetField<int>(sim_pdg_field_id_);
      std::cout << "mother pdg " << mother_pdg << std::endl;

      if(mother_pdg == 3122)
        is_good = true;
    }
  }
  return is_good;
}

std::vector<float> InConverter::GetCovMatrixCbm(const AnalysisTree::Track& track) const
{
  const double tx = track.GetField<float>(par_field_id_+3); 
  const double ty = track.GetField<float>(par_field_id_+4); 
  const double qp = track.GetField<float>(par_field_id_+5);
  const Int_t q = track.GetField<int>(q_field_id_);
  
  //calculate covariance matrix
  const double t = sqrt(1.f + tx*tx + ty*ty);
  const double t3 = t*t*t;
  const double dpxdtx = q/qp*(1.f+ty*ty)/t3;
  const double dpxdty = -q/qp*tx*ty/t3;
  const double dpxdqp = -q/(qp*qp)*tx/t;
  const double dpydtx = -q/qp*tx*ty/t3;
  const double dpydty = q/qp*(1.f+tx*tx)/t3;
  const double dpydqp = -q/(qp*qp)*ty/t;
  const double dpzdtx = -q/qp*tx/t3;
  const double dpzdty = -q/qp*ty/t3;
  const double dpzdqp = -q/(qp*qp*t3);
  
  const double F[6][5] = { {1.f, 0.f, 0.f,    0.f,    0.f},
  {0.f, 1.f, 0.f,    0.f,    0.f},
  {0.f, 0.f, 0.f,    0.f,    0.f},
  {0.f, 0.f, dpxdtx, dpxdty, dpxdqp},
  {0.f, 0.f, dpydtx, dpydty, dpydqp},
  {0.f, 0.f, dpzdtx, dpzdty, dpzdqp} };
    
  double VFT[5][6];
  for(int i=0; i<5; i++)
    for(int j=0; j<6; j++)
    {
      VFT[i][j] = 0;
      for(int k=0; k<5; k++)
      {
        if (k<=i)
          VFT[i][j] += track.GetField<float>(cov_field_id_ + k + i*(i+1)/2) * F[j][k];   //parameters->GetCovariance(i,k) * F[j][k];
        else
          VFT[i][j] += track.GetField<float>(cov_field_id_ + i + k*(k+1)/2) * F[j][k];   //parameters->GetCovariance(i,k) * F[j][k];
      }
    }

    
  std::vector<float> cov(21, 0);
  for(int i=0, l=0; i<6; i++)
    for(int j=0; j<=i; j++, l++)
    {
      cov[l] = 0;
      for(int k=0; k<5; k++)
      {
        cov[l] += F[i][k] * VFT[k][j];
      }
    }
  return cov;
}

std::vector<float> InConverter::GetCovMatrixShine(const AnalysisTree::Track& track) const
{
  std::vector<float> cov(21, 0.);
  
  for (int iCov=0; iCov<21; ++iCov)
    cov[iCov] = track.GetField<float>(cov_field_id_+iCov);
  
  return cov;
}
