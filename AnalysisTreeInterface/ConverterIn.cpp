#include "ConverterIn.h"

#include "Interface/InputContainer.h"
#include "KFSimple/Constants.hpp"

#include "AnalysisTree/TaskManager.hpp"
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Matching.hpp>

using namespace AnalysisTree;

void ConverterIn::FillParticle(const AnalysisTree::Track& rec_particle, InputContainer& input_info) const {
  std::vector<float> mf(kNumberOfFieldPars, 0.f);
  for (int iF = 0; iF < kNumberOfFieldPars; iF++) {
    mf[iF] = rec_particle.GetField<float>(mf_field_id_ + iF);
  }
  auto cov_matrix = is_shine_ ? GetCovMatrixShine(rec_particle) : GetCovMatrixCbm(rec_particle);
  std::vector<float> par(kNumberOfTrackPars, 0.f);

  par[kX] = rec_particle.GetField<float>(par_field_id_);
  par[kY] = rec_particle.GetField<float>(par_field_id_ + 1);
  par[kZ] = rec_particle.GetField<float>(par_field_id_ + 2);
  par[kPx] = rec_particle.GetPx();
  par[kPy] = rec_particle.GetPy();
  par[kPz] = rec_particle.GetPz();

  const int pdg = rec_particle.GetField<int>(pdg_field_id_);//TODO
                                                            //  const int pdg = rec_particle.GetPid();

  input_info.AddTrack(par, cov_matrix, mf, rec_particle.GetField<int>(q_field_id_), pdg, rec_particle.GetId(), rec_particle.GetField<int>(nhits_field_id_));
}

void ConverterIn::Init() {
  auto* chain = AnalysisTree::TaskManager::GetInstance()->GetChain();

  rec_event_header_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(rec_event_header_name_));
  sim_event_header_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(sim_event_header_name_));
  kf_tracks_ = ANALYSISTREE_UTILS_GET<AnalysisTree::TrackDetector*>(chain->GetPointerToBranch(kf_tracks_name_));
  sim_tracks_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch(sim_tracks_name_));
  kf2sim_tracks_ = chain->GetMatchPointers().find(config_->GetMatchName(kf_tracks_name_, sim_tracks_name_))->second;

  const auto& branch_conf_kftr = config_->GetBranchConfig(kf_tracks_name_);
  q_field_id_ = branch_conf_kftr.GetFieldId("q");
  par_field_id_ = branch_conf_kftr.GetFieldId("x");   // par0
  mf_field_id_ = branch_conf_kftr.GetFieldId("cx0");  // magnetic field par0
  cov_field_id_ = branch_conf_kftr.GetFieldId("cov1");// cov matrix 0
  passcuts_field_id_ = branch_conf_kftr.GetFieldId("pass_cuts");
  pdg_field_id_ = branch_conf_kftr.GetFieldId("mc_pdg");
  nhits_field_id_ = branch_conf_kftr.GetFieldId("nhits");

  const auto& branch_conf_simtr = config_->GetBranchConfig(sim_tracks_name_);
  mother_id_field_id_ = branch_conf_simtr.GetFieldId("mother_id");
  sim_pdg_field_id_ = branch_conf_simtr.GetFieldId("pdg");

  if (track_cuts_) {
    track_cuts_->Init(*config_);
  }
}

InputContainer ConverterIn::CreateInputContainer() const {

  InputContainer input_container;
  const int n_tracks = kf_tracks_->GetNumberOfChannels();
  //  std::cout << " Ntracks = " << n_tracks << std::endl;
  input_container.SetDecay(decay_);
  input_container.SetCuts(cuts_);
  input_container.SetPV(rec_event_header_->GetVertexX(), rec_event_header_->GetVertexY(), rec_event_header_->GetVertexZ());

  int n_good_tracks{0};
  input_container.Reserve(n_tracks);
  for (int i_track = 0; i_track < n_tracks; ++i_track) {
    const auto& rec_track = kf_tracks_->GetChannel(i_track);
    if (!IsGoodTrack(rec_track)) continue;
    FillParticle(rec_track, input_container);
    n_good_tracks++;
  }
  //  std::cout << "Good tracks = " << n_good_tracks << "\n" << std::endl;
  return input_container;
}

std::vector<float> ConverterIn::GetCovMatrixCbm(const AnalysisTree::Track& particle) const {
  const auto tx = particle.GetField<float>(par_field_id_ + 3);
  const auto ty = particle.GetField<float>(par_field_id_ + 4);
  const auto qp = particle.GetField<float>(par_field_id_ + 5);
  const auto q = particle.GetField<int>(q_field_id_);

  //calculate covariance matrix
  const auto t = sqrt(1.f + tx * tx + ty * ty);
  const auto t3 = t * t * t;
  const auto dpxdtx = q / qp * (1.f + ty * ty) / t3;
  const auto dpxdty = -q / qp * tx * ty / t3;
  const auto dpxdqp = -q / (qp * qp) * tx / t;
  const auto dpydtx = -q / qp * tx * ty / t3;
  const auto dpydty = q / qp * (1.f + tx * tx) / t3;
  const auto dpydqp = -q / (qp * qp) * ty / t;
  const auto dpzdtx = -q / qp * tx / t3;
  const auto dpzdty = -q / qp * ty / t3;
  const auto dpzdqp = -q / (qp * qp * t3);

  const float F[kNumberOfTrackPars][5] = {{1.f, 0.f, 0.f, 0.f, 0.f},
                                          {0.f, 1.f, 0.f, 0.f, 0.f},
                                          {0.f, 0.f, 0.f, 0.f, 0.f},
                                          {0.f, 0.f, dpxdtx, dpxdty, dpxdqp},
                                          {0.f, 0.f, dpydtx, dpydty, dpydqp},
                                          {0.f, 0.f, dpzdtx, dpzdty, dpzdqp}};

  float VFT[5][kNumberOfTrackPars];
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < kNumberOfTrackPars; j++) {
      VFT[i][j] = 0;
      for (int k = 0; k < 5; k++) {
        VFT[i][j] += particle.GetField<float>(cov_field_id_ + std::min(i, k) + std::max(i, k) * (std::max(i, k) + 1) / 2) * F[j][k];
      }
    }
  }

  std::vector<float> cov(21, 0);
  for (int i = 0, l = 0; i < kNumberOfTrackPars; i++) {
    for (int j = 0; j <= i; j++, l++) {
      cov.at(l) = 0;
      for (int k = 0; k < 5; k++) {
        cov.at(l) += F[i][k] * VFT[k][j];
      }
    }
  }
  return cov;
}

std::vector<float> ConverterIn::GetCovMatrixShine(const AnalysisTree::Track& particle) const {
  std::vector<float> cov(21, 0.);
  for (int iCov = 0; iCov < 21; ++iCov) {
    cov[iCov] = particle.GetField<float>(cov_field_id_ + iCov);
  }
  return cov;
}

bool ConverterIn::IsGoodTrack(const AnalysisTree::Track& rec_track) const {
  return track_cuts_ ? track_cuts_->Apply(rec_track) : true;
}
