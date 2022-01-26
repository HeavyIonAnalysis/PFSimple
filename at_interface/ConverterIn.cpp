#include "ConverterIn.hpp"

#include "Constants.hpp"
#include "InputContainer.hpp"

#include "AnalysisTree/TaskManager.hpp"
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Matching.hpp>

using namespace AnalysisTree;

void ConverterIn::FillParticle(const AnalysisTree::BranchChannel& rec_particle) {
  
  std::vector<float> mf(NumberOfFieldPars, 0.f);
  for(int i=0; i<NumberOfFieldPars; i++)
    mf.at(i) = rec_particle[mf_field_.at(i)];

  auto cov_matrix = is_shine_ ? GetCovMatrixShine(rec_particle) : GetCovMatrixCbm(rec_particle);
  
  std::vector<float> par(kNumberOfTrackPars, 0.f);
  par.at(kX) = rec_particle[x_field_];
  par.at(kY) = rec_particle[y_field_];
  par.at(kZ) = rec_particle[z_field_];
  par.at(kPx) = rec_particle[px_field_];
  par.at(kPy) = rec_particle[py_field_];
  par.at(kPz) = rec_particle[pz_field_];

  const int q = rec_particle[q_field_];

  int pdg = -999;
  if (pid_mode_ == 0) {
    pdg = rec_particle[q_field_];
    container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[nhits_field_]);
  } else if (pid_mode_ == 1) {
    pdg = rec_particle[mc_pdg_field_];
    container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[nhits_field_]);
  } else {
    if (rec_particle[prob_p_field_] == -1)// needs to be done to exclude tracks with no TOF id (tracks with no TOF id have the same pid than negative background)
      return;

    if (pid_mode_ == 2 && pid_purity_.at(0) == 0.5) {
      const int pdg = rec_particle[rec_pdg_field_] * q;    // Be careful if use electrons and muons, because then the sign of pdg and charge do not match
      container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[nhits_field_]);

    } else {

      std::vector<float> pdg_prob;
      pdg_prob.push_back(rec_particle[prob_p_field_]); 
      pdg_prob.push_back(rec_particle[prob_pi_field_]);
      pdg_prob.push_back(rec_particle[prob_K_field_]); 
      pdg_prob.push_back(rec_particle[prob_d_field_]);
      pdg_prob.push_back(rec_particle[prob_bg_field_]);     

      if (pid_mode_ == 2) {
        if (*std::max_element(pdg_prob.begin(), pdg_prob.end()) < pid_purity_.at(0))
          return;
        auto it_prob = std::max_element(pdg_prob.begin(), pdg_prob.end());
        int ipid = std::distance(pdg_prob.begin(), it_prob);
        pdg = pid_codes_rec[ipid] * q;
        container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[nhits_field_]);
      }

      if (pid_mode_ == 3) {
        for (size_t ipid = 0; ipid < pid_codes_rec.size(); ipid++)
          if (pdg_prob[ipid] >= pid_purity_.at(ipid)) {
            pdg = pid_codes_rec[ipid] * q;
            container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[nhits_field_]);
          }
      }
    }
  }
}
void ConverterIn::Init() {
  auto* chain = AnalysisTree::TaskManager::GetInstance()->GetChain();

  if (pid_mode_ > 1) kf_tracks_name_ = "RecTracks";

  rec_event_header_ = chain->GetBranch(rec_event_header_name_);
  sim_event_header_ = chain->GetBranch(sim_event_header_name_);
  kf_tracks_ = chain->GetBranch(kf_tracks_name_);
  sim_tracks_ = chain->GetBranch(sim_tracks_name_);
  kf2sim_tracks_ = chain->GetMatching(kf_tracks_name_, sim_tracks_name_);
  
  for(auto& mf_comp : {"cx0", "cx1", "cx2", "cy0", "cy1", "cy2", "cz0", "cz1", "cz2", "z0"})
    mf_field_.push_back(kf_tracks_.GetField(mf_comp));
    
  x_field_ = kf_tracks_.GetField("x");
  y_field_ = kf_tracks_.GetField("y");  
  z_field_ = kf_tracks_.GetField("z");  
  px_field_ = kf_tracks_.GetField("px");  
  py_field_ = kf_tracks_.GetField("py");  
  pz_field_ = kf_tracks_.GetField("pz");  
    
  q_field_ = kf_tracks_.GetField("q");  
  mc_pdg_field_ = kf_tracks_.GetField("mc_pdg");  
  nhits_field_ = kf_tracks_.GetField("nhits");
    
  if(pid_mode_>1) {
    rec_pdg_field_ = kf_tracks_.GetField("pid");  
    prob_p_field_ = kf_tracks_.GetField("prob_p");  
    prob_pi_field_ = kf_tracks_.GetField("prob_pi");  
    prob_K_field_ = kf_tracks_.GetField("prob_K");  
    prob_d_field_ = kf_tracks_.GetField("prob_d");  
    prob_bg_field_ = kf_tracks_.GetField("prob_bg"); }
  
  vtx_x_field_ = rec_event_header_.GetField("vtx_x");
  vtx_y_field_ = rec_event_header_.GetField("vtx_y");
  vtx_z_field_ = rec_event_header_.GetField("vtx_z");
  
  tx_field_ = kf_tracks_.GetField("tx");
  ty_field_ = kf_tracks_.GetField("ty");
  qp_field_ = kf_tracks_.GetField("qp");
  
  int Ncov = is_shine_ ? 21 : 15;
  
  for(int i=0; i<Ncov; i++)
    cov_field_.push_back(kf_tracks_.GetField(("cov" + std::to_string(i+1)).c_str()));
  
  mother_id_field_ = sim_tracks_.GetField("mother_id");
  sim_pdg_field_ = sim_tracks_.GetField("pid");

  if (track_cuts_) {
    track_cuts_->Init(*config_);
  }
}

void ConverterIn::Exec() {

  container_ = InputContainer();
  const int n_tracks = kf_tracks_.size();

  container_.SetPV(rec_event_header_[0][vtx_x_field_],
                   rec_event_header_[0][vtx_y_field_],
                   rec_event_header_[0][vtx_z_field_]);

  int n_good_tracks{0};
  container_.Reserve(n_tracks);
  for (int i_track = 0; i_track < n_tracks; ++i_track) {
    const auto& rec_track = kf_tracks_[i_track];
    if (!IsGoodTrack(rec_track)) continue;
    if (!CheckMotherPdgs(rec_track)) continue;
    FillParticle(rec_track);
    n_good_tracks++;
  }
}

std::vector<float> ConverterIn::GetCovMatrixCbm(const AnalysisTree::BranchChannel& particle) const {
  const float tx = particle[tx_field_];
  const float ty = particle[ty_field_];
  const float qp = particle[qp_field_];
  const float q =  particle[q_field_];

  //calculate covariance matrix
  const float t = sqrt(1.f + tx * tx + ty * ty);
  const float t3 = t * t * t;
  const float dpxdtx = q / qp * (1.f + ty * ty) / t3;
  const float dpxdty = -q / qp * tx * ty / t3;
  const float dpxdqp = -q / (qp * qp) * tx / t;
  const float dpydtx = -q / qp * tx * ty / t3;
  const float dpydty = q / qp * (1.f + tx * tx) / t3;
  const float dpydqp = -q / (qp * qp) * ty / t;
  const float dpzdtx = -q / qp * tx / t3;
  const float dpzdty = -q / qp * ty / t3;
  const float dpzdqp = -q / (qp * qp * t3);

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
        VFT[i][j] += particle[cov_field_.at(std::min(i, k) + std::max(i, k) * (std::max(i, k) + 1) / 2)] * F[j][k];
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

std::vector<float> ConverterIn::GetCovMatrixShine(const AnalysisTree::BranchChannel& particle) const {
  std::vector<float> cov(21, 0.);
  for (int iCov = 0; iCov < 21; ++iCov) {
    cov[iCov] = particle[cov_field_.at(iCov)];
  }
  return cov;
}

bool ConverterIn::IsGoodTrack(const AnalysisTree::BranchChannel& rec_track) const {
  return track_cuts_ ? track_cuts_->Apply(rec_track) : true;
}

bool ConverterIn::CheckMotherPdgs(const AnalysisTree::BranchChannel& rec_track) const {
  if (mother_pdgs_to_be_considered_.size() == 0)
    return true;

//   if (!sim_tracks_ || !kf2sim_tracks_) {
//     std::cout << "No MC info available!\n";
//     assert(false);
//   }

  const int sim_id = kf2sim_tracks_->GetMatch(rec_track.GetId());
  if (sim_id < 0)
    return false;

  const AnalysisTree::BranchChannel& sim_track = sim_tracks_[sim_id];
  const int mother_id = sim_track[mother_id_field_];
  if (mother_id < 0)
    return false;

  const int mother_pdg = sim_tracks_[mother_id][sim_pdg_field_];

  bool ok = false;

  for (auto& good_mother_pdgs : mother_pdgs_to_be_considered_)
    if (mother_pdg == good_mother_pdgs) {
      ok = true;
      break;
    }

  return ok;
}
