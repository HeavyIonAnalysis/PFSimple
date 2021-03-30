#include "ConverterOut.h"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"
#include "AnalysisTree/TaskManager.hpp"

void CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) {

  particle.SetMomentum(kf_particle.GetPx(), kf_particle.GetPy(), kf_particle.GetPz());
  particle.SetMass(kf_particle.GetMass());

  //  float mass, masserr;
  //  particle.GetMass(mass, masserr);
  //  lambdarec.SetField(particle.DaughterIds()[0], daughter1_id_field_id_);
  //  lambdarec.SetField(particle.DaughterIds()[1], daughter2_id_field_id_);
  //  if (decay_.GetNdaughters() == 3)
  //    lambdarec.SetField(particle.DaughterIds()[2], daughter3_id_field_id_);
  //
  //  lambdarec.SetField(particle.X(), x_field_id_);
  //  lambdarec.SetField(particle.Y(), y_field_id_);
  //  lambdarec.SetField(particle.Z(), z_field_id_);
  //  lambdarec.SetMomentum(particle.GetPx(), particle.GetPy(), particle.GetPz());
  //
  //  lambdarec.SetField(particle.GetErrPx(), px_err_field_id_);
  //  lambdarec.SetField(particle.GetErrPy(), py_err_field_id_);
  //  lambdarec.SetField(particle.GetErrPz(), pz_err_field_id_);
  //
  //  lambdarec.SetMass(mass);
  //  lambdarec.SetField(masserr, mass_err_field_id_);
  //  lambdarec.SetPid(particle.GetPDG());
  //
  //  lambdarec.SetField(candidate.GetChi2PrimPos(), chi2primpos_field_id_);
  //  lambdarec.SetField(candidate.GetChi2PrimNeg(), chi2primneg_field_id_);
  //  lambdarec.SetField(candidate.GetDistance(), distance_field_id_);
  //  lambdarec.SetField(candidate.GetCosineDaughterPos(), cosinepos_field_id_);
  //  lambdarec.SetField(candidate.GetCosineDaughterNeg(), cosineneg_field_id_);
  //  lambdarec.SetField(candidate.GetChi2Geo(), chi2geo_field_id_);
  //  lambdarec.SetField(candidate.GetL(), l_field_id_);
  //  lambdarec.SetField(candidate.GetLdL(), ldl_field_id_);
  //  lambdarec.SetField(candidate.GetIsFromPV(), isfrompv_field_id_);
  //  lambdarec.SetField(candidate.GetCosineTopo(), cosinetopo_field_id_);
  //  lambdarec.SetField(candidate.GetChi2Topo(), chi2topo_field_id_);
  //  lambdarec.SetField(candidate.GetNHitsPos(), nhits_pos_field_id_);
  //  lambdarec.SetField(candidate.GetNHitsNeg(), nhits_neg_field_id_);
}

void ConverterOut::Exec() {
  lambda_reco_->ClearChannels();
  lambda_sim_->ClearChannels();
  lambda_reco2sim_->Clear();

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  events_->SetField(float(sim_events_->GetField<float>(b_field_id_)), 0);//TODO

  for (const auto& candidate : candidates_) {
    auto& lambdarec = lambda_reco_->AddChannel(out_config->GetBranchConfig(lambda_reco_->GetId()));
    CopyParticle(candidate, lambdarec);
  }

  //  MatchWithMc();  //TODO
}

void ConverterOut::Init() {

  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();

  sim_events_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(sim_events_name_));
  mc_particles_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch(mc_particles_name_));
  rec_tracks_ = ANALYSISTREE_UTILS_GET<AnalysisTree::TrackDetector*>(chain->GetPointerToBranch(rec_tracks_name_));
  rec_to_mc_ = chain->GetMatchPointers().find(config_->GetMatchName(rec_tracks_name_, mc_particles_name_))->second;

  std::string out_branch_event = "Events";
  std::string out_branch = decay_.GetNameMother() + std::string("Candidates");
  std::string out_branch_sim = decay_.GetNameMother() + std::string("Simulated");
  std::string out_branch_reco2sim = out_branch + "2" + out_branch_sim;

  AnalysisTree::BranchConfig EventBranch(out_branch_event, AnalysisTree::DetType::kEventHeader);
  EventBranch.AddField<float>("b");

  AnalysisTree::BranchConfig out_particles(out_branch, AnalysisTree::DetType::kParticle);

  out_particles.AddField<float>("chi2_prim");
  out_particles.AddField<float>("distance");
  out_particles.AddField<float>("cosineneg");
  out_particles.AddField<float>("chi2_geo");
  out_particles.AddField<float>("l");
  out_particles.AddField<float>("l_over_dl");
  out_particles.AddField<float>("cosinetopo");
  out_particles.AddField<int>("nhits");

  out_particles.AddField<bool>("isfrompv");

  out_particles.AddFields<float>({"x", "y", "z"});
  out_particles.AddFields<int>({"daughter1_id", "daughter2_id", "daughter3_id"});
  out_particles.AddFields<float>({"pT_err", "phi_err", "eta_err", "mass_err"});

  if (mc_particles_) {
    out_particles.AddField<bool>("is_signal");
  }
  AnalysisTree::BranchConfig LambdaSimBranch(out_branch_sim, AnalysisTree::DetType::kParticle);

  man->AddBranch(out_branch_event, events_, EventBranch);
  man->AddBranch(out_branch, lambda_reco_, out_particles);
  man->AddBranch(out_branch_sim, lambda_sim_, LambdaSimBranch);
  man->AddMatching(out_branch, out_branch_sim, lambda_reco2sim_);

  events_->Init(EventBranch);
  InitIndexes();
}

void ConverterOut::MatchWithMc() {

  for (auto& lambdarec : *lambda_reco_) {

    const int simtrackid1 = rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter1_id_field_id_));
    const int simtrackid2 = rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter2_id_field_id_));

    bool is_signal = false;
    int mother_id = -999;
    if (!(simtrackid1 < 0 || simtrackid2 < 0)) {
      const AnalysisTree::Particle& simtrack1 = mc_particles_->GetChannel(simtrackid1);
      const AnalysisTree::Particle& simtrack2 = mc_particles_->GetChannel(simtrackid2);

      if (simtrack1.GetField<int>(mother_id_field_id_) == simtrack2.GetField<int>(mother_id_field_id_)) {

        mother_id = simtrack1.GetField<int>(mother_id_field_id_);
        if (mother_id < 0) continue;

        const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

        if (decay_.GetNdaughters() == 2) is_signal = simtrackmother.GetPid() == decay_.GetPdgMother();

        if (decay_.GetNdaughters() == 3) {
          const int simtrackid3 = rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter3_id_field_id_));
          if (simtrackid3 >= 0) {
            const AnalysisTree::Particle& simtrack3 = mc_particles_->GetChannel(simtrackid3);
            if (simtrack1.GetField<int>(mother_id_field_id_) == simtrack3.GetField<int>(mother_id_field_id_))
              is_signal = simtrackmother.GetPid() == decay_.GetPdgMother();
          }
        }
      }
    }

    lambdarec.SetField(is_signal, is_signal_field_id_);
    auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

    if (is_signal) {
      const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

      auto& lambdasim = lambda_sim_->AddChannel(out_config->GetBranchConfig(lambda_sim_->GetId()));

      lambdasim.SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
      lambdasim.SetMass(simtrackmother.GetMass());
      lambdasim.SetPid(simtrackmother.GetPid());
      lambda_reco2sim_->AddMatch(lambdarec.GetId(), lambdasim.GetId());
    }
  }
}

void ConverterOut::InitIndexes() {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  const auto& out_branch = out_config->GetBranchConfig(lambda_reco_->GetId());

  auto branch_conf_sim_event = config_->GetBranchConfig(sim_events_name_);
  b_field_id_ = branch_conf_sim_event.GetFieldId("b");

  x_field_id_ = out_branch.GetFieldId("x");
  daughter1_id_field_id_ = out_branch.GetFieldId("daughter1_id");
  px_err_field_id_ = out_branch.GetFieldId("pT_err");

  //  if (mc_particles_) {
  //    auto branch_conf_sim = config_->GetBranchConfig(mc_particles_name_);
  //    mother_id_field_id_ = branch_conf_sim.GetFieldId("mother_id");
  //    is_signal_field_id_ = out_branch.GetFieldId("is_signal");
  //  }

  //  chi2primpos_field_id_ = out_branch.GetFieldId("chi2primpos");
  //  chi2primneg_field_id_ = out_branch.GetFieldId("chi2primneg");
  //  distance_field_id_ = out_branch.GetFieldId("distance");
  //  cosinepos_field_id_ = out_branch.GetFieldId("cosinepos");
  //  cosineneg_field_id_ = out_branch.GetFieldId("cosineneg");
  //  chi2geo_field_id_ = out_branch.GetFieldId("chi2geo");
  //  l_field_id_ = out_branch.GetFieldId("l");
  //  ldl_field_id_ = out_branch.GetFieldId("ldl");
  //  isfrompv_field_id_ = out_branch.GetFieldId("isfrompv");
  //  cosinetopo_field_id_ = out_branch.GetFieldId("cosinetopo");
  //  chi2topo_field_id_ = out_branch.GetFieldId("chi2topo");
  //  nhits_pos_field_id_ = out_branch.GetFieldId("nhitspos");
  //  nhits_neg_field_id_ = out_branch.GetFieldId("nhitsneg");
}
