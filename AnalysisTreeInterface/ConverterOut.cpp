#include "ConverterOut.h"

#include "TTree.h"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"

void ConverterOut::Exec() {
  lambda_reco_->ClearChannels();
  lambda_sim_->ClearChannels();
  lambda_reco2sim_->Clear();

  events_->SetField(float(sim_events_->GetField<float>(b_field_id_)), out_config_->GetBranchConfig(events_->GetId()).GetFieldId("b"));

  for (const auto& candidate : canditates_) {
    auto* lambdarec = lambda_reco_->AddChannel();
    lambdarec->Init(out_config_->GetBranchConfig(lambda_reco_->GetId()));

    const KFParticle& particle = candidate.GetParticle();
    float mass, masserr;
    particle.GetMass(mass, masserr);
    lambdarec->SetField(particle.DaughterIds()[0], daughter1_id_field_id_);
    lambdarec->SetField(particle.DaughterIds()[1], daughter2_id_field_id_);
    if (decay_.GetNdaughters() == 3)
      lambdarec->SetField(particle.DaughterIds()[2], daughter3_id_field_id_);

    lambdarec->SetField(particle.X(), x_field_id_);
    lambdarec->SetField(particle.Y(), y_field_id_);
    lambdarec->SetField(particle.Z(), z_field_id_);
    lambdarec->SetMomentum(particle.GetPx(), particle.GetPy(), particle.GetPz());

    lambdarec->SetField(particle.GetErrPx(), px_err_field_id_);
    lambdarec->SetField(particle.GetErrPy(), py_err_field_id_);
    lambdarec->SetField(particle.GetErrPz(), pz_err_field_id_);

    lambdarec->SetMass(mass);
    lambdarec->SetField(masserr, mass_err_field_id_);
    lambdarec->SetPid(particle.GetPDG());

    lambdarec->SetField(candidate.GetChi2PrimPos(), chi2primpos_field_id_);
    lambdarec->SetField(candidate.GetChi2PrimNeg(), chi2primneg_field_id_);
    lambdarec->SetField(candidate.GetDistance(), distance_field_id_);
    lambdarec->SetField(candidate.GetCosineDaughterPos(), cosinepos_field_id_);
    lambdarec->SetField(candidate.GetCosineDaughterNeg(), cosineneg_field_id_);
    lambdarec->SetField(candidate.GetChi2Geo(), chi2geo_field_id_);
    lambdarec->SetField(candidate.GetL(), l_field_id_);
    lambdarec->SetField(candidate.GetLdL(), ldl_field_id_);
    lambdarec->SetField(candidate.GetIsFromPV(), isfrompv_field_id_);
    lambdarec->SetField(candidate.GetCosineTopo(), cosinetopo_field_id_);
    lambdarec->SetField(candidate.GetChi2Topo(), chi2topo_field_id_);
    lambdarec->SetField(candidate.GetNHitsPos(), nhits_pos_field_id_);
    lambdarec->SetField(candidate.GetNHitsNeg(), nhits_neg_field_id_);

    if (decay_.GetNdaughters() == 3) {
      lambdarec->SetField(candidate.GetChi2PrimThird(), chi2primthird_field_id_);
      lambdarec->SetField(candidate.GetDistanceThird(), distancethird_field_id_);
      lambdarec->SetField(candidate.GetCosineDaughterThird(), cosinethird_field_id_);
      lambdarec->SetField(candidate.GetChi2GeoThree(), chi2geothree_field_id_);
      lambdarec->SetField(candidate.GetCosineTopoThree(), cosinetopothree_field_id_);
      lambdarec->SetField(candidate.GetChi2TopoThree(), chi2topothree_field_id_);
      lambdarec->SetField(candidate.GetNHitsThird(), nhits_third_field_id_);
    }
  }

  MatchWithMc();
  out_tree_->Fill();
}

void ConverterOut::Init() {
  assert(out_config_ && out_tree_);
  std::string out_branch_event = "Events";
  out_branch_ = decay_.GetNameMother() + std::string("Candidates");
  std::string out_branch_sim = decay_.GetNameMother() + std::string("Simulated");
  std::string out_branch_reco2sim = decay_.GetNameMother() + std::string("Candidates2") + decay_.GetNameMother() + std::string("Simulated");

  if (!in_branches_.empty()) {
    sim_events_ = (AnalysisTree::EventHeader*) branches.find(in_branches_[2])->second;
    mc_particles_ = (AnalysisTree::Particles*) branches.find(in_branches_[0])->second;
    rec_tracks_ = (AnalysisTree::TrackDetector*) branches.find(in_branches_[1])->second;
    const auto& match = config_->GetMatchName(in_branches_[1], in_branches_[0]);
    rec_to_mc_ = (AnalysisTree::Matching*) branches.find(match)->second;
  }

  AnalysisTree::BranchConfig EventBranch(out_branch_event, AnalysisTree::DetType::kEventHeader);
  EventBranch.AddField<float>("b");

  AnalysisTree::BranchConfig LambdaRecoBranch(out_branch_, AnalysisTree::DetType::kParticle);

  LambdaRecoBranch.AddField<float>("chi2primpos");
  LambdaRecoBranch.AddField<float>("chi2primneg");
  LambdaRecoBranch.AddField<float>("distance");
  LambdaRecoBranch.AddField<float>("cosinepos");
  LambdaRecoBranch.AddField<float>("cosineneg");
  LambdaRecoBranch.AddField<float>("chi2geo");
  LambdaRecoBranch.AddField<float>("l");
  LambdaRecoBranch.AddField<float>("ldl");
  LambdaRecoBranch.AddField<int>("isfrompv");
  LambdaRecoBranch.AddField<float>("cosinetopo");
  LambdaRecoBranch.AddField<float>("chi2topo");
  LambdaRecoBranch.AddField<int>("nhitspos");
  LambdaRecoBranch.AddField<int>("nhitsneg");

  if (decay_.GetNdaughters() == 3) {
    LambdaRecoBranch.AddField<float>("chi2primthird");
    LambdaRecoBranch.AddField<float>("distancethird");
    LambdaRecoBranch.AddField<float>("cosinethird");
    LambdaRecoBranch.AddField<float>("chi2geothree");
    LambdaRecoBranch.AddField<float>("cosinetopothree");
    LambdaRecoBranch.AddField<float>("chi2topothree");
    LambdaRecoBranch.AddField<int>("nhitsthird");
  }

  LambdaRecoBranch.AddField<float>("x");
  LambdaRecoBranch.AddField<float>("y");
  LambdaRecoBranch.AddField<float>("z");
  LambdaRecoBranch.AddField<int>("daughter1id");
  LambdaRecoBranch.AddField<int>("daughter2id");
  LambdaRecoBranch.AddField<int>("daughter3id");
  LambdaRecoBranch.AddField<float>("pxerr");
  LambdaRecoBranch.AddField<float>("pyerr");
  LambdaRecoBranch.AddField<float>("pzerr");
  LambdaRecoBranch.AddField<float>("masserr");

  if (mc_particles_) {
    LambdaRecoBranch.AddField<bool>("is_signal");
  }

  out_config_->AddBranchConfig(LambdaRecoBranch);
  lambda_reco_ = new AnalysisTree::Particles(out_config_->GetLastId());

  out_config_->AddBranchConfig(EventBranch);
  events_ = new AnalysisTree::EventHeader(out_config_->GetLastId());

  AnalysisTree::BranchConfig LambdaSimBranch(out_branch_sim, AnalysisTree::DetType::kParticle);
  out_config_->AddBranchConfig(LambdaSimBranch);
  lambda_sim_ = new AnalysisTree::Particles(out_config_->GetLastId());

  lambda_reco2sim_ = new AnalysisTree::Matching(out_config_->GetBranchConfig(out_branch_).GetId(), out_config_->GetBranchConfig(out_branch_sim).GetId());

  out_config_->AddMatch(lambda_reco2sim_);

  out_tree_->Branch(out_branch_event.c_str(), "AnalysisTree::EventHeader", &events_);
  out_tree_->Branch(out_branch_.c_str(), "AnalysisTree::Particles", &lambda_reco_);
  out_tree_->Branch(out_branch_sim.c_str(), "AnalysisTree::Particles", &lambda_sim_);
  out_tree_->Branch(out_branch_reco2sim.c_str(), "AnalysisTree::Matching", &lambda_reco2sim_);
  events_->Init(out_config_->GetBranchConfig(events_->GetId()));
  InitIndexes();
}

void ConverterOut::MatchWithMc() {

  //  for(auto& lambdarec : *lambda_reco_->Channels()) {
  for (int i_ch = 0; i_ch < lambda_reco_->GetNumberOfChannels(); i_ch++) {
    AnalysisTree::Particle& lambdarec = lambda_reco_->GetChannel(i_ch);

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
          if (!(simtrackid3 < 0)) {
            const AnalysisTree::Particle& simtrack3 = mc_particles_->GetChannel(simtrackid3);
            if (simtrack1.GetField<int>(mother_id_field_id_) == simtrack3.GetField<int>(mother_id_field_id_))
              is_signal = simtrackmother.GetPid() == decay_.GetPdgMother();
          }
        }
      }
    }

    lambdarec.SetField(is_signal, is_signal_field_id_);

    if (is_signal) {
      const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

      auto* lambdasim = lambda_sim_->AddChannel();
      lambdasim->Init(out_config_->GetBranchConfig(lambda_sim_->GetId()));
      lambdasim->SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
      lambdasim->SetMass(simtrackmother.GetMass());
      lambdasim->SetPid(simtrackmother.GetPid());
      lambda_reco2sim_->AddMatch(lambdarec.GetId(), lambdasim->GetId());
    }
  }
}

void ConverterOut::InitIndexes() {

  const auto& out_branch = out_config_->GetBranchConfig(lambda_reco_->GetId());

  auto branch_conf_sim_event = config_->GetBranchConfig(in_branches_[2]);
  b_field_id_ = branch_conf_sim_event.GetFieldId("b");

  x_field_id_ = out_branch.GetFieldId("x");
  y_field_id_ = out_branch.GetFieldId("y");
  z_field_id_ = out_branch.GetFieldId("z");

  daughter1_id_field_id_ = out_branch.GetFieldId("daughter1id");
  daughter2_id_field_id_ = out_branch.GetFieldId("daughter2id");
  daughter3_id_field_id_ = out_branch.GetFieldId("daughter3id");

  px_err_field_id_ = out_branch.GetFieldId("pxerr");
  py_err_field_id_ = out_branch.GetFieldId("pyerr");
  pz_err_field_id_ = out_branch.GetFieldId("pzerr");
  mass_err_field_id_ = out_branch.GetFieldId("masserr");

  if (mc_particles_) {
    auto branch_conf_sim = config_->GetBranchConfig(in_branches_[0]);
    mother_id_field_id_ = branch_conf_sim.GetFieldId("mother_id");
    is_signal_field_id_ = out_branch.GetFieldId("is_signal");
  }

  chi2primpos_field_id_ = out_branch.GetFieldId("chi2primpos");
  chi2primneg_field_id_ = out_branch.GetFieldId("chi2primneg");
  distance_field_id_ = out_branch.GetFieldId("distance");
  cosinepos_field_id_ = out_branch.GetFieldId("cosinepos");
  cosineneg_field_id_ = out_branch.GetFieldId("cosineneg");
  chi2geo_field_id_ = out_branch.GetFieldId("chi2geo");
  l_field_id_ = out_branch.GetFieldId("l");
  ldl_field_id_ = out_branch.GetFieldId("ldl");
  isfrompv_field_id_ = out_branch.GetFieldId("isfrompv");
  cosinetopo_field_id_ = out_branch.GetFieldId("cosinetopo");
  chi2topo_field_id_ = out_branch.GetFieldId("chi2topo");
  nhits_pos_field_id_ = out_branch.GetFieldId("nhitspos");
  nhits_neg_field_id_ = out_branch.GetFieldId("nhitsneg");

  chi2primthird_field_id_ = out_branch.GetFieldId("chi2primthird");
  distancethird_field_id_ = out_branch.GetFieldId("distancethird");
  cosinethird_field_id_ = out_branch.GetFieldId("cosinethird");
  chi2geothree_field_id_ = out_branch.GetFieldId("chi2geothree");
  cosinetopothree_field_id_ = out_branch.GetFieldId("cosinetopothree");
  chi2topothree_field_id_ = out_branch.GetFieldId("chi2topothree");
  nhits_third_field_id_ = out_branch.GetFieldId("nhitsthird");
}
