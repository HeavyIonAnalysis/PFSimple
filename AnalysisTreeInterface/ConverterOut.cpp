#include "ConverterOut.h"

#include "TTree.h"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"

void ConverterOut::Exec()
{
  lambda_reco_->ClearChannels();
  lambda_sim_->ClearChannels();
  lambda_reco2sim_->Clear();

  for(const auto& candidate : canditates_)
  {
    auto* lambdarec = lambda_reco_->AddChannel();
    lambdarec->Init(out_config_->GetBranchConfig(lambda_reco_->GetId()));

    const KFParticle& particle = candidate.GetParticle();
    float mass, masserr;
    particle.GetMass(mass, masserr);
    lambdarec->SetField(particle.DaughterIds()[0], daughter1_id_field_id_);
    lambdarec->SetField(particle.DaughterIds()[1], daughter2_id_field_id_);

    lambdarec->SetField(particle.X(), x_field_id_);
    lambdarec->SetField(particle.Y(), y_field_id_);
    lambdarec->SetField(particle.Z(), z_field_id_);
    lambdarec->SetMomentum(particle.GetPx(), particle.GetPy(), particle.GetPz());

    lambdarec->SetPid( particle.GetPDG());
    lambdarec->SetMass(mass);

    lambdarec->SetField( candidate.GetChi2PrimPos(), chi2primpos_field_id_);
    lambdarec->SetField( candidate.GetChi2PrimNeg(), chi2primneg_field_id_);
    lambdarec->SetField( candidate.GetDistance(), distance_field_id_);
    lambdarec->SetField( candidate.GetCosineDaughterPos(), cosinepos_field_id_);
    lambdarec->SetField( candidate.GetCosineDaughterNeg(), cosineneg_field_id_);
    lambdarec->SetField( candidate.GetChi2Geo(), chi2geo_field_id_);
    lambdarec->SetField( candidate.GetL(), l_field_id_);
    lambdarec->SetField( candidate.GetLdL(), ldl_field_id_);
    lambdarec->SetField( candidate.GetIsFromPV(), isfrompv_field_id_);
    lambdarec->SetField( candidate.GetCosineTopo(), cosinetopo_field_id_);
    lambdarec->SetField( candidate.GetChi2Topo(), chi2topo_field_id_);
    lambdarec->SetField( candidate.GetNHitsPos(), nhits_pos_field_id_);
    lambdarec->SetField( candidate.GetNHitsNeg(), nhits_neg_field_id_);
  }
  MatchWithMc();
  out_tree_->Fill();
}


void ConverterOut::Init(std::map<std::string, void*>& branches)
{
  assert(out_config_ && out_tree_);
  out_branch_ = "LambdaCandidates";

  if(!in_branches_.empty()){
    mc_particles_ = (AnalysisTree::Particles*) branches.find(in_branches_[0])->second;
    rec_tracks_ = (AnalysisTree::TrackDetector*) branches.find(in_branches_[1])->second;
    const auto& match = config_->GetMatchName(in_branches_[1], in_branches_[0]);
    rec_to_mc_ = (AnalysisTree::Matching*) branches.find(match)->second;
  }

  AnalysisTree::BranchConfig LambdaRecoBranch(out_branch_, AnalysisTree::DetType::kParticle);

  LambdaRecoBranch.AddField<float>("chi2primpos");
  LambdaRecoBranch.AddField<float>("chi2primneg");
  LambdaRecoBranch.AddField<float>("distance");
  LambdaRecoBranch.AddField<float>("cosinepos");
  LambdaRecoBranch.AddField<float>("cosineneg");
  LambdaRecoBranch.AddField<float>("chi2geo");
  LambdaRecoBranch.AddField<float>("l");
  LambdaRecoBranch.AddField<float>("ldl");
  LambdaRecoBranch.AddField<int>  ("isfrompv");
  LambdaRecoBranch.AddField<float>("cosinetopo");
  LambdaRecoBranch.AddField<float>("chi2topo");
  LambdaRecoBranch.AddField<int>  ("nhitspos");
  LambdaRecoBranch.AddField<int>  ("nhitsneg");

  LambdaRecoBranch.AddField<float>("x");
  LambdaRecoBranch.AddField<float>("y");
  LambdaRecoBranch.AddField<float>("z");
  LambdaRecoBranch.AddField<int>("daughter1id");
  LambdaRecoBranch.AddField<int>("daughter2id");
  if(mc_particles_) {
    LambdaRecoBranch.AddField<bool>("is_signal");
  }

  out_config_->AddBranchConfig( LambdaRecoBranch );
  lambda_reco_ = new AnalysisTree::Particles(out_config_->GetLastId());
  
  AnalysisTree::BranchConfig LambdaSimBranch("LambdaSimulated", AnalysisTree::DetType::kParticle);
  out_config_->AddBranchConfig(LambdaSimBranch);
  lambda_sim_ = new AnalysisTree::Particles(out_config_->GetLastId());
  
  lambda_reco2sim_ = new AnalysisTree::Matching(out_config_->GetBranchConfig(out_branch_).GetId(), out_config_->GetBranchConfig("LambdaSimulated").GetId());

  out_config_->AddMatch(lambda_reco2sim_);
  
  out_tree_ -> Branch(out_branch_.c_str(), "AnalysisTree::Particles", &lambda_reco_);
  out_tree_ -> Branch("LambdaSimulated", "AnalysisTree::Particles", &lambda_sim_);
  out_tree_ -> Branch("LambdaCandidates2LambdaSimulated", "AnalysisTree::Matching", &lambda_reco2sim_);

  InitIndexes();
}

void ConverterOut::MatchWithMc(){

//  for(auto& lambdarec : *lambda_reco_->Channels()) {
  for(int i_ch=0; i_ch<lambda_reco_->GetNumberOfChannels(); i_ch++){
    AnalysisTree::Particle& lambdarec = lambda_reco_->GetChannel(i_ch);

    const int simtrackid1 = rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter1_id_field_id_));
    const int simtrackid2 = rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter2_id_field_id_));
    bool is_signal = false;
    int mother_id = -999;
    if(!(simtrackid1 < 0 || simtrackid2 < 0))
    {
      const AnalysisTree::Particle& simtrack1 = mc_particles_->GetChannel(simtrackid1);
      const AnalysisTree::Particle& simtrack2 = mc_particles_->GetChannel(simtrackid2);
      if(simtrack1.GetField<int>(mother_id_field_id_) == simtrack2.GetField<int>(mother_id_field_id_))
      {
        mother_id = simtrack1.GetField<int>(mother_id_field_id_);
        if(mother_id < 0) continue;
        const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

        is_signal = simtrackmother.GetPid() == pdg_lambda;
      }
    }
    lambdarec.SetField( is_signal, is_signal_field_id_);

    if(is_signal) {
      const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);
      
      auto* lambdasim = lambda_sim_->AddChannel();
      lambdasim->Init(out_config_->GetBranchConfig(lambda_sim_->GetId()));
      lambdasim->SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
      lambdasim->SetPid(simtrackmother.GetPid());
      lambdasim->SetMass(simtrackmother.GetMass());      
      lambda_reco2sim_->AddMatch(lambdarec.GetId(), lambdasim->GetId());
    }

  }
}


void ConverterOut::InitIndexes(){

  const auto& out_branch = out_config_->GetBranchConfig( lambda_reco_->GetId());

  x_field_id_ = out_branch.GetFieldId("x");
  y_field_id_ = out_branch.GetFieldId("y");
  z_field_id_ = out_branch.GetFieldId("z");

  daughter1_id_field_id_ = out_branch.GetFieldId("daughter1id");
  daughter2_id_field_id_ = out_branch.GetFieldId("daughter2id");
  
  if(mc_particles_){
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
}