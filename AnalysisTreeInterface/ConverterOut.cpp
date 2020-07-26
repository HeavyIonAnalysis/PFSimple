#include "ConverterOut.h"

#include "TTree.h"

#include "AnalysisTree/DataHeader.hpp"

void ConverterOut::Exec()
{
  lambda_reco_->ClearChannels();

  for(const auto& candidate : canditates_)
  {
    auto* lambdarec = lambda_reco_->AddChannel();
    lambdarec->Init(out_config_->GetBranchConfig(lambda_reco_->GetId()));

    const KFParticle& particle = candidate.GetParticle();
    float mass, masserr;
    particle.GetMass(mass, masserr);
//    lambdarec->SetField(mass, mass_field_id_);
    lambdarec->SetField(particle.DaughterIds()[0], daughter1_id_field_id_);
    lambdarec->SetField(particle.DaughterIds()[1], daughter2_id_field_id_);
//     std::cout << "FC\t" << particle.DaughterIds()[0] << "\t" << particle.DaughterIds()[1] << "\n";

    lambdarec->SetField(particle.X(), x_field_id_);
    lambdarec->SetField(particle.Y(), y_field_id_);
    lambdarec->SetField(particle.Z(), z_field_id_);
    lambdarec->SetMomentum(particle.GetPx(), particle.GetPy(), particle.GetPz());
//    lambdarec->SetField( particle.GetRapidity(), rap_lab_field_id_);
    lambdarec->SetField( particle.GetRapidity() - data_header_->GetBeamRapidity(), rap_cm_field_id_);
//    lambdarec->SetField( particle.GetPDG(), pdg_field_id_w_ );

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

  LambdaRecoBranch.AddField<float>("x");
  LambdaRecoBranch.AddField<float>("y");
  LambdaRecoBranch.AddField<float>("z");
  LambdaRecoBranch.AddField<float>("rapidity_cm");
  LambdaRecoBranch.AddField<int>("daughter1id");
  LambdaRecoBranch.AddField<int>("daughter2id");
  if(mc_particles_) {
    LambdaRecoBranch.AddField<bool>("is_signal");
  }

  out_config_->AddBranchConfig( LambdaRecoBranch );
  lambda_reco_ = new AnalysisTree::Particles(out_config_->GetLastId());
  out_tree_ -> Branch(out_branch_.c_str(), "AnalysisTree::Particles", &lambda_reco_);

  InitIndexes();
}

void ConverterOut::MatchWithMc(){

//  std::cout << mc_particles_->GetNumberOfChannels() << " " << rec_tracks_->GetNumberOfChannels() << std::endl;
//  std::cout << rec_to_mc_->GetMatch(0) << std::endl;

}


void ConverterOut::InitIndexes(){

  const auto& out_branch = out_config_->GetBranchConfig( lambda_reco_->GetId());

  x_field_id_ = out_branch.GetFieldId("x");
  y_field_id_ = out_branch.GetFieldId("y");
  z_field_id_ = out_branch.GetFieldId("z");
//  mass_field_id_ = out_branch.GetFieldId("invmass");
//  rap_lab_field_id_ = out_branch.GetFieldId("rapidity_lab");
  rap_cm_field_id_ = out_branch.GetFieldId("rapidity_cm");
//  pdg_field_id_w_ = out_branch.GetFieldId("pdg");
  daughter1_id_field_id_ = out_branch.GetFieldId("daughter1id");
  daughter2_id_field_id_ = out_branch.GetFieldId("daughter2id");

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
}