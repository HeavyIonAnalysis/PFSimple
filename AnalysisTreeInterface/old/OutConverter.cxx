#include "OutConverter.h"

void OutConverter::InitAT()
{
  out_file_ = TFile::Open(out_file_name_.c_str(), "recreate");

//  ***** Lambda Candidates *******    

  AnalysisTree::BranchConfig LambdaRecoBranch("LambdaCandidates", AnalysisTree::DetType::kParticle);
  
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
  LambdaRecoBranch.AddField<float>("invmass");
  LambdaRecoBranch.AddField<float>("rapidity_lab");
  LambdaRecoBranch.AddField<float>("rapidity_cm");
  LambdaRecoBranch.AddField<int>("pdg");
  LambdaRecoBranch.AddField<int>("daughter1id");
  LambdaRecoBranch.AddField<int>("daughter2id");
  
  out_config_.AddBranchConfig( LambdaRecoBranch );
  lambda_reco_ = new AnalysisTree::Particles(out_config_.GetLastId());
  
  out_tree_ = new TTree("aTree", "AnalysisTree Lambda");
  out_tree_ -> Branch("LambdaCandidates", "AnalysisTree::Particles", &lambda_reco_);
  out_config_.Write("Configuration");
  
  x_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("x");
  y_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("y");
  z_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("z");
  mass_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("invmass");
  rap_lab_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("rapidity_lab");
  rap_cm_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("rapidity_cm");
  pdg_field_id_w_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("pdg");
  daughter1_id_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("daughter1id");
  daughter2_id_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("daughter2id");
  
  chi2primpos_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("chi2primpos");
  chi2primneg_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("chi2primneg");
  distance_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("distance");
  cosinepos_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("cosinepos");
  cosineneg_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("cosineneg");
  chi2geo_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("chi2geo");
  l_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("l");
  ldl_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("ldl");
  isfrompv_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("isfrompv");
  cosinetopo_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("cosinetopo");
  chi2topo_field_id_ = out_config_.GetBranchConfig( lambda_reco_->GetId() ).GetFieldId("chi2topo");
}

void OutConverter::WriteCandidates(const std::vector<OutputContainer>& canditates)
{
  lambda_reco_->ClearChannels();
  
  for(const auto& candidate : canditates)
  {
    auto* lambdarec = lambda_reco_->AddChannel();
    lambdarec->Init(out_config_.GetBranchConfig(lambda_reco_->GetId()));
    
    const KFParticle& particle = candidate.GetParticle();
    float mass, masserr;
    particle.GetMass(mass, masserr);
    lambdarec->SetField(mass, mass_field_id_);
    lambdarec->SetField(particle.DaughterIds()[0], daughter1_id_field_id_);
    lambdarec->SetField(particle.DaughterIds()[1], daughter2_id_field_id_);
//     std::cout << "FC\t" << particle.DaughterIds()[0] << "\t" << particle.DaughterIds()[1] << "\n";
    
    lambdarec->SetField(particle.X(), x_field_id_);
    lambdarec->SetField(particle.Y(), y_field_id_);
    lambdarec->SetField(particle.Z(), z_field_id_);
    lambdarec->SetMomentum(particle.GetPx(), particle.GetPy(), particle.GetPz());
    lambdarec->SetField( particle.GetRapidity(), rap_lab_field_id_);
    lambdarec->SetField( particle.GetRapidity() - RapidityShift(12.), rap_cm_field_id_);
    lambdarec->SetField( particle.GetPDG(), pdg_field_id_w_ );
    
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
  out_tree_->Fill();
}

void OutConverter::Finish()
{
  out_tree_->Write();  
  out_file_->Close();
}

float OutConverter::RapidityShift(float pbeam)
{
  const float nucleonmass = 0.939;
  return TMath::Log(((sqrt(nucleonmass*nucleonmass + pbeam*pbeam) + pbeam) / (sqrt(nucleonmass*nucleonmass + pbeam*pbeam) - pbeam)))/4.;    
}