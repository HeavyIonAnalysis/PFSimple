#include "ConverterOut.hpp"

#include "PFSimpleTask.hpp"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <algorithm>

void ConverterOut::CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) const {

  particle.SetMomentum(kf_particle.GetPx(), kf_particle.GetPy(), kf_particle.GetPz());
  particle.SetMass(kf_particle.GetMass());
  particle.SetPid(kf_particle.GetPdg());

  particle.SetField(kf_particle.GetX(), x_field_id_);
  particle.SetField(kf_particle.GetY(), x_field_id_ + 1);
  particle.SetField(kf_particle.GetZ(), x_field_id_ + 2);
  particle.SetField(kf_particle.GetXError(), x_field_id_ + 3);
  particle.SetField(kf_particle.GetYError(), x_field_id_ + 4);
  particle.SetField(kf_particle.GetZError(), x_field_id_ + 5);

  particle.SetField(kf_particle.GetPtError(), pt_err_field_id_);
  particle.SetField(kf_particle.GetPhiError(), pt_err_field_id_ + 1);
  particle.SetField(kf_particle.GetEtaError(), pt_err_field_id_ + 2);
  particle.SetField(kf_particle.GetMassError(), pt_err_field_id_ + 3);

  for (int i = 0; i < decay_.GetNDaughters(); ++i) {
    particle.SetField(kf_particle.GetChi2Prim(i), chi2prim_field_id_ + i);
    particle.SetField(kf_particle.GetCos(i), cosine_field_id_ + i);
    particle.SetField(kf_particle.GetDaughterIds().at(i), daughter_id_field_id_ + i);
  }

  particle.SetField(kf_particle.GetDistance(), distance_field_id_);

  if (decay_.GetNDaughters() == 3) {
    particle.SetField(kf_particle.GetDistanceToSV(), distance_field_id_ + 1);
    for (int i = 0; i < decay_.GetNDaughters(); ++i) {
      particle.SetField(kf_particle.GetChi2Geo(i + 1), chi2geo_sm_field_id_ + i);
      particle.SetField(kf_particle.GetChi2Topo(i + 1), chi2topo_sm_field_id_ + i);
      particle.SetField(kf_particle.GetCosineTopo(i + 1), cosine_topo_sm_field_id_ + i);
    }
  }

  particle.SetField(kf_particle.GetChi2Geo(0), chi2geo_field_id_);
  particle.SetField(kf_particle.GetL(), chi2geo_field_id_ + 1);
  particle.SetField(kf_particle.GetLdL(), chi2geo_field_id_ + 2);
  particle.SetField(kf_particle.GetChi2Topo(0), chi2geo_field_id_ + 3);
  particle.SetField(kf_particle.GetCosineTopo(0), chi2geo_field_id_ + 4);
}

void ConverterOut::Exec() {

  candidates_ = pfsimple_task_->GetSimpleFinder()->GetCandidates();

  lambda_reco_->ClearChannels();
  lambda_sim_->ClearChannels();
  lambda_reco2sim_->Clear();

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  events_->SetField(float(sim_events_->GetField<float>(b_field_id_)), 0);//TODO

  const auto& br_conf = out_config->GetBranchConfig(lambda_reco_->GetId());

  for (const auto& candidate : candidates_) {

    AnalysisTree::Particle particle(lambda_reco_->GetNumberOfChannels(), br_conf);
    CopyParticle(candidate, particle);
    if (mc_particles_) {
      MatchWithMc(particle);
    }

    bool is_write = true;
    if (output_cuts_) {
      is_write = output_cuts_->Apply(particle);
    }

    if (is_write) {
      auto& lambdarec = lambda_reco_->AddChannel(br_conf);
      lambdarec = particle;
    }
  }
  delete pfsimple_task_->GetSimpleFinder();
}

void ConverterOut::Init() {

  this->SetInputBranchNames({sim_events_name_, rec_tracks_name_, mc_particles_name_});

  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();

  sim_events_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(sim_events_name_));
  mc_particles_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch(mc_particles_name_));
  rec_to_mc_ = chain->GetMatchPointers().find(config_->GetMatchName(rec_tracks_name_, mc_particles_name_))->second;

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  std::string out_branch_event = "Events";
  std::string out_branch = std::string("Candidates");
  std::string out_branch_sim = std::string("Simulated");
  std::string out_branch_reco2sim = out_branch + "2" + out_branch_sim;

  AnalysisTree::BranchConfig EventBranch(out_branch_event, AnalysisTree::DetType::kEventHeader);
  EventBranch.AddField<float>("b", "Impact parameter, fm");

  AnalysisTree::BranchConfig out_particles(out_branch, AnalysisTree::DetType::kParticle);
  out_particles.AddFields<float>({"x", "y", "z", "x_error", "y_error", "z_error"}, "Position and its error, cm");
  out_particles.AddFields<float>({"pT_err", "phi_err", "eta_err", "mass_err"}, "Momentum error");

  if (decay_.GetNDaughters() == 3) {
    out_particles.AddFields<int>({"daughter1_id", "daughter2_id", "daughter3_id"}, "");
    out_particles.AddFields<float>({"chi2_prim_first", "chi2_prim_second", "chi2_prim_third"}, "");
    out_particles.AddFields<float>({"distance", "distance_sv"}, "Distance between the particles, cm");
    out_particles.AddFields<float>({"cosine_first", "cosine_second", "cosine_third"}, "Cos between mother and daughter particle");
    out_particles.AddFields<float>({"chi2_geo_sm1", "chi2_geo_sm2", "chi2_geo_sm3"}, "");
    out_particles.AddFields<float>({"chi2_topo_sm1", "chi2_topo_sm2", "chi2_topo_sm3"}, "");
    out_particles.AddFields<float>({"cosine_topo_sm1", "cosine_topo_sm2", "cosine_topo_sm3"}, "");
  } else if (decay_.GetNDaughters() == 2) {
    out_particles.AddFields<int>({"daughter1_id", "daughter2_id"}, "");
    out_particles.AddFields<float>({"chi2_prim_first", "chi2_prim_second"}, "");
    out_particles.AddField<float>("distance", "Distance between the particles, cm");
    out_particles.AddFields<float>({"cosine_first", "cosine_second"}, "Cos between mother and daughter particle");
  }

  out_particles.AddFields<float>({"chi2_geo", "l", "l_over_dl", "chi2_topo", "cosine_topo"}, "");

  AnalysisTree::BranchConfig LambdaSimBranch(out_branch_sim, AnalysisTree::DetType::kParticle);

  if (mc_particles_) {
    out_particles.AddField<int>("generation", "");
    LambdaSimBranch.AddField<int>("geant_process_id", "");
  }

  man->AddBranch(events_, EventBranch);
  man->AddBranch(lambda_reco_, out_particles);
  man->AddBranch(lambda_sim_, LambdaSimBranch);
  man->AddMatching(out_branch, out_branch_sim, lambda_reco2sim_);

  if (output_cuts_)
    output_cuts_->Init(*out_config);

  events_->Init(EventBranch);
  InitIndexes();
}

int ConverterOut::GetMothersSimId(AnalysisTree::Particle& lambdarec) {
  std::vector<int> daughter_sim_id;
  for (int i = 0; i < decay_.GetNDaughters(); i++)
    daughter_sim_id.push_back(rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_ + i)));

  if (*std::min_element(daughter_sim_id.begin(), daughter_sim_id.end()) < 0)// at least one daughter has no matching with mc
    return -1;

  std::vector<int> mother_sim_id;
  for (int i = 0; i < decay_.GetNDaughters(); i++)
    mother_sim_id.push_back(mc_particles_->GetChannel(daughter_sim_id.at(i)).GetField<int>(mother_id_field_id_));

  if (*std::min_element(mother_sim_id.begin(), mother_sim_id.end()) != *std::max_element(mother_sim_id.begin(), mother_sim_id.end()))// daughters belong to not the same mother
    return -1;

  if (mother_sim_id.at(0) < 0)// mother has negative id
    return -1;

  if (mc_particles_->GetChannel(mother_sim_id.at(0)).GetPid() != lambdarec.GetPid())// mother has not PDG which was supposed
    return -1;

  return mother_sim_id.at(0);
}

int ConverterOut::DetermineGeneration(int mother_sim_id) {
  int generation = 0;
  int older_id = mother_sim_id;
  while (older_id >= 0) {
    const auto& simtrackolder = mc_particles_->GetChannel(older_id);
    older_id = simtrackolder.GetField<int>(mother_id_field_id_);
    generation++;
  }

  return generation;
}

void ConverterOut::MatchWithMc(AnalysisTree::Particle& lambdarec) {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  int mother_id = GetMothersSimId(lambdarec);
  int generation = DetermineGeneration(mother_id);
  if(is_detailed_bg_ && generation==0) generation = DetermineBGType(lambdarec);
  lambdarec.SetField(generation, generation_field_id_);

  if (generation < 1) return;

  const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

  auto& lambdasim = lambda_sim_->AddChannel(out_config->GetBranchConfig(lambda_sim_->GetId()));

  lambdasim.SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
  lambdasim.SetMass(simtrackmother.GetMass());
  lambdasim.SetPid(simtrackmother.GetPid());
  lambdasim.SetField(simtrackmother.GetField<int>(g4process_field_id_), g4process_field_id_w_);
  lambda_reco2sim_->AddMatch(lambdarec.GetId(), lambdasim.GetId());
}

void ConverterOut::InitIndexes() {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  const auto& out_branch_reco = out_config->GetBranchConfig(lambda_reco_->GetId());
  const auto& out_branch_sim = out_config->GetBranchConfig(lambda_sim_->GetId());

  auto branch_conf_sim_event = config_->GetBranchConfig(sim_events_name_);
  b_field_id_ = branch_conf_sim_event.GetFieldId("b");

  x_field_id_ = out_branch_reco.GetFieldId("x");
  daughter_id_field_id_ = out_branch_reco.GetFieldId("daughter1_id");
  pt_err_field_id_ = out_branch_reco.GetFieldId("pT_err");

  if (mc_particles_) {
    auto branch_conf_sim = config_->GetBranchConfig(mc_particles_name_);
    mother_id_field_id_ = branch_conf_sim.GetFieldId("mother_id");
    g4process_field_id_ = branch_conf_sim.GetFieldId("geant_process_id");
    generation_field_id_ = out_branch_reco.GetFieldId("generation");
    g4process_field_id_w_ = out_branch_sim.GetFieldId("geant_process_id");
  }

  chi2prim_field_id_ = out_branch_reco.GetFieldId("chi2_prim_first");
  distance_field_id_ = out_branch_reco.GetFieldId("distance");
  cosine_field_id_ = out_branch_reco.GetFieldId("cosine_first");

  chi2geo_sm_field_id_ = out_branch_reco.GetFieldId("chi2_geo_sm1");
  chi2topo_sm_field_id_ = out_branch_reco.GetFieldId("chi2_topo_sm1");
  cosine_topo_sm_field_id_ = out_branch_reco.GetFieldId("cosine_topo_sm1");

  chi2geo_field_id_ = out_branch_reco.GetFieldId("chi2_geo");
}

std::pair<int, int> ConverterOut::DetermineDaughtersMCStatus(const int daughter_rec_id, const Pdg_t mother_expected_pdg) const {
  const int daughter_sim_id = rec_to_mc_->GetMatch(daughter_rec_id);
  if(daughter_sim_id < 0) return std::make_pair(1, -999); // no match to MC

  auto& daughter_sim = mc_particles_->GetChannel(daughter_sim_id);
  const int mother_sim_id = daughter_sim.GetField<int>(mother_id_field_id_);
  if(mother_sim_id < 0) return std::make_pair(2, -999); // daughter is primary

  const int geant_process = daughter_sim.GetField<int>(g4process_field_id_);
  auto& mother_sim = mc_particles_->GetChannel(mother_sim_id);
  auto mother_pdg = mother_sim.GetPid();

  int daughter_status{-999};

  if(geant_process != 4 && mother_pdg != mother_expected_pdg) daughter_status = 3; // daughter not from decay, mother's pdg unexpected
  if(geant_process != 4 && mother_pdg == mother_expected_pdg) daughter_status = 4; // daughter not from decay, mother's pdg expected
  if(geant_process == 4 && mother_pdg != mother_expected_pdg) daughter_status = 5; // daughter from decay, mother's pdg unexpected
  if(geant_process == 4 && mother_pdg == mother_expected_pdg) daughter_status = 6; // daughter from decay, mother's pdg expected

  return std::make_pair(daughter_status, mother_sim_id);
}

int ConverterOut::DetermineMotherMCStatus(const int mid1, const int mid2) {
  if(mid1 == -999 || mid2 == -999) return 0;
  if(mid1 == mid2) return 1;
  else             return 2;
}

int ConverterOut::DetermineBGType(AnalysisTree::Particle& particle) {
  std::vector<std::pair<int, int>> daughters_statuses;
  for (int i = 0; i < decay_.GetNDaughters(); i++) {
    auto daughter_rec_id = particle.GetField<int>(daughter_id_field_id_ + i);
    daughters_statuses.emplace_back(DetermineDaughtersMCStatus(daughter_rec_id, particle.GetPid()));
  }

  int result{0};
  int decimal{1};
  for(auto& ds : daughters_statuses) {
    result += decimal * ds.first;
    decimal *= 10;
  }

  std::vector<int> common_mother_statuses;
  common_mother_statuses.emplace_back(DetermineMotherMCStatus(daughters_statuses.at(0).second, daughters_statuses.at(1).second));

  if(daughters_statuses.size() == 3) {
    common_mother_statuses.emplace_back(DetermineMotherMCStatus(daughters_statuses.at(1).second, daughters_statuses.at(2).second));
    common_mother_statuses.emplace_back(DetermineMotherMCStatus(daughters_statuses.at(0).second, daughters_statuses.at(2).second));
  }

  decimal = 1000;
  for(auto& cms : common_mother_statuses) {
    result += decimal * cms;
    decimal *= 10;
  }

  return -result;
}
