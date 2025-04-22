#include "ConverterOut.hpp"

#include "PFSimpleTask.hpp"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <algorithm>

void ConverterOut::CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle, const int Ndaughters) const {
  
  //particle.SetIsAllowedSetMassAndChargeExplicitly(true); //Uncomment when most recent AnalysisTree is installed in cbmroot
  
  particle.SetMomentum(kf_particle.GetPx(), kf_particle.GetPy(), kf_particle.GetPz());
  particle.SetMass(kf_particle.GetMass());
  //particle.SetCharge(kf_particle.GetCharge());  //Uncomment when most recent AnalysisTree is installed in cbmroot
  particle.SetPid(kf_particle.GetPdg());
  particle.SetField(kf_particle.GetId(), particle_id_field_id_);
  
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
  
  for (int i = 0; i < Ndaughters; ++i) {
    particle.SetField(kf_particle.GetChi2Prim(i), chi2prim_field_id_ + i);
    particle.SetField(kf_particle.GetCos(i), cosine_field_id_ + i);
    particle.SetField(kf_particle.GetDaughterIds().at(i), daughter_id_field_id_ + i);
    if(Ndaughters == 3) {
      particle.SetField(kf_particle.GetChi2Geo(i + 1), chi2geo_sm_field_id_ + i);
      particle.SetField(kf_particle.GetCosOpen(i + 1), cosopen_sm_field_id_ + i);
      particle.SetField(kf_particle.GetChi2Topo(i + 1), chi2topo_sm_field_id_ + i);
      particle.SetField(kf_particle.GetCosineTopo(i + 1), cosine_topo_sm_field_id_ + i);
    }
  }
  
  if(Ndaughters_max_ == 3 && Ndaughters < 3) {
    particle.SetField(-999.f, chi2prim_field_id_ + 2);
    particle.SetField(-9.f, cosine_field_id_ + 2);
    particle.SetField(-999, daughter_id_field_id_ + 2);

    for (int i = 0; i < Ndaughters_max_; ++i) {
      particle.SetField(-999.f, chi2geo_sm_field_id_ + i);
      particle.SetField(-9.f, cosopen_sm_field_id_ + i);
      particle.SetField(-999.f, chi2topo_sm_field_id_ + i);
      particle.SetField(-9.f, cosine_topo_sm_field_id_ + i);
    }
  }
  
  particle.SetField(kf_particle.GetDistance(), distance_field_id_);

  if(Ndaughters_max_ == 3) {
    if(Ndaughters == 3)
      particle.SetField(kf_particle.GetDistanceToSV(), distance_field_id_ + 1);
    else
      particle.SetField(-999.f, distance_field_id_ + 1);
  }
  
  particle.SetField(kf_particle.GetChi2Geo(0), chi2geo_field_id_);
  particle.SetField(kf_particle.GetCosOpen(0), chi2geo_field_id_ + 1);
  particle.SetField(kf_particle.GetL(), chi2geo_field_id_ + 2);
  particle.SetField(kf_particle.GetLdL(), chi2geo_field_id_ + 3);
  particle.SetField(kf_particle.GetChi2Topo(0), chi2geo_field_id_ + 4);
  particle.SetField(kf_particle.GetCosineTopo(0), chi2geo_field_id_ + 5);
  particle.SetField(kf_particle.GetChi2PrimMother(), chi2prim_mother_field_id_);
  particle.SetField(kf_particle.GetInvMassDiscr(), invmass_discr_field_id_); 
}

void ConverterOut::Exec() {

  candidates_ = pfsimple_task_->GetSimpleFinder()->GetCandidates();

  particle_reco_->ClearChannels();

  if (is_write_mc_ == true) {
    particle_sim_->ClearChannels();
    particle_reco2sim_->Clear();
  }
  
  rec2sim_daughters_.clear();

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  events_->SetField(float(sim_events_->GetField<float>(b_field_id_)), 0);//TODO

  const auto& br_conf = out_config->GetBranchConfig(particle_reco_->GetId());

  for (const auto& candidate : candidates_) {

    const int Ndaughters = candidate.GetDaughterIds().size();

    AnalysisTree::Particle particle(particle_reco_->GetNumberOfChannels(), br_conf);
    CopyParticle(candidate, particle, Ndaughters);
    if (is_write_mc_ == true) {
      MatchWithMc(candidate, particle, Ndaughters);
    }
    
    bool is_write = true;
    if (output_cuts_) {
      is_write = output_cuts_->Apply(particle);
    }

    if (is_write) {
      auto& particlerec = particle_reco_->AddChannel(br_conf);
      particlerec = particle;
    }
  }
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
  out_particles.AddField<int>("particle_id", "");
  out_particles.AddFields<float>({"x", "y", "z", "x_error", "y_error", "z_error"}, "Position and its error, cm");
  out_particles.AddFields<float>({"pT_err", "phi_err", "eta_err", "mass_err"}, "Momentum error");

  std::vector<int> ndaughters;
  for (const auto& decay : decays_)
    ndaughters.push_back(decay.GetNDaughters());
  Ndaughters_max_ = *std::max_element(ndaughters.begin(), ndaughters.end());
  
  if (Ndaughters_max_ == 3) {
    out_particles.AddFields<int>({"daughter1_id", "daughter2_id", "daughter3_id"}, "");
    out_particles.AddFields<float>({"chi2_prim_first", "chi2_prim_second", "chi2_prim_third"}, "");
    out_particles.AddFields<float>({"distance", "distance_sv"}, "Distance between the particles / between third particle and SV, cm");
    out_particles.AddFields<float>({"cosine_first", "cosine_second", "cosine_third"}, "Cos between mother and daughter particle");
    out_particles.AddFields<float>({"chi2_geo_sm1", "chi2_geo_sm2", "chi2_geo_sm3"}, "");
    out_particles.AddFields<float>({"cosopen_sm1", "cosopen_sm2", "cosopen_sm3"}, "");
    out_particles.AddFields<float>({"chi2_topo_sm1", "chi2_topo_sm2", "chi2_topo_sm3"}, "");
    out_particles.AddFields<float>({"cosine_topo_sm1", "cosine_topo_sm2", "cosine_topo_sm3"}, "");
  } else if (Ndaughters_max_ == 2) {
    out_particles.AddFields<int>({"daughter1_id", "daughter2_id"}, "");
    out_particles.AddFields<float>({"chi2_prim_first", "chi2_prim_second"}, "");
    out_particles.AddFields<float>({"distance", "distance_third"}, "Distance between the particles, cm");
    out_particles.AddFields<float>({"cosine_first", "cosine_second"}, "Cos between mother and daughter particle");
  }

  out_particles.AddFields<float>({"chi2_geo", "cosopen", "l", "l_over_dl", "chi2_topo", "cosine_topo"}, "");
  out_particles.AddField<float>("chi2_prim_mother", "chi2 of the mother to the primary vertex (PV)");
  out_particles.AddField<float>("invmass_discr", "Discrepancy in mass of the V0 candidate invariant mass from the PDG value in terms of characteristic sigma");                

  
  AnalysisTree::BranchConfig out_particles_sim(out_branch_sim, AnalysisTree::DetType::kParticle);

  if (!mc_particles_)
     is_write_mc_ = false;

  for (const auto& decay : decays_) 
    if (decay.GetIsDoNotWriteMother() == true) {
      is_write_mc_ = false;
      break;
    }
    
  if (is_write_mc_ == true) {
    out_particles.AddField<int>("generation", "");
    out_particles_sim.AddField<int>("geant_process_id", "");
    out_particles_sim.AddField<int>("mother_id", "particle mother's id in SimParticles branch");   // particle mother's id in SimParticles branch
    out_particles_sim.AddField<int>("mc_id", "");
  }

  man->AddBranch(events_, EventBranch);
  man->AddBranch(particle_reco_, out_particles);
  
  if (is_write_mc_ == true) {
    man->AddBranch(particle_sim_, out_particles_sim);
    man->AddMatching(out_branch, out_branch_sim, particle_reco2sim_);
  }
  
  if (output_cuts_)
    output_cuts_->Init(*out_config);

  events_->Init(EventBranch);
  InitIndexes();
}

int ConverterOut::GetMothersSimId(const OutputContainer& candidate, AnalysisTree::Particle& particlerec, const int Ndaughters) const {
  /** Returnes mother's id in SimParticles branch if:
   ** 1) all reconstructed daughters have matching with simulated particles
   ** 2) all of them belong to the same mother
   ** 3) mother has PDG code as expected according to decay reconstruction hypothesis
   **/
  std::vector<int> daughter_sim_id;
  for (int i = 0; i < Ndaughters; i++) {
    if (candidate.GetDaughterGenerations().at(i) == 0)
      daughter_sim_id.push_back(rec_to_mc_->GetMatch(particlerec.GetField<int>(daughter_id_field_id_ + i)));
    else {
      if (rec2sim_daughters_.find(particlerec.GetField<int>(daughter_id_field_id_+ i)) != rec2sim_daughters_.end()) {
	int daughter_sim_id_pf = rec2sim_daughters_.find(particlerec.GetField<int>(daughter_id_field_id_ + i))->second;
	daughter_sim_id.push_back(particle_sim_->GetChannel(daughter_sim_id_pf).GetField<int>(mc_id_field_id_));
      }
      else
	return -1;
    }
  }

  if (*std::min_element(daughter_sim_id.begin(), daughter_sim_id.end()) < 0)// at least one daughter has no matching with mc
    return -1;
  
  std::vector<int> mother_sim_id;
  for (int i = 0; i < Ndaughters; i++) 
      mother_sim_id.push_back(mc_particles_->GetChannel(daughter_sim_id.at(i)).GetField<int>(mother_id_field_id_));
    
  if (*std::min_element(mother_sim_id.begin(), mother_sim_id.end()) != *std::max_element(mother_sim_id.begin(), mother_sim_id.end()))// daughters belong to not the same mother
	return -1;

  if (mother_sim_id.at(0) < 0)// mother has negative id
    return -1;

  if (mc_particles_->GetChannel(mother_sim_id.at(0)).GetPid() != particlerec.GetPid())// mother has not PDG which was supposed
    return -1;

  return mother_sim_id.at(0);
}

int ConverterOut::DetermineGeneration(int mother_sim_id) const {
  int generation = 0;
  int older_id = mother_sim_id;
  while (older_id >= 0) {
    const auto& simtrackolder = mc_particles_->GetChannel(older_id);
    older_id = simtrackolder.GetField<int>(mother_id_field_id_);
    generation++;
  }

  return generation;
}

void ConverterOut::MatchWithMc(const OutputContainer& candidate, AnalysisTree::Particle& particlerec, const int Ndaughters) {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  int mother_id = GetMothersSimId(candidate, particlerec, Ndaughters);
  int generation = DetermineGeneration(mother_id);
  if(is_detailed_bg_ && generation==0) generation = DetermineBGType(candidate, particlerec, Ndaughters);
  particlerec.SetField(generation, generation_field_id_);
 
  if (generation < 1) return;

  const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

  auto& particlesim = particle_sim_->AddChannel(out_config->GetBranchConfig(particle_sim_->GetId()));
  //particlesim.SetIsAllowedSetMassAndChargeExplicitly(true); //Uncomment when most recent AnalysisTree is installed in cbmroot
  particlesim.SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
  particlesim.SetMass(simtrackmother.GetMass());
  //particlesim.SetCharge(simtrackmother.GetCharge());  //Uncomment when most recent AnalysisTree is installed in cbmroot
  particlesim.SetPid(simtrackmother.GetPid());
  particlesim.SetField(simtrackmother.GetField<int>(mother_id_field_id_), mother_id_field_id_w_);
  particlesim.SetField(simtrackmother.GetField<int>(g4process_field_id_), g4process_field_id_w_);
  int id = simtrackmother.GetId();
  particlesim.SetField(id, mc_id_field_id_);
  particle_reco2sim_->AddMatch(particlerec.GetId(), particlesim.GetId());
  rec2sim_daughters_[particlerec.GetField<int>(particle_id_field_id_)] = particlesim.GetId();
}

void ConverterOut::InitIndexes() {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  const auto& out_branch_reco = out_config->GetBranchConfig(particle_reco_->GetId());
  const auto& out_branch_sim = out_config->GetBranchConfig(particle_sim_->GetId());

  auto branch_conf_sim_event = config_->GetBranchConfig(sim_events_name_);
  b_field_id_ = branch_conf_sim_event.GetFieldId("b");

  x_field_id_ = out_branch_reco.GetFieldId("x");
  daughter_id_field_id_ = out_branch_reco.GetFieldId("daughter1_id");
  particle_id_field_id_ = out_branch_reco.GetFieldId("particle_id");
  pt_err_field_id_ = out_branch_reco.GetFieldId("pT_err");

  if (is_write_mc_ == true) {
    auto branch_conf_sim  = config_->GetBranchConfig(mc_particles_name_);
    auto out_branch_sim   = out_config->GetBranchConfig(particle_sim_->GetId());
    mother_id_field_id_   = branch_conf_sim.GetFieldId("mother_id");
    g4process_field_id_   = branch_conf_sim.GetFieldId("geant_process_id");
    generation_field_id_  = out_branch_reco.GetFieldId("generation");
    g4process_field_id_w_ = out_branch_sim.GetFieldId("geant_process_id");
    mother_id_field_id_w_ = out_branch_sim.GetFieldId("mother_id");
    mc_id_field_id_       = out_branch_sim.GetFieldId("mc_id");
  }

  chi2prim_field_id_ = out_branch_reco.GetFieldId("chi2_prim_first");
  distance_field_id_ = out_branch_reco.GetFieldId("distance");
  cosine_field_id_ = out_branch_reco.GetFieldId("cosine_first");

  chi2geo_sm_field_id_ = out_branch_reco.GetFieldId("chi2_geo_sm1");
  chi2topo_sm_field_id_ = out_branch_reco.GetFieldId("chi2_topo_sm1");
  cosine_topo_sm_field_id_ = out_branch_reco.GetFieldId("cosine_topo_sm1");

  chi2geo_field_id_ = out_branch_reco.GetFieldId("chi2_geo");
  invmass_discr_field_id_ = out_branch_reco.GetFieldId("invmass_discr");
}

std::pair<int, int> ConverterOut::DetermineDaughtersMCStatus(const int daughter_rec_id, const Pdg_t mother_expected_pdg, const int daughter_gen) const {

  int daughter_sim_id;
  if (daughter_gen == 0)
    daughter_sim_id = rec_to_mc_->GetMatch(daughter_rec_id);
  else 
    if (rec2sim_daughters_.find(daughter_rec_id) != rec2sim_daughters_.end()) {
      int daughter_sim_id_pf = rec2sim_daughters_.find(daughter_rec_id)->second;
      daughter_sim_id = particle_sim_->GetChannel(daughter_sim_id_pf).GetField<int>(mc_id_field_id_);
    }
    else daughter_sim_id = -1;

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

int ConverterOut::DetermineBGType(const OutputContainer& candidate, AnalysisTree::Particle& particle, const int Ndaughters) {
  std::vector<std::pair<int, int>> daughters_statuses;
  for (int i = 0; i < Ndaughters; i++) {
    auto daughter_rec_id = particle.GetField<int>(daughter_id_field_id_ + i);
    int daughter_gen = candidate.GetDaughterGenerations().at(i);
    daughters_statuses.emplace_back(DetermineDaughtersMCStatus(daughter_rec_id, particle.GetPid(), daughter_gen));
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

