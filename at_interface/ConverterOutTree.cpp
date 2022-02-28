#include "ConverterOutTree.hpp"

#include "PFSimpleTask.hpp"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include "TString.h"
#include <algorithm>

void ConverterOutTree::CopyParticle(const OutputContainer& kf_particle) {

  pid_ = kf_particle.GetPdg();
  mass_ = kf_particle.GetMass();
  mass_err_ = kf_particle.GetMassError();

  px_ = kf_particle.GetPx();
  py_ = kf_particle.GetPy();
  pz_ = kf_particle.GetPz();
  x_ = kf_particle.GetX();
  y_ = kf_particle.GetY();
  z_ = kf_particle.GetZ();
  x_err_ = kf_particle.GetXError();
  y_err_ = kf_particle.GetYError();
  z_err_ = kf_particle.GetZError();
  pt_err_ = kf_particle.GetPtError();
  phi_err_ = kf_particle.GetPhiError();
  eta_err_ = kf_particle.GetEtaError();

  daughter_id_1_ = kf_particle.GetDaughterIds().at(0);
  chi2prim_1_ = kf_particle.GetChi2Prim(0);
  cos_1_ = kf_particle.GetCos(0);
  daughter_id_2_ = kf_particle.GetDaughterIds().at(1);
  chi2prim_2_ = kf_particle.GetChi2Prim(1);
  cos_2_ = kf_particle.GetCos(1);

  if (decay_.GetNDaughters() == 3) {
    daughter_id_3_ = kf_particle.GetDaughterIds().at(2);
    chi2prim_3_ = kf_particle.GetChi2Prim(2);
    cos_3_ = kf_particle.GetCos(2);
  }
  distance_ = kf_particle.GetDistance();
  distance_sv_ = kf_particle.GetDistanceToSV();

  if (decay_.GetNDaughters() == 3) {
    chi2geo_sm_1_ = kf_particle.GetChi2Geo(1);
    chi2topo_sm_1_ = kf_particle.GetChi2Topo(1);
    costopo_sm_1_ = kf_particle.GetCosineTopo(1);
    chi2geo_sm_2_ = kf_particle.GetChi2Geo(2);
    chi2topo_sm_2_ = kf_particle.GetChi2Topo(2);
    costopo_sm_2_ = kf_particle.GetCosineTopo(2);
    chi2geo_sm_3_ = kf_particle.GetChi2Geo(3);
    chi2topo_sm_3_ = kf_particle.GetChi2Topo(3);
    costopo_sm_3_ = kf_particle.GetCosineTopo(3);
  }
  chi2geo_ = kf_particle.GetChi2Geo(0);
  chi2topo_ = kf_particle.GetChi2Topo(0);
  costopo_ = kf_particle.GetCosineTopo(0);
  L_ = kf_particle.GetL();
  LdL_ = kf_particle.GetLdL();
  distance_pv_line_ = kf_particle.GetDistanceToPVLine();
}

void ConverterOutTree::Exec() {

  b_ = sim_events_->GetField<float>(config_->GetBranchConfig(sim_events_name_).GetFieldId("b"));
  out_events_->Fill();

  candidates_ = pfsimple_task_->GetSimpleFinder()->GetCandidates();

  for (const auto& candidate : candidates_) {
    id_rec_++;
    CopyParticle(candidate);
    if (mc_particles_)
      MatchWithMc();
    out_reco_->Fill();
  }
}

void ConverterOutTree::Init() {
  
  this->SetInputBranchNames({sim_events_name_, rec_tracks_name_, mc_particles_name_});

  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();

  sim_events_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(sim_events_name_));
  mc_particles_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch(mc_particles_name_));
  rec_to_mc_ = chain->GetMatchPointers().find(config_->GetMatchName(rec_tracks_name_, mc_particles_name_))->second;

  out_file_ = TFile::Open(out_file_name_.c_str(), "RECREATE");

  TDirectory* dir_events = out_file_->mkdir("Events");
  dir_events->cd();
  out_events_ = new TTree("Events", "Events");
  out_events_->Branch("b", &b_, "b_/F");
  TDirectory* dir = out_file_->mkdir(decay_.GetName().c_str());
  dir->cd();
  out_reco_ = new TTree("Signal", "Signal");
  out_reco_->Branch("id", &id_rec_, "id_rec_/I");
  out_reco_->Branch("pid", &pid_, "pid_/I");
  out_reco_->Branch("mass", &mass_, "mass_/F");
  out_reco_->Branch("mass_err", &mass_err_, "mass_err_/F");
  out_reco_->Branch("x", &x_, "x_/F");
  out_reco_->Branch("y", &y_, "y_/F");
  out_reco_->Branch("z", &z_, "z_/F");
  out_reco_->Branch("x_err", &x_err_, "x_err_/F");
  out_reco_->Branch("y_err", &y_err_, "y_err_/F");
  out_reco_->Branch("z_err", &z_err_, "z_err_/F");
  out_reco_->Branch("px", &px_, "px_/F");
  out_reco_->Branch("py", &py_, "py_/F");
  out_reco_->Branch("pz", &pz_, "pz_/F");
  out_reco_->Branch("pt_err", &pt_err_, "pt_err_/F");
  out_reco_->Branch("phi_err", &phi_err_, "phi_err_/F");
  out_reco_->Branch("eta_err", &eta_err_, "eta_err_/F");

  out_reco_->Branch("daughter_id_1", &daughter_id_1_, "daughter_id_1_/I");
  out_reco_->Branch("daughter_id_2", &daughter_id_2_, "daughter_id_2_/I");
  if (decay_.GetNDaughters() == 3) out_reco_->Branch("daughter_id_3", &daughter_id_3_, "daughter_id_3_/I");
  out_reco_->Branch("chi2prim_1", &chi2prim_1_, "chi2prim_1_/F");
  out_reco_->Branch("chi2prim_2", &chi2prim_2_, "chi2prim_2_/F");
  if (decay_.GetNDaughters() == 3) out_reco_->Branch("chi2prim_3", &chi2prim_3_, "chi2prim_3_/F");
  out_reco_->Branch("cos_1", &cos_1_, "cos_1_/F");
  out_reco_->Branch("cos_2", &cos_2_, "cos_2_/F");
  if (decay_.GetNDaughters() == 3) out_reco_->Branch("cos_2", &cos_2_, "cos_2_/F");

  out_reco_->Branch("distance", &distance_, "distance_/F");
  out_reco_->Branch("distance_sv", &distance_sv_, "distance_sv_/F");
  if (decay_.GetNDaughters() == 3) {
    out_reco_->Branch("chi2geo_sm_1", &chi2geo_sm_1_, "chi2geo_sm_1_/F");
    out_reco_->Branch("chi2geo_sm_2", &chi2geo_sm_2_, "chi2geo_sm_2_/F");
    out_reco_->Branch("chi2geo_sm_3", &chi2geo_sm_3_, "chi2geo_sm_3_/F");
    out_reco_->Branch("chi2topo_sm_1", &chi2topo_sm_1_, "chi2topo_sm_1_/F");
    out_reco_->Branch("chi2topo_sm_2", &chi2topo_sm_2_, "chi2topo_sm_2_/F");
    out_reco_->Branch("chi2topo_sm_3", &chi2topo_sm_3_, "chi2topo_sm_3_/F");
    out_reco_->Branch("costopo_sm_1", &costopo_sm_1_, "costopo_sm_1_/F");
    out_reco_->Branch("costopo_sm_2", &costopo_sm_2_, "costopo_sm_2_/F");
    out_reco_->Branch("costopo_sm_3", &costopo_sm_3_, "costopo_sm_3_/F");
  }
  out_reco_->Branch("chi2geo", &chi2geo_, "chi2geo_/F");
  out_reco_->Branch("chi2topo", &chi2topo_, "chi2topo_/F");
  out_reco_->Branch("costopo", &costopo_, "costopo_/F");
  out_reco_->Branch("L", &L_, "L_/F");
  out_reco_->Branch("LdL", &LdL_, "LdL_/F");
  out_reco_->Branch("distance_pv_line", &distance_pv_line_, "distance_pv_line_/F");

  if (mc_particles_) {
    out_reco_->Branch("generation", &generation_, "generation_/I");
    out_reco_->Branch("is_mc", &is_mc_, "is_mc_/I");

    out_sim_ = new TTree("MC", "MC");
    out_sim_->Branch("id", &id_mc_, "id_mc_/I");
    out_sim_->Branch("pid", &pid_, "pid_/I");
    out_sim_->Branch("mass_mc", &mass_mc_, "mass_mc_/F");
    out_sim_->Branch("px_mc", &px_mc_, "px_mc_/F");
    out_sim_->Branch("py_mc", &py_mc_, "py_mc_/F");
    out_sim_->Branch("pz_mc", &pz_mc_, "pz_mc_/F");
    out_sim_->Branch("g4process", &g4process_, "g4process_/I");

    out_match_ = new TTree("Match", "Match");
    out_match_->Branch("id_rec", &id_rec_, "id_rec_/I");
    out_match_->Branch("id_mc", &id_mc_, "id_mc_/I");

    mother_id_field_id_ = config_->GetBranchConfig(mc_particles_name_).GetFieldId("mother_id");
    g4process_field_id_ = config_->GetBranchConfig(mc_particles_name_).GetFieldId("geant_process_id");
  }
  std::cout << "finish init out" << std::endl;
}

int ConverterOutTree::GetMothersSimId() {
  std::vector<int> daughter_sim_id;
  daughter_sim_id.push_back(rec_to_mc_->GetMatch(daughter_id_1_));
  daughter_sim_id.push_back(rec_to_mc_->GetMatch(daughter_id_2_));
  if (decay_.GetNDaughters() == 3) daughter_sim_id.push_back(rec_to_mc_->GetMatch(daughter_id_3_));

  if (*std::min_element(daughter_sim_id.begin(), daughter_sim_id.end()) < 0)// at least one daughter has no matching with mc
    return -1;

  std::vector<int> mother_sim_id;
  for (int i = 0; i < decay_.GetNDaughters(); i++)
    mother_sim_id.push_back(mc_particles_->GetChannel(daughter_sim_id.at(i)).GetField<int>(mother_id_field_id_));

  if (*std::min_element(mother_sim_id.begin(), mother_sim_id.end()) != *std::max_element(mother_sim_id.begin(), mother_sim_id.end()))// daughters belong to not the same mother
    return -1;

  if (mother_sim_id.at(0) < 0)// mother has negative id
    return -1;

  if (mc_particles_->GetChannel(mother_sim_id.at(0)).GetPid() != pid_)// mother has not PDG which was supposed
    return -1;

  return mother_sim_id.at(0);
}

int ConverterOutTree::DetermineGeneration(int mother_sim_id) {
  int generation = 0;
  int older_id = mother_sim_id;
  while (older_id >= 0) {
    const auto& simtrackolder = mc_particles_->GetChannel(older_id);
    older_id = simtrackolder.GetField<int>(mother_id_field_id_);
    generation++;
  }

  return generation;
}

void ConverterOutTree::MatchWithMc() {
  is_mc_ = 0;
  int mother_id = GetMothersSimId();
  generation_ = DetermineGeneration(mother_id);
  if (generation_ < 1) return;
  is_mc_ = 1;
  id_mc_++;
  out_match_->Fill();

  const auto& particle_mc = mc_particles_->GetChannel(mother_id);
  pid_ = particle_mc.GetPid();
  mass_mc_ = particle_mc.GetMass();
  px_mc_ = particle_mc.GetPx();
  py_mc_ = particle_mc.GetPy();
  pz_mc_ = particle_mc.GetPz();
  g4process_ = particle_mc.GetField<int>(g4process_field_id_);
  out_sim_->Fill();
}

void ConverterOutTree::Finish() {
  out_file_->cd();
  out_file_->Write();
  out_file_->Close();
}
