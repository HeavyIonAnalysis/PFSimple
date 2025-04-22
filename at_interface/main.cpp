#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

#include "AnalysisTree/PlainTreeFiller.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

int main(int argc, char** argv) {

  if (argc < 4) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt dataset parfile.txt\n";
    return EXIT_FAILURE;
  }

  Int_t ndecays = argc - 3;
  
  const std::string& filename   = argv[1];
  std::string outfilename       = std::string(argv[2]) + std::string(".PFSimpleOutput.root");
  std::string outfilename_plain = std::string(argv[2]) + std::string(".PFSimplePlainTree.root");
 
  const Int_t npdgs = 5;
  const Int_t ndaughters_max = 3;

  char name_mother [5];
  Int_t pdg_mother;
  Float_t massPdg, massPdgSigma;
  Int_t ndaughters;
  std::array<Int_t, 3> pdg_1, pdg_2, pdg_3;
  std::array<Float_t, ndaughters_max> chi2prim, cos;
  Float_t dist, distSV, chi2geo, cosopen, chi2topo, costopo, LdL, decaylength, distPVline, invmass;
  std::array<Float_t, ndaughters_max>  chi2geoSM, cosopenSM, chi2topoSM, costopoSM;
  Int_t apply_mass_constr, transport_to_pv, do_not_write_mother;
  Int_t pid_mode;
  Float_t pid_purity;
  std::array<Float_t, npdgs> purity_pdg;
  char atree_name_c [7]; char rec_tracks_name_c [11];
  std::string atree_name, rec_tracks_name;
  Int_t nevents;
  Int_t write_detailed_bg, make_plain_tree;

  std::vector<Decay> decays;
    
  for (int idecay = 0; idecay < ndecays; idecay++) {
      
    TString inputFileInfo = argv[3+idecay];
    
    FILE *inputInfo = fopen(inputFileInfo, "r");
    fscanf(inputInfo, "%*[^\n]%*c");
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_mother);
    fscanf(inputInfo, "%s %*[^\n]%*c", name_mother);
    fscanf(inputInfo, "%f %*[^\n]%*c", &massPdg);
    fscanf(inputInfo, "%f %*[^\n]%*c", &massPdgSigma);

    fscanf(inputInfo, "%i %*[^\n]%*c", &ndaughters);
    
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_1.at(0));
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_1.at(1));
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_1.at(2));
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_2.at(0));
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_2.at(1));
    fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_2.at(2));
    if (ndaughters == 3) {
      fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_3.at(0));
      fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_3.at(1));
      fscanf(inputInfo, "%i %*[^\n]%*c", &pdg_3.at(2));
    }
    fscanf(inputInfo, "%*[^\n]%*c");
    fscanf(inputInfo, "%*[^\n]%*c");
  
    fscanf(inputInfo, "%f %*[^\n]%*c", &chi2prim.at(0));
    fscanf(inputInfo, "%f %*[^\n]%*c", &chi2prim.at(1));
    if (ndaughters == 3)
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2prim.at(2));
    fscanf(inputInfo, "%f %*[^\n]%*c", &cos.at(0));
    fscanf(inputInfo, "%f %*[^\n]%*c", &cos.at(1));
    if (ndaughters == 3)
      fscanf(inputInfo, "%f %*[^\n]%*c", &cos.at(2));
    fscanf(inputInfo, "%*[^\n]%*c");
    fscanf(inputInfo, "%*[^\n]%*c");
  
    fscanf(inputInfo, "%f %*[^\n]%*c", &dist);
    fscanf(inputInfo, "%f %*[^\n]%*c", &distSV);
    if (ndaughters == 3) {
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2geoSM.at(0));
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2geoSM.at(1));
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2geoSM.at(2));
    }
    fscanf(inputInfo, "%f %*[^\n]%*c", &chi2geo);
    if (ndaughters == 3) {
      fscanf(inputInfo, "%f %*[^\n]%*c", &cosopenSM.at(0));
      fscanf(inputInfo, "%f %*[^\n]%*c", &cosopenSM.at(1));
      fscanf(inputInfo, "%f %*[^\n]%*c", &cosopenSM.at(2));
    }
    fscanf(inputInfo, "%f %*[^\n]%*c", &cosopen);
    fscanf(inputInfo, "%*[^\n]%*c");
    fscanf(inputInfo, "%*[^\n]%*c");

    if (ndaughters == 3) {
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2topoSM.at(0));
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2topoSM.at(1));
      fscanf(inputInfo, "%f %*[^\n]%*c", &chi2topoSM.at(2));
      fscanf(inputInfo, "%f %*[^\n]%*c", &costopoSM.at(0));
      fscanf(inputInfo, "%f %*[^\n]%*c", &costopoSM.at(1));
      fscanf(inputInfo, "%f %*[^\n]%*c", &costopoSM.at(2));
      fscanf(inputInfo, "%*[^\n]%*c");
      fscanf(inputInfo, "%*[^\n]%*c");
    }
    
    fscanf(inputInfo, "%f %*[^\n]%*c", &chi2topo);
    fscanf(inputInfo, "%f %*[^\n]%*c", &costopo);
    fscanf(inputInfo, "%f %*[^\n]%*c", &LdL);
    fscanf(inputInfo, "%f %*[^\n]%*c", &decaylength);
    fscanf(inputInfo, "%f %*[^\n]%*c", &distPVline);
    fscanf(inputInfo, "%f %*[^\n]%*c", &invmass);
    fscanf(inputInfo, "%*[^\n]%*c");
    fscanf(inputInfo, "%*[^\n]%*c");

    fscanf(inputInfo, "%i %*[^\n]%*c", &apply_mass_constr);
    fscanf(inputInfo, "%i %*[^\n]%*c", &transport_to_pv);
    fscanf(inputInfo, "%i %*[^\n]%*c", &do_not_write_mother);
    fscanf(inputInfo, "%*[^\n]%*c");
    fscanf(inputInfo, "%*[^\n]%*c");

    if (idecay == 0) {

      fscanf(inputInfo, "%i %*[^\n]%*c", &pid_mode);
      fscanf(inputInfo, "%*[^\n]%*c");
      fscanf(inputInfo, "%*[^\n]%*c");

      fscanf(inputInfo, "%f %*[^\n]%*c", &pid_purity);
      fscanf(inputInfo, "%f %*[^\n]%*c", &purity_pdg.at(0));
      fscanf(inputInfo, "%f %*[^\n]%*c", &purity_pdg.at(1));
      fscanf(inputInfo, "%f %*[^\n]%*c", &purity_pdg.at(2));
      fscanf(inputInfo, "%f %*[^\n]%*c", &purity_pdg.at(3));
      fscanf(inputInfo, "%f %*[^\n]%*c", &purity_pdg.at(4));
      fscanf(inputInfo, "%*[^\n]%*c");
      fscanf(inputInfo, "%*[^\n]%*c");

      fscanf(inputInfo, "%s %*[^\n]%*c", atree_name_c);
      fscanf(inputInfo, "%s %*[^\n]%*c", rec_tracks_name_c);
      fscanf(inputInfo, "%i %*[^\n]%*c", &write_detailed_bg);
      fscanf(inputInfo, "%i %*[^\n]%*c", &nevents);
      fscanf(inputInfo, "%i %*[^\n]%*c", &make_plain_tree);
	
      atree_name = atree_name_c;
      rec_tracks_name = rec_tracks_name_c;
    }

    fclose(inputInfo);
 
    std::vector<Daughter> daughters;
      
    std::vector<Pdg_t> pdg_1_vec;
    if (pdg_1.at(1) != 0 || pdg_1.at(2) != 0)	{
      pdg_1_vec.push_back(pdg_1.at(0));
      for (int i = 1; i < 3; i++)
	if (pdg_1.at(i) != 0) pdg_1_vec.push_back(pdg_1.at(i));
      Daughter daughter_1(pdg_1_vec.at(0), pdg_1_vec);
      daughters.push_back(daughter_1);
    }
    else
      daughters.push_back(pdg_1.at(0));
      
    std::vector<Pdg_t> pdg_2_vec;
    if (pdg_2.at(1) != 0 || pdg_2.at(2) != 0)	{
      pdg_2_vec.push_back(pdg_2.at(0));
      for (int i = 1; i < 3; i++)
	if (pdg_2.at(i) != 0) pdg_2_vec.push_back(pdg_2.at(i));
      Daughter daughter_2(pdg_2_vec.at(0), pdg_2_vec);
      daughters.push_back(daughter_2);
    }
    else
      daughters.push_back(pdg_2.at(0));
    
    if (ndaughters == 3) {
      std::vector<Pdg_t> pdg_3_vec;
      if (pdg_3.at(1) != 0 || pdg_3.at(2) != 0)	{
	pdg_3_vec.push_back(pdg_3.at(0));
	for (int i = 1; i < 3; i++)
	  if (pdg_3.at(i) != 0) pdg_3_vec.push_back(pdg_3.at(i));
	Daughter daughter_3(pdg_3_vec.at(0), pdg_3_vec);
	daughters.push_back(daughter_3);
      }
      else
	daughters.push_back(pdg_3.at(0));
    }

    for (size_t idaughter = 0; idaughter < daughters.size(); ++idaughter) {
      daughters.at(idaughter).CancelCuts();
      if (chi2prim.at(idaughter) != -1) daughters.at(idaughter).SetCutChi2Prim(chi2prim.at(idaughter));
      if (cos.at(idaughter)  != -1) daughters.at(idaughter).SetCutCos(cos.at(idaughter));
    }
    Mother mother(pdg_mother);
    if (massPdg      != -1) mother.SetMassPdg(massPdg);
    if (massPdgSigma != -1) mother.SetMassPdgSigma(massPdgSigma);
    
    mother.CancelCuts();

    // Set cut values
    if (dist        != -1) mother.SetCutDistance(dist);
    if (ndaughters == 3)
      if (distSV      != -1) mother.SetCutDistanceToSV(distSV);
    if (chi2geo     != -1) mother.SetCutChi2Geo(chi2geo);
    if (cosopen     != -1) mother.SetCutCosOpen(cosopen);
    if (chi2topo    != -1) mother.SetCutChi2Topo(chi2topo);
    if (costopo     != -1) mother.SetCutCosTopo(costopo);
    if (LdL         != -1) mother.SetCutLdL(LdL);
    if (decaylength != -1) mother.SetCutDecayLength(decaylength);
    if (distPVline  != -1) mother.SetCutDistancePVLine(distPVline);
    if (invmass     != -1) mother.SetCutInvMass(invmass);

    if (ndaughters == 3) {
      std::vector<Float_t> chi2geoSM_vec;
      for (int i = 0; i < ndaughters; i++)
	if (chi2geoSM.at(i) != -1)
	  chi2geoSM_vec.push_back(chi2geoSM.at(i));
	else break;
      if (chi2geoSM_vec.size() > 0)
	mother.SetCutChi2GeoSM({chi2geoSM_vec});

      std::vector<Float_t> cosopenSM_vec;
      for (int i = 0; i < ndaughters; i++)
	if (cosopenSM.at(i) != -1)
	  cosopenSM_vec.push_back(cosopenSM.at(i));
	else break;
      if (cosopenSM_vec.size() > 0)
	mother.SetCutCosOpenSM({cosopenSM_vec});

      std::vector<Float_t> chi2topoSM_vec;
      for (int i = 0; i < ndaughters; i++)
	if (chi2topoSM.at(i) != -1)
	  chi2topoSM_vec.push_back(chi2topoSM.at(i));
	else break;
      if (chi2topoSM_vec.size() > 0)
	mother.SetCutChi2TopoSM({chi2topoSM_vec});

      std::vector<Float_t> costopoSM_vec;
      for (int i = 0; i < ndaughters; i++)
	if (costopoSM.at(i) != -1)
	  costopoSM_vec.push_back(costopoSM.at(i));
	else break;
      if (costopoSM_vec.size() > 0)
	mother.SetCutCosTopoSM({costopoSM_vec});
    }
    
    Decay decay(name_mother, mother, {daughters});
    if (apply_mass_constr   == 1) decay.SetIsApplyMassConstraint();
    if (transport_to_pv     == 1) decay.SetIsTransportToPV();
    if (do_not_write_mother == 1) decay.SetIsDoNotWriteMother();
    
    decays.push_back(decay);
  }

  auto* man = TaskManager::GetInstance();

  auto* in_converter = new ConverterIn();
  in_converter->SetRecEventHeaderName("RecEventHeader");
  in_converter->SetRecTracksName(rec_tracks_name);
  in_converter->SetSimTracksName("SimParticles");
  in_converter->SetTrackCuts(new Cuts("Cut to reproduce KFPF", {EqualsCut((rec_tracks_name + ".pass_cuts").c_str(), 1)}));

  in_converter->SetPidMode(pid_mode);
  if (pid_purity       != -1) in_converter->SetPidPurity(pid_purity);
  if (purity_pdg.at(0) != -1) in_converter->SetPidPurityProton(purity_pdg.at(0)); // in pid-mode 4 pdg-specific purity possible
  if (purity_pdg.at(1) != -1) in_converter->SetPidPurityPion(purity_pdg.at(1));
  if (purity_pdg.at(2) != -1) in_converter->SetPidPurityKaon(purity_pdg.at(2));
  if (purity_pdg.at(3) != -1) in_converter->SetPidPurityDeuteron(purity_pdg.at(3));
  if (purity_pdg.at(4) != -1) in_converter->SetPidPurityBG(purity_pdg.at(4));
 
  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({decays});

  auto* out_converter = new ConverterOut();
  man->SetOutputName(outfilename, "pTree");
  out_converter->SetSimEventHeaderName("SimEventHeader");
  out_converter->SetRecTracksName(rec_tracks_name);
  out_converter->SetSimTracksName("SimParticles");
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetDecays(decays);
  out_converter->SetIsWriteDetailedBG(write_detailed_bg);
  
  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);
  man->SetVerbosityPeriod(100);
  
  man->Init({filename}, {atree_name});
  man->Run(nevents);// -1 = all events
  man->Finish();
  man->ClearTasks();

  if (make_plain_tree) {
    std::ofstream filelist;
    filelist.open("filelist.txt");
    filelist << outfilename;
    filelist << "\n";
    filelist.close();

    auto* tree_task = new PlainTreeFiller();
    tree_task->SetOutputName(outfilename_plain, "plain_tree");
    std::string branchname_rec = "Candidates";
    tree_task->SetInputBranchNames({branchname_rec});
    tree_task->AddBranch(branchname_rec);
    //tree_task->SetIsPrependLeavesWithBranchName(false); //Uncomment when most recent AnalysisTree is installed in cbmroot

    man->AddTask(tree_task);

    man->Init({"filelist.txt"}, {"pTree"});
    man->Run(nevents);// -1 = all events
    man->Finish();
  }

  return EXIT_SUCCESS;
}
