/**
 ** @class DecayContainer
 ** @brief Container with decay parameters.
 ** @authors Susanne Glaessel, Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov
 **
 ** A decay is characterized by: \n
 ** name of the mother, mother pdg, number of daughters, positive daughter pdg, negative daughter pdg, third daughter pdg, mass of the mother, sigma of mass. \n
 ** The user can either set the decay manually in main.cxx or select from a decay in the list giving the mother pdg (list will be extended). \n
 ** Optional: The user can set additional particle pids for a daughter that will be considered as candidates for the reconstruction of the mother. \n
 ** In this case, a daughter is represented by a vector<int> with: \n
 ** 1st position: daughter pid, 2nd position: =1: nuclei pids (pid > 1000000000) will be considered as daughter candidates, =0: no nuclei pids considered, 3rd, 4th, etc. & positions: alternative pids \n
   eg.: decay.SetPdgDaughterPos({2212,1,211,-1})  \n
**/

#ifndef DecayContainer_H
#define DecayContainer_H

#include <string>
#include "TObject.h"
#include <iostream>
#include <map>

//DecayInfo: Class to define variables of the decay list
class DecayList
{
 public:
 DecayList():name_mother_("nullptr"),pdg_mother_(0),ndaughters_(0),pdg_daughter_pos_(0),pdg_daughter_neg_(0),pdg_daughter_third_(0),mass_(0.0), mass_sigma_(0.001), use_alt_pos_(true) {};
  // Constructor with all parameters
 DecayList(std::string name_mother, int pdg_mother, int ndaughters, int pdg_daughter_pos, int pdg_daughter_neg, int pdg_daughter_third, float mass, float mass_sigma, bool use_alt_pos):
  name_mother_(name_mother),pdg_mother_(pdg_mother),ndaughters_(ndaughters),pdg_daughter_pos_(pdg_daughter_pos),pdg_daughter_neg_(pdg_daughter_neg),pdg_daughter_third_(pdg_daughter_third),mass_(mass), mass_sigma_(mass_sigma), use_alt_pos_(use_alt_pos) {};
  ~DecayList() {};

  //Getters
  std::string NameMother()       const { return name_mother_; }
  int         PdgMother()        const { return pdg_mother_; }
  int         NDaughters()       const { return ndaughters_; }
  int         PdgDaughterPos()   const { return pdg_daughter_pos_;}
  int         PdgDaughterNeg()   const { return pdg_daughter_neg_; }
  int         PdgDaughterThird() const { return pdg_daughter_third_; }
  float       Mass()             const { return mass_; }
  float       MassSigma()        const { return mass_sigma_; }
  bool        UseAltPos()        const { return use_alt_pos_; }

 private:
  std::string name_mother_;
  int         pdg_mother_;
  int         ndaughters_;
  int         pdg_daughter_pos_;
  int         pdg_daughter_neg_;
  int         pdg_daughter_third_;
  float       mass_;
  float       mass_sigma_;
  bool        use_alt_pos_;
};

class DecayContainer :public TObject
{
 public:
  DecayContainer()
  {
    DecayList decaylist_[] = 
      {
	//         name     pdg mother   ndaughters   pdg daughter pos,  pdg daughter neg,   pdg daughter third  mass mother   mass sigma mother   use all pos
	DecayList("Lambda", 3122,        2,           2212,              -211,               0,                  1.115683,      2.7e-3,            true),    // Lambda -> p + pi-
	DecayList("H3L",    3004,        3,           2212,              -211,               1000010020,         2.993876,      3.057e-3,          false)   // H3L -> p + pi- + d
      };

    for (int idecay = 0; idecay < ndecays_; idecay++)  
      {
	list_pdg_mother_[idecay] = decaylist_[idecay].PdgMother() ;
	list_name_mother_[idecay] = decaylist_[idecay].NameMother();
	list_ndaughters_[idecay] = decaylist_[idecay].NDaughters();
	list_pdg_daughter_[0][idecay] = decaylist_[idecay].PdgDaughterPos() ;
	list_pdg_daughter_[1][idecay] = decaylist_[idecay].PdgDaughterNeg() ;
	if (list_ndaughters_[idecay] == 3) 
	  list_pdg_daughter_[2][idecay] = decaylist_[idecay].PdgDaughterThird() ;
	list_mass_[idecay] = decaylist_[idecay].Mass();
	list_mass_sigma_[idecay] = decaylist_[idecay].MassSigma();
	list_use_alt_pos_[idecay] = decaylist_[idecay].UseAltPos();
	list_pdg2index_[decaylist_[idecay].PdgMother()] = idecay;
      }
  }
  
  virtual ~DecayContainer() = default;

  // Set decay from list based on mother pdg
  bool SetDecay(const int pdg_mother)
  {
    std::map<int, int>::iterator pdg2index;
    pdg2index=list_pdg2index_.find(pdg_mother);
    int idecay;
    if (pdg2index != list_pdg2index_.end())
      {
	idecay = pdg2index->second;
        pdg_mother_ = pdg_mother;
	name_mother_ = list_name_mother_[idecay];
	ndaughters_ = list_ndaughters_[idecay];
	pdg_daughter_pos_.at(0) = list_pdg_daughter_[0][idecay];
	pdg_daughter_neg_.at(0) = list_pdg_daughter_[1][idecay];
	if (ndaughters_ == 3) 
	  pdg_daughter_third_.at(0) = list_pdg_daughter_[2][idecay];
	mass_ = list_mass_[idecay];
	mass_sigma_ = list_mass_sigma_[idecay];
	if (list_use_alt_pos_[idecay] == true) {
	  pdg_daughter_pos_.resize(4);
	  pdg_daughter_pos_.at(1) = 1;
	  pdg_daughter_pos_.at(2) = 211;
	  pdg_daughter_pos_.at(3) = -1;
	 }
	return true;
      }
    else
      {
	std::cout<<"Decay of particle "<<pdg_mother<<" is not listed!\n  --> Please set decay manually."<<std::endl;
	return false;
      }
  };

  //  decay parameters setters
  void SetNameMother(std::string value) {name_mother_ = value;};
  void SetPdgMother(int value) {pdg_mother_ = value;};
  void SetNdaughters(int value) {ndaughters_ = value;};
  void SetPdgDaughterPos(int value) {pdg_daughter_pos_.at(0) = value;};
  void SetPdgDaughterPos(std::vector<int> value) {pdg_daughter_pos_ = value;};
  void SetPdgDaughterNeg(int value) {pdg_daughter_neg_.at(0) = value;};
  void SetPdgDaughterNeg(std::vector<int> value) {pdg_daughter_neg_ = value;};
  void SetPdgDaughterThird(int value) {pdg_daughter_third_.at(0) = value;};
  void SetPdgDaughterThird(std::vector<int> value) {pdg_daughter_third_ = value;};
  
  //  decay parameters getters
  std::string GetNameMother() const {return name_mother_;};
  int GetPdgMother() const {return pdg_mother_;};
  int GetNdaughters() const {return ndaughters_;};
  int GetPdgDaughterPos() const {return pdg_daughter_pos_.at(0);};
  std::vector<int> GetPdgDaughterPosCandidates() const {return pdg_daughter_pos_;};
  int GetPdgDaughterNeg() const {return pdg_daughter_neg_.at(0);};
   std::vector<int> GetPdgDaughterNegCandidates() const {return pdg_daughter_neg_;};
  int GetPdgDaughterThird() const {return pdg_daughter_third_.at(0);};
   std::vector<int> GetPdgDaughterThirdCandidates() const {return pdg_daughter_third_;};
  
 protected:
  // Decay parameters with their default values
  std::string name_mother_{"Lambda"};                       ///< Name of the mother
  int pdg_mother_{3122};                                    ///< Mother pdg
  int ndaughters_{2};                                       ///< Number of daughters              
  std::vector<int> pdg_daughter_pos_{2212};                 ///< Pid of first positive daughter (1st position),
                                                            ///< optional: additional pids that will be considered as candidates for the reconstruction of the mother:
                                                            ///< 2nd position: = 1, light nuclei as optional candidates, =0: no light nuclei as optional candidates}
                                                            ///< 3rd, 4th, etc. position: optional pdgs for positive daughter,
  std::vector<int> pdg_daughter_neg_{-211};                 ///< Pid of first negative daughter (1st position) (for 2nd, 3rd, etc, see comment above)
  std::vector<int> pdg_daughter_third_{0};                  ///< Pid of third (pos or neg) daughter (1st position) (for 2nd, 3rd, etc,, see comment above)
  float mass_{1.115683};                                    ///< Mean mass of mother particle
  float mass_sigma_{2.7e-3};                                ///< Sigma mass of mother particle

 private:
  static const int ndecays_=2;
  std::string list_name_mother_[ndecays_];
  int list_pdg_mother_[ndecays_];
  int list_ndaughters_[ndecays_];
  int list_pdg_daughter_[3][ndecays_];
  float list_mass_[ndecays_];
  float list_mass_sigma_[ndecays_];
  bool list_use_alt_pos_[ndecays_];
  std::map<int, int> list_pdg2index_;   
	
};
#endif//DecayContainer_H
