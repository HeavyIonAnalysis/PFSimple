/**
 ** @class OutputContainer
 ** @brief Container with output information about reconstructed particles and geometrical decay parameters (quantities to be cut in order to select particles)
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov, Susanne Glaessel
 **
 ** Each particle candidate is characterized with set of geometrical decay parameters. Depending on the
 ** value of each parameter the candidate is saved or rejected.
 ** In order to save the reconstructed particle, the KFParticle object is used. It contains all
 ** information about the particle (mass, momentum etc), and access to this information is
 ** possible via KFParticle methods.
 **/

#ifndef OutputContainer_H
#define OutputContainer_H

#include "KFParticle.h"

#include "Constants.hpp"

class OutputContainer {
 public:
  OutputContainer() = default;
  explicit OutputContainer(const KFParticle& particle) : px_(particle.GetPx()),
                                                         py_(particle.GetPy()),
                                                         pz_(particle.GetPz()),
                                                         mass_(particle.GetMass()),
                                                         pt_error_(particle.GetErrPt()),
                                                         phi_error_(particle.GetErrPhi()),
                                                         eta_error_(particle.GetErrEta()),
                                                         mass_error_(particle.GetErrMass()),
                                                         x_(particle.GetX()),
                                                         y_(particle.GetY()),
                                                         z_(particle.GetZ()),
                                                         x_error_(particle.GetErrX()),
                                                         y_error_(particle.GetErrY()),
                                                         z_error_(particle.GetErrZ())
                                                         {}

  virtual ~OutputContainer() = default;

  const std::vector<int>& GetDaughterIds() const { return daughter_ids_; }
  float GetPx() const { return px_; }
  float GetPy() const { return py_; }
  float GetPz() const { return pz_; }
  float GetMass() const { return mass_; }
  float GetPtError() const { return pt_error_; }
  float GetPhiError() const { return phi_error_; }
  float GetEtaError() const { return eta_error_; }
  float GetMassError() const { return mass_error_; }

  float GetChi2Prim(int i) const { return values_.chi2_prim[i]; }
  float GetChi2Geo() const { return values_.chi2_geo; }
  float GetChi2Topo() const { return values_.chi2_topo; }
  float GetDistance() const { return values_.distance; }
  float GetL() const { return values_.l; }
  float GetLdL() const { return values_.l_over_dl; }
//  float GetCosineTopo() const { return values_.; }
//  int GetNHits() const { return n_hits_; }
  int GetId() const { return id_; }
  bool IsFromPV() const { return values_.is_from_PV; }

//  void SetNHits(int n_hits) { n_hits_ = n_hits; }
  void SetId(int id) { id_ = id; }

  void SetSelectionValues(const SelectionValues& v){ values_ = v; }

 protected:
  int id_{-1};
  std::vector<int> daughter_ids_{};

  float px_{-1.};
  float py_{-1.};
  float pz_{-1.};
  float mass_{-1.};

  float pt_error_{-1.};
  float phi_error_{-1.};
  float eta_error_{-1.};
  float mass_error_{-1.};

  float x_{-1.};
  float y_{-1.};
  float z_{-1.};
  float x_error_{-1.};
  float y_error_{-1.};
  float z_error_{-1.};

//  int n_hits_{-1};

  SelectionValues values_{};

};

#endif// OutputContainer_H