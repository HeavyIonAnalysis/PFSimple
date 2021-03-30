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
                                                         mass_error_(particle.GetErrMass()) {}

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

  float GetChi2Prim() const { return chi_2_prim_; }
  float GetChi2Geo() const { return chi2_geo_; }
  float GetChi2Topo() const { return chi2_topo_; }
  float GetDistance() const { return distance_; }
  float GetL() const { return l_; }
  float GetLdL() const { return ldl_; }
  float GetCosineTopo() const { return cosine_topo_; }
  int GetNHits() const { return n_hits_; }
  int GetId() const { return id_; }
  bool IsFromPV() const { return is_from_pv_; }

  void SetChi2Prim(float chi_2_prim) { chi_2_prim_ = chi_2_prim; }
  void SetChi2Geo(float chi_2_geo) { chi2_geo_ = chi_2_geo; }
  void SetChi2Topo(float chi_2_topo) { chi2_topo_ = chi_2_topo; }
  void SetDistance(float distance) { distance_ = distance; }
  void SetL(float l) { l_ = l; }
  void SetLdL(float ldl) { ldl_ = ldl; }
  void SetCosineTopo(float cosine_topo) { cosine_topo_ = cosine_topo; }
  void SetNHits(int n_hits) { n_hits_ = n_hits; }
  void SetId(int id) { id_ = id; }
  void SetIsFromPV(bool is_from_pv) { is_from_pv_ = is_from_pv; }

 protected:
  std::vector<int> daughter_ids_{};

  float px_{-1.};
  float py_{-1.};
  float pz_{-1.};
  float mass_{-1.};

  float pt_error_{-1.};
  float phi_error_{-1.};
  float eta_error_{-1.};
  float mass_error_{-1.};

  float chi_2_prim_{-1.}; ///< \f$\chi^2\f$ of the particle to the primary vertex (PV)
  float chi2_geo_{-1.};   ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
  float chi2_topo_{-1.};  ///< \f$\chi^2\f$ of the mother's track to the PV
  float distance_{-1.};   ///< Distance between daughter tracks in their closest approach
  float l_{0.};           ///< Lenght of interpolated track from secondary to primary vertex
  float ldl_{0.};         ///< Distance between primary and secondary vertices divided by error
  float cosine_topo_{-1.};///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV
  int n_hits_{-1};
  int id_{-1};

  bool is_from_pv_{true};
};

#endif// OutputContainer_H