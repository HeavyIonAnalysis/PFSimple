#ifndef Constants_H
#define Constants_H

#include <array>
#include <limits>

typedef long long Pdg_t;

struct SelectionValues {
  SelectionValues() = default;

  std::array<float, 3> chi2_prim {{-1.f, -1.f, -1.f}};  ///< \f$\chi^2\f$ of the particle to the primary vertex (PV)
  std::array<float, 3> cos{{-1.f, -1.f, -1.f}};         ///< cosine of angle between daughter track and mother particle
  std::array<float, 2> distance{{-1.f, -1.f}};          ///< Distance between daughter tracks in their closest approach
  float l{-1.f};                                        ///< Lenght of interpolated track from secondary to primary vertex
  float l_over_dl{-1.f};                                ///< Distance between primary and secondary vertices divided by error
  float chi2_geo{-1.f};                                 ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
  float chi2_topo{-1.f};                                ///< \f$\chi^2\f$ of the mother's track to the PV
  float cos_topo{-1.f};                                 ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV
  float invmassdisc{-1.f};                              ///< Discrepancy of the V0 candidate invariant mass from the PDG value in terms of characteristic sigma (hardcoded)
  float chi2_prim_mother{-1.f};                        ///< \f$\chi^2\f$ of the mother to the primary vertex (PV)


  bool is_from_PV{false};
};

enum eParams : short {
  kX = 0,
  kY,
  kZ,
  kPx,
  kPy,
  kPz,
  kNumberOfTrackPars
};

constexpr float huge_value = std::numeric_limits<float>::max();
constexpr float lambda_mass = 1.115683;
constexpr float lambda_mass_sigma = 1.5e-3;

constexpr int NumberOfCovElements = 21;
constexpr int NumberOfFieldPars = 10;

constexpr int NumberOfPids = 5;

constexpr std::array<int, NumberOfPids> pid_codes_rec = { 
	     2212, 211, 321, 1000010020, 1
};
#endif// Constants_H