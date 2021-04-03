#ifndef Constants_H
#define Constants_H

#include <limits>

typedef long long Pdg_t;

struct SelectionValues {
  SelectionValues() = default;

  float chi2_prim[3]{-1.f, -1.f, -1.f};
  float distance{std::numeric_limits<float>::max()};
  float l{-1.f};
  float l_over_dl{-1.f};
  float chi2_geo{-1.f};
  float chi2_topo{-1.f};

  bool is_from_PV{false};
};

//float chi_2_prim_{-1.}; ///< \f$\chi^2\f$ of the particle to the primary vertex (PV)
//float chi2_geo_{-1.};   ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
//float chi2_topo_{-1.};  ///< \f$\chi^2\f$ of the mother's track to the PV
//float distance_{-1.};   ///< Distance between daughter tracks in their closest approach
//float l_{0.};           ///< Lenght of interpolated track from secondary to primary vertex
//float ldl_{0.};         ///< Distance between primary and secondary vertices divided by error
//float cosine_topo_{-1.};///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV

enum eParams : short {
  kX = 0,
  kY,
  kZ,
  kPx,
  kPy,
  kPz,
  kNumberOfTrackPars
};

constexpr int NumberOfCovElements = 21;
constexpr int NumberOfFieldPars = 10;

#endif// Constants_H