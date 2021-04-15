#ifndef Constants_H
#define Constants_H

#include <limits>

typedef long long Pdg_t;

struct SelectionValues {
  SelectionValues() = default;

  float chi2_prim[3]{-1.f, -1.f, -1.f};///< \f$\chi^2\f$ of the particle to the primary vertex (PV)
  float cos[3]{-1.f, -1.f, -1.f};      ///< cosine of angle between daughter track and mother particle
  float distance[2]{-1.f, -1.f};       ///< Distance between daughter tracks in their closest approach
  float l{-1.f};                       ///< Lenght of interpolated track from secondary to primary vertex
  float l_over_dl{-1.f};               ///< Distance between primary and secondary vertices divided by error
  float chi2_geo{-1.f};                ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
  float chi2_topo{-1.f};               ///< \f$\chi^2\f$ of the mother's track to the PV
  float cos_topo{-1.f};                ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV

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

constexpr int NumberOfCovElements = 21;
constexpr int NumberOfFieldPars = 10;

#endif// Constants_H