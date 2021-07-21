#ifndef Constants_H
#define Constants_H

#include <array>
#include <limits>

typedef long long Pdg_t;

struct SelectionValues {
  SelectionValues() = default;

  std::array<float, 3> chi2_prim {{-1.f, -1.f, -1.f}};        ///< \f$\chi^2\f$ of the particle to the primary vertex (PV)
  std::array<float, 3> cos{{-1.f, -1.f, -1.f}};               ///< cosine of angle between daughter track and mother particle
  std::array<float, 2> distance{{-1.f, -1.f}};                ///< Distance between daughter tracks in their closest approach
  std::array<float, 3> distance_sv{{-1.f, -1.f, -1.f}};       ///< Distance between daughter track and SV
  float l{-1.f};                                              ///< Lenght of interpolated track from secondary to primary vertex
  float l_over_dl{-1.f};                                      ///< Distance between primary and secondary vertices divided by error
  float distance_pv{-1.f};                                    ///< Distance between secondary vertex and primary vertex line
  std::array<float, 4> chi2_geo{{-1.f, -1.f, -1.f, -1.f}};    ///< \f$\chi^2\f$ of daughters' tracks in their closest approach (prim & sec mothers)
  std::array<float, 4> chi2_topo{{-1.f, -1.f, -1.f,-1.f}};    ///< \f$\chi^2\f$ of the mother's track to the PV (prim & sec mothers)
  std::array<float, 4> cos_topo{{-1.f, -1.f, -1.f,-1.f}};     ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV (prim & sec mothers)

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

constexpr int NumberOfCovElements = 21;
constexpr int NumberOfFieldPars = 10;

#endif// Constants_H
