#ifndef Constants_H
#define Constants_H

#include <array>
#include <limits>

typedef long long Pdg_t;

struct SelectionValues {
  SelectionValues() = default;

  std::array<float, 3> chi2_prim{{-1.f, -1.f, -1.f}};      ///< \f$\chi^2\f$ of the particle to the primary vertex (PV)
  std::array<float, 3> cos{{-1.f, -1.f, -1.f}};            ///< cosine of angle between daughter track and mother particle
  float distance{-1.f};                                    ///< Distance between daughter tracks in their closest approach
  float distance_sv{-1.f};                                 ///< Distance between daughter track and SV
  float l{-1.f};                                           ///< Lenght of interpolated track from secondary to primary vertex
  float l_over_dl{-1.f};                                   ///< Lenght of interpolated track from secondary to primary vertex divided by error
  float distance_pv{-1.f};                                 ///< Distance between secondary vertex and primary vertex line
  std::array<float, 4> chi2_geo{{-1.f, -1.f, -1.f, -1.f}}; ///< \f$\chi^2\f$ of daughters' tracks in their closest approach (prim & sec mothers)
  std::array<float, 4> cos_open{{-1.f, -1.f, -1.f, -1.f}}; ///< cosine of angle between daughter tracks (for mother of 3: cos of 2 daughters with max angle)
  std::array<float, 4> chi2_topo{{-1.f, -1.f, -1.f, -1.f}};///< \f$\chi^2\f$ of the mother's track to the PV (prim & sec mothers)
  std::array<float, 4> cos_topo{{-1.f, -1.f, -1.f, -1.f}}; ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV (prim & sec mothers)

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
constexpr int NumberOfPids = (kNumberOfTrackPars + 1) * kNumberOfTrackPars / 2;

constexpr std::array<int, NumberOfPids> pid_codes_rec = {
    2212, 211, 321, 1000010020, 1};
#endif// Constants_H
