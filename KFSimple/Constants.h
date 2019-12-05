#ifndef Constants_H
#define Constants_H

  enum eTrackTypes
  {
    kSecPos = 0,
    kSecNeg,
    kPrimPos,
    kPrimNeg,
    kNumberOfTrackTypes
  };
  
  constexpr int pdg_proton = 2212;
  constexpr int pdg_pionMinus = -211;
  constexpr int pdg_lambda = 3122;
  
  constexpr float mass_lambda = 1.115683;
  constexpr float sigma_lambda = 2.7e-3;

#endif // Constants_H