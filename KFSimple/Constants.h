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

  enum eField : short
  {
    kCx0 = 0,
    kCx1,
    kCx2,
    kCy0,
    kCy1,
    kCy2,
    kCz0,
    kCz1,
    kCz2,
    kZ0,
    kNumberOfFieldPars
  };

  enum eParams : short
  {
    kX = 0,
    kY,
    kZ,
    kPx,
    kPy,
    kPz,
    kNumberOfTrackPars
  };


  constexpr int pdg_proton = 2212;
  constexpr int pdg_pionMinus = -211;
  constexpr int pdg_lambda = 3122;
  
  constexpr float mass_lambda = 1.115683;
  constexpr float sigma_lambda = 2.7e-3;

#endif // Constants_H