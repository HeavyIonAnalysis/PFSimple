#include "SimpleFinder.h"

#include "KFParticleSIMD.h"
#include "KFParticle.h"


void SimpleFinder::Init(const KFPTrackVector& tracks, const KFVertex& pv)
{
  tracks_ = tracks;
  prim_vx_ = pv;
}

void SimpleFinder::Init(const InputContainer& input)
{
  KFPTrackVector track_tmp;
  const std::vector<KFParticle>& tracks = input.GetTracks();
  track_tmp.Resize(tracks.size());

  for(int iTr=0; iTr<tracks.size(); iTr++)
  {
    for(Int_t iP=0; iP<6; iP++)
      track_tmp.SetParameter(tracks[iTr].GetParameter(iP), iP, iTr);
    for(Int_t iC=0; iC<21; iC++)
      track_tmp.SetCovariance(tracks[iTr].GetCovariance(iC), iC, iTr);
    for(Int_t iF=0; iF<10; iF++)
      track_tmp.SetFieldCoefficient(tracks[iTr].GetFieldCoeff()[iF], iF, iTr);
    track_tmp.SetPDG(tracks[iTr].GetPDG(), iTr);
    track_tmp.SetQ(tracks[iTr].GetQ(), iTr);
    track_tmp.SetPVIndex(-1, iTr);
    track_tmp.SetId(tracks[iTr].Id(), iTr);
  }

  Init(track_tmp, input.GetVertex());
  SetDecay(input.GetDecay());
  SetCuts(input.GetCuts());
}

void SimpleFinder::SortTracks()
{
  /**
   *   * Sorts tracks' indices into 4 groups:\n
   * 1) secondary positive\n
   * 2) secondary negative\n
   * 3) primary positive\n
   * 4) primary negative\n
   */
  const int Size = tracks_.Size();
  
  for(int iTr = 0; iTr < Size; iTr++)
  {
    if(tracks_.PVIndex()[iTr] < 0)          // secondary
    {
      if(tracks_.Q()[iTr] > 0)              // secondary positive
        trIndex_[kSecPos].push_back(iTr);
      else                                  // secondary negative
        trIndex_[kSecNeg].push_back(iTr);
    }
    else                                    // primary
    {
      if(tracks_.Q()[iTr] > 0)              // secondary positive
        trIndex_[kPrimPos].push_back(iTr);
      else                                  // secondary negative
        trIndex_[kPrimNeg].push_back(iTr);
    }
  }
}

/*float SimpleFinder::CalculateChiToPrimaryVertex(const KFPTrack &track, const int pid) const
{
  // Scalar version. Not optimal (was written earlier and was not improved after SIMD'ization of KFParticle package)
  KFParticle tmpPart(track, pid);
  const float point[3] = {prim_vx_.X(), prim_vx_.Y(), prim_vx_.Z()};
  tmpPart.TransportToPoint(point);
  
  return tmpPart.GetDeviationFromVertex(prim_vx_);
}*/




/** Functions for the reconstruction of the mother from two daughters **/

float SimpleFinder::CalculateChiToPrimaryVertex(const KFPTrack &track, int pid) const
{
  // SIMD'ized version
  KFParticle tmpPart(track, pid);
  KFParticleSIMD tmpPartSIMD(tmpPart);
  
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);  
  const float_v point[3] = {prim_vx_Simd.X(), prim_vx_Simd.Y(), prim_vx_Simd.Z()};
  tmpPartSIMD.TransportToPoint(point);

  float_v chi2vec = tmpPartSIMD.GetDeviationFromVertex(prim_vx_Simd);
  
  return chi2vec[0];
}

void SimpleFinder::CalculateParamsInPCA(const KFPTrack &track1, int pid1, const KFPTrack &track2, int pid2, std::array<float, 8> &pars1, std::array<float, 8> &pars2) const
{
  KFParticle particle1(track1, pid1);
  KFParticle particle2(track2, pid2);
  KFParticleSIMD particleSIMD1(particle1);    // the same particle is copied to each SIMD element
  KFParticleSIMD particleSIMD2(particle2);
  
  float_v dS[2];
  float_v params1[8], params2[8];
  particleSIMD1.GetDStoParticleFast(particleSIMD2, dS);
  particleSIMD1.TransportFast(dS[0], params1);
  particleSIMD2.TransportFast(dS[1], params2);
  
  float_v parbuf;
  for(int i=0; i<8; i++)
  {
    parbuf = params1[i];
    pars1.at(i) = parbuf[0];
    parbuf = params2[i];
    pars2.at(i) = parbuf[0];    
  }
}

float SimpleFinder::CalculateDistanceBetweenParticles(const std::array<float, 8> &pars1, const std::array<float, 8> &pars2) const
{
  const float dx = pars1.at(0) - pars2.at(0);
  const float dy = pars1.at(1) - pars2.at(1);
  const float dz = pars1.at(2) - pars2.at(2);
  const float dr = std::sqrt(dx*dx+dy*dy+dz*dz);
  
  return dr;
}

float SimpleFinder::CalculateCosMomentumSum(const std::array<float, 8> &pars1, const std::array<float, 8> &pars2) const
{
  /**
   * Find cosine bitween daughter1 and mother momenta
   */
  const std::array<float, 3> P1 = {pars1.at(3), pars1.at(4), pars1.at(5)};
  const std::array<float, 3> P2 = {pars2.at(3), pars2.at(4), pars2.at(5)};
  const std::array<float, 3> PSum = {P1.at(0)+P2.at(0), P1.at(1)+P2.at(1), P1.at(2)+P2.at(2)};
  
  return (P1.at(0)*PSum.at(0) + P1.at(1)*PSum.at(1) + P1.at(2)*PSum.at(2)) /
         (sqrt(P1.at(0)*P1.at(0) + P1.at(1)*P1.at(1) + P1.at(2)*P1.at(2)) * std::sqrt(PSum.at(0)*PSum.at(0) + PSum.at(1)*PSum.at(1) + PSum.at(2)*PSum.at(2)));
}

KFParticleSIMD SimpleFinder::ConstructMother(const KFPTrack &track1, const int pid1, const KFPTrack &track2, const int pid2) const
{
  KFParticle particle1(track1, pid1);
  KFParticle particle2(track2, pid2);
  particle1.SetId(track1.Id());
  particle2.SetId(track2.Id());
  KFParticleSIMD particleSIMD1(particle1);    // the same particle is copied to each SIMD element
  KFParticleSIMD particleSIMD2(particle2);
  
  float_v ds[2] = {0.f,0.f};
  float_v dsdr[4][6];

  particleSIMD1.GetDStoParticle( particleSIMD2, ds, dsdr );
  particleSIMD1.TransportToDS(ds[0], dsdr[0]);
  particleSIMD2.TransportToDS(ds[1], dsdr[3]);
  const KFParticleSIMD* vDaughtersPointer[2] = {&particleSIMD1, &particleSIMD2};
  
  KFParticleSIMD mother;
  mother.Construct(vDaughtersPointer, 2, nullptr);
  
  return mother;
}

float SimpleFinder::CalculateChi2Geo(const KFParticleSIMD& mother) const
{
 float_v chi2 = mother.Chi2()/simd_cast<float_v>(mother.NDF());
 
 return chi2[0];
}

void SimpleFinder::CalculateMotherProperties(const KFParticleSIMD& mother, float &l, float &ldl, int &isFromPV) const
{
  float_v l_Simd;
  float_v dl_Simd;
  float_m isFromPV_Simd;
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);

  mother.GetDistanceToVertexLine(prim_vx_Simd, l_Simd, dl_Simd, &isFromPV_Simd);
  
  KFParticleSIMD motherTopo = mother;
  motherTopo.SetProductionVertex(prim_vx_Simd);
  motherTopo.KFParticleBaseSIMD::GetDecayLength(l_Simd, dl_Simd);
  
  l = l_Simd[0];
  ldl = l_Simd[0]/dl_Simd[0];
  isFromPV = isFromPV_Simd[0];
}

float SimpleFinder::CalculateCosTopo(const KFParticleSIMD& mother) const
{
  const float x_mother = mother.GetX()[0];
  const float y_mother = mother.GetY()[0];
  const float z_mother = mother.GetZ()[0];

  const float px_mother = mother.GetPx()[0];
  const float py_mother = mother.GetPy()[0];
  const float pz_mother = mother.GetPz()[0];

  const float delta_x = x_mother - prim_vx_.GetX();
  const float delta_y = y_mother - prim_vx_.GetY();
  const float delta_z = z_mother - prim_vx_.GetZ();

  const float sp = delta_x*px_mother + delta_y*py_mother + delta_z*pz_mother;
  const float norm = std::sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z) * std::sqrt(px_mother*px_mother + py_mother*py_mother + pz_mother*pz_mother);
  
  return sp/norm;  
}

float SimpleFinder::CalculateChi2Topo(const KFParticleSIMD& mother) const
{
  KFParticleSIMD motherTopo = mother;
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);
  motherTopo.SetProductionVertex(prim_vx_Simd);
  const float_v& chi2 = motherTopo.GetChi2()/simd_cast<float_v>(motherTopo.GetNDF());
  
  return chi2[0];
}

void SimpleFinder::SaveParticle(const OutputContainer& Mother)
{
  vec_mass_.push_back(mass_);
  vec_mother_.push_back(Mother);
}



/** Additional functions for adding a third daughter to the mother**/

void SimpleFinder::CalculateCoordinatesSecondaryVertex(const std::array<float, 8> &pars1, const std::array<float, 8> &pars2, std::array<float_v, 3> &sv) const
{

  //** Calculate coordinates of the secondary vertex of the first two daughters
  sv.at(0) = (pars1.at(0) + pars2.at(0))/2;
  sv.at(1) = (pars1.at(1) + pars2.at(1))/2;
  sv.at(2) = (pars1.at(2) + pars2.at(2))/2;
  
}

void SimpleFinder::CalculateParamsInSecondaryVertex(const KFParticleSIMD &particleSIMD1, const std::array<float_v, 3> sec_vx, std::array<float, 8> &pars1) const
{  

  //** Calculate parameters of a daughter in the PCA to the secondary vertex
  
  float_v dS;
  float_v dsdr[6];
  float_v params1[8];
  float_v xyz[3];

  for (int i=0; i<3; i++) xyz[i]=sec_vx.at(i);
       
  dS = particleSIMD1.GetDStoPoint(xyz,dsdr);

  particleSIMD1.TransportFast(dS,params1);
  
  float_v parbuf;
  for(int i=0; i<8; i++)
  {
    parbuf = params1[i];
    pars1.at(i) = parbuf[0];  
  }
}

float SimpleFinder::CalculateDistanceToSecondaryVertex(const std::array<float, 8> &pars1, std::array<float_v, 3> &sec_vx) const
{

  //** Calculate the distance of a daughter and the secondary vertex in the PCA

  std::array<float, 3> sv;
  float_v parbuf;

  for(int i=0; i<3; i++)
  {
    parbuf = sec_vx.at(i);
    sv.at(i) = parbuf[0];  
  }

  const float dx = pars1.at(0) - sv.at(0);
  const float dy = pars1.at(1) - sv.at(1);
  const float dz = pars1.at(2) - sv.at(2);
  const float dr = sqrt(dx*dx+dy*dy+dz*dz);
  
  return dr;
  
}

float SimpleFinder::CalculateCosMomentumSumThird(const std::array<float, 8> &pars1, const std::array<float, 8> &pars2, const std::array<float, 8> &pars3) const
{
  /**
   * Find cosine between the momenta of the third daughter and the mother 
   */
  const std::array<float, 3> P1 = {pars1.at(3), pars1.at(4), pars1.at(5)};
  const std::array<float, 3> P2 = {pars2.at(3), pars2.at(4), pars2.at(5)};
  const std::array<float, 3> P3 = {pars3.at(3), pars3.at(4), pars3.at(5)};
  const std::array<float, 3> PSum = {P1.at(0)+P2.at(0)+P3.at(0), P1.at(1)+P2.at(1)+P3.at(1), P1.at(2)+P2.at(2)+P3.at(2)};
  
  return (P3.at(0)*PSum.at(0) + P3.at(1)*PSum.at(1) + P3.at(2)*PSum.at(2)) /
         (sqrt(P3.at(0)*P3.at(0) + P3.at(1)*P3.at(1) + P3.at(2)*P3.at(2)) * sqrt(PSum.at(0)*PSum.at(0) + PSum.at(1)*PSum.at(1) + PSum.at(2)*PSum.at(2)));
}

KFParticleSIMD SimpleFinder::ConstructMotherThree(KFParticleSIMD &particleSIMD1, KFParticleSIMD &particleSIMD2, KFParticleSIMD &particleSIMD3,const std::array<float_v, 3> sec_vx) const
{

  //** Construct the mother from three daughters
  
  float_v ds[3] = {0.f,0.f,0.f};
  float_v dsdr1[6], dsdr2[6], dsdr3[6];
  float_v xyz[3];

  for (int i=0; i<3; i++) xyz[i]=sec_vx.at(i);

  ds[0] = particleSIMD1.GetDStoPoint(xyz,dsdr1);
  ds[1] = particleSIMD2.GetDStoPoint(xyz,dsdr2);
  ds[2] = particleSIMD3.GetDStoPoint(xyz,dsdr3);
  
  particleSIMD1.TransportToDS(ds[0], dsdr1);
  particleSIMD2.TransportToDS(ds[1], dsdr2);
  particleSIMD3.TransportToDS(ds[2], dsdr3);
  
  const KFParticleSIMD* vDaughtersPointer[3] = {&particleSIMD1, &particleSIMD2, &particleSIMD3};
  
  KFParticleSIMD mother;
  //mother.SetConstructMethod(construct_method);
  //if(construct_method == 2) mother.SetMassHypo(mass_H3L);
  mother.Construct(vDaughtersPointer, 3, nullptr);
  
  return mother;
}

void SimpleFinder::FindParticles()
{
  /*
   * The main function which performs the mother-candidate selection algorithm for two or three daugthers.
   */

  // Reconstruction of mother candidates with two daugthers / with the first two of three daughters
  
  int nSecPoses = trIndex_[kSecPos].size();
  int nSecNegs  = trIndex_[kSecNeg].size();
  
  int Ntwo = 0;
  int Nthree = 0;
  
  for(int iSecPos=0; iSecPos<nSecPoses; iSecPos++)
  {
    for(int iSecNeg=0; iSecNeg<nSecNegs; iSecNeg++)
    {

      KFPTrack trackPos;
      tracks_.GetTrack(trackPos, trIndex_[kSecPos][iSecPos]);
      int pidPos = tracks_.PDG()[trIndex_[kSecPos][iSecPos]];
      std::vector<int> pidCandidates;
      pidCandidates.clear();
      pidCandidates = decay_.GetPdgDaughterPosCandidates();
      
      if (pidCandidates.size() > 1 && pidCandidates.at(1) == 1 && pidPos > 1000000000)
      	pidPos = pidCandidates.at(0);

      if (pidCandidates.size() > 2)
	for (int icandidate=2; icandidate < pidCandidates.size() ;icandidate++)
	  if (pidPos == pidCandidates.at(icandidate))
	    pidPos = pidCandidates.at(0);
      
      if(pidPos != pidCandidates.at(0)) continue;
      
      KFPTrack trackNeg;
      tracks_.GetTrack(trackNeg, trIndex_[kSecNeg][iSecNeg]);
      int pidNeg = tracks_.PDG()[trIndex_[kSecNeg][iSecNeg]];
      
      pidCandidates.clear();
      pidCandidates = decay_.GetPdgDaughterNegCandidates();
      
      if (pidCandidates.size() > 1 && pidCandidates.at(1) == 1 && pidNeg < -1000000000)
      	pidNeg = pidCandidates.at(0);

      if (pidCandidates.size() > 2)
	for (int icandidate=2; icandidate < pidCandidates.size() ;icandidate++)
	  if (pidNeg == pidCandidates.at(icandidate))
	    pidNeg = pidCandidates.at(0);
         
      if(pidNeg != pidCandidates.at(0)) continue;
      
      OutputContainer mother;

      mother.SetNHitsPos(tracks_.NPixelHits()[trIndex_[kSecPos][iSecPos]]);
      mother.SetNHitsNeg(tracks_.NPixelHits()[trIndex_[kSecNeg][iSecNeg]]);
            
      mother.SetChi2PrimPos(CalculateChiToPrimaryVertex(trackPos, pidPos));
      if(((cuts_.GetIsApplyCutChi2PrimPos())&&(mother.GetChi2PrimPos() <= cuts_.GetCutChi2PrimPos())) || mother.GetChi2PrimPos()!=mother.GetChi2PrimPos()) continue;
      
      mother.SetChi2PrimNeg(CalculateChiToPrimaryVertex(trackNeg, pidNeg));
      if(((cuts_.GetIsApplyCutChi2PrimNeg())&&(mother.GetChi2PrimNeg() <= cuts_.GetCutChi2PrimNeg())) || mother.GetChi2PrimNeg()!=mother.GetChi2PrimNeg()) continue;
                  
      std::array<float, 8> pars1{};
      std::array<float, 8> pars2{};
      CalculateParamsInPCA(trackNeg, pidNeg, trackPos, pidPos, pars1, pars2);
      
      mother.SetDistance(CalculateDistanceBetweenParticles(pars1, pars2));
      if(((cuts_.GetIsApplyCutDistance())&&(mother.GetDistance() >= cuts_.GetCutDistance())) || mother.GetDistance()!=mother.GetDistance()) continue;

      mother.SetCosineDaughterPos(CalculateCosMomentumSum(pars2, pars1));
      mother.SetCosineDaughterNeg(CalculateCosMomentumSum(pars1, pars2));
      
      if(((cuts_.GetIsApplyCutCosineDaughterPos())&&(mother.GetCosineDaughterPos() < cuts_.GetCutCosineDaughterPos())) || ((cuts_.GetIsApplyCutCosineDaughterNeg())&&(mother.GetCosineDaughterNeg() < cuts_.GetCutCosineDaughterNeg()))
	 || mother.GetCosineDaughterPos()!=mother.GetCosineDaughterPos() || mother.GetCosineDaughterNeg()!=mother.GetCosineDaughterNeg()) continue;
            
      KFParticleSIMD motherTwoKF = ConstructMother(trackNeg, pidNeg, trackPos, pidPos);
      
      mother.SetChi2Geo(CalculateChi2Geo(motherTwoKF));
      
      if(!finite(mother.GetChi2Geo()) || mother.GetChi2Geo() < 0.) continue;
      if((cuts_.GetIsApplyCutChi2Geo())&&(mother.GetChi2Geo() >= cuts_.GetCutChi2Geo())) continue;
      
      float l, ldl;
      int isfrompv = -1;
      CalculateMotherProperties(motherTwoKF, l, ldl, isfrompv);
      mother.SetL(l);
      mother.SetLdL(ldl);
      mother.SetIsFromPV(isfrompv);

      if(((cuts_.GetIsApplyCutLUp())&&(mother.GetL() >= cuts_.GetCutLUp())) || mother.GetL()!=mother.GetL()) continue;
      if(((cuts_.GetIsApplyCutLdL())&&(mother.GetLdL() <= cuts_.GetCutLdL())) || mother.GetLdL()!=mother.GetLdL()) continue;
      if((cuts_.GetIsApplyCutLDown())&&(mother.GetL() <= cuts_.GetCutLDown())) continue;
      //if((cuts_.GetIsApplyCutIsFromPV())&&(mother.GetIsFromPV() == cuts_.GetCutIsFromPV())) continue;
      
      mother.SetCosineTopo(CalculateCosTopo(motherTwoKF));

      if(((cuts_.GetIsApplyCutCosineTopoDown())&&(mother.GetCosineTopo() <= cuts_.GetCutCosineTopoDown()))|| mother.GetCosineTopo()!=mother.GetCosineTopo()) continue;
      if((cuts_.GetIsApplyCutCosineTopoUp())&&(mother.GetCosineTopo() >= cuts_.GetCutCosineTopoUp())) continue;
      
      mother.SetChi2Topo(CalculateChi2Topo(motherTwoKF));
      
      if (decay_.GetNdaughters() == 2)
	if((cuts_.GetIsApplyCutChi2Topo() && mother.GetChi2Topo() > cuts_.GetCutChi2Topo()) || mother.GetChi2Topo()!=mother.GetChi2Topo()) continue;

      if (decay_.GetNdaughters() == 3)
	 if((cuts_.GetIsApplyCutChi2Topo() && mother.GetChi2Topo() < cuts_.GetCutChi2Topo()) || mother.GetChi2Topo()!=mother.GetChi2Topo()) continue;

      KFParticle particle;
      float mass_err; // unused
      
      if (decay_.GetNdaughters() == 2)
	{
	
	  motherTwoKF.GetKFParticle(particle, 0);     
	  particle.GetMass(mass_, mass_err);     
	  particle.SetPDG(decay_.GetPdgMother());	
	  mother.SetParticle(particle);
	  SaveParticle(mother);
	  Ntwo++;
	}
      
      if (decay_.GetNdaughters() == 3)
	
	//** Add the third daughter to the mother candidates

	{
	  
	  KFParticle particlePos(trackPos,pidPos);
	  KFParticle particleNeg (trackNeg, pidNeg);
	  particlePos.SetId(trackPos.Id());
	  particleNeg.SetId(trackNeg.Id());
	  KFParticleSIMD particleSIMDPos(particlePos);
	  KFParticleSIMD particleSIMDNeg(particleNeg);

	  for(int iThird=0; iThird<nSecPoses; iThird++)
	    {
	      KFPTrack trackThird;
	      int pidThird;
	      if (decay_.GetPdgDaughterThird() > 0) {
		tracks_.GetTrack(trackThird, trIndex_[kSecPos][iThird]);
		pidThird = tracks_.PDG()[trIndex_[kSecPos][iThird]];
		mother.SetNHitsThird(tracks_.NPixelHits()[trIndex_[kSecPos][iThird]]);
	      }
	      if (decay_.GetPdgDaughterThird() < 0) {
		tracks_.GetTrack(trackThird, trIndex_[kSecNeg][iThird]);
		pidThird = tracks_.PDG()[trIndex_[kSecNeg][iThird]];
		mother.SetNHitsThird(tracks_.NPixelHits()[trIndex_[kSecNeg][iThird]]);
	      }

	      pidCandidates.clear();
	      pidCandidates = decay_.GetPdgDaughterThirdCandidates();
	      
	      if (pidCandidates.size() > 1 && pidCandidates.at(1) == 1 && ((pidCandidates.at(0) > 0 && pidThird > 1000000000) || (pidCandidates.at(0) < 0 && pidThird < -1000000000)))
		pidThird = pidCandidates.at(0);

	      if (pidCandidates.size() > 2)
		for (int icandidate=2; icandidate < pidCandidates.size() ;icandidate++)
		  if (pidThird == pidCandidates.at(icandidate))
		    pidThird = pidCandidates.at(0);

		 
	      if(pidThird != pidCandidates.at(0)) continue;
			      
	      mother.SetChi2PrimThird(CalculateChiToPrimaryVertex(trackThird, pidThird));

	      if(((cuts_.GetIsApplyCutChi2PrimThird()) && (mother.GetChi2PrimThird() <= cuts_.GetCutChi2PrimThird())) || mother.GetChi2PrimThird()!=mother.GetChi2PrimThird()) continue;
	    
	      KFParticle particleThird(trackThird,pidThird);
	      particleThird.SetId(trackThird.Id());
	      KFParticleSIMD particleSIMDThird(particleThird);

	      std::array<float_v, 3> sec_vx;
	      CalculateCoordinatesSecondaryVertex(pars1,pars2,sec_vx);

	      std::array<float, 8> pars3{};
	      CalculateParamsInSecondaryVertex(particleSIMDThird, sec_vx, pars3);

	      mother.SetDistanceThird(CalculateDistanceToSecondaryVertex(pars3, sec_vx));
	       
	      if(((cuts_.GetIsApplyCutDistanceThird()) && (mother.GetDistanceThird() >= cuts_.GetCutDistanceThird())) || mother.GetDistanceThird()!=mother.GetDistanceThird()) continue;

	      mother.SetCosineDaughterThird(CalculateCosMomentumSumThird(pars1, pars2, pars3));

	      if(((cuts_.GetIsApplyCutCosineDaughterThird()) &&(mother.GetCosineDaughterThird() < cuts_.GetCutCosineDaughterThird())) || mother.GetCosineDaughterThird()!=mother.GetCosineDaughterThird()) continue;

	      KFParticleSIMD motherThreeKF = ConstructMotherThree(particleSIMDPos, particleSIMDNeg, particleSIMDThird,sec_vx);

	      mother.SetChi2GeoThree(CalculateChi2Geo(motherThreeKF));
	      
	      if(!finite(mother.GetChi2GeoThree()) || mother.GetChi2GeoThree() < 0.) continue;
	      if((cuts_.GetIsApplyCutChi2GeoThree())&&(mother.GetChi2GeoThree() >= cuts_.GetCutChi2GeoThree())) continue;
      
	      l=-1;
	      ldl=-1;
	      isfrompv = -1;
	      CalculateMotherProperties(motherThreeKF, l, ldl, isfrompv);
	      mother.SetL(l);
	      mother.SetLdL(ldl);
	      mother.SetIsFromPV(isfrompv);

	      if(((cuts_.GetIsApplyCutLThreeDown()) && (mother.GetL() <= cuts_.GetCutLThreeDown())) || mother.GetL()!=mother.GetL()) continue;
	      if((cuts_.GetIsApplyCutLThreeUp())&&(mother.GetL() >= cuts_.GetCutLThreeUp())) continue;
	      if((cuts_.GetIsApplyCutIsFromPVThree())&&(mother.GetIsFromPV() == cuts_.GetCutIsFromPVThree())) continue;
	      
	      mother.SetCosineTopoThree(CalculateCosTopo(motherThreeKF));
	      
	      if(((cuts_.GetIsApplyCutCosineTopoThreeDown())&&(mother.GetCosineTopoThree() <= cuts_.GetCutCosineTopoThreeDown())) || mother.GetCosineTopoThree() != mother.GetCosineTopoThree()) continue;
	      if((cuts_.GetIsApplyCutCosineTopoThreeUp())&&(mother.GetCosineTopoThree() >= cuts_.GetCutCosineTopoThreeUp())) continue;
 
	      mother.SetChi2TopoThree(CalculateChi2Topo(motherThreeKF));
	      
	      if((cuts_.GetIsApplyCutChi2TopoThree() && mother.GetChi2TopoThree() > cuts_.GetCutChi2TopoThree()) || mother.GetChi2TopoThree()!=mother.GetChi2TopoThree()) continue;
	      
	      motherThreeKF.GetKFParticle(particle, 0);
      
	      particle.GetMass(mass_, mass_err);
	      
	      particle.SetPDG(decay_.GetPdgMother());
	      
	      Nthree++;
	      mother.SetParticle(particle);
	      SaveParticle(mother);
	}  
      
      }
    }
  }
  
  //std::cout <<"N="<< N << ", H="<< H << std::endl;
  
}
