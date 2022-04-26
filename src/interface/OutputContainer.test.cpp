#ifndef PFSIMPLE_OUTPUTCONTAINER_TEST_HPP_
#define PFSIMPLE_OUTPUTCONTAINER_TEST_HPP_

#include <gtest/gtest.h>

#include "OutputContainer.hpp"

namespace {

TEST(OutputContainer, Basics) {
  
  KFParticle particle;
  particle.AddDaughterId(27);
  particle.AddDaughterId(39);
  particle.X() = 1.f;
  particle.Y() = 2.f;
  particle.Z() = 3.f;
  particle.Px() = 4.f;
  particle.Py() = 5.f;
  particle.Pz() = 6.f;
  particle.E() = 10.f;
  particle.SetPDG(3122);
  
  SelectionValues values;
  values.chi2_prim = {1.f, 2.f, 3.f};
  values.cos = {4.f, 5.f, 6.f};
  values.distance = 7.f;    
  values.distance_sv = 8.f; 
  values.l = 9.f;           
  values.l_over_dl = 10.f;   
  values.distance_pv = 11.f;
  values.chi2_geo = {12.f, 13.f, 14.f, 15.f};
  values.chi2_topo = {16.f, 17.f, 18.f, 19.f};
  values.cos_topo = {20.f, 21.f, 22.f, 23.f};
  values.is_from_PV = true; 

  OutputContainer outputcontainer = OutputContainer(particle);
  outputcontainer.SetSelectionValues(values);
  
  EXPECT_EQ(outputcontainer.GetDaughterIds().size(), 2);
  EXPECT_EQ(outputcontainer.GetDaughterIds().at(0), 27);
  EXPECT_EQ(outputcontainer.GetDaughterIds().at(1), 39);
  EXPECT_FLOAT_EQ(outputcontainer.GetX(), 1.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetY(), 2.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetZ(), 3.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetPx(), 4.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetPy(), 5.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetPz(), 6.f);
  EXPECT_NEAR(outputcontainer.GetMass(), std::sqrt(10*10 - 4*4 - 5*5 - 6*6), 0.001);
  EXPECT_EQ(outputcontainer.GetPdg(), 3122);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Prim(0), 1.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Prim(1), 2.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Prim(2), 3.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCos(0), 4.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCos(1), 5.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCos(2), 6.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetDistance(), 7.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetDistanceToSV(), 8.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetL(), 9.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetLdL(), 10.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetDistanceToPVLine(), 11.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Geo(0), 12.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Geo(1), 13.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Geo(2), 14.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Geo(3), 15.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Topo(0), 16.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Topo(1), 17.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Topo(2), 18.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetChi2Topo(3), 19.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCosineTopo(0), 20.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCosineTopo(1), 21.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCosineTopo(2), 22.f);
  EXPECT_FLOAT_EQ(outputcontainer.GetCosineTopo(3), 23.f);
}

}
#endif // PFSIMPLE_OUTPUTCONTAINER_TEST_HPP_