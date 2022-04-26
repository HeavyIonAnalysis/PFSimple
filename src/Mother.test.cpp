#ifndef PFSIMPLE_MOTHER_TEST_HPP_
#define PFSIMPLE_MOTHER_TEST_HPP_

#include <gtest/gtest.h>

#include "Mother.hpp"

namespace {

TEST(Mother, Basics) {
  
  Mother mother_1(3122);
  mother_1.SetCutDistance(1.f);
  mother_1.SetCutDistanceToSV(2.f);
  mother_1.SetCutChi2Geo(3.f);
  mother_1.SetCutChi2GeoSM({4.f, 5.f, 6.f});
  mother_1.SetCutLdL(7.f);
  mother_1.SetCutDecayLength(8.f);
  mother_1.SetCutDistancePVLine(9.f);
  mother_1.SetCutChi2Topo(10.f);
  mother_1.SetCutChi2TopoSM({11.f, 12.f, 13.f});
  mother_1.SetCutCosTopo(14.f);
  mother_1.SetCutCosTopoSM({15.f, 16.f, 17.f});
  
  EXPECT_EQ(mother_1.GetPdg(), 3122);
  EXPECT_FLOAT_EQ(mother_1.GetCutDistance(), 1.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutDistanceToSV(), 2.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Geo().at(0), 3.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Geo().at(1), 4.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Geo().at(2), 5.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Geo().at(3), 6.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutLdL(), 7.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutDecayLength(), 8.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutDistancePVLine(), 9.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Topo().at(0), 10.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Topo().at(1), 11.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Topo().at(2), 12.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutChi2Topo().at(3), 13.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutCosTopo().at(0), 14.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutCosTopo().at(1), 15.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutCosTopo().at(2), 16.f);
  EXPECT_FLOAT_EQ(mother_1.GetCutCosTopo().at(3), 17.f);
  
  
  Mother mother_2(310);
  mother_2.CancelCuts();
  
  EXPECT_GT(mother_2.GetCutDistance(), 1e9);
  EXPECT_GT(mother_2.GetCutDistanceToSV(), 1e9);
  EXPECT_GT(mother_2.GetCutChi2Geo().at(0), 1e9);
  EXPECT_GT(mother_2.GetCutChi2Geo().at(1), 1e9);
  EXPECT_GT(mother_2.GetCutChi2Geo().at(2), 1e9);
  EXPECT_GT(mother_2.GetCutChi2Geo().at(3), 1e9);
  EXPECT_LT(mother_2.GetCutLdL(), -1e9);
  EXPECT_LT(mother_2.GetCutDecayLength(), -1e9);
  EXPECT_LT(mother_2.GetCutDistancePVLine(), -1e9);
  EXPECT_GT(mother_2.GetCutChi2Topo().at(0), 1e9);
  EXPECT_LT(mother_2.GetCutChi2Topo().at(1), -1e9);
  EXPECT_LT(mother_2.GetCutChi2Topo().at(2), -1e9);
  EXPECT_LT(mother_2.GetCutChi2Topo().at(3), -1e9);
  EXPECT_LT(mother_2.GetCutCosTopo().at(0), -1e9);
  EXPECT_GT(mother_2.GetCutCosTopo().at(1), 1e9);
  EXPECT_GT(mother_2.GetCutCosTopo().at(2), 1e9);
  EXPECT_GT(mother_2.GetCutCosTopo().at(3), 1e9);
  
}

}
#endif // PFSIMPLE_MOTHER_TEST_HPP_