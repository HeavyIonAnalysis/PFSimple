#ifndef PFSIMPLE_DAUGHTER_TEST_HPP_
#define PFSIMPLE_DAUGHTER_TEST_HPP_

#include <gtest/gtest.h>

#include "Daughter.hpp"

namespace {

TEST(Daughter, Basics) {
  
  Daughter daughter_1(211);
  daughter_1.SetCutChi2Prim(15.f);
  daughter_1.SetCutCos(0.99f);
  
  EXPECT_EQ(daughter_1.GetPdgHypo(), 211);
  EXPECT_EQ(daughter_1.GetPids().size(), 1);
  EXPECT_EQ(daughter_1.GetPids().at(0), 211);
  EXPECT_FLOAT_EQ(daughter_1.GetCutChi2Prim(), 15.f);
  EXPECT_FLOAT_EQ(daughter_1.GetCutCos(), 0.99f);
  
  
  Daughter daughter_2(211, {2112, 1});
  daughter_2.CancelCuts();
  
  EXPECT_EQ(daughter_2.GetPdgHypo(), 211);
  EXPECT_EQ(daughter_2.GetPids().size(), 2);
  EXPECT_EQ(daughter_2.GetPids().at(0), 2112);
  EXPECT_EQ(daughter_2.GetPids().at(1), 1);
  EXPECT_LT(daughter_2.GetCutChi2Prim(), -1e9);
  EXPECT_LT(daughter_2.GetCutCos(), -1e9);

}

}
#endif // PFSIMPLE_DAUGHTER_TEST_HPP_