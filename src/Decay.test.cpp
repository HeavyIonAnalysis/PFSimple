#ifndef PFSIMPLE_DECAY_TEST_HPP_
#define PFSIMPLE_DECAY_TEST_HPP_

#include <gtest/gtest.h>

#include "Decay.hpp"

namespace {

TEST(Decay, Basics) {
  
  Mother lambda(3122);
  Daughter proton(2212);
  Daughter pion(-211);
  Decay decay("lambda", lambda, {pion, proton});
  
  EXPECT_STREQ(decay.GetName().c_str(), "lambda");
  EXPECT_EQ(decay.GetMother().GetPdg(), 3122);
  EXPECT_EQ(decay.GetDaughters().size(), 2);
  EXPECT_EQ(decay.GetDaughters().at(0).GetPdgHypo(), -211);
  EXPECT_EQ(decay.GetDaughters().at(1).GetPdgHypo(), 2212);
  EXPECT_EQ(decay.GetDaughters().at(0).GetId(), 0);
  EXPECT_EQ(decay.GetDaughters().at(1).GetId(), 1);
  EXPECT_THROW(decay.SetDaughters({pion, proton}), std::runtime_error);
  
}

}
#endif // PFSIMPLE_DECAY_TEST_HPP_