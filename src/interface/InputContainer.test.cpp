#ifndef PFSIMPLE_INPUTCONTAINER_TEST_HPP_
#define PFSIMPLE_INPUTCONTAINER_TEST_HPP_

#include <gtest/gtest.h>

#include "InputContainer.hpp"

#include "Constants.hpp"

namespace {

TEST(InputContainer, Basics) {

  InputContainer inputcontainer = InputContainer();

  EXPECT_EQ(inputcontainer.GetTracks().size(), 0);

  const int n_tracks = 10;
  inputcontainer.Reserve(n_tracks);

  EXPECT_EQ(inputcontainer.GetTracks().size(), 0);

  inputcontainer.SetPV(1.f, 2.f, 3.f);

  const float mf_value = 2.f;
  const float cov_value = 3.5f;
  const float par_value = 0.5f;
  const int q_value = 3;
  const int id_value = 21;
  const int pdg_value = 11;

  for (int i_track = 0; i_track < n_tracks; i_track++) {

    std::vector<float> mf(NumberOfFieldPars, 0.f);
    for (int i = 0; i < NumberOfFieldPars; i++)
      mf.at(i) = mf_value * i_track * i;

    std::vector<float> cov(NumberOfCovElements, 0.f);
    for (int i = 0; i < NumberOfCovElements; i++)
      cov.at(i) = cov_value * i_track * i;

    std::vector<float> par(kNumberOfTrackPars, 0.f);
    for (int i = 0; i < kNumberOfTrackPars; i++)
      par.at(i) = par_value * i_track * i;

    const int q = q_value * i_track;
    const int id = id_value * i_track;
    const int pdg = pdg_value * i_track;

    inputcontainer.AddTrack(par, cov, mf, q, pdg, id);
  }

  EXPECT_EQ(inputcontainer.GetTracks().size(), n_tracks);
  EXPECT_FLOAT_EQ(inputcontainer.GetVertex().X(), 1.f);
  EXPECT_FLOAT_EQ(inputcontainer.GetVertex().Y(), 2.f);
  EXPECT_FLOAT_EQ(inputcontainer.GetVertex().Z(), 3.f);
  EXPECT_NEAR(inputcontainer.GetTracks().at(0).X(), par_value * 0 * 0, 0.01);
  EXPECT_NEAR(inputcontainer.GetTracks().at(1).Z(), par_value * 1 * 2, 0.01);
  EXPECT_NEAR(inputcontainer.GetTracks().at(2).GetP(), std::sqrt(par_value * 2 * 3 * par_value * 2 * 3 + par_value * 2 * 4 * par_value * 2 * 4 + par_value * 2 * 5 * par_value * 2 * 5), 0.01);
  EXPECT_NEAR(inputcontainer.GetTracks().at(3).GetCovariance(6), cov_value * 3 * 6, 0.01);
  EXPECT_NEAR(inputcontainer.GetTracks().at(4).GetFieldCoeff()[5], mf_value * 4 * 5, 0.01);
  EXPECT_EQ(inputcontainer.GetTracks().at(5).Q(), q_value * 5);
  EXPECT_EQ(inputcontainer.GetTracks().at(6).Id(), id_value * 6);
  EXPECT_EQ(inputcontainer.GetTracks().at(7).GetPDG(), pdg_value * 7);

  inputcontainer.Clear();
  EXPECT_EQ(inputcontainer.GetTracks().size(), 0);
}

}// namespace
#endif// PFSIMPLE_INPUTCONTAINER_TEST_HPP_