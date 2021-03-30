// Flow Vector Correction Framework
//
// Copyright (C) 2020  Lukas Kreis Ilya Selyuzhenkov
// Contact: l.kreis@gsi.de; ilya.selyuzhenkov@gmail.com
// For a full list of contributors please see docs/Credits
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <gtest/gtest.h>

#include "AxesConfiguration.hpp"
#include "QVectorFunctor.hpp"

/**
 * Check the the constructor and the definition of a Q-vector with 1 subevents.
 * input 4 phi angles (45 deg, 135 deg, 215 deg, 305 deg)
 *       4 weights (1, 1, 1, 1)
 *       4 etas    (25, 25, 75, 75)
 *       subevents are binned in eta with the bin edges (0, 50, 100)
 *       expectation of the first harmonic Q-vector in the only subevent:
 *       x_1^{ev_1} = 0, y_1^{ev_1} = 0
 *
 */
TEST(QVectorFunctor, BasicIntegrated) {
  auto q_action = Qn::MakeQVectorFunctor("test", {1, 2, 3, 4});
  int si = 4;
  ROOT::RVec<double> phis(si);
  ROOT::RVec<double> weights(si);
  ROOT::RVec<double> etas(si);
  double angle = 1.;
  for (int i = 0.; i < si; ++i) {
    phis[i] = angle / 8 * 2 * TMath::Pi();
    weights[i] = (double)1;
    if (i < 2)
      etas[i] = (double)25.;
    else
      etas[i] = (double)75.;
    angle = angle + 2.;
  }
  auto q = q_action(phis, weights);
  EXPECT_NEAR(q.At(0).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(0).y(1), 0., 1e-5);
}

/**
 * Check the the constructor and the definition of a Q-vector with 2 subevents.
 * input 4 phi angles (45 deg, 135 deg, 215 deg, 305 deg)
 *       4 weights (1, 1, 1, 1)
 *       4 etas    (25, 25, 75, 75)
 *       subevents are binned in eta with the bin edges (0, 50, 100)
 *       expectation of the first harmonic Q-vector in the first subevent:
 *       x_1^{ev_1} = 0, y_1^{ev_1} = sqrt(2)
 *       expectation of the the first harmonic Q-vector in the second subevent:
 *       x_1^{ev_2} = 0, y_1^{ev_2} = -sqrt(2)
 *
 */
TEST(QVectorFunctor, Basic1D) {
  auto axes = Qn::MakeAxes(Qn::AxisD("t1", 2, 0, 100));
  Qn::QVectorFunctor<decltype(axes), std::tuple<ROOT::RVec<double>>> q_action(
      "test", {1, 2, 3, 4}, axes);
  int si = 4;
  ROOT::RVec<double> phis(si);
  ROOT::RVec<double> weights(si);
  ROOT::RVec<double> etas(si);
  double angle = 1.;
  for (int i = 0.; i < si; ++i) {
    phis[i] = angle / 8 * 2 * TMath::Pi();
    weights[i] = (double)1;
    if (i < 2)
      etas[i] = (double)25.;
    else
      etas[i] = (double)75.;
    angle = angle + 2.;
  }
  auto q = q_action(phis, weights, etas);
  EXPECT_NEAR(q.At(0).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(0).y(1), std::sqrt(2), 1e-5);
  EXPECT_NEAR(q.At(1).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(1).y(1), -std::sqrt(2), 1e-5);
}

/**
 * Check the the constructor and the definition of a Q-vector with 2 subevents.
 * input 4 phi angles (45 deg, 135 deg, 215 deg, 305 deg)
 *       4 weights (1, 1, 1, 1)
 *       4 etas    (25, 25, 75, 75)
 *       subevents are binned in eta with the bin edges (0, 50, 100)
 *       expectation of the first harmonic Q-vector in the first subevent:
 *       x_1^{ev_1} = 0, y_1^{ev_1} = sqrt(2)
 *       expectation of the first harmonic Q-vector in the second subevent:
 *       x_1^{ev_2} = 0, y_1^{ev_2} = -sqrt(2)
 *       expectation of the first harmonic Q-vector in the third subevent:
 *       x_1^{ev_1} = 0, y_1^{ev_1} = sqrt(2)
 *       expectation of the first harmonic Q-vector in the fourth subevent:
 *       x_1^{ev_2} = 0, y_1^{ev_2} = -sqrt(2)
 *
 */
TEST(QVectorFunctor, Basic2D) {
  auto axes = Qn::MakeAxes(Qn::AxisD("eta", 2, 2, 4), Qn::AxisD("pt", 2, 0, 2));
  auto q_action = Qn::MakeQVectorFunctor("test", {1, 2, 3, 4}, axes);
  int si = 8;
  ROOT::RVec<double> phis(si);
  ROOT::RVec<double> weights(si);
  ROOT::RVec<double> etas(si);
  ROOT::RVec<double> pts(si);
  double angle = 1.;
  for (int i = 0.; i < si; ++i) {
    if (angle > 7) angle = 1;
    phis[i] = angle / 8 * 2 * TMath::Pi();
    weights[i] = (double)1;
    angle = angle + 2.;
    if (i == 0 || i == 1) {
      etas[i] = 2.;
      pts[i] = 0.;
    }
    if (i == 2 || i == 3) {
      etas[i] = 2.;
      pts[i] = 1.;
    }
    if (i == 4 || i == 5) {
      etas[i] = 3.;
      pts[i] = 0.;
    }
    if (i == 6 || i == 7) {
      etas[i] = 3.;
      pts[i] = 1.;
    }
  }
  auto q = q_action(phis, weights, etas, pts);
  EXPECT_NEAR(q.At(0).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(0).y(1), std::sqrt(2), 1e-5);
  EXPECT_NEAR(q.At(1).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(1).y(1), -std::sqrt(2), 1e-5);
  EXPECT_NEAR(q.At(2).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(2).y(1), std::sqrt(2), 1e-5);
  EXPECT_NEAR(q.At(3).x(1), 0., 1e-5);
  EXPECT_NEAR(q.At(3).y(1), -std::sqrt(2), 1e-5);
}

/**
 * Check the the constructor and the definition of a Q-vector with 2 subevents
 * using the RDataFrame. input 4 phi angles (45 deg, 135 deg, 215 deg, 305 deg)
 *       4 weights (1, 1, 1, 1)
 *       4 etas    (25, 25, 75, 75)
 *       subevents are binned in eta with the bin edges (0, 50, 100)
 *       expectation of the first harmonic Q-vector in the first subevent:
 *       x_1^{ev_1} = 0, y_1^{ev_1} = sqrt(2)
 *       expectation of the first harmonic Q-vector in the second subevent:
 *       x_1^{ev_2} = 0, y_1^{ev_2} = -sqrt(2)
 *       All events in the DataFrame should give the same result.
 */
TEST(QVectorFunctor, RDataFrame) {
  ROOT::RDataFrame df(8);
  std::size_t phi_size = 4;
  auto df1 = df.Define("phis",
                       [phi_size]() {
                         ROOT::RVec<double> phis(phi_size);
                         double angle = 1.;
                         for (int i = 0.; i < phi_size; ++i) {
                           phis[i] = angle / 8 * 2 * TMath::Pi();
                           angle = angle + 2.;
                         }
                         return phis;
                       },
                       {})
                 .Define("weights",
                         [phi_size]() {
                           ROOT::RVec<double> weights(phi_size);
                           for (int i = 0.; i < phi_size; ++i) {
                             weights[i] = 1;
                           }
                           return weights;
                         },
                         {})
                 .Define("etas",
                         [phi_size]() {
                           ROOT::RVec<double> etas(phi_size);
                           for (int i = 0.; i < phi_size; ++i) {
                             if (i < phi_size / 2) {
                               etas[i] = (double)25.;
                             } else {
                               etas[i] = (double)75.;
                             }
                           }
                           return etas;
                         },
                         {});
  auto axes = Qn::MakeAxes(Qn::AxisD("t1", 2, 0, 100));
  auto q_action = Qn::MakeQVectorFunctor("test", {1, 2, 3, 4}, axes);
  auto df2 =
      df1.Define(q_action.GetName(), q_action, {"phis", "weights", "etas"});
  auto qs = df2.Take<Qn::DataContainerQVector>("test");
  for (int i = 0; i < 8; ++i) {
    auto q = qs.GetValue().at(i);
    EXPECT_NEAR(q.At(0).x(1), 0., 1e-5);
    EXPECT_NEAR(q.At(0).y(1), std::sqrt(2), 1e-5);
    EXPECT_NEAR(q.At(1).x(1), 0., 1e-5);
    EXPECT_NEAR(q.At(1).y(1), -std::sqrt(2), 1e-5);
  }
}
