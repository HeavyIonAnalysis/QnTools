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
#include "QVectorNormalizationFunctor.hpp"

/**
 * Check the the constructor and the definition of a Q-vector with 2 subevents
 * using the RDataFrame. input 4 phi angles (45 deg, 135 deg, 215 deg, 305 deg)
 * All qvectors are normalized by the number of contributors (2).
 *       4 weights (1, 1, 1, 1)
 *       4 etas    (25, 25, 75, 75)
 *       subevents are binned in eta with the bin edges (0, 50, 100)
 *       expectation of the first harmonic Q-vector in the first subevent:
 *       x_1^{ev_1} = 0, y_1^{ev_1} = sqrt(2)/2
 *       expectation of the first harmonic Q-vector in the second subevent:
 *       x_1^{ev_2} = 0, y_1^{ev_2} = -sqrt(2)/2
 *       All events in the DataFrame should give the same result.
 */
TEST(QVectorNormalizationFunctor, RDataFrame) {
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
auto q_normalization = Qn::QVectorNormalizationFunctor(Qn::QVectorNormalizationFunctor::Normalization::M);
auto df2 =
    df1.Define(q_action.GetName(), q_action, {"phis", "weights", "etas"})
       .Define(q_normalization.GetName(q_action.GetName()), q_normalization, {q_action.GetName()});
auto qs = df2.Take<Qn::DataContainerQVector>("test");
auto qsnorm = df2.Take<Qn::DataContainerQVector>("test_norm_M");

for (int i = 0; i < 8; ++i) {
auto q = qs.GetValue().at(i);
auto qnorm = qsnorm.GetValue().at(i);
EXPECT_NEAR(qnorm.At(0).x(1), 0., 1e-5);
EXPECT_NEAR(qnorm.At(0).y(1), std::sqrt(2)/2, 1e-5);
EXPECT_NEAR(qnorm.At(1).x(1), 0., 1e-5);
EXPECT_NEAR(qnorm.At(1).y(1), -std::sqrt(2)/2, 1e-5);
}
}
