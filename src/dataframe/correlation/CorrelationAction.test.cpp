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
#include <ROOT/RDataFrame.hxx>
#include <gtest/gtest.h>

#include "AverageHelper.hpp"
#include "AxesConfiguration.hpp"
#include "CorrelationAction.hpp"
#include "CorrelationFunctions.hpp"
#include "ReSampleFunctor.hpp"

/**
 * Basic test of the correlation
 * input: 1D Qvector with first harmonic Q vector of (0., 0.70710678)
 *        Event axes with 2 bins
 *        no event-to-event fluctuation
 *        correlation function: c = Qy_1 + Qx_2
 * expectation: mean of all events should be 0.70710678
 *              number of events in the first class should be 91
 *              number of events in the second class should be 9
 */
TEST(CorrelationAction, BasicIntegratedQ) {
  ROOT::RDataFrame df(100);
  auto event_axes = Qn::MakeAxes(Qn::AxisD("event", 2, 0, 2));
  auto q = Qn::DataContainerQVector();
  q.At(0).ActivateHarmonic(1);
  q.At(0).InitializeHarmonics();
  auto df1 = df.DefineSlotEntry("event",
                                [](unsigned int, ULong64_t entry) {
                                  double a = 0;
                                  if (entry > 90) a = 1.;
                                  return a;
                                },
                                {})
                 .Define("q",
                         [q](double event) -> Qn::DataContainerQVector {
                           auto ret = q;
                           for (int i = 0; i < 100; ++i) {
                             ret.At(0).Add(3. / 4 * TMath::Pi(), 1.);
                             ret.At(0).Add(1. / 4 * TMath::Pi(), 1.);
                           }
                           ret.At(0).CheckQuality();
                           ret.At(0) =
                               ret.At(0).Normal(Qn::QVector::Normalization::M);
                           return ret;
                         },
                         {"event"});
  auto dfs = Qn::Correlation::Resample(df1, 10);

  auto f_test = [](const Qn::QVector &a) { return a.x(1) + a.y(1); };
  auto f_weight = [](const Qn::QVector &a) { return a.sumweights()*(a.sumweights()-1); };
  std::vector<Qn::DataContainerQVector> initialization{q};
  auto correlation = Qn::Correlation::MakeCorrelationAction(
      "test", f_test, f_weight, Qn::Correlation::UseWeights::Yes, {"q"}, event_axes, 0);
  auto avg = Qn::MakeAverageHelper(correlation);
  auto res = avg.SetInitializationWithInitializationObject(&initialization)
                 .BookMe(dfs);
  auto data_container = res->GetDataContainer();
  EXPECT_NEAR(data_container.At(0).GetStatistics().Mean(), 0.70710678, 1e-5);
  EXPECT_EQ(data_container.At(0).GetStatistics().N(), 91);
  EXPECT_EQ(data_container.At(1).GetStatistics().N(), 9);


  Qn::QVector p;
  Qn::QVector r;
  r.ActivateHarmonic(2);
  r.InitializeHarmonics();
  p.ActivateHarmonic(2);
  p.InitializeHarmonics();

  p.Add(0.7,1.);
  p.Add(0.8,1.);
  r.Add(0.5,1.);
  r.Add(0.6,1.);

  auto func = Qn::Correlation::TwoParticle::d2(2);
  auto result = func(p,r);
}
