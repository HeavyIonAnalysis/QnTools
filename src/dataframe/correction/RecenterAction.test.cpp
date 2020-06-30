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
#include "RecenterAction.hpp"

TEST(RecenterActionUnitTest, Constructor) { ASSERT_EQ(1, 1); }

TEST(RecenterActionUnitTest, RDataFrame) {
  ROOT::RDataFrame df(100);
  std::size_t phi_size = 4;
  Qn::DataContainerQVector qvec_proto;
  auto event_axis = Qn::AxisD("event", 2, 0, 2);
  for (auto &q : qvec_proto) {
    q.ActivateHarmonic(1);
    q.InitializeHarmonics();
  }
  auto df0 = df.DefineSlotEntry("event",
                                [](unsigned int, ULong64_t entry) {
                                  double a = 0;
                                  if (entry > 90) a = 1.;
                                  return a;
                                },
                                {});
  auto df1 = df0.Define("q",
                        [qvec_proto](double event) -> Qn::DataContainerQVector {
                          auto ret = qvec_proto;
                          for (int i = 0; i < 100; ++i) {
                            ret.At(0).Add(0, 1.);
                          }
                          ret.At(0).CheckQuality();
                          ret.At(0).Normal(Qn::QVector::Normalization::M);
                          return ret;
                        },
                        {"event"});
  auto axes = Qn::MakeAxes(event_axis);
  auto res = Qn::Correction::MakeRecenterAction("test", axes, "q");
  TTreeReader reader;
  auto rec = Qn::AverageHelper(res)
                 .SetInitializationWithInitializationObject(&qvec_proto)
                 .BookMe(df1);
  auto resv = rec.GetValue();
  auto dfc = resv.ApplyCorrection(df1);
  auto qs = dfc.Take<Qn::DataContainerQVector>("q").GetValue();
  auto qsc = dfc.Take<Qn::DataContainerQVector>(resv.GetName()).GetValue();
}