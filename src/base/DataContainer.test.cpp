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

#include <gtest/gtest.h>
#include "DataContainer.hpp"
#include "Stats.hpp"
#include "Axis.hpp"

TEST(ContainerUnitTest, Projection) {
  std::mt19937 gen{0};
  std::normal_distribution<> dist{2.,1.};
  Qn::DataContainerStats test2dobs{{Qn::AxisD("a1",10,0,10), Qn::AxisD("a2",10,0,10)}};
  Qn::DataContainerStats test2dref{{Qn::AxisD("a1",10,0,10), Qn::AxisD("a2",10,0,10)}};
  Qn::DataContainerStats test1d{{Qn::AxisD("a1",10,0,10)}};
  Qn::DataContainerStats test0d;
  for (auto &bin : test2dobs) { bin.SetWeights(Qn::Stats::Weights::OBSERVABLE); bin.ResetBits(Qn::Stats::CORRELATEDERRORS); }
  for (auto &bin : test2dref) { bin.SetWeights(Qn::Stats::Weights::REFERENCE); bin.ResetBits(Qn::Stats::CORRELATEDERRORS); }
  for (auto &bin : test1d) { bin.SetWeights(Qn::Stats::Weights::REFERENCE); bin.ResetBits(Qn::Stats::CORRELATEDERRORS); }
  for (auto &bin : test0d) { bin.SetWeights(Qn::Stats::Weights::REFERENCE); bin.ResetBits(Qn::Stats::CORRELATEDERRORS); }
  for (unsigned long i = 0; i < 10; ++i) {
    for (unsigned long j = 0; j < 10; ++j) {
      for (int n =0; n < 10000; ++n) {
        auto rndm = dist(gen);
        test0d.At(0).Fill(rndm,1.0,std::vector<size_t>{});
        test1d.At(i).Fill(rndm,1.0,std::vector<size_t>{});
        test2dobs.At({i,j}).Fill(rndm,1.0,std::vector<size_t>{});
        test2dref.At({i,j}).Fill(rndm,1.0,std::vector<size_t>{});
      }
    }
  }
  EXPECT_NEAR(test0d.At(0).Mean(), test1d.Projection().At(0).Mean(),0.0001);
  EXPECT_NEAR(test0d.At(0).Mean(), test2dref.Projection().At(0).Mean(),0.0001);
  EXPECT_NEAR((test2dref * (1./test0d.At(0).Mean())).At(0).Mean(), 1.0, test2dref.At(0).MeanError());
  EXPECT_NEAR((test2dref / test1d).At(0).Mean(), 1.0, test2dref.At(0).MeanError());
  EXPECT_FALSE(std::isnan((test2dref / test1d).Projection().At(0).Mean()));
  EXPECT_NEAR((test2dobs / test1d).Projection().At(0).Mean(), 1.0, test2dobs.At(0).MeanError());
}