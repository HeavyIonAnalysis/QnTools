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

#include <random>

#include <gtest/gtest.h>

#include "Statistics.hpp"

class StatisticsUnitTest : public ::testing::Test {
 protected:
  void SetUp() override {
    using namespace Qn;
    std::normal_distribution<> dist{2., 0.5};
    std::mt19937 gen{0};
    auto halfev = nevents_ / 2;
    for (size_t i = 0; i < nevents_; ++i) {
      auto val = dist(gen);
      statsum_.Fill(val, 1.);
      if (i < halfev) {
        stat1_.Fill(val, 1.);
      } else {
        stat2_.Fill(val, 1.);
      }
    }
  }

  Qn::Statistics stat1_;
  Qn::Statistics stat2_;
  Qn::Statistics statempty_;
  Qn::Statistics statsum_;
  long nevents_ = 1000;
  int nsamples_ = 100;
};

TEST_F(StatisticsUnitTest, Filling) {
  using namespace Qn;
  EXPECT_NEAR(statsum_.Mean(), 2., 0.5/std::sqrt(nevents_));
  EXPECT_EQ(statsum_.N(), nevents_);
  EXPECT_NEAR(statsum_.StandardErrorOfMean(), 0.5/std::sqrt(nevents_), 0.5/std::sqrt(nevents_)*0.1);
  EXPECT_NEAR(statsum_.Variance(), 0.25, 0.25*0.1);
}

TEST_F(StatisticsUnitTest, StatBinMerging) {
  using namespace Qn;
  auto merged = Qn::Merge(stat1_, stat2_);
  auto manual_merge = (stat1_.Mean()* stat1_.SumWeights()
      +stat2_.Mean()* stat2_.SumWeights())
      /(stat1_.SumWeights()+ stat2_.SumWeights());
  EXPECT_FLOAT_EQ(statsum_.Mean(), merged.Mean());
  EXPECT_FLOAT_EQ(statsum_.Mean(), manual_merge);
  EXPECT_FLOAT_EQ(statsum_.Variance(), merged.Variance());
}

TEST_F(StatisticsUnitTest, StatBinMergingWithEmptyLHS) {
  using namespace Qn;
  auto merged = Qn::Merge(statempty_, stat2_);
  auto manual_merge = (statempty_.Mean()* statempty_.SumWeights()
      +stat2_.Mean()* stat2_.SumWeights())
      /(statempty_.SumWeights()+ stat2_.SumWeights());
  EXPECT_FLOAT_EQ(stat2_.Mean(), merged.Mean());
  EXPECT_FLOAT_EQ(stat2_.Mean(), manual_merge);
  EXPECT_FLOAT_EQ(stat2_.Variance(), merged.Variance());
}

TEST_F(StatisticsUnitTest, StatBinMergingWithEmptyRHS) {
  using namespace Qn;
  auto merged = Qn::Merge(stat2_, statempty_);
  auto manual_merge = (statempty_.Mean()* statempty_.SumWeights()
      +stat2_.Mean()* stat2_.SumWeights())
      /(statempty_.SumWeights()+ stat2_.SumWeights());
  EXPECT_FLOAT_EQ(stat2_.Mean(), merged.Mean());
  EXPECT_FLOAT_EQ(stat2_.Mean(), manual_merge);
  EXPECT_FLOAT_EQ(stat2_.Variance(), merged.Variance());
}

TEST_F(StatisticsUnitTest, StatBinMergingWithBothEmpty) {
  using namespace Qn;
  auto merged = Qn::Merge(statempty_, statempty_);
  EXPECT_FLOAT_EQ(0., merged.Mean());
  EXPECT_FLOAT_EQ(0., merged.Variance());
}

