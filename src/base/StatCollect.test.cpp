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
#include "StatCollect.hpp"

class StatCollectUnitTest : public ::testing::Test {
 protected:
  void SetUp() override {
    using namespace Qn;
    std::normal_distribution<> dist{2., 0.5};
    std::poisson_distribution<> poisson{1};
    std::vector<size_t> samples(nsamples_);
    stat1_.SetNumberOfSamples(nsamples_);
    stat2_.SetNumberOfSamples(nsamples_);
    statsum_.SetNumberOfSamples(nsamples_);
    statempty_.SetNumberOfSamples(nsamples_);
    std::mt19937 gen{0};
    auto halfev = nevents_ / 2;
    for (size_t i = 0; i < nevents_; ++i) {
      for (auto &sample : samples) {
        sample = poisson(gen);
      }
      auto val = dist(gen);
      statsum_.Fill(val, 1., samples);
      if (i < halfev) {
        stat1_.Fill(val, 1., samples);
      } else {
        stat2_.Fill(val, 1., samples);
      }
    }
  }

  Qn::StatCollect stat1_;
  Qn::StatCollect stat2_;
  Qn::StatCollect statempty_;
  Qn::StatCollect statsum_;
  long nevents_ = 1000;
  int nsamples_ = 100;
};

TEST_F(StatCollectUnitTest, Filling) {
  using namespace Qn;
  auto statistics = statsum_.GetStatistics();
  auto bootstrap = statsum_.GetBootStrap();
  EXPECT_NEAR(statistics.Mean(), 2., 0.5/std::sqrt(nevents_));
  EXPECT_EQ(statistics.N(), nevents_);
  EXPECT_NEAR(statistics.StandardErrorOfMean(), 0.5/std::sqrt(nevents_), 0.5/std::sqrt(nevents_)*0.1);
  EXPECT_NEAR(statistics.Variance(), 0.25, 0.25*0.1);
  auto means = bootstrap.GetMeans();
  auto weights = bootstrap.GetWeights();
  Statistics bs_stat;
  for (int i = 0; i < means.size(); ++i) { bs_stat.Fill(means[i], weights[i]); }
  EXPECT_NEAR(std::sqrt(bs_stat.Variance()), 0.5/std::sqrt(nevents_), 0.5/std::sqrt(nevents_)*0.10);
}

TEST_F(StatCollectUnitTest, StatBinMerging) {
  using namespace Qn;
  auto merged = Qn::Merge(stat1_, stat2_);
  auto manual_merge = (stat1_.GetStatistics().Mean()* stat1_.GetStatistics().SumWeights()
      +stat2_.GetStatistics().Mean()* stat2_.GetStatistics().SumWeights())
      /(stat1_.GetStatistics().SumWeights()+ stat2_.GetStatistics().SumWeights());
  EXPECT_FLOAT_EQ(statsum_.GetStatistics().Mean(), merged.GetStatistics().Mean());
  EXPECT_FLOAT_EQ(statsum_.GetStatistics().Mean(), manual_merge);
  EXPECT_FLOAT_EQ(statsum_.GetStatistics().Variance(), merged.GetStatistics().Variance());
}

TEST_F(StatCollectUnitTest, StatBinMergingWithEmptyLHS) {
  using namespace Qn;
  auto merged = Qn::Merge(statempty_, stat2_);
  auto manual_merge = (statempty_.GetStatistics().Mean()* statempty_.GetStatistics().SumWeights()
      +stat2_.GetStatistics().Mean()* stat2_.GetStatistics().SumWeights())
      /(statempty_.GetStatistics().SumWeights()+ stat2_.GetStatistics().SumWeights());
  EXPECT_FLOAT_EQ(stat2_.GetStatistics().Mean(), merged.GetStatistics().Mean());
  EXPECT_FLOAT_EQ(stat2_.GetStatistics().Mean(), manual_merge);
  EXPECT_FLOAT_EQ(stat2_.GetStatistics().Variance(), merged.GetStatistics().Variance());
}

TEST_F(StatCollectUnitTest, StatBinMergingWithEmptyRHS) {
  using namespace Qn;
  auto merged = Qn::Merge(stat2_, statempty_);
  auto manual_merge = (statempty_.GetStatistics().Mean()* statempty_.GetStatistics().SumWeights()
      +stat2_.GetStatistics().Mean()* stat2_.GetStatistics().SumWeights())
      /(statempty_.GetStatistics().SumWeights()+ stat2_.GetStatistics().SumWeights());
  EXPECT_FLOAT_EQ(stat2_.GetStatistics().Mean(), merged.GetStatistics().Mean());
  EXPECT_FLOAT_EQ(stat2_.GetStatistics().Mean(), manual_merge);
  EXPECT_FLOAT_EQ(stat2_.GetStatistics().Variance(), merged.GetStatistics().Variance());
}

TEST_F(StatCollectUnitTest, StatBinMergingWithBothEmpty) {
  using namespace Qn;
  auto merged = Qn::Merge(statempty_, statempty_);
  EXPECT_FLOAT_EQ(0., merged.GetStatistics().Mean());
  EXPECT_FLOAT_EQ(0., merged.GetStatistics().Variance());
}

