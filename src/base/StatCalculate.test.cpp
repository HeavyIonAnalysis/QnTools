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

#include "StatCalculate.hpp"
#include "StatCollect.hpp"

class StatCalculateUnitTest : public ::testing::Test {
 protected:
  void SetUp() override {
    using namespace Qn;
    std::normal_distribution<> dist{2., 0.5};
    std::poisson_distribution<> poisson{1};
    std::mt19937 gen{0};
    std::vector<size_t> samples(nsamples_);
    Qn::StatCollect stats1;
    Qn::StatCollect stats2;
    Qn::StatCollect statssum;
    Qn::StatCollect statsempty;
    InitStats(stats1);
    InitStats(stats2);
    InitStats(statssum);
    InitStats(statsempty);
    auto ev = 1000l;
    auto halfev = ev / 2;
    for (size_t i = 0; i < ev; ++i) {
      for (auto &sample : samples) {
        sample = poisson(gen);
      }
      auto val = dist(gen);
      statssum.Fill(val, 1., samples);
      if (i < halfev) {
        stats1.Fill(val, 1., samples);
      } else {
        stats2.Fill(val, 1., samples);
      }
    }
    stat1_     = Qn::StatCalculate(stats1);
    stat2_     = Qn::StatCalculate(stats2);
    statsum_   = Qn::StatCalculate(statssum);
    statempty_ = Qn::StatCalculate(statsempty);
  }
  void InitStats(Qn::StatCollect &stat) {
    stat.SetNumberOfSamples(nsamples_);
    stat.SetWeightType(Qn::StatCollect::WeightType::OBSERVABLE);
  }

  Qn::StatCalculate stat1_;
  Qn::StatCalculate stat2_;
  Qn::StatCalculate statempty_;
  Qn::StatCalculate statsum_;
  double bootstrap_precision = 0.05;
  int nsamples_ = 100;
};

TEST_F(StatCalculateUnitTest, MergingWithEmptyLHS) {
  using namespace Qn;
  auto merged = Qn::Merge(statempty_, stat2_);
  auto manual_merge = (statempty_.Mean()* statempty_.SumWeights()
      +stat2_.Mean()* stat2_.SumWeights())
      /(statempty_.SumWeights()+ stat2_.SumWeights());
  EXPECT_FLOAT_EQ(stat2_.Mean(), merged.Mean());
  EXPECT_FLOAT_EQ(stat2_.Mean(), manual_merge);
  EXPECT_NEAR(stat2_.Variance(), merged.Variance(), stat2_.Variance()*0.01);
  // bootstrap vs propagated errors
  auto prop_error = merged.StdDevOfMeanFromPropagation();
  auto bs_error = merged.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);}

TEST_F(StatCalculateUnitTest, StatBinMergingBothEmpty) {
  using namespace Qn;
  auto merged = Qn::Merge(statempty_, statempty_);
  EXPECT_FLOAT_EQ(0., merged.Mean());
  EXPECT_FLOAT_EQ(0., merged.Variance());
  // bootstrap vs propagated errors
  auto prop_error = merged.StdDevOfMeanFromPropagation();
  auto bs_error = merged.StdDevOfMeanFromBootstrapVariance();
  EXPECT_FLOAT_EQ(prop_error, 0.);
  EXPECT_FLOAT_EQ(bs_error, 0.);
}

TEST_F(StatCalculateUnitTest, MergingWithEmptyRHS) {
  using namespace Qn;
  auto merged = Qn::Merge(stat2_, statempty_);
  auto manual_merge = (statempty_.Mean()* stat1_.SumWeights()
      +stat2_.Mean()* stat2_.SumWeights())
      /(statempty_.SumWeights()+ stat2_.SumWeights());
  EXPECT_FLOAT_EQ(stat2_.Mean(), merged.Mean());
  EXPECT_FLOAT_EQ(stat2_.Mean(), manual_merge);
  EXPECT_NEAR(stat2_.Variance(), merged.Variance(), stat2_.Variance()*0.01);
  // bootstrap vs propagated errors
  auto prop_error = merged.StdDevOfMeanFromPropagation();
  auto bs_error = merged.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(StatCalculateUnitTest, Merging) {
  using namespace Qn;
  auto merged = Qn::Merge(stat1_, stat2_);
  auto manual_merge = (stat1_.Mean()* stat1_.SumWeights()
      +stat2_.Mean()* stat2_.SumWeights())
      /(stat1_.SumWeights()+ stat2_.SumWeights());
  EXPECT_FLOAT_EQ(statsum_.Mean(), merged.Mean());
  EXPECT_FLOAT_EQ(statsum_.Mean(), manual_merge);
  EXPECT_NEAR(statsum_.Variance(), merged.Variance(),statsum_.Variance()*0.01);
  // bootstrap vs propagated errors
  auto prop_error = merged.StdDevOfMeanFromPropagation();
  auto bs_error = merged.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(StatCalculateUnitTest, Addition) {
  using namespace Qn;
  auto sum = stat1_ + stat2_;
  auto man_sum = stat1_.Mean() + stat2_.Mean();
  auto man_sum_v = stat1_.Variance() + stat2_.Variance();
  EXPECT_FLOAT_EQ(sum.Mean(), man_sum);
  EXPECT_FLOAT_EQ(sum.Variance(), man_sum_v);
  // bootstrap vs propagated errors
  auto prop_error = sum.StdDevOfMeanFromPropagation();
  auto bs_error = sum.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);}

TEST_F(StatCalculateUnitTest, Subtraction) {
  auto diff = stat1_ - stat2_;
  auto man_diff = stat1_.Mean() - stat2_.Mean();
  auto man_diff_v = stat1_.Variance() + stat2_.Variance();
  EXPECT_FLOAT_EQ(diff.Mean(), man_diff);
  EXPECT_FLOAT_EQ(diff.Variance(), man_diff_v);
  // bootstrap vs propagated errors
  auto prop_error = diff.StdDevOfMeanFromPropagation();
  auto bs_error = diff.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(StatCalculateUnitTest, Multiplication) {
  auto mult = stat1_ * stat2_;
  auto man_mult = stat1_.Mean() * stat2_.Mean();
  auto man_mult_v =   stat2_.Mean() * stat2_.Mean() * stat1_.Variance()
      + stat1_.Mean() * stat1_.Mean() * stat2_.Variance();
  EXPECT_FLOAT_EQ(mult.Mean(), man_mult);
  EXPECT_FLOAT_EQ(mult.Variance(), man_mult_v);
  // bootstrap vs propagated errors
  auto prop_error = mult.StdDevOfMeanFromPropagation();
  auto bs_error = mult.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(StatCalculateUnitTest, Division) {
  auto div = stat1_ / stat2_;
  auto man_div = stat1_.Mean() / stat2_.Mean();
  auto man_div_v =   stat1_.Variance() / (stat2_.Mean() * stat2_.Mean())
      +   stat1_.Mean() * stat1_.Mean() * stat2_.Variance()
          / (stat2_.Mean() * stat2_.Mean() * stat2_.Mean() * stat2_.Mean());
  EXPECT_FLOAT_EQ(div.Mean(), man_div);
  EXPECT_FLOAT_EQ(div.Variance(), man_div_v);
  // bootstrap vs propagated errors
  auto prop_error = div.StdDevOfMeanFromPropagation();
  auto bs_error = div.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);}

TEST_F(StatCalculateUnitTest, Pow) {
  auto pow = Qn::Pow(stat1_, 2.);
  auto man_pow = std::pow(stat1_.Mean(), 2.);
  auto man_pow_v = 4 * stat1_.Mean() * stat1_.Mean() * stat1_.Variance();
  EXPECT_FLOAT_EQ(pow.Mean(), man_pow);
  EXPECT_FLOAT_EQ(pow.Variance(), man_pow_v);
  // bootstrap vs propagated errors
  auto prop_error = pow.StdDevOfMeanFromPropagation();
  auto bs_error = pow.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);}

TEST_F(StatCalculateUnitTest, Scaling) {
  auto scale1 = statsum_ * (1./2);
  auto scale2 = (1./2) * statsum_;
  auto scale3 = statsum_ / 2;
  auto man_scale = statsum_.Mean() / 2;
  auto man_scale_v =   statsum_.Variance() / 4.;
  EXPECT_FLOAT_EQ(scale1.Mean(), scale2.Mean());
  EXPECT_FLOAT_EQ(scale1.Mean(), scale3.Mean());
  EXPECT_FLOAT_EQ(scale1.Mean(), man_scale);
  EXPECT_FLOAT_EQ(scale1.Variance(), man_scale_v);
  auto prop_error = scale1.StdDevOfMeanFromPropagation();
  auto bs_error = scale1.StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(StatCalculateUnitTest, ConvergenceOfSampleVariance) {
  using namespace Qn;
  Qn::StatCollect stats_largest;
  std::normal_distribution<> dist{2., 0.5};
  std::mt19937 gen{0};
  auto nsamples = 100;
  std::vector<size_t> samples(nsamples);
  std::poisson_distribution<> poisson{1};
  stats_largest.SetNumberOfSamples(nsamples);
  stats_largest.SetWeightType(Stat::WeightType::OBSERVABLE);
  auto stats_small = stats_largest;
  auto stats_larger = stats_largest;
  for (size_t i = 0; i < 1000l; ++i) {
    for (auto &sample : samples) {
      sample = poisson(gen);
    }
    auto val = dist(gen);
    stats_largest.Fill(val, 1., samples);
    if (i < 10) stats_small.Fill(val, 1., samples);
    if (i < 100) stats_larger.Fill(val, 1., samples);
  }
  auto statbin_largest = Qn::StatCalculate(stats_largest);
  auto statbin_larger = Qn::StatCalculate(stats_larger);
  auto statbin_small = Qn::StatCalculate(stats_small);

  EXPECT_TRUE(statbin_largest.StdDevOfMeanFromPropagation() < statbin_larger.StdDevOfMeanFromPropagation()
                  && statbin_larger.StdDevOfMeanFromPropagation() < statbin_small.StdDevOfMeanFromPropagation());
  EXPECT_TRUE(statbin_largest.StdDevOfMeanFromBootstrapVariance() < statbin_larger.StdDevOfMeanFromBootstrapVariance()
                  && statbin_larger.StdDevOfMeanFromBootstrapVariance() < statbin_small.StdDevOfMeanFromBootstrapVariance());
  EXPECT_TRUE(statbin_largest.StdDevOfMeanFromBootstrapPercentile() < statbin_larger.StdDevOfMeanFromBootstrapPercentile()
                  && statbin_larger.StdDevOfMeanFromBootstrapPercentile() < statbin_small.StdDevOfMeanFromBootstrapPercentile());
}

TEST_F(StatCalculateUnitTest, ConvergenceOfBootstrap) {
  using namespace Qn;
  Qn::StatCollect stats_largest;
  Qn::StatCollect stats_larger;
  Qn::StatCollect stats_small;
  std::normal_distribution<> dist{2., 0.5};
  std::mt19937 gen{0};
  auto nsamples_largest = 500;
  auto nsamples_larger = 50;
  auto nsamples_small = 5;
  std::vector<size_t> samples_largest(nsamples_largest);
  std::vector<size_t> samples_larger(nsamples_larger);
  std::vector<size_t> samples_small(nsamples_small);
  std::poisson_distribution<> poisson{1};
  stats_largest.SetNumberOfSamples(nsamples_largest);
  stats_larger.SetNumberOfSamples(nsamples_larger);
  stats_small.SetNumberOfSamples(nsamples_small);
  for (size_t i = 0; i < 1000l; ++i) {
    for (int is = 0; is < nsamples_largest; ++is) {
      auto sample  = poisson(gen);
      samples_largest[is] = sample;
      if (is < nsamples_larger) samples_larger[is] = sample;
      if (is < nsamples_small) samples_small[is] = sample;
    }
    auto val = dist(gen);
    stats_largest.Fill(val, 1., samples_largest);
    stats_larger.Fill(val, 1., samples_larger);
    stats_small.Fill(val, 1., samples_small);
  }
  auto statbin_largest = Qn::StatCalculate(stats_largest);
  auto statbin_larger = Qn::StatCalculate(stats_larger);
  auto statbin_small = Qn::StatCalculate(stats_small);

  auto diff_small   = std::abs(statbin_largest.StdDevOfMeanFromPropagation() - statbin_largest.StdDevOfMeanFromBootstrapVariance());
  auto diff_larger  = std::abs(statbin_larger.StdDevOfMeanFromPropagation() - statbin_larger.StdDevOfMeanFromBootstrapVariance());
  auto diff_largest = std::abs(statbin_small.StdDevOfMeanFromPropagation() - statbin_small.StdDevOfMeanFromBootstrapVariance());

  EXPECT_TRUE(diff_small < diff_larger && diff_larger < diff_largest);
}
