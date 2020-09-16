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

#include "TCanvas.h"
#include "TFile.h"

#include "Axis.hpp"
#include "DataContainer.hpp"
#include "StatCollect.hpp"

class DataContainerUnitTest : public ::testing::Test {
 protected:
  void SetUp() override {
    using namespace Qn;
    stats_all.AddAxes({AxisD("Centrality", 4, 0, 100)});
    stats_nd.AddAxes({AxisD("a1", 2, 0, 100),AxisD("a2", 2, 0, 100)});
    std::normal_distribution<> dist{2., 0.5};
    std::poisson_distribution<> poisson{1};
    std::mt19937 gen{0};
    size_t nsamples = 100;
    std::vector<size_t> samples(nsamples);
    for (auto & stat : stats_all) {
      stat.SetNumberOfSamples(nsamples);
      stat.SetWeightType(Stat::WeightType::OBSERVABLE);
    }
    for (auto & stat : stats_nd) {
      stat.SetNumberOfSamples(nsamples);
      stat.SetWeightType(Stat::WeightType::OBSERVABLE);
    }
    stats_part_1 = stats_all;
    stats_part_2 = stats_all;
    auto ev = 10000l;
    auto halfev = ev / 2;
    for (size_t i = 0; i < ev; ++i) {
      for (auto &sample : samples) {
        sample = poisson(gen);
      }
      for (int ibin = 0; ibin < stats_all.size(); ++ibin) {
        auto val = dist(gen);
        stats_all[ibin].Fill(val, 1., samples);
        if (i < halfev) {
          stats_part_1[ibin].Fill(val, 1., samples);
        } else {
          stats_part_2[ibin].Fill(val, 1., samples);
        }
      }
      for (auto & bin : stats_nd) {
        auto val = dist(gen);
        bin.Fill(val, 1., samples);
      }
    }
    statbin_all = Qn::BinnedStatistics(stats_all);
    statbin_nd = Qn::BinnedStatistics(stats_nd);
    statbin_part_1 = Qn::BinnedStatistics(stats_part_1);
    statbin_part_2 = Qn::BinnedStatistics(stats_part_2);
  }

  double bootstrap_precision = 0.05;
  Qn::DataContainerStatCollect stats_nd;
  Qn::DataContainerStatCollect stats_all;
  Qn::DataContainerStatCollect stats_part_1;
  Qn::DataContainerStatCollect stats_part_2;
  Qn::BinnedStatistics statbin_nd;
  Qn::BinnedStatistics statbin_all;
  Qn::BinnedStatistics statbin_part_1;
  Qn::BinnedStatistics statbin_part_2;
};

TEST_F(DataContainerUnitTest, StatBinMerging) {
  using namespace Qn;
  Qn::BinnedStatistics merged_statbin(stats_part_1);
  auto list = new TList();
  list->Add(&statbin_part_2);
  auto manual_merge = (statbin_part_1[0].Mean()* statbin_part_1[0].SumWeights()
                      +statbin_part_2[0].Mean()* statbin_part_2[0].SumWeights())
                      /(statbin_part_1[0].SumWeights()+ statbin_part_2[0].SumWeights());
  merged_statbin.Merge(list);
  EXPECT_FLOAT_EQ(statbin_all[0].Mean(), merged_statbin[0].Mean());
  EXPECT_NE(statbin_all[0].Mean(), statbin_part_1[0].Mean());
  EXPECT_NE(statbin_all[0].Mean(), statbin_part_2[0].Mean());
  EXPECT_FLOAT_EQ(statbin_all[0].Mean(), manual_merge);
  EXPECT_NEAR(statbin_all[0].Variance(), merged_statbin[0].Variance(), statbin_all[0].Variance()*0.01);
  // bootstrap vs propagated errors
  auto prop_error = merged_statbin[0].StdDevOfMeanFromPropagation();
  auto bs_error = merged_statbin[0].StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(DataContainerUnitTest, StatBinAddition) {
  using namespace Qn;
  auto sum = statbin_part_1 + statbin_part_2;
  auto man_sum = statbin_part_1[0].Mean() + statbin_part_2[0].Mean();
  auto man_sum_v = statbin_part_1[0].Variance() + statbin_part_2[0].Variance();
  EXPECT_FLOAT_EQ(sum[0].Mean(), man_sum);
  EXPECT_FLOAT_EQ(sum[0].Variance(), man_sum_v);
  // bootstrap vs propagated errors
  auto prop_error = sum[0].StdDevOfMeanFromPropagation();
  auto bs_error = sum[0].StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(DataContainerUnitTest, StatBinSubtraction) {
  auto diff = statbin_part_1 - statbin_part_2;
  auto man_diff = statbin_part_1[0].Mean() - statbin_part_2[0].Mean();
  auto man_diff_v = statbin_part_1[0].Variance() + statbin_part_2[0].Variance();
  EXPECT_FLOAT_EQ(diff[0].Mean(), man_diff);
  EXPECT_FLOAT_EQ(diff[0].Variance(), man_diff_v);
  // bootstrap vs propagated errors
  auto prop_error = diff[0].StdDevOfMeanFromPropagation();
  auto bs_error = diff[0].StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(DataContainerUnitTest, StatBinMultiplication) {  
  auto mult = statbin_part_1 * statbin_part_2;
  auto man_mult = statbin_part_1[0].Mean() * statbin_part_2[0].Mean();
  auto man_mult_v =   statbin_part_2[0].Mean() * statbin_part_2[0].Mean() * statbin_part_1[0].Variance()
      + statbin_part_1[0].Mean() * statbin_part_1[0].Mean() * statbin_part_2[0].Variance();
  EXPECT_FLOAT_EQ(mult[0].Mean(), man_mult);
  EXPECT_FLOAT_EQ(mult[0].Variance(), man_mult_v);
  // bootstrap vs propagated errors
  auto prop_error = mult[0].StdDevOfMeanFromPropagation();
  auto bs_error = mult[0].StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(DataContainerUnitTest, StatBinDivision) {
  auto div = statbin_part_1 / statbin_part_2;
  auto man_div = statbin_part_1[0].Mean() / statbin_part_2[0].Mean();
  auto man_div_v =   statbin_part_1[0].Variance() / (statbin_part_2[0].Mean() * statbin_part_2[0].Mean())
      +   statbin_part_1[0].Mean() * statbin_part_1[0].Mean() * statbin_part_2[0].Variance()
          / (statbin_part_2[0].Mean() * statbin_part_2[0].Mean() * statbin_part_2[0].Mean() * statbin_part_2[0].Mean());
  EXPECT_FLOAT_EQ(div[0].Mean(), man_div);
  EXPECT_FLOAT_EQ(div[0].Variance(), man_div_v);
  // bootstrap vs propagated errors
  auto prop_error = div[0].StdDevOfMeanFromPropagation();
  auto bs_error = div[0].StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(DataContainerUnitTest, StatBinPow) {
  auto pow = Qn::Pow(statbin_part_1, 2.);
  auto man_pow = std::pow(statbin_part_1[0].Mean(), 2.);
  auto man_pow_v = 4 * statbin_part_1[0].Mean() * statbin_part_1[0].Mean() * statbin_part_1[0].Variance();
  EXPECT_FLOAT_EQ(pow[0].Mean(), man_pow);
  EXPECT_FLOAT_EQ(pow[0].Variance(), man_pow_v);
  // bootstrap vs propagated errors
  auto prop_error = pow[0].StdDevOfMeanFromPropagation();
  auto bs_error = pow[0].StdDevOfMeanFromBootstrapVariance();
  auto is_equal = std::abs(prop_error - bs_error) / prop_error < bootstrap_precision;
  auto is_larger = bs_error > prop_error;
  EXPECT_TRUE(is_equal || is_larger);
}

TEST_F(DataContainerUnitTest, StatBinScaling) {
  auto scale1 = statbin_part_1 * (1./2);
  auto scale2 = (1./2) * statbin_part_1;
  auto scale3 = statbin_part_1 / 2;
  auto man_scale = statbin_part_1[0].Mean() / 2;
  auto man_scale_v =   statbin_part_1[0].Variance() / 4.;
  EXPECT_FLOAT_EQ(scale1[0].Mean(), scale2[0].Mean());
  EXPECT_FLOAT_EQ(scale1[0].Mean(), scale3[0].Mean());
  EXPECT_FLOAT_EQ(scale1[0].Mean(), man_scale);
  EXPECT_FLOAT_EQ(scale1[0].Variance(), man_scale_v);
}

TEST_F(DataContainerUnitTest, RebinAfterDivision) {
  auto div = statbin_part_1 / statbin_part_2;
  auto rebin = div.Rebin(Qn::AxisD("Centrality", 2, 0, 100));
  auto div_rebin_m_0 = rebin[0].Mean();
  auto div_m_0 = div[0].Mean();
  auto div_m_1 = div[1].Mean();
  auto div_w_0 = div[0].SumWeights();
  auto div_w_1 = div[1].SumWeights();
  EXPECT_FLOAT_EQ(div_rebin_m_0, (div_m_0*div_w_0 + div_m_1*div_w_1)/(div_w_0+div_w_1));
}

TEST_F(DataContainerUnitTest, ConvergenceOfSampleVariance) {
  using namespace Qn;
  Qn::DataContainerStatCollect stats_largest;
  std::normal_distribution<> dist{2., 0.5};
  std::mt19937 gen{0};
  std::vector<size_t> samples(10, 1);
  for (auto & stat : stats_all) {
    stat.SetNumberOfSamples(10);
    stat.SetWeightType(Stat::WeightType::OBSERVABLE);
  }
  auto stats_small = stats_largest;
  auto stats_larger = stats_largest;
  for (size_t i = 0; i < 1000l; ++i) {
      auto val = dist(gen);
      stats_largest[0].Fill(val, 1., samples);
      if (i < 10) stats_small[0].Fill(val, 1., samples);
      if (i < 100) stats_larger[0].Fill(val, 1., samples);
  }
  auto statbin_largest = Qn::BinnedStatistics(stats_largest);
  auto statbin_larger = Qn::BinnedStatistics(stats_larger);
  auto statbin_small = Qn::BinnedStatistics(stats_small);

  EXPECT_TRUE(statbin_largest[0].StdDevOfMeanFromPropagation() < statbin_larger[0].StdDevOfMeanFromPropagation()
  && statbin_larger[0].StdDevOfMeanFromPropagation() < statbin_small[0].StdDevOfMeanFromPropagation());
}

TEST_F(DataContainerUnitTest, BootstrapErrorComparison) {
  auto test = statbin_all;
  auto prop_error = test[0].StdDevOfMeanFromPropagation();
  auto bs_error = test[0].StdDevOfMeanFromBootstrapVariance();
  auto bs_error_percentile = test[0].StdDevOfMeanFromBootstrapVariance();
  EXPECT_NEAR(prop_error, bs_error, prop_error*0.05);
  EXPECT_NEAR(bs_error_percentile, bs_error, prop_error*0.05);
}

TEST_F(DataContainerUnitTest, ErrorAfterMultiplicationAndRebin) {
  auto mult = statbin_all * statbin_all;
  auto rebin = mult.Rebin(Qn::AxisD("Centrality", 2, 0, 100));
  auto graph = Qn::ToTGraph(mult);
  auto graph_rebin = Qn::ToTGraph(rebin);
  rebin.SetErrors(Qn::StatCalculate::ErrorType::BOOTSTRAP);
  auto graph_rebin_bs = Qn::ToTGraph(rebin);
  auto canvas = new TCanvas("c1","c1",800, 600);
  graph->Draw("APL");
  graph_rebin->Draw("PL");
  graph_rebin->SetLineColor(kRed);
  canvas->SaveAs("ErrorsAfterRebin.png");
  auto file = new TFile("test.root", "RECREATE");
  rebin.Write("test");
  stats_all.Write("test2");
  file->Close();
}

TEST_F(DataContainerUnitTest, StatBinToTGraphConversion) {
  auto graph = Qn::ToTGraph(statbin_all);
  auto canvas = new TCanvas("c1","c1",800, 600);
  graph->Draw("APL");
  canvas->SaveAs("StatBinToTGraph.png");
}

TEST_F(DataContainerUnitTest, MultiDimProjection0d) {
  auto projection = statbin_nd.Projection();
  auto manual_sumwx = 0.;
  auto manual_sumw2 = 0.;
  auto manual_sumw = 0.;
  for (int i = 0; i < statbin_nd.size(); ++i) {
    manual_sumwx += statbin_nd[i].Mean() * statbin_nd[i].SumWeights();
    manual_sumw += statbin_nd[i].SumWeights();
    manual_sumw2 += statbin_nd[i].SumWeights2();
  }
  double manual_variance = 0;
  double t1 = 0;
  for (int i = 0; i < statbin_nd.size(); ++i) {
    t1 += statbin_nd[i].Variance()*(statbin_nd[i].SumWeights()) +
          statbin_nd[i].Mean()*statbin_nd[i].Mean()*
              statbin_nd[i].SumWeights();
  }
  auto manual_mean = manual_sumwx / manual_sumw;
  manual_variance = 1/(manual_sumw) * (t1 - manual_sumw*manual_mean*manual_mean);
  EXPECT_FLOAT_EQ(projection[0].Mean(), manual_mean);
  EXPECT_FLOAT_EQ(projection[0].Variance(), manual_variance);
}

TEST_F(DataContainerUnitTest, MultiDimProjection1d) {
  auto projection = statbin_nd.Projection({"a1"});
  auto manual_sumwx = 0.;
  auto manual_sumw = 0.;
  for (int i = 0; i < 2; ++i) {
    manual_sumwx += statbin_nd[i].Mean() * statbin_nd[i].SumWeights();
    manual_sumw += statbin_nd[i].SumWeights();
  }
  double manual_variance = 0;
  double t1 = 0;
  for (int i = 0; i < 2; ++i) {
    t1 += statbin_nd[i].Variance()*(statbin_nd[i].SumWeights()) +
        statbin_nd[i].Mean()*statbin_nd[i].Mean()*
              statbin_nd[i].SumWeights();
  }
  auto manual_mean = manual_sumwx / manual_sumw;
  manual_variance = 1/(manual_sumw) * (t1 - manual_sumw*manual_mean*manual_mean);
  EXPECT_FLOAT_EQ(projection[0].Mean(), manual_mean);
  EXPECT_FLOAT_EQ(projection[0].Variance(), manual_variance);
}

TEST_F(DataContainerUnitTest, EmptyBin) {
  using namespace Qn;
  Qn::DataContainerStatCollect stats_one_empty_bin;
  stats_one_empty_bin.AddAxes({AxisD("Centrality", 10, 0, 100)});
  std::normal_distribution<> dist{2., 0.5};
  std::poisson_distribution<> poisson{1};
  std::mt19937 gen{0};
  size_t nsamples = 100;
  std::vector<size_t> samples(nsamples);
  for (auto & stat : stats_one_empty_bin) {
    stat.SetNumberOfSamples(nsamples);
    stat.SetWeightType(Stat::WeightType::REFERENCE);
  }
  for (size_t i = 0; i < 1000l; ++i) {
    for (auto &sample : samples) {
      sample = poisson(gen);
    }
    for (int ibin = 0; ibin < stats_one_empty_bin.size(); ++ibin) {
      auto val = dist(gen);
      if (ibin == 1 || ibin == 0) {

      } else {
        stats_one_empty_bin[ibin].Fill(val, 1., samples);
      }
    }
  }
  auto statbin_one_empty = Qn::BinnedStatistics(stats_one_empty_bin);

  auto projection = statbin_one_empty.Projection();
  EXPECT_FALSE(std::isnan(projection[0].Mean()));
  EXPECT_FALSE(std::isnan(projection[0].StandardErrorOfMean()));
  EXPECT_NEAR(projection[0].Mean(), 2., 0.05);
}

TEST_F(DataContainerUnitTest, RebinReference) {
  using namespace Qn;
  DataContainerStatCollect r({AxisD("Centrality", 10, 0, 100)});
  std::normal_distribution<> dist{2.,1.};
  std::mt19937 gen{0};
  for (auto &s: r) {
    std::vector<size_t> multiplicities(10, 1);
    s.SetNumberOfSamples(10);
    s.SetWeightType(Stat::WeightType::REFERENCE);
    for (size_t i = 0; i < 100l; ++i) {
      s.Fill(dist(gen), 1., multiplicities);
    }
  }
  auto r_statbin = DataContainer<StatCalculate>(r);
  auto r_statbin_rebinned = r_statbin.Rebin(AxisD("Centrality", 5, 0, 100));
  for (StatCalculate &s : r_statbin_rebinned) {
    EXPECT_NEAR(s.Mean(), 2., 1.);
  }
}

TEST_F(DataContainerUnitTest, RebinProjection) {
  using namespace Qn;
  DataContainerStatCollect uQ({
                            AxisD("Centrality", 10, 0, 100),
                            AxisD("Pt", 20, 0., 2.),
                            AxisD("Y", 80, -4., 4.),
                        });
  DataContainerStatCollect r({AxisD("Centrality", 10, 0, 100)});
  std::normal_distribution<> dist{2.,1.};
  std::mt19937 gen{0};
  for (auto &s : uQ) {
    std::vector<size_t> multiplicities(10, 1);
    s.SetNumberOfSamples(10);
    s.SetWeightType(Stat::WeightType::OBSERVABLE);
    for (size_t i = 0; i < 100l; ++i) {
      s.Fill(dist(gen), 1., multiplicities);
    }
  }
  for (auto &s: r) {
    std::vector<size_t> multiplicities(10, 1);
    s.SetNumberOfSamples(10);
    s.SetWeightType(Stat::WeightType::REFERENCE);
    for (size_t i = 0; i < 100l; ++i) {
      s.Fill(dist(gen), 1., multiplicities);
    }
  }

  auto v1_stat_bin = Qn::DataContainerStatCalculate(uQ)/Qn::DataContainerStatCalculate(r);

  EXPECT_NO_THROW(v1_stat_bin.Rebin(AxisD("Centrality", 1, 0, 10)));
  EXPECT_NO_THROW(v1_stat_bin.Rebin(AxisD("Pt", 1, 0, 0.1)));
  EXPECT_NO_THROW(v1_stat_bin.Rebin(AxisD("Y", 1, 0, 0.1)));
  /* Merging two bins */
  EXPECT_NO_THROW(v1_stat_bin.Rebin(AxisD("Centrality", 1, 0, 20)));
  EXPECT_NO_THROW(v1_stat_bin.Rebin(AxisD("Pt", 1, 0, 0.2)));
  EXPECT_NO_THROW(v1_stat_bin.Rebin(AxisD("Y", 1, 0, 0.2)));
  EXPECT_NO_THROW(v1_stat_bin.Projection({"Centrality"}));
  EXPECT_NO_THROW(v1_stat_bin.Projection({"Pt"}));
  EXPECT_NO_THROW(v1_stat_bin.Projection({"Y"}));
}
