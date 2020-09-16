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

#ifndef QNTOOLS_STATBIN_H_
#define QNTOOLS_STATBIN_H_

#include <cmath>
#include <vector>

#include "Rtypes.h"
#include "Stat.hpp"
#include "StatCollect.hpp"
#include "Statistics.hpp"

namespace Qn {

class StatCalculate : public Stat {
 public:

  /// Default constructor.
  StatCalculate() = default;

  /// Construct StatCalculate from Stats.
  explicit StatCalculate(StatCollect &stats) {
    mean_     = stats.GetStatistics().Mean();
    variance_ = stats.GetStatistics().Variance();
    sum_weight_ = stats.GetStatistics().SumWeights();
    sum_weight2_ = stats.GetStatistics().SumWeights2();
    const auto& samples = stats.GetBootStrap();
    sample_means_ = samples.GetMeans();
    sample_weights_ = samples.GetWeights();
    if (stats.GetWeightType() == Stat::WeightType::OBSERVABLE) {
      SetWeightType(WeightType::OBSERVABLE);
    } else {
      SetWeightType(WeightType::REFERENCE);
    }
  }

  virtual ~StatCalculate() = default;

  /// Returns the sample Variance from error propagation
  [[nodiscard]] double VarianceFromPropagation() const { return variance_; }

  /// Returns variance of the sample mean from error propagation
  [[nodiscard]] double VarianceOfMeanFromPropagation() const {
    if (sum_weight_ > 0.) {
      return VarianceFromPropagation() * NEffectiveInverse();
    } else {
      return 0.;
    }
  }

  /// Returns variance of the sample mean from bootstrapping using the variance statistic.
  [[nodiscard]] double VarianceOfMeanFromBootstrap() const {
    Statistics stats;
    for (int i = 0; i < sample_means_.size(); ++i) {
      stats.Fill(sample_means_[i], sample_weights_[i]);
    }
    return stats.Variance();
  }

  /// Returns standard error of the sample mean from bootstrapping using the percentile method.
  [[nodiscard]] double StdDevOfMeanFromBootstrapPercentile() const {
    auto means = sample_means_;
    const double alpha = 0.3173;
    std::sort(means.begin(), means.end());
    auto lower = means.at(std::floor(means.size() * alpha / 2.));
    auto upper = means.at(std::ceil(means.size() * (1 - alpha / 2.)));
    return (upper - lower) / 2;
  }

  /// Returns standard error from error propagation
  [[nodiscard]] double StdDevFromPropagation() const {
    return std::sqrt(VarianceFromPropagation());
  }

  /// Returns the standard error of the mean from error propagation
  [[nodiscard]] double StdDevOfMeanFromPropagation() const {
    return std::sqrt(VarianceOfMeanFromPropagation());
  }

  /// Returns the standard error of the mean from bootstrapping using the variance statistic.
  [[nodiscard]] double StdDevOfMeanFromBootstrapVariance() const {
    return std::sqrt(VarianceOfMeanFromBootstrap());
  }
  /// Returns the sample mean.
  [[nodiscard]] double Mean() const { return mean_; }

  /// Returns the sample variance.
  [[nodiscard]] double Variance() const { return variance_; }

  /// Returns the sum of weights.
  [[nodiscard]] double SumWeights() const { return sum_weight_; }

  /// Returns the sum of weights squared.
  [[nodiscard]] double SumWeights2() const { return sum_weight2_; }

  /// Retunrs the effective number of events.
  /// See For more information:
  /// https://en.wikipedia.org/w/index.php?title=Effective_sample_size&oldid=976981114#Weighted_samples
  [[nodiscard]] double NEffectiveInverse() const { return SumWeights2() / (SumWeights() * SumWeights()); }

  /// Returns standard error of the mean using the chosen method.
  [[nodiscard]] double StandardErrorOfMean() const {
    if (GetErrorType() == ErrorType::PROPAGATION) return StdDevOfMeanFromPropagation();
    else                                          return StdDevOfMeanFromBootstrapVariance();
  }

  /// Merge using pooled statistics
  /// For merging of pooled variance see:
  /// https://en.wikipedia.org/w/index.php?title=Pooled_variance&oldid=953883248#Sample-based_statistics
  friend StatCalculate Merge(const StatCalculate &, const StatCalculate &);
  friend StatCalculate MergeBins(const StatCalculate &, const StatCalculate &);

  /// Implements algebra of random variables
  /// For algebra of random variables used here see:
  /// https://en.wikipedia.org/w/index.php?title=Algebra_of_random_variables&oldid=945632545
  friend StatCalculate operator+(const StatCalculate &, const StatCalculate &);
  friend StatCalculate operator-(const StatCalculate &, const StatCalculate &);
  friend StatCalculate operator*(const StatCalculate &, const StatCalculate &);
  friend StatCalculate operator/(const StatCalculate &, const StatCalculate &);
  friend StatCalculate operator/(const StatCalculate &, double);
  friend StatCalculate operator*(const StatCalculate &, double);
  friend StatCalculate operator*(double, const StatCalculate &);
  friend StatCalculate Pow(const StatCalculate &, double);
  friend StatCalculate Sqrt(const StatCalculate &);

 private:
  std::vector<double> sample_means_; /// means of bootstrap samples
  std::vector<double> sample_weights_; /// weights of bootstrap samples
  double mean_ = 0.; /// mean of the distribution
  double variance_ = 0.; /// variance of the distribution
  double sum_weight_ = 0.; /// sum of weights
  double sum_weight2_ = 0.; /// sum of {weights}^{2}

  void CopySettings(const StatCalculate &other);
  void CopySettings(const StatCalculate &lhs, const StatCalculate &rhs);

  /// \cond CLASSIMP
 ClassDef(StatCalculate, 2);
  /// \endcond
};

StatCalculate Merge(const StatCalculate &lhs, const StatCalculate &rhs);

StatCalculate operator+(const StatCalculate &lhs, const StatCalculate &rhs);
StatCalculate operator-(const StatCalculate &lhs, const StatCalculate &rhs);

StatCalculate operator*(const StatCalculate &lhs, const StatCalculate &rhs);
StatCalculate operator/(const StatCalculate &num, const StatCalculate &den);

StatCalculate operator*(const StatCalculate &operand, double scale);
StatCalculate operator*(double operand, const StatCalculate &rhs);
StatCalculate operator/(const StatCalculate &operand, double scale);

StatCalculate Pow(const StatCalculate &base, double exp);
StatCalculate Sqrt(const StatCalculate &operand);

}

#endif  // QNTOOLS_INCLUDE_BASE_QNTOOLS_STATBIN_H_
