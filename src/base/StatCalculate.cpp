// Qn Tools
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

#include "StatCalculate.hpp"

namespace Qn {

/**
 * Default constructor
 */
StatCalculate::~StatCalculate() = default;

/**
 * Merges two StatBins. Returns the variance and mean from the pooled sample.
 * The sum of weights and sum of weights^2 are seperately accumulated to be able
 * to calculate the effective number of events, which is required to calculate
 * the standard error of the mean.
 * @param lhs
 * @param rhs
 * @return merged result.
 */
StatCalculate Merge(const StatCalculate &lhs, const StatCalculate &rhs) {
  StatCalculate merged;
  // mean and variance
  merged.SelectWeightandErrors(lhs, rhs);
  merged.sum_weight_ = lhs.SumWeights() + rhs.SumWeights();
  merged.sum_weight2_ = lhs.SumWeights2() + rhs.SumWeights2();
  if (merged.SumWeights() > 0 ) {
    merged.mean_ =
        (lhs.Mean() * lhs.SumWeights() + rhs.Mean() * rhs.SumWeights()) /
        merged.SumWeights();
    merged.variance_ =
        ((lhs.SumWeights()) * lhs.Variance() +
         (rhs.SumWeights()) * rhs.Variance() +
                        lhs.SumWeights() * lhs.Mean() * lhs.Mean() +
                        rhs.SumWeights() * rhs.Mean() * rhs.Mean() -
                        merged.SumWeights() * merged.Mean() * merged.Mean()
        ) / (merged.SumWeights());
  } else {
    merged.mean_ = 0.;
    merged.variance_ = 0.;
    merged.sum_weight_ = 0.;
    merged.sum_weight2_ = 0.;
  }
  // Bootstrap samples
  if (lhs.sum_weight_ == 0.) {
    merged.sample_means_ = rhs.sample_means_;
    merged.sample_weights_ = rhs.sample_weights_;
  } else if (rhs.sum_weight_ == 0.) {
    merged.sample_means_ = lhs.sample_means_;
    merged.sample_weights_ = lhs.sample_weights_;
  } else {
    for (int i = 0; i < rhs.sample_means_.size(); ++i) {
      merged.sample_weights_.push_back(lhs.sample_weights_[i] +
          rhs.sample_weights_[i]);
      merged.sample_means_.push_back(
          (   lhs.sample_means_[i] * lhs.sample_weights_[i]
              + rhs.sample_means_[i] * rhs.sample_weights_[i])
              / merged.sample_weights_[i]);
    }
  }
  return merged;
}

/**
 * Sum of two StatBins
 * @param lhs
 * @param rhs
 * @return left operand + right operand
 */
StatCalculate operator+(const StatCalculate &lhs, const StatCalculate &rhs) {
  StatCalculate sum(lhs, rhs);
  // mean and variance
  sum.mean_ = lhs.Mean() + rhs.Mean();
  sum.variance_ = lhs.Variance() + rhs.Variance();
  // Bootstrap samples
  for (int i = 0; i < lhs.sample_means_.size(); ++i) {
    sum.sample_weights_[i] = lhs.sample_weights_[i];
    sum.sample_means_[i] = lhs.sample_means_[i] + rhs.sample_means_[i];
  }
  return sum;
}

/**
 * Difference of two StatBins
 * @param lhs
 * @param rhs
 * @return left operand - right operand
 */
StatCalculate operator-(const StatCalculate &lhs, const StatCalculate &rhs) {
  StatCalculate difference(lhs, rhs);
  // mean and variance
  difference.mean_ = lhs.Mean() - rhs.Mean();
  difference.variance_ = lhs.Variance() + rhs.Variance();
  // Bootstrap samples
  for (int i = 0; i < lhs.sample_means_.size(); ++i) {
    difference.sample_weights_[i] = lhs.sample_weights_[i];
    difference.sample_means_[i] = lhs.sample_means_[i] - rhs.sample_means_[i];
  }
  return difference;
}

/**
 * Returns product of two StatBins
 * @param lhs left operand
 * @param rhs right operand
 * @return left operand * right operand
 */
StatCalculate operator*(const StatCalculate &lhs, const StatCalculate &rhs) {
  StatCalculate product(lhs, rhs);
  // mean and variance
  product.mean_ = lhs.Mean() * rhs.Mean();
  product.variance_ =   lhs.Variance() * rhs.Mean() * rhs.Mean()
      + rhs.Variance() * lhs.Mean() * lhs.Mean();
  // Bootstrap samples
  for (int i = 0; i < lhs.sample_means_.size(); ++i) {
    product.sample_weights_[i] = lhs.sample_weights_[i];
    product.sample_means_[i] = lhs.sample_means_[i] * rhs.sample_means_[i];
  }
  return product;
}

/**
 * Returns ratio of two StatBins
 * @param num numerator
 * @param den denominator
 * @return numerator / denominator
 */
StatCalculate operator/(const StatCalculate &num, const StatCalculate &den) {
  StatCalculate ratio(num, den);
  // mean and variance
  ratio.mean_ = num.Mean() / den.Mean();
  ratio.variance_ =   num.Variance() / (den.Mean() * den.Mean())
      + num.Mean() * num.Mean() * den.Variance() / std::pow(den.Mean(), 4);
  // Bootstrap samples
  for (int i = 0; i < num.sample_means_.size(); ++i) {
    ratio.sample_weights_[i] = num.sample_weights_[i];
    ratio.sample_means_[i] = num.sample_means_[i] / den.sample_means_[i];
  }
  return ratio;
}

/**
 * Returns Scaled StatCalculate StatCalculate * scale.
 * @param operand
 * @param scale
 * @return scale * StatCalculate.
 */
StatCalculate operator*(const StatCalculate &operand, double scale) {
  StatCalculate scaled(operand);
  // mean and variance
  scaled.mean_ = operand.Mean() * scale;
  scaled.variance_ = operand.Variance() * scale * scale;
  // Bootstrap samples
  for (int i = 0; i < operand.sample_means_.size(); ++i) {
    scaled.sample_weights_[i] = operand.sample_weights_[i];
    scaled.sample_means_[i] = operand.sample_means_[i] * scale;
  }
  return scaled;
}

/**
 * Returns Scaled StatCalculate StatCalculate * scale.
 * @param operand
 * @param scale
 * @return StatCalculate * scale.
 */
StatCalculate operator*(double scale, const StatCalculate &operand) {
  return operand * scale;
}

/**
 * Returns Scaled StatCalculate StatCalculate / scale.
 * @param operand
 * @param scale
 * @return StatCalculate / scale.
 */
StatCalculate operator/(const StatCalculate &operand, double scale) {
  return operand * (1./scale);
}

/**
 * Returns StatCalculate^{exponent}.
 * @param base StatCalculate base of the operation.
 * @param exp Exponent of the operation.
 * @return base^{exp}
 */
StatCalculate Pow(const StatCalculate &base, double exp) {
  StatCalculate result(base);
  // mean and variance
  result.mean_ = std::pow(base.Mean(), exp);
  result.variance_ =   exp / base.Mean() * result.Mean()
      * exp / base.Mean() *
      result.Mean()
      * base.Variance();
  // Bootstrap samples
  for (int i = 0; i < base.sample_means_.size(); ++i) {
    result.sample_weights_[i] = base.sample_weights_[i];
    result.sample_means_[i] = std::pow(base.sample_means_[i],exp);
  }
  return result;
}

/**
 * Specialization of Pow for Sqrt.
 * @param operand
 * @return sqrt(StatCalculate)
 */
StatCalculate Sqrt(const StatCalculate &operand) {
  return Pow(operand,(1./2));
}


/**
 * Returns OBSERVABLE type StatCalculate
 * @param lhs left candidate
 * @param rhs right candidate
 * @return returns OBSERVABLE type StatCalculate prefers lhs over rhs in case none is found.
 */
StatCalculate StatCalculate::PreferObservable(const StatCalculate &lhs, const StatCalculate &rhs) {
  if ( rhs.GetWeightType() == StatCalculate::WeightType::OBSERVABLE) {
    return rhs;
  } else {
    return lhs;
  }
}


}
