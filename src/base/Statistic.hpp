// Qn Tools
//
// Copyright (C) 2019  Lukas Kreis, Ilya Selyuzhenkov
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

#ifndef FLOW_STATISTIC_HPP
#define FLOW_STATISTIC_HPP

#include <algorithm>
#include <cmath>

namespace Qn {

class Statistic {
 public:
  void Fill(double value, double weight) {
    if (weight == 0) return;
    ++n_entries_;
    min_ = std::min(min_, value);
    max_ = std::max(max_, value);
    auto old_weights = sum_weights_;
    auto newvalue = value * weight;
    sum_weights_ += weight;
    auto num = weight * sum_values_ - old_weights * newvalue;
    if (old_weights != 0) sum_sq_ += num * num / (sum_weights_ * weight * old_weights);
    sum_weights2_ += weight * weight;
    sum_values_ += newvalue;
  }
  double SumSq() const { return sum_sq_; }
  double Neff() const { return sum_weights2_ > 0 ? sum_weights_ * sum_weights_ / sum_weights2_ : 0.; }
  double N() const { return n_entries_; }
  double Variance() const { return sum_weights_ > 0 ? (n_entries_ > 1 ? (sum_sq_ / (sum_weights_ - 1)) : 0) : -1; }
  double Sigma() const {
    double var = Variance();
    return var > 0 ? std::sqrt(var) : -1;
  }
  double SumWeights() const { return sum_weights_; }
  double Mean() const { return sum_weights_ > 0 ? sum_values_ / sum_weights_ : 0.0; }
  double MeanError() const { return Neff() > 0 ? std::sqrt(Variance() / Neff()) : 0.0; }
  double Entries() const { return n_entries_; }
  double Min() const { return min_; }
  double Max() const { return max_; }
  friend Statistic Merge(const Statistic &lhs, const Statistic &rhs);
  friend Statistic MergeBins(const Statistic &lhs, const Statistic &rhs);

 private:
  double sum_values_ = 0;
  double sum_sq_ = 0;
  double sum_weights_ = 0.;
  double sum_weights2_ = 0.;
  double n_entries_ = 0;
  double min_ = std::numeric_limits<double>::max();
  double max_ = std::numeric_limits<double>::min();
};

inline Statistic Merge(const Statistic &lhs, const Statistic &rhs) {
  Statistic result = lhs;
  result.sum_weights_ += rhs.sum_weights_;
  result.n_entries_ += rhs.n_entries_;
  result.sum_weights2_ += rhs.sum_weights2_;
  result.sum_values_ += rhs.sum_values_;
  result.max_ = std::max(result.max_, rhs.max_);
  result.min_ = std::min(result.min_, rhs.min_);
  double num = rhs.sum_weights_ * lhs.sum_values_ - lhs.sum_weights_ * rhs.sum_values_;
  result.sum_sq_ = lhs.sum_sq_ + rhs.sum_sq_;
  if (lhs.sum_weights_ != 0. && rhs.sum_weights_ != 0. && result.sum_weights_ != 0.) {
    result.sum_sq_ += (num * num) / (lhs.sum_weights_ * rhs.sum_weights_ * result.sum_weights_);
  }
  return result;
}

inline Statistic MergeBins(const Statistic &lhs, const Statistic &rhs) {
  Statistic result = lhs;
  result.sum_weights_ += rhs.sum_weights_;
  result.n_entries_ += rhs.n_entries_;
  result.sum_weights2_ += rhs.sum_weights2_;
  result.sum_values_ += rhs.sum_values_;
  result.max_ = std::max(result.max_, rhs.max_);
  result.min_ = std::min(result.min_, rhs.min_);
  double num = rhs.sum_weights_ * lhs.sum_values_ - lhs.sum_weights_ * rhs.sum_values_;
  result.sum_sq_ = lhs.sum_sq_ + rhs.sum_sq_;
  if (lhs.sum_weights_ != 0. && rhs.sum_weights_ != 0. && result.sum_weights_ != 0.) {
    result.sum_sq_ += (num * num) / (lhs.sum_weights_ * rhs.sum_weights_ * result.sum_weights_);
  }
  return result;
}

}// namespace Qn

#endif