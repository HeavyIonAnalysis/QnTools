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

#include <vector>
#include <algorithm>
#include <cmath>

#include "Rtypes.h"


namespace Qn {

class Statistics {
 public:
  /**
   * Fill variable and weight to the statistics.
   * @param value value of event
   * @param weight weight of event
   */
  void Fill(double value, double weight);

  /// Returns second moment of distribution.
  [[nodiscard]] double SumSq() const { return sum_sq_; }
  /// Returns effective number of events.
  [[nodiscard]] double Neff() const { return sum_weights2_ > 0 ? sum_weights_ * sum_weights_ / sum_weights2_ : 0.; }
  /// Returns number of events.
  [[nodiscard]] double N() const { return n_entries_; }
  /// Returns variance of distribution.
  [[nodiscard]] double Variance() const { return sum_weights_ > 0. ? (n_entries_ > 1 ? (sum_sq_ / (sum_weights_ - 1)) : 0) : 0.; }
  /// Returns standard deviation of distribution.
  [[nodiscard]] double StandardDeviation() const { return Variance() > 0 ? std::sqrt(Variance()) : -1; }
  /// Returns sum of weights.
  [[nodiscard]] double SumWeights() const { return sum_weights_; }
  /// Returns sum of weights * weights.
  [[nodiscard]] double SumWeights2() const { return sum_weights2_; }
  /// Returns mean of the distribution.
  [[nodiscard]] double Mean() const { return sum_weights_ > 0 ? sum_values_ / sum_weights_ : 0.0; }
  /// Return error of the mean.
  [[nodiscard]] double StandardErrorOfMean() const { return Neff() > 0 ? std::sqrt(Variance() / Neff()) : 0.0; }
  /// Returns number of entries in the distribution.
  [[nodiscard]] double Entries() const { return n_entries_; }
  /// Returns the maximum value of the distribution.
  [[nodiscard]] double Min() const { return min_; }
  /// Returns the minimum value of the distribution.
  [[nodiscard]] double Max() const { return max_; }

  friend Statistics Merge(const Statistics &lhs, const Statistics &rhs);
  friend Statistics MergeBins(const Statistics &lhs, const Statistics &rhs);

 private:
  double sum_values_ = 0;
  double sum_sq_ = 0;
  double sum_weights_ = 0.;
  double sum_weights2_ = 0.;
  double n_entries_ = 0;
  double min_ = std::numeric_limits<double>::max();
  double max_ = std::numeric_limits<double>::min();

  /// \cond CLASSIMP
  ClassDef(Statistics, 2);
  /// \endcond
};

Statistics Merge(const Statistics &lhs, const Statistics &rhs);

}// namespace Qn

#endif
