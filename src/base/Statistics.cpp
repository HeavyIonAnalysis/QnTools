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

#include "Statistics.hpp"

namespace Qn {

void Statistics::Fill(double value, double weight) {
  if (weight == 0) return;
  ++n_entries_;
  min_ = std::min(min_, value);
  max_ = std::max(max_, value);
  auto old_weights = sum_weights_;
  auto newvalue = value * weight;
  sum_weights_ += weight;
  auto num = weight * sum_values_ - old_weights * newvalue;
  if (old_weights != 0)
    sum_sq_ += num * num / (sum_weights_ * weight * old_weights);
  sum_weights2_ += weight * weight;
  sum_values_ += newvalue;
}

Statistics Merge(const Statistics &lhs, const Statistics &rhs) {
  Statistics result = lhs;
  result.sum_weights_ += rhs.sum_weights_;
  result.n_entries_ += rhs.n_entries_;
  result.sum_weights2_ += rhs.sum_weights2_;
  result.sum_values_ += rhs.sum_values_;
  result.max_ = std::max(result.max_, rhs.max_);
  result.min_ = std::min(result.min_, rhs.min_);
  double num =
      rhs.sum_weights_ * lhs.sum_values_ - lhs.sum_weights_ * rhs.sum_values_;
  result.sum_sq_ = lhs.sum_sq_ + rhs.sum_sq_;
  if (lhs.sum_weights_ != 0. && rhs.sum_weights_ != 0. &&
      result.sum_weights_ != 0.) {
    result.sum_sq_ += (num * num) / (lhs.sum_weights_ * rhs.sum_weights_ *
                                     result.sum_weights_);
  }
  return result;
}

} // namespace Qn
