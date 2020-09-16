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

#include "BootStrap.hpp"

namespace Qn {

BootStrap Merge(const BootStrap &lhs, const BootStrap &rhs) {
  BootStrap result;
  // Bootstrap samples
  if (lhs.vector_sum_weights_.empty()) {
    result.vector_sum_values_ = rhs.vector_sum_values_;
    result.vector_sum_weights_ = rhs.vector_sum_weights_;
  } else if (rhs.vector_sum_weights_.empty()) {
    result.vector_sum_values_ = lhs.vector_sum_values_;
    result.vector_sum_weights_ = lhs.vector_sum_weights_;
  } else {
    for (int i = 0; i < rhs.vector_sum_values_.size(); ++i) {
      result.vector_sum_weights_.push_back(lhs.vector_sum_weights_[i] +
                                           rhs.vector_sum_weights_[i]);
      result.vector_sum_values_.push_back(lhs.vector_sum_values_[i] +
                                          rhs.vector_sum_values_[i]);
    }
  }
  return result;
}

BootStrap MergeBins(const BootStrap &lhs, const BootStrap &rhs) {
  return Merge(lhs, rhs);
}

} // namespace Qn
