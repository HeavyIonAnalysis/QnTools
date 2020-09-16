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

#ifndef QNTOOLS_BOOTSTRAP_HPP_
#define QNTOOLS_BOOTSTRAP_HPP_

#include <vector>
#include <algorithm>
#include <cmath>

#include "Rtypes.h"

namespace Qn {

class BootStrap {
 public:
  template <typename SAMPLES>
  void Fill(const double value, const double weight,
            SAMPLES &&sample_multiplicities_) {
    // Fill Bootstrapping samples
    using ArrayValueType =
        typename std::decay<decltype(std::declval<SAMPLES &>()[0])>::type;
    for (std::size_t i = 0; i < vector_sum_values_.size(); ++i) {
      for (ArrayValueType j = 0; j < sample_multiplicities_[i]; ++j) {
        vector_sum_values_[i] += value * weight;
        vector_sum_weights_[i] += weight;
      }
    }
  }

  void SetNumberOfSamples(int n_samples_) {
    vector_sum_values_.resize(n_samples_);
    vector_sum_weights_.resize(n_samples_);
  }

  [[nodiscard]] std::vector<double> GetMeans() const {
    size_t size = vector_sum_values_.size();
    std::vector<double> means_(size);
    for (int i = 0; i < size; ++i) {
      means_[i] = vector_sum_values_[i] / vector_sum_weights_[i];
    }
    return means_;
  }

  [[nodiscard]] std::vector<double> GetWeights() const { return vector_sum_weights_; }

  friend BootStrap Merge(const BootStrap &lhs, const BootStrap &rhs);
  friend BootStrap MergeBins(const BootStrap &lhs, const BootStrap &rhs);

 private:
  std::vector<double> vector_sum_values_;
  std::vector<double> vector_sum_weights_;

  /// \cond CLASSIMP
 ClassDef(BootStrap, 2);
  /// \endcond
};

BootStrap Merge(const BootStrap &lhs, const BootStrap &rhs);
BootStrap MergeBins(const BootStrap &lhs, const BootStrap &rhs);

}
#endif  // QNTOOLS_BOOTSTRAP_HPP_
