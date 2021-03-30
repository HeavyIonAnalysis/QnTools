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

#ifndef QNTOOLS_QVECTORNORMALIZATIONFUNCTOR_HPP_
#define QNTOOLS_QVECTORNORMALIZATIONFUNCTOR_HPP_

#include "DataContainer.hpp"
#include "QVector.hpp"

namespace Qn {
class QVectorNormalizationFunctor {
 public:
  using Normalization = typename Qn::QVector::Normalization;

  explicit QVectorNormalizationFunctor(Normalization normalization)
      : normalization_(normalization) {}

  [[nodiscard]] std::string GetName(const std::string& name) const {
    return name + "_norm_" + Qn::QVector::GetNormalizationName(normalization_);
  }

  Qn::DataContainerQVector operator()(
      const Qn::DataContainerQVector& q_vector) {
    auto ret = q_vector;
    for (unsigned int i = 0; i < q_vector.size(); ++i) {
      ret[i] = q_vector[i].Normal(normalization_);
    }
    return ret;
  }

 private:
  Normalization normalization_;
};

}  // namespace Qn
#endif  // QNTOOLS_QVECTORNORMALIZATIONFUNCTOR_HPP_
