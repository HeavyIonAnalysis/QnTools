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

#ifndef QNTOOLS_STAT_HPP_
#define QNTOOLS_STAT_HPP_

namespace Qn {

class Stat {
 public:
  /// WeightType used for merging.
  /// OBSERVABLE types are preferred over REFERENCE type weights when merging.
  enum class WeightType {
    REFERENCE,  /// lower priority: event weight not preferred in merging
    OBSERVABLE  /// higher priority: event weight preferred in merging
  };
  /// Type of Error used for plotting the results.
  enum class ErrorType {
    BOOTSTRAP,   /// Using uncertainties bootstrapping
    PROPAGATION  /// Using uncertainties from error propagation
  };

  virtual ~Stat() = default;

  /// Sets the weight type
  void SetWeightType(WeightType type) { type_of_weight_ = type; }

  /// Gets the weight type
  [[nodiscard]] auto GetWeightType() const { return type_of_weight_; }

  /// Sets the error type.
  void SetErrorType(ErrorType type) { type_of_error_ = type; }

  /// Returns the error type.
  [[nodiscard]] ErrorType GetErrorType() const { return type_of_error_; }

  /// Selects the weight and error with the priority
  /// lhs > rhs and observable > reference.
  void SelectWeightandErrors(const Stat &lhs, const Stat &rhs) {
    if (lhs.GetWeightType() == Stat::WeightType::OBSERVABLE) {
      type_of_weight_ = Stat::WeightType::OBSERVABLE;
      type_of_error_ = lhs.type_of_error_;
    } else if (rhs.GetWeightType() == Stat::WeightType::OBSERVABLE){
      type_of_weight_ = Stat::WeightType::OBSERVABLE;
      type_of_error_ = rhs.type_of_error_;
    } else {
      type_of_weight_ = Stat::WeightType::REFERENCE;
      type_of_error_ = lhs.type_of_error_;
    }
  }

 protected:
  /// Type of weight used for merging
  WeightType type_of_weight_ = WeightType::REFERENCE;
  /// Type of error used for plotting
  ErrorType type_of_error_ = ErrorType::BOOTSTRAP;

  /// \cond CLASSIMP
  ClassDef(Stat, 2);
  /// \endcond
};

}
#endif // QNTOOLS_STAT_HPP_
