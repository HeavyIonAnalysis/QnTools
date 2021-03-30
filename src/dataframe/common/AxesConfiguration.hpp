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
#ifndef QNTOOLS_AXESCONFIGURATION_H_
#define QNTOOLS_AXESCONFIGURATION_H_

#include <array>
#include <string>
#include <tuple>
#include <vector>

#include "Axis.hpp"
#include "TemplateFunctions.hpp"

namespace Qn {
namespace Impl {

/**
 * template class StrideCalculator
 * @tparam position recursion step
 */
template <std::size_t position>
struct StrideCalculator {
  template <typename Stride, typename Axes>
  static Stride Calculate(Stride &stride, const Axes &axes) {
    stride[position] = stride[position + 1] * std::get<position>(axes).size();
    return StrideCalculator<position - 1>::Calculate(stride, axes);
  }
};

/**
 * base template class StrideCalculator
 */
template <>
struct StrideCalculator<0> {
  template <typename Stride, typename Axes>
  static Stride Calculate(Stride &stride, const Axes &axes) {
    stride[0] = stride[1] * std::get<0>(axes).size();
    return stride;
  }
};

}  // namespace Impl

/**
 * Template class for the configuration of the event axes. Calculates the
 * binning of a multidimensional tensor in row-major using the specified
 * binning. Does not store the bin data itself.
 * @tparam Axes type and number of the axes.
 */
template <typename... Axes>
class AxesConfiguration {
 public:
  static constexpr auto kDimension =
      sizeof...(Axes);  /// Dimensionality of the matrix
  using AxisTuple = typename std::tuple<Axes...>;
  using AxisType = typename std::tuple_element<0, std::tuple<Axes...>>::type;
  using AxisValueType = typename AxisType::ValueType;
  using AxisValueTypeTuple =
      typename TemplateFunctions::TupleOf<sizeof...(Axes), AxisValueType>;

  /**
   * Constructor
   * @param axes Axes going into the construction of the AxisConfiguration
   * Calculates the stride of the axis configuration. The Axes and the stride
   * are both const.
   */
  explicit AxesConfiguration(Axes... axes)
      : axes_(axes...), stride_(CalculateStride()) {}

  /**
   * Returns the axes inform of a vector. Can be used to configure the
   * RDataFrame during runtime.
   * @return vector of axes.
   */
  std::vector<AxisType> GetVector() const {
    return TemplateFunctions::ToVector(axes_);
  }

  /**
   * Get the total size
   */
  [[nodiscard]] std::size_t GetSize() const { return stride_[0]; }
  /**
   * Returns the linear bin offset for the axes configuration from passed
   * coordinates.
   * @tparam Coordinates type of the coordinates
   * @param coordinates coordinates for which the corresponding bin is searched.
   * @return linear index of the bin.
   */
  template <typename... Coordinates>
  long GetLinearIndex(Coordinates &&... coordinates) const {
    return FindBin(std::forward<Coordinates>(coordinates)...);
  }

  /**
   * Returns a mapping between the linear index and the corresponding upper and
   * lower bin edges.
   * @return mapping between the linear index and the corresponding upper and
   * lower bin edges
   */
  auto GetBinEdgesIndexMap() const {
    using edgest =
        typename Qn::TemplateFunctions::TupleOf<kDimension,
                                              typename AxisType::ValueType>;
    std::vector<std::pair<edgest, edgest>> index_map(GetSize());
    for (int i = 0; i < GetSize(); ++i) {
      AxisValueTypeTuple lower;
      AxisValueTypeTuple higher;
      auto ndimbins = GetNDimBins(i);
      auto caller = [](auto &lower, auto &higher, auto axis, auto bin) {
        lower = axis.GetLowerBinEdge(bin);
        higher = axis.GetUpperBinEdge(bin);
      };
      Qn::TemplateFunctions::TupleForEach(caller, lower, higher, axes_, ndimbins);
      index_map.at(i) = std::pair(lower, higher);
    }
    return index_map;
  }

  /**
   * Returns a vector of names of the axes.
   * @return vector of names.
   */
  [[nodiscard]] std::vector<std::string> GetNames() const {
    const auto axes_vector = GetVector();
    std::vector<std::string> axes_names_;
    for (const auto &axis : axes_vector) {
      axes_names_.push_back(axis.Name());
    }
    return axes_names_;
  }

  /**
   * Returns the Tuple of axes. Used to optimize the binning.
   * @return the tuple of axes.
   */
  const std::tuple<Axes...> &GetAxes() const { return axes_; }

  /**
   * Returns the Tuple of axes. Used to optimize the binning.
   * @return the tuple of axes.
   */
  std::tuple<Axes...> &GetAxes() { return axes_; }

 private:
  /**
   * Recursive function to calculate the bins.
   * Recursion
   * @tparam FirstCoordinate type of the first coordinates.
   * @tparam Rest types of the rest of coordinates.
   * @param first_coordinate first coordinate
   * @param rest rest of coordinates
   * @return bin
   */
  template <typename FirstCoordinate, typename... Rest>
  long FindBin(FirstCoordinate &&first_coordinate, Rest &&... rest) const {
    constexpr std::size_t position = kDimension - sizeof...(rest) - 1;
    auto bin = stride_[position + 1] *
               std::get<position>(axes_).FindBin(
                   std::forward<FirstCoordinate>(first_coordinate));
    auto rest_bin = FindBin(std::forward<Rest>(rest)...);
    if (rest_bin < 0) return -1;
    return bin + rest_bin;
  }

  /**
   * Recursive function to calculate the bins.
   * base case
   * @tparam FirstCoordinate type of the first coordinates.
   * @param first_coordinate first coordinate
   * @return bin
   */
  template <typename FirstCoordinate>
  long FindBin(FirstCoordinate &&first_coordinate) const {
    constexpr std::size_t position = kDimension - 1;
    return stride_[position + 1] *
           std::get<position>(axes_).FindBin(
               std::forward<FirstCoordinate>(first_coordinate));
  }

  /**
   * Calculates the n dimensional bin from a linear_index
   * @param linear_index linear bin index
   * @return tuple of length n with index of axis in each entry
   */
  auto GetNDimBins(int linear_index) const {
    Qn::TemplateFunctions::TupleOf<kDimension, std::size_t> tuple;
    int dimension = 1;
    auto caller = [&](std::size_t &t) {
      std::size_t index = linear_index / stride_[dimension];
      linear_index = linear_index % stride_[dimension];
      ++dimension;
      t = index;
    };
    Qn::TemplateFunctions::TupleForEach(caller, tuple);
    return tuple;
  }

  /**
   * Calculates the strides.
   * @return return the array of strides.
   */
  std::array<std::size_t, kDimension + 1> CalculateStride() {
    std::array<std::size_t, kDimension + 1> temp;
    temp[kDimension] = 1;
    return Impl::StrideCalculator<kDimension - 1>::Calculate(temp, axes_);
  }

  std::tuple<Axes...> axes_;  /// tuple of event axes.
  std::array<std::size_t, kDimension + 1>
      stride_;  /// array of strides used for linear bin calculation.
};

/**
 * Wrapper function to create the axes configuration.
 * @tparam Axes types of the axes.
 * @param axes Axes used for the construction of the configuration.
 * @return The configuration.
 */
template <typename... Axes>
auto MakeAxes(Axes... axes) {
  return AxesConfiguration<Axes...>(axes...);
}

}  // namespace Qn
#endif  // QNTOOLS_AXESCONFIGURATION_H_
