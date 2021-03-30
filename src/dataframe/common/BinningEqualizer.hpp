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
#ifndef QNTOOLS_EQUALBINNINGHELPER_H_
#define QNTOOLS_EQUALBINNINGHELPER_H_

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <Math/Interpolator.h>
#include <ROOT/RResultPtr.hxx>

#include "Axis.hpp"

namespace Qn {

/**
 * Template class to facilitate equal number of events per bin of the correction
 * axes. Currently only performs bin optimization in one dimension.
 * @tparam AxesConfig type of the AxesConfiguration
 * @tparam DataFrame  type of the RDataFrame
 */
template <typename AxesConfig, typename DataFrame>
class BinningEqualizer {
 public:
  using AxisTuple =
      typename AxesConfig::AxisTuple;  /// Tuple containing the axes.

  /**
   * Constructor
   * @param df RDataFrame containing the input information
   * @param to_equalize Vector of names of the event axes which are supposed to
   * be optimized.
   */
  BinningEqualizer(DataFrame df, const std::vector<std::string> &to_equalize)
      : dataframe_(df) {
    for (const auto &name : to_equalize)
      value_map_.emplace(name, df.template Take<double>(name));
  }

  /**
   * Performs the optimization and returns the AxesConfiguration with the
   * optimized binning.
   * @return AxesConfiguration with the optimized binning.
   */
  void operator()(AxesConfig &config) { Equalize(config.GetAxes()); }

 private:
  DataFrame dataframe_;  /// Input data
  std::map<std::string, ROOT::RDF::RResultPtr<std::vector<double>>>
      value_map_;  /// map of vectors of event parameters

  /**
   * Recursion base case.
   */
  template <std::size_t I = 0>
  inline typename std::enable_if<I == std::tuple_size<AxisTuple>{}, void>::type
  Equalize(AxisTuple &axis_tuple) {
    (void)axis_tuple;
  }

  /**
   * Recursion.
   * Applies the equalization of number of entries algorithm on all specified
   * axes. Modifes the axes in axis_tuple.
   * @tparam I recursion step
   */
  template <std::size_t I = 0>
      inline typename std::enable_if <
      I<std::tuple_size<AxisTuple>{}, void>::type Equalize(
          AxisTuple &axis_tuple) {
    auto &axis = std::get<I>(axis_tuple);
    if (value_map_.find(axis.Name()) != value_map_.end()) {
      axis = OptimizeBins(value_map_[axis.Name()].GetValue(), axis);
    }
    Equalize<I + 1>(axis_tuple);
  }

  /**
   * Equalizes the bin of the specified axis reusing the axis borders and the
   * number of bins.
   * @param data input data used to perform the optimization
   * @param axis Axis to be optimized.
   * @return Optimized axis
   */
  Qn::AxisD OptimizeBins(std::vector<double> data, const Qn::AxisD &axis) {
    auto nbins = axis.size();
    auto low = axis.GetFirstBinEdge();
    auto high = axis.GetLastBinEdge();
    std::sort(data.begin(), data.end());
    if (low < data.front()) {
      data.insert(data.begin(), low);
    }
    if (high > data.back()) {
      data.push_back(high);
    }
    auto npoints = data.size();
    auto position_low = std::distance(
        data.begin(), std::lower_bound(data.begin(), data.end(), low));
    auto position_high = std::distance(
        data.begin(), std::lower_bound(data.begin(), data.end(), high));
    std::vector<double> x(npoints);
    std::iota(x.begin(), x.end(), 0);
    ROOT::Math::Interpolator interpolator(
        x, data, ROOT::Math::Interpolation::Type::kLINEAR);
    auto evalat = LinearSpacing(position_low, position_high, nbins + 1);
    std::vector<double> bin_edges(nbins + 1);
    for (unsigned int ibin = 0; ibin < nbins + 1; ++ibin) {
      bin_edges[ibin] = interpolator.Eval(evalat[ibin]);
    }
    bin_edges.front() = low;
    bin_edges.back() = high;
    return {axis.Name(), bin_edges};
  }

  /**
   * Calculates linearly spaced floating point numbers between two floating
   * point values.
   * @param start start of the interval
   * @param end end of the interval
   * @param num number of points
   * @return vector of linearly spaced number between start and end.
   */
  std::vector<double> LinearSpacing(double start, double end, int num) {
    int partitions = num - 1;
    std::vector<double> positions;
    double length = (end - start) / partitions;
    positions.push_back(start);
    for (int i = 1; i < num - 1; i++) {
      positions.push_back(start + i * length);
    }
    positions.push_back(end);
    return positions;
  }
};

/**
 * Function which returns a axesconfig with optimized bins
 * @tparam AxesConfig type of the AxisConfiguration
 * @tparam DataFrame type of the RDataFrame
 * @param config axis configuration which is supposed to be optimized
 * @param dataframe RDataFrame which is contains the data used for optimization.
 * @param to_equalize vector of names of the axes which are supposed to be
 * optimized
 * @return Axis configuration with optimized axes.
 */
template <typename AxesConfig, typename DataFrame>
auto EqualizeBinning(AxesConfig &config, DataFrame dataframe,
                     const std::vector<std::string> &to_equalize) {
  Qn::BinningEqualizer<AxesConfig, DataFrame>(dataframe, to_equalize)(config);
}

}  // namespace Qn
#endif  // QNTOOLS_EQUALBINNINGHELPER_H_
