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

#ifndef QNTOOLS_QVECTORACTION_HPP_
#define QNTOOLS_QVECTORACTION_HPP_

#include <string>
#include <string_view>
#include <tuple>

#include <ROOT/RVec.hxx>

#include "AxesConfiguration.hpp"
#include "Axis.hpp"
#include "DataContainer.hpp"

namespace Qn {
/**
 * predeclaration of DefineQAction class. Used to make variadic templates from a
 * passed tuple.
 * @tparam AxesConfig
 * @tparam EventParameters
 */
template <typename AxesConfig, typename DetectorParameters>
class QVectorHelper;

/**
 * Defines a Q vector with a variable number of harmonics from a vector of
 * azimuthal angles and a vector of weights Subevents can be defined with an
 * AxesConfiguration. For Integrated Q vectors
 * @tparam AxesConfig
 * @tparam DetectorParameters
 */
template <typename AxesConfig, typename... DetectorParameters>
class QVectorHelper<AxesConfig, std::tuple<DetectorParameters...>> {
 public:
  using AxisValueType = typename AxesConfig::AxisValueType;
  using AxisValueTypeTuple = typename AxesConfig::AxisValueTypeTuple;
  using BinEdgesMapping =
      typename std::vector<std::pair<AxisValueTypeTuple, AxisValueTypeTuple>>;
  static constexpr std::size_t NumberOfParameters =
      sizeof...(DetectorParameters);
  QVectorHelper(std::string_view name, const std::vector<int> &harmonics,
                AxesConfig axes_config)
      : name_(name), axes_(axes_config), q_vector_(axes_config.GetVector()) {
    bin_edges_mapping_ = axes_.GetBinEdgesIndexMap();
    for (auto &q : q_vector_) {
      for (auto &harmonic : harmonics) {
        q.ActivateHarmonic(harmonic);
      }
      q.InitializeHarmonics();
    }
  }
  /**
   * Returns the name of the output Q-vector.
   * @return name of the output Q-vector
   */
  [[nodiscard]] std::string GetName() const { return name_; }

  Qn::DataContainerQVector operator()(const ROOT::RVec<double> &phis,
                                      const ROOT::RVec<double> &weights,
                                      DetectorParameters... parameters) {
    Qn::DataContainerQVector ret = q_vector_;
    auto rest_phis = phis;
    auto rest_weights = weights;
    auto rest_parameters_tuple = std::tuple(parameters...);
    for (int ibin = 0; ibin < axes_.GetSize(); ++ibin) {
      auto [lower_edges, upper_edges] = bin_edges_mapping_[ibin];
      Qn::TemplateHelpers::TupleOf<NumberOfParameters, ROOT::RVec<double>>
          selections_tuple;
      auto selector = [](auto &res, auto lower, auto upper, auto value) {
        res = value >= lower && value < upper;
      };
      Qn::TemplateHelpers::TupleForEach(selector, selections_tuple, lower_edges,
                                        upper_edges, rest_parameters_tuple);
      auto is_selected = std::apply([&](auto... cond) { return (cond && ...); },
                                    selections_tuple);
      auto selected_phis = rest_phis[is_selected];
      auto selected_weights = rest_weights[is_selected];
      auto &bin = ret.At(ibin);
      for (int i = 0; i < selected_phis.size(); ++i) {
        bin.Add(selected_phis[i], selected_weights[i]);
      }
      bin.CheckQuality();
      rest_phis = rest_phis[!is_selected];
      rest_weights = rest_weights[!is_selected];
      rest_parameters_tuple = std::apply(
          [&is_selected](auto... t) { return std::tuple(t[!is_selected]...); },
          rest_parameters_tuple);
    }
    return ret;
  }

  Qn::DataContainerQVector* GetPrototype() {return &q_vector_;}

 private:
  std::string name_;
  AxesConfig axes_;
  BinEdgesMapping bin_edges_mapping_;
  Qn::DataContainerQVector q_vector_;
};

/**
 * Template specialization for integrated Q-Vectors.
 * AxesConfiguration<AxisD> is a dummy type and will not be used in the
 * calculation The empty tuple will be ignored.
 */
template <>
inline Qn::DataContainerQVector
QVectorHelper<AxesConfiguration<AxisD>, typename std::tuple<>>::operator()(
    const ROOT::RVec<double> &phis, const ROOT::RVec<double> &weights) {
  Qn::DataContainerQVector ret = q_vector_;
  auto &outputbin = ret.At(0);
  for (int i = 0; i < phis.size(); ++i) {
    outputbin.Add(phis[i], weights[i]);
  }
  outputbin.CheckQuality();
  return ret;
}

template <typename ParameterAxes>
auto MakeQVectorHelper(std::string_view name, const std::vector<int> &harmonics,
                       ParameterAxes axes_config) {
  using Parameters = typename Qn::TemplateHelpers::TupleOf<
      ParameterAxes::kDimension,
      ROOT::RVec<typename ParameterAxes::AxisValueType>>;
  return QVectorHelper<ParameterAxes, Parameters>(name, harmonics, axes_config);
}

inline auto MakeQVectorHelper(std::string_view name,
                              const std::vector<int> &harmonics) {
  auto axes_null = MakeAxes(AxisD("integrated", 1, 0, 1));
  return QVectorHelper<decltype(axes_null), std::tuple<>>(name, harmonics,
                                                          axes_null);
}
}  // namespace Qn

#endif  // QNTOOLS_QVECTORACTION_HPP+