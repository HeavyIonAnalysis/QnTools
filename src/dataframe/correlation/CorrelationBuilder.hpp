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

#ifndef QNTOOLS_CORRELATIONVECTOR_HPP_
#define QNTOOLS_CORRELATIONVECTOR_HPP_

#include <array>
#include <string>
#include <vector>

#include <ROOT/RResultPtr.hxx>

#include "AverageHelper.hpp"
#include "AxesConfiguration.hpp"
#include "CorrelationAction.hpp"
#include "TTreeReader.h"

namespace Qn::Correlation {
template <typename EventAxes, typename DataFrame>
class CorrelationBuilder {
 public:
  using ResultPtr =
      ROOT::RDF::RResultPtr<Qn::Correlation::CorrelationActionBase>;
  using InitializationObjects = std::vector<Qn::DataContainerQVector>;

  /**
   * Add correlations to the correlation manager
   * @param df Dataframe
   * @param n_samples number of samples
   * @param axes correlation axes
   */
  CorrelationBuilder(DataFrame *df, int n_samples, EventAxes axes)
      : df_(df), n_samples_(n_samples), event_axes_(axes) {}

  /**
   * Add external TTreeReader needs to be called before Adding correlations
   * using the external TTreeReader
   * @param reader external TTreeReader pointer
   */
  void AddExternalReader(TTreeReader *reader) { reader_ = reader; }

  /**
   * Add Correlations using the internal TTreeReader
   * @tparam N
   * @tparam Function
   * @tparam WF
   * @param name
   * @param function
   * @param weight_function
   * @param use_weights
   * @param input_names
   * @param input_reader_names
   * @return
   */
  template <std::size_t N, typename Function, typename WF>
  auto AddCorrelationWithInternalReader(
      const std::string &name, Function &&function, WF &&weight_function,
      Qn::Correlation::UseWeights use_weights,
      std::array<std::string, N> input_names,
      std::array<std::string, N> input_reader_names) {
    auto correlation_name = MakeName(name);
    auto correlation = Qn::Correlation::MakeCorrelationAction(
        correlation_name, function, weight_function, use_weights, input_names,
        event_axes_, n_samples_);
    correlation.SetReaderInputNames(input_reader_names);
    result_ptrs_.push_back(Qn::MakeAverageHelper(correlation).BookMe(*df_));
  }

  /**
   * Add Correlations using the initialization object
   * @tparam N
   * @tparam Function
   * @tparam WF
   * @param name
   * @param function
   * @param weight_function
   * @param use_weights
   * @param input_names
   * @param initializations
   * @return
   */
  template <std::size_t N, typename Function, typename WF>
  auto AddCorrelationWithInitializationObject(
      const std::string &name, Function &&function, WF &&weight_function,
      Qn::Correlation::UseWeights use_weights,
      std::array<std::string, N> input_names,
      InitializationObjects initializations) {
    auto correlation_name = MakeName(name);
    auto correlation = Qn::Correlation::MakeCorrelationAction(
        correlation_name, function, weight_function, use_weights, input_names,
        event_axes_, n_samples_);
    initializations_.emplace_back(
        std::make_unique<InitializationObjects>(initializations));
    result_ptrs_.push_back(Qn::MakeAverageHelper(correlation)
                               .SetInitializationWithInitializationObject(
                                   initializations_.back().get())
                               .BookMe(*df_));
  }

  /**
   * Add Correlations using the external TTreeReader
   * @tparam N
   * @tparam Function
   * @tparam WF
   * @param name
   * @param function
   * @param weight_function
   * @param use_weights
   * @param input_names
   * @param initializations
   * @return
   */
  template <std::size_t N, typename Function, typename WF>
  auto AddCorrelationWithExternalReader(
      const std::string &name, Function &&function, WF &&weight_function,
      Qn::Correlation::UseWeights use_weights,
      const std::array<std::string, N> &input_names,
      const std::array<std::string, N> &input_reader_names) {
    if (!reader_)
      throw std::runtime_error(
          "External TTreeReader not specified before calling"
          "AddCorrelationWithExternalReader for " +
          name);
    auto correlation_name = MakeName(name);
    auto correlation = Qn::Correlation::MakeCorrelationAction(
        correlation_name, function, weight_function, use_weights, input_names,
        event_axes_, n_samples_);
    correlation.SetReaderInputNames(input_reader_names);
    result_ptrs_.push_back(
        Qn::MakeAverageHelper(correlation)
            .SetInitializationWithExternalTTreeReader(reader_)
            .BookMe(*df_));
  }

  /**
   * Get the results. Results are not yet evaluated.
   * @return std::vector<ROOT::RResultPtr>
   */
  auto GetResults() { return result_ptrs_; }

 private:
  TTreeReader *reader_ =
      nullptr;            /// pointer to external TTreeReader not owning
  DataFrame *df_;         /// pointer to the dataframe not owning
  EventAxes event_axes_;  /// correlation axes
  std::vector<ResultPtr> result_ptrs_;  /// vector of result pointers
  std::vector<std::unique_ptr<InitializationObjects>>
      initializations_;  /// all initialization objects to keep them in scope.
  int n_samples_;  /// number of bootstrap samples needs to be equal to the
                   /// preconfigured samples in the dataframe

  auto MakeName(const std::string &name) {
    auto correlation = std::string{name};
    auto axesv = event_axes_.GetVector();
    for (const auto &axis : axesv) correlation += axis.Name();
    return correlation;
  }
};

}  // namespace Qn::Correlation
#endif
