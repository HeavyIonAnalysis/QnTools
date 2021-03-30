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
#ifndef QNTOOLS_CORRELATIONACTION_H
#define QNTOOLS_CORRELATIONACTION_H

#include <algorithm>
#include <array>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <ROOT/RVec.hxx>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "AverageHelper.hpp"
#include "DataContainer.hpp"
#include "QVector.hpp"
#include "Stat.hpp"
#include "StatCollect.hpp"
#include "TemplateFunctions.hpp"

namespace Qn::Correlation {

/**
 * CorrelationActionBase
 */
class CorrelationActionBase {
 public:
  using InitializationObject = std::vector<Qn::DataContainerQVector>;
  /**
   * Constructor
   * @param name name of the correlation
   */
  explicit CorrelationActionBase(std::string_view name) : action_name_(name) {}
  /**
   * Returns the name of the output Q-vector.
   * @return name of the output Q-vector
   */
  [[nodiscard]] std::string GetName() const { return action_name_; }

  /**
   * Writes the result of the correlation to file. Triggers the evaluation of
   * the event loop
   */
  void Write() const { correlation_.Write(GetName().data()); }

  /**
   * Gives a const reference to the result of the correlation. Triggers the
   * evaluation of the event loop.
   * @return returns const reference.
   */
  [[nodiscard]] const Qn::DataContainerStatCollect &GetDataContainer() const {
    return correlation_;
  }

 protected:
  std::string action_name_;             /// Name of the action
  Qn::DataContainerStatCollect correlation_;  /// Result data container.
  std::vector<std::string> input_from_reader_names_; /// names of branches used for initialization
};

/**
 * predeclaration of CorrelationAction class. Used to make variadic templates
 * from a passed tuple.
 * @tparam InputDataContainers
 * @tparam AxesConfig
 * @tparam EventParameters
 */
template <typename Function, typename WeightFunction,
          typename InputDataContainers, typename AxesConfig,
          typename EventParameters>
class CorrelationAction;

/**
 * CorrelationAction used to calculate correlations in one event.
 * @tparam Function
 * @tparam InputDataContainers
 * @tparam AxesConfig
 * @tparam EventParameters
 */
template <typename Function, typename WeightFunction,
          typename... InputDataContainers, typename AxesConfig,
          typename... EventParameters>
class CorrelationAction<Function, WeightFunction,
                        std::tuple<InputDataContainers...>, AxesConfig,
                        std::tuple<EventParameters...>>
    : public CorrelationActionBase {
 private:
  constexpr static std::size_t NumberOfInputs = sizeof...(InputDataContainers);
  using DataContainerRef = std::reference_wrapper<const DataContainerQVector>;
 public:
  using InitializationObject = CorrelationActionBase::InitializationObject;
  /**
   * Constructor
   * @param function correlation function of type (const Qn::QVector &q...) ->
   * double
   * @param input_names Array of Q-vectors going into the correlation
   * @param weights Array of weights of the Q-vectors
   * @param correlation_name Name of the correlation
   * @param event_axes Event axes for the binning of the correlation.
   * @param n_samples Number of samples used for the bootstrapping
   */
  CorrelationAction(std::string_view correlation_name, Function function,
                    WeightFunction weight_function, bool use_weights,
                    const std::array<std::string, NumberOfInputs> &input_names,
                    AxesConfig event_axes, unsigned int n_samples)
      : CorrelationActionBase(correlation_name),
        n_samples_(n_samples),
        input_names_(input_names),
        event_axes_(event_axes),
        function_(function),
        weight_function_(weight_function),
        use_weights_(use_weights) {}

  template <std::size_t N>
  auto SetReaderInputNames(const std::array<std::string, N> &names) {
    input_from_reader_names_ = names;
  }

 private:
  friend class AverageHelper<CorrelationAction>;  /// Helper friend

  unsigned int stride_ = 1;     /// Offset of due to the non-event axes.
  unsigned int n_samples_ = 1;  /// Number of samples used in the ReSampler.
  std::array<std::string, NumberOfInputs>
      input_names_;                 /// Names of the input Q-vectors.
 std::array<std::string, NumberOfInputs>
      input_from_reader_names_;     /// Names of the input Q-vectors.
  AxesConfig event_axes_;           /// Configuration of the event axes.
  Function function_;               /// correlation function.
  WeightFunction weight_function_;  /// weight function.
  bool use_weights_;                /// determines if weight is used.

  /**
   * Returns the name of the columns used in the correction step. This includes
   * both the input Q-vector and the name of the axes parameters. This function
   * is required by the AverageHelper.
   * @return returns a vector of the column names.
   */
  [[nodiscard]] std::vector<std::string> GetColumnNames() const {
    auto columns = std::vector<std::string>{};
    std::copy(std::begin(input_names_), std::end(input_names_),
              std::back_inserter(columns));
    columns.emplace_back("samples");
    const auto event_axes_names = event_axes_.GetNames();
    std::copy(std::begin(event_axes_names), std::end(event_axes_names),
              std::back_inserter(columns));
    return columns;
  }

  /**
   * Checks if at least one Q-vector is not a reference Q-vector.
   * @return true if at least one of them is not a reference Q-vector.
   */
  [[nodiscard]] inline bool IsObservable() const { return use_weights_; }

  /**
   * Initializes the CorrelationAction using the initialization object.
   * @param object The TTreeReader gives access to the input TTree.
   */
  void Initialize(const InitializationObject &object) {
    correlation_.AddAxes(event_axes_.GetVector());
    for (std::size_t i = 0; i < object.size(); ++i) {
      const auto &i_data = object[i];
      if (!i_data.IsIntegrated()) AddAxes(object, i);
    }
    stride_ = correlation_.size() / event_axes_.GetSize();
    for (auto &bin : correlation_) {
      bin.SetNumberOfSamples(n_samples_);
      if (IsObservable())
        bin.SetWeightType(Qn::Stat::WeightType::OBSERVABLE);
      else
        bin.SetWeightType(Qn::Stat::WeightType::REFERENCE);
    }
  }

  void TryNextEventInReader(TTreeReader &reader,
                            std::vector<TTreeReaderValue<DataContainerQVector>> input_data,
                            Long64_t i_event,
                            Long64_t n_events) {
    using namespace std::literals::string_literals;
    if (i_event > n_events) {
      throw std::out_of_range("ERROR: Tried 10% of the events."
                              "Could not find a good event to initialize."
                              "Try to initialize using the InitializationObject.");
    }
    try {
      reader.Next();
      auto initialization_object = std::vector<Qn::DataContainerQVector>{};
      // Move valid entries to the initialization Object
      std::transform(std::begin(input_data), std::end(input_data),
                     std::back_inserter(initialization_object),
                     [](auto &reader_value) {
                       if (reader_value.GetSetupStatus() < 0)
                         throw std::runtime_error("Q-Vector branch "s +
                             reader_value.GetBranchName() +
                             " not found.");
                       else
                         return *reader_value;
                     });
      Initialize(initialization_object);
    } catch (std::out_of_range &) {
      TryNextEventInReader(reader, input_data, i_event+1, n_events);
    }
  }

  /**
   * Initializes the CorrelationAction using the input Q-vectors in the input
   * TTree. Throws a runtime error, when a entry of the input_names is not found
   * in the TTree.
   * @param reader The TTreeReader gives access to the input TTree.
   */
  void Initialize(TTreeReader &reader) {
    using namespace std::literals::string_literals;
    reader.Restart();
    auto input_data = std::vector<TTreeReaderValue<DataContainerQVector>>{};
    // Get a vector of TTreeReaderValues.
    if (input_from_reader_names_.empty()) input_from_reader_names_ = input_names_;
    std::transform(
        std::begin(input_from_reader_names_), std::end(input_from_reader_names_),
        std::back_inserter(input_data), [&reader](const std::string &name) {
          return TTreeReaderValue<DataContainerQVector>{reader, name.data()};
        });
    // Read in the first event in the TTree
    reader.Next();
    auto initialization_object = std::vector<Qn::DataContainerQVector>{};
    // Move valid entries to the initialization Object
    std::transform(std::begin(input_data), std::end(input_data),
                   std::back_inserter(initialization_object),
                   [](auto &reader_value) {
                     if (reader_value.GetSetupStatus() < 0)
                       throw std::runtime_error("Q-Vector branch "s +
                                                reader_value.GetBranchName() +
                                                " not found.");
                     else
                       return *reader_value;
                   });
    // Initialize
    try {
      Initialize(initialization_object);
    } catch (std::out_of_range&) {
      TryNextEventInReader(reader, input_data, 0, 0.1*reader.GetEntries());
    }
  }

  /**
   * Calculates the result of the correlation function.
   * @param input_q input Q-vectors
   * @param sample_ids frequency map of the resamples.
   * @param event_parameters event parameters for event classification.
   */
  void CalculateAction(InputDataContainers... input_q,
                       const ROOT::RVec<ULong64_t> &sample_ids,
                       EventParameters... event_parameters) {
    long initial_offset =
        event_axes_.GetLinearIndex(event_parameters...) * stride_;
    if (initial_offset < 0) return;
    std::array<const Qn::QVector *, NumberOfInputs> q_vectors;
    const std::array<const DataContainerRef, NumberOfInputs> input_array = {
        std::cref(input_q)...};
    LoopOverBins(initial_offset, q_vectors, input_array, sample_ids, 0);
  }

  /**
   * Merges the correction histogram after the collection of statistics is
   * complete. This function is required by the AverageHelper class.
   * @param results the other results which are to be merged with this one.
   */
  void Merge(const std::vector<std::shared_ptr<CorrelationAction>> &results) {
    TList correlation_list;
    for (auto &result : results) {
      correlation_list.Add(&result->correlation_);
    }
    correlation_.Merge(&correlation_list);
  }

  /**
   * Called during booking process.
   * Allows to hide template parameters.
   * @tparam DataFrame type of dataframe
   * @param df dataframe to which the action is booked.
   * @param helper helper
   * @return Result pointer to the result of the action
   */
  template <typename DataFrame>
  auto BookMe(DataFrame df, Qn::AverageHelper<CorrelationAction> helper) {
    return df.template Book<InputDataContainers..., ROOT::RVec<ULong64_t>,
                            EventParameters...>(std::move(helper),
                                                GetColumnNames());
  }

  /**
   * Copy one state to the next.
   * Used during multi threading.
   * Does not copy function_.
   * @param other another state.
   */
  void CopyInitializedState(CorrelationAction &other) {
    stride_ = other.stride_;
    n_samples_ = other.n_samples_;
    action_name_ = other.action_name_;
    input_names_ = other.input_names_;
    use_weights_ = other.use_weights_;
    event_axes_ = other.event_axes_;
    correlation_ = other.correlation_;
  }

  /**
   * Loops recursivly over all bins of the correlation result.
   * @param out_bin Initial offset based on the current event parameters.
   * @param q_array Holder for the Q-vectors entering the calculation of the
   * current bin.
   * @param input_array Holder of the DataContainers for the current event.
   * @param sample_ids Frequency map of the resamples of the current event.
   * @param step Recursion step.
   */
  void LoopOverBins(
      long &out_bin, std::array<const Qn::QVector *, NumberOfInputs> &q_array,
      const std::array<const DataContainerRef, NumberOfInputs> &input_array,
      const ROOT::RVec<ULong64_t> &sample_ids, std::size_t step) {
    if (step + 1 == NumberOfInputs) {  /// base case
      for (const auto &bin : input_array[step].get()) {
        if (bin.n() < 1.) {
          ++out_bin;
          continue;
        }
        q_array[step] = &bin;
        correlation_[out_bin].Fill(
            TemplateFunctions::Call(function_, q_array),
            TemplateFunctions::Call(weight_function_, q_array), sample_ids);
        ++out_bin;
      }
    } else {  /// recursion
      for (const auto &bin : input_array[step].get()) {
        if (bin.n() < 1.) {
          ++out_bin;
          continue;
        }
        q_array[step] = &bin;
        LoopOverBins(out_bin, q_array, input_array, sample_ids, step + 1);
      }
    }
  }

  /**
   * Adds axes of the Q-vectors to the correlation DataContainer.
   * If similar axes are found the axis is renamed to make sure the name is
   * unique.
   * @param data_containers Datacontainers going into the correlation
   * @param i Current position of the data container.
   */
  void AddAxes(const std::vector<DataContainerQVector> &data_containers,
               std::size_t i) {
    for (auto axis : data_containers[i].GetAxes()) {
      std::string name = axis.Name();
      for (std::size_t j = 0; j < NumberOfInputs; ++j) {
        if (i == j) continue;
        auto &other = data_containers[j];
        if (other.IsIntegrated()) continue;
        for (const auto &other_axis : other.GetAxes()) {
          if (axis != other_axis) continue;
          auto this_name = std::string(input_names_[i]);
          if (this_name == input_names_[j])
            name = this_name + "_" + std::to_string(i) + "_" + axis.Name();
          else
            name = this_name + "_" + axis.Name();
        }
      }
      axis.SetName(name);
      correlation_.AddAxis(axis);
    }
  }
};

using namespace TemplateFunctions;
enum class UseWeights {
  Yes,
  No
};
/**
 * Helper function to create a correlation without specifying the template
 * parameters directly.
 * @tparam Function type of the correlation function.
 * @tparam AxesConfig Type of the axis configuration.
 * @param function correlation function
 * @param input_names Array of input Q-vector names.
 * @param weights Array of weights.
 * @param correlation_name Name of the correlation
 * @param event_axes Event axes used in the correlation.
 * @param n_samples Number of samples used in the Resampling step.
 * @return result of the correlation.
 */
template <typename Function, typename WeightFunction, typename AxesConfig>
auto MakeCorrelationAction(
    std::string_view correlation_name, Function function,
    WeightFunction weight_function, UseWeights use_weights,
    const std::array<std::string, FunctionTraits<Function>::Arity> input_names,
    AxesConfig event_axes, unsigned int n_samples) {
  using DataContainerTuple =
      TupleOf<FunctionTraits<Function>::Arity, Qn::DataContainerQVector>;
  using EventParameterTuple = typename AxesConfig::AxisValueTypeTuple;
  bool usew = use_weights==UseWeights::Yes ? true : false;
  return CorrelationAction<Function, WeightFunction, DataContainerTuple,
                           AxesConfig, EventParameterTuple>(
      correlation_name, function, weight_function, usew, input_names,
      event_axes, n_samples);
}

}  // namespace Qn::Correlation
#endif  // QNTOOLS_CORRELATIONHELPER_H_
