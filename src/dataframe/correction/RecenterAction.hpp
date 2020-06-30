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
#ifndef QNTOOLS_RECENTERACTION_H_
#define QNTOOLS_RECENTERACTION_H_

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "AverageHelper.hpp"
#include "DataContainer.hpp"
#include "QVector.hpp"
#include "TemplateHelpers.hpp"

namespace Qn::Correction {

/**
 * predeclaration of RecenterAction class. Used to make variadic templates from
 * a passed tuple.
 * @tparam AxesConfig
 * @tparam EventParameters
 */
template <typename AxesConfig, typename EventParameters>
class RecenterAction;

/**
 * Recentering class to be used with the AverageHelper class in a DataFrame
 * analysis to implement the recentering procedure. Recentering works using the
 * formula: \f{eqnarray*}{
 * x' &=& x - \langle x \rangle \\
 * y' &=& y - \langle y \rangle
 * \f}
 * where the averages are built in the subevents specified in the
 * Qn::AxesConfiguration. It is possible to enable width equalization.
 * \f{eqnarray*}{
 * x' &=& \frac{x - \langle x \rangle}{\sigma_{x}} \\
 * y' &=& \frac{y - \langle y \rangle}{\sigma{y}}
 * \f}
 * @tparam AxesConfig
 * @tparam EventParameters
 */
template <typename AxesConfig, typename... EventParameters>
class RecenterAction<AxesConfig, std::tuple<EventParameters...>> {
 public:
  /**
   * Constructor
   * @param correction_name name of the correction step
   * @param previous_q_name name of the input Q vector in the input TTree.
   * @param axes_configuration Qn::AxesConfiguration determining the sub samples
   * used for corrections.
   */
  RecenterAction(std::string correction_name, AxesConfig axes_configuration,
                 std::string base_q_name, std::string previous_correction_name)
      : correction_name_(std::move(correction_name)),
        previous_q_name_(std::move(previous_correction_name)),
        base_q_name_(std::move(base_q_name)),
        event_axes_(axes_configuration) {}

  /**
   * Enables the width equalization.
   * @return returns it self.
   */
  RecenterAction EnableWidthEqualization() {
    use_width_equalization_ = true;
    return *this;
  }

  /**
   * Sets the minumum number of contributors required for the correction to be
   * applied in this bin.
   * @return returns it self.
   */
  RecenterAction SetMinimumNumberOfEntries(const unsigned int minimum_entries) {
    min_entries_ = minimum_entries;
    return *this;
  }

  /**
   * Returns the name of the output Q-vector.
   * @return name of the output Q-vector
   */
  [[nodiscard]] std::string GetName() const {
    return AddCorrectionStepToName(correction_name_);
  }

  /**
   * Applies the correction on the input Q-vector and returns a corrected
   * Q-vector. Uses the information saved in a previous iteration over the data.
   * This function is used by the Define functionality of the RDataFrame.
   * @param input_q input to the correction step.
   * @param coordinates event parameters determining the correction histogram
   * bin.
   * @return corrected Q-vector
   */
  Qn::DataContainerQVector operator()(const Qn::DataContainerQVector &input_q,
                                      EventParameters... coordinates) {
    Qn::DataContainerQVector corrected_q(input_q);
    const auto event_bin = event_axes_.GetLinearIndex(coordinates...) * stride_;
    if (event_bin < 0) return corrected_q;
    for (std::size_t ibin = 0; ibin < stride_; ++ibin) {
      const auto correction_bin = event_bin + ibin;
      if (x_[0].At(correction_bin).Entries() < min_entries_) continue;
      auto x_width = 1.;
      auto y_width = 1.;
      corrected_q[ibin].UpdateCorrectionStep(
          Qn::QVector::CorrectionStep::RECENTERED);
      for (std::size_t i_harmonic = 0; i_harmonic < harmonics_vector_.size();
           ++i_harmonic) {
        if (use_width_equalization_) {
          x_width = x_[i_harmonic][correction_bin].Sigma();
          y_width = y_[i_harmonic][correction_bin].Sigma();
        }
        const auto harmonic = harmonics_vector_[i_harmonic];
        corrected_q[ibin].SetQ(harmonic,
                               (input_q[ibin].x(harmonic) -
                                x_[i_harmonic][correction_bin].Mean()) /
                                   x_width,
                               (input_q[ibin].y(harmonic) -
                                y_[i_harmonic][correction_bin].Mean()) /
                                   y_width);
      }
    }
    return corrected_q;
  }

  /**
   * Wrapper function to easily apply the correction on the Q-vector
   * @tparam DataFrame a RDataFrame type
   * @param df the RDataframe which contains the input Q-vector.
   * @return returns a RDataframe which has a column of the corrected Q-vector.
   */
  template <typename DataFrame>
  auto ApplyCorrection(DataFrame &df) const {
    auto columns = GetColumnNames();
    return df.Define(GetName(), *this, columns);
  }

  /**
   * Loads the correction histograms from a previous iteration from the file.
   * @param file file which contains the correction histograms.
   * @param reader reader which wraps the input tree. Needed to perform the
   * initialization.
   */
  bool LoadCorrectionFromFile(TDirectory *dir, TTreeReader &reader) {
    if (!dir) {
      std::cout << "Rerunning the correction step." << std::endl;
      return false;
    }
    Initialize(reader);
    if (!dir->FindKey(GetNameIncludingAxesAndHarmonics().data())) {
      std::cout << "correction " << GetName() << ": Not found in the file "
                << dir->GetName() << ". ";
      std::cout << "Rerunning the correction step." << std::endl;
      Reset();
      return false;
    }
    for (std::size_t i_harmonic = 0; i_harmonic < harmonics_vector_.size();
         ++i_harmonic) {
      auto harmonic = harmonics_vector_[i_harmonic];
      auto x = dynamic_cast<Qn::DataContainerStatistic *>(
          dir->Get((GetNameIncludingAxesAndHarmonics() + "/X_" +
                    std::to_string(harmonic))
                       .data()));
      auto y = dynamic_cast<Qn::DataContainerStatistic *>(
          dir->Get((GetNameIncludingAxesAndHarmonics() + "/Y_" +
                    std::to_string(harmonic))
                       .data()));
      if (!(x && y)) {
        std::cout << "correction " << GetName() << ": Harmonic" << harmonic;
        std::cout << " not found in the file " << dir->GetName() << ". ";
        std::cout << "Rerunning the correction step." << std::endl;
        Reset();
        return false;
      }
      auto read_axes = x->GetAxes();
      auto configured_axes = x_.at(i_harmonic).GetAxes();
      std::vector<bool> result;
      std::transform(
          std::begin(read_axes), std::end(read_axes),
          std::begin(configured_axes), std::back_inserter(result),
          [](const Qn::AxisD &a, const Qn::AxisD &b) { return a == b; });
      if (!std::all_of(result.begin(), result.end(),
                       [](bool x) { return x; })) {
        std::cout << "correction " << GetName()
                  << ": Axes not equal to the ones found in the file ";
        std::cout << dir->GetName() << ". ";
        std::cout << "Rerunning the correction step." << std::endl;
        Reset();
        return false;
      }
      x_.at(i_harmonic) = *x;
      y_.at(i_harmonic) = *y;
    }
    std::cout << "correction " << GetName() << ": Found in the file "
              << dir->GetName() << "." << std::endl;
    return true;
  }

  /**
   * Writes the correction histograms to file.
   * @param directory directory in which the correction histograms are supposed
   * to be written.
   */
  void Write(TDirectory *directory) {
    using namespace std::literals::string_literals;
    directory->cd();
    auto unique_name = GetNameIncludingAxesAndHarmonics();
    directory->mkdir(unique_name.data());
    directory->cd(unique_name.data());
    for (std::size_t i_harmonic = 0; i_harmonic < harmonics_vector_.size();
         ++i_harmonic) {
      x_.at(i_harmonic)
          .Write(
              ("X_"s + std::to_string(harmonics_vector_[i_harmonic])).data());
      y_.at(i_harmonic)
          .Write(
              ("Y_"s + std::to_string(harmonics_vector_[i_harmonic])).data());
    }
    CreateBinOccupancyHisto().Write("EntriesPerBin");
    directory->cd("../");
  }

 private:
  friend class AverageHelper<RecenterAction>;  /// Helper friend
  bool use_width_equalization_ =
      false;  /// Switch for applying the width equalization procedure.
  unsigned int min_entries_ =
      10;  /// Number of minimum entries in a bin required to apply corrections.
  unsigned int stride_ = 1.;  /// stride of the differential Q-vector.
  std::vector<unsigned int> harmonics_vector_;  /// vector of enabled harmonics.
  std::string correction_name_;                 /// name of the correction step.
  std::string previous_q_name_;                 /// name of the input Q-vector.
  std::string
      base_q_name_;  /// name of the input Q-vector used for Initialization.
  AxesConfig event_axes_;  /// event axes used to classify the events in classes
                           /// for the correction step.
  std::vector<Qn::DataContainerStatistic>
      x_;  /// x component of correction histograms.
  std::vector<Qn::DataContainerStatistic>
      y_;  /// y component correction histograms.

  /**
   * Initializes the correction step using the information inside the input
   * tree.
   * @param reader reader wrapping the input Q-vector tree. This function is
   * required by the AverageHelper.
   */
  void Initialize(TTreeReader &reader) {
    using namespace std::literals::string_literals;
    reader.Restart();
    TTreeReaderValue<DataContainerQVector> input_data(reader,
                                                      base_q_name_.data());
    reader.Next();
    if (input_data.GetSetupStatus() < 0) {
      throw std::runtime_error(
          "The Q-Vector entry "s + input_data.GetBranchName() +
          " in the tree is not valid. Cannot setup the recentering");
    }
    auto input_q = input_data->At(0);
    input_q.InitializeHarmonics();
    int i_harmonic = 0;
    int harmonic = input_q.GetFirstHarmonic();
    Reset();
    while (harmonic != -1) {
      x_.emplace_back();
      y_.emplace_back();
      x_[i_harmonic].AddAxes(event_axes_.GetVector());
      y_[i_harmonic].AddAxes(event_axes_.GetVector());
      if (!input_data->IsIntegrated()) {
        x_[i_harmonic].AddAxes(input_data->GetAxes());
        y_[i_harmonic].AddAxes(input_data->GetAxes());
      }
      harmonics_vector_.push_back(harmonic);
      harmonic = input_q.GetNextHarmonic(harmonic);
      ++i_harmonic;
    }
    stride_ = input_data->size();
  }

  /**
   * Calculates the corrections and saves them in the x_ and y_ member
   * variables. This function is required by the AverageHelper.
   * @param input input Q-vector, which is supposed to be corrected.
   * @param event_parameters event parameters determining the correction
   * histogram bin.
   */
  void CalculateAction(const Qn::DataContainerQVector &input,
                       EventParameters... event_parameters) {
    const auto event_bin =
        event_axes_.GetLinearIndex(event_parameters...) * stride_;
    if (event_bin < 0) return;
    for (std::size_t i = 0; i < harmonics_vector_.size(); ++i) {
      for (std::size_t ibin = 0; ibin < stride_; ++ibin) {
        const auto output_bin = event_bin + ibin;
        if (input[ibin].sumweights() > 0.) {
          x_[i].At(output_bin).Fill(input[ibin].x(harmonics_vector_[i]), 1.);
          y_[i].At(output_bin).Fill(input[ibin].y(harmonics_vector_[i]), 1.);
        }
      }
    }
  }

  /**
   * Merges the correction histogram after the collection of statistics is
   * complete. This function is required by the AverageHelper class.
   * @param results the other results which are to be merged with this one.
   */
  void Merge(const std::vector<std::shared_ptr<RecenterAction>> &results) {
    for (const auto &result : results) {
      for (std::size_t i_harmonic = 0; i_harmonic < harmonics_vector_.size();
           ++i_harmonic) {
        TList xs;
        TList ys;
        xs.Add(&result->x_[i_harmonic]);
        ys.Add(&result->y_[i_harmonic]);
        x_[i_harmonic].Merge(&xs);
        y_[i_harmonic].Merge(&ys);
      }
    }
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
  auto BookMe(DataFrame df, Qn::AverageHelper<RecenterAction> helper) {
    return df.template Book<DataContainerQVector, EventParameters...>(
        std::move(helper), GetColumnNames());
  }

  /**
   * Copy one state to the next.
   * Used during multi threading.
   * @param other another state.
   */
  void CopyInitializedState(const RecenterAction &other) {
    stride_ = other.stride_;
    harmonics_vector_ = other.harmonics_vector_;
    event_axes_ = other.event_axes_;
    x_ = other.x_;
    y_ = other.y_;
  }

  [[nodiscard]] std::string AddCorrectionStepToName(
      std::string correction_step) const {
    if (correction_step.empty())
      return base_q_name_;
    else
      return base_q_name_ + "_" + correction_step;
  }

  /**
   * Returns the name of the columns used in the correction step. This includes
   * both the input Q-vector and the name of the axes parameters. This function
   * is required by the AverageHelper.
   * @return returns a vector of the column names.
   */
  [[nodiscard]] std::vector<std::string> GetColumnNames() const {
    std::vector<std::string> columns;
    columns.emplace_back(AddCorrectionStepToName(previous_q_name_));
    const auto event_axes_names = event_axes_.GetNames();
    std::copy(std::begin(event_axes_names), std::end(event_axes_names),
              std::back_inserter(columns));
    return columns;
  }

  /**
   * Resets the correction.
   */
  void Reset() {
    x_.clear();
    y_.clear();
    harmonics_vector_.clear();
    stride_ = 1;
  }

  /**
   * Gives a unique name to the Corrections step for persisting it to the file.
   */
  [[nodiscard]] std::string GetNameIncludingAxesAndHarmonics() const {
    std::stringstream name;
    name << GetName();
    name << "_h";
    for (const auto &harmonic : harmonics_vector_) name << harmonic;
    for (const auto &axis : event_axes_.GetVector()) name << axis.ShortName();
    return name.str();
  }

  TH1D CreateBinOccupancyHisto() {
    auto min_max = std::minmax_element(
        x_[0].begin(), x_[0].end(), [](const Statistic &a, const Statistic &b) {
          return a.Entries() < b.Entries();
        });
    auto min_entries = min_max.first->Entries();
    auto max_entries = min_max.second->Entries();
    auto difference = (max_entries - min_entries) * 0.05;
    max_entries = max_entries + difference;
    min_entries = min_entries - difference;
    if (min_entries < 1.) min_entries = 1.;
    TH1D histo_bin_occupancy("bin_occupancy",
                             "occupancy per bin; Counts; Entries in a bin", 100,
                             min_entries, max_entries);
    for (auto &bin : x_[0]) {
      auto n = bin.Entries();
      histo_bin_occupancy.Fill(n);
    }
    return histo_bin_occupancy;
  }
};

/**
 * Function which creates the Recentering correction step
 * @tparam EventAxes Qn::AxesConfiguration
 * @param correction_name name of the correction step
 * @param sub_event_name name of the input Q vector in the input TTree.
 * @param axes_configuration Qn::AxesConfiguration determining the sub samples
 * used for corrections.
 * @return RecenterAction
 */
template <typename EventAxes>
auto Recentering(const std::string &correction_name,
                 EventAxes axes_configuration, const std::string &base_q_name,
                 const std::string &previous_correction = "") {
  return RecenterAction<EventAxes, typename EventAxes::AxisValueTypeTuple>{
      correction_name, axes_configuration, base_q_name, previous_correction};
}

/**
 * Convenience function to apply the corrections of multiple Recentering
 * procedures. Ends the recursion.
 * @tparam DataFrame type of a RDataFrame
 * @tparam Last Type of the last remaining correction
 * @param df RDataFrame which contains the input Q-vectors.
 * @param last The last remaining correction step.
 * @return RDataFrame, which contains the the corrected Q-vectors.
 */
template <typename DataFrame, typename Last>
auto ApplyCorrections(DataFrame df, Last last) {
  return last.ApplyCorrection(df);
}

/**
 * Convenience function to apply the corrections of multiple Recentering
 * procedures. Ends the recursion.
 * @tparam DataFrame type of a RDataFrame
 * @tparam First Type of the first correction step.
 * @tparam Rest Types of the remaining correction steps
 * @param df RDataFrame which contains the input Q-vectors.
 * @param first the first correction step to be applied.
 * @param rest rest of the correction steps.
 * @return RDataFrame, which contains the the corrected Q-vectors.
 */
template <typename DataFrame, typename First, typename... Rest>
inline auto ApplyCorrections(DataFrame df, First first, Rest... rest) {
  return ApplyCorrections(first.ApplyCorrection(df), rest...);
}

/**
 * Convenience function to apply the corrections of multiple Recentering
 * procedures.
 * @tparam DataFrame type of a RDataFrame
 * @param df RDataFrame which contains the input Q-vectors.
 * @param first the first correction step to be applied.
 * @param rest rest of the correction steps.
 * @return RDataFrame, which contains the the corrected Q-vectors.
 */
template <typename DataFrame, typename VectorOfCorrections>
inline auto ApplyCorrectionsVector(DataFrame df,
                                   VectorOfCorrections corrections) {
  ROOT::RDF::RNode corrected_df(df);
  for (auto &result : corrections) {
    corrected_df =
        Qn::TemplateHelpers::DereferenceRResultPtr(result).ApplyCorrection(
            corrected_df);
  }
  return corrected_df;
}

}  // namespace Qn::Correction
#endif  // QNTOOLS_RECENTERACTION_H_
