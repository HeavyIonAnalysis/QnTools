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

#ifndef FLOW_CORRECTIONMANAGER_H
#define FLOW_CORRECTIONMANAGER_H

#include <string>
#include <map>
#include <utility>

#include <TFile.h>

#include "Detector.hpp"
#include "InputVariableManager.hpp"
#include "Cuts.hpp"
#include "CorrectionProfile3DCorrelations.hpp"
#include "CorrectionProfileCorrelationComponents.hpp"
#include "SubEventChannels.hpp"
#include "SubEvent.hpp"
#include "SubEventTracks.hpp"
#include "Recentering.hpp"
#include "TwistAndRescale.hpp"

#include "GainEqualization.hpp"
#include "Alignment.hpp"
#include "DataContainer.hpp"
#include "RunList.hpp"
#include "DetectorList.hpp"

namespace Qn {
class CorrectionManager {
  using CutCallBack =  std::function<std::unique_ptr<ICut>(Qn::InputVariableManager *)>;
 public:
  CorrectionManager() = default;
  virtual ~CorrectionManager() = default;
  /**
   * Add a variable to the variable manager
   * @param name Name of the variable
   * @param id Id of the variable inside the array used to pass the data into the framework.
   * @param length Lenght of the variable inside the array e.g. number of channels.
   */
  void AddVariable(const std::string &name, const int id, const int length) {
    variable_manager_.CreateVariable(name, id, length);
  }

  /**
   * Adds a axis used for correction.
   * @param axis Axis used for correction. The name of the axis corresponds to the name of a variable.
   */
  void AddCorrectionAxis(const Qn::AxisD& axis) { correction_axes_.Add(axis); }

  /**
   * Adds a correction step based on Q vectors to the specified detector
   * @tparam CORRECTION
   * @param name Name of the detector
   * @param correction the preconfigured correction step
   */
  template<typename CORRECTION>
  void AddCorrectionOnQnVector(const std::string &name, CORRECTION &&correction) {
    detectors_.FindDetector(name).AddCorrectionOnQnVector(std::forward<CORRECTION>(correction));
  }

  void SetChannelGroups(const std::string &name, const std::vector<int> &channel_groups) {
    detectors_.FindDetector(name).SetChannelScheme(channel_groups);
  }

  /**
   * Adds a correction step based on the input data to the specified detector
   * @tparam CORRECTION
   * @param name Name of the detector
   * @param correction the preconfigured correction step
   */
  template<typename CORRECTION>
  void AddCorrectionOnInputData(const std::string &name, CORRECTION &&correction) {
    detectors_.FindDetector(name).AddCorrectionOnInputData(std::forward<CORRECTION>(correction));
  }

  /**
   * Adds a event variable, which will be saved to the output tree.
   * @param name Name corresponds to a variable defined in the variable manager.
   */
  void AddEventVariable(const std::string &name) { variable_manager_.RegisterOutputF(name); }

  /**
   * Adds a event variable, which will be saved to the output tree.
   * Remember to add them as a variable first.
   * @param name Name corresponds to a variable defined in the variable manager.
   */
  void AddEventVariableInt(const std::string &name) { variable_manager_.RegisterOutputL(name); }

  /**
   * Adds a detector to the correction framework.
   * @param name Name of the detector
   * @param type Type of the detector: Channel or Track
   * @param phi_name Name of the variable which saves the angular information e.g. phi angle of a channel or particle
   * @param weight_name Name of the variable which saves a weight. "Ones" is used for a weight of 1.
   * @param axes Axes used for differential corrections. Names correspond to the name of variables.
   * @param harmonics activated harmonics in a bitset.
   */
  void AddDetector(std::string name,
                   Qn::DetectorType type,
                   const std::string &phi_name,
                   const std::string &weight_name,
                   const std::vector<Qn::AxisD> &axes,
                   const std::bitset<Qn::QVector::kmaxharmonics> &harmonics,
                   QVector::Normalization norm = QVector::Normalization::M,
                   const std::string &radial_offset_name = "Ones") {
    auto phi = variable_manager_.FindVariable(phi_name);
    auto weight = variable_manager_.FindVariable(weight_name);
    auto radial_offset = variable_manager_.FindVariable(radial_offset_name);
    detectors_.AddDetector(name, type, phi, weight, radial_offset, axes, harmonics, norm);
  }

  /**
   * @brief Add a detector to the correction framework. Old signature with
   * static polymorphism
   */
  template <size_t N>
  void AddDetector(std::string name,
                   Qn::DetectorType type,
                   const std::string &phi_name,
                   const std::string &weight_name,
                   const std::vector<Qn::AxisD> &axes,
                   int const(&harmonics_array)[N],
                   QVector::Normalization norm = QVector::Normalization::M,
                   const std::string &radial_offset_name = "Ones"
                   ) {
    std::bitset<Qn::QVector::kmaxharmonics> harmonics_bitset;
    for (size_t i = 0; i < N; ++i) {
      harmonics_bitset[harmonics_array[i] - 1] = true;
    }
    AddDetector(name, type, phi_name, weight_name, axes, harmonics_bitset, norm, radial_offset_name);
    return;
  }

  template<std::size_t N, typename Function>
  void AddCutOnDetector(const std::string &detector_name,
                        const char *const (&variable_names)[N],
                        Function cut_function,
                        const std::string &cut_description) {
    std::array<std::string, N> variable_names_arr;
    std::copy(std::begin(variable_names), std::end(variable_names), std::begin(variable_names_arr));
    bool is_channel_wise = variable_manager_.FindVariable(variable_names_arr[0]).size() > 1;
    detectors_.AddCut(detector_name,
                      CallBacks::MakeCut(variable_names_arr, cut_function, cut_description),
                      is_channel_wise);
  }

  template<typename Function>
  void AddCutOnDetector(const std::string &detector_name,
                        const std::vector<std::string> &variable_names,
                        Function cut_function,
                        const std::string &cut_description) {
    bool is_channel_wise = variable_manager_.FindVariable(variable_names[0]).size() > 1;
    detectors_.AddCut(detector_name,
                      CallBacks::MakeCut(variable_names, cut_function, cut_description),
                      is_channel_wise);
  }

  template<std::size_t N, typename FUNCTION>
  void AddEventCut(const char *const (&variable_names)[N],
                   FUNCTION cut_function,
                   const std::string &cut_description) {
    event_cuts_.AddCut(CallBacks::MakeCut(variable_names, cut_function, cut_description));
  }

  /**
   * @brief Adds a one dimensional event histogram
   * @param axes axis of the histogram. Name corresponds to the axis.
   * @param weight Name of the weights used when filling. Standard is "Ones" (1).
   */
  void AddEventHisto1D(const Qn::AxisD &axis, const std::string &weight = "Ones") {
    event_histograms_.Add("Event", axis, weight);
  }

  /**
   * @brief Adds a two n event histogram
   * @param axes axes of the histogram. Name corresponds to the axes.
   * @param weight Name of the weights used when filling. Standard is "Ones" (1).
   */
  void AddEventHisto2D(const std::vector<Qn::AxisD> &axes, const std::string &weight = "Ones") {
    event_histograms_.Add("Event", axes, weight);
  }

  void AddEventHisto2DArray(const std::vector<Qn::AxisD> &axes,
                            const Qn::AxisD &axis,
                            const std::string &weight = "Ones") {
    event_histograms_.Add("Event", axes, weight, axis);
  }

  /**
  * @brief Adds a one dimensional histogram to a detector.
  * @param Name name of the detector
  * @param axes axis of the histogram. Name corresponds to the axis.
  * @param weight Name of the weights used when filling. Standard is "Ones" (1).
  */
  void AddHisto1D(const std::string &detector, const Qn::AxisD &axis, const std::string &weight = "Ones") {
    detectors_.FindDetector(detector).AddHistogram(detector, axis, weight);
  }

  /**
  * Adds a two dimensional histogram to a detector.
  * @param Name name of the detector
  * @param axes axis of the histogram. Name corresponds to the axis.
  * @param weight Name of the weights used when filling. Standard is "Ones" (1).
  */
  void AddHisto2D(const std::string &detector,
                  const std::vector<Qn::AxisD> &axes,
                  const std::string &weight = "Ones") {
    detectors_.FindDetector(detector).AddHistogram(detector, axes, weight);
  }

  void SetOutputQVectors(const std::string &name, const std::vector<Qn::QVector::CorrectionStep> &steps) {
    auto &detector = detectors_.FindDetector(name);
    for (auto &step : steps) {
      detector.SetOutputQVector(step);
    }
  }
  void SetFillOutputTree(bool tree) { fill_output_tree_ = tree; }
  void SetFillCalibrationQA(bool calibration) { fill_qa_histos_ = calibration; }
  void SetFillValidationQA(bool validation) { fill_validation_qa_histos_ = validation; }
  void SetCurrentRunName(const std::string &name);
  void SetCalibrationInputFileName(const std::string &file_name) { correction_input_file_name_ = file_name; }
  void SetCalibrationInputFile(TFile *file) { correction_input_file_.reset(file); }

  /**
   * @brief Set output tree.
   * Lifetime of the tree is managed by the user.
   * @param tree non-owning pointer to the tree
   */
  void ConnectOutputTree(TTree *tree) { if (fill_output_tree_) out_tree_ = tree; }

  /**
   * @brief Initializes the correction framework
   * @param in_calibration_file_ non-owning pointer to the calibration file.
   * Lifetime of the file has to be managed by the user.
   */
  void InitializeOnNode();

  bool ProcessEvent();

  double *GetVariableContainer() { return variable_manager_.GetVariableContainer(); }

  inline void FillTrackingDetectors() { if (event_passed_cuts_) detectors_.FillTracking(); }
  inline void FillChannelDetectors() { if (event_passed_cuts_) detectors_.FillChannel(); }

  void ProcessCorrections();
  /**
   * @brief Resets the correction framework. To be called before a new event is processed.
   */
  void Reset();

  /**
 * @brief Finalizes the correction framework. To be called after all events are processed.
 */
  void Finalize();

  void CreateReport() {
    detectors_.CreateReport();
  }

  /**
   * @brief Get the list containing the calibration histograms.
   * @return A pointer of the list to which the calibration histograms will be saved.
   */
  TList *GetCorrectionList() { return correction_output.get(); }

  /**
   * @brief Get the list containing the calibration QA histograms.
   * @return A pointer of the list to which the calibration QA histograms will be saved.
   */
  TList *GetCorrectionQAList() { return correction_qa_histos_.get(); }

 private:
  void InitializeCorrections();
  void AttachQAHistograms();
  static constexpr auto kCorrectionListName = "CorrectionHistograms";
  bool fill_qa_histos_ = true; ///< Flag for filling QA histograms
  bool fill_validation_qa_histos_ = true; ///< Flag for filling calibration bin validation histograms
  bool fill_output_tree_ = false; ///< Flag for filling the output tree
  bool event_passed_cuts_ = false; ///< variable holding status if an event passed the cuts.
  RunList runs_; ///< list of processed runs
  DetectorList detectors_; ///< list of detectors
  InputVariableManager variable_manager_; ///< manager of the variables
  std::string correction_input_file_name_; ///< name of the calibration input file
  std::unique_ptr<TList> correction_input_;      //!<! the list of the input calibration histograms
  std::unique_ptr<TList> correction_output;      //!<! the list of the support histograms
  std::unique_ptr<TList> correction_qa_histos_;  //!<! the list of QA histograms
  std::unique_ptr<TFile> correction_input_file_; //!<! input calibration file
  CorrectionAxisSet correction_axes_; /// CorrectionCalculator correction axes
  CorrectionCuts event_cuts_; ///< Pointer to the event cuts
  QAHistograms event_histograms_; ///< event QA histograms
  TTree *out_tree_ = nullptr;  //!<! Tree of Qn Vectors and event variables. Lifetime is managed by the user.
 /// \cond CLASSIMP
 ClassDef(CorrectionManager, 1);
 /// \endcond
};
}

#endif //FLOW_CORRECTIONMANAGER_H
