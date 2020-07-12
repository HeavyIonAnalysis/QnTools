#ifndef QN_DETECTOR_H
#define QN_DETECTOR_H

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

#include <utility>
#include <memory>
#include <utility>

#include "ROOT/RMakeUnique.hxx"

#include "CorrectionAxisSet.hpp"
#include "SubEventChannels.hpp"
#include "SubEventTracks.hpp"
#include "InputVariableManager.hpp"
#include "DataContainer.hpp"
#include "QVector.hpp"
#include "QAHistogram.hpp"
#include "CorrectionCuts.hpp"

namespace Qn {
class DetectorList;
class CorrectionManager;
/**
 * Enumerator class used to determine the type of detector
 */
enum class DetectorType {
  TRACK,
  CHANNEL
};

class Detector {
 public:
  Detector() = default;
  Detector(std::string name,
           DetectorType type,
           std::vector<AxisD> axes,
           InputVariable phi,
           InputVariable weight,
           InputVariable radial_offset,
           std::bitset<Qn::QVector::kmaxharmonics> harmonics,
           QVector::Normalization norm);
  virtual ~Detector() = default;
  Detector(Detector &&detector) = default;
  Detector &operator=(Detector &&detector) = default;
  Detector(const Detector &other);
  /**
   * @brief Clears data before filling new event.
   */
  void ClearData() {
    for (auto &bin : sub_events_) {
      bin->Clear();
    }
  }
  /**
   * @brief Adds a cut to the detector
   * @param cut unique pointer to the cut.
   */
  void AddCut(const CorrectionCut::CallBack &callback, bool is_channel_wise) {
    if (is_channel_wise) {
      cuts_.AddCut(callback);
    } else {
      int_cuts_.AddCut(callback);
    }
  }
  /**
   * Adds a qa histogram to the detector
   * @tparam Args
   * @param args arguments of the QAHistogram constructor
   */
  template<typename... Args>
  void AddHistogram(Args &&... args) {
    histograms_.Add(std::forward<Args>(args)...);
  }
  /**
   * Adds a correction of type CorrectionOnQnVector to the detector
   * @tparam CORRECTION
   * @param correction preconfigured correction step which is added to the detector
   */
  template<typename CORRECTION>
  void AddCorrectionOnQnVector(const CORRECTION &correction) {
    auto corr = new CORRECTION(correction);
    correction_on_q_vector.Add(corr);
  }
  /**
   * Adds a correction of type CorrectionOnInputData to the detector
   * @tparam CORRECTION
   * @param correction preconfigured correction step which is added to the detector
   */
  template<typename CORRECTION>
  void AddCorrectionOnInputData(const CORRECTION &correction) {
    auto corr = new CORRECTION(correction);
    correction_on_input_data.Add(corr);
  }
  /**
 * Configure which corrected Q-vectors are propagated to the output.
 * @param step
 */
  void SetOutputQVector(QVector::CorrectionStep step) { output_tree_q_vectors_.emplace_back(step); }

  void Initialize(DetectorList &detectors, InputVariableManager &var, CorrectionAxisSet &correction_axis);
  void FillData();
  void FillReport() {
    int_cuts_.FillReport();
    cuts_.FillReport();
  }
  void CreateSupportQVectors() { for (auto &ev : sub_events_) { ev->CreateSupportQVectors(); }}

  void AttachCorrectionInputs(TList *list) {
    for (auto &ev : sub_events_) {
      ev->AttachCorrectionInput(list);
    }
  }

  void AfterInputAttachAction() {
    for (auto &ev : sub_events_) {
      ev->AfterInputAttachAction();
    }
  }

  void CreateCorrectionHistograms() {
    for (auto &ev : sub_events_) {
      ev->CreateCorrectionHistograms();
    }
  }
  void CopyToOutputList(TList* list) {
    for (auto &ev : sub_events_) {
      ev->CopyToOutputList(list);
    }
  }


  bool IsIntegrated() const { return sub_events_.IsIntegrated(); }
  void ProcessCorrections();
  void IncludeQnVectors();
  void AttachToTree(TTree *tree);
  void ActivateHarmonic(unsigned int i) {
    harmonics_bits_.set(i - 1);
    for (auto &ev : sub_events_) {
      ev->ActivateHarmonic(i);
    }
  }

  void SetChannelScheme(std::vector<int> channel_groups) {
    channel_groups_ = channel_groups;
  }

  DetectorList *GetDetectors() const { return detectors_; }
  Qn::QVector::Normalization GetNormalizationMethod() const { return q_vector_normalization_method_; }
  std::string GetName() const { return name_; }
  std::string GetBinName(unsigned int id) const { return sub_events_.GetBinDescription(id); }
  SubEvent *GetSubEvent(unsigned int ibin) { return sub_events_.At(ibin).get(); }
  TList *CreateQAHistogramList(bool fill_qa, bool fill_validation);

  DataContainerQVector *GetQVector(QVector::CorrectionStep step) { return q_vectors_.at(step).get(); }

 private:
  InputVariable phi_; /// variable holding the azimuthal angle
  InputVariable weight_; /// variable holding the weight which is used for the calculation of the Q vector.
  InputVariable radial_offset_; /// variable holding the radial offset
  std::string name_; /// name of  the detector
  DetectorType type_; /// type of detector
  int nchannels_ = 0; /// number of channels in case of channel detector
  std::bitset<Qn::QVector::kmaxharmonics> harmonics_bits_; /// bitset of all activated harmonics
  Qn::QVector::Normalization q_vector_normalization_method_ = Qn::QVector::Normalization::NONE;
  std::vector<InputVariable> input_variables_; //!<! variables used for the binning of the Q vector.
  std::vector<float> coordinates_;  //!<!  vector holding the temporary coordinates of one track or channel.
  std::map<QVector::CorrectionStep, std::unique_ptr<DataContainerQVector>> q_vectors_; //!<! output qvectors
  std::vector<QVector::CorrectionStep> output_tree_q_vectors_; /// Holds correction steps used for the output
  CorrectionCuts cuts_; /// per channel selection  cuts
  CorrectionCuts int_cuts_; /// integrated selection cuts
  QAHistograms histograms_; /// QA histograms of the detector
  std::vector<Qn::AxisD> axes_; /// Holds axes till they are used to configure the subevents
  Qn::DataContainer<std::unique_ptr<SubEvent>, AxisD> sub_events_; //!<! SubEvents of the detector
  Qn::DetectorList *detectors_ = nullptr; /// Pointer to the list of detectors
  TObjArray correction_on_q_vector; /// Holds the correction steps till they are used to configure the sub events
  TObjArray correction_on_input_data; /// Holds the correction steps till they are used to configure the sub events

  std::vector<int> channel_groups_; /// for gain equalization


  /// \cond CLASSIMP
 ClassDef(Detector, 2);
  /// \endcond
};

}

#endif //FLOW_DETECTOR_H
