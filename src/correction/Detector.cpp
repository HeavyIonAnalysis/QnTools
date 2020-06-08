
#include "Detector.hpp"
#include "CorrectionManager.hpp"

namespace Qn {

Detector::Detector(std::string name,
                   DetectorType type,
                   std::vector<AxisD> axes,
                   InputVariable phi,
                   InputVariable weight,
                   InputVariable radial_offset,
                   std::bitset<Qn::QVector::kmaxharmonics> harmonics,
                   QVector::Normalization norm) : phi_(phi),
                                                  weight_(weight),
                                                  radial_offset_(radial_offset),
                                                  name_(std::move(name)),
                                                  type_(type),
                                                  nchannels_(phi_.size()),
                                                  harmonics_bits_(harmonics),
                                                  q_vector_normalization_method_(norm),
                                                  axes_(axes),
                                                  correction_on_q_vector(),
                                                  correction_on_input_data() {
}

Detector::Detector(const Detector &other) : phi_(other.phi_),
                                            weight_(other.weight_),
                                            radial_offset_(other.radial_offset_),
                                            name_(other.name_),
                                            type_(other.type_),
                                            nchannels_(other.nchannels_),
                                            harmonics_bits_(other.harmonics_bits_),
                                            q_vector_normalization_method_(other.q_vector_normalization_method_),
                                            output_tree_q_vectors_(other.output_tree_q_vectors_),
                                            cuts_(other.cuts_),
                                            int_cuts_(other.int_cuts_),
                                            histograms_(other.histograms_),
                                            axes_(other.axes_),
                                            correction_on_q_vector(other.correction_on_q_vector),
                                            correction_on_input_data(other.correction_on_input_data) {
}

/**
 * Initializes the detectors before computing the corrections
 * @param detectors the list of other detectors (needed for some corrections)
 * @param var variable manager needed to configure the data inputs
 * @param correction_axis the set of axes used for calculating the corrections for classes
 */
void Detector::Initialize(DetectorList &detectors, InputVariableManager &var, CorrectionAxisSet &correction_axis) {
  sub_events_.AddAxes(axes_);
  detectors_ = &detectors;
  var.InitVariable(phi_);
  var.InitVariable(weight_);
  var.InitVariable(radial_offset_);
  int ibin = 0;
  for (auto &event : sub_events_) {
    if (type_ == DetectorType::CHANNEL) {
      event = std::make_unique<SubEventChannels>(ibin, &correction_axis, nchannels_, harmonics_bits_);
      if (channel_groups_.empty()) {
        for (int i = 0; i < nchannels_; ++i) {
          channel_groups_.push_back(0);
        }
      }
      event->SetChannelsScheme(channel_groups_);
    } else if (type_ == DetectorType::TRACK) {
      event = std::make_unique<SubEventTracks>(ibin, &correction_axis, harmonics_bits_);
    }
    event->SetDetector(this);
    for (int i = 0; i < correction_on_input_data.GetEntriesFast(); ++i) {
      event->AddCorrectionOnInputData(dynamic_cast<CorrectionOnInputData *>(correction_on_input_data.At(i))->MakeCopy());
    }
    for (int i = 0; i < correction_on_q_vector.GetEntriesFast(); ++i) {
      event->AddCorrectionOnQnVector(dynamic_cast<CorrectionOnQnVector *>(correction_on_q_vector.At(i))->MakeCopy());
    }
    ++ibin;
  }
  if (!sub_events_.IsIntegrated()) {
    for (const auto &axis : sub_events_.GetAxes()) {
      input_variables_.push_back(var.FindVariable(axis.Name()));
    }
    coordinates_.resize(input_variables_.size());
  }
  // Initialize the cuts
  cuts_.Initialize(var);
  int_cuts_.Initialize(var);
  histograms_.Initialize(var);
}

TList *Detector::CreateQAHistogramList(bool fill_qa, bool fill_validation) {
  auto list = new TList();
  list->SetName(name_.data());
  // Add detector qa
  auto detector_qa = new TList();
  detector_qa->SetOwner(true);
  detector_qa->SetName("detector_QA");
  int_cuts_.CreateCutReport(name_ + ":", 1);
  cuts_.CreateCutReport(name_ + ":", phi_.size());
  int_cuts_.AddToList(detector_qa);
  cuts_.AddToList(detector_qa);
  histograms_.AddToList(detector_qa);
  list->Add(detector_qa);
  // Add qvector qa
  if (fill_qa) {
    auto qvector_qa = new TList();
    qvector_qa->SetName("q_vector_QA");
    qvector_qa->SetOwner(true);
    for (auto &ev : sub_events_) ev->AttachQAHistograms(qvector_qa);
    list->Add(qvector_qa);
  }
  // Add qvector bin validation qa
  if (fill_validation) {
    auto valid_qa = new TList();
    valid_qa->SetName("bin_validation_QA");
    valid_qa->SetOwner(true);
    for (auto &ev : sub_events_) { ev->AttachNveQAHistograms(valid_qa); }
    list->Add(valid_qa);
  }
  return list;
}

void Detector::ProcessCorrections() {
  for (auto &ev : sub_events_) { ev->ProcessCorrections(); }
  for (auto &ev : sub_events_) { ev->ProcessDataCollection(); }
  // passes the corrected Q-vectors to the output container.
  for (auto &pair_step_qvector : q_vectors_) {
    for (unsigned int i = 0; i < sub_events_.size(); ++i) {
      try {
        (*pair_step_qvector.second)[i] = *sub_events_[i]->GetQVector(pair_step_qvector.first);
      } catch (std::out_of_range &) {
        throw std::out_of_range(name_ + " bin " + std::to_string(i) + " correctionstep: " + kCorrectionStepNamesArray[pair_step_qvector.first] + "not found.");
      }
    }
  }
}

void Detector::AttachToTree(TTree *tree) {
  for (const auto &qvec : q_vectors_) {
    auto is_output_variable = std::find(output_tree_q_vectors_.begin(), output_tree_q_vectors_.end(), qvec.first);
    if (is_output_variable != output_tree_q_vectors_.end()) {
      auto suffix = kCorrectionStepNamesArray[qvec.first];
      auto name = name_ + "_" + suffix;
      tree->Branch(name.data(), qvec.second.get());
    }
  }
}

void Detector::IncludeQnVectors() {
  for (auto &ev : sub_events_) { ev->IncludeQnVectors(); }
  // Adds DataContainerQVector for each of the active correction steps.
  auto correction_steps = sub_events_[0]->GetCorrectionSteps();
  for (auto correction_step : correction_steps) {
    if (sub_events_.IsIntegrated()) {
      q_vectors_.emplace(correction_step, std::make_unique<DataContainerQVector>());
    } else {
      q_vectors_.emplace(correction_step, std::make_unique<DataContainerQVector>(sub_events_.GetAxes()));
    }
  }
}

void Detector::FillData() {
  if (!int_cuts_.CheckCuts(0)) return;
  histograms_.Fill();
  for (std::size_t channel = 0; channel < phi_.size(); ++channel) {
    if (!cuts_.CheckCuts(channel)) continue;
    /// Integrated case (detector only has one bin)
    if (input_variables_.empty()) {
      sub_events_.At(0)->AddDataVector(channel, phi_[channel], weight_[channel], radial_offset_[channel]);
      /// differential case (detector has more than one bin)
    } else {
      for (std::size_t coordinate = 0; coordinate < input_variables_.size(); ++coordinate) {
        coordinates_.at(coordinate) = input_variables_[coordinate][channel];
      }
      const auto ibin = sub_events_.FindBin(coordinates_);
      if (ibin > -1) {
        sub_events_.At(ibin)->AddDataVector(channel, phi_[channel], weight_[channel], radial_offset_[channel]);
      }
    }
  }
}

}// namespace Qn