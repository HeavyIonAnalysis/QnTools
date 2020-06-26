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

#include "CorrectionManager.hpp"
#include "TList.h"

namespace Qn {

void CorrectionManager::InitializeCorrections() {
  // Connects the correction histogram list
  if (!correction_input_file_) {
    correction_input_file_ = std::make_unique<TFile>(correction_input_file_name_.data(), "READ");
  }
  if (correction_input_file_ && !correction_input_file_->IsZombie()) {
    auto input = dynamic_cast<TList *>(correction_input_file_->FindObjectAny(kCorrectionListName));
    correction_input_.reset(input);
    if (input) correction_input_->SetOwner(true);
  }
  // Prepares the correctionsteps
  detectors_.CreateSupportQVectors();
  correction_output = std::make_unique<TList>();
  correction_output->SetName(kCorrectionListName);
  correction_output->SetOwner(true);
}

void CorrectionManager::SetCurrentRunName(const std::string &name) {
  runs_.SetCurrentRun(name);
  TList *current_output = nullptr;
  if (!runs_.empty()) {
    current_output = new TList();
    current_output->SetName(runs_.GetCurrent().data());
    current_output->SetOwner(true);
    correction_output->Add(current_output);
    detectors_.CreateCorrectionHistograms();
  }
  if (correction_input_) {
    auto current_run = (TList *) correction_input_->FindObject(runs_.GetCurrent().data());
    if (current_run) {
      detectors_.AttachCorrectionInput(current_run);
    }
  }
  detectors_.CopyToOutputList(current_output);
  detectors_.IncludeQnVectors();
  if (fill_output_tree_ && out_tree_) {
    detectors_.SetOutputTree(out_tree_);
    variable_manager_.SetOutputTree(out_tree_);
  }
  detectors_.CreateReport();
}

void CorrectionManager::AttachQAHistograms() {
  correction_qa_histos_ = std::make_unique<TList>();
  correction_qa_histos_->SetName("QA_histograms");
  correction_qa_histos_->SetOwner(true);
  detectors_.AttachQAHistograms(correction_qa_histos_.get(),
                                fill_qa_histos_,
                                fill_validation_qa_histos_);
  auto event_qa_list = new TList();
  event_qa_list->SetOwner(true);
  event_qa_list->SetName("event_QA");
  event_cuts_.CreateCutReport("Cut_Report:");
  event_cuts_.AddToList(event_qa_list);
  event_histograms_.AddToList(event_qa_list);
  correction_qa_histos_->Add(event_qa_list);
}

void CorrectionManager::InitializeOnNode() {
  variable_manager_.Initialize();
  correction_axes_.Initialize(variable_manager_);
  event_histograms_.Initialize(variable_manager_);
  detectors_.Initialize(detectors_, variable_manager_, correction_axes_);
  event_cuts_.Initialize(variable_manager_);
  InitializeCorrections();
  AttachQAHistograms();
}

bool CorrectionManager::ProcessEvent() {
  event_passed_cuts_ = event_cuts_.CheckCuts(0);
  return event_passed_cuts_;
}

void CorrectionManager::ProcessCorrections() {
  if (event_passed_cuts_) {
    event_cuts_.FillReport();
    variable_manager_.UpdateOutVariables();
    event_histograms_.Fill();
    detectors_.ProcessCorrections();
    detectors_.FillReport();
    if (fill_output_tree_) out_tree_->Fill();
  }
}

void CorrectionManager::Reset() {
  event_passed_cuts_ = false;
  detectors_.ResetDetectors();
}

void CorrectionManager::Finalize() {
  auto calibration_list = (TList *) correction_output->FindObject(runs_.GetCurrent().data());
  if (calibration_list) {
    correction_output->Add(calibration_list->Clone("all"));
  }
}

}// namespace Qn