// Qn Tools
//
// Copyright (C) 2019  Lukas Kreis, Ilya Selyuzhenkov
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

#ifndef FLOW_DETECTORLIST_H
#define FLOW_DETECTORLIST_H

#include "Detector.hpp"
namespace Qn {
class DetectorList {
 public:
  DetectorList() = default;
  virtual ~DetectorList() = default;
  void AddDetector(std::string name,
                   Qn::DetectorType type,
                   const InputVariable& phi,
                   const InputVariable& weight,
                   const InputVariable& radial_offset,
                   const std::vector<Qn::AxisD> &axes,
                   std::bitset<Qn::QVector::kmaxharmonics> harmonics,
                   QVector::Normalization norm) {
    auto find = [&name](const Detector &a) { return name==a.GetName(); };
    if (type==Qn::DetectorType::CHANNEL) {
      if (std::find_if(channel_detectors_.begin(), channel_detectors_.end(), find)==channel_detectors_.end()) {
        channel_detectors_.emplace_back(name, type, axes, phi, weight, radial_offset, harmonics, norm);

      }
    } else if (type==Qn::DetectorType::TRACK) {
      if (std::find_if(tracking_detectors_.begin(), tracking_detectors_.end(), find)==tracking_detectors_.end()) {
        tracking_detectors_.emplace_back(name, type, axes, phi, weight, radial_offset, harmonics, norm);
      }
    }
  }

  void AddCut(const std::string &name, const CorrectionCut::CallBack& cut, bool is_channel_wise) {
    auto &det = FindDetector(name);
    det.AddCut(cut, is_channel_wise);
  }

  Detector &FindDetector(const std::string name) {
    auto find = [&name](const Detector &a) { return a.GetName()==name; };
    auto det_ch = std::find_if(channel_detectors_.begin(), channel_detectors_.end(), find);
    auto det_trk = std::find_if(tracking_detectors_.begin(), tracking_detectors_.end(), find);
    if (det_ch!=channel_detectors_.end()) {
      return *det_ch;
    } else if (det_trk!=tracking_detectors_.end()) {
      return *det_trk;
    } else {
      throw std::runtime_error("Detector" + name + "not found.");
    }
  }

  void SetOutputTree(TTree *output_tree) {
    if (output_tree) {
      for (auto &detector : tracking_detectors_) {
        detector.AttachToTree(output_tree);
      }
      for (auto &detector : channel_detectors_) {
        detector.AttachToTree(output_tree);
      }
    }
  }

  void Initialize(DetectorList &detectors, InputVariableManager &var, CorrectionAxisSet &axes) {
    for (auto &detector : channel_detectors_) {
      all_detectors_.push_back(&detector);
      detector.Initialize(detectors, var, axes);
    }
    for (auto &detector : tracking_detectors_) {
      all_detectors_.push_back(&detector);
      detector.Initialize(detectors, var, axes);
    }
  }

  void FillTracking() {
    for (auto &dp : tracking_detectors_) {
      dp.FillData();
    }
  }

  void FillChannel() {
    for (auto &dp : channel_detectors_) {
      dp.FillData();
    }
  }

  void ProcessCorrections() {
    for (auto &d : all_detectors_) {
      d->ProcessCorrections();
    }
  }

  void CreateSupportQVectors() {
    for (auto &d : all_detectors_) {
      d->CreateSupportQVectors();
    }
  }

  void FillReport() {
    for (auto &d : all_detectors_) {
      d->FillReport();
    }
  }

  void CreateCorrectionHistograms() {
    for (auto &d : all_detectors_) {
      d->CreateCorrectionHistograms();
    }
  }

  void CopyToOutputList(TList* list) {
    for (auto &d : all_detectors_) {
      d->CopyToOutputList(list);
    }
  }

  void AttachCorrectionInput(TList *list) {
    for (auto &d : all_detectors_) {
      d->AttachCorrectionInputs(list);
    }
    for (auto &d : all_detectors_) {
      d->AfterInputAttachAction();
    }
  }

  void AttachQAHistograms(TList *list, bool fill_qa, bool fill_validation) {
    for (auto &detector: all_detectors_) {
      auto detector_list = detector->CreateQAHistogramList(fill_qa, fill_validation);
      list->Add(detector_list);
    }
  }

  void ResetDetectors() {
    for (auto &d : all_detectors_) {
      d->ClearData();
    }
  }

  void IncludeQnVectors() {
    for (auto &d : all_detectors_) {
      d->IncludeQnVectors();
    }
  }

  void CreateReport() {
    auto iteration = CalculateProgress(all_detectors_);
    std::cout << "iteration " << iteration.first << " of " << iteration.second << std::endl;
    for (auto &d : all_detectors_) {
      auto corrections = d->GetSubEvent(0)->ReportOnCorrections();
      std::cout << d->GetName() << std::endl;
      for (const auto &step : corrections) {
        std::cout << step.first << " : ";
        if (step.second.first) {
          std::cout << "collecting ";
        }
        if (step.second.second) {
          std::cout << "applying";
        }
        if (!step.second.first && !step.second.second) {
          std::cout << "waiting";
        }
        std::cout << std::endl;
      }
    }
  }

 private:

  std::pair<int, int> CalculateProgress(const std::vector<Detector *> &detectors) {
    int remaining_iterations_global = 0;
    int total_iterations_global = 0;
    for (const auto &d : detectors) {
      auto corrections = d->GetSubEvent(0)->ReportOnCorrections();
      int performed_iterations = 0;
      int total_iterations = corrections.size();
      for (auto &c : corrections) { if (c.second.second) performed_iterations++; }
      if (total_iterations - performed_iterations > remaining_iterations_global) {
        remaining_iterations_global = total_iterations - performed_iterations;
      }
      if (total_iterations_global < total_iterations) total_iterations_global = total_iterations;
    }
    return {total_iterations_global - remaining_iterations_global + 1, total_iterations_global + 1};
  }

  std::vector<Detector> tracking_detectors_; ///< vector of tracking detectors
  std::vector<Detector> channel_detectors_; ///< vector of channel detectors
  std::vector<Detector *> all_detectors_; ///!<! storing pointers to all detectors

  /// \cond CLASSIMP
 ClassDef(DetectorList, 1);
  /// \endcond
};
}
#endif