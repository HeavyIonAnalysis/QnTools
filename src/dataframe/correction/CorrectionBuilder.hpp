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

#ifndef QNTOOLS_RECENTERVECTOR_HPP_
#define QNTOOLS_RECENTERVECTOR_HPP_

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <TFile.h>

#include "CorrectionAction.hpp"

namespace Qn::Correction {
template <typename EventAxes>
class CorrectionBuilder {
 public:
  using CorrectionAction = typename Qn::Correction::CorrectionAction<
      EventAxes, typename EventAxes::AxisValueTypeTuple>;
  using ResultPtr = typename ROOT::RDF::RResultPtr<CorrectionAction>;
  using InitializationObject =
      typename CorrectionAction::InitializationObject;
  explicit CorrectionBuilder(EventAxes event_axes) : event_axes_(event_axes) {}

  void SetMinimumEntries(int min_entries) { min_entries_ = min_entries; }

  void AddCorrection(std::string qvector, std::string this_step,
                     std::string previous_step) {
    corrections_.push_back(
        MakeRecenterAction(this_step, event_axes_, qvector, previous_step));
  }

  void EnableDiagonalization(std::string name) {
    for (auto &correction : corrections_) {
      if (correction.GetInputName() == name) {
        correction.EnableDiagonalization();
      }
    }
  }

  void EnableRescaling(std::string name) {
    for (auto &correction : corrections_) {
      if (correction.GetInputName() == name) {
        correction.EnableRescaling();
      }
    }
  }

  void EnableWidthEqualization(std::string name) {
    for (auto &correction : corrections_) {
      if (correction.GetInputName() == name) {
        correction.EnableWidthEqualization();
      }
    }
  }


  void AddCorrection(std::string qvector, std::string this_step,
                     std::string previous_step,
                     InitializationObject initialization) {
    corrections_.push_back(
        MakeRecenterAction(this_step, event_axes_, qvector, previous_step));
    initializations_.emplace(corrections_.size()-1, std::make_unique<InitializationObject>(initialization));
  }

  template <typename DataFrame>
  auto Apply(DataFrame &df, TFile *correction_file, TTreeReader *reader = nullptr) {
    auto dir = correction_file;
    /// Disable calculation of correction factors, if the result was found in
    /// the corrections input file. The correction factors will not to be
    /// recalculated.
    int iq = -1;
    corrections_.erase(
        std::remove_if(
            std::begin(corrections_), std::end(corrections_),
            [this, &iq, &dir, reader, &df](auto &correction) {
              ++iq;
              if (initializations_[iq]) {
                if (!correction.LoadCorrectionFromFile(dir, *initializations_[iq])) {
                  result_ptrs_.push_back(
                      Qn::MakeAverageHelper(correction.SetMinimumNumberOfEntries(min_entries_))
                          .SetInitializationWithInitializationObject(
                              initializations_[iq].get())
                          .BookMe(df));
                  return true;
                } else {
                  return false;
                }
              } else {
                if (!correction.LoadCorrectionFromFile(dir, *reader)) {
                  result_ptrs_.push_back(
                      Qn::MakeAverageHelper(correction).BookMe(df));
                  return true;
                } else {
                  return false;
                }
              }
            }),
        corrections_.end());
    auto corrected_temp =
        Qn::Correction::ApplyCorrectionsVector(df, result_ptrs_);
    auto corrected =
        Qn::Correction::ApplyCorrectionsVector(corrected_temp, corrections_);
    return corrected;
  }

  void Write(TDirectory* dir) {
    for (auto & result_ptr : result_ptrs_) result_ptr->Write(dir);
  }

  const CorrectionAction & GetRecenterCorrection(const std::string &name) {
    for (auto correction : result_ptrs_) {
      std::cout << correction->GetName() << std::endl;
      if (correction->GetName() == name) {
        return correction.GetValue();
      }
    }
    throw std::out_of_range(name+ " not found in corrections");
  }

 private:
  std::map<int, std::unique_ptr<InitializationObject>> initializations_;
  std::vector<CorrectionAction> corrections_;
  std::vector<ResultPtr> result_ptrs_;
  int min_entries_ = 0;
  EventAxes event_axes_;

};
}  // namespace Qn::Correction

#endif  // QNTOOLS_RECENTERVECTOR_HPP_
