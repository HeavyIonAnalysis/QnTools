#ifndef QNCORRECTIONS_CORRECTIONSETONQNVECTOR_H
#define QNCORRECTIONS_CORRECTIONSETONQNVECTOR_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsCorrectionsSetOnQvector.h
/// \brief Set of corrections on Qn vector support within Q vector correction framework
///

#include <list>
#include <set>

#include "CorrectionBase.hpp"
#include "CorrectionOnQnVector.hpp"
#include "CorrectionOnInputData.hpp"

namespace Qn {
/// \class QnCorrectionsCorrectionsSetOnQvector
/// \brief Encapsulate the set of corrections to apply on Q vectors
///
/// Order matters so, the list must be built with the order in which
/// corrections should be applied.
///
/// The correction steps are own by the object instance so they will
/// be destroyed with it.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 05, 2016

template<typename T>
class CorrectionsSet {
 public:
  CorrectionsSet() = default;
  virtual ~CorrectionsSet() = default;
  typename std::list<std::unique_ptr<T>>::iterator begin() { return list_.begin(); }
  typename std::list<std::unique_ptr<T>>::iterator end() { return list_.end(); }
  typename std::list<std::unique_ptr<T>>::const_iterator begin() const { return list_.begin(); }
  typename std::list<std::unique_ptr<T>>::const_iterator end() const { return list_.end(); }

  void AddCorrection(T *correction) {
    if (list_.empty()) {
      list_.emplace_front(correction);
    } else if (*correction < *list_.front()) {
      list_.emplace_front(correction);
    } else if (*list_.back() < *correction) {
      list_.emplace_back(correction);
    } else {
      for (auto it = list_.begin(); it!=list_.end(); ++it) {
        if (*correction < **it) {
          list_.insert(it, std::unique_ptr<T>(correction));
          break;
        }
      }
    }
  }

  /// Fill the global list of correction steps
/// \param correctionlist (partial) global list of corrections ordered by correction key
  void FillOverallCorrectionsList(std::set<CorrectionBase *> &set) const {
    for (auto &entry : list_) {
      set.emplace(entry.get());
    }
  }

/// Gets the correction on Qn vector previous to the one passed as argument
/// \param correction the correction to find the previous one
/// \return the previous correction, NULL if none
  const CorrectionOnQnVector *GetPrevious(const T *correction) const {
    if (correction==nullptr) return nullptr;
    if (list_.empty()) return nullptr;
    if (list_.front()->GetName()==correction->GetName()) return nullptr;
    if (list_.size()==1) return nullptr;
    for (auto it = list_.begin(); it!=list_.end(); ++it) {
      auto nextit = std::next(it, 1);
      if ((*nextit)->GetName()==correction->GetName()) return (*it).get();
    }
    return nullptr;
  }

/// Check if a concrete correction step is being applied on this detector configuration
/// It is not enough having the correction step configured or collecting data. To
/// get an affirmative answer the correction step must be being applied.
/// Transfer the order to each of the Qn correction steps.
/// \param step the name of the correction step
/// \return TRUE if the correction step is being applied
  Bool_t IsCorrectionStepBeingApplied(const char *step) const {
    for (auto &entry : list_) {
      if (std::string(entry->GetName()).find(step)) {
        return entry->IsBeingApplied();
      }
    }
    return false;
  }

  bool IsLastStepApplied() const {
    if (!list_.empty()) {
      return (*list_.end())->IsBeingApplied();
    }
    return false;
  }

  std::map<std::string, std::pair<bool, bool>> ReportOnUsage() const {
    std::map<std::string, std::pair<bool, bool>> report;
    for (const auto &correction : list_) {
      report.emplace(correction->GetName(), correction->ReportUsage());
    }
    return report;
  }

  void CreateCorrectionHistograms() {
    for (auto &correction : list_) {
      correction->CreateCorrectionHistograms();
    }
  }

  void AttachInputs(TList *input_list) {
    T *previous_correction = nullptr;
    for (auto &correction : list_) {
      if (previous_correction) {
        if (previous_correction->GetState()==CorrectionBase::State::APPLYCOLLECT) {
          correction->Enable();
        }
      }
      if (correction->GetState()==CorrectionBase::State::CALIBRATION) {
        correction->AttachInput(input_list);
      }
      previous_correction = correction.get();
    }
  }

  void CopyToOutputList(TList *output_list) {
    for (auto &correction : list_) {
      correction->CopyToOutputList(output_list);
    }
  }

  void EnableFirstCorrection() {
    if (!list_.empty()) (*list_.begin())->Enable();
  }

  bool Empty() {
    return list_.empty();
  }

 private:
  std::list<std::unique_ptr<T>> list_;
/// \cond CLASSIMP
 ClassDef(CorrectionsSet, 3);
/// \endcond
};
// typedefs for convenience
using CorrectionsSetOnQvector = CorrectionsSet<CorrectionOnQnVector>;
using CorrectionsSetOnInputData = CorrectionsSet<CorrectionOnInputData>;
}
#endif // QNCORRECTIONS_CORRECTIONSETONQNVECTOR_H
