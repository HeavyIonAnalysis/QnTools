#ifndef QNCORRECTIONS_DETECTORCONFTRACKS_H
#define QNCORRECTIONS_DETECTORCONFTRACKS_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsDetectorConfigurationTracks.h
/// \brief Track detector configuration class for Q vector correction framework
///

#include <iostream>

#include "CorrectionDataVector.hpp"
#include "SubEvent.hpp"
namespace Qn {
class CorrectionProfileComponents;

/// \class QnCorrectionsDetectorConfigurationTracks
/// \brief Track detector configuration within Q vector correction framework
///
/// A track detector within the Q vector correction framework is defined
/// as one for which its data vectors only involve azimuthal angles and a
/// potential weight. Apart from that no other input data calibration is
/// available.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016

class SubEventTracks : public SubEvent {
 public:
  friend class CorrectionBase;
  friend class SubEvent;
  SubEventTracks() = default;
  SubEventTracks(unsigned int bin_id,
                 const CorrectionAxisSet *eventClassesVariables,
                 std::bitset<QVector::kmaxharmonics> harmonics);

  /// Get if the detector configuration is own by a tracking detector
  /// \return TRUE, this is tracking detector configuration
  virtual Bool_t GetIsTrackingDetector() const { return kTRUE; }

  virtual void CreateSupportQVectors();
  virtual void CreateCorrectionHistograms();
  virtual void CopyToOutputList(TList* list);
  virtual void AttachQAHistograms(TList *list);
  virtual void AttachNveQAHistograms(TList *list);
  virtual void AttachCorrectionInput(TList *list);
  virtual void AfterInputAttachAction();

  virtual Bool_t ProcessCorrections();
  virtual Bool_t ProcessDataCollection();

  virtual void IncludeQnVectors();
  virtual void FillOverallInputCorrectionStepList(std::set<CorrectionBase *> &set) const;
  virtual void FillOverallQnVectorCorrectionStepList(std::set<CorrectionBase *> &set) const;
  virtual void Clear();
  virtual std::map<std::string, Report> ReportOnCorrections() const {
    return fQnVectorCorrections.ReportOnUsage();
  }
  
 private:
  /* QA section */
  void FillQAHistograms();

/// \cond CLASSIMP
 ClassDef(SubEventTracks, 2);
/// \endcond
};

/// Clean the configuration to accept a new event
///
/// Transfers the order to the Q vector correction steps and
/// cleans the own Q vector and the input data vector bank
/// for accepting the next event.
inline void SubEventTracks::Clear() {
  /* transfer the order to the Q vector corrections */
  for (auto &correction : fQnVectorCorrections) {
    correction->ClearCorrectionStep();
  }
  /* clean the own Q vectors */
  fPlainQnVector.Reset();
  fPlainQ2nVector.Reset();
  fCorrectedQnVector.Reset();
  fCorrectedQ2nVector.Reset();
  /* and now clear the the input data bank */
  fDataVectorBank.clear();
}

/// Ask for processing corrections for the involved detector configuration
///
/// The request is transmitted to the Q vector correction steps.
/// The first not applied correction step breaks the loop and kFALSE is returned
/// \return kTRUE if all correction steps were applied
inline Bool_t SubEventTracks::ProcessCorrections() {
  /* first we build the Q vector with the chosen calibration */
  BuildQnVector();
  /* then we transfer the request to the Q vector correction steps */
  /* the loop is broken when a correction step has not been applied */
  for (auto &correction : fQnVectorCorrections) {
    if (correction->ProcessCorrections())
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}

/// Ask for processing corrections data collection for the involved detector configuration
/// Fill own QA histogram information and then
/// the request is transmitted to the Q vector correction steps.
/// The first not applied correction step breaks the loop and kFALSE is returned
/// \return kTRUE if all correction steps were applied
inline Bool_t SubEventTracks::ProcessDataCollection() {
  FillQAHistograms();
  /* we transfer the request to the Q vector correction steps */
  /* the loop is broken when a correction step has not been applied */
  for (auto &correction : fQnVectorCorrections) {
    if (correction->ProcessDataCollection())
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}
}
#endif // QNCORRECTIONS_DETECTORCONFTRACKS_H
