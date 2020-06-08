/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/

/// \file QnCorrectionsDetectorConfigurationTracks.cxx
/// \brief Implementation of the track detector configuration class
#include <utility>

#include "CorrectionProfileComponents.hpp"
#include "SubEventTracks.hpp"

/// \cond CLASSIMP
ClassImp(Qn::SubEventTracks);
/// \endcond
namespace Qn {

/// Normal constructor
/// Allocates the data vector bank.
/// \param name the name of the detector configuration
/// \param eventClassesVariables the set of event classes variables
/// \param nNoOfHarmonics the number of harmonics that must be handled
/// \param harmonicMap an optional ordered array with the harmonic numbers
SubEventTracks::SubEventTracks(unsigned int bin_id,
                               const CorrectionAxisSet *eventClassesVariables,
                               std::bitset<QVector::kmaxharmonics> harmonics) : SubEvent(bin_id, eventClassesVariables, harmonics) {
}

/// Asks for support data structures creation
///
/// The input data vector bank is allocated and the request is
/// transmitted to the Q vector corrections.
void SubEventTracks::CreateSupportQVectors() {

  /* this is executed in the remote node so, allocate the data bank */
  fDataVectorBank.reserve(Qn::SubEvent::INITIALSIZE);
  for (auto &correction : fQnVectorCorrections) {
    correction->CreateSupportQVectors();
  }
}

/// Asks for support histograms creation
///
/// The request is transmitted to the Q vector corrections.
///
/// A new histograms list is created for the detector configuration and incorporated
/// to the passed list. Then the new list is passed to the corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void SubEventTracks::CreateCorrectionHistograms() {
  fQnVectorCorrections.CreateCorrectionHistograms();
  fQnVectorCorrections.EnableFirstCorrection();
}

void SubEventTracks::CopyToOutputList(TList *list) {
  auto correction_list = new TList();
  correction_list->SetName(GetName().data());
  correction_list->SetOwner(kTRUE);
  fQnVectorCorrections.CopyToOutputList(correction_list);
  /* if list is empty delete it if not incorporate it */
  if (!correction_list->IsEmpty()) {
    list->Add(correction_list);
  } else {
    delete correction_list;
  }
}

/// Asks for QA histograms creation
///
/// The request is transmitted to the Q vector corrections.
///
/// A new histograms list is created for the detector and incorporated
/// to the passed list. Then the new list is passed to the corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void SubEventTracks::AttachQAHistograms(TList *list) {
  TList *detectorConfigurationList;
  if (GetName().empty()) {
    detectorConfigurationList = list;
  } else {
    detectorConfigurationList = new TList();
    detectorConfigurationList->SetName(GetName().data());
    detectorConfigurationList->SetOwner(kTRUE);
  }
  /* the own QA average Qn vector components histogram */
  auto name = std::string(szQAQnAverageHistogramName) + " " + GetName();
  fQAQnAverageHistogram = std::make_unique<CorrectionProfileComponents>(name, GetEventClassVariablesSet());
  /* get information about the configured harmonics to pass it for histogram creation */
  Int_t nNoOfHarmonics = this->GetNoOfHarmonics();
  auto harmonicsMap = new Int_t[nNoOfHarmonics];
  this->GetHarmonicMap(harmonicsMap);
  fQAQnAverageHistogram->CreateComponentsProfileHistograms(detectorConfigurationList, nNoOfHarmonics, harmonicsMap);
  delete[] harmonicsMap;
  for (auto &correction : fQnVectorCorrections) {
    correction->AttachQAHistograms(detectorConfigurationList);
  }
  /* if list is empty delete it if not incorporate it */
  if (!detectorConfigurationList->IsEmpty()) {
    list->Add(detectorConfigurationList);
  } else {
    delete detectorConfigurationList;
  }
}

/// Asks for non validated entries QA histograms creation
///
/// The request is transmitted to the Q vector corrections.
///
/// A new histograms list is created for the detector and incorporated
/// to the passed list. Then the new list is passed to the corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void SubEventTracks::AttachNveQAHistograms(TList *list) {
  auto detectorConfigurationList = new TList();
  detectorConfigurationList->SetName(GetName().data());
  detectorConfigurationList->SetOwner(kTRUE);
  for (auto &correction : fQnVectorCorrections) {
    correction->AttachNveQAHistograms(detectorConfigurationList);
  }
  /* if list is empty delete it if not incorporate it */
  if (!detectorConfigurationList->IsEmpty()) {
    list->Add(detectorConfigurationList);
  } else {
    delete detectorConfigurationList;
  }
}

/// Asks for attaching the needed input information to the correction steps
///
/// The detector list is extracted from the passed list and then
/// the request is transmitted to the Q vector corrections with the found list.
/// \param list list where the input information should be found
/// \return kTRUE if everything went OK
void SubEventTracks::AttachCorrectionInput(TList *list) {
  auto detectorConfigurationList = (TList *) list->FindObject(GetName().data());
  if (detectorConfigurationList) {
    fQnVectorCorrections.EnableFirstCorrection();
    fQnVectorCorrections.AttachInputs(detectorConfigurationList);
  }
}

/// Perform after calibration histograms attach actions
/// It is used to inform the different correction step that
/// all conditions for running the network are in place so
/// it is time to check if their requirements are satisfied
///
/// The request is transmitted to the Q vector corrections
void SubEventTracks::AfterInputAttachAction() {
  /* now propagate it to Q vector corrections */
  for (auto &correction : fQnVectorCorrections) {
    correction->AfterInputAttachAction();
  }
}

/// Fills the QA plain Qn vector average components histogram
/// \param variableContainer pointer to the variable content bank
void SubEventTracks::FillQAHistograms() {
  if (fQAQnAverageHistogram) {
    Int_t harmonic = fPlainQnVector.GetFirstHarmonic();
    while (harmonic != -1) {
      fQAQnAverageHistogram->FillX(harmonic, fPlainQnVector.x(harmonic));
      fQAQnAverageHistogram->FillY(harmonic, fPlainQnVector.y(harmonic));
      harmonic = fPlainQnVector.GetNextHarmonic(harmonic);
    }
  }
}

/// Include the the list of Qn vector associated to the detector configuration
/// into the passed list
///
/// A new list is created for the detector configuration and incorporated
/// to the passed list.
///
/// Always includes first the fully corrected Qn vector,
/// and then includes the plain Qn vector and asks to the different correction
/// steps to include their partially corrected Qn vectors.
/// The check if we are already there is because it could be late information
/// about the process name and then the correction histograms could still not
/// be attached and the constructed list does not contain the final Qn vectors.
/// \param list list where the corrected Qn vector should be added
void SubEventTracks::IncludeQnVectors() {
  qvectors_.emplace(QVector::CorrectionStep::PLAIN, &fPlainQnVector);
  for (auto &correction : fQnVectorCorrections) {
    correction->IncludeCorrectedQnVector(qvectors_);
  }
}

/// Include only one instance of each input correction step
/// in execution order
///
/// There are not input correction so we do nothing
/// \param list list where the correction steps should be incorporated
void SubEventTracks::FillOverallInputCorrectionStepList(std::set<CorrectionBase *> &set) const {
  (void) set;
}

/// Include only one instance of each Qn vector correction step
/// in execution order
///
/// The request is transmitted to the set of Qn vector corrections
/// \param list list where the correction steps should be incorporated
void SubEventTracks::FillOverallQnVectorCorrectionStepList(std::set<CorrectionBase *> &set) const {
  fQnVectorCorrections.FillOverallCorrectionsList(set);
}

/// Provide information about assigned corrections
///
/// We create three list which items they own, incorporate info from the
/// correction steps and add them to the passed list
/// \param steps list for incorporating the list of assigned correction steps
/// \param calib list for incorporating the list of steps in calibrating status
/// \param apply list for incorporating the list of steps in applying status

//void SubEventTracks::ReportOnCorrections(std::vector<std::string> steps,
//                                         std::vector<std::string> calib,
//                                         std::vector<std::string> apply) const {
//  auto mysteps = new TList();
//  mysteps->SetOwner(kTRUE);
//  mysteps->SetName(GetName().data());
//  auto mycalib = new TList();
//  mycalib->SetOwner(kTRUE);
//  mycalib->SetName(GetName().data());
//  auto myapply = new TList();
//  myapply->SetOwner(kTRUE);
//  myapply->SetName(GetName().data());
//  /* incorporate Qn vector corrections */
//  Bool_t keepIncorporating = kTRUE;
//  for (auto &correction : fQnVectorCorrections) {
//    mysteps->Add(new TObjString(correction->GetName()));
//    /* incorporate additional info if the step will be reached */
//    if (keepIncorporating) {
//      Bool_t keep = correction->ReportUsage(mycalib, myapply);
//      keepIncorporating = keepIncorporating && keep;
//    }
//  }
//  steps->Add(mysteps);
//  calib->Add(mycalib);
//  apply->Add(myapply);
//}

}// namespace Qn