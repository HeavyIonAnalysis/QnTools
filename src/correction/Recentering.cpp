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

/// \file QnCorrectionsQnVectorRecentering.cxx
/// \brief Implementation of procedures for Qn vector recentering.
#include <memory>

#include "CorrectionAxisSet.hpp"
#include "CorrectionHistogramSparse.hpp"
#include "CorrectionProfileComponents.hpp"
#include "ROOT/RMakeUnique.hxx"
#include "Recentering.hpp"
#include "SubEvent.hpp"

/// \cond CLASSIMP
ClassImp(Qn::Recentering);
/// \endcond
namespace Qn {
const Int_t Recentering::fDefaultMinNoOfEntries = 2;
const char *Recentering::szCorrectionName = "Recentering and width equalization";
const char *Recentering::szSupportHistogramName = "Qn";
const char *Recentering::szCorrectedQnVectorName = "rec";
const char *Recentering::szQANotValidatedHistogramName = "Rec NvE";
const char *Recentering::szQAQnAverageHistogramName = "Rec Qn avg ";

/// Default constructor
/// Passes to the base class the identity data for the recentering and width equalization correction step
Recentering::Recentering() : CorrectionOnQnVector(szCorrectionName, CorrectionOnQnVector::Step::kRecentering) {
  fApplyWidthEqualization = kFALSE;
  fMinNoOfEntriesToValidate = fDefaultMinNoOfEntries;
}

/// Asks for support data structures creation
///
/// Creates the recentered Qn vector
void Recentering::CreateSupportQVectors() {
  fInputQnVector = fSubEvent->GetPreviousCorrectedQnVector(this);
  fCorrectedQnVector = std::make_unique<QVector>(fSubEvent->GetHarmonics(),
                                                 QVector::CorrectionStep::RECENTERED,
                                                 fInputQnVector->GetNorm());
}

/// Asks for support histograms creation
///
/// Allocates the histogram objects and creates the calibration histograms.
/// The histograms are constructed with standard deviation error calculation
/// for the proper behavior of optional gain equalization step.
///
/// Process concurrency requires Calibration Histograms creation for all c
/// concurrent processes but not for Input Histograms so, we delete previously
/// allocated ones.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void Recentering::CreateCorrectionHistograms() {
  auto hname = std::string(szSupportHistogramName) + "_" + fSubEvent->GetName();
  fInputHistograms = std::make_unique<CorrectionProfileComponents>(hname, fSubEvent->GetEventClassVariablesSet(),
                                                                   CorrectionHistogramBase::ErrorMode::SPREAD);
  fInputHistograms->SetNoOfEntriesThreshold(fMinNoOfEntriesToValidate);
  fCalibrationHistograms = std::make_unique<CorrectionProfileComponents>(hname, fSubEvent->GetEventClassVariablesSet(),
                                                                         CorrectionHistogramBase::ErrorMode::SPREAD);
  /* get information about the configured harmonics to pass it for histogram creation */
  Int_t nNoOfHarmonics = fSubEvent->GetNoOfHarmonics();
  auto harmonicsMap = new Int_t[nNoOfHarmonics];
  fSubEvent->GetHarmonicMap(harmonicsMap);
  fCalibrationHistograms->CreateComponentsProfileHistograms(&output_histograms, nNoOfHarmonics, harmonicsMap);
  delete[] harmonicsMap;
}

/// Attaches the needed input information to the correction step
/// \param list list where the inputs should be found
/// \return kTRUE if everything went OK
void Recentering::AttachInput(TList *list) {
  if (fInputHistograms->AttachHistograms(list)) {
    fState = State::APPLYCOLLECT;
  }
}

/// Asks for QA histograms creation
///
/// Allocates the histogram objects and creates the QA histograms.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void Recentering::AttachQAHistograms(TList *list) {
  auto hname = std::string(szQAQnAverageHistogramName) + "_" + fSubEvent->GetName();
  fQAQnAverageHistogram = std::make_unique<CorrectionProfileComponents>(hname, fSubEvent->GetEventClassVariablesSet());
  /* get information about the configured harmonics to pass it for histogram creation */
  auto nNoOfHarmonics = fSubEvent->GetNoOfHarmonics();
  auto harmonicsMap = new Int_t[nNoOfHarmonics];
  fSubEvent->GetHarmonicMap(harmonicsMap);
  fQAQnAverageHistogram->CreateComponentsProfileHistograms(list, nNoOfHarmonics, harmonicsMap);
  delete[] harmonicsMap;
}

/// Asks for non validated entries QA histograms creation
///
/// Allocates the histogram objects and creates the non validated entries QA histograms.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void Recentering::AttachNveQAHistograms(TList *list) {
  auto hname = std::string(szQANotValidatedHistogramName) + "_" + fSubEvent->GetName();
  fQANotValidatedBin = std::make_unique<CorrectionHistogramSparse>(hname, fSubEvent->GetEventClassVariablesSet());
  fQANotValidatedBin->CreateHistogram(list);
}

/// Processes the correction step
///
/// Pure virtual function
/// \return kTRUE if the correction step was applied
bool Recentering::ProcessCorrections() {
  int harmonic;
  bool applied = false;
  switch (fState) {
    case State::CALIBRATION:
      /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
      /* we have not perform any correction yet */
      break;
    case State::APPLYCOLLECT:
      /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
      /* and proceed to ... */
    case State::APPLY: /* apply the correction if the current Qn vector is good enough */
      if (fSubEvent->GetCurrentQnVector()->IsGoodQuality()) {
        /* we get the properties of the current Qn vector but its name */
        fCorrectedQnVector->CopyNumberOfContributors(*fSubEvent->GetCurrentQnVector());
        harmonic = fSubEvent->GetCurrentQnVector()->GetFirstHarmonic();
        /* let's check the correction histograms */
        Long64_t bin = fInputHistograms->GetBin();
        if (fInputHistograms->BinContentValidated(bin)) {
          /* correction information validated */
          while (harmonic != -1) {
            Float_t widthX = 1.0;
            Float_t widthY = 1.0;
            if (fApplyWidthEqualization) {
              widthX = fInputHistograms->GetXBinError(harmonic, bin);
              widthY = fInputHistograms->GetYBinError(harmonic, bin);
            }
            fCorrectedQnVector->SetX(harmonic, (fSubEvent->GetCurrentQnVector()->x(harmonic) - fInputHistograms->GetXBinContent(harmonic, bin)) / widthX);
            fCorrectedQnVector->SetY(harmonic, (fSubEvent->GetCurrentQnVector()->y(harmonic) - fInputHistograms->GetYBinContent(harmonic, bin)) / widthY);
            harmonic = fSubEvent->GetCurrentQnVector()->GetNextHarmonic(harmonic);
          }
        } /* correction information not validated, we leave the Q vector untouched */
        else {
          if (fQANotValidatedBin) fQANotValidatedBin->Fill(1.0);
        }
      } else {
        /* not done! input vector with bad quality */
        fCorrectedQnVector->SetGood(kFALSE);
      }
      /* and update the current Qn vector */
      fSubEvent->UpdateCurrentQnVector(*fCorrectedQnVector);
      applied = true;
      break;
    case State::PASSIVE:
      /* we are in passive state waiting for proper conditions, no corrections applied */
      break;
  }
  return applied;
}

/// Processes the correction step data collection
///
/// Pure virtual function
/// \return kTRUE if the correction step was applied
bool Recentering::ProcessDataCollection() {
  int harmonic;
  bool applied = false;
  switch (fState) {
    case State::CALIBRATION:
      /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
      if (fInputQnVector->IsGoodQuality()) {

        harmonic = fInputQnVector->GetFirstHarmonic();
        while (harmonic != -1) {
          fCalibrationHistograms->FillX(harmonic, fInputQnVector->x(harmonic));
          fCalibrationHistograms->FillY(harmonic, fInputQnVector->y(harmonic));
          harmonic = fInputQnVector->GetNextHarmonic(harmonic);
        }
      }
      /* we have not perform any correction yet */
      break;
    case State::APPLYCOLLECT:
      /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
      if (fInputQnVector->IsGoodQuality()) {
        harmonic = fInputQnVector->GetFirstHarmonic();
        while (harmonic != -1) {
          fCalibrationHistograms->FillX(harmonic, fInputQnVector->x(harmonic));
          fCalibrationHistograms->FillY(harmonic, fInputQnVector->y(harmonic));
          harmonic = fInputQnVector->GetNextHarmonic(harmonic);
        }
      }
      /* and proceed to ... */
      /* FALLTHRU */
    case State::APPLY: /* apply the correction if the current Qn vector is good enough */
      /* provide QA info if required */
      if (fQAQnAverageHistogram) {
        harmonic = fCorrectedQnVector->GetFirstHarmonic();
        while (harmonic != -1) {
          fQAQnAverageHistogram->FillX(harmonic, fCorrectedQnVector->x(harmonic));
          fQAQnAverageHistogram->FillY(harmonic, fCorrectedQnVector->y(harmonic));
          harmonic = fCorrectedQnVector->GetNextHarmonic(harmonic);
        }
      }
      applied = true;
      break;
    case State::PASSIVE:
      /* we are in passive state waiting for proper conditions, no corrections applied */
      break;
  }
  return applied;
}

/// Clean the correction to accept a new event
void Recentering::ClearCorrectionStep() {
  fCorrectedQnVector->Reset();
}

}// namespace Qn