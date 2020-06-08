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

/// \file QnCorrectionsQnVectorAlignment.cxx
/// \brief Implementation of procedures for Qn vector alignment correction.
#include "Alignment.hpp"
#include "CorrectionAxisSet.hpp"
#include "DetectorList.hpp"
#include "ROOT/RMakeUnique.hxx"

/// \cond CLASSIMP
ClassImp(Qn::Alignment);
/// \endcond
namespace Qn {
/// Set the detector configuration used as reference for alignment
/// The detector configuration name is stored for further use.
/// \param name the name of the reference detector configuration
void Alignment::SetReferenceConfigurationForAlignment(const char *name) {
  //  QnCorrectionsInfo(Form("Reference name: %s, attached to detector configuration: %s",
  //                         name,
  //                         ((fDetector) ? "yes" : "no")));
  fDetectorForAlignmentName = name;
  /* we could be in different situations of framework attachment */
  /* so, we do nothing for the time being */
}

/// Asks for support data structures creation
///
/// Locates the reference detector configuration for alignment if its name has been previously stored
/// Creates the recentered Qn vector
void Alignment::CreateSupportQVectors() {
  /* now, definitely, we should have the reference detector configurations */
  if (!fDetectorForAlignmentName.empty()) {
    auto &aligndetector = fSubEvent->GetDetector()->GetDetectors()->FindDetector(fDetectorForAlignmentName);
    if (!aligndetector.IsIntegrated()) throw std::logic_error("Alignment detector needs to be integrated.");
    fDetectorForAlignment = aligndetector.GetSubEvent(0);
  } else {
    throw std::logic_error("Correctionstep not configured. Please add detector for alignment.");
  }
  /* make sure the alignment harmonic processing is active */
  fSubEvent->GetDetector()->ActivateHarmonic(fHarmonicForAlignment);
  /* in both configurations */
  fDetectorForAlignment->GetDetector()->ActivateHarmonic(fHarmonicForAlignment);
  fInputQnVector = fSubEvent->GetPreviousCorrectedQnVector(this);
  /* and now create the corrected Qn vector */
  fCorrectedQnVector = std::make_unique<QVector>(fSubEvent->GetHarmonics(), QVector::CorrectionStep::ALIGNED,
                                                 fInputQnVector->GetNorm());
}

/// Asks for support histograms creation
///
/// Allocates the histogram objects and creates the calibration histograms.
///
/// Process concurrency requires Calibration Histograms creation for all
/// concurrent processes but not for Input Histograms so, we delete previously
/// allocated ones.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void Alignment::CreateCorrectionHistograms() {
  auto hname = std::string(szSupportHistogramName) + "_" + fSubEvent->GetName() + "#times"
      + fDetectorForAlignment->GetName();
  fInputHistograms =
      std::make_unique<CorrectionProfileCorrelationComponents>(hname, fSubEvent->GetEventClassVariablesSet());
  fInputHistograms->SetNoOfEntriesThreshold(fMinNoOfEntriesToValidate);
  fCalibrationHistograms =
      std::make_unique<CorrectionProfileCorrelationComponents>(hname, fSubEvent->GetEventClassVariablesSet());
  fCalibrationHistograms->CreateCorrelationComponentsProfileHistograms(&output_histograms);
}

/// Attaches the needed input information to the correction step
/// \param list list where the inputs should be found
/// \return kTRUE if everything went OK
void Alignment::AttachInput(TList *list) {
  if (fInputHistograms->AttachHistograms(list)) {
    fState = State::APPLYCOLLECT;
  }
}

/// Asks for QA histograms creation
///
/// Allocates the histogram objects and creates the QA histograms.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void Alignment::AttachQAHistograms(TList *list) {
  auto hname = std::string(szQAQnAverageHistogramName) + "_" + fSubEvent->GetName();
  fQAQnAverageHistogram =
      std::make_unique<CorrectionProfileComponents>(hname, fSubEvent->GetEventClassVariablesSet());
  /* get information about the configured harmonics to pass it for histogram creation */
  Int_t nNoOfHarmonics = fSubEvent->GetNoOfHarmonics();
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
void Alignment::AttachNveQAHistograms(TList *list) {
  auto hname = std::string(szQANotValidatedHistogramName) + "_" + fSubEvent->GetName();
  fQANotValidatedBin =
      std::make_unique<CorrectionHistogramSparse>(hname, fSubEvent->GetEventClassVariablesSet());
  fQANotValidatedBin->CreateHistogram(list);
}

/// Processes the correction step
///
/// Apply the correction step
/// \return kTRUE if the correction step was applied
bool Alignment::ProcessCorrections() {
  bool applied = false;
  switch (fState) {
    case State::CALIBRATION:
      /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
      /* we have not perform any correction yet */
      break;
    case State::APPLYCOLLECT:
      /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
      /* and proceed to ... */
      /* FALLTHRU */
    case State::APPLY: /* apply the correction if the current Qn vector is good enough */
      /* logging */
      //      QnCorrectionsInfo(Form("Alignment process in detector %s with reference %s: applying correction.",
      //                             fDetector->GetName(),
      //                             fDetectorConfigurationForAlignment->GetName()));
      if (fSubEvent->GetCurrentQnVector()->IsGoodQuality()) {
        /* we get the properties of the current Qn vector but its name */
        fCorrectedQnVector->CopyNumberOfContributors(*fSubEvent->GetCurrentQnVector());
        /* let's check the correction histograms */
        Long64_t bin = fInputHistograms->GetBin();
        if (fInputHistograms->BinContentValidated(bin)) {
          /* the bin content is validated so, apply the correction */
          Double_t XX = fInputHistograms->GetXXBinContent(bin);
          Double_t YY = fInputHistograms->GetYYBinContent(bin);
          Double_t XY = fInputHistograms->GetXYBinContent(bin);
          Double_t YX = fInputHistograms->GetYXBinContent(bin);
          Double_t eXY = fInputHistograms->GetXYBinError(bin);
          Double_t eYX = fInputHistograms->GetYXBinError(bin);
          Double_t deltaPhi = -TMath::ATan2((XY - YX), (XX + YY)) * (1.0 / fHarmonicForAlignment);
          /* significant correction? */
          if (!(TMath::Sqrt((XY - YX) * (XY - YX) / (eXY * eXY + eYX * eYX)) < 2.0)) {
            Int_t harmonic = fSubEvent->GetCurrentQnVector()->GetFirstHarmonic();
            while (harmonic != -1) {
              fCorrectedQnVector->SetX(harmonic,
                                       fSubEvent->GetCurrentQnVector()->x(harmonic)
                                               * TMath::Cos(((Double_t) harmonic) * deltaPhi)
                                           + fSubEvent->GetCurrentQnVector()->y(harmonic)
                                               * TMath::Sin(((Double_t) harmonic) * deltaPhi));
              fCorrectedQnVector->SetY(harmonic,
                                       fSubEvent->GetCurrentQnVector()->y(harmonic)
                                               * TMath::Cos(((Double_t) harmonic) * deltaPhi)
                                           - fSubEvent->GetCurrentQnVector()->x(harmonic)
                                               * TMath::Sin(((Double_t) harmonic) * deltaPhi));
              harmonic = fSubEvent->GetCurrentQnVector()->GetNextHarmonic(harmonic);
            }
          } /* if the correction is not significant we leave the Q vector untouched */
        }   /* if the correction bin is not validated we leave the Q vector untouched */
        else {
          if (fQANotValidatedBin) fQANotValidatedBin->Fill(1.0);
        }
      } else {
        /* not done! input Q vector with bad quality */
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
/// Collect data for the correction step.
/// \return kTRUE if the correction step was applied
Bool_t Alignment::ProcessDataCollection() {
  bool applied = false;
  switch (fState) {
    case State::CALIBRATION:
      /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
      if ((fInputQnVector->IsGoodQuality()) && (fDetectorForAlignment->GetCurrentQnVector()->IsGoodQuality())) {
        fCalibrationHistograms->FillXX(fInputQnVector->x(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->x(
                                           fHarmonicForAlignment));
        fCalibrationHistograms->FillXY(fInputQnVector->x(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->y(
                                           fHarmonicForAlignment));
        fCalibrationHistograms->FillYX(fInputQnVector->y(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->x(
                                           fHarmonicForAlignment));
        fCalibrationHistograms->FillYY(fInputQnVector->y(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->y(
                                           fHarmonicForAlignment));
      }
      /* we have not perform any correction yet */
      break;
    case State::APPLYCOLLECT:
      /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
      if ((fInputQnVector->IsGoodQuality()) && (fDetectorForAlignment->GetCurrentQnVector()->IsGoodQuality())) {
        fCalibrationHistograms->FillXX(fInputQnVector->x(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->x(
                                           fHarmonicForAlignment));
        fCalibrationHistograms->FillXY(fInputQnVector->x(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->y(
                                           fHarmonicForAlignment));
        fCalibrationHistograms->FillYX(fInputQnVector->y(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->x(
                                           fHarmonicForAlignment));
        fCalibrationHistograms->FillYY(fInputQnVector->y(fHarmonicForAlignment)
                                       * fDetectorForAlignment->GetCurrentQnVector()->y(
                                           fHarmonicForAlignment));
      }
      /* and proceed to ... */
      /* FALLTHRU */
    case State::APPLY: /* apply the correction if the current Qn vector is good enough */
      /* provide QA info if required */
      if (fQAQnAverageHistogram) {
        Int_t harmonic = fCorrectedQnVector->GetFirstHarmonic();
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
void Alignment::ClearCorrectionStep() {
  fCorrectedQnVector->Reset();
}

}// namespace Qn
