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

/// \file QnCorrectionsCorrectionDetectorBase.cxx
/// \brief Implementation of the base detector configuration class within Q vector correction framework

#include "SubEvent.hpp"
#include "Detector.hpp"

/// \cond CLASSIMP
ClassImp(Qn::SubEvent);
/// \endcond
namespace Qn {
const char *SubEvent::szPlainQnVectorName = "plain";
const char *SubEvent::szQAQnAverageHistogramName = "Plain Qn avg ";

/// Incorporates the passed correction to the set of Q vector corrections
/// \param correctionOnQn the correction to add
void SubEvent::AddCorrectionOnQnVector(CorrectionOnQnVector *correctionOnQn) {
  correctionOnQn->SetOwner(this);
  fQnVectorCorrections.AddCorrection(correctionOnQn);
}

/// Incorporates the passed correction to the set of input data corrections
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param correctionOnInputData the correction to add
void SubEvent::AddCorrectionOnInputData(CorrectionOnInputData *correctionOnInputData) {
  (void) correctionOnInputData;
}

/// Get the corrected Qn vector from the step previous to the one given
/// If not previous step the plain Qn vector is returned.
/// The user is not able to modify it.
/// \param correctionOnQn the correction to find its predecessor corrected Qn vector
/// \return the corrected Qn vector from the correction step predecessor or the plain Qn vector
const QVector *SubEvent::GetPreviousCorrectedQnVector(CorrectionOnQnVector *correctionOnQn) const {
  if (fQnVectorCorrections.GetPrevious(correctionOnQn)) {
    return fQnVectorCorrections.GetPrevious(correctionOnQn)->GetCorrectedQnVector();
  } else {
    return &fPlainQnVector;
  }
}

/// Check if a concrete correction step is being applied on this detector configuration
/// It is not enough having the correction step configured or collecting data. To
/// get an affirmative answer the correction step must be being applied.
/// Transfers the request to the set of Qn vector corrections.
/// \param step the name of the correction step
/// \return TRUE if the correction step is being applied
Bool_t SubEvent::IsCorrectionStepBeingApplied(const char *step) const {
  return fQnVectorCorrections.IsCorrectionStepBeingApplied(step);
}

/// Activate the processing for the passed harmonic
/// \param harmonic the desired harmonic number to activate
void SubEvent::ActivateHarmonic(Int_t harmonic) {
  fPlainQnVector.ActivateHarmonic(harmonic);
  fCorrectedQnVector.ActivateHarmonic(harmonic);
  fPlainQ2nVector.ActivateHarmonic(harmonic);
  fCorrectedQ2nVector.ActivateHarmonic(harmonic);
  fTempQnVector.ActivateHarmonic(harmonic);
  fTempQ2nVector.ActivateHarmonic(harmonic);
}

std::string SubEvent::GetName() const {
  if (fDetector->GetBinName(binid_).empty()) {
    return fDetector->GetName();
  } else {
    return fDetector->GetName() + "_" + fDetector->GetBinName(binid_);
  }
}

void SubEvent::BuildQnVector() {
  fPlainQnVector.SetNormalization(QVector::Normalization::NONE);
  fPlainQ2nVector.SetNormalization(QVector::Normalization::NONE);
  for (const auto &dataVector : fDataVectorBank) {
    fPlainQnVector.Add(dataVector.Phi(), dataVector.RadialOffset(), dataVector.EqualizedWeight());
    fPlainQ2nVector.Add(dataVector.Phi(), dataVector.RadialOffset(), dataVector.EqualizedWeight());
  }
  /* check the quality of the Qn vector */
  fPlainQnVector.CheckQuality();
  fPlainQ2nVector.CheckQuality();
  fPlainQnVector = fPlainQnVector.Normal(fDetector->GetNormalizationMethod());
  fPlainQ2nVector = fPlainQ2nVector.Normal(fDetector->GetNormalizationMethod());
  fCorrectedQnVector = fPlainQnVector;
  fCorrectedQ2nVector = fPlainQ2nVector;
}

}// namespace Qn
