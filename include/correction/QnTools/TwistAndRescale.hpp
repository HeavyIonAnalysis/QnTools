#ifndef QNCORRECTIONS_QNVECTORTWISTANDRESCALE_H
#define QNCORRECTIONS_QNVECTORTWISTANDRESCALE_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsQnVectorTwistAndRescale.h
/// \brief Definition of the class that implements twist and rescale on Q vectors.
///
/// Twist and rescale are applied on ongoing Q vector from the involved detector
/// configuration. It is supposed that such a Q vector has already been built after input data
/// equalization and is with the proper normalization for the twist and rescaling corrections being applied.
/// Both twist and rescale application are user configurable.
///
/// Twist correction is applied according to:
/// \f[
///        Q_{n,(x,y)}' = \frac{Q_{n,(x,y)} - \Lambda^{s(-,+)}_{2n} Q_{n,(y,x)}}{1 - \Lambda^{s-}_{2n}\Lambda^{s+}_{2n}}
/// \f]
///
/// The rescaling correction is applied according to:
/// \f[
///     Q''_{n(x,y)} = \frac{Q'_{n(x,y)}}{A_{2n}^{(+,-)}}
/// \f]
/// Parameters \f$ A_{2n}^{(+,-)} \f$ and \f$ \Lambda^{s(+,-)}_{2n} \f$ are extracted from data according to the
/// method selected for the twist and rescale correction. Currently two methods are supported: the double harmonic
/// method and the correlations method.
///
/// For the double harmonic method
/// \f[
///     A^{\pm}_{2n} = 1 \pm \langle X_{2n} \rangle
/// \f]
/// and
/// \f[
///     \Lambda ^{\pm}_{2n} = \frac{\langle Y_{2n} \rangle}{A^{\pm}_{2n}}
/// \f]
/// For the correlations method two additional subdetectors need to be configured. Let these two subdetectors be
/// denominated \f$ B \f$ and \f$ C \f$ being \f$ A \f$ the current subdetector on which the twist and rescaling
/// correction will be applied. \f$ B \f$ must be a tracking subdetector and must have been twist corrected for the
/// method to work properly. Then
///
/// \f[
/// A^{A+}_{2n} = \frac{\sqrt{2 \langle X^{A}_{n} X^{C}_{n} \rangle}\langle X^{A}_{n} X^{B}_{n} \rangle}
/// {\sqrt{\langle X^{A}_{n} X^{B}_{n} \rangle \langle X^{B}_{n} X^{C}_{n} \rangle + \langle X^{A}_{n} Y^{B}_{n} \rangle \langle X^{B}_{n} Y^{C}_{n} \rangle}},
/// \f]
///
/// \f[
/// A^{A-}_{2n} = \frac{\sqrt{2 \langle X^{A}_{n} X^{C}_{n} \rangle}\langle Y^{A}_{n} Y^{B}_{n} \rangle}
/// {\sqrt{\langle X^{A}_{n} X^{B}_{n} \rangle \langle X^{B}_{n} X^{C}_{n} \rangle + \langle X^{A}_{n} Y^{B}_{n} \rangle \langle X^{B}_{n} Y^{C}_{n} \rangle}},
/// \f]
///
/// \f[
///  \Lambda^{A,+}_{2n} =  \frac{\langle X^{A}_{n} Y^{B}_{n} \rangle}{\langle X^{A}_{n} X^{B}_{n} \rangle}
/// \f]
/// and
/// \f[
///  \Lambda^{A,-}_{2n} =  \frac{\langle X^{A}_{n} Y^{B}_{n} \rangle}{\langle Y^{A}_{n} Y^{B}_{n} \rangle}
/// \f]
///
/// Twist and rescale are only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper correction histograms that will
/// provide the required averages per event class.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the correction histograms.
///
/// Data collection for twist and rescale correction parameters building in the double harmonic method is
/// performed on plain Qn vector and by using only the harmonics n such as its double 2n is within the range
/// of the harmonics handled by the configuration.
///
/// Data collection for twist and rescale correction parameters building in the correlations method is
/// performed on the highest corrected Qn vector on the involved detectors.
///
/// Corrections are performed for the harmonics for which there are data collection support.

/* harmonic multiplier */

#include "CorrectionHistogramSparse.hpp"
#include "CorrectionOnQnVector.hpp"
#include "CorrectionProfile3DCorrelations.hpp"
#include "CorrectionProfileComponents.hpp"
namespace Qn {
/// \class QnCorrectionsQnVectorTwistAndRescale
/// \brief Encapsulates twist and rescale on Q vector
///
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jun 22, 2016
///
/// Twist and rescale are applied on ongoing Q vector from the involved detector
/// configuration.
///
/// Twist correction is applied according to:
/// \f[
///        Q_{n,(x,y)}' = \frac{Q_{n,(x,y)} - \Lambda^{s(-,+)}_{2n} Q_{n,(y,x)}}{1 - \Lambda^{s-}_{2n}\Lambda^{s+}_{2n}}
/// \f]
///
/// The rescaling correction is applied according to:
/// \f[
///     Q''_{n(x,y)} = \frac{Q'_{n(x,y)}}{A_{2n}^{(+,-)}}
/// \f]
/// Parameters \f$ A_{2n}^{(+,-)} \f$ and \f$ \Lambda^{s(+,-)}_{2n} \f$ are extracted from data according to the
/// method selected for the twist and rescale correction. Currently two methods are supported: the double harmonic
/// method and the correlations method.
///
/// Twist and rescale are only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper correction histograms that will
/// provide the required averages per event class.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the correction histograms.
///
/// Data collection for twist and rescale correction parameters building in the double harmonic method is
/// performed on plain Qn vector and by using only the harmonics n such as its double 2n is within the range
/// of the harmonics handled by the configuration.
///
/// Data collection for twist and rescale correction parameters building in the correlations method is
/// performed on the highest corrected Qn vector on the involved detectors.
///
/// Correction are performed for the harmonics for which there are data collection support.
///

class CorrectionHistogramSparse;

class TwistAndRescale : public CorrectionOnQnVector {
 public:
  /// \enum QnTwistAndRescaleMethod
  /// \brief The class of the id of the supported twist and rescale methods
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  ///
  enum class Method {
    DOUBLE_HARMONIC,///< \f$ A^{\pm}_{2n} = 1 \pm \langle X_{2n} \rangle, \qquad \Lambda ^{\pm}_{2n} = \frac{\langle Y_{2n} \rangle}{A^{\pm}_{2n}} \f$
    CORRELATIONS    ///< \f$ A^{A+}_{2n} = \frac{\sqrt{2 \langle X^{A}_{n} X^{C}_{n} \rangle}\langle X^{A}_{n} X^{B}_{n} \rangle}
    /// {\sqrt{\langle X^{A}_{n} X^{B}_{n} \rangle \langle X^{B}_{n} X^{C}_{n} \rangle + \langle X^{A}_{n} Y^{B}_{n} \rangle \langle X^{B}_{n} Y^{C}_{n} \rangle}}, \\ \,
    /// A^{A-}_{2n} = \frac{\sqrt{2 \langle X^{A}_{n} X^{C}_{n} \rangle}\langle Y^{A}_{n} Y^{B}_{n} \rangle}
    /// {\sqrt{\langle X^{A}_{n} X^{B}_{n} \rangle \langle X^{B}_{n} X^{C}_{n} \rangle + \langle X^{A}_{n} Y^{B}_{n} \rangle \langle X^{B}_{n} Y^{C}_{n} \rangle}}, \\ \,
    ///  \Lambda^{A,+}_{2n} =  \frac{\langle X^{A}_{n} Y^{B}_{n} \rangle}{\langle X^{A}_{n} X^{B}_{n} \rangle}, \quad
    ///  \Lambda^{A,-}_{2n} =  \frac{\langle X^{A}_{n} Y^{B}_{n} \rangle}{\langle Y^{A}_{n} Y^{B}_{n} \rangle} \f$
  };

  TwistAndRescale();
  virtual ~TwistAndRescale() = default;

  TwistAndRescale(const TwistAndRescale &other) : CorrectionOnQnVector(other),
                                                  fTwistAndRescaleMethod(other.fTwistAndRescaleMethod),
                                                  fApplyTwist(other.fApplyTwist),
                                                  fApplyRescale(other.fApplyRescale),
                                                  fBDetectorConfigurationName(other.fBDetectorConfigurationName),
                                                  fSubEventB(other.fSubEventB),
                                                  fCDetectorConfigurationName(other.fCDetectorConfigurationName),
                                                  fSubEventC(other.fSubEventC),
                                                  fMinNoOfEntriesToValidate(other.fMinNoOfEntriesToValidate) {
  }

  virtual CorrectionOnQnVector *MakeCopy() const { return dynamic_cast<CorrectionOnQnVector *>(new TwistAndRescale(*this)); }

  /// Sets the method for extracting twist and rescale correction parameters
  /// \param method the chosen method
  void SetTwistAndRescaleMethod(Method method) { fTwistAndRescaleMethod = method; }
  /// Controls if twist step shall be applied
  /// \param apply kTRUE for applying the twist step
  void SetApplyTwist(Bool_t apply) { fApplyTwist = apply; }
  /// Controls if rescale step shall be applied
  /// \param apply kTRUE for applying the rescale step
  void SetApplyRescale(Bool_t apply) { fApplyRescale = apply; }
  void SetReferenceConfigurationsForTwistAndRescale(const char *nameB, const char *nameC);
  /// Set the minimum number of entries for calibration histogram bin content validation
  /// \param nNoOfEntries the number of entries threshold
  void SetNoOfEntriesThreshold(Int_t nNoOfEntries) { fMinNoOfEntriesToValidate = nNoOfEntries; }
  virtual void AttachInput(TList *list);
  virtual void AfterInputAttachAction();
  virtual void CreateSupportQVectors();
  virtual void CreateCorrectionHistograms();
  virtual void AttachQAHistograms(TList *list);
  virtual void AttachNveQAHistograms(TList *list);
  virtual Bool_t ProcessCorrections();
  virtual Bool_t ProcessDataCollection();
  virtual void ClearCorrectionStep();
  virtual void IncludeCorrectedQnVector(std::map<QVector::CorrectionStep, QVector *> &qvectors) const;
  virtual void IncludeCorrectionStep(std::vector<QVector::CorrectionStep> &steps) {
    if (fApplyRescale) {
      steps.push_back(fRescaleCorrectedQnVector->GetCorrectionStep());
    }
    if (fApplyTwist) {
      steps.push_back(fTwistCorrectedQnVector->GetCorrectionStep());
    }
  }

 private:
  using State = Qn::CorrectionBase::State;
  static constexpr const unsigned int szPriority =
      CorrectionOnQnVector::Step::kTwistAndRescale;///< the key of the correction step for ordering purpose
  static const Int_t fDefaultMinNoOfEntries;       ///< the minimum number of entries for bin content validation
  static const Double_t fMaxThreshold;             ///< highest absolute value for meaningful results
  static const char *szTwistCorrectionName;        ///< the name of the twist correction step
  static const char *szRescaleCorrectionName;      ///< the name of the rescale correction step
  static const char *
      szDoubleHarmonicSupportHistogramName;///< the name and title for double harmonic method support histograms
  static const char
      *szCorrelationsSupportHistogramName;///< the name and title for correlations method support histograms
  static const char
      *szTwistCorrectedQnVectorName;///< the name of the Qn vector after applying the twist correction
  static const char
      *szRescaleCorrectedQnVectorName;             ///< the name of the Qn vector after applying the rescale correction
  static const char *szQANotValidatedHistogramName;///< the name and title for bin not validated QA histograms
  static const char
      *szQATwistQnAverageHistogramName;///< the name and title for after twist Qn components average QA histograms
  static const char *
      szQARescaleQnAverageHistogramName;///< the name and title for after rescale Qn components average QA histograms
  std::unique_ptr<CorrectionProfileComponents>
      fDoubleHarmonicInputHistograms;//!<! the histogram with calibration information for the double harmonic method
  std::unique_ptr<CorrectionProfileComponents>
      fDoubleHarmonicCalibrationHistograms;//!<! the histogram for building calibration information for the doubel harmonic method
  std::unique_ptr<CorrectionProfile3DCorrelations>
      fCorrelationsInputHistograms;//!<! the histogram with calibration information for the correlations method
  std::unique_ptr<CorrectionProfile3DCorrelations>
      fCorrelationsCalibrationHistograms;//!<! the histogram for building calibration information for the correlations method
  std::unique_ptr<CorrectionHistogramSparse>
      fQANotValidatedBin;//!<! the histogram with non validated bin information
  std::unique_ptr<CorrectionProfileComponents>
      fQATwistQnAverageHistogram;//!<! the after twist correction step average Qn components QA histogram
  std::unique_ptr<CorrectionProfileComponents>
      fQARescaleQnAverageHistogram;//!<! the after rescale correction step average Qn components QA histogram

  Method fTwistAndRescaleMethod;          ///< the chosen method for extracting twist and rescale correction parameters
  Bool_t fApplyTwist;                     ///< apply the twist step
  Bool_t fApplyRescale;                   ///< apply the rescale step
  std::string fBDetectorConfigurationName;///< the name of the B detector configuration
  SubEvent *fSubEventB;                   ///< pointer to the B detector configuration
  std::string fCDetectorConfigurationName;///< the name of the C detector configuration
  SubEvent *fSubEventC;                   ///< pointer to the C detector configuration
  Int_t fMinNoOfEntriesToValidate;        ///< number of entries for bin content validation threshold

  std::unique_ptr<QVector> fTwistCorrectedQnVector;  //!<! twisted Qn vector
  std::unique_ptr<QVector> fRescaleCorrectedQnVector;//!<! rescaled Qn vector

  /// \cond CLASSIMP
  ClassDef(TwistAndRescale, 3);
  /// \endcond
};
}// namespace Qn
#endif// QNCORRECTIONS_QNVECTORTWISTANDRESCALE_H
