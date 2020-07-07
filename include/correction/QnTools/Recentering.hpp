#ifndef QNCORRECTIONS_QNVECTORRECENTERING_H
#define QNCORRECTIONS_QNVECTORRECENTERING_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsQnVectorRecentering.h
/// \brief Definition of the class that implements recentering and width equalization on Q vectors.
///
/// Recentering and weight equalization is applied on ongoing Q vector from the involved detector
/// configuration. It is supposed that such a Q vector has already been built after input data
/// equalization and is with proper normalization for the recentering correction being applied.
/// The recentering correction is always applied while the with equalization is user configurable.
///
/// The recentering is applied according to:
/// \f[
///        Q' = Q - {\langle Q \rangle}
/// \f]
/// where  \f$\langle Q \rangle\f$ is an average over events in a given event class
/// \f[
///        \langle Q \rangle = \frac{1}{\mbox{N}_{ev}} \sum_{i}^{\mbox{N}_{ev}} Q_i
/// \f]
///
/// The recentering and width equalization is applied according to:
/// \f[
///     Q' = \frac{Q- \langle Q \rangle}{\sigma_Q}
/// \f]
/// where \f$ \sigma_Q \f$ is the standard deviation of the \f$ Q \f$ values
/// for the considered event class.
/// \f[
///        \sigma_Q = \sqrt{
///          \frac{1}{\mbox{N}_{ev}} \sum_{i}^{\mbox{N}_{ev}} Q^2_i -
///          \frac{1}{\mbox{N}^2_{ev}} \left(\sum_{i}^{\mbox{N}_{ev}} Q_i \right)^2}
/// \f]
///
/// Recentering (and width equalization) is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper correction histograms that will
/// provide the required averages per event class.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the correction histograms.
///
/// Correction and data collecting during calibration is performed for all harmonics
/// defined within the involved detector configuration

#include "CorrectionOnQnVector.hpp"
namespace Qn {
/// \class QnCorrectionsQnVectorRecentering
/// \brief Encapsulates recentering and width equalization on Q vector
///
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Apr 01, 2016
///
/// Recentering and weight equalization is applied on ongoing Q vector from the involved detector
/// configuration.
/// The recentering correction is always applied while the with equalization is user configurable.
///
/// The recentering is applied according to:
/// \f[
///        Q' = Q - {\langle Q \rangle}
/// \f]
/// where  \f$\langle Q \rangle\f$ is an average over events in a given event class
/// \f[
///        \langle Q \rangle = \frac{1}{\mbox{N}_{ev}} \sum_{i}^{\mbox{N}_{ev}} Q_i
/// \f]
///
/// The recentering and width equalization is applied according to:
/// \f[
///     Q' = \frac{Q- \langle Q \rangle}{\sigma_Q}
/// \f]
/// where \f$ \sigma_Q \f$ is the standard deviation of the \f$ Q \f$ values
/// for the considered event class.
/// \f[
///        \sigma_Q = \sqrt{
///          \frac{1}{\mbox{N}_{ev}} \sum_{i}^{\mbox{N}_{ev}} Q^2_i -
///          \frac{1}{\mbox{N}^2_{ev}} \left(\sum_{i}^{\mbox{N}_{ev}} Q_i \right)^2}
/// \f]
///
/// Recentering (and width equalization) is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper correction histograms that will
/// provide the required averages per event class.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the correction histograms.
///
/// Correction and data collecting during calibration is performed for all harmonics
/// defined within the involved detector configuration

class CorrectionHistogramSparse;

class Recentering : public CorrectionOnQnVector {
 public:
  Recentering();
  virtual ~Recentering() = default;
  Recentering(const Recentering &other) : CorrectionOnQnVector(other),
                                          fApplyWidthEqualization(other.fApplyWidthEqualization),
                                          fMinNoOfEntriesToValidate(other.fMinNoOfEntriesToValidate) {
  }

  virtual CorrectionOnQnVector *MakeCopy() const { return dynamic_cast<CorrectionOnQnVector *>(new Recentering(*this)); }

  /// Controls if width equalization step shall be additionally applied
  /// \param apply kTRUE for applying the width equalization step
  void SetApplyWidthEqualization(Bool_t apply) { fApplyWidthEqualization = apply; }
  /// Set the minimum number of entries for calibration histogram bin content validation
  /// \param nNoOfEntries the number of entries threshold
  void SetNoOfEntriesThreshold(Int_t nNoOfEntries) { fMinNoOfEntriesToValidate = nNoOfEntries; }
  /// Informs when the detector configuration has been attached to the framework manager
  /// Basically this allows interaction between the different framework sections at configuration time
  /// No action for Qn vector recentering
  virtual void AttachInput(TList *list);
  /// Perform after calibration histograms attach actions
  /// It is used to inform the different correction step that
  /// all conditions for running the network are in place so
  /// it is time to check if their requirements are satisfied
  /// Does nothing for the time being
  virtual void AfterInputAttachAction() {}
  virtual void CreateSupportQVectors();
  virtual void CreateCorrectionHistograms();
  virtual void AttachQAHistograms(TList *list);
  virtual void AttachNveQAHistograms(TList *list);
  virtual Bool_t ProcessCorrections();
  virtual Bool_t ProcessDataCollection();
  virtual void ClearCorrectionStep();

 private:
  using State = Qn::CorrectionBase::State;
  static constexpr const unsigned int
      szPriority = CorrectionOnQnVector::Step::kRecentering; ///< the key of the correction step for ordering purpose
  static const Int_t fDefaultMinNoOfEntries;         ///< the minimum number of entries for bin content validation
  static const char *szCorrectionName;               ///< the name of the correction step
  static const char *szSupportHistogramName;         ///< the name and title for support histograms
  static const char *szCorrectedQnVectorName;        ///< the name of the Qn vector after applying the correction
  static const char *szQANotValidatedHistogramName;  ///< the name and title for bin not validated QA histograms
  static const char *szQAQnAverageHistogramName;     ///< the name and title for Qn components average QA histograms
  std::unique_ptr<CorrectionProfileComponents> fInputHistograms;       //!<! the histogram with calibration information
  std::unique_ptr<CorrectionProfileComponents>
      fCalibrationHistograms; //!<! the histogram for building calibration information
  std::unique_ptr<CorrectionHistogramSparse>
      fQANotValidatedBin;     //!<! the histogram with non validated bin information
  std::unique_ptr<CorrectionProfileComponents>
      fQAQnAverageHistogram;  //!<! the after correction step average Qn components QA histogram
  Bool_t fApplyWidthEqualization;               ///< apply the width equalization step
  Int_t fMinNoOfEntriesToValidate;              ///< number of entries for bin content validation threshold

/// \cond CLASSIMP
 ClassDef(Recentering, 3);
/// \endcond
};
}
#endif // QNCORRECTIONS_QNVECTORRECENTERING_H
