#ifndef QNCORRECTIONS_INPUTGAINEQUALIZATION_H
#define QNCORRECTIONS_INPUTGAINEQUALIZATION_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsInputGainEqualization.h
/// \brief Definition of the class that implements gain equalization of individual channels.
///
/// Gain equalization is applied on raw data from the involved detector.
/// Two procedures are implemented: average gain equalization and width equalization.
///
/// The average gain equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///        \mbox{M}'_{c,i} = \frac{\mbox{M}_{c,i}}{\langle\mbox{M}_{c}\rangle}
/// \f]
/// where  \f$\langle\mbox{M}_{c}\rangle\f$ is an average over events in a given event class
/// \f[
///        \langle\mbox{M}_{c}\rangle = \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}
/// \f]
///
/// The width equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///     \mbox{M}'_{c,i} = \mbox{A} + \mbox{B} \frac{\mbox{M}_{c,i} - \langle\mbox{M}_{c}\rangle}
///                               {\sigma_{{M}_{c}}}
/// \f]
/// with A and B are parameters that are the same for all channels and
/// \f$\sigma_{{M}_{c}}\f$ is the standard deviation of the signals of the channel \f$c\f$
/// for the considered event class
/// \f[
///        \sigma_{{M}_{c}} = \sqrt{
///          \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}^2_{c,i} -
///          \frac{1}{\mbox{N}^2_{ev}} \left(\Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}\right)^2}
/// \f]
///
/// The gain equalization is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper calibration histograms that will
/// provide the required averages per event class and channel.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the calibration histograms.

#include "CorrectionOnInputData.hpp"
#include "CorrectionProfileChannelizedIngress.hpp"
#include "CorrectionProfileChannelized.hpp"
#include "CorrectionHistogramChannelizedSparse.hpp"

namespace Qn {

/// \class QnCorrectionsInputGainEqualization
/// \brief Encapsulates the gain equalization on input data correction step
///
/// Gain equalization is applied on raw data from the involved detector.
/// Two procedures are implemented: average gain equalization and width equalization.
///
/// The average gain equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///        \mbox{M}'_{c,i} = \frac{\mbox{M}_{c,i}}{\langle\mbox{M}_{c}\rangle}
/// \f]
/// where  \f$\langle\mbox{M}_{c}\rangle\f$ is an average over events in a given event class
/// \f[
///        \langle\mbox{M}_{c}\rangle = \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}
/// \f]
///
/// The width equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///     \mbox{M}'_{c,i} = \mbox{A} + \mbox{B} \frac{\mbox{M}_{c,i} - \langle\mbox{M}_{c}\rangle}
///                               {\sigma_{{M}_{c}}}
/// \f]
/// with A and B are parameters that are the same for all channels and
/// \f$\sigma_{{M}_{c}}\f$ is the standard deviation of the signals of the channel \f$c\f$
/// for the considered event class
/// \f[
///        \sigma_{{M}_{c}} = \sqrt{
///          \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}^2_{c,i} -
///          \frac{1}{\mbox{N}^2_{ev}} \left(\Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}\right)^2}
/// \f]
/// At the class level A is known as the shift and B is known as the scale.
/// The class also provides support for group equalization where a group weight can be
/// extracted from the channels that conform a group or could be passed as hard coded group
/// weights at detector configuration definition.
///
/// The gain equalization is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper calibration histograms that will
/// provide the required averages per event class and channel.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the calibration histograms.

class GainEqualization : public CorrectionOnInputData {
 public:
  /// \enum QnGainEqualizationMethod
  /// \brief The class of the id of the supported gain equalization methods
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  ///
  enum class Method {
    NONE,         ///< \f$ \mbox{M'} = \mbox{M}\f$
    AVERAGE,    ///< \f$ \mbox{M}' = \frac{\mbox{M}}{\langle\mbox{M}\rangle} \f$
    WIDTH,      ///< \f$ \mbox{M}' = \mbox{A} + \mbox{B} \frac{\mbox{M} - \langle\mbox{M} \rangle}{\sigma_{{M}}} \f$
  };

  GainEqualization();
  ~GainEqualization() = default;
  GainEqualization(const GainEqualization &other) : CorrectionOnInputData(other),
                                                    fEqualizationMethod(other.fEqualizationMethod),
                                                    fShift(other.fShift),
                                                    fScale(other.fScale),
                                                    fUseChannelGroupsWeights(other.fUseChannelGroupsWeights),
                                                    fHardCodedWeights(other.fHardCodedWeights),
                                                    fMinNoOfEntriesToValidate(other.fMinNoOfEntriesToValidate) {
  }

  virtual CorrectionOnInputData *MakeCopy() const { return dynamic_cast<CorrectionOnInputData *>(new GainEqualization(*this)); }

  /// Stores the passed equalization method
  /// \param method the desired equalization method
  void SetEqualizationMethod(Method method) { fEqualizationMethod = method; }

  /// Set the shift (A) width equalization parameter
  /// \param shift the shift parameter value
  void SetShift(Float_t shift) { fShift = shift; }
  /// Set the scale (B) equalization parameter
  /// \param scale the scale parameter value
  void SetScale(Float_t scale) { fScale = scale; }
  /// Enable or disable the group weights extracted from channel multiplicity
  /// \param enable kTRUE / kFALSE for enable / disable it
  void SetUseChannelGroupsWeights(Bool_t enable) { fUseChannelGroupsWeights = enable; }
  /// Set the minimum number of entries for calibration histogram bin content validation
  /// \param nNoOfEntries the number of entries threshold
  void SetNoOfEntriesThreshold(Int_t nNoOfEntries) { fMinNoOfEntriesToValidate = nNoOfEntries; }

  /// Informs when the detector configuration has been attached to the framework manager
  /// Basically this allows interaction between the different framework sections at configuration time
  /// No action for input gain equalization
  virtual void AttachedToFrameworkManager() {}
  virtual void AttachInput(TList *list);
  virtual void CreateSupportQVectors();
  virtual void CreateCorrectionHistograms();
  virtual void AttachQAHistograms(TList *list);
  virtual void AttachNveQAHistograms(TList *list);

  virtual Bool_t ProcessCorrections();
  virtual Bool_t ProcessDataCollection();
  /// Clean the correction to accept a new event
  /// Does nothing for the time being
  virtual void ClearCorrectionStep() {}

 private:
  using State = Qn::CorrectionBase::State;
  static constexpr const unsigned int szPriority =
      CorrectionOnInputData::Priority::kGainEqualization; ///< the key of the correction step for ordering purpose
  static constexpr const Float_t
      fMinimumSignificantValue = 1e-6;     ///< the minimum value that will be considered as meaningful for processing
  static constexpr const char *szCorrectionName = "Gain equalization"; ///< the name of the correction step
  static constexpr const char *szSupportHistogramName = "Multiplicity"; ///< the name and title for support histograms
  static constexpr const char *szQAHistogramName = "QA Multiplicity"; ///< the name and title for QA histograms
  static constexpr const char
      *szQANotValidatedHistogramName = "GE NvE";  ///< the name and title for bin not validated QA histograms
  std::unique_ptr<CorrectionProfileChannelizedIngress>
      fInputHistograms; //!<! the histogram with calibration information
  std::unique_ptr<CorrectionProfileChannelized>
      fCalibrationHistograms; //!<! the histogram for building calibration information
  std::unique_ptr<CorrectionProfileChannelized>
      fQAMultiplicityBefore;  //!<! the channel multiplicity histogram before gain equalization
  std::unique_ptr<CorrectionProfileChannelized>
      fQAMultiplicityAfter;   //!<! the channel multiplicity histogram after gain equalization
  std::unique_ptr<CorrectionHistogramChannelizedSparse>
      fQANotValidatedBin;    //!<! the histogram with non validated bin information
  Method fEqualizationMethod = Method::NONE; ///< the selected equalization method
  Float_t fShift = 0.0;                               ///< the shift (A) parameter for width equalization
  Float_t fScale = 1.0;                               ///< the scale (B) parameter for width equalization
  Bool_t fUseChannelGroupsWeights = false;              ///< use group weights extracted from channel multiplicity
  const Float_t
      *fHardCodedWeights = nullptr;             //!<! group hard coded weights stored in the detector configuration
  Int_t fMinNoOfEntriesToValidate = 2;              ///< number of entries for bin content validation threshold

/// \cond CLASSIMP
 ClassDef(GainEqualization, 2);
/// \endcond
};
}
#endif // QNCORRECTIONS_INPUTGAINEQUALIZATION_H
