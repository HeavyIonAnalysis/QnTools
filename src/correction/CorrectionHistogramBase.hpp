#ifndef QNCORRECTIONS_HISTOGRAMSBASE_H
#define QNCORRECTIONS_HISTOGRAMSBASE_H

/// \file QnCorrectionsHistogramBase.h
/// \brief Multidimensional profile histograms base class for the Q vector correction framework

#include <THn.h>
#include "CorrectionAxisSet.hpp"
namespace Qn {
/// \class QnCorrectionsHistogramBase
/// \brief Base class for the Q vector correction histograms
///
/// Basically stores the set of variables that identify
/// the different event classes the involved histograms
/// are storing information about. It also stores (in its
/// parent) the base name and base title for the different
/// histograms it will encapsulate.
///
/// The passed at construction option parameter is the option for the
/// the computation of the  errors in the descendant profiles. Possible
/// values for the options are:
///
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
///
/// The encapsulated bin axes values provide an efficient
/// runtime storage for computing bin numbers.
///
/// Provides the interface for the whole set of histogram
/// classes providing error information that helps debugging.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 11, 2016
class CorrectionHistogramBase {
  /// \typedef QnCorrectionHistogramErrorMode
  /// \brief The type of bin errors supported by the framework histograms.
 public:
  enum class ErrorMode{
    MEAN,  ///< the bin errors are the standard error on the mean
    SPREAD ///< the bin errors are the standard deviation
  };
 public:
  CorrectionHistogramBase() = default;
  CorrectionHistogramBase(std::string name, std::string title, const CorrectionAxisSet &ecvs, ErrorMode mode);
  CorrectionHistogramBase(std::string name, std::string title, const CorrectionAxisSet &ecvs);
  virtual ~CorrectionHistogramBase();
  /// Set the minimum number of entries needed to validate the bin content
  /// \param nNoOfEntries the number of entries threshold
  virtual void SetNoOfEntriesThreshold(Int_t nNoOfEntries) { fMinNoOfEntriesToValidate = nNoOfEntries; }
  std::string GetName() const { return fName;}
  std::string GetTitle() const { return fTitle;}

 protected:
  void FillBinAxesValues(Int_t chgrpId = -1);
  THnF *DivideTHnF(THnF *values, THnI *entries, THnC *valid = nullptr);
  void CopyTHnF(THnF *hDest, THnF *hSource, Int_t *binsArray);
  void CopyTHnFDimension(THnF *hDest, THnF *hSource, Int_t *binsArray, Int_t dimension);

  std::string fName;
  std::string fTitle;
  CorrectionAxisSet fEventClassVariables;  //!<! The variables set that determines the event classes
  Double_t *fBinAxesValues = nullptr;                                  //!<! Runtime place holder for computing bin number
  ErrorMode fErrorMode = ErrorMode::MEAN;                 //!<! The error type for the current instance
  Int_t fMinNoOfEntriesToValidate = nDefaultMinNoOfEntriesValidated;     ///< the minimum number of entries for validating a bin content
  /// \cond CLASSIMP
 ClassDef(CorrectionHistogramBase, 2);
  /// \endcond
  static const char *szChannelAxisTitle;                 ///< The title for the channel extra axis
  static const char *szGroupAxisTitle;                   ///< The title for the channel group extra axis
  static const char *szGroupHistoPrefix;                 ///< The prefix for the name of the group histograms
  static const char *szEntriesHistoSuffix;               ///< The suffix for the name of the entries histograms
  static const char *szXComponentSuffix;                 ///< The suffix for the name of X component histograms
  static const char *szYComponentSuffix;                 ///< The suffix for the name of Y component histograms
  static const char *szXXCorrelationComponentSuffix;     ///< The suffix for the name of XX correlation component histograms
  static const char *szXYCorrelationComponentSuffix;     ///< The suffix for the name of XY correlation component histograms
  static const char *szYXCorrelationComponentSuffix;     ///< The suffix for the name of YX correlation component histograms
  static const char *szYYCorrelationComponentSuffix;     ///< The suffix for the name of YY correlation component histograms
  static const Int_t nMaxHarmonicNumberSupported;        ///< The maximum external harmonic number the framework support
  static const UInt_t harmonicNumberMask[];              ///< Mask for each external harmonic number
  static const UInt_t correlationXXmask;                 ///< Maks for XX correlation component
  static const UInt_t correlationXYmask;                 ///< Maks for XY correlation component
  static const UInt_t correlationYXmask;                 ///< Maks for YX correlation component
  static const UInt_t correlationYYmask;                 ///< Maks for YY correlation component
  constexpr static const Int_t nDefaultMinNoOfEntriesValidated =2; ///< The default minimum number of entries for validating a bin content
};

/// Fills the axes values for the current passed variable container
///
/// Core of the GetBin members. Stores the current values of the involved
/// variables in the internal place holder. Space is prepared for potential
/// channel or group id.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param chgrpId additional optional channel or group Id
inline void CorrectionHistogramBase::FillBinAxesValues(Int_t chgrpId) {
  unsigned int ivar = 0;
    for (const auto &var : fEventClassVariables) {
    fBinAxesValues[ivar] = var.GetValue();//variableContainer[var.GetId()];
    ++ivar;
  }
  fBinAxesValues[fEventClassVariables.GetSize()] = chgrpId;
}

}
#endif
