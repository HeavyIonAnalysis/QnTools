/// file QnCorrectionsHistogramBase.cxx
/// \brief Implementation of the multidimensional profile base class
#include "CorrectionHistogramBase.hpp"
#include "TList.h"
#include <utility>

namespace Qn {
const char *CorrectionHistogramBase::szChannelAxisTitle = "Channel number";
const char *CorrectionHistogramBase::szGroupAxisTitle = "Channels group";
const char *CorrectionHistogramBase::szGroupHistoPrefix = "Group";
const char *CorrectionHistogramBase::szEntriesHistoSuffix = "_entries";
const char *CorrectionHistogramBase::szXComponentSuffix = "X";
const char *CorrectionHistogramBase::szYComponentSuffix = "Y";
const char *CorrectionHistogramBase::szXXCorrelationComponentSuffix = "XX";
const char *CorrectionHistogramBase::szXYCorrelationComponentSuffix = "XY";
const char *CorrectionHistogramBase::szYXCorrelationComponentSuffix = "YX";
const char *CorrectionHistogramBase::szYYCorrelationComponentSuffix = "YY";
const Int_t CorrectionHistogramBase::nMaxHarmonicNumberSupported = 15;
const UInt_t CorrectionHistogramBase::harmonicNumberMask[] =
    {0x0000, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080,
     0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000};
const UInt_t CorrectionHistogramBase::correlationXXmask = 0x0001;
const UInt_t CorrectionHistogramBase::correlationXYmask = 0x0002;
const UInt_t CorrectionHistogramBase::correlationYXmask = 0x0004;
const UInt_t CorrectionHistogramBase::correlationYYmask = 0x0008;

/// Default destructor
///
/// restores the taken memory for the bin axes values bank
CorrectionHistogramBase::~CorrectionHistogramBase() {
  delete[] fBinAxesValues;
}

/// Normal constructor
///
/// Basically stores the set of variables that identify
/// the different event classes the involved histograms
/// are storing information about
///
/// This is the base class. For simplicity and consistency
/// we leave open an extra variable storage that will be
/// used by the channelized histograms.
///
/// \param name base for the name of the histograms
/// \param title base for the title of the histograms
/// \param ecvs the event classes variables set
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
CorrectionHistogramBase::CorrectionHistogramBase(std::string name,
                                                 std::string title,
                                                 const CorrectionAxisSet &ecvs,
                                                 ErrorMode mode) : fName(std::move(name)),
                                                                   fTitle(std::move(title)),
                                                                   fEventClassVariables(ecvs),
                                                                   fBinAxesValues(new Double_t[fEventClassVariables.GetSize() + 1]),
                                                                   fErrorMode(mode) {}

CorrectionHistogramBase::CorrectionHistogramBase(std::string name,
                                                 std::string title,
                                                 const CorrectionAxisSet &ecvs) : fName(std::move(name)),
                                                                                  fTitle(std::move(title)),
                                                                                  fEventClassVariables(ecvs),
                                                                                  fBinAxesValues(new Double_t[fEventClassVariables.GetSize() + 1]) {}

/// Divide two THn histograms
///
/// Creates a value / error multidimensional histogram from
/// a values and entries multidimensional histograms.
/// The validation histogram is filled according to entries threshold value.
/// \param hValues the values multidimensional histogram
/// \param hEntries the entries multidimensional histogram
/// \param hValid optional multidimensional histogram where validation information is stored
/// \return the values / error multidimensional histogram
THnF *CorrectionHistogramBase::DivideTHnF(THnF *hValues, THnI *hEntries, THnC *hValid) {
  THnF *hResult = (THnF *) THn::CreateHn(hValues->GetName(), hValues->GetTitle(), hValues);
  Double_t value, error2;
  Int_t nEntries;
  Int_t nNotValidatedBins = 0;
  for (Long64_t bin = 0; bin < hResult->GetNbins(); bin++) {
    value = hValues->GetBinContent(bin);
    nEntries = Int_t(hEntries->GetBinContent(bin));
    error2 = hValues->GetBinError2(bin);
    if (nEntries < fMinNoOfEntriesToValidate) {
      /* bin content not validated */
      hResult->SetBinContent(bin, 0.0);
      hResult->SetBinError(bin, 0.0);
      if (hValid) hValid->SetBinContent(bin, 0.0);
      if (value != 0.0) {
        nNotValidatedBins++;
      }
    } else {
      Double_t average = value / nEntries;
      Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
      switch (fErrorMode) {
        case ErrorMode::MEAN:
          /* standard error on the mean of the bin values */
          hResult->SetBinContent(bin, average);
          hResult->SetBinError(bin, serror / TMath::Sqrt(nEntries));
          break;
        case ErrorMode::SPREAD:
          /* standard deviation of the bin values */
          hResult->SetBinContent(bin, average);
          hResult->SetBinError(bin, serror);
          break;
      }
      if (hValid != nullptr) hValid->SetBinContent(bin, 1.0);
    }
    hResult->SetEntries(hValues->GetEntries());
  }
  return hResult;
}

/// Starts the copy of two THnF histograms.
///
/// Source should not have channel/group structure (axis)
/// while dest should. The channel/group involved in dest
/// should be stored in binsArray[nVariables].
/// \param hDest the histogram that will receive the copy
/// \param hSource the histogram to copy
/// \param binsArray the array to build the bin numbers on each dimension
void CorrectionHistogramBase::CopyTHnF(THnF *hDest, THnF *hSource, Int_t *binsArray) {
  CopyTHnFDimension(hDest, hSource, binsArray, 0);
}

/// Process a dimension in the copy of two THnF histograms process.
///
/// Source should not have channel/group structure (axis)
/// while dest should. The channel/group involved in dest
/// should be stored in binsArray[nVariables].
/// It is called in a recursive way until the number of
/// event class variables (dimensions in the source) is
/// exhausted then, the corresponding bins are copied.
/// \param hDest the histogram that will receive the copy
/// \param hSource the histogram to copy
/// \param binsArray the array to build the bin numbers on each dimension
/// \param dimension the current dimension being handled
void CorrectionHistogramBase::CopyTHnFDimension(THnF *hDest, THnF *hSource, Int_t *binsArray, Int_t dimension) {
  /* are all variables settled */
  if ((unsigned long) dimension < fEventClassVariables.GetSize()) {
    /* no then, scan this dimension and move to the next one */
    for (Long64_t bin = 0; bin < hSource->GetAxis(dimension)->GetNbins(); bin++) {
      binsArray[dimension] = bin + 1;
      CopyTHnFDimension(hDest, hSource, binsArray, dimension + 1);
    }
  } else {
    /* all variables have a new bin configuration, let's do the actual copy */
    Double_t value = hSource->GetBinContent(binsArray);
    Double_t error = hSource->GetBinError(binsArray);
    hDest->SetBinContent(binsArray, value);
    hDest->SetBinError(binsArray, error);
  }
}
}// namespace Qn
