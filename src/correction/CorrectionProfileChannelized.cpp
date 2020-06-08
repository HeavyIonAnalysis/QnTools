/// \file QnCorrectionsProfileChannelized.cxx
/// \brief Implementation of multidimensional channelized profile class

#include "TList.h"

#include "CorrectionAxisSet.hpp"
#include "CorrectionProfileChannelized.hpp"

/// \cond CLASSIMP
ClassImp(Qn::CorrectionProfileChannelized);
/// \endcond
namespace Qn {
/// Normal constructor
///
/// Stores the set of variables that identify the
/// different event classes passing them to its parent
/// and prepares the object for actual histogram
/// creation or attachment
///
/// \param name base for the name of the histograms
/// \param title base for the title of the histograms
/// \param ecvs the event classes variables set
/// \param nNoOfChannels the number of channels associated
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
CorrectionProfileChannelized::CorrectionProfileChannelized(std::string name,
                                                           std::string title,
                                                           const CorrectionAxisSet &ecvs,
                                                           Int_t nNoOfChannels,
                                                           ErrorMode mode) : CorrectionHistogramBase(name, title, ecvs, mode),
                                                                             fNoOfChannels(nNoOfChannels) {}

CorrectionProfileChannelized::CorrectionProfileChannelized(std::string name,
                                                           const CorrectionAxisSet &ecvs,
                                                           Int_t nNoOfChannels,
                                                           ErrorMode mode) : CorrectionHistogramBase(name, name, ecvs, mode),
                                                                             fNoOfChannels(nNoOfChannels) {}

/// Default destructor
/// Releases the memory taken
CorrectionProfileChannelized::~CorrectionProfileChannelized() {
  delete[] fUsedChannel;
  delete[] fChannelGroup;
  delete[] fChannelMap;
}

/// Creates the support histograms for the profile function
///
/// Based in the event classes variables set in the parent class
/// and the channel information passed as parameters
/// the values and entries multidimensional histograms are
/// created.
///
/// Both histograms are added to the passed histogram list
///
/// The actual number of channels is stored and a mask from
/// external channel number to histogram channel number. If
/// bUsedChannel is nullptr all channels
/// within fNoOfChannels are assigned to this profile.
/// If nChannelGroup is nullptr all channels assigned to this
/// profile are allocated to the same group.
/// \param histogramList list where the histograms have to be added
/// \param bUsedChannel array of booleans one per each channel
/// \param nChannelGroup array of group number for each channel
/// \return true if properly created
Bool_t CorrectionProfileChannelized::CreateProfileHistograms(TList *histogramList,
                                                             const Bool_t *bUsedChannel,
                                                             const Int_t *nChannelGroup) {
  /* let's build the histograms names and titles */
  TString histoName = GetName();
  TString histoTitle = GetTitle();
  TString entriesHistoName = GetName();
  entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle();
  entriesHistoTitle += szEntriesHistoSuffix;
  /* we open space for channel variable as well */
  Int_t nVariables = fEventClassVariables.GetSize();
  Double_t *minvals = new Double_t[nVariables + 1];
  Double_t *maxvals = new Double_t[nVariables + 1];
  Int_t *nbins = new Int_t[nVariables + 1];
  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins, minvals, maxvals);
  /* lets consider now the channel information */
  fUsedChannel = new Bool_t[fNoOfChannels];
  fChannelGroup = new Int_t[fNoOfChannels];
  fChannelMap = new Int_t[fNoOfChannels];
  fActualNoOfChannels = 0;
  for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
    if (bUsedChannel) {
      fUsedChannel[ixChannel] = bUsedChannel[ixChannel];
    } else {
      fUsedChannel[ixChannel] = kTRUE;
    }
    if (nChannelGroup) {
      fChannelGroup[ixChannel] = nChannelGroup[ixChannel];
    } else {
      fChannelGroup[ixChannel] = 0;
    }

    if (fUsedChannel[ixChannel]) {
      fChannelMap[ixChannel] = fActualNoOfChannels;
      fActualNoOfChannels++;
    }
  }
  /* There will be a wrong external view of the channel number especially */
  /* manifested when there are holes in the channel assignment */
  /* so, lets complete the dimension information */
  /* WARNING: be aware that ROOT does not keep label information when projecting THn */
  minvals[nVariables] = -0.5;
  maxvals[nVariables] = -0.5 + fActualNoOfChannels;
  nbins[nVariables] = fActualNoOfChannels;
  /* create the values and entries multidimensional histograms */
  fValues = new THnF((const char *) histoName, (const char *) histoTitle, nVariables + 1, nbins, minvals, maxvals);
  fEntries = new THnI((const char *) entriesHistoName,
                      (const char *) entriesHistoTitle,
                      nVariables + 1,
                      nbins,
                      minvals,
                      maxvals);
  /* now let's set the proper binning and label on each axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fEntries->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fEntries->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
  }
  /* and now the channel axis */
  fValues->GetAxis(nVariables)->SetTitle(szChannelAxisTitle);
  fEntries->GetAxis(nVariables)->SetTitle(szChannelAxisTitle);
  /* and now set the proper channel labels if needed */
  if (fActualNoOfChannels != fNoOfChannels) {
    for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
      if (fUsedChannel[ixChannel]) {
        fValues->GetAxis(nVariables)->SetBinLabel(fChannelMap[ixChannel] + 1, Form("%d", ixChannel));
        fEntries->GetAxis(nVariables)->SetBinLabel(fChannelMap[ixChannel] + 1, Form("%d", ixChannel));
      }
    }
  }
  fValues->Sumw2();
  histogramList->Add(fValues);
  histogramList->Add(fEntries);
  delete[] minvals;
  delete[] maxvals;
  delete[] nbins;
  return kTRUE;
}

/// Get the bin number for the current variable content and passed channel
///
/// The bin number identifies the event class the current
/// variable content points to under the passed channel.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \return the associated bin to the current variables content
Long64_t CorrectionProfileChannelized::GetBin(Int_t nChannel) {
  FillBinAxesValues(fChannelMap[nChannel]);
  /* store the channel number */
  return fEntries->GetBin(fBinAxesValues);
}

/// Check the validity of the content of the passed bin
/// If the number of entries is lower
/// than the minimum number of entries to validate it
/// the bin content is not considered valid and kFALSE is returned,
/// otherwise kTRUE is returned
/// \param bin the bin to check its content validity
/// \return kTRUE if the content is valid kFALSE otherwise
Bool_t CorrectionProfileChannelized::BinContentValidated(Long64_t bin) {
  auto nEntries = Int_t(fEntries->GetBinContent(bin));
  return nEntries >= fMinNoOfEntriesToValidate;
}

/// Get the bin content for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// is requested. If the bin content is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileChannelized::GetBinContent(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fValues->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content error
Float_t CorrectionProfileChannelized::GetBinError(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fValues->GetBinContent(bin);
    Float_t error2 = fValues->GetBinError2(bin);
    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN:
        return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD:
        return serror;
      default:
        return 0.0;
    }
  }
}

/// Fills the histogram
///
/// The involved bin is computed according to the current variables
/// content and the passed external channel number. The bin is then
/// increased by the given weight and the entries also increased properly.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \param weight the increment in the bin content
void CorrectionProfileChannelized::Fill(Int_t nChannel, Float_t weight) {
  Double_t nEntries = fValues->GetEntries();
  FillBinAxesValues(fChannelMap[nChannel]);
  fValues->Fill(fBinAxesValues, weight);
  fValues->SetEntries(nEntries + 1);
  fEntries->Fill(fBinAxesValues, 1.0);
}
}// namespace Qn
