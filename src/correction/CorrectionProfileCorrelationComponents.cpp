/// \file QnCorrectionsProfileCorrelationComponents.cxx
/// \brief Implementation of the multidimensional correlation components based set of profiles

#include "TList.h"

#include "CorrectionAxisSet.hpp"
#include "CorrectionProfileCorrelationComponents.hpp"

/// \cond CLASSIMP
ClassImp(Qn::CorrectionProfileCorrelationComponents);
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
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
CorrectionProfileCorrelationComponents::CorrectionProfileCorrelationComponents(std::string name,
                                                                               std::string title,
                                                                               const CorrectionAxisSet &ecvs,
                                                                               ErrorMode mode) : CorrectionHistogramBase(name, title, ecvs, mode) {}

CorrectionProfileCorrelationComponents::CorrectionProfileCorrelationComponents(std::string name,
                                                                               const CorrectionAxisSet &ecvs,
                                                                               ErrorMode mode) : CorrectionHistogramBase(name, name, ecvs, mode) {}

/// Creates the XX, XY, YX, YY correlation components support histograms
/// for the profile function
///
/// Based in the event classes variables set in the parent class
/// the values and entries multidimensional histograms are
/// created.
///
/// Four values histograms are created, XX,XY, YX and YY.
/// The fully filled condition is computed and stored
///
/// The whole set of histograms are added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \return true if properly created
Bool_t CorrectionProfileCorrelationComponents::CreateCorrelationComponentsProfileHistograms(TList *histogramList) {
  /* let's build the histograms names and titles */
  TString histoXXName = GetName();
  histoXXName += szXXCorrelationComponentSuffix;
  TString histoXYName = GetName();
  histoXYName += szXYCorrelationComponentSuffix;
  TString histoYXName = GetName();
  histoYXName += szYXCorrelationComponentSuffix;
  TString histoYYName = GetName();
  histoYYName += szYYCorrelationComponentSuffix;
  TString histoXXTitle = GetTitle();
  histoXXTitle += szXXCorrelationComponentSuffix;
  TString histoXYTitle = GetTitle();
  histoXYTitle += szXYCorrelationComponentSuffix;
  TString histoYXTitle = GetTitle();
  histoYXTitle += szYXCorrelationComponentSuffix;
  TString histoYYTitle = GetTitle();
  histoYYTitle += szYYCorrelationComponentSuffix;
  TString entriesHistoName = GetName();
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle();
  entriesHistoTitle += szXXCorrelationComponentSuffix;
  entriesHistoTitle += szXYCorrelationComponentSuffix;
  entriesHistoTitle += szYXCorrelationComponentSuffix;
  entriesHistoTitle += szYYCorrelationComponentSuffix;
  entriesHistoTitle += szEntriesHistoSuffix;
  /* now prepare the construction of the histograms */
  Int_t nVariables = fEventClassVariables.GetSize();
  auto *minvals = new Double_t[nVariables];
  auto *maxvals = new Double_t[nVariables];
  auto *nbins = new Int_t[nVariables];
  TString sVariableLabels = "";
  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins, minvals, maxvals);
  /* create the values multidimensional histograms */
  fXXValues = new THnF((const char *) histoXXName, (const char *) histoXXTitle,
                       nVariables, nbins, minvals, maxvals);
  fXYValues = new THnF((const char *) histoXYName, (const char *) histoXYTitle,
                       nVariables, nbins, minvals, maxvals);
  fYXValues = new THnF((const char *) histoYXName, (const char *) histoYXTitle,
                       nVariables, nbins, minvals, maxvals);
  fYYValues = new THnF((const char *) histoYYName, (const char *) histoYYTitle,
                       nVariables, nbins, minvals, maxvals);
  /* now let's set the proper binning and label on each axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fXXValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fXXValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fXXValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fXXValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fXYValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fXYValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fXYValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fXYValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fYXValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fYXValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fYXValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fYXValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fYYValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fYYValues->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fYYValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fYYValues->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
  }
  /* ask for square sum accumulation */
  fXXValues->Sumw2();
  fXYValues->Sumw2();
  fYXValues->Sumw2();
  fYYValues->Sumw2();
  /* and finally add the histograms to the list */
  histogramList->Add(fXXValues);
  histogramList->Add(fXYValues);
  histogramList->Add(fYXValues);
  histogramList->Add(fYYValues);
  /* and store the fully filled condition */
  fXXXYYXYYFillMask = 0x0000;
  fFullFilled = correlationXXmask | correlationXYmask | correlationYXmask | correlationYYmask;
  /* create the entries multidimensional histogram */
  fEntries =
      new THnI((const char *) entriesHistoName, (const char *) entriesHistoTitle, nVariables, nbins, minvals, maxvals);
  /* now let's set the proper binning and label on each entries histogram axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fEntries->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fEntries->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
  }
  /* and finally add the entries histogram to the list */
  histogramList->Add(fEntries);
  delete[] minvals;
  delete[] maxvals;
  delete[] nbins;
  return kTRUE;
}

/// Attaches existing histograms as the support histograms for XX, XY, YX, YY
/// correlation component of the profile function for different harmonics
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// The harmonic map is inferred from the found histograms within the list
/// that match the naming scheme.
///
/// \param histogramList list where the histograms have to be located
/// \return true if properly attached else false
Bool_t CorrectionProfileCorrelationComponents::AttachHistograms(TList *histogramList) {
  /* let's build the histograms names */
  TString histoXXName = GetName();
  histoXXName += szXXCorrelationComponentSuffix;
  TString histoXYName = GetName();
  histoXYName += szXYCorrelationComponentSuffix;
  TString histoYXName = GetName();
  histoYXName += szYXCorrelationComponentSuffix;
  TString histoYYName = GetName();
  histoYYName += szYYCorrelationComponentSuffix;
  TString entriesHistoName = GetName();
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  /* initialize. Remember we don't own the histograms */
  fEntries = nullptr;
  fXXValues = nullptr;
  fXYValues = nullptr;
  fYXValues = nullptr;
  fYYValues = nullptr;
  fXXXYYXYYFillMask = 0x0000;
  fFullFilled = 0x0000;
  fEntries = (THnI *) histogramList->FindObject((const char *) entriesHistoName);
  if (fEntries != nullptr && fEntries->GetEntries() != 0) {
    /* search the values multidimensional histograms */
    fXXValues = (THnF *) histogramList->FindObject((const char *) histoXXName);
    fXYValues = (THnF *) histogramList->FindObject((const char *) histoXYName);
    fYXValues = (THnF *) histogramList->FindObject((const char *) histoYXName);
    fYYValues = (THnF *) histogramList->FindObject((const char *) histoYYName);
    /* and update the fully filled condition whether applicable */
    if ((fXXValues != nullptr) && (fXYValues != nullptr) && (fYXValues != nullptr) && (fYYValues != nullptr))
      fFullFilled = correlationXXmask | correlationXYmask | correlationYXmask | correlationYYmask;
  } else
    return kFALSE;
  /* check that we actually got something */
  return fFullFilled != 0x0000;
}

/// Get the bin number for the current variable content
///
/// The bin number identifies the event class the current
/// variable content points to.
///
/// \param variableContainer the current variables content addressed by var Id
/// \return the associated bin to the current variables content
Long64_t CorrectionProfileCorrelationComponents::GetBin() {
  FillBinAxesValues();
  return fEntries->GetBin(fBinAxesValues);
}

/// Check the validity of the content of the passed bin
/// If the number of entries is lower
/// than the minimum number of entries to validate it
/// the bin content is not considered valid and kFALSE is returned,
/// otherwise kTRUE is returned
/// \param bin the bin to check its content validity
/// \return kTRUE if the content is valid kFALSE otherwise
Bool_t CorrectionProfileCorrelationComponents::BinContentValidated(Long64_t bin) {
  auto nEntries = Int_t(fEntries->GetBinContent(bin));
  return nEntries >= fMinNoOfEntriesToValidate;
}

/// Get the XX correlation component bin content.
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileCorrelationComponents::GetXXBinContent(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXXValues->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XY correlation component bin content.
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileCorrelationComponents::GetXYBinContent(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXYValues->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YX correlation component bin content.
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileCorrelationComponents::GetYXBinContent(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYXValues->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YY correlation component bin content.
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileCorrelationComponents::GetYYBinContent(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYYValues->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XX correlation component bin content error for the passed bin number.
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfileCorrelationComponents::GetXXBinError(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXXValues->GetBinContent(bin);
    Float_t error2 = fXXValues->GetBinError2(bin);

    Double_t average = values / Float_t(nEntries);
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / Float_t(nEntries) - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN:
        /* standard error on the mean of the bin values */
        return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD:
        /* standard deviation of the bin values */
        return serror;
      default:
        return 0.0;
    }
  }
}

/// Get the XY correlation component bin content error for the passed bin number.
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfileCorrelationComponents::GetXYBinError(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXYValues->GetBinContent(bin);
    Float_t error2 = fXYValues->GetBinError2(bin);

    Double_t average = values / Float_t(nEntries);
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / Float_t(nEntries) - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN:
        /* standard error on the mean of the bin values */
        return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD:
        /* standard deviation of the bin values */
        return serror;
      default:
        return 0.0;
    }
  }
}

/// Get the YX correlation component bin content error for the passed bin number.
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfileCorrelationComponents::GetYXBinError(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYXValues->GetBinContent(bin);
    Float_t error2 = fYXValues->GetBinError2(bin);

    Double_t average = values / Float_t(nEntries);
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / Float_t(nEntries) - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN:
        /* standard error on the mean of the bin values */
        return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD:
        /* standard deviation of the bin values */
        return serror;
      default:
        return 0.0;
    }
  }
}

/// Get the YY correlation component bin content error for the passed bin number.
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin is not validated zero is returned.
///
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfileCorrelationComponents::GetYYBinError(Long64_t bin) {
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYYValues->GetBinContent(bin);
    Float_t error2 = fYYValues->GetBinError2(bin);

    Double_t average = values / Float_t(nEntries);
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / Float_t(nEntries) - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN:
        /* standard error on the mean of the bin values */
        return serror / TMath::Sqrt(nEntries);
        break;
      case ErrorMode::SPREAD:
        /* standard deviation of the bin values */
        return serror;
        break;
      default: return 0.0;
    }
  }
}

/// Fills the XX correlation component.
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills without a previous entries update.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void CorrectionProfileCorrelationComponents::FillXX(Float_t weight) {
  /* first the sanity checks */
  if (!(fXXXYYXYYFillMask & correlationXXmask)) {
    /* now it's safe to continue */
    /* keep total entries in fValues updated */
    Double_t nEntries = fXXValues->GetEntries();
    FillBinAxesValues();
    fXXValues->Fill(fBinAxesValues, weight);
    fXXValues->SetEntries(nEntries + 1);
    /* update fill mask */
    fXXXYYXYYFillMask |= correlationXXmask;
    /* now check if time for updating entries histogram */
    if (fXXXYYXYYFillMask != fFullFilled) return;
    /* update entries and reset the masks */
    fEntries->Fill(fBinAxesValues, 1.0);
    fXXXYYXYYFillMask = 0x0000;
  }
}

/// Fills the XY correlation component.
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills without a previous entries update.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void CorrectionProfileCorrelationComponents::FillXY(Float_t weight) {
  /* first the sanity checks */
  if (!(fXXXYYXYYFillMask & correlationXYmask)) {
    /* now it's safe to continue */
    /* keep total entries in fValues updated */
    Double_t nEntries = fXYValues->GetEntries();
    FillBinAxesValues();
    fXYValues->Fill(fBinAxesValues, weight);
    fXYValues->SetEntries(nEntries + 1);
    /* update fill mask */
    fXXXYYXYYFillMask |= correlationXYmask;
    /* now check if time for updating entries histogram */
    if (fXXXYYXYYFillMask != fFullFilled) return;
    /* update entries and reset the masks */
    fEntries->Fill(fBinAxesValues, 1.0);
    fXXXYYXYYFillMask = 0x0000;
  }
}

/// Fills the YX correlation component.
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills without a previous entries update.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void CorrectionProfileCorrelationComponents::FillYX(Float_t weight) {
  /* first the sanity checks */
  if (!(fXXXYYXYYFillMask & correlationYXmask)) {
    /* now it's safe to continue */
    /* keep total entries in fValues updated */
    Double_t nEntries = fYXValues->GetEntries();

    FillBinAxesValues();
    fYXValues->Fill(fBinAxesValues, weight);
    fYXValues->SetEntries(nEntries + 1);

    /* update fill mask */
    fXXXYYXYYFillMask |= correlationYXmask;

    /* now check if time for updating entries histogram */
    if (fXXXYYXYYFillMask != fFullFilled) return;
    /* update entries and reset the masks */
    fEntries->Fill(fBinAxesValues, 1.0);
    fXXXYYXYYFillMask = 0x0000;
  }
}

/// Fills the YY correlation component.
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries count is only updated if the whole set for the four components
/// has been already filled. A check is done for detecting consecutive
/// fills without a previous entries update.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void CorrectionProfileCorrelationComponents::FillYY(Float_t weight) {
  /* first the sanity checks */
  if (!(fXXXYYXYYFillMask & correlationYYmask)) {
    /* now it's safe to continue */
    /* keep total entries in fValues updated */
    Double_t nEntries = fYYValues->GetEntries();
    FillBinAxesValues();
    fYYValues->Fill(fBinAxesValues, weight);
    fYYValues->SetEntries(nEntries + 1);
    /* update harmonic fill mask */
    fXXXYYXYYFillMask |= correlationYYmask;
    /* now check if time for updating entries histogram */
    if (fXXXYYXYYFillMask != fFullFilled) return;
    /* update entries and reset the masks */
    fEntries->Fill(fBinAxesValues, 1.0);
    fXXXYYXYYFillMask = 0x0000;
  }
}
}// namespace Qn
