/// \file QnCorrectionsProfileComponents.cxx
/// \brief Implementation of the multidimensional component based set of profiles

#include "TList.h"

#include "CorrectionAxisSet.hpp"
#include "CorrectionProfileComponents.hpp"

/// \cond CLASSIMP
ClassImp(Qn::CorrectionProfileComponents);
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
CorrectionProfileComponents::CorrectionProfileComponents(std::string name,
                                                         std::string title,
                                                         const CorrectionAxisSet &ecvs,
                                                         ErrorMode mode) : CorrectionHistogramBase(name, title, ecvs, mode) {
}

CorrectionProfileComponents::CorrectionProfileComponents(std::string name,
                                                         const CorrectionAxisSet &ecvs,
                                                         ErrorMode mode) : CorrectionHistogramBase(name, name, ecvs, mode) {
}
/// Default destructor
///
/// Returns the only taken memory, the harmonic histograms storage,
/// the own histograms and other members are not own at destruction time
CorrectionProfileComponents::~CorrectionProfileComponents() {
  delete[] fXValues;
  delete[] fYValues;
}

/// Creates the X, Y components support histograms for the profile function
///
/// Based in the event classes variables set in the parent class
/// the values and entries multidimensional histograms are
/// created.
///
/// For each harmonic number two values histograms are created, X and Y.
/// The histograms are organized to support external harmonic number.
/// By default the external harmonic number is always considered to
/// start by one. If no map is passed as parameter the external harmonic
/// numbers are considered as: 1, 2, ..., nNoOfHarmonic.
/// If the user wants a different assignment he has to provide an
/// ordered map, for instance: four harmonics with external harmonic numbers
/// 2, 4, 6 and 8 will require nNoOfHarmonics = 4 and harmonicMap = [2, 4, 6, 8].
/// The fully filled condition is computed and stored
///
/// The whole set of histograms are added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \param nNoOfHarmonics the desired number of harmonics
/// \param harmonicMap ordered array with the external number of the harmonics
/// \return true if properly created
Bool_t CorrectionProfileComponents::CreateComponentsProfileHistograms(TList *histogramList,
                                                                      Int_t nNoOfHarmonics,
                                                                      Int_t *harmonicMap) {
  /* let's build the histograms names and titles */
  TString histoXName = GetName();
  histoXName += szXComponentSuffix;
  TString histoYName = GetName();
  histoYName += szYComponentSuffix;
  TString histoXTitle = GetTitle();
  histoXTitle += szXComponentSuffix;
  TString histoYTitle = GetTitle();
  histoYTitle += szYComponentSuffix;
  TString entriesHistoName = GetName();
  entriesHistoName += szXComponentSuffix;
  entriesHistoName += szYComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle();
  entriesHistoTitle += szXComponentSuffix;
  entriesHistoTitle += szYComponentSuffix;
  entriesHistoTitle += szEntriesHistoSuffix;

  /* check whether within the supported harmonic range */
  Int_t nHigherHarmonic = nNoOfHarmonics;
  if (harmonicMap) {
    nHigherHarmonic = harmonicMap[nNoOfHarmonics - 1];
  }
  if (nMaxHarmonicNumberSupported < nHigherHarmonic) {
  }

  /* let's support the external harmonic number map */
  /* external harmonic number will always start from one */
  Int_t nNumberOfSlots = 1;
  if (harmonicMap) {
    /* the highest harmonic number within the map if any */
    nNumberOfSlots += harmonicMap[nNoOfHarmonics - 1];
  } else {
    nNumberOfSlots += nNoOfHarmonics;
  }
  /* now allocate the slots for the values histograms */
  fXValues = new THnF *[nNumberOfSlots];
  fYValues = new THnF *[nNumberOfSlots];
  /* and initiallize them */
  for (Int_t i = 0; i < nNumberOfSlots; i++) {
    fXValues[i] = nullptr;
    fYValues[i] = nullptr;
  }
  /* now prepare the construction of the histograms */
  Int_t nVariables = fEventClassVariables.GetSize();
  auto *minvals = new Double_t[nVariables];
  auto *maxvals = new Double_t[nVariables];
  auto *nbins = new Int_t[nVariables];
  TString sVariableLabels = "";
  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins, minvals, maxvals);
  /* create the values multidimensional histograms for each harmonic */
  Int_t currentHarmonic = 0;
  for (Int_t i = 0; i < nNoOfHarmonics; i++) {
    if (harmonicMap) {
      currentHarmonic = harmonicMap[i];
    } else {
      currentHarmonic++;
    }
    fXValues[currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoXName, currentHarmonic),
                                         Form("%s h%d", (const char *) histoXTitle, currentHarmonic),
                                         nVariables, nbins, minvals, maxvals);
    fYValues[currentHarmonic] = new THnF(Form("%s_h%d", (const char *) histoYName, currentHarmonic),
                                         Form("%s h%d", (const char *) histoYTitle, currentHarmonic),
                                         nVariables, nbins, minvals, maxvals);
    /* now let's set the proper binning and label on each axis */
    for (Int_t var = 0; var < nVariables; var++) {
      fXValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                   fEventClassVariables[var].GetBins());
      fXValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
      fXValues[currentHarmonic]->GetAxis(var)->SetName(fEventClassVariables[var].GetLabel().data());
      fYValues[currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                   fEventClassVariables[var].GetBins());
      fYValues[currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
      fYValues[currentHarmonic]->GetAxis(var)->SetName(fEventClassVariables[var].GetLabel().data());
    }
    /* ask for square sum accumulation */
    fXValues[currentHarmonic]->Sumw2();
    fYValues[currentHarmonic]->Sumw2();
    /* and finally add the histograms to the list */
    histogramList->Add(fXValues[currentHarmonic]);
    histogramList->Add(fYValues[currentHarmonic]);
    /* and update the fully filled condition */
    fFullFilled |= harmonicNumberMask[currentHarmonic];
  }
  /* create the entries multidimensional histogram */
  fEntries =
      new THnI((const char *) entriesHistoName, (const char *) entriesHistoTitle, nVariables, nbins, minvals, maxvals);
  /* now let's set the proper binning and label on each entries histogram axis */
  for (Int_t var = 0; var < nVariables; var++) {
    fEntries->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(), fEventClassVariables[var].GetBins());
    fEntries->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
    fEntries->GetAxis(var)->SetName(fEventClassVariables[var].GetLabel().data());
  }
  /* and finally add the entries histogram to the list */
  histogramList->Add(fEntries);
  delete[] minvals;
  delete[] maxvals;
  delete[] nbins;
  return kTRUE;
}

/// Attaches existing histograms as the support histograms for X, Y, component
/// of the profile function for different harmonics
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// The harmonic map is inferred from the found histograms within the list
/// that match the naming scheme.
///
/// \param histogramList list where the histograms have to be located
/// \return true if properly attached else false
Bool_t CorrectionProfileComponents::AttachHistograms(TList *histogramList) {
  /* let's build the histograms names */
  TString histoXName = GetName();
  histoXName += szXComponentSuffix;
  TString histoYName = GetName();
  histoYName += szYComponentSuffix;
  TString entriesHistoName = GetName();
  entriesHistoName += szXComponentSuffix;
  entriesHistoName += szYComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  /* initialize. Remember we don't own the histograms */
  fEntries = nullptr;
  if (fXValues != nullptr) {
    delete[] fXValues;
    fXValues = nullptr;
  }
  if (fYValues != nullptr) {
    delete[] fYValues;
    fYValues = nullptr;
  }
  fXharmonicFillMask = 0x0000;
  fYharmonicFillMask = 0x0000;
  fFullFilled = 0x0000;
  fEntries = (THnI *) histogramList->FindObject((const char *) entriesHistoName);
  if (fEntries != nullptr) {
    /* allocate enough space for the supported harmonic numbers */
    fXValues = new THnF *[nMaxHarmonicNumberSupported + 1];
    fYValues = new THnF *[nMaxHarmonicNumberSupported + 1];
    /* search the multidimensional histograms for each harmonic */
    Int_t currentHarmonic = 0;
    for (Int_t i = 0; i < nMaxHarmonicNumberSupported; i++) {
      currentHarmonic++;
      fXValues[currentHarmonic] =
          (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoXName, currentHarmonic));
      fYValues[currentHarmonic] =
          (THnF *) histogramList->FindObject(Form("%s_h%d", (const char *) histoYName, currentHarmonic));
      /* and update the fully filled condition whether applicable */
      if ((fXValues[currentHarmonic] != nullptr) && (fYValues[currentHarmonic] != nullptr))
        fFullFilled |= harmonicNumberMask[currentHarmonic];
    }
  } else {
    return kFALSE;
  }
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
Long64_t CorrectionProfileComponents::GetBin() {
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
Bool_t CorrectionProfileComponents::BinContentValidated(Long64_t bin) {
  return Int_t(fEntries->GetBinContent(bin)) >= fMinNoOfEntriesToValidate;
}

/// Get the X component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileComponents::GetXBinContent(Int_t harmonic, Long64_t bin) {
  /* sanity check */
  if (!fXValues[harmonic]) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    return fXValues[harmonic]->GetBinContent(bin) / fEntries->GetBinContent(bin);
  }
}

/// Get the Y component bin content for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfileComponents::GetYBinContent(Int_t harmonic, Long64_t bin) {
  /* sanity check */
  if (!fYValues[harmonic]) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    return fYValues[harmonic]->GetBinContent(bin) / fEntries->GetBinContent(bin);
  }
}

/// Get the X component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfileComponents::GetXBinError(Int_t harmonic, Long64_t bin) {
  /* sanity check */
  if (!fXValues[harmonic]) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXValues[harmonic]->GetBinContent(bin);
    Float_t error2 = fXValues[harmonic]->GetBinError2(bin);
    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
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

/// Get the Y component bin content error for the passed bin number
/// for the corresponding harmonic
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfileComponents::GetYBinError(Int_t harmonic, Long64_t bin) {
  /* sanity check */
  if (fYValues[harmonic] == nullptr) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYValues[harmonic]->GetBinContent(bin);
    Float_t error2 = fYValues[harmonic]->GetBinError2(bin);
    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
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

/// Fills the X component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries is only updated if the whole set for both components
/// has been already filled. A check is done for detecting consecutive
/// fills for certain harmonic without a previous entries update.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void CorrectionProfileComponents::FillX(Int_t harmonic, Float_t weight) {
  /* first the sanity checks */
  if (fXValues[harmonic] == nullptr) {
    return;
  }
  if (fXharmonicFillMask & harmonicNumberMask[harmonic]) {
  }
  /* now it's safe to continue */
  /* keep total entries in fValues updated */
  Double_t nEntries = fXValues[harmonic]->GetEntries();
  FillBinAxesValues();
  fXValues[harmonic]->Fill(fBinAxesValues, weight);
  fXValues[harmonic]->SetEntries(nEntries + 1);
  /* update harmonic fill mask */
  fXharmonicFillMask |= harmonicNumberMask[harmonic];
  /* now check if time for updating entries histogram */
  if (fXharmonicFillMask != fFullFilled) return;
  if (fYharmonicFillMask != fFullFilled) return;
  /* update entries and reset the masks */
  fEntries->Fill(fBinAxesValues, 1.0);
  fXharmonicFillMask = 0x0000;
  fYharmonicFillMask = 0x0000;
}

/// Fills the Y component for the corresponding harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the given weight.
/// The entries is only updated if the whole set for both components
/// has been already filled. A check is done for detecting consecutive
/// fills for certain harmonic without a previous entries update.
///
/// \param harmonic the interested external harmonic number
/// \param variableContainer the current variables content addressed by var Id
/// \param weight the increment in the bin content
void CorrectionProfileComponents::FillY(Int_t harmonic, Float_t weight) {
  /* first the sanity checks */
  if (fYValues[harmonic] == nullptr) {
    return;
  }
  if (fYharmonicFillMask & harmonicNumberMask[harmonic]) {
    return;
  }
  /* now it's safe to continue */
  /* keep total entries in fValues updated */
  Double_t nEntries = fYValues[harmonic]->GetEntries();
  FillBinAxesValues();
  fYValues[harmonic]->Fill(fBinAxesValues, weight);
  fYValues[harmonic]->SetEntries(nEntries + 1);
  /* update harmonic fill mask */
  fYharmonicFillMask |= harmonicNumberMask[harmonic];
  /* now check if time for updating entries histogram */
  if (fYharmonicFillMask != fFullFilled) return;
  if (fXharmonicFillMask != fFullFilled) return;
  /* update entries and reset the masks */
  fEntries->Fill(fBinAxesValues, 1.0);
  fXharmonicFillMask = 0x0000;
  fYharmonicFillMask = 0x0000;
}
}// namespace Qn
