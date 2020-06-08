/// \file QnCorrectionsProfile3DCorrelations.cxx
/// \brief Implementation of the multidimensional correlation components based set of profiles
/// for three Qn vectors with harmonic support

#include "TList.h"

#include "CorrectionAxisSet.hpp"
#include "QVector.hpp"

#include "CorrectionProfile3DCorrelations.hpp"

/// \cond CLASSIMP
ClassImp(Qn::CorrectionProfile3DCorrelations);
/// \endcond
namespace Qn {
///< the number of Qn supported
#define CORRELATIONSNOOFQNVECTORS 3
/// Normal constructor
///
/// Stores the set of variables that identify the
/// different event classes passing them to its parent
/// and prepares the object for actual histogram
/// creation or attachment
///
/// \param name base for the name of the histograms
/// \param title base for the title of the histograms
/// \param nameA A detector name
/// \param nameB B detector name
/// \param nameC C detector name
/// \param ecvs the event classes variables set
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
CorrectionProfile3DCorrelations::CorrectionProfile3DCorrelations(std::string name,
                                                                 std::string title,
                                                                 std::string nameA,
                                                                 std::string nameB,
                                                                 std::string nameC,
                                                                 const CorrectionAxisSet &ecvs,
                                                                 ErrorMode mode) : CorrectionHistogramBase(name, title, ecvs, mode), fNameA(nameA), fNameB(nameB), fNameC(nameC) {
}

CorrectionProfile3DCorrelations::CorrectionProfile3DCorrelations(std::string name,
                                                                 std::string nameA,
                                                                 std::string nameB,
                                                                 std::string nameC,
                                                                 const CorrectionAxisSet &ecvs,
                                                                 ErrorMode mode) : CorrectionHistogramBase(name, name, ecvs, mode), fNameA(nameA), fNameB(nameB), fNameC(nameC) {
}

/// Default destructor
///
/// Returns the only taken memory, the harmonic histograms storage,
/// the own histograms and other members are not own at destruction time
CorrectionProfile3DCorrelations::~CorrectionProfile3DCorrelations() {
  if (fXXValues) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      delete[] fXXValues[ixComb];
    }
    delete[] fXXValues;
  }
  if (fXYValues) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      delete[] fXYValues[ixComb];
    }
    delete[] fXYValues;
  }
  if (fYXValues) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      delete[] fYXValues[ixComb];
    }
    delete[] fYXValues;
  }
  if (fYYValues) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      delete[] fYYValues[ixComb];
    }
    delete[] fYYValues;
  }
}

/// Creates the XX, XY, YX, YY correlation components support histograms
/// for the profile function for each combination of Qn vectors
///
/// Based in the event classes variables set in the parent class
/// the values and entries multidimensional histograms are
/// created.
///
/// For each Qn vector comibination, three,
/// and for each harmonic number four values histograms are created, XX,
/// XY, YX and YY. The histograms are organized to support external harmonic
/// number. By default the external harmonic number is always considered to
/// start by one. If no map is passed as parameter the external harmonic
/// numbers are considered as: 1, 2, ..., nNoOfHarmonic.
/// If the user wants a different assignment he has to provide an
/// ordered map, for instance: four harmonics with external harmonic numbers
/// 2, 4, 6 and 8 will require nNoOfHarmonics = 4 and harmonicMap = [2, 4, 6, 8].
/// The fully filled condition is computed and stored
///
/// The potential situation where the Qn vector has an harmonic multiplier
/// is properly supported
///
/// The whole set of histograms are added to the passed histogram list
///
/// \param histogramList list where the histograms have to be added
/// \param nNoOfHarmonics the desired number of harmonics
/// \param nHarmonicMultiplier the multiplier for the harmonic number
/// \param harmonicMap ordered array with the external number of the harmonics
/// \return true if properly created
Bool_t CorrectionProfile3DCorrelations::CreateCorrelationComponentsProfileHistograms(TList *histogramList,
                                                                                     Int_t nNoOfHarmonics,
                                                                                     Int_t nHarmonicMultiplier,
                                                                                     Int_t *harmonicMap) {
  /* store the potential harmonic multiplier */
  fHarmonicMultiplier = nHarmonicMultiplier;
  /* for now on, everything is handled as if the multiplier were m=1 so, we only consider n */

  /* check whether within the supported harmonic range */
  Int_t nHigherHarmonic = nNoOfHarmonics;
  if (harmonicMap != nullptr) {
    nHigherHarmonic = harmonicMap[nNoOfHarmonics - 1];
  }
  if (nMaxHarmonicNumberSupported < nHigherHarmonic) {
    return false;
  }

  /* let's support the external harmonic number map */
  /* external harmonic number will always start from one */
  Int_t nNumberOfSlots = 1;
  if (harmonicMap != nullptr) {
    /* the highest harmonic number within the map if any */
    nNumberOfSlots += harmonicMap[nNoOfHarmonics - 1];
  } else {
    nNumberOfSlots += nNoOfHarmonics;
  }

  /* now allocate the slots for the values histograms for each Qn vector correlation combination */
  fXXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  fXYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  fYXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  fYYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
  for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
    fXXValues[ixComb] = new THnF *[nNumberOfSlots];
    fXYValues[ixComb] = new THnF *[nNumberOfSlots];
    fYXValues[ixComb] = new THnF *[nNumberOfSlots];
    fYYValues[ixComb] = new THnF *[nNumberOfSlots];
    /* and initiallize them */
    for (Int_t i = 0; i < nNumberOfSlots; i++) {
      fXXValues[ixComb][i] = nullptr;
      fXYValues[ixComb][i] = nullptr;
      fYXValues[ixComb][i] = nullptr;
      fYYValues[ixComb][i] = nullptr;
    }
  }
  /* now prepare the construction of the histograms */
  Int_t nVariables = fEventClassVariables.GetSize();
  Double_t *minvals = new Double_t[nVariables];
  Double_t *maxvals = new Double_t[nVariables];
  Int_t *nbins = new Int_t[nVariables];
  TString sVariableLabels = "";
  /* get the multidimensional structure */
  fEventClassVariables.GetMultidimensionalConfiguration(nbins, minvals, maxvals);
  /* create the values multidimensional histograms for each Qn vector correlation combination and harmonic */
  const char *combNames[CORRELATIONSNOOFQNVECTORS] = {fNameA.data(), fNameB.data(), fNameC.data()};
  for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
    Int_t currentHarmonic = 0;
    for (Int_t i = 0; i < nNoOfHarmonics; i++) {
      if (harmonicMap != nullptr) {
        currentHarmonic = harmonicMap[i];
      } else {
        currentHarmonic++;
      }
      /* let's build the histograms names and titles */
      TString
          BaseName = GetName() + " " + combNames[ixComb] + "x" + combNames[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS];
      TString
          BaseTitle = GetTitle() + " " + combNames[ixComb] + "x" + combNames[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS];
      TString histoXXName = BaseName;
      histoXXName += szXXCorrelationComponentSuffix;
      TString histoXYName = BaseName;
      histoXYName += szXYCorrelationComponentSuffix;
      TString histoYXName = BaseName;
      histoYXName += szYXCorrelationComponentSuffix;
      TString histoYYName = BaseName;
      histoYYName += szYYCorrelationComponentSuffix;
      TString histoXXTitle = BaseTitle;
      histoXXTitle += szXXCorrelationComponentSuffix;
      TString histoXYTitle = BaseTitle;
      histoXYTitle += szXYCorrelationComponentSuffix;
      TString histoYXTitle = BaseTitle;
      histoYXTitle += szYXCorrelationComponentSuffix;
      TString histoYYTitle = BaseTitle;
      histoYYTitle += szYYCorrelationComponentSuffix;
      fXXValues[ixComb][currentHarmonic] =
          new THnF(Form("%s_h%d", (const char *) histoXXName, currentHarmonic * fHarmonicMultiplier),
                   Form("%s h%d", (const char *) histoXXTitle, currentHarmonic),
                   nVariables, nbins, minvals, maxvals);
      fXYValues[ixComb][currentHarmonic] =
          new THnF(Form("%s_h%d", (const char *) histoXYName, currentHarmonic * fHarmonicMultiplier),
                   Form("%s h%d", (const char *) histoXYTitle, currentHarmonic),
                   nVariables, nbins, minvals, maxvals);
      fYXValues[ixComb][currentHarmonic] =
          new THnF(Form("%s_h%d", (const char *) histoYXName, currentHarmonic * fHarmonicMultiplier),
                   Form("%s h%d", (const char *) histoYXTitle, currentHarmonic),
                   nVariables, nbins, minvals, maxvals);
      fYYValues[ixComb][currentHarmonic] =
          new THnF(Form("%s_h%d", (const char *) histoYYName, currentHarmonic * fHarmonicMultiplier),
                   Form("%s h%d", (const char *) histoYYTitle, currentHarmonic),
                   nVariables, nbins, minvals, maxvals);
      /* now let's set the proper binning and label on each axis */
      for (Int_t var = 0; var < nVariables; var++) {
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fXXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fXYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fYXValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->Set(fEventClassVariables[var].GetNBins(),
                                                              fEventClassVariables[var].GetBins());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
        fYYValues[ixComb][currentHarmonic]->GetAxis(var)->SetTitle(fEventClassVariables[var].GetLabel().data());
      }
      /* ask for square sum accumulation */
      fXXValues[ixComb][currentHarmonic]->Sumw2();
      fXYValues[ixComb][currentHarmonic]->Sumw2();
      fYXValues[ixComb][currentHarmonic]->Sumw2();
      fYYValues[ixComb][currentHarmonic]->Sumw2();
      /* and finally add the histograms to the list */
      histogramList->Add(fXXValues[ixComb][currentHarmonic]);
      histogramList->Add(fXYValues[ixComb][currentHarmonic]);
      histogramList->Add(fYXValues[ixComb][currentHarmonic]);
      histogramList->Add(fYYValues[ixComb][currentHarmonic]);
    }
  }
  /* now the entries histogram name and title */
  TString entriesHistoName = GetName();
  entriesHistoName = entriesHistoName + fNameA + fNameB + fNameC;
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  TString entriesHistoTitle = GetTitle();
  entriesHistoTitle = entriesHistoTitle + fNameA + fNameB + fNameC;
  entriesHistoTitle += szXXCorrelationComponentSuffix;
  entriesHistoTitle += szXYCorrelationComponentSuffix;
  entriesHistoTitle += szYXCorrelationComponentSuffix;
  entriesHistoTitle += szYYCorrelationComponentSuffix;
  entriesHistoTitle += szEntriesHistoSuffix;
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
/// and for the different Qn vector correlation combinations.
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// The harmonic map is inferred from the found histograms within the list
/// that match the naming scheme.
///
/// The potential situation where the Qn vectors have an harmonic multiplier
/// is properly supported
///
/// \param histogramList list where the histograms have to be located
/// \return true if properly attached else false
Bool_t CorrectionProfile3DCorrelations::AttachHistograms(TList *histogramList) {
  /* initialize. Remember we don't own the histograms */
  fEntries = nullptr;
  if (fXXValues != nullptr) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fXXValues[ixComb] != nullptr)
        delete[] fXXValues[ixComb];
    }
    delete[] fXXValues;
  }
  if (fXYValues != nullptr) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fXYValues[ixComb] != nullptr)
        delete[] fXYValues[ixComb];
    }
    delete[] fXYValues;
  }
  if (fYXValues != nullptr) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fYXValues[ixComb] != nullptr)
        delete[] fYXValues[ixComb];
    }
    delete[] fYXValues;
  }
  if (fYYValues != nullptr) {
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      if (fYYValues[ixComb] != nullptr)
        delete[] fYYValues[ixComb];
    }
    delete[] fYYValues;
  }
  /* let's build the entries histogram name */
  TString entriesHistoName = GetName();
  entriesHistoName = entriesHistoName + fNameA + fNameB + fNameC;
  entriesHistoName += szXXCorrelationComponentSuffix;
  entriesHistoName += szXYCorrelationComponentSuffix;
  entriesHistoName += szYXCorrelationComponentSuffix;
  entriesHistoName += szYYCorrelationComponentSuffix;
  entriesHistoName += szEntriesHistoSuffix;
  UInt_t harmonicFilledMask = 0x0000;
  fEntries = (THnI *) histogramList->FindObject((const char *) entriesHistoName);
  if (fEntries && fEntries->GetEntries() != 0) {
    /* allocate enough space for the supported harmonic numbers */
    fXXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    fXYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    fYXValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    fYYValues = new THnF **[CORRELATIONSNOOFQNVECTORS];
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      fXXValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
      fXYValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
      fYXValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
      fYYValues[ixComb] = new THnF *[nMaxHarmonicNumberSupported + 1];
    }
    /* search the multidimensional histograms for each harmonic and Qn vecto correlation combination */
    const char *combNames[CORRELATIONSNOOFQNVECTORS] = {fNameA.data(), fNameB.data(), fNameC.data()};
    for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
      Int_t currentHarmonic = 0;
      /* let's build the histograms names */
      TString
          BaseName = GetName() + " " + combNames[ixComb] + "x" + combNames[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS];
      TString
          BaseTitle = GetTitle() + " " + combNames[ixComb] + "x" + combNames[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS];
      TString histoXXName = BaseName;
      histoXXName += szXXCorrelationComponentSuffix;
      TString histoXYName = BaseName;
      histoXYName += szXYCorrelationComponentSuffix;
      TString histoYXName = BaseName;
      histoYXName += szYXCorrelationComponentSuffix;
      TString histoYYName = BaseName;
      histoYYName += szYYCorrelationComponentSuffix;
      for (Int_t i = 0; i < nMaxHarmonicNumberSupported; i++) {
        currentHarmonic++;
        fXXValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d",
                                                                                     (const char *) histoXXName,
                                                                                     currentHarmonic
                                                                                         * fHarmonicMultiplier));
        fXYValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d",
                                                                                     (const char *) histoXYName,
                                                                                     currentHarmonic
                                                                                         * fHarmonicMultiplier));
        fYXValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d",
                                                                                     (const char *) histoYXName,
                                                                                     currentHarmonic
                                                                                         * fHarmonicMultiplier));
        fYYValues[ixComb][currentHarmonic] = (THnF *) histogramList->FindObject(Form("%s_h%d",
                                                                                     (const char *) histoYYName,
                                                                                     currentHarmonic
                                                                                         * fHarmonicMultiplier));

        /* update the correcto condition */
        if ((fXXValues[ixComb][currentHarmonic] != nullptr) && (fXYValues[ixComb][currentHarmonic] != nullptr)
            && (fYXValues[ixComb][currentHarmonic] != nullptr) && (fYYValues[ixComb][currentHarmonic] != nullptr))
          harmonicFilledMask |= harmonicNumberMask[currentHarmonic];
      }
    }
  } else {
    return kFALSE;
  }

  /* check that we actually got something */
  return harmonicFilledMask != 0x0000;
}

/// Get the bin number for the current variable content
///
/// The bin number identifies the event class the current
/// variable content points to.
///
/// \param variableContainer the current variables content addressed by var Id
/// \return the associated bin to the current variables content
Long64_t CorrectionProfile3DCorrelations::GetBin() {
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
Bool_t CorrectionProfile3DCorrelations::BinContentValidated(Long64_t bin) {
  auto nEntries = Int_t(fEntries->GetBinContent(bin));

  return nEntries >= fMinNoOfEntriesToValidate;
}

/// Get the XX correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfile3DCorrelations::GetXXBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (fXXValues[ixComb][harmonic] == nullptr) {
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXXValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XY correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfile3DCorrelations::GetXYBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (!fXYValues[ixComb][harmonic]) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fXYValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YX correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfile3DCorrelations::GetYXBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (!fYXValues[ixComb][harmonic]) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYXValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the YY correlation component bin content for the passed bin number
/// for the corresponding harmonic and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin number content
Float_t CorrectionProfile3DCorrelations::GetYYBinContent(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (fYYValues[ixComb][harmonic] == nullptr) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    return fYYValues[ixComb][harmonic]->GetBinContent(bin) / Float_t(nEntries);
  }
}

/// Get the XX correlation component bin content error for the passed bin number
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfile3DCorrelations::GetXXBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (fXXValues[ixComb][harmonic] == nullptr) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXXValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fXXValues[ixComb][harmonic]->GetBinError2(bin);

    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN: return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD: return serror;
      default: return 0.0;
    }
  }
}

/// Get the XY correlation component bin content error for the passed bin number
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfile3DCorrelations::GetXYBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (fXYValues[ixComb][harmonic] == nullptr) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    Int_t nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fXYValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fXYValues[ixComb][harmonic]->GetBinError2(bin);

    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN: return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD: return serror;
      default: return 0.0;
    }
  }
}

/// Get the YX correlation component bin content error for the passed bin number
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfile3DCorrelations::GetYXBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;

  /* sanity check */
  if (fYXValues[ixComb][harmonic] == nullptr) {
    return 0.0;
  }

  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYXValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fYXValues[ixComb][harmonic]->GetBinError2(bin);

    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN: return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD: return serror;
      default: return 0.0;
    }
  }
}

/// Get the YY correlation component bin content error for the passed bin number
/// for the corresponding harmonic  and Qn vector combination
///
/// The bin number identifies a desired event class whose content is
/// error is requested. If the bin content is not validated zero is returned.
///
/// \param comb the name of the desired Qn vector combination: "AB", "BC" or "AC"
/// \param harmonic the interested external harmonic number
/// \param bin the interested bin number
/// \return the bin content error
Float_t CorrectionProfile3DCorrelations::GetYYBinError(const char *comb, Int_t harmonic, Long64_t bin) {
  TString szComb = comb;
  Int_t ixComb = -1;
  if (szComb.EqualTo("AB"))
    ixComb = 0;
  else if (szComb.EqualTo("BC"))
    ixComb = 1;
  else if (szComb.EqualTo("AC"))
    ixComb = 2;
  else
    return 0.0;
  /* sanity check */
  if (fYYValues[ixComb][harmonic] == nullptr) {
    return 0.0;
  }
  if (!BinContentValidated(bin)) {
    return 0.0;
  } else {
    auto nEntries = Int_t(fEntries->GetBinContent(bin));
    Float_t values = fYYValues[ixComb][harmonic]->GetBinContent(bin);
    Float_t error2 = fYYValues[ixComb][harmonic]->GetBinError2(bin);
    Double_t average = values / nEntries;
    Double_t serror = TMath::Sqrt(TMath::Abs(error2 / nEntries - average * average));
    switch (fErrorMode) {
      case ErrorMode::MEAN: return serror / TMath::Sqrt(nEntries);
      case ErrorMode::SPREAD: return serror;
      default: return 0.0;
    }
  }
}

/// Fills the correlation component for the different Qn vector correlation combinations
/// and for all handled harmonic histogram
///
/// The involved bin is computed according to the current variables
/// content. The bin is then increased by the corresponding values.
/// The entries count is updated accordingly.
///
/// It is considered that the three Qn vectors have the same harmonic
/// structure including the harmonic multiplier. If this is not the case
/// and that situation should be supported this member must be modified.
/// \param QnA A Qn vector
/// \param QnB B Qn vector
/// \param QnC C Qn vector
/// \param variableContainer the current variables content addressed by var Id
void CorrectionProfile3DCorrelations::Fill(const QVector *QnA,
                                           const QVector *QnB,
                                           const QVector *QnC) {
  /* first the sanity checks */
  if (!((QnA->IsGoodQuality()) && (QnB->IsGoodQuality()) && (QnC->IsGoodQuality()))) return;
  if ((QnA->GetHarmonicMultiplier() != QnB->GetHarmonicMultiplier())
      || (QnA->GetHarmonicMultiplier() != QnC->GetHarmonicMultiplier())) {
    return;
  }
  /* let's get the axis information */
  FillBinAxesValues();
  /* consider all combinations */
  const QVector *combQn[CORRELATIONSNOOFQNVECTORS] = {QnA, QnB, QnC};
  for (Int_t ixComb = 0; ixComb < CORRELATIONSNOOFQNVECTORS; ixComb++) {
    /* and all harmonics */
    Int_t nCurrentHarmonic = QnA->GetFirstHarmonic();
    while (nCurrentHarmonic != -1) {
      /* first the sanity checks */
      if (fXXValues[ixComb][nCurrentHarmonic] == nullptr) {
        return;
      }
      /* keep total entries in fValues updated */
      Double_t nXXEntries = fXXValues[ixComb][nCurrentHarmonic]->GetEntries();
      Double_t nXYEntries = fXYValues[ixComb][nCurrentHarmonic]->GetEntries();
      Double_t nYXEntries = fYXValues[ixComb][nCurrentHarmonic]->GetEntries();
      Double_t nYYEntries = fYYValues[ixComb][nCurrentHarmonic]->GetEntries();
      fXXValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues,
                                                combQn[ixComb]->x(nCurrentHarmonic) * combQn[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS]->x(nCurrentHarmonic));
      fXYValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues,
                                                combQn[ixComb]->x(nCurrentHarmonic) * combQn[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS]->y(nCurrentHarmonic));
      fYXValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues,
                                                combQn[ixComb]->y(nCurrentHarmonic) * combQn[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS]->x(nCurrentHarmonic));
      fYYValues[ixComb][nCurrentHarmonic]->Fill(fBinAxesValues,
                                                combQn[ixComb]->y(nCurrentHarmonic) * combQn[(ixComb + 1) % CORRELATIONSNOOFQNVECTORS]->y(nCurrentHarmonic));
      fXXValues[ixComb][nCurrentHarmonic]->SetEntries(nXXEntries + 1);
      fXYValues[ixComb][nCurrentHarmonic]->SetEntries(nXYEntries + 1);
      fYXValues[ixComb][nCurrentHarmonic]->SetEntries(nYXEntries + 1);
      fYYValues[ixComb][nCurrentHarmonic]->SetEntries(nYYEntries + 1);
      nCurrentHarmonic = QnA->GetNextHarmonic(nCurrentHarmonic);
    }
  }
  /* update the profile entries */
  fEntries->Fill(fBinAxesValues, 1.0);
}
}// namespace Qn