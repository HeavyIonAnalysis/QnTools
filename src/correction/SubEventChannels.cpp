#include <utility>

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

/// \file QnCorrectionsDetectorConfigurationChannels.cxx
/// \brief Implementation of the channel detector configuration class

#include "CorrectionProfileComponents.hpp"
#include "SubEventChannels.hpp"

#include "Detector.hpp"

/// \cond CLASSIMP
ClassImp(Qn::SubEventChannels);
/// \endcond
namespace Qn {
const char *SubEventChannels::szRawQnVectorName = "raw";
const char *SubEventChannels::szQAMultiplicityHistoName = "Multiplicity";

/// Normal constructor
/// Allocates the data vector bank.
/// \param name the name of the detector configuration
/// \param eventClassesVariables the set of event classes variables
/// \param nNoOfChannels the number of channels of the associated detector
/// \param nNoOfHarmonics the number of harmonics that must be handled
/// \param harmonicMap an optional ordered array with the harmonic numbers
SubEventChannels::SubEventChannels(unsigned int bin_id,
                                   const CorrectionAxisSet *eventClassesVariables,
                                   Int_t nNoOfChannels, std::bitset<QVector::kmaxharmonics> harmonics) : SubEvent(bin_id, eventClassesVariables, harmonics),
                                                                                                         fRawQnVector(harmonics, QVector::CorrectionStep::RAW),
                                                                                                         fInputDataCorrections() {
  fNoOfChannels = nNoOfChannels;
}

/// Default destructor
/// Releases the memory taken
SubEventChannels::~SubEventChannels() {
  delete[] fUsedChannel;
  delete[] fChannelMap;
  delete[] fChannelGroup;
  delete[] fHardCodedGroupWeights;
}

void SubEventChannels::SetChannelsScheme(std::vector<int> channel_groups) {
  /* TODO: there should be smart procedures on how to improve the channels scan for actual data */
  fUsedChannel = new Bool_t[fNoOfChannels];
  fChannelMap = new Int_t[fNoOfChannels];
  fChannelGroup = new Int_t[fNoOfChannels];
  Int_t nMinGroup = 0xFFFF;
  Int_t nMaxGroup = 0x0000;
  Int_t intChannelNo = 0;
  for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
    fUsedChannel[ixChannel] = kTRUE;
    if (fUsedChannel[ixChannel]) {
      fChannelMap[ixChannel] = intChannelNo;
      intChannelNo++;
      fChannelGroup[ixChannel] = channel_groups[ixChannel];
      /* update min max group number */
      if (channel_groups[ixChannel] < nMinGroup)
        nMinGroup = channel_groups[ixChannel];
      if (nMaxGroup < channel_groups[ixChannel])
        nMaxGroup = channel_groups[ixChannel];
    }
  }
}

/// Asks for support data structures creation
///
/// The input data vector bank is allocated and the request is
/// transmitted to the input data corrections and then to the Q vector corrections.
void SubEventChannels::CreateSupportQVectors() {
  /* this is executed in the remote node so, allocate the data bank */
  fDataVectorBank.reserve(Qn::SubEvent::INITIALSIZE);
  for (auto &correction : fInputDataCorrections) {
    correction->CreateSupportQVectors();
  }
  for (auto &correction : fQnVectorCorrections) {
    correction->CreateSupportQVectors();
  }
}

/// Asks for support histograms creation
///
/// A new histograms list is created for the detector and incorporated
/// to the passed list. Then the new list is passed first to the input data corrections
/// and then to the Q vector corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void SubEventChannels::CreateCorrectionHistograms() {
  fInputDataCorrections.CreateCorrectionHistograms();
  fQnVectorCorrections.CreateCorrectionHistograms();
  if (!fInputDataCorrections.Empty()) {
    fInputDataCorrections.EnableFirstCorrection();
  } else {
    fQnVectorCorrections.EnableFirstCorrection();
  }
  //  /* if list is empty delete it if not incorporate it */
  //  if (!correction_list->IsEmpty()) {
  //    list->Add(correction_list);
  //  } else {
  //    delete correction_list;
  //  }
}

void SubEventChannels::CopyToOutputList(TList *list) {
  auto correction_list = new TList();
  correction_list->SetName(GetName().data());
  correction_list->SetOwner(kTRUE);
  fInputDataCorrections.CopyToOutputList(correction_list);
  fQnVectorCorrections.CopyToOutputList(correction_list);
  /* if list is empty delete it if not incorporate it */
  if (!correction_list->IsEmpty()) {
    list->Add(correction_list);
  } else {
    delete correction_list;
  }
}

/// Asks for QA histograms creation
///
/// A new histograms list is created for the detector and incorporated
/// to the passed list. The own QA histograms are then created and incorporated
/// to the new list. Then the new list is passed first to the input data corrections
/// and then to the Q vector corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void SubEventChannels::AttachQAHistograms(TList *list) {
  TList *detectorConfigurationList;
  if (GetName().empty()) {
    detectorConfigurationList = list;
  } else {
    detectorConfigurationList = new TList();
    detectorConfigurationList->SetName(GetName().data());
    detectorConfigurationList->SetOwner(kTRUE);
  }
  /* first create our own QA histograms */
  std::string beforeName{GetName() + szQAMultiplicityHistoName + "Before"};
  std::string beforeTitle{GetName() + " " + szQAMultiplicityHistoName + " before input equalization"};
  std::string afterName{GetName() + szQAMultiplicityHistoName + "After"};
  std::string afterTitle{GetName() + " " + szQAMultiplicityHistoName + " after input equalization"};
  /* let's pick the centrality variable and its binning */
  Int_t ixVarId = -1;
  int ivar = 0;
  for (const auto &var : *fEventClassVariables) {
    if (var.GetId() != fQACentralityVarId) {
      ++ivar;
      continue;
    } else {
      ixVarId = ivar;
      break;
    }
  }
  /* let's get the effective number of channels */
  Int_t nNoOfChannels = 0;
  for (Int_t i = 0; i < fNoOfChannels; i++)
    if (fUsedChannel[i])
      nNoOfChannels++;
  if (ixVarId != -1) {
    fQAMultiplicityBefore3D = new TH3F(
        beforeName.data(),
        beforeTitle.data(),
        fEventClassVariables->At(ixVarId).GetNBins(),
        fEventClassVariables->At(ixVarId).GetLowerEdge(),
        fEventClassVariables->At(ixVarId).GetUpperEdge(),
        nNoOfChannels,
        0.0,
        nNoOfChannels,
        fQAnBinsMultiplicity,
        fQAMultiplicityMin,
        fQAMultiplicityMax);
    fQAMultiplicityAfter3D = new TH3F(
        afterName.data(),
        afterTitle.data(),
        fEventClassVariables->At(ixVarId).GetNBins(),
        fEventClassVariables->At(ixVarId).GetLowerEdge(),
        fEventClassVariables->At(ixVarId).GetUpperEdge(),
        nNoOfChannels,
        0.0,
        nNoOfChannels,
        fQAnBinsMultiplicity,
        fQAMultiplicityMin,
        fQAMultiplicityMax);
    /* now set the proper labels and titles */
    fQAMultiplicityBefore3D->GetXaxis()->SetTitle(fEventClassVariables->At(ixVarId).GetLabel().data());
    fQAMultiplicityBefore3D->GetYaxis()->SetTitle("channel");
    fQAMultiplicityBefore3D->GetZaxis()->SetTitle("M");
    fQAMultiplicityAfter3D->GetXaxis()->SetTitle(fEventClassVariables->At(ixVarId).GetLabel().data());
    fQAMultiplicityAfter3D->GetYaxis()->SetTitle("channel");
    fQAMultiplicityAfter3D->GetZaxis()->SetTitle("M");
    if (fNoOfChannels != nNoOfChannels) {
      Int_t bin = 1;
      for (Int_t i = 0; i < fNoOfChannels; i++)
        if (fUsedChannel[i]) {
          fQAMultiplicityBefore3D->GetYaxis()->SetBinLabel(bin, Form("%d", i));
          fQAMultiplicityAfter3D->GetYaxis()->SetBinLabel(bin, Form("%d", i));
          bin++;
        }
    }
    detectorConfigurationList->Add(fQAMultiplicityBefore3D);
    detectorConfigurationList->Add(fQAMultiplicityAfter3D);
  }
  /* now propagate it to the input data corrections */
  for (auto &correction : fInputDataCorrections) {
    correction->AttachQAHistograms(detectorConfigurationList);
  }
  /* the own QA average Qn vector components histogram */
  const auto qaname = std::string(szQAQnAverageHistogramName) + " " + GetName();
  fQAQnAverageHistogram = std::make_unique<CorrectionProfileComponents>(qaname, qaname, GetEventClassVariablesSet());
  /* get information about the configured harmonics to pass it for histogram creation */
  Int_t nNoOfHarmonics = this->GetNoOfHarmonics();
  auto harmonicsMap = new Int_t[nNoOfHarmonics];
  this->GetHarmonicMap(harmonicsMap);
  fQAQnAverageHistogram->CreateComponentsProfileHistograms(detectorConfigurationList, nNoOfHarmonics, harmonicsMap);
  delete[] harmonicsMap;
  /* if everything right propagate it to Q vector corrections */
  for (auto &correction : fQnVectorCorrections) {
    correction->AttachQAHistograms(detectorConfigurationList);
  }
  /* now incorporate the list to the passed one */
  if (!detectorConfigurationList->IsEmpty()) {
    list->Add(detectorConfigurationList);
  } else {
    delete detectorConfigurationList;
  }
}

/// Asks for non validated entries QA histograms creation
///
/// A new histograms list is created for the detector and incorporated
/// to the passed list. Then the new list is passed first to the input data corrections
/// and then to the Q vector corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
void SubEventChannels::AttachNveQAHistograms(TList *list) {
  auto detectorConfigurationList = new TList();
  detectorConfigurationList->SetName(GetName().data());
  detectorConfigurationList->SetOwner(kTRUE);
  for (auto &correction : fInputDataCorrections) {
    correction->AttachNveQAHistograms(detectorConfigurationList);
  }
  for (auto &correction : fQnVectorCorrections) {
    correction->AttachNveQAHistograms(detectorConfigurationList);
  }
  /* now incorporate the list to the passed one */
  if (!detectorConfigurationList->IsEmpty()) {
    list->Add(detectorConfigurationList);
  } else {
    delete detectorConfigurationList;
  }
}

/// Asks for attaching the needed input information to the correction steps
///
/// The detector list is extracted from the passed list and then
/// the request is transmitted to the input data corrections
/// and then propagated to the Q vector corrections
/// \param list list where the input information should be found
/// \return kTRUE if everything went OK
void SubEventChannels::AttachCorrectionInput(TList *list) {
  auto input_list = (TList *) list->FindObject(GetName().data());
  if (input_list) {
    if (!fInputDataCorrections.Empty()) {
      fInputDataCorrections.EnableFirstCorrection();
      fInputDataCorrections.AttachInputs(input_list);
    }
    if (fInputDataCorrections.Empty() || fInputDataCorrections.IsLastStepApplied()) {
      fQnVectorCorrections.EnableFirstCorrection();
      fQnVectorCorrections.AttachInputs(input_list);
    }
  }
}

/// Perform after calibration histograms attach actions
/// It is used to inform the different correction step that
/// all conditions for running the network are in place so
/// it is time to check if their requirements are satisfied
///
/// The request is transmitted to the input data corrections
/// and then propagated to the Q vector corrections
void SubEventChannels::AfterInputAttachAction() {
  for (auto &correction : fInputDataCorrections) {
    correction->AfterInputAttachAction();
  }
  for (auto &correction : fQnVectorCorrections) {
    correction->AfterInputAttachAction();
  }
}

/// Incorporates the passed correction to the set of input data corrections
/// \param correctionOnInputData the correction to add
void SubEventChannels::AddCorrectionOnInputData(CorrectionOnInputData *correctionOnInputData) {
  correctionOnInputData->SetOwner(this);
  fInputDataCorrections.AddCorrection(correctionOnInputData);
}

/// Fills the QA multiplicity histograms before and after input equalization
/// and the plain Qn vector average components histogram
/// \param variableContainer pointer to the variable content bank
void SubEventChannels::FillQAHistograms() {
  if (fQAMultiplicityBefore3D && fQAMultiplicityAfter3D) {
    for (const auto &dataVector : fDataVectorBank) {
      fQAMultiplicityBefore3D->Fill(fEventClassVariables->At(fQACentralityVarId).GetValue(),
                                    fChannelMap[dataVector.GetId()],
                                    dataVector.Weight());
      fQAMultiplicityAfter3D->Fill(fEventClassVariables->At(fQACentralityVarId).GetValue(),
                                   fChannelMap[dataVector.GetId()],
                                   dataVector.EqualizedWeight());
    }
  }
  if (fQAQnAverageHistogram) {
    Int_t harmonic = fPlainQnVector.GetFirstHarmonic();
    while (harmonic != -1) {
      fQAQnAverageHistogram->FillX(harmonic, fPlainQnVector.x(harmonic));
      fQAQnAverageHistogram->FillY(harmonic, fPlainQnVector.y(harmonic));
      harmonic = fPlainQnVector.GetNextHarmonic(harmonic);
    }
  }
}

/// Include the the list of Qn vector associated to the detector configuration
/// into the passed list
///
/// A new list is created for the detector configuration and incorporated
/// to the passed list.
///
/// Always includes first the fully corrected Qn vector,
/// and then includes the raw Qn vector and the plain Qn vector and then
/// asks to the different correction
/// steps to include their partially corrected Qn vectors.
///
/// The check if we are already there is because it could be late information
/// about the process name and then the correction histograms could still not
/// be attached and the constructed list does not contain the final Qn vectors.
/// \param list list where the corrected Qn vector should be added
inline void SubEventChannels::IncludeQnVectors() {
  qvectors_.emplace(QVector::CorrectionStep::RAW, &fRawQnVector);
  qvectors_.emplace(QVector::CorrectionStep::PLAIN, &fPlainQnVector);
  for (auto &correction : fQnVectorCorrections) {
    correction->IncludeCorrectedQnVector(qvectors_);
  }
}

/// Include only one instance of each input correction step
/// in execution order
///
/// The request is transmitted to the set of Qn vector corrections
/// \param list list where the correction steps should be incorporated
void SubEventChannels::FillOverallInputCorrectionStepList(std::set<CorrectionBase *> &set) const {
  fInputDataCorrections.FillOverallCorrectionsList(set);
}

/// Include only one instance of each Qn vector correction step
/// in execution order
///
/// The request is transmitted to the set of Qn vector corrections
/// \param list list where the correction steps should be incorporated
void SubEventChannels::FillOverallQnVectorCorrectionStepList(std::set<CorrectionBase *> &set) const {
  fQnVectorCorrections.FillOverallCorrectionsList(set);
}

/// Builds raw Qn vector before Q vector corrections and before input
/// data corrections but considering the chosen calibration method.
/// This is a channelized configuration so this Q vector will NOT be
/// the one to be used for subsequent Q vector corrections.
void SubEventChannels::BuildRawQnVector() {
  fRawQnVector.SetNormalization(QVector::Normalization::NONE);
  for (const auto &dataVector : fDataVectorBank) {
    fRawQnVector.Add(dataVector.Phi(), dataVector.RadialOffset(), dataVector.Weight());
  }
  fRawQnVector.CheckQuality();
  fRawQnVector.Normal(fDetector->GetNormalizationMethod());
}

/// Ask for processing corrections for the involved detector configuration
///
/// The request is transmitted to the incoming data correction steps
/// and then to Q vector correction steps.
/// The first not applied correction step breaks the loop and kFALSE is returned
/// \return kTRUE if all correction steps were applied
Bool_t SubEventChannels::ProcessCorrections() {
  /* first we build the raw Q vector with the chosen calibration */
  BuildRawQnVector();
  /* then we transfer the request to the input data correction steps */
  for (auto &correction : fInputDataCorrections) {
    if (correction->ProcessCorrections()) {
      continue;
    } else {
      return kFALSE;
    }
  }
  /* input corrections were applied so let's build the Q vector with the chosen calibration */
  BuildQnVector();
  /* now let's propagate it to Q vector corrections */
  for (auto &correction : fQnVectorCorrections) {
    if (correction->ProcessCorrections())
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}

/// Ask for processing corrections data collection for the involved detector configuration
///
/// The request is transmitted to the incoming data correction steps
/// and then to Q vector correction steps.
/// The first not applied correction step should break the loop after collecting the data and kFALSE is returned
/// \return kTRUE if all correction steps were applied
Bool_t SubEventChannels::ProcessDataCollection() {
  /* we transfer the request to the input data correction steps */
  for (auto &correction : fInputDataCorrections) {
    if (correction->ProcessDataCollection()) {
      continue;
    } else {
      return kFALSE;
    }
  }
  /* check whether QA histograms must be filled */
  FillQAHistograms();
  /* now let's propagate it to Q vector corrections */
  for (auto &correction : fQnVectorCorrections) {
    if (correction->ProcessDataCollection())
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}

/// Clean the configuration to accept a new event
///
/// Transfers the order to the Q vector correction steps then
/// to the input data correction steps and finally
/// cleans the own Q vector and the input data vector bank
/// for accepting the next event.
void SubEventChannels::Clear() {
  for (auto &correction : fQnVectorCorrections) {
    correction->ClearCorrectionStep();
  }
  for (auto &correction : fInputDataCorrections) {
    correction->ClearCorrectionStep();
  }
  /* clean the raw Q vector */
  fRawQnVector.Reset();
  /* clean the own Q vector */
  fPlainQnVector.Reset();
  fPlainQ2nVector.Reset();
  fCorrectedQnVector.Reset();
  fCorrectedQ2nVector.Reset();
  /* and now clear the the input data bank */
  fDataVectorBank.clear();
}

}// namespace Qn
