#ifndef QNCORRECTIONS_DETECTORCONFCHANNEL_H
#define QNCORRECTIONS_DETECTORCONFCHANNEL_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsDetectorConfigurationChannels.h
/// \brief Channel detector configuration class for Q vector correction framework
///

#include "CorrectionsSet.hpp"
#include "SubEvent.hpp"
#include "CorrectionDataVector.hpp"

namespace Qn {
class CorrectionProfileComponents;

/// \class QnCorrectionsDetectorConfigurationChannels
/// \brief Channel detector configuration within Q vector correction framework
///
/// A channel detector within the Q vector correction framework is defined
/// as one for which its data vectors involve azimuthal angles and channels
/// susceptible of weighting and / or grouping and / or calibration, etc.
///
/// According to that, the proper channelized data vector is used and an extra
/// Q vector builder is incorporated.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016

class SubEventChannels : public SubEvent {
 public:
  friend class CorrectionBase;
  friend class SubEvent;
  SubEventChannels() = default;
  SubEventChannels(unsigned int bin_id,
                   const CorrectionAxisSet *axes, Int_t nNoOfChannels, std::bitset<QVector::kmaxharmonics> harmonics);
  virtual ~SubEventChannels();
  SubEventChannels(const SubEventChannels &) = delete;
  SubEventChannels &operator=(const SubEventChannels &) = delete;

  /// Gets the number of channels
  /// \return the number of channels of the associated detector
  Int_t GetNoOfChannels() { return fNoOfChannels; }
  /// Gets the used channels mask
  /// \return the used channels mask
  const Bool_t *GetUsedChannelsMask() const { return fUsedChannel; }
  /// Gets the channels groups
  /// \return the group associated to each channel
  const Int_t *GetChannelsGroups() const { return fChannelGroup; }
  /// Gets the hard coded group weights
  /// \return the groups hard coded weights
  const Float_t *GetHardCodedGroupWeights() const { return fHardCodedGroupWeights; }
  /// Get if the detector configuration is own by a tracking detector
  /// \return FALSE, this is a hit / channel detector configuration
  Bool_t GetIsTrackingDetector() const { return kFALSE; }

  virtual void SetChannelsScheme(std::vector<int> channel_groups);

  /* QA section */
  /// Sets the variable id used for centrality in QA histograms.
  /// It must be one of the event class variables.
  /// \param id id for the variable used for centrality
  void SetQACentralityVar(Int_t id) { fQACentralityVarId = id; }
  /// Sets the characteristics of the multiplicity axis in QA histograms
  /// \param nbins the number of bins
  /// \param min minimum multiplicity value
  /// \param max maximum multiplicity value
  void SetQAMultiplicityAxis(Int_t nbins, Float_t min, Float_t max) {
    fQAnBinsMultiplicity = nbins;
    fQAMultiplicityMin = min;
    fQAMultiplicityMax = max;
  }

  void BuildRawQnVector();

  virtual void CreateSupportQVectors();
  virtual void CreateCorrectionHistograms();
  virtual void CopyToOutputList(TList* list);
  virtual void AttachQAHistograms(TList *list);
  virtual void AttachNveQAHistograms(TList *list);

  /// Activate the processing for the passed harmonic
  /// \param harmonic the desired harmonic number to activate
  virtual void ActivateHarmonic(Int_t harmonic) {
    SubEvent::ActivateHarmonic(harmonic);
    fRawQnVector.ActivateHarmonic(harmonic);
  }
  virtual void AttachCorrectionInput(TList *list);
  virtual void AfterInputAttachAction();
  virtual Bool_t ProcessCorrections();
  virtual Bool_t ProcessDataCollection();
  virtual void AddCorrectionOnInputData(CorrectionOnInputData *correctionOnInputData);
  virtual void IncludeQnVectors();
  virtual void FillOverallInputCorrectionStepList(std::set<CorrectionBase *> &set) const;
  virtual void FillOverallQnVectorCorrectionStepList(std::set<CorrectionBase *> &set) const;
  virtual void Clear();

  virtual std::map<std::string, Report> ReportOnCorrections() const {
    auto report = fInputDataCorrections.ReportOnUsage();
    auto qn_correction_report = fQnVectorCorrections.ReportOnUsage();
    report.insert(qn_correction_report.begin(), qn_correction_report.end());
    return report;
  }

 private:
  static const char
      *szRawQnVectorName;   ///< the name of the raw Qn vector from raw data without input data corrections
  QVector fRawQnVector;                       ///< Q vector from input data before pre-processing
  Int_t fNoOfChannels = 0;                    ///< The number of channels associated
  Bool_t *fUsedChannel =
      nullptr;                   //[fNoOfChannels]   /// array, which of the detector channels is used for this configuration
  Int_t *fChannelMap =
      nullptr;                     //[fNoOfChannels]   /// array, mapping external to internal channel id.
  Int_t
      *fChannelGroup = nullptr;                   //[fNoOfChannels]   /// array, the group to which the channel pertains
  Float_t *fHardCodedGroupWeights = nullptr;         //[fNoOfChannels]  /// array, group hard coded weight
  CorrectionsSetOnInputData fInputDataCorrections; ///< set of corrections to apply on input data vectors

  /* QA section */
  void FillQAHistograms();
  static const char *szQAMultiplicityHistoName; ///< QA multiplicity histograms name
  Int_t fQACentralityVarId = -1;   ///< the id of the variable used for centrality in QA histograms
  Int_t fQAnBinsMultiplicity = 100; ///< number of bins for multiplicity in QA histograms
  Float_t fQAMultiplicityMin = 0.0; ///< minimum multiplicity value
  Float_t fQAMultiplicityMax = 1000.0; ///< maximum multiplicity value
  TH3F *fQAMultiplicityBefore3D = nullptr; //!<! 3D channel multiplicity histogram before input equalization
  TH3F *fQAMultiplicityAfter3D = nullptr;  //!<! 3D channel multiplicity histogram after input equalization

/// \cond CLASSIMP
 ClassDef(SubEventChannels, 2);
/// \endcond
};

}
#endif // QNCORRECTIONS_DETECTORCONFCHANNEL_H
