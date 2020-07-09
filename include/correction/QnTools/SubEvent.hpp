#ifndef QNCORRECTIONS_DETECTORCONFIGBASE_H
#define QNCORRECTIONS_DETECTORCONFIGBASE_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsDetectorConfigurationBase.h
/// \brief The base of a concrete detector configuration (sub-detector) within Q vector correction framework
///

#include <map>

#include "TObject.h"
#include "TList.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TH3.h"
#include "TObjString.h"

#include "CorrectionsSet.hpp"
#include "CorrectionAxisSet.hpp"
#include "QVector.hpp"
#include "CorrectionDataVector.hpp"
#include "CorrectionProfileComponents.hpp"

namespace Qn {
class Detector;
/// \class QnCorrectionsDetectorConfigurationBase
/// \brief The base of a concrete detector configuration within Q vector correction framework
///
/// The detector configuration shapes a detector with a concrete
/// set of cuts to make it the target of a Q vector correction process.
///
/// It receives the data input stream and build the corresponding Q vector associated
/// to it for each processing request.
///
/// As such, it incorporates the set of corrections to carry on the input data
/// and the set of corrections to perform on the produced Qn vector. It always stores
/// the plain Qn vector produced after potential input data corrections and the
/// Qn vector that incorporates the latest Qn vector correction step.
///
/// It also incorporates the equivalent support for Q2n vectors which could be the
/// seed for future Q(m,n) support.
///
/// It receives at construction time the set of event classes variables and the
/// detector reference. The reference of the detector should only be obtained at
/// creation time and the detector configuration object, once created, does not
/// allow its modification.
///
/// The class is a base class for further refined detector configurations.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016
class SubEvent {
 public:
  friend class CorrectionBase;
  SubEvent() = default;
  SubEvent(unsigned int bin_id,
           const CorrectionAxisSet *eventClassesVariables,
           std::bitset<QVector::kmaxharmonics> harmonics) :
      binid_(bin_id),
      fPlainQnVector(harmonics, QVector::CorrectionStep::PLAIN),
      fPlainQ2nVector(harmonics, QVector::CorrectionStep::PLAIN),
      fCorrectedQnVector(harmonics, QVector::CorrectionStep::PLAIN),
      fCorrectedQ2nVector(harmonics, QVector::CorrectionStep::PLAIN),
      fTempQnVector(harmonics, QVector::CorrectionStep::PLAIN),
      fTempQ2nVector(harmonics, QVector::CorrectionStep::PLAIN),
      fQnVectorCorrections() {
    fEventClassVariables = eventClassesVariables;
    fPlainQ2nVector.SetHarmonicMultiplier(2);
    fCorrectedQ2nVector.SetHarmonicMultiplier(2);
    fTempQ2nVector.SetHarmonicMultiplier(2);
  }

  virtual ~SubEvent() = default;
  SubEvent(const SubEvent &) = delete;
  SubEvent &operator=(const SubEvent &) = delete;

  std::string GetName() const;
  /// Get the input data bank.
  /// Makes it available for input corrections steps.
  /// \return pointer to the input data bank
  std::vector<Qn::CorrectionDataVector> &GetInputDataBank() { return fDataVectorBank; }
  /// Get the event class variables set
  /// Makes it available for corrections steps
  /// \return pointer to the event class variables set
  const CorrectionAxisSet &GetEventClassVariablesSet() { return *fEventClassVariables; }
  /// Get the current Qn vector
  /// Makes it available for subsequent correction steps.
  /// It could have already supported previous correction steps
  /// \return pointer to the current Qn vector instance
  QVector *GetCurrentQnVector() { return &fCorrectedQnVector; }
  const QVector *GetPreviousCorrectedQnVector(CorrectionOnQnVector *correctionOnQn) const;
  Bool_t IsCorrectionStepBeingApplied(const char *step) const;
  /// Get the current Q2n vector
  /// Makes it available for subsequent correction steps.
  /// It could have already supported previous correction steps
  /// \return pointer to the current Q2n vector instance
  QVector *GetCurrentQ2nVector() { return &fCorrectedQ2nVector; }
  /// Get the plain Qn vector
  /// Makes it available for correction steps which need it.
  /// \return pointer to the plain Qn vector instance
  QVector *GetPlainQnVector() { return &fPlainQnVector; }
  /// Get the plain Q2n vector
  /// Makes it available for correction steps which need it.
  /// \return pointer to the plain Qn vector instance
  const QVector &GetPlainQ2nVector() const { return fPlainQ2nVector; }
  /// Update the current Qn vector
  /// Update towards what is the latest values of the Qn vector after executing a
  /// correction step to make it available to further steps.
  /// \param newQnVector the new values for the Qn vector
  void UpdateCurrentQnVector(const QVector &newQnVector) { fCorrectedQnVector = newQnVector; }
  /// Update the current Q2n vector
  /// Update towards what is the latest values of the Q2n vector after executing a
  /// correction step to make it available to further steps.
  /// \param newQ2nVector the new values for the Q2n vector
  void UpdateCurrentQ2nVector(const QVector &newQ2nVector) { fCorrectedQ2nVector = newQ2nVector; }
  /// Get the number of harmonics handled by the detector configuration
  /// \return the number of handled harmonics
  Int_t GetNoOfHarmonics() const { return fCorrectedQnVector.GetNoOfHarmonics(); }
  /// Get the harmonics map handled by the detector configuration
  /// \param store pointer to the memory for storing the harmonics map
  void GetHarmonicMap(Int_t *store) const { fCorrectedQnVector.GetHarmonicsMap(store); }
  std::bitset<QVector::kmaxharmonics> GetHarmonics() const { return fCorrectedQnVector.GetHarmonics(); }
  /// Get the pointer to the framework manager
  /// \return the stored pointer to the corrections framework
  Detector *GetDetector() const { return fDetector; }
  void SetDetector(Detector *det) { fDetector = det; }
  /// Get if the detector configuration is own by a tracking detector
  /// Pure virtual function
  /// \return TRUE if it is a tracking detector configuration
  virtual Bool_t GetIsTrackingDetector() const = 0;
  /// Asks for support data structures creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  virtual void CreateSupportQVectors() = 0;

  /// Asks for support histograms creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual void CreateCorrectionHistograms() = 0;

  virtual void CopyToOutputList(TList* list) = 0;

  /// Asks for QA histograms creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual void AttachQAHistograms(TList *list) = 0;

  /// Asks for non validated entries QA histograms creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual void AttachNveQAHistograms(TList *list) = 0;

  /// Asks for attaching the needed input information to the correction steps
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the input information should be found
  /// \return kTRUE if everything went OK
  virtual void AttachCorrectionInput(TList *list) = 0;
  /// Perform after calibration histograms attach actions
  /// It is used to inform the different correction step that
  /// all conditions for running the network are in place so
  /// it is time to check if their requirements are satisfied
  ///
  /// Pure virtual function
  virtual void AfterInputAttachAction() = 0;
  /// Ask for processing corrections for the involved detector configuration
  ///
  /// Pure virtual function.
  /// The request is transmitted to the correction steps
  /// \return kTRUE if everything went OK
  virtual Bool_t ProcessCorrections() = 0;
  /// Ask for processing corrections data collection for the involved detector configuration
  ///
  /// Pure virtual function.
  /// The request is transmitted to the correction steps
  /// \return kTRUE if everything went OK
  virtual Bool_t ProcessDataCollection() = 0;
  virtual void ActivateHarmonic(Int_t harmonic);
  virtual void AddCorrectionOnQnVector(CorrectionOnQnVector *correctionOnQn);
  virtual void AddCorrectionOnInputData(CorrectionOnInputData *correctionOnInputData);

  /// Builds Qn vectors before Q vector corrections but
  /// considering the chosen calibration method.
  /// Remember, this configuration does not have a channelized
  /// approach so, the built Q vectors are the ones to be used for
  /// subsequent corrections.
  void BuildQnVector();
//  /// Include the list of associated Qn vectors into the passed list
//  ///
//  /// Pure virtual function
//  /// \param list list where the Qn vectors list should be added
  virtual void IncludeQnVectors() = 0;
  /// Include only one instance of each input correction step
  /// in execution order
  ///
  /// Pure virtual function
  /// \param list list where the correction steps should be incorporated
  virtual void FillOverallInputCorrectionStepList(std::set<CorrectionBase *> &set) const = 0;
  /// Include only one instance of each Qn vector correction step
  /// in execution order
  ///
  /// Pure virtual function
  /// \param list list where the correction steps should be incorporated
  virtual void FillOverallQnVectorCorrectionStepList(std::set<CorrectionBase *> &set) const = 0;
  /// Provide information about assigned corrections

  using Report = std::pair<bool, bool>;

  virtual std::map<std::string, Report> ReportOnCorrections() const = 0;

  /**
   * Adds a data vector to the sub event.
   * @tparam Args type of the arguments of the CorrectionDataVector constructor.
   * @param args arguments of the CorrectionDataVector constructor.
   */
  template<typename... Args>
  void AddDataVector(Args &&... args) {
    fDataVectorBank.emplace_back(std::forward<Args>(args)...);
  }
  /// Clean the configuration to accept a new event
  /// Pure virtual function
  virtual void Clear() = 0;

  virtual void SetChannelsScheme(std::vector<int> channel_groups) {
    (void) channel_groups;
  }

  QVector *GetQVector(Qn::QVector::CorrectionStep step) { return qvectors_.at(step); }

  Qn::QVector::CorrectionStep GetLatestCorrectionStep() {
    Qn::QVector::CorrectionStep latest = Qn::QVector::CorrectionStep::PLAIN;
    for (auto &correction_step : fQnVectorCorrections) {
      if (correction_step->GetState() >= Qn::CorrectionBase::State::APPLY) {
        auto temp = correction_step->GetCorrectedQnVector()->GetCorrectionStep();
        if (latest < temp) latest = temp;
      }
    }
    return latest;
  }

  std::vector<QVector::CorrectionStep> GetCorrectionSteps() {
    std::vector<QVector::CorrectionStep> steps;
    steps.emplace_back(QVector::CorrectionStep::PLAIN);
    for (auto &correction_step : fQnVectorCorrections) {
      if (correction_step->GetState() >= CorrectionBase::State::APPLY) {
        correction_step->IncludeCorrectionStep(steps);
      }
    }
    return steps;
  }

 protected:
  unsigned int binid_;
  Detector *fDetector = nullptr;
  std::vector<Qn::CorrectionDataVector> fDataVectorBank; //!<! input data for the current process / event
  QVector fPlainQnVector;      ///< Qn vector from the post processed input data
  QVector fPlainQ2nVector;     ///< Q2n vector from the post processed input data
  QVector fCorrectedQnVector;  ///< Qn vector after subsequent correction steps
  QVector fCorrectedQ2nVector; ///< Q2n vector after subsequent correction steps
  QVector fTempQnVector;  ///< temporary Qn vector for efficient Q vector building
  QVector fTempQ2nVector; ///< temporary Qn vector for efficient Q vector building
  std::map<QVector::CorrectionStep, QVector *> qvectors_;
  CorrectionsSetOnQvector fQnVectorCorrections; ///< set of corrections to apply on Q vectors
  const CorrectionAxisSet *fEventClassVariables = nullptr; //-> /// set of variables that define event classes
  std::unique_ptr<CorrectionProfileComponents> fQAQnAverageHistogram; //!<! the plain average Qn components QA histogram
  static const char *szPlainQnVectorName; ///< the name of the Qn plain, not corrected Qn vectors
  static const char *szQAQnAverageHistogramName; ///< name and title for plain Qn vector average QA histograms
  static constexpr unsigned int INITIALSIZE = 4; ///< The default initial size of data vectors banks

/// \cond CLASSIMP
 ClassDef(SubEvent, 3);
/// \endcond
};
}
#endif // QNCORRECTIONS_DETECTORCONFIGBASE_H
