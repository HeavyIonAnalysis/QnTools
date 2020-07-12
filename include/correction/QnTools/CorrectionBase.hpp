#ifndef QNCORRECTIONS_CORRECTIONSTEPBASE_H
#define QNCORRECTIONS_CORRECTIONSTEPBASE_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsCorrectionStepBase.h
/// \brief Base class for the support of the different correction steps within Q vector correction framework
///

#include "TObject.h"
#include "TList.h"

namespace Qn {
class SubEvent;
class QVector;


/// \class QnCorrectionsCorrectionStepBase
/// \brief Base class for correction steps
///
/// Each correction has a name and a key. The name identifies it
/// in an open way while the key is used to codify its position
/// in an ordered list of consecutive corrections.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 05, 2016

class CorrectionBase : public TObject {
 public:
  /// \typedef QnCorrectionStepStatus
  /// \brief The class of the id of the correction steps states
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  ///
  /// When referring as "data being collected" means that the needed data
  /// for producing new correction parameters are being collected.

  enum class State {
    PASSIVE,             ///< the correction step is waiting for external conditions fulfillment
    CALIBRATION,         ///< the correction step is in calibration mode collecting data
    APPLY,               ///< the correction step is being applied
    APPLYCOLLECT,        ///< the correction step is being applied and data are being collected
  };

  friend class SubEvent;
  CorrectionBase() = default;
  virtual ~CorrectionBase() = default;
  CorrectionBase(const char *name, unsigned int prio) : fPriority(prio), fName(name) {
    output_histograms.SetOwner(true);
  }
  virtual const char *GetName() const { return fName.data(); }
  CorrectionBase(const CorrectionBase &other) : TObject(other) {
    fPriority = other.fPriority;
    fName = other.fName;
    fState = other.fState;
    fSubEvent = other.fSubEvent;
  }

  /// Attaches the needed input information to the correction step
  ///
  /// Pure virtual function
  /// \param list list where the inputs should be found
  /// \return kTRUE if everything went OK
  virtual void AttachInput(TList *list) { (void) list; }
  /// Perform after calibration histograms attach actions
  /// It is used to inform the different correction step that
  /// all conditions for running the network are in place so
  /// it is time to check if their requirements are satisfied
  ///
  /// Pure virtual function
  virtual void AfterInputAttachAction() {}
  /// Asks for support data structures creation
  ///
  /// Pure virtual function
  virtual void CreateSupportQVectors() {}
  /// Asks for support histograms creation
  ///
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual void CreateCorrectionHistograms() {}
  /// Asks for QA histograms creation
  ///
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK

  virtual void AttachQAHistograms(TList *list) { (void) list; }
  /// Asks for non validated entries QA histograms creation
  ///
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual void AttachNveQAHistograms(TList *list) { (void) list; }
  /// Processes the correction step
  ///
  /// Pure virtual function
  /// \return kTRUE if everything went OK
  virtual bool ProcessCorrections() { return false; }
  /// Processes the correction step data collection
  ///
  /// Pure virtual function
  /// \return kTRUE if everything went OK
  virtual bool ProcessDataCollection() { return false; }
  /// Clean the correction to accept a new event
  /// Pure virtual function
  virtual void ClearCorrectionStep() {}
  /// Reports if the correction step is being applied
  /// Pure virutal function
  /// \return TRUE if the correction step is being applied
  virtual Bool_t IsBeingApplied() const { return false; }
  /// Report on correction usage
  /// Pure virtual function
  /// Correction step should incorporate its name in calibration
  /// list if it is producing information calibration in the ongoing
  /// step and in the apply list if it is applying correction in
  /// the ongoing step.
  /// \param calibrationList list containing the correction steps producing calibration information
  /// \param applyList list containing the correction steps applying corrections
  /// \return kTRUE if the correction step is being applied
  virtual std::pair<bool, bool> ReportUsage() {
    bool applying = false;
    bool collecting = false;
    switch (fState) {
      case State::CALIBRATION:collecting = true;
        applying = false;
        break;
      case State::APPLYCOLLECT:collecting = true;
        /* FALLTHRU */
      case State::APPLY:applying = true;
        break;
      case State::PASSIVE:collecting = false;
        applying = false;
        break;
    }
    return std::make_pair(collecting, applying);
  }
  State GetState() { return fState; }
  void Enable() {fState = State::CALIBRATION;}

  void CopyToOutputList(TList *list) {
    if (fState == State::PASSIVE) return;
    output_histograms.SetOwner(false);
    list->AddAll(&output_histograms);
    output_histograms.Clear();
  }
 protected:
/// Stores the detector configuration owner
/// \param subevent the detector configuration owner
  void SetOwner(SubEvent *subevent) { fSubEvent = subevent; }
  unsigned int fPriority = 0; ///< the correction key that codifies order information
  std::string fName;
  State fState = State::PASSIVE; ///< the state in which the correction step is
  SubEvent *fSubEvent = nullptr; ///< pointer to the detector configuration owner
  TList output_histograms;
  friend bool operator<(const CorrectionBase &lh, const CorrectionBase &rh);

/// \cond CLASSIMP
 ClassDef(CorrectionBase, 2);
/// \endcond
};
inline bool operator<(const Qn::CorrectionBase &lh, const Qn::CorrectionBase &rh) {
  return lh.fPriority < rh.fPriority;
}

}
#endif // QNCORRECTIONS_CORRECTIONSTEPBASE_H
