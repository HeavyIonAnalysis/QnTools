#ifndef QNCORRECTIONS_CORRECTIONONINPUTDATA_H
#define QNCORRECTIONS_CORRECTIONONINPUTDATA_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsCorrectionOnInputData.h
/// \brief Correction steps on input data support within Q vector correction framework
///

#include <TNamed.h>
#include <TList.h>

#include "CorrectionBase.hpp"

namespace Qn {

/// \class QnCorrectionsCorrectionOnInputData
/// \brief Base class for correction steps applied to input data
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 05, 2016

class CorrectionOnInputData : public CorrectionBase {
 public:
  enum Priority {
    kGainEqualization
  };
  friend class SubEventChannels;
  CorrectionOnInputData() = default;
  CorrectionOnInputData(const char *name, unsigned int prio) : CorrectionBase(name, prio) {}
  virtual ~CorrectionOnInputData() = default;
  CorrectionOnInputData(const CorrectionOnInputData &other) : CorrectionBase(other) {}
  virtual CorrectionOnInputData *MakeCopy() const { return new CorrectionOnInputData(*this); }
  /// Reports if the correction step is being applied
  /// \return FALSE, input data correction step dont make use of this service, yet
  virtual Bool_t IsBeingApplied() const {
    bool applied = false;
    switch (fState) {
      case State::CALIBRATION: break;
      case State::APPLYCOLLECT:
        /* FALLTHRU */
      case State::APPLY: applied = true;
        break;
      case State::PASSIVE: break;
    }
    return applied;
  }
/// \cond CLASSIMP
 ClassDef(CorrectionOnInputData, 1);
/// \endcond
};
}

#endif // QNCORRECTIONS_CORRECTIONONINPUTDATA_H
