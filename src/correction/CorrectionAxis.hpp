#ifndef QNCORRECTIONS_EVENTCLASSVAR_H
#define QNCORRECTIONS_EVENTCLASSVAR_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsEventClassVariable.h
/// \brief Class that models variables used for defining an event class within the Q vector correction framework

/// \class QnCorrectionsEventClassVariable
/// \brief One variable used for defining an event class
///
/// Class defining one variable and its associated binning allowing
/// its use for the definition of event classes within the Q vector
/// correction framework.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 4, 2016


#include "Axis.hpp"
#include "InputVariable.hpp"
#include "InputVariableManager.hpp"

namespace Qn {
class CorrectionAxis {
 public:
  CorrectionAxis() = default;
  CorrectionAxis(const AxisD &axis) : variable_(), axis_(axis) {}

  virtual ~CorrectionAxis() = default;
  void Initialize(const InputVariableManager &man) {
    variable_ = man.FindVariable(axis_.Name());
  }
  /// Gets the variable unique Id
  int GetId() const { return variable_.GetID(); }
  /// Gets the value of the variable
  double GetValue() const { return *variable_.begin(); }
  /// Gets the variable name / label
  std::string GetLabel() const { return axis_.Name(); }
  /// Gets the number of bins
  Int_t GetNBins() const { return axis_.GetNBins(); }
  /// Gets the actual bins edges array
  const double *GetBins() const { return axis_.GetPtr(); }
  /// Gets the lower edge for the passed bin number
  /// \param bin bin number starting from one
  double GetBinLowerEdge(const unsigned long bin) const { return axis_.GetLowerBinEdge(bin); }
  /// Gets the upper edge for the passed bin number
  /// \param bin bin number starting from one
  double GetBinUpperEdge(const unsigned long bin) const { return axis_.GetUpperBinEdge(bin); }
  /// Gets the lowest variable value considered
  double GetLowerEdge() const { return axis_.GetFirstBinEdge(); }
  /// Gets the highest variabel value considered
  double GetUpperEdge() const { return axis_.GetLastBinEdge(); }
 private:
  InputVariable variable_;
  AxisD axis_;
/// \cond CLASSIMP
 ClassDef(CorrectionAxis, 1);
/// \endcond
};
}
#endif /* QNCORRECTIONS_EVENTCLASSVAR_H */
