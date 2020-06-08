#ifndef QNCORRECTIONS_EVENTCLASSVARSET_H
#define QNCORRECTIONS_EVENTCLASSVARSET_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file QnCorrectionsEventClassVariablesSet.h
/// \brief Class that models the set of variables that define an event class for the Q vector correction framework

#include "CorrectionAxis.hpp"
#include "InputVariableManager.hpp"
namespace Qn {
class CorrectionAxisSet {
  using size_type = std::vector<CorrectionAxis>::size_type;
 public:
  CorrectionAxisSet() = default;
  explicit CorrectionAxisSet(int initialsize) : axes_(initialsize) {}
  virtual ~CorrectionAxisSet() = default;
  std::vector<CorrectionAxis>::iterator begin() { return axes_.begin(); }
  std::vector<CorrectionAxis>::iterator end() { return axes_.end(); }
  std::vector<CorrectionAxis>::const_iterator begin() const { return axes_.begin(); }
  std::vector<CorrectionAxis>::const_iterator end() const { return axes_.end(); }
  CorrectionAxis &operator[](size_type i) { return axes_.at(i); }
  const CorrectionAxis &operator[](size_type i) const { return axes_.at(i); }
  CorrectionAxis &At(unsigned int i) { return axes_.at(i); }
  const CorrectionAxis &At(unsigned int i) const { return axes_.at(i); }
  template<typename... ARGS>
  void Add(ARGS &&... args) { axes_.emplace_back(std::forward<ARGS>(args)...); }
  void Initialize(const InputVariableManager &variable_manager) {
    for (auto &axis : axes_) {
      axis.Initialize(variable_manager);
    }
  }
  std::vector<CorrectionAxis>::size_type GetSize() const { return axes_.size(); }
  void GetMultidimensionalConfiguration(int *nbins, double *minvals, double *maxvals) const {
    unsigned int i = 0;
    for (auto &axis : axes_) {
      nbins[i] = axis.GetNBins();
      minvals[i] = axis.GetLowerEdge();
      maxvals[i] = axis.GetUpperEdge();
      ++i;
    }
  }
 private:
  std::vector<CorrectionAxis> axes_;
/// \cond CLASSIMP
 ClassDef(CorrectionAxisSet, 1);
/// \endcond
};
}
#endif /* QNCORRECTIONS_EVENTCLASSVARSET_H */
