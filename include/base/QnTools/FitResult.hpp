//
// Created by eugene on 17/05/2021.
//

#ifndef QNTOOLS_INCLUDE_BASE_QNTOOLS_FITRESULT_HPP
#define QNTOOLS_INCLUDE_BASE_QNTOOLS_FITRESULT_HPP

#include <string>

#include "Rtypes.h"
#include "TFitResultPtr.h"

namespace Qn {

class FitResult {
 public:
  unsigned int NPar() const;
  std::string ParName(unsigned int i) const;

  double Parameter(unsigned int i) const;
  double Value(unsigned int i) const;
  double ParErrorFromBootStrap(unsigned int i) const;

 private:
  TFitResultPtr mean_fit_result_;
  
  std::vector<TFitResultPtr> sample_fit_results_;
  std::vector<double> sample_weights_;

  ClassDef(Qn::FitResult, 0);
};

}

#endif  // QNTOOLS_INCLUDE_BASE_QNTOOLS_FITRESULT_HPP
