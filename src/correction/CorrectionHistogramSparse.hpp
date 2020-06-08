#ifndef QNCORRECTIONS_HISTOGRAMSPARSE_H
#define QNCORRECTIONS_HISTOGRAMSPARSE_H

/// \file QnCorrectionsHistogramSparse.h
/// \brief Single multidimensional sparse histograms

#include "CorrectionHistogramBase.hpp"
#include <THnSparse.h>
namespace Qn {
/// \class QnCorrectionsHistogramSparse
/// \brief Single sparse histogram class for the Q vector correction histograms
///
/// Encapsulates a multi dimensional sparse histogram. Each dimension
/// corresponds to one of the event classes variables so,
/// the number of dimensions matches the number of variables within
/// the set passed in the constructor.
///
/// The involved histograms can be created on the fly when needed,
/// and included in a provided list. They are not destroyed because
/// the are not own by the class but by the involved list.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jun 18, 2016
class CorrectionHistogramSparse : public CorrectionHistogramBase {
 public:
  CorrectionHistogramSparse() = default;
  CorrectionHistogramSparse(const std::string name, const std::string title, const CorrectionAxisSet &ecvs);
  CorrectionHistogramSparse(const std::string name, const CorrectionAxisSet &ecvs);
  virtual ~CorrectionHistogramSparse() = default;
  Bool_t CreateHistogram(TList *histogramList);
  Long64_t GetBin();
  Bool_t BinContentValidated(Long64_t bin);
  Float_t GetBinContent(Long64_t bin);
  Float_t GetBinError(Long64_t bin);
  virtual void Fill(Float_t weight);
 private:
  THnSparseF *fValues = nullptr; //!<! Cumulates values for each of the event classes
  /// \cond CLASSIMP
 ClassDef(CorrectionHistogramSparse, 1);
  /// \endcond
};
}
#endif
