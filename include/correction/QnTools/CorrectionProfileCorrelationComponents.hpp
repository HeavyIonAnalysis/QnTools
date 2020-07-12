#ifndef QNCORRECTIONS_PROFILECORRCOMP_H
#define QNCORRECTIONS_PROFILECORRCOMP_H

/// \file QnCorrectionsProfileCorrelationComponents.h
/// \brief Correlation components based set of profiles for the Q vector correction framework

#include "CorrectionHistogramBase.hpp"
namespace Qn {
/// \class QnCorrectionsProfileCorrelationComponents
/// \brief Base class for the correlation components based set of profiles
///
/// Provides profile histograms for storing component, XX, XY, YX, YY, based
/// information.
///
/// As any histogram derived from QnCorrectionsHistogramBase the set
/// of variables that identify the different event classes has to
/// be passed to the constructor.
/// Of course,  the base name and base title for the different
/// histograms has also to be provided.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date May 08, 2016
class CorrectionProfileCorrelationComponents : public CorrectionHistogramBase {
 public:
  CorrectionProfileCorrelationComponents() = default;
  CorrectionProfileCorrelationComponents(const std::string name,
                                         const CorrectionAxisSet &ecvs,
                                         ErrorMode mode = ErrorMode::MEAN);
  CorrectionProfileCorrelationComponents(std::string name,
                                         std::string title,
                                         const CorrectionAxisSet &ecvs,
                                         ErrorMode mode = ErrorMode::MEAN);
  virtual ~CorrectionProfileCorrelationComponents() = default;
  Bool_t CreateCorrelationComponentsProfileHistograms(TList *histogramList);
  Bool_t AttachHistograms(TList *histogramList);
  Long64_t GetBin();
  Bool_t BinContentValidated(Long64_t bin);
  Float_t GetXXBinContent(Long64_t bin);
  Float_t GetXYBinContent(Long64_t bin);
  Float_t GetYXBinContent(Long64_t bin);
  Float_t GetYYBinContent(Long64_t bin);
  Float_t GetXXBinError(Long64_t bin);
  Float_t GetXYBinError(Long64_t bin);
  Float_t GetYXBinError(Long64_t bin);
  Float_t GetYYBinError(Long64_t bin);
  void FillXX(Float_t weight);
  void FillXY(Float_t weight);
  void FillYX(Float_t weight);
  void FillYY(Float_t weight);
 private:
  THnF *fXXValues = nullptr;           //!<! XX component histogram
  THnF *fXYValues = nullptr;           //!<! XY component histogram
  THnF *fYXValues = nullptr;           //!<! YX component histogram
  THnF *fYYValues = nullptr;           //!<! YY component histogram
  UInt_t fXXXYYXYYFillMask = 0x0000;   //!<! keeps track of component filled values
  UInt_t fFullFilled = 0x0000;         //!<! mask for the fully filled condition
  THnI *fEntries = nullptr;            //!<! Cumulates the number on each of the event classes
  /// \cond CLASSIMP
 ClassDef(CorrectionProfileCorrelationComponents, 1);
  /// \endcond
};
}
#endif
