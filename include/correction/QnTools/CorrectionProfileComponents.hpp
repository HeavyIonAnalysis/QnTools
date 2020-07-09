#ifndef QNCORRECTIONS_PROFILECOMP_H
#define QNCORRECTIONS_PROFILECOMP_H

/// \file QnCorrectionsProfileComponents.h
/// \brief Component based set of profiles for the Q vector correction framework

#include "CorrectionHistogramBase.hpp"
namespace Qn {
/// \class QnCorrectionsProfileComponents
/// \brief Base class for the components based set of profiles
///
/// Provides profile histograms for storing component, X, Y, based
/// information for a set of predefined harmonics defined at creation
/// time. The user can select the harmonic addressing procedure so that
/// it will possible to ask for just one harmonic support and assign
/// to it any desired number.
///
/// As any histogram derived from QnCorrectionsHistogramBase the set
/// of variables that identify the different event classes has to
/// be passed to the constructor. Of course,  the base name and base
/// title for the different histograms has to be also provided.
/// At creation time the required number of harmonics and an optional
/// expected harmonic numbering scheme has to be passed.
///
/// The harmonic map passed should contain an ordered array with
/// as many items as requested harmonics that provides the external
/// number to be used for request the corresponding harmonic.
/// Requesting five harmonics without maps is equivalent to pass
/// {1,2,3,4,5} as map. Requesting just support for the harmonic
/// four will require a map {4}.
///
/// When you fill the histograms care is taken for not increasing
/// the number of entries until all components for the whole set of
/// harmonics have been filled. If you try to fill twice a harmonic
/// component before the whole set is filled you will get an execution
/// error because you are doing something that shall be corrected
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 15, 2016
class CorrectionProfileComponents : public CorrectionHistogramBase {
 public:
  CorrectionProfileComponents() = default;
  CorrectionProfileComponents(std::string name,
                              std::string title,
                              const CorrectionAxisSet &ecvs,
                              ErrorMode mode = ErrorMode::MEAN);
  CorrectionProfileComponents(std::string name, const CorrectionAxisSet &ecvs, ErrorMode mode = ErrorMode::MEAN);
  virtual ~CorrectionProfileComponents();
  Bool_t CreateComponentsProfileHistograms(TList *histogramList, Int_t nNoOfHarmonics, Int_t *harmonicMap = NULL);
  Bool_t AttachHistograms(TList *histogramList);
  Long64_t GetBin();
  Bool_t BinContentValidated(Long64_t bin);
  Float_t GetXBinContent(Int_t harmonic, Long64_t bin);
  Float_t GetYBinContent(Int_t harmonic, Long64_t bin);
  Float_t GetXBinError(Int_t harmonic, Long64_t bin);
  Float_t GetYBinError(Int_t harmonic, Long64_t bin);
  void FillX(Int_t harmonic, Float_t weight);
  void FillY(Int_t harmonic, Float_t weight);
 private:
  THnF **fXValues = nullptr;            //!<! X component histogram for each requested harmonic
  THnF **fYValues = nullptr;            //!<! Y component histogram for each requested harmonic
  UInt_t fXharmonicFillMask = 0x0000;  //!<! keeps track of harmonic X component filled values
  UInt_t fYharmonicFillMask = 0x0000;  //!<! keeps track of harmonic Y component filled values
  UInt_t fFullFilled = 0x0000;         //!<! mask for the fully filled condition
  THnI *fEntries = nullptr;            //!<! Cumulates the number on each of the event classes
  /// \cond CLASSIMP
 ClassDef(CorrectionProfileComponents, 1);
  /// \endcond
};
}
#endif
