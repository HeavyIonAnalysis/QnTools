#ifndef QNCORRECTIONS_PROFILE3DCORRELATIONS_H
#define QNCORRECTIONS_PROFILE3DCORRELATIONS_H

/// \file QnCorrectionsProfile3DCorrelations.h
/// \brief Three detector correlation components based set of profiles with harmonic support for the Q vector correction framework

#include "CorrectionHistogramBase.hpp"
namespace Qn {
class QVector;

/// \class QnCorrectionsProfile3DCorrelations
/// \brief Base class for three detectors correlation components based set of profiles with harmonic support
///
/// Provides profile histograms for storing component, XX, XY, YX, YY, based
/// information for a set of predefined harmonics defined at creation
/// time and for a set of three different subdetectors generically identified
/// as A, B and C. The user can select the harmonic addressing procedure so that
/// it will possible to ask for just one harmonic support and assign
/// to it any desired number.
///
/// As any histogram derived from QnCorrectionsHistogramBase the set
/// of variables that identify the different event classes has to
/// be passed to the constructor together with the required number of
/// harmonics and an optional harmonic expected numbering scheme.
/// Of course,  the base name and base title for the different
/// histograms has also to be provided.
///
/// The harmonic map passed should contain an ordered array with
/// as many items as requested harmonics that provides the external
/// number to be used for request the corresponding harmonic.
/// Requesting five harmonics without maps is equivalent to pass
/// {1,2,3,4,5} as map. Requesting just support for the harmonic
/// four will require a map {4}.
///
/// Externally the harmonic number is addressed as usual. An additional
/// harmonic multiplier field allows to handle mxn harmonics. n is always
/// the external harmonic required internally it si handled as well as
/// n but all the information manipulated is really associated to mxn.
/// Only in the histograms name it appears the proper mxn harmonic to
/// not confuse the external user which browse the histograms.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 19, 2016
class CorrectionProfile3DCorrelations : public CorrectionHistogramBase {
 public:
  CorrectionProfile3DCorrelations() = default;
  CorrectionProfile3DCorrelations(
      std::string name,
      std::string title,
      std::string nameA,
      std::string nameB,
      std::string nameC,
      const CorrectionAxisSet &ecvs,
      ErrorMode mode = ErrorMode::MEAN);
  CorrectionProfile3DCorrelations(
      std::string name,
      std::string nameA,
      std::string nameB,
      std::string nameC,
      const CorrectionAxisSet &ecvs,
      ErrorMode mode = ErrorMode::MEAN);
  virtual ~CorrectionProfile3DCorrelations();
  Bool_t CreateCorrelationComponentsProfileHistograms(TList *histogramList,
                                                      Int_t nNoOfHarmonics,
                                                      Int_t nHarmonicMultiplier = 1,
                                                      Int_t *harmonicMap = NULL);
  virtual Bool_t AttachHistograms(TList *histogramList);
  virtual Long64_t GetBin();
  virtual Bool_t BinContentValidated(Long64_t bin);
  virtual Float_t GetXXBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinContent(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetXXBinError(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetXYBinError(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYXBinError(const char *comb, Int_t harmonic, Long64_t bin);
  virtual Float_t GetYYBinError(const char *comb, Int_t harmonic, Long64_t bin);
  void Fill(const QVector *QnA, const QVector *QnB,
            const QVector *QnC);
 private:
  THnF ***fXXValues = nullptr;            //!<! XX component histogram for each requested harmonic
  THnF ***fXYValues = nullptr;            //!<! XY component histogram for each requested harmonic
  THnF ***fYXValues = nullptr;            //!<! YX component histogram for each requested harmonic
  THnF ***fYYValues = nullptr;            //!<! YY component histogram for each requested harmonic
  THnI *fEntries = nullptr;             //!<! Cumulates the number on each of the event classes
  std::string fNameA;               ///< the name of the A detector
  std::string fNameB;               ///< the name of the B detector
  std::string fNameC;               ///< the name of the C detector
  Int_t fHarmonicMultiplier = 1;    ///< the multiplier for the harmonic number
  /// \cond CLASSIMP
 ClassDef(CorrectionProfile3DCorrelations, 1);
  /// \endcond
};
}
#endif
