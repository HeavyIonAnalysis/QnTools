#ifndef QNCORRECTIONS_PROFILECHANNEL_H
#define QNCORRECTIONS_PROFILECHANNEL_H

/// \file QnCorrectionsProfileChannelized.h
/// \brief Channelized profile class for the Q vector correction framework

#include "CorrectionHistogramBase.hpp"
namespace Qn {
/// \class QnCorrectionsProfileChannelized
/// \brief Channelized profile class for the Q vector correction histograms
///
/// Encapsulates a multidimensional profile. Each dimension
/// corresponds to one of the event classes variables so,
/// the number of dimensions matches the number of variables within
/// the set passed in the constructor. Additionally incorporates an
/// extra dimension to consider the channel number
///
/// The involved histograms can be created on the fly when needed,
/// and included in a provided list. They cannot be attached to existing
/// ones. In any case they are not destroyed because
/// they are not own by the class but by the involved list.
///
/// Storage efficiency reasons dictate that channels were stored in
/// consecutive order although externally to the class everything is
/// handled with the actual external channel number. But if the
/// histograms stored in a file are draw the channels will appear
/// as enumerated form 0 to the number of active channels handled
/// by the detector configuration that is associated to the histogram
/// and as such by the own histogram.
///
/// GetBinContent (once the intended bin is obtained by mean
/// of GetBin) returns in the profile way
/// \f[
///    \frac{\Sigma \mbox{fValues(bin)}}{\mbox{fEntries(bin)}}
/// \f]
/// while depending on the option passed at construction GetBinError returns
/// * "" the standard error on the mean of the values in the interested bin
/// \f[
///    \frac{\sqrt{\frac{\Sigma \mbox{fValues}^2\mbox{(bin)}}{\mbox{fEntries(bin)}}
///          - \left(\frac{\Sigma \mbox{fValues(bin)}}{\mbox{fEntries(bin)}}\right)^2}}
///         {\sqrt{\mbox{fEntries(bin)}}}
/// \f]
/// * "s" the standard deviation of the values
/// in the interested bin
/// \f[
///    \sqrt{\frac{\Sigma \mbox{fValues}^2\mbox{(bin)}}{\mbox{fEntries(bin)}}
///          - \left(\frac{\Sigma \mbox{fValues(bin)}}{\mbox{fEntries(bin)}}\right)^2}
/// \f]
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 11, 2016
class CorrectionProfileChannelized : public CorrectionHistogramBase {
 public:
  CorrectionProfileChannelized() = default;
  CorrectionProfileChannelized(std::string name,
                               std::string title,
                               const CorrectionAxisSet &ecvs,
                               Int_t nNoOfChannels,
                               ErrorMode mode = ErrorMode::MEAN);
  CorrectionProfileChannelized(std::string name,
                               const CorrectionAxisSet &ecvs,
                               Int_t nNoOfChannels,
                               ErrorMode mode = ErrorMode::MEAN);

  virtual ~CorrectionProfileChannelized();
  Bool_t CreateProfileHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup);
  Long64_t GetBin(Int_t nChannel);
  Bool_t BinContentValidated(Long64_t bin);
  Float_t GetBinContent(Long64_t bin);
  Float_t GetBinError(Long64_t bin);
  void Fill(Int_t nChannel, Float_t weight);
 private:
  THnF *fValues = nullptr;              //!<! Cumulates values for each of the event classes
  THnI *fEntries = nullptr;             //!<! Cumulates the number on each of the event classes
  Bool_t *fUsedChannel = nullptr;       //!<! array, which of the detector channels is used for this configuration
  Int_t *fChannelGroup = nullptr;       //!<! array, the group to which the channel pertains
  Int_t fNoOfChannels = 0;        //!<! The number of channels associated to the whole detector
  Int_t fActualNoOfChannels = 0;  //!<! The actual number of channels handled by the histogram
  Int_t *fChannelMap = nullptr;         //!<! array, the map from histo to detector channel number

  /// \cond CLASSIMP
 ClassDef(CorrectionProfileChannelized, 1);
  /// \endcond
};
}
#endif
