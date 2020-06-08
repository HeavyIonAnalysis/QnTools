#ifndef QNCORRECTIONS_PROFILECHANNELINGRESS_H
#define QNCORRECTIONS_PROFILECHANNELINGRESS_H

/// \file QnCorrectionsProfileChannelizedIngress.h
/// \brief Ingress channelized profile for the Q vector correction framework

#include "CorrectionHistogramBase.hpp"
namespace Qn {
/// \class QnCorrectionsProfileChannelizedIngress
/// \brief Ingress channelized profile class for the Q vector correction histograms
///
/// Encapsulates a multidimensional profile. Each dimension
/// corresponds to one of the event classes variables so,
/// the number of dimensions matches the number of variables within
/// the set passed in the constructor. Additionally incorporates an
/// extra dimension to consider the channel number
///
/// The involved histograms can only be attached to existing ones
/// from a given list. The histograms to attach are the ones created from
/// QnCorrectionsProfileChannelized class. The difference with that class
/// is that now, the histograms are the source for applying corrections
/// and as such should be stored in a more efficient way and the
/// implicit group information be recovered.
///
/// Now the class will own histograms that must be deleted when the
/// class gets destroyed. The behavior regarding the channel handling
/// matches the one from QnCorrectionsProfileChannelized class.
///
/// The entries histogram disappears and so, GetBinContent and GetBinError
/// return in the standard way. Additionally, if applicable, the group
/// histogram is created and GetGrpBinContent and GetGrpBinError are
/// functional.
///
/// The profile as such cannot be filled. It should be considered as a
/// read only profile.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 23, 2016
class CorrectionProfileChannelizedIngress : public CorrectionHistogramBase {
 public:
  CorrectionProfileChannelizedIngress();
  CorrectionProfileChannelizedIngress(std::string name, std::string title,
                                      const CorrectionAxisSet &ecvs,
                                      Int_t nNoOfChannels,
                                      ErrorMode mode);
  CorrectionProfileChannelizedIngress(std::string name,
                                      const CorrectionAxisSet &ecvs,
                                      Int_t nNoOfChannels,
                                      ErrorMode mode);
  virtual ~CorrectionProfileChannelizedIngress();
  Bool_t AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup);
  Long64_t GetBin(Int_t nChannel);
  Long64_t GetGrpBin(Int_t nChannel);
  Bool_t BinContentValidated(Long64_t bin);
  Float_t GetBinContent(Long64_t bin);
  Float_t GetGrpBinContent(Long64_t bin);
  Float_t GetBinError(Long64_t bin);
  Float_t GetGrpBinError(Long64_t bin);
 private:
  THnF *fValues = nullptr;              //!<! the values and errors on each event class and channel
  THnF *fGroupValues = nullptr;         //!<! the values and errors on each event class and group
  THnC *fValidated = nullptr;            //!<! bin content validated flag
  Bool_t *fUsedChannel = nullptr;       //!<! array, which of the detector channels are used for this configuration
  Int_t *fChannelGroup = nullptr;        //!<! array, the group to which the channel pertains
  Int_t fNoOfChannels = 0;        //!<! The number of channels associated to the whole detector
  Int_t fActualNoOfChannels = 0;  //!<! The actual number of channels handled by the histogram
  Int_t *fChannelMap = nullptr;          //!<! array, the map from histo to detector channel number
  Bool_t fUseGroups = false;          //!<! the groups structures must be used
  Bool_t *fUsedGroup = nullptr;         //!<! array, which of the detector groups are used for this configuration
  Int_t fNoOfGroups = 0;          //!<! the number of groups associated with the whole detector
  Int_t fActualNoOfGroups = 0;    //!<! The actual number of groups handled by the histogram
  Int_t *fGroupMap = nullptr;           //!<! array, the map from histo to detector channel group number

  /// \cond CLASSIMP
 ClassDef(CorrectionProfileChannelizedIngress, 2);
  /// \endcond
};
}
#endif
