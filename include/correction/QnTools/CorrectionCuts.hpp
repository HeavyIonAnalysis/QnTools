#include <utility>

// Qn Tools
//
// Copyright (C) 2020  Lukas Kreis Ilya Selyuzhenkov
// Contact: l.kreis@gsi.de; ilya.selyuzhenkov@gmail.com
// For a full list of contributors please see docs/Credits
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef FLOW_CORRECTIONCUTS_H
#define FLOW_CORRECTIONCUTS_H

#include <array>
#include <vector>
#include <functional>

#include "ROOT/RMakeUnique.hxx"
#include "ROOT/RIntegerSequence.hxx"

#include "InputVariableManager.hpp"
#include "QAHistogram.hpp"
#include "Cuts.hpp"

namespace Qn {

class CorrectionCut {
 public:
  using CallBack = std::function<std::unique_ptr<Qn::ICut>(const Qn::InputVariableManager &)>;
  explicit CorrectionCut(CallBack callback) : callback_(std::move(callback)) {}
  ~CorrectionCut() = default;
  CorrectionCut(CorrectionCut &&) = default;
  void Initialize(const Qn::InputVariableManager &var) {
    cut_ = callback_(var);
  }
  std::string Name() const {
    return cut_->Name();
  }
  bool Check(unsigned int i) const {
    return cut_->Check(i);
  }
  bool Check() const {
    return cut_->Check();
  }
  CorrectionCut::CallBack GetCallBack() const { return callback_; }
 private:
  std::unique_ptr<ICut> cut_;
  CorrectionCut::CallBack callback_;
};

/**
 * Manages cuts class and allows checking if the current variables passes the cut
 */
class CorrectionCuts {
 public:
  CorrectionCuts() = default;

  CorrectionCuts(const CorrectionCuts &cuts) {
    n_channels_ = cuts.n_channels_;
    var_values_ = cuts.var_values_;
    cut_number = cuts.cut_number;
    cut_weight_ = cuts.cut_weight_;
    cut_channel_ = cuts.cut_channel_;
    for (auto &cut : cuts.cuts_) {
      cuts_.emplace_back(cut.GetCallBack());
    }
    if (cuts.report_) {
      CreateCutReport(cuts.report_->GetName(), n_channels_);
    }
  }

  CorrectionCuts(CorrectionCuts &&cuts) {
    n_channels_ = cuts.n_channels_;
    var_values_ = cuts.var_values_;
    cut_number = cuts.cut_number;
    cut_weight_ = cuts.cut_weight_;
    cut_channel_ = cuts.cut_channel_;
    cuts_ = std::move(cuts.cuts_);
    report_ = std::move(cuts.report_);
  }
  CorrectionCuts &operator=(CorrectionCuts &&cuts)  noexcept {
    n_channels_ = cuts.n_channels_;
    var_values_ = cuts.var_values_;
    cut_number = cuts.cut_number;
    cut_weight_ = cuts.cut_weight_;
    cut_channel_ = cuts.cut_channel_;
    cuts_ = std::move(cuts.cuts_);
    report_ = std::move(cuts.report_);
    return *this;
  };
  
  virtual ~CorrectionCuts() { delete[] var_values_; }
//  /**
//   * @brief Adds a cut to the manager.
//   * @param cut pointer to the cut.
//   */
  void AddCut(const CorrectionCut::CallBack &callback) {
    cuts_.emplace_back(callback);
  }

  void Initialize(const Qn::InputVariableManager &var) {
    for (auto &cut : cuts_) {
      cut.Initialize(var);
    }
  }

  /**
   * Checks if the current variables pass the cuts
   * Creates entries in the cut report
   * @param i offset of the variable in case it has a length longer than 1
   * @return Returns true if the cut was passed.
   */
  inline bool CheckCuts(std::size_t i) {
    int icut = 1;
    if (cuts_.empty()) return true;
    ++*((cut_weight_).begin() + i);
    bool passed = true;
    for (auto &cut : cuts_) {
      bool ipass = cut.Check(i) && passed;
      if (ipass) {
        ++*cut_weight_.at(i + n_channels_*icut);
      }
      passed = ipass;
      ++icut;
    }
    return passed;
  }

  /**
   * @brief Fills the cut report.
   */
  void FillReport() {
    if (report_) {
      report_->Fill();
      // reset binvalues container before filling next event.
      auto offset = n_channels_*(cuts_.size() + 1);
      for (std::size_t i = 0; i < n_channels_; ++i) {
        for (std::size_t j = 0; j < (cuts_.size() + 1); ++j) {
          var_values_[2*offset + i + n_channels_*j] = 0;
        }
      }
    }
  }

  /**
   * @brief Initializes the cut report histogram
   * @param report_name name of the histogram.
   * @param n_channels number of channels.
   */
  void CreateCutReport(const std::string &report_name, std::size_t n_channels = 1) {
    if (!cuts_.empty()) {
      n_channels_ = n_channels;
      auto offset = n_channels_*(cuts_.size() + 1);
      cut_number = InputVariable(0, offset);
      cut_channel_ = InputVariable(offset, offset);
      cut_weight_ = InputVariable(2*offset, offset);
      var_values_ = new double[3*offset];
      cut_weight_.values_container_ = var_values_;
      cut_number.values_container_ = var_values_;
      cut_channel_.values_container_ = var_values_;
      for (std::size_t i = 0; i < n_channels_; ++i) {
        for (std::size_t j = 0; j < (cuts_.size() + 1); ++j) {
          var_values_[i + n_channels_*j] = j;
          var_values_[offset + i + n_channels_*j] = i;
          var_values_[2*offset + i + n_channels_*j] = 0;
        }
      }
      if (n_channels_==1) {
        std::string name = report_name + "Cut_Report";
        std::string title(";cuts;entries");
        auto nbins = cuts_.size() + 1;
        float low = 0.;
        float high = cuts_.size() + 1;
        auto histo = new TH1F(name.data(), title.data(), nbins, low, high);
        int icut = 2;
        histo->GetXaxis()->SetBinLabel(1, "all");
        for (auto &cut : cuts_) {
          histo->GetXaxis()->SetBinLabel(icut, cut.Name().data());
          ++icut;
        }
        std::array<InputVariable, 2> arr = {{cut_number, cut_weight_}};
        report_ = std::make_unique<QAHisto1DPtr>(arr, histo);
      } else {
        std::string name = report_name + "Cut_Report";
        std::string title(";cuts;channels");
        auto x_nbins = cuts_.size() + 1;
        auto y_nbins = n_channels_;
        float low = 0.;
        float x_high = cuts_.size() + 1;
        float y_high = n_channels_;
        auto histo = new TH2F(name.data(), title.data(), x_nbins, low, x_high, y_nbins, low, y_high);
        histo->GetXaxis()->SetBinLabel(1, "all");
        int icut = 2;
        for (auto &cut : cuts_) {
          histo->GetXaxis()->SetBinLabel(icut, cut.Name().data());
          ++icut;
        }
        std::array<InputVariable, 3> arr = {{cut_number, cut_channel_, cut_weight_}};
        report_ = std::make_unique<QAHisto2DPtr>(arr, histo);
      }
    }
  }

  /**
   * @brief Adds the cut report to the list.
   * Lifetime of the list and the histogram has to be managed by the user.
   * @param list list containing output histograms.
   */
  void AddToList(TList *list) {
    if (report_) report_->AddToList(list);
  }

 private:
  std::size_t n_channels_ = 0; /// number of channels is zero in case of no report
  double *var_values_ = nullptr; /// pointer to the values which are filled to the histogram.
  InputVariable cut_number; /// Variable of saving cut number
  InputVariable cut_weight_; /// Variable saving a weight used for filling the cut histogram
  InputVariable cut_channel_; /// Variable saving the channel number
  std::vector<CorrectionCut> cuts_; /// vector of cuts which are applied
  std::unique_ptr<Impl::QAHistoBase> report_; //!<! histogram of the cut report.
};

namespace CallBacks {

namespace Details {


template <size_t N>
inline auto FindVariables(const std::array<std::string, N> &names, const Qn::InputVariableManager &manager) {
  std::array<InputVariable, N> result;
  std::transform(std::begin(names), std::end(names), std::begin(result),
                 [&manager] (const std::string& name) { return manager.FindVariable(name); });
  return result;
}

inline auto FindVariables(const std::vector<std::string> &names, const Qn::InputVariableManager &manager) {
  std::vector<InputVariable> result(names.size());
  std::transform(std::begin(names), std::end(names), std::begin(result),
                 [&manager] (const std::string& name) {return manager.FindVariable(name); });
  return result;
}

}

template<typename VariableNames, typename Function>
CorrectionCut::CallBack MakeCut(VariableNames names, Function lambda, const std::string &cut_description) {
  return CorrectionCut::CallBack{[names, lambda, cut_description](const Qn::InputVariableManager &vm) {
    auto variables = Details::FindVariables(names, vm);
    return MakeUniqueCut(variables, lambda, cut_description);
  }};
}
}

}

#endif //FLOW_CORRECTIONCUTS_H
