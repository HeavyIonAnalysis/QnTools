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

#ifndef FLOW_QAHISTOGRAM_H
#define FLOW_QAHISTOGRAM_H

#include <utility>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TH3F.h"
#include "TList.h"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/RIntegerSequence.hxx"
#include "TROOT.h"

#include "Axis.hpp"
#include "InputVariableManager.hpp"

namespace Qn {
namespace Impl {
/**
 * Base class of a QA histogram
 */
struct QAHistoBase {
  virtual ~QAHistoBase() = default;
  virtual void Fill() = 0;
  virtual void AddToList(TList *) = 0;
  static inline void FillHistogram(std::unique_ptr<QAHistoBase> &histo) { histo->Fill(); }
  virtual std::string GetName() const = 0;

};
/**
 * Wrapper for a ROOT histogram, which allows it to be filled by the correction manager.
 * @tparam HISTO type of histogram.
 * @tparam N number of dimensions
 * @tparam VAR Type of the variable
 */
template<typename HISTO, int N, typename VAR>
class QAHisto : public QAHistoBase {
 public:
  QAHisto(std::array<VAR, N> vec, HISTO histo, std::unique_ptr<Qn::AxisD> axis, const Qn::InputVariable &axisvar) :
      vars_(std::move(vec)),
      axis_(std::move(axis)),
      axisvar_(axisvar) {
    if (axis_) {
      for (std::size_t i = 0; i < axis_->size(); ++i) {
        auto histname = histo->GetName();
        auto binname = std::string(histname) + "_" + axis_->GetBinName(i);
        auto binhisto = dynamic_cast<HISTO>(histo->Clone(binname.data()));
        histo_.push_back(binhisto);
      }
      name_ = histo->GetName();
      delete histo;
    } else {
      histo_.push_back(histo);
    }
  }
  QAHisto(std::array<VAR, N> vec, HISTO histo) :
      vars_(std::move(vec)) {
    histo_.push_back(histo);
  }

  /**
   * Implementation of the fill function
   * @tparam VARS type of array
   * @tparam I index sequence
   * @param vars Array of variables used for filling the histogram
   */
  template<typename VARS, std::size_t... I>
  void FillImpl(const VARS &variables, std::index_sequence<I...>) {
    auto bin = 0;
    if (axis_) {
      bin = axis_->FindBin(*axisvar_.begin());
      if (bin > -1) histo_.at(bin)->FillN(variables[0].size(), (variables[I].begin())...);
    } else {
      histo_.at(bin)->FillN(variables[0].size(), (variables[I].begin())...);
    }
  }
  /**
   * Fill function.
   */
  void Fill() override {
    return FillImpl(vars_, std::make_index_sequence<N>{});
  };
  /**
   * Add the histogram to the list.
   * @param list pointer to the list. Lifetime of the histogram hast to be managed by the list.
   */
  void AddToList(TList *list) override {
    if (axis_) {
      auto dir = new ::TList();
      dir->SetName(name_.data());
      for (auto &histo : histo_) {
        dir->Add(histo);
      }
      list->Add(dir);
    } else {
      for (auto &histo : histo_) {
        list->Add(histo);
      }
    }
  }

  std::string GetName() const override { return name_; }
 private:
  std::array<VAR, N> vars_; /// Array of variables to be filled in the histogram.
  std::vector<HISTO> histo_; /// Histogram (e.g. TH1, TH2) which support the filling with FillN(...).
  std::unique_ptr<Qn::AxisD> axis_ = nullptr; // Creates a histogram for each bin of the axis
  VAR axisvar_;              /// input variable associated with the axis
  std::string name_;         /// name of the QA histogram
};
}
/// specializations used in the framework
using QAHisto1DPtr = Impl::QAHisto<TH1F *, 2, Qn::InputVariable>;
using QAHisto2DPtr = Impl::QAHisto<TH2F *, 3, Qn::InputVariable>;


/**
 * HistogramWrapper for the delayed construction of the histogram.
 */
class QAHistogram {
 public:
  enum class Type {
    kUndefined,
    kOneDim,
    kTwoDim,
    kTwoDimArray
  };
  QAHistogram() = default;
  virtual ~QAHistogram() = default;
  QAHistogram(QAHistogram &&) = default;
  QAHistogram &operator=(QAHistogram &&) = default;
  QAHistogram(const QAHistogram &other);
  QAHistogram(std::string, const AxisD&, std::string);
  QAHistogram(std::string, std::vector<AxisD>, std::string);
  QAHistogram(std::string, std::vector<AxisD>, std::string, const Qn::AxisD&);
  void Initialize(InputVariableManager &var);
  void Fill() { histogram_->Fill(); }
  void AddToList(TList *list) { histogram_->AddToList(list); }
  std::string GetName() const {return name_;}
 private:
  std::unique_ptr<Impl::QAHistoBase> histogram_; //!<! pointer to the QAHisto<>
  std::string name_; // name of the histogram
  std::vector<AxisD> axes_; /// vector of histogram axes.
  std::string weight_; /// name of the weight
  Type type_ = Type::kUndefined; /// style type of histogram
  AxisD histoaxis_; /// additional axis used for array-of-2d-histograms style.

  std::unique_ptr<Impl::QAHistoBase> MakeHisto1D(InputVariableManager &var);
  std::unique_ptr<Impl::QAHistoBase> MakeHisto2D(InputVariableManager &var);
  std::unique_ptr<Impl::QAHistoBase> MakeHisto2DArray(InputVariableManager &var);
  /// \cond CLASSIMP
 ClassDef(QAHistogram, 1);
/// \endcond
};

class QAHistograms {
 public:
  QAHistograms() = default;
  virtual ~QAHistograms() = default;
  template<typename... Args>
  void Add(Args &&... args) {
    histograms_.emplace_back(std::forward<Args>(args)...);
  }
  void AddToList(TList *list) {
    for (auto &histo : histograms_) {
      histo.AddToList(list);
    }
  }
  void Initialize(InputVariableManager &var) {
    for (auto &histo : histograms_) {
      histo.Initialize(var);
    }
  }
  void Fill() {
    for (auto &histo : histograms_) {
      histo.Fill();
    }
  }
 private:
  std::vector<QAHistogram> histograms_;
  /// \cond CLASSIMP
 ClassDef(QAHistograms, 1);
/// \endcond
};
}

#endif //FLOW_QAHISTOGRAM_H
