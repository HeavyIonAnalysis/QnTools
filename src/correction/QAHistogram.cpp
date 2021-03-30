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

#include "QAHistogram.hpp"

namespace Qn {

QAHistogram::QAHistogram(std::string name, const AxisD &axis, std::string weight) : name_(std::move(name)),
                                                                                    weight_(std::move(weight)),
                                                                                    type_(Type::kOneDim) {
  axes_.push_back(axis);
}

QAHistogram::QAHistogram(std::string name, std::vector<AxisD> axes, std::string weight) : name_(std::move(name)),
                                                                                          axes_(std::move(axes)),
                                                                                          weight_(std::move(weight)),
                                                                                          type_(Type::kTwoDim) {
}

QAHistogram::QAHistogram(std::string name, std::vector<AxisD> axes, std::string weight, const Qn::AxisD &histoaxis) : name_(std::move(name)),
                                                                                                                      axes_(std::move(axes)),
                                                                                                                      weight_(std::move(weight)),
                                                                                                                      type_(Type::kTwoDimArray),
                                                                                                                      histoaxis_(histoaxis) {
}


/**
 * Copy Constructor Only use before Initialize is called.
 * @param other QAHistogram which is supposed to be copied.
 */
QAHistogram::QAHistogram(const QAHistogram &other) {
  histogram_ = nullptr;
  name_ = other.name_;
  axes_ = other.axes_;
  weight_ = other.weight_;
  type_ = other.type_;
  histoaxis_ = other.histoaxis_;
}
/**
 * Initialized the QA histogram
 * @param var InputVariableManager used to link to the input variables of the qa histogram
 */
void QAHistogram::Initialize(InputVariableManager &var) {
  ROOT::GetROOT()->cd();
  if (type_ == Type::kOneDim) {
    histogram_ = MakeHisto1D(var);
  } else if (type_ == Type::kTwoDim) {
    histogram_ = MakeHisto2D(var);
  } else if (type_ == Type::kTwoDimArray) {
    histogram_ = MakeHisto2DArray(var);
  }
}
/**
 * Helper function to create the QA histograms
 * @param name name of the histogram
 * @param axes axes used for the histogram
 * @param weightname name of the weight
 * @return unique pointer to the histogram
 */
std::unique_ptr<Impl::QAHistoBase> QAHistogram::MakeHisto1D(InputVariableManager &var) {
  auto axis = axes_[0];
  auto hist_name = name_ + "_A1_" + axis.Name() + "_W_" + weight_;
  auto axisname = std::string(";") + axis.Name();
  auto size = static_cast<const int>(axis.size());
  try {
    var.FindVariable(axis.Name());
  } catch (std::out_of_range &) {
    std::cout << "QAHistogram " << name_ << ": Variable " << axis.Name()
              << " not found. Creating new channel variable."
              << std::endl;
    var.CreateChannelVariable(axis.Name(), axis.size());
  }
  std::array<InputVariable, 2>
      arr = {{var.FindVariable(axis.Name()), var.FindVariable(weight_)}};

  /// 1. double* with bin edges is copied with std::memmove
  /// original array remains untouched
  /// 2. size-1 since last element in the array is the up edge of the last bin
  return std::make_unique<QAHisto1DPtr>(arr, new TH1F(hist_name.data(),
                                                      axisname.data(),
                                                      size - 1,
                                                      (Double_t*) axis.GetPtr()));
}

/**
 * Helper function to create the QA histograms
 * @param name name of the histogram
 * @param axes axes used for the histogram
 * @param weightname name of the weight
 * @return unique pointer to the histogram
 */
std::unique_ptr<Impl::QAHistoBase> QAHistogram::MakeHisto2D(InputVariableManager &var) {
  auto hist_name = name_ + "_A1_" + axes_[0].Name() + "_A2_" + axes_[1].Name() + "_W_" + weight_;
  auto axisname = std::string(";") + axes_[0].Name() + std::string(";") + axes_[1].Name();
  auto size_x = static_cast<const int>(axes_[0].size());
  auto size_y = static_cast<const int>(axes_[1].size());
  for (const auto &axis : axes_) {
    try {
      var.FindVariable(axis.Name());
    } catch (std::out_of_range &) {
      std::cout << "QAHistogram " << name_ << ": Variable " << axis.Name()
                << " not found. Creating new channel variable." << std::endl;
      var.CreateChannelVariable(axis.Name(), axis.size());
    }
  }
  auto upper_edge_x = axes_[0].GetLastBinEdge();
  auto lower_edge_x = axes_[0].GetFirstBinEdge();
  auto upper_edge_y = axes_[1].GetLastBinEdge();
  auto lower_edge_y = axes_[1].GetFirstBinEdge();
  std::array<InputVariable, 3>
      arr = {{var.FindVariable(axes_[0].Name()), var.FindVariable(axes_[1].Name()),
              var.FindVariable(weight_)}};
  auto histo = new TH2F(hist_name.data(), axisname.data(),
                        size_x, lower_edge_x, upper_edge_x,
                        size_y, lower_edge_y, upper_edge_y);
  return std::make_unique<QAHisto2DPtr>(arr, histo);
}

/**
 * Helper function to create the QA histograms
 * @param name name of the histogram
 * @param axes axes used for the histogram
 * @param weightname name of the weight
 * @return unique pointer to the histogram
 */
std::unique_ptr<Impl::QAHistoBase> QAHistogram::MakeHisto2DArray(InputVariableManager &var) {
  auto hist_name = name_ + "_A1_" + axes_[0].Name() + "_A2_" + axes_[1].Name() + "_W_" + weight_;
  auto axisname = std::string(";") + axes_[0].Name() + std::string(";") + axes_[1].Name();
  auto size_x = static_cast<const int>(axes_[0].size());
  auto size_y = static_cast<const int>(axes_[1].size());
  for (const auto &axis : axes_) {
    try {
      var.FindVariable(axis.Name());
    } catch (std::out_of_range &) {
      std::cout << "QAHistogram " << name_ << ": Variable " << axis.Name()
                << " not found. Creating new channel variable." << std::endl;
      var.CreateChannelVariable(axis.Name(), axis.size());
    }
  }
  auto upper_edge_x = axes_[0].GetLastBinEdge();
  auto lower_edge_x = axes_[0].GetFirstBinEdge();
  auto upper_edge_y = axes_[1].GetLastBinEdge();
  auto lower_edge_y = axes_[1].GetFirstBinEdge();
  std::array<InputVariable, 3>
      arr = {{var.FindVariable(axes_[0].Name()), var.FindVariable(axes_[1].Name()),
              var.FindVariable(weight_)}};
  auto histo = new TH2F(hist_name.data(),
                        axisname.data(),
                        size_x,
                        lower_edge_x,
                        upper_edge_x,
                        size_y,
                        lower_edge_y,
                        upper_edge_y);
  auto haxis = std::make_unique<Qn::AxisD>(histoaxis_);
  auto haxisvar = var.FindVariable(histoaxis_.Name());
  return std::make_unique<QAHisto2DPtr>(arr, histo, std::move(haxis), haxisvar);
}
}// namespace Qn