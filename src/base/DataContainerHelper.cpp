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

#include <iostream>
#include <iterator>

#include "DataContainer.hpp"
#include "DataContainerHelper.hpp"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TText.h"

namespace Qn {

TGraphAsymmErrors *DataContainerHelper::ToTGraphShifted(const DataContainerStatistic &data,
                                                        int i,
                                                        int maxi, Errors drawerrors) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout << "Cannot draw as Graph. Use Projection() to make it one dimensional." << std::endl;
    return nullptr;
  }
  auto graph = new TGraphAsymmErrors();
  unsigned int ibin = 0;
  for (const auto &bin : data) {
    auto tbin = bin;
    if (tbin.Entries() == 0) continue;
    auto y = tbin.Mean();
    auto ylo = tbin.MeanError();
    auto yhi = tbin.MeanError();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(ibin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(ibin);
    auto x = xlo + ((xhi - xlo) * static_cast<double>(i) / maxi);
    double exl = 0;
    double exh = 0;
    if (drawerrors == Errors::XandY) {
      exl = x - xlo;
      exh = xhi - x;
    }
    graph->SetPoint(ibin, x, y);
    graph->SetPointError(ibin, exl, exh, ylo, yhi);
    graph->SetMarkerStyle(kFullCircle);
    ibin++;
  }
  return graph;
}

TGraphAsymmErrors *DataContainerHelper::ToTGraphShifted(const DataContainerStats &data,
                                                        int i,
                                                        int maxi, Errors drawerrors) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout << "Cannot draw as Graph. Use Projection() to make it one dimensional." << std::endl;
    return nullptr;
  }
  auto graph = new TGraphAsymmErrors();
  unsigned int ibin = 0;
  for (const auto &bin : data) {
    auto tbin = bin;
    if (tbin.N() == 0 && tbin.GetState() != Stats::State::MEAN_ERROR) continue;
    if (tbin.GetState() != Stats::State::MEAN_ERROR) tbin.CalculateMeanAndError();
    auto y = tbin.Mean();
    auto ylo = tbin.LowerMeanError();
    auto yhi = tbin.UpperMeanError();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(ibin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(ibin);
    auto x = xlo + ((xhi - xlo) * static_cast<double>(i) / maxi);
    double exl = 0;
    double exh = 0;
    if (drawerrors == Errors::XandY) {
      exl = x - xlo;
      exh = xhi - x;
    }
    if (!std::isnan(y) && !std::isnan(ylo) && !std::isnan(yhi)) {
      graph->SetPoint(ibin, x, y);
      graph->SetPointError(ibin, exl, exh, ylo, yhi);
      graph->SetMarkerStyle(kFullCircle);
      ibin++;
    }
  }
  return graph;
}

TGraphAsymmErrors *DataContainerHelper::ToTGraph(const DataContainerStats &data,
                                                 Errors drawerrors) {
  return ToTGraphShifted(data, 1, 2, drawerrors);
}

TGraphAsymmErrors *DataContainerHelper::ToTGraph(const DataContainerStatistic &data,
                                                 Errors drawerrors) {
  return ToTGraphShifted(data, 1, 2, drawerrors);
}

/**
 * Create a TMultiGraph containing profile graphs for each bin of the specified axis.
 * @param data datacontainer to be graphed.
 * @param axisname name of the axis which is used for the multigraph.
 * @return graph containing a profile graph for each bin.
 */
TMultiGraph *DataContainerHelper::ToTMultiGraph(const DataContainerStats &data,
                                                const std::string &axisname,
                                                Errors drawerrors) {
  auto multigraph = new TMultiGraph();
  if (data.GetAxes().size() != 2) {
    std::cout << "Data Container dimension has wrong dimension " << data.GetAxes().size() << "\n";
    return nullptr;
  }
  AxisD axis;
  try {
    axis = data.GetAxis(axisname);
  } catch (std::exception &) {
    throw std::logic_error("axis not found");
  }
  for (unsigned int ibin = 0; ibin < axis.size(); ++ibin) {
    auto subdata = data.Select({axisname, {axis.GetLowerBinEdge(ibin), axis.GetUpperBinEdge(ibin)}});
    auto subgraph = ToTGraphShifted(subdata, ibin, axis.size(), drawerrors);
    subgraph->SetTitle(Form("%.2f - %.2f", axis.GetLowerBinEdge(ibin), axis.GetUpperBinEdge(ibin)));
    subgraph->SetMarkerStyle(kFullCircle);
    multigraph->Add(subgraph);
  }
  return multigraph;
}

void DataContainerHelper::StatisticBrowse(DataContainerStatistic *data, TBrowser *b) {
  using DrawErrorGraph = Internal::ProjectionDrawable<TGraphAsymmErrors *>;
  if (!data->list_) data->list_ = new TList();
  data->list_->SetOwner(true);
  for (auto &axis : data->axes_) {
    TGraphAsymmErrors *graph;
    if (data->dimension_ > 1) {
      graph = DataContainerHelper::ToTGraph(data->Projection({axis.Name()}), Errors::Yonly);
    } else {
      graph = DataContainerHelper::ToTGraph(*data, Errors::Yonly);
    }
    graph->SetName(axis.Name().data());
    graph->SetTitle(axis.Name().data());
    graph->GetXaxis()->SetTitle(axis.Name().data());
    auto *drawable = new DrawErrorGraph(graph);
    data->list_->Add(drawable);
  }
  for (int i = 0; i < data->list_->GetSize(); ++i) {
    b->Add(data->list_->At(i));
  }
}

void DataContainerHelper::StatsBrowse(DataContainerStats *data, TBrowser *b) {
  using DrawErrorGraph = Internal::ProjectionDrawable<TGraphAsymmErrors *>;
  using DrawMultiGraph = Internal::ProjectionDrawable<TMultiGraph *>;
  if (!data->list_) data->list_ = new TList();
  data->list_->SetOwner(true);
  for (auto &axis : data->axes_) {
    TGraphAsymmErrors *graph;
    if (data->dimension_ > 1) {
      graph = DataContainerHelper::ToTGraph(data->Projection({axis.Name()}), Errors::Yonly);
    } else {
      graph = DataContainerHelper::ToTGraph(*data, Errors::Yonly);
    }
    graph->SetName(axis.Name().data());
    graph->SetTitle(axis.Name().data());
    graph->GetXaxis()->SetTitle(axis.Name().data());
    auto *drawable = new DrawErrorGraph(graph);
    data->list_->Add(drawable);
  }
  if (data->dimension_ > 1) {
    auto list2d = new TList();
    for (auto iaxis = std::begin(data->axes_); iaxis < std::end(data->axes_); ++iaxis) {
      for (auto jaxis = std::begin(data->axes_); jaxis < std::end(data->axes_); ++jaxis) {
        if (iaxis != jaxis) {
          auto proj = data->Projection({iaxis->Name(), jaxis->Name()});
          auto mgraph = DataContainerHelper::ToTMultiGraph(proj, iaxis->Name(), Errors::Yonly);
          auto name = (jaxis->Name() + ":" + iaxis->Name());
          mgraph->SetName(name.data());
          mgraph->SetTitle(name.data());
          mgraph->GetXaxis()->SetTitle(jaxis->Name().data());
          mgraph->GetYaxis()->SetTitle("Correlation");
          auto *drawable = new DrawMultiGraph(mgraph);
          list2d->Add(drawable);
        }
      }
    }
    list2d->SetName("2D");
    list2d->SetOwner(true);
    data->list_->Add(list2d);
  }
  for (int i = 0; i < data->list_->GetSize(); ++i) {
    b->Add(data->list_->At(i));
  }
}

void DataContainerHelper::NDraw(DataContainerStats &data,
                                std::string option,
                                const std::string &axis_name = {}) {
  option = (option.empty()) ? std::string("ALP PMC PLC Z") : option;
  Errors err = Errors::Yonly;
  auto foundoption = option.find("XErrors");
  if (foundoption != std::string::npos) {
    err = Errors::XandY;
    option.erase(foundoption, std::string("XErrors").size());
  }
  if (data.axes_.size() == 1) {
    auto graph = DataContainerHelper::ToTGraph(data, err);
    graph->Draw(option.data());
  } else if (data.axes_.size() == 2) {
    try {
      auto mgraph = DataContainerHelper::ToTMultiGraph(data, axis_name, err);
      mgraph->Draw(option.data());
    } catch (std::exception &) {
      std::string axes_in_dc = "Axis \"" + axis_name + "\" not found.\n";
      axes_in_dc += "Available axes:";
      axes_in_dc += " \"" + data.axes_.at(0).Name() + "\"";
      for (size_t i = 1; i < data.axes_.size(); ++i) {
        axes_in_dc += ", \"" + data.axes_.at(i).Name() + "\"";
      }
      std::cout << axes_in_dc << std::endl;
    }
  } else {
    std::cout << "Not drawn because the DataContainer has dimension: " << data.dimension_ << std::endl;
    std::cout << "Only DataContainers with dimension <=2 can be drawn." << std::endl;
  }
}

TGraph *DataContainerHelper::ToBootstrapScatterGraph(const Qn::DataContainer<Stats, AxisD> &data) {
  auto scatter_graph = new TGraph;
  std::size_t i_bin = 0;
  for (const auto &bin : data) {
    if (bin.N() == 0 && bin.GetState() != Stats::State::MEAN_ERROR) continue;
    auto tbin = bin;
    if (tbin.GetState() != Stats::State::MEAN_ERROR) tbin.CalculateMeanAndError();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(i_bin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(i_bin);
    auto offset = (xhi + xlo) / 2;
    tbin.GetReSamples().ScatterGraph(*scatter_graph, offset, (xhi - xlo) / 2);
    ++i_bin;
  }
  return scatter_graph;
}

TGraph *DataContainerHelper::ToErrorComparisonGraph(const Qn::DataContainer<Stats, AxisD> &data) {
  auto error_ratio = new TGraph;
  std::size_t i_bin = 0;
  for (const auto &bin : data) {
    if (bin.N() == 0 && bin.GetState() != Stats::State::MEAN_ERROR) continue;
    auto tbin = bin;
    if (tbin.GetState() != Stats::State::MEAN_ERROR) tbin.CalculateMeanAndError();
    auto ratio = tbin.RatioOfErrors();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(i_bin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(i_bin);
    auto xmean = (xhi + xlo) / 2;
    error_ratio->SetPoint(error_ratio->GetN(), xmean, ratio);
    ++i_bin;
  }
  return error_ratio;
}

TCanvas *UncertaintyComparison(const Qn::DataContainer<Stats, AxisD> &data) {
  auto axes = data.GetAxes();
  auto size = axes.size();
  auto canvas =
      new TCanvas("UncertaintyComparison", "ratio of bootstrap uncertainty to statistical uncertainty", 600, 150 * size);
  std::size_t i_axis = 1;
  float stepsize = 0.9 / size;
  auto label = new TText(0.25, 0.925, "ratio of bootstrap to statistical uncertainty");
  label->SetNDC(true);
  label->SetTextColor(kGray + 1);
  label->SetTextFont(43);
  label->SetTextSize(16);
  label->Draw();
  for (auto &axis : axes) {
    canvas->cd();
    auto pad = new TPad(axis.Name().data(), "", 0.0, stepsize * (i_axis - 1), 1.0, stepsize * (i_axis));
    pad->Draw();
    pad->cd();
    TGraph *graph = nullptr;
    if (axes.size() > 1) {
      graph = ToErrorComparisonGraph(data.Projection({axis.Name()}));
    } else {
      graph = ToErrorComparisonGraph(data);
    }
    auto y_max = TMath::MaxElement(graph->GetN(), graph->GetY());
    auto y_min = TMath::MinElement(graph->GetN(), graph->GetY());
    auto y_range = y_max - y_min;
    auto c_y_max = y_max + y_range * 0.3;
    double c_y_min;
    if (y_min > 1.0) {
      c_y_min = 1.0 - y_range * 0.4;
      y_min = 1.0;
    } else {
      c_y_min = y_min - y_range * 0.4;
    }

    auto x_max = axis.GetLastBinEdge();
    auto x_min = axis.GetFirstBinEdge();
    auto x_range = x_max - x_min;
    x_max = x_max + x_range * 0.1;
    x_min = x_min - x_range * 0.1;
    pad->Range(x_min, c_y_min, x_max, c_y_max);
    std::string axis_name = std::string(" ;") + axis.Name();

    auto x_axis = new TGaxis(axis.GetFirstBinEdge(),
                             1.0, axis.GetLastBinEdge(),
                             1.0,
                             axis.GetFirstBinEdge(),
                             axis.GetLastBinEdge(),
                             5,
                             "-+");
    x_axis->SetLabelFont(43);
    x_axis->SetLabelSize(16);
    x_axis->SetTitle(axis.Name().data());
    x_axis->SetTitleSize(16);
    x_axis->SetTitleFont(43);
    x_axis->SetLineColor(kGray + 1);
    x_axis->SetTitleColor(kGray + 1);
    x_axis->SetLabelColor(kGray + 1);
    x_axis->SetTitleOffset(4.0);
    x_axis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    x_axis->CenterTitle();
    x_axis->Draw();

    auto y_axis = new TGaxis(axis.GetFirstBinEdge(), y_min, axis.GetFirstBinEdge(), y_max, y_min, y_max, 5, "+S");
    y_axis->SetLabelFont(43);
    y_axis->SetLabelSize(16);
    //    y_axis->SetTitle("ratio");
    y_axis->SetTitleSize(16);
    y_axis->SetTitleFont(43);
    y_axis->SetTickSize(0.01);
    y_axis->SetLineColor(kGray + 1);
    y_axis->SetTitleColor(kGray + 1);
    y_axis->SetLabelColor(kGray + 1);
    y_axis->SetLabelOffset(-0.01);
    y_axis->SetTitleOffset(-0.8);
    y_axis->CenterTitle();
    y_axis->Draw();

    graph->SetLineWidth(2);
    graph->Draw("L");
    ++i_axis;
  }
  canvas->SaveAs("bserrorscomparison.pdf");
  return canvas;
}

}// namespace Qn