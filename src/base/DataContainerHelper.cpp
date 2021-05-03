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

#include "DataContainerHelper.hpp"

#include <iostream>
#include <iterator>

#include "DataContainer.hpp"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TText.h"
#include "TH2.h"

namespace Qn {

TGraphErrors *DataContainerHelper::ToTGraph(const Qn::DataContainer<Statistics, AxisD> &data, DrawErrors drawerrors) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout
        << "Cannot draw as Graph. Use Projection() to make it one dimensional."
        << std::endl;
    return nullptr;
  }
  auto graph = new TGraphErrors();
  unsigned int ibin = 0;
  for (const auto &bin : data) {
    if (bin.SumWeights() <= 0) {
      ibin++;
      continue;
    }
    auto y = bin.Mean();
    auto ey = bin.StandardErrorOfMean();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(ibin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(ibin);
    auto xhalfwidth = (xhi - xlo)/2.;
    auto x = xlo + xhalfwidth;
    double ex = 0;
    if (drawerrors == Errors::XandY) { ex = xhalfwidth; }
    graph->SetPoint(graph->GetN(), x, y); /* append point to TGraph */
    graph->SetPointError(graph->GetN() - 1, ex, ey); /* edit last point */
    graph->SetMarkerStyle(kFullCircle);
    ibin++;
  }
  return graph;
}

TGraphErrors *DataContainerHelper::ToTGraph(const DataContainerStatCollect &data, DrawErrors drawerrors) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout
        << "Cannot draw as Graph. Use Projection() to make it one dimensional."
        << std::endl;
    return nullptr;
  }
  auto graph = new TGraphErrors();
  unsigned int ibin = 0;
  for (const auto &bin : data) {
    const auto &stats = bin.GetStatistics();
    if (stats.SumWeights() <= 0) {
      ibin++;
      continue;
    }
    auto y = stats.Mean();
    auto ey = stats.StandardErrorOfMean();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(ibin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(ibin);
    auto xhalfwidth = (xhi - xlo)/2.;
    auto x = xlo + xhalfwidth;
    double ex = 0;
    if (drawerrors == Errors::XandY) { ex = xhalfwidth; }
    graph->SetPoint(graph->GetN(), x, y);
    graph->SetPointError(graph->GetN() - 1, ex, ey);
    graph->SetMarkerStyle(kFullCircle);
    ibin++;
  }
  return graph;
}

TGraphErrors *DataContainerHelper::ToTGraph(const DataContainer<StatCalculate, AxisD> &data, DrawErrors drawerror) {
  if (data.GetAxes().size() > 1) {
    std::cout << "Data container has more than one dimension. " << std::endl;
    std::cout << "Cannot draw as Graph. Use Projection() to make it one dimensional." << std::endl;
    return nullptr;
  }
  auto graph = new TGraphErrors();
  unsigned int ibin = 0;
  for (const auto &bin : data) {
    if (bin.SumWeights() <= 0.) {
      ibin++;
      continue;
    }
    auto y = bin.Mean();
    auto ey = bin.StandardErrorOfMean();
    auto xhi = data.GetAxes().front().GetUpperBinEdge(ibin);
    auto xlo = data.GetAxes().front().GetLowerBinEdge(ibin);
    auto xhalfwidth = (xhi - xlo)/2.;
    auto x = xlo + xhalfwidth;
    double ex = 0;
    if (drawerror == Errors::XandY) { ex = xhalfwidth; }
    graph->SetPoint(graph->GetN(), x, y);
    graph->SetPointError(graph->GetN() - 1, ex, ey);
    graph->SetMarkerStyle(kFullCircle);
    ibin++;
  }
  return graph;
}

void DataContainerHelper::Browse(DataContainerStatistic *data, TBrowser *b) {
  using DrawErrorGraph = Internal::ProjectionDrawable<TGraphErrors *>;
  if (!data->list_) data->list_ = new TList();
  data->list_->SetOwner(true);
  for (auto &axis : data->axes_) {
    TGraphErrors *graph;
    if (data->dimension_ > 1) {
      graph = DataContainerHelper::ToTGraph(data->Projection({axis.Name()}),
                                            Errors::Yonly);
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
    for (auto iaxis = std::begin(data->axes_); iaxis < std::end(data->axes_);
         ++iaxis) {
      for (auto jaxis = std::begin(data->axes_); jaxis < std::end(data->axes_);
           ++jaxis) {
        if (iaxis != jaxis) {
          auto proj = data->Projection({iaxis->Name(), jaxis->Name()});
          auto axes = proj.GetAxes();
          auto axis1 = axes.at(0);
          auto axis2 = axes.at(1);
          auto name = (jaxis->Name() + ":" + iaxis->Name());
          auto histo2d = new TH2F(
              name.c_str(), name.c_str(), axis1.size(), axis1.GetFirstBinEdge(),
              axis1.GetLastBinEdge(), axis2.size(), axis2.GetFirstBinEdge(),
              axis2.GetLastBinEdge());
          for (int i = 0; i < proj.size(); ++i) {
            auto value = proj.At(i).Mean();
            auto bins = proj.GetIndex(i);
            histo2d->SetBinContent(bins[0]+1, bins[1]+1, value);
          }
          list2d->Add(histo2d);
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

void DataContainerHelper::Browse(DataContainerStatCalculate *data, TBrowser *b) {
  using DrawErrorGraph = Internal::ProjectionDrawable<TGraphErrors *>;
  if (!data->list_) data->list_ = new TList();
  data->list_->SetOwner(true);
  for (auto &axis : data->axes_) {
    TGraphErrors *graph;
    if (data->dimension_ > 1) {
      graph = DataContainerHelper::ToTGraph(data->Projection({axis.Name()}),
                                            Errors::Yonly);
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

void DataContainerHelper::Browse(DataContainerStatCollect *data, TBrowser *b) {
  using DrawErrorGraph = Internal::ProjectionDrawable<TGraphErrors *>;
  if (!data->list_) data->list_ = new TList();
  data->list_->SetOwner(true);
  for (auto &axis : data->axes_) {
    TGraphErrors *graph;
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

}// namespace Qn
