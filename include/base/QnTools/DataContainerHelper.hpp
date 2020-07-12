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

#ifndef FLOW_DATACONTAINERHELPER_H
#define FLOW_DATACONTAINERHELPER_H

#include "TBrowser.h"
#include "TEnv.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "Axis.hpp"
#include "ReSamples.hpp"
#include "Stats.hpp"

namespace Qn {
//Forward declaration of DataContainer
template<typename T, typename AxisType>
class DataContainer;

inline float MergeBins(const float &a, const float &b) { return a + b; }

namespace Internal {

/**
 * @brief       Template class to facilitate the drawing of Projections of DataContainer.
 * @attention   Takes ownership of graph and will delete it at the end of the lifetime.
 * @tparam T    Type of Graph to draw.
 */
template<typename T>
struct ProjectionDrawable : public TNamed {
  ProjectionDrawable() = default;
  ~ProjectionDrawable() override { delete graph; }
  explicit ProjectionDrawable(T projection) : graph(std::move(projection)) {
    SetName(projection->GetName());
    SetTitle(projection->GetTitle());
  }
  void Browse(TBrowser *b) override {
    TString opt = gEnv->GetValue("TGraph.BrowseOption", "");
    if (opt.IsNull()) {
      opt = b ? b->GetDrawOption() : "AlP PLC PMC Z";
      opt = (opt == "") ? "ALP PMC PLC Z" : opt.Data();
    }
    graph->Draw(opt);
    gPad->Update();
  }
  T graph = nullptr;
};
}// namespace Internal

/**
 * Helper class for DataContainer used for visualization.
 */
class DataContainerHelper {
 public:
  enum class Errors { Yonly,
                      XandY };

  static TGraph *ToErrorComparisonGraph(const Qn::DataContainer<Stats, AxisD> &data);
  static TGraph *ToBootstrapScatterGraph(const Qn::DataContainer<Stats, AxisD> &data);
  static TGraphAsymmErrors *ToTGraph(const Qn::DataContainer<Stats, AxisD> &data, Errors x = Errors::Yonly);
  static TGraphAsymmErrors *ToTGraph(const Qn::DataContainer<Statistic, AxisD> &data, Errors x = Errors::Yonly);
  static TGraphAsymmErrors *ToTGraphShifted(const Qn::DataContainer<Stats, AxisD> &data,
                                            int i,
                                            int max,
                                            Errors x = Errors::Yonly);
  static TGraphAsymmErrors *ToTGraphShifted(const Qn::DataContainer<Statistic, AxisD> &data,
                                            int i,
                                            int max,
                                            Errors x = Errors::Yonly);
  static TMultiGraph *ToTMultiGraph(const Qn::DataContainer<Stats, AxisD> &data,
                                    const std::string &axisname,
                                    Errors x = Errors::Yonly);

 private:
  friend Qn::DataContainer<Stats, AxisD>;
  friend Qn::DataContainer<Statistic, AxisD>;
  static void StatsBrowse(Qn::DataContainer<Stats, AxisD> *data, TBrowser *b);
  static void StatisticBrowse(Qn::DataContainer<Statistic, AxisD> *data, TBrowser *b);
  static void NDraw(Qn::DataContainer<Stats, AxisD> &data, std::string option, const std::string &axis_name);
};

using Errors = DataContainerHelper::Errors;
constexpr auto ToErrorComparisonGraph = &DataContainerHelper::ToErrorComparisonGraph;
constexpr auto ToBootstrapScatterGraph = &DataContainerHelper::ToBootstrapScatterGraph;
constexpr auto ToTMultiGraph = &DataContainerHelper::ToTMultiGraph;

inline TGraphAsymmErrors *ToTGraph(const Qn::DataContainer<Statistic, AxisD> &data, Qn::Errors x = Errors::Yonly) {
  return DataContainerHelper::ToTGraph(data, x);
}

inline TGraphAsymmErrors *ToTGraph(const Qn::DataContainer<Stats, AxisD> &data, Qn::Errors x = Errors::Yonly) {
  return DataContainerHelper::ToTGraph(data, x);
}

TCanvas *UncertaintyComparison(const Qn::DataContainer<Stats, AxisD> &data);

}// namespace Qn

#endif//FLOW_DATACONTAINERHELPER_H
