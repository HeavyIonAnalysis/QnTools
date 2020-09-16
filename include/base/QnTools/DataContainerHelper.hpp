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

#include "TVirtualPad.h"
#include "TBrowser.h"
#include "TEnv.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "Axis.hpp"
#include "StatCalculate.hpp"
#include "StatCollect.hpp"

namespace Qn {
//Forward declaration of DataContainer
template<typename T, typename AxisType>
class DataContainer;

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
  enum class Errors { Yonly, XandY };

  static TGraphErrors *ToTGraph(const Qn::DataContainer<Statistics, AxisD> &data, Errors x = Errors::Yonly);
  static TGraphErrors *ToTGraph(const Qn::DataContainer<StatCollect, AxisD> &data, Errors x = Errors::Yonly);
  static TGraphErrors *ToTGraph(const Qn::DataContainer<StatCalculate, AxisD> &data, Errors x = Errors::Yonly);

 private:
  friend Qn::DataContainer<StatCollect, AxisD>;
  friend Qn::DataContainer<Statistics, AxisD>;
  friend Qn::DataContainer<StatCalculate, AxisD>;

  static void Browse(Qn::DataContainer<Statistics, AxisD> *data, TBrowser *b);
  static void Browse(Qn::DataContainer<StatCollect, AxisD> *data, TBrowser *b);
  static void Browse(Qn::DataContainer<StatCalculate, AxisD> *data, TBrowser *b);

  template <typename DataContainer>
  static void ProjectandDraw(DataContainer &data, std::string option, const std::string &axis_name) {
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
    } else if (data.axes_.size() > 1 && !axis_name.empty()) {
      auto projected = data.Projection({axis_name});
      auto graph = DataContainerHelper::ToTGraph(projected, err);
      graph->Draw(option.data());
    } else {
      std::cout << "Not drawn because the DataContainer has dimension: "
                << data.dimension_ << std::endl;
      std::cout << "Only DataContainers with dimension <=2 can be drawn."
                << std::endl;
    }
  }

};

using DrawErrors = DataContainerHelper::Errors;

inline TGraphErrors *ToTGraph(const Qn::DataContainer<Statistics, AxisD> &data, Qn::DrawErrors x = DrawErrors::Yonly) {
  return DataContainerHelper::ToTGraph(data, x);
}

inline TGraphErrors *ToTGraph(const Qn::DataContainer<StatCollect, AxisD> &data, Qn::DrawErrors x = DrawErrors::Yonly) {
  return DataContainerHelper::ToTGraph(data, x);
}

inline TGraphErrors *ToTGraph(DataContainer<StatCalculate, AxisD> &data, Qn::DrawErrors x = DrawErrors::Yonly) {
  return DataContainerHelper::ToTGraph(data, x);
}

}// namespace Qn

#endif//FLOW_DATACONTAINERHELPER_H
