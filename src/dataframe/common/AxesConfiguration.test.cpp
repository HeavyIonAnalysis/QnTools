// Flow Vector Correction Framework
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

#include <gtest/gtest.h>

#include "AxesConfiguration.hpp"
#include "Axis.hpp"

TEST(AxesConfiguration, GetLowerBinEdges) {
  auto axes =
      Qn::MakeAxes(Qn::AxisD("test1", 2, 0, 2), Qn::AxisD("test2", 2, 0, 2));
  auto nbins = axes.GetBinEdgesIndexMap();
  auto low = nbins[0].first;
  auto high = nbins[0].first;
}