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

#include "AverageHelper.hpp"
#include "AxesConfiguration.hpp"
#include "QVectorHelper.hpp"
#include "RecenterAction.hpp"

/**
 * Testing if exceptions are thrown in case of double initialization.
 */
TEST(AverageHelper, OnlyOneInitialization) {
  Qn::DataContainerQVector qvec_proto;
  auto event_axis = Qn::AxisD("event", 1, 0, 2);
  auto axes = Qn::MakeAxes(event_axis);
  auto res = Qn::Correction::MakeRecenterAction("test", axes, "q");
  TTreeReader reader;
  try {
    auto rec = Qn::AverageHelper(res)
                   .SetInitializationWithInitializationObject(&qvec_proto)
                   .SetInitializationWithExternalTTreeReader(&reader);
  } catch (std::runtime_error &e) {
    EXPECT_EQ(e.what(),
              std::string("Error! SetInitializationWithInitializationObject "
                          "has already been called."));
  }
  try {
    auto rec = Qn::AverageHelper(res)
                   .SetInitializationWithExternalTTreeReader(&reader)
                   .SetInitializationWithInitializationObject(&qvec_proto);
  } catch (std::runtime_error &e) {
    EXPECT_EQ(e.what(),
              std::string("Error! SetInitializationWithExternalTTreeReader has "
                          "already been called."));
  }
}