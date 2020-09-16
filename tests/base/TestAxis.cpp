//
// Created by eugene on 07/09/2020.
//

#include <gtest/gtest.h>

#include "Axis.hpp"

using namespace Qn;

namespace {

TEST(Axis,FindBin) {
  AxisD axis("Centrality", 10, 0, 10);

  EXPECT_EQ(axis.FindBin(0.5), 0);
  EXPECT_EQ(axis.FindBin(0), 0);
  EXPECT_EQ(axis.FindBin(1), 1);
  EXPECT_EQ(axis.FindBin(-1), -1);
  EXPECT_EQ(axis.FindBin(9.9), 9);
  EXPECT_EQ(axis.FindBin(10), -1); // Is it expected?
  EXPECT_EQ(axis.FindBin(11), -1);
}

}
