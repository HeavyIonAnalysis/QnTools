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

#include "QVector.hpp"

TEST(QVectorUnitTest, Constructor) {
  std::bitset<Qn::QVector::kmaxharmonics> bits("11110000", 8, '0', '1');
  Qn::QVector q_vector(bits,Qn::QVector::PLAIN, Qn::QVector::Normalization::M);
  EXPECT_EQ(q_vector.GetNorm(), Qn::QVector::Normalization::M);
  EXPECT_EQ(q_vector.GetCorrectionStep(), Qn::QVector::PLAIN);
  EXPECT_EQ(q_vector.GetHarmonics(), bits);
}