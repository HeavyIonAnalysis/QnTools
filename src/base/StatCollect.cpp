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

#include "StatCollect.hpp"

namespace Qn {

StatCollect Merge(const StatCollect& lhs, const StatCollect& rhs) {
  StatCollect merged;
  merged.statistics_ = Merge(lhs.statistics_, rhs.statistics_);
  merged.boot_strap_ = Merge(lhs.boot_strap_, rhs.boot_strap_);
  merged.SelectWeightandErrors(lhs, rhs);
  return merged;
}

} // namespace Qn
