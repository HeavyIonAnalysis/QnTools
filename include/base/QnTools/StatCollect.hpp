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

#ifndef FLOW_STATS_H
#define FLOW_STATS_H

#include "BootStrap.hpp"
#include "Statistics.hpp"
#include "Stat.hpp"
#include "TDigest.hpp"

namespace Qn {
class StatCollect : public Stat {
 public:
  virtual ~StatCollect() = default;

  template <typename SAMPLES>
  void Fill(const double value, const double weight,
            SAMPLES &&sample_multiplicities_) {
    statistics_.Fill(value, weight);
    boot_strap_.Fill(value, weight, sample_multiplicities_);
    tdigest_.add(value, weight);
  }

  friend StatCollect Merge(const StatCollect &lhs, const StatCollect &rhs);
  friend StatCollect MergeBins(const StatCollect &lhs, const StatCollect &rhs);

  [[nodiscard]] auto &GetStatistics() const { return statistics_; }
  [[nodiscard]] auto &GetBootStrap() const { return boot_strap_; }
  [[nodiscard]] auto &GetTDigest() const {return tdigest_; }

  void SetNumberOfSamples(int n_samples) { boot_strap_.SetNumberOfSamples(n_samples); }

 private:
  Statistics statistics_;
  BootStrap boot_strap_;
  Qn::TDigest tdigest_;

  /// \cond CLASSIMP
 ClassDef(StatCollect, 3);
  /// \endcond
};

StatCollect Merge(const StatCollect &lhs, const StatCollect &rhs);

}// namespace Qn

#endif//FLOW_STATS_H
