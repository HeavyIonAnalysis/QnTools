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

#ifndef QNTOOLS_SRC_DATAFRAME_TOYMC_TOYMCTRACKINGDETECTORFUNCTOR_HPP_
#define QNTOOLS_SRC_DATAFRAME_TOYMC_TOYMCTRACKINGDETECTORFUNCTOR_HPP_

#include <vector>
#include <random>
#include "ROOT/RVec.hxx"

#include "Axis.hpp"

#ifndef ROOT_THREAD_POOL_SIZE
# if ROOT_VERSION_CODE >= ROOT_VERSION(6, 22, 0)
#   define ROOT_THREAD_POOL_SIZE (ROOT::GetThreadPoolSize())
# else
#   define ROOT_THREAD_POOL_SIZE (ROOT::GetImplicitMTPoolSize())
# endif
#endif

namespace Qn::ToyMC {

class TrackingDetectorFunctor {
  static constexpr double kPi = M_PI;
  static constexpr double k2Pi = 2*kPi;
 public:
  explicit TrackingDetectorFunctor(std::function<double(double)> efficiency_function) :
  random_engines_(ROOT_THREAD_POOL_SIZE+1, std::mt19937_64{std::random_device()()}),
  efficiency_function_(efficiency_function) {}

  ROOT::RVec<std::size_t> operator()(unsigned int slot, ROOT::RVec<double> &phis) {
    ROOT::RVec<std::size_t> is_reconstructed;
    auto uniform = std::uniform_real_distribution<>{0.,1.};
    for (auto iphi = 0u; iphi < phis.size(); ++iphi) {
     auto detection_threshold = efficiency_function_(phis[iphi]);
     auto detection = uniform(random_engines_.at(slot));
     if (detection < detection_threshold) {
       is_reconstructed.emplace_back(iphi);
     }
    }
    return is_reconstructed;
 }

 private:
  std::function<double(double)> efficiency_function_;
  std::vector<std::mt19937_64> random_engines_;
};

}
#endif  // QNTOOLS_SRC_DATAFRAME_TOYMC_TOYMCTRACKINGDETECTORFUNCTOR_HPP_
