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

#ifndef QNTOOLS_TOYMC_TRACKINGDETECTORFUNCTOR_HPP_
#define QNTOOLS_TOYMC_TRACKINGDETECTORFUNCTOR_HPP_

#include <vector>
#include <random>
#include "ROOT/RVec.hxx"

#include "Axis.hpp"

namespace Qn::ToyMC {

class ChannelDetectorFunctor {
  static constexpr double kPi = M_PI;
  static constexpr double k2Pi = 2*kPi;
 public:
  ChannelDetectorFunctor(Qn::AxisD phi_bins, std::vector<double> efficiency) :
    phi_bins_(phi_bins),
    efficiency_(efficiency) {}

  ROOT::RVec<double> operator()(ROOT::RVec<double> &phis) {
    ROOT::RVec<double> weights(phi_bins_.GetNBins());
    ROOT::RVec<double> multiplicities(phi_bins_.GetNBins());
    for (const double &phi : phis) {
      auto bin = phi_bins_.FindBin(phi);
      multiplicities.at(bin) += 1.;
    }
    for (auto ibin = 0u; ibin < phi_bins_.GetNBins(); ++ibin) {
      weights.at(ibin) = efficiency_.at(ibin) * multiplicities.at(ibin);
    }
    return weights;
  }

  auto GetPhiBins() {
    ROOT::RVec<double> phis(phi_bins_.GetNBins());
    for (auto iphi = 0u; iphi < phi_bins_.GetNBins(); ++iphi) {
      phis[iphi] = phi_bins_.GetBinCenter(iphi);
    }
    return [phis](){return phis;};
  }

 private:
  Qn::AxisD phi_bins_;
  std::vector<double> efficiency_;
};

}
#endif  // QNTOOLS_CHANNELDETECTORFUNCTOR_HPP_
