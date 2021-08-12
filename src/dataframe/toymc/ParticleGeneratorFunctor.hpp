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

#ifndef QNTOOLS_TOYMC_PARTICLEGENERATOR_HPP_
#define QNTOOLS_TOYMC_PARTICLEGENERATOR_HPP_

#include <array>
#include <cmath>
#include <random>

#ifndef ROOT_THREAD_POOL_SIZE
# if ROOT_VERSION_CODE >= ROOT_VERSION(6, 22, 0)
#   define ROOT_THREAD_POOL_SIZE (ROOT::GetThreadPoolSize())
# else
#   define ROOT_THREAD_POOL_SIZE (ROOT::GetImplicitMTPoolSize())
# endif
#endif

namespace Qn::ToyMC {
class ParticleGenerator {
  static constexpr double kPi = M_PI;
  static constexpr double k2Pi = 2*kPi;

 public:
  explicit ParticleGenerator(std::vector<double> harmonics, int nphi_slices = 100):
    nphi_slices_(nphi_slices),
        random_engines_(ROOT_THREAD_POOL_SIZE + 1,std::mt19937_64{std::random_device()()}),
        vns_(harmonics) {}

  ROOT::RVec<double> operator()(unsigned int slot, double psi, int multiplicity) {
    auto phis = ROOT::RVec<double>(multiplicity);
    auto phi_distribution = std::piecewise_linear_distribution<>{
        nphi_slices_, 0, 2 * kPi, [&](const double x) { return PhiPdf(slot, x, psi); }};
    for (auto &phi : phis) {
      phi = phi_distribution(random_engines_.at(slot));
    }
    return phis;
  }

 private:
  std::vector<std::mt19937_64> random_engines_;
  std::size_t nphi_slices_;
  std::vector<double> vns_;

  double PhiPdf(double slot, double phi, double psi) {
    double value = 1.;
//    auto vn_fluctuations = std::uniform_real_distribution<>();
    for (unsigned int n = 1; n < vns_.size() + 1; ++n) {
      value += 2 * vns_[n - 1] * std::cos(n * (phi - psi));
    }
    return value;
  }
};

inline ROOT::RVec<double> MakeUnitWeights(ROOT::RVec<double>& phi) {
  return ROOT::RVec<double>(phi.size(), 1.);
}

inline auto MakeUnitWeightsLambda() {
  return [](ROOT::RVec<double>& phi){return ROOT::RVec<double>(phi.size(), 1.);};
}


}
#endif  // QNTOOLS_TOYMC_PARTICLEGENERATOR_HPP_
