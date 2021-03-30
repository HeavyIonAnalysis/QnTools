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
#ifndef QNTOOLS_RESAMPLEHELPER_H_
#define QNTOOLS_RESAMPLEHELPER_H_

#include <algorithm>
#include <random>
#include <vector>

#include <ROOT/RVec.hxx>
#include <RtypesCore.h>

namespace Qn::Correlation {

/**
 * Class for creating samples in the DataFrame using
 */
class ReSampleFunctor {
 public:
  /**
   * Constructor
   * @param df Dataframe to which the samples are added.
   * @param n number of samples which are to be used for error estimation
   */
  explicit ReSampleFunctor(std::size_t n = 100)
      : n_(n), generator_(std::random_device{}()), poisson_(1) {}

  /**
   * Defines the tiems a single event enters in the different samples.
   * This function is called for each event during the event loop.
   * @return returns the vector of multiplicity for the different samples.
   */
  ROOT::RVec<ULong64_t> operator()() {
    ROOT::RVec<ULong64_t> vec(n_);
    for (auto &entry : vec) entry = poisson_(generator_);
    return vec;
  }

 private:
  const std::size_t n_;                  /// Number of samples
  std::mt19937 generator_;               /// Random number generator
  std::poisson_distribution<> poisson_;  /// distribution of events per sample.
};

/**
 * Helper function to add the samples to the RDataFrame.
 * The Event loop is lazily executed as soon as the "samples" branch information
 * is consumed by a correlation.
 * @tparam DataFrame type of the RDataFrame.
 * @param df RDataFrame wrapping the input data.
 * @param n Number of samples.
 * @return The resulting RDataFrame with the sample definition.
 */
template <typename DataFrame>
auto Resample(DataFrame df, std::size_t n) {
  return df.template Define("samples", ReSampleFunctor(n), {});
}

}  // namespace Qn::Correlation
#endif  // RESAMPLEHELPER_H_
