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
#ifndef CHANNELDETECTOR_H_
#define CHANNELDETECTOR_H_

#include <algorithm>
#include <functional>
#include <utility>

#include "Axis.hpp"

class ChannelDetector {
 public:
  ChannelDetector(const Qn::AxisD &axis, std::function<bool(double)> cuts) : n_channels_(axis.size()),
                                                                             axis_(axis),
                                                                             multiplicity_rec_(axis.size()),
                                                                             multiplicity_true_(axis.size()),
                                                                             efficiencies_(axis.size()),
                                                                             cuts_(std::move(cuts)) {
    for (std::size_t i = 0; i < n_channels_; ++i) {
      auto bin_center = axis_.GetLowerBinEdge(i) + (axis_.GetUpperBinEdge(i) - axis_.GetLowerBinEdge(i)) / 2;
      phi_.push_back(bin_center);
      efficiencies_.at(i) = 1.;
    }
  }

  void Reset() {
    for (auto &mult : multiplicity_rec_) mult = 0.;
    for (auto &mult : multiplicity_true_) mult = 0.;
  }

  void Detect(double phi) {
    auto bin = axis_.FindBin(phi);
    if (cuts_(phi)) {
      multiplicity_rec_[bin] = multiplicity_rec_[bin] + 1;
    }
    multiplicity_true_[bin] = multiplicity_true_[bin] + 1;
  }

  void FillDataRec(double *values_array, std::size_t position_phi, std::size_t position_weights) {
    for (std::size_t i = 0; i < n_channels_; ++i) {
      values_array[position_phi + i] = phi_.at(i);
      values_array[position_weights + i] = multiplicity_rec_.at(i) * efficiencies_.at(i);
    }
  }

  void FillDataTruth(double *values_array, std::size_t position_phi, std::size_t position_weights) {
    for (std::size_t i = 0; i < n_channels_; ++i) {
      values_array[position_phi + i] = phi_.at(i);
      values_array[position_weights + i] = multiplicity_true_.at(i);
    }
  }

  void SetChannelEfficencies(std::vector<double> efficiencies) {
    if (efficiencies.size() == axis_.size()) {
      efficiencies_ = efficiencies;
    } else {
      throw std::runtime_error("efficiencies must be the same size as the channel number");
    }
  }

  std::size_t NumberOfChannels() const { return n_channels_; }

  void Print() {
    std::cout << "Detector Multiplicity" << std::endl;
    auto max = *std::max_element(multiplicity_rec_.begin(), multiplicity_rec_.end());
    for (const auto &x : multiplicity_rec_) {
      auto nstars = x * 20. / max;
      for (unsigned int i = 0; i < nstars; ++i) {
        std::cout << "*";
      }
      std::cout << std::endl;
    }
  }

  std::string Name() const { return axis_.Name(); }

 private:
  std::size_t n_channels_;
  Qn::AxisD axis_;
  std::vector<double> phi_;
  std::vector<double> multiplicity_rec_;
  std::vector<double> multiplicity_true_;
  std::vector<double> efficiencies_;
  std::function<bool(double)> cuts_;
};

#endif
