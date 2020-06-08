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
#ifndef TRACKINGDETECTOR_H_
#define TRACKINGDETECTOR_H_

template<typename RandomEngine>
class TrackingDetector {
 public:
  TrackingDetector(const Qn::AxisD &axis, std::function<bool(double)> cuts) : axis_(axis),
                                                                              efficiencies_(axis.size()),
                                                                              cuts_(std::move(cuts)) {
    for (std::size_t i = 0; i < axis.size(); ++i) {
      efficiencies_.at(i) = 1.;
    }
  }

  bool Detect(double phi) {
    phi_ = phi;
    return cuts_(phi);
  }

  void FillDataRec(RandomEngine &engine, double *values_array, std::size_t position_phi, std::size_t position_weight) {
    auto bin = axis_.FindBin(phi_);
    if (bin > -1) {
      values_array[position_phi] = phi_;
      if (detection_efficiency_(engine) < efficiencies_[bin]) {
        values_array[position_weight] = 1.;
      } else {
        values_array[position_weight] = 0.;
      }
    }
  }

  std::string Name() const { return axis_.Name(); }

  void FillDataTruth(double *values_array, std::size_t position_phi, std::size_t position_weight) {
    values_array[position_phi] = phi_;
    values_array[position_weight] = 1.0;
  }

  void SetChannelEfficencies(std::vector<double> efficiencies) {
    if (efficiencies.size() == axis_.size()) {
      efficiencies_ = efficiencies;
    } else {
      throw std::runtime_error("efficiencies must be the same size as the channel number");
    }
  }

 private:
  std::uniform_real_distribution<> detection_efficiency_{0, 1};
  Qn::AxisD axis_;
  double phi_;
  std::vector<double> efficiencies_;
  std::function<bool(double)> cuts_;
};

#endif
