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

#ifndef FLOW_RESAMPLES_H
#define FLOW_RESAMPLES_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "TGraphAsymmErrors.h"
#include "TH2.h"

#include "Statistic.hpp"

namespace Qn {

struct ConfidenceInterval {
  double lower_limit;
  double upper_limit;
  double Uncertainty() const {
    return (upper_limit - lower_limit) / 2.;
  }
};

ConfidenceInterval ConfidenceIntervalPercentile(std::vector<double>);
ConfidenceInterval ConfidenceIntervalPivot(std::vector<double>, double);
ConfidenceInterval ConfidenceIntervalNormal(std::vector<double>, double);

class ReSamples {
  using size_type = std::size_t;

 public:
  using ValueType = double;

  enum class CIMethod {
    percentile,
    pivot,
    normal
  };

  ReSamples() = default;
  explicit ReSamples(unsigned int size) : statistics_(size), means_(size), weights_(size) {}
  virtual ~ReSamples() = default;
  ReSamples(const ReSamples &sample) {
    statistics_ = sample.statistics_;
    means_ = sample.means_;
    weights_ = sample.weights_;
  }
  ReSamples(ReSamples &&sample) = default;
  ReSamples &operator=(const ReSamples &sample) = default;

  void SetNumberOfSamples(unsigned int i) {
    statistics_.resize(i);
    means_.resize(i);
    weights_.resize(i);
  }
  size_type size() const { return means_.size(); }

  const ValueType &GetSampleMean(int i) const { return means_.at(i); }

  std::vector<double> GetMeans() const { return means_; }

  ConfidenceInterval GetConfidenceInterval(const double mean, CIMethod method) const {
    return ConstructConfidenceInterval(means_, mean, method);
  }

  template<typename SAMPLES>
  void Fill(const double value, const double weight, SAMPLES &&sample_multiplicities_) {
    using ArrayValueType = std::decay_t<decltype(std::declval<SAMPLES &>()[0])>;
    for (std::size_t i = 0; i < sample_multiplicities_.size(); ++i) {
      for (ArrayValueType j = 0; j < sample_multiplicities_[i]; ++j) {
        statistics_[i].Fill(value, weight);
      }
    }
  }

  void CalculateMeans() {
    if (!using_means_) {
      unsigned int i = 0;
      for (auto &statistic : statistics_) {
        double mean = statistic.Mean();
        double weight = statistic.SumWeights();
        means_.at(i) = mean;
        weights_.at(i) = weight;
        ++i;
        using_means_ = true;
      }
    }
  }

  std::pair<TGraph *, TGraph *> CIvsNSamples(double mean, CIMethod method, unsigned int nsteps = 10) const;

  void ScatterGraph(TGraph &graph, double offset, double width) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> position(offset - width, offset + width);
    for (auto mean : means_) {
      graph.SetPoint(graph.GetN(), position(gen), mean);
    }
  }

  void FillHistogram(TH2 &histogram, double x_value) const {
    for (const auto &mean : means_) {
      histogram.Fill(x_value, mean);
    }
  }

  static ReSamples Addition(const ReSamples &, const ReSamples &);

  static ReSamples Subtraction(const ReSamples &, const ReSamples &);

  static ReSamples Division(const ReSamples &, const ReSamples &);

  static ReSamples Multiplication(const ReSamples &, const ReSamples &);

  static ReSamples Scaling(const ReSamples &, double);

  static ReSamples Sqrt(const ReSamples &);

  static ReSamples Abs(const ReSamples &);

  static ReSamples PowSqrt(const ReSamples &, unsigned int);

  static ReSamples Merge(const ReSamples &, const ReSamples &, bool);

  static ReSamples MergeStatistics(const ReSamples &, const ReSamples &);

  static ReSamples Concatenate(const ReSamples &, const ReSamples &);

 private:
  ConfidenceInterval ConfidenceIntervalNSamplesMethod(const double mean,
                                                      const size_type nsamples,
                                                      CIMethod method = CIMethod::pivot) const {
    std::vector<ValueType> bsmeans = means_;
    bsmeans.resize(nsamples);
    return ConstructConfidenceInterval(bsmeans, mean, method);
  }

  inline ConfidenceInterval ConstructConfidenceInterval(std::vector<ValueType> means,
                                                        const double mean,
                                                        CIMethod method) const {
    ConfidenceInterval interval{};
    switch (method) {
      case CIMethod::percentile:
        interval = ConfidenceIntervalPercentile(means);
        break;
      case CIMethod::pivot:
        interval = ConfidenceIntervalPivot(means, mean);
        break;
      case CIMethod::normal:
        interval = ConfidenceIntervalNormal(means, mean);
        break;
    }
    return interval;
  }

  bool using_means_ = false;
  std::vector<Statistic> statistics_;
  std::vector<ValueType> means_;
  std::vector<ValueType> weights_;

  /// \cond CLASSIMP
  ClassDef(ReSamples, 2);
  /// \endcond
};
}// namespace Qn

#endif//FLOW_SAMPLE_H
