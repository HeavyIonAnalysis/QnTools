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

#include <TLegend.h>
#include <bitset>
#include <iostream>
#include <vector>

#include "Rtypes.h"

#include "ReSamples.hpp"
#include "Statistic.hpp"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"

namespace Qn {

class Stats {
 public:
  enum Settings {
    CORRELATEDERRORS = BIT(16),
    CONCATENATE_SUBSAMPLES = BIT(17),
    ASYMMERRORS = BIT(18)
  };

  enum class State {
    MOMENTS,
    MEAN_ERROR
  };

  enum class Weights {
    REFERENCE,
    OBSERVABLE
  };

  using size_type = std::size_t;

  Stats() = default;
  Stats(double mean, double error, double weight) : state_(State::MEAN_ERROR),
                                                    mean_(mean),
                                                    error_(error),
                                                    weight_(weight) {
  }
  virtual ~Stats() = default;

  double N() const { return statistic_.N(); }

  double Neff() const { return statistic_.Neff(); }

  double SumWeights() const { return statistic_.SumWeights(); }

  bool IsObservable() const { return Qn::Stats::Weights::OBSERVABLE == weights_flag; }

  double RatioOfErrors() const {
    auto bootstrap_error = resamples_.GetConfidenceInterval(mean_, ReSamples::CIMethod::normal).Uncertainty();
    double error = 0;
    if (state_ == State::MOMENTS) {
      error = statistic_.MeanError();
    } else if (state_ == State::MEAN_ERROR) {
      error = error_;
    }
    std::cout << bootstrap_error << " " << error << std::endl;
    return bootstrap_error / error;
  }

  double Weight() const {
    if (state_ == Qn::Stats::State::MEAN_ERROR) {
      return weight_;
    } else {
      return statistic_.SumWeights();
    }
  }

  double Mean() const {
    double mean = 0;
    switch (state_) {
      case State::MOMENTS:
        mean = statistic_.Mean();
        break;
      case State::MEAN_ERROR:
        mean = mean_;
        break;
    }
    return mean;
  }

  double LowerMeanError() const {
    double lower_error = 0;
    if (bits_ & Settings::ASYMMERRORS) {
      double mean = mean_;
      if (state_ != State::MEAN_ERROR) { mean = statistic_.Mean(); }
      auto ci = resamples_.GetConfidenceInterval(mean, ReSamples::CIMethod::pivot);
      lower_error = mean - ci.lower_limit;
    } else {
      lower_error = MeanError();
    }
    return lower_error;
  }

  double UpperMeanError() const {
    double upper_error = 0;
    if (bits_ & Settings::ASYMMERRORS) {
      double mean = mean_;
      if (state_ != State::MEAN_ERROR) { mean = statistic_.Mean(); }
      auto ci = resamples_.GetConfidenceInterval(mean, ReSamples::CIMethod::pivot);
      upper_error = ci.upper_limit - mean;
    } else {
      upper_error = MeanError();
    }
    return upper_error;
  }

  double MeanError() const {
    double error = 0;
    if (bits_ & Settings::CORRELATEDERRORS) {
      error = MeanErrorBoot();
    } else {
      switch (state_) {
        case State::MOMENTS:
          error = statistic_.MeanError();
          break;
        case State::MEAN_ERROR:
          error = error_;
          break;
      }
    }
    return error;
  }

  double MeanErrorStat() const {
    double error = 0;
    switch (state_) {
      case State::MOMENTS:
        error = statistic_.MeanError();
        break;
      case State::MEAN_ERROR:
        error = error_;
        break;
    }
    return error;
  }

  double MeanErrorBoot() const {
    if (state_ != State::MEAN_ERROR) {
      auto t_samples = resamples_;
      t_samples.CalculateMeans();
      return t_samples.GetConfidenceInterval(statistic_.Mean(), ReSamples::CIMethod::normal).Uncertainty();
    } else {
      return resamples_.GetConfidenceInterval(mean_, ReSamples::CIMethod::normal).Uncertainty();
    }
  }

  void CalculateMeanAndError() {
    if (state_ != State::MEAN_ERROR) {
      state_ = State::MEAN_ERROR;
      mean_ = statistic_.Mean();
      error_ = statistic_.MeanError();
      weight_ = statistic_.SumWeights();
      resamples_.CalculateMeans();
    }
  }

  friend Stats Merge(const Stats &, const Stats &);
  friend Stats MergeBins(const Stats &, const Stats &);
  friend Stats operator+(const Stats &, const Stats &);
  friend Stats operator-(const Stats &, const Stats &);
  friend Stats operator*(const Stats &, const Stats &);
  friend Stats operator*(const Stats &, double);
  friend Stats operator*(const Stats &, std::pair<double, double>);
  friend Stats operator/(const Stats &, double);
  friend Stats operator*(double, const Stats &);
  friend Stats operator/(const Stats &, const Stats &);
  friend Stats Sqrt(const Stats &);
  friend Stats Abs(const Stats &);
  friend Stats PowSqrt(const Stats &, unsigned int);
  friend Stats OllitraultExtrapolation( const Stats &, unsigned int );

  template<typename SAMPLES>
  inline void Fill(const double value, const double weight, SAMPLES &&samples) {
    if (!std::isnan(value)) {
      resamples_.Fill(value, weight, std::forward<SAMPLES>(samples));
      statistic_.Fill(value, weight);
    }
  }

  void SetNumberOfReSamples(size_type nsamples) {
    resamples_.SetNumberOfSamples(nsamples);
  }

  void SetWeights(Weights weights) { weights_flag = weights; }
  State GetState() const { return state_; }

  bool TestBit(unsigned int bit) const { return static_cast<bool>(bits_ & bit); }
  void SetBits(unsigned int bits) { bits_ = bits; }
  void ResetBits(unsigned int bits) { bits_ &= ~(bits & 0x00ffffff); }

  const Statistic &GetStatistic() const { return statistic_; }

  size_type GetNSamples() const { return resamples_.size(); }
  const ReSamples &GetReSamples() const { return resamples_; }

  TCanvas *CIvsNSamples(const int nsteps = 10) const;

 private:
  ReSamples resamples_;                                /// resamples used for error calculation
  Statistic statistic_;                                /// Used in the state of MOMENTS
  unsigned int bits_ = 0 | Qn::Stats::CORRELATEDERRORS;// configuration bits
  State state_ = State::MOMENTS;                       /// state of the calculations
  Weights weights_flag = Weights::REFERENCE;           /// weighting of the Stats
  bool mergeable_ = true;                              /// flag is false if an operation leading to undefined weights has been performed
  /// Used in the state of MEAN_ERROR
  double mean_ = 0.;  /// mean
  double error_ = 0.; /// uncertainty
  double weight_ = 0.;/// relative weight for rebinning

  /// \cond CLASSIMP
  ClassDef(Stats, 4);
  /// \endcond
};

Stats MergeBins(const Stats &, const Stats &);
Stats Merge(const Stats &, const Stats &);
Stats operator+(const Stats &, const Stats &);
Stats operator-(const Stats &, const Stats &);
Stats operator*(const Stats &, const Stats &);
Stats operator*(const Stats &, double);
Stats operator*(const Stats &, std::pair<double, double>);
Stats operator/(const Stats &, double);
Stats operator*(double, const Stats &);
Stats operator/(const Stats &, const Stats &);
Stats Abs(const Stats &);
Stats Sqrt(const Stats &);
Stats PowSqrt(const Stats &, unsigned int);
// fuction for extrapolation RND-sub resolution to full event
Stats OllitraultExtrapolation( const Stats &, unsigned int );
}// namespace Qn

#endif//FLOW_STATS_H
