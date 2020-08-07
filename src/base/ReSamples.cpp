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

#include "ReSamples.hpp"

namespace Qn {

ConfidenceInterval ConfidenceIntervalPercentile(std::vector<double> means) {
  if (means.empty()) { throw std::out_of_range("vector of samples is empty."); }
  const double alpha = 0.3173;
  std::sort(means.begin(), means.end());
  auto lowerpos = std::round(means.size() * alpha / 2.);
  auto upperpos = std::floor(means.size() * (1 - alpha / 2.));
  return ConfidenceInterval{means.at(lowerpos), means.at(upperpos)};
}

ConfidenceInterval ConfidenceIntervalPivot(std::vector<double> means, const double real_mean) {
  if (means.empty()) { throw std::out_of_range("vector of samples is empty."); }
  const double alpha = 0.3173;
  std::sort(means.begin(), means.end());
  auto upperpos = std::round(means.size() * alpha / 2.);
  auto lowerpos = std::floor(means.size() * (1 - alpha / 2.));
  return ConfidenceInterval{2 * real_mean - means.at(lowerpos), 2 * real_mean - means.at(upperpos)};
}

ConfidenceInterval ConfidenceIntervalNormal(std::vector<double> means, const double real_mean) {
  if (means.empty()) { throw std::out_of_range("vector of samples is empty."); }
  Statistic stats;
  for (const auto &mean : means) {
    stats.Fill(mean, 1.0);
  }
  auto stddev = std::sqrt(stats.Variance());
  return ConfidenceInterval{real_mean - stddev, real_mean + stddev};
}

ReSamples ReSamples::Addition(const ReSamples &a, const ReSamples &b) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    result.means_[i] += b.means_[i];
  }
  return result;
}

ReSamples ReSamples::Subtraction(const ReSamples &a, const ReSamples &b) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    result.means_[i] -= b.means_[i];
  }
  return result;
}
ReSamples ReSamples::Division(const ReSamples &a, const ReSamples &b) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    result.means_[i] /= b.means_[i];
  }
  return result;
}

ReSamples ReSamples::Multiplication(const ReSamples &a, const ReSamples &b) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    result.means_[i] *= b.means_[i];
  }
  return result;
}

ReSamples ReSamples::Scaling(const ReSamples &a, const double scale) {
  ReSamples result(a);
  for (auto &mean : result.means_) {
    mean *= scale;
  }
  return result;
}

ReSamples ReSamples::Sqrt(const ReSamples &a) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    auto root = sqrt(fabs(a.means_[i]));
    result.means_[i] = std::signbit(a.means_[i]) ? -1. * root : root;
  }
  return result;
}

ReSamples ReSamples::Abs(const ReSamples &a) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    result.means_[i] = fabs(a.means_[i]);
  }
  return result;
}

ReSamples ReSamples::PowSqrt(const ReSamples &a, unsigned int k) {
  ReSamples result(a);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    auto res = std::pow(fabs(a.means_[i]), 1. / k);
    result.means_[i] = std::signbit(a.means_[i]) ? -1. * res : res;
  }
  return result;
}

ReSamples ReSamples::OllitraultExtrapolation(const ReSamples &a, unsigned int k){
  ReSamples result(a);
  // factorial lambda-function for Bessel's
  auto fact = [](int x) {
    double result = 1;
    for (int i = 1; i <= x; ++i)
      result *= i;
    return result;
  };
  // Modified Bessel's function for resolution calculation
  auto I = [fact](double nu, double z) {
    double result = 0;
    for (int i = 0; i < 10; ++i)
      result +=
          pow(z / 2.0, 2 * i + nu) / (fact(i) * tgamma(i + nu + 1));
    return result;
  };
  // Resolution as the function of chi = vn/sqrt(M)
  auto R = [I](double chi, double m) {
    double chi2_2 = chi * chi / 2;
    double result = sqrt( M_PI_2 * chi * exp(-chi2_2) *
        (I((m -1.0)/2, chi2_2) + I( (m +1.0)/2, chi2_2) ) );
    return result;
  };
  // firstly need to solve equation R(chi) = mean to find chi
  auto f = [R, k](double chi, double value) { return R(chi, k)-sqrt(value); };
  // introducing dichotomy method to solve the equation and find resolution for found chi
  auto dichotomy = [f, R]( double value ) {
    double a = 0.0;
    double b = 3.0;
    int i = 0;
    while (fabs(a - b) > pow(10, -6)) {
      double c = (a + b) / 2;
      double fc = f(c, value);
      if (fc == 0) break;
      double fa = f(a, value);
      double fb = f(b, value);
      if (fa * fc < 0.0) {
        b = c;
        i++;
        continue;
      }
      if (fb * fc < 0.0) {
        a = c;
        i++;
        continue;
      }
    }
    double chi = (a + b) / 2;
    double res = R(sqrt(2) * chi, 2);
    return res;
  };
  for (size_t i = 0; i < result.means_.size(); ++i) {
    auto res = dichotomy(a.means_.at(i));
    result.means_[i] = res;
  }
  return result;
}

ReSamples ReSamples::Merge(const ReSamples &a, const ReSamples &b, bool merge_weights) {
  ReSamples result(b);
  for (size_t i = 0; i < result.means_.size(); ++i) {
    ValueType a_mean = 0.;
    ValueType a_weight = 0.;
    if (i < a.size()) {
      a_mean = a.means_[i];
      a_weight = a.weights_[i];
    }
    auto b_mean = b.means_[i];
    auto b_weight = b.weights_[i];
    result.means_[i] = (a_mean * a_weight + b_mean * b_weight) / (a_weight + b_weight);
    if (merge_weights) {
      result.weights_[i] += b.weights_[i];
    }
  }
  return result;
}

ReSamples ReSamples::MergeStatistics(const ReSamples &a, const ReSamples &b) {
  ReSamples result(b);
  for (size_t i = 0; i < result.statistics_.size(); ++i) {
    Statistic stat_a;
    if (i < a.size()) stat_a = a.statistics_[i];
    auto stat_b = b.statistics_[i];
    result.statistics_[i] = Qn::Merge(stat_a, stat_b);
  }
  return result;
}

ReSamples ReSamples::Concatenate(const ReSamples &a, const ReSamples &b) {
  ReSamples result(a);
  result.means_.insert(result.means_.end(), b.means_.begin(), b.means_.end());
  result.weights_.insert(result.weights_.end(), b.weights_.begin(), b.weights_.end());
  result.statistics_.insert(result.statistics_.end(), b.statistics_.begin(), b.statistics_.end());
  return result;
}

std::pair<TGraph *, TGraph *> ReSamples::CIvsNSamples(double mean,
                                                      ReSamples::CIMethod method,
                                                      unsigned int nsteps) const {
  auto graphup = new TGraph(nsteps);
  auto graphlo = new TGraph(nsteps);
  const size_t step = means_.size() / nsteps;
  for (unsigned int i = 0; i < nsteps; ++i) {
    auto ci = ConfidenceIntervalNSamplesMethod(mean, step * (i + 1), method);
    graphlo->SetPoint(i, (i + 1) * step, ci.lower_limit);
    graphup->SetPoint(i, (i + 1) * step, ci.upper_limit);
  }
  return std::make_pair(graphlo, graphup);
}

}// namespace Qn