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
#ifndef QNTOOLS_CORRELATIONFUNCTIONS_H_
#define QNTOOLS_CORRELATIONFUNCTIONS_H_

#include <cmath>
#include <complex>

namespace Qn::Correlation::TwoParticle {
inline auto xx(unsigned int h_a, unsigned int h_b) {
  return [h_a, h_b](const Qn::QVector &a, const Qn::QVector &b) {
    return a.x(h_a) * b.x(h_b);
  };
}
inline auto yy(unsigned int h_a, unsigned int h_b) {
  return [h_a, h_b](const Qn::QVector &a, const Qn::QVector &b) {
    return a.y(h_a) * b.y(h_b);
  };
}
inline auto yx(unsigned int h_a, unsigned int h_b) {
  return [h_a, h_b](const Qn::QVector &a, const Qn::QVector &b) {
    return a.y(h_a) * b.x(h_b);
  };
}
inline auto xy(unsigned int h_a, unsigned int h_b) {
  return [h_a, h_b](const Qn::QVector &a, const Qn::QVector &b) {
    return a.x(h_a) * b.y(h_b);
  };
}
inline auto ScalarProduct(unsigned int h_u, unsigned int h_Q) {
  return [h_u, h_Q](const Qn::QVector &u, const Qn::QVector &Q) {
    return u.x(h_u) * Q.x(h_Q) + u.y(h_u) * Q.y(h_Q);
  };
}

inline auto c2(unsigned int h_u) {
  return [h_u](const Qn::QVector &u) {
    double ret = 0.;
    auto m = u.sumweights();
    if (m < 2.) {
      ret = NAN;
    } else {
      auto Q = u.DeNormal();
      ret = (ScalarProduct(Q, Q, h_u) - m) / (m * (m - 1.));
    }
    return ret;
  };
}

inline auto d2(unsigned int h_u) {
  return [h_u](const Qn::QVector &p, const Qn::QVector &r) {
    double ret = 0.;
    auto m = p.sumweights();
    if (m < 2.) {
      return ret = NAN;
    } else {
      auto P = p.DeNormal();
      auto R = r.DeNormal();
      std::complex<double> p0{P.sumweights(), 0};
      std::complex<double> r0{R.sumweights(), 0};
      std::complex<double> p1{P.x(h_u), P.y(h_u)};
      std::complex<double> r1{R.x(h_u), R.y(h_u)};
      return (p1 * std::conj(r1) - p0).real() / (p0 * r0 - p0).real();
    }
  };
}

inline auto n2() {
  return [](const Qn::QVector &q) {
    auto m = q.sumweights();
    return m * (m - 1);
  };
}

inline auto nd2() {
  return [](const Qn::QVector &p, const Qn::QVector &r) {
    auto mp = p.sumweights();
    auto mr = r.sumweights();
    return mp * mr - mp;
  };
}

}  // namespace Qn::Correlation::TwoParticle

/**
 * mixed harmonics
 */
namespace Qn::Correlation::MixedHarmonics {
inline auto xxx(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.x(h_a) * Qb.x(h_b) * Qc.x(h_c);
  };
}
inline auto xyy(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.x(h_a) * Qb.y(h_b) * Qc.y(h_c);
  };
}
inline auto yxy(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.y(h_a) * Qb.x(h_b) * Qc.y(h_c);
  };
}
inline auto yyx(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.y(h_a) * Qb.y(h_b) * Qc.x(h_c);
  };
}
inline auto yyy(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.y(h_a) * Qb.y(h_b) * Qc.y(h_c);
  };
}
inline auto xyx(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.x(h_a) * Qb.y(h_b) * Qc.x(h_c);
  };
}
inline auto yxx(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.y(h_a) * Qb.x(h_b) * Qc.x(h_c);
  };
}
inline auto xxy(unsigned int h_a, unsigned int h_b, unsigned int h_c) {
  return [h_a, h_b, h_c](const Qn::QVector &u, const Qn::QVector &Qb,
                         const Qn::QVector &Qc) {
    return u.x(h_a) * Qb.x(h_b) * Qc.y(h_c);
  };
}
}  // namespace Qn::Correlation::MixedHarmonics

/**
 * Four particle correlation functions
 */
namespace Qn::Correlation::FourParticle {

inline auto n4() {
  return [](const Qn::QVector &q) {
    auto m = q.sumweights();
    return m * (m - 1.) * (m - 2.) * (m - 3.);
  };
}

inline auto nd4() {
  return [](const Qn::QVector &p, const Qn::QVector &r) {
    auto mp = p.sumweights();
    auto mr = r.sumweights();
    return (mp * mr - 3. * mp) * (mr - 1.) * (mr - 2.);
  };
}

inline auto d4(unsigned int h_u) {
  return [h_u](const Qn::QVector &p, const Qn::QVector &r) {
    const auto R = r.DeNormal();
    const auto RM = r.sumweights();
    const auto P = p.DeNormal();
    const auto PM = p.sumweights();
    const auto Q = p.DeNormal();
    const auto QM = p.sumweights();
    if (PM < 4.) return std::nan("");
    const auto r0 = std::complex<double>{RM, 0.};
    const auto q0 = std::complex<double>{QM, 0.};
    const auto p1 = std::complex<double>{P.x(h_u), P.y(h_u)};
    const auto r1 = std::complex<double>{R.x(h_u), R.y(h_u)};
    const auto r2 = std::complex<double>{R.x(2 * h_u), R.y(2 * h_u)};
    const auto q1 = std::complex<double>{P.x(h_u), P.y(h_u)};
    const auto q2 = std::complex<double>{P.x(2 * h_u), P.y(2 * h_u)};

    auto c4 = p1 * r1 * std::conj(r1) * std::conj(r1) -
              q2 * std::conj(r1) * std::conj(r1) - p1 * r1 * std::conj(r2) -
              2. * r1 * q0 * std::conj(r1) - 2. * p1 * r0 * std::conj(r1) +
              4. * q1 * std::conj(r1) + 2. * r0 * q0 + 2. * r1 * std::conj(q1) +
              2. * p1 * std::conj(r1) + q2 * std::conj(r2) - 6. * q0;
    return c4.real() / ((q0 * r0 - 3. * q0) * (r0 - 1.) * (r0 - 2.)).real();
  };
};

inline auto c4(unsigned int h_u) {
  return [h_u](const Qn::QVector &u) {
    float ret = 0.;
    auto Q = u.DeNormal();
    auto M = u.sumweights();
    if (M < 4.) {
      ret = NAN;
    } else {
      auto x = Q.x(h_u);
      auto y = Q.y(h_u);
      auto x2 = Q.x(2 * h_u);
      auto y2 = Q.y(2 * h_u);
      auto Q_mag = std::sqrt(x * x + y * y);
      auto Q_2n_mag = std::sqrt(x2 * x2 + y2 * y2);
      auto real = x2 * x * x - x2 * y * y + y2 * y * x + y2 * x * y;
      auto term_1 =
          (Q_mag * Q_mag * Q_mag * Q_mag + Q_2n_mag * Q_2n_mag - 2 * real) /
          (M * (M - 1) * (M - 2) * (M - 3));
      auto term_2 = (2 * (M - 2) * Q_mag * Q_mag - M * (M - 3)) /
                    (M * (M - 1) * (M - 2) * (M - 3));
      ret = term_1 - 2 * term_2;
    }
    return ret;
  };
}

}  // namespace Qn::Correlation::FourParticle
#endif  // QNTOOLS_CORRELATIONFUNCTIONS_H_
