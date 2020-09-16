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

#ifndef FLOW_QVECTOR_H
#define FLOW_QVECTOR_H

#include <array>
#include <bitset>
#include <cmath>
#include <stdexcept>// std::out_of_range

#include "Rtypes.h"

namespace Qn {
/**
 * Struct of a Q-vector of a single harmonic with x and y component
 */
struct QVec {
  QVec() = default;
  QVec(float x, float y) : x(x), y(y) {}
  float x{0.};
  float y{0.};
};

/**
 * Implements vector addition
 * @param a input Q-vector
 * @param b input Q-vector
 * @return sum of Q-vectors
 */
inline QVec operator+(const QVec a, const QVec b) { return {a.x + b.x, a.y + b.y}; }
/**
 * Implements vector subtraction
 * @param a input Q-vector
 * @param b input Q-vector
 * @return difference of Q-vectors
 */
inline QVec operator-(const QVec a, const QVec b) { return {a.x - b.x, a.y - b.y}; }
/**
 * Implements scalar division
 * @param a input Q-vector
 * @param s scalar divisor
 * @return divided Q-vector
 */
inline QVec operator/(const QVec a, const float s) { return {a.x / s, a.y / s}; }
/**
 * Implements scalar multiplication
 * @param a input Q-vector
 * @param s input multiplier
 * @return multiplied Q-vector
 */
inline QVec operator*(const QVec a, const float s) { return {a.x * s, a.y * s}; }
/**
 * Returns Euclidean norm of the Q-vector
 * @param a input Q-vector
 * @return norm of the Q-vector
 */
inline float norm(QVec a) { return std::sqrt(a.x * a.x + a.y * a.y); }

/**
 * @class QVector
 * @brief It carries the information of multiple harmonics of a specified normalization and correction step.
 */
class QVector {
 public:
  static constexpr int kmaxharmonics = 8;
  static constexpr float kminimumweight = 1e-6;
  static constexpr float kPi = 3.14159265358979323846;
  /**
   * @enum available correction steps in the order of application
   */
  enum CorrectionStep {
    RAW,
    PLAIN,
    RECENTERED,
    TWIST,
    RESCALED,
    ALIGNED
  };

  /**
   * @enum available normalizations
   */
  enum class Normalization {
    NONE,    ///< \f$ \mbox{Q'} = \mbox{Q}\f$
    SQRT_M,  ///< \f$ \mbox{Q'} = \frac{\mbox{Q}}{\sqrt{\mbox{M}}} \f$
    M,       ///< \f$ \mbox{Q'} = \frac{\mbox{Q}}{\mbox{M}} \f$
    MAGNITUDE///< \f$ \mbox{Q'} = \frac{\mbox{Q}}{|\mbox{Q}|} \f$
  };

  static std::string GetNormalizationName(Normalization normal) {
    switch (normal){
      case Normalization::NONE: return "None";
      case Normalization::M: return "M";
      case Normalization::SQRT_M: return "SqrtM";
      case Normalization::MAGNITUDE: return "Magnitude";
      default: return "None";
    }
  }

  QVector() = default;

  virtual ~QVector() = default;

  /**
   * Constructor
   * @param bits set of activated harmonics
   * @param step current correction step
   */
  QVector(std::bitset<kmaxharmonics> bits, CorrectionStep step) : correction_step_(step),
                                                                  bits_(bits) {
    q_.resize(bits.count());
    maximum_harmonic_ = highestharmonic();
  }

  /**
 * Constructor
 * @param bits set of activated harmonics
 * @param step current correction step
 * @param norm normalization of Q-vector
 */
  QVector(std::bitset<kmaxharmonics> bits, CorrectionStep step, Qn::QVector::Normalization norm) : norm_(norm),
                                                                                                   correction_step_(step),
                                                                                                   bits_(bits) {
    q_.resize(bits.count());
    maximum_harmonic_ = highestharmonic();
  }

  /**
   * Copy constructor
   * @param rhs Q-vector to be copied
   */
  QVector(const QVector &rhs) = default;

  /**
   * Copies the harmonics of the other Q-vector
   * @param other Q-vector whose harmonics are copied to the current Q-vector.
   */
  void CopyHarmonics(const QVector &other) {
    this->bits_ = other.bits_;
    this->q_.resize(other.q_.size());
  }

  void InitializeHarmonics() {
    maximum_harmonic_ = highestharmonic();
  }

  /**
   * Reset the Q-vector. Is called before the Q-vector is filled in every event.
   */
  void Reset() {
    n_ = 0;
    sum_weights_ = 0.;
    quality_ = false;
    for (auto &q : q_) { q = QVec(); }
  }

  /**
   * Get a c-style array of all activated harmonics. Used to create the correction histograms
   * @param store Pointer to the first element of the array.
   *        Array needs to have the exact size given by the number of harmonics.
   */
  void GetHarmonicsMap(Int_t *store) const {
    unsigned int iharmonics = 0;
    for (unsigned char h = 1; h <= maximum_harmonic_; h++) {
      if (bits_.test(h - 1)) {
        store[iharmonics] = h;
        iharmonics++;
      }
    }
  }

  /**
   * Returns the set of activated harmonics
   * @return std::bitset with all activaed harmonics set to 1 and the others set to 0.
   */
  std::bitset<kmaxharmonics> GetHarmonics() const { return bits_; }

  /**
   * Returns the number of harmonics
   * @return total count of activated harmonics
   */
  unsigned int GetNoOfHarmonics() const { return bits_.count(); }

  /**
   * Sets the harmonic multiplier. Used for some correction steps
   * @param mult harmonic multiplier
   */
  void SetHarmonicMultiplier(const unsigned char mult) { harmonic_multiplier_ = mult; }

  /**
   * Get the harmonic multiplier
   * @return the harmonic multiplier of the current Q-vector.
   */
  unsigned char GetHarmonicMultiplier() const { return harmonic_multiplier_; }

  /**
   * Sets the quality to good.
   * @param quality set to true if the Q-vector satisfies the conditions and can be used in the further processing steps.
   */
  void SetGood(bool quality) { quality_ = quality; }

  /**
   * Checks the quality of the Q-vector.
   * @return true if the Q-vector can be used in the further processing steps.
   */
  bool IsGoodQuality() const { return quality_; }

  /**
   * Sets the quality of the Q-vector. Is called after the q-vector creation is complete.
   * Q-vector needs to have at minimum 1 contributor.
   */
  void CheckQuality() { quality_ = 0 < n_; }

  /**
   * Activate a harmonic in the qvector.
   * @param i specified harmonic. needs to be smaller than the maximum kmaxharmonics.
   */
  void ActivateHarmonic(const unsigned int i) {
    bits_.set(i - 1);
    q_.resize(static_cast<size_t>(bits_.count()));
  }

  /**
   * Set the normalization of the Q-vector
   * @param norm normalization
   */
  void SetNormalization(const Normalization norm) { norm_ = norm; }

  /**
   * Update correction step
   */
  void UpdateCorrectionStep(const CorrectionStep step) {
    correction_step_ = step;
  }

  /**
   * Get the correction step
   * @return correction step of the Q-vector
   */
  CorrectionStep GetCorrectionStep() const { return correction_step_; }

  /**
   * Returns x-component of Q-vector of the i-th harmonic.
   * Throws exception, when the harmonic is out of the range.
   * @param i harmonic i of the Q-vector
   * @return x-component
   */
  inline float x(const unsigned int i) const {
    if (bits_.test(i - 1)) {
      auto position = std::bitset<kmaxharmonics>(bits_ & std::bitset<kmaxharmonics>((1UL << (i)) - 1)).count() - 1;
      return q_[position].x;
    } else {
      throw std::out_of_range("harmonic not in range.");
    }
  }
  /**
   * Returns y-component of Q-vector of the i-th harmonic
   * Throws exception, when the harmonic is out of the range.
   * @param i harmonic i of the Q-vector
   * @return y-component
   */
  inline float y(const unsigned int i) const {
    if (bits_.test(i - 1)) {
      auto position = std::bitset<kmaxharmonics>(bits_ & std::bitset<kmaxharmonics>((1UL << (i)) - 1)).count() - 1;
      return q_[position].y;
    } else {
      throw std::out_of_range("harmonic not in range.");
    }
  }

  /**
   * Sets the x-component of the Q-vector of the i-th harmonic.
   * @param i harmonic i of the Q-vector
   * @param x new x component.
   */
  inline void SetX(const unsigned int i, double x) {
    auto position = std::bitset<kmaxharmonics>(bits_ & std::bitset<kmaxharmonics>((1UL << (i)) - 1)).count() - 1;
    q_[position].x = x;
  }

  /**
   * Sets the y-component of the Q-vector of the i-th harmonic.
   * @param i harmonic i of the Q-vector
   * @param y new y component.
   */
  inline void SetY(const unsigned int i, double y) {
    auto position = std::bitset<kmaxharmonics>(bits_ & std::bitset<kmaxharmonics>((1UL << (i)) - 1)).count() - 1;
    q_[position].y = y;
  }

  /**
 * Sets the x-component of the Q-vector of the i-th harmonic.
 * @param i harmonic i of the Q-vector
 * @param x new x component.
 */
  inline void SetQ(const unsigned int i, double x, double y) {
    auto position = std::bitset<kmaxharmonics>(bits_ & std::bitset<kmaxharmonics>((1UL << (i)) - 1)).count() - 1;
    q_[position].x = x;
    q_[position].y = y;
  }

  /**
   * Copy the number of contributors of the other Q-vector.
   * This is used during the construction of the corrected Q-vector.
   * @param other Q-vector whose properties are being copied.
   */
  void CopyNumberOfContributors(const QVector &other) {
    norm_ = other.norm_;
    quality_ = other.quality_;
    n_ = other.n_;
    sum_weights_ = other.sum_weights_;
  }

  /**
   * Gets the first harmonic.
   * @return first harmonic number. Returns 0 if none are found.
   */
  int GetFirstHarmonic() const {
    for (Int_t h = 1; h < maximum_harmonic_ + 1; h++) {
      if (bits_.test(h - 1)) {
        return h;
      }
    }
    return -1;
  }

  /**
   * Gets the next harmonic using the current one as input.
   * @param harmonic current harmonic.
   * @return returns the next harmonic. Returns -1 if none are found.
   */
  int GetNextHarmonic(const unsigned char harmonic) const {
    for (unsigned char h = harmonic + 1; h < maximum_harmonic_ + 1; h++) {
      if (bits_.test(h - 1)) {
        return h;
      }
    }
    return -1;
  }

  /**
   * Returns the Magnitude or euclidean norm of the Q-vector.
   * @param i i-th harmonic of the Q-vector.
   * @return euclidean norm of the Q-vector.
   */
  inline float mag(const unsigned int i) const { return sqrt(x(i) * x(i) + y(i) * y(i)); }
  /**
   * Returns the sum of weights of the Q-Vector.
   * @return Sum of weights of the Q-Vector.
   */
  inline float sumweights() const { return sum_weights_; }
  /**
   * Returns the number of contributors of the Q-Vector.
   * @return number of contributors of the Q-Vector.
   */
  inline float n() const { return n_; }
  /**
   * Returns the plane angle of the i-th harmonic.
   * @param i i-th harmonic
   * @return event plane angle
   */
  inline float psi(const unsigned int i) const { return std::atan2(y(i), x(i)) + kPi; }
  /**
   * Get the normalization method. of the Q-vector
   * @return normalization method
   */
  inline Normalization GetNorm() const { return norm_; }
  /**
   * Qvector addition operator for all harmonics.
   * @param a input Q-vector
   * @param b input Q-vector
   * @return  sum of Q-vectors
   */
  friend QVector operator+(QVector a, QVector b);

  QVector Normal(const Normalization norm) const;

  QVector DeNormal() const;

  /**
   * Adds a new data vector to the qvector.
   * @param phi angle of the particle or channel.
   * @param weight weight of the particle or channel e.g. channel multiplicity.
   */
  inline void Add(const double phi, const double weight) {
    if (weight < kminimumweight) return;
    unsigned int pos = 0;
    for (unsigned int h = 1; h <= maximum_harmonic_; ++h) {
      if (bits_.test(h - 1)) {
        q_[pos].x += (weight * std::cos(h * harmonic_multiplier_ * phi));
        q_[pos].y += (weight * std::sin(h * harmonic_multiplier_ * phi));
        ++pos;
      }
    }
    sum_weights_ += weight;
    n_ += 1;
  }

  /**
 * Adds a new data vector to the qvector.
 * @param phi angle of the particle or channel.
 * @param offset offset of the phi channel.
 * @param weight weight of the particle or channel e.g. channel multiplicity.
 */
  inline void Add(const double phi, const double offset, const double weight) {
    if (weight < kminimumweight) return;
    unsigned int pos = 0;
    for (unsigned int h = 1; h <= maximum_harmonic_; ++h) {
      if (bits_.test(h - 1)) {
        q_[pos].x += (weight * std::cos(h * harmonic_multiplier_ * phi) * offset);
        q_[pos].y += (weight * std::sin(h * harmonic_multiplier_ * phi) * offset);
        ++pos;
      }
    }
    sum_weights_ += weight;
    n_ += 1;
  }

  /**
   * Returns the highest harmonic configured in the Q-vector.
   * @return highest harmonic.
   */
  inline unsigned int highestharmonic() const {
    unsigned char val = bits_.to_ulong();
    if (val == 0) return 255;
    if (val == 1) return 0 + 1;
    unsigned char ret = 0;
    while (val > 1) {
      val >>= 1;
      ret++;
    }
    return (ret + 1);
  }

 private:
  Normalization norm_ = Normalization::NONE;            ///< normalization method
  CorrectionStep correction_step_ = CorrectionStep::RAW;///< correction step defined by enumerator
  int n_ = 0;                                           ///< number of data vectors contributing to the q vector
  float sum_weights_ = 0.0;                             ///< sum of weights
  std::bitset<kmaxharmonics> bits_{};                   ///< Bitset for keeping track of the harmonics
  std::vector<QVec> q_;                                 ///< array of qvectors for the different harmonics
  /**
   * Data members only used during the construction and correction of the Q-vectors.
   * They are not saved to the root file, as they are not used to read the data.
   */
  bool quality_ = false;                 //!<! quality of the Q-vector (only used during construction)
  unsigned char maximum_harmonic_ = 0;   //!<! maximum harmonic
  unsigned char harmonic_multiplier_ = 1;//!<! harmonic multiplier (used for some correction steps)

  /// \cond CLASSIMP
  ClassDef(QVector, 12);
  /// \endcond
};

inline double ScalarProduct(QVector a, QVector b, unsigned int harmonic) {
  return a.x(harmonic) * b.x(harmonic) + a.y(harmonic) * b.y(harmonic);
}

static constexpr std::array<const char *, 6> kCorrectionStepNamesArray = {
    "RAW",
    "PLAIN",
    "RECENTERED",
    "TWIST",
    "RESCALED",
    "ALIGNED"};

}// namespace Qn
#endif//FLOW_QVECTOR_H
