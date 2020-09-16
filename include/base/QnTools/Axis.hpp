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

#ifndef FLOW_QNAXIS_H
#define FLOW_QNAXIS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include "Rtypes.h"

namespace Qn {
/**
 * @brief Parameter axis with variable bin widths
 *
 * Basic axis implementation
 */
template<typename T>
class Axis {
 public:
  using ValueType = T;

  Axis() = default;         ///< default constructor
  virtual ~Axis() = default;///< default destructor

  /**
   * Constructor for variable or fixed bin width.
   * @param name name of axis.
   * @param bin_edges vector of bin edges. starting with lowest bin edge and ending with uppermost bin edge.
   */
  Axis(std::string name, std::vector<T> bin_edges)
      : name_(std::move(name)), bin_edges_(std::move(bin_edges)) {}

  /**
   * Constructor for fixed bin width. Calculates bin width automatically and sets bin edges.
   * @param name name of axis
   * @param nbins number of bins
   * @param lowbin lowest bin edge
   * @param upbin uppermost bin edge
   */
  Axis(std::string name, const int nbins, const T lowbin, const T upbin)
      : name_(std::move(name)) {
    for (int i = 0; i < nbins + 1; ++i) {
      T bin_width = (upbin - lowbin) / (T) nbins;
      bin_edges_.push_back(lowbin + i * bin_width);
    }
  }

  Axis(const Axis<T> &axis) : name_(axis.name_), bin_edges_(axis.bin_edges_) {}
  bool operator==(const Axis &axis) const {
    auto same_bins = false;
    if (size() == axis.size()) {
      std::vector<bool> result;
      std::transform(std::begin(bin_edges_),
                     std::end(bin_edges_),
                     std::begin(axis.bin_edges_),
                     std::back_inserter(result),
                     [](const T &a, const T &b) {
                       double epsilon = std::fabs(a*1e-3);
                       return std::fabs(a-b) < epsilon;
                     });
      same_bins = std::all_of(std::begin(result), std::end(result), [](bool a) { return a; });
    }
    return name_ == axis.name_ && same_bins;
  }

  bool operator!=(const Axis &axis) const {
    return !operator==(axis);
  }

  typedef typename std::vector<T>::const_iterator citerator;
  typedef typename std::vector<T>::iterator iterator;
  citerator begin() const { return bin_edges_.cbegin(); }///< iterator for external use
  citerator end() const { return bin_edges_.cend(); }    ///< iterator for external use
  iterator begin() { return bin_edges_.begin(); }        ///< iterator for external use
  iterator end() { return bin_edges_.end(); }            ///< iterator for external use
  /**
   * Set Name of axis.
   * @param name name of axis
   */
  void SetName(const std::string name) { name_ = name; }
  /**
   * Returns Name of axis.
   * @return name of axis
   */
  std::string Name() const { return name_; }

  /**
   * Returns a short name for the axis.
   * @return short name
   */
  std::string ShortName() const {
    auto is_vowel = [](const char p_char) {
      constexpr char vowels[] = {'a', 'e', 'i', 'o', 'u', 'A', 'E', 'I', 'O', 'U'};
      return std::find(std::begin(vowels), std::end(vowels), p_char) != std::end(vowels);
    };
    std::string t_name = name_;
    t_name.erase(std::remove_if(t_name.begin(), t_name.end(), is_vowel), t_name.end());
    std::for_each(t_name.begin(), t_name.end(), [](char &c) { c = ::tolower(c); });
    std::ostringstream axislimits;
    axislimits.precision(2);
    axislimits << std::scientific;
    axislimits << "_" << t_name << "_" << size() << "_" << GetFirstBinEdge() << "_" << GetLastBinEdge() << "_";
    return axislimits.str();
  }

  typename std::vector<T>::size_type GetNBins() const { return bin_edges_.size() - 1; }
  const T *GetPtr() const { return bin_edges_.data(); }

  /**
   * Finds bin index for a given value
   * if value is smaller than lowest bin return -1.
   * @param value for finding corresponding bin
   * @return bin index
   */
  inline long FindBin(const T value) const {
    long bin = 0;
    if (value < *bin_edges_.begin()) {
      bin = -1;
    } else {
      auto lb = std::lower_bound(bin_edges_.begin(), bin_edges_.end(), value);
      if (lb == bin_edges_.end()) {
        bin = -1;
      } else if (lb == bin_edges_.begin() || *lb == value)
        bin = (lb - bin_edges_.begin());
      else
        bin = (lb - bin_edges_.begin()) - 1;
    }
    if (bin >= (long) bin_edges_.size() - 1 || bin < 0)
      bin = -1;
    return bin;
  };

  std::string GetBinName(unsigned int i) const {
    auto lower = std::to_string(GetLowerBinEdge(i));
    lower = lower.erase(lower.find_last_not_of('0') + 1, std::string::npos);
    lower = lower.erase(lower.find_last_not_of('.') + 1, std::string::npos);
    return name_ + ":" + lower;
  }
  /**
   * Returns number of bins.
   * @return number of bins.
   */
  constexpr typename std::vector<T>::size_type size() const { return bin_edges_.size() - 1; }
  /**
   * Gets lower bin edge
   * @param bin Index of bin of interest
   * @return lower edge of bin of interest
   */
  inline T GetLowerBinEdge(const unsigned long bin) const { return bin_edges_.at(bin); }
  /**
   * Gets upper bin edge
   * @param bin Index of bin of interest
   * @return upper edge of bin of interest
   */
  inline T GetUpperBinEdge(const unsigned long bin) const { return bin_edges_.at(bin + 1); }
  /**
   * Gets bin center
   * @param bin Index of bin of interest
   * @return upper edge of bin of interest
   */
  inline T GetBinCenter(const unsigned long bin) const { return (bin_edges_.at(bin + 1) + bin_edges_.at(bin)) / 2.; }
  /**
   * Gets lower bin edge
   * @param bin Index of bin of interest
   * @return lower edge of bin of interest
   */
  inline T GetFirstBinEdge() const { return bin_edges_.front(); }
  /**
   * Gets upper bin edge
   * @param bin Index of bin of interest
   * @return upper edge of bin of interest
   */
  inline T GetLastBinEdge() const { return bin_edges_.back(); }

  void Print() const {
    std::cout << "OBJ: Qn::Axis " << name_ << "\n";
    std::cout << "number of Bins:" << GetNBins() << "\n";
    std::cout << "bin edges: ";
    for (unsigned int i = 0; i < bin_edges_.size() - 1; ++i) {
      std::cout << bin_edges_[i] << ", ";
    }
    std::cout << bin_edges_.back() << "\n";
  }

 private:
  std::string name_;
  std::vector<T> bin_edges_;

  /// \cond CLASSIMP
  ClassDef(Axis, 5);
  /// \endcond
};
using AxisF = Axis<float>;
using AxisD = Axis<double>;
}// namespace Qn
#endif//FLOW_QNAXIS_H
