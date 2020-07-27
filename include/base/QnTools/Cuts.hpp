#include <utility>

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

#ifndef QNCUTS_H
#define QNCUTS_H

#include "ROOT/RIntegerSequence.hxx"
#include "ROOT/RMakeUnique.hxx"
#include <string>

namespace Qn {

/**
 * Base class of the Cut.
 */
struct ICut {
  virtual ~ICut() = default;
  virtual bool Check() const = 0;
  virtual bool Check(unsigned int) const = 0;
  virtual std::string Name() const = 0;
};

/**
 * Template class of a cut applied in the correction step.
 * Number of dimensions is determined by the signature of the cut function
 * and by the size of the variables array passed to the constructor.
 * @tparam T Type of variable
 */
template<typename VAR, typename... T>
class Cut : public ICut {
 public:
  Cut(VAR (&arr)[sizeof...(T)], std::function<bool(T...)> lambda, std::string name)
      : lambda_(lambda), name_(std::move(name)) {
    int i = 0;
    for (auto &a : arr) {
      variables_.at(i) = std::move(a);
      ++i;
    }
  }
  bool Check(const unsigned int i_channel) const override {
    return CheckImpl(i_channel, std::make_index_sequence<sizeof...(T)>{});
  }
  bool Check() const override { return CheckImpl(0, std::make_index_sequence<sizeof...(T)>{}); }

  std::string Name() const override { return name_; }

 private:
  /**
   * Implements the evaluation of the cut.
   * @tparam I index sequence
   * @param i offset from the variable id.
   * @return true if the cut is passed.
   */
  template<std::size_t... I>
  bool CheckImpl(const unsigned int i, std::index_sequence<I...>) const {
    return lambda_(*(variables_[I].Get() + i)...);
  }
  std::array<VAR, sizeof...(T)> variables_;/// array of the variables used in the cut.
  std::function<bool(T...)> lambda_;       /// function used to evaluate the cut.
  std::string name_;
};

namespace Details {
template<typename T, std::size_t>
using CutDataType = T &;
template<typename T, typename VAR, std::size_t N, typename FUNC, std::size_t... Is>
std::unique_ptr<Cut<VAR, CutDataType<T, Is>...>> MakeUniqueCutImpl(std::index_sequence<Is...>,
                                                                   VAR (&arr)[N],
                                                                   FUNC &&func, std::string name) {
  return std::make_unique<Cut<VAR, CutDataType<T, Is>...>>(arr, std::forward<FUNC>(func), name);
}
}// namespace Details

/**
 * Function which create a unique_ptr of a cut of n variables.
 * @tparam N length of the variable array.
 * @tparam FUNC type of the cut function
 * @param arr array of variables
 * @param func cut function
 * @return Returns a unique pointer to the cut.
 */
template<typename T, std::size_t N, typename FUNC, typename VAR>
std::unique_ptr<ICut> MakeUniqueCut(VAR (&arr)[N], FUNC &&func, std::string name) {
  return Details::MakeUniqueCutImpl<T>(std::make_index_sequence<N>{}, arr, std::forward<FUNC>(func), name);
}
}// namespace Qn

#endif