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
 * Interface class of the Cut.
 */
struct ICut {
  virtual ~ICut() = default;
  virtual bool Check() const = 0;
  virtual bool Check(unsigned int) const = 0;
  virtual std::string Name() const = 0;
};

namespace Details {

template <typename T>
using remove_cvref_t = typename std::remove_cv<std::remove_reference_t<T>>::type;

template <typename Function, typename ... Arg>
constexpr bool is_static_cut_v = std::is_same_v<bool, std::invoke_result_t<Function&&, Arg...>>;

template <typename Function, typename T>
constexpr bool is_dynamic_cut_v = std::is_same_v<bool, std::invoke_result_t<Function&&,const std::vector<remove_cvref_t<T>>>>;

template <typename Variable>
using variable_get_t = std::remove_pointer_t<decltype(std::declval<const Variable>().Get())>;

template<typename T, std::size_t>
using enumerated_arg = T &;


}// namespace Details

/**
 * Template class of a cut applied in the correction step.
 * Number of dimensions is determined by the signature of the cut function
 * and by the size of the variables array passed to the constructor.
 * @tparam T Type of variable
 */
template<typename Variable, typename... T>
class StaticCut : public ICut {
 public:
  StaticCut() = delete;
  StaticCut(Variable (&arr)[sizeof...(T)], std::function<bool(T...)> lambda, std::string name)
      : lambda_(lambda), name_(std::move(name)) {
    int i = 0;
    for (auto &a : arr) {
      variables_.at(i) = std::move(a);
      ++i;
    }
  }

  StaticCut(
      const std::array<Variable, sizeof...(T)> &args,
      std::function<bool (T...)> lambda,
      std::string name
      ) : variables_(args), lambda_(lambda), name_(std::move(name)) {}


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

  std::array<Variable, sizeof...(T)> variables_;/// array of the variables used in the cut.
  std::function<bool(T...)> lambda_;       /// function used to evaluate the cut.
  std::string name_;
};


template <typename Variable>
class DynamicCut : public ICut {
 public:
  typedef typename Details::remove_cvref_t<Details::variable_get_t<Variable>> lambda_arg_ele_t;
  typedef typename std::vector<lambda_arg_ele_t> lambda_arg_t;
  typedef typename std::function<bool (const lambda_arg_t&)> lambda_t;

  DynamicCut(const std::vector<Variable> &variables,
             lambda_t lambda,
             std::string name)
      : variables_(variables), lambda_(lambda), name_(std::move(name)) {
  }
  bool Check() const override { return CheckImpl(0); }
  bool Check(unsigned int i) const override { return CheckImpl(i); }
  std::string Name() const override { return name_; }

 private:
  bool CheckImpl(unsigned int channel_id) const {
    lambda_arg_t lambda_args(variables_.size());
    std::transform(std::begin(variables_), std::end(variables_),
                     std::begin(lambda_args), [this,channel_id] (const Variable& v) {
      return *(v.Get() + channel_id);
    });

    return lambda_(lambda_args);
  }


  std::vector<Variable> variables_;
  lambda_t lambda_;
  std::string name_;
};

namespace Details {

template<typename Variable, std::size_t N, typename Function, std::size_t... Is>
typename std::enable_if_t<is_static_cut_v<Function, enumerated_arg<variable_get_t<Variable>, Is>...>,std::unique_ptr<ICut>>
MakeStaticCutImpl(std::index_sequence<Is...>, const std::array<Variable, N> &arr, Function &&func, std::string name) {
  return std::make_unique<StaticCut<Variable, enumerated_arg<variable_get_t<Variable>, Is>...>>(arr, std::forward<Function>(func), name);
}


template<typename Variable, typename Function>
typename std::enable_if_t<is_dynamic_cut_v<Function, variable_get_t<Variable>>,std::unique_ptr<ICut>>
MakeDynamicCutImpl(const std::vector<Variable> &arr, Function &&func, std::string name) {
  return std::make_unique<DynamicCut<Variable>>(arr, func, name);
}


}

/**
 * Function which create a unique_ptr of a cut of n variables.
 * @tparam N length of the variable array.
 * @tparam Function type of the cut function
 * @param arr array of variables
 * @param func cut function
 * @return Returns a unique pointer to the cut.
 */
template<std::size_t N, typename Function, typename Variable>
std::unique_ptr<ICut> MakeUniqueCut(const std::array<Variable, N> &arr, Function &&func, std::string name) {
  return Details::MakeStaticCutImpl(std::make_index_sequence<N>{}, arr, std::forward<Function>(func), name);
}

template <typename Function, typename Variable>
std::unique_ptr<ICut> MakeUniqueCut(const std::vector<Variable> &arr, Function &&func, std::string name) {
  return Details::MakeDynamicCutImpl(arr, std::forward<Function>(func), name);
}


}// namespace Qn



#endif