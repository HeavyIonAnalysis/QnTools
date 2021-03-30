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
#ifndef TEMPLATEHELPERS_H
#define TEMPLATEHELPERS_H

#include <algorithm>
#include <tuple>
#include <utility>

#include <ROOT/RResultPtr.hxx>

namespace Qn::TemplateFunctions {

namespace Impl {
template <typename T>
std::enable_if_t<std::is_pointer<T>::value, std::remove_pointer_t<T> &>
Dereference(T &t) {
  return *t;
}
template <typename T>
std::enable_if_t<!std::is_pointer<T>::value, T &> Dereference(T &t) {
  return t;
}
template <typename Function, typename Tuple, size_t... I>
auto Call(Function f, Tuple t, std::index_sequence<I...>) {
  return f(Dereference(std::get<I>(t))...);
}
}  // namespace Impl
template <typename Function, typename Tuple>
auto Call(Function f, Tuple t) {
  static constexpr auto size = std::tuple_size<Tuple>::value;
  return Impl::Call(f, t, std::make_index_sequence<size>{});
}

namespace Impl {
template <std::size_t Index, typename Function, typename... Tuples>
constexpr void invoke_at(Function &&func, Tuples &&... tuples) {
  func(std::get<Index>(std::forward<Tuples>(tuples))...);
}

template <std::size_t... Indices, typename Function, typename... Tuples>
constexpr void apply_sequence(Function &&func, std::index_sequence<Indices...>,
                              Tuples &&... tuples) {
  (((void)invoke_at<Indices>(std::forward<Function>(func),
                             std::forward<Tuples>(tuples)...),
    ...));
}
}  // namespace Impl

template <typename Function, typename... Tuples>
constexpr void TupleForEach(Function &&func, Tuples &&... tuples) {
  static_assert(sizeof...(tuples) > 0,
                "Must be called with at least one tuple argument");

  constexpr auto min_length =
      std::min({std::tuple_size_v<std::remove_reference_t<Tuples>>...});
  if constexpr (min_length != 0) {
    Impl::apply_sequence(std::forward<Function>(func),
                         std::make_index_sequence<min_length>{},
                         std::forward<Tuples>(tuples)...);
  } else {
    func();
  }
}

namespace Impl {
template <size_t I, typename T>
struct TupleN {
  template <typename... Args>
  using type = typename TupleN<I - 1, T>::template type<T, Args...>;
};
template <typename T>
struct TupleN<0, T> {
  template <typename... Args>
  using type = std::tuple<Args...>;
};
}  // namespace Impl

template <size_t I, typename T>
using TupleOf = typename Impl::TupleN<I, T>::template type<>;

template <typename T>
struct FunctionTraits : public FunctionTraits<decltype(&T::operator())> {};
// For generic types, directly use the result of the signature of its
// 'operator()'
template <typename ClassType, typename ReturnType, typename... Args>
struct FunctionTraits<ReturnType (ClassType::*)(Args...) const> {
  // we specialize for pointers to member function
  enum { Arity = sizeof...(Args) };
  // arity is the number of arguments.
  typedef ReturnType result_type;

  template <size_t i>
  struct Arg {
    typedef
        typename std::tuple_element<i, std::tuple<Args..., void>>::type type;
    // the i-th argument is equivalent to the i-th tuple element of a tuple
    // composed of those arguments.
  };
};

namespace Impl {

template <typename check_type, typename first_type, typename... more_types>
struct IsEqualToCheckType {
  static const bool value =
      std::is_same<check_type, first_type>::value &&
      IsEqualToCheckType<check_type, more_types...>::value;
};

template <typename check_type, typename first_type>
struct IsEqualToCheckType<check_type, first_type>
    : std::is_same<check_type, first_type> {};

/**
 * Implements construction of the vector from a tuple
 * @tparam first_type
 * @tparam tuple_type
 * @tparam index
 * @param t
 * @return
 */
template <typename first_type, typename tuple_type, size_t... index>
auto ToVectorHelper(const tuple_type &t, std::index_sequence<index...>) {
  return std::vector<first_type>{std::get<index>(t)...};
}
}  // namespace Impl
/**
 * Constructs a vector from a tuple of all of the same types
 * @tparam first_type
 * @tparam others
 * @param t
 * @return
 */
template <typename first_type, typename... others>
auto ToVector(const std::tuple<first_type, others...> &t) {
  static_assert(
      Impl::IsEqualToCheckType<first_type, first_type, others...>::value,
      "types need to be equal");
  typedef typename std::remove_reference<decltype(t)>::type tuple_type;
  constexpr auto s = std::tuple_size<tuple_type>::value;
  return Impl::ToVectorHelper<first_type, tuple_type>(
      t, std::make_index_sequence<s>{});
}

template <typename T>
struct IsRResultPtr : std::false_type {};
template <typename T>
struct IsRResultPtr<ROOT::RDF::RResultPtr<T>> : std::true_type {};

template <typename T>
std::enable_if_t<IsRResultPtr<T>::value, typename T::Value_t>
DereferenceRResultPtr(T t) {
  return *t;
}

template <typename T>
std::enable_if_t<!IsRResultPtr<T>::value, T> DereferenceRResultPtr(T t) {
  return t;
}

}  // namespace Qn::TemplateHelpers
#endif  // QNTOOLS_TEMPLATEHELPERS_H
