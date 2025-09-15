#ifndef VEXCL_UTIL_HPP
#define VEXCL_UTIL_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   vexcl/util.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL general utilities.
 */

#if defined(_MSC_VER)
#  if defined(min) || defined(max)
#    error Please define NOMINMAX macro globally in your project
#  endif
#  if defined(_VARIADIC_MAX) && (_VARIADIC_MAX < 10)
#    error Please define _VARIADIC_MAX=10 or greater in your project
#  endif
#endif

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <type_traits>

#include <boost/config.hpp>

#ifdef BOOST_NO_VARIADIC_TEMPLATES
#  include <boost/proto/proto.hpp>
#  include <boost/preprocessor/repetition.hpp>
#  ifndef VEXCL_MAX_ARITY
#    define VEXCL_MAX_ARITY BOOST_PROTO_MAX_ARITY
#  endif
#endif

namespace vex {

/// Check run-time condition.
/** Throws std::runtime_error if condition is false */
template <class Condition, class Message>
inline void precondition(const Condition &condition, const Message &fail_message) {
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4800)
#endif
    if (!condition) throw std::runtime_error(fail_message);
#ifdef _MSC_VER
#  pragma warning(pop)
#endif
}

/// Return next power of 2.
inline size_t nextpow2(size_t x) {
    --x;
    x |= x >> 1U;
    x |= x >> 2U;
    x |= x >> 4U;
    x |= x >> 8U;
    x |= x >> 16U;
    return ++x;
}

/// Align n to the next multiple of m.
inline size_t alignup(size_t n, size_t m = 16U) {
    return (n + m - 1) / m * m;
}

template <class T>
struct is_tuple : std::false_type {};


#ifndef BOOST_NO_VARIADIC_TEMPLATES

template <class... Elem>
struct is_tuple < std::tuple<Elem...> > : std::true_type {};

#else

#define VEXCL_IS_TUPLE(z, n, unused)                                           \
  template <BOOST_PP_ENUM_PARAMS(n, class Elem)>                               \
  struct is_tuple<                                                             \
      std::tuple<BOOST_PP_ENUM_PARAMS(n, Elem)> > : std::true_type {           \
  };

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_IS_TUPLE, ~)

#undef VEXCL_IS_TUPLE

#endif

#ifndef BOOST_NO_VARIADIC_TEMPLATES
/// Create std::array from arguments
template <class T, class... Tail>
std::array<T, 1 + sizeof...(Tail)>
make_array(T t, Tail... tail) {
    std::array<T, 1 + sizeof...(Tail)> a = {{t, static_cast<T>(tail)...}};
    return a;
}
#else

#define VEXCL_INIT_ARRAY(z, n, data) static_cast<T0>(t ## n)
#define VEXCL_MAKE_ARRAY(z, n, data)                                           \
  template <BOOST_PP_ENUM_PARAMS(n, class T)>                                  \
  std::array<T0, n> make_array(BOOST_PP_ENUM_BINARY_PARAMS(n, T, t)) {         \
    std::array<T0, n> a = { { BOOST_PP_ENUM(n, VEXCL_INIT_ARRAY, ~) } };       \
    return a;                                                                  \
  }

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_MAKE_ARRAY, ~)

#undef VEXCL_MAKE_ARRAY
#undef VEXCL_INIT_ARRAY

#endif

struct column_owner {
    const std::vector<size_t> &part;

    column_owner(const std::vector<size_t> &part) : part(part) {}

    size_t operator()(size_t c) const {
        return std::upper_bound(part.begin(), part.end(), c)
            - part.begin() - 1;
    }
};

inline const char* getenv(const char *name, const char *defval = NULL) {
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4996)
#endif

    const char *val = ::getenv(name);

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

    return val ? val : (defval ? defval : NULL);
}

} // namespace vex

#endif
