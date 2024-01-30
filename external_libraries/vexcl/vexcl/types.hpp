#ifndef VEXCL_TYPES_HPP
#define VEXCL_TYPES_HPP

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
 * \file   vexcl/types.hpp
 * \author Pascal Germroth <pascal@ensieve.org>
 * \brief  Support for using native C++ and OpenCL types in expressions.
 */

#include <string>
#include <type_traits>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/io/ios_state.hpp>

// This should be available for all backends
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl_platform.h>
#else
#include <CL/cl_platform.h>
#endif

typedef unsigned int  uint;
typedef unsigned char uchar;

namespace vex {

    /// Get the corresponding scalar type for a CL vector (or scalar) type.
    /** \code cl_scalar_of<cl_float4>::type == cl_float \endcode */
    template <class T>
    struct cl_scalar_of {};

    /// Get the corresponding vector type for a CL scalar type.
    /** \code cl_vector_of<cl_float, 4>::type == cl_float4 \endcode */
    template <class T, int dim>
    struct cl_vector_of {};

    /// Get the number of values in a CL vector (or scalar) type.
    /** \code cl_vector_length<cl_float4>::value == 4 \endcode */
    template <class T>
    struct cl_vector_length {};

} // namespace vex

#define VEXCL_BIN_OP(name, len, op)                                            \
  inline cl_##name##len &operator op## =(cl_##name##len & a,                   \
                                         const cl_##name##len & b) {           \
    for (size_t i = 0; i < len; i++) {                                         \
      a.s[i] op## = b.s[i];                                                    \
    }                                                                          \
    return a;                                                                  \
  }                                                                            \
  inline cl_##name##len operator op(const cl_##name##len & a,                  \
                                    const cl_##name##len & b) {                \
    cl_##name##len res = a;                                                    \
    return res op## = b;                                                       \
  }

// `scalar OP vector` acts like `(vector_t)(scalar) OP vector` in OpenCl:
// all components are set to the scalar value.
#define VEXCL_BIN_SCALAR_OP(name, len, op)                                     \
  inline cl_##name##len &operator op## =(cl_##name##len & a,                   \
                                         const cl_##name & b) {                \
    for (size_t i = 0; i < len; i++) {                                         \
      a.s[i] op## = b;                                                         \
    }                                                                          \
    return a;                                                                  \
  }                                                                            \
  inline cl_##name##len operator op(const cl_##name##len & a,                  \
                                    const cl_##name & b) {                     \
    cl_##name##len res = a;                                                    \
    return res op## = b;                                                       \
  }                                                                            \
  inline cl_##name##len operator op(const cl_##name & a,                       \
                                    const cl_##name##len & b) {                \
    cl_##name##len res = b;                                                    \
    return res op## = a;                                                       \
  }

#define VEXCL_VEC_TYPE(name, len)                                              \
  VEXCL_BIN_OP(name, len, +)                                                   \
  VEXCL_BIN_OP(name, len, -)                                                   \
  VEXCL_BIN_OP(name, len, *)                                                   \
  VEXCL_BIN_OP(name, len, /)                                                   \
  VEXCL_BIN_SCALAR_OP(name, len, +)                                            \
  VEXCL_BIN_SCALAR_OP(name, len, -)                                            \
  VEXCL_BIN_SCALAR_OP(name, len, *)                                            \
  VEXCL_BIN_SCALAR_OP(name, len, /)                                            \
  inline cl_##name##len operator-(const cl_##name##len & a) {                  \
    cl_##name##len res;                                                        \
    for (size_t i = 0; i < len; i++) {                                         \
      res.s[i] = -a.s[i];                                                      \
    }                                                                          \
    return res;                                                                \
  }                                                                            \
  inline std::ostream &operator<<(std::ostream & os,                           \
                                  const cl_##name##len & value) {              \
    boost::io::ios_all_saver stream_state(os);                                 \
    os << "{";                                                                 \
    for (std::size_t i = 0; i < len; i++) {                                    \
      if (i != 0) { os << ','; }                                               \
      os << std::setw(13) << std::scientific << value.s[i];                    \
    }                                                                          \
    return os << '}';                                                          \
  }                                                                            \
  namespace vex {                                                              \
  template <> struct cl_scalar_of<cl_##name##len> {                            \
    typedef cl_##name type;                                                    \
  };                                                                           \
  template <> struct cl_vector_of<cl_##name, len> {                            \
    typedef cl_##name##len type;                                               \
  };                                                                           \
  template <>                                                                  \
  struct cl_vector_length<cl_##name##len>                                      \
      : std::integral_constant<unsigned, len> { };                             \
  }

#define VEXCL_TYPES(name)                                                      \
  VEXCL_VEC_TYPE(name, 2)                                                      \
  VEXCL_VEC_TYPE(name, 4)                                                      \
  VEXCL_VEC_TYPE(name, 8)                                                      \
  VEXCL_VEC_TYPE(name, 16)                                                     \
  namespace vex {                                                              \
  template <> struct cl_scalar_of<cl_##name> {                                 \
    typedef cl_##name type;                                                    \
  };                                                                           \
  template <> struct cl_vector_of<cl_##name, 1> {                              \
    typedef cl_##name type;                                                    \
  };                                                                           \
  template <>                                                                  \
  struct cl_vector_length<cl_##name> : std::integral_constant<unsigned, 1> {}; \
  }

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4146)
#endif
VEXCL_TYPES(float)
VEXCL_TYPES(double)
VEXCL_TYPES(char)
VEXCL_TYPES(uchar)
VEXCL_TYPES(short)
VEXCL_TYPES(ushort)
VEXCL_TYPES(int)
VEXCL_TYPES(uint)
VEXCL_TYPES(long)
VEXCL_TYPES(ulong)
#ifdef _MSC_VER
#  pragma warning(pop)
#endif


#undef VEXCL_BIN_OP
#undef VEXCL_BIN_SCALAR_OP
#undef VEXCL_VEC_TYPE
#undef VEXCL_TYPES


namespace vex {

/// Convert each element of the vector to another type.
template<class To, class From>
inline To cl_convert(const From &val) {
    const size_t n = cl_vector_length<To>::value;
    static_assert(n == cl_vector_length<From>::value, "Vectors must be same length.");
    To out;
    for(size_t i = 0 ; i != n ; i++)
        out.s[i] = val.s[i];
    return out;
}

/// Declares a type as CL native, allows using it as a literal.
template <class T> struct is_cl_native : std::false_type {};

/// Convert typename to string.
template <class T, class Enable = void>
struct type_name_impl;

template <class T>
inline std::string type_name() {
    return type_name_impl<T>::get();
}

template<typename T>
struct type_name_impl<T&>
{
    static std::string get() {
        return type_name_impl<T>::get() + " &";
    }
};

#define VEXCL_STRINGIFY(name)                                                  \
  template<> struct type_name_impl<cl_##name> {                                \
    static std::string get() { return #name; }                                 \
  };

#define VEXCL_NATIVE(name)                                                     \
  template<> struct is_cl_native<cl_##name> : std::true_type { };

// enable use of OpenCL vector types as literals
#define VEXCL_VEC_TYPE(name, len)                                              \
  template <> struct type_name_impl<cl_##name##len> {                          \
    static std::string get() { return #name #len; }                            \
  };                                                                           \
  template <> struct is_cl_native<cl_##name##len> : std::true_type { };

#define VEXCL_TYPES(name)                                                      \
  VEXCL_STRINGIFY(name)                                                        \
  VEXCL_NATIVE(name)                                                           \
  VEXCL_VEC_TYPE(name, 2)                                                      \
  VEXCL_VEC_TYPE(name, 4)                                                      \
  VEXCL_VEC_TYPE(name, 8)                                                      \
  VEXCL_VEC_TYPE(name, 16)

VEXCL_TYPES(float)
VEXCL_TYPES(double)
VEXCL_TYPES(char)
VEXCL_TYPES(uchar)
VEXCL_TYPES(short)
VEXCL_TYPES(ushort)
VEXCL_TYPES(int)
VEXCL_TYPES(uint)
VEXCL_TYPES(long)
VEXCL_TYPES(ulong)

#undef VEXCL_TYPES
#undef VEXCL_VEC_TYPE
#undef VEXCL_STRINGIFY

// One can not pass these to the kernel, but the overloads are useful
template <> struct type_name_impl<bool> {
    static std::string get() { return "bool"; }
};
template <> struct type_name_impl<void> {
    static std::string get() { return "void"; }
};

// char and cl_char are different types. Hence, special handling is required:
template <> struct type_name_impl<char> {
    static std::string get() { return "char"; }
};
template <> struct is_cl_native<char> : std::true_type {};
template <> struct cl_vector_length<char> : std::integral_constant<unsigned, 1> {};
template <> struct cl_scalar_of<char> { typedef char type; };

#if defined(__APPLE__)
template <> struct type_name_impl<size_t>
    : public type_name_impl<
        boost::mpl::if_c<
            sizeof(std::size_t) == sizeof(uint),
            cl_uint, cl_ulong
        >::type
    >
{};

template <> struct type_name_impl<ptrdiff_t>
    : public type_name_impl<
        boost::mpl::if_c<
            sizeof(std::size_t) == sizeof(uint),
            cl_int, cl_long
        >::type
    >
{};

template <> struct is_cl_native<size_t>    : std::true_type {};
template <> struct is_cl_native<ptrdiff_t> : std::true_type {};

template <> struct cl_vector_length<size_t>    : std::integral_constant<unsigned, 1> {};
template <> struct cl_vector_length<ptrdiff_t> : std::integral_constant<unsigned, 1> {};

template <> struct cl_scalar_of<size_t>       { typedef size_t    type; };
template <> struct cl_vector_of<size_t, 1>    { typedef size_t    type; };
template <> struct cl_scalar_of<ptrdiff_t>    { typedef ptrdiff_t type; };
template <> struct cl_vector_of<ptrdiff_t, 1> { typedef ptrdiff_t type; };
#endif

template <class T, class Enable = void>
struct is_cl_scalar : std::false_type {};

template <class T>
struct is_cl_scalar<T,
    typename std::enable_if<is_cl_native<T>::value && (cl_vector_length<T>::value == 1)>::type
    > : std::true_type
{};

template <class T, class Enable = void>
struct is_cl_vector : std::false_type {};

template <class T>
struct is_cl_vector<T,
    typename std::enable_if<is_cl_native<T>::value && (cl_vector_length<T>::value > 1)>::type
    > : std::true_type
{};

template <unsigned I, class Enable = void>
struct cl_fit_vec_size { };

template <> struct cl_fit_vec_size<1> : std::integral_constant<int, 1> {};
template <> struct cl_fit_vec_size<2> : std::integral_constant<int, 2> {};

template <unsigned I>
struct cl_fit_vec_size<I, typename std::enable_if<(I>2) && (I<=4)>::type>
: std::integral_constant<unsigned, 4> {};

template <unsigned I>
struct cl_fit_vec_size<I, typename std::enable_if<(I>4) && (I<=8)>::type>
: std::integral_constant<unsigned, 8> {};

template <unsigned I>
struct cl_fit_vec_size<I, typename std::enable_if<(I>8) && (I<=16)>::type>
: std::integral_constant<unsigned, 16> {};

template <unsigned I>
struct cl_fit_vec_size<I, typename std::enable_if<(I<1) || (I>16)>::type>
{
    static_assert(I >= 1 && I <= 16, "Unsupported vector size");
};

}

#endif
