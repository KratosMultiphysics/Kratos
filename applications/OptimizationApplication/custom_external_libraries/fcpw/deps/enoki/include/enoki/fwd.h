/*
    enoki/fwd.h -- Preprocessor definitions and forward declarations

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#if defined(_MSC_VER)
#  if !defined(_USE_MATH_DEFINES)
#    define _USE_MATH_DEFINES
#  endif
#endif

#include <cstddef>
#include <cstring>
#include <type_traits>

#if defined(_MSC_VER)
#  define ENOKI_NOINLINE               __declspec(noinline)
#  define ENOKI_INLINE                 __forceinline
#  define ENOKI_INLINE_LAMBDA
#  define ENOKI_PURE
#  define ENOKI_MALLOC                 __declspec(restrict)
#  define ENOKI_MAY_ALIAS
#  define ENOKI_ASSUME_ALIGNED(x, s)   x
#  define ENOKI_UNROLL
#  define ENOKI_NOUNROLL
#  define ENOKI_IVDEP                  __pragma(loop(ivdep))
#  define ENOKI_PACK
#  define ENOKI_LIKELY(x)              x
#  define ENOKI_UNLIKELY(x)            x
#  define ENOKI_REGCALL
#  define ENOKI_IMPORT                 __declspec(dllimport)
#  define ENOKI_EXPORT                 __declspec(dllexport)
#else
#  define ENOKI_NOINLINE               __attribute__ ((noinline))
#  define ENOKI_INLINE                 __attribute__ ((always_inline)) inline
#  define ENOKI_INLINE_LAMBDA          __attribute__ ((always_inline))
#  define ENOKI_PURE                   __attribute__ ((const,nothrow))
#  define ENOKI_MALLOC                 __attribute__ ((malloc))
#  define ENOKI_ASSUME_ALIGNED(x, s)   __builtin_assume_aligned(x, s)
#  define ENOKI_LIKELY(x)              __builtin_expect(!!(x), 1)
#  define ENOKI_UNLIKELY(x)            __builtin_expect(!!(x), 0)
#  define ENOKI_PACK                   __attribute__ ((packed))
#  if defined(__clang__)
#    define ENOKI_UNROLL               _Pragma("unroll")
#    define ENOKI_NOUNROLL             _Pragma("nounroll")
#    define ENOKI_IVDEP
#    define ENOKI_MAY_ALIAS            __attribute__ ((may_alias))
#    define ENOKI_REGCALL              __attribute__ ((regcall))
#  elif defined(__INTEL_COMPILER)
#    define ENOKI_MAY_ALIAS
#    define ENOKI_UNROLL               _Pragma("unroll")
#    define ENOKI_NOUNROLL             _Pragma("nounroll")
#    define ENOKI_IVDEP                _Pragma("ivdep")
#    define ENOKI_REGCALL              __attribute__ ((regcall))
#  else
#    define ENOKI_MAY_ALIAS            __attribute__ ((may_alias))
#    define ENOKI_UNROLL
#    define ENOKI_NOUNROLL
#    if defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9))
#      define ENOKI_IVDEP              _Pragma("GCC ivdep")
#    else
#      define ENOKI_IVDEP
#    endif
#    define ENOKI_REGCALL
#  endif
#  define ENOKI_IMPORT
#  define ENOKI_EXPORT                 __attribute__ ((visibility("default")))
#endif

#define ENOKI_MARK_USED(x) (void) x

#if !defined(NAMESPACE_BEGIN)
#  define NAMESPACE_BEGIN(name) namespace name {
#endif

#if !defined(NAMESPACE_END)
#  define NAMESPACE_END(name) }
#endif

#define ENOKI_VERSION_MAJOR 0
#define ENOKI_VERSION_MINOR 1
#define ENOKI_VERSION_PATCH 0

#define ENOKI_STRINGIFY(x) #x
#define ENOKI_TOSTRING(x)  ENOKI_STRINGIFY(x)
#define ENOKI_VERSION                                                          \
    (ENOKI_TOSTRING(ENOKI_VERSION_MAJOR) "."                                   \
     ENOKI_TOSTRING(ENOKI_VERSION_MINOR) "."                                   \
     ENOKI_TOSTRING(ENOKI_VERSION_PATCH))

#if defined(__clang__) && defined(__apple_build_version__)
#  if __clang_major__ < 10
#    error Enoki requires a very recent version of AppleClang (XCode >= 10.0)
#  endif
#elif defined(__clang__)
#  if __clang_major__ < 7 && !defined(EMSCRIPTEN)
#    error Enoki requires a very recent version of Clang/LLVM (>= 7.0)
#  endif
#elif defined(__GNUC__)
#  if (__GNUC__ < 8) || (__GNUC__ == 8 && __GNUC_MINOR__ < 2)
#    error Enoki requires a very recent version of GCC (>= 8.2)
#  endif
#endif

#if defined(__x86_64__) || defined(_M_X64)
#  define ENOKI_X86_64 1
#endif

#if (defined(__i386__) || defined(_M_IX86)) && !defined(ENOKI_X86_64)
#  define ENOKI_X86_32 1
#endif

#if defined(__aarch64__)
#  define ENOKI_ARM_64 1
#elif defined(__arm__)
#  define ENOKI_ARM_32 1
#endif

#if (defined(_MSC_VER) && defined(ENOKI_X86_32)) && !defined(ENOKI_DISABLE_VECTORIZATION)
// Enoki does not support vectorization on 32-bit Windows due to various
// platform limitations (unaligned stack, calling conventions don't allow
// passing vector registers, etc.).
# define ENOKI_DISABLE_VECTORIZATION 1
#endif

# if !defined(ENOKI_DISABLE_VECTORIZATION)
#  if defined(__AVX512F__)
#    define ENOKI_X86_AVX512F 1
#  endif
#  if defined(__AVX512CD__)
#    define ENOKI_X86_AVX512CD 1
#  endif
#  if defined(__AVX512DQ__)
#    define ENOKI_X86_AVX512DQ 1
#  endif
#  if defined(__AVX512VL__)
#    define ENOKI_X86_AVX512VL 1
#  endif
#  if defined(__AVX512BW__)
#    define ENOKI_X86_AVX512BW 1
#  endif
#  if defined(__AVX512PF__)
#    define ENOKI_X86_AVX512PF 1
#  endif
#  if defined(__AVX512ER__)
#    define ENOKI_X86_AVX512ER 1
#  endif
#  if defined(__AVX512VBMI__)
#    define ENOKI_X86_AVX512VBMI 1
#  endif
#  if defined(__AVX512VPOPCNTDQ__)
#    define ENOKI_X86_AVX512VPOPCNTDQ 1
#  endif
#  if defined(__AVX2__)
#    define ENOKI_X86_AVX2 1
#  endif
#  if defined(__FMA__)
#    define ENOKI_X86_FMA 1
#  endif
#  if defined(__F16C__)
#    define ENOKI_X86_F16C 1
#  endif
#  if defined(__AVX__)
#    define ENOKI_X86_AVX 1
#  endif
#  if defined(__SSE4_2__)
#    define ENOKI_X86_SSE42 1
#  endif
#  if defined(__ARM_NEON)
#    define ENOKI_ARM_NEON
#  endif
#  if defined(__ARM_FEATURE_FMA)
#    define ENOKI_ARM_FMA
#  endif
#endif

/* Fix missing/inconsistent preprocessor flags */
#if defined(ENOKI_X86_AVX512F) && !defined(ENOKI_X86_AVX2)
#  define ENOKI_X86_AVX2
#endif

#if defined(ENOKI_X86_AVX2) && !defined(ENOKI_X86_F16C)
#  define ENOKI_X86_F16C
#endif

#if defined(ENOKI_X86_AVX2) && !defined(ENOKI_X86_FMA)
#  define ENOKI_X86_FMA
#endif

#if defined(ENOKI_X86_AVX2) && !defined(ENOKI_X86_AVX)
#  define ENOKI_X86_AVX
#endif

#if defined(ENOKI_X86_AVX) && !defined(ENOKI_X86_SSE42)
#  define ENOKI_X86_SSE42
#endif

/* The following macro is used by the test suite to detect
   unimplemented methods in vectorized backends */

#if !defined(ENOKI_TRACK_SCALAR)
#  define ENOKI_TRACK_SCALAR(reason)
#endif

#if defined(ENOKI_ALLOC_VERBOSE)
#  define ENOKI_TRACK_ALLOC(ptr, size)                                         \
      printf("Enoki: %p: alloc(%llu)\n", (ptr), (unsigned long long) (size));
#  define ENOKI_TRACK_DEALLOC(ptr, size)                                       \
      printf("Enoki: %p: dealloc(%llu)\n", (ptr), (unsigned long long) (size));
#endif

#if !defined(ENOKI_TRACK_ALLOC)
#  define ENOKI_TRACK_ALLOC(ptr, size)
#endif

#if !defined(ENOKI_TRACK_DEALLOC)
#  define ENOKI_TRACK_DEALLOC(ptr, size)
#endif

#define ENOKI_CHKSCALAR(reason)                                                \
    if (std::is_arithmetic_v<std::decay_t<Value>>) {                           \
        ENOKI_TRACK_SCALAR(reason)                                             \
    }

#if !defined(ENOKI_APPROX_DEFAULT)
#  define ENOKI_APPROX_DEFAULT 1
#endif

NAMESPACE_BEGIN(enoki)

using ssize_t = std::make_signed_t<size_t>;

/// Maximum hardware-supported packet size in bytes
#if defined(ENOKI_X86_AVX512F)
    static constexpr size_t max_packet_size = 64;
#elif defined(ENOKI_X86_AVX)
    static constexpr size_t max_packet_size = 32;
#elif defined(ENOKI_X86_SSE42) || defined(ENOKI_ARM_NEON)
    static constexpr size_t max_packet_size = 16;
#else
    static constexpr size_t max_packet_size = 4;
#endif

constexpr size_t array_default_size = max_packet_size / 4;

/// Base class of all arrays
template <typename Value_, typename Derived_> struct ArrayBase;

/// Base class of all statically sized arrays
template <typename Value_, size_t Size_, bool IsMask_, typename Derived_>
struct StaticArrayBase;

/// Generic array class, which broadcasts from the outer to inner dimensions
template <typename Value_, size_t Size_ = array_default_size>
struct Array;

/// Generic array class, which broadcasts from the inner to outer dimensions
template <typename Value_, size_t Size_ = array_default_size>
struct Packet;

/// Generic mask class, which broadcasts from the outer to inner dimensions
template <typename Value_, size_t Size_ = array_default_size>
struct Mask;

/// Generic mask class, which broadcasts from the inner to outer dimensions
template <typename Value_, size_t Size_ = array_default_size>
struct PacketMask;

/// Dynamically sized array
template <typename Packet_> struct DynamicArray;
template <typename Packet_> struct DynamicMask;

/// Reverse-mode autodiff array
template <typename Value> struct DiffArray;

template <typename Value_, size_t Size_>
struct Matrix;

template <typename Value_>
struct Complex;

template <typename Value_>
struct Quaternion;

/// Helper class for custom data structures
template <typename T, typename = int>
struct struct_support;

template <typename Value>
struct CUDAArray;

template <typename T> class cuda_host_allocator;
template <typename T> class cuda_managed_allocator;

extern ENOKI_IMPORT void* cuda_host_malloc(size_t);
extern ENOKI_IMPORT void cuda_host_free(void *);

/// Half-precision floating point value
struct half;

template <typename T> struct MaskBit;

namespace detail {
    struct reinterpret_flag { };
}

template <typename T, bool UseIntrinsic = false, typename = int>
struct divisor;
template <typename T>
struct divisor_ext;

/// Reinterpret the binary represesentation of a data type
template<typename T, typename U> ENOKI_INLINE T memcpy_cast(const U &val) {
    static_assert(sizeof(T) == sizeof(U), "memcpy_cast: sizes did not match!");
    T result;
    std::memcpy(&result, &val, sizeof(T));
    return result;
}

NAMESPACE_END(enoki)
