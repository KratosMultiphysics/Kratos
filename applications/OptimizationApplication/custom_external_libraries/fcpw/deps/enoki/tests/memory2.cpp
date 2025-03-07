/*
    tests/memory2.cpp -- tests for more complex memory operations

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/matrix.h>
#include <enoki/dynamic.h>

#if defined(_MSC_VER)
#  define NOMINMAX
#  include <windows.h>
// MSVC sometimes reorders scatter/gather operations and plain array reads/writes.
// This macro is necessary to stop the behavior.
#  define MEMORY_BARRIER() MemoryBarrier()
#else
#  define MEMORY_BARRIER()
#endif

ENOKI_TEST_ALL(test01_extract) {
    auto idx = arange<T>();
    for (size_t i = 0; i < Size; ++i)
        assert(extract(idx, eq(idx, Value(i))) == Value(i));
}

ENOKI_TEST_ALL(test02_compress) {
    alignas(alignof(T)) Value tmp[T::ActualSize];
    auto value = arange<T>();
    Value *tmp2 = tmp;
    compress(tmp2, value, value >= Value(2));
    for (int i = 0; i < int(Size) - 2; ++i)
        assert(tmp[i] == Value(2 + i));
    assert(int(tmp2 - tmp) == std::max(0, int(Size) - 2));
}

ENOKI_TEST_ALL(test03_transform) {
    Value tmp[T::ActualSize] = { 0 };
    auto index = arange<uint_array_t<T>>();
    auto index2 = uint_array_t<T>(0u);

    transform<T>(tmp, index,
                 [](auto &value, auto & /*mask*/) { value += Value(1); });

    transform<T>(tmp, index,
                 [](auto &value, auto & /* mask */) { value += Value(1); },
                 mask_t<T>(false));

    transform<T>(
        tmp, index2,
        [](auto &value, auto &amount, auto & /* mask */) { value += amount; },
        T(2));

    transform<T>(
        tmp, index2,
        [](auto &value, auto &amount, auto & /* mask */) { value += amount; },
        T(2), mask_t<T>(false));

    assert(tmp[0] == 2*Size + 1);
    for (size_t i = 1; i < Size; ++i) {
        assert(tmp[i] == 1);
    }
}

#if defined(_MSC_VER)
  template <typename T, std::enable_if_t<T::Size == 8, int> = 0>
#else
  template <typename T, std::enable_if_t<T::Size != 31, int> = 0>
#endif
void test04_nested_gather_packed_impl() {
    using Value = value_t<T>;
    using UInt32P = Packet<uint32_t, T::Size>;
    using Vector4 = Array<Value, 4>;
    using Matrix4 = Matrix<Value, 4>;
    using Matrix4P = Matrix<T, 4>;
    using Array44 = Array<Vector4, 4>;

    Matrix4 x[32];
    for (size_t i = 0; i<32; ++i)
        for (size_t k = 0; k<Matrix4::Size; ++k)
            for (size_t j = 0; j<Matrix4::Size; ++j)
                x[i /* slice */ ][k /* col */ ][j /* row */ ] = Value(i*100+k*10+j);

    for (size_t i = 0; i<32; ++i) {
        /* Direct load */
        Matrix4 q = gather<Matrix4>(x, i);
        Matrix4 q2 = Array44(Matrix4(
            0, 10, 20, 30,
            1, 11, 21, 31,
            2, 12, 22, 32,
            3, 13, 23, 33
        )) + Value(100*i);
        assert(q == q2);
    }

    /* Variant 2 */
    for (size_t i = 0; i < 32; ++i)
        assert(gather<Matrix4>(x, i) == x[i]);

    /* Nested gather + slice */
    auto q = gather<Matrix4P>(x, arange<UInt32P>());
    for (size_t i = 0; i < T::Size; ++i)
        assert(slice(q, i) == x[i]);

    /* Masked nested gather */
    q = gather<Matrix4P>(x, arange<UInt32P>(), arange<UInt32P>() < 2u);
    for (size_t i = 0; i < T::Size; ++i) {
        if (i < 2u)
            assert(slice(q, i) == x[i]);
        else
            assert(slice(q, i) == zero<Matrix4>());
    }

    /* Nested scatter */
    q = gather<Matrix4P>(x, arange<UInt32P>());
    scatter(x, q + full<Matrix4>(Value(1000)), arange<UInt32P>());
    scatter(x, q + full<Matrix4>(Value(2000)), arange<UInt32P>(), arange<UInt32P>() < 2u);
    for (size_t i = 0; i < T::Size; ++i) {
        if (i < 2u)
            assert(slice(q, i) + full<Matrix4>(Value(2000)) == x[i]);
        else
            assert(slice(q, i) + full<Matrix4>(Value(1000)) == x[i]);
    }

    /* Nested prefetch -- doesn't do anything here, let's just make sure it compiles */
    prefetch<Matrix4P>(x, arange<UInt32P>());
    prefetch<Matrix4P>(x, arange<UInt32P>(), arange<UInt32P>() < 2u);

    /* Nested gather + slice for dynamic arrays */
    using UInt32X   = DynamicArray<UInt32P>;
    using TX        = DynamicArray<T>;
    using Matrix4X  = Matrix<TX, 4>;
    auto q2 = gather<Matrix4X>(x, arange<UInt32X>(2));
    q2 = q2 + full<Matrix4>(Value(1000));
    scatter(x, q2, arange<UInt32X>(2));

    for (size_t i = 0; i < T::Size; ++i) {
        if (i < 2u)
            assert(slice(q, i) + full<Matrix4>(Value(3000)) == x[i]);
        else
            assert(slice(q, i) + full<Matrix4>(Value(1000)) == x[i]);
    }
}

#if defined(_MSC_VER)
  template <typename T, std::enable_if_t<T::Size != 8, int> = 0>
#else
  template <typename T, std::enable_if_t<T::Size == 31, int> = 0>
#endif
void test04_nested_gather_packed_impl() { }

ENOKI_TEST_ALL(test04_nested_gather_packed) {
    test04_nested_gather_packed_impl<T>();
}

#if defined(_MSC_VER)
  template <typename T, std::enable_if_t<T::Size == 8, int> = 0>
#else
  template <typename T, std::enable_if_t<T::Size != 31, int> = 0>
#endif
void test05_nested_gather_nonpacked_impl() {
    using Value = value_t<T>;
    using UInt32P = Packet<uint32_t, T::Size>;
    using Vector3 = Array<Value, 3>;
    using Matrix3 = Matrix<Value, 3>;
    using Matrix3P = Matrix<T, 3>;
    using Array33 = Array<Vector3, 3>;

    Matrix3 x[32];
    for (size_t i = 0; i<32; ++i)
        for (size_t k = 0; k<Matrix3::Size; ++k)
            for (size_t j = 0; j<Matrix3::Size; ++j)
                x[i /* slice */ ][k /* col */ ][j /* row */ ] = Value(i*100+k*10+j);

    for (size_t i = 0; i<32; ++i) {
        /* Direct load */
        Matrix3 q = gather<Matrix3, 0, false>(x, i);
        Matrix3 q2 = Array33(Matrix3(
            0, 10, 20,
            1, 11, 21,
            2, 12, 22
        )) + Value(100*i);
        assert(q == q2);
    }

    /* Variant 2 */
    for (size_t i = 0; i < 32; ++i)
        assert((gather<Matrix3, 0, false>(x, i)) == x[i]);

    /* Nested gather + slice */
    auto q = gather<Matrix3P, 0, false>(x, arange<UInt32P>());
    for (size_t i = 0; i < T::Size; ++i)
        assert(slice(q, i) == x[i]);

    /* Masked nested gather */
    q = gather<Matrix3P, 0, false>(x, arange<UInt32P>(), arange<UInt32P>() < 2u);
    for (size_t i = 0; i < T::Size; ++i) {
        if (i < 2u)
            assert(slice(q, i) == x[i]);
        else
            assert(slice(q, i) == zero<Matrix3>());
    }

    /* Nested scatter */
    q = gather<Matrix3P, 0, false>(x, arange<UInt32P>());
    scatter<0, false>(x, q + full<Matrix3>(Value(1000)), arange<UInt32P>());
    scatter<0, false>(x, q + full<Matrix3>(Value(2000)), arange<UInt32P>(), arange<UInt32P>() < 2u);
    for (size_t i = 0; i < T::Size; ++i) {
        if (i < 2u)
            assert(slice(q, i) + full<Matrix3>(Value(2000)) == x[i]);
        else
            assert(slice(q, i) + full<Matrix3>(Value(1000)) == x[i]);
    }

    /* Nested prefetch -- doesn't do anything here, let's just make sure it compiles */
    prefetch<Matrix3P, false, 2, 0, false>(x, arange<UInt32P>());
    prefetch<Matrix3P, false, 2, 0, false>(x, arange<UInt32P>(), arange<UInt32P>() < 2u);

    /* Nested gather + slice for dynamic arrays */
    using UInt32X   = DynamicArray<UInt32P>;
    using TX        = DynamicArray<T>;
    using Matrix3X  = Matrix<TX, 3>;
    auto q2 = gather<Matrix3X, 0, false>(x, arange<UInt32X>(2));
    q2 = q2 + full<Matrix3>(Value(1000));
    scatter<0, false>(x, q2, arange<UInt32X>(2));

    for (size_t i = 0; i < T::Size; ++i) {
        if (i < 2u)
            assert(slice(q, i) + full<Matrix3>(Value(3000)) == x[i]);
        else
            assert(slice(q, i) + full<Matrix3>(Value(1000)) == x[i]);
    }
}

#if defined(_MSC_VER)
  template <typename T, std::enable_if_t<T::Size != 8, int> = 0>
#else
  template <typename T, std::enable_if_t<T::Size == 31, int> = 0>
#endif
void test05_nested_gather_nonpacked_impl() { }

ENOKI_TEST_ALL(test05_nested_gather_nonpacked) {
    test05_nested_gather_nonpacked_impl<T>();
}

ENOKI_TEST_ALL(test06_range) {
    alignas(alignof(T)) Value mem[Size*10];
    for (size_t i = 0; i < Size*10; ++i)
        mem[i] = 1;
    using Index = uint_array_t<T>;
    MEMORY_BARRIER();
    T sum = zero<T>();
    for (auto pair : range<Index>((10*Size)/3))
        sum += gather<T>(mem, pair.first, pair.second);
    assert(((10*Size)/3) == hsum(sum));
}

ENOKI_TEST_ALL(test07_range_2d) {
    alignas(alignof(T)) Value mem[4*5*6];
    for (size_t i = 0; i < 4*5*6; ++i)
        mem[i] = 0;
    using Index3 = Array<uint_array_t<T>, 3>;
    for (auto pair : range<Index3>(4u, 5u, 6u)) {
        auto index = pair.first[0] +
                     pair.first[1] * 4u + pair.first[2] * 20u;
        scatter(mem, T(1), index, pair.second);
    }
    MEMORY_BARRIER();
    for (size_t i = 0; i < 4*5*6; ++i) {
        assert(mem[i] == 1);
    }
}

ENOKI_TEST_ALL(test08_nested_gather_strides) {
    using Vector4x = Array<Value, 4>;
    using Vector4xP = Array<T, 4>;
    using Int = int_array_t<T>;

    Value x[100];
    memset(x, 0, sizeof(x));
    for (int i = 0; i < 8; ++i) {
        x[i] = (Value) i;
        x[i + 60] = (Value) (i + 10);
    }
    MEMORY_BARRIER();

    assert((gather<Vector4x>(x, 0) == Vector4x(0, 1, 2, 3)));
    assert((gather<Vector4x>(x, 1) == Vector4x(4, 5, 6, 7)));
    assert((gather<Vector4x, 60*sizeof(Value)>(x, 1) == Vector4x(10, 11, 12, 13)));
    assert((gather<Vector4xP>(x, Int(0)) == Vector4xP(0, 1, 2, 3)));
    assert((gather<Vector4xP>(x, Int(1)) == Vector4xP(4, 5, 6, 7)));
    assert((gather<Vector4xP, 60*sizeof(Value)>(x, Int(1)) == Vector4xP(10, 11, 12, 13)));
}

ENOKI_TEST_INT(test09_gather_mask) {
    if ((T::Size & (T::Size - 1)) != 0)
        return;
    using Scalar = scalar_t<T>;
    using TX = DynamicArray<T>;
    using MaskP = mask_t<T>;
    using MaskX = mask_t<TX>;

    MaskX p = eq(arange<TX>(50) & (Scalar) 1, (Scalar) 0);
    MaskP result = gather<MaskP>(p, arange<T>() + (Scalar) 1);
    MaskP target = eq((arange<T>() + (Scalar) 1) & (Scalar) 1, (Scalar) 0);
    assert(target == result);
}
