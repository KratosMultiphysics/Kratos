/*
    tests/basic.cpp -- tests linear algebra-related routines

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/matrix.h>
#include <enoki/transform.h>
#include <enoki/random.h>
#include <enoki/dynamic.h>

template <typename T, std::enable_if_t<std::is_signed<typename T::Value>::value, int> = 0>
void test01_dot_signed() {
    using Value = typename T::Value;
    constexpr size_t Size = T::Size;
    Value expected1 = Value((Size * (2 * Size - 1) * (Size - 1)) / 6);
    Value expected2 = Value((Size * (Size - 1)) / 2);
    T value = arange<T>();
    T value2(Value(1));
    assert(dot(value, -value) == -expected1);
    assert(dot(value, -value2) == -expected2);
    assert(abs_dot(value, -value) == expected1);
    assert(abs_dot(value, -value2) == expected2);
}

template <typename T, std::enable_if_t<!std::is_signed<typename T::Value>::value, int> = 0>
void test01_dot_signed() { }

ENOKI_TEST_ALL(test01_dot) {
    Value expected1 = Value((Size * (2 * Size - 1) * (Size - 1)) / 6);
    Value expected2 = Value((Size * (Size - 1)) / 2);
    T value = arange<T>();
    T value2(Value(1));
    assert(dot(value, value)  == expected1);
    assert(dot(value, value2) == expected2);
    test01_dot_signed<T>();
}

template <typename T> void test02_vecops_float() {
    using Value = value_t<T>;

    /* Extra tests for horizontal operations */
    T v(1.f, 2.f, 3.f);
    assert(v.x() == 1);
    assert(v.y() == 2);
    assert(v.z() == 3);
    assert(v.x() == Value(1));
    assert(v.y() == Value(2));
    assert(v.z() == Value(3));
    assert(std::abs(norm(v) - 3.74165738677394f) < 1e-5f);
    assert(std::abs(squared_norm(v) - std::pow(3.74165738677394f, 2.f)) < 1e-5f);
    assert(hsum(abs(normalize(v) - T(0.26726f, 0.53452f, 0.80178f))) < 1e-5f);
    assert(hsum(v) == 6.f);
    assert(hprod(v) == 6.f);
    assert(hmax(v) == 3.f);
    if (T::ActualSize == 4)
        v.coeff(3) = -4;
    assert(hmin(v) == 1.f);
    assert(dot(v, v) == 14.f);

    v = T(0);
    if (T::ActualSize == 4)
        v.coeff(3) = -1;
    assert(!any(v < Value(0)));

    v = T(-1);
    if (T::ActualSize == 4)
        v.coeff(3) = 0;
    assert(all(v < Value(0)));

    /* Test cross product */
    assert(cross(T(1, 2, 3), T(1, 0, 0)) == T(0, 3, -2));
    assert(cross(T(1, 2, 3), T(0, 1, 0)) == T(-3, 0, 1));
    assert(cross(T(1, 2, 3), T(0, 0, 1)) == T(2, -1, 0));
    assert(cross(T(1, 1, 1), T(1, 2, 3)) == T(1, -2, 1));

    v.x() = 3; v.y() = 4; v.z() = 5;
    if (T::ActualSize == 4)
        v.coeff(3) = -1;
    assert(v == T(3, 4, 5));
    assert(v != T(3, 2, 5));
}

ENOKI_TEST(test02_vecops_float) {
    test02_vecops_float<Array<float, 3>>();
}

ENOKI_TEST(test02_vecops_double) {
    test02_vecops_float<Array<double, 3>>();
}

ENOKI_TEST(array_float_04_transpose) {
    using T  = Array<float, 4>;
    using T2 = Array<T, 4>;

    assert(
        transpose(T2(T(1, 2, 3, 4), T(5, 6, 7, 8), T(9, 10, 11, 12), T(13, 14, 15, 16))) ==
        T2(T(1, 5, 9, 13), T(2, 6, 10, 14), T(3, 7, 11, 15), T(4, 8, 12, 16)));
}

ENOKI_TEST(array_double_04_transpose) {
    using T  = Array<double, 4>;
    using T2 = Array<T, 4>;

    assert(
        transpose(T2(T(1, 2, 3, 4), T(5, 6, 7, 8), T(9, 10, 11, 12), T(13, 14, 15, 16))) ==
        T2(T(1, 5, 9, 13), T(2, 6, 10, 14), T(3, 7, 11, 15), T(4, 8, 12, 16)));
}

ENOKI_TEST(array_float_05_outer_product) {
    using Vector4f = Array<float, 4>;
    using Vector3f = Array<float, 3>;

    assert(to_string(full<Vector3f>(Vector4f(1, 2, 3, 4)) * Vector3f(0, 1, 0)) == "[[0, 1, 0],\n [0, 2, 0],\n [0, 3, 0],\n [0, 4, 0]]");
    assert(to_string(full<float>(Vector4f(1, 2, 3, 4)) * 3.f) == "[3, 6, 9, 12]");
    assert(to_string(full<Vector3f>(3.f) * Vector3f(1, 2, 3)) == "[3, 6, 9]");
}

ENOKI_TEST(array_float_06_head_tail) {
    using T  = Array<float, 4>;

    auto t = T(1, 2, 3, 4);
    assert(to_string(head<2>(t)) == "[1, 2]");
    assert(to_string(tail<2>(t)) == "[3, 4]");
}

ENOKI_TEST(array_float_02_test07_matrix) {
    using M2f = Matrix<float, 2>;
    using V2f = Array<float, 2>;

    auto a = M2f(1, 2, 3, 4);
    auto b = V2f(1, 1);
    assert(to_string(a) == "[[1, 2],\n [3, 4]]");
    assert(to_string(transpose(a)) == "[[1, 3],\n [2, 4]]");
    assert(a(0, 1) == 2);
    assert(to_string(a.coeff(0)) == "[1, 3]");
    assert(to_string(a.coeff(1)) == "[2, 4]");
    assert(to_string(a*a) == "[[7, 10],\n [15, 22]]");
    assert(to_string(a*a) == "[[7, 10],\n [15, 22]]");
    assert(to_string(a*b) == "[3, 7]");

    assert(to_string(M2f::from_rows(V2f(1, 2), V2f(3, 4))) == "[[1, 2],\n [3, 4]]");
    assert(to_string(M2f::from_cols(V2f(1, 2), V2f(3, 4))) == "[[1, 3],\n [2, 4]]");

    /* Ensure that matrix-vector multiplication broadcasts compile */ {
        using Float = float;

        using FloatP = Packet<Float, 8>;

        using Vector4f = Array<Float, 4>;
        using Vector4fP = Array<FloatP, 4>;

        using Matrix4f = Matrix<Float, 4>;

        /* Variant 1 */ {
            Matrix4f M;
            Vector4f x;
            Vector4f y = M * x;
            (void) y;
        }

        /* Variant 2 */ {
            Matrix4f M;
            Vector4fP x;
            Vector4fP y = M * x;
            (void) y;
        }
    }
}

template <typename T> void test_concat() {
    assert((concat(Array<T, 3>((T) 1, (T) 2, (T) 3), (T) 4) ==
            Array<T, 4>((T) 1, (T) 2, (T) 3, (T) 4)));
}

ENOKI_TEST(array_float_04_concat) { test_concat<float>(); }
ENOKI_TEST(array_uint32_04_concat) { test_concat<uint32_t>(); }
ENOKI_TEST(array_double_04_concat) { test_concat<double>(); }
ENOKI_TEST(array_uint64_t_04_concat) { test_concat<uint64_t>(); }

template <typename Type, size_t Size_>
struct Vector : enoki::StaticArrayImpl<Type, Size_, false, Vector<Type, Size_>> {

    using Base = enoki::StaticArrayImpl<Type, Size_, false, Vector<Type, Size_>>;
    ENOKI_ARRAY_IMPORT(Base, Vector)

    using ArrayType = Vector;
    using MaskType = Mask<Type, Size_>;

    /// Helper alias used to transition between vector types
    template <typename T> using ReplaceValue = Vector<T, Size_>;
};


template <typename T> void test_bcast() {
    using Packet = enoki::Packet<T, 4>;
    using Vector4 = Vector<T, 4>;
    using Vector4P = Vector<Packet, 4>;

    assert(hsum(Vector4P(Vector4(T(1), T(2), T(3), T(4)))) == Vector4(T(10), T(10), T(10), T(10)));
    assert(hsum(Vector4P(Packet(T(1), T(2), T(3), T(4)))) == Vector4(T(4), T(8), T(12), T(16)));

    using Array4s = Array<size_t, 4>;
    assert(count(mask_t<Vector4P>(Vector4(T(1), T(2), T(3), T(4)) < value_t<T>(3))) == Array4s(2u, 2u, 2u, 2u));
    assert(count(mask_t<Vector4P>(Packet(T(1), T(2), T(3), T(4)) < value_t<T>(3))) == Array4s(4u, 4u, 0u, 0u));
}

ENOKI_TEST(array_float_04_bcast) { test_bcast<float>(); }
ENOKI_TEST(array_uint32_04_bcast) { test_bcast<uint32_t>(); }
ENOKI_TEST(array_double_04_bcast) { test_bcast<double>(); }
ENOKI_TEST(array_uint64_t_04_bcast) { test_bcast<uint64_t>(); }


ENOKI_TEST(transform_decompose) {
    using Matrix4f = Matrix<float, 4>;
    using Matrix3f = Matrix<float, 3>;
    using Vector3f = Array<float, 3>;
    using Quaternion4f = Quaternion<float>;

    Matrix4f A(
         0.99646652, -0.13867289,  0.14220636,  1.,
         0.07042366,  1.99456394, -0.06498755,  2.,
        -0.04577128,  0.04984837,  2.9959228 ,  3.,
         0.        ,  0.        ,  0.        ,  1.
    );
    Matrix3f s;
    Quaternion4f q;
    Vector3f t;
    std::tie(s, q, t) = transform_decompose(A);
    auto result2 = transform_compose(s, q, t);
    assert(frob(A - result2) < 1e-6f);
}

ENOKI_TEST(transform_compose_inverse) {
    using Matrix = enoki::Matrix<float, 4>;
    auto rng = PCG32<float>();

    for (int k = 0; k<10; ++k) {
        Matrix m = zero<Matrix>();
        for (size_t i=0; i<3; ++i)
            for (size_t j=0; j<4; ++j)
                m(i, j) = rng.next_float32();
        m(3, 3) = 1.f;

        auto x = enoki::transform_decompose(m);
        auto m2 = enoki::transform_compose(std::get<0>(x), std::get<1>(x), std::get<2>(x));
        auto m3 = enoki::transform_compose_inverse(std::get<0>(x), std::get<1>(x), std::get<2>(x));

        auto diff1 = frob(m-m2)/frob(m);
        auto diff2 = frob(inverse(m)-m3)/frob(inverse(m));

        assert(diff1 < 1e-6f);
        assert(diff2 < 1e-6f);
    }
}

ENOKI_TEST(test_unit_angle) {
    using Vector3f = Array<float, 3>;
    using Random = PCG32<Vector3f>;

    Random rng;

    for (int i = 0; i < 30; ++i) {
        Vector3f a = normalize(rng.next_float32() * 2.f - 1.f);
        Vector3f b = normalize(rng.next_float32() * 2.f - 1.f);

        assert(std::abs(acos(dot(a, b)) - unit_angle(a, b)) < 1e-6f);
    }
}

ENOKI_TEST(full) {
    using Vector4f  = Array<float, 4>;
    using MyMatrix = Array<Vector4f, 4>;
    MyMatrix result = full<MyMatrix>(10.f);
    assert(to_string(result) == "[[10, 10, 10, 10],\n [10, 10, 10, 10],\n [10, 10, 10, 10],\n [10, 10, 10, 10]]");
    result = full<Vector4f>(Vector4f(1, 2, 3, 4));
    assert(to_string(result) == "[[1, 1, 1, 1],\n [2, 2, 2, 2],\n [3, 3, 3, 3],\n [4, 4, 4, 4]]");
    result = MyMatrix(Vector4f(1, 2, 3, 4));
    assert(to_string(result) == "[[1, 2, 3, 4],\n [1, 2, 3, 4],\n [1, 2, 3, 4],\n [1, 2, 3, 4]]");

    Matrix<float, 4> result2 = full<Matrix<float, 4>>(10.f);
    assert(to_string(result2) == "[[10, 10, 10, 10],\n [10, 10, 10, 10],\n [10, 10, 10, 10],\n [10, 10, 10, 10]]");
}

ENOKI_TEST(masked_assignment) {
    using FloatP = Packet<float>;
    using Matrix4f = Matrix<float, 4>;
    using Matrix4fP = Matrix<FloatP, 4>;

    Matrix4fP z = identity<Matrix4f>();
    masked(z, true) *= 2.f;
    masked(z, z > 0.f) *= 2.f;
    masked(z, eq(arange<FloatP>(), 0.f)) *= 2.f;
    assert(z.coeff(0, 0, 0) == 8.f);

    if (FloatP::Size > 1) {
        assert(z.coeff(0, 0, 1) == 4.f);
        assert(z.coeff(1, 0, 1) == 0.f);
    }
}
