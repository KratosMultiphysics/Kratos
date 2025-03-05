/*
    tests/complex.cpp -- tests complex numbers and quaternions

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/complex.h>
#include <enoki/quaternion.h>
#include <enoki/transform.h>
#include <enoki/dynamic.h>

using Cf = Complex<double>;
using Qf = Quaternion<double>;
using V3 = Array<double, 3>;

ENOKI_TEST(test00_complex_str) {
    assert(to_string(Cf(1.0)) == "1 + 0i");
    assert(to_string(Cf(1.0, 2.0)) == "1 + 2i");
    assert(to_string(conj(Cf(1.0, 2.0))) == "1 - 2i");
}

ENOKI_TEST(test01_quat_str) {
    assert(to_string(Qf(2.0, 3.0, 4.0, 1.0)) == "1 + 2i + 3j + 4k");
    assert(to_string(conj(Qf(2.0, 3.0, 4.0, 1.0))) == "1 - 2i - 3j - 4k");
}

ENOKI_TEST(test02_complex_mult) {
    assert(to_string(Cf(1.0, 2.0) * Cf(-4.f, 3.0)) == "-10 - 5i");
}

ENOKI_TEST(test03_quat_mult) {
    assert(to_string(Qf(1.0, 2.0, -1.0, 3.0) * Qf(-4.f, 3.0, 5.f, -2.f)) == "-3 - 1i + 4j + 28k");
}

ENOKI_TEST(test04_complex_rcp) {
    auto a = Cf(2.f, 1.f);
    assert(to_string(rcp(a)) == "0.4 - 0.2i");
    assert(abs(a/a - Cf(1.f)) < 1e-5);
}

ENOKI_TEST(test05_quat_rcp) {
    auto a = normalize(Qf(1.f, 2.f, 3.f, 4.f));
    auto b = rcp(a);
    auto c = a*b - Qf(1.f);
    assert(abs(c) < 1e-5);
    assert(abs(a/a - Qf(1.f)) < 1e-5);
}

ENOKI_TEST(test06_complex_decomp) {
    assert(real(Cf(1, 2)) == 1);
    assert(imag(Cf(1, 2)) == 2);
}

ENOKI_TEST(test07_quat_decomp) {
    assert(real(Qf(1, 2, 3, 4)) == 4);
    assert(imag(Qf(1, 2, 3, 4)) == V3(1, 2, 3));
}

ENOKI_TEST(test08_complex_exp) {
    assert(abs(exp(Cf(1, 2)) - Cf(-1.1312, 2.47173)) < 1e-5);
}

ENOKI_TEST(test09_quat_exp) {
    assert(abs(exp(Qf(2, 3, 4, 1)) - Qf(-0.78956, -1.18434, -1.57912, 1.69392)) < 1e-5);
}

ENOKI_TEST(test10_complex_log) {
    assert(abs(log(Cf(1, 2)) - Cf(0.804719, 1.10715)) < 1e-5);
}

ENOKI_TEST(test11_quat_log) {
    assert(abs(log(Qf(2, 3, 4, 1)) - Qf(0.51519, 0.772785, 1.03038, 1.7006)) < 1e-5);
}

ENOKI_TEST(test12_complex_sqrt) {
    assert(abs(sqrt(Cf(1, 2)) - Cf(1.27202, 0.786151)) < 1e-6);
}

ENOKI_TEST(test13_quat_sqrt) {
    assert(abs(sqrt(Qf(2, 3, 4, 1)) - Qf(0.555675, 0.833512, 1.11135, 1.79961)) < 1e-5);
}

ENOKI_TEST(test14_complex_sin_cos_tan) {
    assert(abs(sin(Cf(1, 2)) - Cf(3.16578, 1.9596)) < 1e-5);
    assert(abs(cos(Cf(1, 2)) - Cf(2.03272, -3.0519)) < 1e-5);
    assert(abs(tan(Cf(1, 2)) - Cf(0.0338128, 1.01479)) < 1e-5);
    auto sc = sincos(Cf(1, 2));
    assert(abs(sc.first - Cf(3.16578, 1.9596)) < 1e-5);
    assert(abs(sc.second - Cf(2.03272, - 3.0519)) < 1e-5);
}

ENOKI_TEST(test15_complex_sinh_cosh_tanh) {
    assert(abs(sinh(Cf(1, 2)) - Cf(-0.489056, 1.40312)) < 1e-5);
    assert(abs(cosh(Cf(1, 2)) - Cf(-0.642148, 1.06861)) < 1e-5);
    assert(abs(tanh(Cf(1, 2)) - Cf(1.16674, -0.243458)) < 1e-5);
    auto sc = sincosh(Cf(1, 2));
    assert(abs(sc.first - Cf(-0.489056, 1.40312)) < 1e-5);
    assert(abs(sc.second - Cf(-0.642148, 1.06861)) < 1e-5);
}

ENOKI_TEST(test16_complex_asin_acos_atan) {
    assert(abs(asin(Cf(1, 2)) - Cf(0.427079, 1.52857)) < 1e-5);
    assert(abs(acos(Cf(1, 2)) - Cf(1.14372, -1.52857)) < 1e-5);
    assert(abs(atan(Cf(1, 2)) - Cf(1.33897, 0.402359)) < 1e-5);
}

ENOKI_TEST(test17_complex_asinh_acosh_atanh) {
    assert(abs(asinh(Cf(1, 2)) - Cf(1.46935, 1.06344)) < 1e-5);
    assert(abs(acosh(Cf(1, 2)) - Cf(1.52857, 1.14372)) < 1e-5);
    assert(abs(atanh(Cf(1, 2)) - Cf(0.173287, 1.1781)) < 1e-5);
}

using FloatP = Packet<float>;
using FloatX = DynamicArray<FloatP>;
using Quaternion4f = Quaternion<float>;
using Quaternion4X = Quaternion<FloatX>;
using Matrix4X  = Matrix<FloatX, 4>;
using Matrix4f  = Matrix<float, 4>;
using Matrix4fP = Matrix<FloatP, 4>;
using Vector3f  = Array<float, 3>;
using Vector4f  = Array<float, 4>;

Matrix4X slerp_matrix(const Quaternion4X &x, const Quaternion4X &y, float t) {
    return vectorize([t](auto &&x, auto &&y) { return quat_to_matrix<Matrix4fP>(slerp(x, y, t)); }, x, y);
};

Quaternion4X to_quat(const Matrix4X &m) {
    return vectorize([](auto &&m) { return matrix_to_quat(m); }, m);
};

ENOKI_TEST(test18_complex_vectorize_scalar) {
    Quaternion4f a = normalize(Quaternion4f(1, 2, 3, 4));
    Quaternion4f b = normalize(Quaternion4f(0, 0, 0, 1));

    Quaternion4X x, y;
    set_slices(x, 1);
    set_slices(y, 1);
    slice(x, 0) = a;
    slice(y, 0) = b;
    auto tmp0 = slerp_matrix(x, y, 0.5f);
    auto tmp1 = to_quat(tmp0);
    Quaternion4f result = slice(tmp1, 0);
    Quaternion4f ref = normalize(a+b);
    assert(abs(result - ref) < 1e-5f);
}

ENOKI_TEST(test19_rotation) {
    auto axis = normalize(Vector3f(1.f, 2.f, 3.f));
    Vector4f input(0.8f, 0.3f, 0.2f, 0.0f);
    float angle = 0.5f;

    auto quat1 = rotate<Quaternion4f>(axis, angle);
    auto r1 = Vector4f(quat1 * Quaternion4f(input) * conj(quat1));

    auto mtx2 = rotate<Matrix4f>(axis, angle);
    auto r2 = mtx2 * input;

    auto mtx1 = quat_to_matrix<Matrix4f>(quat1);
    auto r3 = mtx1 * input;

    auto quat2 = matrix_to_quat(mtx2);
    auto r4 = Vector4f(quat2 * Quaternion4f(input) * conj(quat2));

    assert(norm(r1-r2) < 1e-6f);
    assert(norm(r1-r3) < 1e-6f);
    assert(norm(r1-r4) < 1e-6f);
}

ENOKI_TEST(test20_sincos_arg) {
    auto result = sincos_arg_diff(Cf(-1.01264771f, 1.1261553f), Cf(-0.70017226f, 1.24072149f));
    assert(abs(result.first - 0.2168644f) < 1e-6f);
    assert(abs(result.second - 0.97620174f) < 1e-6f);

    result = sincos_arg_diff(Cf(-0.08012004f, 0.86251237f), Cf(-1.22284338f, 0.86829703f));
    assert(abs(result.first + 0.75831358f) < 1e-6f);
    assert(abs(result.second - 0.65188996f) < 1e-6f);
}
