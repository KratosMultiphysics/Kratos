/*
    tests/basic.cpp -- tests basic floating point operations

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"

ENOKI_TEST_FLOAT(test01_div_fp) {
    auto sample = test::sample_values<Value>();

    test::validate_binary<T>(sample,
        [](const T &a, const T &b) -> T { return a / b; },
        [](Value a, Value b) -> Value { return a / b; },
#if !defined(ENOKI_ARM_32)
        0.f
#else
        1e-6f
#endif
    );

    test::validate_binary<T>(sample,
        [](const T &a, const T &b) -> T { T x(a); x /= b; return x; },
        [](Value a, Value b) -> Value { return a / b; },
#if !defined(ENOKI_ARM_32)
        0.f
#else
        1e-6f
#endif
    );

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return a / 3.f; },
        [](Value a) -> Value { return a / 3.f; }, 1e-6f
    );

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return a / Value(3); },
        [](Value a) -> Value { return a / Value(3); }, 1e-6f
    );

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return Value(3) / a; },
        [](Value a) -> Value { return Value(3) / a; }, 1e-6f
    );

#if !defined(ENOKI_X86_AVX512F)
    /* In AVX512 mode, the approximate reciprocal function is
       considerably more accurate and this test fails */
    if (std::is_same<Value, float>::value && has_sse42) {
        // Make sure that division optimization is used
        assert(T (123.f) / 123.f != T (1.f));
    }
#endif
}

ENOKI_TEST_FLOAT(test02_ceil) {
    auto sample = test::sample_values<Value>();

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return ceil(a); },
        [](Value a) -> Value { return std::ceil(a); }
    );

    Array<T, 4> x(3.4f); Array<T&, 4> y(x);
    assert(ceil(x) == ceil(y));
}

ENOKI_TEST_FLOAT(test03_floor) {
    auto sample = test::sample_values<Value>();

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return floor(a); },
        [](Value a) -> Value { return std::floor(a); }
    );

    Array<T, 4> x(3.4f); Array<T&, 4> y(x);
    assert(floor(x) == floor(y));
}

ENOKI_TEST_FLOAT(test04_round) {
    auto sample = test::sample_values<Value>();

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return round(a); },
        [](Value a) -> Value { return std::rint(a); }
    );

    Array<T, 4> x(3.4f); Array<T&, 4> y(x);
    assert(round(x) == round(y));
}

ENOKI_TEST_FLOAT(test05_trunc) {
    auto sample = test::sample_values<Value>();

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return trunc(a); },
        [](Value a) -> Value { return std::trunc(a); }
    );

    Array<T, 4> x(3.4f); Array<T&, 4> y(x);
    assert(trunc(x) == trunc(y));
}

ENOKI_TEST_FLOAT(test06_sqrt) {
    auto sample = test::sample_values<Value>();

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return sqrt(a); },
        [](Value a) -> Value { return std::sqrt(a); },
#if !defined(ENOKI_ARM_32)
        0.f
#else
        1e-6f
#endif
    );

    Array<T, 4> x(3.4f); Array<T&, 4> y(x);
    assert(sqrt(x) == sqrt(y));
}

ENOKI_TEST_FLOAT(test07_rsqrt) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return rsqrt(a); },
        [](double a) { return 1 / std::sqrt(a); },
        Value(1e-6), Value(1024), 3
    );

    test::probe_accuracy<T>(
        [](const T &a) -> T {
            T result;
            for (size_t i = 0; i < Size; ++i)
               result.coeff(i) = rsqrt(a.coeff(i));
            return result;
        },
        [](double a) { return 1 / std::sqrt(a); },
        Value(1e-6), Value(1024), 3
    );
}

ENOKI_TEST_FLOAT(test08_rcp) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return rcp(a); },
        [](double a) { return 1 / a; },
        Value(1e-6), Value(1024), 2
    );

    test::probe_accuracy<T>(
        [](const T &a) -> T {
            T result;
            for (size_t i = 0; i < Size; ++i)
               result.coeff(i) = rcp(a.coeff(i));
            return result;
        },
        [](double a) { return 1 / a; },
        Value(1e-6), Value(1024), 2
    );
}

ENOKI_TEST_FLOAT(test09_sign) {
    auto sample = test::sample_values<Value>();

    test::validate_unary<T>(sample,
        [](const T &a) -> T { return enoki::sign(a); },
        [](Value a) -> Value { return std::copysign(1.f, a); }
    );

    Array<T, 4> x(3.4f); Array<T&, 4> y(x);
    assert(sign(x) == sign(y));
}

ENOKI_TEST_FLOAT(test10_isinf) {
    auto sample = test::sample_values<Value>();

    using enoki::isinf;
    test::validate_unary<T>(sample,
        [](const T &a) -> T { return select(enoki::isinf(a), T(1), T(0)); },
        [](Value a) -> Value { return Value(std::isinf(a) ? 1 : 0); }
    );
}

ENOKI_TEST_FLOAT(test11_isnan) {
    auto sample = test::sample_values<Value>();

    using enoki::isnan;
    test::validate_unary<T>(sample,
        [](const T &a) -> T { return select(enoki::isnan(a), T(1), T(0)); },
        [](Value a) -> Value { return Value(std::isnan(a) ? 1 : 0); }
    );
}

ENOKI_TEST_FLOAT(test12_isfinite) {
    auto sample = test::sample_values<Value>();

    using enoki::isfinite;
    test::validate_unary<T>(sample,
        [](const T &a) -> T { return select(enoki::isfinite(a), T(1), T(0)); },
        [](Value a) -> Value { return Value(std::isfinite(a) ? 1 : 0); }
    );
}

ENOKI_TEST_FLOAT(test13_nan_initialization) {
    T x;
    for (size_t i = 0; i < Size; ++i)
        assert(std::isnan(x[i]));
}

ENOKI_TEST(test14_half) {
    using T = Array<float, 4>;
    using THalf = Array<half, 4>;

    for (uint32_t i = 0; i < 0xFFFF; ++i) {
        uint16_t data[8]  = { (uint16_t) i };
        uint16_t data3[8] = { (uint16_t) i };

        float f1 = T(load<THalf>((const half *) data))[0];
        float f2 = (float) half::from_binary(data[0]);

        bool both_nan = std::isnan(f1) && std::isnan(f2);

        assert((memcpy_cast<uint32_t>(f1) ==
                memcpy_cast<uint32_t>(f2)) || both_nan);

        half data2[8];
        store(data2, THalf(T(f1)));
        data3[0] = half(f2).value;

        assert((data2[0].value == data3[0]) || both_nan);
    }
}

ENOKI_TEST_FLOAT(test16_hypot) {
    auto sample = test::sample_values<Value>();

    test::validate_binary<T>(sample,
                             [](const T &a, const T &b) -> T {
                                 return enoki::hypot(a, b);
                             },
                             [](Value a, Value b) -> Value {
                                 if (std::isnan(a) || std::isnan(b))
                                     return std::numeric_limits<Value>::quiet_NaN();
                                 else if (!std::isfinite(a) || !std::isfinite(b))
                                     return std::numeric_limits<Value>::infinity();
                                 else
                                    return std::hypot(a, b);
                             },
                             1e-6f);
}

ENOKI_TEST_FLOAT(test17_next_float) {
    Value inf = std::numeric_limits<Value>::infinity();
    Value nan = std::numeric_limits<Value>::quiet_NaN();
    Value zero = Value(0), one = Value(1.f);

    assert(next_float(T( zero))  == T(std::nextafter(zero, inf)));
    assert(next_float(T(-zero)) == T(std::nextafter(-zero, inf)));
    assert(next_float(T( one)) == T(std::nextafter( one, inf)));
    assert(next_float(T(-one)) == T(std::nextafter(-one, inf)));
    assert(next_float(T( inf))  == inf);
    assert(next_float(T(-inf)) == -inf);
    assert(all(enoki::isnan(next_float(T(nan)))));

    assert(prev_float(T( zero))  == T(std::nextafter(zero, -inf)));
    assert(prev_float(T(-zero)) == T(std::nextafter(-zero, -inf)));
    assert(prev_float(T( one)) == T(std::nextafter( one, -inf)));
    assert(prev_float(T(-one)) == T(std::nextafter(-one, -inf)));
    assert(prev_float(T( inf))  == inf);
    assert(prev_float(T(-inf)) == -inf);
    assert(all(enoki::isnan(prev_float(T(nan)))));
}

ENOKI_TEST_FLOAT(test18_fmod) {
    T a = Value(5.1);
    T b = Value(3.0);
    T c = Value(2.1);

    assert(abs(enoki::fmod( a,  b) - c)[0] < 1e-12f);
    assert(abs(enoki::fmod(-a,  b) + c)[0] < 1e-12f);
    assert(abs(enoki::fmod( a, -b) - c)[0] < 1e-12f);
    assert(abs(enoki::fmod(-a, -b) + c)[0] < 1e-12f);
}

ENOKI_TEST_FLOAT(test19_ceil2int) {
    T a = Value(-5.1);
    using Int = int_array_t<T>;
    assert(floor2int<Int>(a) == Int(-6));
    assert(ceil2int<Int>(a) == Int(-5));
}

ENOKI_TEST_FLOAT(test20_cbrt) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return cbrt(a); },
        [](double a) { return std::cbrt(a); },
        Value(-10), Value(10),
        3
    );
}
