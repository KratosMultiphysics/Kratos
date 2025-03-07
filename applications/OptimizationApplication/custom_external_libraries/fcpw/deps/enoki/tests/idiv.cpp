/*
    tests/idiv.cpp -- tests integer division by constants

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <random>

#define ITERATIONS 1000000

ENOKI_TEST(test01_idiv_u64) {
    std::mt19937_64 mt;

    for (uint64_t i = 2; i < ITERATIONS; ++i) {
        uint64_t x = (uint64_t) mt(), y = (uint64_t) mt(), z = y / x;

        divisor<uint64_t> precomp(x);
        uint64_t q = precomp(y);
        if (q != z)
            std::cout << y << " / " << x << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);

        divisor<uint64_t> precomp2(i);
        q = precomp2(y);
        z = y / i;
        if (q != z)
            std::cout << y << " / " << i << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);
    }
}

ENOKI_TEST(test02_idiv_u32) {
    std::mt19937 mt;

    for (uint32_t i = 2; i < ITERATIONS; ++i) {
        uint32_t x = (uint32_t) mt(), y = (uint32_t) mt(), z = y / x;

        divisor<uint32_t> precomp(x);
        uint32_t q = precomp(y);
        if (q != z)
            std::cout << y << " / " << x << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);

        divisor<uint32_t> precomp2(i);
        q = precomp2(y);
        z = y / i;
        if (q != z)
            std::cout << y << " / " << i << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);
    }
}

ENOKI_TEST(test03_idiv_s64) {
    std::mt19937_64 mt;

    for (uint64_t i = 2; i < ITERATIONS; ++i) {
        int64_t x = (int64_t) mt(), y = (int64_t) mt(), z = y / x;

        divisor<int64_t> precomp(x);
        int64_t q = precomp(y);
        if (q != z)
            std::cout << y << " / " << x << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);

        divisor<int64_t> precomp2((int64_t) i);
        q = precomp2(y);
        z = y / (int64_t) i;
        if (q != z)
            std::cout << y << " / " << i << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);

        divisor<int64_t> precomp3(-(int64_t) i);
        q = precomp3(y);
        z = y / -(int64_t) i;
        if (q != z)
            std::cout << y << " / " << i << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);
    }
}

ENOKI_TEST(test03_idiv_s32) {
    std::mt19937 mt;

    for (uint32_t i = 2; i < ITERATIONS; ++i) {
        int32_t x = (int32_t) mt(), y = (int32_t) mt(), z = y / x;

        divisor<int32_t> precomp(x);
        int32_t q = precomp(y);
        if (q != z)
            std::cout << y << " / " << x << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);

        divisor<int32_t> precomp2((int32_t) i);
        q = precomp2(y);
        z = y / (int32_t) i;
        if (q != z)
            std::cout << y << " / " << i << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);

        divisor<int32_t> precomp3(-(int32_t) i);
        q = precomp3(y);
        z = y / -(int32_t) i;
        if (q != z)
            std::cout << y << " / " << i << " = " << q << " vs " << z
                      << std::endl;
        assert(q == z);
    }
}

ENOKI_TEST_INT(test04_idiv_vector) {
    std::mt19937_64 mt;
    for (Value i = 2; i < 1000; ++i) {
        Value x = (Value) mt(), y = (Value) mt();
        assert((T(y) / x)[0] == y / x);
        assert((T(y) / i)[0] == y / i);
        assert((T(y) % x)[0] == y % x);
        assert((T(y) % i)[0] == y % i);
    }
}
