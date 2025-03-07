/*
    tests/morton.cpp -- tests Morton/Z-order curve encoding and decoding

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/morton.h>

ENOKI_TEST(test01_morton_u32_2d_scalar) {
    using T = uint32_t;
    using T2 = enoki::Array<T, 2>;

    T2 value = T2(123u, 456u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test02_morton_u32_2d_array) {
    using T = enoki::Array<uint32_t>;
    using T2 = enoki::Array<T, 2>;

    T2 value = T2(123u, 456u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test03_morton_u32_3d_scalar) {
    using T = uint32_t;
    using T2 = enoki::Array<T, 3>;

    T2 value = T2(123u, 456u, 789u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test04_morton_u32_3d_array) {
    using T = enoki::Array<uint32_t>;
    using T2 = enoki::Array<T, 3>;

    T2 value = T2(123u, 456u, 789u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test05_morton_u64_2d_scalar) {
    using T = uint64_t;
    using T2 = enoki::Array<T, 2>;

    T2 value = T2(123u, 456u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test06_morton_u64_2d_array) {
    using T = enoki::Array<uint64_t>;
    using T2 = enoki::Array<T, 2>;

    T2 value = T2(123u, 456u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test07_morton_u64_3d_scalar) {
    using T = uint64_t;
    using T2 = enoki::Array<T, 3>;

    T2 value = T2(123u, 456u, 789u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}

ENOKI_TEST(test08_morton_u64_3d_array) {
    using T = enoki::Array<uint64_t>;
    using T2 = enoki::Array<T, 3>;

    T2 value = T2(123u, 456u, 789u);
    T value2 = morton_encode(value);
    T2 value3 = morton_decode<T2>(value2);

    assert(value == value3);
}
