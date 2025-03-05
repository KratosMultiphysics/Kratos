/*
    tests/conv.cpp -- tests value and mask conversion routines

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#if defined(__GNUG__)
#  pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "test.h"

template <typename T, typename Value2> void convtest() {
    using T2 = replace_scalar_t<T, Value2>;
    auto value1 = arange<T>();
    auto value2 = T2(value1);
    auto value3 = T(value2);
    assert(value1 == value3);
}

template <typename T, typename Value2> void masktest() {
    using Value = typename T::Value;
    using T2 = replace_scalar_t<T, Value2>;
    for (size_t i = 0; i < T::Size; ++i) {
        mask_t<T> mask = eq(arange<T>() - T(Value(i)), T(0));
        mask_t<T2> mask2(mask);
        T2 result = select(mask2, T2(Value2(1)), T2(Value2(0)));
        Value2 out[T::Size];
        store_unaligned(out, result);
        for (size_t j = 0; j < T::Size; ++j)
            assert(out[j] == ((j == i) ? Value2(1) : Value2(0)));
    }
}

ENOKI_TEST_ALL(test01_conv_int32_t)  { convtest<T, int32_t>();  }
ENOKI_TEST_ALL(test02_conv_uint32_t) { convtest<T, uint32_t>(); }
ENOKI_TEST_ALL(test03_conv_int64_t)  { convtest<T, int64_t>();  }
ENOKI_TEST_ALL(test04_conv_uint64_t) { convtest<T, uint64_t>(); }
ENOKI_TEST_ALL(test05_conv_half)     { convtest<T, half>();     }
ENOKI_TEST_ALL(test06_conv_float)    { convtest<T, float>();    }
ENOKI_TEST_ALL(test07_conv_double)   { convtest<T, double>();   }

ENOKI_TEST_ALL(test08_mask_int32_t)  { masktest<T, int32_t>();  }
ENOKI_TEST_ALL(test09_mask_uint32_t) { masktest<T, uint32_t>(); }
ENOKI_TEST_ALL(test10_mask_int64_t)  { masktest<T, int64_t>();  }
ENOKI_TEST_ALL(test11_mask_uint64_t) { masktest<T, uint64_t>(); }
ENOKI_TEST_ALL(test12_mask_float)    { masktest<T, float>();    }
ENOKI_TEST_ALL(test13_mask_double)   { masktest<T, double>();   }
ENOKI_TEST_ALL(test14_mask_half)     { masktest<T, half>();     }

ENOKI_TEST_ALL(test15_bool_conv) {
    for (size_t i = 0; i < T::Size; ++i) {
        mask_t<T> mask = eq(arange<T>() - T(Value(i)), T(0));
        bool_array_t<T> mask3(mask);
        mask_t<T> mask4(mask3);
        T result  = select(mask, T(Value(1)), T(Value(0)));
        T result2 = select(mask4, T(Value(1)), T(Value(0)));

        assert(result == result2);
    }
}
