/*
    tests/custom.cpp -- tests operations involving custom data structures

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/array.h>
#include <enoki/dynamic.h>

template <typename Value_> struct Custom {
    using Value = Value_;
    using FloatVector = Array<Value, 3>;
    using DoubleVector = float64_array_t<FloatVector>;
    using IntVector = int32_array_t<Value>;

    FloatVector o;
    DoubleVector d;
    IntVector i = 0;

    template <typename T>
    bool operator==(const Custom<T> &other) const {
        return other.o == o && other.d == d && other.i == i;
    }

    template <typename T>
    bool operator!=(const Custom<T> &other) const {
        return !operator==(other);
    }

    ENOKI_STRUCT(Custom, o, d, i)
};

ENOKI_STRUCT_SUPPORT(Custom, o, d, i)

ENOKI_TEST(test01_mask_slice_custom) {
    using FloatP = Packet<float>;
    using Vector3f = Array<float, 3>;
    using Vector3d = Array<double, 3>;
    using Custom3f = Custom<float>;
    using Custom3fP = Custom<FloatP>;

    Custom3fP x = zero<Custom3fP>();
    Custom3fP y;
    y.o = Vector3f(1, 2, 3);
    y.d = Vector3d(4, 5, 6);
    auto mask = arange<FloatP>() > 0.f;

    masked(x, mask) = y;

    assert((slice(x, 0) == Custom3f(Vector3f(0, 0, 0), Vector3f(0, 0, 0), 0)));
    if (FloatP::Size > 1)
        assert((slice(x, 1) == Custom3f(Vector3f(1, 2, 3), Vector3f(4, 5, 6), 0)));
}

ENOKI_TEST(test02_mask_slice_custom_scalar) {
    using Custom3f = Custom<float>;
    using Vector3f = Array<float, 3>;

    Custom3f x = zero<Custom3f>();
    Custom3f y(Vector3f(1, 2, 3), Vector3f(4, 5, 6), 0);
    Custom3f z = zero<Custom3f>();
    masked(z, true) = y;

    assert(y != x);
    assert(y == y);
}

struct Test { };

template <typename T> struct TrickyStruct {
    using Ptr = replace_scalar_t<T, Test *>;
    using Mask = mask_t<T>;

    Ptr ptr;
    Mask mask;

    ENOKI_STRUCT(TrickyStruct, ptr, mask);
};

ENOKI_STRUCT_SUPPORT(TrickyStruct, ptr, mask);

ENOKI_TEST(test03_tricky) {
    using FloatP = Packet<float>;
    using FloatX = DynamicArray<FloatP>;
    using Tricky = TrickyStruct<float>;
    using TrickyP = TrickyStruct<FloatP>;
    using TrickyX = TrickyStruct<FloatX>;

    TrickyP x;
    for (size_t i = 0; i<FloatP::Size; ++i)
        slice(x, i) = Tricky((Test *) (0xdeadbeef + i), (i & 1) != 0);

    for (size_t i = 0; i<FloatP::Size; ++i) {
        assert(x.mask.coeff(i) == ((i & 1) != 0));
        assert(x.ptr.coeff(i) == (Test *) (0xdeadbeef + i));
    }

    TrickyX xl;
    set_slices(xl, 3);
    xl.ptr = (Test *) 0x1337c0d3;
    auto pslice = slice_ptr(xl, 1);
    static_assert(std::is_same_v<decltype(pslice.ptr), Test**>);
    *pslice.ptr = (Test *) 0xdeadbeef;
    assert(slice(xl.ptr, 1) == (Test *) 0xdeadbeef);
}

ENOKI_TEST(test04_gather_custom_struct) {
    using FloatP = Packet<float>;
    using UInt32P = Packet<uint32_t>;
    using FloatX = DynamicArray<FloatP>;
    using UInt32X = DynamicArray<UInt32P>;
    using Custom3f = Custom<float>;
    using Custom3fP = Custom<FloatP>;
    using Custom3fX = Custom<FloatX>;

    Custom3fX z;
    z.o.x() = arange<FloatX>(20) + 1.f;
    z.o.y() = arange<FloatX>(20) + 100.f;
    z.o.z() = arange<FloatX>(20) + 1000.f;
    z.d.x() = arange<FloatX>(20) + 2.f;
    z.d.y() = arange<FloatX>(20) + 200.f;
    z.d.z() = arange<FloatX>(20) + 2000.f;
    z.i = arange<UInt32X>(20) + 1234u;

    Custom3fP p = gather<Custom3fP>(z, arange<UInt32P>() + 1u);

    assert(p.o.x() == arange<FloatP>() + 2.f);
    assert(p.o.y() == arange<FloatP>() + 101.f);
    assert(p.o.z() == arange<FloatP>() + 1001.f);

    assert(p.d.x() == arange<FloatP>() + 3.f);
    assert(p.d.y() == arange<FloatP>() + 201.f);
    assert(p.d.z() == arange<FloatP>() + 2001.f);

    assert(p.i == arange<UInt32P>() + 1235u);

    Custom3f s = gather<Custom3f>(z, 1u);

    assert(s.o.x() == 2.f);
    assert(s.o.y() == 101.f);
    assert(s.o.z() == 1001.f);

    assert(s.d.x() == 3.f);
    assert(s.d.y() == 201.f);
    assert(s.d.z() == 2001.f);

    assert(s.i == 1235u);
}
