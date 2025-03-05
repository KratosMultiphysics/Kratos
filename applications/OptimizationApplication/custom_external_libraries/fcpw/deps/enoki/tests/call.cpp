/*
    tests/complex.cpp -- tests vectorized function calls

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/dynamic.h>
#include "ray.h"
#include <enoki/stl.h>

struct Test;
struct TestChild;

using Int32P = Array<int>;
using TestP = Array<const Test*, Int32P::Size>;
using ChildP = Array<const TestChild*, Int32P::Size>;

using TestX  = DynamicArray<TestP>;
using Int32X = DynamicArray<Int32P>;

using TestPMask = mask_t<TestP>;
using TestXMask = mask_t<TestX>;

using FloatP    = Packet<float>;
using Vector3f  = Array<float, 3>;
using Vector3fP = Array<FloatP, 3>;
using Ray3fP    = Ray<Vector3fP>;


struct Test {
    ENOKI_CALL_SUPPORT_FRIEND()

    Test(int32_t value) : value(value) { }
    virtual ~Test() { }

    // Vectorized function (accepts a mask, which is ignored here)
    virtual Int32P func1(Int32P i, TestPMask /* unused */) const { return i + value; }
    virtual Int32X func1(Int32X i, TestXMask /* unused */) const { return i + value; }

    // Vectorized function (accepts a mask, which is ignored here)
    virtual void func2(Int32P &i, TestPMask mask) const { i[mask] += value; }

    bool func3() const { return value == 20; }

    std::pair<Int32P, Int32P> func4(TestPMask) const { return std::make_pair(value, value+1); }

    Ray3fP make_ray(TestPMask) const { return Ray3fP(Vector3f(1, 1, 1), Vector3f(1, 2, 3));}

protected:
    int32_t value;
};

struct TestChild : public Test {
    TestChild() : Test(42) { }

    bool is_child() const { return value == 42; }
};

// Allow Enoki arrays containing pointers to transparently forward function
// calls (with the appropriate masks).
ENOKI_CALL_SUPPORT_BEGIN(Test)
ENOKI_CALL_SUPPORT_METHOD(func1)
ENOKI_CALL_SUPPORT_METHOD(func2)
ENOKI_CALL_SUPPORT_METHOD(func3)
ENOKI_CALL_SUPPORT_GETTER(get_value, value)
ENOKI_CALL_SUPPORT_METHOD(func4)
ENOKI_CALL_SUPPORT_METHOD(make_ray)
ENOKI_CALL_SUPPORT_END(Test)

ENOKI_CALL_SUPPORT_BEGIN(TestChild)
ENOKI_CALL_SUPPORT_METHOD(is_child)
ENOKI_CALL_SUPPORT_END(TestChild)


ENOKI_TEST(test01_call) {
    size_t offset = std::min((size_t) 2, TestP::Size-1);
    Test *a = new Test(10);
    Test *b = new Test(20);

    TestP pointers(a);
    pointers.coeff(offset) = b;

    Int32P index = arange<Int32P>();
    Int32P ref = arange<Int32P>() + 10;
    if (offset < Int32P::Size)
        ref.coeff(offset) += 10;

    Int32P result = pointers->func1(index);
    assert(result == ref);

    Int32P ref2 = 10;
    if (offset < Int32P::Size)
        ref2.coeff(offset) += 10;

    assert(pointers->get_value() == ref2);

    std::pair<Int32P, Int32P> result2 = pointers->func4();
    assert(result2.first == ref2);
    assert(result2.second == ref2+1);

    TestX pointers_x;
    Int32X index_x;
    set_slices(pointers_x, TestP::Size);
    set_slices(index_x, TestP::Size);
    packet(pointers_x, 0) = pointers;
    packet(index_x, 0) = index;
    Int32X result_x = pointers_x->func1(index_x);
    assert(packet(result_x, 0) == ref);

    pointers->func2(index);
    assert(index == ref);

    auto mask = mask_t<TestP>(pointers->func3());
    assert(mask == eq(pointers, b));

    /* The following should not crash */
    pointers.coeff(0) = nullptr;
    pointers->func3();

    delete a;
    delete b;
}


ENOKI_TEST(test02_reinterpret_pointer_array) {
    using Mask = mask_t<ChildP>;
    Test *a = new Test(1);
    Test *b = new TestChild();

    TestP objects(b);
    objects[std::min((size_t) 2, TestP::Size-1)] = a;

    auto children = reinterpret_array<ChildP>(objects);
    // is_child returns an Array of bools, need to cast to a mask type for the
    // comparison to be correct.
    assert(all(Mask(children->is_child()) == eq(objects, b)));

    delete a;
    delete b;
}

ENOKI_TEST(test03_call_with_structure) {
    Test *a = new Test(1);
    TestP objects(a);
    Vector3fP t = objects->make_ray()(1);
    assert(all_nested(eq(t, Vector3f(2, 3, 4))));
    delete a;
}
