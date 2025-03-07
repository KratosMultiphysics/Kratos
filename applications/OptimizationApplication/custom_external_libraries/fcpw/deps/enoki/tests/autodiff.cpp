/*
    tests/autodiff.cpp -- tests the reverse-mode automatic differentation layer

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include "test.h"
#include <enoki/dynamic.h>
#include <enoki/autodiff.h>
#include <enoki/color.h>

using Float  = float;
using FloatP = Packet<Float>;
using FloatX = DynamicArray<FloatP>;
using FloatD = DiffArray<FloatX>;

using UInt32  = uint32_t;
using UInt32P = Packet<UInt32>;
using UInt32X = DynamicArray<UInt32P>;
using UInt32D = DiffArray<UInt32X>;

using Vector2f  = Array<Float, 2>;
using Vector2fP = Array<FloatP, 2>;
using Vector2fX = Array<FloatX, 2>;
using Vector2fD = Array<FloatD, 2>;

using Vector3f  = Array<Float, 3>;
using Vector3fP = Array<FloatP, 3>;
using Vector3fX = Array<FloatX, 3>;
using Vector3fD = Array<FloatD, 3>;

template <typename T> void my_backward(T &t, bool clear_graph = true) {
    FloatD::simplify_graph_();
    backward(t, clear_graph);
}

template <typename T> void my_forward(T &t, bool clear_graph = true) {
    FloatD::simplify_graph_();
    forward(t, clear_graph);
}

ENOKI_TEST(test00_identity) {
    FloatD x = 2.f;
    set_requires_gradient(x);
    my_backward(x);
    assert(gradient(x)[0] == 1.f);
}

ENOKI_TEST(test00_identity_fwd) {
    FloatD x = 2.f;
    set_requires_gradient(x);
    my_forward(x);
    assert(gradient(x)[0] == 1.f);
}

ENOKI_TEST(test01_add) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x + y;
    my_backward(z);
    assert(gradient(x)[0] == 1.f);
    assert(gradient(y)[0] == 1.f);
}

ENOKI_TEST(test01_add_fwd) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x + y;
    my_forward(x);
    assert(gradient(z)[0] == 1.f);
}

ENOKI_TEST(test02_sub) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x - y;
    my_backward(z);
    assert(gradient(x)[0] == 1.f);
    assert(gradient(y)[0] == -1.f);
}

ENOKI_TEST(test02_sub_fwd) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x - y;
    my_forward(x, false);
    assert(gradient(z)[0] == 1.f);
    my_forward(y);
    assert(gradient(z)[0] == -1.f);
}

ENOKI_TEST(test03_mul) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x * y;
    my_backward(z);
    assert(gradient(x)[0] == 3.f);
    assert(gradient(y)[0] == 2.f);
}

ENOKI_TEST(test03_mul_fwd) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x * y;
    my_forward(x, false);
    assert(gradient(z)[0] == 3.f);
    my_forward(y);
    assert(gradient(z)[0] == 2.f);
}

ENOKI_TEST(test04_div) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x / y;
    my_backward(z);
    assert(std::abs(gradient(x)[0] - 1.f / 3.f) < 1e-6f);
    assert(std::abs(gradient(y)[0] + 2.f / 9.f) < 1e-6f);
}

ENOKI_TEST(test04_div_fwd) {
    FloatD x = 2.f, y = 3.f;
    set_requires_gradient(x);
    set_requires_gradient(y);
    FloatD z = x / y;
    my_forward(x, false);
    assert(std::abs(gradient(z)[0] - 1.f / 3.f) < 1e-6f);
    my_forward(y);
    assert(std::abs(gradient(z)[0] + 2.f / 9.f) < 1e-6f);
}

ENOKI_TEST(test05_hsum_0) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = hsum(x*x);
    my_backward(y);
    assert(y.size() == 1 && allclose(y.coeff(0), 95.f/27.f));
    assert(allclose(gradient(x), 2.f * detach(x)));
}

ENOKI_TEST(test05_hsum_0_fwd) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = hsum(x*x);
    my_forward(x);
    assert(y.size() == 1 && allclose(y.coeff(0), 95.f/27.f));
    assert(allclose(gradient(y), 10));
}

ENOKI_TEST(test05_hsum_1) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 11);
    set_requires_gradient(x);
    FloatD z = hsum(hsum(x)*x);
    my_backward(z);
    assert(gradient(x) == 11.f);
}

ENOKI_TEST(test05_hsum_1_fwd) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = hsum(hsum(x)*x);
    my_forward(x);
    assert(allclose(gradient(y), 100));
}

ENOKI_TEST(test05_hsum_2) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 11);
    set_requires_gradient(x);
    FloatD z = hsum(hsum(x*x)*x*x);
    my_backward(z);

    assert(allclose(gradient(x),
        FloatX(0.f, 1.54f, 3.08f, 4.62f, 6.16f, 7.7f, 9.24f, 10.78f, 12.32f, 13.86f, 15.4f)));
}

ENOKI_TEST(test05_hsum_2_fwd) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = hsum(hsum(x*x)*hsum(x*x));
    my_forward(x);
    assert(allclose(gradient(y), 1900.f/27.f));
}

ENOKI_TEST(test06_hprod) {
    FloatD x = linspace<FloatD>(1.f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = hprod(x);
    my_backward(y);
    assert(allclose(gradient(x), hprod(detach(x)) / detach(x)) &&
           y.size() == 1 && allclose(y.coeff(0), 45.5402f));
}

ENOKI_TEST(test07_sqrt) {
    FloatD x = linspace<FloatD>(1.f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = sqrt(x);
    my_backward(y);
    assert(allclose(y, sqrt(x)));
    assert(allclose(gradient(x), .5f*rsqrt(detach(x))));
}

ENOKI_TEST(test08_rsqrt) {
    FloatD x = linspace<FloatD>(1.f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = rsqrt(x);
    my_backward(y);
    assert(allclose(y, rsqrt(x)));
    assert(allclose(gradient(x), -.5f * pow(detach(x), -3.f / 2.f)));
}

ENOKI_TEST(test09_exp) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = exp(x*x);
    my_backward(y);
    assert(allclose(y, exp(sqr(detach(x)))));
    assert(allclose(gradient(x), 2.f * detach(x) * exp(sqr(detach(x)))));
}

ENOKI_TEST(test10_log) {
    FloatD x = linspace<FloatD>(1.f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = log(x*x);
    my_backward(y);
    assert(allclose(y, log(sqr(detach(x)))));
    assert(allclose(gradient(x), 2.f / detach(x)));
}

ENOKI_TEST(test11_sin) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = sin(x*x);
    my_backward(y);
    assert(allclose(y, sin(sqr(detach(x)))));
    assert(allclose(gradient(x), 2*detach(x)*cos(sqr(detach(x)))));
}

ENOKI_TEST(test12_cos) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = cos(x*x);
    my_backward(y);
    assert(allclose(y, cos(sqr(detach(x)))));
    assert(allclose(gradient(x), -2*detach(x)*sin(sqr(detach(x)))));
}

ENOKI_TEST(test13_tan) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = tan(x*x);
    my_backward(y);
    assert(allclose(y, tan(sqr(detach(x)))));
    assert(allclose(gradient(x), 2*detach(x)*sqr(sec(sqr(detach(x))))));
}

ENOKI_TEST(test14_csc) {
    FloatD x = linspace<FloatD>(1.f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = csc(x*x);
    my_backward(y);
    assert(allclose(y, csc(sqr(detach(x)))));
    assert(allclose(gradient(x), -2.f * detach(x) * cot(sqr(detach(x))) *
                                     csc(sqr(detach(x))), 1e-4f, 1e-4f));
}

ENOKI_TEST(test15_sec) {
    FloatD x = linspace<FloatD>(1.f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = cot(x*x);
    my_backward(y);
    assert(allclose(y, cot(sqr(detach(x)))));
    assert(allclose(gradient(x), -2 * detach(x) * sqr(csc(sqr(detach(x))))));
}

ENOKI_TEST(test17_asin) {
    FloatD x = linspace<FloatD>(-.8f, .8f, 10);
    set_requires_gradient(x);
    FloatD y = asin(x*x);
    my_backward(y);
    assert(allclose(y, asin(sqr(detach(x)))));
    assert(allclose(gradient(x), 2*detach(x) / sqrt(1-sqr(sqr(detach(x))))));
}

ENOKI_TEST(test18_acos) {
    FloatD x = linspace<FloatD>(-.8f, .8f, 10);
    set_requires_gradient(x);
    FloatD y = acos(x*x);
    my_backward(y);
    assert(allclose(y, acos(sqr(detach(x)))));
    assert(allclose(gradient(x), -2*detach(x) / sqrt(1-sqr(sqr(detach(x))))));
}

ENOKI_TEST(test19_atan) {
    FloatD x = linspace<FloatD>(-.8f, .8f, 10);
    set_requires_gradient(x);
    FloatD y = atan(x*x);
    my_backward(y);
    assert(allclose(y, atan(sqr(detach(x)))));
    assert(allclose(gradient(x), 2.f*detach(x) / (1.f+sqr(sqr(detach(x))))));
}

ENOKI_TEST(test20_sinh) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = sinh(x*x);
    my_backward(y);
    assert(allclose(y, sinh(sqr(detach(x)))));
    assert(allclose(gradient(x), 2*detach(x) * cosh(sqr(detach(x)))));
}

ENOKI_TEST(test21_cosh) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = cosh(x*x);
    my_backward(y);
    assert(allclose(y, cosh(sqr(detach(x)))));
    assert(allclose(gradient(x), 2*detach(x) * sinh(sqr(detach(x)))));
}

ENOKI_TEST(test22_tanh) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = tanh(x*x);
    my_backward(y);
    assert(allclose(y, tanh(sqr(detach(x)))));
    assert(allclose(gradient(x), 2*detach(x) * sqr(sech(sqr(detach(x))))));
}

ENOKI_TEST(test23_csch) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = csch(x*x);
    my_backward(y);
    assert(allclose(y, csch(sqr(detach(x)))));
    assert(allclose(gradient(x), -2*detach(x) * csch(sqr(detach(x))) * coth(sqr(detach(x)))));
}

ENOKI_TEST(test24_sech) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = sech(x*x);
    my_backward(y);
    assert(allclose(y, sech(sqr(detach(x)))));
    assert(allclose(gradient(x), -2.f*detach(x) * sech(sqr(detach(x))) * tanh(sqr(detach(x)))));
}

ENOKI_TEST(test25_coth) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = asinh(x*x);
    my_backward(y);
    assert(allclose(y, asinh(sqr(detach(x)))));
    assert(allclose(gradient(x), 2.f*detach(x) * rsqrt(1.f + sqr(sqr(detach(x))))));
}

ENOKI_TEST(test26_acosh) {
    FloatD x = linspace<FloatD>(1.01f, 2.f, 10);
    set_requires_gradient(x);
    FloatD y = acosh(x*x);
    my_backward(y);
    assert(allclose(y, acosh(sqr(detach(x)))));
    assert(allclose(gradient(x), 2.f*detach(x) * rsqrt(sqr(sqr(detach(x))) - 1.f)));
}

ENOKI_TEST(test27_atanh) {
    FloatD x = linspace<FloatD>(-.99f, .99f, 10);
    set_requires_gradient(x);
    FloatD y = atanh(x*x);
    my_backward(y);
    assert(allclose(y, atanh(sqr(detach(x)))));
    assert(allclose(gradient(x), -2.f*detach(x) * rcp(sqr(sqr(detach(x))) - 1.f)));
}

ENOKI_TEST(test28_linear_to_srgb) {
    FloatD x = linspace<FloatD>(0.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = linear_to_srgb(x);

    std::cout << graphviz(y) << std::endl;
    my_backward(y);
    /// from mathematica
    FloatX ref{ 12.92f,   1.58374f,  1.05702f, 0.834376f, 0.705474f,
                0.61937f, 0.556879f, 0.50899f, 0.470847f, 0.439583f };
    assert(hmax(abs(gradient(x) - ref)) < 1e-5f);
}

ENOKI_TEST(test29_scatter_add) {
    UInt32D idx1 = arange<UInt32D>(5);
    UInt32D idx2 = arange<UInt32D>(4)+3u;

    FloatD x = linspace<FloatD>(0, 1, 5);
    FloatD y = linspace<FloatD>(1, 2, 4);

    set_requires_gradient(x);
    set_requires_gradient(y);
    set_label(x, "x");
    set_label(y, "y");

    FloatD buf = zero<FloatD>(10);
    scatter_add(buf, x, idx1);
    scatter_add(buf, y, idx2);

    FloatD ref_buf { 0.0000f, 0.2500f, 0.5000f, 1.7500f, 2.3333f,
                     1.6667f, 2.0000f, 0.0000f, 0.0000f, 0.0000f };

    assert(allclose(ref_buf, buf, 1e-4f, 1e-4f));

    FloatD s = dot(buf, buf);
    std::cout << graphviz(s) << std::endl;

    my_backward(s);

    FloatD ref_x {0.0000f, 0.5000f, 1.0000f, 3.5000f, 4.6667f};
    FloatD ref_y {3.5000f, 4.6667f, 3.3333f, 4.0000f};

    assert(allclose(gradient(y), detach(ref_y), 1e-4f, 1e-4f));
    assert(allclose(gradient(x), detach(ref_x), 1e-4f, 1e-4f));
}

ENOKI_TEST(test30_scatter) {
    UInt32D idx1 = arange<UInt32D>(5);
    UInt32D idx2 = arange<UInt32D>(4)+3u;

    FloatD x = linspace<FloatD>(0, 1, 5);
    FloatD y = linspace<FloatD>(1, 2, 4);

    set_requires_gradient(x);
    set_requires_gradient(y);
    set_label(x, "x");
    set_label(y, "y");

    FloatD buf = zero<FloatD>(10);
    scatter(buf, x, idx1);
    if constexpr (is_cuda_array_v<FloatD>)
        cuda_eval();
    scatter(buf, y, idx2);

    FloatD ref_buf{ 0.0000f, 0.2500f, 0.5000f, 1.0000f, 1.3333f,
                    1.6667f, 2.0000f, 0.0000f, 0.0000f, 0.0000f };

    assert(allclose(ref_buf, buf, 1e-4f, 1e-4f));

    FloatD s = dot(buf, buf);
    std::cout << graphviz(s) << std::endl;

    my_backward(s);

    FloatD ref_x{ 0.0000f, 0.5000f, 1.0000f, 0.0000f, 0.0000f };
    FloatD ref_y{ 2.0000f, 2.6667f, 3.3333f, 4.0000f };

    assert(allclose(gradient(y), detach(ref_y), 1e-4f, 1e-4f));
    assert(allclose(gradient(x), detach(ref_x), 1e-4f, 1e-4f));
}

template <typename Vector2>
Vector2 square_to_uniform_disk_concentric(Vector2 sample) {
    using Value  = value_t<Vector2>;
    using Mask   = mask_t<Value>;
    using Scalar = scalar_t<Value>;

    Value x = fmsub(Scalar(2), sample.x(), Scalar(1)),
          y = fmsub(Scalar(2), sample.y(), Scalar(1));

    Mask is_zero         = eq(x, zero<Value>()) &&
                           eq(y, zero<Value>()),
         quadrant_1_or_3 = abs(x) < abs(y);

    Value r  = select(quadrant_1_or_3, y, x),
          rp = select(quadrant_1_or_3, x, y);

    Value phi = rp / r * Scalar(.25f * Float(M_PI));
    masked(phi, quadrant_1_or_3) = Scalar(.5f * Float(M_PI)) - phi;
    masked(phi, is_zero) = zero<Value>();

    auto [sin_phi, cos_phi] = sincos(phi);

    return Vector2(r * cos_phi, r * sin_phi);
}


/// Warp a uniformly distributed square sample to a Beckmann distribution
template <typename Vector2,
          typename Value   = expr_t<value_t<Vector2>>,
          typename Vector3 = Array<Value, 3>>
Vector3 square_to_beckmann(const Vector2 &sample, const value_t<Vector2> &alpha) {
    Vector2 p = square_to_uniform_disk_concentric(sample);
    Value r2 = squared_norm(p);

    Value tan_theta_m_sqr = -alpha * alpha * log(1 - r2);
    Value cos_theta_m = rsqrt(1 + tan_theta_m_sqr);
    p *= safe_sqrt((1 - cos_theta_m * cos_theta_m) / r2);

    return Vector3(p.x(), p.y(), cos_theta_m);
}

ENOKI_TEST(test32_sample_disk) {
    Vector2f x(.2f, .3f);

    Vector2fD y(x);
    set_requires_gradient(y);

    Vector3fD z = square_to_beckmann(y, .4f);
    Vector3f z_ref (-0.223574f, -0.12908f, 0.966102f);
    assert(allclose(detach(z), z_ref));

    auto sum = hsum(z);
    std::cout << graphviz(sum) << std::endl;
    my_backward(sum);
    Float eps = 1e-3f;
    Float dx = hsum(square_to_beckmann(x + Vector2f(eps, 0.f), .4f) -
                    square_to_beckmann(x - Vector2f(eps, 0.f), .4f)) /
                        (2 * eps);
    Float dy = hsum(square_to_beckmann(x + Vector2f(0.f, eps), .4f) -
                    square_to_beckmann(x - Vector2f(0.f, eps), .4f)) /
                        (2 * eps);

    assert(allclose(gradient(y), Vector2f(dx, dy), 1e-3f, 1e-3f));
}

ENOKI_TEST(test33_bcast) {
    FloatD x(5.f);
    FloatD y = arange<FloatD>(10);

    set_requires_gradient(x);
    set_requires_gradient(y);
    set_label(x, "x");
    set_label(y, "y");
    FloatD z = hsum(sqr(sin(x)*cos(y)));
    my_backward(z);

    assert(allclose(gradient(x), -2.8803, 1e-4f, 1e-4f));
    assert(allclose(gradient(y),
                    FloatX(-0.0000, -0.8361, 0.6959, 0.2569, -0.9098, 0.5002,
                           0.4934, -0.9109, 0.2647, 0.6906), 1e-4f, 1e-4f));
}

ENOKI_TEST(test34_gradient_descent) {
    FloatD x = zero<FloatD>(10);
    set_label(x, "x");
    float loss_f = 0.f;
    for (size_t i = 0; i < 10; ++i) {
        set_requires_gradient(x);
        FloatD loss = norm(x - linspace<FloatD>(0, 1, 10));
        my_backward(loss);
        x = detach(x) - gradient(x)*2e-1f;
        loss_f = detach(loss)[0];
    }
    assert(loss_f < 1e-1f);
}

struct Function {
    virtual FloatD eval(const FloatD &x) const = 0;
    virtual ~Function() = default;
};

struct Square : Function {
    FloatD eval(const FloatD &x) const override {
        return x*x;
    }
};

struct Reciprocal : Function {
    FloatD eval(const FloatD &x) const override {
        return 1.f / x;
    }
};

using FunctionP = Packet<Function *, FloatP::Size>;
using FunctionX = DynamicArray<FunctionP>;
using FunctionD = DiffArray<FunctionX>;

ENOKI_CALL_SUPPORT_BEGIN(Function)
ENOKI_CALL_SUPPORT_METHOD(eval)
ENOKI_CALL_SUPPORT_END(Function)

ENOKI_TEST(test35_call) {
    Function *square = new Square();
    Function *reciprocal = new Reciprocal();
    FunctionD f = select(arange<UInt32D>(10) < 4, FunctionD(square),
                         FunctionD(reciprocal));

    FloatD x = linspace<FloatD>(1, 2, 10);
    set_requires_gradient(x);

    FloatD out = f->eval(x);
    my_backward(out);

    FloatX ref_gradient{ 2.f, 2.22222f, 2.44444f, 2.66667f, -0.47929f,
                         -0.413265f, -0.36f, -0.316406f, -0.280277f, -0.25f };
    assert(allclose(ref_gradient, gradient(x), 1e-4f, 1e-4f));

    delete reciprocal;
    delete square;
}

ENOKI_TEST(test36_gather) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = gather<FloatD>(x*x, UInt32D(1, 2, 3));
    FloatD z = hsum(y);
    my_backward(z);
    FloatX ref_gradient { 0.f, -1.55556, -1.11111, -0.666667, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
    assert(allclose(ref_gradient, gradient(x), 1e-4f, 1e-4f));
}

ENOKI_TEST(test36_gather_fwd) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 10);
    set_requires_gradient(x);
    FloatD y = gather<FloatD>(x*x, UInt32D(1, 2, 3));
    my_forward(x);
    FloatX ref_gradient { -1.55556f, -1.11111f, -0.666667f };
    assert(allclose(ref_gradient, gradient(y), 1e-4f, 1e-4f));
}

ENOKI_TEST(test37_scatter_fwd) {
    FloatD x = linspace<FloatD>(-1.f, 1.f, 5);
    set_requires_gradient(x);
    FloatD y = zero<FloatD>(10);
    scatter(y, x*x, arange<UInt32D>(5) + 2);
    my_forward(x);
    FloatX ref_gradient { 0.f, 0.f, -2.f, -1.f, 0.f, 1.f, 2.f, 0.f, 0.f, 0.f };
    assert(allclose(ref_gradient, gradient(y), 1e-4f, 1e-4f));
}
