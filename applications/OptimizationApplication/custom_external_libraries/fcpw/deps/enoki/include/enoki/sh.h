/*
    enoki/matrix.h -- Real spherical harmonics evaluation routines

    The generated code is based on the paper `Efficient Spherical Harmonic
    Evaluation, Journal of Computer Graphics Techniques (JCGT), vol. 2, no. 2,
    84-90, 2013 by Peter-Pike Sloan

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include "array.h"

NAMESPACE_BEGIN(enoki)

template <typename Array>
void sh_eval(const Array &d, size_t order, value_t<expr_t<Array>> *out) {
    switch (order) {
        case 0: sh_eval_0(d, out); break;
        case 1: sh_eval_1(d, out); break;
        case 2: sh_eval_2(d, out); break;
        case 3: sh_eval_3(d, out); break;
        case 4: sh_eval_4(d, out); break;
        case 5: sh_eval_5(d, out); break;
        case 6: sh_eval_6(d, out); break;
        case 7: sh_eval_7(d, out); break;
        case 8: sh_eval_8(d, out); break;
        case 9: sh_eval_9(d, out); break;
        default: throw std::runtime_error("sh_eval(): order too high!");
    }
}

template <typename Array>
void sh_eval_0(const Array &, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    store(out + 0, Value(Scalar(0.28209479177387814)));
}

template <typename Array>
void sh_eval_1(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z();
    Value c0, s0, tmp_a;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
}

template <typename Array>
void sh_eval_2(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_c = Scalar(0.546274215296039478);
    store(out + 8, tmp_c * c1);
    store(out + 4, tmp_c * s1);
}

template <typename Array>
void sh_eval_3(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_c = Scalar(-0.590043589926643519);
    store(out + 15, tmp_c * c0);
    store(out + 9, tmp_c * s0);
}

template <typename Array>
void sh_eval_4(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    store(out + 20, fmadd(z * Scalar(1.98431348329844304), load<Value>(out + 12), load<Value>(out + 6) * Scalar(-1.00623058987490532)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-4.6833258049010249), Scalar(2.00713963067186763));
    store(out + 21, tmp_a * c0);
    store(out + 19, tmp_a * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(3.31161143515146028), Scalar(-0.473087347878779985));
    store(out + 22, tmp_c * c1);
    store(out + 18, tmp_c * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.590043589926643519);
    store(out + 15, tmp_a * c0);
    store(out + 9, tmp_a * s0);
    tmp_b = z * Scalar(-1.77013076977993067);
    store(out + 23, tmp_b * c0);
    store(out + 17, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_c = Scalar(0.625835735449176256);
    store(out + 24, tmp_c * c1);
    store(out + 16, tmp_c * s1);
}

template <typename Array>
void sh_eval_5(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    store(out + 20, fmadd(z * Scalar(1.98431348329844304), load<Value>(out + 12), load<Value>(out + 6) * Scalar(-1.00623058987490532)));
    store(out + 30, fmadd(z * Scalar(1.98997487421323993), load<Value>(out + 20), load<Value>(out + 12) * Scalar(-1.00285307284481395)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-4.6833258049010249), Scalar(2.00713963067186763));
    store(out + 21, tmp_a * c0);
    store(out + 19, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.03100960115899021), tmp_a, tmp_c * Scalar(-0.991031208965114985));
    store(out + 31, tmp_b * c0);
    store(out + 29, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(3.31161143515146028), Scalar(-0.473087347878779985));
    store(out + 22, tmp_c * c1);
    store(out + 18, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(7.19030517745998665), Scalar(-2.39676839248666207));
    store(out + 32, tmp_a * c1);
    store(out + 28, tmp_a * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.590043589926643519);
    store(out + 15, tmp_a * c0);
    store(out + 9, tmp_a * s0);
    tmp_b = z * Scalar(-1.77013076977993067);
    store(out + 23, tmp_b * c0);
    store(out + 17, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-4.40314469491725369), Scalar(0.48923829943525049));
    store(out + 33, tmp_c * c0);
    store(out + 27, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.625835735449176256);
    store(out + 24, tmp_a * c1);
    store(out + 16, tmp_a * s1);
    tmp_b = z * Scalar(2.07566231488104114);
    store(out + 34, tmp_b * c1);
    store(out + 26, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_c = Scalar(-0.656382056840170258);
    store(out + 35, tmp_c * c0);
    store(out + 25, tmp_c * s0);
}

template <typename Array>
void sh_eval_6(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    store(out + 20, fmadd(z * Scalar(1.98431348329844304), load<Value>(out + 12), load<Value>(out + 6) * Scalar(-1.00623058987490532)));
    store(out + 30, fmadd(z * Scalar(1.98997487421323993), load<Value>(out + 20), load<Value>(out + 12) * Scalar(-1.00285307284481395)));
    store(out + 42, fmadd(z * Scalar(1.99304345718356646), load<Value>(out + 30), load<Value>(out + 20) * Scalar(-1.00154202096221923)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-4.6833258049010249), Scalar(2.00713963067186763));
    store(out + 21, tmp_a * c0);
    store(out + 19, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.03100960115899021), tmp_a, tmp_c * Scalar(-0.991031208965114985));
    store(out + 31, tmp_b * c0);
    store(out + 29, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.02131498923702768), tmp_b, tmp_a * Scalar(-0.995226703056238504));
    store(out + 43, tmp_c * c0);
    store(out + 41, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(3.31161143515146028), Scalar(-0.473087347878779985));
    store(out + 22, tmp_c * c1);
    store(out + 18, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(7.19030517745998665), Scalar(-2.39676839248666207));
    store(out + 32, tmp_a * c1);
    store(out + 28, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.11394181566096995), tmp_a, tmp_c * Scalar(-0.973610120462326756));
    store(out + 44, tmp_b * c1);
    store(out + 40, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.590043589926643519);
    store(out + 15, tmp_a * c0);
    store(out + 9, tmp_a * s0);
    tmp_b = z * Scalar(-1.77013076977993067);
    store(out + 23, tmp_b * c0);
    store(out + 17, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-4.40314469491725369), Scalar(0.48923829943525049));
    store(out + 33, tmp_c * c0);
    store(out + 27, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-10.1332578546641603), Scalar(2.76361577854477058));
    store(out + 45, tmp_a * c0);
    store(out + 39, tmp_a * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.625835735449176256);
    store(out + 24, tmp_a * c1);
    store(out + 16, tmp_a * s1);
    tmp_b = z * Scalar(2.07566231488104114);
    store(out + 34, tmp_b * c1);
    store(out + 26, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(5.55021390801596581), Scalar(-0.504564900728724064));
    store(out + 46, tmp_c * c1);
    store(out + 38, tmp_c * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.656382056840170258);
    store(out + 35, tmp_a * c0);
    store(out + 25, tmp_a * s0);
    tmp_b = z * Scalar(-2.3666191622317525);
    store(out + 47, tmp_b * c0);
    store(out + 37, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_c = Scalar(0.683184105191914415);
    store(out + 48, tmp_c * c1);
    store(out + 36, tmp_c * s1);
}

template <typename Array>
void sh_eval_7(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    store(out + 20, fmadd(z * Scalar(1.98431348329844304), load<Value>(out + 12), load<Value>(out + 6) * Scalar(-1.00623058987490532)));
    store(out + 30, fmadd(z * Scalar(1.98997487421323993), load<Value>(out + 20), load<Value>(out + 12) * Scalar(-1.00285307284481395)));
    store(out + 42, fmadd(z * Scalar(1.99304345718356646), load<Value>(out + 30), load<Value>(out + 20) * Scalar(-1.00154202096221923)));
    store(out + 56, fmadd(z * Scalar(1.99489143482413467), load<Value>(out + 42), load<Value>(out + 30) * Scalar(-1.00092721392195827)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-4.6833258049010249), Scalar(2.00713963067186763));
    store(out + 21, tmp_a * c0);
    store(out + 19, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.03100960115899021), tmp_a, tmp_c * Scalar(-0.991031208965114985));
    store(out + 31, tmp_b * c0);
    store(out + 29, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.02131498923702768), tmp_b, tmp_a * Scalar(-0.995226703056238504));
    store(out + 43, tmp_c * c0);
    store(out + 41, tmp_c * s0);
    tmp_a = fmadd(z * Scalar(2.01556443707463773), tmp_c, tmp_b * Scalar(-0.99715504402183186));
    store(out + 57, tmp_a * c0);
    store(out + 55, tmp_a * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(3.31161143515146028), Scalar(-0.473087347878779985));
    store(out + 22, tmp_c * c1);
    store(out + 18, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(7.19030517745998665), Scalar(-2.39676839248666207));
    store(out + 32, tmp_a * c1);
    store(out + 28, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.11394181566096995), tmp_a, tmp_c * Scalar(-0.973610120462326756));
    store(out + 44, tmp_b * c1);
    store(out + 40, tmp_b * s1);
    tmp_c = fmadd(z * Scalar(2.08166599946613307), tmp_b, tmp_a * Scalar(-0.984731927834661791));
    store(out + 58, tmp_c * c1);
    store(out + 54, tmp_c * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.590043589926643519);
    store(out + 15, tmp_a * c0);
    store(out + 9, tmp_a * s0);
    tmp_b = z * Scalar(-1.77013076977993067);
    store(out + 23, tmp_b * c0);
    store(out + 17, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-4.40314469491725369), Scalar(0.48923829943525049));
    store(out + 33, tmp_c * c0);
    store(out + 27, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-10.1332578546641603), Scalar(2.76361577854477058));
    store(out + 45, tmp_a * c0);
    store(out + 39, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.20794021658196149), tmp_a, tmp_c * Scalar(-0.95940322360024699));
    store(out + 59, tmp_b * c0);
    store(out + 53, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.625835735449176256);
    store(out + 24, tmp_a * c1);
    store(out + 16, tmp_a * s1);
    tmp_b = z * Scalar(2.07566231488104114);
    store(out + 34, tmp_b * c1);
    store(out + 26, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(5.55021390801596581), Scalar(-0.504564900728724064));
    store(out + 46, tmp_c * c1);
    store(out + 38, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(13.4918050467267694), Scalar(-3.11349347232156193));
    store(out + 60, tmp_a * c1);
    store(out + 52, tmp_a * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.656382056840170258);
    store(out + 35, tmp_a * c0);
    store(out + 25, tmp_a * s0);
    tmp_b = z * Scalar(-2.3666191622317525);
    store(out + 47, tmp_b * c0);
    store(out + 37, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-6.7459025233633847), Scalar(0.518915578720260395));
    store(out + 61, tmp_c * c0);
    store(out + 51, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.683184105191914415);
    store(out + 48, tmp_a * c1);
    store(out + 36, tmp_a * s1);
    tmp_b = z * Scalar(2.64596066180190048);
    store(out + 62, tmp_b * c1);
    store(out + 50, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_c = Scalar(-0.707162732524596271);
    store(out + 63, tmp_c * c0);
    store(out + 49, tmp_c * s0);
}

template <typename Array>
void sh_eval_8(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    store(out + 20, fmadd(z * Scalar(1.98431348329844304), load<Value>(out + 12), load<Value>(out + 6) * Scalar(-1.00623058987490532)));
    store(out + 30, fmadd(z * Scalar(1.98997487421323993), load<Value>(out + 20), load<Value>(out + 12) * Scalar(-1.00285307284481395)));
    store(out + 42, fmadd(z * Scalar(1.99304345718356646), load<Value>(out + 30), load<Value>(out + 20) * Scalar(-1.00154202096221923)));
    store(out + 56, fmadd(z * Scalar(1.99489143482413467), load<Value>(out + 42), load<Value>(out + 30) * Scalar(-1.00092721392195827)));
    store(out + 72, fmadd(z * Scalar(1.9960899278339137), load<Value>(out + 56), load<Value>(out + 42) * Scalar(-1.00060078106951478)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-4.6833258049010249), Scalar(2.00713963067186763));
    store(out + 21, tmp_a * c0);
    store(out + 19, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.03100960115899021), tmp_a, tmp_c * Scalar(-0.991031208965114985));
    store(out + 31, tmp_b * c0);
    store(out + 29, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.02131498923702768), tmp_b, tmp_a * Scalar(-0.995226703056238504));
    store(out + 43, tmp_c * c0);
    store(out + 41, tmp_c * s0);
    tmp_a = fmadd(z * Scalar(2.01556443707463773), tmp_c, tmp_b * Scalar(-0.99715504402183186));
    store(out + 57, tmp_a * c0);
    store(out + 55, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.01186954040739119), tmp_a, tmp_c * Scalar(-0.998166817890174474));
    store(out + 73, tmp_b * c0);
    store(out + 71, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(3.31161143515146028), Scalar(-0.473087347878779985));
    store(out + 22, tmp_c * c1);
    store(out + 18, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(7.19030517745998665), Scalar(-2.39676839248666207));
    store(out + 32, tmp_a * c1);
    store(out + 28, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.11394181566096995), tmp_a, tmp_c * Scalar(-0.973610120462326756));
    store(out + 44, tmp_b * c1);
    store(out + 40, tmp_b * s1);
    tmp_c = fmadd(z * Scalar(2.08166599946613307), tmp_b, tmp_a * Scalar(-0.984731927834661791));
    store(out + 58, tmp_c * c1);
    store(out + 54, tmp_c * s1);
    tmp_a = fmadd(z * Scalar(2.06155281280883029), tmp_c, tmp_b * Scalar(-0.990337937660287326));
    store(out + 74, tmp_a * c1);
    store(out + 70, tmp_a * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.590043589926643519);
    store(out + 15, tmp_a * c0);
    store(out + 9, tmp_a * s0);
    tmp_b = z * Scalar(-1.77013076977993067);
    store(out + 23, tmp_b * c0);
    store(out + 17, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-4.40314469491725369), Scalar(0.48923829943525049));
    store(out + 33, tmp_c * c0);
    store(out + 27, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-10.1332578546641603), Scalar(2.76361577854477058));
    store(out + 45, tmp_a * c0);
    store(out + 39, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.20794021658196149), tmp_a, tmp_c * Scalar(-0.95940322360024699));
    store(out + 59, tmp_b * c0);
    store(out + 53, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.15322168769582012), tmp_b, tmp_a * Scalar(-0.975217386560017774));
    store(out + 75, tmp_c * c0);
    store(out + 69, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.625835735449176256);
    store(out + 24, tmp_a * c1);
    store(out + 16, tmp_a * s1);
    tmp_b = z * Scalar(2.07566231488104114);
    store(out + 34, tmp_b * c1);
    store(out + 26, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(5.55021390801596581), Scalar(-0.504564900728724064));
    store(out + 46, tmp_c * c1);
    store(out + 38, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(13.4918050467267694), Scalar(-3.11349347232156193));
    store(out + 60, tmp_a * c1);
    store(out + 52, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.30488611432322132), tmp_a, tmp_c * Scalar(-0.948176387355465389));
    store(out + 76, tmp_b * c1);
    store(out + 68, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.656382056840170258);
    store(out + 35, tmp_a * c0);
    store(out + 25, tmp_a * s0);
    tmp_b = z * Scalar(-2.3666191622317525);
    store(out + 47, tmp_b * c0);
    store(out + 37, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-6.7459025233633847), Scalar(0.518915578720260395));
    store(out + 61, tmp_c * c0);
    store(out + 51, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-17.2495531104905417), Scalar(3.44991062209810817));
    store(out + 77, tmp_a * c0);
    store(out + 67, tmp_a * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.683184105191914415);
    store(out + 48, tmp_a * c1);
    store(out + 36, tmp_a * s1);
    tmp_b = z * Scalar(2.64596066180190048);
    store(out + 62, tmp_b * c1);
    store(out + 50, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(7.98499149089313942), Scalar(-0.532332766059542606));
    store(out + 78, tmp_c * c1);
    store(out + 66, tmp_c * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.707162732524596271);
    store(out + 63, tmp_a * c0);
    store(out + 49, tmp_a * s0);
    tmp_b = z * Scalar(-2.91570664069931995);
    store(out + 79, tmp_b * c0);
    store(out + 65, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_c = Scalar(0.728926660174829988);
    store(out + 80, tmp_c * c1);
    store(out + 64, tmp_c * s1);
}

template <typename Array>
void sh_eval_9(const Array &d, value_t<expr_t<Array>> *out) {
    static_assert(array_size_v<Array> == 3, "The parameter 'd' should be a 3D vector.");

    using Value = value_t<expr_t<Array>>;
    using Scalar = scalar_t<Value>;

    Value x = d.x(), y = d.y(), z = d.z(), z2 = z * z;
    Value c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    store(out + 0, Value(Scalar(0.28209479177387814)));
    store(out + 2, z * Scalar(0.488602511902919923));
    store(out + 6, fmadd(z2, Scalar(0.94617469575756008), Scalar(-0.315391565252520045)));
    store(out + 12, z * fmadd(z2, Scalar(1.865881662950577), Scalar(-1.1195289977703462)));
    store(out + 20, fmadd(z * Scalar(1.98431348329844304), load<Value>(out + 12), load<Value>(out + 6) * Scalar(-1.00623058987490532)));
    store(out + 30, fmadd(z * Scalar(1.98997487421323993), load<Value>(out + 20), load<Value>(out + 12) * Scalar(-1.00285307284481395)));
    store(out + 42, fmadd(z * Scalar(1.99304345718356646), load<Value>(out + 30), load<Value>(out + 20) * Scalar(-1.00154202096221923)));
    store(out + 56, fmadd(z * Scalar(1.99489143482413467), load<Value>(out + 42), load<Value>(out + 30) * Scalar(-1.00092721392195827)));
    store(out + 72, fmadd(z * Scalar(1.9960899278339137), load<Value>(out + 56), load<Value>(out + 42) * Scalar(-1.00060078106951478)));
    store(out + 90, fmadd(z * Scalar(1.99691119506793657), load<Value>(out + 72), load<Value>(out + 56) * Scalar(-1.0004114379931337)));
    c0 = x;
    s0 = y;

    tmp_a = Scalar(-0.488602511902919978);
    store(out + 3, tmp_a * c0);
    store(out + 1, tmp_a * s0);
    tmp_b = z * Scalar(-1.09254843059207896);
    store(out + 7, tmp_b * c0);
    store(out + 5, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-2.28522899732232876), Scalar(0.457045799464465774));
    store(out + 13, tmp_c * c0);
    store(out + 11, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-4.6833258049010249), Scalar(2.00713963067186763));
    store(out + 21, tmp_a * c0);
    store(out + 19, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.03100960115899021), tmp_a, tmp_c * Scalar(-0.991031208965114985));
    store(out + 31, tmp_b * c0);
    store(out + 29, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.02131498923702768), tmp_b, tmp_a * Scalar(-0.995226703056238504));
    store(out + 43, tmp_c * c0);
    store(out + 41, tmp_c * s0);
    tmp_a = fmadd(z * Scalar(2.01556443707463773), tmp_c, tmp_b * Scalar(-0.99715504402183186));
    store(out + 57, tmp_a * c0);
    store(out + 55, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.01186954040739119), tmp_a, tmp_c * Scalar(-0.998166817890174474));
    store(out + 73, tmp_b * c0);
    store(out + 71, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.00935312974101166), tmp_b, tmp_a * Scalar(-0.998749217771908837));
    store(out + 91, tmp_c * c0);
    store(out + 89, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.546274215296039478);
    store(out + 8, tmp_a * c1);
    store(out + 4, tmp_a * s1);
    tmp_b = z * Scalar(1.44530572132027735);
    store(out + 14, tmp_b * c1);
    store(out + 10, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(3.31161143515146028), Scalar(-0.473087347878779985));
    store(out + 22, tmp_c * c1);
    store(out + 18, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(7.19030517745998665), Scalar(-2.39676839248666207));
    store(out + 32, tmp_a * c1);
    store(out + 28, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.11394181566096995), tmp_a, tmp_c * Scalar(-0.973610120462326756));
    store(out + 44, tmp_b * c1);
    store(out + 40, tmp_b * s1);
    tmp_c = fmadd(z * Scalar(2.08166599946613307), tmp_b, tmp_a * Scalar(-0.984731927834661791));
    store(out + 58, tmp_c * c1);
    store(out + 54, tmp_c * s1);
    tmp_a = fmadd(z * Scalar(2.06155281280883029), tmp_c, tmp_b * Scalar(-0.990337937660287326));
    store(out + 74, tmp_a * c1);
    store(out + 70, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.04812235835781919), tmp_a, tmp_c * Scalar(-0.993485272670404207));
    store(out + 92, tmp_b * c1);
    store(out + 88, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.590043589926643519);
    store(out + 15, tmp_a * c0);
    store(out + 9, tmp_a * s0);
    tmp_b = z * Scalar(-1.77013076977993067);
    store(out + 23, tmp_b * c0);
    store(out + 17, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-4.40314469491725369), Scalar(0.48923829943525049));
    store(out + 33, tmp_c * c0);
    store(out + 27, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-10.1332578546641603), Scalar(2.76361577854477058));
    store(out + 45, tmp_a * c0);
    store(out + 39, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.20794021658196149), tmp_a, tmp_c * Scalar(-0.95940322360024699));
    store(out + 59, tmp_b * c0);
    store(out + 53, tmp_b * s0);
    tmp_c = fmadd(z * Scalar(2.15322168769582012), tmp_b, tmp_a * Scalar(-0.975217386560017774));
    store(out + 75, tmp_c * c0);
    store(out + 69, tmp_c * s0);
    tmp_a = fmadd(z * Scalar(2.11804417118980526), tmp_c, tmp_b * Scalar(-0.983662844979209416));
    store(out + 93, tmp_a * c0);
    store(out + 87, tmp_a * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.625835735449176256);
    store(out + 24, tmp_a * c1);
    store(out + 16, tmp_a * s1);
    tmp_b = z * Scalar(2.07566231488104114);
    store(out + 34, tmp_b * c1);
    store(out + 26, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(5.55021390801596581), Scalar(-0.504564900728724064));
    store(out + 46, tmp_c * c1);
    store(out + 38, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(13.4918050467267694), Scalar(-3.11349347232156193));
    store(out + 60, tmp_a * c1);
    store(out + 52, tmp_a * s1);
    tmp_b = fmadd(z * Scalar(2.30488611432322132), tmp_a, tmp_c * Scalar(-0.948176387355465389));
    store(out + 76, tmp_b * c1);
    store(out + 68, tmp_b * s1);
    tmp_c = fmadd(z * Scalar(2.22917715070623501), tmp_b, tmp_a * Scalar(-0.967152839723182112));
    store(out + 94, tmp_c * c1);
    store(out + 86, tmp_c * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.656382056840170258);
    store(out + 35, tmp_a * c0);
    store(out + 25, tmp_a * s0);
    tmp_b = z * Scalar(-2.3666191622317525);
    store(out + 47, tmp_b * c0);
    store(out + 37, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-6.7459025233633847), Scalar(0.518915578720260395));
    store(out + 61, tmp_c * c0);
    store(out + 51, tmp_c * s0);
    tmp_a = z * fmadd(z2, Scalar(-17.2495531104905417), Scalar(3.44991062209810817));
    store(out + 77, tmp_a * c0);
    store(out + 67, tmp_a * s0);
    tmp_b = fmadd(z * Scalar(2.40163634692206163), tmp_a, tmp_c * Scalar(-0.939224604204370817));
    store(out + 95, tmp_b * c0);
    store(out + 85, tmp_b * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.683184105191914415);
    store(out + 48, tmp_a * c1);
    store(out + 36, tmp_a * s1);
    tmp_b = z * Scalar(2.64596066180190048);
    store(out + 62, tmp_b * c1);
    store(out + 50, tmp_b * s1);
    tmp_c = fmadd(z2, Scalar(7.98499149089313942), Scalar(-0.532332766059542606));
    store(out + 78, tmp_c * c1);
    store(out + 66, tmp_c * s1);
    tmp_a = z * fmadd(z2, Scalar(21.3928901909086377), Scalar(-3.77521591604270101));
    store(out + 96, tmp_a * c1);
    store(out + 84, tmp_a * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_a = Scalar(-0.707162732524596271);
    store(out + 63, tmp_a * c0);
    store(out + 49, tmp_a * s0);
    tmp_b = z * Scalar(-2.91570664069931995);
    store(out + 79, tmp_b * c0);
    store(out + 65, tmp_b * s0);
    tmp_c = fmadd(z2, Scalar(-9.26339318284890467), Scalar(0.544905481344053255));
    store(out + 97, tmp_c * c0);
    store(out + 83, tmp_c * s0);
    c1 = fmsub(x, c0, y * s0);
    s1 = fmadd(x, s0, y * c0);

    tmp_a = Scalar(0.728926660174829988);
    store(out + 80, tmp_a * c1);
    store(out + 64, tmp_a * s1);
    tmp_b = z * Scalar(3.17731764895469793);
    store(out + 98, tmp_b * c1);
    store(out + 82, tmp_b * s1);
    c0 = fmsub(x, c1, y * s1);
    s0 = fmadd(x, s1, y * c1);

    tmp_c = Scalar(-0.74890095185318839);
    store(out + 99, tmp_c * c0);
    store(out + 81, tmp_c * s0);
}

NAMESPACE_END(enoki)
