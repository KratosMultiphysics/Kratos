/*
    enoki/color.h -- Color space transformations (only sRGB so far)

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array.h>

NAMESPACE_BEGIN(enoki)

template <typename T> expr_t<T> linear_to_srgb(const T &x) {
    using Value  = expr_t<T>;
    using Mask   = mask_t<Value>;
    using Scalar = scalar_t<Value>;
    constexpr bool Single = std::is_same_v<Scalar, float>;

    Value r = Scalar(12.92);
    Mask large_mask = x > Scalar(0.0031308);

    if (ENOKI_LIKELY(any(large_mask))) {
        Value y = sqrt(x), p, q;

        if constexpr (Single) {
            p = poly5(y, -0.0016829072605308378, 0.03453868659826638,
                      0.7642611304733891, 2.0041169284241644,
                      0.7551545191665577, -0.016202083165206348);
            q = poly5(y, 4.178892964897981e-7, -0.00004375359692957097,
                      0.03467195408529984, 0.6085338522168684,
                      1.8970238036421054, 1.);
        } else {
            p = poly10(y, -3.7113872202050023e-6, -0.00021805827098915798,
                       0.002531335520959116, 0.2263810267005674,
                       3.0477578489880823, 15.374469584296442,
                       32.44669922192121, 27.901125077137042, 8.450947414259522,
                       0.5838023820686707, -0.0031151377052754843);
            q = poly10(y, 2.2380622409188757e-11, -8.387527630781522e-9,
                       0.00007045228641004039, 0.007244514696840552,
                       0.21749170309546628, 2.575446652731678,
                       13.297981743005433, 30.50364355650628, 29.70548706952188,
                       10.723011300050162, 1.);
        }

        masked(r, large_mask) = p / q;
    }

    return r * x;
}

template <typename T> expr_t<T> srgb_to_linear(const T &x) {
    using Value  = expr_t<T>;
    using Mask   = mask_t<Value>;
    using Scalar = scalar_t<Value>;
    constexpr bool Single = std::is_same_v<Scalar, float>;

    Value r = Scalar(1.0 / 12.92);
    Mask large_mask = x > Scalar(0.04045);

    if (ENOKI_LIKELY(any(large_mask))) {
        Value p, q;

        if constexpr (Single) {
            p = poly4(x, -0.0163933279112946, -0.7386328024653209,
                      -11.199318357635072, -47.46726633009393,
                      -36.04572663838034);
            q = poly4(x, -0.004261480793199332, -19.140923959601675,
                      -59.096406619244426, -18.225745396846637, 1.);
        } else {
            p = poly9(x, -0.008042950896814532, -0.5489744177844188,
                      -14.786385491859248, -200.19589605282445,
                      -1446.951694673217, -5548.704065887224,
                      -10782.158977031822, -9735.250875334352,
                      -3483.4445569178347, -342.62884098034357);
            q = poly9(x, -2.2132610916769585e-8, -9.646075249097724,
                      -237.47722999429413, -2013.8039726540235,
                      -7349.477378676199, -11916.470977597566,
                      -8059.219012060384, -1884.7738197074218,
                      -84.8098437770271, 1.);
        }

        masked(r, large_mask) = p / q;
    }

    return r * x;
}

NAMESPACE_END(enoki)
