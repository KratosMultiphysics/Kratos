#include "test.h"
#include <enoki/color.h>

ENOKI_TEST_FLOAT(test01_linear_to_srgb) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return linear_to_srgb(a); },
        [](double value) {
            auto branch1 = 12.92 * value;
            auto branch2 = 1.055 * std::pow(value, 1.0 / 2.4) - 0.055;
            return select(value <= 0.0031308, branch1, branch2);
        },
        Value(0), Value(1),
        60,
        false
    );
}

ENOKI_TEST_FLOAT(test02_srgb_to_linear) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return srgb_to_linear(a); },
        [](double value) {
            auto branch1 = (1.0 / 12.92) * value;
            auto branch2 = std::pow((value + 0.055) * (1.0 / 1.055), 2.4);
            return select(value <= 0.04045, branch1, branch2);
        },
        Value(0), Value(1),
        60,
        false
    );
}

