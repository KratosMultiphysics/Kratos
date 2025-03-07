#include "test.h"

ENOKI_TEST_FLOAT(test01_sinh) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return sinh(a); },
        [](double a) { return std::sinh(a); },
        Value(-10), Value(10),
        8
    );

    Array<T, 4> x((Value) 1);
    Array<T&, 4> y(x);
    assert(sinh(x) == sinh(y));
}

ENOKI_TEST_FLOAT(test02_cosh) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return cosh(a); },
        [](double a) { return std::cosh(a); },
        Value(-10), Value(10),
        8
    );

    Array<T, 4> x((Value) 1);
    Array<T&, 4> y(x);
    assert(cosh(x) == cosh(y));
}

ENOKI_TEST_FLOAT(test03_sincosh_sin) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return sincosh(a).first; },
        [](double a) { return std::sinh(a); },
        Value(-10), Value(10),
        8
    );

    Array<T, 4> x((Value) 1), s, c;
    Array<T&, 4> y(x);
    auto result = sincosh(y);
#if !defined(_WIN32)
    assert(result.first == sinh(y) && result.second == cosh(y));
#else
    assert(all_nested(abs(result.first - sinh(y)) < T(1e-6f)) &&
           all_nested(abs(result.second - cosh(y)) < T(1e-6f)));
#endif
}

ENOKI_TEST_FLOAT(test04_sincosh_cos) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return sincosh(a).second; },
        [](double a) { return std::cosh(a); },
        Value(-10), Value(10),
        8
    );
}

ENOKI_TEST_FLOAT(test05_tanh) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return tanh(a); },
        [](double a) { return std::tanh(a); },
        Value(-10), Value(10),
        7
    );

    Array<T, 4> x((Value) 1);
    Array<T&, 4> y(x);
    assert(tanh(x) == tanh(y));
}

ENOKI_TEST_FLOAT(test06_csch) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return csch(a); },
        [](double a) { return 1/std::sinh(a); },
        Value(-10), Value(10),
        8
    );
}

ENOKI_TEST_FLOAT(test07_sech) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return sech(a); },
        [](double a) { return 1/std::cosh(a); },
        Value(-10), Value(10),
        9
    );
}

ENOKI_TEST_FLOAT(test08_coth) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return coth(a); },
        [](double a) { return 1/std::tanh(a); },
        Value(-10), Value(10),
        8
    );
}

ENOKI_TEST_FLOAT(test09_asinh) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return asinh(a); },
        [](double a) { return std::asinh(a); },
        Value(-30), Value(30),
        3
    );
    Array<T, 4> x((Value) 2);
    Array<T&, 4> y(x);
    assert(asinh(x) == asinh(y));
}

ENOKI_TEST_FLOAT(test11_acosh) {
    if (Size == 2 && has_avx512er)
        return; /// Skip for KNL, Clang 7 generates an unsupported SKX+ instruction :(
    test::probe_accuracy<T>(
        [](const T &a) -> T { return acosh(a); },
        [](double a) { return std::acosh(a); },
        Value(1), Value(10),
        5
    );
    Array<T, 4> x((Value) 2);
    Array<T&, 4> y(x);
    assert(acosh(x) == acosh(y));
}

ENOKI_TEST_FLOAT(test12_atanh) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return atanh(a); },
        [](double a) { return std::atanh(a); },
        Value(-1 + 0.001), Value(1 - 0.001),
        3
    );
    Array<T, 4> x((Value) 0.5);
    Array<T&, 4> y(x);
    assert(atanh(x) == atanh(y));
}

