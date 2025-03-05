#include "test.h"

ENOKI_TEST_FLOAT(test01_ldexp) {
    const Value inf = std::numeric_limits<Value>::infinity();
    const Value nan = std::numeric_limits<Value>::quiet_NaN();

    for (int i = -10; i < 10; ++i) {
        for (int j = -10; j < 10; ++j) {
            T f = T(std::ldexp(Value(i), j));
            T f2 = ldexp(T(Value(i)), T(Value(j)));
            assert(f == f2);
        }
    }

    assert(T(ldexp(T(inf), T(Value(2)))) == T(inf));
    assert(T(ldexp(T(-inf), T(Value(2)))) == T(-inf));
    assert(all(enoki::isnan(ldexp(T(nan), T(Value(2))))));
}

// AVX512F frexp() uses slightly different conventions
// It is used by log() where this is not a problem though
ENOKI_TEST_FLOAT(test02_frexp) {
    const Value inf = std::numeric_limits<Value>::infinity();
    const Value nan = std::numeric_limits<Value>::quiet_NaN();
    using int_array_t = enoki::int_array_t<T>;
    using Int = typename int_array_t::Value;

    for (int i = -10; i < 10; ++i) {
        if (i == 0)
            continue;
        int e;
        Value f = std::frexp(Value(i), &e);
        T e2, f2;
        std::tie(f2, e2) = frexp(T(Value(i)));
        assert(T(Value(e)) == e2 + 1.f);
        assert(T(f) == f2);
    }

    T e, f;

    std::tie(f, e) = frexp(T(inf));
    assert((std::isinf(f[0]) && !std::isinf(e[0])) ||
           (std::isinf(e[0]) && !std::isinf(f[0])));
    assert(!std::isnan(f[0]) && !std::isnan(e[0]));
    assert(f[0] > 0);

    std::tie(f, e) = frexp(T(-inf));
    assert((std::isinf(f[0]) && !std::isinf(e[0])) ||
           (std::isinf(e[0]) && !std::isinf(f[0])));
    assert(!std::isnan(f[0]) && !std::isnan(e[0]));
    assert(f[0] < 0);

    if (!has_avx512f) {
        std::tie(f, e) = frexp(T(+0.f));
        assert((reinterpret_array<int_array_t>(f) == int_array_t(memcpy_cast<Int>(Value(+0.f)))));

        std::tie(f, e) = frexp(T(-0.f));
        assert((reinterpret_array<int_array_t>(f) == int_array_t(memcpy_cast<Int>(Value(-0.f)))));
    }

    std::tie(f, e) = frexp(T(nan));
    assert(std::isnan(f[0]) || std::isnan(e[0]));
}

ENOKI_TEST_FLOAT(test03_exp) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return exp(a); },
        [](double a) { return std::exp(a); },
        Value(-20), Value(30),
#if defined(ENOKI_X86_AVX512ER)
        27
#else
        3
#endif
    );

    Array<T, 4> x((Value) M_PI);
    Array<T&, 4> y(x);
    assert(exp(x) == exp(y));
}

ENOKI_TEST_FLOAT(test04_log) {
    test::probe_accuracy<T>(
        [](const T &a) -> T { return log(a); },
        [](double a) { return std::log(a); },
        Value(0), Value(2e30),
        2
    );

    Array<T, 4> x((Value) M_PI);
    Array<T&, 4> y(x);
    assert(log(x) == log(y));
}

ENOKI_TEST_FLOAT(test05_pow) {
    assert(T(abs(pow(T(Value(M_PI)), T(Value(-2))) -
               T(Value(0.101321183642338))))[0] < 1e-6f);
}
