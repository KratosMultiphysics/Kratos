/*
    test/test.h -- Rudimentary test runner framework

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2018 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#if defined(NDEBUG)
#  undef NDEBUG
#endif

namespace test {
    /* Global parameters and statistics */
    static bool detailed = false;
    static int nonvectorized_count = 0;
    static int alloc_count = 0;
    static int realloc_count = 0;
    static int dealloc_count = 0;
}

//#define ENOKI_TRACK_SCALAR(reason) { printf("[%s] ", reason); ++test::nonvectorized_count; }
#define ENOKI_TRACK_SCALAR(reason) { ++test::nonvectorized_count; }
#define ENOKI_TRACK_ALLOC(ptr, size) { ++test::alloc_count; }
#define ENOKI_TRACK_DEALLOC(ptr, size) { ++test::dealloc_count; }
#define ENOKI_TEST(name) void name(); static test::Test name##_test{#name, &name}; void name()

#include <enoki/array.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <random>
#include <sstream>

/// Generic string conversion routine
template <typename T> inline std::string to_string(const T& value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

using namespace enoki;

NAMESPACE_BEGIN(test)

bool replace(std::string &str, const std::string &from, const std::string &to) {
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

class Test {
public:
    using TestStorage = std::vector<std::pair<std::string, void (*)()>>;

    Test(const char *name, void (*test)()) {
        if (!registered_tests)
            registered_tests = new TestStorage();
        registered_tests->push_back(std::make_pair(name, test));
    }

    static void run_all() {
        if (!registered_tests)
            return;

        for (auto &item : *registered_tests) {
            std::string &name = item.first;
            replace(name, "uint32_t", "uint32");
            replace(name, "int32_t", "int32");
            replace(name, "uint64_t", "uint64");
            replace(name, "int64_t", "int64");
            replace(name, "float", "float32");
            replace(name, "half", "float16");
            replace(name, "double", "float64");
        }

        std::sort(
            registered_tests->begin(), registered_tests->end(),
            [](const auto &t1, const auto &t2) { return t1.first < t2.first; });

        std::string last_prefix;
        size_t nonvectorized_tests = 0, total_tests = 0;
        for (auto& test : *registered_tests) {
            std::string name = test.first, prefix = "";
            if (name.substr(0, 6) == "array_") {
                name = name.substr(6);
                auto offset = name.find("_");
                if (offset == std::string::npos)
                    throw std::out_of_range("Invalid test name");
                std::string type = name.substr(0, offset);
                name = name.substr(offset+1);
                offset = name.find("_");
                if (offset == std::string::npos)
                    throw std::out_of_range("Invalid test name");
                std::string length = name.substr(0, offset);
                while (!length.empty() && length[0] == '0')
                    length = length.substr(1);
                name = name.substr(offset+1);
                prefix = "Array<" + type + ", " + length + ">";
            }

            if (prefix != last_prefix) {
                if (!last_prefix.empty()) {
                    std::cout << std::endl
                              << "Done. (" << total_tests - nonvectorized_tests << "/"
                              << total_tests << " tests were vectorized)" << std::endl
                              << std::endl;
                    nonvectorized_tests = total_tests = 0;
                }

                if (prefix.empty())
                    std::cout << "Remaining tests .." << std::endl;
                else
                    std::cout << "Testing \"" << prefix << "\" .." << std::endl;
            }
            last_prefix = prefix;

            std::cout << "    " << name << ".. ";
            std::cout.flush();
            int nvi = test::nonvectorized_count;
            test.second();
            std::cout << "ok.";
            if (test::nonvectorized_count != nvi) {
                std::cout << " [non-vectorized]";
                nonvectorized_tests++;
            }
            std::cout << std::endl;
            total_tests++;
        }

        if (total_tests > 0) {
            std::cout << std::endl
                      << "Done. (" << total_tests - nonvectorized_tests << "/"
                      << total_tests << " tests were vectorized)" << std::endl
                      << std::endl;
        }
    }

private:
    static TestStorage *registered_tests;
};

Test::TestStorage *Test::registered_tests = nullptr;

/// Compute the difference in ULPs between two values
template <typename T> T ulpdiff(T ref, T val) {
    const T eps = std::numeric_limits<T>::epsilon() / 2;

    /* Express mantissa wrt. same exponent */
    int e_ref, e_val;
    T m_ref = std::frexp(ref, &e_ref);
    T m_val = std::frexp(val, &e_val);

    T diff;
    if (e_ref == e_val)
        diff = m_ref - m_val;
    else
        diff = m_ref - std::ldexp(m_val, e_val - e_ref);

    return std::abs(diff) / eps;
}

template <typename Value>
void assert_close(Value result, Value result_ref, Value eps) {
    using UInt = std::conditional_t<sizeof(Value) == 4, uint32_t, uint64_t>;
    using Float = std::conditional_t<sizeof(Value) == 4, float, double>;
    bool fail = false;

    if (eps == 0 || !std::isfinite(Float(result_ref))) {
        fail = memcpy_cast<UInt>(result) != memcpy_cast<UInt>(result_ref)
            && !(std::isnan(Float(result_ref)) && std::isnan(Float(result)));
    } else {
        fail = !(std::abs(Float(result - result_ref)) /
                  std::max(Float(1e-6), std::abs(Float(result_ref))) < Float(eps));
    }

    if (fail) {
        std::cerr << std::endl << "    * Failure: got "
            << result << " (0x" << std::hex << memcpy_cast<UInt>(result) << std::dec << "), expected "
            << result_ref << " (0x" << std::hex << memcpy_cast<UInt>(result_ref) << ")!" << std::endl;
        abort();
    }
}

std::string print_float(double value) {
    double exp = std::log10(std::abs(value));
    exp += std::copysign(1.f, exp);
    double mant = value * std::pow(10.0, -double(int(exp)));
    std::ostringstream oss;
    oss.precision(2);
    oss << mant << " \\cdot 10^{" << int(exp) << "}";
    return oss.str();
}


template <typename T>
void probe_accuracy(T (*func)(const T &), double (*ref)(double),
                    const value_t<T> &min, const value_t<T> &max,
                    const value_t<T> &max_ulp_assert = 0,
                    bool check_special = true) {
    using Value = value_t<T>;

    Value max_err = 0, avg_err = 0;
    Value max_err_rel = 0, avg_err_rel = 0;
    Value avg_ulp = 0, max_ulp = 0;
    Value max_ulp_pos = -1;
    Value max_err_pos = -1;

    auto update_errors = [&](Value x, Value v, Value v_ref) {
        Value err = std::abs(v - v_ref);
        Value err_rel = err / std::abs(v_ref);
        if (v_ref == 0 && v == 0)
            err_rel = 0;

        Value ulp = ulpdiff(v_ref, v);
        if (ulp > max_ulp)
            max_ulp_pos = x;
        if (err > max_err)
            max_err_pos = x;

        avg_err += err;
        max_err = std::max(max_err, err);
        avg_err_rel += err_rel;
        max_err_rel = std::max(max_err_rel, err_rel);
        avg_ulp += ulp;
        max_ulp = std::max(max_ulp, ulp);
    };

    T value, result_ref;

    size_t idx = 0;
    size_t nTries = test::detailed ? 1000000 : 10000;
    for (size_t i = 0; i < nTries; i++) {
        value[idx] = min + (Value(i) / Value(nTries - 1)) * (max - min);
        result_ref[idx] = Value(ref(double(value[idx])));
        idx++;

        if (idx == value.size()) {
            auto&& result = func(value);
            for (size_t k = 0; k < idx; ++k)
                update_errors(value[k], result[k], result_ref[k]);
            idx = 0;
        }
    }

    if (idx > 0) {
        auto&& result = func(value);
        for (size_t k = 0; k < idx; ++k)
            update_errors(value[k], result[k], result_ref[k]);
    }

    if (check_special) {
        /// Test function behavior for +/- inf and NaN-valued inputs
        const Value inf = std::numeric_limits<Value>::infinity();
        const Value nan = std::numeric_limits<Value>::quiet_NaN();
        assert_close(func(T(inf))[0], (Value) ref(double(inf)), Value(1e-6));
        assert_close(func(T(-inf))[0], (Value) ref(double(-inf)), Value(1e-6));
        assert_close(func(T(nan))[0], (Value) ref(double(nan)), Value(1e-6));
    }

    bool success = max_ulp <= max_ulp_assert || max_ulp_assert == 0;

    if (test::detailed || !success) {
#if 1
        std::cout << "(in [" << min << ", " << max << "]):" << std::endl
              << "     * avg abs. err = " << avg_err / float(nTries) << std::endl
              << "     * avg rel. err = " << avg_err_rel / float(nTries) << std::endl
              << "        -> in ULPs  = " << avg_ulp / float(nTries) << std::endl
              << "     * max abs. err = " << max_err << std::endl
              << "       (at x=" << max_err_pos << ")" << std::endl
              << "     * max rel. err = " << max_err_rel << std::endl
              << "       -> in ULPs   = " << max_ulp << std::endl
              << "       (at x=" << max_ulp_pos << ")" << std::endl;
#else
        std::cout << "* - :math:`name(x)`" << std::endl;
        std::cout << "  - :math:`" << min << " < x < " << max << "`" << std::endl;
        // abs.error
        std::cout << "  - :math:`" << print_float(avg_err / float(nTries)) << "`"<< std::endl;
        std::cout << "  - :math:`" << print_float(max_err) << "`"<< std::endl;
        std::cout.precision(2);
        std::cout << "  - :math:`" << print_float(avg_err_rel / float(nTries))
                        << "\\,(" << avg_ulp / float(nTries) << "\\,\\mathrm{ulp})"<< "`"<< std::endl;
        std::cout << "  - :math:`" << print_float(max_err_rel) << "\\,(" << max_ulp<< "\\,\\mathrm{ulp})"<< "`"<< std::endl;

#if 0

    * - ``sin(x)``
      - :math:`|x|<8192`
      - :math:`6.61\cdot 10^{-9}`
      - :math:`1.37\cdot 10^{-8}\,(0.16)`
      - :math:`5.96\cdot 10^{-8}`
      - :math:`1.176\cdot 10^{-6}\,(19)`
#endif
#endif
    }
    assert(success);
}

template <typename T>
void validate_unary(const std::vector<value_t<T>> &args,
                    T (*func)(const T &),
                    value_t<T> (*ref)(value_t<T>),
                    const value_t<T> &eps = 0) {
    T value, result_ref;
    size_t idx = 0;

    for (size_t i = 0; i < args.size(); ++i) {
        value_t<T> arg_i = args[i];
        value[idx] = arg_i;
        result_ref[idx] = ref(arg_i);
        idx++;

        if (idx == value.size()) {
            auto&& result = func(value);
            for (size_t k = 0; k < idx; ++k)
                assert_close(result[k], result_ref[k], eps);
            idx = 0;
        }
    }

    if (idx > 0) {
        auto&& result = func(value);
        for (size_t k = 0; k < idx; ++k)
            assert_close(result[k], result_ref[k], eps);
    }
}


template <typename T>
void validate_binary(const std::vector<value_t<T>> &args,
                     T (*func)(const T &, const T &),
                     value_t<T> (*ref)(value_t<T>, value_t<T>),
                     value_t<T> eps = 0) {
    using Value = value_t<T>;
    T value1, value2, result_ref;
    size_t idx = 0;

    for (size_t i=0; i<args.size(); ++i) {
        for (size_t j=0; j<args.size(); ++j) {
            Value arg_i = args[i], arg_j = args[j];
            value1[idx] = arg_i; value2[idx] = arg_j;
            result_ref[idx] = ref(arg_i, arg_j);
            idx++;

            if (idx == value1.size()) {
                auto&& result = func(value1, value2);
                for (size_t k = 0; k < idx; ++k)
                    assert_close(result[k], result_ref[k], eps);
                idx = 0;
            }
        }
    }

    if (idx > 0) {
        auto&& result = func(value1, value2);
        for (size_t k = 0; k < idx; ++k)
            assert_close(result[k], result_ref[k], eps);
    }
}

/// Generate an array containing values that are used to test various operations
template <typename Value> std::vector<Value> sample_values(bool has_nan = true) {
    std::vector<Value> args;
    if (std::is_integral<Value>::value) {
        args = { Value(0), Value(1),    Value(5),
                 Value(7), Value(1234), Value(1000000) };
        if (std::is_signed<Value>::value) {
            args.push_back(Value(-1));
            args.push_back(Value(-5));
            args.push_back(Value(-7));
            args.push_back(Value(-1234));
            args.push_back(Value(-1000000));
        }
    } else if (std::is_floating_point<Value>::value) {
        args = { Value(0), Value(0.5), Value(0.6), Value(1), Value(2),
                 Value(3), Value(M_PI), Value(-0), Value(-0.5), Value(-0.6),
                 Value(-1), Value(-2), Value(-3),
                 Value(std::numeric_limits<float>::infinity()),
                 Value(-std::numeric_limits<float>::infinity())
        };
        if (has_nan)
            args.push_back(Value(std::numeric_limits<float>::quiet_NaN()));
    }
    return args;
}

template <typename T>
void validate_ternary(const std::vector<value_t<T>> &args,
                      T (*func)(const T &, const T &, const T &),
                      value_t<T> (*refFunc)(value_t<T>,
                                                    value_t<T>,
                                                    value_t<T>),
                      const value_t<T> &eps = 0) {
    T value1, value2, value3, result_ref;
    size_t idx = 0;

    for (size_t i = 0; i < args.size(); ++i) {
        for (size_t j = 0; j < args.size(); ++j) {
            for (size_t k = 0; k < args.size(); ++k) {
                value_t<T> arg_i = args[i], arg_j = args[j],
                                   arg_k = args[k];
                value1[idx] = arg_i;
                value2[idx] = arg_j;
                value3[idx] = arg_k;
                result_ref[idx] = refFunc(arg_i, arg_j, arg_k);
                idx++;

                if (idx == value1.size()) {
                    auto &&result = func(value1, value2, value3);
                    for (size_t k = 0; k < idx; ++k)
                        assert_close(result[k], result_ref[k], eps);
                    idx = 0;
                }
            }
        }
    }

    if (idx > 0) {
        auto &&result = func(value1, value2, value3);
        for (size_t k = 0; k < idx; ++k)
            assert_close(result[k], result_ref[k], eps);
    }
}

template <typename T>
void validate_horizontal(const std::vector<value_t<T>> &args,
                         value_t<T> (*func)(const T &),
                         value_t<T> (*refFunc)(const T &),
                         const value_t<T> &eps = 0) {
    std::mt19937 gen;
    std::uniform_int_distribution<> dis(0, (int) args.size()-1);
    T value;

    for (int i=0; i<1000; ++i) {
        for (size_t i=0; i<value.size(); ++i)
            value[i] = args[(size_t) dis(gen)];

        assert_close(func(value), refFunc(value), eps);
    }
}

NAMESPACE_END(test)

#if defined(_MSC_VER)
    /* Don't build the large 31x/32x tests to reduce compilation time on AppVeyor */
#define ENOKI_TEST_HELPER(name, type)                                           \
    ENOKI_TEST(array_##type##_01##_##name) { name<type, 1>();  }                \
    ENOKI_TEST(array_##type##_02##_##name) { name<type, 2>();  }                \
    ENOKI_TEST(array_##type##_03##_##name) { name<type, 3>();  }                \
    ENOKI_TEST(array_##type##_04##_##name) { name<type, 4>();  }                \
    ENOKI_TEST(array_##type##_08##_##name) { name<type, 8>();  }                \
    ENOKI_TEST(array_##type##_16##_##name) { name<type, 16>(); }
#elif defined(ENOKI_ARM_32) || defined(ENOKI_ARM_64)
#define ENOKI_TEST_HELPER(name, type)                                           \
    ENOKI_TEST(array_##type##_01##_##name) { name<type, 1>();  }                \
    ENOKI_TEST(array_##type##_02##_##name) { name<type, 2>();  }                \
    ENOKI_TEST(array_##type##_03##_##name) { name<type, 3>();  }                \
    ENOKI_TEST(array_##type##_04##_##name) { name<type, 4>();  }                \
    ENOKI_TEST(array_##type##_08##_##name) { name<type, 8>();  }
#else
#define ENOKI_TEST_HELPER(name, type)                                           \
    ENOKI_TEST(array_##type##_01##_##name) { name<type, 1>();  }                \
    ENOKI_TEST(array_##type##_02##_##name) { name<type, 2>();  }                \
    ENOKI_TEST(array_##type##_03##_##name) { name<type, 3>();  }                \
    ENOKI_TEST(array_##type##_04##_##name) { name<type, 4>();  }                \
    ENOKI_TEST(array_##type##_08##_##name) { name<type, 8>();  }                \
    ENOKI_TEST(array_##type##_16##_##name) { name<type, 16>(); }                \
    ENOKI_TEST(array_##type##_31##_##name) { name<type, 31>(); }                \
    ENOKI_TEST(array_##type##_32##_##name) { name<type, 32>(); }
#endif

#define ENOKI_TEST_TYPE(name, type)                                             \
    template <typename Value, size_t Size,                                      \
              typename T = enoki::Array<Value, Size>>                           \
    void name##_##type();                                                       \
    ENOKI_TEST_HELPER(name##_##type, type)                                      \
    template <typename Value, size_t Size, typename T>                          \
    void name##_##type()

#define ENOKI_TEST_FLOAT(name)                                                  \
    template <typename Value, size_t Size,                                      \
              typename T = enoki::Array<Value, Size>>                           \
    void name();                                                                \
    ENOKI_TEST_HELPER(name, float)                                              \
    ENOKI_TEST_HELPER(name, double)                                             \
    template <typename Value, size_t Size, typename T>                          \
    void name()

#define ENOKI_TEST_INT(name)                                                    \
    template <typename Value, size_t Size,                                      \
              typename T = enoki::Array<Value, Size>>                           \
    void name();                                                                \
    ENOKI_TEST_HELPER(name, int32_t)                                            \
    ENOKI_TEST_HELPER(name, uint32_t)                                           \
    ENOKI_TEST_HELPER(name, int64_t)                                            \
    ENOKI_TEST_HELPER(name, uint64_t)                                           \
    template <typename Value, size_t Size, typename T>                          \
    void name()

#define ENOKI_TEST_ALL(name)                                                    \
    template <typename Value, size_t Size,                                      \
              typename T = enoki::Array<Value, Size>>                           \
    void name();                                                                \
    ENOKI_TEST_HELPER(name, float)                                              \
    ENOKI_TEST_HELPER(name, double)                                             \
    ENOKI_TEST_HELPER(name, int32_t)                                            \
    ENOKI_TEST_HELPER(name, uint32_t)                                           \
    ENOKI_TEST_HELPER(name, int64_t)                                            \
    ENOKI_TEST_HELPER(name, uint64_t)                                           \
    template <typename Value, size_t Size, typename T>                          \
    void name()

int main(int argc, char** argv) {
    std::cout << "=== Enoki test suite (version " << ENOKI_VERSION << ") ===" << std::endl;
    std::cout << "Enabled compiler features: ";

    if (has_avx512f)         std::cout << "avx512f ";
    if (has_avx512cd)        std::cout << "avx512cd ";
    if (has_avx512dq)        std::cout << "avx512dq ";
    if (has_avx512vl)        std::cout << "avx512vl ";
    if (has_avx512bw)        std::cout << "avx512bw ";
    if (has_avx512pf)        std::cout << "avx512pf ";
    if (has_avx512er)        std::cout << "avx512er ";
    if (has_avx512vpopcntdq) std::cout << "avx512vpopcntdq ";
    if (has_avx2)            std::cout << "avx2 ";
    if (has_avx)             std::cout << "avx ";
    if (has_fma)             std::cout << "fma ";
    if (has_f16c)            std::cout << "f16c ";
    if (has_sse42)           std::cout << "sse4.2 ";
    if (has_x86_64)          std::cout << "x86_64 ";
    if (has_x86_32)          std::cout << "x86 ";
    if (has_neon)            std::cout << "neon ";
    if (has_arm_32)          std::cout << "arm ";
    if (has_arm_64)          std::cout << "aarch64 ";

    /* Turn on verbose mode if requested */
    if (argc == 2 && strcmp(argv[1], "-v") == 0)
        test::detailed = true;

    (void) test::alloc_count;
    (void) test::realloc_count;
    (void) test::dealloc_count;

    std::cout << std::endl;
    std::cout << std::endl;
    test::Test::run_all();
    return 0;
}
