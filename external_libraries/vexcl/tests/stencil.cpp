#define BOOST_TEST_MODULE StencilConvolution
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/stencil.hpp>
#include "context_setup.hpp"

struct index {
    size_t n;
    index(size_t n) : n(n) {}

    size_t operator()(size_t i, long shift) const {
        return std::min<size_t>(n - 1, std::max<long>(0, static_cast<long>(i) + shift));
    }
};

BOOST_AUTO_TEST_CASE(stencil_convolution)
{
    const size_t n = 1024;

    std::vector<double> s = random_vector<double>(rand() % 64 + 1);

    int center = rand() % s.size();

    vex::stencil<double> S(ctx, s, center);

    std::vector<double> x = random_vector<double>(n);
    std::generate(x.begin(), x.end(), [](){ return (double)rand() / RAND_MAX; });

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = 1;
    Y += X * S;

    index idx(n);

    check_sample(Y, [&](size_t i, double a) {
        double sum = 1;
        size_t j = 0;
        int k = -center;
        for(; j < s.size(); k++, j++)
            sum += s[j] * x[idx(i, k)];
        BOOST_CHECK_CLOSE(a, sum, 1e-8);
    });

    Y = 42 * (X * S);

    check_sample(Y, [&](size_t i, double a) {
        double sum = 0;
        size_t j = 0;
        int k = -center;
        for(; j < s.size(); k++, j++)
            sum += s[j] * x[idx(i, k)];
        BOOST_CHECK_CLOSE(a, 42 * sum, 1e-8);
    });
}

#if BOOST_VERSION >= 105000
// Boost upto v1.49 segfaults on this test
BOOST_AUTO_TEST_CASE(two_stencils)
{
    const size_t n = 32;
    std::vector<double> s(5, 1);
    vex::stencil<double> S(ctx, s, 3);
    vex::vector<double> X(ctx, n);
    vex::vector<double> Y(ctx, n);

    X = 0;
    Y = X * S + X * S;

    BOOST_CHECK(Y[ 0] == 0);
    BOOST_CHECK(Y[16] == 0);
    BOOST_CHECK(Y[31] == 0);
}
#endif

BOOST_AUTO_TEST_CASE(small_vector)
{
    const size_t n = 128;

    std::vector<double> s = random_vector<double>(rand() % 64 + 1);

    int center = rand() % s.size();

    vex::stencil<double> S(ctx, s, center);

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = 1;
    Y += X * S;

    index idx(n);

    check_sample(Y, [&](size_t i, double a) {
        double sum = 1;
        size_t j = 0;
        int k = -center;
        for(; j < s.size(); k++, j++)
            sum += s[j] * x[idx(i, k)];
        BOOST_CHECK_CLOSE(a, sum, 1e-8);
    });
}

BOOST_AUTO_TEST_CASE(multivector)
{
    typedef std::array<double, 2> elem_t;
    const size_t n = 1024;

    std::vector<double> s = random_vector<double>(rand() % 64 + 1);
    int center = rand() % s.size();

    vex::stencil<double> S(ctx, s.begin(), s.end(), center);

    std::vector<double> x = random_vector<double>(2 * n);

    vex::multivector<double,2> X(ctx, x);
    vex::multivector<double,2> Y(ctx, n);

    Y = 1;
    Y += X * S;

    index idx(n);

    check_sample(Y, [&](size_t i, elem_t a) {
        double sum[2] = {1, 1};
        size_t j = 0;
        int k = -center;
        for(; j < s.size(); k++, j++) {
            sum[0] += s[j] * x[0 + idx(i, k)];
            sum[1] += s[j] * x[n + idx(i, k)];
        }

        BOOST_CHECK_CLOSE(a[0], sum[0], 1e-8);
        BOOST_CHECK_CLOSE(a[1], sum[1], 1e-8);
    });

    Y = 42 * (X * S);

    check_sample(Y, [&](size_t i, elem_t a) {
        double sum[2] = {0, 0};
        size_t j = 0;
        int k = -center;
        for(; j < s.size(); k++, j++) {
            sum[0] += s[j] * x[0 + idx(i, k)];
            sum[1] += s[j] * x[n + idx(i, k)];
        }

        BOOST_CHECK_CLOSE(a[0], 42 * sum[0], 1e-8);
        BOOST_CHECK_CLOSE(a[1], 42 * sum[1], 1e-8);
    });
}

BOOST_AUTO_TEST_CASE(big_stencil)
{
    const size_t n = 1 << 16;

    std::vector<double> s = random_vector<double>(2048);
    int center = rand() % s.size();

    vex::stencil<double> S(ctx, s, center);

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    index idx(n);

    Y = X * S;
    check_sample(Y, [&](size_t i, double a) {
        double sum = 0;
        size_t j = 0;
        int k = -center;
        for(; j < s.size(); k++, j++)
            sum += s[j] * x[idx(i, k)];
        BOOST_CHECK_CLOSE(a, sum, 1e-8);
    });
}

BOOST_AUTO_TEST_CASE(user_defined_stencil)
{
    const size_t n = 1024;

    VEX_STENCIL_OPERATOR(oscillate,
            double, 3, 1,  "return sin(X[1] - X[0]) + sin(X[0] - X[-1]);",
            ctx);

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = oscillate(X);

    index idx(n);

    check_sample(Y, [&](size_t i, double a) {
        size_t left  = idx(i, -1);
        size_t right = idx(i, +1);
        double s = sin(x[right] - x[i]) + sin(x[i] - x[left]);
        BOOST_CHECK_CLOSE(a, s, 1e-8);
    });

#if BOOST_VERSION >= 105000
    Y = 41 * oscillate(X) + oscillate(X);

    check_sample(Y, [&](size_t i, double a) {
        size_t left  = idx(i, -1);
        size_t right = idx(i, +1);
        double s = sin(x[right] - x[i]) + sin(x[i] - x[left]);
        BOOST_CHECK_CLOSE(a, 42 * s, 1e-8);
    });
#endif
}

BOOST_AUTO_TEST_SUITE_END()
