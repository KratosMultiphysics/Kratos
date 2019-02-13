#define BOOST_TEST_MODULE MultivectorArithmetics
#include <boost/test/unit_test.hpp>
#include <vexcl/constants.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(arithmetics)
{
    typedef std::array<double, 4> elem_t;

    const size_t n = 1024;

    vex::multivector<double, 4> x(ctx, n);
    vex::multivector<double, 4> y(ctx, random_vector<double>(n * 4));
    vex::multivector<double, 4> z(ctx, random_vector<double>(n * 4));

    vex::Reductor<double,vex::MIN> min(ctx);
    vex::Reductor<double,vex::MAX> max(ctx);

    elem_t v = {{0, 1, 2, 3}};
    x = v;

    BOOST_CHECK(min(x) == v);
    BOOST_CHECK(max(x) == v);

    x = std::make_tuple(1, 2, 3, 4) * y + z;

    check_sample(x, y, z, [](size_t, elem_t a, elem_t b, elem_t c) {
            for(size_t i = 0; i < 4; ++i)
                BOOST_CHECK_CLOSE(a[i], (i + 1) * b[i] + c[i], 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(multivector_multiexpressions)
{
    typedef std::array<double, 2> elem_t;

    const size_t n = 1024;

    vex::multivector<double, 2> x(ctx, n);
    vex::multivector<double, 2> y(ctx, random_vector<double>(n * 2));

    x = std::tie(
            sin( y(0) ) + cos( y(1) ),
            cos( y(0) ) + sin( y(1) )
            );

    check_sample(x, y, [](size_t, elem_t a, elem_t b) {
            BOOST_CHECK_CLOSE(a[0], sin(b[0]) + cos(b[1]), 1e-8);
            BOOST_CHECK_CLOSE(a[1], cos(b[0]) + sin(b[1]), 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(tied_vectors)
{
    const size_t n = 1024;

    vex::vector<double> X(ctx, random_vector<double>(n));
    vex::vector<double> Y(ctx, random_vector<double>(n));

    vex::vector<double> A(ctx, n);
    vex::vector<double> B(ctx, n);

    vex::tie(A, B) = std::tie(X + Y, X - Y);

    check_sample(A, X, Y, [](size_t, double a, double x, double y) {
            BOOST_CHECK(a == x + y);
            });

    check_sample(B, X, Y, [](size_t, double b, double x, double y) {
            BOOST_CHECK(b == x - y);
            });
}

BOOST_AUTO_TEST_CASE(builtin_functions)
{
    typedef std::array<double, 2> elem_t;
    const size_t n = 1024;

    vex::multivector<double, 2> x(ctx, random_vector<double>(n * 2));
    vex::multivector<double, 2> y(ctx, n);

    y = pow(sin(x), 2.0) + pow(cos(x), 2.0);

    check_sample(y, [](size_t, elem_t a) {
            BOOST_CHECK_CLOSE(a[0], 1, 1e-8);
            BOOST_CHECK_CLOSE(a[1], 1, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(user_defined_functions)
{
    typedef std::array<double, 2> elem_t;

    const size_t n = 1024;
    const size_t m = 2;

    vex::multivector<double, m> x(ctx, n);
    vex::multivector<double, m> y(ctx, n);

    elem_t v1 = {{1, 2}};
    elem_t v2 = {{2, 1}};

    VEX_FUNCTION(size_t, greater, (double, x)(double, y), return x > y;);

    x = v1;
    y = v2;

    x = greater(x, y);

    check_sample(x, [](size_t, elem_t a) {
            BOOST_CHECK(a[0] == 0);
            BOOST_CHECK(a[1] == 1);
            });
}

BOOST_AUTO_TEST_CASE(reduction)
{
    typedef std::array<double, 2> elem_t;

    const size_t n = 1024;

    std::vector<double> x = random_vector<double>(n);
    std::vector<double> y = random_vector<double>(n);

    vex::multivector<double, 2> m(ctx, n);
    copy(x, m(0));
    copy(y, m(1));

    vex::Reductor<double, vex::SUM> sum(ctx);
    vex::Reductor<double, vex::MIN> min(ctx);
    vex::Reductor<double, vex::MAX> max(ctx);

    elem_t summ = sum(m);
    elem_t minm = min(m);
    elem_t maxm = max(m);

    BOOST_CHECK_CLOSE(summ[0], std::accumulate(x.begin(), x.end(), 0.0), 1e-6);
    BOOST_CHECK_CLOSE(summ[1], std::accumulate(y.begin(), y.end(), 0.0), 1e-6);

    BOOST_CHECK_CLOSE(minm[0], *std::min_element(x.begin(), x.end()), 1e-12);
    BOOST_CHECK_CLOSE(minm[1], *std::min_element(y.begin(), y.end()), 1e-12);

    BOOST_CHECK_CLOSE(maxm[0], *std::max_element(x.begin(), x.end()), 1e-12);
    BOOST_CHECK_CLOSE(maxm[1], *std::max_element(y.begin(), y.end()), 1e-12);
}

BOOST_AUTO_TEST_CASE(element_index)
{
    typedef std::array<double, 2> elem_t;

    const size_t N = 1024;

    vex::multivector<double, 2> x(ctx, N);

    x = 0.5 * vex::element_index();

    check_sample(x, [](size_t idx, elem_t a) {
            BOOST_CHECK_CLOSE(a[0], 0.5 * idx, 1e-12);
            BOOST_CHECK_CLOSE(a[1], 0.5 * idx, 1e-12);
            });

    x = std::tie(
            sin(0.5 * vex::element_index()),
            cos(0.5 * vex::element_index())
            );

    check_sample(x, [](size_t idx, elem_t a) {
            BOOST_CHECK_CLOSE(a[0], sin(0.5 * idx), 1e-6);
            BOOST_CHECK_CLOSE(a[1], cos(0.5 * idx), 1e-6);
            });
}

BOOST_AUTO_TEST_CASE(compound_assignment)
{
    const size_t n = 1024;
    const size_t m = 2;

    typedef std::array<double, m> elem_t;

    vex::multivector<double, m> x(ctx, n);
    vex::multivector<double, m> y(ctx, random_vector<double>(n * m));

    x = 0;

    x += sin(2 * y);

    check_sample(x, y, [&](size_t, elem_t a, elem_t b) {
            for(size_t i = 0; i < m; ++i)
                BOOST_CHECK_CLOSE(a[i], sin(2 * b[i]), 1e-8);
            });

    x = 0;
    x -= sin(2 * y);

    check_sample(x, y, [&](size_t, elem_t a, elem_t b) {
            for(size_t i = 0; i < m; ++i)
                BOOST_CHECK_CLOSE(a[i], -sin(2 * b[i]), 1e-8);
            });

    x = 1;
    x *= std::tie(y(1), sin(y(0)));

    check_sample(x, y, [](size_t, elem_t a, elem_t b) {
            BOOST_CHECK_CLOSE(a[0], b[1], 1e-8);
            BOOST_CHECK_CLOSE(a[1], sin(b[0]), 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(integral_constants)
{
    typedef std::array<double, 4> elem_t;

    const size_t n = 1024;

    vex::multivector<double, 4> x(ctx, n);

    x = std::integral_constant<int, 42>();
    check_sample(x, [](size_t, elem_t a) {
            for(size_t i = 0; i < 4; ++i) BOOST_CHECK_EQUAL(a[i], 42);
            });

    x = sin( vex::constants::e() * vex::element_index() );
    check_sample(x, [](size_t idx, elem_t a) {
            for(size_t i = 0; i < 4; ++i)
                BOOST_CHECK_CLOSE(a[i], sin(boost::math::constants::e<double>() * idx), 1e-8);
            });
}

#if (VEXCL_CHECK_SIZES > 0)
BOOST_AUTO_TEST_CASE(expression_size_check)
{
    vex::multivector<int, 2> x(ctx, 16);
    vex::multivector<int, 2> y(ctx, 32);

    BOOST_CHECK_THROW(x = y, std::runtime_error);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
