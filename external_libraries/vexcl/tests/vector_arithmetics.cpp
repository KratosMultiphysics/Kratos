#define BOOST_TEST_MODULE VectorArithmetics
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/test/unit_test.hpp>
#include <vexcl/constants.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/tagged_terminal.hpp>
#include <vexcl/function.hpp>
#ifndef VEXCL_BACKEND_CUDA
#include <vexcl/vector_view.hpp>
#include <vexcl/constant_address_space.hpp>
#endif
#include <boost/math/constants/constants.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(access_element)
{
    vex::vector<double> x(ctx, 5);
    x = 5;
    x[2] = 2;
    x.at(3) = 3;


    BOOST_CHECK_EQUAL(x[1], 5.0);
    BOOST_CHECK_EQUAL(x.at(2), 2.0);
    BOOST_CHECK_EQUAL(x[3], 3.0);
    BOOST_CHECK_THROW(x.at(5), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(assign_expression)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);
    vex::vector<double> y(ctx, N);
    vex::vector<double> z(ctx, N);

    y = 42;
    z = 67;
    x = 5 * sin(y) + z;

    check_sample(x, [](size_t, double a) {
            BOOST_CHECK_CLOSE(a, 5 * sin(42.0) + 67, 1e-12);
            });
}

BOOST_AUTO_TEST_CASE(compound_assignment)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, n);

    x = 0;
    x += 1;

    check_sample(x, [](size_t, double a) { BOOST_CHECK(a == 1); });

    x -= 2;

    check_sample(x, [](size_t, double a) { BOOST_CHECK(a == -1); });
}

BOOST_AUTO_TEST_CASE(reduce_expression)
{
    namespace acc = boost::accumulators;

    const size_t N = 1024;

    std::vector<double> x = random_vector<double>(N);
    std::transform(x.begin(), x.end(), x.begin(), [](double v){ return (v - 0.5) * 1e8; });

    vex::vector<double> X(ctx, x);

    vex::Reductor<double,vex::SUM      > sum(ctx);
    vex::Reductor<double,vex::MIN      > min(ctx);
    vex::Reductor<double,vex::MAX      > max(ctx);
    vex::Reductor<double,vex::SUM_Kahan> csum(ctx);

    acc::accumulator_set< double, acc::stats< acc::tag::sum_kahan > > stat;
    std::for_each(x.begin(), x.end(), std::ref(stat));

    BOOST_CHECK_CLOSE(sum(X),  acc::sum_kahan(stat), 1e-8);
    BOOST_CHECK_CLOSE(csum(X), acc::sum_kahan(stat), 1e-8);

    BOOST_CHECK_EQUAL(min(X), *std::min_element(x.begin(), x.end()));
    BOOST_CHECK_EQUAL(max(X), *std::max_element(x.begin(), x.end()));

#ifndef BOOST_NO_VARIADIC_TEMPLATES
    vex::Reductor<double,vex::MIN_MAX  > minmax(ctx);
    auto mm = minmax(X);
    BOOST_CHECK_EQUAL(mm.s[0], *std::min_element(x.begin(), x.end()));
    BOOST_CHECK_EQUAL(mm.s[1], *std::max_element(x.begin(), x.end()));
#endif

    BOOST_CHECK_EQUAL( max( fabs(X - X) ), 0.0);
}

BOOST_AUTO_TEST_CASE(builtin_functions)
{
    const size_t N = 1024;
    std::vector<double> x = random_vector<double>(N);
    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, N);

    Y = pow(sin(X), 2.0) + pow(cos(X), 2.0);

    check_sample(Y, [](size_t, double a) { BOOST_CHECK_CLOSE(a, 1, 1e-8); });
}

BOOST_AUTO_TEST_CASE(user_defined_functions)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);
    vex::vector<double> y(ctx, N);

    x = 1;
    y = 2;

    VEX_FUNCTION(size_t, greater, (double, x)(double, y), return x > y;);

    vex::Reductor<size_t,vex::SUM> sum(ctx);

    BOOST_CHECK( sum( greater(x, y) ) == 0U);
    BOOST_CHECK( sum( greater(y, x) ) == N);
}

BOOST_AUTO_TEST_CASE(user_defined_functions_same_signature)
{
    const size_t N = 1024;
    vex::vector<double> x(ctx, N);

    x = 1;

    VEX_FUNCTION(double, times2, (double, x), return x * 2;);
    VEX_FUNCTION(double, times4, (double, x), return x * 4;);

    vex::Reductor<size_t,vex::SUM> sum(ctx);

    BOOST_CHECK( sum( times2(x) ) == 2 * N );
    BOOST_CHECK( sum( times4(x) ) == 4 * N );
}

BOOST_AUTO_TEST_CASE(element_index)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);

    x = sin(0.5 * vex::element_index());

    check_sample(x, [](size_t idx, double a) { BOOST_CHECK_CLOSE(a, sin(0.5 * idx), 1e-6); });
}

#if !defined(VEXCL_BACKEND_CUDA) && !defined(VEXCL_BACKEND_JIT) && !defined(__APPLE__)
BOOST_AUTO_TEST_CASE(vector_values)
{
    const size_t N = 1024;

    VEX_FUNCTION(cl_int4, make_int4, (int, x), return (int4)(x, x, x, x););

    cl_int4 c = {{1, 2, 3, 4}};

    vex::vector<cl_int4> X(ctx, N);
    X = c * (make_int4(5 + vex::element_index()));

    check_sample(X, [c](size_t idx, cl_int4 v) {
            for(int j = 0; j < 4; ++j)
                BOOST_CHECK(v.s[j] == c.s[j] * (5 + static_cast<int>(idx)));
            });
}
#endif

BOOST_AUTO_TEST_CASE(nested_functions)
{
    const size_t N = 1024;

    VEX_FUNCTION(int, f, (int, x), return 2 * x;);
    VEX_FUNCTION(int, g, (int, x), return 3 * x;);

    vex::vector<int> x(ctx, N);

    x = 1;
    x = f(f(x));
    check_sample(x, [](size_t, int a) { BOOST_CHECK(a == 4); });

    x = 1;
    x = g(f(x));
    check_sample(x, [](size_t, int a) { BOOST_CHECK(a == 6); });
}

BOOST_AUTO_TEST_CASE(custom_header)
{
    const size_t n = 1024;

    vex::vector<int> x(ctx, n);

    vex::push_program_header(ctx, "#define THE_ANSWER 42\n");

    VEX_FUNCTION(int, answer, (int, x), return x * THE_ANSWER;);

    x = answer(1);

    check_sample(x, [](size_t, int a) {
            BOOST_CHECK(a == 42);
            });

    vex::pop_program_header(ctx);
}

// Have to define these outside of the following test case scope.
// In Visual C++ is not types defined in enclosing function scope
// are not able to refeence each other.
VEX_FUNCTION(double, sin2, (double, x), return pow(sin(x), 2.0););
VEX_FUNCTION(double, cos2, (double, x), return pow(cos(x), 2.0););

BOOST_AUTO_TEST_CASE(function_with_preamble)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, random_vector<double>(n));
    vex::vector<double> y(ctx, n);

    VEX_FUNCTION_D(double, one, (double, x), (sin2)(cos2),
            return sin2(x) + cos2(x);
            );

    y = one(x);

    check_sample(y, [](size_t, double a) {
            BOOST_CHECK_CLOSE(a, 1.0, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(ternary_operator)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, random_vector<double>(n));
    vex::vector<double> y(ctx, n);

    y = if_else(x > 0.5, sin(x), cos(x));

    check_sample(x, y, [&](size_t, double X, double Y) {
            BOOST_CHECK_CLOSE(Y, X > 0.5 ? sin(X) : cos(X), 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(assign_to_ternary_operator)
{
    const size_t n = 32;

    vex::vector<double> x(ctx, random_vector<double>(n));
    vex::vector<double> y(ctx, n);
    vex::vector<double> z(ctx, n);

    y = 0;
    z = 0;

    vex::tie( *if_else(x < 0.5, &y, &z) ) = 42;

    check_sample(x, y, z, [&](size_t, double X, double Y, double Z) {
            BOOST_CHECK_EQUAL(X < 0.5 ? Y : Z, 42);
            BOOST_CHECK_EQUAL(X < 0.5 ? Z : Y, 0);
            });
}

BOOST_AUTO_TEST_CASE(combine_expressions)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, n);

    auto alpha  = vex::constants::two_pi() * vex::tag<1>(vex::element_index());
    auto sine   = sin(alpha);
    auto cosine = cos(alpha);

    x = pow(sine, 2.0) + pow(cosine, 2.0);

    check_sample(x, [](size_t, double v) { BOOST_CHECK_CLOSE(v, 1.0, 1e-8); });
}

BOOST_AUTO_TEST_CASE(constants)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, n);

    x = std::integral_constant<int, 42>();
    check_sample(x, [](size_t, double v) { BOOST_CHECK_EQUAL(v, 42); });

    x = vex::constants::pi();
    check_sample(x, [](size_t, double v) { BOOST_CHECK_CLOSE(v, boost::math::constants::pi<double>(), 1e-8); });
}

#ifndef VEXCL_BACKEND_CUDA
BOOST_AUTO_TEST_CASE(constant_address_space)
{
    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    const size_t N = 1024;

    int host_data[] = {1, 2, 3, 4};

    vex::vector<int> x(queue, 4, host_data, vex::backend::MEM_READ_ONLY);
    vex::vector<int> y(queue, N);

    y = vex::permutation(vex::element_index() % 4)( vex::constant(x) );

    check_sample(y, [&](size_t idx, int v) {
            BOOST_CHECK_EQUAL(v, host_data[idx % 4]);
            });
}
#endif

#if (VEXCL_CHECK_SIZES > 0)
BOOST_AUTO_TEST_CASE(expression_size_check)
{
    vex::vector<int> x(ctx, 16);
    vex::vector<int> y(ctx, 32);

    BOOST_CHECK_THROW(x = y, std::runtime_error);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
