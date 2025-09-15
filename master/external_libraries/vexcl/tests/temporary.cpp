#define BOOST_TEST_MODULE TemporaryTerminal
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/temporary.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(temporary)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, random_vector<double>(n));
    vex::vector<double> y(ctx, n);

    VEX_FUNCTION(double, sqr, (double, x), return x * x;);

    {
        // Deduce temporary type
        auto s = vex::make_temp<1>( sqr(x) + 25 );
        y = s * (x + s);

        check_sample(y, [&](size_t idx, double v) {
                double X = x[idx];
                double S = X * X + 25;
                BOOST_CHECK_CLOSE(v, S * (X + S), 1e-8);
                });
    }

    {
        // Provide temporary type
        auto s = vex::make_temp<1, double>( sqr(x) + 25 );
        y = s * (x + s);

        check_sample(y, [&](size_t idx, double v) {
                double X = x[idx];
                double S = X * X + 25;
                BOOST_CHECK_CLOSE(v, S * (X + S), 1e-8);
                });
    }
}

BOOST_AUTO_TEST_CASE(nested_temporary)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, random_vector<double>(n));
    vex::vector<double> y(ctx, n);

    auto t1 = vex::make_temp<1>( log(x) );
    auto t2 = vex::make_temp<2>( t1 + sin(x) );

    y = t1 * t2;

    check_sample(y, [&](size_t idx, double v) {
            double X = x[idx];
            double T1 = log(X);
            double T2 = T1 + sin(X);
            BOOST_CHECK_CLOSE(v, T1 * T2, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(reduce_temporary)
{
    const size_t n = 1024;

    vex::vector<double> x(ctx, random_vector<double>(n));

    auto t1 = vex::make_temp<1>( pow(sin(x), 2) );
    auto t2 = vex::make_temp<2>( pow(cos(x), 2) );

    vex::Reductor<double, vex::SUM> sum(ctx);

    BOOST_CHECK_CLOSE(sum(10 * (t1 + t2)), 10.0 * n, 1e-6);
}

BOOST_AUTO_TEST_CASE(multiexpression_temporary)
{
    typedef std::array<double, 2> elem_t;
    const size_t n = 1024;

    vex::vector<double> x(ctx, random_vector<double>(n));

    vex::multivector<double, 2> y(ctx,n);

    auto tmp = vex::make_temp<1>( sin(x) );

    y = std::tie(tmp, sqrt(1 - tmp * tmp));

    check_sample(y, [&](size_t idx, elem_t v){
            double X = x[idx];
            BOOST_CHECK_CLOSE(v[0], sin(X), 1e-8);
            BOOST_CHECK_CLOSE(v[1], cos(X), 1e-8);
            });
}

#if (BOOST_VERSION >= 105200) && !defined(_MSC_VER)
BOOST_AUTO_TEST_CASE(multivector_temporary)
{
    typedef std::array<double, 2> elem_t;
    const size_t n = 1024;

    vex::multivector<double, 2> x(ctx, random_vector<double>(2 * n));
    vex::multivector<double, 2> y(ctx,n);

    auto tmp = vex::make_temp<1, double>( tan(x) );

    y = tmp * tmp;

    check_sample(y, [&](size_t idx, elem_t v){
            elem_t X = x[idx];

            BOOST_CHECK_CLOSE(v[0], pow(tan(X[0]), 2.0), 1e-8);
            BOOST_CHECK_CLOSE(v[1], pow(tan(X[1]), 2.0), 1e-8);
            });
}
#endif

BOOST_AUTO_TEST_SUITE_END()

