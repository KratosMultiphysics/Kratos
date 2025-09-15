#define BOOST_TEST_MODULE RandomNumbers
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/random.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/tagged_terminal.hpp>
#include <vexcl/temporary.hpp>
#include <vexcl/function.hpp>
#include <boost/math/constants/constants.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(random_numbers)
{
    const size_t N = 1 << 20;

    vex::Reductor<size_t, vex::SUM> sumi(ctx);
    vex::Reductor<double, vex::SUM> sumd(ctx);

    vex::Random<cl_int> rand0;
    vex::vector<cl_uint> x0(ctx, N);
    x0 = rand0(vex::element_index(), std::rand());

    vex::Random<cl_float4> rand1;
    vex::vector<cl_float4> x1(ctx, N);
    x1 = rand1(vex::element_index(), std::rand());

    vex::Random<cl_double4> rand2;
    vex::vector<cl_double4> x2(ctx, N);
    x2 = rand2(vex::element_index(), std::rand());

    vex::Random<cl_double> rand3;
    vex::vector<cl_double> x3(ctx, N);
    x3 = rand3(vex::element_index(), std::rand());

    // X in [0,1]
    BOOST_CHECK(sumi(x3 > 1) == 0);
    BOOST_CHECK(sumi(x3 < 0) == 0);

    // mean = 0.5
    BOOST_CHECK(std::abs((sumd(x3) / N) - 0.5) < 1e-2);

    vex::RandomNormal<cl_double> rand4;
    vex::vector<cl_double> x4(ctx, N);
    x4 = rand4(vex::element_index(), std::rand());

    // E(X ~ N(0,s)) = 0
    BOOST_CHECK(std::abs(sumd(x4)/N) < 1e-2);

    // E(abs(X) ~ N(0,s)) = sqrt(2/M_PI) * s
    BOOST_CHECK(std::abs(sumd(fabs(x4))/N - std::sqrt(2 / boost::math::constants::pi<double>())) < 1e-2);

    vex::Random<cl_double, vex::random::threefry> rand5;
    vex::vector<cl_double> x5(ctx, N);
    x5 = rand5(vex::element_index(), std::rand());

    BOOST_CHECK(std::abs(sumd(x5)/N - 0.5) < 1e-2);
}

BOOST_AUTO_TEST_CASE(monte_carlo_pi)
{
    vex::Random<double, vex::random::threefry> rnd;

    vex::Reductor<size_t, vex::SUM> sum(ctx);

    const size_t n = 1 << 20;

    auto i = vex::tag<0>(vex::element_index(0, n));

    auto x = vex::make_temp<1>(rnd(i, std::rand()));
    auto y = vex::make_temp<2>(rnd(i, std::rand()));

    double pi = 4.0 * sum( (x * x + y * y) < 1 ) / n;

    BOOST_CHECK_CLOSE(pi, boost::math::constants::pi<double>(), 0.5);
}

BOOST_AUTO_TEST_SUITE_END()

