#define BOOST_TEST_MODULE Constants
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/constants.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(constants) {
    const size_t N = 1024;

    vex::vector<int> x(ctx, N);

    VEX_CONSTANT(forty_two, 42);

    x = forty_two();

    BOOST_CHECK_EQUAL(forty_two, 42);

    check_sample(x, [=](size_t, int v) { BOOST_CHECK_EQUAL(v, forty_two); });
}

BOOST_AUTO_TEST_SUITE_END()

