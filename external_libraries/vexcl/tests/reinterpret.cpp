#define BOOST_TEST_MODULE ReinterpretVector
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(reinterpret_vector)
{
    const size_t n = 1024;
    vex::vector<cl_int2> x(ctx, n);

    x.reinterpret<int>() = vex::element_index();

    check_sample(x, [&](size_t i, cl_int2 v) {
            BOOST_CHECK_EQUAL(i*2+0, v.s[0]);
            BOOST_CHECK_EQUAL(i*2+1, v.s[1]);
            });
}

BOOST_AUTO_TEST_SUITE_END()
