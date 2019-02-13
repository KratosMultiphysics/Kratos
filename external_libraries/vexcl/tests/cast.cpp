#define BOOST_TEST_MODULE VectorArithmetics
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/cast.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(casted_expession)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);

    x = vex::cast<double>(5);

    check_sample(x, [](size_t, double a) { BOOST_CHECK_EQUAL(a, 5); });
}

#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
BOOST_AUTO_TEST_CASE(convert_functions)
{
    const size_t N = 1024;

    vex::vector<cl_int2> x(ctx, N);

    union {
        cl_float2 f;
        cl_int2   i;
    } y;

    y.f = {{4.2f, 8.4f}};

    x = vex::convert_int2(y.f);

    check_sample(x, [y](size_t, cl_int2 a) {
            BOOST_CHECK_EQUAL(a.s[0], static_cast<int>(y.f.s[0]));
            BOOST_CHECK_EQUAL(a.s[1], static_cast<int>(y.f.s[1]));
            });

    x = vex::as_int2(y.f);

    check_sample(x, [y](size_t, cl_int2 a) {
            BOOST_CHECK_EQUAL(a.s[0], y.i.s[0]);
            BOOST_CHECK_EQUAL(a.s[1], y.i.s[1]);
            });

}
#endif

BOOST_AUTO_TEST_SUITE_END()

