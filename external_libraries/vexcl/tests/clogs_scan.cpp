#define BOOST_TEST_MODULE ClogsScan
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/external/clogs.hpp>
#include "context_setup.hpp"

// Test that SFINAE is working correctly for the type introspection
static_assert(vex::clogs::is_scannable<int>::value, "clogs_is_scannable not working");
static_assert(vex::clogs::is_scannable<cl_int4>::value, "clogs_is_scannable not working");
static_assert(!vex::clogs::is_scannable<bool>::value, "clogs_is_scannable not working");
static_assert(!vex::clogs::is_scannable<float>::value, "clogs_is_scannable not working");
static_assert(!vex::clogs::is_scannable<int(int)>::value, "clogs_is_scannable not working");

BOOST_AUTO_TEST_CASE(clogs_scan_scalar)
{
    const size_t n = 1024;

    std::vector<cl_int> x = random_vector<cl_int>(n);

    vex::vector<cl_int> X(ctx, x);
    vex::vector<cl_int> Y(ctx, n);

    vex::clogs::exclusive_scan(X, Y, 1000000000);

    cl_int sum = 1000000000;
    for (size_t i = 0; i < n; i++) {
        cl_int next = sum + x[i];
        x[i] = sum;
        sum = next;
    }

    check_sample(Y, [&](size_t idx, cl_int a) {
        BOOST_CHECK_EQUAL(a, x[idx]);
    });
}

// Test with a vector type to ensure that the type inference works
BOOST_AUTO_TEST_CASE(clogs_scan_vector)
{
    const size_t n = 1234;
    cl_uint4 init = {{1, 2, 3, 4}};

    std::vector<cl_uint4> x = random_vector<cl_uint4>(n);

    vex::vector<cl_uint4> X(ctx, x);
    vex::vector<cl_uint4> Y(ctx, n);

    vex::clogs::exclusive_scan(X, Y, init);

    cl_uint4 sum = init;
    for (size_t i = 0; i < n; i++) {
        cl_uint4 next = sum + x[i];
        x[i] = sum;
        sum = next;
    }

    check_sample(Y, [&](size_t idx, cl_uint4 a) {
        for (int i = 0; i < 4; i++)
            BOOST_CHECK_EQUAL(a.s[i], x[idx].s[i]);
    });
}

BOOST_AUTO_TEST_SUITE_END()
