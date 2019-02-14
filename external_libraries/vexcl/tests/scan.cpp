#define BOOST_TEST_MODULE Scan
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/scan.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(inclusive)
{
    const size_t n = 1000 * 1000;

    std::vector<int> x = random_vector<int>(n);
    vex::vector<int> X(ctx, x);

    vex::inclusive_scan(X, X);

    std::partial_sum(x.begin(), x.end(), x.begin());

    check_sample(X, [&](size_t idx, int v) {
            BOOST_CHECK_EQUAL(v, x[idx]);
            });
}

BOOST_AUTO_TEST_CASE(exclusive)
{
    const size_t n = 1000 * 1000;

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);

    vex::exclusive_scan(X, X);

    std::partial_sum(x.begin(), x.end(), x.begin());
    std::rotate(x.rbegin(), x.rbegin() + 1, x.rend());
    x[0] = 0;

    check_sample(X, [&](size_t idx, double v) {
            BOOST_CHECK_CLOSE(v, x[idx], 1e-8f);
            });
}

BOOST_AUTO_TEST_SUITE_END()
