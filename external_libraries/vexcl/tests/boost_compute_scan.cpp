#define BOOST_TEST_MODULE BoostComputeScan
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/external/boost_compute.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(boost_compute_inclusive_scan)
{
    const size_t n = 1024;

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    vex::compute::inclusive_scan(X, Y);

    std::partial_sum(x.begin(), x.end(), x.begin());

    check_sample(Y, [&](size_t idx, double a) {
            BOOST_CHECK_CLOSE(a, x[idx], 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(boost_compute_exclusive_scan)
{
    const size_t n = 1024;

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    vex::compute::exclusive_scan(X, Y);

    double sum = 0;
    for(auto xi = x.begin(); xi != x.end(); ++xi) {
        double next = sum + *xi;
        *xi = sum;
        sum = next;
    }

    check_sample(Y, [&](size_t idx, double a) {
            BOOST_CHECK_CLOSE(a, x[idx], 1e-8);
            });
}

BOOST_AUTO_TEST_SUITE_END()
