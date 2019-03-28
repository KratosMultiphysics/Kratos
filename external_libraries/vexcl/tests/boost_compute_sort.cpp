#define BOOST_TEST_MODULE BoostComputeSort
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/random.hpp>
#include <vexcl/external/boost_compute.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(boost_compute_sort)
{
    const size_t n = 1024;

    vex::vector<float> X(ctx, n);

    vex::Random<float, vex::random::philox> rnd;
    X = rnd(vex::element_index(), std::rand());

    vex::compute::sort(X);

    std::vector<float> x(n);
    copy(X, x);

    BOOST_CHECK(std::is_sorted(x.begin(), x.end()));
}

BOOST_AUTO_TEST_SUITE_END()
