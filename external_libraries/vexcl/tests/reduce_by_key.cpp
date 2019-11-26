#define BOOST_TEST_MODULE ReduceByKey
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/reduce_by_key.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(rbk)
{
    const int n = 1024 * 1024;

    std::vector<int>    x = random_vector<int>   (n);
    std::vector<double> y = random_vector<double>(n);

    std::sort(x.begin(), x.end());

    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    vex::vector<int>    ikeys(queue, x);
    vex::vector<double> ivals(queue, y);

    vex::vector<int>    okeys;
    vex::vector<double> ovals;

    int num_keys = vex::reduce_by_key(ikeys, ivals, okeys, ovals);

    std::vector<int> ux = x;
    ux.erase( std::unique(ux.begin(), ux.end()), ux.end() );

    BOOST_CHECK_EQUAL(ux.size(),    num_keys);
    BOOST_CHECK_EQUAL(okeys.size(), num_keys);
    BOOST_CHECK_EQUAL(ovals.size(), num_keys);

    check_sample(okeys, ovals, [&](size_t, int key, double dev_sum) {
        double host_sum = std::accumulate(
                y.begin() + (std::lower_bound(x.begin(), x.end(), key) - x.begin()),
                y.begin() + (std::upper_bound(x.begin(), x.end(), key) - x.begin()),
                0.0);
        BOOST_CHECK_CLOSE(dev_sum, host_sum, 1e-8);
        });
}

BOOST_AUTO_TEST_CASE(rbk_tuple)
{
    const int n = 1000 * 1000;

    std::vector<cl_int>  k1(n);
    std::vector<cl_long> k2(n);

    {
        std::vector<cl_int>  k1s = random_vector<cl_int> (n);
        std::vector<cl_long> k2s = random_vector<cl_long>(n);

        std::vector<size_t> idx(n);
        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) {
                return std::make_tuple(k1s[i], k2s[i]) < std::make_tuple(k1s[j], k2s[j]);
                });

        for(int i = 0; i < n; ++i) {
            k1[i] = k1s[idx[i]];
            k2[i] = k2s[idx[i]];
        }
    }

    std::vector<double> y = random_vector<double>(n);

    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    vex::vector<cl_int>  ikey1(queue, k1);
    vex::vector<cl_long> ikey2(queue, k2);
    vex::vector<double>  ivals(queue, y);

    vex::vector<cl_int>  okey1;
    vex::vector<cl_long> okey2;
    vex::vector<double>  ovals;

    VEX_FUNCTION(bool, equal, (cl_int, a1)(cl_long, a2)(cl_int, b1)(cl_long, b2),
            return (a1 == b1) && (a2 == b2);
            );

    VEX_FUNCTION(double, plus, (double, x)(double, y),
            return x + y;
            );

    int num_keys = vex::reduce_by_key(
            boost::fusion::vector_tie(ikey1, ikey2), ivals,
            boost::fusion::vector_tie(okey1, okey2), ovals,
            equal, plus
            );

    std::vector<cl_int>  rk1;
    std::vector<cl_long> rk2;
    std::vector<double>  rsum;

    rk1.reserve(num_keys);
    rk2.reserve(num_keys);
    rsum.reserve(num_keys);

    rk1.push_back(k1[0]);
    rk2.push_back(k2[0]);
    rsum.push_back(y[0]);

    for(int i = 1; i < n; ++i) {
        if (std::make_tuple(k1[i-1], k2[i-1]) == std::make_tuple(k1[i], k2[i])) {
            rsum.back() += y[i];
        } else {
            rk1.push_back(k1[i]);
            rk2.push_back(k2[i]);
            rsum.push_back(y[i]);
        }
    }

    BOOST_CHECK_EQUAL(rsum.size(),  num_keys);
    BOOST_CHECK_EQUAL(okey1.size(), num_keys);
    BOOST_CHECK_EQUAL(okey2.size(), num_keys);
    BOOST_CHECK_EQUAL(ovals.size(), num_keys);

    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    check_sample(okey1, okey2, ovals, [&](size_t i, cl_int key1, cl_long key2, double dsum) {
        BOOST_CHECK_EQUAL(key1, rk1[i]);
        BOOST_CHECK_EQUAL(key2, rk2[i]);
        BOOST_CHECK_CLOSE(dsum, rsum[i], 1e-8);
        });
}

BOOST_AUTO_TEST_SUITE_END()
