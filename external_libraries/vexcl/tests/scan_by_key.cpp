#define BOOST_TEST_MODULE ReduceByKey
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/scan_by_key.hpp>
#include "context_setup.hpp"

bool nvidia_cl(const vex::Context &ctx) {
#if defined(VEXCL_BACKEND_CUDA) || defined(VEXCL_BACKEND_JIT)
    return false;
#else
    return vex::Filter::Platform("NVIDIA")(ctx.device(0));
#endif
}

BOOST_AUTO_TEST_CASE(sbk)
{
    // NVIDIA OpenCL compiler crashes on scan_by_key kernels.
    if (nvidia_cl(ctx)) return;

    const int n = 1000;

    std::vector<int> x = random_vector<int>(n);
    std::vector<int> y = random_vector<int>(n);

    std::sort(x.begin(), x.end());

    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    vex::vector<int> ikeys(queue, x);
    vex::vector<int> ivals(queue, y);
    vex::vector<int> ovals(queue, n);

    vex::inclusive_scan_by_key(ikeys, ivals, ovals);

    check_sample(ovals, [&](size_t i, int v) {
            if (i == 0)
                BOOST_CHECK_EQUAL(v, y[i]);
            else if (x[i-1] == x[i])
                BOOST_CHECK_EQUAL(
                    y[i],
                    static_cast<int>(ovals[i]) - static_cast<int>(ovals[i-1])
                    );
            else
                BOOST_CHECK_EQUAL(v, y[i]);
            });

    vex::exclusive_scan_by_key(ikeys, ivals, ovals);

    check_sample(ovals, [&](size_t i, int v) {
            if (i == 0)
                BOOST_CHECK_EQUAL(v, 0);
            else if (x[i-1] == x[i])
                BOOST_CHECK_EQUAL(
                    y[i-1],
                    static_cast<int>(ovals[i]) - static_cast<int>(ovals[i-1])
                    );
            else
                BOOST_CHECK_EQUAL(v, 0);
            });
}

BOOST_AUTO_TEST_CASE(sbk_tuple)
{
    // NVIDIA OpenCL compiler crashes on scan_by_key kernels.
    if (nvidia_cl(ctx)) return;

    const int n = 1000;

    std::vector<int> x1 = random_vector<int>(n);
    std::vector<int> x2 = random_vector<int>(n);
    std::vector<int> y  = random_vector<int>(n);

    std::sort(x1.begin(), x1.end());
    std::sort(x2.begin(), x2.end());

    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    vex::vector<int> ikey1(queue, x1);
    vex::vector<int> ikey2(queue, x2);
    vex::vector<int> ivals(queue, y);
    vex::vector<int> ovals(queue, n);

    VEX_FUNCTION(bool, equal, (int, a1)(int, a2)(int, b1)(int, b2),
            return a1 == b1 && a2 == b2;
            );

    VEX_FUNCTION(int, plus, (int, x)(int, y),
            return x + y;
            );

    vex::inclusive_scan_by_key(
            boost::fusion::vector_tie(ikey1, ikey2), ivals, ovals,
            equal, plus);

    check_sample(ovals, [&](size_t i, int v) {
            if (i == 0)
                BOOST_CHECK_EQUAL(v, y[i]);
            else if (std::make_tuple(x1[i-1], x2[i-1]) == std::make_tuple(x1[i], x2[i]))
                BOOST_CHECK_EQUAL(
                    y[i],
                    static_cast<int>(ovals[i]) - static_cast<int>(ovals[i-1])
                    );
            else
                BOOST_CHECK_EQUAL(v, y[i]);
            });

    vex::exclusive_scan_by_key(
            boost::fusion::vector_tie(ikey1, ikey2), ivals, ovals,
            equal, plus);

    check_sample(ovals, [&](size_t i, int v) {
            if (i == 0)
                BOOST_CHECK_EQUAL(v, 0);
            else if (std::make_tuple(x1[i-1], x2[i-1]) == std::make_tuple(x1[i], x2[i]))
                BOOST_CHECK_EQUAL(
                    y[i-1],
                    static_cast<int>(ovals[i]) - static_cast<int>(ovals[i-1])
                    );
            else
                BOOST_CHECK_EQUAL(v, 0);
            });
}

BOOST_AUTO_TEST_SUITE_END()
