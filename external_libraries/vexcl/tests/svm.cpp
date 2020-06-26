#define BOOST_TEST_MODULE SVM
#include <boost/test/unit_test.hpp>

#include <vexcl/vector.hpp>

#if defined(CL_VERSION_2_0) || defined(VEXCL_BACKEND_CUDA)
#include <vexcl/svm_vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(svm)
{
    const int n = 1024;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

#if !defined(VEXCL_BACKEND_CUDA)
    if (!vex::Filter::CLVersion(2,0)(ctx.device(0)))
        return;

    vex::push_compile_options(ctx.queue(0), "-cl-std=CL2.0");
#endif

    vex::svm_vector<int> x(queue[0], n);

    {
        auto p = x.map(vex::backend::MAP_WRITE);
        for(int i = 0; i < n; ++i) p[i] = i;
    }

    vex::vector<int> y(queue, n);
    y = x * 2;

    check_sample(y, [&](size_t idx, int v){ BOOST_CHECK_EQUAL(2 * idx, v); });

    x = y / 2;

    {
        auto p = x.map(vex::backend::MAP_READ);
        for(int i = 0; i < n; ++i)
            BOOST_CHECK_EQUAL(p[i], i);
    }

    vex::svm_vector<int> z(queue[0], n);
    z = x;

    {
        auto p = z.map(vex::backend::MAP_READ);
        for(int i = 0; i < n; ++i)
            BOOST_CHECK_EQUAL(p[i], i);
    }

    VEX_FUNCTION(int, dbl, (size_t, i)(int*, x),
            return x[i] * 2;
            );

    y = dbl(vex::element_index(), vex::raw_pointer(x));
    check_sample(y, [&](size_t idx, int v){ BOOST_CHECK_EQUAL(2 * idx, v); });
}

BOOST_AUTO_TEST_SUITE_END()

#else

BOOST_AUTO_TEST_CASE(empty)
{
}

#endif
