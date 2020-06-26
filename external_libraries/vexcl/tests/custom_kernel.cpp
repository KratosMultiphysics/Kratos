#define BOOST_TEST_MODULE CustomKernel
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(custom_kernel)
{
    const cl_ulong n = 1024;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    vex::vector<int> x(queue, n);

    // Single kernel per program
    {
        vex::backend::source_generator src(queue[0]);

        src.begin_kernel("zeros");
        src.begin_kernel_parameters();
        src.parameter<size_t>("n");
        src.parameter<int*>("x");
        src.end_kernel_parameters();
        src.grid_stride_loop("idx", "n").open("{");
        src.new_line() << "x[idx] = 0;";
        src.close("}");
        src.end_kernel();

        vex::backend::kernel zeros(queue[0], src.str(), "zeros");

#ifdef BOOST_NO_VARIADIC_TEMPLATES
        zeros.push_arg(n);
        zeros.push_arg(x(0));
        zeros(queue[0]);
#else
        zeros(queue[0], n, x(0));
#endif

        check_sample(x, [](size_t, int v) { BOOST_CHECK_EQUAL(v, 0); });
    }

    // A couple of kernels per program
    {
        vex::backend::source_generator src(queue[0]);

        src.begin_kernel("ones");
        src.begin_kernel_parameters();
        src.parameter<size_t>("n");
        src.parameter<int*>("x");
        src.end_kernel_parameters();
        src.grid_stride_loop("idx", "n").open("{");
        src.new_line() << "x[idx] = 1;";
        src.close("}");
        src.end_kernel();

        src.begin_kernel("twos");
        src.begin_kernel_parameters();
        src.parameter<size_t>("n");
        src.parameter<int*>("x");
        src.end_kernel_parameters();
        src.grid_stride_loop("idx", "n").open("{");
        src.new_line() << "x[idx] = 2;";
        src.close("}");
        src.end_kernel();

        auto program = vex::backend::build_sources(queue[0], src.str());

        vex::backend::kernel ones(queue[0], program, "ones");
        vex::backend::kernel twos(queue[0], program, "twos");

#ifdef BOOST_NO_VARIADIC_TEMPLATES
        ones.push_arg(n);
        ones.push_arg(x(0));
        ones(queue[0]);
#else
        ones(queue[0], n, x(0));
#endif

        check_sample(x, [](size_t, int v) { BOOST_CHECK_EQUAL(v, 1); });

#ifdef BOOST_NO_VARIADIC_TEMPLATES
        twos.push_arg(n);
        twos.push_arg(x(0));
        twos(queue[0]);
#else
        twos(queue[0], n, x(0));
#endif

        check_sample(x, [](size_t, int v) { BOOST_CHECK_EQUAL(v, 2); });
    }
}

BOOST_AUTO_TEST_SUITE_END()
