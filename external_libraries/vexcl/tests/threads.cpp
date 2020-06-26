#define BOOST_TEST_MODULE Sort
#include <boost/thread.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/reductor.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(threads)
{
    const size_t n = 1024 * 1024;

    auto run = [n](vex::backend::command_queue queue, cl_long *s) {
        std::vector<vex::backend::command_queue> q(1, queue);
        vex::vector<int> x(q, n);
        x = 1;
        vex::Reductor<cl_long, vex::SUM> sum(q);
        *s = sum(x);
    };

    boost::ptr_vector< boost::thread > threads;
    std::vector< cl_long       > results(ctx.size(), 0);

    for(unsigned d = 0; d < ctx.size(); ++d) {
        threads.push_back( new boost::thread(run, ctx.queue(d), &results[d]) );
    }

    cl_long sum = 0;
    for(unsigned d = 0; d < ctx.size(); ++d) {
        threads[d].join();
        sum += results[d];
    }

    BOOST_CHECK_EQUAL(sum, n * ctx.size());
}

BOOST_AUTO_TEST_SUITE_END()
