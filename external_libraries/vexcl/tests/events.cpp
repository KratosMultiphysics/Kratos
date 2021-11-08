#define BOOST_TEST_MODULE Let
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/enqueue.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(compute_overlap)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q1(1, ctx.queue(0));
    std::vector<vex::command_queue> q2(1, vex::backend::duplicate_queue(ctx.queue(0)));

    vex::vector<int> x(q1, n);
    vex::vector<int> y(q2, n);

    vex::Reductor<size_t> count(q2);

    x = 1;
    q1[0].finish();

    x = 2;

    vex::backend::enqueue_barrier(q2[0], {vex::backend::enqueue_marker(q1[0])});

    y = x;

    BOOST_CHECK_EQUAL(count(y != 2), 0);
}

BOOST_AUTO_TEST_CASE(assignment_queue)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q1(1, ctx.queue(0));
    std::vector<vex::command_queue> q2(1, vex::backend::duplicate_queue(ctx.queue(0)));

    vex::vector<int> x(q1, n);
    vex::vector<int> y(q1, n);

    vex::Reductor<size_t> count(q2);

    x = 1;
    q1[0].finish();

    x = 2;

    vex::backend::enqueue_barrier(q2[0], {vex::backend::enqueue_marker(q1[0])});

    enqueue(q2, y) = x;

    BOOST_CHECK_EQUAL(count(y != 2), 0);
}

BOOST_AUTO_TEST_CASE(transfer_overlap)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q1(1, ctx.queue(0));
    std::vector<vex::command_queue> q2(1, vex::backend::duplicate_queue(ctx.queue(0)));

    vex::vector<int> xd(q1, n);
    std::vector<int> xh(n);

    xd = 1;
    q1[0].finish();

    xd = 2;

    vex::backend::enqueue_barrier(q2[0], {vex::backend::enqueue_marker(q1[0])});

    vex::copy(q2, xd, xh);

    BOOST_CHECK_EQUAL(std::count(xh.begin(), xh.end(), 2), n);
}

BOOST_AUTO_TEST_CASE(enqueue_multiexpression)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q1(1, ctx.queue(0));
    std::vector<vex::command_queue> q2(1, vex::backend::duplicate_queue(ctx.queue(0)));

    vex::multivector<int, 2> x(q1, n);
    vex::multivector<int, 2> y(q1, n);

    vex::Reductor<size_t> count(q2);

    x = std::make_tuple(1, 2);
    q1[0].finish();

    x = std::make_tuple(3, 4);

    vex::backend::enqueue_barrier(q2[0], {vex::backend::enqueue_marker(q1[0])});

    enqueue(q2, y) = x * 2;

    BOOST_CHECK_EQUAL(count(y(0) != 6), 0);
    BOOST_CHECK_EQUAL(count(y(1) != 8), 0);
}

BOOST_AUTO_TEST_SUITE_END()

