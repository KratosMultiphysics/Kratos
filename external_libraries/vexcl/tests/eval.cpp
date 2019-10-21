#define BOOST_TEST_MODULE Eval
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/function.hpp>
#include <vexcl/vector_view.hpp>
#include <vexcl/eval.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(eval_atomic) {
    const size_t M = 16;
    const size_t C = 64;
    const size_t N = M * C;

    std::vector<vex::command_queue> q(1, ctx.queue(0));

    vex::vector<int> x(q, N);
    vex::vector<int> y(q, M);

    y = 0;
    x = vex::element_index() % M;

    vex::eval(atomic_add(&vex::permutation(x)(y), 1));
    check_sample(y, [=](size_t, int v) { BOOST_CHECK_EQUAL(v, C); });

    vex::eval(atomic_sub(&vex::permutation(x)(y), 1));
    check_sample(y, [=](size_t, int v) { BOOST_CHECK_EQUAL(v, 0); });
}

BOOST_AUTO_TEST_SUITE_END()
