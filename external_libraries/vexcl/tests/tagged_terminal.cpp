#define BOOST_TEST_MODULE TaggedTerminal
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/tagged_terminal.hpp>
#include <vexcl/temporary.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/vector_view.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(tagged_terminal)
{
    using vex::tag;

    const size_t n = 1024;

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y = tag<1>(X) * tag<1>(X);

    check_sample(Y, [&](size_t idx, double v) {
            BOOST_CHECK_CLOSE(v, x[idx] * x[idx], 1e-8); });

    vex::Reductor<double, vex::SUM> sum(ctx);

    BOOST_CHECK_CLOSE(
            sum(tag<3>(2) * tag<1>(X) * tag<1>(X) + tag<3>(2) * tag<2>(Y) * tag<2>(Y)),
            std::accumulate(x.begin(), x.end(), 0.0, [](double sum, double v) -> double {
                double y = v * v;
                return sum + 2 * (y + y * y);
                }),
            1e-6
            );

    tag<1>(X) = 1;
    check_sample(X, [&](size_t, double v) { BOOST_CHECK_CLOSE(v, 1, 1e-8); });

    tag<1>(X) = tag<1>(X) + vex::element_index();
    check_sample(X, [&](size_t idx, double v) { BOOST_CHECK_CLOSE(v, 1.0 + idx, 1e-8); });
}

BOOST_AUTO_TEST_CASE(tagged_slice)
{
    using vex::tag;
    using vex::range;

    const size_t n = 1024;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(queue, x);

    vex::slicer<1> slice(&n);

    vex::vector<double> Y = tag<1>(41) * tag<2>(X) + tag<3>(slice[range()](X));

    check_sample(Y, [&](size_t idx, double v) {
            BOOST_CHECK_CLOSE(v, 42 * x[idx], 1e-8); });
}

BOOST_AUTO_TEST_CASE(temporary_inside_tag)
{
    using vex::tag;

    const size_t n = 1024;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(queue, x);

    auto i = vex::make_temp<1>(vex::element_index() + 1);
    auto reverse = vex::permutation(n - i);
    vex::vector<double> Y(queue, n);

    Y = vex::tag<1>(reverse(X));

    check_sample(Y, [&](size_t idx, double v) {
            BOOST_CHECK_EQUAL(v, x[n - (idx + 1)]); });
}

BOOST_AUTO_TEST_CASE(tied_tagged_terminals)
{
    using vex::tag;

    const size_t n = 1024;

    vex::multivector<double, 2> X(ctx, random_vector<double>(2 * n));
    vex::multivector<double, 2> Y(ctx, n);

#define X0 tag<0>(X(0))
#define X1 tag<1>(X(1))

    Y = std::tie(X1, X0);

    check_sample(Y(1), X(0), [&](size_t, double y, double x) {
            BOOST_CHECK_EQUAL(y, x); });

    check_sample(Y(0), X(1), [&](size_t, double y, double x) {
            BOOST_CHECK_EQUAL(y, x); });
}

BOOST_AUTO_TEST_SUITE_END()
