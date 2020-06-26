#define BOOST_TEST_MODULE VectorCreate
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(empty)
{
    vex::vector<double> x;

    BOOST_CHECK(0U == x.size());
    BOOST_CHECK(x.end() - x.begin() == 0);
}

BOOST_AUTO_TEST_CASE(size)
{
    const size_t N = 1024;
    vex::vector<double> x(ctx, N);

    BOOST_CHECK(x.size() == N);
    BOOST_CHECK(x.end() == x.begin() + N);
}

BOOST_AUTO_TEST_CASE(std_vector)
{
    const size_t N = 1024;

    std::vector<double> x = random_vector<double>(N);
    vex::vector<double> X(ctx, x);

    BOOST_CHECK(X.size() == x.size());

    std::vector<double> y(N);
    copy(X, y);

    check_sample(x, y, [](size_t, double a, double b) { BOOST_CHECK(a == b); });
}

BOOST_AUTO_TEST_CASE(host_pointer)
{
    const size_t N = 1024;

    std::vector<double> x = random_vector<double>(N);
    vex::vector<double> X(ctx, N, x.data());

    BOOST_CHECK(X.size() == x.size());

    std::vector<double> y(N);
    copy(X, y);

    check_sample(x, y, [](size_t, double a, double b) { BOOST_CHECK(a == b); });
}

#ifndef VEXCL_NO_COPY_CONSTRUCTORS
BOOST_AUTO_TEST_CASE(copy_constructor)
{
    const size_t N = 1024;

    vex::vector<double> x1;
    vex::vector<double> x2(x1);

    BOOST_CHECK(x1.size() == 0U);
    BOOST_CHECK(x1.size() == x2.size());

    vex::vector<double> y1(ctx, random_vector<double>(N));
    vex::vector<double> y2(y1);

    BOOST_CHECK(y1.size() == N);
    BOOST_CHECK(y1.size() == y2.size());

    check_sample(y1, y2, [](size_t, double a, double b) { BOOST_CHECK(a == b); });
}
#endif

BOOST_AUTO_TEST_CASE(move_constructor)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);
    x = 42;

    vex::vector<double> y = std::move(x);

    BOOST_CHECK(x.size() == 0U);
    BOOST_CHECK(y.size() == N);

    check_sample(y, [](size_t, double a) { BOOST_CHECK(a == 42); });
}

BOOST_AUTO_TEST_CASE(move_assign)
{
    const size_t N = 1024;
    std::vector<double> x = random_vector<double>(N);
    vex::vector<double> X(ctx, x);

    vex::vector<double> Y;
    Y = std::move(X);

    BOOST_CHECK(Y.size() == x.size());

    check_sample(Y, x, [](size_t, double a, double b) { BOOST_CHECK(a == b); });
}

BOOST_AUTO_TEST_CASE(vector_swap)
{
    const size_t N = 1024;
    const size_t M = 512;

    vex::vector<double> x(ctx, N);
    vex::vector<double> y(ctx, M);

    x = 42;
    y = 67;

    swap(x, y);

    BOOST_CHECK(y.size() == N);
    BOOST_CHECK(x.size() == M);

    check_sample(x, [](size_t, double a) { BOOST_CHECK(a == 67); });
    check_sample(y, [](size_t, double a) { BOOST_CHECK(a == 42); });
}

BOOST_AUTO_TEST_CASE(vector_resize_to_std_vector)
{
    const size_t N = 1024;

    std::vector<double> x = random_vector<double>(N);
    vex::vector<double> X;

    X.resize(ctx, x);
    BOOST_CHECK(X.size() == x.size());

    check_sample(X, x, [](size_t, double a, double b) { BOOST_CHECK(a == b); });
}

BOOST_AUTO_TEST_CASE(vector_resize_to_vex_vector)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);
    x = 42;

    vex::vector<double> y;
    y.resize(x);

    BOOST_CHECK(y.size() == x.size());

    check_sample(x, y, [](size_t, double a, double b) { BOOST_CHECK(a == b); });
}

BOOST_AUTO_TEST_CASE(stl_container_of_vex_vector)
{
    typedef vex::device_vector<unsigned>::raw_type raw_mem;
    const size_t N = 1024;
    const size_t M = 16 + generator<size_t>::get() ;

    std::vector< vex::vector<unsigned> > x;

    std::vector< raw_mem > bufs;

    for(size_t i = 0; i < M; ++i) {
        x.push_back( vex::vector<unsigned>(ctx, random_vector<unsigned>(N)) );
        x.back() = i;
        bufs.push_back( x.back()(0).raw() );
    }

    for(size_t i = 0; i < M; ++i)
        x[i] -= i;

    for(size_t i = 0; i < M; ++i) {
        BOOST_CHECK_EQUAL(N, x[i].size());
        BOOST_CHECK_EQUAL(bufs[i], x[i](0).raw());
    }
}

BOOST_AUTO_TEST_CASE(initialize_with_expression)
{
    const size_t n = 1024;

    std::vector<double> x = random_vector<double>(n);

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y = sin(X);

    check_sample(Y, [&](size_t idx, double v) { BOOST_CHECK_CLOSE(v, sin(x[idx]), 1e-8); });
}

BOOST_AUTO_TEST_CASE(some_devices_are_empty)
{
    vex::vector<double> x(ctx, 1);
    x = 0;
    BOOST_CHECK(x[0] == 0);
}

BOOST_AUTO_TEST_CASE(expression_properties)
{
    const size_t n = 16;
    vex::vector<int> x(ctx, n);

    std::vector<vex::backend::command_queue> q;
    size_t s;

    std::tie(q, s) = vex::expression_properties(2 * x);
    BOOST_CHECK_EQUAL(q.size(), ctx.size());
    BOOST_CHECK_EQUAL(s, n);

    std::tie(q, s) = vex::expression_properties(std::make_tuple(2 * x, x - 1));
    BOOST_CHECK_EQUAL(q.size(), ctx.size());
    BOOST_CHECK_EQUAL(s, n);
}

BOOST_AUTO_TEST_SUITE_END()
