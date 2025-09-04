#define BOOST_TEST_MODULE VectorCopy
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/gather.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(iterate_over_vector)
{
    const size_t N = 1024;

    vex::vector<double> x(ctx, N);
    x = 42;

    BOOST_CHECK(42 == *std::min_element(x.begin(), x.end()));
}

BOOST_AUTO_TEST_CASE(element_access)
{
    const size_t N = 1024;
    vex::vector<double> x(ctx, N);

    for(size_t i = 0; i < N; i++)
        x[i] = 42;

    check_sample(x, [](size_t, double a) { BOOST_CHECK(a == 42); });
}

BOOST_AUTO_TEST_CASE(copy_to_std_vector)
{
    const size_t N = 1024;
    vex::vector<double> X(ctx, N);
    std::vector<double> x(N);

    X = 42;
    copy(X, x);
    check_sample(x, [](size_t, double a) { BOOST_CHECK(a == 42); });

    X = 67;
    vex::copy(X.begin(), X.end(), x.begin());
    check_sample(x, [](size_t, double a) { BOOST_CHECK(a == 67); });
}

BOOST_AUTO_TEST_CASE(copy_from_std_vector)
{
    const size_t N = 1024;

    std::vector<double> x = random_vector<double>(N);
    vex::vector<double> X(ctx, N);

    copy(x, X);
    check_sample(X, x, [](size_t, double a, double b) { BOOST_CHECK(a == b); });

    std::fill(x.begin(), x.end(), 42);
    vex::copy(x.begin(), x.end(), X.begin());
    check_sample(X, [](size_t, double a) { BOOST_CHECK(a == 42); });
}

BOOST_AUTO_TEST_CASE(map_buffer)
{
    const size_t N = 1 << 20;
    vex::vector<size_t> x(ctx, N);

    for(unsigned d = 0; d < ctx.size(); ++d) {
        auto ptr = x.map(d);
        for(size_t i = 0; i < x.part_size(d); ++i)
            ptr[i] = i + x.part_start(d);
    }

    check_sample(x, [](size_t idx, size_t a) { BOOST_CHECK(a == idx); });
}

BOOST_AUTO_TEST_CASE(gather)
{
    const size_t n = 1 << 20;
    const size_t m = 100;

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);

    std::vector<size_t> i(m);
    std::generate(i.begin(), i.end(), [n](){ return rand() % n; });

    for (int sorted = 0; sorted < 2; ++sorted) {
        if (sorted) {
            std::sort(i.begin(), i.end());
            i.resize( std::unique(i.begin(), i.end()) - i.begin() );
        }

        std::vector<double>  data(i.size());
        vex::gather  get(ctx, x.size(), i);
        vex::scatter put(ctx, x.size(), i);

        get(X, data);

        for(size_t p = 0; p < i.size(); ++p)
            BOOST_CHECK(data[p] == x[i[p]]);

        vex::vector<double> Y(ctx, n);
        Y = 0;
        put(data, Y);

        for(size_t p = 0; p < i.size(); ++p)
            BOOST_CHECK(Y[i[p]] == x[i[p]]);
    }
}

BOOST_AUTO_TEST_CASE(std_sort_vex_vector)
{
    const size_t n = 1 << 10;

    vex::vector<double> x(ctx, random_vector<double>(n));

    std::sort(x.begin(), x.end());
    BOOST_CHECK(std::is_sorted(x.begin(), x.end()));
}

BOOST_AUTO_TEST_SUITE_END()

