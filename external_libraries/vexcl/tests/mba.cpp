#define BOOST_TEST_MODULE VectorArithmetics
#include <boost/test/unit_test.hpp>
#include <vexcl/mba.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/temporary.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

template <typename T>
inline std::array<T, 2> make_array(T x,  T y) {
    std::array<T, 2> p = {{x, y}};
    return p;
}

BOOST_AUTO_TEST_CASE(mba)
{
    std::vector< std::array<double,2> > p;

    p.push_back(make_array<double>(0.0, 0.0));
    p.push_back(make_array<double>(0.0, 1.0));
    p.push_back(make_array<double>(1.0, 0.0));
    p.push_back(make_array<double>(1.0, 1.0));
    p.push_back(make_array<double>(0.4, 0.4));
    p.push_back(make_array<double>(0.6, 0.6));

    std::vector<double> v;

    v.push_back( 0.2);
    v.push_back( 0.0);
    v.push_back( 0.0);
    v.push_back(-0.2);
    v.push_back(-1.0);
    v.push_back( 1.0);

    vex::mba<2> cloud(ctx,
            make_array<double>(-0.01, -0.01),
            make_array<double>( 1.01,  1.01),
            p, v, make_array<size_t>(2, 2)
            );

    const size_t n = 11;
    vex::vector<double> x(ctx, n);
    vex::vector<double> z(ctx, n);

    x = 1.0 * vex::element_index() / (n - 1.0);
    z = sin(cloud(x, 1.0 * vex::element_index() / (n - 1.0)));

    BOOST_CHECK_CLOSE(static_cast<double>(z[ 0]), sin( 0.2), 1e-6);
    BOOST_CHECK_CLOSE(static_cast<double>(z[ 4]), sin(-1.0), 1e-6);
    BOOST_CHECK_CLOSE(static_cast<double>(z[ 6]), sin( 1.0), 1e-6);
    BOOST_CHECK_CLOSE(static_cast<double>(z[10]), sin(-0.2), 1e-6);

    auto t = vex::make_temp<1>(1.0 * vex::element_index() / (n - 1.0));

    z = cloud(t, t);

    BOOST_CHECK_CLOSE(static_cast<double>(z[ 0]),  0.2, 1e-6);
    BOOST_CHECK_CLOSE(static_cast<double>(z[ 4]), -1.0, 1e-6);
    BOOST_CHECK_CLOSE(static_cast<double>(z[ 6]),  1.0, 1e-6);
    BOOST_CHECK_CLOSE(static_cast<double>(z[10]), -0.2, 1e-6);

}

BOOST_AUTO_TEST_SUITE_END()
