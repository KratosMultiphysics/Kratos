#define BOOST_TEST_MODULE ExpressionTypeDeduction
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/mba.hpp>
#ifdef RESOLVED_70
#include <vexcl/spmat.hpp>
#endif
#include <vexcl/tagged_terminal.hpp>
#include <vexcl/temporary.hpp>
#include <vexcl/cast.hpp>
#include <vexcl/vector_view.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

template <typename T>
inline std::array<T, 2> make_array(T x,  T y) {
    std::array<T, 2> p = {{x, y}};
    return p;
}

template <class Result, class Expr>
void check(const Expr &expr) {
    typedef typename vex::detail::return_type<Expr>::type Deduced;

    boost::proto::display_expr( boost::proto::as_child(expr) );
    std::cout << vex::type_name<Deduced>() << std::endl << std::endl;

    BOOST_CHECK( (std::is_same<Deduced, Result>::value) );
}

BOOST_AUTO_TEST_CASE(terminals)
{
    const size_t n = 1024;

    vex::vector<double>     x;
    vex::vector<int>        y;
    vex::vector<cl_double2> z;

    check<int>       (5);
    check<double>    (4.2);
    check<double>    (x);
    check<int>       (y);
    check<cl_double2>(z);
    check<size_t>(vex::element_index());

    {
        vex::mba<2, float> *surf= 0;
        check<float>( (*surf)(x, y) );
    }

#ifdef RESOLVED_70
    {
        vex::SpMatCCSR<double> *A = 0;
        check<double>( (*A) * x );
    }

    {
        std::vector<cl::CommandQueue> q1(1, ctx.queue(0));
        vex::vector<double> x1(q1, n);
        vex::SpMat<double> *A = 0;
        check<double>( vex::make_inline( (*A) * x1 ) );
    }
#endif

    check<double>( vex::tag<1>(x) );

    {
        auto tmp = vex::make_temp<1>(x * y);
        check<double>( tmp );
    }

    {
        std::vector<vex::command_queue> q1(1, ctx.queue(0));
        vex::vector<int> x1(q1, n);
        vex::slicer<1> slice(vex::extents[n]);
        check<int>( slice[vex::range(0, 2, n)](x1) );
    }
}

BOOST_AUTO_TEST_CASE(logical_expr)
{
    vex::vector<double> x;
    vex::vector<int> y;

    check<cl_long>(x < y);
    check<cl_long>(5 > pow(x, 2.0 * y));
    check<cl_long>(!x);
}

BOOST_AUTO_TEST_CASE(nary_expr)
{
    vex::vector<double>     x;
    vex::vector<cl_double2> x2;
    vex::vector<int>        y;
    vex::vector<cl_int2>    y2;

    check<double>    (x + y);
    check<double>    (x + 2 * y);
    check<int>       (-y);
    check<cl_double2>(x * x2);
    check<cl_double2>(y * x2);
    check<cl_double2>(x * y2);
}

BOOST_AUTO_TEST_CASE(user_functions)
{
    vex::vector<double> x;
    vex::vector<int> y;

    VEX_FUNCTION(double, f1, (double, x),            return 42;);
    VEX_FUNCTION(int,    f2, (double, x)(double, y), return 42;);

    check<double>( f1(x) );
    check<int>   ( f2(x, y) );
    check<int>   ( f2(x + y, x - y) );
}

BOOST_AUTO_TEST_CASE(ternary_operator)
{
    vex::vector<double> x;
    vex::vector<int> y;

    check<int>   (  if_else(x < 0,  1,  y) );
    check<double>( *if_else(x < 0, &x, &y) );
}

BOOST_AUTO_TEST_CASE(builtin_functions)
{
    vex::vector<double> x;
    vex::vector<int> y;

    check<double>( cos(x) - sin(y) );
    check<double>( pow(x, 2.0 * y) );
}

BOOST_AUTO_TEST_CASE(reduced_view)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q1(1, ctx.queue(0));
    vex::vector<double> x(q1, n);

    vex::slicer<2> s(vex::extents[32][32]);

    check<double>( vex::reduce<vex::SUM>(s[vex::_](x), 1) );
}

BOOST_AUTO_TEST_CASE(casted_terminals)
{
    check<double>(vex::cast<double>(42));
    check<int>(vex::cast<int>(4.2));
}

BOOST_AUTO_TEST_SUITE_END()
