#define BOOST_TEST_MODULE KernelGenerator
#include <boost/test/unit_test.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/generator.hpp>
#include <vexcl/tagged_terminal.hpp>
#include "context_setup.hpp"

template <class state_type>
state_type sys_func(const state_type &x) {
    return sin(x);
}

template <class state_type, class SysFunction>
void runge_kutta_2(SysFunction sys, state_type &x, double dt) {
    state_type k1 = dt * sys(x);
    state_type k2 = dt * sys(x + 0.5 * k1);

    x += k2;
}

BOOST_AUTO_TEST_CASE(kernel_generator)
{
    typedef vex::symbolic<double> sym_state;

    const size_t n  = 1024;
    const double dt = 0.01;

    std::ostringstream body;
    vex::generator::set_recorder(body);

    sym_state sym_x(sym_state::VectorParameter);

    // Record expression sequence.
    runge_kutta_2(sys_func<sym_state>, sym_x, dt);

    // Build kernel.
    auto kernel = vex::generator::build_kernel(
            ctx, "rk2_stepper", body.str(), sym_x);

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);

    for(int i = 0; i < 100; i++) kernel(X);

    check_sample(X, [&](size_t idx, double a) {
            double s = x[idx];
            for(int i = 0; i < 100; i++)
                runge_kutta_2(sys_func<double>, s, dt);

            BOOST_CHECK_CLOSE(a, s, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(kernel_generator_with_user_function)
{
    typedef vex::symbolic<double> sym_state;

    const size_t n  = 1024;

    std::ostringstream body;
    vex::generator::set_recorder(body);

    sym_state sym_x(sym_state::VectorParameter, sym_state::Const);
    sym_state sym_y(sym_state::VectorParameter);

    VEX_FUNCTION(double, sin2, (double, x),
            double s = sin(x);
            return s * s;
            );

    sym_y = sin2(sym_x);

    auto kernel = vex::generator::build_kernel(
            ctx, "test_sin2", body.str(), sym_x, sym_y);

    vex::vector<double> X(ctx, random_vector<double>(n));
    vex::vector<double> Y(ctx, n);

    for(int i = 0; i < 100; i++) kernel(X, Y);

    check_sample(X, Y, [&](size_t, double x, double y) {
            BOOST_CHECK_CLOSE(y, sin(x) * sin(x), 1e-8);
            });
}
BOOST_AUTO_TEST_CASE(function_generator)
{
    typedef vex::symbolic<double> sym_state;

    const size_t n  = 1024;
    const double dt = 0.01;

    std::ostringstream body;
    vex::generator::set_recorder(body);

    sym_state sym_x(sym_state::VectorParameter);

    // Record expression sequence.
    runge_kutta_2(sys_func<sym_state>, sym_x, dt);

    // Build function.
    // Body string has to be static:
    static std::string function_body = vex::generator::make_function(
            body.str(), sym_x, sym_x);

    VEX_FUNCTION_S(double, rk2, (double, prm1), function_body);

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);

    for(int i = 0; i < 100; i++) {
        X = rk2(X);
    }

    check_sample(X, [&](size_t idx, double a) {
            double s = x[idx];
            for(int i = 0; i < 100; i++)
                runge_kutta_2(sys_func<double>, s, dt);

            BOOST_CHECK_CLOSE(a, s, 1e-8);
            });
}

struct rk2_stepper {
    double dt;

    rk2_stepper(double dt) : dt(dt) {}

    template <class State>
    State operator()(const State &x) const {
        State new_x = x;
        runge_kutta_2(sys_func<State>, new_x, dt);
        return new_x;
    }
};

BOOST_AUTO_TEST_CASE(function_adapter)
{
    const size_t n  = 1024;
    const double dt = 0.01;

    rk2_stepper step(dt);

    auto rk2 = vex::generator::make_function<double(double)>(step);

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);

    for(int i = 0; i < 100; i++) {
        X = rk2(X);
    }

    check_sample(X, [&](size_t idx, double a) {
            double s = x[idx];
            for(int i = 0; i < 100; i++) s = step(s);

            BOOST_CHECK_CLOSE(a, s, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(function_adapter_and_phoenix_lambda)
{
    using namespace boost::phoenix::arg_names;

    const size_t n  = 1024;

    auto squared_radius = vex::generator::make_function<double(double, double)>(
            arg1 * arg1 + arg2 * arg2);

    vex::vector<double> X(ctx, random_vector<double>(n));
    vex::vector<double> Y(ctx, random_vector<double>(n));

    vex::vector<double> Z = squared_radius(X, Y);

    check_sample(X, Y, Z, [&](size_t, double x, double y, double z) {
            BOOST_CHECK_CLOSE(z, x * x + y * y, 1e-8);
            });
}

/*
An alternative variant, which does not use the generator facility.
Intermediate subexpression are captured with help of 'auto' keyword, and
are combined into larger expression.

Note how vex::tag<>() facilitates reuse of kernel parameters.
*/
BOOST_AUTO_TEST_CASE(lazy_evaluation)
{
    const size_t n  = 1024;
    const double dt = 0.01;

    auto rk2 = [](vex::vector<double> &x, double dt) {
        auto X  = vex::tag<1>(x);
        auto DT = vex::tag<2>(dt);

        auto k1 = DT * sin(X);
        auto x1 = X + 0.5 * k1;

        auto k2 = DT * sin(x1);

        x = X + k2;
    };

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);

    for(int i = 0; i < 100; i++) {
        rk2(X, dt);
    }

    check_sample(X, [&](size_t idx, double a) {
            double s = x[idx];
            for(int i = 0; i < 100; i++)
                runge_kutta_2(sys_func<double>, s, dt);

            BOOST_CHECK_CLOSE(a, s, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(element_index)
{
    const size_t n  = 1024;
    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    typedef vex::symbolic<int> sym_vector;

    std::ostringstream body;
    vex::generator::set_recorder(body);

    sym_vector sym_x(sym_vector::VectorParameter);

    sym_x = vex::generator::index();
    auto kernel = vex::generator::build_kernel(queue, "element_index", body.str(), sym_x);

    vex::vector<int> x(queue, n);

    kernel(x);

    check_sample(x, [&](size_t idx, int i) { BOOST_CHECK_EQUAL(i, idx); });
}

BOOST_AUTO_TEST_SUITE_END()
