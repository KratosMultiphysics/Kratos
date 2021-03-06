Converting generic C++ algorithms to OpenCL/CUDA
================================================

CUDA and OpenCL differ in their handling of compute kernels compilation. In
NVIDIA's framework the compute kernels are compiled to PTX code together with
the host program. In OpenCL the compute kernels are compiled at runtime from
high-level C-like sources, adding an overhead which is particularly noticeable
for smaller sized problems. This distinction leads to higher initialization
cost of OpenCL programs, but at the same time it allows one to generate better
optimized kernels for the problem at hand. VexCL exploits this possibility with
help of its kernel generator mechanism. Moreover, VexCL's CUDA backend uses the
same technique to generate and compile CUDA kernels at runtime.

An instance of :cpp:class:`vex::symbolic\<T>` dumps to an output stream any
arithmetic operations it is being subjected to. For example, this code snippet:

.. code-block:: cpp

    vex::generator::set_recorder(std::cout);
    vex::symbolic<double> x = 6, y = 7;
    x = sin(x * y);

results in the following output::

    double var1 = 6;
    double var2 = 7;
    var1 = sin( ( var1 * var2 ) );

Kernel generator
----------------

The symbolic type allows one to record a sequence of arithmetic operations made
by a generic C++ algorithm. To illustrate the idea, consider the generic
implementation of a 4th order Runge-Kutta ODE stepper:

.. code-block:: cpp

    template <class state_type, class SysFunction>
    void runge_kutta_4(SysFunction sys, state_type &x, double dt) {
        state_type k1 = dt * sys(x);
        state_type k2 = dt * sys(x + 0.5 * k1);
        state_type k3 = dt * sys(x + 0.5 * k2);
        state_type k4 = dt * sys(x + k3);

        x += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

This function takes a system function ``sys``, state variable ``x``, and advances
``x`` by the time step ``dt``. For example, to model the equation
:math:`\mbox{d}x/\mbox{d}t = \sin(x)`,
one has to provide the following system function:

.. code-block:: cpp

    template <class state_type>
    state_type sys_func(const state_type &x) {
        return sin(x);
    }

The following code snippet makes one hundred RK4 iterations for a single
``double`` value on a CPU:

.. code-block:: cpp

    double x = 1, dt = 0.01;

    for(int step = 0; step < 100; ++step)
        runge_kutta_4(sys_func<double>, x, dt);

Let's now generate the kernel for a single RK4 step and apply the kernel to a
:cpp:class:`vex::vector\<double>` (by doing this we essentially simultaneously
solve a large number of identical ODEs with different initial conditions).

.. code-block:: cpp

    // Set recorder for expression sequence.
    std::ostringstream body;
    vex::generator::set_recorder(body);

    // Create symbolic variable.
    typedef vex::symbolic<double> sym_state;
    sym_state sym_x(sym_state::VectorParameter);

    // Record expression sequience for a single RK4 step.
    double dt = 0.01;
    runge_kutta_4(sys_func<sym_state>, sym_x, dt);

    // Build kernel from the recorded sequence.
    auto kernel = vex::generator::build_kernel(ctx, "rk4_stepper", body.str(), sym_x);

    // Create initial state.
    const size_t n = 1024 * 1024;
    vex::vector<double> x(ctx, n);
    x = 10.0 * vex::element_index() / n;

    // Make 100 RK4 steps.
    for(int i = 0; i < 100; i++) kernel(x);

This approach has some obvious restrictions. Namely, the C++ code has to be
embarrassingly parallel and is not allowed to contain any branching or
data-dependent loops. Nevertheless, the kernel generation facility may save
a substantial amount of both human and machine time when applicable.

.. doxygenclass:: vex::symbolic
    :members:

.. doxygenfunction:: vex::generator::set_recorder
.. doxygenfunction:: vex::generator::build_kernel


Function generator
------------------

VexCL also provides a user-defined function generator which takes a function
signature and generic function object, and returns custom VexCL function ready
to be used in vector expressions. Let's rewrite the above example using an
autogenerated function for a Runge-Kutta stepper. First, we need to implement
generic functor:

.. code-block:: cpp

    struct rk4_stepper {
        double dt;

        rk4_stepper(double dt) : dt(dt) {}

        template <class state_type>
        state_type operator()(const state_type &x) const {
            state_type new_x = x;
            runge_kutta_4(sys_func<state_type>, new_x, dt);
            return new_x;
        }
    };

Now we can generate and apply the custom function:

.. code-block:: cpp

    double dt = 0.01;
    rk4_stepper stepper(dt);

    // Generate custom VexCL function:
    auto rk4 = vex::generator::make_function<double(double)>(stepper);

    // Create initial state.
    const size_t n = 1024 * 1024;
    vex::vector<double> x(ctx, n);
    x = 10.0 * vex::element_index() / n;

    // Use the function to advance initial state:
    for(int i = 0; i < 100; i++) x = rk4(x);

Note that both ``runge_kutta_4()`` and ``rk4_stepper`` may be reused for
the host-side computations.

It is very easy to generate a VexCL function from a Boost.Phoenix_ lambda
expression (since Boost.Phoenix_ lambdas are themselves generic functors):

.. code-block:: cpp

    using namespace boost::phoenix::arg_names;
    using vex::generator::make_function;

    auto squared_radius = make_function<double(double, double)>(arg1 * arg1 + arg2 * arg2);

    Z = squared_radius(X, Y);

.. _Boost.Phoenix: http://www.boost.org/doc/libs/release/libs/phoenix/doc/html/index.html

.. doxygenfunction:: vex::generator::make_function(Functor&&)
