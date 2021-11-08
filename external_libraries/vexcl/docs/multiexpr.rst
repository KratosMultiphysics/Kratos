Multivectors and multiexpressions
=================================

The :cpp:class:`vex::multivector\<T,N>` class allows to store several equally
sized device vectors and perform computations on each component in sync.  Each
operation is delegated to the underlying vectors, but usually results in the
launch of a single fused kernel. Expressions may include values of
``std::array<T,N>`` where ``N`` is equal to the number of multivector
components, or appropriately sized tuples. Each component gets the
corresponding element of either the array or the tuple when the expression is
applied. Similarly, :cpp:func:`vex::multivector::operator[]` or reduction of a
multivector returns an instance of ``std::array<T,N>``.
:cpp:func:`vex::multivector::operator()()` allows to access individual components
of a multivector.

Some examples:

.. code-block:: cpp

    VEX_FUNCTION(bool, between, (double, a)(double, b)(double, c),
        return a <= b && b <= c;
        );

    vex::Reductor<double, vex::SUM> sum(ctx);
    vex::SpMat<double> A(ctx, ... );
    std::array<double, 2> v = {6.0, 7.0};

    vex::multivector<double, 2> X(ctx, N), Y(ctx, N);

    // ...

    X = sin(v * Y + 1);             // X(k) = sin(v[k] * Y(k) + 1);
    v = sum( between(0, X, Y) );    // v[k] = sum( between( 0, X(k), Y(k) ) );
    X = A * Y;                      // X(k) = A * Y(k);

Some operations can not be expressed with simple multivector arithmetic. For
example, an operation of two dimensional rotation mixes components in the right
hand side expressions:

.. math::

    y_0 = x_0 \cos(\alpha) - x_1 \sin(\alpha),\\
    y_1 = x_0 \sin(\alpha) + x_1 \cos(\alpha).

This may in principle be implemented as:

.. code-block:: cpp

    double alpha;
    vex::multivector<double, 2> X(ctx, N), Y(ctx, N);

    Y(0) = X(0) * cos(alpha) - X(1) * sin(alpha);
    Y(1) = X(0) * sin(alpha) + X(1) * cos(alpha);

But this would result in two kernel launches instead of single fused launch.
VexCL allows one to assign a tuple of expressions to a multivector, which will
lead to the launch of a single fused kernel:

.. code-block:: cpp

    Y = std::make_tuple(
        X(0) * cos(alpha) - X(1) * sin(alpha),
        X(0) * sin(alpha) + X(1) * cos(alpha) );

:cpp:func:`vex::tie` function even allows to get rid of multivectors completely
and fuse several vector expressions into a single kernel. We can rewrite the
above examples with just :cpp:class:`vex::vector\<double>` instances:

.. code-block:: cpp

    vex::vector<double> x0(ctx, N), x1(ctx,N), y0(ctx, N), y1(ctx, N);

    vex::tie(y0, y1) = std::make_tuple(
        x0 * cos(alpha) - x1 * sin(alpha),
        x0 * sin(alpha) + x1 * cos(alpha) );

.. doxygenclass:: vex::multivector
    :members:

.. doxygenfunction:: vex::tie
