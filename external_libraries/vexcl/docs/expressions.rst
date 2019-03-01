Vector expressions
==================

VexCL allows the use of convenient and intuitive notation for vector
operations. In order to be used in the same expression, all participating
vectors have to be `compatible`:

- Have same size;
- Span same set of compute devices.

If these conditions are satisfied, then vectors may be combined with rich set
of available expressions. Vector expressions are processed in parallel across
all devices they were allocated on. Each vector expression results in the
launch of a single compute kernel. The kernel is automatically generated and
compiled the first time the expression is encountered in the program, and is
submitted to command queues associated with the vector that is being assigned
to.

VexCL will dump the sources of the generated kernels to stdout if either the
:c:macro:`VEXCL_SHOW_KERNELS` preprocessor macro is defined, or there exists
``VEXCL_SHOW_KERNELS`` environment variable. For example, the expression:

.. code-block:: cpp

    X = 2 * Y - sin(Z);

will lead to the launch of the following compute kernel:

.. code-block:: cpp

    kernel void vexcl_vector_kernel(
        ulong n,
        global double * prm_1,
        int prm_2,
        global double * prm_3,
        global double * prm_4
    )
    {
        for(size_t idx = get_global_id(0); idx < n; idx += get_global_size(0)) {
            prm_1[idx] = ( ( prm_2 * prm_3[idx] ) - sin( prm_4[idx] ) );
        }
    }

Here and in the rest of examples ``X``, ``Y``, and ``Z`` are compatible
instances of :cpp:class:`vex::vector\<double>`; it is also assumed that OpenCL
backend is selected.

VexCL is able to cache the compiled kernels offline. The compiled binaries are
stored in ``$HOME/.vexcl`` on Linux and MacOSX, and in ``%APPDATA%\vexcl`` on
Windows systems. In order to enable this functionality for OpenCL-based
backends, the user has to define the :c:macro:`VEXCL_CACHE_KERNELS` prprocessor
macro. NVIDIA OpenCL implementation does the caching already, but on AMD or
Intel platforms this may lead to dramatic decrease of program initialization
time (e.g. VexCL tests take around 20 seconds to complete without kernel
caches, and 2 seconds when caches are available). In case of the CUDA backend
the offline caching is always enabled.

.. c:macro:: VEXCL_SHOW_KERNELS

    When defined, VexCL will dump source code of the generated kernels to
    stdout. Same effect may be achieved by exporting an environment variable
    with the same name.

.. c:macro:: VEXCL_CACHE_KERNELS

    When defined, VexCL will use offline cache to store the compiled kernels.
    The first time a kernel is compiled on the system, its binaries are saved
    to the cache folder (``$HOME/.vexcl`` on Unix-like systems;
    ``%APPDATA%\vexcl`` on Windows). Next time the program is run, the binaries
    will be obtained from the cache, thus speeding up the program startup.

Builtin operations
------------------

VexCL expressions may combine device vectors and scalars with arithmetic,
logic, or bitwise operators as well as with builtin OpenCL/CUDA functions. If
some builtin operator or function is unavailable, it should be considered a
bug.  Please do not hesitate to open an issue in this case.

.. code-block:: cpp

    Z = sqrt(2 * X) + pow(cos(Y), 2.0);

Constants
---------

As you have seen above, ``2`` in the expression ``2 * Y - sin(Z)`` is passed to
the generated compute kernel as an ``int`` parameter (``prm_2``). Sometimes
this is desired behaviour, because the same kernel will be reused for the
expressions ``42 * Z - sin(Y)`` or ``a * Y - sin(Y)`` (where ``a`` is an
integer variable). But this may lead to a slight overhead if an expression
involves true constant that will always have same value. The
:c:macro:`VEX_CONSTANT` macro allows one to define such constants for use in
vector expressions. Compare the generated kernel for the following example with
the kernel above:

.. code-block:: cpp

    VEX_CONSTANT(two, 2);
    X = two() * Y - sin(Z);

.. code-block:: cpp

    kernel void vexcl_vector_kernel(
        ulong n,
        global double * prm_1,
        global double * prm_3,
        global double * prm_4
    )
    {
        for(ulong idx = get_global_id(0); idx < n; idx += get_global_size(0)) {
            prm_1[idx] = ( ( ( 2 ) * prm_3[idx] ) - sin( prm_4[idx] ) );
        }
    }

VexCL provides some predefined constants in the ``vex::constants`` namespace
that correspond to `boost::math::constants`_ (e.g.
:cpp:func:`vex::constants::pi`).

.. doxygendefine:: VEX_CONSTANT

.. _`boost::math::constants`: http://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/constants.html

Element indices
---------------

The function :cpp:func:`vex::element_index` allows one to use the index
of each vector element inside vector expressions. The numbering is continuous
across all compute devices and starts with an optional ``offset``.

.. code-block:: cpp

    // Linear function:
    double x0 = 0.0, dx = 1.0 / (N - 1);
    X = x0 + dx * vex::element_index();

    // Single period of sine function:
    Y = sin(vex::constants::two_pi() * vex::element_index() / N);

.. doxygenfunction:: vex::element_index

User-defined functions
----------------------

Users may define custom functions for use in vector expressions. One has to
define the function signature and the function body. The body may contain any
number of lines of valid OpenCL or CUDA code, depending on the selected
backend. The most convenient way to define a function is via the
:c:macro:`VEX_FUNCTION` macro:

.. code-block:: cpp

    VEX_FUNCTION(double, squared_radius, (double, x)(double, y),
        return x * x + y * y;
        );
    Z = sqrt(squared_radius(X, Y));

The first macro parameter here defines the function return type, the second
parameter is the function name, the third parameter defines function arguments
in form of a preprocessor sequence. Each element of the sequence is a tuple of
argument type and name. The rest of the macro is the function body (compare
this with how functions are defined in C/C++).  The resulting
``squared_radius`` function object is stateless; only its type is used for
kernel generation.  Hence, it is safe to define commonly used functions at the
global scope.

Note that any valid vector expression may be passed as a function parameter,
including nested function calls:

.. code-block:: cpp

    Z = squared_radius(sin(X + Y), cos(X - Y));


Another version of the macro takes the function body directly as a string:

.. code-block:: cpp

    VEX_FUNCTION_S(double, squared_radius, (double, x)(double, y),
        "return x * x + y * y;"
        );
    Z = sqrt(squared_radius(X, Y));

In case the function that is being defined calls other custom function inside
its body, one can use the version of the :c:macro:`VEX_FUNCTION` macro that
takes sequence of parent function names as the fourth parameter. This way the
kernel generator will know to include the function definitions into the kernel
source:

.. code-block:: cpp

    VEX_FUNCTION(double, bar, (double, x),
            double s = sin(x);
            return s * s;
            );
    VEX_FUNCTION(double, baz, (double, x),
            double c = cos(x);
            return c * c;
            );
    VEX_FUNCTION_D(double, foo, (double, x)(double, y), (bar)(baz),
            return bar(x - y) * baz(x + y);
            );


Similarly to :c:macro:`VEX_FUNCTION_S`, there is a version called
:c:macro:`VEX_FUNCTION_DS` (or symmetrical :c:macro:`VEX_FUNCTION_SD`) that
takes the function body as a string parameter.

Custom functions may be used not only for convenience, but also for performance
reasons. The above example with ``squared_radius`` could in principle be
rewritten as:

.. code-block:: cpp

    Z = sqrt(X * X + Y * Y);

The drawback of this version is that ``X`` and ``Y`` will be passed to the
kernel and read `twice` (see the next section for an explanation).

.. doxygendefine:: VEX_FUNCTION
.. doxygendefine:: VEX_FUNCTION_S
.. doxygendefine:: VEX_FUNCTION_D
.. doxygendefine:: VEX_FUNCTION_SD
.. doxygendefine:: VEX_STRINGIZE_SOURCE

Tagged terminals
----------------

The last example in the previous section is ineffective because the compiler
cannot tell if any two terminals in an expression tree are actually referring
to the same data. But programmers often have this information. VexCL allows one
to pass this knowledge to compiler by tagging terminals with unique tags.  By
doing this, the programmer guarantees that any two terminals with matching tags
are referencing the same data.

Below is a more effective variant of the above example:

.. code-block:: cpp

    using vex::tag;
    Z = sqrt(tag<1>(X) * tag<1>(X) + tag<2>(Y) * tag<2>(Y));

Here, the generated kernel will have one parameter for each of the vectors
``X`` and ``Y``:

.. code-block:: cpp

    kernel void vexcl_vector_kernel(
      ulong n,
      global double * prm_1,
      global double * prm_tag_1_1,
      global double * prm_tag_2_1
    )
    {
      for(ulong idx = get_global_id(0); idx < n; idx += get_global_size(0)) {
	prm_1[idx] = sqrt( ( ( prm_tag_1_1[idx] * prm_tag_1_1[idx] )
				+ ( prm_tag_2_1[idx] * prm_tag_2_1[idx] ) ) );
      }
    }

.. doxygenfunction:: vex::tag

Temporary values
----------------

Some expressions may have several occurences of the same subexpression.
Unfortunately, VexCL is not able to determine these cases without the
programmer's help. For example, let us consider the following expression:

.. code-block:: cpp

    Y = log(X) * (log(X) + Z);

Here, ``log(X)`` would be computed twice. One could tag vector ``X`` as in:

.. code-block:: cpp

    auto x = vex::tag<1>(X);
    Y = log(x) * (log(x) + Z);

and hope that the backend compiler is smart enough to reuse result of
``log(x)``.  In fact, most modern compilers will in this simple case. But in
harder cases it is possible to explicitly tell VexCL to store the result of a
subexpression in a local variable and reuse it. The
:cpp:func:`vex::make_temp\<size_t>` function template serves this purpose:

.. code-block:: cpp

    auto tmp1 = vex::make_temp<1>( sin(X) );
    auto tmp2 = vex::make_temp<2>( cos(X) );
    Y = (tmp1 - tmp2) * (tmp1 + tmp2);

This will result in the following kernel:

.. code-block:: cpp

    kernel void vexcl_vector_kernel(
      ulong n,
      global double * prm_1,
      global double * prm_2_1,
      global double * prm_3_1
    )
    {
      for(ulong idx = get_global_id(0); idx < n; idx += get_global_size(0))
      {
	double temp_1 = sin( prm_2_1[idx] );
	double temp_2 = cos( prm_3_1[idx] );
	prm_1[idx] = ( ( temp_1 - temp_2 ) * ( temp_1 + temp_2 ) );
      }
    }

Any valid vector or multivector expression (but not additive expressions, such
as sparse matrix-vector products) may be wrapped into a
:cpp:func:`vex::make_temp` call.

.. doxygenfunction:: vex::make_temp


Raw pointers [#sd]_
-------------------

Most of the expressions in VexCL are element-wise. That is, the user describes
what needs to be done on an element-by-element basis, and has no access to
neighboring elements. :cpp:func:`vex::raw_pointer` allows to use pointer
arithmetic with either :cpp:class:`vex::vector\<T>` or
:cpp:class:`vex::svm_vector\<T>`.

.. rubric:: The :math:`N`-body problem

Let us consider the :math:`N`-body problem as an example. The :math:`N`-body
problem considers :math:`N` point masses, :math:`m_i,\;i=1,2,\ldots,N` in three
dimensional space :math:`\mathbb{R}^3` moving under the influence of mutual
gravitational attraction. Each mass :math:`m_i` has a position vector
:math:`\vec q_i`. `Newton's law of gravity`_ says that the gravitational force
felt on mass :math:`m_i` by a single mass :math:`m_j` is given by

.. math::

    \vec F_{ij} = \frac{G m_i m_j (\vec q_j - \vec q_i)}
    {||\vec q_j - \vec q_i||^3},

where :math:`G` is the `gravitational constant`_ and
:math:`||\vec q_j - \vec q_i||` is the the distance between :math:`\vec q_i`
and :math:`\vec q_j`.

We can find the total force acting on mass :math:`m_i` by summing over all masses:

.. math::

    \vec F_i = \sum\limits_{j=1,j \neq i}^N \frac{G m_i m_j (\vec q_j - \vec q_i)}
    {||\vec q_j - \vec q_i||^3}.

.. _`Newton's law of gravity`: https://en.wikipedia.org/wiki/Newton%27s_law_of_gravity
.. _`gravitational constant`: https://en.wikipedia.org/wiki/Gravitational_constant

In VexCL, we can encode the formula above with the following custom function:

.. code-block:: cpp

    vex::vector<double>     m(ctx, n);
    vex::vector<cl_double3> q(ctx, n), f(ctx, n);

    VEX_FUNCTION(cl_double3, force, (size_t, n)(size_t, i)(double*, m)(cl_double3*, q),
        const double G = 6.674e-11;

        double3 sum = {0.0, 0.0, 0.0};
        double  m_i = m[i];
        double3 q_i = q[i];

        for(size_t j = 0; j < n; ++j) {
            if (j == i) continue;

            double  m_j = m[j];
            double3 d   = q[j] - q_i;
            double  r   = length(d);

            sum += G * m_i * m_j * d / (r * r * r);
        }
        return sum;
        );

    f = force(n, vex::element_index(), vex::raw_pointer(m), vex::raw_pointer(q));

The function takes number of elements ``n``, index of the current element
``i``, and pointers to arrays of point masses ``m`` and positions ``q``. It
returns the force acting on the current point. Note that we use host-side types
(``cl_double3``) in declaration of function return type and parameter types,
and we use OpenCL types (``double3``) inside the function body.

.. rubric:: Constant address space

In the OpenCL-based backends VexCL allows one to use constant cache on GPUs in
order to speed up the read-only access to small vectors. Usually around 64Kb of
constant cache per compute unit is available. Vectors wrapped in
:cpp:func:`vex::constant` will be decorated with the ``constant`` keyword
instead of the usual ``global`` one. For example, the following expression:

.. code-block:: cpp

    x = 2 * vex::constant(y);

will result in the OpenCL kernel below:

.. code-block:: cpp

    kernel void vexcl_vector_kernel(
      ulong n,
      global int * prm_1,
      int prm_2,
      constant int * prm_3
    )
    {
      for(ulong idx = get_global_id(0); idx < n; idx += get_global_size(0)) {
        prm_1[idx] = ( prm_2 * prm_3[idx] );
      }
    }

In cases where access to arbitrary vector elements is required,
:cpp:func:`vex::constant_pointer` may be used similarly to
:cpp:func:`vex::raw_pointer`. The extracted pointer will be decorated with the
``constant`` keyword.

.. doxygenfunction:: vex::raw_pointer(const vector<T>&)
.. doxygenfunction:: vex::raw_pointer(const svm_vector<T>&)
.. doxygenfunction:: vex::constant
.. doxygenfunction:: vex::constant_pointer

.. _random-number-generation:

Random number generation
------------------------

VexCL provides a counter-based random number generators from Random123_
suite, in which  N-th random number is obtained by applying a stateless mixing
function to N instead of the conventional approach of using N iterations of a
stateful transformation. This technique is easily parallelizable and is well
suited for use in GPGPU applications.

.. _Random123: http://www.deshawresearch.com/resources_random123.html

For integral types, the generated values span the complete range; for floating
point types, the generated values lie in the interval [0,1].

In order to use a random number sequence in a vector expression, the user has
to declare an instance of either :cpp:class:`vex::Random` or
:cpp:class:`vex::RandomNormal` class template as in the following example:

.. code-block:: cpp

    vex::Random<double, vex::random::threefry> rnd;

    // X will contain random numbers from [-1, 1]:
    X = 2 * rnd(vex::element_index(), std::rand()) - 1;

Note that :cpp:func:`vex::element_index` function here provides the random
number generator with a sequence position N, and ``std::rand()`` is used to
obtain a seed for this specific sequence.

.. rubric:: Monte Carlo :math:`\pi`

Here is a more interesting example of using random numbers to estimate the
value of :math:`\pi`. In order to do this we remember that area of a circle
with radius :math:`r` is equal to :math:`\pi r^2`. A square of the same
'radius' has area of :math:`(2r)^2`. Then we can write

.. math::

    \frac{\text{area of circle}}{\text{area of square}} =
    \frac{\pi r^2}{(2r)^2} = \frac{\pi}{4},

.. math::

    \pi = 4 \frac{\text{area of circle}}{\text{area of square}}

We can estimate the last fraction in the formula above with the Monte-Carlo
method. If we generate a lot of random points in a square, then ratio of circle
area over square area will be approximately equal to the ratio of points in the
circle over all points. This is illustrated by the following figure:

.. plot::
    :align: center

    from pylab import *

    chameleon1 = ( 51.0/255,  17.0/255,   1.0/255)
    chameleon2 = (255.0/255, 195.0/255,   0.0/255)
    chameleon3 = (199.0/255,  10.0/255,  30.0/255)
    chameleon4 = (107.0/255,   6.0/255,   5.0/255)

    figure(figsize=(4,4))

    x = np.random.uniform(0, 1, 3000)
    y = np.random.uniform(0, 1, 3000)

    i = (x * x + y * y) < 1.0

    plot(x[i==0], y[i==0], 'o', markersize=2, markeredgecolor=chameleon2, markerfacecolor=chameleon2)
    plot(x[i==1], y[i==1], 'o', markersize=2, markeredgecolor=chameleon3, markerfacecolor=chameleon3)

    x = np.linspace(0, 1, 100)
    y = np.sqrt(1 - x * x)

    plot(x, y, '-', color=chameleon4, linewidth=2)

    xticks([0, 0.5, 1])
    yticks([   0.5, 1])

In VexCL we can compute the estimate with a single expression, that will
generate single compute-bound kernel:

.. code-block:: cpp

    vex::Random<cl_double2> rnd;
    vex::Reductor<size_t> sum(ctx);

    double pi = sum(length(rnd(vex::element_index(0, n), std::rand())) < 1) * 4.0 / n;

Here we generate ``n`` random 2D points and use the builtin OpenCL function
``length`` to see which points are located withing the circle. Then we use the
``sum`` functor to count the points within the circle and finally multiply the
number with ``4.0/n`` to get the estimated value of :math:`\pi`.

.. doxygenstruct:: vex::Random
    :members:

.. doxygenstruct:: vex::RandomNormal
    :members:

Permutations [#sd]_
-------------------

:cpp:func:`vex::permutation` allows the use of a permuted vector in a vector
expression. The function accepts a vector expression that returns integral
values (indices).  The following example reverses `X` and assigns it to `Y` in
two different ways:

.. code-block:: cpp

    Y = vex::permutation(n - 1 - vex::element_index())(X);

    // Permutation expressions are writable!
    vex::permutation(n - 1 - vex::element_index())(Y) = X;

.. doxygenfunction:: vex::permutation

Slicing [#sd]_
--------------

An instance of the :cpp:class:`vex::slicer\<NDim>` class allows one to
conveniently access sub-blocks of multi-dimensional arrays that are stored in
:cpp:class:`vex::vector\<T>` in row-major order. The constructor of the class
accepts the dimensions of the array to be sliced. The following example
extracts every other element from interval ``[100, 200)`` of a
one-dimensional vector ``X``:

.. code-block:: cpp

    vex::vector<double> X(ctx, n);
    vex::vector<double> Y(ctx, 50);

    vex::slicer<1> slice(vex::extents[n]);

    Y = slice[vex::range(100, 2, 200)](X);

And the example below shows how to work with a two-dimensional matrix:

.. code-block:: cpp

    using vex::range, vex::_; // vex::_ is a shortcut for an empty range

    vex::vector<double> X(ctx, n * n); // n-by-n matrix stored in row-major order.
    vex::vector<double> Y(ctx, n);

    // vex::extents is a helper object similar to boost::multi_array::extents.
    vex::slicer<2> slice(vex::extents[n][n]);

    Y = slice[42](X);    // Put 42-nd row of X into Y.
    Y = slice[_][42](X); // Put 42-nd column of X into Y.

    slice[_][10](X) = Y; // Put Y into 10-th column of X.

    // Assign sub-block [10,20)x[30,40) of X to Z:
    vex::vector<double> Z = slice[range(10, 20)][range(30, 40)](X);
    assert(Z.size() == 100);

.. doxygenstruct:: vex::slicer
    :members:

.. doxygenvariable:: vex::extents
.. doxygenstruct:: vex::range
    :members:

.. doxygenvariable:: vex::_

Reducing multidimensional expressions [#sd]_
--------------------------------------------

:cpp:func:`vex::reduce` function allows one to reduce a multidimensional
expression along its one or more dimensions. The result is again a vector
expression, which may be used in other expressions. The supported reduction
operations are :cpp:class:`vex::SUM`, :cpp:class:`vex::MIN`, and
:cpp:class:`vex::MAX`. The function takes three arguments: the shape of the
expression to reduce (with the slowest changing dimension in the front), the
expression to reduce, and the dimension(s) to reduce along. The latter are
specified as indices into the shape array. Both the shape and indices are
specified as static arrays of integers, but :cpp:var:`vex::extents` object
may be used for convenience.

In the following example we find maximum absolute value of each row in a
two-dimensional matrix and assign the result to a vector:

.. code-block:: cpp

    vex::vector<double> A(ctx, N * M);
    vex::vector<double> x(ctx, N);

    x = vex::reduce<vex::MAX>(vex::extents[N][M], fabs(A), vex::extents[1]);

It is also possible to use :cpp:class:`vex::slicer` instance to provide
information about the expression shape:

.. code-block:: cpp

    vex::slicer<2> Adim(vex::extents[N][M]);
    x = vex::reduce<vex::MAX>(Adim[_](fabs(A)), vex::extents[1]);

.. doxygenfunction:: vex::reduce(const SlicedExpr&, const ReduceDims&)
.. doxygenfunction:: vex::reduce(const ExprShape&, const Expr&, const ReduceDims&)

Reshaping [#sd]_
----------------

:cpp:func:`vex::reshape` function is a powerful primitive that
allows one to conveniently manipulate multidimensional data. It takes three
arguments -- an arbitrary vector expression to reshape, the dimensions
``dst_dims`` of the final result (with the slowest changing dimension in the
front), and the native dimensions of the expression, which are specified as
indices into ``dst_dims``. The function returns a vector expression. The
dimensions may be conveniently specified with the help of
:cpp:var:`vex::extents` object.

Here is an example that shows how a two-dimensional matrix of size
:math:`N \times M` could be transposed:

.. code-block:: cpp

    vex::vector<double> A(ctx, N * M);
    vex::vector<double> B = vex::reshape(A,
                                vex::extents[M][N], // new shape
                                vex::extents[1][0]  // A is shaped as [N][M]
                                );

If the source expression lacks some of the destination dimensions, then those
will be introduced by replicating the available data. For example, to make a
two-dimensional matrix from a one-dimensional vector by copying the vector to
each row of the matrix, one could do the following:

.. code-block:: cpp

    vex::vector<double> x(ctx, N);
    vex::vector<double> y(ctx, M);
    vex::vector<double> A(ctx, M * N);

    // Copy x into each row of A:
    A = vex::reshape(x, vex::extents[M][N], vex::extents[1]);
    // Now, copy y into each column of A:
    A = vex::reshape(y, vex::extents[M][N], vex::extents[0]);

Here is a more realistic example of a dense matrix-matrix multiplication.
Elements of a matrix product :math:`C = AB` are defined as
:math:`C_{ij} = \sum_k A_{ik} B_{kj}`. Let's assume that matrix :math:`A` has
shape :math:`N \times L`, and matrix :math:`B` is shaped as :math:`L \times M`.
Then matrix :math:`C` has dimensions :math:`N \times M`. In order to implement
the multiplication we extend matrices :math:`A` and :math:`B` to the
shape of :math:`N \times L \times M`, multiply the resulting expressions
elementwise, and reduce the product along the middle dimension (:math:`L`):

.. code-block:: cpp

    vex::vector<double> A(ctx, N * L);
    vex::vector<double> B(ctx, L * M);
    vex::vector<double> C(ctx, N * M);

    C = vex::reduce<vex::SUM>(
            vex::extents[N][L][M],
            vex::reshape(A, vex::extents[N][L][M], vex::extents[0][1]) *
            vex::reshape(B, vex::extents[N][L][M], vex::extents[1][2]),
            1
            );

This of course would not be as efficient as a carefully crafted custom
implementation or a call to a vendor BLAS function. Also, this particular
operation is more efficiently done with tensor product function described in
the next section.

.. doxygenfunction:: vex::reshape

Tensor product [#sd]_
---------------------

Given two tensors (arrays of dimension greater than or equal to one), ``A`` and
``B``, and a list of axes pairs (where each pair represents corresponding axes
from each of the two tensors), the tensor product operation sums the products
of ``A``'s and ``B``'s elements over the given axes. In VexCL this is
implemented as :cpp:func:`vex::tensordot` operation (compare with python's
numpy.tensordot_).

For example, the above matrix-matrix product may be implemented much more
efficiently with :cpp:func:`vex::tensordot`:

.. code-block:: cpp

    using vex::_;

    vex::slicer<2> Adim(vex::extents[N][M]);
    vex::slicer<2> Bdim(vex::extents[M][L]);

    C = vex::tensordot(Adim[_](A), Bdim[_](B), vex::axes_pairs(1, 0));

Here instances of :cpp:class:`vex::slicer` class are used to provide shape
information for the ``A`` and ``B`` vectors.

.. _numpy.tensordot: http://docs.scipy.org/doc/numpy/reference/generated/numpy.tensordot.html

.. doxygenfunction:: vex::tensordot
.. doxygenfunction:: vex::axes_pairs

Scattered data interpolation with multilevel B-Splines
------------------------------------------------------

VexCL provides an implementation of the MBA algorithm based on paper by Lee,
Wolberg, and Shin [LeWS97]_. This is a fast algorithm for scattered
N-dimensional data interpolation and approximation.  Multilevel B-splines are
used to compute a C2-continuously differentiable surface through a set of
irregularly spaced points. The algorithm makes use of a coarse-to-fine
hierarchy of control lattices to generate a sequence of bicubic B-spline
functions whose sum approaches the desired interpolation function. Large
performance gains are realized by using B-spline refinement to reduce the sum
of these functions into one equivalent B-spline function.  High-fidelity
reconstruction is possible from a selected set of sparse and irregular samples.

The algorithm is setup on a CPU. After that, it may be used in vector
expressions. Here is an example in 2D:

.. code-block:: cpp

    // Coordinates of data points:
    std::vector< std::array<double,2> > coords = {
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.4, 0.4},
        {0.6, 0.6}
    };

    // Data values:
    std::vector<double> values = {
        0.2, 0.0, 0.0, -0.2, -1.0, 1.0
    };

    // Bounding box:
    std::array<double, 2> xmin = {-0.01, -0.01};
    std::array<double, 2> xmax = { 1.01,  1.01};

    // Initial grid size:
    std::array<size_t, 2> grid = {5, 5};

    // Algorithm setup.
    vex::mba<2> surf(ctx, xmin, xmax, coords, values, grid);

    // x and y are coordinates of arbitrary 2D points
    // (here the points are placed on a regular grid):
    vex::vector<double> x(ctx, n*n), y(ctx, n*n), z(ctx, n*n);

    auto I = vex::element_index() % n;
    auto J = vex::element_index() / n;
    vex::tie(x, y) = std::make_tuple(h * I, h * J);

    // Get interpolated values:
    z = surf(x, y);

.. [LeWS97] S. Lee, G. Wolberg, and S. Y. Shin. \
            Scattered data interpolation with multilevel B-Splines. \
            IEEE Transactions on Visualization and Computer Graphics, 3:228–244, 1997

.. doxygenclass:: vex::mba
    :members:

Fast Fourier Transform [#sd]_
-----------------------------

VexCL provides an implementation of the Fast Fourier Transform (FFT) that
accepts arbitrary vector expressions as input, allows one to perform
multidimensional transforms (of any number of dimensions), and supports
arbitrary sized vectors:

.. code-block:: cpp

    vex::FFT<double, cl_double2> fft(ctx, n);
    vex::FFT<cl_double2, double> ifft(ctx, n, vex::fft::inverse);

    vex::vector<double> rhs(ctx, n), u(ctx, n), K(ctx, n);

    // Solve Poisson equation with FFT:
    u = ifft( K * fft(rhs) );

.. doxygenstruct:: vex::FFT
    :members:

.. doxygenenum:: vex::fft::direction

.. [#sd]

    This operation involves access to arbitrary elements of its subexpressions
    and may lead to unpredictable device-to-device communication. Hence, it is
    restricted to single-device expressions. That is, only vectors that are
    located on a single device are allowed to participate in this operation.

