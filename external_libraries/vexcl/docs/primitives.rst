Parallel primitives and algorithms
==================================

Reductions
----------

An instance of :cpp:class:`vex::Reductor\<T, ReduceOP=vex::SUM>` allows one to
reduce an arbitrary vector expression to a single value of type T. Supported
reduction operations are :cpp:class:`vex::SUM`, :cpp:class:`vex::MIN`, and
:cpp:class:`vex::MAX`. Reductor objects are steteful -- they keep small
temporary buffers on compute devices and receive a list of command queues at
construction.

In the following example an inner product of two vectors is computed:

.. code-block:: cpp

    vex::Reductor<double, vex::SUM> sum(ctx);
    double s = sum(x * y);

Also, see :ref:`random-number-generation` for an example of estimating value of
:math:`\pi` with the Monte Carlo method.

Reduce operations may be combined with the
:cpp:class:`vex::CombineReductors` class. This way several
reduction operations will be fused into single compute kernel. The operations
should return the same scalar type, and the result of the combined reduction
operation will be appropriately sized OpenCL/CUDA vector type.

In the following example minimum and maximum values of the vector are computed
at the same time:

.. code-block:: cpp

    vex::Reductor<double, vex::CombineReductors<vex::MIN, vex::MAX>> minmax(ctx);
    cl_double2 m = minmax(x);
    std::cout << "min(x) = " << m.s[0] << std::endl;
    std::cout << "max(x) = " << m.s[1] << std::endl;

In fact, the operation is so common, that VexCL provides a convenience typedef
:cpp:class:`vex::MIN_MAX`.

.. doxygenclass:: vex::Reductor
    :members:

.. doxygenstruct:: vex::SUM
.. doxygenstruct:: vex::MIN
.. doxygenstruct:: vex::MAX
.. doxygenstruct:: vex::CombineReductors
.. doxygentypedef:: vex::MIN_MAX

Sparse matrix-vector products
-----------------------------

One of the most common operations in linear algebra is the matrix-vector
product. An instance of :cpp:class:`vex::SpMat` class holds a representation of
a sparse matrix. Its constructor accepts a sparse matrix in common CRS_ format.
In the example below a `vex::SpMat` is constructed from an Eigen_ sparse
matrix:

.. _CRS: http://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR_or_CRS.29
.. _Eigen: http://eigen.tuxfamily.org/

.. code-block:: cpp

    Eigen::SparseMatrix<double, Eigen::RowMajor, int> E;

    vex::SpMat<double, int> A(ctx, E.rows(), E.cols(),
        E.outerIndexPtr(), E.innerIndexPtr(), E.valuePtr());


Matrix-vector products may be used in vector expressions. The only restriction
is that the expressions have to be additive. This is due to the fact that the
operation involves inter-device communication for multi-device contexts.

.. code-block:: cpp

    // Compute residual value for a system of linear equations:
    Z = Y - A * X;

This restriction may be lifted for single-device contexts. In this case VexCL
does not need to worry about inter-device communication. Hence, it is possible
to inline matrix-vector product into a normal vector expression with the help of
:cpp:func:`vex::make_inline`:

.. code-block:: cpp

    residual = sum(Y - vex::make_inline(A * X));
    Z = sin(vex::make_inline(A * X));

.. doxygenclass:: vex::SpMat
    :members:

.. doxygenfunction:: vex::make_inline(const MVProdExpr&)

Sort, scan, reduce-by-key algorithms
------------------------------------

VexCL provides several standalone parallel primitives that may not be used as
part of a vector expression. These are :cpp:func:`vex::inclusive_scan_by_key`,
:cpp:func:`vex::exclusive_scan_by_key`, :cpp:func:`vex::sort`,
:cpp:func:`vex::sort_by_key`, :cpp:func:`vex::reduce_by_key`. All of these
functions take VexCL vectors both as input and output parameters.

Sort and scan functions take an optional function object used for comparison
and summing of elements. The functor should provide the same interface as, e.g.
``std::less`` for sorting or ``std::plus`` for summing; additionally, it should
provide a VexCL function for device-side operations.

Here is an example of such an object comparing integer elements in such a way
that even elements precede odd ones:

.. code-block:: cpp

    template <typename T>
    struct even_first {
        #define BODY                        \
            char bit1 = 1 & a;              \
            char bit2 = 1 & b;              \
            if (bit1 == bit2) return a < b; \
            return bit1 < bit2;

        // Device version.
        VEX_FUNCTION(bool, device, (int, a)(int, b), BODY);

        // Host version.
        bool operator()(int a, int b) const { BODY }

        #undef BODY
    };

Same functor could be created with the help of :c:macro:`VEX_DUAL_FUNCTOR`
macro, which takes return type, sequence of arguments (similar to the
:c:macro:`VEX_FUNCTION`), and the body of the functor:

.. code-block:: cpp

    template <typename T>
    struct even_first {
        VEX_DUAL_FUNCTOR(bool, (T, a)(T, b),
            char bit1 = 1 & a;
            char bit2 = 1 & b;
            if (bit1 == bit2) return a < b;
            return bit1 < bit2;
        )
    };

Note that VexCL already provides :cpp:class:`vex::less\<T>`,
:cpp:class:`vex::less_equal\<T>`, :cpp:class:`vex::greater\<T>`,
:cpp:class:`vex::greater_equal\<T>`, and :cpp:class:`vex::plus\<T>`.

The need to provide both host-side and device-side parts of the functor comes
from the fact that multidevice vectors are first sorted partially on each of
the compute devices and then merged on the host.

Sorting algorithms may also take tuples of keys/values (in fact, any
Boost.Fusion_ sequence will do).  One will have to explicitly specify the
comparison functor in this case. Both host and device variants of the
comparison functor should take ``2n`` arguments, where ``n`` is the number of
keys.  The first ``n`` arguments correspond to the left set of keys, and the
second ``n`` arguments correspond to the right set of keys. Here is an example
that sorts values by a tuple of two keys:

.. code-block:: cpp

    vex::vector<int>    keys1(ctx, n);
    vex::vector<float>  keys2(ctx, n);
    vex::vector<double> vals (ctx, n);

    struct {
        VEX_FUNCTION(bool, device, (int, a1)(float, a2)(int, b1)(float, b2),
                return (a1 == b1) ? (a2 < b2) : (a1 < b1);
                );
        bool operator()(int a1, float a2, int b1, float b2) const {
            return std::make_tuple(a1, a2) < std::make_tuple(b1, b2);
        }
    } comp;

    vex::sort_by_key(std::tie(keys1, keys2), vals, comp);

.. _Boost.Fusion: http://www.boost.org/doc/libs/release/libs/fusion/doc/html/index.html

.. doxygenfunction:: vex::inclusive_scan(vector<T> const&, vector<T>&, T, Oper)
.. doxygenfunction:: vex::exclusive_scan(vector<T> const&, vector<T>&, T, Oper)
.. doxygenfunction:: vex::inclusive_scan_by_key(K&&, const vector<V>&, vector<V>&, Comp, Oper, V)
.. doxygenfunction:: vex::exclusive_scan_by_key(K&&, const vector<V>&, vector<V>&, Comp, Oper, V)
.. doxygenfunction:: vex::sort(K&&, Comp)
.. doxygenfunction:: vex::sort_by_key(K&&, V&&, Comp)
.. doxygendefine:: VEX_DUAL_FUNCTOR
.. doxygenstruct:: vex::less
.. doxygenstruct:: vex::less_equal
.. doxygenstruct:: vex::greater
.. doxygenstruct:: vex::greater_equal
.. doxygenstruct:: vex::plus
