.. cpp:namespace:: enoki

Spherical harmonics
===================

Enoki can efficiently evaluate real spherical harmonics basis functions for
both scalar and vector arguments. To use this feature, include the following
header file:

.. code-block:: cpp

    #include <enoki/sh.h>

The evaluation routines rely on efficient pre-generated branch-free code
processed using aggressive constant folding and common subexpression
elimination passes. Evaluation routines are provided up to order 10.

The generated code is based on the paper `Efficient Spherical Harmonic
Evaluation <http://jcgt.org/published/0002/02/06/>`_, *Journal of Computer
Graphics Techniques (JCGT)*, vol. 2, no. 2, 84-90, 2013 by `Peter-Pike Sloan
<http://www.ppsloan.org/publications/>`_.

.. note::

    The directions provided to ``sh_eval_*`` must be normalized 3D vectors
    (i.e. using Cartesian instead of spherical coordinates).

    The Mathematica equivalent of the real spherical harmonic basis implemented
    in :file:`enoki/sh.h` is given by the following definition:

    .. code-block:: wolfram-language

        SphericalHarmonicQ[l_, m_, d_] := Block[{θ, ϕ},
          θ = ArcCos[d[[3]]];
          ϕ = ArcTan[d[[1]], d[[2]]];
          Piecewise[{
            {SphericalHarmonicY[l, m, θ, ϕ], m == 0},
            {Sqrt[2] * Re[SphericalHarmonicY[l,  m, θ, ϕ]], m > 0},
            {Sqrt[2] * Im[SphericalHarmonicY[l, -m, θ, ϕ]], m < 0}
          }]
        ]

Usage
-----

The following example shows how to evaluate the spherical harmonics basis up to
and including order 2 producing a total of 9 function evaluations.

.. code-block:: cpp

    using Vector3f = Array<float, 3>;
    Vector3f d = normalize(Vector3f(1, 2, 3));

    float coeffs[9];
    sh_eval(d, 2, coeffs);

    // Prints: [0.282095, -0.261169, 0.391754, -0.130585, 0.156078, -0.468235, 0.292864, -0.234118, -0.117059]
    std::cout << load<Array<float, 9>>(coeffs) << std::endl;


Reference
---------

.. cpp:function:: template <typename Array> void sh_eval(const Array &d, size_t order, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order ``order``. The output array must have room for ``(order + 1)*(order +
    1)`` entries. This function dispatches to one of the ``sh_eval_*``
    implementations and throws an exception if ``order > 9``.

.. cpp:function:: template <typename Array> void sh_eval_0(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 0. The output array must have room for ``1`` entry.

.. cpp:function:: template <typename Array> void sh_eval_1(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 1. The output array must have room for ``4`` entries.

.. cpp:function:: template <typename Array> void sh_eval_2(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 2. The output array must have room for ``9`` entries.

.. cpp:function:: template <typename Array> void sh_eval_3(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 3. The output array must have room for ``16`` entries.

.. cpp:function:: template <typename Array> void sh_eval_4(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 4. The output array must have room for ``25`` entries.

.. cpp:function:: template <typename Array> void sh_eval_5(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 5. The output array must have room for ``36`` entries.

.. cpp:function:: template <typename Array> void sh_eval_6(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 6. The output array must have room for ``49`` entries.

.. cpp:function:: template <typename Array> void sh_eval_7(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 7. The output array must have room for ``64`` entries.

.. cpp:function:: template <typename Array> void sh_eval_8(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 8. The output array must have room for ``81`` entries.

.. cpp:function:: template <typename Array> void sh_eval_9(const Array &d, expr_t<value_t<Array>> *out)

    Evaluates the real spherical harmonics basis functions up to and including
    order 9. The output array must have room for ``100`` entries.
