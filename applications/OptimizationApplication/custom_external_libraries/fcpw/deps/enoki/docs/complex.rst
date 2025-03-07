.. cpp:namespace:: enoki

Complex numbers
===============

Enoki provides a vectorizable type for complex number arithmetic analogous to
``std::complex<T>``. To use it, include the following header file:

.. code-block:: cpp

    #include <enoki/complex.h>

Usage
-----

The following example shows how to define and perform basic arithmetic using
:cpp:class:`enoki::Complex` vectorized over 4-wide packets.

.. code-block:: cpp

    /* Declare underlying packet type, could just be 'float' for scalar arithmetic */
    using FloatP   = Packet<float, 4>;

    /* Define complex number type */
    using ComplexP = Complex<FloatP>;

    const ComplexP I(0.f, 1.f);
    const ComplexP z(0.2f, 0.3f);

    /* Two different ways of computing the tangent function */
    ComplexP t0 = (exp(I * z) - exp(-I * z)) / (I * (exp(I * z) + exp(-I * z)));
    ComplexP t1 = tan(z);

    std::cout << "t0 = " << t0 << std::endl << std::endl;
    std::cout << "t1 = " << t1 << std::endl;

    /* Prints

        t0 = [0.184863 + 0.302229i,
         0.184863 + 0.302229i,
         0.184863 + 0.302229i,
         0.184863 + 0.302229i]

        t1 = [0.184863 + 0.302229i,
         0.184863 + 0.302229i,
         0.184863 + 0.302229i,
         0.184863 + 0.302229i]
     */

Reference
---------

.. cpp:class:: template <typename Type> Complex : StaticArrayImpl<Type, 2>

    The class :cpp:class:`enoki::Complex` is a 2D Enoki array whose components
    are of type ``Type``. Various arithmetic operators (e.g. multiplication)
    and transcendental functions are overloaded so that they provide the
    correct behavior for complex-valued inputs.

    .. cpp:function:: Complex(Type real, Type imag)

        Creates a :cpp:class:`enoki::Complex` instance from the given real and
        imaginary inputs.

    .. cpp:function:: Complex(Type f)

        Creates a real-valued :cpp:class:`enoki::Complex` instance from ``f``.
        This constructor effectively changes the broadcasting behavior of
        non-complex inputs---for instance, the snippet

        .. code-block:: cpp

            auto value_a = zero<Array<float, 2>>();
            auto value_c = zero<Complex<float>>();

            value_a += 1.f; value_c += 1.f;

            std::cout << "value_a = "<< value_a << ", value_c = " << value_c << std::endl;

        prints ``value_a = [1, 1], value_c = 1 + 0i``, which is the desired
        behavior for complex numbers. For standard Enoki arrays, the number
        ``1.f`` is broadcast to both components.

Elementary operations
*********************

.. cpp:function:: template <typename T> T real(Complex<T> z)

    Extracts the real part of ``z``.

.. cpp:function:: template <typename T> T imag(Complex<T> z)

    Extracts the imaginary part of ``z``.

.. cpp:function:: template <typename T> Complex<T> arg(Complex<T> z)

    Evaluates the complex argument of ``z``.

.. cpp:function:: template <typename T> Complex<T> abs(Complex<T> z)

    Compute the absolute value of ``z``.

.. cpp:function:: template <typename T> Complex<T> sqrt(Complex<T> z)

    Compute the square root of ``z``.

.. cpp:function:: template <typename T> Complex<T> conj(Complex<T> z)

    Evaluates the complex conjugate of ``z``.

.. cpp:function:: template <typename T> Complex<T> rcp(Complex<T> z)

    Evaluates the complex reciprocal of ``z``.

Arithmetic operators
********************

Only a few arithmetic operators need to be overridden to support complex
arithmetic. The rest are automatically provided by Enoki's existing operators
and broadcasting rules.

.. cpp:function:: template <typename T> Complex<T> operator*(Complex<T> z0, Complex<T> z1)

    Evaluates the complex product of ``z1`` and ``z2``.

.. cpp:function:: template <typename T> Complex<T> operator/(Complex<T> z0, Complex<T> z1)

    Evaluates the complex division of ``z1`` and ``z2``.

Stream operators
****************

.. cpp:function:: std::ostream& operator<<(std::ostream &os, const Complex<T> &z)

    Sends the complex number ``z`` to the stream ``os`` using the format
    ``1 + 2i``.


Exponential, logarithm, and power function
******************************************

.. cpp:function:: template <typename T> Complex<T> exp(Complex<T> z)

    Evaluates the complex exponential of ``z``.

.. cpp:function:: template <typename T> Complex<T> log(Complex<T> z)

    Evaluates the complex logarithm of ``z``.

.. cpp:function:: template <typename T> Complex<T> pow(Complex<T> z0, Complex<T> z1)

    Evaluates the complex power of ``z0`` raised to the ``z1``.

Trigonometric functions
***********************

.. cpp:function:: template <typename T> Complex<T> sin(Complex<T> z)

    Evaluates the complex sine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> cos(Complex<T> z)

    Evaluates the complex cosine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> tan(Complex<T> z)

    Evaluates the complex tangent function for ``z``.

.. cpp:function:: template <typename T> std::pair<Complex<T>, Complex<T>> sincos(Complex<T> z)

    Jointly evaluates the complex sine and cosine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> asin(Complex<T> z)

    Evaluates the complex arc sine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> acos(Complex<T> z)

    Evaluates the complex arc cosine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> atan(Complex<T> z)

    Evaluates the complex arc tangent function for ``z``.

Hyperbolic functions
********************

.. cpp:function:: template <typename T> Complex<T> sinh(Complex<T> z)

    Evaluates the complex hyperbolic sine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> cosh(Complex<T> z)

    Evaluates the complex hyperbolic cosine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> tanh(Complex<T> z)

    Evaluates the complex hyperbolic tangent function for ``z``.

.. cpp:function:: template <typename T> std::pair<Complex<T>, Complex<T>> sincosh(Complex<T> z)

    Jointly evaluates the complex hyperbolic sine and cosine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> asinh(Complex<T> z)

    Evaluates the complex hyperbolic arc sine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> acosh(Complex<T> z)

    Evaluates the complex hyperbolic arc cosine function for ``z``.

.. cpp:function:: template <typename T> Complex<T> atanh(Complex<T> z)

    Evaluates the complex hyperbolic arc tangent function for ``z``.

Miscellaneous functions
***********************

.. cpp:function:: std::pair<T, T> sincos_arg_diff(const Complex<T> &z1, const Complex<T> &z2)

   Efficiently evaluates ``sin(arg(z1) - arg(z2))`` and ``cos(arg(z1) - arg(z2))``.

.. cpp:function:: template <typename T> Complex<T> sqrtz(T x)

    Compute the complex square root of a real-valued argument ``x`` (which may
    be negative). This is considerably more efficient than the general complex
    square root above.
