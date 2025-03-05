.. cpp:namespace:: enoki

Half-precision floats
=====================

Enoki provides a compact implementation of a 16 bit *half-precision* floating
point type that is compatible with the FP16 format on GPUs and high dynamic
range image libraries such as OpenEXR. To use this feature, include the
following header:

.. code-block:: cpp

    #include <enoki/half.h>


Current processors don't natively implement half precision arithmetic, hence
mathematical operations involving this type always involve a
``half``:math:`\to` ``float``:math:`\to` ``half`` roundtrip. For this reason,
it is unwise to rely on it for expensive parts of a computation.

The main reason for including a dedicated half precision type in Enoki is that
it provides an ideal storage format for floating point data that does not
require the full accuracy of the single precision representation, which leads
to an immediate storage savings of :math:`2\times`.

.. note::

    If supported by the target architecture, Enoki uses the *F16C* instruction
    set to perform efficient vectorized conversion between half and single
    precision variables (however, this only affects conversion and no other
    arithmetic operations). ARM NEON also provides native conversion
    instructions.

Usage
-----

The following example shows how to use the :cpp:class:`enoki::half` type in a
typical use case.

.. code-block:: cpp

    using Color4f = Array<float, 4>;
    using Color4h = Array<half, 4>;

    uint8_t *image_ptr = ...;

    Color4f pixel(load<Color4h>(image_ptr)); // <- conversion vectorized using F16C

    /* ... update 'pixel' using single-precision arithmetic ... */

    store(image_ptr, Color4h(pixel)); // <- conversion vectorized using F16C

Reference
---------

.. cpp:class:: half

    A :cpp:class:`half` instance encodes a sign bit, an exponent width of 5
    bits, and 10 explicitly stored mantissa bits.

    All standard mathematical operators are overloaded and implemented using
    the processor's floating point unit after a conversion to a IEEE754 single
    precision. The result of the operation is then converted back to half
    precision.

    .. cpp:var:: uint16_t value

        Stores the represented half precision value as an unsigned 16-bit integer.

    .. cpp:function:: half(float value)

        Constructs a half-precision value from the given single precision
        argument.

    .. cpp:function:: operator float() const

        Implicit ``half`` to ``float`` conversion operator.

    .. cpp:function:: static half from_binary(uint16_t value)

        Reinterpret a 16-bit unsigned integer as a half-precision variable.

    .. cpp:function:: half operator+(half h) const

        Addition operator.

    .. cpp:function:: half& operator+=(half h)

        Addition compound assignment operator.

    .. cpp:function:: half operator-() const

        Unary minus operator

    .. cpp:function:: half operator*(half h) const

        Multiplication operator.

    .. cpp:function:: half& operator*=(half h)

        Multiplication compound assignment operator.

    .. cpp:function:: half operator/(half h) const

        Division operator.

    .. cpp:function:: half& operator/=(half h)

        Division compound assignment operator.

    .. cpp:function:: bool operator<(half h) const

        Less-than comparison operator.

    .. cpp:function:: bool operator<=(half h) const

        Less-than-or-equal comparison operator.

    .. cpp:function:: bool operator>(half h) const

        Greater-than comparison operator.

    .. cpp:function:: bool operator>=(half h) const

        Greater-than-or-equal comparison operator.

    .. cpp:function:: bool operator==(half h) const

        Equality operator.

    .. cpp:function:: bool operator!=(half h) const

        Inequality operator.

    .. cpp:function:: friend std::ostream& operator<<(std::ostream &os, const half &h)

        Stream insertion operator.
