.. _random:
.. cpp:namespace:: enoki

Random number generation
========================

Enoki ships with a vectorized implementation of the `PCG32 random number
generator <http://www.pcg-random.org/>`_ developed by `Melissa O'Neill
<https://www.cs.hmc.edu/~oneill>`_. To use it, include the following header:

.. code-block:: cpp

    #include <enoki/random.h>

The following reference is based on the original `PCG32 documentation
<http://www.pcg-random.org/using-pcg-c.html>`_.

Usage
-----

The :cpp:class:`enoki::PCG32` class takes a single template parameter ``T``
that denotes the "shape" of the output. This can be any scalar type like
``uint32_t``, in which case the implementation generates scalar variates:

.. code-block:: cpp

    /* Scalar RNG */
    using RNG_1x = PCG32<uint32_t>;

    RNG_1x my_rng;
    float value = my_rng.next_float32();

Alternatively, it can be an Enoki array, in which case the implementation
produces arrays of variates.

.. code-block:: cpp

    using FloatP = Packet<float, 16>;

    /* Vector RNG -- generates 16 independent variates at once */
    using RNG_16x = PCG32<FloatP>;

    RNG_16x my_rng;
    FloatP value = my_rng.next_float32();

PCG32 is *fast*: on a Skylake i7-6920HQ processor, the vectorized
implementation provided here generates around 1.4 billion single precision
variates per second.

Reference
---------

.. cpp:class:: template <typename T> PCG32

    This class implements the PCG32 pseudorandom number generator. It has a
    period of :math:`2^{64}` and supports :math:`2^{63}` separate *streams*.
    Each stream produces a different unique sequence of random numbers, which
    is particularly useful in the context of vectorized computations.

Member types
************

.. cpp:namespace:: template <typename T> enoki::PCG32

.. cpp:member:: constexpr size_t Size

    Denotes the SIMD width of the random number generator (i.e. how many
    pseudorandom variates are generated at the same time)

.. cpp:type:: Int64 = int64_array_t<T>

    Type alias for a signed 64-bit integer (or an array thereof).

.. cpp:type:: UInt64 = uint64_array_t<T>

    Type alias for a unsigned 64-bit integer (or an array thereof).

.. cpp:type:: UInt32 = uint32_array_t<T>

    Type alias for a unsigned 32-bit integer (or an array thereof).

.. cpp:type:: Float32 = float32_array_t<T>

    Type alias for a single precision float (or an array thereof).

.. cpp:type:: Float64 = float64_array_t<T>

    Type alias for a double precision float (or an array thereof).

Member variables
****************

.. cpp:member:: UInt64 state

    Stores the RNG state.  All values are possible.

.. cpp:member:: UInt64 inc

    Controls which RNG sequence (stream) is selected. Must *always* be odd,
    which is ensured by the constructor and :cpp:func:`seed()` method.

Constructors
************

.. cpp:function:: PCG32(const UInt64 &initstate = PCG32_DEFAULT_STATE, \
                        const UInt64 &initseq = PCG32_DEFAULT_STREAM + arange<UInt64>())

     Seeds the PCG32 with the default state. When ``T`` is an array, every
     entry by default uses a different stream index, which yields an
     uncorrelated and non-overlapping set of sequences.

Methods
*******

.. cpp:function:: void seed(const UInt64 &initstate, const UInt64 &initseq)

    This function initializes (a.k.a. "seeds") the random number generator, a
    required initialization step before the generator can be used. The provided
    arguments are defined as follows:

    - ``initstate`` is the starting state for the RNG. Any 64-bit value is
      permissible.

    - ``initseq`` selects the output sequence for the RNG. Any 64-bit value is
      permissible, although only the low 63 bits are used.

    For this generator, there are :math:`2^{63}` possible sequences of
    pseudorandom numbers. Each sequence is entirely distinct and has a period
    of :math:`2^{64}`. The ``initseq`` argument selects which stream is used.
    The ``initstate`` argument specifies the location within the :math:`2^{64}`
    period.

    Calling :cpp:func:`PCG32::seed` with the same arguments produces the same
    output, allowing programs to use random number sequences repeatably.

.. cpp:function:: UInt32 next_uint32(const mask_t<UInt64> &mask = true)

    Generate a uniformly distributed unsigned 32-bit random number (i.e.
    :math:`x`, where :math:`0\le x< 2^{32}`)

    If a mask parameter is provided, only the pseudorandom number generators
    of active SIMD lanes are advanced.

.. cpp:function:: UInt64 next_uint64(const mask_t<UInt64> &mask = true)

    Generate a uniformly distributed unsigned 64-bit random number (i.e.
    :math:`x`, where :math:`0\le x< 2^{64}`)

    If a mask parameter is provided, only the pseudorandom number generators
    of active SIMD lanes are advanced.

    .. note::

        This function performs two internal calls to :cpp:func:`next_uint32()`.

.. cpp:function:: UInt32 next_uint32_bound(uint32_t bound, const mask_t<UInt64> &mask = true)

    Generate a uniformly distributed unsigned 32-bit random number less
    than ``bound`` (i.e. :math:`x`, where :math:`0\le x<` ``bound``)

    If a mask parameter is provided, only the pseudorandom number generators
    of active SIMD lanes are advanced.

    .. note::

        This may involve multiple internal calls to
        :cpp:func:`next_uint32()`, in which case the RNG advances by
        several steps. This is only relevant when using the
        :cpp:func:`advance()` or :cpp:func:`operator-()` method.

.. cpp:function:: UInt64 next_uint64_bound(uint64_t bound, const mask_t<UInt64> &mask = true)

    Generate a uniformly distributed unsigned 64-bit random number less
    than ``bound`` (i.e. :math:`x`, where :math:`0\le x<` ``bound``)

    If a mask parameter is provided, only the pseudorandom number generators of
    active SIMD lanes are advanced.

    .. note::

        This may involve multiple internal calls to
        :cpp:func:`next_uint64()`, in which case the RNG advances by
        several steps. This is only relevant when using the
        :cpp:func:`advance()` or :cpp:func:`operator-()` method.

.. cpp:function:: Float32 next_float32(const mask_t<UInt64> &mask = true)

    Generate a single precision floating point value on the interval :math:`[0, 1)`

    If a mask parameter is provided, only the pseudorandom number generators of
    active SIMD lanes are advanced.

.. cpp:function:: Float64 next_float64(const mask_t<UInt64> &mask = true)

    Generate a double precision floating point value on the interval :math:`[0, 1)`

    If a mask parameter is provided, only the pseudorandom number generators of
    active SIMD lanes are advanced.

    .. warning::

        Since the underlying random number generator produces 32 bit
        output, only the first 32 mantissa bits will be filled (however,
        the resolution is still finer than in :cpp:func:`next_float32`,
        which only uses 23 mantissa bits)

.. cpp:function:: void advance(const Int64 &delta)

    This operation provides jump-ahead; it advances the RNG by ``delta`` steps,
    doing so in :math:`\log(\texttt{delta})` time. Because of the periodic
    nature of generation, advancing by :math:`2^{64}-d` (i.e., passing
    :math:`-d`) is equivalent to backstepping the generator by :math:`d` steps.

.. cpp:function:: Int64 operator-(const PCG32 &other)

    Compute the distance between two PCG32 pseudorandom number generators

.. cpp:function:: bool operator==(const PCG32 &other)

    Equality operator

.. cpp:function:: bool operator!=(const PCG32 &other)

    Inequality operator

Macros
******

The following macros are defined in :file:`enoki/random.h`:

.. cpp:var:: uint64_t PCG32_DEFAULT_STATE = 0x853c49e6748fea9bULL

    Default initialization passed to :cpp:func:`PCG32::seed`.

.. cpp:var:: uint64_t PCG32_DEFAULT_STREAM = 0xda3e39cb94b95bdbULL

    Default stream index passed to :cpp:func:`PCG32::seed`.
