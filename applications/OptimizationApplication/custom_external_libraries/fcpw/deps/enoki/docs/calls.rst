.. cpp:namespace:: enoki
.. _calls:

Method calls
============

Method calls and virtual method calls are important building blocks of modern
object-oriented C++ applications. When vectorization enters the picture, it is
not immediately clear how they should be dealt with. This section introduces
Enoki's method call vectorization support, focusing on a hypothetical
``Sensor`` class that decodes a measurement performed by a sensor.

Note that the examples will refer to CPU SIMD-style vectorization, but
everything in this section also applies to other kinds of Enoki arrays (GPU
arrays, differentiable arrays).

Suppose that the interface of the ``Sensor`` class originally looks as follows:

.. code-block:: cpp

    class Sensor {
    public:
        /// Decode a measurement based on the sensor's response curve
        virtual float decode(float input) = 0;

        /// Return sensor's serial number
        virtual uint32_t serial_number() = 0;
    };

It is trivial to add a second method that takes vector inputs, like so:

.. code-block:: cpp
    :emphasize-lines: 9

    using FloatP = Packet<float, 8>;
    using MaskP  = mask_t<FloatP>;

    class Sensor {
    public:
        /// Scalar version
        virtual float decode(float input) = 0;

        /// Vector version
        virtual FloatP decode(FloatP input) = 0;

        /// Return sensor's serial number
        virtual uint32_t serial_number() = 0;
    };

This will work fine if there is just a single ``Sensor`` instance. But what if
there are many of them, e.g. when each ``FloatP`` array of measurements also
comes with a ``SensorP`` structure whose entries reference the sensor that
produced the measurement?

.. code-block:: cpp

    class Sensor;
    using SensorP = Array<Sensor *, 8>;

Ideally, we'd still be able to write the following code, but this sort of thing
is clearly not supported by standard C++.

.. code-block:: cpp

    SensorP sensor = ...;
    FloatP data = ...;

    data = sensor->decode(data);

Enoki provides a support layer that can handle such vectorized method calls. It
performs as many method calls as there are unique instances in the ``sensor``
array, and an optional mask is forwarded to the callee indicating the
associated active SIMD lanes. Null pointers in the ``data`` array are legal and
are considered as masked entries. The return value of masked entries is always
zero (or a zero-filled array/structure, depending on the method's return type).

The :c:macro:`ENOKI_CALL_SUPPORT_METHOD` macro is required to support the
above syntax. This generates the Enoki support layer that intercepts and
carries out the function call:

.. code-block:: cpp
    :emphasize-lines: 7, 13, 14, 15, 17

    class Sensor {
    public:
        // Scalar version
        virtual float decode(float input) = 0;

        // Vector version with optional mask argument
        virtual FloatP decode(FloatP input, MaskP mask) = 0;

        /// Return sensor's serial number
        virtual uint32_t serial_number() = 0;
    };

    ENOKI_CALL_SUPPORT_BEGIN(Sensor)
    ENOKI_CALL_SUPPORT_METHOD(decode)
    ENOKI_CALL_SUPPORT_METHOD(serial_number)
    /// .. potentially other methods ..
    ENOKI_CALL_SUPPORT_END(Sensor)

The macro supports functions taking an arbitrary number of arguments but
assumes that results are provided to the caller via the return value only
(i.e. no writing to arguments passed by reference). The mask, if present,
must be the last argument of the function.

Here is a hypothetical implementation of the ``Sensor`` interface:

.. code-block:: cpp

    class MySensor : Sensor {
    public:
        /// Vector version
        virtual FloatP decode(FloatP input, MaskP active) override {
            /// Keep track of invalid samples
            n_invalid += count(isnan(input) && active);

            /* Transform e.g. from log domain. */
            return log(input);
        }

        /// Return sensor's serial number
        uint32_t serial_number() { return 363436u; }

        // ...

        size_t n_invalid = 0;
    };

With this interface, the following vectorized expressions are now valid:

.. code-block:: cpp

    SensorP sensor = ...;
    FloatP data = ...;

    /* Unmasked version */
    data = sensor->decode(data);

    /* Masked version */
    auto mask = sensor->serial_number() > 1000;
    data = sensor->decode(data, mask);

Note how both functions with scalar and vector return values are vectorized
automatically.

The implementation of vector method calls depends on the array type and
hardware capabilities.

- On machines with the AVX512 instruction set, the ``vpextractq`` instruction
  is used to efficiently extract the unique set of instance pointers.

- The CUDA backend performs a parallel radix sort and run-length encoding of
  the pointer array using NVIDIA's `CUB library
  <https://nvlabs.github.io/cub/>`_ to obtain the list of unique pointers and
  indices referring to them. It then gathers the argument values corresponding
  to a particular pointer, evaluates the function, and then scatters the result
  into an output array.

- In all other cases, the unique elements are found using a linear sweep.

Supporting scalar *getter* functions
************************************

The above way of vectorizing a scalar *getter* function may involve multiple
virtual method calls, which is not particularly efficient when the invoked
function is very simple (e.g. a *getter*). Enoki provides an alternative macro
:c:macro:`ENOKI_CALL_SUPPORT_GETTER` that turns any such attribute lookup into
a *gather* operation. The macro takes the getter name and field name as
arguments. The macro :c:macro:`ENOKI_CALL_SUPPORT_FRIEND` is needed if the
field in question is a private member.

.. code-block:: cpp
    :emphasize-lines: 2, 14

    class Sensor {
        ENOKI_CALL_SUPPORT_FRIEND()
    public:
        /// ...

        /// Return sensor's serial number
        uint32_t serial_number() { return m_serial_number; }

    private:
        uint32_t m_serial_number;
    };

    ENOKI_CALL_SUPPORT_BEGIN(Sensor)
    ENOKI_CALL_SUPPORT_GETTER(serial_number, m_serial_number)
    ENOKI_CALL_SUPPORT_END(Sensor)

The usage is identical to before, i.e.:

.. code-block:: cpp

    using UInt32P = Packet<uint32_t, 8>;

    SensorP sensor = ...;
    UInt32P serial = sensor->serial_number();

Note that this trick even works for GPU arrays! In this case, the GPU will
directly fetch the value of the ``m_serial_number`` field from the CPU via
shared memory. However, this only works when the ``Sensor`` instance has been
allocated in *host-pinned* address space that will be reachable on the GPU. To
do so, add the :c:macro:`ENOKI_PINNED_OPERATOR_NEW` annotation that will
override the ``new`` and ``delete`` operator to ensure that this is always the
case for `Sensor` instances.

.. code-block:: cpp

    class Sensor {
        ENOKI_CALL_SUPPORT_FRIEND()
        ENOKI_PINNED_OPERATOR_NEW(UInt32P)
    public:
        // ...
    };

