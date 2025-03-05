.. cpp:namespace:: enoki

Quaternions
===========

Enoki provides a vectorizable type for quaternion arithmetic.
To use it, include the following header:

.. code-block:: cpp

    #include <enoki/quaternion.h>

Usage
-----

The following example shows how to define and perform basic arithmetic using
:cpp:class:`enoki::Quaternion` vectorized over 4-wide packets.

.. code-block:: cpp

    /* Declare underlying packet type, could just be 'float' for scalar arithmetic */
    using FloatP = Packet<float, 4>;

    /* Define vectorized quaternion type */
    using QuaternionP = Quaternion<FloatP>;

    QuaternionP a = QuaternionP(1.f, 0.f, 0.f, 0.f);
    QuaternionP b = QuaternionP(0.f, 1.f, 0.f, 0.f);

    /* Compute several rotations that interpolate between 'a' and 'b' */
    FloatP t = linspace<FloatP>(0.f, 1.f);
    QuaternionP c = slerp(a, b, t);

    std::cout << "Interpolated quaternions:" << std::endl;
    std::cout << c << std::endl << std::endl;

    /* Turn into a 4x4 homogeneous coordinate rotation matrix packet */
    using MatrixP = Matrix<FloatP, 4>;
    MatrixP c_rot = quat_to_matrix<Matrix4P>(c);

    std::cout << "Rotation matrices:" << std::endl;
    std::cout << c_rot << std::endl << std::endl;

    /* Round trip: turn the rotation matrices back into rotation quaternions */
    QuaternionP c2 = matrix_to_quat(c_rot);

    if (hsum(abs(c-c2)) < 1e-6f)
        std::cout << "Test passed." << std::endl;
    else
        std::cout << "Test failed." << std::endl;

    /* Prints:

        Interpolated quaternions:
        [0 + 1i + 0j + 0k,
         0 + 0.866025i + 0.5j + 0k,
         0 + 0.5i + 0.866025j + 0k,
         0 - 4.37114e-08i + 1j + 0k]

        Rotation matrices:
        [[[1, 0, 0, 0],
          [0, -1, 0, 0],
          [0, 0, -1, 0],
          [0, 0, 0, 1]],
         [[0.5, 0.866025, 0, 0],
          [0.866025, -0.5, 0, 0],
          [0, 0, -1, 0],
          [0, 0, 0, 1]],
         [[-0.5, 0.866025, 0, 0],
          [0.866025, 0.5, 0, 0],
          [0, 0, -1, 0],
          [0, 0, 0, 1]],
         [[-1, -8.74228e-08, 0, 0],
          [-8.74228e-08, 1, 0, 0],
          [-0, 0, -1, 0],
          [0, 0, 0, 1]]]

         Test passed.
    */

Reference
---------

.. cpp:class:: template <typename Type> Quaternion : StaticArrayImpl<Type, 4>

    The class :cpp:class:`enoki::Quaternion` is a 4D Enoki array whose
    components are of type ``Type``. Various arithmetic operators (e.g.
    multiplication) and transcendental functions are overloaded so that they
    provide the correct behavior for quaternion-valued inputs.

    .. cpp:function:: Quaternion(Type x, Type y, Type z, Type w)

        Initialize a new :cpp:class:`enoki::Quaternion` instance with the value
        :math:`x\mathbf{i} + y\mathbf{j} + z\mathbf{k} + w`

        .. warning::

            Note the different order convention compared to
            :cpp:func:`Complex::Complex`.

    .. cpp:function:: Quaternion(Array<Type, 3> imag, Type real)

        Creates a :cpp:class:`enoki::Quaternion` instance from the given
        imaginary and real inputs.

        .. warning::

            Note the different order convention compared to
            :cpp:func:`Complex::Complex`.

    .. cpp:function:: Quaternion(Type f)

        Creates a real-valued :cpp:class:`enoki::Quaternion` instance from ``f``.
        This constructor effectively changes the broadcasting behavior of
        non-quaternion inputs---for instance, the snippet

        .. code-block:: cpp

            auto value_a = zero<Array<float, 4>>();
            auto value_q = zero<Quaternion<float>>();

            value_a += 1.f; value_q += 1.f;

            std::cout << "value_a = "<< value_a << ", value_q = " << value_q << std::endl;

        prints ``value_a = [1, 1, 1, 1], value_q = 1 + 0i + 0j + 0k``, which is
        the desired behavior for quaternions. For standard Enoki arrays, the
        number ``1.f`` is broadcast to all four components.

Elementary operations
*********************

.. cpp:function:: template <typename Quat> Quat identity()

    Returns the identity quaternion.

.. cpp:function:: template <typename T> T real(Quaternion<T> q)

    Extracts the real part of ``q``.

.. cpp:function:: template <typename T> Array<T, 3> imag(Quaternion<T> q)

    Extracts the imaginary part of ``q``.

.. cpp:function:: template <typename T> T abs(Quaternion<T> q)

    Compute the absolute value of ``q``.

.. cpp:function:: template <typename T> T sqrt(Quaternion<T> q)

    Compute the square root of ``q``.

.. cpp:function:: template <typename T> Quaternion<T> conj(Quaternion<T> q)

    Evaluates the quaternion conjugate of ``q``.

.. cpp:function:: template <typename T> Quaternion<T> rcp(Quaternion<T> q)

    Evaluates the quaternion reciprocal of ``q``.

Arithmetic operators
********************

Only a few arithmetic operators need to be overridden to support quaternion
arithmetic. The rest are automatically provided by Enoki's existing operators
and broadcasting rules.

.. cpp:function:: template <typename T> Quaternion<T> operator*(Quaternion<T> q0, Quaternion<T> q1)

    Evaluates the quaternion product of ``q1`` and ``z2``.

.. cpp:function:: template <typename T> Quaternion<T> operator/(Quaternion<T> q0, Quaternion<T> q1)

    Evaluates the quaternion division of ``q1`` and ``z2``.

Stream operators
****************

.. cpp:function:: std::ostream& operator<<(std::ostream &os, const Quaternion<T> &z)

    Sends the quaternion number ``q`` to the stream ``os`` using the format
    ``1 + 2i + 3j + 4k``.

Exponential, logarithm, and power function
******************************************

.. cpp:function:: template <typename T> Quaternion<T> exp(Quaternion<T> q)

    Evaluates the quaternion exponential of ``q``.

.. cpp:function:: template <typename T> Quaternion<T> log(Quaternion<T> q)

    Evaluates the quaternion logarithm of ``q``.

.. cpp:function:: template <typename T> Quaternion<T> pow(Quaternion<T> q0, Quaternion<T> q1)

    Evaluates the quaternion power of ``q0`` raised to the ``q1``.

Operations for rotation-related computations
********************************************

.. cpp:function:: template <typename T, typename Float> Quaternion<T> slerp(Quaternion<T> q0, Quaternion<T> q1, Float t)

    Performs a spherical linear interpolation between two rotation quaternions,
    where ``slerp(q0, q1, 0.f) == q0`` and ``slerp(q0, q1, 1.f) == q1``.

.. cpp:function:: template <typename Matrix, typename T> Matrix quat_to_matrix(Quaternion<T> q)

    Converts a rotation quaternion into a :math:`3\times 3` or :math:`4\times4`
    homogeneous coordinate transformation matrix (depending on the
    ``Matrix`` template argument).

.. cpp:function:: template <typename T, size_t Size> Quaternion<T> matrix_to_quat(MatrixP<T, Size> q)

    Converts a :math:`3\times 3` or :math:`4\times 4` homogeneous containing a
    pure rotation into a rotation quaternion.

.. cpp:function:: template <typename Quat, typename Vector3, typename Float> Quat rotate(Vector3 v, Float angle)

    Constructs a rotation quaternion, which rotates by ``angle`` radians
    around the axis ``v``. The function requires ``v`` to be normalized.
