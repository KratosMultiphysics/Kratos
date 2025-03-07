.. cpp:namespace:: enoki

Standard Template Library
=========================

When Enoki extracts packets or slices through custom data structures, it also
handles STL data structures including ``std::array``,  and ``std::pair``,  and
``std::tuple``. Please review the section on :ref:`dynamic arrays <dynamic>`
for general details on vectorizing over dynamic arrays and working with slices.

To use this feature, include the following header file:

.. code-block:: cpp

    #include <enoki/stl.h>

Usage
-----

Consider the following example, where a function returns a ``std::tuple``
containing a 3D position and a mask specifying whether the computation was
successful. When the :file:`enoki/stl.h` header file is included, Enoki's
dynamic vectorization machinery can be applied to vectorize such functions over
arbitrarily large inputs.

.. code-block:: cpp
    :emphasize-lines: 2,3,4,5,6,30,36

    /// Return value of 'my_function'
    template <typename T>
    using Return = std::tuple<
        Array<T, 3>,
        mask_t<T>
    >;

    template <typename T> Return<T> my_function(T theta, T phi) {
        /* Turn spherical -> cartesian coordinates */
        Array<T, 3> pos(
            sin(theta) * cos(phi),
            sin(theta) * sin(phi),
            cos(theta)
        );

        /* Only points on the top hemisphere are 'valid' */
        return std::make_pair(pos, pos.z() > 0);
    }

    /// Packet of floats
    using FloatP  = Packet<float>;

    /// Arbitrarily large sequence of floats
    using FloatX  = DynamicArray<FloatP>;

    /// Tuple containing a packet of results
    using ReturnP = Return<FloatP>;

    /// Tuple containing dynamic arrays with arbitrarily many results
    using ReturnX = Return<FloatX>;

    int main(int argc, char *argv[]) {
        FloatX theta = linspace<FloatX>(-10.f, 10.f, 10);
        FloatX phi = linspace<FloatX>(0.f, 60.f, 10);

        ReturnX result = vectorize(my_function<FloatP>, theta, phi);

        /* Prints:
            [[0.544021, 0, -0.839072],
             [-0.924676, -0.373065, 0.0761302],
             [0.478888, 0.461548, 0.746753],
             [0.0777672, 0.173978, -0.981674],
             [-0.0330365, -0.895583, 0.443666],
             [-0.304446, 0.842896, 0.443666],
             [0.127097, -0.141995, -0.981674],
             [0.596782, -0.293616, 0.746753],
             [-0.994388, 0.0734624, 0.0761293],
             [0.518133, 0.165823, -0.839072]]
        */

        std::cout << std::get<0>(result) << std::endl;

        /* Prints:
            [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
        */
        std::cout << std::get<1>(result) << std::endl;
    }
