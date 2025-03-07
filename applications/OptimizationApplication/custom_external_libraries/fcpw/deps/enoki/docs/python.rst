.. cpp:namespace:: enoki

Python bindings
===============

Enoki provides support for `pybind11 <https://github.com/pybind/pybind11/>`_, a
lightweight header-only binding library that is used to expose C++ types in
Python and vice versa.

To use this feature, include the following header file in the extension module:

.. code-block:: cpp

    #include <enoki/python.h>

Usage
-----

The example below shows how to create bindings for a simple vector computation
that converts spherical to Cartesian coordinates. A CMake build system file is
provided at the :ref:`bottom <py-build>` of this page.

Extension module
****************

.. code-block:: cpp

    #include <enoki/python.h>

    /* Import pybind11 and Enoki namespaces */
    namespace py = pybind11;
    using namespace enoki;
    using namespace py::literals; // Enables the ""_a parameter tags used below

    /* Define a packet type used for vectorization */
    using FloatP    = Packet<float>;

    /* Dynamic type for arbitrary-length arrays */
    using FloatX    = DynamicArray<FloatP>;

    /* Various flavors of 3D vectors */
    using Vector3f  = Array<float, 3>;
    using Vector3fP = Array<FloatP, 3>;
    using Vector3fX = Array<FloatX, 3>;

    /* The function we want to expose in Python */
    template <typename Float>
    Array<Float, 3> sph_to_cartesian(Float theta, Float phi) {
        auto sc_theta = sincos(theta);
        auto sc_phi   = sincos(phi);

        return {
            sc_theta.first * sc_phi.second,
            sc_theta.first * sc_phi.first,
            sc_theta.second
        };
    }

    /* The function below is called when the extension module is loaded. It performs a
       sequence of m.def(...) calls which define functions in the module namespace 'm' */
    PYBIND11_MODULE(pybind11_test /* <- name of extension module */, m) {
        m.doc() = "Enoki & pybind11 test plugin"; // Set a docstring

        /* 1. Bind the scalar version of the function */
        m.def(
              /* Name of the function in the Python extension module */
              "sph_to_cartesian",

              /* Function that should be exposed */
              sph_to_cartesian<float>,

              /* Function docstring */
              "Convert from spherical to cartesian coordinates [scalar version]",

              /* Parameter names for function signature in docstring */
              "theta"_a, "phi"_a
        );

        /* 2. Bind the packet version of the function */
        m.def("sph_to_cartesian",
               /* The only differnce is the FloatP template argument */
               sph_to_cartesian<FloatP>,
              "Convert from spherical to cartesian coordinates [packet version]",
              "theta"_a, "phi"_a);

        /* 3. Bind dynamic version of the function */
        m.def("sph_to_cartesian",
               /* Note the use of 'vectorize_wrapper', which is described below */
               vectorize_wrapper(sph_to_cartesian<FloatP>),
              "Convert from spherical to cartesian coordinates [dynamic version]",
              "theta"_a, "phi"_a);
    }

pybind11 infers the necessary binding code from the type of the function
provided to the ``def()`` calls. Including the :file:`enoki/python.h` header is
all it takes to make the pybind11 library fully Enoki-aware---arbitrarily
nested dynamic and static arrays will be converted automatically.

In practice, one would usually skip the packet version since it is subsumed by
the dynamic case.

Using the extension from Python
*******************************

The following iteractive session shows how to load the extension module and
query its automatically generated help page.

.. code-block:: pycon

    Python 3.5.2 |Anaconda 4.2.0 (x86_64)| (default, Jul  2 2016, 17:52:12)
    [GCC 4.2.1 Compatible Apple LLVM 4.2 (clang-425.0.28)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.

    >>> import pybind11_test
    >>> help(pybind11_test)

    Help on module pybind11_test

    NAME
        pybind11_test - Enoki & pybind11 test plugin

    FUNCTIONS
        sph_to_cartesian(...)
            sph_to_cartesian(*args, **kwargs)
            Overloaded function.

            1. sph_to_cartesian(theta: float, phi: float)
                   -> numpy.ndarray[dtype=float32, shape=(3)]

            Convert from spherical to cartesian coordinates [scalar version]

            2. sph_to_cartesian(theta: numpy.ndarray[dtype=float32, shape=(8)],
                                phi: numpy.ndarray[dtype=float32, shape=(8)])
                   -> numpy.ndarray[dtype=float32, shape=(8, 3)]

            Convert from spherical to cartesian coordinates [packet version]

            3. sph_to_cartesian(theta: numpy.ndarray[dtype=float32, shape=(n)],
                                phi: numpy.ndarray[dtype=float32, shape=(n)])
                   -> numpy.ndarray[dtype=float32, shape=(n, 3)]

            Convert from spherical to cartesian coordinates [dynamic version]

    FILE
        /Users/wjakob/pybind11_test/pybind11_test.cpython-35m-darwin.so

As can be seen, the help describes all three overloads along with the name and shape of their input arguments.
Let's try calling one of them:

.. code-block:: python

    >>> from pybind11_test import sph_to_cartesian
    >>> sph_to_cartesian(theta=1, phi=2)
    array([-0.35017547,  0.76514739,  0.54030228], dtype=float32)

Note how the returned Enoki array was automatically converted into a NumPy array.

Let's now call the dynamic version of the function. We will use ``np.linspace``
to generate inputs, which actually have an *incorrect* ``dtype`` of
``np.float64``. The binding layer detects this and automatically creates a
temporary single precision input array before performing the function call.

.. code-block:: python

    >>> import numpy as np
    >>> sph_to_cartesian(theta=np.linspace(0.0, 1.0, 10),
    ...                  phi=np.linspace(1.0, 2.0, 10))
    array([[ 0.        ,  0.        ,  1.        ],
           [ 0.04919485,  0.09937215,  0.99383354],
           [ 0.07527862,  0.20714317,  0.9754101 ],
           [ 0.07696848,  0.31801295,  0.9449569 ],
           [ 0.05418137,  0.42652887,  0.90284967],
           [ 0.00803789,  0.52735412,  0.84960753],
           [-0.05919253,  0.61553025,  0.7858873 ],
           [-0.14420365,  0.68672061,  0.71247464],
           [-0.24281444,  0.73742425,  0.63027501],
           [-0.35017547,  0.76514739,  0.54030228]], dtype=float32)

.. _py-build:

Build system
************

The following ``CMakeLists.txt`` file can be used to build the module on
various platforms.

.. code-block:: cmake

    cmake_minimum_required (VERSION 2.8.12)
    project(pybind11_test CXX)
    include(CheckCXXCompilerFlag)

    # Set a default build configuration (Release)
    if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
      message(STATUS "Setting build type to 'Release' as none was specified.")
      set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo")
    endif()

    # Enable C++14 support
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|Intel")
      CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)
      if (HAS_CPP14_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
      else()
        message(FATAL_ERROR "Unsupported compiler -- C++14 support is needed!")
      endif()
    endif()

    # Assumes that pybind11 is located in the 'pybind11' subdirectory
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)

    # Assumes that enoki is located in the 'enoki' subdirectory
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/enoki)

    # Enable some helpful vectorization-related compiler flags
    enoki_set_compile_flags()
    enoki_set_native_flags()

    include_directories(enoki/include pybind11/include)

    # Compile our pybind11 module
    pybind11_add_module(pybind11_test pybind11_test.cpp)

Reference
---------

Please refer to pybind11's extensive `documentation
<http://pybind11.readthedocs.io/en/master/?badge=master>`_. for details on
using it in general. The :file:`enoki/python.h` API only provides one public
function:

.. cpp:function:: template <typename Func> auto vectorize_wrapper(Func func)

    "Converts" a function that takes a set of packets and structures of packets
    as inputs into a new function that processes dynamic versions of these
    parameters. Non-array arguments are not transformed. For instance, it would
    turn the following hypothetical signature

    .. code-block:: cpp

        FloatP my_func(Array<FloatP, 3> position, GPSRecord2<FloatP> record, int scalar);

    into

    .. code-block:: cpp

        FloatX my_func(Array<FloatX, 3> position, GPSRecord2<FloatX> record, int scalar);

    where

    .. code-block:: cpp

        using FloatX = DynamicArray<FloatP>;

    This is handy because a one-liner like
    ``vectorize_wrapper(sph_to_cartesian<FloatP>)`` in the above example is all
    it takes to take a packet version of a function and expose a dynamic
    version that can process arbitrarily large NumPy arrays.
