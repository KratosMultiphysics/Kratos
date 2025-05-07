Building VexCL programs with CMake
==================================

In order to build a VexCL program with the `CMake`_ build system you need just
a couple of lines in your ``CmakeLists.txt``:

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.24.0)
    project(example)

    find_package(VexCL)

    add_executable(example example.cpp)
    target_link_libraries(example VexCL::OpenCL)

VexCL provides interface targets for the backends supported on the current
system. Possible choices are ``VexCL::OpenCL`` for the OpenCL backend,
``VexCL::Compute`` for Boost.Compute, ``VexCL::CUDA`` for CUDA, and
``VexCL::JIT`` for the just-in-time compiled OpenMP kernels.
The targets will take care of the appropriate compiler and linker flags for the
selected backend.

If you are interested in generating all possible backends, you can use
``vexcl_add_executables(example example.cpp)``, which will generate up to
four different versions of the same program, with ``_cl``, ``_comp``, ``_cuda``,
and ``_jit`` appended, depending on what backends were discovered. An interface
target is available for you to add dependencies to all targets at once:

.. code-block:: cmake

    vexcl_add_executables(example example.cpp)

    target_link_libraries(example INTERFACE MyDependenices)
    target_link_libraries(example_cl OpenCLOnlyDependency)


``find_package(VexCL)`` may be used when VexCL was installed system wide. If
that is not the case, you can just copy the VexCL into a subdirectory of your
project (or use git submodules) and replace the line with

.. code-block:: cmake

    add_subdirectory(vexcl)

.. _`CMake`: https://cmake.org/
