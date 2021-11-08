Initialization
==============

Selecting backend
-----------------

VexCL provides the following backends:

- **OpenCL**, built on top of `Khronos C++ API`_. The backend is selected when
  :c:macro:`VEXCL_BACKEND_OPENCL` macro is defined, or by default. Link with
  ``libOpenCL.so`` on unix-like systems or with ``OpenCL.dll`` on Windows.
- **Boost.Compute**. The backend is also based on OpenCL, but uses core
  functionality of the Boost.Compute_ library instead of somewhat outdated
  Khronos C++ API. The additional advantage is the increased interoperability
  between VexCL and Boost.Compute_. The backend is selected when
  :c:macro:`VEXCL_BACKEND_COMPUTE` macro is defined. Link with
  ``libOpenCL.so``/``OpenCL.dll`` and make sure that Boost.Compute_ headers are
  in the include path.
- **CUDA**, uses the NVIDIA CUDA technology. The backend is selected when
  :c:macro:`VEXCL_BACKEND_CUDA` macro is defined. Link with
  ``libcuda.so``/``cuda.dll``. For the CUDA backend to work, CUDA Toolkit has
  to be installed, and NVIDIA CUDA compiler driver nvcc has to be in executable
  PATH and usable at runtime.

Whatever backend is selected, you will need to link to Boost.System_ and
Boost.Filesystem_ libraries. Some systems may also require linking to
Boost.Thread_ and Boost.Date_Time_. All of those are distributed with
Boost_ libraries collection.

.. _`Khronos C++ API`: https://www.khronos.org/registry/cl
.. _Boost.Compute: https://github.com/boostorg/compute
.. _Boost.System: http://www.boost.org/doc/libs/release/libs/system/doc/index.html
.. _Boost.Filesystem: http://www.boost.org/doc/libs/release/libs/filesystem/doc/index.html
.. _Boost.Thread: http://www.boost.org/doc/libs/release/libs/thread/doc/index.html
.. _Boost.Date_Time: http://www.boost.org/doc/libs/release/libs/date_time/doc/index.html
.. _Boost: http://www.boost.org

Context initialization
----------------------

VexCL transparently works with multiple compute devices that are present in the
system. A VexCL context is initialized with a device filter, which is just a
functor that takes a const reference to a :cpp:class:`vex::backend::device`
instance and returns a boolean value. Several standard filters are provided
(see below), and one can easily add a custom functor. Filters may be combined
with logical operators. All compute devices that satisfy the provided filter
are added to the created context. In the example below all GPU devices that
support double precision arithmetic are selected:

.. code-block:: cpp

    #include <iostream>
    #include <stdexcept>
    #include <vexcl/vexcl.hpp>

    int main() {
        vex::Context ctx( vex::Filter::GPU && vex::Filter::DoublePrecision );

        if (!ctx) throw std::runtime_error("No devices available.");

        // Print out list of selected devices:
        std::cout << ctx << std::endl;
    }

One of the most convenient filters is :cpp:class:`vex::Filter::Env` which
selects compute devices based on environment variables. It allows to switch
the compute device without the need to recompile the program.

Each stateful object in VexCL, like :cpp:class:`vex::vector\<T>`, takes an STL
vector of :cpp:class:`vex::backend::command_queue` instances. The
:cpp:class:`vex::Context` class is just a convenient way to initialize and hold
the command queues. Since it provides the corresponding type conversion
operator, it also may be used directly for object initialization:

.. code-block:: cpp

    vex::vector<double> x(ctx, n);

But the users are not required to actually create a :cpp:class:`vex::Context`
instance. They may just use the command queues initialized elsewhere. In the
following example the Boost.Compute_ is used as a backend and takes care of
initializing the OpenCL context:

.. code-block:: cpp

    #include <iostream>
    #include <boost/compute.hpp>

    #define VEXCL_BACKEND_COMPUTE
    #include <vexcl/vexcl.hpp>

    int main() {
        boost::compute::command_queue bcq = boost::compute::system::default_queue();

        // Use Boost.Compute queue to allocate VexCL vectors:
        vex::vector<int> x({bcq}, 16);
    }

Device filters
--------------

Common filters
^^^^^^^^^^^^^^

These filters are supported for all backends:

- :cpp:class:`vex::Filter::Any`. Selects all available devices.
- :cpp:class:`vex::Filter::DoublePrecision`. Selects devices that support
  double precision arithmetics.
- :cpp:class:`vex::Filter::Count(n)`. Selects first ``n`` devices that are
  passed through the filter.  This filter should be the last in the filter
  chain. This will assure that it will be applied only to devices which passed
  all other filters.  Otherwise, you could get less devices than planned (every
  time this filter is applied, its internal counter is decremented).
- :cpp:class:`vex::Filter::Position(n)`. Selects single device at the given
  position.
- :cpp:class:`vex::Filter::Env`. Selects devices with respect to environment
  variables. Recognized variables are:

  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_DEVICE`      | Name of the device or its substring.         |
  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_MAX_DEVICES` | Maximum number of devices to select. The     |
  |                            | effect is similar to the                     |
  |                            | ``vex::Filter::Count`` filter above.         |
  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_POSITION`    | Single device with the specified position    |
  |                            | in the list of available devices. The effect |
  |                            | is similar to the ``vex::Filter::Position``  |
  |                            | filter above.                                |
  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_PLATFORM`    | OpenCL platform name or its substring.       |
  |                            | Only supported for OpenCL-based backends.    |
  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_VENDOR`      | OpenCL device vendor name or its substring.  |
  |                            | Only supported for OpenCL-based backends.    |
  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_TYPE`        | OpenCL device type. Possible values are      |
  |                            | ``CPU``, ``GPU``, ``ACCELERATOR``.           |
  |                            | Only supported for OpenCL-based backends.    |
  +----------------------------+----------------------------------------------+
  | :c:macro:`OCL_EXTENSION`   | OpenCL device supporting the specified       |
  |                            | extension.                                   |
  |                            | Only supported for OpenCL-based backends.    |
  +----------------------------+----------------------------------------------+

- :cpp:class:`vex::Filter::Exclusive(filter)`. This is a filter wrapper that
  allows to obtain exclusive access to compute devices. This may be helpful if
  several compute devices are present in the system and several processes are
  trying to grab a single device. The exclusivity is only guaranteed between
  processes that use the ``Exclusive`` filter wrapper.

OpenCL-specific filters
^^^^^^^^^^^^^^^^^^^^^^^

These filters are only available for OpenCL and Boost.Compute_ backends:

- :cpp:class:`vex::Filter::CLVersion(major,minor)`. Selects devices that
  support the specified version of OpenCL standard.
- :cpp:class:`vex::Filter::Extension(string)`. Selects devices that provide the
  specified extension.
- :cpp:class:`vex::Filter::GLSharing`. Selects devices that support OpenGL
  sharing extension.
  This is a shortcut for ``vex::Filter::Extension("cl_khr_gl_sharing")``.
- :cpp:class:`vex::Filter::Type(cl_device_type)`. Selects devices with the
  specified device type. The device type is a bit mask.
- :cpp:class:`vex::Filter::GPU`. Selects GPU devices.
  This is a shortcut for ``vex::Filter::Type(CL_DEVICE_TYPE_GPU)``.
- :cpp:class:`vex::Filter::CPU`. Selects CPU devices.
  This is a shortcut for ``vex::Filter::Type(CL_DEVICE_TYPE_CPU)``.
- :cpp:class:`vex::Filter::Accelerator`. Selects Accelerator devices.
  This is a shortcut for ``vex::Filter::Type(CL_DEVICE_TYPE_ACCELERATOR)``.

Custom filters
^^^^^^^^^^^^^^

In case more complex functionality is required than provided by the builtin
filters, the users may introduce their own functors:

.. code-block:: cpp

    // Select a GPU with more than 4GiB of global memory:
    vex::Context ctx(vex::Filter::GPU &&
                     [](const vex::backend::device &d) {
                         size_t GiB = 1024 * 1024 * 1024;
                         return d.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() >= 4 * GiB;
                     });

Reference
---------

.. doxygenclass:: vex::Context
    :members:

.. cpp:function:: std::vector<vex::backend::device> vex::backend::device_list<DevFilter>(DevFilter &&filter)

    Returns vector of compute devices satisfying the given criteria without
    trying to initialize the contexts on the devices.


