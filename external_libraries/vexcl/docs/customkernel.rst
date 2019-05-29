Custom kernels
==============

As `Kozma Prutkov`_ repeatedly said, "One cannot embrace the unembraceable". So
in order to be usable, VexCL has to support custom kernels.
``vex::backend::kernel`` is a thin wrapper around a compute kernel for each of
the contexts. Its constructor takes a command queue and the kernel source code,
and its function call operator submits the kernel to the specified command
queue.  The following example builds and launches a custom kernel for the a
context with a single device:

.. code-block:: cpp

    // Compile the kernel. This can be done once per program lifetime.
    // If offline kernel cache is enabled, it will be used for custom kernels as well.
    vex::backend::kernel dummy(ctx.queue(0), VEX_STRINGIZE_SOURCE(
        kernel void dummy(ulong n, global int *x) {
            for(size_t i = get_global_id(0); i < n; i += get_global_size(0))
                x[i] = 42;
        }),
        "dummy");

    vex::vector<int> x(ctx, n);

    // Apply the kernel to the vector partition located on the first device:
    dummy(ctx.queue(0), static_cast<cl_ulong>(n), x(0));

In case there are several devices in the context, you will need to create an
instance of the kernel for each of the devices.
:cpp:func:`vex::vector::operator()()` returns vector partition located on the
given device.  If the result depends on the neighboring points, one has to keep
in mind that these points are possibly located on a different compute device.
In this case the exchange of these halo points has to be addressed manually.

.. code-block:: cpp

    std::vector<vex::backend::kernel> kernel;

    // Compile and store the kernels for the later use.
    for(uint d = 0; d < ctx.size(); d++) {
        kernel.emplace_back(ctx.queue(d), VEX_STRINGIZE_SOURCE(
            kernel void dummy(ulong n, global float *x) {
                for(size_t i = get_global_id(0); i < n; i += get_global_size(0))
                    x[i] = 4.2;
            }),
            "dummy");
    }

    // Apply the kernels to the vector partitions on each device.
    for(uint d = 0; d < ctx.size(); d++)
        kernel[d](ctx.queue(d), static_cast<cl_ulong>(x.part_size(d)), x(d));

.. _`Kozma Prutkov`: http://en.wikipedia.org/wiki/Kozma_Prutkov
