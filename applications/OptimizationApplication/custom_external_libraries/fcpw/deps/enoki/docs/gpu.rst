.. cpp:namespace:: enoki
.. _gpu:

GPU Arrays
==========

We will now switch gears to GPU arrays, whose operation is markedly different
from that of the other previously discussed Enoki array types.

The first major change when working with GPU arrays is that Enoki is no longer
a pure header file library. A compilation step becomes necessary, which
produces shared libraries against which applications must be linked. The reason
for this is to minimize template bloat: GPU arrays involve a certain amount of
additional machinery, and it would be wasteful to have to include the
underlying implementation in every file of an application that relies on GPU
arrays.

Compiling Enoki produces up to three libraries that are potentially of interest:

1. ``libenoki-cuda.so``: a just-in-time compiler that is used to realize the
   Enoki backend of the ``CUDAArray<T>`` array type discussed in this section.

2. ``libenoki-autodiff.so``: a library for maintaining a computation graph for
   automatic differentiation discussed in the next section.

3. ``enoki.cpyhon-37m-x86_64-linux-gnu.so`` (platform-dependent filename): a
   Python binding library that provides interoperability with Enoki's GPU
   arrays and differentiable GPU arrays.

Enter the following CMake command to compile all of them:

.. code-block:: bash

   cd <path-to-enoki>
   mkdir build
   cd build
   cmake -DENOKI_CUDA=ON -DENOKI_AUTODIFF=ON -DENOKI_PYTHON=ON ..
   make

For educational reasons, it is instructive to compile in Enoki in debug mode,
which enables a number of log messages that we will refer to in the remainder
of this section. Use the following CMake command to do so:

.. code-block:: bash

   cmake -DCMAKE_BUILD_TYPE=Debug -DENOKI_CUDA=ON -DENOKI_AUTODIFF=ON -DENOKI_PYTHON=ON ..

Using GPU Arrays in Python
--------------------------

We find it easiest to introduce Enoki's GPU arrays from within an interactive
Python interpreter and postpone the discussion of the C++ interface to the end
of this section. We'll start by importing the Enoki extension into an
interactive Python session and set the log level to a high value via
:cpp:func:`cuda_set_log_level`, which will be helpful in the subsequent
discussion.

.. code-block:: python

   >>> from enoki import *
   >>> cuda_set_log_level(4)

.. note::

    The first time that Enoki is imported on a new machine, it will trigger a
    kernel pre-compilation step that takes a few seconds.

The Enoki python bindings expose a number of types in the `enoki.cuda` submodule that
correspond to GPU-resident arrays. The following example initializes such an array with a
constant followed by a simple addition operation.

.. code-block:: python

   >>> from enoki.cuda import Float32 as FloatC
   >>> a = FloatC(1)
   cuda_trace_append(10): mov.$t1 $r1, 0f3f800000

   >>> a = a + a
   cuda_trace_append(11 <- 10, 10): add.rn.ftz.$t1 $r1, $r2, $r3

Observe the two ``cuda_trace_append`` log messages, which begin to reveal the
mechanics underlying the GPU backend: neither of these two operations has
actually occurred at this point. Instead, Enoki simply queued this computation
for later execution using an assembly-like intermediate language named `NVIDIA
PTX <https://docs.nvidia.com/cuda/parallel-thread-execution/index.html>`_.

For instance, the first line ``cuda_trace_append(10): mov.$t1 $r1, 0f3f800000``
indicates the creation of a new variable with index ``10``, and the associated
line of PTX will eventually be used to initialize this variable with a binary
constant representing the floating point value ``1.0``. The next
``cuda_trace_append`` command introduces a new variable ``11`` that will record
the result of the addition, while keeping track of the dependence on the
original variable ``10``, etc. More complex numerical operations (e.g. a
hyperbolic tangent) result in a longer sequence of steps that are similarly
enqueued:

.. code-block:: python

   >>> a = tanh(a)
   cuda_trace_append(12 <- 11): abs.ftz.$t1 $r1, $r2
   cuda_trace_append(13): mov.$t1 $r1, 0f3f200000
   ... 25 lines skipped ...
   cuda_trace_append(39 <- 38, 37): sub.rn.ftz.$t1 $r1, $r2, $r3
   cuda_trace_append(40 <- 39, 29, 14): selp.$t1 $r1, $r2, $r3, $r4

Eventually, numerical evaluation can no longer be postponed, e.g. when we try
to print the array contents:

.. code-block:: python

   >>> print(a)
   cuda_eval(): launching kernel (n=1, in=0, out=1, ops=31)
   .... many lines skipped ...
   cuda_jit_run(): cache miss, jit: 541 us, ptx compilation: 43.534, 10 registers
   [0.964028]

At this point, Enoki's JIT backend compiles and launches a kernel that contains
all of the computation queued thus far.

.. container:: toggle

   .. container:: header

      **Show/Hide the resulting PTX code**

   .. code-block:: bash

      .version 6.3
      .target sm_75
      .address_size 64

      .visible .entry enoki_8a163272(.param .u64 ptr,
                                     .param .u32 size) {
          .reg.b8 %b<41>;
          .reg.b16 %w<41>;
          .reg.b32 %r<41>;
          .reg.b64 %rd<41>;
          .reg.f32 %f<41>;
          .reg.f64 %d<41>;
          .reg.pred %p<41>;


          // Grid-stride loop setup
          ld.param.u64 %rd0, [ptr];
          ld.param.u32 %r1, [size];
          mov.u32 %r4, %tid.x;
          mov.u32 %r5, %ctaid.x;
          mov.u32 %r6, %ntid.x;
          mad.lo.u32 %r2, %r5, %r6, %r4;
          setp.ge.u32 %p0, %r2, %r1;
          @%p0 bra L0;

          mov.u32 %r7, %nctaid.x;
          mul.lo.u32 %r3, %r6, %r7;

      L1:
          // Loop body

          mov.f32 %f10, 0f3f800000;
          add.rn.ftz.f32 %f11, %f10, %f10;
          mul.rn.ftz.f32 %f12, %f11, %f11;
          mul.rn.ftz.f32 %f13, %f12, %f12;
          mul.rn.ftz.f32 %f14, %f13, %f13;
          mov.f32 %f15, 0fbbbaf0ea;
          mul.rn.ftz.f32 %f16, %f15, %f14;
          mov.f32 %f17, 0f3e088393;
          mov.f32 %f18, 0fbeaaaa99;
          fma.rn.ftz.f32 %f19, %f12, %f17, %f18;
          add.rn.ftz.f32 %f20, %f19, %f16;
          mov.f32 %f21, 0f3ca9134e;
          mov.f32 %f22, 0fbd5c1e2d;
          fma.rn.ftz.f32 %f23, %f12, %f21, %f22;
          fma.rn.ftz.f32 %f24, %f13, %f23, %f20;
          mul.rn.ftz.f32 %f25, %f12, %f11;
          fma.rn.ftz.f32 %f26, %f24, %f25, %f11;
          add.rn.ftz.f32 %f27, %f11, %f11;
          mov.f32 %f28, 0f3fb8aa3b;
          mul.rn.ftz.f32 %f29, %f28, %f27;
          ex2.approx.ftz.f32 %f30, %f29;
          mov.f32 %f31, 0f3f800000;
          add.rn.ftz.f32 %f32, %f30, %f31;
          rcp.approx.ftz.f32 %f33, %f32;
          add.rn.ftz.f32 %f34, %f33, %f33;
          mov.f32 %f35, 0f3f800000;
          sub.rn.ftz.f32 %f36, %f35, %f34;
          abs.ftz.f32 %f37, %f11;
          mov.f32 %f38, 0f3f200000;
          setp.ge.f32 %p39, %f37, %f38;
          selp.f32 %f40, %f36, %f26, %p39;

          // Store register %f40
          ldu.global.u64 %rd8, [%rd0 + 0];
          st.global.f32 [%rd8], %f40;

          add.u32     %r2, %r2, %r3;
          setp.ge.u32 %p0, %r2, %r1;
          @!%p0 bra L1;

      L0:
          ret;
      }

Internally, Enoki hands the PTX code over to CUDA's runtime compiler (`NVRTC
<https://docs.nvidia.com/cuda/nvrtc/index.html>`_), which performs a second
pass that translates from PTX to the native GPU instruction set *SASS*.

.. container:: toggle

    .. container:: header

        **Show/Hide the resulting SASS code**

    .. code-block:: bash

        enoki_8a163272:
            MOV R1, c[0x0][0x28];
            S2R R0, SR_TID.X;
            S2R R3, SR_CTAID.X;
            IMAD R0, R3, c[0x0][0x0], R0;
            ISETP.GE.U32.AND P0, PT, R0, c[0x0][0x168], PT;
        @P0 EXIT;
            BSSY B0, `(.L_2);
            ULDC.64 UR4, c[0x0][0x160];
        .L_3:
             LDG.E.64.SYS R2, [UR4];
             MOV R5, 0x3f76ca83;
             MOV R7, c[0x0][0x0];
             IMAD R0, R7, c[0x0][0xc], R0;
             ISETP.GE.U32.AND P0, PT, R0, c[0x0][0x168], PT;
             STG.E.SYS [R2], R5;
        @!P0 BRA `(.L_3);
             BSYNC B0;
        .L_2:
             EXIT ;
        .L_4:
             BRA `(.L_4);

This second phase is a full-fledged optimizing compiler with constant
propagation and common subexpression elimination. You can observe this in the
previous example because the second snippet is *much smaller*---in fact, almost
all of the computation was optimized away and replaced by a simple constant
(:math:`\tanh(2)\approx 0.964028`).

Enoki's approach is motivated by efficiency considerations: most array
operations are individually very simple and do not involve a sufficient amount
of computation to outweigh overheads related to memory accesses and GPU kernel
launches. Enoki therefore accumulates larger amounts of work (potentially
hundreds of thousands of individual operations) before creating and launching
an optimized GPU kernel. Once evaluated, array contents can be accessed without
triggering further computation:

.. code-block:: python

    >>> print(a)
    [0.964028]

Kernel caching
--------------

GPU kernel compilation consists of two steps: the first generates a PTX kernel
from the individual operations---this is essentially just string concatenation
and tends to be very fast (541 µs in the above example, most of which is caused
by printing assembly code onto the console due to the high log level).

The second step (``ptx compilation``) that converts the PTX intermediate
representation into concrete machine code that can be executed on the installed
graphics card is orders of magnitude slower (43 ms in the above example) but
only needs to happen once: whenever the same computation occurs again (e.g. in
subsequent iterations of an optimization algorithm), the previously generated
kernel is reused:

.. code-block:: python
    :emphasize-lines: 7

    >>> b = FloatC(1)
    >>> b = b + b
    >>> b = tanh(b)
    >>> print(b)
    cuda_eval(): launching kernel (n=1, in=0, out=1, ops=31)
    .... many lines skipped ...
    cuda_jit_run(): cache hit, jit backend: 550 us
    [0.964028]

A more complex example
----------------------

We now turn to a more complex example: computing the three-dimensional volume
of a sphere using Monte Carlo integration. To do so, we create a random number
generator RNG that will generate 1 million samples:

.. code-block:: python

    >>> from enoki.cuda import PCG32 as PCG32C, UInt64 as UInt64C
    >>> rng = PCG32C(UInt64C.arange(1000000))

Here, *PCG32* refers to a linear congruential generator from the section on
:ref:`random number generation <random>`. We use it to sample three random
number vectors from the RNG and create a dynamic array of 3D vectors
(``Vector3fC``).

.. code-block:: python

    >>> from enoki.cuda import Vector3f as Vector3fC
    >>> v = Vector3fC([rng.next_float32() * 2 - 1 for _ in range(3)])

Finally, we compute a mask that determines which of the uniformly distributed
vectors on the set :math:`[-1, 1]^3` lie within the unit sphere:

.. code-block:: python

    >>> inside = norm(v) < 1

At this point, seeding of the random number generator and subsequent sampling
steps touching its internal state have produced over a hundred different
operations generating various intermediate results along with the output
variable of interest.

To understand the specifics of this process, we assign a label to this variable
and enter the command :cpp:func:`cuda_whos`, which is analogous to ``whos`` in
IPython and MATLAB and generates a listing of all variables that are currently
registered (with the JIT compiler, in this case).

.. code-block:: python

    >>> set_label(inside, 'inside')
    >>> cuda_whos()

      ID        Type   E/I Refs   Size        Memory     Ready    Label
      =================================================================
      10        u32    0 / 1      1000000     3.8147 MiB  [ ]
      11        u64    0 / 1      1000000     7.6294 MiB  [ ]
      ... 126 lines skipped ...
      178       f32    0 / 1      1           4 B         [ ]
      179       msk    1 / 0      1000000     976.56 KiB  [ ]     inside
      =================================================================

      Memory usage (ready)     : 0 B
      Memory usage (scheduled) : 0 B + 20.027 MiB = 20.027 MiB
      Memory savings           : 350.95 MiB


The resulting output lists variables of many types (single precision floating
point values, 32/64 bit unsigned integers, masks, etc..), of which the last one
corresponds to the ``inside`` variable named above.

Note how each variable lists two *reference counts* (in the column ``E/I
refs``): the first (*external*) specifies how many times the variable is
referenced from an external application like the interactive Python prompt,
while the second (*internal*) counts how many times it is referenced as part of
queued arithmetic expressions. Variables with zero references in both categories
are automatically purged from the list.

Most of the variables are only referenced *internally*---these correspond to
temporaries created during a computation. Because they can no longer be
"reached" through external references, it would be impossible to ask the system
for the contents of such a temporary variable. Enoki relies on this observation
to perform an important optimization: rather than storing temporaries in
global GPU memory, their contents can be represented using cheap temporary GPU
registers. This yields significant storage and memory traffic savings: over 350
MiB of storage can be elided in the last example, leaving only roughly 20 MiB
of required storage.

In fact, these numbers can still change: we have not actually executed the
computation yet, and Enoki currently conservatively assumes that we plan to
continue using the random number generator ``rng`` and list of 3D vectors ``v``
later on. If we instruct Python to garbage-collect these two variables, the
required storage drops to less than a megabyte:

.. code-block:: python
   :emphasize-lines: 14

   >>> del v, rng
   >>> cuda_whos()

     ID        Type   E/I Refs   Size        Memory     Ready    Label
     =================================================================
     10        u32    0 / 1      1000000     3.8147 MiB  [ ]
     11        u64    0 / 1      1000000     7.6294 MiB  [ ]
     ... 126 lines skipped ...
     178       f32    0 / 1      1           4 B         [ ]
     179       msk    1 / 0      1000000     976.56 KiB  [ ]     inside
     =================================================================

     Memory usage (ready)     : 0 B
     Memory usage (scheduled) : 0 B + 976.56 KiB = 976.56 KiB
     Memory savings           : 324.25 MiB


Finally, we can "peek" into the ``inside`` array to compute the fraction of
points that lie within the sphere, which approximates the expected value
:math:`\frac{4}{3\cdot 2^3}\pi\approx0.523599`.

.. code-block:: python

   >>> count(inside) / len(inside)
   ... many lines skipped ...
   0.523946


Manually triggering JIT compilation
-----------------------------------

It is sometimes desirable to manually force Enoki's JIT compiler to generate a
kernel containing the computation queued thus far. For instance, rather than
compiling a long-running iterative algorithm into a single huge kernel, a
single kernel per iteration may be preferable. This can be accomplished by
explicitly invoking the :cpp:func:`cuda_eval` function periodically. An example:

.. code-block:: python

    >>> from enoki.cuda import UInt32 as UInt32C
    >>> a = UInt32C.arange(1234)

    >>> cuda_eval()
    cuda_eval(): launching kernel (n=1234, in=0, out=1, ops=1)

    >>> cuda_whos()

      ID        Type   E/I Refs   Size        Memory     Ready    Label
      =================================================================
      10        u32    1 / 0      1234        4.8203 KiB  [x]
      =================================================================

      Memory usage (ready)     : 4.8203 KiB
      Memory usage (scheduled) : 4.8203 KiB + 0 B = 4.8203 KiB
      Memory savings           : 0 B

The array is now marked "ready", which means that its contents were evaluated
and reside in GPU memory at an address that can be queried via the ``data``
field.

.. code-block:: python

    >>> a.data
    140427428626432


Actually, that is not entirely accurate: kernels are always launched
*asynchronously*, which means that the function :cpp:func:`cuda_eval` may have
returned before the GPU finished executing the kernel. Nonetheless, is
perfectly safe to begin using the variable immediately as asynchronous
communication with the GPU still observes a linear ordering guarantee.

In very rare cases (e.g. kernel benchmarking), it may be desirable to wait
until all currently running kernels have terminated. For this, invoke
:cpp:func:`cuda_sync` following :cpp:func:`cuda_eval`.

Parallelization and horizontal operations
-----------------------------------------

Recall the difference between :ref:`vertical <vertical>` and :ref:`horizontal
<horizontal>` operations: vertical operations are applied independently to each
element of a vector, while horizontal ones combine the different elements of a
vector. Enoki's GPU arrays are designed to operate very efficiently when
working with vertical operations that can be parallelized over the entire chip.

Horizontal operations (e.g. :cpp:func:`hsum`, :cpp:func:`all`,
:cpp:func:`count`, etc.) are best avoided whenever possible, because they
require that all prior computation has finished. In other words: each time
Enoki encounters a horizontal operation involving an unevaluated array, it
triggers a call to :cpp:func:`cuda_eval`. That said, horizontal reductions are
executed in parallel using NVIDIA's `CUB <https://nvlabs.github.io/cub/>`_
library, which is a highly performant implementation of these primitives.

Interfacing with NumPy
----------------------

Enoki GPU arrays support bidirectional conversion from/to NumPy arrays, which
will of course involve some communication between the CPU and GPU:

.. code-block:: python

   >>> x = FloatC.linspace(0, 1, 5)

   >>> # Enoki -> NumPy
   >>> y = Vector3fC(x, x*2, x*3).numpy()
   cuda_eval(): launching kernel (n=5, in=1, out=6, ops=36)

   >>> print(y)
   array([[0.  , 0.  , 0.  ],
          [0.25, 0.5 , 0.75],
          [0.5 , 1.  , 1.5 ],
          [0.75, 1.5 , 2.25],
          [1.  , 2.  , 3.  ]], dtype=float32)

   >>> # NumPy -> Enoki
   >>> Vector3fC(y)
   cuda_eval(): launching kernel (n=5, in=1, out=3, ops=27)
   [[0, 0, 0],
    [0.25, 0.5, 0.75],
    [0.5, 1, 1.5],
    [0.75, 1.5, 2.25],
    [1, 2, 3]]

Interfacing with PyTorch
------------------------

`PyTorch <https://pytorch.org/>`_ GPU tensors are supported as well. In this
case, copying occurs on the GPU (but is still necessary, as the two frameworks
use different memory layouts for tensors).

.. code-block:: python

   >>> x = FloatC.linspace(0, 1, 5)

   >>> # Enoki -> PyTorch
   >>> y = Vector3fC(x, x*2, x*3).torch()
   cuda_eval(): launching kernel (n=5, in=2, out=5, ops=31)

   >>> y
   tensor([[0.0000, 0.0000, 0.0000],
           [0.2500, 0.5000, 0.7500],
           [0.5000, 1.0000, 1.5000],
           [0.7500, 1.5000, 2.2500],
           [1.0000, 2.0000, 3.0000]], device='cuda:0')

   >>> # PyTorch -> Enoki
   >>> Vector3fC(y)
   cuda_eval(): launching kernel (n=5, in=1, out=3, ops=27)
   [[0, 0, 0],
    [0.25, 0.5, 0.75],
    [0.5, 1, 1.5],
    [0.75, 1.5, 2.25],
    [1, 2, 3]]

Note how the ``.numpy()`` and ``.torch()`` function calls both triggered a
mandatory kernel launch to ensure that that the array contents were ready
before returning a representation in the other framework. This can be wasteful
when converting many variables at an interface between two frameworks. For this
reason, both ``.numpy()`` and ``.torch()`` functions take an optional ``eval``
argument that is set to ``True`` by default. Passing ``False`` causes the
operation to return an uninitialized NumPy or PyTorch array, while at the same
time scheduling Enoki code that will eventually fill this memory with valid
contents the next time that :cpp:func:`cuda_eval` is triggered. An example is
shown below. This feature is to be used with caution.

.. code-block:: python

   >>> x = FloatC.linspace(0, 1, 5)

   >>> y = Vector3fC(x, x*2, x*3).numpy(False)

   >>> y
   array([[0., 0., 0.],
          [0., 0., 0.],
          [0., 0., 0.],
          [0., 0., 0.],
          [0., 0., 0.]], dtype=float32)

   >>> cuda_eval()
   cuda_eval(): launching kernel (n=5, in=1, out=4, ops=36)

   >>> y
   array([[0.  , 0.  , 0.  ],
          [0.25, 0.5 , 0.75],
          [0.5 , 1.  , 1.5 ],
          [0.75, 1.5 , 2.25],
          [1.  , 2.  , 3.  ]], dtype=float32)

Scatter/gather operations
-------------------------

The GPU backend also supports scatter and gather operations
involving GPU arrays as targets.

.. code-block:: python

    >>> a = FloatC.zero(10)
    >>> b = UInt32C.arange(5)
    >>> scatter(target=a, source=FloatC(b), index=b*2)
    >>> a
    cuda_eval(): launching kernel (n=5, in=1, out=2, ops=9)
    [0, 0, 1, 0, 2, 0, 3, 0, 4, 0]

Note that gathering from an unevaluated Enoki array is not guaranteed to be a
vertical operation, hence it triggers a call to :cpp:func:`cuda_eval`.

Caching memory allocations
--------------------------

Similar to the `PyTorch memory allocator
<https://pytorch.org/docs/stable/notes/cuda.html#cuda-memory-management>`_,
Enoki uses a caching scheme to avoid very costly device synchronizations when
releasing memory. This means that freeing a large GPU variable doesn't cause
the associated memory region to become available for use by the operating
system or other frameworks like Tensorflow or PyTorch. Use the function
:cpp:func:`cuda_malloc_trim` to fully purge all unused memory. The function is
only relevant when working with other frameworks and does not need to be called
to free up memory for use by Enoki itself.

C++ interface
-------------

Everything demonstrated in the above sections can be directly applied to
C++ programs as well. To use the associated (dynamic array) type :cpp:class:`CUDAArray`,
include the header

.. code-block:: cpp

    #include <enoki/cuda.h>

Furthermore, applications must be linked against the ``cuda`` and
``enoki-cuda`` libraries. The following snippet contains a C++ translation of
the Monte Carlo integration Python example shown earlier.

.. code-block:: cpp

    #include <enoki/cuda.h>
    #include <enoki/random.h>

    using namespace enoki;

    using FloatC    = CUDAArray<float>;
    using Vector3fC = Array<FloatC, 3>;
    using PCG32C    = PCG32<FloatC>;
    using MaskC     = mask_t<FloatC>;

    int main(int argc, char **argv) {
        PCG32C rng(PCG32_DEFAULT_STATE, arange<FloatC>(1000000));

        Vector3fC v(
            rng.next_float32() * 2.f - 1.f,
            rng.next_float32() * 2.f - 1.f,
            rng.next_float32() * 2.f - 1.f
        );

        MaskC inside = norm(v) < 1.f;

        std::cout << count(inside) / (float) inside.size() << std::endl;
    }

.. _horizontal_ops_on_gpu:

Suggestions regarding horizontal operations
-------------------------------------------

When vectorizing code, we may sometimes want to skip an expensive computation
when it is not actually needed by any elements in the array being processed.
This is usually done with the :cpp:func:`any` function and yields good
performance in when targeting the *CPU* (e.g. with the AVX512 backend). An
example:

.. code-block:: cpp

    auto condition = variable > 1.f;
    if (any(condition))
        result[condition] = /* expensive-to-evaluate expression */;

However, recall the discussion earlier in this section, which explained how
horizontal operations tend to be fairly expensive in conjunction with the GPU
backend because they flush the JIT compiler. This effectively breaks up the
program into smaller kernels, increasing memory traffic and missing potential
optimization opportunities. Arrays processed by the GPU backend tend to be much
larger, and from a probabilistic viewpoint it is often likely that the
:cpp:func:`any` function call will in any case evaluate to ``true``. For these
reasons, skipping test and always evaluating the expression often leads to
better performance on the GPU.

Enoki provides alternative horizontal reductions of masks named
:cpp:func:`any_or`, :cpp:func:`all_or`, :cpp:func:`none_or` that do exactly
this: they skip evaluation when compiling for GPU targets and simply return the
supplied template argument. For other targets, they behave as usual. With this
change, the example looks as follows:

.. code-block:: cpp

    auto condition = variable > 1.f;
    if (any_or<true>(condition))
        result[condition] = /* expensive-to-evaluate expression */;


Differences between Enoki and existing frameworks
-------------------------------------------------
Enoki was designed as a numerical foundation for differentiable physical
simulations, specifically the `Mitsuba renderer
<https://github.com/mitsuba-renderer/mitsuba2>`_, though it is significantly
more general and should be a trusty tool for a variety of simulation and
optimization problems.

Its GPU and Autodiff backends are related to well-known frameworks like
`TensorFlow <https://www.tensorflow.org/>`_ and `PyTorch
<https://pytorch.org/>`_ that have become standard tools for training and
evaluating neural networks. In the following, we outline the main differences
between these frameworks and Enoki.

Both PyTorch and Tensorflow provide two main operational modes: *eager mode*
directly evaluates arithmetic operations on the GPU, which yields excellent
performance in conjunction with arithmetically intensive operations like
convolutions and large matrix-vector multiplications, both of which are
building blocks of neural networks. When evaluating typical simulation code
that mainly consists of much simpler arithmetic (e.g. additions,
multiplications, etc.), the resulting memory traffic and scheduling overheads
induce severe bottlenecks. An early prototype of Enoki provided a
``TorchArray<T>`` type that carried out operations using PyTorch's eager mode,
and the low performance of this combination eventually motivated us to develop
the technique based on JIT compilation introduced in the previous section.

The second operational mode requires an up-front specification of the complete
computation graph to generate a single optimized GPU kernel (e.g. via XLA in
TensorFlow and ``jit.trace`` in PyTorch). This is feasible for neural networks,
whose graph specification is very regular and typically only consists of a few
hundred operations. Simulation code, on the other hand, involves much larger
graphs, whose structure is *unpredictable*: program execution often involves
randomness, which could cause jumps to almost any part of the system. The full
computation graph would simply be the entire codebase (potentially on the order
of hundreds of thousands lines of code), which is of course far too big.

Enoki's approach could be interpreted as a middle ground between the two
extremes discussed above. Graphs are created on the fly during a simulation,
and can be several orders of magnitude larger compared to typical neural
networks. They consist mostly of unstructured and comparably simple arithmetic
that is lazily fused into optimized CUDA kernels. Since our system works
without an up-front specification of the full computation graph, it must
support features like dynamic indirection via virtual function calls that can
simultaneously branch to multiple different implementations. The details of
this are described in the section on :ref:`function calls <calls>`.

Note that that there are of of course many use cases where PyTorch, Tensorflow,
etc. are vastly superior to Enoki, and it is often a good idea to combine the
two in such cases (e.g. to feed the output of a differentiable simulation into
a neural network).

One last related framework is `ArrayFire
<https://github.com/arrayfire/arrayfire>`_, which provides a JIT compiler that
lazily fuses instructions similar to our ``CUDAArray<T>`` type. ArrayFire
targets a higher-level language (C), but appears to be limited to fairly small
kernels (100 operations by default), and does not support a mechanism for
automatic differentiation. In contrast, Enoki emits an intermediate
representation (PTX) and fuses instructions into comparatively larger kernels
that often exceed 100K instructions.
