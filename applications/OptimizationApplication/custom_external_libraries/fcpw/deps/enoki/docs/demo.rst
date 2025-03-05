.. _demo:

A quick demonstration
=====================

This section contains a quick and somewhat contrived demonstration that
illustrates a number of basic Enoki features. The remainder of the
documentation contains a more systematic introduction of this functionality.

Consider the function

.. math::

   f(x)=\begin{cases}
      12.92x, &x \le 0.0031308\\
      1.055x^{1/2.4} -0.055, &\mathrm{otherwise,}
   \end{cases}


which converts a linear color value into one that can be displayed on
a modern display following the `sRGB <https://en.wikipedia.org/wiki/SRGB>`_
standard---this is known as gamma correction. A standard C++ implementation
of this function might look as follows:

.. code-block:: cpp

   float srgb_gamma(float x) {
      if (x <= 0.0031308f)
         return x * 12.92f;
      else
         return 1.055f * std::pow(x, 1.f / 2.4f) - 0.055f;
   }

A Enoki implementation of the same computation first replaces ``float`` by a
generic ``Value`` type. The second difference is that scalar conditionals are
replaced by generalized expressions involving masks.

.. code-block:: cpp

   template <typename Value> Value srgb_gamma(Value x) {
      return enoki::select(
         x <= 0.0031308f,
         x * 12.92f,
         1.055f * enoki::pow(x, 1.f / 2.4f) - 0.055f
      );
   }

Vectorization
-------------

The simple generalization from normal C++ code to an Enoki function template
enables a number of interesting applications. For instance, the function
automatically extends to cases where ``Value`` is a color type with three
components, in which case the arithmetic operations recursively thread through
the array.

.. code-block:: cpp

   using Color3f = enoki::Array<float, 3>;

   Color3f input = /* ... */;
   Color3f output = srgb_gamma(input);

Arrays can be nested arbitrarily: the following snippet declares a 16-wide
``FloatP`` "packet" type (hence the "P" suffix) and uses it to construct a
a new type storing 16 colors that will all be processed in parallel.

.. code-block:: cpp

   using FloatP   = enoki::Array<float, 16>;
   using Color3fP = enoki::Array<FloatP, 3>;

   Color3fP input = /* ... */;
   Color3fP output = srgb_gamma(input);

If this code is compiled on a machine supporting the SSE4.2, AVX, AVX2, or
AVX512 instructions set extensions, these vector instructions will be leveraged
to carry out the computation more efficiently.


Execution on the GPU
--------------------

Vectorization is not restricted to the CPU---for instance, the following type
declarations create a special array that is resident in GPU memory. In this mode
of operation, Enoki relies on an internal just-in-time compiler to generate
efficient CUDA kernels on the fly.

.. code-block:: cpp

   using FloatC   = enoki::CUDAArray<float>;
   using Color3fC = enoki::Array<FloatC, 3>;

   Color3fC input = /* ... */;
   Color3fC output = srgb_gamma(input);

Enoki's ``CUDAArray<T>`` type applies an important optimization that leads to
significantly improved performance: in contrast to the previous examples, the
function call ``srgb_gamma(input)`` now merely records the sequence of
computations that is needed to determine the value of ``output`` but does not
yet execute it.

Eventually, this evaluation can no longer be postponed (e.g. when we try to
access or print the array contents). At this point, Enoki's JIT backend
compiles and executes a kernel that contains all queued computations using
NVIDIA's PTX intermediate representation. All of this happens automatically: in
particular, no CUDA-specific rewrite (e.g. to ``nvcc`` compatible kernels) of
the program is necessary!

Automatic differentiation
-------------------------

Enoki can also apply transparent forward or reverse-mode automatic
differentiation to a program using a special ``enoki::DiffArray<T>`` array that
wraps a number type or another Enoki array ``T``.

For instance, the following example computes the gradient of a loss function
that measures L2 distance from a given gamma-corrected color value. Both primal
and gradient-related computations involve GPU-resident arrays, and the
resulting computation is queued up as in the previously example using Enoki's
just-in-time compiler.

.. code-block:: cpp

   using FloatC   = enoki::CUDAArray<float>;
   using FloatD   = enoki::DiffArray<FloatC>;
   using Color3fD = enoki::Array<FloatD, 3>;

   Color3fD input = /* ... */;
   enoki::set_requires_gradient(input);

   Color3fD output = srgb_gamma(input);

   FloatD loss = enoki::norm(output - Color3fD(.1f, .2f, .3f));
   enoki::backward(loss);

   std::cout << enoki::gradient(input) << std::endl;

The scalar case
---------------

All Enoki functions also accept non-array arguments, hence the original scalar
implementation remains available:

.. code-block:: cpp

   float input = /* ... */;
   float output = srgb_gamma(input);

Python bindings
---------------

Modern C++ systems often strive to provide fine-grained Python bindings to
facilitate rapid prototyping and interoperability with other software. Enoki is
designed to work with the widely used `pybind11
<https://github.com/pybind/pybind11>`_ library (itself based on template
metaprogramming) to facilitate this. Exposing an Enoki function on the Python
side is usually a 1-liner, even for the "fancy" GPU+autodiff variants, as in the
following example:

.. code-block:: cpp

   /// Create python bindings with 2 overloads (here, 'm' is a py::module)
   m.def("srgb_gamma", &srgb_gamma<float>);
   m.def("srgb_gamma", &srgb_gamma<Color3fD>);


Summary
-------

In summary: Enoki, along with a generalized template implementation of a
computation enables several powerful transformations:

1. A simple type substitution yields an equivalent vectorized computation that
   leverages vector instructions on modern processor architectures.

2. Symbolic execution of the computation using a a just-in-time compiler
   yields efficient kernels that run on NVIDIA GPUs.

3. Further type transformations enable tracking of derivatives through
   a calculation, either on the CPU or the GPU.

4. The above transformations can all be deduced from the type of the resulting
   functions. This is an ideal fit for metaprogramming-based libraries like
   `pybind11 <https://github.com/pybind/pybind11>`_ which inspect the
   type of a function to generate high-quality binding code.

There are many missing pieces that weren't discussed in this basic overview:
how to handle more complex control flow, types and data structures, virtual
method calls, and so on. The remainder of this documentation provides a more
systematic overview of these topics.
