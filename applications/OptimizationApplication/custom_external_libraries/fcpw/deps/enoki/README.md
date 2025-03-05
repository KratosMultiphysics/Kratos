<p align="center"><img src="https://github.com/mitsuba-renderer/enoki/raw/master/docs/enoki-logo.png" alt="Enoki logo" width="300"/></p>

# Enoki — structured vectorization and differentiation on modern processor architectures

| Documentation   | Linux             | Windows             |
|      :---:      |       :---:       |        :---:        |
| [![docs][1]][2] | [![rgl-ci][3]][4] | [![appveyor][5]][6] |


[1]: https://readthedocs.org/projects/enoki/badge/?version=master
[2]: http://enoki.readthedocs.org/en/master
[3]: https://rgl-ci.epfl.ch/app/rest/builds/buildType(id:Enoki_Build)/statusIcon.svg
[4]: https://rgl-ci.epfl.ch/viewType.html?buildTypeId=Enoki_Build&guest=1
[5]: https://ci.appveyor.com/api/projects/status/68db7e5es7el1btd/branch/master?svg=true
[6]: https://ci.appveyor.com/project/wjakob/enoki/branch/master

## Introduction

**Enoki** is a C++ template library that enables automatic transformations of
numerical code, for instance to create a "wide" vectorized variant of an
algorithm that runs on the CPU or GPU, or to compute gradients via transparent
forward/reverse-mode automatic differentation.

The core parts of the library are implemented as a set of header files with no
dependencies other than a sufficiently C++17-capable compiler (GCC >= 8.2,
Clang >= 7.0, Visual Studio >= 2017). Enoki code reduces to efficient SIMD
instructions available on modern CPUs and GPUs—in particular, Enoki supports:

* **Intel**: AVX512, AVX2, AVX, and SSE4.2,
* **ARM**: NEON/VFPV4 on armv7-a, Advanced SIMD on 64-bit armv8-a,
* **NVIDIA**: CUDA via a *Parallel Thread Execution* (PTX) just-in-time compiler.
* **Fallback**: a scalar fallback mode ensures that programs still run even
  if none of the above are available.

Deploying a program on top of Enoki usually serves three goals:

1. Enoki ships with a convenient library of special functions and data
   structures that facilitate implementation of numerical code (vectors,
   matrices, complex numbers, quaternions, etc.).

2. Programs built using these can be instantiated as *wide* versions that
   process many arguments at once (either on the CPU or the GPU).

   Enoki is also *structured* in the sense that it handles complex programs
   with custom data structures, lambda functions, loadable modules, virtual
   method calls, and many other modern C++ features.

3. If derivatives are desired (e.g. for stochastic gradient descent), Enoki
   performs transparent forward or reverse-mode automatic differentiation of
   the entire program.

Finally, Enoki can do all of the above simultaneously: if desired, it can
compile the same source code to multiple different implementations (e.g.
scalar, AVX512, and CUDA+autodiff).

### Motivation

The development of this library was prompted by the author's frustration
with the current vectorization landscape:

1. Auto-vectorization in state-of-the-art compilers is inherently local. A
   computation whose call graph spans separate compilation units (e.g. multiple
   shared libraries) simply can't be vectorized.

2. Data structures must be converted into a *Structure of Arrays* (SoA) layout
   to be eligible for vectorization.

   <p align="center">
       <img src="https://github.com/mitsuba-renderer/enoki/raw/master/docs/intro-01.png" alt="SoA layout" width="400"/>
   </p>

   This is analogous to performing a matrix transpose of an application's
   entire memory layout—an intrusive change that is likely to touch almost
   every line of code.

3. Parts of the application likely have to be rewritten using [intrinsic
   instructions](https://software.intel.com/sites/landingpage/IntrinsicsGuide),
   which is going to look something like this:

   <p align="center">
       <img src="https://github.com/mitsuba-renderer/enoki/raw/master/docs/intro-02.png" alt="intrinsics" width="400"/>
   </p>


   Intrinsics-heavy code is challenging to read and modify once written, and it
   is inherently non-portable. CUDA provides a nice language environment
   for programming GPUs but does nothing to help with the other requirements
   (vectorization on CPUs, automatic differentiation).

4. Vectorized transcendental functions (*exp*, *cos*, *erf*, ..) are not widely
   available. Intel, AMD, and CUDA provide proprietary implementations, but many
   compilers don't include them by default.

5. It is desirable to retain both scalar and vector versions of an algorithm,
   but ensuring their consistency throughout the development cycle becomes a
   maintenance nightmare.

6. *Domain-specific languages* (DSLs) for vectorization such as
   [ISPC](https://ispc.github.io) address many of the above issues but assume
   that the main computation underlying an application can be condensed into a
   compact kernel that is implementable using the limited language subset of
   the DSL (e.g. plain C in the case of ISPC).

   This is not the case for complex applications, where the "kernel" may be
   spread out over many separate modules involving high-level language features
   such as functional or object-oriented programming.

### What Enoki does differently

Enoki addresses these issues and provides a *complete* solution for vectorizing
and differentiating modern C++ applications with nontrivial control flow and
data structures, dynamic memory allocation, virtual method calls, and vector
calls across module boundaries. It has the following design goals:

1. **Unobtrusive**. Only minor modifications are necessary to convert existing
   C++ code into its Enoki-vectorized equivalent, which remains readable and
   maintainable.

2. **No code duplication**. It is generally desirable to provide both scalar
   and vectorized versions of an API, e.g. for debugging, and to preserve
   compatibility with legacy code. Enoki code extensively relies on class and
   function templates to achieve this goal without any code duplication—the
   same code template can be leveraged to create scalar, CPU SIMD, and GPU
   implementations, and each variant can provide gradients via automatic
   differentiation if desired.

3. **Custom data structures**. Enoki can also vectorize custom data
   structures. All the hard work (e.g. conversion to SoA format) is handled by
   the C++17 type system.

4. **Function calls**. Vectorized calls to functions in other compilation units
   (e.g. a dynamically loaded plugin) are possible. Enoki can even vectorize
   method or virtual method calls (e.g. ``instance->my_function(arg1, arg2,
   ...);`` when ``instance`` turns out to be an array containing many different
   instances).

5. **Mathematical library**. Enoki includes an extensive mathematical support
   library with complex numbers, matrices, quaternions, and related operations
   (determinants, matrix, inversion, etc.). A set of transcendental and special
   functions supports real, complex, and quaternion-valued arguments in single
   and double-precision using polynomial or rational polynomial approximations,
   generally with an average error of <1/2 ULP on their full domain.
   These include exponentials, logarithms, and trigonometric and hyperbolic
   functions, as well as their inverses. Enoki also provides real-valued
   versions of error function variants, Bessel functions, the Gamma function,
   and various elliptic integrals.

   <p align="center">
       <img src="https://github.com/mitsuba-renderer/enoki/raw/master/docs/intro-03.png" alt="Transcendentals" width="720"/>
   </p>

   Importantly, all of this functionality is realized using the abstractions of
   Enoki, which means that it transparently composes with vectorization,
   the JIT compiler for generating CUDA kernels, automatic differentiation, etc.

6. **Portability**. When creating vectorized CPU code, Enoki supports arbitrary
   array sizes that don't necessarily match what is supported by the underlying
   hardware (e.g. 16 x single precision on a machine, whose SSE vector only has
   hardware support for 4 x single precision operands). The library uses
   template metaprogramming techniques to efficiently map array expressions
   onto the available hardware resources. This greatly simplifies development
   because it's enough to write a single implementation of a numerical
   algorithm that can then be deployed on any target architecture. There are
   non-vectorized fallbacks for everything, thus programs will run even on
   unsupported architectures (albeit without the performance benefits of
   vectorization).

7. **Modular architecture**. Enoki is split into two major components: the
   front-end provides various high-level array operations, while the back-end
   provides the basic ingredients that are needed to realize these operations
   using the SIMD instruction set(s) supported by the target architecture.

   The CPU vector back-ends e.g. make heavy use of SIMD intrinsics to
   ensure that compilers generate efficient machine code. The
   intrinsics are contained in separate back-end header files (e.g.
   ``array_avx.h`` for AVX intrinsics), which provide rudimentary
   arithmetic and bit-level operations. Fancier operations (e.g.
   *atan2*) use the back-ends as an abstract interface to the hardware,
   which means that it's simple to support other instruction sets such
   as a hypothetical future AVX1024 or even an entirely different
   architecture (e.g. a DSP chip) by just adding a new back-end.

8. **License**. Enoki is available under a non-viral open source license
   (3-clause BSD).

## Cloning

Enoki depends on two other repositories
([pybind11](https://github.com/pybind/pybind11) and
[cub](https://nvlabs.github.io/cub)) that are required when using certain
optional features, specifically differentiable GPU arrays with Python bindings.

To fetch the entire project including these dependencies, clone the project
using the ``--recursive`` flag as follows:

```bash
$ git clone --recursive https://github.com/mitsuba-renderer/enoki
```

## Documentation

An extensive set of tutorials and reference documentation are available at
[readthedocs.org](http://enoki.readthedocs.org/en/master).

## About

This project was created by [Wenzel Jakob](http://rgl.epfl.ch/people/wjakob).
It is named after [Enokitake](https://en.wikipedia.org/wiki/Enokitake), a type
of mushroom with many long and parallel stalks reminiscent of data flow in
vectorized arithmetic.

Enoki is the numerical foundation of version 2 of the [Mitsuba
renderer](https://github.com/mitsuba-renderer/mitsuba2), though it is
significantly more general and should be a trusty tool for a variety of
simulation and optimization problems.

When using Enoki in academic projects, please cite

```bibtex
@misc{Enoki,
   author = {Wenzel Jakob},
   year = {2019},
   note = {https://github.com/mitsuba-renderer/enoki},
   title = {Enoki: structured vectorization and differentiation on modern processor architectures}
}
```
