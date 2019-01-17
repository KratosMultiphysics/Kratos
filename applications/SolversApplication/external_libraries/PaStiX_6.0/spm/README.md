# SPM: SParse Matrix package

[![pipeline status](https://gitlab.inria.fr/solverstack/spm/badges/master/pipeline.svg)](https://gitlab.inria.fr/solverstack/spm/pipelines) [![coverage report](https://gitlab.inria.fr/solverstack/spm/badges/master/coverage.svg)](https://sonarqube.bordeaux.inria.fr/sonarqube/dashboard?id=hiepacs%3Aspm%3Agitlab%3Amaster)

SPM (SParse Matrix package) is a scientific library that provides
basic operation coverage to manipulate sparse matrices in CSC, CSR,
and IJV format.
The functionalities covered are:

* sparse matrix -by- dense matrix products
* sparse matrix -by- vector products
* norm computations
* matrix and vector scaling
* sort functions
* graph symmetrization and merge of duplicate entries
* random generators for right hand sides
* check routines for linear solvers
* In-place format conversion routines
* drivers to read matrices from files (MatrixMarket, Harwell-Boeing/RSA, IJV, ...)
* Laplacian generators for stencils
* multi-dof and variadic dof as input

This package is for now available in sequential for shared memory
system, and will be further developed to handle distributed matrices
over MPI processes.

Python and Fortran90 wrappers are included in the package.

## Get SPM

To use last development state of SPM, please clone the master
branch. Note that SPM contains a `git submodule` **morse_cmake**.
To get sources please use these commands:

    # if git version >= 1.9
      git clone --recursive git@gitlab.inria.fr:solverstack/spm.git
      cd spm
    # else
      git clone git@gitlab.inria.fr:solverstack/spm.git
      cd spm
      git submodule init
      git submodule update

## Documentation

The documentation will be soon available as Doxygen pages.

## Installation


### Build and install with CMake

SPM can be built using [CMake](https://cmake.org/). This
installation requires to have some library dependencies already
installed on the system:

* BLAS (MKL, OpenBlas, ...) and CBLAS (sequential version required)
* LAPACK and LAPACKE
* Python and Fortran compiler for wrappers and examples (optional)

For instance, on debian-like systems, dependencies can be installed with the following command:

      sudo apt-get install cmake gcc gfortran libopenblas-dev liblapacke-dev python-numpy

The main options to configure the SPM configuration build are:

* Classic cmake options:
  * CMAKE_BUILD_TYPE: Debug, RelWithDebInfo, Release, MinSizeRel; we recommend to use the Release, or RelWithDebInfo, for performance.
  * CMAKE_INSTALL_PREFIX: Specify the prefix directory to install the library
  * BUILD_SHARED_LIBS=[OFF]: Enable the shared libraries build. This option needs to be enabled for the Python wrapper.
* Integer type:
  * SPM_INT64[=ON]: Enable/disable int64_t for integer arrays.
* Documentation:
  * BUILD_DOCUMENTATION[=OFF] to enable the Doxygen documentation generation

## Get involved!

### Reporting an issue

We strongly recommend all users to use the issue tracker to report any
problems with the software, or for any feature request. We will try
our best to answer them in a short time frame.

### Contributions

https://gitlab.inria.fr/solverstack/spm/blob/master/CONTRIBUTING.md

### Authors

The following people contribute or contributed to the development of SPM:

* Mathieu Faverge
* Matthieu Kuhn
* Xavier Lacoste
* Grégoire Pichon
* Florent Pruvost
* Pierre Ramet
* Theophile Terraz

If we forgot your name, please let us know that we can fix that mistake.

### Licence

<https://gitlab.inria.fr/solverstack/spm/blob/master/LICENCE>
