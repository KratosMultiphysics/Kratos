## External Solvers Application

External Solvers application is the place where interface to third party solvers are located.

In the following list you can find which solvers are currently supported by kratos and links to the library projects or algortihms in which are bases:

### Direct Solvers
* __SuperLUSolver__: Based on [SuperLu](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) project.
* __PastixSolver__: Based on [Pastix](http://pastix.gforge.inria.fr/files/README-txt.html) project.

### Iterative Solvers:
* __SuperLUIterativeSolver__: Based on [SuperLu](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) project.
* __GMRESSolver__: Based in the [GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method) method.


## Installing Blas and Lapack
Blas and Lapack are required to compile the `ExternalSolversApplication`.
After installing they have to be added to the `configure` file using the CMake Flags `BLAS_LIBRARIES` and `LAPACK_LIBRARIES`.
Typical installation paths are:

### Windows
~~~
-DLAPACK_LIBRARIES="C:\CompiledLibs\blas_x64\liblapack.lib"
-DBLAS_LIBRARIES="C:\CompiledLibs\blas_x64\libblas.lib"
~~~

### MacOS
~~~
-DLAPACK_LIBRARIES="/usr/lib/liblapack.dylib"
-DBLAS_LIBRARIES="/usr/lib/libblas.dylib"
~~~
