Test MORSE CMake modules
========================

Procedure to test the `find_package` modules hosted in `modules/find`

```
git clone https://gitlab.inria.fr/solverstack/morse_cmake.git
cd morse_cmake/modules/find/tests
mkdir build
cd build
# to test all the packages, e.g.
cmake .. -DENABLE_CTEST=ON
# to test some only, e.g.
cmake .. -DPACKAGES="HWLOC;STARPU;PASTIX"
# to enable some components in some find, e.g. here we look for QUARK linked with HWLOC and PASTIX with PARSEC and STARPU
cmake .. -DENABLE_CTEST=ON -DQUARK_COMPONENTS="HWLOC" -DPASTIX_COMPONENTS="PARSEC;STARPU"
```
