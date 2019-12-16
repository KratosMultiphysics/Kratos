## CoSimIO

_CoSimIO_ for exchanging data between different solvers.

The (**C++**) interface is defined in **`co_sim_io.h`**.
Interfaces for **C** (`c/co_sim_c_io.h`) and **Fortran** (`fortran/co_sim_io.f90`) also exist, which internally use the C++ interface.

### Usage

Examples for the usage of the _CoSimIO_ in solvers can be found in `CoSimulationApplication/custom_helpers/dummy_solver/`

**C++**
Include `co_sim_io.h` in your code and use it.

**C**
1. Include `c/co_sim_c_io.h` in your code and use it
2. Compile `c/co_sim_c_io.h` into a shared library (see `c/CMakeLists.txt`) and link against it

**Fortran**
1. Include `fortran/co_sim_io.f90` in your code and use it
2. Compile `c/co_sim_c_io.cpp` into a shared library (see `c/CMakeLists.txt`) and link against it. The `ISO_C_BINDING` module of fortran 2003 is used to connect the C-library with Fortran

### Notes on implementation:
- Dependency on Kratos: There is **NO** dependency on Kratos. The _CoSimIO_ can be used completely without including or linking against Kratos
- _CoSimIO_ is **header-only**, no linking is required
- The only dependency is **C++11**
- Different means of communication / data-exchange are available. Communication through files is the default, others (e.g. through Sockets or MPI) can be enabled at compile time (and selected at run time). This can introduce other dependencies such as boost or MPI. Except MPI, all these dependencies are header-only.