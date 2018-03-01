# EigenSolversApplication

The *EigenSolversApplication* is a thin wrapper for the [Eigen linear algebra library](http://eigen.tuxfamily.org/index.php?title=Main_Page).

## Direct sparse solvers

The application provides the following direct sparse solvers:

| Python class      | solver_type          | Matrix kind | Dependencies |
|-------------------|----------------------|:-----------:|:------------:|
| SparseLUSolver    | `eigen_sparse_lu`    | Square      | None         |
| SparseQRSolver    | *not available*      | Rectangular | None         |
| PardisoLLTSolver  | `eigen_pardiso_llt`  | SPD         | Intel® MKL   |
| PardisoLDLTSolver | `eigen_pardiso_ldlt` | SPD         | Intel® MKL   |
| PardisoLUSolver   | `eigen_pardiso_lu`   | Square      | Intel® MKL   |

SPD = symmetric positive definite

**Example**:

```json
{
    "solver_type": "eigen_sparse_lu"
}
```

## Generalized eigensystem solver

The application provides a generalized eigensystem solver for sparse matrices. It gives the eigenvalues and eigenvectors for the smallest eigenvalues. MKL routines are used automatically if they are available.

**Example:**

```json
{
    "solver_type": "eigen_eigensystem",
    "number_of_eigenvalues": 1,
    "max_iteration": 1000,
    "tolerance": 1e-6,
    "echo_level": 1
}
```

## Build instructions

1. Download and unpack [Eigen](http://eigen.tuxfamily.org/)

2. Set the required definitions for cmake

    **Windows:** in `configure.bat`

    ```batch
    -DEIGEN_SOLVERS_APPLICATION=ON ^
    -DEIGEN_ROOT="<path to eigen>" ^
    ```

    **Linux:** in `configure.sh`

    ```bash
    -DEIGEN_SOLVERS_APPLICATION=ON \
    -DEIGEN_ROOT="<path to eigen>" \
    ```

    > **Hint:** The `EIGEN_ROOT` directory should contain a file called `README.md`.

3. Build Kratos

4. Setup the `ProjectParameters.json`

    ```json
    "linear_solver_settings": {
        "solver_type" : "eigen_sparse_lu"
    }
    ```

5. Run the simulation

## Enable MKL (optional)

In case you have installed [MKL](https://software.intel.com/en-us/mkl), you can also use the Pardiso solvers.

1. Run the MKL setup script before building Kratos:

    **Windows:**

    ```batch
    call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" intel64 lp64
    ```

    **Linux:**

    ```batch
    source ~/intel/mkl/bin/mklvars.sh intel64 lp64
    ```

2. Add the following flag to CMake to your configure script:

    ```batch
    -DUSE_EIGEN_MKL=ON ^
    ```

3. Build Kratos

4. Usage:

    **Windows:**

    Copy the required MKL libraries to the Kratos `lib`

    ```batch
    C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\mkl\mkl_core.dll
    C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\mkl\mkl_rt.dll
    C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\mkl\mkl_intel_thread.dll
    C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\mkl\mkl_def.dll
    C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\compiler\libiomp5md.dll
    ```

    or add the folders to your `PATH`/`LD_LIBRARY_PATH` variable.

    **Linux:**

    Set the environment before using MKL
    ```batch
    source ~/intel/mkl/bin/mklvars.sh intel64 lp64
    ```