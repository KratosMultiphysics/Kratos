# EigenSolversApplication

The *EigenSolversApplication* is a thin wrapper for the [Eigen linear algebra library](http://eigen.tuxfamily.org/index.php?title=Main_Page).

It provides the following direct sparse solvers:

| solver_type          | Matrix kind | Dependencies |
|----------------------|:-----------:|:------------:|
| `eigen_sparse_lu`    | Square      | None         |
| `eigen_pardiso_llt`  | SPD         | Intel® MKL   |
| `eigen_pardiso_ldlt` | SPD         | Intel® MKL   |
| `eigen_pardiso_lu`   | Square      | Intel® MKL   |

SPD = symmetric positive definite

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
        "solver_type" : "Eigen_SparseLU"
    }
    ```

5. Run the simulation

## Enable MKL (optional)

In case you have installed [MKL](https://software.intel.com/en-us/mkl), you can also use the Pardiso solvers.

1. Run the MKL setup script

    **Windows:** before cmake in `configure.bat`

    ```batch
    call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" intel64 lp64
    ```
    
    **Linux:** before cmake in `configure.sh`

    ```batch
    source "~/intel/mkl/bin/mklvars.sh" intel64 lp64
    ```

2. Build Kratos

3. Copy the required MKL library to the Kratos `lib`
    
    **Windows:**

    ```
    C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64_win\mkl\mkl_rt.dll
    ```

    **Linux:**

    ```
    ~/intel/mkl/lib/intel64/libmkl_rt.so
    ```

    or add the folder to your `PATH`/`LD_LIBRARY_PATH` variable.
