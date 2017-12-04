# EigenSolversApplication

The *EigenSolversApplication* is a thin wrapper for the [Eigen linear algebra library](http://eigen.tuxfamily.org/index.php?title=Main_Page).

It provides the following direct sparse solvers:

| solver_type         | Matrix kind | Dependencies |
|---------------------|:-----------:|:------------:|
| `Eigen_SparseLU`    | Square      | None         |
| `Eigen_PardisoLLT`  | SPD         | Intel® MKL   |
| `Eigen_PardisoLDLT` | SPD         | Intel® MKL   |
| `Eigen_PardisoLU`   | Square      | Intel® MKL   |

SPD = symmetric positive definite

## Windows instructions

1. Download an unpack [Eigen](http://bitbucket.org/eigen/eigen/get/3.3.4.zip)

2. Add the following definitions to your Kratos `configure.bat`:

    ```batch
    -DEIGEN_SOLVERS_APPLICATION=ON ^
    -DEIGEN_ROOT="<path to eigen>" ^
    ```

    > **Hint:** The `EIGEN_ROOT` directory should contain a file called `README.md`.

3. Build Kratos

4. Setup your `ProjectParameters.json`

    ```json
    "linear_solver_settings": {
        "solver_type" : "Eigen_SparseLU"
    }
    ```

5. Run your simulation

## Enable MKL (optional)

In case your have installed [MKL](https://software.intel.com/en-us/mkl), you can also use the Pardiso solvers by running the MKL setup script at the beginning of your `configure.bat`.

For example:

```batch
call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" intel64
```