# LinearSolversApplication

The *LinearSolversApplication* is a thin wrapper for the [Eigen linear algebra library](http://eigen.tuxfamily.org/index.php?title=Main_Page).

## Direct sparse solvers

The application provides the following direct sparse solvers:

| Python class             | solver_type            | Matrix kind | Domain   | Dependencies |
|--------------------------|------------------------|:-----------:|:--------:|:------------:|
| SparseLUSolver           | `sparse_lu`            | Square      | Real     | None         |
| SparseQRSolver           | `sparse_qr`            | Rectangular | Real     | None         |
| SparseCGSolver           | `sparse_cg`            | SPD*        | Real     | None         |
| PardisoLLTSolver         | `pardiso_llt`          | SPD*        | Real     | Intel® MKL   |
| PardisoLDLTSolver        | `pardiso_ldlt`         | SPD*        | Real     | Intel® MKL   |
| PardisoLUSolver          | `pardiso_lu`           | Square      | Real     | Intel® MKL   |
| ComplexSparseLUSolver    | `sparse_lu_complex`    | Square      | Complex  | None         |
| ComplexPardisoLLTSolver  | `pardiso_llt_complex`  | SPD*        | Complex  | Intel® MKL   |
| ComplexPardisoLDLTSolver | `pardiso_ldlt_complex` | SPD*        | Complex  | Intel® MKL   |
| ComplexPardisoLUSolver   | `pardiso_lu_complex`   | Square      | Complex  | Intel® MKL   |

*SPD = Symmetric Positive Definite

**Example**:

```json
{
    "solver_type": "eigen_sparse_lu"
}
```

## Direct dense solvers

The application provides the following direct solvers for dense systems of equations:

| Python class                            | solver_type                            | Matrix requirements | Domain  | Dependencies |
| --------------------------------------- | -------------------------------------- | :-----------------: | :-----: | :----------: |
| DenseColPivHouseholderQRSolver**        | `dense_col_piv_householder_qr`         |        None         |  Real   |     None     |
| DenseHouseholderQRSolver**              | `dense_householder_qr`                 |        None         |  Real   |     None     |
| DenseLLTSolver**                        | `dense_llt`                            |        SPD*         |  Real   |     None     |
| DensePartialPivLUSolver**               | `dense_partial_piv_lu`                 |     Invertible      |  Real   |     None     |
| ComplexDenseColPivHouseholderQRSolver   | `complex_dense_col_piv_householder_qr` |        None         | Complex |     None     |
| ComplexDenseHouseholderQRSolver         | `complex_dense_householder_qr`         |        None         | Complex |     None     |
| ComplexDensePartialPivLUSolver          | `complex_dense_partial_piv_lu`         |     Invertible      | Complex |     None     |

*SPD = Symmetric Positive Definite

**Can also be used to solve equation systems with multiple right hand sides.

## Generalized eigensystem solvers

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
If the application is compiled with MKL, [FEAST 4.0](http://www.ecs.umass.edu/~polizzi/feast/) can be used to solve the generalized eigenvalue problem for real and complex systems (symmetric or unsymmetric). The cmake switch `USE_EIGEN_FEAST` must be set to `ON` with
```batch
-DUSE_EIGEN_FEAST=ON \
```

**Example:**
```json
{
    "solver_type": "feast",
    "symmetric": true,
    "number_of_eigenvalues": 3,
    "search_lowest_eigenvalues": true,
    "e_min" : 0.0,
    "e_max" : 0.2
}
```

## Build instructions

1. Set the required definitions for cmake

    As any other app:

    **Windows:** in `configure.bat`

    ```batch
    set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\LinearSolversApplication;
    ```

    **Linux:** in `configure.sh`

    ```bash
    add_app ${KRATOS_APP_DIR}/LinearSolversApplication
    ```

2. Build Kratos

3. Setup the `ProjectParameters.json`

    ```json
    "linear_solver_settings": {
        "solver_type" : "LinearSolversApplication.sparse_lu"
    }
    ```

4. Run the simulation

## Enable MKL (optional)

In case you have installed [MKL](https://software.intel.com/en-us/mkl) (see below), you can also use the Pardiso solvers.

1. Run the MKL setup script before building Kratos:

    **Windows:**

    ```batch
    call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" intel64 lp64
    ```

    **Linux:**

    ```bash
    source /opt/intel/mkl/bin/mklvars.sh intel64 lp64
    ```

2. Add the following flag to CMake to your configure script:

    **Windows:**

    ```batch
    -DUSE_EIGEN_MKL=ON ^
    ```

    **Linux:**

    ```bash
    -DUSE_EIGEN_MKL=ON \
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
    ```bash
    source /opt/intel/mkl/bin/mklvars.sh intel64 lp64
    ```

## Install MKL on Ubuntu with apt

Intel MKL can be installed with apt on Ubuntu. A guide can be found in https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo.
For example to install the MKL 2020 version

```bash
sudo bash
# <type your user password when prompted.  this will put you in a root shell>
# cd to /tmp where this shell has write permission
cd /tmp
# now get the key:
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# now install that key
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# now remove the public key file exit the root shell
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# add the repository for mkl only (other packages are not needed here)
sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
apt-get update -y
# install specific version of intel-mkl
apt-get install -y intel-mkl-2020.0-088
exit
```

To enable the MKL environment (needs to be done before build/run) use

```bash
source /opt/intel/mkl/bin/mklvars.sh intel64 lp64
```
