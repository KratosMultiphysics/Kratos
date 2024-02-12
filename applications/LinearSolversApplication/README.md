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

The application provides the following generalized eigensystem `Ax=λBx` solver for sparse matrices.

| Python class                           | solver_type                | Matrix kind A | Matrix kind B | Domain   | Dependencies |
|----------------------------------------|----------------------------|:-------------:|:-------------:|:--------:|:------------:|
| EigensystemSolver                      | `eigen_eigensystem`        | Symmetric     | SPD*          | Real     | None         |
| SpectraSymGEigsShiftSolver             | `spectra_sym_g_eigs_shift` | Symmetric     | SPD*          | Real     | None         |
| FEASTGeneralEigensystemSolver**        | `feast`                    | General       | General       | Real     | Intel® MKL   |
| ComplexFEASTGeneralEigensystemSolver** | `feast_complex`            | General       | General       | Complex  | Intel® MKL   |

*SPD = Symmetric Positive Definite
**A special version for symmetric matrices can be triggered in the solver settings.

`EigensystemSolver` and `SpectraSymGEigsShiftSolver` compute the smallest eigenvalues and corresponding eigenvectors of the system. MKL routines are used automatically if they are available.

`SpectraSymGEigsShiftSolver` interfaces a solver from the [Spectra library](https://spectralib.org/), and has a shift mode that can be used to compute the smallest eigenvalues > `shift`.

**Example:**

```json
{
    "solver_type": "spectra_sym_g_eigs_shift",
    "number_of_eigenvalues": 3,
    "max_iteration": 1000,
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

    ```console
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
    call "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\env\vars.bat" intel64 lp64
    ```

    **Linux:**

    ```console
    source /opt/intel/oneapi/setvars.sh intel64
    ```

2. Add the following flag to CMake to your configure script:

    **Windows:**

    ```batch
    -DUSE_EIGEN_MKL=ON ^
    ```

    **Linux:**

    ```console
    -DUSE_EIGEN_MKL=ON \
    ```

3. Build Kratos

4. Usage:

    **Windows:**

    ```batch
    call "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\env\vars.bat" intel64 lp64
    ```

    **Linux:**

    Set the environment before using MKL

    ```console
    source /opt/intel/oneapi/setvars.sh intel64
    ```

## Install MKL on Ubuntu with apt

Intel MKL can be installed with apt on Ubuntu. A guide can be found in [here](https://neelravi.com/post/intel-oneapi-install/).
For example to install the MKL 2022 version

```console
sudo bash
# <type your user password when prompted.  this will put you in a root shell>
# If they are not installed, you can install using the following command:
sudo apt update
sudo apt -y install cmake pkg-config build-essential
# use wget to fetch the Intel repository public key
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# add to your apt sources keyring so that archives signed with this key will be trusted.
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# remove the public key
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# Configure apt client to use Intel repository
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
# Install all MKL related dependencies. You can install full HPC with: sudo apt install intel-hpckit
sudo apt install intel-oneapi-mkl-devel
# Exit
exit
```

To enable the MKL environment (needs to be done before build/run) use

```console
source /opt/intel/oneapi/setvars.sh intel64
```
