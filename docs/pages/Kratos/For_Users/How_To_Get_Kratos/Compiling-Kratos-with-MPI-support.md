---
title: Compiling Kratos with MPI support
keywords: 
tags: [Compiling-Kratos-with-MPI-support.md]
sidebar: kratos_for_users
summary: 
---

# Compiling Kratos with MPI Support

Kratos supports distributed-memory parallelism through **MPI**, requiring the following dependencies:

- **[METIS Application](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MetisApplication#readme)**: Provides mesh partitioning capabilities through the [`METIS`](https://github.com/KarypisLab/METIS) library.
- **[Trilinos Application](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/TrilinosApplication#readme)**: Integrates the [`Trilinos`](https://github.com/trilinos/Trilinos) libraries to manage distributed-memory matrices, vectors, and linear solvers.

The following sections describe how to obtain and compile these dependencies and how to configure Kratos for MPI execution.

---

## 1. Installing Dependencies

Dependencies can be obtained either through your Linux distribution’s package manager or by compiling from source.

### Option A: Installing from Packages (Recommended)

On Ubuntu/Debian-based distributions:

- **Install METIS and Trilinos (all packages)**:
  ```bash
  sudo apt install libmetis-dev trilinos-all-dev
  ```

- **Or, minimal Trilinos packages required by Kratos**:
  ```bash
  sudo apt install \
      libtrilinos-amesos-dev \
      libtrilinos-amesos2-dev \
      libtrilinos-aztecoo-dev \
      libtrilinos-epetra-dev \
      libtrilinos-epetraext-dev \
      libtrilinos-ifpack-dev \
      libtrilinos-ml-dev \
      libtrilinos-teuchos-dev \
      libtrilinos-tpetra-dev \
      libtrilinos-kokkos-dev \
      libtrilinos-kokkos-kernels-dev \
      libtrilinos-shylu-dev
  ```

### Option B: Compiling Dependencies from Source

If packages are unavailable, outdated, or if administrative privileges are lacking, build dependencies manually.

**Recommended directory structure:**

```bash
export ROOT=${HOME}/Projects
```

Tested environment: `cmake v3.24`, `gcc v10.2`

---

### Step-by-Step Compilation of Dependencies

#### Step 1: Compiling GKlib

`GKlib` is required by `METIS`:

```bash
# Clone GKlib
cd ${ROOT}
git clone --depth 1 https://github.com/KarypisLab/GKlib.git
cd GKlib
git checkout a7f8172703cf6e999dd0710eb279bba513da4fec

# Build GKlib
cmake \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_INSTALL_PREFIX=${ROOT}/GKlib/build \
  -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=199309L" \
  .
make install
```

---

#### Step 2: Compiling METIS

Patch `METIS` CMake configuration as follows:

```bash
cd ${ROOT}
git clone --depth 1 https://github.com/KarypisLab/METIS.git 
cd METIS
git checkout 94c03a6e2d1860128c2d0675cbbb86ad4f261256

# Patch CMakeLists.txt
cd libmetis
cat <<EOF > CMakeLists.txt
include_directories(.)

file(GLOB metis_sources *.c)

add_library(metis \${METIS_LIBRARY_TYPE} \${metis_sources})
target_link_libraries(metis PUBLIC GKlib)

if(METIS_INSTALL)
  install(TARGETS metis
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)
endif()
EOF

cd ..
make config shared=1 prefix=${ROOT}/METIS/build gklib_path=${ROOT}/GKlib/build
make install
```

---

#### Step 3: Compiling Trilinos

Build Trilinos only with the packages necessary for Kratos. Ensure that BLAS, LAPACK, and MPI (OpenMPI recommended) are available.

```bash
cd ${ROOT}
git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
git checkout 9c03d9b

module load intel/oneapi2021/mkl/2021.1.1  # Provides LAPACK/BLAS
module load openmpi  # MPI support

mkdir -p build
cd build

cmake \
  -D CMAKE_INSTALL_PREFIX=${ROOT}/Trilinos/build \
  -D CMAKE_BUILD_TYPE="RelWithDebInfo" \
  -D TPL_ENABLE_MPI=ON \
  -D BUILD_SHARED_LIBS=ON \
  -D Trilinos_ENABLE_Epetra=ON \
  -D Trilinos_ENABLE_EpetraExt=ON \
  -D Trilinos_ENABLE_Triutils=ON \
  -D Trilinos_ENABLE_AztecOO=ON \
  -D Trilinos_ENABLE_Ifpack=ON \
  -D Trilinos_ENABLE_ML=ON \
  -D Trilinos_ENABLE_Amesos=ON \
  -D Trilinos_ENABLE_Amesos2=ON \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -D TPL_ENABLE_MKL=OFF \
  -D BLAS_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
  -D BLAS_LIBRARY_NAMES="mkl_rt" \
  -D LAPACK_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
  -D LAPACK_LIBRARY_NAMES="mkl_rt" \
  ..
make -j $(nproc) install
```

---

## 2. Compiling Kratos with MPI

Once dependencies are set up, configure and compile Kratos:

```bash
#!/bin/bash

# Define application selection function
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Load required modules
module purge
module load boost python git intel/oneapi2021/mkl/2021.1.1 openmpi gcc/10.2.0

# Set environment variables
export KRATOS_BRANCH="mpi-gcc"
export KRATOS_SOURCE="$( cd "$(dirname "$0")"; pwd -P )/.."
export KRATOS_BIN="${KRATOS_SOURCE}/bin/${KRATOS_BRANCH}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build/${KRATOS_BRANCH}"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"

# Applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/FluidDynamicApplication
add_app ${KRATOS_APP_DIR}/TrilinosApplication
add_app ${KRATOS_APP_DIR}/MetisApplication

# Clean previous build caches
rm -rf "${KRATOS_BUILD}/cmake_install.cmake" "${KRATOS_BUILD}/CMakeCache.txt" "${KRATOS_BUILD}/CMakeFiles"

# Configure Kratos
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}" -GNinja \
  -DUSE_MPI=ON \
  -DKRATOS_SHARED_MEMORY_PARALLELIZATION="OpenMP" \
  -DKRATOS_BUILD_TESTING=OFF \
  -DINSTALL_RUNKRATOS=OFF \
  -DCMAKE_INSTALL_PREFIX="${KRATOS_BIN}" \
  -DTRILINOS_ROOT="${ROOT}/Trilinos/build" \
  -DMETIS_ROOT_DIR="${ROOT}/METIS/build"

# Build and install
cmake --build "${KRATOS_BUILD}" --target install -- -j $(nproc)

# Create symbolic link for current build
ln -sfn ${KRATOS_BRANCH} ${KRATOS_SOURCE}/bin/current
```

---

## 3. Setting Up Runtime Environment

Ensure that `LD_LIBRARY_PATH` includes all compiled libraries:

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOT}/Trilinos/build/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOT}/GKlib/build/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOT}/METIS/build/lib
```

---

## 4. Running Kratos with MPI

To run Kratos simulations in parallel, edit `ProjectParameters.json`:

```json
{
  "problem_data": {
    ...
    "parallel_type": "MPI",
    ...
  },
  ...
}
```

Now, you are ready to execute Kratos with MPI:

```bash
mpirun -np <number_of_processes> python3 run_simulation.py
```