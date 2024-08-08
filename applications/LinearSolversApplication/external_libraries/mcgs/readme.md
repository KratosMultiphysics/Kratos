# MCGS - Multicolor Gauss-Seidel Smoother

MCGS is a lightweight library for performing parallelized Gauss-Sidel smoothing, focusing on sparse systems with imbalanced topologies. Implementations for related tasks such as graph coloring and reordering are included as well.

<p style="display: flex; justify-content: center; align-items: center;">
<img src=".github/assets/matrix.png" width=300/>
<img src=".github/assets/right_arrow.png" width=100 style="object-fit: scale-down;">
<img src=".github/assets/matrix_reordered.png" width=300/>
</p>

## Features

- Gauss-Seidel iterations in parallel (OpenMP) or serial
- SOR (successive over-relaxation) in parallel (OpenMP) or serial
- Graph coloring in parallel (loose implementation of `doi:10.1006/jagm.2000.1097`)
- Matrix/vector reordering and reverse reordering

## Usage

Typical workflow

- construct your linear system $A x = b$
- construct and adaptor for $A$ (see [mcgs::CSRAdaptor](structmcgs_1_1CSRAdaptor.html))
- compute a coloring of $A$ (see [mcgs::color](namespacemcgs.html#ad660f970843b8c8edea18c6e9291f6e5))
- construct a partition of $A$ with respect to the coloring (see [mcgs::makePartition](namespacemcgs.html#adbeb4189f3eadcb713e803cf94aa38cf))
- reorder the system with respect to the coloring (see [mcgs::reorder](namespacemcgs.html#a5291808c16a69190ac0bb31a1f3ee81d))
- perform Gauss-Seidel iterations **using the reordered partition** (see [mcgs::solve](namespacemcgs.html#ae862fac411e001950f012872f6ac7e0c))
- *optional: restore the original order of your system* (see [mcgs::revertReorder](namespacemcgs.html#aa5b1a78cfa8d230b2100320dde50f3c7))
- deallocate partitions (see [mcgs::destroyPartition](namespacemcgs.html#ad619ded9f67d8a9f379ad7e4b759d854))


### C++ Example Snippet

```cpp
#include "mcgs/mcgs.hpp"

...

// Any CSR matrix will do but for the sake of familiarity, let's assume you're using Eigen.
Eigen::SparseMatrix<double,Eigen::RowMajor> A;  // <== left hand side matrix
Eigen::Matrix<double,Eigen::Dynamic,1> b;       // <== right hand side vector

// Construct an adaptor for your matrix.
mcgs::CSRAdaptor</*index type=*/ int, /*value type=*/double> adaptor;
adaptor.rowCount        = A.rows();
adaptor.columnCount     = A.cols();
adaptor.entryCount      = A.nonZeros();
adaptor.pRowExtents     = A.outerIndexPtr();
adaptor.pColumnIndices  = A.innerIndexPtr();
adaptor.pEntries        = A.innerNonZeroPtr();

// Color the rows of your matrix.
std::vector<unsigned> colors(adaptor.rowCount);
mcgs::color(colors.data(),
            adaptor,
            mcgs::ColorSettings {});

// Construct a partition for your matrix with respect to the coloring,
// and reorder the system accordingly. Note that this mutates your original matrix!
auto pPartition = mcgs::makePartition(colors.data(), adaptor.rowCount);
auto pReorderedPartition = mcgs::reorder(A.rows(), A.cols(), A.nonZeros(),
                                         A.outerIndexPtr(), A.innerIndexPtr(), A.innerNonZeroPtr(),
                                         b.data());

// Do 10 Gauss-Seidel iterations.
std::vector<double> x(adaptor.columnCount);
mcgs::SolveSettings</*index type=*/int, /*value type=*/double> settings;
settings.maxIterations = 10;
settings.parallelization = mcgs::Parallelization::RowWise; // <== default parallelization strategy, check out the other ones as well.
mcgs::solve(x.data(), adaptor, b.data(), pReorderedPartition, settings);

// Optional: if you need to recover your original system,
//           you need to undo the reordering.
// See mcgs::revertReorder

// Cleanup
mcgs::destroyPartition(pPartition);
mcgs::destroyPartition(pReorderedPartition);
```

## Requirements

- C++ compiler with full C++17 support (GCC or Clang are tested)
- CMake version 3.15 or later
- [optional] OpenMP 2.0 or later for shared memory parallelization

## Installation

MCGS is written in C++ and uses CMake as a build system. Building produces a single shared library and matching header.

### Build and Install
  ```bash
  cmake <path-to-mcgs-source>                       \
        -B<path-to-build-dir>                       \
        -DCMAKE_INSTALL_PREFIX=<your-install-dir>   \
        -DCMAKE_BUILD_TYPE=Release                  \
        -DMCGS_SHARED_MEMORY_PARALLELISM=OpenMP

  cmake <path-to-build-dir> --build --target install
  ```

### Build options

- CMake options for `MCGS_SHARED_MEMORY_PARALLELISM`
  - `None`: no parallelization
  - `OpenMP`: use OpenMP for shared memory parallelization

### Include MCGS in a CMake project

```cmake
find_package(MCGS REQUIRED)
target_link_libraries(<your_target> PRIVATE mcgs)
```
