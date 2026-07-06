//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <complex>

// External includes

// Project includes
#include "spaces/ublas_space.h"
#ifdef KRATOS_USE_EIGEN_BACKEND
#include "spaces/eigen_space.h"
#endif

namespace Kratos {

///@name Type Definitions
///@{
// Default linear-algebra spaces, resolved at configure time through the
// KRATOS_LINEAR_ALGEBRA_BACKEND CMake option ("ublas", the default, or
// "eigen" which defines KRATOS_USE_EIGEN_BACKEND).
//
// Only the real *sparse* system space switches backend:
// - The dense/local space stays uBLAS in both modes, because the
//   Element/Condition virtual interfaces hardcode Kratos::Matrix/Vector
//   (includes/element.h) and the schemes pass TDenseSpace::MatrixType
//   directly into CalculateLocalSystem.
// - The complex spaces stay uBLAS in both modes: they are only used by the
//   eigenvalue-related solvers (e.g. skyline_lu_complex), which gain nothing
//   from the switch.
//
// NOTE: the backend define changes the meaning of these aliases and hence the
// mangled names of everything instantiated with them. The define is set
// globally by the root CMakeLists.txt; never mix binaries compiled with
// different KRATOS_LINEAR_ALGEBRA_BACKEND values.
///@{

#ifdef KRATOS_USE_EIGEN_BACKEND
namespace Internals {
// Real sparse spaces switch to Eigen; complex ones stay uBLAS (see note above)
template<class TDataType>
struct DefaultSparseSpaceSelector { using Type = TEigenSparseSpace<TDataType>; };
template<class TDataType>
struct DefaultSparseSpaceSelector<std::complex<TDataType>> { using Type = TUblasSparseSpace<std::complex<TDataType>>; };
} // namespace Internals

template<class TDataType>
using TDefaultSparseSpace = typename Internals::DefaultSparseSpaceSelector<TDataType>::Type;
#else
template<class TDataType>
using TDefaultSparseSpace = TUblasSparseSpace<TDataType>;
#endif

/// The dense/local space does not switch (see note above)
template<class TDataType>
using TDefaultDenseSpace = TUblasDenseSpace<TDataType>;

using DefaultSparseSpaceType = TDefaultSparseSpace<double>;
using DefaultLocalSpaceType = TDefaultDenseSpace<double>;

/// Resolve to uBLAS in both modes (see note above)
using DefaultComplexSparseSpaceType = TDefaultSparseSpace<std::complex<double>>;
using DefaultComplexLocalSpaceType = TDefaultDenseSpace<std::complex<double>>;

///@}

} // namespace Kratos
