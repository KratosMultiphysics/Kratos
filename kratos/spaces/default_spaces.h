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
// Under the Eigen backend both the real sparse system space AND the real
// dense/local space switch to Eigen (matching the Eigen-backed
// Kratos::Matrix/Vector aliases the Element/Condition virtual interfaces are
// written against, see includes/ublas_interface.h). The complex spaces stay
// uBLAS in both modes: they are only used by the eigenvalue-related solvers
// (e.g. skyline_lu_complex), and there is no Eigen complex space.
//
// NOTE: the backend define changes the meaning of these aliases and hence the
// mangled names of everything instantiated with them. The define is set
// globally by the root CMakeLists.txt; never mix binaries compiled with
// different KRATOS_LINEAR_ALGEBRA_BACKEND values.
///@{

#ifdef KRATOS_USE_EIGEN_BACKEND
namespace Internals {
// Real spaces switch to Eigen; complex ones stay uBLAS (see note above)
template<class TDataType>
struct DefaultSparseSpaceSelector { using Type = TEigenSparseSpace<TDataType>; };
template<class TDataType>
struct DefaultSparseSpaceSelector<std::complex<TDataType>> { using Type = TUblasSparseSpace<std::complex<TDataType>>; };

template<class TDataType>
struct DefaultDenseSpaceSelector { using Type = TEigenDenseSpace<TDataType>; };
template<class TDataType>
struct DefaultDenseSpaceSelector<std::complex<TDataType>> { using Type = TUblasDenseSpace<std::complex<TDataType>>; };
} // namespace Internals

template<class TDataType>
using TDefaultSparseSpace = typename Internals::DefaultSparseSpaceSelector<TDataType>::Type;

template<class TDataType>
using TDefaultDenseSpace = typename Internals::DefaultDenseSpaceSelector<TDataType>::Type;
#else
template<class TDataType>
using TDefaultSparseSpace = TUblasSparseSpace<TDataType>;

template<class TDataType>
using TDefaultDenseSpace = TUblasDenseSpace<TDataType>;
#endif

using DefaultSparseSpaceType = TDefaultSparseSpace<double>;
using DefaultLocalSpaceType = TDefaultDenseSpace<double>;

/// Resolve to uBLAS in both modes (see note above)
using DefaultComplexSparseSpaceType = TDefaultSparseSpace<std::complex<double>>;
using DefaultComplexLocalSpaceType = TDefaultDenseSpace<std::complex<double>>;

///@}

} // namespace Kratos
