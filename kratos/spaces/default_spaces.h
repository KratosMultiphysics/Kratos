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
#ifdef KRATOS_USE_EIGEN_BACKEND
#include "spaces/eigen_space.h"
#else
#include "spaces/ublas_space.h"
#endif

namespace Kratos {

///@name Type Definitions
///@{
// Default linear-algebra spaces, resolved at configure time through the
// KRATOS_LINEAR_ALGEBRA_BACKEND CMake option ("ublas", the default, or
// "eigen" which defines KRATOS_USE_EIGEN_BACKEND). Under the Eigen backend
// every space (real and complex, sparse and dense) is Eigen-backed, matching
// the Eigen-backed Kratos::Matrix/Vector/CompressedMatrix aliases of
// includes/ublas_interface.h.
//
// NOTE: the backend define changes the meaning of these aliases and hence the
// mangled names of everything instantiated with them. The define is set
// globally by the root CMakeLists.txt; never mix binaries compiled with
// different KRATOS_LINEAR_ALGEBRA_BACKEND values.
///@{

#ifdef KRATOS_USE_EIGEN_BACKEND
template<class TDataType>
using TDefaultSparseSpace = TEigenSparseSpace<TDataType>;

template<class TDataType>
using TDefaultDenseSpace = TEigenDenseSpace<TDataType>;
#else
template<class TDataType>
using TDefaultSparseSpace = TUblasSparseSpace<TDataType>;

template<class TDataType>
using TDefaultDenseSpace = TUblasDenseSpace<TDataType>;
#endif

using DefaultSparseSpaceType = TDefaultSparseSpace<double>;
using DefaultLocalSpaceType = TDefaultDenseSpace<double>;

using DefaultComplexSparseSpaceType = TDefaultSparseSpace<std::complex<double>>;
using DefaultComplexLocalSpaceType = TDefaultDenseSpace<std::complex<double>>;

///@}

} // namespace Kratos
