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

namespace Kratos {

///@name Type Definitions
///@{

template<class TDataType>
using TDefaultSparseSpace = TUblasSparseSpace<TDataType>;

template<class TDataType>
using TDefaultDenseSpace = TUblasDenseSpace<TDataType>;

using DefaultSparseSpaceType = TDefaultSparseSpace<double>;
using DefaultLocalSpaceType = TDefaultDenseSpace<double>;

/// Resolve to uBLAS in both modes (see note above)
using DefaultComplexSparseSpaceType = TDefaultSparseSpace<std::complex<double>>;
using DefaultComplexLocalSpaceType = TDefaultDenseSpace<std::complex<double>>;

///@}

} // namespace Kratos
