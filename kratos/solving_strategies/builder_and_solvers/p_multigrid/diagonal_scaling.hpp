//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "includes/kratos_parameters.h" // Parameters
#include "includes/process_info.h" // ProcessInfo


namespace Kratos {


enum class DiagonalScaling
{
    None        = 0,
    AbsMax      = 1,
    Norm        = 2,
    Constant    = 3
}; // enum class DiagonalScaling


DiagonalScaling ParseDiagonalScaling(Parameters Settings);


template <class TSparse>
typename TSparse::DataType
GetDiagonalScaleFactor(const typename TSparse::MatrixType& rMatrix,
                       const DiagonalScaling ScalingStrategy,
                       [[maybe_unused]] const ProcessInfo& rProcessInfo);


} // namespace Kratos
