//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <array>
#include <type_traits>

// External includes

// Project includes

// Application includes
#include "optimization_application_variables.h"


namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes>
class HelmholtzVectorSurfaceDataContainer
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    static constexpr IndexType NumberOfNodes = TNumNodes;

    static constexpr auto TargetVariablesList = (TDim == 2) ?
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_X, &HELMHOLTZ_VECTOR_Y} :
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_X, &HELMHOLTZ_VECTOR_Y, &HELMHOLTZ_VECTOR_Z};

    static constexpr auto SourceVariablesList = (TDim == 2) ?
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_SOURCE_X, &HELMHOLTZ_VECTOR_SOURCE_Y} :
                                                    std::array<const Variable<double>*, TDim>{&HELMHOLTZ_VECTOR_SOURCE_X, &HELMHOLTZ_VECTOR_SOURCE_Y, &HELMHOLTZ_VECTOR_SOURCE_Z};

    ///@}

};

}