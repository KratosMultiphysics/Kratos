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

// External includes

// Project includes
// #include "arrat"


// Application includes
#include "optimization_application_variables.h"


namespace Kratos
{

template<unsigned int TNumNodes>
class HelmholtzScalarSurfaceDataContainer
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using DataType = double;

    static constexpr IndexType NumberOfNodes = TNumNodes;

    static constexpr std::array<const Variable<double>*, 1> TargetVariablesList = {&HELMHOLTZ_SCALAR};

    static constexpr std::array<const Variable<double>*, 1> SourceVariablesList = {&HELMHOLTZ_SCALAR_SOURCE};

    ///@}

};

}