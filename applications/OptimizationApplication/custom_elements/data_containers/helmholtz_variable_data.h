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

// Application includes
#include "optimization_application_variables.h"

namespace Kratos {

template <unsigned int TDataDimension>
struct HelmholtzVariableData;

template<>
struct HelmholtzVariableData<1> {
    static constexpr auto TargetVariablesList = std::array<const Variable<double>*, 1>{&HELMHOLTZ_SCALAR};
    static constexpr auto SourceVariablesList = std::array<const Variable<double>*, 1>{&HELMHOLTZ_SCALAR_SOURCE};
};

template<>
struct HelmholtzVariableData<2> {
    static constexpr auto TargetVariablesList = std::array<const Variable<double>*, 2>{&HELMHOLTZ_VECTOR_X, &HELMHOLTZ_VECTOR_Y};
    static constexpr auto SourceVariablesList = std::array<const Variable<double>*, 2>{&HELMHOLTZ_VECTOR_SOURCE_X, &HELMHOLTZ_VECTOR_SOURCE_Y};
};

template<>
struct HelmholtzVariableData<3> {
    static constexpr auto TargetVariablesList = std::array<const Variable<double>*, 3>{&HELMHOLTZ_VECTOR_X, &HELMHOLTZ_VECTOR_Y, &HELMHOLTZ_VECTOR_Z};
    static constexpr auto SourceVariablesList = std::array<const Variable<double>*, 3>{&HELMHOLTZ_VECTOR_SOURCE_X, &HELMHOLTZ_VECTOR_SOURCE_Y, &HELMHOLTZ_VECTOR_SOURCE_Z};
};

} // namespace Kratos