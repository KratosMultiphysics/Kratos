//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLUID_ADJOINT_UTILITIES_H)
#define KRATOS_FLUID_ADJOINT_UTILITIES_H

// System includes
#include <array>

// External includes

// Project includes
#include "containers/variable.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template<unsigned int TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAdjointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static Operations
    ///@{

    template <class TDataType>
    static const Variable<double>& GetRelevantVariable(
        const IndexType DirectionIndex,
        const Variable<TDataType>& rVariable,
        const std::array<const Variable<double>*, 3>& rAllVariableComponents);

    template <class TDataType, std::size_t TGradientVariableTotalDimensionality>
    static std::array<const Variable<double>*, TDim> GetRelevantGradientVariableComponentList(
        const IndexType DirectionIndex,
        const Variable<TDataType>& rVariable,
        const std::array<const Variable<double>*, TGradientVariableTotalDimensionality>& rAllGradientVariableComponents);

    ///@}
};

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_UTILITIES_H
