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

#pragma once

// System includes
#include <functional>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes

namespace Kratos
{
///@name Classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAdjointTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using ConditionType = ModelPart::ConditionType;

    using ElementType = ModelPart::ElementType;

    ///@}
    ///@name Static Operations
    ///@{

    template<class TDataType>
    static TDataType CalculateRelaxedVariableRate(
        const double BossakAlpha,
        const Variable<TDataType>& rVariable,
        const NodeType& rNode);

    template<class TDataType>
    static void RunAdjointEntityDerivativesTest(
        ModelPart& rPrimalModelPart,
        ModelPart& rAdjointModelPart,
        const std::function<void(ModelPart&)>& rUpdateModelPart,
        const Variable<TDataType>& rVariable,
        const std::function<void(Matrix&, ConditionType&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
        const IndexType EquationOffset,
        const IndexType DerivativeOffset,
        const double Delta,
        const double Tolerance);

    template<class TDataType>
    static void RunAdjointEntityDerivativesTest(
        ModelPart& rPrimalModelPart,
        ModelPart& rAdjointModelPart,
        const std::function<void(ModelPart&)>& rUpdateModelPart,
        const Variable<TDataType>& rVariable,
        const std::function<void(Matrix&, ElementType&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
        const IndexType EquationOffset,
        const IndexType DerivativeOffset,
        const double Delta,
        const double Tolerance);

    ///@}
};

///@}

} // namespace Kratos