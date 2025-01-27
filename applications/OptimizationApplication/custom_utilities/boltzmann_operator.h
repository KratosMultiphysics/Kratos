//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class BoltzmannOperator
{
public:
    ///@name Life cycle
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(BoltzmannOperator);

    BoltzmannOperator(
        const double Beta,
        const double ScalingFactor);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue() const;

    typename ContainerExpression<TContainerType>::Pointer CalculateGradient() const;

    void Update(typename ContainerExpression<TContainerType>::Pointer pInputContainerExpression);

    ///@}
private:

    ///@name Private member variables
    ///@{

    double mBeta;

    double mScalingFactor;

    double mNumerator;

    double mDenominator;

    typename ContainerExpression<TContainerType>::Pointer mpContainerExpression;

    ///@}
};

///@}
} // namespace Kratos