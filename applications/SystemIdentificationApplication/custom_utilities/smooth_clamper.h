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
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class SmoothClamper
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SmoothClamper);

    ///@}
    ///@name Life cycle
    ///@{

    SmoothClamper(
        const double Min,
        const double Max);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clamps the given input expression smoothly.
     * @details The input value (x) can be in [-infty, +infty] range. But the output (y) will be in
     *          [Min, Max] range.
     * @param rX                                    Input expression in input space.
     * @return ContainerExpression<TContainerType>  Clamped output expression in clamped space.
     */
    ContainerExpression<TContainerType> ProjectForward(const ContainerExpression<TContainerType>& rX) const;

    /**
     * @brief Compute the derivative dy/dx at given X for clamping.
     * @details Computes the clamped derivative where clamped value is Y and input is X, and the
     *          derivative is dY/dX
     * @param rX                                    Input expression in input space.
     * @return ContainerExpression<TContainerType>  Output expression with dY/dX.
     */
    ContainerExpression<TContainerType> CalculateForwardProjectionGradient(const ContainerExpression<TContainerType>& rX) const;

    /**
     * @brief Computes the X given Y.
     * @param rY                                    Input expression in clamped space.
     * @return ContainerExpression<TContainerType>  Output expression in input space.
     */
    ContainerExpression<TContainerType> ProjectBackward(const ContainerExpression<TContainerType>& rY) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mMin;

    const double mMax;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/