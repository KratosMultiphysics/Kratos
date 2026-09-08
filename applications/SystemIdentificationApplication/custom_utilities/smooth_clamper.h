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
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SmoothClamper
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
     * @brief Clamps the given value smoothly.
     * @details The input value (x) can be in [-infty, +infty] range. But the output (y) will be in
     *          [Min, Max] range.
     *
     * @param X            Value to be clamped.
     * @return double       Clamped value.
     */
    double ProjectForward(const double X) const;

    /**
     * @brief Compute the derivative dy/dx at given X for clamping.
     * @details Computes the clamped derivative where clamped value is Y and input is X, and the
     *          derivative is dY/dX
     *
     * @param X            Input value in input space.
     * @return double       Output with dY/dX.
     */
    double CalculateForwardProjectionGradient(const double X) const;

    /**
     * @brief Computes the X given Y.
     * @param Y        Input value in clamped space.
     * @return double  Output value in input space.
     */
    double ProjectBackward(const double Y) const;

    /**
     * @brief Clamps the given input tensor adaptor smoothly.
     * @details The input value (x) can be in [-infty, +infty] range. But the output (y) will be in
     *          [Min, Max] range.
     * @param rX                                    Input tensor adaptor in input space.
     * @return TensorAdaptor<double>::Pointer       Clamped output tensor adaptor in clamped space.
     */
    TensorAdaptor<double>::Pointer ProjectForward(const TensorAdaptor<double>& rX) const;

    /**
     * @brief Compute the derivative dy/dx at given X for clamping.
     * @details Computes the clamped derivative where clamped value is Y and input is X, and the
     *          derivative is dY/dX
     * @param rX                                    Input tensor adaptor in input space.
     * @return TensorAdaptor<double>::Pointer       Output tensor adaptor with dY/dX.
     */
    TensorAdaptor<double>::Pointer CalculateForwardProjectionGradient(const TensorAdaptor<double>& rX) const;

    /**
     * @brief Computes the X given Y.
     * @param rY                                    Input tensor adaptor in clamped space.
     * @return TensorAdaptor<double>::Pointer       Output tensor adaptor in input space.
     */
    TensorAdaptor<double>::Pointer ProjectBackward(const TensorAdaptor<double>& rY) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mMin;

    const double mMax;

    const double mDelta;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/