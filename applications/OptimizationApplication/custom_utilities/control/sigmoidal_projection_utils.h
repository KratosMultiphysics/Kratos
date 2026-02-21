//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) SigmoidalProjectionUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    static TensorAdaptor<double>::Pointer ProjectForward(
        const TensorAdaptor<double>& rInputTensorAdaptor,
        const std::vector<double>& rXValues,
        const std::vector<double>& rYValues,
        const double Beta,
        const int PenaltyFactor);

    static TensorAdaptor<double>::Pointer ProjectBackward(
        const TensorAdaptor<double>& rInputTensorAdaptor,
        const std::vector<double>& rXValues,
        const std::vector<double>& rYValues,
        const double Beta,
        const int PenaltyFactor);

    static TensorAdaptor<double>::Pointer CalculateForwardProjectionGradient(
        const TensorAdaptor<double>& rInputTensorAdaptor,
        const std::vector<double>& rXValues,
        const std::vector<double>& rYValues,
        const double Beta,
        const int PenaltyFactor);

    ///@}
};

///@}
}