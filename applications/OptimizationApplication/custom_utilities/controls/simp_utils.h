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

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder_base.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) SimpUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType>
    static void ProjectContainerVariableDataHolderForward(
        ContainerVariableDataHolderBase<TContainerType>& rY,
        const ContainerVariableDataHolderBase<TContainerType>& rX,
        const Vector& rYLimits,
        const double Beta,
        const double PenaltyFactor = 1);

    template<class TContainerType>
    static void ProjectContainerVariableDataHolderBackward(
        ContainerVariableDataHolderBase<TContainerType>& rX,
        const ContainerVariableDataHolderBase<TContainerType>& rY,
        const Vector& rYLimits,
        const double Beta,
        const double PenaltyFactor = 1);

    template<class TContainerType>
    static void ProjectContainerVariableDataHolderDerivative(
        ContainerVariableDataHolderBase<TContainerType>& rDerivative,
        const ContainerVariableDataHolderBase<TContainerType>& rX,
        const Vector& rYLimits,
        const double Beta,
        const double PenaltyFactor = 1);

    ///@}
private:
    ///@name Private operations
    ///@{

    static double ProjectForward(
        const double X,
        const Vector& rYLimits,
        const double Beta,
        const double PenaltyFactor = 1);

    static double ProjectBackward(
        const double Y,
        const Vector& rYLimits,
        const double Beta,
        const double PenaltyFactor = 1);

    static double ProjectionDerivative(
        const double X,
        const Vector& rYLimits,
        const double Beta,
        const double PenaltyFactor = 1);

    ///@}
};

///@}
}