//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class SensorSensitivityBoltzmannOperatorResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorSensitivityBoltzmannOperatorResponseUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorSensitivityBoltzmannOperatorResponseUtils(
        ModelPart& rSensorModelPart,
        const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rSensitivityDistributions,
        const double Beta);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue();

    ContainerExpression<ModelPart::NodesContainerType> CalculateGradient() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart* mpSensorModelPart;

    double mBeta;

    double mNumerator;

    double mDenominator;

    const std::vector<typename ContainerExpression<TContainerType>::Pointer> mSensorSensitivities;

    Vector mElementRedundancy;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/