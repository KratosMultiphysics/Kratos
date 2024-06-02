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

// External includes

// Project includes
#include "includes/define.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class SensorDistanceVarianceUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorDistanceVarianceUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorDistanceVarianceUtils(ModelPart& rSensorModelPart);

    ///@}
    ///@name Public operations
    ///@{

    void Initialize();

    double CalculateValue();

    ContainerExpression<ModelPart::NodesContainerType> CalculateGradient() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart* mpSensorModelPart;

    double mMean;

    Matrix mDistances;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/