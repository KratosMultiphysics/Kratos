//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_sensors/sensor_specification.h"

namespace Kratos {
///@name Kratos Classes
///@{

class AdjointSensor: public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointSensor);

    using BaseType = AdjointResponseFunction;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointSensor() = default;

    /// Destructor.
    virtual ~AdjointSensor() = default;

    ///@}
    ///@name Operations
    ///@{

    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_ERROR << "Calling base class AdjointSensor::CalculateValue. Please implement it in the derrived class.";
    }

    virtual void SetSensorSpecification(SensorSpecification& rSensorSpecification)
    {
        KRATOS_ERROR << "Calling base class AdjointSensor::SetSensorSpecification. Please implement it in the derrived class.";
    }

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/