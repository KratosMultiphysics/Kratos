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
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "sensor_specification.h"

namespace Kratos {
///@name Kratos Classes
///@{

class NodalSensorSpecification: public SensorSpecification
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = SensorSpecification;

    using IndexType = BaseType::IndexType;

    using NodeType = ModelPart::NodeType;

    KRATOS_CLASS_POINTER_DEFINITION(NodalSensorSpecification);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalSensorSpecification(
        const std::string& rName,
        const IndexType NewId,
        const double SensorWeight,
        const NodeType::Pointer pNode);

    /// Destructor.
    ~NodalSensorSpecification() = default;

    ///@}
    ///@name Public Operations
    ///@{

    Point GetLocation() const override;

    NodeType::Pointer GetNode() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    NodeType::Pointer mpNode;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/