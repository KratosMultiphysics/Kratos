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

// System includes

// External includes

// Project includes

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "nodal_sensor_specification.h"

namespace Kratos {

NodalSensorSpecification::NodalSensorSpecification(
    const std::string& rName,
    const IndexType NewId,
    const double SensorWeight,
    const NodeType::Pointer pNode)
    : BaseType(rName, NewId),
      mpNode(pNode)
{
    KRATOS_TRY

    this->SetValue(SENSOR_NODE_ID, static_cast<int>(pNode->Id()));
    this->SetValue(SENSOR_WEIGHT, SensorWeight);

    KRATOS_CATCH("");
}

Point NodalSensorSpecification::GetLocation() const
{
    return *mpNode;
}

NodalSensorSpecification::NodeType::Pointer NodalSensorSpecification::GetNode() const
{
    return mpNode;
}

} /* namespace Kratos.*/