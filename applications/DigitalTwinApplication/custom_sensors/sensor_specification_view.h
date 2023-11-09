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
#include <string>

// External includes

// Project includes
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "sensor_specification.h"

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class SensorSpecificationView
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorSpecificationView);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensorSpecificationView(
        SensorSpecification::Pointer pSensorSpecification,
        const std::string& rExpressionName);

    ///@}
    ///@name Public operations
    ///@{

    SensorSpecification::Pointer GetSensorSpecification() const;

    typename ContainerExpression<TContainerType>::Pointer GetContainerExpression() const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const;

    void PrintInfo(std::ostream& rOStream) const;

    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    SensorSpecification::Pointer mpSensorSpecification;

    typename ContainerExpression<TContainerType>::Pointer mpContainerExpression;

    ///@}
};

/// output stream functions
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const SensorSpecificationView<TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@} // Kratos Classes

} /* namespace Kratos.*/