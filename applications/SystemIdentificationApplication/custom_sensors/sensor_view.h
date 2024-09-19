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
#include <string>
#include <vector>

// External includes

// Project includes
#include "expression/container_expression.h"
#include "includes/define.h"

// Application includes
#include "sensor.h"

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorView
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorView);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensorView(
        Sensor::Pointer pSensor,
        const std::string& rExpressionName);

    ///@}
    ///@name Public operations
    ///@{

    Sensor::Pointer GetSensor() const;

    typename ContainerExpression<TContainerType>::Pointer GetContainerExpression() const;

    std::string GetExpressionName() const;

    void  AddAuxiliaryExpression(
        const std::string& rSuffix,
        typename ContainerExpression<TContainerType>::Pointer pContainerExpression);

    typename ContainerExpression<TContainerType>::Pointer GetAuxiliaryExpression(const std::string& rSuffix) const;

    std::vector<std::string> GetAuxiliarySuffixes() const;

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

    Sensor::Pointer mpSensor;

    const std::string mExpressionName;

    typename ContainerExpression<TContainerType>::Pointer mpContainerExpression;

    ///@}
};

/// output stream functions
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const SensorView<TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@} // Kratos Classes

} /* namespace Kratos.*/