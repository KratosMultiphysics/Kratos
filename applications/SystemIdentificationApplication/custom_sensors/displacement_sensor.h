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

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/element.h"

// Application includes
#include "sensor.h"


namespace Kratos
{
///@addtogroup SystemIdentificationApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) DisplacementSensor : public Sensor
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = Sensor;

    KRATOS_CLASS_POINTER_DEFINITION(DisplacementSensor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DisplacementSensor(
        const std::string& rName,
        Node::Pointer pNode,
        const array_1d<double, 3>& rDirection,
        const Element& rElement,
        const double Weight);

    /// Destructor.
    ~DisplacementSensor() override = default;

    ///@}
    ///@name Static operations
    ///@{

    static Sensor::Pointer Create(
        ModelPart& rDomainModelPart,
        ModelPart& rSensorModelPart,
        const IndexType Id,
        Parameters SensorParameters);

    static Parameters GetDefaultParameters();

    ///@}
    ///@name Operations
    ///@{

    Parameters GetSensorParameters() const override;

    double CalculateValue(ModelPart& rModelPart) override;

    void CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<double>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<double>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    const IndexType mElementId;

    array_1d<double, 3> mDirection;

    Vector mNs;

    ///@}
    ///@name Private operations
    ///@{

    void SetVectorToZero(
        Vector& rVector,
        const IndexType Size);

    ///@}
};

///@} // Kratos Classes

///@} //Digital Twin Application group

} /* namespace Kratos.*/
