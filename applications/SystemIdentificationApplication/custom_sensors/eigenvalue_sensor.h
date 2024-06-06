//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Talhah Ansari
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/element.h"

// Application includes
#include "sensor.h"


namespace Kratos
{
///@addtogroup SystemIdentificationApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(DIGITAL_TWIN_APPLICATION) EigenvalueSensor : public Sensor
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = Sensor;

    KRATOS_CLASS_POINTER_DEFINITION(EigenvalueSensor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigenvalueSensor(
        const std::string& rName,
        const Point& rLocation,
        const double Weight);

    /// Destructor.
    ~EigenvalueSensor() override = default;

    ///@}
    ///@name Operations
    ///@{

    static Parameters GetDefaultParameters();

    const Parameters GetSensorParameters() const override;

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

    Point mLocalPoint;

    ///@}
    ///@name Private operations
    ///@{

    void SetVectorToZero(
        Vector& rVector,
        const IndexType Size);

    void DetermineEigenvectorOfElement(
        ModelPart::ElementType& rElement, 
        const int eigenfrequency_id, 
        Vector& rEigenvectorOfElement, 
        const ProcessInfo& CurrentProcessInfo);

    // void CalculateLeftHandSideDerivative(
    //     Element& rElement,
    //     const Matrix& rLHS,
    //     const double& rPertubationSize,
    //     Matrix& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo);

    // void CalculateMassMatrixDerivative(
    //     Element& rElement,
    //     const Matrix& rMassMatrix,
    //     const double& rPerturbationSize,
    //     Matrix& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo);

    void CalculateElementContributionToPartialSensitivity(
        Element& rAdjointElement,
        const std::string& rVariableName,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo);

    ///@}
};

///@} // Kratos Classes

///@} //Digital Twin Application group

} /* namespace Kratos.*/
