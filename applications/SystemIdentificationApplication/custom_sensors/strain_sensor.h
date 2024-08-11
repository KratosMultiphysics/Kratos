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
#include "includes/element.h"

// Application includes
#include "sensor.h"


namespace Kratos
{
///@addtogroup SystemIdentificationApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) StrainSensor : public Sensor
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = Sensor;

    KRATOS_CLASS_POINTER_DEFINITION(StrainSensor);

    ///@}
    ///@name Enums
    ///@{

    enum StrainType
    {
        STRAIN_XX = 0,
        STRAIN_YY = 4,
        STRAIN_ZZ = 8,
        STRAIN_XY = 1,
        STRAIN_XZ = 2,
        STRAIN_YZ = 5
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    StrainSensor(
        const std::string& rName,
        Node::Pointer pNode,
        const Variable<Matrix>& rStrainVariable,
        const StrainType& rStrainType,
        const Element& rElement,
        const double Weight);

    /// Destructor.
    ~StrainSensor() override = default;

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

    const IndexType mElementId;

    StrainType mStrainType;

    Point mLocalPoint;

    const Variable<Matrix>& mrStrainVariable;

    ///@}
    ///@name Private operations
    ///@{

    void SetVectorToZero(
        Vector& rVector,
        const IndexType Size);

    double CalculateStrainDirectionalSensitivity(
        const double Perturbation,
        const Variable<double>& rPerturbationVariable,
        ModelPart::NodeType& rNode,
        ModelPart::ElementType& rElement,
        std::vector<Matrix>& rPerturbedStrains,
        const std::vector<Matrix>& rRefStrains,
        const ProcessInfo& rProcessInfo) const;

    ///@}
};

///@} // Kratos Classes

///@} //Digital Twin Application group

} /* namespace Kratos.*/
