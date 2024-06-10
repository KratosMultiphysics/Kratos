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
#include <vector>

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

class KRATOS_API(DIGITAL_TWIN_APPLICATION) EigenvectorSensor : public Sensor
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = Sensor;

    KRATOS_CLASS_POINTER_DEFINITION(EigenvectorSensor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigenvectorSensor(
        const std::string& rName,
        const Point& rLocation,
        const double Weight,
        const Vector& rSensorValueVector);

    /// Destructor.
    ~EigenvectorSensor() override = default;

    ///@}
    ///@name Operations
    ///@{

    static Parameters GetDefaultParameters();

    const Parameters GetSensorParameters() const override;

        /**
     * @brief Set the Sensor Value
     *
     * @param Value         Value to be set.
     */
    void SetSensorValueVector(const Vector Value) override;

    /**
     * @brief Get the Sensor value
     *
     * @return Vector       Value of the sensor.
     */
    Vector GetSensorValueVector() const override;



    /**
     * @brief Calculate the current value of the sensor using the given model part.
     *
     * @param rModelPart        Model part to calculate the sensor value.
     * @return Vector        Calculated sensor value.
     */
    Vector CalculateValueVector(ModelPart& rModelPart) override;

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

    void CalculateGlobalPartialSensitivity(
        Element& rAdjointElement,
        const Variable<double>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo,
        ModelPart& rModelPart) override;

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

    Vector mSensorValueVector;

    std::vector<Vector> mEigenVectorDerivative;


    //array_1d<double> mSensorValueVector;

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
        const ProcessInfo& rProcessInfo,
        ModelPart& rModelPart);

    ///@}
};

///@} // Kratos Classes

///@} //Digital Twin Application group

} /* namespace Kratos.*/
