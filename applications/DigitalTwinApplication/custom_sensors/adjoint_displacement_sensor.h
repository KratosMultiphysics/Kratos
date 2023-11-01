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
#include "containers/array_1d.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "adjoint_sensor.h"


namespace Kratos
{
///@addtogroup DigitalTwinApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(DIGITAL_TWIN_APPLICATION) AdjointDisplacementSensor : public AdjointSensor
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(AdjointDisplacementSensor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointDisplacementSensor(
        Model& rModel,
        Parameters SensorSettings);

    /// Destructor.
    virtual ~AdjointDisplacementSensor() override = default;

    ///@}
    ///@name Operations
    ///@{

    void SetSensorSpecification(SensorSpecification& rSensorSpecification) override;

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

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;

    array_1d<double, 3> mDirection;

    double mWeight;

    int mNodeId;

    int mElementId;

    int mNodeIndex;

    bool mIsNodeAvailable;

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
