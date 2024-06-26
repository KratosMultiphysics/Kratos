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
#include "includes/element.h"
#include "includes/condition.h"
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "sensor.h"

namespace Kratos
{
///@addtogroup SystemIdentificationApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(DIGITAL_TWIN_APPLICATION) MeasurementResidualResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = AdjointResponseFunction;

    KRATOS_CLASS_POINTER_DEFINITION(MeasurementResidualResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MeasurementResidualResponseFunction();

    /// Destructor.
    ~MeasurementResidualResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void AddSensor(Sensor::Pointer pSensor);

    void Clear();

    std::vector<Sensor::Pointer>& GetSensorsList();

    void Initialize() override;

    void InitializeSolutionStep() override;

    void FinalizeSolutionStep() override;

    void Finalize();

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

    std::string Info() const;

    void PrintInfo(std::ostream& rOStream) const;

    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<Sensor::Pointer> mpSensorsList;

    std::vector<Vector> mResponseGradientList;

    ///@}
    ///@name Private operations
    ///@{

    template<class TCalculationType, class... TArgs>
    void CalculateDerivative(
        Vector& rResponseGradient,
        const Matrix& rResidualGradient,
        TArgs&&... rArgs);

    ///@}
};

///@} // Kratos Classes

/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const MeasurementResidualResponseFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@} //Digital Twin Application group

} /* namespace Kratos.*/
