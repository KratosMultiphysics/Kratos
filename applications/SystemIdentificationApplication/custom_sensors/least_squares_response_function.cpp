//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya,
//                   Ihar Antonau
//                   Jonas Hillenbrand
//

// System includes
#include <iostream>
#include <cmath>
#include <sstream>

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "least_squares_response_function.h"

namespace Kratos
{

namespace LeastSquaresResponseFunctionUtilities
{

struct GradientCalculation {
    template <class TEntityType>
    static inline void Calculate(
        Sensor& rSensor,
        Vector& rResponseGradient,
        const Matrix& rResidualGradient,
        const TEntityType& rEntity,
        const ProcessInfo& rProcessInfo)
    {
        rSensor.CalculateGradient(rEntity, rResidualGradient, rResponseGradient, rProcessInfo);
    }
};

struct FirstDerivativesGradientCalculation {
    template <class TEntityType>
    static inline void Calculate(
        Sensor& rSensor,
        Vector& rResponseGradient,
        const Matrix& rResidualGradient,
        const TEntityType& rEntity,
        const ProcessInfo& rProcessInfo)
    {
        rSensor.CalculateFirstDerivativesGradient(rEntity, rResidualGradient, rResponseGradient, rProcessInfo);
    }
};

struct SecondDerivativesGradientCalculation {
    template <class TEntityType>
    static inline void Calculate(
        Sensor& rSensor,
        Vector& rResponseGradient,
        const Matrix& rResidualGradient,
        const TEntityType& rEntity,
        const ProcessInfo& rProcessInfo)
    {
        rSensor.CalculateSecondDerivativesGradient(rEntity, rResidualGradient, rResponseGradient, rProcessInfo);
    }
};

struct PartialSensitivity
{
    template <class TEntityType, class TVariableType>
    static inline void Calculate(
        Sensor& rSensor,
        Vector& rSensitivityGradient,
        const Matrix& rSensitivityMatrix,
        TEntityType& rEntity,
        const TVariableType& rVariable,
        const ProcessInfo& rProcessInfo)
    {
        rSensor.CalculatePartialSensitivity(rEntity, rVariable, rSensitivityMatrix, rSensitivityGradient, rProcessInfo);
    }
};

} // namespace LeastSquaresResponseFunctionUtilities

LeastSquaresResponseFunction::LeastSquaresResponseFunction()
{
    mResponseGradientList.resize(ParallelUtilities::GetNumThreads());
}

void LeastSquaresResponseFunction::AddSensor(Sensor::Pointer pSensor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(pSensor->Has(SENSOR_MEASURED_VALUE))
        << pSensor->GetName() << " does not have the SENSOR_MEASURED_VALUE defined.\n";

    mpSensorsList.push_back(pSensor);

    KRATOS_CATCH("");
}

void LeastSquaresResponseFunction::Clear()
{
    mpSensorsList.clear();
}

std::vector<Sensor::Pointer>& LeastSquaresResponseFunction::GetSensorsList()
{
    return mpSensorsList;
}

void LeastSquaresResponseFunction::Initialize()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->Initialize();
    }
}

void LeastSquaresResponseFunction::InitializeSolutionStep()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->InitializeSolutionStep();
    }
}

void LeastSquaresResponseFunction::FinalizeSolutionStep()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->FinalizeSolutionStep();
    }
}

void LeastSquaresResponseFunction::Finalize()
{

}

double LeastSquaresResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY

    double sum = 0.0;

    for (auto& p_sensor : mpSensorsList) {
        const double sensor_value = p_sensor->CalculateValue(rModelPart);
        const double current_sensor_error = sensor_value - p_sensor->GetValue(SENSOR_MEASURED_VALUE);

        // TODO: Why do we store the sensor_value?
        p_sensor->SetSensorValue(sensor_value);

        sum += 0.5 * std::pow(current_sensor_error, 2) * p_sensor->GetWeight();
    }

    return sum;

    KRATOS_CATCH("");
}

template<class TCalculationType, class... TArgs>
void LeastSquaresResponseFunction::CalculateDerivative(
    Vector& rResponseGradient,
    const Matrix& rResidualGradient,
    TArgs&&... rArgs)
{
    KRATOS_TRY

    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    auto& local_sensor_response_gradient = mResponseGradientList[OpenMPUtils::ThisThread()];

    for (auto& p_sensor : mpSensorsList) {
        TCalculationType::Calculate(*p_sensor, local_sensor_response_gradient, rResidualGradient, rArgs...);
        noalias(rResponseGradient) += p_sensor->GetValue(SENSOR_ERROR) * p_sensor->GetWeight() * local_sensor_response_gradient;
    }

    KRATOS_CATCH("");
}

void LeastSquaresResponseFunction::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::GradientCalculation>(rResponseGradient, rResidualGradient, rAdjointElement, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::GradientCalculation>(rResponseGradient, rResidualGradient, rAdjointCondition, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualFirstDerivativesGradient,
    Vector& rResponseFirstDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::FirstDerivativesGradientCalculation>(rResponseFirstDerivativesGradient, rResidualFirstDerivativesGradient, rAdjointElement, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualFirstDerivativesGradient,
    Vector& rResponseFirstDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::FirstDerivativesGradientCalculation>(rResponseFirstDerivativesGradient, rResidualFirstDerivativesGradient, rAdjointCondition, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualSecondDerivativesGradient,
    Vector& rResponseSecondDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::SecondDerivativesGradientCalculation>(rResponseSecondDerivativesGradient, rResidualSecondDerivativesGradient, rAdjointElement, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualSecondDerivativesGradient,
    Vector& rResponseSecondDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::SecondDerivativesGradientCalculation>(rResponseSecondDerivativesGradient, rResidualSecondDerivativesGradient, rAdjointCondition, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointElement, rVariable, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointCondition, rVariable, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointElement, rVariable, rProcessInfo);
}

void LeastSquaresResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<LeastSquaresResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointCondition, rVariable, rProcessInfo);
}

std::string LeastSquaresResponseFunction::Info() const
{
    std::stringstream msg;
    msg << "LeastSquaresResponseFunction with " << mpSensorsList.size() << " sensors.";
    return msg.str();
}

void LeastSquaresResponseFunction::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info() << std::endl;
}

void LeastSquaresResponseFunction::PrintData(std::ostream& rOStream) const
{
    PrintInfo(rOStream);
    for (const auto& p_sensor : mpSensorsList) {
        rOStream << "\t" << p_sensor->GetName() << std::endl;
    }
}

} /* namespace Kratos.*/
