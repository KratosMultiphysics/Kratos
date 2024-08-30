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

// System includes
#include <cmath>
#include <sstream>

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "measurement_residual_response_function.h"

namespace Kratos
{

namespace MeasurementResidualResponseFunctionUtilities
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

} // namespace MeasurementResidualResponseFunctionUtilities

MeasurementResidualResponseFunction::MeasurementResidualResponseFunction(const double PCoefficient)
    : mPCoefficient(PCoefficient)
    {
    mResponseGradientList.resize(ParallelUtilities::GetNumThreads());
}

void MeasurementResidualResponseFunction::AddSensor(Sensor::Pointer pSensor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(pSensor->Has(SENSOR_MEASURED_VALUE))
        << pSensor->GetName() << " does not have the SENSOR_MEASURED_VALUE defined.\n";

    mpSensorsList.push_back(pSensor);

    KRATOS_CATCH("");
}

void MeasurementResidualResponseFunction::Clear()
{
    mpSensorsList.clear();
}

std::vector<Sensor::Pointer>& MeasurementResidualResponseFunction::GetSensorsList()
{
    return mpSensorsList;
}

void MeasurementResidualResponseFunction::Initialize()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->Initialize();
    }
}

void MeasurementResidualResponseFunction::InitializeSolutionStep()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->InitializeSolutionStep();
    }
}

void MeasurementResidualResponseFunction::FinalizeSolutionStep()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->FinalizeSolutionStep();
    }
}

void MeasurementResidualResponseFunction::Finalize()
{

}

double MeasurementResidualResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY

    double sum_B_p = 0.0;
    double value = 0.0;
    for (auto& p_sensor : mpSensorsList) {
        const double sensor_value = p_sensor->CalculateValue(rModelPart);
        p_sensor->SetSensorValue(sensor_value);
        const double current_sensor_error = sensor_value - p_sensor->GetValue(SENSOR_MEASURED_VALUE);
        p_sensor->SetValue(SENSOR_ERROR, current_sensor_error);
        sum_B_p += ( std::pow( 0.5 * pow(current_sensor_error, 2) * p_sensor->GetWeight(), mPCoefficient ) );
        value += p_sensor->GetWeight() * std::pow(sensor_value - p_sensor->GetValue(SENSOR_MEASURED_VALUE), 2) * 0.5;
    }
    mC1 = std::pow( sum_B_p, 1/mPCoefficient - 1 ) / std::pow(2, mPCoefficient - 1);  

    KRATOS_WATCH(value);
    KRATOS_WATCH(std::pow(sum_B_p, 1 / mPCoefficient));

    return std::pow(sum_B_p, 1 / mPCoefficient);


    KRATOS_CATCH("");
}

template<class TCalculationType, class... TArgs>
void MeasurementResidualResponseFunction::CalculateDerivative(
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
        noalias(rResponseGradient) += mC1 * (std::pow(p_sensor->GetWeight(), mPCoefficient) * std::pow(p_sensor->GetValue(SENSOR_ERROR), mPCoefficient * 2 - 1 ) ) * local_sensor_response_gradient ;      
    }

    KRATOS_CATCH("");
}

void MeasurementResidualResponseFunction::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::GradientCalculation>(rResponseGradient, rResidualGradient, rAdjointElement, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::GradientCalculation>(rResponseGradient, rResidualGradient, rAdjointCondition, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualFirstDerivativesGradient,
    Vector& rResponseFirstDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::FirstDerivativesGradientCalculation>(rResponseFirstDerivativesGradient, rResidualFirstDerivativesGradient, rAdjointElement, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualFirstDerivativesGradient,
    Vector& rResponseFirstDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::FirstDerivativesGradientCalculation>(rResponseFirstDerivativesGradient, rResidualFirstDerivativesGradient, rAdjointCondition, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualSecondDerivativesGradient,
    Vector& rResponseSecondDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::SecondDerivativesGradientCalculation>(rResponseSecondDerivativesGradient, rResidualSecondDerivativesGradient, rAdjointElement, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualSecondDerivativesGradient,
    Vector& rResponseSecondDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::SecondDerivativesGradientCalculation>(rResponseSecondDerivativesGradient, rResidualSecondDerivativesGradient, rAdjointCondition, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointElement, rVariable, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointCondition, rVariable, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointElement, rVariable, rProcessInfo);
}

void MeasurementResidualResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointCondition, rVariable, rProcessInfo);
}

std::string MeasurementResidualResponseFunction::Info() const
{
    std::stringstream msg;
    msg << "MeasurementResidualResponseFunction with " << mpSensorsList.size() << " sensors.";
    return msg.str();
}

void MeasurementResidualResponseFunction::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info() << std::endl;
}

void MeasurementResidualResponseFunction::PrintData(std::ostream& rOStream) const
{
    PrintInfo(rOStream);
    for (const auto& p_sensor : mpSensorsList) {
        rOStream << "\t" << p_sensor->GetName() << std::endl;
    }
}

} /* namespace Kratos.*/
