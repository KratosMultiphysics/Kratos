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

// System includes
#include <cmath>
#include <sstream>

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "measurement_residual_p_norm_response_function.h"

namespace Kratos
{

namespace MeasurementResidualPNormResponseFunctionUtilities
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

} // namespace MeasurementResidualPNormResponseFunctionUtilities

MeasurementResidualPNormResponseFunction::MeasurementResidualPNormResponseFunction(const double PCoefficient)
    : mPCoefficient(PCoefficient)
{
    mResponseGradientList.resize(ParallelUtilities::GetNumThreads());
}

void MeasurementResidualPNormResponseFunction::AddSensor(Sensor::Pointer pSensor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(pSensor->Has(SENSOR_MEASURED_VALUE))
        << pSensor->GetName() << " does not have the SENSOR_MEASURED_VALUE defined.\n";

    mpSensorsList.push_back(pSensor);

    KRATOS_CATCH("");
}

void MeasurementResidualPNormResponseFunction::Clear()
{
    mpSensorsList.clear();
}

std::vector<Sensor::Pointer>& MeasurementResidualPNormResponseFunction::GetSensorsList()
{
    return mpSensorsList;
}

void MeasurementResidualPNormResponseFunction::Initialize()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->Initialize();
    }
}

void MeasurementResidualPNormResponseFunction::InitializeSolutionStep()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->InitializeSolutionStep();
    }
}

void MeasurementResidualPNormResponseFunction::FinalizeSolutionStep()
{
    for (auto& p_sensor : mpSensorsList) {
        p_sensor->FinalizeSolutionStep();
    }
}

void MeasurementResidualPNormResponseFunction::Finalize()
{

}

double MeasurementResidualPNormResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY

    double value = 0.0;
    for (auto& p_sensor : mpSensorsList) {
        const double sensor_value = p_sensor->CalculateValue(rModelPart);
        p_sensor->SetSensorValue(sensor_value);
        const double current_sensor_error_square  = std::pow(sensor_value - p_sensor->GetValue(SENSOR_MEASURED_VALUE), 2);
        p_sensor->SetValue(SENSOR_ERROR_SQUARE, current_sensor_error_square);
        value += std::pow(p_sensor->GetValue(SENSOR_ERROR_SQUARE), mPCoefficient);
    }

    return std::pow(value, 1 / mPCoefficient);

    KRATOS_CATCH("");
}

template<class TCalculationType, class... TArgs>
void MeasurementResidualPNormResponseFunction::CalculateDerivative(
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
    double temp = 0.0;
    for (auto& p_sensor : mpSensorsList) {
        temp += ( std::pow( p_sensor->GetValue(SENSOR_ERROR_SQUARE), mPCoefficient ) );
    }
    double c1 = 1 / mPCoefficient * std::pow( temp, 1/mPCoefficient - 1 );

    for (auto& p_sensor : mpSensorsList) {
        TCalculationType::Calculate(*p_sensor, local_sensor_response_gradient, rResidualGradient, rArgs...);
        double error = p_sensor->GetSensorValue() - p_sensor->GetValue(SENSOR_MEASURED_VALUE);
        noalias(rResponseGradient) += c1 * local_sensor_response_gradient * 2.0 * error * mPCoefficient * std::pow( p_sensor->GetValue(SENSOR_ERROR_SQUARE), mPCoefficient - 1 ) ;
    }

    KRATOS_CATCH("");
}

void MeasurementResidualPNormResponseFunction::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::GradientCalculation>(rResponseGradient, rResidualGradient, rAdjointElement, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::GradientCalculation>(rResponseGradient, rResidualGradient, rAdjointCondition, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualFirstDerivativesGradient,
    Vector& rResponseFirstDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::FirstDerivativesGradientCalculation>(rResponseFirstDerivativesGradient, rResidualFirstDerivativesGradient, rAdjointElement, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualFirstDerivativesGradient,
    Vector& rResponseFirstDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::FirstDerivativesGradientCalculation>(rResponseFirstDerivativesGradient, rResidualFirstDerivativesGradient, rAdjointCondition, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualSecondDerivativesGradient,
    Vector& rResponseSecondDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::SecondDerivativesGradientCalculation>(rResponseSecondDerivativesGradient, rResidualSecondDerivativesGradient, rAdjointElement, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualSecondDerivativesGradient,
    Vector& rResponseSecondDerivativesGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::SecondDerivativesGradientCalculation>(rResponseSecondDerivativesGradient, rResidualSecondDerivativesGradient, rAdjointCondition, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointElement, rVariable, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointCondition, rVariable, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointElement, rVariable, rProcessInfo);
}

void MeasurementResidualPNormResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDerivative<MeasurementResidualPNormResponseFunctionUtilities::PartialSensitivity>(rSensitivityGradient, rSensitivityMatrix, rAdjointCondition, rVariable, rProcessInfo);
}

std::string MeasurementResidualPNormResponseFunction::Info() const
{
    std::stringstream msg;
    msg << "MeasurementResidualPNormResponseFunction with " << mpSensorsList.size() << " sensors.";
    return msg.str();
}

void MeasurementResidualPNormResponseFunction::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info() << std::endl;
}

void MeasurementResidualPNormResponseFunction::PrintData(std::ostream& rOStream) const
{
    PrintInfo(rOStream);
    for (const auto& p_sensor : mpSensorsList) {
        rOStream << "\t" << p_sensor->GetName() << std::endl;
    }
}

} /* namespace Kratos.*/
