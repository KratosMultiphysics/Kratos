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

// External includes

// Project includes

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "adjoint_displacement_sensor.h"

namespace Kratos {

/// Constructor.
AdjointDisplacementSensor::AdjointDisplacementSensor(
    Model& rModel,
    Parameters SensorSettings)
    : mrModel(rModel),
      mElementId(0),
      mIsElementAvailable(false)
{
    KRATOS_TRY;

    Parameters defaut_parameters(R"({
            "model_part_name"        : "PLEASE_SPECIFY_A_MODEL_PART_NAME",
            "perturbation_size"      : 1e-8,
            "adapt_perturbation_size": false
        })");

    SensorSettings.ValidateAndAssignDefaults(defaut_parameters);

    mModelPartName = SensorSettings["model_part_name"].GetString();

    auto& r_process_info = mrModel.GetModelPart(mModelPartName).GetProcessInfo();
    r_process_info[PERTURBATION_SIZE] = SensorSettings["perturbation_size"].GetDouble();
    r_process_info[ADAPT_PERTURBATION_SIZE] = SensorSettings["adapt_perturbation_size"].GetBool();

    KRATOS_CATCH("");
}

void AdjointDisplacementSensor::SetSensorSpecification(SensorSpecification& rSensorSpecification)
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_communicator = r_model_part.GetCommunicator();
    const auto& r_local_mesh = r_communicator.LocalMesh();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();

    KRATOS_ERROR_IF_NOT(rSensorSpecification.Has(SENSOR_ELEMENT_ID))
        << "No SENSOR_ELEMENT_ID found in the sensor specification.";

    KRATOS_ERROR_IF_NOT(rSensorSpecification.Has(SENSOR_DIRECTION))
        << "No SENSOR_DIRECTION found in the sensor specification.";

    KRATOS_ERROR_IF_NOT(rSensorSpecification.Has(SENSOR_WEIGHT))
        << "No SENSOR_WEIGHT found in the sensor specification.";

    mElementId = rSensorSpecification.GetValue(SENSOR_ELEMENT_ID);
    mDirection = rSensorSpecification.GetValue(SENSOR_DIRECTION);
    mWeight = rSensorSpecification.GetValue(SENSOR_WEIGHT);

    mIsElementAvailable = r_local_mesh.HasElement(mElementId);

    KRATOS_ERROR_IF_NOT(r_data_communicator.OrReduceAll(mIsElementAvailable))
        << "Element id " << mElementId << " not found in " << r_model_part.FullName() << ".\n";

    if (mIsElementAvailable) {
        const auto& r_element = r_local_mesh.GetElement(mElementId);
        r_element.GetGeometry().PointLocalCoordinates(mLocalCoordinates, rSensorSpecification.GetLocation());
    }

    KRATOS_CATCH("");
}

double AdjointDisplacementSensor::CalculateValue(ModelPart& rModelPart)
{
    double directional_displacement = 0.0;
    if (mIsElementAvailable) {
        const auto& r_element = rModelPart.GetElement(mElementId);
        const auto& r_geometry = r_element.GetGeometry();

        Vector Ns;
        r_geometry.ShapeFunctionsValues(Ns, mLocalCoordinates);
        array_1d<double, 3> displacement = ZeroVector(3);
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            displacement += r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT) * Ns[i];
        }

        directional_displacement = inner_prod(displacement, mDirection) / mWeight;
    }
    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(directional_displacement);
}

void AdjointDisplacementSensor::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    if (static_cast<int>(rAdjointElement.Id()) == mElementId) {
        const auto& r_geometry = rAdjointElement.GetGeometry();
        const IndexType block_size = rResidualGradient.size1() / r_geometry.size();

        Vector Ns;
        r_geometry.ShapeFunctionsValues(Ns, mLocalCoordinates);
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            rResponseGradient[i * block_size] = Ns[i] * mDirection[0];
            rResponseGradient[i * block_size + 1] = Ns[i] * mDirection[1];
            rResponseGradient[i * block_size + 2] = Ns[i] * mDirection[2];
        }
    }

    rResponseGradient /= mWeight;

    KRATOS_CATCH("");
}

void AdjointDisplacementSensor::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void AdjointDisplacementSensor::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void AdjointDisplacementSensor::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void AdjointDisplacementSensor::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void AdjointDisplacementSensor::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void AdjointDisplacementSensor::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void AdjointDisplacementSensor::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void AdjointDisplacementSensor::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void AdjointDisplacementSensor::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void AdjointDisplacementSensor::SetVectorToZero(
    Vector& rVector,
    const IndexType Size)
{
    if (rVector.size() != Size) {
        rVector.resize(Size, false);
    }

    rVector.clear();
}

}; // namespace Kratos
