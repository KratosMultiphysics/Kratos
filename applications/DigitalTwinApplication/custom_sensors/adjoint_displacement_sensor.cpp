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
      mNodeId(0),
      mElementId(0),
      mNodeIndex(0),
      mIsNodeAvailable(false)
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

    KRATOS_ERROR_IF_NOT(rSensorSpecification.Has(SENSOR_NODE_ID))
        << "No SENSOR_NODE_ID found in the sensor specification.";

    KRATOS_ERROR_IF_NOT(rSensorSpecification.Has(SENSOR_DIRECTION))
        << "No SENSOR_DIRECTION found in the sensor specification.";

    KRATOS_ERROR_IF_NOT(rSensorSpecification.Has(SENSOR_WEIGHT))
        << "No SENSOR_WEIGHT found in the sensor specification.";

    mNodeId = rSensorSpecification.GetValue(SENSOR_NODE_ID);
    mIsNodeAvailable = r_local_mesh.HasNode(mNodeId);

    KRATOS_ERROR_IF_NOT(r_data_communicator.OrReduceAll(mIsNodeAvailable))
        << "Node id " << mNodeId << " not found in " << r_model_part.FullName() << ".\n";

    mWeight = rSensorSpecification.GetValue(SENSOR_WEIGHT);
    mDirection = rSensorSpecification.GetValue(SENSOR_DIRECTION);

    mElementId = 0;
    mNodeIndex = 0;

    if (!rSensorSpecification.Has(SENSOR_ELEMENT_ID)){
        if (mIsNodeAvailable) {
            bool found_node = false;
            for (const auto& r_element : r_local_mesh.Elements()) {
                const auto& r_geometry = r_element.GetGeometry();
                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    if (static_cast<int>(r_geometry[i].Id()) == mNodeId) {
                        mElementId = r_element.Id();
                        mNodeIndex = i;
                        found_node = true;
                        break;
                    }
                }

                if (found_node) {
                    break;
                }
            }

            KRATOS_ERROR_IF(!found_node)
                << "The node " << mNodeId << " is not associated with any elements in "
                << r_model_part.FullName() << ".\n";
        }

        rSensorSpecification.SetValue(SENSOR_ELEMENT_ID, r_data_communicator.SumAll(mElementId));
        rSensorSpecification.SetValue(SENSOR_ELEMENT_NODE_INDEX, r_data_communicator.SumAll(mNodeIndex));
    } else {
        mElementId = rSensorSpecification.GetValue(SENSOR_ELEMENT_ID);
        mNodeIndex = rSensorSpecification.GetValue(SENSOR_ELEMENT_NODE_INDEX);
    }

    KRATOS_CATCH("");
}

double AdjointDisplacementSensor::CalculateValue(ModelPart& rModelPart)
{
    double displacement = 0.0;
    if (mIsNodeAvailable) {
        const auto& r_node = rModelPart.GetNode(mNodeId);
        displacement = inner_prod(r_node.FastGetSolutionStepValue(DISPLACEMENT), mDirection) / mWeight;
    }
    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(displacement);
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
        const IndexType block_start_index = mNodeIndex * block_size;

        rResponseGradient[block_start_index] = mDirection[0];
        rResponseGradient[block_start_index + 1] = mDirection[1];
        rResponseGradient[block_start_index + 2] = mDirection[2];
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
