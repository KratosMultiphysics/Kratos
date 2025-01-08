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

// External includes

// Project includes

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "strain_sensor.h"

namespace Kratos {

/// Constructor.
StrainSensor::StrainSensor(
    const std::string& rName,
    const Point& rLocation,
    const Variable<Matrix>& rStrainVariable,
    const StrainType& rStrainType,
    const Element& rElement,
    const double Weight)
    : BaseType(rName, rLocation, Weight),
      mElementId(rElement.Id()),
      mStrainType(rStrainType),
      mrStrainVariable(rStrainVariable)
{
    KRATOS_ERROR_IF_NOT(rElement.GetGeometry().IsInside(this->GetLocation(), mLocalPoint))
        << "The point " << this->GetLocation() << " is not inside or on the boundary of the geometry of element with id "
        << mElementId << ".";

    // if the element is of type 2d
    if (rElement.GetGeometry().WorkingSpaceDimension() == 2) {
        switch (mStrainType) {
            case StrainType::STRAIN_XX:
                mStrainIndex = 0;
                break;
            case StrainType::STRAIN_YY:
                mStrainIndex = 3;
                break;
            case StrainType::STRAIN_XY:
                mStrainIndex = 1;
                break;
            default:
                KRATOS_ERROR << "The element with id = " << rElement.Id() << " is of 2d type, hence only xx, yy, xy strains are allowed.";
        }
    } else if (rElement.GetGeometry().WorkingSpaceDimension() == 3) {
        mStrainIndex = mStrainType;
    }

    this->SetValue(SENSOR_ELEMENT_ID, static_cast<int>(mElementId));
}

Parameters StrainSensor::GetDefaultParameters()
{
    return Parameters(R"(
    {
        "type"           : "strain_sensor",
        "name"           : "",
        "value"          : 0.0,
        "location"       : [0.0, 0.0, 0.0],
        "strain_type"    : "strain_xx",
        "strain_variable": "SHELL_STRAIN",
        "weight"         : 0.0,
        "variable_data": {}
    })" );
}

const Parameters StrainSensor::GetSensorParameters() const
{
    Parameters parameters = Parameters(R"(
    {
        "type"           : "strain_sensor",
        "name"           : "",
        "value"          : 0.0,
        "location"       : [0.0, 0.0, 0.0],
        "strain_type"    : "strain_xx",
        "strain_variable": "SHELL_STRAIN",
        "weight"         : 0.0
    })" );
    parameters["name"].SetString(this->GetName());
    parameters["value"].SetDouble(this->GetSensorValue());
    parameters["location"].SetVector(this->GetLocation());
    parameters["weight"].SetDouble(this->GetWeight());

    switch (mStrainType) {
        case StrainType::STRAIN_XX:
            parameters["strain_type"].SetString("strain_xx");
            break;
        case StrainType::STRAIN_YY:
            parameters["strain_type"].SetString("strain_yy");
            break;
        case StrainType::STRAIN_ZZ:
            parameters["strain_type"].SetString("strain_zz");
            break;
        case StrainType::STRAIN_XY:
            parameters["strain_type"].SetString("strain_xy");
            break;
        case StrainType::STRAIN_XZ:
            parameters["strain_type"].SetString("strain_xz");
            break;
        case StrainType::STRAIN_YZ:
            parameters["strain_type"].SetString("strain_yz");
            break;
    };

    return parameters;
}

double StrainSensor::CalculateValue(ModelPart& rModelPart)
{
    double directional_strain = 0.0;
    if (rModelPart.HasElement(mElementId)) {
        auto& r_element = rModelPart.GetElement(mElementId);

        std::vector<Matrix> strains;
        r_element.CalculateOnIntegrationPoints(mrStrainVariable, strains, rModelPart.GetProcessInfo());

        for (const auto& strain : strains) {
            KRATOS_ERROR_IF(strain.data().size() <= mStrainIndex)
                << "The size of the strain " << strain.data().size() << " does not contain the index = " << mStrainIndex << ".";
            directional_strain += *(strain.data().begin() + mStrainIndex);
        }

        directional_strain /= strains.size();
    }

    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(directional_strain);
}

void StrainSensor::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();

    if (rAdjointElement.Id() == mElementId) {
        auto& r_element = const_cast<Element&>(rAdjointElement);
        std::vector<Matrix> perturbed_strains, ref_strains;
        r_element.CalculateOnIntegrationPoints(mrStrainVariable, ref_strains, rProcessInfo);

        const double delta = rProcessInfo[PERTURBATION_SIZE];
        auto& r_geometry = r_element.GetGeometry();

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
            auto& r_node = r_geometry[i_node];

            if (r_node.HasDofFor(ADJOINT_DISPLACEMENT_X)) {
                rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                    delta, DISPLACEMENT_X, r_node, r_element,
                    perturbed_strains, ref_strains, rProcessInfo);
            }

            if (r_node.HasDofFor(ADJOINT_DISPLACEMENT_Y)) {
                rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                    delta, DISPLACEMENT_Y, r_node, r_element,
                    perturbed_strains, ref_strains, rProcessInfo);
            }

            if (r_node.HasDofFor(ADJOINT_DISPLACEMENT_Z)) {
                rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                    delta, DISPLACEMENT_Z, r_node, r_element,
                    perturbed_strains, ref_strains, rProcessInfo);
            }

            if (r_node.HasDofFor(ADJOINT_ROTATION_X)) {
                rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                    delta, ROTATION_X, r_node, r_element,
                    perturbed_strains, ref_strains, rProcessInfo);
            }

            if (r_node.HasDofFor(ADJOINT_ROTATION_Y)) {
                rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                    delta, ROTATION_Y, r_node, r_element,
                    perturbed_strains, ref_strains, rProcessInfo);
            }

            if (r_node.HasDofFor(ADJOINT_ROTATION_Z)) {
                rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                    delta, ROTATION_Z, r_node, r_element,
                    perturbed_strains, ref_strains, rProcessInfo);
            }
        }
    }

    KRATOS_CATCH("");
}

void StrainSensor::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void StrainSensor::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void StrainSensor::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void StrainSensor::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void StrainSensor::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void StrainSensor::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void StrainSensor::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void StrainSensor::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void StrainSensor::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

std::string StrainSensor::Info() const
{
    std::stringstream msg;
    msg << "StrainSensor " << this->GetName();
    return msg.str();
}

void StrainSensor::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info() << std::endl;
}

void StrainSensor::PrintData(std::ostream& rOStream) const
{
    PrintInfo(rOStream);
    rOStream << "    Location: " << this->GetLocation() << std::endl;
    rOStream << "    Value: " << this->GetSensorValue() << std::endl;
    rOStream << "    Weight: " << this->GetWeight() << std::endl;
    rOStream << "    Element Id: " << mElementId << std::endl;
    switch (mStrainType) {
        case StrainType::STRAIN_XX:
            rOStream << "    Direction: STRAIN_XX";
            break;
        case StrainType::STRAIN_YY:
            rOStream << "    Direction: STRAIN_YY";
            break;
        case StrainType::STRAIN_ZZ:
            rOStream << "    Direction: STRAIN_ZZ";
            break;
        case StrainType::STRAIN_XY:
            rOStream << "    Direction: STRAIN_XY";
            break;
        case StrainType::STRAIN_XZ:
            rOStream << "    Direction: STRAIN_XZ";
            break;
        case StrainType::STRAIN_YZ:
            rOStream << "    Direction: STRAIN_YZ";
            break;
    }
    DataValueContainer::PrintData(rOStream);
}

void StrainSensor::SetVectorToZero(
    Vector& rVector,
    const IndexType Size)
{
    if (rVector.size() != Size) {
        rVector.resize(Size, false);
    }

    rVector.clear();
}

double StrainSensor::CalculateStrainDirectionalSensitivity(
    const double Perturbation,
    const Variable<double>& rPerturbationVariable,
    ModelPart::NodeType& rNode,
    ModelPart::ElementType& rElement,
    std::vector<Matrix>& rPerturbedStrains,
    const std::vector<Matrix>& rRefStrains,
    const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    rNode.FastGetSolutionStepValue(rPerturbationVariable) += Perturbation;
    rElement.CalculateOnIntegrationPoints(mrStrainVariable, rPerturbedStrains, rProcessInfo);
    rNode.FastGetSolutionStepValue(rPerturbationVariable) -= Perturbation;

    double strain_sensitivity = 0.0;
    for (IndexType j = 0; j < rPerturbedStrains.size(); ++j) {
        strain_sensitivity += (*(rPerturbedStrains[j].data().begin() + mStrainIndex) -
                               *(rRefStrains[j].data().begin() + mStrainIndex)) /
                              Perturbation;
    }

    return strain_sensitivity / rPerturbedStrains.size();

    KRATOS_CATCH("");
}

}; // namespace Kratos
