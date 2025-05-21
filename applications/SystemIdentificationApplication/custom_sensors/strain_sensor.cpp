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
#include "includes/kratos_components.h"
#include "utilities/brute_force_point_locator.h"

// Application includes
#include "custom_utilities/sensor_utils.h"
#include "system_identification_application_variables.h"

// Include base h
#include "strain_sensor.h"

namespace Kratos {

/// Constructor.
StrainSensor::StrainSensor(
    const std::string& rName,
    Node::Pointer pNode,
    const Variable<Matrix>& rStrainVariable,
    const StrainType& rStrainType,
    const Element& rElement,
    const double Weight)
    : BaseType(rName, pNode, Weight),
      mElementId(rElement.Id()),
      mStrainType(rStrainType),
      mrStrainVariable(rStrainVariable)
{
    KRATOS_ERROR_IF_NOT(rElement.GetGeometry().IsInside(*(this->GetNode()), mLocalPoint))
        << "The point " << this->GetNode()->Coordinates() << " is not inside or on the boundary of the geometry of element with id "
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
    } else {
        KRATOS_ERROR << "Unsupported working space dimension = "
                     << rElement.GetGeometry().WorkingSpaceDimension()
                     << " in element with id = " << rElement.Id() << ".";
    }

    this->GetNode()->SetValue(SENSOR_ELEMENT_ID, static_cast<int>(mElementId));
}

Sensor::Pointer StrainSensor::Create(
    ModelPart& rDomainModelPart,
    ModelPart& rSensorModelPart,
    const IndexType Id,
    Parameters SensorParameters)
{
    KRATOS_TRY

    SensorParameters.ValidateAndAssignDefaults(StrainSensor::GetDefaultParameters());

    const auto& location = SensorParameters["location"].GetVector();
    KRATOS_ERROR_IF_NOT(location.size() == 3)
        << "Location of the sensor \"" << SensorParameters["name"].GetString()
        << "\" should have 3 components. [ location = " << location << " ].\n";

    Point loc(location[0], location[1], location[2]);

    Vector dummy_shape_functions;

    const auto element_id = BruteForcePointLocator(rDomainModelPart).FindElement(loc, dummy_shape_functions);
    const auto& r_element = rDomainModelPart.GetElement(element_id);

    auto p_node = rSensorModelPart.CreateNewNode(Id, location[0], location[1], location[2]);

    const auto& strain_type_str = SensorParameters["strain_type"].GetString();

    StrainType strain_type;

    if (strain_type_str == "strain_xx") {
        strain_type = StrainType::STRAIN_XX;
    } else if (strain_type_str == "strain_yy") {
        strain_type = StrainType::STRAIN_YY;
    } else if (strain_type_str == "strain_zz") {
        strain_type = StrainType::STRAIN_ZZ;
    } else if (strain_type_str == "strain_xy") {
        strain_type = StrainType::STRAIN_XY;
    } else if (strain_type_str == "strain_xz") {
        strain_type = StrainType::STRAIN_XZ;
    } else if (strain_type_str == "strain_yz") {
        strain_type = StrainType::STRAIN_YZ;
    } else {
        KRATOS_ERROR << "Unsupported strain type = \""
                     << strain_type_str << "\". Followings are supported:"
                     << "\n\tstrain_xx"
                     << "\n\tstrain_yy"
                     << "\n\tstrain_zz"
                     << "\n\tstrain_xy"
                     << "\n\tstrain_xz"
                     << "\n\tstrain_yz";
    }

    auto p_sensor = Kratos::make_shared<StrainSensor>(
        SensorParameters["name"].GetString(),
        p_node,
        KratosComponents<Variable<Matrix>>::Get(SensorParameters["strain_variable"].GetString()),
        strain_type,
        r_element,
        SensorParameters["weight"].GetDouble()
    );

    SensorUtils::ReadVariableData(p_sensor->GetNode()->GetData(), SensorParameters["variable_data"]);

    return p_sensor;

    KRATOS_CATCH("");
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

Parameters StrainSensor::GetSensorParameters() const
{
    auto parameters = BaseType::GetSensorParameters();
    parameters.AddString("type", "strain_sensor");

    switch (mStrainType) {
        case StrainType::STRAIN_XX:
            parameters.AddString("strain_type", "strain_xx");
            break;
        case StrainType::STRAIN_YY:
            parameters.AddString("strain_type", "strain_yy");
            break;
        case StrainType::STRAIN_ZZ:
            parameters.AddString("strain_type", "strain_zz");
            break;
        case StrainType::STRAIN_XY:
            parameters.AddString("strain_type", "strain_xy");
            break;
        case StrainType::STRAIN_XZ:
            parameters.AddString("strain_type", "strain_xz");
            break;
        case StrainType::STRAIN_YZ:
            parameters.AddString("strain_type", "strain_yz");
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

        Element::DofsVectorType elemental_dofs;
        r_element.GetDofList(elemental_dofs, rProcessInfo);

        const double delta = rProcessInfo[PERTURBATION_SIZE];
        auto& r_geometry = r_element.GetGeometry();

        // now only keep the dofs of the first node since we only need the variable
        // type to be checked.
        const IndexType block_size = elemental_dofs.size() / r_geometry.size();
        elemental_dofs.erase(elemental_dofs.begin() + block_size, elemental_dofs.end());

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
            auto& r_node = r_geometry[i_node];

            for (const auto& p_dof : elemental_dofs) {
                if (p_dof->GetVariable() == ADJOINT_DISPLACEMENT_X) {
                    rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                        delta, DISPLACEMENT_X, r_node, r_element,
                        perturbed_strains, ref_strains, rProcessInfo);
                } else if (p_dof->GetVariable() == ADJOINT_DISPLACEMENT_Y) {
                    rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                        delta, DISPLACEMENT_Y, r_node, r_element,
                        perturbed_strains, ref_strains, rProcessInfo);
                } else if (p_dof->GetVariable() == ADJOINT_DISPLACEMENT_Z) {
                    rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                        delta, DISPLACEMENT_Z, r_node, r_element,
                        perturbed_strains, ref_strains, rProcessInfo);
                } else if (p_dof->GetVariable() == ADJOINT_ROTATION_X) {
                    rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                        delta, ROTATION_X, r_node, r_element,
                        perturbed_strains, ref_strains, rProcessInfo);
                } else if (p_dof->GetVariable() == ADJOINT_ROTATION_Y) {
                    rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                        delta, ROTATION_Y, r_node, r_element,
                        perturbed_strains, ref_strains, rProcessInfo);
                } else if (p_dof->GetVariable() == ADJOINT_ROTATION_Z) {
                    rResponseGradient[local_index++] = this->CalculateStrainDirectionalSensitivity(
                        delta, ROTATION_Z, r_node, r_element,
                        perturbed_strains, ref_strains, rProcessInfo);
                } else {
                    KRATOS_ERROR << "Unsupported dof " << p_dof->GetVariable().Name();
                }
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
    rOStream << "    Element Id: " << mElementId << std::endl;
    switch (mStrainType) {
        case StrainType::STRAIN_XX:
            rOStream << "    Direction: STRAIN_XX" << std::endl;
            break;
        case StrainType::STRAIN_YY:
            rOStream << "    Direction: STRAIN_YY" << std::endl;
            break;
        case StrainType::STRAIN_ZZ:
            rOStream << "    Direction: STRAIN_ZZ" << std::endl;
            break;
        case StrainType::STRAIN_XY:
            rOStream << "    Direction: STRAIN_XY" << std::endl;
            break;
        case StrainType::STRAIN_XZ:
            rOStream << "    Direction: STRAIN_XZ" << std::endl;
            break;
        case StrainType::STRAIN_YZ:
            rOStream << "    Direction: STRAIN_YZ" << std::endl;
            break;
    }
    Sensor::PrintData(rOStream);
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
