// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_element_strain_response_function.h"

namespace Kratos
{
AdjointElementStrainResponseFunction::AdjointElementStrainResponseFunction(
    ModelPart& rModelPart,
    Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
{
    mElementId = ResponseSettings["element_id"].GetInt();

    mWeight = ResponseSettings["weight"].GetDouble();

    const auto& strain_type = ResponseSettings["strain_type"].GetString();
    if (strain_type == "x") {
        mStrainType = StrainType::STRAIN_X;
    } else if (strain_type == "y") {
        mStrainType = StrainType::STRAIN_Y;
    } else if (strain_type == "z") {
        mStrainType = StrainType::STRAIN_Z;
    } else {
        KRATOS_ERROR << "Unsupported strain type.";
    }
}

double AdjointElementStrainResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY;

    auto& r_element = mrModelPart.GetElement(mElementId);
    std::vector<Matrix> strains;
    r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, strains, mrModelPart.GetProcessInfo());

    double result = 0.0;
    for (const auto& strain : strains) {
        if (mStrainType == StrainType::STRAIN_X) {
            result += strain(0, 0);
        } else if (mStrainType == StrainType::STRAIN_Y) {
            result += strain(1, 1);
        } else if (mStrainType == StrainType::STRAIN_Z) {
            result += strain(2, 2);
        }
    }

    return result / (mWeight * strains.size());

    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculateGradient(
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

    if (rAdjointElement.Id() == mElementId) {
        auto& r_element = mrModelPart.GetElement(mElementId);
        std::vector<Matrix> perturbed_strains, ref_strains;
        Matrix strain_sensitivities(3, 3);
        r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, ref_strains, rProcessInfo);

        DofsVectorType dofs_of_element;
        r_element.GetDofList(dofs_of_element, rProcessInfo);
        const double delta = rProcessInfo[PERTURBATION_SIZE];

        auto& r_geometry = r_element.GetGeometry();

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
            auto& r_node = r_geometry[i_node];

            // if (r_node.HasDofFor(ADJOINT_DISPLACEMENT_X)) {
                r_node.FastGetSolutionStepValue(DISPLACEMENT_X) += delta;
                r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, perturbed_strains, rProcessInfo);
                r_node.FastGetSolutionStepValue(DISPLACEMENT_X) -= delta;
                strain_sensitivities.clear();
                for (IndexType j = 0; j < perturbed_strains.size(); ++j) {
                    noalias(strain_sensitivities) += (perturbed_strains[j] - ref_strains[j]) / delta;
                }
                if (mStrainType == StrainType::STRAIN_X) {
                    rResponseGradient[local_index++] = strain_sensitivities(0, 0);
                } else if (mStrainType == StrainType::STRAIN_Y) {
                    rResponseGradient[local_index++] = strain_sensitivities(1, 1);
                } else if (mStrainType == StrainType::STRAIN_Z) {
                    rResponseGradient[local_index++] = strain_sensitivities(2, 2);
                }
            // }
            // if (r_node.HasDofFor(ADJOINT_DISPLACEMENT_Y)) {
                r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) += delta;
                r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, perturbed_strains, rProcessInfo);
                r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) -= delta;
                strain_sensitivities.clear();
                for (IndexType j = 0; j < perturbed_strains.size(); ++j) {
                    noalias(strain_sensitivities) += (perturbed_strains[j] - ref_strains[j]) / delta;
                }
                if (mStrainType == StrainType::STRAIN_X) {
                    rResponseGradient[local_index++] = strain_sensitivities(0, 0);
                } else if (mStrainType == StrainType::STRAIN_Y) {
                    rResponseGradient[local_index++] = strain_sensitivities(1, 1);
                } else if (mStrainType == StrainType::STRAIN_Z) {
                    rResponseGradient[local_index++] = strain_sensitivities(2, 2);
                }
            // }
            // if (r_node.HasDofFor(ADJOINT_DISPLACEMENT_Z)) {
                r_node.FastGetSolutionStepValue(DISPLACEMENT_Z) += delta;
                r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, perturbed_strains, rProcessInfo);
                r_node.FastGetSolutionStepValue(DISPLACEMENT_Z) -= delta;
                strain_sensitivities.clear();
                for (IndexType j = 0; j < perturbed_strains.size(); ++j) {
                    noalias(strain_sensitivities) += (perturbed_strains[j] - ref_strains[j]) / delta;
                }
                if (mStrainType == StrainType::STRAIN_X) {
                    rResponseGradient[local_index++] = strain_sensitivities(0, 0);
                } else if (mStrainType == StrainType::STRAIN_Y) {
                    rResponseGradient[local_index++] = strain_sensitivities(1, 1);
                } else if (mStrainType == StrainType::STRAIN_Z) {
                    rResponseGradient[local_index++] = strain_sensitivities(2, 2);
                }
            // }
            // if (r_node.HasDofFor(ADJOINT_ROTATION_X)) {
                r_node.FastGetSolutionStepValue(ROTATION_X) += delta;
                r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, perturbed_strains, rProcessInfo);
                r_node.FastGetSolutionStepValue(ROTATION_X) -= delta;
                strain_sensitivities.clear();
                for (IndexType j = 0; j < perturbed_strains.size(); ++j) {
                    noalias(strain_sensitivities) += (perturbed_strains[j] - ref_strains[j]) / delta;
                }
                if (mStrainType == StrainType::STRAIN_X) {
                    rResponseGradient[local_index++] = strain_sensitivities(0, 0);
                } else if (mStrainType == StrainType::STRAIN_Y) {
                    rResponseGradient[local_index++] = strain_sensitivities(1, 1);
                } else if (mStrainType == StrainType::STRAIN_Z) {
                    rResponseGradient[local_index++] = strain_sensitivities(2, 2);
                }
            // }
            // if (r_node.HasDofFor(ADJOINT_ROTATION_Y)) {
                r_node.FastGetSolutionStepValue(ROTATION_Y) += delta;
                r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, perturbed_strains, rProcessInfo);
                r_node.FastGetSolutionStepValue(ROTATION_Y) -= delta;
                strain_sensitivities.clear();
                for (IndexType j = 0; j < perturbed_strains.size(); ++j) {
                    noalias(strain_sensitivities) += (perturbed_strains[j] - ref_strains[j]) / delta;
                }
                if (mStrainType == StrainType::STRAIN_X) {
                    rResponseGradient[local_index++] = strain_sensitivities(0, 0);
                } else if (mStrainType == StrainType::STRAIN_Y) {
                    rResponseGradient[local_index++] = strain_sensitivities(1, 1);
                } else if (mStrainType == StrainType::STRAIN_Z) {
                    rResponseGradient[local_index++] = strain_sensitivities(2, 2);
                }
            // }
            // if (r_node.HasDofFor(ADJOINT_ROTATION_Z)) {
                r_node.FastGetSolutionStepValue(ROTATION_Z) += delta;
                r_element.CalculateOnIntegrationPoints(SHELL_STRAIN, perturbed_strains, rProcessInfo);
                r_node.FastGetSolutionStepValue(ROTATION_Z) -= delta;
                strain_sensitivities.clear();
                for (IndexType j = 0; j < perturbed_strains.size(); ++j) {
                    noalias(strain_sensitivities) += (perturbed_strains[j] - ref_strains[j]) / delta;
                }
                if (mStrainType == StrainType::STRAIN_X) {
                    rResponseGradient[local_index++] = strain_sensitivities(0, 0);
                } else if (mStrainType == StrainType::STRAIN_Y) {
                    rResponseGradient[local_index++] = strain_sensitivities(1, 1);
                } else if (mStrainType == StrainType::STRAIN_Z) {
                    rResponseGradient[local_index++] = strain_sensitivities(2, 2);
                }
            // }
        }

        rResponseGradient /= (mWeight * ref_strains.size());
    }

    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    KRATOS_CATCH("");
}

void AdjointElementStrainResponseFunction::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    KRATOS_CATCH("");
}

} // namespace Kratos.


