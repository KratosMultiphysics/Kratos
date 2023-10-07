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
#include "adjoint_nodal_disp_response_function.h"

namespace Kratos
{
AdjointNodalDispResponseFunction::AdjointNodalDispResponseFunction(
    ModelPart& rModelPart,
    Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
{
    mNodeId = ResponseSettings["node_id"].GetInt();

    mWeight = ResponseSettings["weight"].GetDouble();

    mDirection = ResponseSettings["direction"].GetVector();

    for (const auto& r_element : mrModelPart.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            if (r_geometry[i].Id() == mNodeId) {
                mElementId = r_element.Id();
                mIndexPos = i;
                break;
            }
        }

        if (mIndexPos >= 0) {
            break;
        }
    }
}

double AdjointNodalDispResponseFunction::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY;

    const auto& r_node = mrModelPart.GetNode(mNodeId);

    return inner_prod(r_node.FastGetSolutionStepValue(DISPLACEMENT), mDirection) / mWeight;

    KRATOS_CATCH("");
}

void AdjointNodalDispResponseFunction::CalculateGradient(
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
        const auto& r_geometry = rAdjointElement.GetGeometry();
        const IndexType block_size = rResidualGradient.size1() / r_geometry.size();
        const IndexType block_start_index = mIndexPos * block_size;

        rResponseGradient[block_start_index] = mDirection[0];
        rResponseGradient[block_start_index + 1] = mDirection[1];
        rResponseGradient[block_start_index + 2] = mDirection[2];
    }

    rResponseGradient /= mWeight;

    KRATOS_CATCH("");
}

void AdjointNodalDispResponseFunction::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointNodalDispResponseFunction::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointNodalDispResponseFunction::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointNodalDispResponseFunction::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;
    rResponseGradient = ZeroVector(rResidualGradient.size1());
    KRATOS_CATCH("");
}

void AdjointNodalDispResponseFunction::CalculatePartialSensitivity(
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

void AdjointNodalDispResponseFunction::CalculatePartialSensitivity(
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

void AdjointNodalDispResponseFunction::CalculatePartialSensitivity(
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

void AdjointNodalDispResponseFunction::CalculatePartialSensitivity(
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


