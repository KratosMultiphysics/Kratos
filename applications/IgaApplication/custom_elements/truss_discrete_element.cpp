/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//
//  Authors: Thomas Oberbichler
*/

// System includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"

// External includes

// Project includes
#include "custom_elements/base_discrete_element.h"
#include "custom_elements/truss_discrete_element.h"

// Application includes
#include "iga_application.h"
#include "iga_application_variables.h"

namespace Kratos
{

constexpr std::size_t TrussDiscreteElement::DofsPerNode()
{
    return 3;
};

std::size_t TrussDiscreteElement::NumberOfNodes() const
{
    return GetGeometry().size();
};

std::size_t TrussDiscreteElement::NumberOfDofs() const
{
    return NumberOfNodes() * DofsPerNode();
};

void TrussDiscreteElement::Initialize()
{
    KRATOS_TRY

    Vector3D base_vector = ZeroVector(3);
    GetBaseVector(base_vector, GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));
    mBaseVector0 = base_vector;

    KRATOS_CATCH("")
}

void TrussDiscreteElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_TRY

    const std::size_t number_of_dofs = NumberOfDofs();

    if (CalculateStiffnessMatrixFlag) {
        if (rLeftHandSideMatrix.size1() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs,
            number_of_dofs);
    }

    if (CalculateResidualVectorFlag == true) {
        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }

        noalias(rRightHandSideVector) = ZeroVector(number_of_dofs);
    }

    // get integration data
    
    const double& integration_weight = GetValue(INTEGRATION_WEIGHT);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    const auto& properties = GetProperties();

    const double E = properties[YOUNG_MODULUS];
    const double A = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];

    // compute base vectors

    Vector3D actual_base_vector = ZeroVector(3);
    GetBaseVector(actual_base_vector, shape_derivatives);

    const double reference_a = norm_2(mBaseVector0);
    const double actual_a = norm_2(actual_base_vector);

    const double actual_aa = actual_a * actual_a;
    const double reference_aa = reference_a * reference_a;

    // green-lagrange strain

    const double e11_membrane = 0.5 * (actual_aa - reference_aa);

    // normal force

    const double s11_membrane = prestress * A + e11_membrane * A * E /
        reference_aa;

    for (size_t r = 0; r < number_of_dofs; r++) {
        const size_t dof_type_r = r % DofsPerNode();
        const size_t shape_index_r = r / DofsPerNode();

        const double epsilon_var_r = actual_base_vector[dof_type_r] *
            shape_derivatives(shape_index_r, 0) / reference_aa;

        for (size_t s = 0; s < number_of_dofs; s++) {
            const size_t dof_type_s = s % DofsPerNode();
            const size_t shape_index_s = s / DofsPerNode();

            const double epsilon_var_s = actual_base_vector[dof_type_s] *
                shape_derivatives(shape_index_s, 0) / reference_aa;

            rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r * epsilon_var_s;

            if (dof_type_r == dof_type_s) {
                const double epsilon_var_rs =
                    shape_derivatives(shape_index_r, 0) *
                    shape_derivatives(shape_index_s, 0) / reference_aa;

                rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs;
            }
        }

        rRightHandSideVector[r] = -s11_membrane * epsilon_var_r;
    }

    rLeftHandSideMatrix *= integration_weight * reference_a;
    rRightHandSideVector *= integration_weight * reference_a;

    KRATOS_CATCH("");
}

} // namespace Kratos
