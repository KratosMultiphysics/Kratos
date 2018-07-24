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

constexpr size_t TrussDiscreteElement::DofsPerNode()
{
    return 3;
};

size_t TrussDiscreteElement::NumberOfNodes() const
{
    return GetGeometry().size();
};

size_t TrussDiscreteElement::NumberOfDofs() const
{
    return NumberOfNodes() * DofsPerNode();
};

void TrussDiscreteElement::Initialize()
{
    KRATOS_TRY

    Vector base_vector = ZeroVector(3);
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

    const size_t number_of_dofs = NumberOfDofs();

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
    
    double& integration_weight = GetValue(INTEGRATION_WEIGHT);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    auto& properties = GetProperties();

    double E = properties[YOUNG_MODULUS];
    double A = properties[CROSS_AREA];
    double prestress = properties[PRESTRESS_CAUCHY];

    // compute base vectors

    Vector actual_base_vector = ZeroVector(3);
    GetBaseVector(actual_base_vector, shape_derivatives);

    double reference_a = norm_2(mBaseVector0);
    double actual_a = norm_2(actual_base_vector);

    // green-lagrange strain

    double e11_membrane = 0.5 * (actual_a * actual_a - reference_a *
        reference_a);

    // normal force

    double s11_membrane = prestress * A + e11_membrane * A * E /
        (reference_a * reference_a);

    // 1st variation of the axial strain

    Vector epsilon_var_1_dof = ZeroVector(number_of_dofs);
    Get1stVariationsAxialStrain(epsilon_var_1_dof, actual_base_vector, 3,
        shape_derivatives);
    epsilon_var_1_dof = epsilon_var_1_dof / (reference_a * reference_a);

    // 2nd variation of the axial strain 

    Matrix epsilon_var_2 = ZeroMatrix(number_of_dofs, number_of_dofs);
    Get2ndVariationsAxialStrain(epsilon_var_2, 3, shape_derivatives);
    epsilon_var_2 = epsilon_var_2 / (reference_a * reference_a);

    for (size_t r = 0; r < number_of_dofs; r++) {
        for (size_t s = 0; s < number_of_dofs; s++) {
            rLeftHandSideMatrix(r, s) = E * A * epsilon_var_1_dof[r] *
                epsilon_var_1_dof[s] + s11_membrane * epsilon_var_2(r, s);
        }
    }

    rRightHandSideVector = -s11_membrane * epsilon_var_1_dof;

    rLeftHandSideMatrix *= integration_weight * reference_a;
    rRightHandSideVector *= integration_weight * reference_a;

    KRATOS_CATCH("");
}

} // namespace Kratos
