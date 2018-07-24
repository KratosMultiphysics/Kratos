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
    
    double& integrationWeight = GetValue(INTEGRATION_WEIGHT);
    Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
    Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    auto& properties = GetProperties();

    double E = properties[YOUNG_MODULUS];
    double A = properties[CROSS_AREA];
    double prestress = properties[PRESTRESS_CAUCHY];

    // compute base vectors

    Vector actualBaseVector = ZeroVector(3);
    GetBaseVector(actualBaseVector, DN_De);

    double referenceA = norm_2(mBaseVector0);
    double actualA = norm_2(actualBaseVector);

    // green-lagrange strain

    double E11_membrane = 0.5 * (actualA * actualA - referenceA * referenceA);

    // normal force

    double S11_membrane = prestress * A + E11_membrane * A * E / (referenceA *
        referenceA);

    // 1st variation of the axial strain

    Vector epsilonVar1Dof = ZeroVector(number_of_dofs);
    Get1stVariationsAxialStrain(epsilonVar1Dof, actualBaseVector, 3, DN_De);
    epsilonVar1Dof = epsilonVar1Dof / (referenceA * referenceA);

    // 2nd variation of the axial strain 

    Matrix epsilonVar2Dof = ZeroMatrix(number_of_dofs, number_of_dofs);
    Get2ndVariationsAxialStrain(epsilonVar2Dof, 3, DN_De);
    epsilonVar2Dof = epsilonVar2Dof / (referenceA * referenceA);

    for (size_t r = 0; r < number_of_dofs; r++) {
        for (size_t s = 0; s < number_of_dofs; s++) {
            rLeftHandSideMatrix(r, s) = E * A * epsilonVar1Dof[r] *
                epsilonVar1Dof[s] + S11_membrane * epsilonVar2Dof(r, s);
        }
    }

    rRightHandSideVector = -S11_membrane * epsilonVar1Dof;

    rLeftHandSideMatrix *= integrationWeight * referenceA;
    rRightHandSideVector *= integrationWeight * referenceA;

    KRATOS_CATCH("");
}

} // namespace Kratos
