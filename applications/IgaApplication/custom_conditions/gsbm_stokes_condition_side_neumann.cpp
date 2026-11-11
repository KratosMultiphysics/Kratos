//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes

// External includes

// Project includes
#include "custom_conditions/gsbm_stokes_condition_side_neumann.h"
#include "custom_elements/gsbm_stokes_element.h"
#include "includes/variables.h"

namespace Kratos
{

void GsbmStokesConditionSideNeumann::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Local system over the dofs of the side point: u_x, u_y, then p
    SizeType number_of_dofs = 0;
    for (IndexType part = 0; part < 3; ++part) {
        number_of_dofs += GetGeometry().GetGeometryPart(part).size();
    }
    if (rLeftHandSideMatrix.size1() != number_of_dofs ||
        rLeftHandSideMatrix.size2() != number_of_dofs) {
        rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs);
    if (rRightHandSideVector.size() != number_of_dofs) {
        rRightHandSideVector.resize(number_of_dofs, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_dofs);

    // Side measure, and the outward normal n
    const GeometryType& r_reference =
        GetGeometry().GetGeometryPart(GsbmStokesElement::VelocityXPart);
    const double weight = r_reference.IntegrationPoints()[0].Weight()
                        * r_reference.DeterminantOfJacobian(0);
    const array_1d<double, 3> normal = r_reference.UnitNormal(0);

    // t_N = -p_D n + NORMAL_STRESS, as in SupportPressureCondition
    Vector traction(2);
    for (IndexType m = 0; m < 2; ++m) {
        traction[m] = -GetValue(PRESSURE) * normal[m];
    }
    if (Has(NORMAL_STRESS)) {
        const Vector& r_normal_stress = GetValue(NORMAL_STRESS);
        KRATOS_ERROR_IF(r_normal_stress.size() != 2) << Info()
            << ": NORMAL_STRESS has size " << r_normal_stress.size() << ", expected 2."
            << std::endl;
        noalias(traction) += r_normal_stress;
    }

    // (v, t_N): component c in the rows of u_c
    IndexType begin = 0;
    for (IndexType c = 0; c < 2; ++c) {
        const Vector N = row(GetGeometry().GetGeometryPart(c).ShapeFunctionsValues(), 0);
        noalias(subrange(rRightHandSideVector, begin, begin + N.size())) +=
            weight * traction[c] * N;
        begin += N.size();
    }

    KRATOS_CATCH("")
}

void GsbmStokesConditionSideNeumann::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmStokesConditionSideNeumann::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmStokesConditionSideNeumann::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rResult.clear();
    for (IndexType part = 0; part < 3; ++part) {
        const auto& r_variable = GsbmStokesElement::GetFieldVariable(part);
        for (const auto& r_node : GetGeometry().GetGeometryPart(part)) {
            rResult.push_back(r_node.GetDof(r_variable).EquationId());
        }
    }
}

void GsbmStokesConditionSideNeumann::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rConditionDofList.clear();
    for (IndexType part = 0; part < 3; ++part) {
        const auto& r_variable = GsbmStokesElement::GetFieldVariable(part);
        for (const auto& r_node : GetGeometry().GetGeometryPart(part)) {
            rConditionDofList.push_back(r_node.pGetDof(r_variable));
        }
    }
}

int GsbmStokesConditionSideNeumann::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(GetGeometry().NumberOfGeometryParts() != 3) << Info()
        << " needs a side point with three parts, the spaces of u_x, u_y and p. Create it "
        << "with \"sides\" in IgaModelerGsbm." << std::endl;

    return 0;

    KRATOS_CATCH("")
}

} // namespace Kratos
