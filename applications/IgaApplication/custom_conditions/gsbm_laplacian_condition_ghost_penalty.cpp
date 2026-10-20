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
#include <algorithm>
#include <cmath>

// External includes

// Project includes
#include "custom_conditions/gsbm_laplacian_condition_ghost_penalty.h"
#include "iga_application_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace
{

using IndexType = GsbmLaplacianConditionGhostPenalty::IndexType;
using SizeType = GsbmLaplacianConditionGhostPenalty::SizeType;
using GeometryType = Geometry<Node>;

/**
 * @brief Derivatives of order Order of a facet side's shape functions along
 * the unit vector rDirection.
 *
 * @param rSide One side of the facet.
 * @param rDirection Unit vector of the derivatives.
 * @param Order Order of the derivatives.
 * @return The directional derivatives of the shape functions.
 */
Vector ComputeDirectionalDerivatives(
    const GeometryType& rSide,
    const array_1d<double, 3>& rDirection,
    const SizeType Order)
{
    // Order-th derivatives: column k is Order - k times in x and k times in y
    const Matrix& r_derivatives = rSide.ShapeFunctionDerivatives(
        Order, 0, rSide.GetDefaultIntegrationMethod()
    );
    Vector directional_derivatives = ZeroVector(rSide.size());

    for (IndexType k = 0; k <= Order; ++k) {
        // D_s^n N = (s . \nabla)^n N = \sum_k binom(n, k) s_x^(n-k) s_y^k
        //           * (column k), with n = Order and s = rDirection
        const double binomial = MathUtils<double>::Factorial(Order)
            / (MathUtils<double>::Factorial(Order - k) * MathUtils<double>::Factorial(k));
        const double factor =
            binomial * std::pow(rDirection[0], Order - k) * std::pow(rDirection[1], k);
        noalias(directional_derivatives) += factor * column(r_derivatives, k);
    }

    return directional_derivatives;
}

} // namespace

void GsbmLaplacianConditionGhostPenalty::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_unknown_variable =
        rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();

    // Both sides of the facet point
    const GeometryType& r_minus = GetGeometry().GetGeometryPart(0);  // (-) side
    const GeometryType& r_plus = GetGeometry().GetGeometryPart(1);   // (+) side
    const SizeType number_of_minus_dofs = r_minus.size();
    const SizeType number_of_dofs = number_of_minus_dofs + r_plus.size();

    // Local system over the dofs of both sides
    if (rLeftHandSideMatrix.size1() != number_of_dofs ||
        rLeftHandSideMatrix.size2() != number_of_dofs) {
        rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs);
    if (rRightHandSideVector.size() != number_of_dofs) {
        rRightHandSideVector.resize(number_of_dofs, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_dofs);

    // Volume fractions \alpha_h on each side
    const double alpha_minus = r_minus.GetValue(VOLUME_FRACTION);
    const double alpha_plus = r_plus.GetValue(VOLUME_FRACTION);

    // \beta_2 = \kappa if {{\alpha_h}} < 1 and |[[\alpha_h]]| < 1, else 0
    const double alpha_average = 0.5 * (alpha_minus + alpha_plus);
    if (alpha_average >= 1.0 || std::abs(alpha_minus - alpha_plus) >= 1.0) {
        return;  // skip early if \beta_2 = 0 (zero contribution)
    }

    // Facet measure, and the normal n- pointing from the - to the + side
    const IndexType point_index = 0;
    const double weight = r_minus.IntegrationPoints()[0].Weight()
                        * r_minus.DeterminantOfJacobian(point_index);
    const array_1d<double, 3> normal = r_minus.UnitNormal(point_index);

    // Characteristic size h: the smallest knot span size of both sides
    const Vector& r_sizes_minus = r_minus.GetValue(KNOT_SPAN_SIZES);
    const Vector& r_sizes_plus = r_plus.GetValue(KNOT_SPAN_SIZES);
    const double h = std::min(
        *std::min_element(r_sizes_minus.begin(), r_sizes_minus.end()),
        *std::min_element(r_sizes_plus.begin(), r_sizes_plus.end())
    );

    // Degree p_n across the facet: only the p_n-th derivative jumps there
    const IndexType normal_direction = std::abs(normal[0]) > std::abs(normal[1]) ? 0 : 1;
    const GeometryType& r_surface = r_minus.GetGeometryParent(0);
    const SizeType p_normal = r_surface.PolynomialDegree(normal_direction);

    // [[D_n^{p_n} v_h]] . n- = D_n^{p_n} v_h- - D_n^{p_n} v_h+
    Vector jump_normal_derivative(number_of_dofs);
    noalias(subrange(jump_normal_derivative, 0, number_of_minus_dofs)) =
        ComputeDirectionalDerivatives(r_minus, normal, p_normal);
    noalias(subrange(jump_normal_derivative, number_of_minus_dofs, number_of_dofs)) =
        -ComputeDirectionalDerivatives(r_plus, normal, p_normal);

    // Ghost penalty: (\kappa h^(2p_n-1) [[D_n^{p_n} v_h]], [[D_n^{p_n} u_h]])
    const double kappa = GetProperties()[GHOST_PENALTY_FACTOR];
    noalias(rLeftHandSideMatrix) += weight * kappa * std::pow(h, 2 * p_normal - 1)
        * outer_prod(jump_normal_derivative, jump_normal_derivative);

    // Residual form
    Vector current_values(number_of_dofs);
    for (IndexType i = 0; i < number_of_minus_dofs; ++i) {
        current_values[i] = r_minus[i].FastGetSolutionStepValue(r_unknown_variable);
    }
    for (IndexType i = number_of_minus_dofs; i < number_of_dofs; ++i) {
        current_values[i] =
            r_plus[i - number_of_minus_dofs].FastGetSolutionStepValue(r_unknown_variable);
    }
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, current_values);

    KRATOS_CATCH("")
}

void GsbmLaplacianConditionGhostPenalty::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmLaplacianConditionGhostPenalty::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmLaplacianConditionGhostPenalty::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_unknown_variable =
        rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();

    rResult.clear();
    for (IndexType side = 0; side < 2; ++side) {
        for (const auto& r_node : GetGeometry().GetGeometryPart(side)) {
            rResult.push_back(r_node.GetDof(r_unknown_variable).EquationId());
        }
    }
}

void GsbmLaplacianConditionGhostPenalty::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_unknown_variable =
        rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]->GetUnknownVariable();

    rConditionDofList.clear();
    for (IndexType side = 0; side < 2; ++side) {
        for (const auto& r_node : GetGeometry().GetGeometryPart(side)) {
            rConditionDofList.push_back(r_node.pGetDof(r_unknown_variable));
        }
    }
}

int GsbmLaplacianConditionGhostPenalty::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(GHOST_PENALTY_FACTOR))
        << "GsbmLaplacianConditionGhostPenalty #" << Id() << ": no GHOST_PENALTY_FACTOR "
        << "(kappa) defined in the properties." << std::endl;
    return 0;
}

} // namespace Kratos
