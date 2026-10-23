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
#include <limits>

// External includes

// Project includes
#include "custom_conditions/gsbm_laplacian_condition_dirichlet.h"
#include "iga_application_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/global_pointer_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

Vector GsbmLaplacianConditionDirichlet::CalculateShiftedShapeFunctions(
    const GeometryType& rSide,
    const Node& rProjectionNode,
    const array_1d<double, 3>& rDistance,
    const SizeType Order) const
{
    Vector shifted_shape_functions = row(rSide.ShapeFunctionsValues(), 0);

    for (IndexType n = 1; n <= Order; ++n) {
        // n-th derivatives: column k is n - k times in x and k times in y
        const Matrix& r_derivatives = rSide.ShapeFunctionDerivatives(
            n, 0, rSide.GetDefaultIntegrationMethod()
        );
        for (IndexType k = 0; k <= n; ++k) {
            // Taylor coefficient of column k, with d = (d_x, d_y) the distance
            // vector: d_x^(n-k) d_y^k / ((n-k)! k!)
            const double factor = std::pow(rDistance[0], n - k) * std::pow(rDistance[1], k)
                / (MathUtils<double>::Factorial(n - k) * MathUtils<double>::Factorial(k));
            noalias(shifted_shape_functions) += factor * column(r_derivatives, k);
        }
    }

    return shifted_shape_functions;
}

double GsbmLaplacianConditionDirichlet::CalculateShiftedBoundaryValue(
    const Node& rProjectionNode,
    const array_1d<double, 3>& rDistance,
    const SizeType Order,
    const Variable<double>& rUnknownVariable) const
{
    return rProjectionNode.GetValue(rUnknownVariable);
}

void GsbmLaplacianConditionDirichlet::CalculateLocalSystem(
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
    if (alpha_minus == alpha_plus) {
        return;  // skip early if jump [[\alpha_h]] = 0 (zero contribution)
    }

    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_NODES)) << Info()
        << " is cut but has no projection node (NEIGHBOUR_NODES). Run "
        << "VolumeFractionGsbmProcess first." << std::endl;

    // Distance d to the closest point
    const Node& r_projection_node = GetValue(NEIGHBOUR_NODES)[0];
    const array_1d<double, 3> distance =
        r_projection_node.Coordinates() - r_minus.Center().Coordinates();

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

    // Degree p of the background patch, and \beta_1 = \chi (p + 1)^2
    const GeometryType& r_surface = r_minus.GetGeometryParent(0);
    const SizeType p = std::max(r_surface.PolynomialDegree(0), r_surface.PolynomialDegree(1));
    const double beta_1 = GetProperties()[PENALTY_FACTOR] * (p + 1) * (p + 1);

    // Value imposed on the shifted trace: the boundary value g_D
    const double boundary_value = CalculateShiftedBoundaryValue(
        r_projection_node, distance, p, r_unknown_variable
    );

    // \gamma = \alpha_h+ / max(\alpha_h+ + \alpha_h-, tol)
    // Average operator: {{.}}_\gamma = \gamma (.)+ + (1 - \gamma) (.)-
    const double gamma =
        alpha_plus / std::max(alpha_plus + alpha_minus, std::numeric_limits<double>::epsilon());

    // [[\alpha_h]] = (\alpha_h- - \alpha_h+) n-
    const double alpha_jump = alpha_minus - alpha_plus;

    // Operators on the local dofs: the - side first, then the + side
    const auto join = [&](const Vector& rMinus, const Vector& rPlus) {
        Vector joined(number_of_dofs);
        noalias(subrange(joined, 0, number_of_minus_dofs)) = rMinus;
        noalias(subrange(joined, number_of_minus_dofs, number_of_dofs)) = rPlus;
        return joined;
    };

    // [[v_\alpha]] . n- = \alpha_h- v_h- - \alpha_h+ v_h+
    const Vector jump_v_alpha = join(
        alpha_minus * row(r_minus.ShapeFunctionsValues(), 0),
        -alpha_plus * row(r_plus.ShapeFunctionsValues(), 0)
    );

    // One-sided \nabla v_h . n-: normal derivatives of a side's shape functions
    const auto compute_normal_derivatives = [&normal](const GeometryType& rSide) {
        const Matrix& r_DN_De = rSide.ShapeFunctionsLocalGradients()[0];
        return Vector(column(r_DN_De, 0) * normal[0] + column(r_DN_De, 1) * normal[1]);
    };

    // {{\nabla v_h}}_\gamma . n-
    const Vector average_normal_derivative = join(
        (1.0 - gamma) * compute_normal_derivatives(r_minus),
        gamma * compute_normal_derivatives(r_plus)
    );

    // {{S(v_h)}}_\gamma, with the shift of order p
    const Vector average_shifted_value = join(
        (1.0 - gamma) * CalculateShiftedShapeFunctions(r_minus, r_projection_node, distance, p),
        gamma * CalculateShiftedShapeFunctions(r_plus, r_projection_node, distance, p)
    );

    // Consistency: -([[v_\alpha]], {{\nabla u_h}}_\gamma)
    noalias(rLeftHandSideMatrix) -=
        weight * outer_prod(jump_v_alpha, average_normal_derivative);

    // Symmetry: -({{\nabla v_h}}_\gamma,
    //             [[\alpha_h]] ({{S(u_h)}}_\gamma - g_D))
    noalias(rLeftHandSideMatrix) -=
        weight * alpha_jump * outer_prod(average_normal_derivative, average_shifted_value);
    noalias(rRightHandSideVector) -=
        weight * alpha_jump * boundary_value * average_normal_derivative;

    // Penalty (stabilization): (\beta_1 / h |[[\alpha_h]]| {{S(v_h)}}_\gamma,
    //                           {{S(u_h)}}_\gamma - g_D)
    const double penalty = weight * beta_1 / h * std::abs(alpha_jump);
    noalias(rLeftHandSideMatrix) +=
        penalty * outer_prod(average_shifted_value, average_shifted_value);
    noalias(rRightHandSideVector) += penalty * boundary_value * average_shifted_value;

    // Residual form, as LaplacianElement: subtract the current solution's part
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

void GsbmLaplacianConditionDirichlet::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmLaplacianConditionDirichlet::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmLaplacianConditionDirichlet::EquationIdVector(
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

void GsbmLaplacianConditionDirichlet::GetDofList(
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

int GsbmLaplacianConditionDirichlet::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << Info() << ": no PENALTY_FACTOR (chi) defined in the properties." << std::endl;
    return 0;
}

} // namespace Kratos
