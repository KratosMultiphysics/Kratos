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
#include <initializer_list>
#include <limits>

// External includes

// Project includes
#include "custom_conditions/gsbm_laplacian_condition_dirichlet.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_surface_geometry.h"
#include "iga_application_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/global_pointer_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace
{

using IndexType = GsbmLaplacianConditionDirichlet::IndexType;
using SizeType = GsbmLaplacianConditionDirichlet::SizeType;
using GeometryType = GsbmLaplacianConditionDirichlet::GeometryType;
using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;

/**
 * @brief Taylor coefficient of a multi-index
 *
 *     c_\boldsymbol\alpha = \mathbf{d}^{\boldsymbol\alpha} / \boldsymbol\alpha!
 *
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param MultiIndex The multi-index \boldsymbol\alpha, one entry per direction.
 * @return The coefficient c_\boldsymbol\alpha.
 */
double ComputeTaylorCoefficient(
    const array_1d<double, 3>& rDistance,
    const std::initializer_list<IndexType> MultiIndex)
{
    double power = 1.0;
    IndexType factorial = 1;

    // \prod_j{ d_j^{\alpha_j} / \alpha_j! }
    for (IndexType j = 0; const IndexType alpha_j : MultiIndex) {
        power *= std::pow(rDistance[j], alpha_j);
        factorial *= MathUtils<double>::Factorial(alpha_j);
        ++j;
    }
    return power / factorial;
}

/**
 * @brief Shape functions of a facet side shifted to the true boundary, from
 *        the side's 1D B-spline bases
 *
 *     S_d^m[N] = \sum_{|\boldsymbol\alpha| \le m} c_\boldsymbol\alpha
 *                \partial^{\boldsymbol\alpha} N(\tilde{x}),
 *                \alpha_u \le p_u, \alpha_v \le p_v
 *
 * m = p_u + p_v is the full Taylor series, the span's polynomial at
 * \tilde{x} + d. The 1D derivatives never exceed the degree, so every order m
 * is available.
 *
 * @param rSide One side of the facet.
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param Order Order m of the shift.
 * @return The shifted shape functions.
 */
Vector ComputeShiftedShapeFunctions(
    const GeometryType& rSide,
    const array_1d<double, 3>& rDistance,
    const SizeType Order)
{
    const auto& r_surface = dynamic_cast<const NurbsSurfaceType&>(rSide.GetGeometryParent(0));
    const SizeType p_u = r_surface.PolynomialDegree(0);
    const SizeType p_v = r_surface.PolynomialDegree(1);

    // 1D bases of the side's span at its point, derivatives up to the degree
    const auto& r_indices = rSide.GetValue(KNOT_SPAN_INDICES);
    const auto& r_point = rSide.IntegrationPoints()[0];
    NurbsCurveShapeFunction shape_functions_u(p_u, p_u);
    shape_functions_u.ComputeBSplineShapeFunctionValuesAtSpan(
        r_surface.KnotsU(), r_indices[0] + p_u - 1, r_point.X()
    );
    NurbsCurveShapeFunction shape_functions_v(p_v, p_v);
    shape_functions_v.ComputeBSplineShapeFunctionValuesAtSpan(
        r_surface.KnotsV(), r_indices[1] + p_v - 1, r_point.Y()
    );

    // c_{(i, j)} N_a^{(i)} N_b^{(j)} for N_{a + b (p_u + 1)}, uv = xy
    Vector shifted_shape_functions = ZeroVector(rSide.size());
    for (IndexType i = 0; i <= std::min(p_u, Order); ++i) {
        for (IndexType j = 0; j <= std::min(p_v, Order - i); ++j) {
            const double factor = ComputeTaylorCoefficient(rDistance, {i, j});
            for (IndexType b = 0; b <= p_v; ++b) {
                for (IndexType a = 0; a <= p_u; ++a) {
                    shifted_shape_functions[a + b * (p_u + 1)] +=
                        factor * shape_functions_u(a, i) * shape_functions_v(b, j);
                }
            }
        }
    }

    return shifted_shape_functions;
}

} // namespace

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

    // Degree p of the background patch, and shift order m: SHIFT_ORDER, p by
    // default
    const GeometryType& r_surface = r_minus.GetGeometryParent(0);
    const SizeType p = std::max(r_surface.PolynomialDegree(0), r_surface.PolynomialDegree(1));
    const SizeType m = GetProperties().Has(SHIFT_ORDER)
        ? static_cast<SizeType>(GetProperties()[SHIFT_ORDER])
        : p;

    // Shifted trace \mathsf S_d^m[N] on each side, and imposed value g_D
    const Vector shifted_minus = ComputeShiftedShapeFunctions(r_minus, distance, m);
    const Vector shifted_plus = ComputeShiftedShapeFunctions(r_plus, distance, m);
    const double boundary_value = r_projection_node.GetValue(r_unknown_variable);

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

    // {{S(v_h)}}_\gamma
    const Vector average_shifted_value = join(
        (1.0 - gamma) * shifted_minus, gamma * shifted_plus
    );

    // Consistency: -([[v_\alpha]], {{\nabla u_h}}_\gamma)
    noalias(rLeftHandSideMatrix) -=
        weight * outer_prod(jump_v_alpha, average_normal_derivative);

    // Adjoint: -\theta ({{\nabla v_h}}_\gamma,
    //                   [[\alpha_h]] ({{S(u_h)}}_\gamma - g_D)),
    // \theta = ADJOINT_FACTOR: -1 skew-symmetric (default), 1 symmetric
    const double theta =
        GetProperties().Has(ADJOINT_FACTOR) ? GetProperties()[ADJOINT_FACTOR] : -1.0;
    noalias(rLeftHandSideMatrix) -= theta * weight * alpha_jump
        * outer_prod(average_normal_derivative, average_shifted_value);
    noalias(rRightHandSideVector) -=
        theta * weight * alpha_jump * boundary_value * average_normal_derivative;

    // Penalty, if \chi = PENALTY_FACTOR > 0 (none by default):
    // (\chi (p + 1)^2 / h |[[\alpha_h]]| {{S(v_h)}}_\gamma,
    //  {{S(u_h)}}_\gamma - g_D)
    const double chi = GetProperties().Has(PENALTY_FACTOR) ? GetProperties()[PENALTY_FACTOR] : 0.0;
    if (chi > 0.0) {
        // Characteristic size h: the smallest knot span size of both sides
        const Vector& r_sizes_minus = r_minus.GetValue(KNOT_SPAN_SIZES);
        const Vector& r_sizes_plus = r_plus.GetValue(KNOT_SPAN_SIZES);
        const double h = std::min(
            *std::min_element(r_sizes_minus.begin(), r_sizes_minus.end()),
            *std::min_element(r_sizes_plus.begin(), r_sizes_plus.end())
        );

        const double penalty = weight * chi * (p + 1) * (p + 1) / h * std::abs(alpha_jump);
        noalias(rLeftHandSideMatrix) +=
            penalty * outer_prod(average_shifted_value, average_shifted_value);
        noalias(rRightHandSideVector) += penalty * boundary_value * average_shifted_value;
    }

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
    KRATOS_ERROR_IF(GetProperties().Has(SHIFT_ORDER) && GetProperties()[SHIFT_ORDER] < 0)
        << Info() << ": SHIFT_ORDER must not be negative." << std::endl;
    return 0;
}

} // namespace Kratos
