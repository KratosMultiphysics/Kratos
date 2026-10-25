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

constexpr IndexType LaplacianOrder = 2;

/// \tau_d = - t^q / q!
double ComputeTau(const array_1d<double, 3>& rDistance)
{
    return - std::pow(norm_2(rDistance), LaplacianOrder)
           / MathUtils<double>::Factorial(LaplacianOrder);
}

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

/// Whether u'_h is added: USE_RESIDUAL_TERMS, off by default
bool UsesSubscale(const Properties& rProperties)
{
    return rProperties.Has(USE_RESIDUAL_TERMS) ? rProperties[USE_RESIDUAL_TERMS] : false;
}

/**
 * @brief Shape functions of a facet side shifted to the true boundary
 *
 *     S_d^m[N] = \sum_{|\boldsymbol\alpha| \le m} c_\boldsymbol\alpha
 *                \partial^{\boldsymbol\alpha} N(\tilde{x})
 *
 * @param rSide One side of the facet.
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param Order Order of the shift.
 * @return The shifted shape functions.
 */
Vector ComputeShiftedShapeFunctions(
    const GeometryType& rSide,
    const array_1d<double, 3>& rDistance,
    const SizeType Order)
{
    Vector shifted_shape_functions = row(rSide.ShapeFunctionsValues(), 0);

    for (IndexType k = 1; k <= Order; ++k) {
        // Derivatives of order k: column i is \partial_x^(k-i) \partial_y^i N
        const Matrix& r_derivatives = rSide.ShapeFunctionDerivatives(
            k, 0, rSide.GetDefaultIntegrationMethod()
        );

        // S_d^m, \boldsymbol\alpha = (k - i, i)
        for (IndexType i = 0; i <= k; ++i) {
            const double factor = ComputeTaylorCoefficient(rDistance, {k - i, i});
            noalias(shifted_shape_functions) += factor * column(r_derivatives, i);
        }
    }

    return shifted_shape_functions;
}

/**
 * @brief Shifted subscale of a facet side's shape functions
 *
 *     u'_h[N] = \tau_d \sum_{k=q}^{m} \binom{k}{q}^{-1}
 *               \sum_{|\boldsymbol\beta| = k-q} c_\boldsymbol\beta
 *               \partial^{\boldsymbol\beta} \Delta N(\tilde{x})
 *
 * @param rSide One side of the facet.
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param Order Order of the shift.
 * @return The shifted subscale of the shape functions.
 */
Vector ComputeShiftedSubscaleShapeFunctions(
    const GeometryType& rSide,
    const array_1d<double, 3>& rDistance,
    const SizeType Order)
{
    const IndexType q = LaplacianOrder;
    const double tau = ComputeTau(rDistance);
    Vector shifted_subscale = ZeroVector(rSide.size());

    // u'_h[N] at \tilde{x}: \tau_d \sum_k (q!/k!) \mathsf D_d^{k-q} \Delta N
    for (IndexType k = q; k <= Order; ++k) {
        // Derivatives of order k: column i is \partial_x^(k-i) \partial_y^i N
        const Matrix& r_derivatives = rSide.ShapeFunctionDerivatives(
            k, 0, rSide.GetDefaultIntegrationMethod()
        );

        // \tau_d \binom{k}{q}^{-1}
        const double weight = tau * MathUtils<double>::Factorial(q)
            * MathUtils<double>::Factorial(k - q) / MathUtils<double>::Factorial(k);

        for (IndexType i = 0; i <= k - q; ++i) {
            // \boldsymbol\beta = (k-q-i, i); \Delta: columns i and i + 2
            const double factor = weight * ComputeTaylorCoefficient(rDistance, {k - q - i, i});
            noalias(shifted_subscale) +=
                factor * (column(r_derivatives, i) + column(r_derivatives, i + 2));
        }
    }

    return shifted_subscale;
}

/**
 * @brief Shifted subscale of the source
 *
 *     u'_h[f] = \tau_d \sum_{k=q}^{m} \binom{k}{q}^{-1}
 *               \sum_{|\boldsymbol\beta| = k-q} c_\boldsymbol\beta
 *               \partial^{\boldsymbol\beta} f(\tilde{x})
 *
 * @param rSourceDerivatives f and its Cartesian derivatives at \tilde{x}: order
 * n at offset n (n + 1) / 2, entry i being \partial_x^(n-i) \partial_y^i f.
 * @param rDistance Distance vector d from the facet point to the skin.
 * @param Order Order of the shift.
 * @return The shifted subscale of the source.
 */
double ComputeShiftedSubscaleSource(
    const Vector& rSourceDerivatives,
    const array_1d<double, 3>& rDistance,
    const SizeType Order)
{
    const IndexType q = LaplacianOrder;
    if (Order < q) {
        return 0.0;  // no Taylor term of order q or higher
    }

    // \partial^{\boldsymbol\beta} f(\tilde{x}), |\boldsymbol\beta| <= m - q
    const SizeType f_size = (Order - q + 1) * (Order - q + 2) / 2;
    KRATOS_ERROR_IF(rSourceDerivatives.size() < f_size) << "GsbmLaplacianConditionDirichlet: "
        << "VOLUME_SOURCE_DERIVATIVES at a cut facet point needs f and its Cartesian "
        << "derivatives of orders 0 to " << Order - q << " (size " << f_size
        << "), but it has size " << rSourceDerivatives.size() << "." << std::endl;

    // u'_h[f] at \tilde{x}: \tau_d \sum_k (q!/k!) \mathsf D_d^{k-q} f
    const double tau = ComputeTau(rDistance);
    double shifted_subscale = 0.0;
    for (IndexType k = q; k <= Order; ++k) {
        // \tau_d \binom{k}{q}^{-1}
        const double weight = tau * MathUtils<double>::Factorial(q)
            * MathUtils<double>::Factorial(k - q) / MathUtils<double>::Factorial(k);
        const IndexType offset = (k - q) * (k - q + 1) / 2;  // |\boldsymbol\beta| = k - q
        for (IndexType i = 0; i <= k - q; ++i) {
            // \boldsymbol\beta = (k - q - i, i)
            const double factor = weight * ComputeTaylorCoefficient(rDistance, {k - q - i, i});
            shifted_subscale += factor * rSourceDerivatives[offset + i];
        }
    }

    return shifted_subscale;
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

    // Shifted trace S(N) on each side, \mathsf S_d^m[N], and imposed value g_D
    Vector shifted_minus = ComputeShiftedShapeFunctions(r_minus, distance, p);
    Vector shifted_plus = ComputeShiftedShapeFunctions(r_plus, distance, p);
    double boundary_value = r_projection_node.GetValue(r_unknown_variable);

    // Shifted subscale: plus u'_h[N] on the trace, minus u'_h[f] on the value
    if (UsesSubscale(GetProperties())) {
        noalias(shifted_minus) += ComputeShiftedSubscaleShapeFunctions(r_minus, distance, p);
        noalias(shifted_plus) += ComputeShiftedSubscaleShapeFunctions(r_plus, distance, p);
        boundary_value -=
            ComputeShiftedSubscaleSource(GetValue(VOLUME_SOURCE_DERIVATIVES), distance, p);
    }

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
        (1.0 - gamma) * shifted_minus, gamma * shifted_plus
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
