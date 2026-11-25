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
#include <array>
#include <cmath>
#include <limits>
#include <map>

// External includes

// Project includes
#include "custom_conditions/gsbm_stokes_condition_dirichlet_multiplier.h"
#include "custom_elements/gsbm_stokes_element.h"
#include "custom_utilities/shifted_boundary_utilities.h"
#include "geometries/nurbs_curve_geometry.h"
#include "iga_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

namespace
{

using MultiplierCurveType = NurbsCurveGeometry<2, PointerVector<Node>>;

/**
 * @brief Viscous traction D B N \nu of the velocity dofs of a facet side, from
 *        the constitutive matrix and the side's velocity gradients
 *
 * @param rD Constitutive matrix D of the Voigt strain [xx, yy, xy].
 * @param rGradientsX Gradients of the shape functions of u_x's space.
 * @param rGradientsY Gradients of the shape functions of u_y's space.
 * @param rNormal Unit normal \nu of the traction.
 * @return Rows x and y of the traction, over the dofs of u_x, then of u_y.
 */
Matrix ComputeMultiplierViscousTraction(
    const Matrix& rD,
    const Matrix& rGradientsX,
    const Matrix& rGradientsY,
    const array_1d<double, 3>& rNormal)
{
    const std::size_t n_x = rGradientsX.size1();
    const std::size_t n_u = n_x + rGradientsY.size1();

    // Voigt strain [xx, yy, xy] of the velocity dofs: u_x's, then u_y's
    Matrix B = ZeroMatrix(3, n_u);
    for (std::size_t i = 0; i < n_x; ++i) {
        B(0, i) = rGradientsX(i, 0);
        B(2, i) = rGradientsX(i, 1);
    }
    for (std::size_t i = 0; i < rGradientsY.size1(); ++i) {
        B(1, n_x + i) = rGradientsY(i, 1);
        B(2, n_x + i) = rGradientsY(i, 0);
    }
    const Matrix stress = prod(rD, B);

    // (D B) \nu, with the Voigt stress [\sigma_xx, \sigma_yy, \sigma_xy]
    Matrix traction(2, n_u);
    for (std::size_t k = 0; k < n_u; ++k) {
        traction(0, k) = stress(0, k) * rNormal[0] + stress(2, k) * rNormal[1];
        traction(1, k) = stress(2, k) * rNormal[0] + stress(1, k) * rNormal[1];
    }

    return traction;
}

/// The multiplier curve of a facet geometry, or an error
const MultiplierCurveType& GetMultiplierCurve(const Geometry<Node>& rGeometry)
{
    const auto p_curve = dynamic_cast<const MultiplierCurveType*>(
        &rGeometry.GetGeometryPart(GsbmStokesConditionDirichletMultiplier::MultiplierPart)
    );
    KRATOS_ERROR_IF(p_curve == nullptr)
        << "GsbmStokesConditionDirichletMultiplier: the third part of the facet is not "
        << "a multiplier curve. Create the facets with \"multiplier_parameters\"." << std::endl;
    return *p_curve;
}

/**
 * @brief Shares c_k = m_k / \sum_j m_j of the periodic multipliers, with
 *        m_k = (t_{k+p+1} - t_k) / (p + 1) the integral of \mu_k over the
 *        curve's parameter
 *
 * @param rCurve Multiplier curve of degree p, with the first p points again.
 * @param NumberOfMultipliers Number n of multipliers.
 * @return c_k for k = 0, ..., n - 1.
 */
Vector ComputeMultiplierShares(
    const MultiplierCurveType& rCurve,
    const std::size_t NumberOfMultipliers)
{
    const Vector& r_knots = rCurve.Knots();
    const std::size_t p = rCurve.PolynomialDegree(0);
    const std::size_t n = NumberOfMultipliers;

    // t_j = knots[j - 1], without the first and last knots, as Kratos stores
    // them, and t_0 = t_n - (t_{n+p} - t_p) by periodicity
    const auto get_knot = [&](const std::size_t J) {
        return J > 0 ? r_knots[J - 1] : r_knots[n - 1] - (r_knots[n + p - 1] - r_knots[p - 1]);
    };

    Vector shares(n);
    for (std::size_t k = 0; k < n; ++k) {
        shares[k] = (get_knot(k + p + 1) - get_knot(k)) / static_cast<double>(p + 1);
    }
    return shares / sum(shares);
}

} // namespace

void GsbmStokesConditionDirichletMultiplier::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}

void GsbmStokesConditionDirichletMultiplier::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw) << Info()
        << " has no constitutive law. Call Initialize first." << std::endl;

    // Both sides of the facet point, each with a part per field: u_x, u_y, p
    const GeometryType& r_minus = GetGeometry().GetGeometryPart(0);  // (-) side
    const GeometryType& r_plus = GetGeometry().GetGeometryPart(1);   // (+) side

    // First local dof of each field: u_x, u_y and p of the - side, then + side
    std::array<std::array<IndexType, 3>, 2> first_dof;
    SizeType number_of_dofs = 0;
    for (IndexType side = 0; side < 2; ++side) {
        for (IndexType part = 0; part < 3; ++part) {
            first_dof[side][part] = number_of_dofs;
            number_of_dofs += GetGeometry().GetGeometryPart(side).GetGeometryPart(part).size();
        }
    }

    // The multipliers after them
    const MultiplierCurveType& r_curve = GetMultiplierCurve(GetGeometry());
    const SizeType number_of_multipliers = NumberOfMultipliers();
    const IndexType first_multiplier = number_of_dofs;
    number_of_dofs += number_of_multipliers;

    // Local system over the dofs of both sides and the multipliers
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

    // Facet measure, and the normal n- pointing from the - to the + side
    const GeometryType& r_reference = r_minus.GetGeometryPart(GsbmStokesElement::VelocityXPart);
    const double weight = r_reference.IntegrationPoints()[0].Weight()
                        * r_reference.DeterminantOfJacobian(0);
    const array_1d<double, 3> normal = r_reference.UnitNormal(0);

    // Distance d, imposed velocity g and unit skin tangent \tau there
    const Node& r_projection_node = GetValue(NEIGHBOUR_NODES)[0];
    const array_1d<double, 3> distance =
        r_projection_node.Coordinates() - r_reference.Center().Coordinates();
    const array_1d<double, 3>& r_imposed_velocity = r_projection_node.GetValue(VELOCITY);
    Vector velocity_data(2);
    velocity_data[0] = r_imposed_velocity[0];
    velocity_data[1] = r_imposed_velocity[1];
    const array_1d<double, 3>& r_skin_tangent = r_projection_node.GetValue(LOCAL_TANGENT);
    const array_1d<double, 3> tangent = r_skin_tangent / norm_2(r_skin_tangent);

    // \mu_k at the closest point's normalized skin parameter t
    KRATOS_ERROR_IF_NOT(r_projection_node.Has(INTEGRATION_COORDINATES)) << Info()
        << ": the projection node has no skin parameter (INTEGRATION_COORDINATES). The "
        << "multipliers need a NURBS skin." << std::endl;
    array_1d<double, 3> skin_parameter = ZeroVector(3);
    skin_parameter[0] = r_projection_node.GetValue(INTEGRATION_COORDINATES)[0];
    std::vector<IndexType> multiplier_ids;
    Vector multiplier_values;
    r_curve.ShapeFunctionsValuesAndCPIndices(skin_parameter, multiplier_ids, multiplier_values);
    std::map<IndexType, IndexType> multiplier_index;
    for (IndexType k = 0; k < number_of_multipliers; ++k) {
        multiplier_index[r_curve[k].Id()] = first_multiplier + k;
    }
    Vector multipliers = ZeroVector(number_of_dofs);
    for (IndexType i = 0; i < multiplier_ids.size(); ++i) {
        multipliers[multiplier_index.at(multiplier_ids[i])] += multiplier_values[i];
    }

    // Largest degree P of the spaces, and shift order, by default P
    SizeType largest_degree = 0;
    for (IndexType part = 0; part < 3; ++part) {
        const GeometryType& r_surface = r_minus.GetGeometryPart(part).GetGeometryParent(0);
        largest_degree = std::max({
            largest_degree, r_surface.PolynomialDegree(0), r_surface.PolynomialDegree(1)
        });
    }
    const SizeType shift_order = GetProperties().Has(SHIFT_ORDER)
        ? static_cast<SizeType>(GetProperties()[SHIFT_ORDER])
        : largest_degree;

    // \gamma = \alpha_h+ / max(\alpha_h+ + \alpha_h-, tol)
    // Average operator: {{.}}_\gamma = \gamma (.)+ + (1 - \gamma) (.)-
    const double gamma =
        alpha_plus / std::max(alpha_plus + alpha_minus, std::numeric_limits<double>::epsilon());

    // [[\alpha_h]] = (\alpha_h- - \alpha_h+) n-
    const double alpha_jump = alpha_minus - alpha_plus;

    // Rows x and y over the local dofs: [[v_\alpha]], {{D B N}}_\gamma n- and
    // {{S(N)}}_\gamma
    Matrix jump_v_alpha = ZeroMatrix(2, number_of_dofs);
    Matrix viscous_traction = ZeroMatrix(2, number_of_dofs);
    Matrix shifted_velocity = ZeroMatrix(2, number_of_dofs);
    // {{S(N_p)}}_\gamma and {{S(N_p)}}_\gamma - {{N_p}}_\gamma
    Vector shifted_pressure_average = ZeroVector(number_of_dofs);
    Vector pressure_increment = ZeroVector(number_of_dofs);

    for (IndexType side = 0; side < 2; ++side) {
        const GeometryType& r_side = GetGeometry().GetGeometryPart(side);
        const std::vector<IndexType>& r_indices = r_side.GetValue(KNOT_SPAN_INDICES);
        const GeometryType& r_velocity_x = r_side.GetGeometryPart(GsbmStokesElement::VelocityXPart);
        const GeometryType& r_velocity_y = r_side.GetGeometryPart(GsbmStokesElement::VelocityYPart);
        const GeometryType& r_pressure = r_side.GetGeometryPart(GsbmStokesElement::PressurePart);
        const double jump_factor = side == 0 ? alpha_minus : -alpha_plus;
        const double average_factor = side == 0 ? 1.0 - gamma : gamma;

        // D B N n-: the traction without the pressure
        const Matrix D = CalculateConstitutiveMatrix(r_side, rCurrentProcessInfo);
        const Matrix traction = ComputeMultiplierViscousTraction(
            D, r_velocity_x.ShapeFunctionsLocalGradients()[0],
            r_velocity_y.ShapeFunctionsLocalGradients()[0], normal
        );
        const IndexType begin = first_dof[side][0];
        const IndexType end = begin + traction.size2();
        noalias(subrange(viscous_traction, 0, 2, begin, end)) = average_factor * traction;

        // N e_c and S(N) e_c for the dofs of each velocity component
        for (IndexType c = 0; c < 2; ++c) {
            const GeometryType& r_velocity = r_side.GetGeometryPart(c);
            const Matrix& r_N = r_velocity.ShapeFunctionsValues();
            const Vector shifted = ShiftedBoundaryUtilities::ComputeShiftedShapeFunctions(
                r_velocity, r_indices, distance, shift_order
            );
            for (IndexType i = 0; i < r_velocity.size(); ++i) {
                jump_v_alpha(c, first_dof[side][c] + i) = jump_factor * r_N(0, i);
                shifted_velocity(c, first_dof[side][c] + i) = average_factor * shifted[i];
            }
        }

        // Taylor increment of the pressure from the facet to the closest point
        const Matrix& r_N_p = r_pressure.ShapeFunctionsValues();
        const Vector shifted_pressure = ShiftedBoundaryUtilities::ComputeShiftedShapeFunctions(
            r_pressure, r_indices, distance, shift_order
        );
        const IndexType first_pressure = first_dof[side][GsbmStokesElement::PressurePart];
        for (IndexType i = 0; i < r_pressure.size(); ++i) {
            shifted_pressure_average[first_pressure + i] = average_factor * shifted_pressure[i];
            pressure_increment[first_pressure + i] =
                average_factor * (shifted_pressure[i] - r_N_p(0, i));
        }
    }

    // [[v_\alpha]] . n-, n- . {{S(N)}}_\gamma and \tau . {{S(N)}}_\gamma
    const Vector jump_normal =
        normal[0] * row(jump_v_alpha, 0) + normal[1] * row(jump_v_alpha, 1);
    const Vector shifted_normal =
        normal[0] * row(shifted_velocity, 0) + normal[1] * row(shifted_velocity, 1);
    const Vector shifted_tangent =
        tangent[0] * row(shifted_velocity, 0) + tangent[1] * row(shifted_velocity, 1);
    const double normal_data = normal[0] * velocity_data[0] + normal[1] * velocity_data[1];
    const double tangent_data = tangent[0] * velocity_data[0] + tangent[1] * velocity_data[1];

    // Viscous consistency: -([[v_\alpha]], {{D B u_h}}_\gamma n-)
    noalias(rLeftHandSideMatrix) -= weight * prod(trans(jump_v_alpha), viscous_traction);

    // Pressure: +([[v_\alpha]] . n-, \lambda - {{S p_h - p_h}}_\gamma)
    noalias(rLeftHandSideMatrix) += weight * outer_prod(jump_normal, multipliers);
    noalias(rLeftHandSideMatrix) -= weight * outer_prod(jump_normal, pressure_increment);

    // Tangential adjoint: -\theta [[\alpha_h]] ({{D B v_h}}_\gamma n- . \tau,
    //                                           ({{S(u_h)}}_\gamma - g) . \tau)
    const double theta =
        GetProperties().Has(ADJOINT_FACTOR) ? GetProperties()[ADJOINT_FACTOR] : -1.0;
    const Vector viscous_tangent =
        tangent[0] * row(viscous_traction, 0) + tangent[1] * row(viscous_traction, 1);
    noalias(rLeftHandSideMatrix) -=
        theta * weight * alpha_jump * outer_prod(viscous_tangent, shifted_tangent);
    noalias(rRightHandSideVector) -= theta * weight * alpha_jump * tangent_data * viscous_tangent;

    // Tangential penalty, if \chi = PENALTY_FACTOR > 0 (none by default):
    // (\chi (P + 1)^2 \mu / h |[[\alpha_h]]| {{S(v_h)}}_\gamma . \tau,
    //  ({{S(u_h)}}_\gamma - g) . \tau), and along n- with PENALTY_ON_NORMAL
    const double chi = GetProperties().Has(PENALTY_FACTOR) ? GetProperties()[PENALTY_FACTOR] : 0.0;
    if (chi > 0.0) {
        // Characteristic size h: the smallest knot span size of both sides
        const Vector& r_sizes_minus = r_minus.GetValue(KNOT_SPAN_SIZES);
        const Vector& r_sizes_plus = r_plus.GetValue(KNOT_SPAN_SIZES);
        const double h = std::min(
            *std::min_element(r_sizes_minus.begin(), r_sizes_minus.end()),
            *std::min_element(r_sizes_plus.begin(), r_sizes_plus.end())
        );

        const double penalty = weight * chi * (largest_degree + 1) * (largest_degree + 1)
            * GetProperties()[DYNAMIC_VISCOSITY] / h * std::abs(alpha_jump);
        noalias(rLeftHandSideMatrix) += penalty * outer_prod(shifted_tangent, shifted_tangent);
        noalias(rRightHandSideVector) += penalty * tangent_data * shifted_tangent;

        // Augmented Lagrangian: the same penalty along n-
        if (GetProperties().Has(PENALTY_ON_NORMAL) && GetProperties()[PENALTY_ON_NORMAL]) {
            noalias(rLeftHandSideMatrix) += penalty * outer_prod(shifted_normal, shifted_normal);
            noalias(rRightHandSideVector) += penalty * normal_data * shifted_normal;
        }
    }

    // Pinned by default, for a divergence-free velocity
    const bool is_pinned =
        !GetProperties().Has(PIN_MULTIPLIERS) || GetProperties()[PIN_MULTIPLIERS];

    // Multiplier rows, mean-free if pinned, with c_k the share of \mu_k:
    // ([[\alpha_h]] (\mu_k - c_k), n- . ({{S(u_h)}}_\gamma - g)) = 0
    Vector tested_multipliers = multipliers;
    if (is_pinned) {
        const Vector shares = ComputeMultiplierShares(r_curve, number_of_multipliers);
        for (IndexType k = 0; k < number_of_multipliers; ++k) {
            tested_multipliers[first_multiplier + k] -= shares[k];
        }
    }
    noalias(rLeftHandSideMatrix) +=
        weight * alpha_jump * outer_prod(tested_multipliers, shifted_normal);
    noalias(rRightHandSideVector) += weight * alpha_jump * normal_data * tested_multipliers;

    // Pinned, they sum to zero, so the last follows from the others. It ties
    // \lambda to the pressure instead, in place of their constant combination,
    // which repeats the continuity equations of the divergence-free velocity:
    // (|[[\alpha_h]]|, \lambda - {{S p_h}}_\gamma) = 0
    if (is_pinned) {
        const IndexType last_multiplier = first_multiplier + number_of_multipliers - 1;
        noalias(row(rLeftHandSideMatrix, last_multiplier)) =
            weight * std::abs(alpha_jump) * (multipliers - shifted_pressure_average);
        rRightHandSideVector[last_multiplier] = 0.0;
    }

    // Residual form
    Vector current_values(number_of_dofs);
    IndexType dof = 0;
    for (IndexType side = 0; side < 2; ++side) {
        for (IndexType part = 0; part < 3; ++part) {
            const auto& r_variable = GsbmStokesElement::GetFieldVariable(part);
            for (const auto& r_node : GetGeometry().GetGeometryPart(side).GetGeometryPart(part)) {
                current_values[dof++] = r_node.FastGetSolutionStepValue(r_variable);
            }
        }
    }
    for (IndexType k = 0; k < number_of_multipliers; ++k) {
        current_values[dof++] = r_curve[k].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
    }
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, current_values);

    KRATOS_CATCH("")
}

void GsbmStokesConditionDirichletMultiplier::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmStokesConditionDirichletMultiplier::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmStokesConditionDirichletMultiplier::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rResult.clear();
    for (IndexType side = 0; side < 2; ++side) {
        for (IndexType part = 0; part < 3; ++part) {
            const auto& r_variable = GsbmStokesElement::GetFieldVariable(part);
            for (const auto& r_node : GetGeometry().GetGeometryPart(side).GetGeometryPart(part)) {
                rResult.push_back(r_node.GetDof(r_variable).EquationId());
            }
        }
    }
    const MultiplierCurveType& r_curve = GetMultiplierCurve(GetGeometry());
    for (IndexType k = 0; k < NumberOfMultipliers(); ++k) {
        rResult.push_back(r_curve[k].GetDof(SCALAR_LAGRANGE_MULTIPLIER).EquationId());
    }
}

void GsbmStokesConditionDirichletMultiplier::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rConditionDofList.clear();
    for (IndexType side = 0; side < 2; ++side) {
        for (IndexType part = 0; part < 3; ++part) {
            const auto& r_variable = GsbmStokesElement::GetFieldVariable(part);
            for (const auto& r_node : GetGeometry().GetGeometryPart(side).GetGeometryPart(part)) {
                rConditionDofList.push_back(r_node.pGetDof(r_variable));
            }
        }
    }
    const MultiplierCurveType& r_curve = GetMultiplierCurve(GetGeometry());
    for (IndexType k = 0; k < NumberOfMultipliers(); ++k) {
        rConditionDofList.push_back(r_curve[k].pGetDof(SCALAR_LAGRANGE_MULTIPLIER));
    }
}

int GsbmStokesConditionDirichletMultiplier::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(GetGeometry().NumberOfGeometryParts() != 3)
        << Info() << " needs three geometry parts: both facet sides and the multiplier curve."
        << std::endl;
    for (IndexType side = 0; side < 2; ++side) {
        KRATOS_ERROR_IF(GetGeometry().GetGeometryPart(side).NumberOfGeometryParts() != 3)
            << Info() << " needs facet sides with three parts, the spaces of u_x, u_y and p."
            << std::endl;
    }
    const MultiplierCurveType& r_curve = GetMultiplierCurve(GetGeometry());
    for (IndexType k = 0; k < NumberOfMultipliers(); ++k) {
        KRATOS_ERROR_IF_NOT(r_curve[k].HasDofFor(SCALAR_LAGRANGE_MULTIPLIER))
            << Info() << ": the multiplier node " << r_curve[k].Id()
            << " has no SCALAR_LAGRANGE_MULTIPLIER dof." << std::endl;
    }
    KRATOS_ERROR_IF(GetProperties().Has(SHIFT_ORDER) && GetProperties()[SHIFT_ORDER] < 0)
        << Info() << ": SHIFT_ORDER must not be negative." << std::endl;
    KRATOS_ERROR_IF(GetProperties().Has(PENALTY_FACTOR) && GetProperties()[PENALTY_FACTOR] > 0.0
        && !GetProperties().Has(DYNAMIC_VISCOSITY)) << Info()
        << ": PENALTY_FACTOR needs DYNAMIC_VISCOSITY in the properties." << std::endl;
    KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW)
        || GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << Info() << ": no CONSTITUTIVE_LAW in the properties." << std::endl;

    const GeometryType& r_velocity_x =
        GetGeometry().GetGeometryPart(0).GetGeometryPart(GsbmStokesElement::VelocityXPart);
    return GetProperties()[CONSTITUTIVE_LAW]->Check(
        GetProperties(), r_velocity_x, rCurrentProcessInfo
    );

    KRATOS_CATCH("")
}

void GsbmStokesConditionDirichletMultiplier::InitializeMaterial()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW)
        || GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << Info() << ": no CONSTITUTIVE_LAW in the properties." << std::endl;

    // The constitutive law of the properties, at the point of u_x's space
    const GeometryType& r_velocity_x =
        GetGeometry().GetGeometryPart(0).GetGeometryPart(GsbmStokesElement::VelocityXPart);
    mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    mpConstitutiveLaw->InitializeMaterial(
        GetProperties(), r_velocity_x, row(r_velocity_x.ShapeFunctionsValues(), 0)
    );

    KRATOS_CATCH("")
}

Matrix GsbmStokesConditionDirichletMultiplier::CalculateConstitutiveMatrix(
    const GeometryType& rSide,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // Current Voigt strain of the side's velocity
    std::array<std::array<double, 2>, 2> velocity_gradients{};
    for (IndexType c = 0; c < 2; ++c) {
        const GeometryType& r_velocity = rSide.GetGeometryPart(c);
        const Matrix& r_DN = r_velocity.ShapeFunctionsLocalGradients()[0];
        const auto& r_variable = GsbmStokesElement::GetFieldVariable(c);
        for (IndexType i = 0; i < r_velocity.size(); ++i) {
            const double value = r_velocity[i].FastGetSolutionStepValue(r_variable);
            velocity_gradients[c][0] += r_DN(i, 0) * value;
            velocity_gradients[c][1] += r_DN(i, 1) * value;
        }
    }
    Vector strain(3);
    strain[0] = velocity_gradients[0][0];
    strain[1] = velocity_gradients[1][1];
    strain[2] = velocity_gradients[0][1] + velocity_gradients[1][0];

    // Stress and D of the constitutive law at that strain
    const GeometryType& r_velocity_x = rSide.GetGeometryPart(GsbmStokesElement::VelocityXPart);
    ConstitutiveLaw::Parameters values(r_velocity_x, GetProperties(), rCurrentProcessInfo);
    Flags& r_options = values.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    Vector stress = ZeroVector(3);
    Matrix D = ZeroMatrix(3, 3);
    values.SetStrainVector(strain);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(D);
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(values);

    return D;
}

GsbmStokesConditionDirichletMultiplier::SizeType
    GsbmStokesConditionDirichletMultiplier::NumberOfMultipliers() const
{
    const MultiplierCurveType& r_curve = GetMultiplierCurve(GetGeometry());
    return r_curve.size() - r_curve.PolynomialDegree(0);
}

} // namespace Kratos
