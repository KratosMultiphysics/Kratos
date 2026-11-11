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

// External includes

// Project includes
#include "custom_conditions/gsbm_laplacian_condition_ghost_penalty.h"
#include "custom_conditions/gsbm_stokes_condition_ghost_penalty.h"
#include "custom_elements/gsbm_stokes_element.h"
#include "custom_utilities/shifted_boundary_utilities.h"
#include "iga_application_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace
{

/// Whether the pressure ghost penalty acts on \nabla . u_h, false by default
bool IsPressureGhostPenaltyOnDivergence(const Properties& rProperties)
{
    return rProperties.Has(PRESSURE_GHOST_PENALTY_ON_DIVERGENCE)
        && rProperties[PRESSURE_GHOST_PENALTY_ON_DIVERGENCE];
}

} // namespace

void GsbmStokesConditionGhostPenalty::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr IndexType pressure_part = GsbmStokesElement::PressurePart;

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

    if (!GsbmLaplacianConditionGhostPenalty::IsPenalized(*this)) {
        return;  // skip early if \beta_2 = 0 (zero contribution)
    }

    // Facet measure, and the normal n- pointing from the - to the + side
    const GeometryType& r_reference = r_minus.GetGeometryPart(GsbmStokesElement::VelocityXPart);
    const double weight = r_reference.IntegrationPoints()[0].Weight()
                        * r_reference.DeterminantOfJacobian(0);
    const array_1d<double, 3> normal = r_reference.UnitNormal(0);

    // Characteristic size h: the smallest knot span size of both sides
    const Vector& r_sizes_minus = r_minus.GetValue(KNOT_SPAN_SIZES);
    const Vector& r_sizes_plus = r_plus.GetValue(KNOT_SPAN_SIZES);
    const double h = std::min(
        *std::min_element(r_sizes_minus.begin(), r_sizes_minus.end()),
        *std::min_element(r_sizes_plus.begin(), r_sizes_plus.end())
    );

    // Degree k of a field's space across the facet: only D_n^k jumps there
    const IndexType normal_direction = std::abs(normal[0]) > std::abs(normal[1]) ? 0 : 1;
    const auto get_degree_across = [&](const IndexType Part) -> SizeType {
        return r_minus.GetGeometryPart(Part).GetGeometryParent(0).PolynomialDegree(
            normal_direction
        );
    };

    // w / ((k!)^2 (2k + 1))
    const auto get_normalized_weight = [&](const SizeType Degree) {
        return weight / (std::pow(MathUtils<double>::Factorial(Degree), 2) * (2 * Degree + 1));
    };

    // Velocity: (\kappa_u \mu h^(2p_c-1) [[D_n^{p_c} v_c]], [[D_n^{p_c} u_c]])
    const double viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    for (IndexType c = 0; c < 2; ++c) {
        const SizeType p_c = get_degree_across(c);

        // [[D_n^{p_c} v_c]] . n- = D_n^{p_c} v_c- - D_n^{p_c} v_c+
        Vector jump_velocity = ZeroVector(number_of_dofs);
        for (IndexType side = 0; side < 2; ++side) {
            const GeometryType& r_side = GetGeometry().GetGeometryPart(side);
            const GeometryType& r_velocity = r_side.GetGeometryPart(c);
            const IndexType begin = first_dof[side][c];
            noalias(subrange(jump_velocity, begin, begin + r_velocity.size())) =
                (side == 0 ? 1.0 : -1.0) * ShiftedBoundaryUtilities::ComputeDirectionalDerivatives(
                    r_velocity, r_side.GetValue(KNOT_SPAN_INDICES), normal, p_c
                );
        }

        const double factor = get_normalized_weight(p_c) * GetProperties()[GHOST_PENALTY_FACTOR]
            * viscosity * std::pow(h, 2 * p_c - 1);
        noalias(rLeftHandSideMatrix) += factor * outer_prod(jump_velocity, jump_velocity);
    }

    // Pressure, with k its degree across the facet: [[D_n^k q]] . n-
    const SizeType k = get_degree_across(pressure_part);
    Vector jump_pressure = ZeroVector(number_of_dofs);
    for (IndexType side = 0; side < 2; ++side) {
        const GeometryType& r_side = GetGeometry().GetGeometryPart(side);
        const GeometryType& r_pressure = r_side.GetGeometryPart(pressure_part);
        const IndexType begin = first_dof[side][pressure_part];
        noalias(subrange(jump_pressure, begin, begin + r_pressure.size())) =
            (side == 0 ? 1.0 : -1.0) * ShiftedBoundaryUtilities::ComputeDirectionalDerivatives(
                r_pressure, r_side.GetValue(KNOT_SPAN_INDICES), normal, k
            );
    }

    // s(a, b) = (\kappa_p h^(2k + 1) [[D_n^k a]], [[D_n^k b]])
    const double pressure_factor = get_normalized_weight(k)
        * GetProperties()[PRESSURE_GHOST_PENALTY_FACTOR] * std::pow(h, 2 * k + 1);

    if (!IsPressureGhostPenaltyOnDivergence(GetProperties())) {
        // s(q, p) / (\mu + \beta_div), as in the paper
        const double grad_div_factor =
            GetProperties().Has(GRAD_DIV_FACTOR) ? GetProperties()[GRAD_DIV_FACTOR] : 0.0;
        noalias(rLeftHandSideMatrix) += pressure_factor / (viscosity + grad_div_factor)
            * outer_prod(jump_pressure, jump_pressure);
    } else {
        // [[D_n^k \nabla . v]] . n- = \sum_c [[D_n^k \partial_c v_c]] . n-
        Vector jump_divergence = ZeroVector(number_of_dofs);
        for (IndexType side = 0; side < 2; ++side) {
            const GeometryType& r_side = GetGeometry().GetGeometryPart(side);
            for (IndexType c = 0; c < 2; ++c) {
                const GeometryType& r_velocity = r_side.GetGeometryPart(c);
                const IndexType begin = first_dof[side][c];
                noalias(subrange(jump_divergence, begin, begin + r_velocity.size())) =
                    (side == 0 ? 1.0 : -1.0)
                    * ShiftedBoundaryUtilities::ComputeDirectionalDerivatives(
                        r_velocity, r_side.GetValue(KNOT_SPAN_INDICES), normal, k, c
                    );
            }
        }

        // +s(q, \nabla . u) and -s(\nabla . v, p), as +(q, \nabla . u) and
        // -(p, \nabla . v) in GsbmStokesElement
        noalias(rLeftHandSideMatrix) +=
            pressure_factor * outer_prod(jump_pressure, jump_divergence);
        noalias(rLeftHandSideMatrix) -=
            pressure_factor * outer_prod(jump_divergence, jump_pressure);
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
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, current_values);

    KRATOS_CATCH("")
}

void GsbmStokesConditionGhostPenalty::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmStokesConditionGhostPenalty::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmStokesConditionGhostPenalty::EquationIdVector(
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
}

void GsbmStokesConditionGhostPenalty::GetDofList(
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
}

int GsbmStokesConditionGhostPenalty::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    for (IndexType side = 0; side < 2; ++side) {
        KRATOS_ERROR_IF(GetGeometry().GetGeometryPart(side).NumberOfGeometryParts() != 3)
            << Info() << " needs facet sides with three parts, the spaces of u_x, u_y and p."
            << std::endl;
    }
    KRATOS_ERROR_IF_NOT(GetProperties().Has(GHOST_PENALTY_FACTOR))
        << Info() << ": no GHOST_PENALTY_FACTOR in the properties." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PRESSURE_GHOST_PENALTY_FACTOR))
        << Info() << ": no PRESSURE_GHOST_PENALTY_FACTOR in the properties." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(DYNAMIC_VISCOSITY))
        << Info() << ": no DYNAMIC_VISCOSITY in the properties." << std::endl;

    // u_x in S^(k_x+1, k_y) and u_y in S^(k_x, k_y+1), with p in S^(k_x, k_y)
    if (IsPressureGhostPenaltyOnDivergence(GetProperties())) {
        const GeometryType& r_side = GetGeometry().GetGeometryPart(0);
        const auto get_degree = [&](const IndexType Part, const IndexType Direction) {
            return r_side.GetGeometryPart(Part).GetGeometryParent(0).PolynomialDegree(Direction);
        };
        for (IndexType c = 0; c < 2; ++c) {
            for (IndexType direction = 0; direction < 2; ++direction) {
                const SizeType expected_degree =
                    get_degree(GsbmStokesElement::PressurePart, direction) + (c == direction);
                KRATOS_ERROR_IF(get_degree(c, direction) != expected_degree) << Info()
                    << ": PRESSURE_GHOST_PENALTY_ON_DIVERGENCE needs divergence-conforming "
                    << "spaces, u_x in S^(k_x+1, k_y) and u_y in S^(k_x, k_y+1) with p in "
                    << "S^(k_x, k_y)." << std::endl;
            }
        }
    }

    return 0;

    KRATOS_CATCH("")
}

} // namespace Kratos
