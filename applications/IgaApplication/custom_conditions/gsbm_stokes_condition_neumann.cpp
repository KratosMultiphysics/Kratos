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
#include <limits>

// External includes

// Project includes
#include "custom_conditions/gsbm_stokes_condition_neumann.h"
#include "custom_elements/gsbm_stokes_element.h"
#include "custom_utilities/shifted_boundary_utilities.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

namespace
{

/**
 * @brief Traction \sigma(N) \nu of the dofs of a facet side, from the
 *        constitutive matrix and the side's velocity gradients and pressure
 *        shape functions
 *
 * @param rD Constitutive matrix D of the Voigt strain [xx, yy, xy].
 * @param rGradientsX Gradients of the shape functions of u_x's space.
 * @param rGradientsY Gradients of the shape functions of u_y's space.
 * @param rPressureValues Shape functions of the pressure space.
 * @param rNormal Unit normal \nu of the traction.
 * @return Rows x and y of the traction, over the side's dofs: u_x, u_y, p.
 */
Matrix ComputeStokesTraction(
    const Matrix& rD,
    const Matrix& rGradientsX,
    const Matrix& rGradientsY,
    const Vector& rPressureValues,
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
    Matrix traction(2, n_u + rPressureValues.size());
    for (std::size_t k = 0; k < n_u; ++k) {
        traction(0, k) = stress(0, k) * rNormal[0] + stress(2, k) * rNormal[1];
        traction(1, k) = stress(2, k) * rNormal[0] + stress(1, k) * rNormal[1];
    }

    // -p \nu
    for (std::size_t i = 0; i < rPressureValues.size(); ++i) {
        traction(0, n_u + i) = -rPressureValues[i] * rNormal[0];
        traction(1, n_u + i) = -rPressureValues[i] * rNormal[1];
    }

    return traction;
}

} // namespace

void GsbmStokesConditionNeumann::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}

void GsbmStokesConditionNeumann::CalculateLocalSystem(
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

    // Facet measure, and the normal n- pointing from the - to the + side
    const GeometryType& r_reference = r_minus.GetGeometryPart(GsbmStokesElement::VelocityXPart);
    const double weight = r_reference.IntegrationPoints()[0].Weight()
                        * r_reference.DeterminantOfJacobian(0);
    const array_1d<double, 3> normal = r_reference.UnitNormal(0);

    // Distance d and outward normal n at the closest point
    const Node& r_projection_node = GetValue(NEIGHBOUR_NODES)[0];
    const array_1d<double, 3> distance =
        r_projection_node.Coordinates() - r_reference.Center().Coordinates();
    const array_1d<double, 3>& r_skin_normal = r_projection_node.GetValue(NORMAL);

    // t_N = -p_D n + NORMAL_STRESS, as in SupportPressureCondition
    Vector traction_data(2);
    for (IndexType m = 0; m < 2; ++m) {
        traction_data[m] = -r_projection_node.GetValue(PRESSURE) * r_skin_normal[m];
    }
    if (r_projection_node.Has(NORMAL_STRESS)) {
        const Vector& r_normal_stress = r_projection_node.GetValue(NORMAL_STRESS);
        KRATOS_ERROR_IF(r_normal_stress.size() != 2) << Info()
            << ": NORMAL_STRESS of the projection node has size " << r_normal_stress.size()
            << ", expected 2." << std::endl;
        noalias(traction_data) += r_normal_stress;
    }

    // Shift order m: SHIFT_ORDER, the largest degree of the spaces by default
    SizeType shift_order = 0;
    if (GetProperties().Has(SHIFT_ORDER)) {
        shift_order = static_cast<SizeType>(GetProperties()[SHIFT_ORDER]);
    } else {
        for (IndexType part = 0; part < 3; ++part) {
            const GeometryType& r_surface = r_minus.GetGeometryPart(part).GetGeometryParent(0);
            shift_order = std::max({
                shift_order, r_surface.PolynomialDegree(0), r_surface.PolynomialDegree(1)
            });
        }
    }

    // \gamma = \alpha_h+ / max(\alpha_h+ + \alpha_h-, tol)
    // Average operator: {{.}}_\gamma = \gamma (.)+ + (1 - \gamma) (.)-
    const double gamma =
        alpha_plus / std::max(alpha_plus + alpha_minus, std::numeric_limits<double>::epsilon());

    // Rows x and y over the local dofs: [[v_\alpha]] = \alpha_h- v_h- -
    // \alpha_h+ v_h+, {{\sigma_h}}_\gamma n- and {{S_d^m[\sigma_h]}}_\gamma n
    Matrix jump_v_alpha = ZeroMatrix(2, number_of_dofs);
    Matrix average_traction = ZeroMatrix(2, number_of_dofs);
    Matrix average_shifted_traction = ZeroMatrix(2, number_of_dofs);

    for (IndexType side = 0; side < 2; ++side) {
        const GeometryType& r_side = GetGeometry().GetGeometryPart(side);
        const std::vector<IndexType>& r_indices = r_side.GetValue(KNOT_SPAN_INDICES);
        const GeometryType& r_velocity_x = r_side.GetGeometryPart(GsbmStokesElement::VelocityXPart);
        const GeometryType& r_velocity_y = r_side.GetGeometryPart(GsbmStokesElement::VelocityYPart);
        const GeometryType& r_pressure = r_side.GetGeometryPart(GsbmStokesElement::PressurePart);
        const double jump_factor = side == 0 ? alpha_minus : -alpha_plus;
        const double average_factor = side == 0 ? 1.0 - gamma : gamma;

        // \sigma(N) n- at \tilde{x}, and S_d^m[\sigma(N)] n at \tilde{x} + d
        const Matrix D = CalculateConstitutiveMatrix(r_side, rCurrentProcessInfo);
        const Vector pressure_values = row(r_pressure.ShapeFunctionsValues(), 0);
        const Matrix traction = ComputeStokesTraction(
            D, r_velocity_x.ShapeFunctionsLocalGradients()[0],
            r_velocity_y.ShapeFunctionsLocalGradients()[0], pressure_values, normal
        );
        const Matrix shifted_traction = ComputeStokesTraction(
            D,
            ShiftedBoundaryUtilities::ComputeShiftedGradients(
                r_velocity_x, r_indices, distance, shift_order
            ),
            ShiftedBoundaryUtilities::ComputeShiftedGradients(
                r_velocity_y, r_indices, distance, shift_order
            ),
            ShiftedBoundaryUtilities::ComputeShiftedShapeFunctions(
                r_pressure, r_indices, distance, shift_order
            ),
            r_skin_normal
        );

        const IndexType begin = first_dof[side][0];
        const IndexType end = begin + traction.size2();
        noalias(subrange(average_traction, 0, 2, begin, end)) = average_factor * traction;
        noalias(subrange(average_shifted_traction, 0, 2, begin, end)) =
            average_factor * shifted_traction;

        // N e_c for the dofs of each velocity component
        for (IndexType c = 0; c < 2; ++c) {
            const Matrix& r_N = r_side.GetGeometryPart(c).ShapeFunctionsValues();
            for (IndexType i = 0; i < r_N.size2(); ++i) {
                jump_v_alpha(c, first_dof[side][c] + i) = jump_factor * r_N(0, i);
            }
        }
    }

    // Consistency: -([[v_\alpha]], {{\sigma_h}}_\gamma n-)
    noalias(rLeftHandSideMatrix) -= weight * prod(trans(jump_v_alpha), average_traction);

    // Shifted traction:
    // ([[v_\alpha]], (n- . n) ({{S_d^m[\sigma_h]}}_\gamma n - t_N))
    const double normal_weight = weight * inner_prod(normal, r_skin_normal);
    noalias(rLeftHandSideMatrix) +=
        normal_weight * prod(trans(jump_v_alpha), average_shifted_traction);
    noalias(rRightHandSideVector) += normal_weight * prod(trans(jump_v_alpha), traction_data);

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

void GsbmStokesConditionNeumann::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmStokesConditionNeumann::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmStokesConditionNeumann::EquationIdVector(
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

void GsbmStokesConditionNeumann::GetDofList(
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

int GsbmStokesConditionNeumann::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    for (IndexType side = 0; side < 2; ++side) {
        KRATOS_ERROR_IF(GetGeometry().GetGeometryPart(side).NumberOfGeometryParts() != 3)
            << Info() << " needs facet sides with three parts, the spaces of u_x, u_y and p."
            << std::endl;
    }
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

void GsbmStokesConditionNeumann::InitializeMaterial()
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

Matrix GsbmStokesConditionNeumann::CalculateConstitutiveMatrix(
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

} // namespace Kratos
