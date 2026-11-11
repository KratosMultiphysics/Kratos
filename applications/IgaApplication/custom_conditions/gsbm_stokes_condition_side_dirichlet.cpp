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

// External includes

// Project includes
#include "custom_conditions/gsbm_stokes_condition_side_dirichlet.h"
#include "custom_elements/gsbm_stokes_element.h"
#include "iga_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

namespace
{

/**
 * @brief Viscous traction D B N n of the velocity dofs of a side point, from
 *        the constitutive matrix and the velocity gradients
 *
 * @param rD Constitutive matrix D of the Voigt strain [xx, yy, xy].
 * @param rGradientsX Gradients of the shape functions of u_x's space.
 * @param rGradientsY Gradients of the shape functions of u_y's space.
 * @param rNormal Unit normal n of the traction.
 * @return Rows x and y of the traction, over the dofs of u_x, then of u_y.
 */
Matrix ComputeSideViscousTraction(
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

    // (D B) n, with the Voigt stress [\sigma_xx, \sigma_yy, \sigma_xy]
    Matrix traction(2, n_u);
    for (std::size_t k = 0; k < n_u; ++k) {
        traction(0, k) = stress(0, k) * rNormal[0] + stress(2, k) * rNormal[1];
        traction(1, k) = stress(2, k) * rNormal[0] + stress(1, k) * rNormal[1];
    }

    return traction;
}

} // namespace

void GsbmStokesConditionSideDirichlet::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}

void GsbmStokesConditionSideDirichlet::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw) << Info()
        << " has no constitutive law. Call Initialize first." << std::endl;

    // First local dof of each field: u_x, u_y, then p
    std::array<IndexType, 3> first_dof;
    SizeType number_of_dofs = 0;
    for (IndexType part = 0; part < 3; ++part) {
        first_dof[part] = number_of_dofs;
        number_of_dofs += GetGeometry().GetGeometryPart(part).size();
    }

    // Local system over the dofs of the side point
    if (rLeftHandSideMatrix.size1() != number_of_dofs ||
        rLeftHandSideMatrix.size2() != number_of_dofs) {
        rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs);
    if (rRightHandSideVector.size() != number_of_dofs) {
        rRightHandSideVector.resize(number_of_dofs, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_dofs);

    const GeometryType& r_velocity_x =
        GetGeometry().GetGeometryPart(GsbmStokesElement::VelocityXPart);
    const GeometryType& r_velocity_y =
        GetGeometry().GetGeometryPart(GsbmStokesElement::VelocityYPart);
    const GeometryType& r_pressure = GetGeometry().GetGeometryPart(GsbmStokesElement::PressurePart);

    // Side measure, and the outward normal n
    const double weight = r_velocity_x.IntegrationPoints()[0].Weight()
                        * r_velocity_x.DeterminantOfJacobian(0);
    const array_1d<double, 3> normal = r_velocity_x.UnitNormal(0);

    // Imposed velocity g of the condition
    const array_1d<double, 3>& r_imposed_velocity = GetValue(VELOCITY);
    Vector velocity_data(2);
    velocity_data[0] = r_imposed_velocity[0];
    velocity_data[1] = r_imposed_velocity[1];

    // Largest degree P of the spaces
    SizeType largest_degree = 0;
    for (IndexType part = 0; part < 3; ++part) {
        const GeometryType& r_surface = GetGeometry().GetGeometryPart(part).GetGeometryParent(0);
        largest_degree = std::max({
            largest_degree, r_surface.PolynomialDegree(0), r_surface.PolynomialDegree(1)
        });
    }

    // Rows x and y over the local dofs: N e_c and D B N n
    Matrix velocity = ZeroMatrix(2, number_of_dofs);
    Matrix viscous_traction = ZeroMatrix(2, number_of_dofs);
    // N_p
    Vector pressure = ZeroVector(number_of_dofs);

    const Matrix D = CalculateConstitutiveMatrix(rCurrentProcessInfo);
    const Matrix traction = ComputeSideViscousTraction(
        D, r_velocity_x.ShapeFunctionsLocalGradients()[0],
        r_velocity_y.ShapeFunctionsLocalGradients()[0], normal
    );
    noalias(subrange(viscous_traction, 0, 2, 0, traction.size2())) = traction;

    for (IndexType c = 0; c < 2; ++c) {
        const Matrix& r_N = GetGeometry().GetGeometryPart(c).ShapeFunctionsValues();
        for (IndexType i = 0; i < r_N.size2(); ++i) {
            velocity(c, first_dof[c] + i) = r_N(0, i);
        }
    }
    const Matrix& r_N_p = r_pressure.ShapeFunctionsValues();
    for (IndexType i = 0; i < r_pressure.size(); ++i) {
        pressure[first_dof[GsbmStokesElement::PressurePart] + i] = r_N_p(0, i);
    }

    // \sigma_h n = D B u_h n - p_h n
    Matrix stress_traction = viscous_traction;
    for (IndexType m = 0; m < 2; ++m) {
        noalias(row(stress_traction, m)) -= normal[m] * pressure;
    }

    // Consistency: -(v_h, \sigma_h n)
    noalias(rLeftHandSideMatrix) -= weight * prod(trans(velocity), stress_traction);

    // Viscous adjoint: -\theta (D B v_h n, u_h - g),
    // \theta = ADJOINT_FACTOR: -1 skew-symmetric (default), 1 symmetric
    const double theta =
        GetProperties().Has(ADJOINT_FACTOR) ? GetProperties()[ADJOINT_FACTOR] : -1.0;
    noalias(rLeftHandSideMatrix) -= theta * weight * prod(trans(viscous_traction), velocity);
    noalias(rRightHandSideVector) -= theta * weight * prod(trans(viscous_traction), velocity_data);

    // q term, unless ADJOINT_ON_PRESSURE is false: -(q_h, n . (u_h - g))
    if (!GetProperties().Has(ADJOINT_ON_PRESSURE) || GetProperties()[ADJOINT_ON_PRESSURE]) {
        const Vector normal_velocity = normal[0] * row(velocity, 0) + normal[1] * row(velocity, 1);
        const double normal_data = normal[0] * velocity_data[0] + normal[1] * velocity_data[1];
        noalias(rLeftHandSideMatrix) -= weight * outer_prod(pressure, normal_velocity);
        noalias(rRightHandSideVector) -= weight * normal_data * pressure;
    }

    // Penalty, if \chi = PENALTY_FACTOR > 0 (none by default):
    // (\chi (P + 1)^2 \mu / h v_h, u_h - g), h the smallest knot span size
    const double chi = GetProperties().Has(PENALTY_FACTOR) ? GetProperties()[PENALTY_FACTOR] : 0.0;
    if (chi > 0.0) {
        const Vector& r_sizes = GetValue(KNOT_SPAN_SIZES);
        const double h = *std::min_element(r_sizes.begin(), r_sizes.end());
        const double penalty = weight * chi * (largest_degree + 1) * (largest_degree + 1)
            * GetProperties()[DYNAMIC_VISCOSITY] / h;
        noalias(rLeftHandSideMatrix) += penalty * prod(trans(velocity), velocity);
        noalias(rRightHandSideVector) += penalty * prod(trans(velocity), velocity_data);
    }

    // Residual form
    Vector current_values(number_of_dofs);
    IndexType dof = 0;
    for (IndexType part = 0; part < 3; ++part) {
        const auto& r_variable = GsbmStokesElement::GetFieldVariable(part);
        for (const auto& r_node : GetGeometry().GetGeometryPart(part)) {
            current_values[dof++] = r_node.FastGetSolutionStepValue(r_variable);
        }
    }
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, current_values);

    KRATOS_CATCH("")
}

void GsbmStokesConditionSideDirichlet::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmStokesConditionSideDirichlet::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmStokesConditionSideDirichlet::EquationIdVector(
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

void GsbmStokesConditionSideDirichlet::GetDofList(
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

int GsbmStokesConditionSideDirichlet::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(GetGeometry().NumberOfGeometryParts() != 3) << Info()
        << " needs a side point with three parts, the spaces of u_x, u_y and p. Create it "
        << "with \"sides\" in IgaModelerGsbm." << std::endl;
    KRATOS_ERROR_IF(GetProperties().Has(PENALTY_FACTOR) && GetProperties()[PENALTY_FACTOR] > 0.0
        && !GetProperties().Has(DYNAMIC_VISCOSITY)) << Info()
        << ": PENALTY_FACTOR needs DYNAMIC_VISCOSITY in the properties." << std::endl;
    KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW)
        || GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << Info() << ": no CONSTITUTIVE_LAW in the properties." << std::endl;

    const GeometryType& r_velocity_x =
        GetGeometry().GetGeometryPart(GsbmStokesElement::VelocityXPart);
    return GetProperties()[CONSTITUTIVE_LAW]->Check(
        GetProperties(), r_velocity_x, rCurrentProcessInfo
    );

    KRATOS_CATCH("")
}

void GsbmStokesConditionSideDirichlet::InitializeMaterial()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW)
        || GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << Info() << ": no CONSTITUTIVE_LAW in the properties." << std::endl;

    // The constitutive law of the properties, at the point of u_x's space
    const GeometryType& r_velocity_x =
        GetGeometry().GetGeometryPart(GsbmStokesElement::VelocityXPart);
    mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    mpConstitutiveLaw->InitializeMaterial(
        GetProperties(), r_velocity_x, row(r_velocity_x.ShapeFunctionsValues(), 0)
    );

    KRATOS_CATCH("")
}

Matrix GsbmStokesConditionSideDirichlet::CalculateConstitutiveMatrix(
    const ProcessInfo& rCurrentProcessInfo) const
{
    // Current Voigt strain of the velocity
    std::array<std::array<double, 2>, 2> velocity_gradients{};
    for (IndexType c = 0; c < 2; ++c) {
        const GeometryType& r_velocity = GetGeometry().GetGeometryPart(c);
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
    const GeometryType& r_velocity_x =
        GetGeometry().GetGeometryPart(GsbmStokesElement::VelocityXPart);
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
