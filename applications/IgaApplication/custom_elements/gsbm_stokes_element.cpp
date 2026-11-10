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
#include "custom_elements/gsbm_stokes_element.h"
#include "iga_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

namespace Kratos
{

namespace
{

/// The factor of an optional term in the properties, 0 if it is not given
double GetOptionalFactor(const Properties& rProperties, const Variable<double>& rVariable)
{
    return rProperties.Has(rVariable) ? rProperties[rVariable] : 0.0;
}

} // namespace

const Variable<double>& GsbmStokesElement::GetFieldVariable(const IndexType Part)
{
    const std::array<const Variable<double>*, 3> variables{&VELOCITY_X, &VELOCITY_Y, &PRESSURE};
    return *variables[Part];
}

void GsbmStokesElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}

void GsbmStokesElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw) << Info()
        << " has no constitutive law. Call Initialize first." << std::endl;

    // One spline space per field: u_x, u_y and p
    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    const GeometryType& r_velocity_y = GetGeometry().GetGeometryPart(VelocityYPart);
    const GeometryType& r_pressure = GetGeometry().GetGeometryPart(PressurePart);

    const SizeType number_of_dofs = r_velocity_x.size() + r_velocity_y.size() + r_pressure.size();

    // Local system over the dofs: u_x, then u_y, then p
    if (rLeftHandSideMatrix.size1() != number_of_dofs ||
        rLeftHandSideMatrix.size2() != number_of_dofs) {
        rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs);
    if (rRightHandSideVector.size() != number_of_dofs) {
        rRightHandSideVector.resize(number_of_dofs, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_dofs);

    // Gauss point weight w, the same in every space
    const double weight = r_velocity_x.IntegrationPoints()[0].Weight()
        * std::abs(r_velocity_x.DeterminantOfJacobian(0));

    // Calculate the B matrix
    Matrix B;
    CalculateB(B);

    // Constitutive law
    ConstitutiveLaw::Parameters values(r_velocity_x, GetProperties(), rCurrentProcessInfo);
    Vector stress_vector = ZeroVector(3);
    Matrix D = ZeroMatrix(3, 3);
    ApplyConstitutiveLaw(B, values, stress_vector, D);

    // Add velocity terms in momentum equation
    AddMomentumTerms(
        rLeftHandSideMatrix, rRightHandSideVector, GetValue(BODY_FORCE), B, D, stress_vector,
        weight
    );

    // Add velocity-pressure terms
    AddContinuityTerms(rLeftHandSideMatrix, rRightHandSideVector, B, weight);

    // Optional PSPG
    const double pspg_factor = GetOptionalFactor(GetProperties(), PSPG_FACTOR);
    const double pspg_tau = pspg_factor > 0.0 ? CalculatePspgTau(pspg_factor) : 0.0;
    if (pspg_tau > 0.0) {
        AddPspgTerms(
            rLeftHandSideMatrix, rRightHandSideVector, GetValue(BODY_FORCE), D, pspg_tau, weight
        );
    }

    // Optional time term
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        AddTimeTerms(
            rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo[BDF_COEFFICIENTS],
            pspg_tau, weight
        );
    }

    // Weighted by the volume fraction \alpha_h, constant on the knot span
    const double alpha = GetValue(VOLUME_FRACTION);
    rLeftHandSideMatrix *= alpha;
    rRightHandSideVector *= alpha;

    // Optional grad-div, not weighted: on the whole background
    const double grad_div_factor = GetOptionalFactor(GetProperties(), GRAD_DIV_FACTOR);
    if (grad_div_factor > 0.0) {
        AddGradDivTerms(rLeftHandSideMatrix, rRightHandSideVector, B, grad_div_factor, weight);
    }

    KRATOS_CATCH("")
}

void GsbmStokesElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector;
    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector, rCurrentProcessInfo);
}

void GsbmStokesElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix;
    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo);
}

void GsbmStokesElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rResult.clear();
    for (IndexType part = 0; part < GetGeometry().NumberOfGeometryParts(); ++part) {
        const auto& r_variable = GetFieldVariable(part);
        for (const auto& r_node : GetGeometry().GetGeometryPart(part)) {
            rResult.push_back(r_node.GetDof(r_variable).EquationId());
        }
    }
}

void GsbmStokesElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    rElementalDofList.clear();
    for (IndexType part = 0; part < GetGeometry().NumberOfGeometryParts(); ++part) {
        const auto& r_variable = GetFieldVariable(part);
        for (const auto& r_node : GetGeometry().GetGeometryPart(part)) {
            rElementalDofList.push_back(r_node.pGetDof(r_variable));
        }
    }
}

int GsbmStokesElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const SizeType number_of_parts = GetGeometry().NumberOfGeometryParts();
    KRATOS_ERROR_IF(number_of_parts != 3) << Info()
        << " needs a geometry with the spaces of u_x, u_y and p (2D only)." << std::endl;
    KRATOS_ERROR_IF_NOT(Has(VOLUME_FRACTION)) << Info()
        << ": VOLUME_FRACTION is not set. Run VolumeFractionGsbmProcess first." << std::endl;
    KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW)
        || GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << Info() << ": no CONSTITUTIVE_LAW in the properties." << std::endl;

    // What the time term needs
    if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS)) {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(DENSITY)) << Info()
            << ": the time term needs DENSITY in the properties." << std::endl;
        KRATOS_ERROR_IF(GetGeometry().GetGeometryPart(VelocityXPart)[0].GetBufferSize()
            < rCurrentProcessInfo[BDF_COEFFICIENTS].size()) << Info()
            << ": the buffer size is smaller than the number of BDF coefficients." << std::endl;
    }

    // What PSPG's \tau needs
    if (GetOptionalFactor(GetProperties(), PSPG_FACTOR) > 0.0) {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(DYNAMIC_VISCOSITY)) << Info()
            << ": PSPG_FACTOR needs DYNAMIC_VISCOSITY in the properties." << std::endl;
        KRATOS_ERROR_IF_NOT(Has(KNOT_SPAN_SIZES)) << Info()
            << ": PSPG_FACTOR needs KNOT_SPAN_SIZES, set by IgaModelerGsbm." << std::endl;
    }

    // Each part's nodes carry its unknown
    for (IndexType part = 0; part < number_of_parts; ++part) {
        const auto& r_variable = GetFieldVariable(part);
        for (const auto& r_node : GetGeometry().GetGeometryPart(part)) {
            if (part == PressurePart) {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);
            } else {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            }
            KRATOS_CHECK_DOF_IN_NODE(r_variable, r_node);
        }
    }

    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    return GetProperties()[CONSTITUTIVE_LAW]->Check(
        GetProperties(), r_velocity_x, rCurrentProcessInfo
    );

    KRATOS_CATCH("")
}

void GsbmStokesElement::InitializeMaterial()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW)
        || GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << Info() << ": no CONSTITUTIVE_LAW in the properties." << std::endl;

    // The constitutive law of the properties, at the Gauss point of u_x's space
    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    mpConstitutiveLaw->InitializeMaterial(
        GetProperties(), r_velocity_x, row(r_velocity_x.ShapeFunctionsValues(), 0)
    );

    KRATOS_CATCH("")
}

void GsbmStokesElement::GetSolutionCoefficientVector(
    Vector& rValues) const
{
    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    const GeometryType& r_velocity_y = GetGeometry().GetGeometryPart(VelocityYPart);
    const SizeType n_x = r_velocity_x.size();
    if (rValues.size() != n_x + r_velocity_y.size()) {
        rValues.resize(n_x + r_velocity_y.size(), false);
    }

    for (IndexType i = 0; i < n_x; ++i) {
        rValues[i] = r_velocity_x[i].FastGetSolutionStepValue(VELOCITY_X);
    }
    for (IndexType i = 0; i < r_velocity_y.size(); ++i) {
        rValues[n_x + i] = r_velocity_y[i].FastGetSolutionStepValue(VELOCITY_Y);
    }
}

void GsbmStokesElement::CalculateB(Matrix& rB) const
{
    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    const GeometryType& r_velocity_y = GetGeometry().GetGeometryPart(VelocityYPart);
    const SizeType n_x = r_velocity_x.size();
    const SizeType n_u = n_x + r_velocity_y.size();

    // uv = xy: the local gradients are physical
    const Matrix& r_DN_x = r_velocity_x.ShapeFunctionsLocalGradients()[0];
    const Matrix& r_DN_y = r_velocity_y.ShapeFunctionsLocalGradients()[0];

    if (rB.size1() != 3 || rB.size2() != n_u) {
        rB.resize(3, n_u, false);
    }
    noalias(rB) = ZeroMatrix(3, n_u);

    // Voigt strain [xx, yy, xy]: u_x's dofs, then u_y's
    for (IndexType i = 0; i < n_x; ++i) {
        rB(0, i) = r_DN_x(i, 0);
        rB(2, i) = r_DN_x(i, 1);
    }
    for (IndexType i = 0; i < r_velocity_y.size(); ++i) {
        rB(1, n_x + i) = r_DN_y(i, 1);
        rB(2, n_x + i) = r_DN_y(i, 0);
    }
}

void GsbmStokesElement::ApplyConstitutiveLaw(
    const Matrix& rB,
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStressVector,
    Matrix& rD) const
{
    // Set constitutive law flags
    Flags& r_options = rValues.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Strain B u of the current velocity
    Vector velocity;
    GetSolutionCoefficientVector(velocity);
    Vector strain = prod(rB, velocity);
    rValues.SetStrainVector(strain);
    rValues.SetStressVector(rStressVector);
    rValues.SetConstitutiveMatrix(rD);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void GsbmStokesElement::AddMomentumTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const array_1d<double, 3>& rBodyForce,
    const Matrix& rB,
    const Matrix& rD,
    const Vector& rStressVector,
    const double Weight) const
{
    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    const GeometryType& r_velocity_y = GetGeometry().GetGeometryPart(VelocityYPart);
    const SizeType n_x = r_velocity_x.size();
    const SizeType n_u = rB.size2();

    // Diffusion term: (D B u, B v)
    noalias(subrange(rLeftHandSideMatrix, 0, n_u, 0, n_u)) +=
        Weight * prod(trans(rB), Matrix(prod(rD, rB)));

    // RHS corresponding term: the internal forces of the current stress
    noalias(subrange(rRightHandSideVector, 0, n_u)) -=
        Weight * prod(trans(rB), rStressVector);

    // Body force: (f, v)
    noalias(subrange(rRightHandSideVector, 0, n_x)) +=
        Weight * rBodyForce[0] * row(r_velocity_x.ShapeFunctionsValues(), 0);
    noalias(subrange(rRightHandSideVector, n_x, n_u)) +=
        Weight * rBodyForce[1] * row(r_velocity_y.ShapeFunctionsValues(), 0);
}

void GsbmStokesElement::AddContinuityTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const Matrix& rB,
    const double Weight) const
{
    const GeometryType& r_pressure = GetGeometry().GetGeometryPart(PressurePart);
    const SizeType n_u = rB.size2();
    const SizeType number_of_dofs = n_u + r_pressure.size();
    const Vector N_p = row(r_pressure.ShapeFunctionsValues(), 0);

    // \nabla . N of the velocity dofs: \partial_x N^x, then \partial_y N^y
    const Vector divergence = row(rB, 0) + row(rB, 1);

    // Current \nabla . u and p
    Vector velocity;
    GetSolutionCoefficientVector(velocity);
    const double div_u_current = inner_prod(divergence, velocity);
    double pressure_current = 0.0;
    for (IndexType i = 0; i < r_pressure.size(); ++i) {
        pressure_current += N_p[i] * r_pressure[i].FastGetSolutionStepValue(PRESSURE);
    }

    // Divergence term +(q, \nabla . u), and gradient term -(p, \nabla . v)
    noalias(subrange(rLeftHandSideMatrix, n_u, number_of_dofs, 0, n_u)) +=
        Weight * outer_prod(N_p, divergence);
    noalias(subrange(rLeftHandSideMatrix, 0, n_u, n_u, number_of_dofs)) -=
        Weight * outer_prod(divergence, N_p);

    // RHS corresponding terms
    noalias(subrange(rRightHandSideVector, n_u, number_of_dofs)) -=
        Weight * div_u_current * N_p;
    noalias(subrange(rRightHandSideVector, 0, n_u)) +=
        Weight * pressure_current * divergence;
}

double GsbmStokesElement::CalculatePspgTau(const double Factor) const
{
    const Vector& r_sizes = GetValue(KNOT_SPAN_SIZES);
    const double h = *std::min_element(r_sizes.begin(), r_sizes.end());
    return Factor * h * h / (4.0 * GetProperties()[DYNAMIC_VISCOSITY]);
}

void GsbmStokesElement::AddPspgTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const array_1d<double, 3>& rBodyForce,
    const Matrix& rD,
    const double Tau,
    const double Weight) const
{
    const GeometryType& r_velocity_x = GetGeometry().GetGeometryPart(VelocityXPart);
    const GeometryType& r_velocity_y = GetGeometry().GetGeometryPart(VelocityYPart);
    const GeometryType& r_pressure = GetGeometry().GetGeometryPart(PressurePart);
    const SizeType n_x = r_velocity_x.size();
    const SizeType n_u = n_x + r_velocity_y.size();

    // Second derivatives of the velocity spaces, columns [xx, xy, yy]
    const Matrix& r_DDN_x = r_velocity_x.ShapeFunctionDerivatives(
        2, 0, r_velocity_x.GetDefaultIntegrationMethod()
    );
    const Matrix& r_DDN_y = r_velocity_y.ShapeFunctionDerivatives(
        2, 0, r_velocity_y.GetDefaultIntegrationMethod()
    );

    // \nabla . \sigma of the velocity dofs, \sigma = D [e_xx, e_yy, 2 e_xy]
    Matrix divergence_stress(2, n_u);
    for (IndexType i = 0; i < n_x; ++i) {
        const double xx = r_DDN_x(i, 0);
        const double xy = r_DDN_x(i, 1);
        const double yy = r_DDN_x(i, 2);
        divergence_stress(0, i) = rD(0, 0) * xx + (rD(0, 2) + rD(2, 0)) * xy + rD(2, 2) * yy;
        divergence_stress(1, i) = rD(2, 0) * xx + (rD(2, 2) + rD(1, 0)) * xy + rD(1, 2) * yy;
    }
    for (IndexType i = 0; i < r_velocity_y.size(); ++i) {
        const double xx = r_DDN_y(i, 0);
        const double xy = r_DDN_y(i, 1);
        const double yy = r_DDN_y(i, 2);
        divergence_stress(0, n_x + i) = rD(0, 2) * xx + (rD(0, 1) + rD(2, 2)) * xy + rD(2, 1) * yy;
        divergence_stress(1, n_x + i) = rD(2, 2) * xx + (rD(2, 1) + rD(1, 2)) * xy + rD(1, 1) * yy;
    }

    // Current momentum residual \nabla p - \nabla . \sigma - f
    Vector velocity;
    GetSolutionCoefficientVector(velocity);
    const Matrix& r_DN_p = r_pressure.ShapeFunctionsLocalGradients()[0];
    array_1d<double, 2> residual;
    for (IndexType m = 0; m < 2; ++m) {
        residual[m] = -inner_prod(row(divergence_stress, m), velocity) - rBodyForce[m];
        for (IndexType j = 0; j < r_pressure.size(); ++j) {
            residual[m] += r_DN_p(j, m) * r_pressure[j].FastGetSolutionStepValue(PRESSURE);
        }
    }

    // (\nabla q, \tau (\nabla p - \nabla . \sigma - f)), in the pressure rows
    const double factor = Weight * Tau;
    for (IndexType i = 0; i < r_pressure.size(); ++i) {
        for (IndexType m = 0; m < 2; ++m) {
            const double row_factor = factor * r_DN_p(i, m);
            for (IndexType j = 0; j < n_u; ++j) {
                rLeftHandSideMatrix(n_u + i, j) -= row_factor * divergence_stress(m, j);
            }
            for (IndexType j = 0; j < r_pressure.size(); ++j) {
                rLeftHandSideMatrix(n_u + i, n_u + j) += row_factor * r_DN_p(j, m);
            }
            rRightHandSideVector[n_u + i] -= row_factor * residual[m];
        }
    }
}

void GsbmStokesElement::AddGradDivTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const Matrix& rB,
    const double Factor,
    const double Weight) const
{
    const SizeType n_u = rB.size2();

    // \nabla . N of the velocity dofs: \partial_x N^x, then \partial_y N^y
    const Vector divergence = row(rB, 0) + row(rB, 1);

    // Current \nabla . u
    Vector velocity;
    GetSolutionCoefficientVector(velocity);
    const double div_u_current = inner_prod(divergence, velocity);

    // (\beta_div \nabla . v, \nabla . u)
    noalias(subrange(rLeftHandSideMatrix, 0, n_u, 0, n_u)) +=
        Factor * Weight * outer_prod(divergence, divergence);
    noalias(subrange(rRightHandSideVector, 0, n_u)) -=
        Factor * Weight * div_u_current * divergence;
}

void GsbmStokesElement::AddTimeTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const Vector& rCoefficients,
    const double PspgTau,
    const double Weight) const
{
    const double density = GetProperties()[DENSITY];
    const GeometryType& r_pressure = GetGeometry().GetGeometryPart(PressurePart);
    const Matrix& r_DN_p = r_pressure.ShapeFunctionsLocalGradients()[0];
    const SizeType n_u = GetGeometry().GetGeometryPart(VelocityXPart).size()
        + GetGeometry().GetGeometryPart(VelocityYPart).size();
    const SizeType number_of_dofs = n_u + r_pressure.size();

    IndexType begin = 0;
    for (IndexType c = 0; c < 2; ++c) {
        const GeometryType& r_velocity = GetGeometry().GetGeometryPart(c);
        const auto& r_variable = GetFieldVariable(c);
        const SizeType n = r_velocity.size();
        const Vector N = row(r_velocity.ShapeFunctionsValues(), 0);

        // \partial_t u_c = \sum_k c_k u_c^{n+1-k} at the Gauss point
        double time_derivative = 0.0;
        for (IndexType k = 0; k < rCoefficients.size(); ++k) {
            for (IndexType i = 0; i < n; ++i) {
                time_derivative += rCoefficients[k] * N[i]
                    * r_velocity[i].FastGetSolutionStepValue(r_variable, k);
            }
        }

        // (v_c, \rho \partial_t u_c)
        noalias(subrange(rLeftHandSideMatrix, begin, begin + n, begin, begin + n)) +=
            Weight * density * rCoefficients[0] * outer_prod(N, N);
        noalias(subrange(rRightHandSideVector, begin, begin + n)) -=
            Weight * density * time_derivative * N;

        // (\partial_c q, \tau \rho \partial_t u_c), PSPG's residual
        if (PspgTau > 0.0) {
            const double pspg_factor = Weight * PspgTau * density;
            noalias(subrange(rLeftHandSideMatrix, n_u, number_of_dofs, begin, begin + n)) +=
                pspg_factor * rCoefficients[0] * outer_prod(column(r_DN_p, c), N);
            noalias(subrange(rRightHandSideVector, n_u, number_of_dofs)) -=
                pspg_factor * time_derivative * column(r_DN_p, c);
        }

        begin += n;
    }
}

} // namespace Kratos
