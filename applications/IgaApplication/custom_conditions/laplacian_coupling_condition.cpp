//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli

// System includes

// External includes

// Project includes
#include "custom_conditions/laplacian_coupling_condition.h"
#include "includes/variables.h"
#include "iga_application_variables.h"

namespace Kratos
{

void LaplacianCouplingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    const double integration_weight = r_integration_points[0].Weight();
    SetValue(INTEGRATION_WEIGHT, integration_weight);

    if (mIsGapSbmCoupling) {
        const auto& r_true = GetGeometry();
        const auto& r_surrogate_B = GetGeometryMirror();

        mDistanceVectorB.resize(3);
        noalias(mDistanceVectorB) = r_true.Center().Coordinates() - r_surrogate_B.Center().Coordinates();
    }
}

void LaplacianCouplingCondition::InitializeMemberVariables()
{
    const auto& r_geometry_patchA = GetGeometry();
    const auto& r_geometry_patchB = GetGeometryMirror();

    // Basis function order of the Taylor expansion is decided by the r_geometry_patchB
    const auto& r_DN_De = r_geometry_patchB.ShapeFunctionsLocalGradients(r_geometry_patchB.GetDefaultIntegrationMethod());
    mDim = r_DN_De[0].size2();
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    mBasisFunctionsOrder *= 2;

    // KRATOS_ERROR_IF(norm_2(r_geometry_patchA.Center()-r_geometry_patchB.Center()) > 1e-12)
    //     << "LaplacianCouplingCondition found non matching geometries." << std::endl;
    if (norm_2(r_geometry_patchA.Center()-r_geometry_patchB.Center()) > 1e-12)
        mIsGapSbmCoupling = true;

    mNormalParameterSpaceA = -r_geometry_patchA.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceA /= MathUtils<double>::Norm(mNormalParameterSpaceA);
    mNormalPhysicalSpaceA = mNormalParameterSpaceA;

    mNormalParameterSpaceB = -r_geometry_patchB.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceB /= MathUtils<double>::Norm(mNormalParameterSpaceB);
    mNormalPhysicalSpaceB = mNormalParameterSpaceB;

    KRATOS_ERROR_IF(std::abs(inner_prod(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB)+1) > 1e-12 && !mIsGapSbmCoupling)
        << "LaplacianCouplingCondition: normals are not opposite." << std::endl;
    
    SetValue(NORMAL, mNormalPhysicalSpaceA);
    SetValue(PROJECTION_NODE_COORDINATES, r_geometry_patchB.Center());
}

void LaplacianCouplingCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void LaplacianCouplingCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    const SizeType total_size = n_patch_A + n_patch_B;

    KRATOS_ERROR_IF(total_size == 0) << "LaplacianCouplingCondition found empty geometry." << std::endl;

    if (rLeftHandSideMatrix.size1() != total_size || rLeftHandSideMatrix.size2() != total_size) {
        rLeftHandSideMatrix.resize(total_size, total_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);

    const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    const auto& integration_points_patch_B = r_patch_B.IntegrationPoints();

    KRATOS_ERROR_IF(integration_points_patch_A.empty() || integration_points_patch_B.empty())
        << "LaplacianCouplingCondition expects at least one integration point." << std::endl;

    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_B =
        r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());

    KRATOS_ERROR_IF(rDN_De_A.size() == 0 || rDN_De_B.size() == 0)
        << "Shape function gradients not available for LaplacianCouplingCondition. Configure the modeler with derivative order >= 2." << std::endl;

    KRATOS_ERROR_IF(rDN_De_A[0].size1() != n_patch_A)
        << "LaplacianCouplingCondition: gradient matrix rows mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << rDN_De_A[0].size1() << "." << std::endl;
    KRATOS_ERROR_IF(rDN_De_B[0].size1() != n_patch_B)
        << "LaplacianCouplingCondition: gradient matrix rows mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << rDN_De_B[0].size1() << "." << std::endl;

    const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    Matrix N_patch_B = r_patch_B.ShapeFunctionsValues();

    if (mIsGapSbmCoupling) {
        Vector N_sum_B_vector(n_patch_B);
        ComputeTaylorExpansionContribution(r_patch_B, mDistanceVectorB, N_sum_B_vector);
        for (IndexType j = 0; j < n_patch_B; ++j) {
            N_patch_B(0, j) = N_sum_B_vector[j];
        }
    }

    KRATOS_ERROR_IF(N_patch_A.size2() != n_patch_A)
        << "LaplacianCouplingCondition: shape function values columns mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << N_patch_A.size2() << "." << std::endl;
    KRATOS_ERROR_IF(N_patch_B.size2() != n_patch_B)
        << "LaplacianCouplingCondition: shape function values columns mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << N_patch_B.size2() << "." << std::endl;

    const unsigned int dim_patch_A = rDN_De_A[0].size2();
    const unsigned int dim_patch_B = rDN_De_B[0].size2();

    KRATOS_ERROR_IF(dim_patch_A < 2 || dim_patch_B < 2)
        << "LaplacianCouplingCondition requires at least 2D parameter space." << std::endl;

    // Jacobians to compute the interface measure
    GeometryType::JacobiansType J_patch_A;
    r_patch_A.Jacobian(J_patch_A, r_patch_A.GetDefaultIntegrationMethod());
    Matrix jacobian_patch_A = ZeroMatrix(3, 3);
    jacobian_patch_A(0, 0) = J_patch_A[0](0, 0);
    jacobian_patch_A(0, 1) = J_patch_A[0](0, 1);
    jacobian_patch_A(1, 0) = J_patch_A[0](1, 0);
    jacobian_patch_A(1, 1) = J_patch_A[0](1, 1);
    jacobian_patch_A(2, 2) = 1.0;
    array_1d<double, 3> tangent_patch_A;
    r_patch_A.Calculate(LOCAL_TANGENT, tangent_patch_A);
    Vector determinant_factor_patch_A = prod(jacobian_patch_A, tangent_patch_A);
    determinant_factor_patch_A[2] = 0.0;
    const double detJ_patch_A = norm_2(determinant_factor_patch_A);
    KRATOS_ERROR_IF(detJ_patch_A <= std::numeric_limits<double>::epsilon())
        << "LaplacianCouplingCondition: degenerate Jacobian for patch A (determinant ~ 0)." << std::endl;
    const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;
    KRATOS_ERROR_IF(std::abs(weight) <= std::numeric_limits<double>::epsilon())
        << "LaplacianCouplingCondition: zero integration weight on patch A." << std::endl;

    double penalty_factor;
    double nitsche_penalty_free = 1.0;
    if (GetProperties().Has(PENALTY_FACTOR)) {
        penalty_factor = GetProperties()[PENALTY_FACTOR];
        penalty_factor = 10000000; 
        if (penalty_factor <= 0.0) {
            penalty_factor = 0.0;
            nitsche_penalty_free = -1.0;
        }
    } else {
        KRATOS_ERROR << "LaplacianCouplingCondition requires PENALTY_FACTOR to be defined in the Properties." << std::endl;
    }

    double characteristic_length;
    if (Has(KNOT_SPAN_SIZES)) {
        const Vector& knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
        if (!knot_span_sizes.empty()) {
            characteristic_length = knot_span_sizes[0];
            if (knot_span_sizes.size() > 1) {
                characteristic_length = std::min(characteristic_length, knot_span_sizes[1]);
            }
        } else {
            KRATOS_ERROR << "KNOT_SPAN_SIZES is empty in LaplacianCouplingCondition" << std::endl;
        }
    } else {
        KRATOS_ERROR << "KNOT_SPAN_SIZES is not defined in LaplacianCouplingCondition" << std::endl;
    }
    const double penalty_over_h = penalty_factor / characteristic_length;

    Matrix DN_patch_A = rDN_De_A[0];
    Matrix DN_patch_B = rDN_De_B[0];

    if (mIsGapSbmCoupling) {
        Matrix grad_H_B_T(3, n_patch_B);
        ComputeGradientTaylorExpansionContribution(r_patch_B, mDistanceVectorB, grad_H_B_T);
        DN_patch_B = trans(grad_H_B_T); // [n_patch_B x 3]
    }

    Vector DNn_A_n(n_patch_A, 0.0);
    Vector DNn_B_n(n_patch_B, 0.0);
    for (IndexType a = 0; a < n_patch_A; ++a) {
        for (unsigned int d = 0; d < dim_patch_A; ++d) {
            DNn_A_n[a] += DN_patch_A(a, d) * mNormalPhysicalSpaceA[d];
        }
    }
    for (IndexType b = 0; b < n_patch_B; ++b) {
        for (unsigned int d = 0; d < dim_patch_B; ++d) {
            DNn_B_n[b] += DN_patch_B(b, d) * mNormalPhysicalSpaceA[d];
        }
    }

    // Integration by parts + coupling condition 
    // Term: -∫ [[w]] · {{∇u}} = -∫ (v_A - v_B) · 0.5(∇u_A dot n_A + ∇u_B dot n_B)
    const IndexType offB = n_patch_A;
    // rows: test v_A with n_A
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const double vA = N_patch_A(0, i);
        // cols: u_A via (∇N_A · n_A)
        for (IndexType j = 0; j < n_patch_A; ++j)
            rLeftHandSideMatrix(i, j) -= weight * 0.5 * vA * DNn_A_n[j];
        // cols: u_B via (∇N_B · n_A)
        for (IndexType j = 0; j < n_patch_B; ++j)
            rLeftHandSideMatrix(i, offB + j) -= weight * 0.5 * vA * DNn_B_n[j];
    }
    // rows: test v_B with n_B
    for (IndexType i = 0; i < n_patch_B; ++i) {
        const double vB = N_patch_B(0, i);
        // cols: u_A via (∇N_A · n_B)
        for (IndexType j = 0; j < n_patch_A; ++j)
            rLeftHandSideMatrix(offB + i, j) += weight * 0.5 * vB * DNn_A_n[j];
        // cols: u_B via (∇N_B · n_B)
        for (IndexType j = 0; j < n_patch_B; ++j)
            rLeftHandSideMatrix(offB + i, offB + j) += weight * 0.5 * vB * DNn_B_n[j];
    }



    // Nitsche stabilization
    // 0.5*(grad v_A n_A + grad v_B n_B) · (u_A - u_B)
    for (IndexType i = 0; i < n_patch_A; ++i) {
        for (IndexType j = 0; j < n_patch_A; ++j) {
            rLeftHandSideMatrix(i, j) -= nitsche_penalty_free * 0.5 * weight * DNn_A_n[i] * N_patch_A(0, j);
        }
        for (IndexType j = 0; j < n_patch_B; ++j) {
            rLeftHandSideMatrix(i, n_patch_A + j) += nitsche_penalty_free * 0.5 * weight * DNn_A_n[i] * N_patch_B(0, j);
        }
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        for (IndexType j = 0; j < n_patch_A; ++j) {
            rLeftHandSideMatrix(n_patch_A + i, j) -= nitsche_penalty_free * 0.5 * weight * DNn_B_n[i] * N_patch_A(0, j);
        }
        for (IndexType j = 0; j < n_patch_B; ++j) {
            rLeftHandSideMatrix(n_patch_A + i, n_patch_A + j) += nitsche_penalty_free * 0.5 * weight * DNn_B_n[i] * N_patch_B(0, j);
        }
    }


    // Penalty term
    if (penalty_over_h != 0.0) {

        for (IndexType i = 0; i < n_patch_A; ++i) {
            for (IndexType j = 0; j < n_patch_A; ++j) {
                rLeftHandSideMatrix(i, j) += penalty_over_h * weight * N_patch_A(0, i) * N_patch_A(0, j);
            }
            for (IndexType j = 0; j < n_patch_B; ++j) {
                rLeftHandSideMatrix(i, n_patch_A + j) -= penalty_over_h * weight * N_patch_A(0, i) * N_patch_B(0, j);
            }
        }

        for (IndexType i = 0; i < n_patch_B; ++i) {
            for (IndexType j = 0; j < n_patch_A; ++j) {
                rLeftHandSideMatrix(n_patch_A + i, j) -= penalty_over_h * weight * N_patch_B(0, i) * N_patch_A(0, j);
            }
            for (IndexType j = 0; j < n_patch_B; ++j) {
                rLeftHandSideMatrix(n_patch_A + i, n_patch_A + j) += penalty_over_h * weight * N_patch_B(0, i) * N_patch_B(0, j);
            }
        }

    }


    KRATOS_CATCH("")
}

void LaplacianCouplingCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType lhs;
    CalculateLeftHandSide(lhs, rCurrentProcessInfo);

    const SizeType size = lhs.size1();
    KRATOS_ERROR_IF(size == 0) << "LaplacianCouplingCondition::CalculateRightHandSide: The left hand side matrix has zero size." << std::endl;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    Vector plus_values;
    GetSolutionCoefficientVectorA(r_unknown_var, plus_values);
    Vector minus_values;
    GetSolutionCoefficientVectorB(r_unknown_var, minus_values);

    Vector values(plus_values.size() + minus_values.size());
    for (IndexType i = 0; i < plus_values.size(); ++i) {
        values[i] = plus_values[i];
    }
    for (IndexType i = 0; i < minus_values.size(); ++i) {
        values[plus_values.size() + i] = minus_values[i];
    }

    if (rRightHandSideVector.size() != size) {
        rRightHandSideVector.resize(size, false);
    }
    noalias(rRightHandSideVector) = -prod(lhs, values);



    // analytical flux 

    // const auto& r_patch_A = GetGeometry();
    // const auto& r_patch_B = GetGeometryMirror();
    // const SizeType n_patch_A = r_patch_A.PointsNumber();
    // const SizeType n_patch_B = r_patch_B.PointsNumber();
    // const SizeType total_size = n_patch_A + n_patch_B;

    // if (rRightHandSideVector.size() != total_size) {
    //     rRightHandSideVector.resize(total_size, false);
    // }
    // noalias(rRightHandSideVector) = ZeroVector(total_size);

    // const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    // const auto& integration_points_patch_B = r_patch_B.IntegrationPoints();

    // KRATOS_ERROR_IF(integration_points_patch_A.empty() || integration_points_patch_B.empty())
    //     << "LaplacianCouplingCondition expects at least one integration point." << std::endl;

    // const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    // const Matrix& N_patch_B = r_patch_B.ShapeFunctionsValues();

    // GeometryType::JacobiansType J_patch_A;
    // r_patch_A.Jacobian(J_patch_A, r_patch_A.GetDefaultIntegrationMethod());
    // Matrix jacobian_patch_A = ZeroMatrix(3, 3);
    // jacobian_patch_A(0, 0) = J_patch_A[0](0, 0);
    // jacobian_patch_A(0, 1) = J_patch_A[0](0, 1);
    // jacobian_patch_A(1, 0) = J_patch_A[0](1, 0);
    // jacobian_patch_A(1, 1) = J_patch_A[0](1, 1);
    // jacobian_patch_A(2, 2) = 1.0;
    // array_1d<double, 3> tangent_patch_A;
    // r_patch_A.Calculate(LOCAL_TANGENT, tangent_patch_A);
    // Vector determinant_factor_patch_A = prod(jacobian_patch_A, tangent_patch_A);
    // determinant_factor_patch_A[2] = 0.0;
    // const double detJ_patch_A = norm_2(determinant_factor_patch_A);
    // const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;

    // array_1d<double, 3> global_coords_A;
    // r_patch_A.GlobalCoordinates(global_coords_A, integration_points_patch_A[0].Coordinates());
    // const double x_A = global_coords_A[0];
    // const double y_A = global_coords_A[1];
    // const double gradx_A = std::cos(x_A) * std::sinh(y_A);
    // const double grady_A = std::sin(x_A) * std::cosh(y_A);
    // const double flux_A = gradx_A * mNormalPhysicalSpaceA[0] + grady_A * mNormalPhysicalSpaceA[1];

    // array_1d<double, 3> global_coords_B;
    // r_patch_B.GlobalCoordinates(global_coords_B, integration_points_patch_B[0].Coordinates());
    // const double x_B = global_coords_B[0];
    // const double y_B = global_coords_B[1];
    // const double gradx_B = std::cos(x_B) * std::sinh(y_B);
    // const double grady_B = std::sin(x_B) * std::cosh(y_B);
    // const double flux_B = gradx_B * mNormalPhysicalSpaceB[0] + grady_B * mNormalPhysicalSpaceB[1];

    // for (IndexType i = 0; i < n_patch_A; ++i) {
    //     rRightHandSideVector[i] += weight * flux_A * N_patch_A(0, i);
    // }
    // for (IndexType i = 0; i < n_patch_B; ++i) {
    //     rRightHandSideVector[n_patch_A + i] += weight * flux_B * N_patch_B(0, i);
    // }
}

void LaplacianCouplingCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometryA = GetGeometry();
    const auto& r_geometryB = GetGeometryMirror();
    const SizeType n_patch_A = r_geometryA.PointsNumber();
    const SizeType n_patch_B = r_geometryB.PointsNumber();
    const SizeType total_size = n_patch_A + n_patch_B;

    if (rResult.size() != total_size) {
        rResult.resize(total_size, false);
    }

    const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = r_settings.GetUnknownVariable();

    SizeType counter = 0;
    for (IndexType i = 0; i < n_patch_A; ++i) {
        rResult[counter++] = r_geometryA[i].GetDof(r_unknown_var).EquationId();
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        rResult[counter++] = r_geometryB[i].GetDof(r_unknown_var).EquationId();
    }
}

void LaplacianCouplingCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();

    const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = r_settings.GetUnknownVariable();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(n_patch_A + n_patch_B);

    for (IndexType i = 0; i < n_patch_A; ++i) {
        rElementalDofList.push_back(r_patch_A[i].pGetDof(r_unknown_var));
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        rElementalDofList.push_back(r_patch_B[i].pGetDof(r_unknown_var));
    }
}

int LaplacianCouplingCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS))
        << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;

    const auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable())
        << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable())
        << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    return Condition::Check(rCurrentProcessInfo);
}

void LaplacianCouplingCondition::GetSolutionCoefficientVectorA(
    const Variable<double>& rUnknown,
    Vector& rValues) const
{
    const auto& r_patch_A = GetGeometry();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    if (rValues.size() != n_patch_A) {
        rValues.resize(n_patch_A, false);
    }
    for (IndexType i = 0; i < n_patch_A; ++i) {
        rValues[i] = r_patch_A[i].FastGetSolutionStepValue(rUnknown);
    }
}

void LaplacianCouplingCondition::GetSolutionCoefficientVectorB(
    const Variable<double>& rUnknown,
    Vector& rValues) const
{
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    if (rValues.size() != n_patch_B) {
        rValues.resize(n_patch_B, false);
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        rValues[i] = r_patch_B[i].FastGetSolutionStepValue(rUnknown);
    }
}

void LaplacianCouplingCondition::ComputeTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Vector& H_sum_vec)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    if (H_sum_vec.size() != number_of_control_points) H_sum_vec.resize(number_of_control_points);
    // derivatives up to order mBasisFunctionsOrder
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term = 0.0;
        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    const double derivative = r_shape_function_derivatives(i, k);
                    H_taylor_term += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
                int countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        const double derivative = r_shape_function_derivatives(i, countDerivativeId);
                        H_taylor_term += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0, i);
    }
}

void LaplacianCouplingCondition::ComputeGradientTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Matrix& grad_H_sum)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(rGeometry.GetDefaultIntegrationMethod());
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points) grad_H_sum.resize(3, number_of_control_points);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term_X = 0.0, H_taylor_term_Y = 0.0, H_taylor_term_Z = 0.0;
        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n - 1];
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k);
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k + 1);
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n - 1];
                IndexType countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        const double derivative = shapeFunctionDerivatives(i, countDerivativeId);
                        if (k_x >= 1) H_taylor_term_X += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x - 1, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        if (k_y >= 1) H_taylor_term_Y += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y - 1, rDistanceVector[2], k_z);
                        if (k_z >= 1) H_taylor_term_Z += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z - 1);
                        countDerivativeId++;
                    }
                }
            }
        }
        grad_H_sum(0, i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1, i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2, i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else
            grad_H_sum(2, i) = 0.0;
    }
}

double LaplacianCouplingCondition::ComputeTaylorTerm(const double derivative, const double dx, const IndexType n_k, const double dy, const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

double LaplacianCouplingCondition::ComputeTaylorTerm3D(const double derivative, const double dx, const IndexType k_x, const double dy, const IndexType k_y, const double dz, const IndexType k_z)
{
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
}

} // namespace Kratos
