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
}

void LaplacianCouplingCondition::InitializeMemberVariables()
{
    const auto& r_geometry_patchA = GetGeometry();
    const auto& r_geometry_patchB = GetGeometryMirror();

    KRATOS_ERROR_IF(norm_2(r_geometry_patchA.Center()-r_geometry_patchB.Center()) > 1e-12)
        << "LaplacianCouplingCondition found non matching geometries." << std::endl;

    mNormalParameterSpaceA = r_geometry_patchA.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceA /= MathUtils<double>::Norm(mNormalParameterSpaceA);
    mNormalPhysicalSpaceA = mNormalParameterSpaceA;

    mNormalParameterSpaceB = r_geometry_patchB.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceB /= MathUtils<double>::Norm(mNormalParameterSpaceB);
    mNormalPhysicalSpaceB = mNormalParameterSpaceB;

    KRATOS_ERROR_IF(std::abs(inner_prod(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB)+1) > 1e-12)
        << "LaplacianCouplingCondition: normals are not opposite." << std::endl;
}

void LaplacianCouplingCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    const SizeType total_size = rLeftHandSideMatrix.size1();
    KRATOS_ERROR_IF(total_size == 0) << "LaplacianCouplingCondition found empty LHS." << std::endl;

    if (rRightHandSideVector.size() != total_size) {
        rRightHandSideVector.resize(total_size, false);
    }

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    Vector plus_values;
    GetSolutionCoefficientVectorA(r_unknown_var, plus_values);
    Vector minus_values;
    GetSolutionCoefficientVectorB(r_unknown_var, minus_values);

    Vector values(total_size);
    for (IndexType i = 0; i < plus_values.size(); ++i) {
        values[i] = plus_values[i];
    }
    for (IndexType i = 0; i < minus_values.size(); ++i) {
        values[plus_values.size() + i] = minus_values[i];
    }

    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, values);
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

    const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    const Matrix& N_patch_B = r_patch_B.ShapeFunctionsValues();

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
    const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;

    double penalty_factor;
    double nitsche_penalty_free = -1.0;
    if (GetProperties().Has(PENALTY_FACTOR)) {
        penalty_factor = -1.0; // GetProperties()[PENALTY_FACTOR];
        // penalty_factor = 10.0; 
        if (penalty_factor <= 0.0) {
            penalty_factor = 0.0;
            nitsche_penalty_free = 1.0;
        }
    } else {
        KRATOS_ERROR << "LaplacianCouplingCondition requires PENALTY_FACTOR to be defined in the Properties." << std::endl;
    }

    double characteristic_length = 1.0;
    if (Has(KNOT_SPAN_SIZES)) {
        const Vector& knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
        if (!knot_span_sizes.empty()) {
            characteristic_length = knot_span_sizes[0];
            if (knot_span_sizes.size() > 1) {
                characteristic_length = std::min(characteristic_length, knot_span_sizes[1]);
            }
        }
    }
    if (characteristic_length <= 0.0) {
        characteristic_length = 1.0;
    }
    const double penalty_over_h = (characteristic_length > 0.0) ? penalty_factor / characteristic_length : 0.0;

    Vector DN_dot_n_patch_A = ZeroVector(n_patch_A);
    Vector DN_dot_n_patch_B = ZeroVector(n_patch_B);
    Matrix DN_patch_A = rDN_De_A[0];
    Matrix DN_patch_B = rDN_De_B[0];

    for (IndexType a = 0; a < n_patch_A; ++a) {
        for (unsigned int d = 0; d < dim_patch_A; ++d) {
            DN_dot_n_patch_A[a] += DN_patch_A(a, d) * mNormalPhysicalSpaceA[d];
        }
    }
    for (IndexType b = 0; b < n_patch_B; ++b) {
        for (unsigned int d = 0; d < dim_patch_B; ++d) {
            DN_dot_n_patch_B[b] += DN_patch_B(b, d) * mNormalPhysicalSpaceB[d];
        }
    }
    // Assembling coupling terms 
    for (IndexType i = 0; i < n_patch_A; ++i) {
        for (IndexType j = 0; j < n_patch_B; ++j) {
            rLeftHandSideMatrix(i, n_patch_A + j) += weight * N_patch_A(0, i) * DN_dot_n_patch_B[j];
        }
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        for (IndexType j = 0; j < n_patch_A; ++j) {
            rLeftHandSideMatrix(n_patch_A + i, j) += weight * N_patch_B(0, i) * DN_dot_n_patch_A[j];
        }
    }


    // Nitsche stabilization
    Vector DNn_A_nA(n_patch_A, 0.0);
    Vector DNn_A_nB(n_patch_A, 0.0);
    Vector DNn_B_nA(n_patch_B, 0.0);
    Vector DNn_B_nB(n_patch_B, 0.0);
    for (IndexType a = 0; a < n_patch_A; ++a) {
        for (unsigned int d = 0; d < dim_patch_A; ++d) {
            const double grad_val = rDN_De_A[0](a, d);
            DNn_A_nA[a] += grad_val * mNormalPhysicalSpaceA[d];
            DNn_A_nB[a] += grad_val * mNormalPhysicalSpaceB[d];
        }
    }
    for (IndexType b = 0; b < n_patch_B; ++b) {
        for (unsigned int d = 0; d < dim_patch_B; ++d) {
            const double grad_val = rDN_De_B[0](b, d);
            DNn_B_nA[b] += grad_val * mNormalPhysicalSpaceA[d];
            DNn_B_nB[b] += grad_val * mNormalPhysicalSpaceB[d];
        }
    }

    // 0.5*(grad v_A + grad v_B) · (u_A n_A + u_B n_B)
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const double grad_test_A_nA = DNn_A_nA[i];
        const double grad_test_A_nB = DNn_A_nB[i];
        for (IndexType j = 0; j < n_patch_A; ++j) {
            rLeftHandSideMatrix(i, j) += nitsche_penalty_free * 0.5 * weight * grad_test_A_nA * N_patch_A(0, j);
        }
        for (IndexType j = 0; j < n_patch_B; ++j) {
            rLeftHandSideMatrix(i, n_patch_A + j) += nitsche_penalty_free * 0.5 * weight * grad_test_A_nB * N_patch_B(0, j);
        }
    }

    for (IndexType i = 0; i < n_patch_B; ++i) {
        const double grad_test_B_nA = DNn_B_nA[i];
        const double grad_test_B_nB = DNn_B_nB[i];
        for (IndexType j = 0; j < n_patch_A; ++j) {
            rLeftHandSideMatrix(n_patch_A + i, j) += nitsche_penalty_free * 0.5 * weight * grad_test_B_nA * N_patch_A(0, j);
        }
        for (IndexType j = 0; j < n_patch_B; ++j) {
            rLeftHandSideMatrix(n_patch_A + i, n_patch_A + j) += nitsche_penalty_free * 0.5 * weight * grad_test_B_nB * N_patch_B(0, j);
        }
    }

    if (penalty_over_h != 0.0) {
        const double dot_nA_nA = MathUtils<double>::Dot(mNormalPhysicalSpaceA, mNormalPhysicalSpaceA);
        const double dot_nA_nB = MathUtils<double>::Dot(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB);
        const double dot_nB_nB = MathUtils<double>::Dot(mNormalPhysicalSpaceB, mNormalPhysicalSpaceB);

        for (IndexType i = 0; i < n_patch_A; ++i) {
            for (IndexType j = 0; j < n_patch_A; ++j) {
                rLeftHandSideMatrix(i, j) += penalty_over_h * weight * N_patch_A(0, i) * N_patch_A(0, j) * dot_nA_nA;
            }
            for (IndexType j = 0; j < n_patch_B; ++j) {
                rLeftHandSideMatrix(i, n_patch_A + j) += penalty_over_h * weight * N_patch_A(0, i) * N_patch_B(0, j) * dot_nA_nB;
            }
        }

        for (IndexType i = 0; i < n_patch_B; ++i) {
            for (IndexType j = 0; j < n_patch_A; ++j) {
                rLeftHandSideMatrix(n_patch_A + i, j) += penalty_over_h * weight * N_patch_B(0, i) * N_patch_A(0, j) * dot_nA_nB;
            }
            for (IndexType j = 0; j < n_patch_B; ++j) {
                rLeftHandSideMatrix(n_patch_A + i, n_patch_A + j) += penalty_over_h * weight * N_patch_B(0, i) * N_patch_B(0, j) * dot_nB_nB;
            }
        }

    }

    
    // // Rubén's penalty 
    // if (penalty_over_h != 0.0) {
    //     const double penalty_weight = penalty_over_h * weight;

    //     // Penalty contribution: (w_A + w_B) * (u_A - u_B)
    //     for (IndexType i = 0; i < n_patch_A; ++i) {
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             rLeftHandSideMatrix(i, j) += penalty_weight * N_patch_A(0, i) * N_patch_A(0, j);
    //         }
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             rLeftHandSideMatrix(i, n_patch_A + j) -= penalty_weight * N_patch_A(0, i) * N_patch_B(0, j);
    //         }
    //     }

    //     for (IndexType i = 0; i < n_patch_B; ++i) {
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             rLeftHandSideMatrix(n_patch_A + i, j) += penalty_weight * N_patch_B(0, i) * N_patch_A(0, j);
    //         }
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             rLeftHandSideMatrix(n_patch_A + i, n_patch_A + j) -= penalty_weight * N_patch_B(0, i) * N_patch_B(0, j);
    //         }
    //     }
    // }

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

} // namespace Kratos
