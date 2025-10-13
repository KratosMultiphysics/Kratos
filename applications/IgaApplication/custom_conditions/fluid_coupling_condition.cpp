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
#include "custom_conditions/fluid_coupling_condition.h"
#include "includes/variables.h"
#include "iga_application_variables.h"

namespace Kratos
{

void FluidCouplingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    const double integration_weight = r_integration_points[0].Weight();
    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void FluidCouplingCondition::InitializeMemberVariables()
{
    const auto& r_geometry_patchA = GetGeometry();
    const auto& r_geometry_patchB = GetGeometryMirror();

    KRATOS_ERROR_IF(norm_2(r_geometry_patchA.Center()-r_geometry_patchB.Center()) > 1e-12)
        << "FluidCouplingCondition found non matching geometries." << std::endl;

    mNormalParameterSpaceA = r_geometry_patchA.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceA /= MathUtils<double>::Norm(mNormalParameterSpaceA);
    mNormalPhysicalSpaceA = mNormalParameterSpaceA;

    mNormalParameterSpaceB = r_geometry_patchB.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceB /= MathUtils<double>::Norm(mNormalParameterSpaceB);
    mNormalPhysicalSpaceB = mNormalParameterSpaceB;

    KRATOS_ERROR_IF(std::abs(inner_prod(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB)+1) > 1e-12)
        << "FluidCouplingCondition: normals are not opposite." << std::endl;
}

void FluidCouplingCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void FluidCouplingCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    // Each node: mDim velocity components + 1 pressure
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A_all =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const unsigned int dim = rDN_De_A_all[0].size2();
    const SizeType dofs_per_node = dim + 1;
    const SizeType total_size = dofs_per_node * (n_patch_A + n_patch_B);

    KRATOS_ERROR_IF(total_size == 0) << "FluidCouplingCondition found empty geometry." << std::endl;

    if (rLeftHandSideMatrix.size1() != total_size || rLeftHandSideMatrix.size2() != total_size) {
        rLeftHandSideMatrix.resize(total_size, total_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);

    const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    const auto& integration_points_patch_B = r_patch_B.IntegrationPoints();

    KRATOS_ERROR_IF(integration_points_patch_A.empty() || integration_points_patch_B.empty())
        << "FluidCouplingCondition expects at least one integration point." << std::endl;

    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_B =
        r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());

    KRATOS_ERROR_IF(rDN_De_A.size() == 0 || rDN_De_B.size() == 0)
        << "Shape function gradients not available for FluidCouplingCondition. Configure the modeler with derivative order >= 2." << std::endl;

    KRATOS_ERROR_IF(rDN_De_A[0].size1() != n_patch_A)
        << "FluidCouplingCondition: gradient matrix rows mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << rDN_De_A[0].size1() << "." << std::endl;
    KRATOS_ERROR_IF(rDN_De_B[0].size1() != n_patch_B)
        << "FluidCouplingCondition: gradient matrix rows mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << rDN_De_B[0].size1() << "." << std::endl;

    const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    const Matrix& N_patch_B = r_patch_B.ShapeFunctionsValues();

    KRATOS_ERROR_IF(N_patch_A.size2() != n_patch_A)
        << "FluidCouplingCondition: shape function values columns mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << N_patch_A.size2() << "." << std::endl;
    KRATOS_ERROR_IF(N_patch_B.size2() != n_patch_B)
        << "FluidCouplingCondition: shape function values columns mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << N_patch_B.size2() << "." << std::endl;

    const unsigned int dim_patch_A = rDN_De_A[0].size2();
    const unsigned int dim_patch_B = rDN_De_B[0].size2();

    KRATOS_ERROR_IF(dim_patch_A < 2 || dim_patch_B < 2)
        << "FluidCouplingCondition requires at least 2D parameter space." << std::endl;

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
        << "FluidCouplingCondition: degenerate Jacobian for patch A (determinant ~ 0)." << std::endl;
    const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;
    KRATOS_ERROR_IF(std::abs(weight) <= std::numeric_limits<double>::epsilon())
        << "FluidCouplingCondition: zero integration weight on patch A." << std::endl;

    double penalty_factor;
    double nitsche_penalty_free = -1.0; // should be = -1
    if (GetProperties().Has(PENALTY_FACTOR)) {
        // penalty_factor = -1.0; // GetProperties()[PENALTY_FACTOR];
        penalty_factor = 1000.0; 
        if (penalty_factor <= 0.0) {
            penalty_factor = 0.0;
            nitsche_penalty_free = -1.0; // does not matter this sign
        }
    } else {
        KRATOS_ERROR << "FluidCouplingCondition requires PENALTY_FACTOR to be defined in the Properties." << std::endl;
    }



    double characteristic_length = 1.0;
    if (Has(KNOT_SPAN_SIZES)) {
        const Vector& knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
        if (!knot_span_sizes.empty()) {
            characteristic_length = knot_span_sizes[0];
            if (knot_span_sizes.size() > 1) {
                characteristic_length = std::min(characteristic_length, knot_span_sizes[1]);
            }
        } else {
            KRATOS_ERROR << "KNOT_SPAN_SIZES is empty in FluidCouplingCondition" << std::endl;
        }
    }
    if (characteristic_length <= 0.0) {
        KRATOS_ERROR << "KNOT_SPAN_SIZES is < 0 in FluidCouplingCondition" << std::endl;
    }
    const double penalty_over_h = penalty_factor / characteristic_length;

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

    // Helper to compute index in block with (node, component, side)
    auto idxA = [dofs_per_node](IndexType node, IndexType comp){ return node * dofs_per_node + comp; };
    auto idxB = [dofs_per_node, n_patch_A](IndexType node, IndexType comp){ return n_patch_A * dofs_per_node + node * dofs_per_node + comp; };

    // Integration by parts style terms, applied per velocity component
    for (IndexType comp = 0; comp < dim; ++comp) {
        // rows: test from A with n_A
        for (IndexType i = 0; i < n_patch_A; ++i) {
            const double vA = N_patch_A(0, i);
            const IndexType row = idxA(i, comp);
            // cols: u_A via (∇N_A · n_A)
            for (IndexType j = 0; j < n_patch_A; ++j) {
                const IndexType col = idxA(j, comp);
                rLeftHandSideMatrix(row, col) -= weight * 0.5 * vA * DNn_A_nA[j];
            }
            // cols: u_B via (∇N_B · n_A)
            for (IndexType j = 0; j < n_patch_B; ++j) {
                const IndexType col = idxB(j, comp);
                rLeftHandSideMatrix(row, col) -= weight * 0.5 * vA * DNn_B_nA[j];
            }
        }
        // rows: test from B with n_B
        for (IndexType i = 0; i < n_patch_B; ++i) {
            const double vB = N_patch_B(0, i);
            const IndexType row = idxB(i, comp);
            // cols: u_A via (∇N_A · n_B)
            for (IndexType j = 0; j < n_patch_A; ++j) {
                const IndexType col = idxA(j, comp);
                rLeftHandSideMatrix(row, col) -= weight * 0.5 * vB * DNn_A_nB[j];
            }
            // cols: u_B via (∇N_B · n_B)
            for (IndexType j = 0; j < n_patch_B; ++j) {
                const IndexType col = idxB(j, comp);
                rLeftHandSideMatrix(row, col) -= weight * 0.5 * vB * DNn_B_nB[j];
            }
        }
    }



    // Nitsche stabilization
    // 0.5*(grad v_A + grad v_B) · (u_A n_A + u_B n_B)
    for (IndexType comp = 0; comp < dim; ++comp) {
        for (IndexType i = 0; i < n_patch_A; ++i) {
            const double grad_test_A_nA = DNn_A_nA[i];
            const double grad_test_A_nB = DNn_A_nB[i];
            const IndexType row = idxA(i, comp);
            for (IndexType j = 0; j < n_patch_A; ++j) {
                const IndexType col = idxA(j, comp);
                rLeftHandSideMatrix(row, col) += nitsche_penalty_free * 0.5 * weight * grad_test_A_nA * N_patch_A(0, j);
            }
            for (IndexType j = 0; j < n_patch_B; ++j) {
                const IndexType col = idxB(j, comp);
                rLeftHandSideMatrix(row, col) += nitsche_penalty_free * 0.5 * weight * grad_test_A_nB * N_patch_B(0, j);
            }
        }
        for (IndexType i = 0; i < n_patch_B; ++i) {
            const double grad_test_B_nA = DNn_B_nA[i];
            const double grad_test_B_nB = DNn_B_nB[i];
            const IndexType row = idxB(i, comp);
            for (IndexType j = 0; j < n_patch_A; ++j) {
                const IndexType col = idxA(j, comp);
                rLeftHandSideMatrix(row, col) += nitsche_penalty_free * 0.5 * weight * grad_test_B_nA * N_patch_A(0, j);
            }
            for (IndexType j = 0; j < n_patch_B; ++j) {
                const IndexType col = idxB(j, comp);
                rLeftHandSideMatrix(row, col) += nitsche_penalty_free * 0.5 * weight * grad_test_B_nB * N_patch_B(0, j);
            }
        }
    }


    // Penalty term
    if (penalty_over_h != 0.0) {
        const double dot_nA_nA = MathUtils<double>::Dot(mNormalPhysicalSpaceA, mNormalPhysicalSpaceA);
        const double dot_nA_nB = MathUtils<double>::Dot(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB);
        const double dot_nB_nB = MathUtils<double>::Dot(mNormalPhysicalSpaceB, mNormalPhysicalSpaceB);

        for (IndexType comp = 0; comp < dim; ++comp) {
            for (IndexType i = 0; i < n_patch_A; ++i) {
                const IndexType rowA = idxA(i, comp);
                for (IndexType j = 0; j < n_patch_A; ++j) {
                    const IndexType colA = idxA(j, comp);
                    rLeftHandSideMatrix(rowA, colA) += penalty_over_h * weight * N_patch_A(0, i) * N_patch_A(0, j) * dot_nA_nA;
                }
                for (IndexType j = 0; j < n_patch_B; ++j) {
                    const IndexType colB = idxB(j, comp);
                    rLeftHandSideMatrix(rowA, colB) += penalty_over_h * weight * N_patch_A(0, i) * N_patch_B(0, j) * dot_nA_nB;
                }
            }

            for (IndexType i = 0; i < n_patch_B; ++i) {
                const IndexType rowB = idxB(i, comp);
                for (IndexType j = 0; j < n_patch_A; ++j) {
                    const IndexType colA = idxA(j, comp);
                    rLeftHandSideMatrix(rowB, colA) += penalty_over_h * weight * N_patch_B(0, i) * N_patch_A(0, j) * dot_nA_nB;
                }
                for (IndexType j = 0; j < n_patch_B; ++j) {
                    const IndexType colB = idxB(j, comp);
                    rLeftHandSideMatrix(rowB, colB) += penalty_over_h * weight * N_patch_B(0, i) * N_patch_B(0, j) * dot_nB_nB;
                }
            }
        }

    }

    KRATOS_CATCH("")
}

void FluidCouplingCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType lhs;
    CalculateLeftHandSide(lhs, rCurrentProcessInfo);

    const SizeType size = lhs.size1();
    KRATOS_ERROR_IF(size == 0) << "FluidCouplingCondition::CalculateRightHandSide: The left hand side matrix has zero size." << std::endl;

    // Build a full unknown vector [u_A, (v_A, w_A,) p_A, u_B, (v_B, w_B,) p_B]
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const SizeType dofs_per_node = dim + 1;

    auto idxA = [dofs_per_node](IndexType node, IndexType comp){ return node * dofs_per_node + comp; };
    auto idxB = [dofs_per_node, n_patch_A](IndexType node, IndexType comp){ return n_patch_A * dofs_per_node + node * dofs_per_node + comp; };

    Vector valsA, valsB;
    GetVelocityCoefficientVectorA(valsA);
    GetVelocityCoefficientVectorB(valsB);

    Vector values_full(size, false);
    noalias(values_full) = ZeroVector(size);

    // Fill velocity entries; pressures remain zero
    for (IndexType i = 0; i < n_patch_A; ++i) {
        for (IndexType c = 0; c < dim; ++c) {
            values_full[idxA(i, c)] = valsA[i * dim + c];
        }
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        for (IndexType c = 0; c < dim; ++c) {
            values_full[idxB(i, c)] = valsB[i * dim + c];
        }
    }

    if (rRightHandSideVector.size() != size) {
        rRightHandSideVector.resize(size, false);
    }
    noalias(rRightHandSideVector) = -prod(lhs, values_full);



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

void FluidCouplingCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometryA = GetGeometry();
    const auto& r_geometryB = GetGeometryMirror();
    const SizeType n_patch_A = r_geometryA.PointsNumber();
    const SizeType n_patch_B = r_geometryB.PointsNumber();

    // Determine velocity dimension from shape gradients (consistent with LHS)
    const auto& rDN_De_A = r_geometryA.ShapeFunctionsLocalGradients(r_geometryA.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();

    const SizeType dofs_per_node = dim + 1; // dim velocities + pressure
    const SizeType total_size = dofs_per_node * (n_patch_A + n_patch_B);

    if (rResult.size() != total_size) {
        rResult.resize(total_size, false);
    }

    SizeType k = 0;
    // Validate DOFs on A side
    for (IndexType i = 0; i < n_patch_A; ++i) {
        KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(VELOCITY_X)) << "Node missing VELOCITY_X DOF in A, node id: " << r_geometryA[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(VELOCITY_Y)) << "Node missing VELOCITY_Y DOF in A, node id: " << r_geometryA[i].Id() << std::endl;
        if (dim == 3) KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(VELOCITY_Z)) << "Node missing VELOCITY_Z DOF in A, node id: " << r_geometryA[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(PRESSURE))   << "Node missing PRESSURE DOF in A, node id: "   << r_geometryA[i].Id() << std::endl;
    }
    // Validate DOFs on B side
    for (IndexType i = 0; i < n_patch_B; ++i) {
        KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(VELOCITY_X)) << "Node missing VELOCITY_X DOF in B, node id: " << r_geometryB[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(VELOCITY_Y)) << "Node missing VELOCITY_Y DOF in B, node id: " << r_geometryB[i].Id() << std::endl;
        if (dim == 3) KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(VELOCITY_Z)) << "Node missing VELOCITY_Z DOF in B, node id: " << r_geometryB[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(PRESSURE))   << "Node missing PRESSURE DOF in B, node id: "   << r_geometryB[i].Id() << std::endl;
    }

    // Fill EIDs and track min/max for debugging
    int min_eid = std::numeric_limits<int>::max();
    int max_eid = std::numeric_limits<int>::lowest();
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const int ex = r_geometryA[i].GetDof(VELOCITY_X).EquationId();   min_eid = std::min(min_eid, ex); max_eid = std::max(max_eid, ex); rResult[k++] = ex;
        const int ey = r_geometryA[i].GetDof(VELOCITY_Y).EquationId();   min_eid = std::min(min_eid, ey); max_eid = std::max(max_eid, ey); rResult[k++] = ey;
        if (dim == 3) {
            const int ez = r_geometryA[i].GetDof(VELOCITY_Z).EquationId(); min_eid = std::min(min_eid, ez); max_eid = std::max(max_eid, ez); rResult[k++] = ez;
        }
        const int ep = r_geometryA[i].GetDof(PRESSURE).EquationId();     min_eid = std::min(min_eid, ep); max_eid = std::max(max_eid, ep); rResult[k++] = ep;
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        const int ex = r_geometryB[i].GetDof(VELOCITY_X).EquationId();   min_eid = std::min(min_eid, ex); max_eid = std::max(max_eid, ex); rResult[k++] = ex;
        const int ey = r_geometryB[i].GetDof(VELOCITY_Y).EquationId();   min_eid = std::min(min_eid, ey); max_eid = std::max(max_eid, ey); rResult[k++] = ey;
        if (dim == 3) {
            const int ez = r_geometryB[i].GetDof(VELOCITY_Z).EquationId(); min_eid = std::min(min_eid, ez); max_eid = std::max(max_eid, ez); rResult[k++] = ez;
        }
        const int ep = r_geometryB[i].GetDof(PRESSURE).EquationId();     min_eid = std::min(min_eid, ep); max_eid = std::max(max_eid, ep); rResult[k++] = ep;
    }

}

void FluidCouplingCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();

    // Determine velocity dimension from shape gradients (consistent with LHS)
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();

    rElementalDofList.resize(0);
    rElementalDofList.reserve((dim + 1) * (n_patch_A + n_patch_B));

    for (IndexType i = 0; i < n_patch_A; ++i) {
        rElementalDofList.push_back(r_patch_A[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(r_patch_A[i].pGetDof(VELOCITY_Y));
        if (dim == 3) rElementalDofList.push_back(r_patch_A[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(r_patch_A[i].pGetDof(PRESSURE));
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        rElementalDofList.push_back(r_patch_B[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(r_patch_B[i].pGetDof(VELOCITY_Y));
        if (dim == 3) rElementalDofList.push_back(r_patch_B[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(r_patch_B[i].pGetDof(PRESSURE));
    }
}

int FluidCouplingCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    return Condition::Check(rCurrentProcessInfo);
}

void FluidCouplingCondition::GetVelocityCoefficientVectorA(
    Vector& rValues) const
{
    const auto& r_patch_A = GetGeometry();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const SizeType nvals = n_patch_A * dim;
    if (rValues.size() != nvals) rValues.resize(nvals, false);
    IndexType k = 0;
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const auto& v = r_patch_A[i].FastGetSolutionStepValue(VELOCITY);
        rValues[k++] = v[0];
        rValues[k++] = v[1];
        if (dim == 3) rValues[k++] = v[2];
    }
}

void FluidCouplingCondition::GetVelocityCoefficientVectorB(
    Vector& rValues) const
{
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    const auto& rDN_De_B = r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_B.empty() ? 2 : rDN_De_B[0].size2();
    const SizeType nvals = n_patch_B * dim;
    if (rValues.size() != nvals) rValues.resize(nvals, false);
    IndexType k = 0;
    for (IndexType i = 0; i < n_patch_B; ++i) {
        const auto& v = r_patch_B[i].FastGetSolutionStepValue(VELOCITY);
        rValues[k++] = v[0];
        rValues[k++] = v[1];
        if (dim == 3) rValues[k++] = v[2];
    }
}

} // namespace Kratos
