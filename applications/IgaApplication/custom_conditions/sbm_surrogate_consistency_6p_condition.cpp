//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

// System includes

// External includes

// Project includes
#include "custom_conditions/sbm_surrogate_consistency_6p_condition.h"

namespace Kratos
{

void SbmSurrogateConsistency6pCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(Has(SBM_SURROGATE_REFERENCE_DIRECTION))
        << "SbmSurrogateConsistency6pCondition " << Id() << ": SBM_SURROGATE_REFERENCE_DIRECTION not set -- "
        << "must be set at condition-creation time (from the SurrogateBoundarySegment's own "
        << "ActiveSpanIndexU/V) before Initialize() runs." << std::endl;

    const Vector& r_reference_direction = GetValue(SBM_SURROGATE_REFERENCE_DIRECTION);
    KRATOS_ERROR_IF(r_reference_direction.size() != 2)
        << "SbmSurrogateConsistency6pCondition " << Id()
        << ": SBM_SURROGATE_REFERENCE_DIRECTION must have size 2 (du, dv)." << std::endl;

    array_1d<double, 3> local_tangent;
    GetGeometry().Calculate(LOCAL_TANGENT, local_tangent);
    const double tangent_norm = std::sqrt(local_tangent[0] * local_tangent[0] + local_tangent[1] * local_tangent[1]);

    const double nu_param_0 = local_tangent[1] / tangent_norm;
    const double nu_param_1 = -local_tangent[0] / tangent_norm;

    const double dot = nu_param_0 * r_reference_direction[0] + nu_param_1 * r_reference_direction[1];
    SetValue(SBM_SURROGATE_CONORMAL_SIGN, (dot > 0.0) ? -1.0 : 1.0);

    KRATOS_CATCH("")
}

void SbmSurrogateConsistency6pCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = 6 * number_of_nodes;

    const double thickness = GetProperties()[THICKNESS];
    const double conormal_sign = GetValue(SBM_SURROGATE_CONORMAL_SIGN);

    Matrix D_local;
    CalculateLocalConstitutiveMatrix(GetProperties(), D_local);

    Matrix K = ZeroMatrix(mat_size, mat_size);

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

    // 2-point Gauss quadrature through the thickness (matches Shell6pElement).
    const double gauss_zeta[2] = { -std::sqrt(1.0 / 3.0), std::sqrt(1.0 / 3.0) };
    const double gauss_weight[2] = { 1.0, 1.0 };

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        KinematicVariables kinematics(3);
        CalculateKinematics(point_number, kinematics, r_geometry);

        Matrix normal_derivatives;
        CalculateNormalVectorDerivatives(point_number, kinematics, normal_derivatives, r_geometry);

        Matrix T;
        CalculateTransformationFromLocalToGlobalCartesian(kinematics, T);
        const Matrix D = prod(T, Matrix(prod(D_local, trans(T))));

        array_1d<double, 3> local_tangent;
        r_geometry.Calculate(LOCAL_TANGENT, local_tangent);
        const array_1d<double, 3> t_vec = local_tangent[0] * kinematics.BaseVector1 + local_tangent[1] * kinematics.BaseVector2;
        const double surface_jacobian = norm_2(t_vec);

        const double integration_weight = integration_points[point_number].Weight();

        for (IndexType gauss_index = 0; gauss_index < 2; ++gauss_index)
        {
            const double zeta = gauss_zeta[gauss_index];

            Matrix jacobian_inv;
            double jacobian_det;
            CalculateThicknessJacobian(zeta, thickness, kinematics, normal_derivatives, jacobian_inv, jacobian_det);

            array_1d<double, 3> unit_conormal;
            double area_scale;
            CalculateLateralConormal(point_number, jacobian_inv, jacobian_det, r_geometry, unit_conormal, area_scale);
            unit_conormal *= conormal_sign;

            Matrix B, N_op;
            CalculateBAndDisplacementOperator(point_number, zeta, thickness, jacobian_inv, normal_derivatives, kinematics, r_geometry, B, N_op);

            // Voigt traction matrix T_n (3x6), Voigt order [xx,yy,zz,xy,yz,xz].
            Matrix Tn = ZeroMatrix(3, 6);
            Tn(0, 0) = unit_conormal[0]; Tn(0, 3) = unit_conormal[1]; Tn(0, 5) = unit_conormal[2];
            Tn(1, 1) = unit_conormal[1]; Tn(1, 3) = unit_conormal[0]; Tn(1, 4) = unit_conormal[2];
            Tn(2, 2) = unit_conormal[2]; Tn(2, 4) = unit_conormal[1]; Tn(2, 5) = unit_conormal[0];

            const Matrix F = prod(Tn, Matrix(prod(D, B)));

            const double dW = integration_weight * surface_jacobian * area_scale * gauss_weight[gauss_index];

            noalias(K) += prod(trans(N_op), F) * (-1.0 * dW);
        }
    }

    if (CalculateStiffnessMatrixFlag) {
        if (rLeftHandSideMatrix.size1() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = K;
    }

    if (CalculateResidualVectorFlag) {
        if (rRightHandSideVector.size() != mat_size) {
            rRightHandSideVector.resize(mat_size, false);
        }

        Vector d_current = ZeroVector(mat_size);
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3>& disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& rot = r_geometry[i].FastGetSolutionStepValue(ROTATION);
            const IndexType index = 6 * i;
            d_current[index] = disp[0]; d_current[index + 1] = disp[1]; d_current[index + 2] = disp[2];
            d_current[index + 3] = rot[0]; d_current[index + 4] = rot[1]; d_current[index + 5] = rot[2];
        }

        noalias(rRightHandSideVector) = -prod(K, d_current);
    }

    KRATOS_CATCH("")
}

int SbmSurrogateConsistency6pCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(THICKNESS))
        << "No THICKNESS defined in property of SbmSurrogateConsistency6pCondition" << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(YOUNG_MODULUS))
        << "No YOUNG_MODULUS defined in property of SbmSurrogateConsistency6pCondition" << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(POISSON_RATIO))
        << "No POISSON_RATIO defined in property of SbmSurrogateConsistency6pCondition" << std::endl;
    KRATOS_ERROR_IF_NOT(Has(SBM_SURROGATE_REFERENCE_DIRECTION))
        << "SbmSurrogateConsistency6pCondition " << Id() << ": SBM_SURROGATE_REFERENCE_DIRECTION not set -- "
        << "must be set at condition-creation time; SBM_SURROGATE_CONORMAL_SIGN is then derived from it "
        << "automatically in Initialize()." << std::endl;
    return 0;
}

void SbmSurrogateConsistency6pCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
    const IndexType number_of_nodes = r_geometry.size();

    if (rResult.size() != 6 * number_of_nodes)
        rResult.resize(6 * number_of_nodes, false);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const IndexType index = i * 6;
        const auto& r_node = r_geometry[i];
        rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
    }

    KRATOS_CATCH("")
}

void SbmSurrogateConsistency6pCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
    const IndexType number_of_nodes = r_geometry.size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(6 * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto& r_node = r_geometry.GetPoint(i);
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
        rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
        rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
    }

    KRATOS_CATCH("")
}

} // Namespace Kratos
