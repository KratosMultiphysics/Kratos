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
#include <unordered_map>

// External includes

// Project includes
#include "custom_conditions/coupling_sbm_extension_operator_6p_condition.h"
#include "custom_utilities/iga_flags.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace {

void BuildExtendedDisplacementOperator(
    const Vector& rWeights,
    double zeta,
    double Thickness,
    const array_1d<double, 3>& rNormalVector,
    Matrix& rNOperatorExtended)
{
    using IndexType = std::size_t;
    const std::size_t n_dof_full = rWeights.size();
    const std::size_t mat_size = 6 * n_dof_full;

    if (rNOperatorExtended.size1() != 3 || rNOperatorExtended.size2() != mat_size)
        rNOperatorExtended.resize(3, mat_size);
    noalias(rNOperatorExtended) = ZeroMatrix(3, mat_size);

    const double normal_x = rNormalVector[0];
    const double normal_y = rNormalVector[1];
    const double normal_z = rNormalVector[2];
    const double half_t_zeta = zeta * (Thickness / 2.0);

    for (IndexType k = 0; k < n_dof_full; ++k) {
        const double w = rWeights[k];
        if (w == 0.0) continue;

        const IndexType index = 6 * k;
        rNOperatorExtended(0, index)     = w;
        rNOperatorExtended(1, index + 1) = w;
        rNOperatorExtended(2, index + 2) = w;

        const double coeff = w * half_t_zeta;
        rNOperatorExtended(0, index + 4) =  coeff * normal_z;
        rNOperatorExtended(0, index + 5) = -coeff * normal_y;
        rNOperatorExtended(1, index + 3) = -coeff * normal_z;
        rNOperatorExtended(1, index + 5) =  coeff * normal_x;
        rNOperatorExtended(2, index + 3) =  coeff * normal_y;
        rNOperatorExtended(2, index + 4) = -coeff * normal_x;
    }
}

} 

CouplingSbmExtensionOperator6pCondition::SizeType CouplingSbmExtensionOperator6pCondition::GetFullDofNodeCount() const
{
    return GetValue(SBM_MLS_ALL_DOF_NODES).size();
}

void CouplingSbmExtensionOperator6pCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const double stabilization_parameter = GetProperties()[NITSCHE_STABILIZATION_FACTOR];
    const bool couple_rotation_x = Is(IgaFlags::FIX_ROTATION_X);
    const bool couple_rotation_y = Is(IgaFlags::FIX_ROTATION_Y);
    const bool couple_rotation_z = Is(IgaFlags::FIX_ROTATION_Z);
    double penalty_rotation = 0.0;
    if (couple_rotation_x || couple_rotation_y || couple_rotation_z)
        penalty_rotation = GetProperties()[PENALTY_ROTATION_FACTOR];

    const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
    const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

    const SizeType number_of_nodes_master = r_geometry_master.size();
    const SizeType number_of_nodes_slave = r_geometry_slave.size();
    const SizeType mat_size_own_master = 6 * number_of_nodes_master;
    const SizeType mat_size_own = 6 * (number_of_nodes_master + number_of_nodes_slave);

    const Matrix& r_master_weights = GetValue(SBM_MLS_MASTER_WEIGHTS);
    const Matrix& r_slave_weights = GetValue(SBM_MLS_SLAVE_WEIGHTS);
    const SizeType n_dof_full = r_master_weights.size2();
    const SizeType mat_size = 6 * n_dof_full;

    KRATOS_ERROR_IF(r_slave_weights.size2() != n_dof_full)
        << "CouplingSbmExtensionOperator6pCondition " << Id()
        << ": mismatched SBM_MLS_MASTER_WEIGHTS/SBM_MLS_SLAVE_WEIGHTS column counts." << std::endl;
    KRATOS_ERROR_IF(mat_size < mat_size_own)
        << "CouplingSbmExtensionOperator6pCondition " << Id()
        << ": stored SBM_MLS_*_WEIGHTS DOF count is smaller than this condition's own edge DOF count "
        << "-- PrecomputeAndStoreCouplingExtensionOperators must be re-run after any geometry change." << std::endl;

    Matrix K = ZeroMatrix(mat_size, mat_size);

    std::unordered_map<IndexType, IndexType> dof_node_id_to_column;
    if (mpMasterShiftSource != nullptr && mpSlaveShiftSource != nullptr) {
        const auto& r_all_dof_nodes_map = GetValue(SBM_MLS_ALL_DOF_NODES);
        for (IndexType k = 0; k < r_all_dof_nodes_map.size(); ++k) {
            dof_node_id_to_column[r_all_dof_nodes_map[k]->Id()] = k;
        }
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

    const Properties& r_props_master = GetProperties().GetSubProperties().front();
    const Properties& r_props_slave = GetProperties().GetSubProperties().back();
    const double thickness_master = r_props_master[THICKNESS];
    const double thickness_slave = r_props_slave[THICKNESS];

    Matrix D_local_master, D_local_slave;
    CalculateLocalConstitutiveMatrix(r_props_master, D_local_master);
    CalculateLocalConstitutiveMatrix(r_props_slave, D_local_slave);

    // 2-point Gauss quadrature through the thickness.
    const double gauss_zeta[2] = { -std::sqrt(1.0 / 3.0), std::sqrt(1.0 / 3.0) };
    const double gauss_weight[2] = { 1.0, 1.0 };

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        KinematicVariables kinematics_master(3);
        KinematicVariables kinematics_slave(3);
        CalculateKinematics(point_number, kinematics_master, PatchType::Master);
        CalculateKinematics(point_number, kinematics_slave, PatchType::Slave);

        Matrix normal_derivatives_master, normal_derivatives_slave;
        CalculateNormalVectorDerivatives(point_number, kinematics_master, normal_derivatives_master, PatchType::Master);
        CalculateNormalVectorDerivatives(point_number, kinematics_slave, normal_derivatives_slave, PatchType::Slave);

        Matrix T_master, T_slave;
        CalculateTransformationFromLocalToGlobalCartesian(kinematics_master, T_master);
        CalculateTransformationFromLocalToGlobalCartesian(kinematics_slave, T_slave);

        const Matrix D_master = prod(T_master, Matrix(prod(D_local_master, trans(T_master))));
        const Matrix D_slave = prod(T_slave, Matrix(prod(D_local_slave, trans(T_slave))));

        array_1d<double, 3> local_tangent_master;
        GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent_master);
        const array_1d<double, 3> t_vec_master =
            local_tangent_master[0] * kinematics_master.BaseVector1 + local_tangent_master[1] * kinematics_master.BaseVector2;
        const double surface_jacobian = norm_2(t_vec_master);

        const double integration_weight = integration_points[point_number].Weight();

        const bool same_director_direction = inner_prod(kinematics_master.NormalVector, kinematics_slave.NormalVector) > 0.0;

        // Shift-source (each patch's own Gamma_tilde_h) kinematics 
        const bool has_shift_sources = (mpMasterShiftSource != nullptr && mpSlaveShiftSource != nullptr);
        const auto* p_geometry_master_shift = has_shift_sources ? &mpMasterShiftSource->GetGeometry() : nullptr;
        const auto* p_geometry_slave_shift = has_shift_sources ? &mpSlaveShiftSource->GetGeometry() : nullptr;

        KinematicVariables kinematics_master_shift(3);
        KinematicVariables kinematics_slave_shift(3);
        Matrix normal_derivatives_master_shift, normal_derivatives_slave_shift;
        double surface_jacobian_master_shift = 0.0, surface_jacobian_slave_shift = 0.0;
        double integration_weight_master_shift = 0.0, integration_weight_slave_shift = 0.0;

        if (has_shift_sources) {
            CalculateKinematics(0, kinematics_master_shift, *p_geometry_master_shift);
            CalculateKinematics(0, kinematics_slave_shift, *p_geometry_slave_shift);
            CalculateNormalVectorDerivatives(0, kinematics_master_shift, normal_derivatives_master_shift, *p_geometry_master_shift);
            CalculateNormalVectorDerivatives(0, kinematics_slave_shift, normal_derivatives_slave_shift, *p_geometry_slave_shift);

            array_1d<double, 3> local_tangent_master_shift, local_tangent_slave_shift;
            p_geometry_master_shift->Calculate(LOCAL_TANGENT, local_tangent_master_shift);
            p_geometry_slave_shift->Calculate(LOCAL_TANGENT, local_tangent_slave_shift);
            const array_1d<double, 3> t_vec_master_shift =
                local_tangent_master_shift[0] * kinematics_master_shift.BaseVector1 + local_tangent_master_shift[1] * kinematics_master_shift.BaseVector2;
            const array_1d<double, 3> t_vec_slave_shift =
                local_tangent_slave_shift[0] * kinematics_slave_shift.BaseVector1 + local_tangent_slave_shift[1] * kinematics_slave_shift.BaseVector2;
            surface_jacobian_master_shift = norm_2(t_vec_master_shift);
            surface_jacobian_slave_shift = norm_2(t_vec_slave_shift);
            integration_weight_master_shift = p_geometry_master_shift->IntegrationPoints()[0].Weight();
            integration_weight_slave_shift = p_geometry_slave_shift->IntegrationPoints()[0].Weight();
        }

        const Vector master_w_row = row(r_master_weights, point_number);
        const Vector slave_w_row = row(r_slave_weights, point_number);

        // Rotation-jump penalty 
        if (couple_rotation_x || couple_rotation_y || couple_rotation_z)
        {
            Matrix H_rot_ext = ZeroMatrix(3, mat_size);
            for (IndexType k = 0; k < n_dof_full; ++k) {
                const double jump_w = master_w_row[k] - slave_w_row[k];
                if (jump_w == 0.0) continue;
                const IndexType index = 6 * k;
                if (couple_rotation_x) H_rot_ext(0, index + 3) = jump_w;
                if (couple_rotation_y) H_rot_ext(1, index + 4) = jump_w;
                if (couple_rotation_z) H_rot_ext(2, index + 5) = jump_w;
            }
            noalias(K) += prod(trans(H_rot_ext), H_rot_ext) * (penalty_rotation * integration_weight * surface_jacobian);
        }

        // Through-thickness SBM flux/consistency/adjoint/penalty.
        for (IndexType gauss_index = 0; gauss_index < 2; ++gauss_index)
        {
            const double zeta = gauss_zeta[gauss_index];
            const double zeta_slave = same_director_direction ? zeta : -zeta;

            Matrix jacobian_inv_master, jacobian_inv_slave;
            double jacobian_det_master, jacobian_det_slave;
            CalculateThicknessJacobian(zeta, thickness_master, kinematics_master, normal_derivatives_master, jacobian_inv_master, jacobian_det_master);
            CalculateThicknessJacobian(zeta_slave, thickness_slave, kinematics_slave, normal_derivatives_slave, jacobian_inv_slave, jacobian_det_slave);

            array_1d<double, 3> unit_conormal_master, unit_conormal_slave;
            double area_scale_master, area_scale_slave;
            CalculateLateralConormal(point_number, jacobian_inv_master, jacobian_det_master, PatchType::Master, unit_conormal_master, area_scale_master);
            CalculateLateralConormal(point_number, jacobian_inv_slave, jacobian_det_slave, PatchType::Slave, unit_conormal_slave, area_scale_slave);

            // traction and jump operators
            Matrix B_master, N_op_master;
            CalculateBAndDisplacementOperator(point_number, zeta, thickness_master, jacobian_inv_master, normal_derivatives_master, kinematics_master, PatchType::Master, B_master, N_op_master);
            Matrix B_slave, N_op_slave;
            CalculateBAndDisplacementOperator(point_number, zeta_slave, thickness_slave, jacobian_inv_slave, normal_derivatives_slave, kinematics_slave, PatchType::Slave, B_slave, N_op_slave);

            Matrix Tn_master = ZeroMatrix(3, 6);
            Tn_master(0, 0) = unit_conormal_master[0]; Tn_master(0, 3) = unit_conormal_master[1]; Tn_master(0, 5) = unit_conormal_master[2];
            Tn_master(1, 1) = unit_conormal_master[1]; Tn_master(1, 3) = unit_conormal_master[0]; Tn_master(1, 4) = unit_conormal_master[2];
            Tn_master(2, 2) = unit_conormal_master[2]; Tn_master(2, 4) = unit_conormal_master[1]; Tn_master(2, 5) = unit_conormal_master[0];

            Matrix Tn_slave = ZeroMatrix(3, 6);
            Tn_slave(0, 0) = unit_conormal_slave[0]; Tn_slave(0, 3) = unit_conormal_slave[1]; Tn_slave(0, 5) = unit_conormal_slave[2];
            Tn_slave(1, 1) = unit_conormal_slave[1]; Tn_slave(1, 3) = unit_conormal_slave[0]; Tn_slave(1, 4) = unit_conormal_slave[2];
            Tn_slave(2, 2) = unit_conormal_slave[2]; Tn_slave(2, 4) = unit_conormal_slave[1]; Tn_slave(2, 5) = unit_conormal_slave[0];

            const Matrix F_master = prod(Tn_master, Matrix(prod(D_master, B_master)));
            const Matrix F_slave = prod(Tn_slave, Matrix(prod(D_slave, B_slave)));

            Matrix F_combined = ZeroMatrix(3, mat_size);
            for (SizeType c = 0; c < mat_size_own_master; ++c) {
                column(F_combined, c) = column(F_master, c);
            }
            for (SizeType c = 0; c < mat_size_own - mat_size_own_master; ++c) {
                column(F_combined, c + mat_size_own_master) = -column(F_slave, c);
            }

            // Average-flux numerator {F(u)} = 0.5*(F_master + F_slave)
            Matrix F_avg = ZeroMatrix(3, mat_size);
            for (SizeType c = 0; c < mat_size_own_master; ++c) {
                column(F_avg, c) = column(F_master, c);
            }
            for (SizeType c = 0; c < mat_size_own - mat_size_own_master; ++c) {
                column(F_avg, c + mat_size_own_master) = column(F_slave, c);
            }

            // Jump of the test function using own representation
            Matrix N_jump_own = ZeroMatrix(3, mat_size);
            for (SizeType c = 0; c < mat_size_own_master; ++c) {
                column(N_jump_own, c) = column(N_op_master, c);
            }
            for (SizeType c = 0; c < mat_size_own - mat_size_own_master; ++c) {
                column(N_jump_own, c + mat_size_own_master) = -column(N_op_slave, c);
            }

            // MLS-extended jump operator
            Matrix N_ext_master, N_ext_slave;
            BuildExtendedDisplacementOperator(master_w_row, zeta, thickness_master, kinematics_master.NormalVector, N_ext_master);
            BuildExtendedDisplacementOperator(slave_w_row, zeta_slave, thickness_slave, kinematics_slave.NormalVector, N_ext_slave);
            const Matrix N_extended = N_ext_master - N_ext_slave;

            const double dW = integration_weight * surface_jacobian * area_scale_master * gauss_weight[gauss_index];

            if (has_shift_sources) {
                Matrix jacobian_inv_master_shift, jacobian_inv_slave_shift;
                double jacobian_det_master_shift, jacobian_det_slave_shift;
                CalculateThicknessJacobian(zeta, thickness_master, kinematics_master_shift, normal_derivatives_master_shift, jacobian_inv_master_shift, jacobian_det_master_shift);
                CalculateThicknessJacobian(zeta_slave, thickness_slave, kinematics_slave_shift, normal_derivatives_slave_shift, jacobian_inv_slave_shift, jacobian_det_slave_shift);

                Matrix B_direct_master, N_direct_master;
                CalculateBAndDisplacementOperator(0, zeta, thickness_master, jacobian_inv_master_shift,
                    normal_derivatives_master_shift, kinematics_master_shift, *p_geometry_master_shift, B_direct_master, N_direct_master);
                Matrix B_direct_slave, N_direct_slave;
                CalculateBAndDisplacementOperator(0, zeta_slave, thickness_slave, jacobian_inv_slave_shift,
                    normal_derivatives_slave_shift, kinematics_slave_shift, *p_geometry_slave_shift, B_direct_slave, N_direct_slave);

                array_1d<double, 3> unit_conormal_master_direct, unit_conormal_slave_direct;
                double area_scale_master_direct, area_scale_slave_direct;
                CalculateLateralConormal(0, jacobian_inv_master_shift, jacobian_det_master_shift, *p_geometry_master_shift, unit_conormal_master_direct, area_scale_master_direct);
                unit_conormal_master_direct *= mpMasterShiftSource->GetValue(SBM_SURROGATE_CONORMAL_SIGN);
                CalculateLateralConormal(0, jacobian_inv_slave_shift, jacobian_det_slave_shift, *p_geometry_slave_shift, unit_conormal_slave_direct, area_scale_slave_direct);
                unit_conormal_slave_direct *= mpSlaveShiftSource->GetValue(SBM_SURROGATE_CONORMAL_SIGN);

                Matrix Tn_master_direct = ZeroMatrix(3, 6);
                Tn_master_direct(0, 0) = unit_conormal_master_direct[0]; Tn_master_direct(0, 3) = unit_conormal_master_direct[1]; Tn_master_direct(0, 5) = unit_conormal_master_direct[2];
                Tn_master_direct(1, 1) = unit_conormal_master_direct[1]; Tn_master_direct(1, 3) = unit_conormal_master_direct[0]; Tn_master_direct(1, 4) = unit_conormal_master_direct[2];
                Tn_master_direct(2, 2) = unit_conormal_master_direct[2]; Tn_master_direct(2, 4) = unit_conormal_master_direct[1]; Tn_master_direct(2, 5) = unit_conormal_master_direct[0];
                Matrix Tn_slave_direct = ZeroMatrix(3, 6);
                Tn_slave_direct(0, 0) = unit_conormal_slave_direct[0]; Tn_slave_direct(0, 3) = unit_conormal_slave_direct[1]; Tn_slave_direct(0, 5) = unit_conormal_slave_direct[2];
                Tn_slave_direct(1, 1) = unit_conormal_slave_direct[1]; Tn_slave_direct(1, 3) = unit_conormal_slave_direct[0]; Tn_slave_direct(1, 4) = unit_conormal_slave_direct[2];
                Tn_slave_direct(2, 2) = unit_conormal_slave_direct[2]; Tn_slave_direct(2, 4) = unit_conormal_slave_direct[1]; Tn_slave_direct(2, 5) = unit_conormal_slave_direct[0];

                const Matrix F_direct_master = prod(Tn_master_direct, Matrix(prod(D_master, B_direct_master)));
                const Matrix F_direct_slave = prod(Tn_slave_direct, Matrix(prod(D_slave, B_direct_slave)));

                Matrix F_direct_master_scattered = ZeroMatrix(3, mat_size);
                for (IndexType i = 0; i < p_geometry_master_shift->PointsNumber(); ++i) {
                    const auto it = dof_node_id_to_column.find((*p_geometry_master_shift)[i].Id());
                    KRATOS_ERROR_IF(it == dof_node_id_to_column.end())
                        << "CouplingSbmExtensionOperator6pCondition " << Id()
                        << ": master shift source node " << (*p_geometry_master_shift)[i].Id()
                        << " not found in this condition's SBM_MLS_ALL_DOF_NODES." << std::endl;
                    for (IndexType d = 0; d < 6; ++d) {
                        column(F_direct_master_scattered, 6 * it->second + d) = column(F_direct_master, 6 * i + d);
                    }
                }
                Matrix F_direct_slave_scattered = ZeroMatrix(3, mat_size);
                for (IndexType i = 0; i < p_geometry_slave_shift->PointsNumber(); ++i) {
                    const auto it = dof_node_id_to_column.find((*p_geometry_slave_shift)[i].Id());
                    KRATOS_ERROR_IF(it == dof_node_id_to_column.end())
                        << "CouplingSbmExtensionOperator6pCondition " << Id()
                        << ": slave shift source node " << (*p_geometry_slave_shift)[i].Id()
                        << " not found in this condition's SBM_MLS_ALL_DOF_NODES." << std::endl;
                    for (IndexType d = 0; d < 6; ++d) {
                        column(F_direct_slave_scattered, 6 * it->second + d) = column(F_direct_slave, 6 * i + d);
                    }
                }

                const double dW_master_direct = integration_weight_master_shift * surface_jacobian_master_shift * area_scale_master_direct * gauss_weight[gauss_index];
                const double dW_slave_direct = integration_weight_slave_shift * surface_jacobian_slave_shift * area_scale_slave_direct * gauss_weight[gauss_index];

                noalias(K) += prod(trans(N_ext_master), F_direct_master_scattered) * dW_master_direct;
                noalias(K) += prod(trans(N_ext_slave), F_direct_slave_scattered) * dW_slave_direct;
            } else {
                noalias(K) += prod(trans(F_combined), N_extended) * (-0.5 * dW);
                noalias(K) += prod(trans(N_jump_own), F_avg) * (-0.5 * dW);
            }

            noalias(K) += prod(trans(N_extended), N_extended) * (stabilization_parameter * dW);
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

        const auto& r_all_dof_nodes = GetValue(SBM_MLS_ALL_DOF_NODES);
        Vector d_current = ZeroVector(mat_size);
        for (IndexType k = 0; k < n_dof_full; ++k) {
            const auto& r_node = *r_all_dof_nodes[k];
            const array_1d<double, 3>& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& rot = r_node.FastGetSolutionStepValue(ROTATION);
            const IndexType index = 6 * k;
            d_current[index] = disp[0]; d_current[index + 1] = disp[1]; d_current[index + 2] = disp[2];
            d_current[index + 3] = rot[0]; d_current[index + 4] = rot[1]; d_current[index + 5] = rot[2];
        }

        noalias(rRightHandSideVector) = -prod(K, d_current);
    }

    KRATOS_CATCH("")
}

int CouplingSbmExtensionOperator6pCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(NITSCHE_STABILIZATION_FACTOR))
        << "No NITSCHE_STABILIZATION_FACTOR defined in property of CouplingSbmExtensionOperator6pCondition" << std::endl;

    KRATOS_ERROR_IF(GetProperties().GetSubProperties().size() < 1)
        << "CouplingSbmExtensionOperator6pCondition needs master/slave material sub-properties (YOUNG_MODULUS, "
        << "POISSON_RATIO, THICKNESS) attached to its Properties" << std::endl;

    KRATOS_ERROR_IF_NOT(Has(SBM_MLS_ALL_DOF_NODES))
        << "CouplingSbmExtensionOperator6pCondition " << Id() << ": SBM_MLS_ALL_DOF_NODES not set -- "
        << "IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingExtensionOperators "
        << "must be called on this condition's ModelPart before the solve." << std::endl;

    KRATOS_ERROR_IF_NOT(Has(SBM_MLS_MASTER_WEIGHTS) && Has(SBM_MLS_SLAVE_WEIGHTS))
        << "CouplingSbmExtensionOperator6pCondition " << Id() << ": SBM_MLS_MASTER_WEIGHTS/SBM_MLS_SLAVE_WEIGHTS not set -- "
        << "IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingExtensionOperators "
        << "must be called on this condition's ModelPart before the solve." << std::endl;

    return 0;
}

void CouplingSbmExtensionOperator6pCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_all_dof_nodes = GetValue(SBM_MLS_ALL_DOF_NODES);
    const SizeType n_dof_full = r_all_dof_nodes.size();

    if (rResult.size() != 6 * n_dof_full)
        rResult.resize(6 * n_dof_full, false);

    for (IndexType k = 0; k < n_dof_full; ++k) {
        const IndexType index = 6 * k;
        const auto& r_node = *r_all_dof_nodes[k];
        rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
    }

    KRATOS_CATCH("")
}

void CouplingSbmExtensionOperator6pCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_all_dof_nodes = GetValue(SBM_MLS_ALL_DOF_NODES);

    rElementalDofList.resize(0);
    rElementalDofList.reserve(6 * r_all_dof_nodes.size());

    for (const auto& p_node : r_all_dof_nodes) {
        rElementalDofList.push_back(p_node->pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(p_node->pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(p_node->pGetDof(DISPLACEMENT_Z));
        rElementalDofList.push_back(p_node->pGetDof(ROTATION_X));
        rElementalDofList.push_back(p_node->pGetDof(ROTATION_Y));
        rElementalDofList.push_back(p_node->pGetDof(ROTATION_Z));
    }

    KRATOS_CATCH("")
}

} // Namespace Kratos
