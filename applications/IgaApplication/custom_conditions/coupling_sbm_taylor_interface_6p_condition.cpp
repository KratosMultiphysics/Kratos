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
#include <cmath>

// External includes

// Project includes
#include "custom_conditions/coupling_sbm_taylor_interface_6p_condition.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace {

using IndexType = CouplingSbmTaylorInterface6pCondition::IndexType;
using SizeType = CouplingSbmTaylorInterface6pCondition::SizeType;

double ComputeTaylorTerm2D(double derivative, double dx, IndexType n_k, double dy, IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k)
         / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

array_1d<double, 3> PhysicalCenter(const Geometry<Node>& rGeometry)
{
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    const SizeType n = rGeometry.PointsNumber();
    array_1d<double, 3> center = ZeroVector(3);
    for (IndexType i = 0; i < n; ++i) {
        center[0] += r_N(0, i) * rGeometry[i].X0();
        center[1] += r_N(0, i) * rGeometry[i].Y0();
        center[2] += r_N(0, i) * rGeometry[i].Z0();
    }
    return center;
}

// Basis-function polynomial order p
SizeType InferCappedTaylorOrder(const Geometry<Node>& rGeometry, GeometryData::IntegrationMethod ThisMethod)
{
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(ThisMethod);
    const SizeType n_cp = r_DN_De[0].size1();
    const SizeType p = static_cast<SizeType>(std::lround(std::sqrt(static_cast<double>(n_cp)))) - 1;
    return std::min(p, static_cast<SizeType>(2));
}

// Taylor-extrapolates rGeometry's own shape function values from its own
// (single) evaluation point out to a target a distance rDistanceVector away.
void ComputeTaylorExpansionValues(
    const Geometry<Node>& rGeometry,
    const Vector& rDistanceVector,
    SizeType TaylorOrder,
    GeometryData::IntegrationMethod ThisMethod,
    Vector& rHSum)
{
    const SizeType n_cp = rGeometry.PointsNumber();
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    const double dx = rDistanceVector[0];
    const double dy = rDistanceVector[1];

    rHSum = ZeroVector(n_cp);

    if (TaylorOrder == 0) {
        for (IndexType i = 0; i < n_cp; ++i) rHSum[i] = r_N(0, i);
        return;
    }

    std::vector<Matrix> shape_function_derivatives(TaylorOrder);
    for (IndexType n = 1; n <= TaylorOrder; ++n) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, ThisMethod);
    }

    for (IndexType i = 0; i < n_cp; ++i) {
        double taylor_term = 0.0;
        for (IndexType n = 1; n <= TaylorOrder; ++n) {
            const Matrix& r_deriv = shape_function_derivatives[n - 1];
            for (IndexType k = 0; k <= n; ++k) {
                const IndexType n_k = n - k;
                taylor_term += ComputeTaylorTerm2D(r_deriv(i, k), dx, n_k, dy, k);
            }
        }
        rHSum[i] = taylor_term + r_N(0, i);
    }
}

// Taylor-shifted GRADIENT (d/dxi, d/deta) of rGeometry's own shape functions
void ComputeTaylorExpansionGradients(
    const Geometry<Node>& rGeometry,
    const Vector& rDistanceVector,
    SizeType TaylorOrder,
    GeometryData::IntegrationMethod ThisMethod,
    Vector& rGradXiSum,
    Vector& rGradEtaSum)
{
    const SizeType n_cp = rGeometry.PointsNumber();
    const double dx = rDistanceVector[0];
    const double dy = rDistanceVector[1];

    rGradXiSum = ZeroVector(n_cp);
    rGradEtaSum = ZeroVector(n_cp);

    if (TaylorOrder == 0) {
        const auto& r_DN_De_local = rGeometry.ShapeFunctionsLocalGradients(ThisMethod)[0];
        for (IndexType i = 0; i < n_cp; ++i) {
            rGradXiSum[i] = r_DN_De_local(i, 0);
            rGradEtaSum[i] = r_DN_De_local(i, 1);
        }
        return;
    }

    std::vector<Matrix> shape_function_derivatives(TaylorOrder);
    for (IndexType n = 1; n <= TaylorOrder; ++n) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, ThisMethod);
    }

    const Matrix& r_deriv1 = shape_function_derivatives[0];

    for (IndexType i = 0; i < n_cp; ++i) {
        double gxi_term = r_deriv1(i, 0);
        double geta_term = r_deriv1(i, 1);
        for (IndexType m = 1; m <= TaylorOrder - 1; ++m) {
            const Matrix& r_deriv_mp1 = shape_function_derivatives[m];
            for (IndexType k = 0; k <= m; ++k) {
                const IndexType m_k = m - k;
                const double coeff = std::pow(dx, m_k) * std::pow(dy, k)
                                    / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(m_k));
                gxi_term += r_deriv_mp1(i, k) * coeff;
                geta_term += r_deriv_mp1(i, k + 1) * coeff;
            }
        }
        rGradXiSum[i] = gxi_term;
        rGradEtaSum[i] = geta_term;
    }
}

// Builds the "displacement-at-zeta" operator from Taylor VALUE
// weights instead of real shape functions. 
void BuildTaylorDisplacementOperator(
    const Vector& rWeights,
    double zeta,
    double Thickness,
    const array_1d<double, 3>& rNormalVector,
    Matrix& rNOperator)
{
    const SizeType n = rWeights.size();
    const SizeType mat_size = 6 * n;
    rNOperator = ZeroMatrix(3, mat_size);

    const double normal_x = rNormalVector[0];
    const double normal_y = rNormalVector[1];
    const double normal_z = rNormalVector[2];
    const double half_t_zeta = zeta * (Thickness / 2.0);

    for (IndexType k = 0; k < n; ++k) {
        const double w = rWeights[k];
        if (w == 0.0) continue;
        const IndexType index = 6 * k;
        rNOperator(0, index)     = w;
        rNOperator(1, index + 1) = w;
        rNOperator(2, index + 2) = w;
        const double coeff = w * half_t_zeta;
        rNOperator(0, index + 4) =  coeff * normal_z;
        rNOperator(0, index + 5) = -coeff * normal_y;
        rNOperator(1, index + 3) = -coeff * normal_z;
        rNOperator(1, index + 5) =  coeff * normal_x;
        rNOperator(2, index + 3) =  coeff * normal_y;
        rNOperator(2, index + 4) = -coeff * normal_x;
    }
}

// Builds the (6 x 6*n) Taylor-shifted B-operator (strain operator) from
// shifted VALUE and GRADIENT weights. 
void BuildTaylorBOperator(
    const Vector& rWeights,
    const Vector& rGradXi,
    const Vector& rGradEta,
    double zeta,
    double Thickness,
    const Matrix& rJacobianInv,
    const Matrix& rNormalVectorDerivatives,
    const array_1d<double, 3>& rNormalVector,
    Matrix& rBOperator)
{
    const SizeType n = rWeights.size();
    const SizeType mat_size = 6 * n;
    rBOperator = ZeroMatrix(6, mat_size);

    const double normal_x = rNormalVector[0];
    const double normal_y = rNormalVector[1];
    const double normal_z = rNormalVector[2];

    const Matrix normal_vector_derivatives_global = prod(rJacobianInv, rNormalVectorDerivatives);
    const double d_normal_x_dx = normal_vector_derivatives_global(0, 0);
    const double d_normal_y_dx = normal_vector_derivatives_global(0, 1);
    const double d_normal_z_dx = normal_vector_derivatives_global(0, 2);
    const double d_normal_x_dy = normal_vector_derivatives_global(1, 0);
    const double d_normal_y_dy = normal_vector_derivatives_global(1, 1);
    const double d_normal_z_dy = normal_vector_derivatives_global(1, 2);
    const double d_normal_x_dz = normal_vector_derivatives_global(2, 0);
    const double d_normal_y_dz = normal_vector_derivatives_global(2, 1);
    const double d_normal_z_dz = normal_vector_derivatives_global(2, 2);

    const double d_zeta_dx = rJacobianInv(0, 2);
    const double d_zeta_dy = rJacobianInv(1, 2);
    const double d_zeta_dz = rJacobianInv(2, 2);

    for (IndexType i = 0; i < n; ++i) {
        const double Ni = rWeights[i];
        const double grad_xi = rGradXi[i];
        const double grad_eta = rGradEta[i];
        if (Ni == 0.0 && grad_xi == 0.0 && grad_eta == 0.0) continue;

        const double dNi_dx = rJacobianInv(0, 0) * grad_xi + rJacobianInv(0, 1) * grad_eta;
        const double dNi_dy = rJacobianInv(1, 0) * grad_xi + rJacobianInv(1, 1) * grad_eta;
        const double dNi_dz = rJacobianInv(2, 0) * grad_xi + rJacobianInv(2, 1) * grad_eta;

        const IndexType index = 6 * i;

        rBOperator(0, index)     = dNi_dx;
        rBOperator(1, index + 1) = dNi_dy;
        rBOperator(2, index + 2) = dNi_dz;
        rBOperator(3, index)     = dNi_dy;
        rBOperator(3, index + 1) = dNi_dx;
        rBOperator(4, index + 1) = dNi_dz;
        rBOperator(4, index + 2) = dNi_dy;
        rBOperator(5, index)     = dNi_dz;
        rBOperator(5, index + 2) = dNi_dx;

        rBOperator(0, index + 4) =  ((dNi_dx * zeta * normal_z) + (Ni * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (Thickness / 2.0);
        rBOperator(0, index + 5) = -((dNi_dx * zeta * normal_y) + (Ni * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (Thickness / 2.0);

        rBOperator(1, index + 3) = -((dNi_dy * zeta * normal_z) + (Ni * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (Thickness / 2.0);
        rBOperator(1, index + 5) =  ((dNi_dy * zeta * normal_x) + (Ni * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (Thickness / 2.0);

        rBOperator(2, index + 3) =  ((dNi_dz * zeta * normal_y) + (Ni * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (Thickness / 2.0);
        rBOperator(2, index + 4) = -((dNi_dz * zeta * normal_x) + (Ni * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (Thickness / 2.0);

        rBOperator(3, index + 3) = -((dNi_dx * zeta * normal_z) + (Ni * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (Thickness / 2.0);
        rBOperator(3, index + 4) =  ((dNi_dy * zeta * normal_z) + (Ni * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (Thickness / 2.0);
        rBOperator(3, index + 5) = (((dNi_dx * zeta * normal_x) + (Ni * (zeta * d_normal_x_dx + d_zeta_dx * normal_x)))
                                   -((dNi_dy * zeta * normal_y) + (Ni * (zeta * d_normal_y_dy + d_zeta_dy * normal_y)))) * (Thickness / 2.0);

        rBOperator(4, index + 3) = (((dNi_dy * zeta * normal_y) + (Ni * (zeta * d_normal_y_dy + d_zeta_dy * normal_y)))
                                   -((dNi_dz * zeta * normal_z) + (Ni * (zeta * d_normal_z_dz + d_zeta_dz * normal_z)))) * (Thickness / 2.0);
        rBOperator(4, index + 4) = -((dNi_dy * zeta * normal_x) + (Ni * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (Thickness / 2.0);
        rBOperator(4, index + 5) =  ((dNi_dz * zeta * normal_x) + (Ni * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (Thickness / 2.0);

        rBOperator(5, index + 3) =  ((dNi_dx * zeta * normal_y) + (Ni * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (Thickness / 2.0);
        rBOperator(5, index + 4) = (((dNi_dz * zeta * normal_z) + (Ni * (zeta * d_normal_z_dz + d_zeta_dz * normal_z)))
                                   -((dNi_dx * zeta * normal_x) + (Ni * (zeta * d_normal_x_dx + d_zeta_dx * normal_x)))) * (Thickness / 2.0);
        rBOperator(5, index + 5) = -((dNi_dz * zeta * normal_y) + (Ni * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (Thickness / 2.0);
    }
}

}

void CouplingSbmTaylorInterface6pCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mpMasterShiftSource == nullptr || mpSlaveShiftSource == nullptr)
        << "CouplingSbmTaylorInterface6pCondition " << Id()
        << ": shift sources not set -- call SetShiftSources() before Initialize()." << std::endl;

    const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
    const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

    const array_1d<double, 3> target_master = PhysicalCenter(r_geometry_master);
    const array_1d<double, 3> target_slave = PhysicalCenter(r_geometry_slave);

    const array_1d<double, 3> center_master_shift = PhysicalCenter(mpMasterShiftSource->GetGeometry());
    const array_1d<double, 3> center_slave_shift = PhysicalCenter(mpSlaveShiftSource->GetGeometry());

    mDistanceVectorMaster.resize(3);
    mDistanceVectorMaster[0] = target_master[0] - center_master_shift[0];
    mDistanceVectorMaster[1] = target_master[1] - center_master_shift[1];
    mDistanceVectorMaster[2] = target_master[2] - center_master_shift[2];

    mDistanceVectorSlave.resize(3);
    mDistanceVectorSlave[0] = target_slave[0] - center_slave_shift[0];
    mDistanceVectorSlave[1] = target_slave[1] - center_slave_shift[1];
    mDistanceVectorSlave[2] = target_slave[2] - center_slave_shift[2];

    mTaylorOrderMaster = InferCappedTaylorOrder(mpMasterShiftSource->GetGeometry(), mpMasterShiftSource->GetIntegrationMethod());
    mTaylorOrderSlave = InferCappedTaylorOrder(mpSlaveShiftSource->GetGeometry(), mpSlaveShiftSource->GetIntegrationMethod());

    KRATOS_CATCH("")
}

void CouplingSbmTaylorInterface6pCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const double stabilization_parameter = GetProperties()[NITSCHE_STABILIZATION_FACTOR];

    const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);

    const auto& r_geometry_master_shift = mpMasterShiftSource->GetGeometry();
    const auto& r_geometry_slave_shift = mpSlaveShiftSource->GetGeometry();
    const SizeType n_master = r_geometry_master_shift.PointsNumber();
    const SizeType n_slave = r_geometry_slave_shift.PointsNumber();
    const SizeType mat_size_master = 6 * n_master;
    const SizeType mat_size = 6 * (n_master + n_slave);

    Matrix K = ZeroMatrix(mat_size, mat_size);

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

    const Properties& r_props_master = GetProperties().GetSubProperties().front();
    const Properties& r_props_slave = GetProperties().GetSubProperties().back();
    const double thickness_master = r_props_master[THICKNESS];
    const double thickness_slave = r_props_slave[THICKNESS];

    Matrix D_local_master, D_local_slave;
    CalculateLocalConstitutiveMatrix(r_props_master, D_local_master);
    CalculateLocalConstitutiveMatrix(r_props_slave, D_local_slave);

    const double gauss_zeta[2] = { -std::sqrt(1.0 / 3.0), std::sqrt(1.0 / 3.0) };
    const double gauss_weight[2] = { 1.0, 1.0 };

    Vector H_sum_master, H_sum_slave;
    ComputeTaylorExpansionValues(r_geometry_master_shift, mDistanceVectorMaster, mTaylorOrderMaster, mpMasterShiftSource->GetIntegrationMethod(), H_sum_master);
    ComputeTaylorExpansionValues(r_geometry_slave_shift, mDistanceVectorSlave, mTaylorOrderSlave, mpSlaveShiftSource->GetIntegrationMethod(), H_sum_slave);

    Vector GradXi_master, GradEta_master, GradXi_slave, GradEta_slave;
    ComputeTaylorExpansionGradients(r_geometry_master_shift, mDistanceVectorMaster, mTaylorOrderMaster, mpMasterShiftSource->GetIntegrationMethod(), GradXi_master, GradEta_master);
    ComputeTaylorExpansionGradients(r_geometry_slave_shift, mDistanceVectorSlave, mTaylorOrderSlave, mpSlaveShiftSource->GetIntegrationMethod(), GradXi_slave, GradEta_slave);

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        // True boundary's own kinematics 
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
        r_geometry_master.Calculate(LOCAL_TANGENT, local_tangent_master);
        const array_1d<double, 3> t_vec_master =
            local_tangent_master[0] * kinematics_master.BaseVector1 + local_tangent_master[1] * kinematics_master.BaseVector2;
        const double surface_jacobian = norm_2(t_vec_master);

        const double integration_weight = integration_points[point_number].Weight();

        const bool same_director_direction = inner_prod(kinematics_master.NormalVector, kinematics_slave.NormalVector) > 0.0;

        // Rotation-jump penalty 
        const bool couple_rotation_x = Is(IgaFlags::FIX_ROTATION_X);
        const bool couple_rotation_y = Is(IgaFlags::FIX_ROTATION_Y);
        const bool couple_rotation_z = Is(IgaFlags::FIX_ROTATION_Z);
        if (couple_rotation_x || couple_rotation_y || couple_rotation_z) {
            const double penalty_rotation = GetProperties()[PENALTY_ROTATION_FACTOR];

            Matrix H_rot = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < n_master; ++i) {
                const IndexType index = 6 * i;
                if (couple_rotation_x) H_rot(0, index + 3) = H_sum_master[i];
                if (couple_rotation_y) H_rot(1, index + 4) = H_sum_master[i];
                if (couple_rotation_z) H_rot(2, index + 5) = H_sum_master[i];
            }
            for (IndexType i = 0; i < n_slave; ++i) {
                const IndexType index = mat_size_master + 6 * i;
                if (couple_rotation_x) H_rot(0, index + 3) = -H_sum_slave[i];
                if (couple_rotation_y) H_rot(1, index + 4) = -H_sum_slave[i];
                if (couple_rotation_z) H_rot(2, index + 5) = -H_sum_slave[i];
            }

            noalias(K) += prod(trans(H_rot), H_rot) * (penalty_rotation * integration_weight * surface_jacobian);
        }

        // Shift-source kinematics 
        KinematicVariables kinematics_master_shift(3);
        KinematicVariables kinematics_slave_shift(3);
        CalculateKinematics(0, kinematics_master_shift, r_geometry_master_shift);
        CalculateKinematics(0, kinematics_slave_shift, r_geometry_slave_shift);
        Matrix normal_derivatives_master_shift, normal_derivatives_slave_shift;
        CalculateNormalVectorDerivatives(0, kinematics_master_shift, normal_derivatives_master_shift, r_geometry_master_shift);
        CalculateNormalVectorDerivatives(0, kinematics_slave_shift, normal_derivatives_slave_shift, r_geometry_slave_shift);

        // Jacobian and quadrature weight
        array_1d<double, 3> local_tangent_master_shift, local_tangent_slave_shift;
        r_geometry_master_shift.Calculate(LOCAL_TANGENT, local_tangent_master_shift);
        r_geometry_slave_shift.Calculate(LOCAL_TANGENT, local_tangent_slave_shift);
        const array_1d<double, 3> t_vec_master_shift =
            local_tangent_master_shift[0] * kinematics_master_shift.BaseVector1 + local_tangent_master_shift[1] * kinematics_master_shift.BaseVector2;
        const array_1d<double, 3> t_vec_slave_shift =
            local_tangent_slave_shift[0] * kinematics_slave_shift.BaseVector1 + local_tangent_slave_shift[1] * kinematics_slave_shift.BaseVector2;
        const double surface_jacobian_master_shift = norm_2(t_vec_master_shift);
        const double surface_jacobian_slave_shift = norm_2(t_vec_slave_shift);
        const double integration_weight_master_shift = r_geometry_master_shift.IntegrationPoints()[0].Weight();
        const double integration_weight_slave_shift = r_geometry_slave_shift.IntegrationPoints()[0].Weight();

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

            Matrix Tn_master = ZeroMatrix(3, 6);
            Tn_master(0, 0) = unit_conormal_master[0]; Tn_master(0, 3) = unit_conormal_master[1]; Tn_master(0, 5) = unit_conormal_master[2];
            Tn_master(1, 1) = unit_conormal_master[1]; Tn_master(1, 3) = unit_conormal_master[0]; Tn_master(1, 4) = unit_conormal_master[2];
            Tn_master(2, 2) = unit_conormal_master[2]; Tn_master(2, 4) = unit_conormal_master[1]; Tn_master(2, 5) = unit_conormal_master[0];

            Matrix Tn_slave = ZeroMatrix(3, 6);
            Tn_slave(0, 0) = unit_conormal_slave[0]; Tn_slave(0, 3) = unit_conormal_slave[1]; Tn_slave(0, 5) = unit_conormal_slave[2];
            Tn_slave(1, 1) = unit_conormal_slave[1]; Tn_slave(1, 3) = unit_conormal_slave[0]; Tn_slave(1, 4) = unit_conormal_slave[2];
            Tn_slave(2, 2) = unit_conormal_slave[2]; Tn_slave(2, 4) = unit_conormal_slave[1]; Tn_slave(2, 5) = unit_conormal_slave[0];

            // Shift-source-frame Jacobians
            Matrix jacobian_inv_master_shift, jacobian_inv_slave_shift;
            double jacobian_det_master_shift, jacobian_det_slave_shift;
            CalculateThicknessJacobian(zeta, thickness_master, kinematics_master_shift, normal_derivatives_master_shift, jacobian_inv_master_shift, jacobian_det_master_shift);
            CalculateThicknessJacobian(zeta_slave, thickness_slave, kinematics_slave_shift, normal_derivatives_slave_shift, jacobian_inv_slave_shift, jacobian_det_slave_shift);

            Matrix sub_master, sub_slave;
            BuildTaylorDisplacementOperator(H_sum_master, zeta, thickness_master, kinematics_master.NormalVector, sub_master);
            BuildTaylorDisplacementOperator(H_sum_slave, zeta_slave, thickness_slave, kinematics_slave.NormalVector, sub_slave);

            Matrix B_master_shifted, B_slave_shifted;
            BuildTaylorBOperator(H_sum_master, GradXi_master, GradEta_master, zeta, thickness_master,
                jacobian_inv_master_shift, normal_derivatives_master_shift, kinematics_master.NormalVector, B_master_shifted);
            BuildTaylorBOperator(H_sum_slave, GradXi_slave, GradEta_slave, zeta_slave, thickness_slave,
                jacobian_inv_slave_shift, normal_derivatives_slave_shift, kinematics_slave.NormalVector, B_slave_shifted);

            const Matrix F_master_shifted = prod(Tn_master, Matrix(prod(D_master, B_master_shifted)));
            const Matrix F_slave_shifted = prod(Tn_slave, Matrix(prod(D_slave, B_slave_shifted)));

            const double dW = integration_weight * surface_jacobian * area_scale_master * gauss_weight[gauss_index];

            Matrix N_combined = ZeroMatrix(3, mat_size);
            for (SizeType c = 0; c < mat_size_master; ++c) column(N_combined, c) = column(sub_master, c);
            for (SizeType c = 0; c < mat_size - mat_size_master; ++c) column(N_combined, c + mat_size_master) = -column(sub_slave, c);

            // true boundary adjoint term
            Matrix B_direct_master, N_direct_master;
            CalculateBAndDisplacementOperator(0, zeta, thickness_master, jacobian_inv_master_shift,
                normal_derivatives_master_shift, kinematics_master_shift, r_geometry_master_shift, B_direct_master, N_direct_master);
            Matrix B_direct_slave, N_direct_slave;
            CalculateBAndDisplacementOperator(0, zeta_slave, thickness_slave, jacobian_inv_slave_shift,
                normal_derivatives_slave_shift, kinematics_slave_shift, r_geometry_slave_shift, B_direct_slave, N_direct_slave);

            array_1d<double, 3> unit_conormal_master_direct, unit_conormal_slave_direct;
            double area_scale_master_direct, area_scale_slave_direct;
            CalculateLateralConormal(0, jacobian_inv_master_shift, jacobian_det_master_shift, r_geometry_master_shift, unit_conormal_master_direct, area_scale_master_direct);
            unit_conormal_master_direct *= mpMasterShiftSource->GetValue(SBM_SURROGATE_CONORMAL_SIGN);
            CalculateLateralConormal(0, jacobian_inv_slave_shift, jacobian_det_slave_shift, r_geometry_slave_shift, unit_conormal_slave_direct, area_scale_slave_direct);
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

            Matrix N_shifted_master_only = ZeroMatrix(3, mat_size);
            Matrix F_direct_master_only = ZeroMatrix(3, mat_size);
            for (SizeType c = 0; c < mat_size_master; ++c) {
                column(N_shifted_master_only, c) = column(sub_master, c);
                column(F_direct_master_only, c) = column(F_direct_master, c);
            }
            Matrix N_shifted_slave_only = ZeroMatrix(3, mat_size);
            Matrix F_direct_slave_only = ZeroMatrix(3, mat_size);
            for (SizeType c = 0; c < mat_size - mat_size_master; ++c) {
                column(N_shifted_slave_only, c + mat_size_master) = column(sub_slave, c);
                column(F_direct_slave_only, c + mat_size_master) = column(F_direct_slave, c);
            }

            const double dW_master_direct = integration_weight_master_shift * surface_jacobian_master_shift * area_scale_master_direct * gauss_weight[gauss_index];
            const double dW_slave_direct = integration_weight_slave_shift * surface_jacobian_slave_shift * area_scale_slave_direct * gauss_weight[gauss_index];

            noalias(K) += prod(trans(N_shifted_master_only), F_direct_master_only) * dW_master_direct;
            noalias(K) += prod(trans(N_shifted_slave_only), F_direct_slave_only) * dW_slave_direct;

            noalias(K) += prod(trans(N_combined), N_combined) * (stabilization_parameter * dW);
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
        for (IndexType i = 0; i < n_master; ++i) {
            const array_1d<double, 3>& disp = r_geometry_master_shift[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& rot = r_geometry_master_shift[i].FastGetSolutionStepValue(ROTATION);
            const IndexType index = 6 * i;
            d_current[index] = disp[0]; d_current[index + 1] = disp[1]; d_current[index + 2] = disp[2];
            d_current[index + 3] = rot[0]; d_current[index + 4] = rot[1]; d_current[index + 5] = rot[2];
        }
        for (IndexType i = 0; i < n_slave; ++i) {
            const array_1d<double, 3>& disp = r_geometry_slave_shift[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& rot = r_geometry_slave_shift[i].FastGetSolutionStepValue(ROTATION);
            const IndexType index = mat_size_master + 6 * i;
            d_current[index] = disp[0]; d_current[index + 1] = disp[1]; d_current[index + 2] = disp[2];
            d_current[index + 3] = rot[0]; d_current[index + 4] = rot[1]; d_current[index + 5] = rot[2];
        }

        noalias(rRightHandSideVector) = -prod(K, d_current);
    }

    KRATOS_CATCH("")
}

int CouplingSbmTaylorInterface6pCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(NITSCHE_STABILIZATION_FACTOR))
        << "No NITSCHE_STABILIZATION_FACTOR defined in property of CouplingSbmTaylorInterface6pCondition" << std::endl;
    KRATOS_ERROR_IF(GetProperties().GetSubProperties().size() < 1)
        << "CouplingSbmTaylorInterface6pCondition " << Id() << ": expected at least one master/slave sub-property "
        << "(front()/back() may coincide when both patches share the same material)." << std::endl;
    KRATOS_ERROR_IF(mpMasterShiftSource == nullptr || mpSlaveShiftSource == nullptr)
        << "CouplingSbmTaylorInterface6pCondition " << Id()
        << ": shift sources not set -- call SetShiftSources() before Check()/Initialize()." << std::endl;
    return 0;
}

void CouplingSbmTaylorInterface6pCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geometry_master_shift = mpMasterShiftSource->GetGeometry();
    const auto& r_geometry_slave_shift = mpSlaveShiftSource->GetGeometry();
    const SizeType n_master = r_geometry_master_shift.size();
    const SizeType n_slave = r_geometry_slave_shift.size();

    if (rResult.size() != 6 * (n_master + n_slave))
        rResult.resize(6 * (n_master + n_slave), false);

    for (IndexType i = 0; i < n_master; ++i) {
        const IndexType index = i * 6;
        const auto& r_node = r_geometry_master_shift[i];
        rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
    }
    for (IndexType i = 0; i < n_slave; ++i) {
        const IndexType index = 6 * n_master + i * 6;
        const auto& r_node = r_geometry_slave_shift[i];
        rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
    }

    KRATOS_CATCH("")
}

void CouplingSbmTaylorInterface6pCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geometry_master_shift = mpMasterShiftSource->GetGeometry();
    const auto& r_geometry_slave_shift = mpSlaveShiftSource->GetGeometry();
    const SizeType n_master = r_geometry_master_shift.size();
    const SizeType n_slave = r_geometry_slave_shift.size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(6 * (n_master + n_slave));

    for (IndexType i = 0; i < n_master; ++i) {
        const auto& r_node = r_geometry_master_shift.GetPoint(i);
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
        rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
        rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
    }
    for (IndexType i = 0; i < n_slave; ++i) {
        const auto& r_node = r_geometry_slave_shift.GetPoint(i);
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
