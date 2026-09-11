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
#include "custom_conditions/coupling_nitsche_6p_condition.h"
#include "utilities/math_utils.h"

namespace Kratos
{
    void CouplingNitsche6pCondition::CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const PatchType& rPatch) const
    {
        const IndexType geometry_part = (rPatch == PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(geometry_part);
        const SizeType number_of_nodes = r_geometry.size();

        const auto& r_shape_functions_gradients = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
        const Matrix& r_DN_De = r_shape_functions_gradients[IntegrationPointIndex];

        KRATOS_ERROR_IF(r_DN_De.size1() != number_of_nodes || r_DN_De.size2() < 2)
            << "CouplingNitsche6pCondition: unexpected 1st-derivative shape function matrix size ("
            << r_DN_De.size1() << "x" << r_DN_De.size2() << "), expected (" << number_of_nodes << "x>=2)" << std::endl;

        noalias(rKinematicVariables.BaseVector1) = ZeroVector(3);
        noalias(rKinematicVariables.BaseVector2) = ZeroVector(3);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3>& r_coordinates = r_geometry[i].GetInitialPosition();
            for (IndexType k = 0; k < 3; ++k) {
                rKinematicVariables.BaseVector1[k] += r_coordinates[k] * r_DN_De(i, 0);
                rKinematicVariables.BaseVector2[k] += r_coordinates[k] * r_DN_De(i, 1);
            }
        }

        MathUtils<double>::CrossProduct(rKinematicVariables.NormalVectorTilde,
            rKinematicVariables.BaseVector1, rKinematicVariables.BaseVector2);

        rKinematicVariables.DifferentialArea = norm_2(rKinematicVariables.NormalVectorTilde);

        noalias(rKinematicVariables.NormalVector) =
            rKinematicVariables.NormalVectorTilde / rKinematicVariables.DifferentialArea;
    }

    void CouplingNitsche6pCondition::CalculateNormalVectorDerivatives(
        IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& rNormalVectorDerivatives,
        const PatchType& rPatch) const
    {
        const IndexType geometry_part = (rPatch == PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(geometry_part);

        const Matrix& r_DDN_DDe = r_geometry.ShapeFunctionDerivatives(
            2, IntegrationPointIndex, r_geometry.GetDefaultIntegrationMethod());

        KRATOS_ERROR_IF(r_DDN_DDe.size1() != r_geometry.size() || r_DDN_DDe.size2() != 3)
            << "CouplingNitsche6pCondition: unexpected 2nd-derivative shape function matrix size ("
            << r_DDN_DDe.size1() << "x" << r_DDN_DDe.size2() << "), expected (" << r_geometry.size() << "x3)" << std::endl;

        const double inv_differential_area = 1.0 / rKinematicVariables.DifferentialArea;
        const double inv_differential_area_cube = 1.0 / std::pow(rKinematicVariables.DifferentialArea, 3);

        array_1d<double, 3> base_vector1_derivative_11 = ZeroVector(3);
        array_1d<double, 3> base_vector1_derivative_12 = ZeroVector(3);
        array_1d<double, 3> base_vector2_derivative_22 = ZeroVector(3);

        for (IndexType i = 0; i < r_geometry.size(); ++i)
        {
            const array_1d<double, 3>& r_coordinates = r_geometry[i].GetInitialPosition();

            base_vector1_derivative_11[0] += r_coordinates[0] * r_DDN_DDe(i, 0);
            base_vector1_derivative_11[1] += r_coordinates[1] * r_DDN_DDe(i, 0);
            base_vector1_derivative_11[2] += r_coordinates[2] * r_DDN_DDe(i, 0);

            base_vector1_derivative_12[0] += r_coordinates[0] * r_DDN_DDe(i, 1);
            base_vector1_derivative_12[1] += r_coordinates[1] * r_DDN_DDe(i, 1);
            base_vector1_derivative_12[2] += r_coordinates[2] * r_DDN_DDe(i, 1);

            base_vector2_derivative_22[0] += r_coordinates[0] * r_DDN_DDe(i, 2);
            base_vector2_derivative_22[1] += r_coordinates[1] * r_DDN_DDe(i, 2);
            base_vector2_derivative_22[2] += r_coordinates[2] * r_DDN_DDe(i, 2);
        }

        array_1d<double, 3> normal_tilde_derivative_1_term1, normal_tilde_derivative_1_term2, normal_tilde_derivative_1;
        array_1d<double, 3> normal_tilde_derivative_2_term1, normal_tilde_derivative_2_term2, normal_tilde_derivative_2;

        MathUtils<double>::CrossProduct(normal_tilde_derivative_1_term1, base_vector1_derivative_11, rKinematicVariables.BaseVector2);
        MathUtils<double>::CrossProduct(normal_tilde_derivative_1_term2, rKinematicVariables.BaseVector1, base_vector1_derivative_12);
        normal_tilde_derivative_1 = normal_tilde_derivative_1_term1 + normal_tilde_derivative_1_term2;

        MathUtils<double>::CrossProduct(normal_tilde_derivative_2_term1, base_vector1_derivative_12, rKinematicVariables.BaseVector2);
        MathUtils<double>::CrossProduct(normal_tilde_derivative_2_term2, rKinematicVariables.BaseVector1, base_vector2_derivative_22);
        normal_tilde_derivative_2 = normal_tilde_derivative_2_term1 + normal_tilde_derivative_2_term2;

        if (rNormalVectorDerivatives.size1() != 3 || rNormalVectorDerivatives.size2() != 3)
            rNormalVectorDerivatives.resize(3, 3);

        for (IndexType j = 0; j < 3; j++)
        {
            rNormalVectorDerivatives(0, j) = normal_tilde_derivative_1[j] * inv_differential_area
                - rKinematicVariables.NormalVectorTilde[j] * inner_prod(rKinematicVariables.NormalVectorTilde, normal_tilde_derivative_1) * inv_differential_area_cube;
            rNormalVectorDerivatives(1, j) = normal_tilde_derivative_2[j] * inv_differential_area
                - rKinematicVariables.NormalVectorTilde[j] * inner_prod(rKinematicVariables.NormalVectorTilde, normal_tilde_derivative_2) * inv_differential_area_cube;
            rNormalVectorDerivatives(2, j) = 0.0;
        }
    }

    void CouplingNitsche6pCondition::CalculateTransformationFromLocalToGlobalCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTransformationMatrix) const
    {
        const double l_a1 = norm_2(rKinematicVariables.BaseVector1);
        array_1d<double, 3> local_base_1 = rKinematicVariables.BaseVector1 / l_a1;
        array_1d<double, 3> local_base_3 = rKinematicVariables.NormalVector;
        array_1d<double, 3> local_base_2 = ZeroVector(3);
        MathUtils<double>::CrossProduct(local_base_2, local_base_3, local_base_1);

        if (rTransformationMatrix.size1() != 6 || rTransformationMatrix.size2() != 6)
            rTransformationMatrix.resize(6, 6);
        noalias(rTransformationMatrix) = ZeroMatrix(6, 6);

        for (IndexType i = 0; i < 3; ++i)
        {
            IndexType j = (i + 1) % 3;

            rTransformationMatrix(i, 0) = local_base_1[i] * local_base_1[i];
            rTransformationMatrix(i, 1) = local_base_2[i] * local_base_2[i];
            rTransformationMatrix(i, 2) = local_base_3[i] * local_base_3[i];
            rTransformationMatrix(i, 3) = 2 * local_base_1[i] * local_base_2[i];
            rTransformationMatrix(i, 4) = 2 * local_base_2[i] * local_base_3[i];
            rTransformationMatrix(i, 5) = 2 * local_base_1[i] * local_base_3[i];

            rTransformationMatrix(i + 3, 0) = local_base_1[i] * local_base_1[j];
            rTransformationMatrix(i + 3, 1) = local_base_2[i] * local_base_2[j];
            rTransformationMatrix(i + 3, 2) = local_base_3[i] * local_base_3[j];
            rTransformationMatrix(i + 3, 3) = (local_base_1[i] * local_base_2[j]) + (local_base_2[i] * local_base_1[j]);
            rTransformationMatrix(i + 3, 4) = (local_base_2[i] * local_base_3[j]) + (local_base_3[i] * local_base_2[j]);
            rTransformationMatrix(i + 3, 5) = (local_base_1[i] * local_base_3[j]) + (local_base_3[i] * local_base_1[j]);
        }
    }

    void CouplingNitsche6pCondition::CalculateLocalConstitutiveMatrix(
        const Properties& rProperties,
        Matrix& rConstitutiveMatrixLocal) const
    {
        if (rConstitutiveMatrixLocal.size1() != 6 || rConstitutiveMatrixLocal.size2() != 6)
            rConstitutiveMatrixLocal.resize(6, 6);
        noalias(rConstitutiveMatrixLocal) = ZeroMatrix(6, 6);

        const double poisson_ratio = rProperties[POISSON_RATIO];
        const double youngs_modulus = rProperties[YOUNG_MODULUS];
        const double lame_lambda = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
        const double shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
        const double shear_correction_factor = 5.0 / 6.0;

        rConstitutiveMatrixLocal(0, 0) = lame_lambda;
        rConstitutiveMatrixLocal(0, 1) = lame_lambda * poisson_ratio;
        rConstitutiveMatrixLocal(1, 0) = lame_lambda * poisson_ratio;
        rConstitutiveMatrixLocal(1, 1) = lame_lambda;
        rConstitutiveMatrixLocal(3, 3) = lame_lambda * (1.0 - poisson_ratio) / 2.0;
        rConstitutiveMatrixLocal(4, 4) = shear_modulus * shear_correction_factor;
        rConstitutiveMatrixLocal(5, 5) = shear_modulus * shear_correction_factor;
    }

    void CouplingNitsche6pCondition::CalculateThicknessJacobian(
        double zeta,
        double Thickness,
        const KinematicVariables& rKinematicVariables,
        const Matrix& rNormalVectorDerivatives,
        Matrix& rJacobianInv,
        double& rJacobianDet) const
    {
        Matrix jacobian = ZeroMatrix(3, 3);
        for (IndexType i = 0; i < 3; ++i) {
            jacobian(0, i) = rKinematicVariables.BaseVector1[i] + (Thickness / 2.0) * zeta * rNormalVectorDerivatives(0, i);
            jacobian(1, i) = rKinematicVariables.BaseVector2[i] + (Thickness / 2.0) * zeta * rNormalVectorDerivatives(1, i);
            jacobian(2, i) = rKinematicVariables.NormalVector[i] * (Thickness / 2.0);
        }

        if (rJacobianInv.size1() != 3 || rJacobianInv.size2() != 3)
            rJacobianInv.resize(3, 3);

        MathUtils<double>::InvertMatrix(jacobian, rJacobianInv, rJacobianDet);
    }

    void CouplingNitsche6pCondition::CalculateLateralConormal(
        IndexType IntegrationPointIndex,
        const Matrix& rJacobianInv,
        double JacobianThicknessDet,
        const PatchType& rPatch,
        array_1d<double, 3>& rUnitConormal,
        double& rAreaScale) const
    {
        const IndexType geometry_part = (rPatch == PatchType::Master) ? 0 : 1;

        array_1d<double, 3> local_tangent;
        GetGeometry().GetGeometryPart(geometry_part).Calculate(LOCAL_TANGENT, local_tangent);

        const double tangent_norm = std::sqrt(local_tangent[0] * local_tangent[0] + local_tangent[1] * local_tangent[1]);

        // Unit parametric conormal
        array_1d<double, 3> nu_param = ZeroVector(3);
        nu_param[0] = local_tangent[1] / tangent_norm;
        nu_param[1] = -local_tangent[0] / tangent_norm;

        // Nanson's formula: n_phys * dA_phys = det(J) * J^-T * nu_param * dA_param
        array_1d<double, 3> n_phys = ZeroVector(3);
        for (IndexType k = 0; k < 3; ++k) {
            double value = 0.0;
            for (IndexType m = 0; m < 3; ++m) {
                value += rJacobianInv(m, k) * nu_param[m];
            }
            n_phys[k] = JacobianThicknessDet * value;
        }

        rAreaScale = norm_2(n_phys);
        noalias(rUnitConormal) = n_phys / rAreaScale;
    }

    void CouplingNitsche6pCondition::CalculateBAndDisplacementOperator(
        IndexType IntegrationPointIndex,
        double zeta,
        double Thickness,
        const Matrix& rJacobianInv,
        const Matrix& rNormalVectorDerivatives,
        const KinematicVariables& rKinematicVariables,
        const PatchType& rPatch,
        Matrix& rBOperator,
        Matrix& rDisplacementOperator) const
    {
        const IndexType geometry_part = (rPatch == PatchType::Master) ? 0 : 1;
        const auto& r_geometry = GetGeometry().GetGeometryPart(geometry_part);
        const IndexType number_of_control_points = r_geometry.size();
        const IndexType mat_size = number_of_control_points * 6;

        const auto& r_N = r_geometry.ShapeFunctionsValues();
        const auto& r_shape_functions_gradients = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
        const Matrix& r_DN_De = r_shape_functions_gradients[IntegrationPointIndex];

        KRATOS_ERROR_IF(r_DN_De.size1() != number_of_control_points || r_DN_De.size2() < 2)
            << "CouplingNitsche6pCondition: unexpected 1st-derivative shape function matrix size ("
            << r_DN_De.size1() << "x" << r_DN_De.size2() << "), expected (" << number_of_control_points << "x>=2)" << std::endl;
        KRATOS_ERROR_IF(rJacobianInv.size1() != 3 || rJacobianInv.size2() != 3)
            << "CouplingNitsche6pCondition: unexpected JacobianInv size (" << rJacobianInv.size1() << "x" << rJacobianInv.size2() << ")" << std::endl;

        Matrix shape_functions_derivatives_local = ZeroMatrix(number_of_control_points, 3);
        column(shape_functions_derivatives_local, 0) = column(r_DN_De, 0);
        column(shape_functions_derivatives_local, 1) = column(r_DN_De, 1);
        const Matrix shape_functions_derivatives_global = trans(prod(rJacobianInv, trans(shape_functions_derivatives_local)));

        const double normal_x = rKinematicVariables.NormalVector[0];
        const double normal_y = rKinematicVariables.NormalVector[1];
        const double normal_z = rKinematicVariables.NormalVector[2];

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

        if (rBOperator.size1() != 6 || rBOperator.size2() != mat_size)
            rBOperator.resize(6, mat_size);
        noalias(rBOperator) = ZeroMatrix(6, mat_size);

        if (rDisplacementOperator.size1() != 3 || rDisplacementOperator.size2() != mat_size)
            rDisplacementOperator.resize(3, mat_size);
        noalias(rDisplacementOperator) = ZeroMatrix(3, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const IndexType index = i * 6;
            const double Ni = r_N(IntegrationPointIndex, i);

            // Displacement DOFs (strain operator)
            rBOperator(0, index)     = shape_functions_derivatives_global(i, 0);
            rBOperator(1, index + 1) = shape_functions_derivatives_global(i, 1);
            rBOperator(2, index + 2) = shape_functions_derivatives_global(i, 2);
            rBOperator(3, index)     = shape_functions_derivatives_global(i, 1);
            rBOperator(3, index + 1) = shape_functions_derivatives_global(i, 0);
            rBOperator(4, index + 1) = shape_functions_derivatives_global(i, 2);
            rBOperator(4, index + 2) = shape_functions_derivatives_global(i, 1);
            rBOperator(5, index)     = shape_functions_derivatives_global(i, 2);
            rBOperator(5, index + 2) = shape_functions_derivatives_global(i, 0);

            // Rotation DOFs (strain operator) -- identical to Shell6pElement::CalculateBOperator
            rBOperator(0, index + 4) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_z) + (Ni * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (Thickness / 2.0);
            rBOperator(0, index + 5) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_y) + (Ni * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (Thickness / 2.0);

            rBOperator(1, index + 3) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_z) + (Ni * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (Thickness / 2.0);
            rBOperator(1, index + 5) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_x) + (Ni * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (Thickness / 2.0);

            rBOperator(2, index + 3) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_y) + (Ni * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (Thickness / 2.0);
            rBOperator(2, index + 4) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_x) + (Ni * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (Thickness / 2.0);

            rBOperator(3, index + 3) = -((shape_functions_derivatives_global(i, 0) * zeta * normal_z) + (Ni * (zeta * d_normal_z_dx + d_zeta_dx * normal_z))) * (Thickness / 2.0);
            rBOperator(3, index + 4) =  ((shape_functions_derivatives_global(i, 1) * zeta * normal_z) + (Ni * (zeta * d_normal_z_dy + d_zeta_dy * normal_z))) * (Thickness / 2.0);
            rBOperator(3, index + 5) = (((shape_functions_derivatives_global(i, 0) * zeta * normal_x) + (Ni * (zeta * d_normal_x_dx + d_zeta_dx * normal_x)))
                                       -((shape_functions_derivatives_global(i, 1) * zeta * normal_y) + (Ni * (zeta * d_normal_y_dy + d_zeta_dy * normal_y)))) * (Thickness / 2.0);

            rBOperator(4, index + 3) = (((shape_functions_derivatives_global(i, 1) * zeta * normal_y) + (Ni * (zeta * d_normal_y_dy + d_zeta_dy * normal_y)))
                                       -((shape_functions_derivatives_global(i, 2) * zeta * normal_z) + (Ni * (zeta * d_normal_z_dz + d_zeta_dz * normal_z)))) * (Thickness / 2.0);
            rBOperator(4, index + 4) = -((shape_functions_derivatives_global(i, 1) * zeta * normal_x) + (Ni * (zeta * d_normal_x_dy + d_zeta_dy * normal_x))) * (Thickness / 2.0);
            rBOperator(4, index + 5) =  ((shape_functions_derivatives_global(i, 2) * zeta * normal_x) + (Ni * (zeta * d_normal_x_dz + d_zeta_dz * normal_x))) * (Thickness / 2.0);

            rBOperator(5, index + 3) =  ((shape_functions_derivatives_global(i, 0) * zeta * normal_y) + (Ni * (zeta * d_normal_y_dx + d_zeta_dx * normal_y))) * (Thickness / 2.0);
            rBOperator(5, index + 4) = (((shape_functions_derivatives_global(i, 2) * zeta * normal_z) + (Ni * (zeta * d_normal_z_dz + d_zeta_dz * normal_z)))
                                       -((shape_functions_derivatives_global(i, 0) * zeta * normal_x) + (Ni * (zeta * d_normal_x_dx + d_zeta_dx * normal_x)))) * (Thickness / 2.0);
            rBOperator(5, index + 5) = -((shape_functions_derivatives_global(i, 2) * zeta * normal_y) + (Ni * (zeta * d_normal_y_dz + d_zeta_dz * normal_y))) * (Thickness / 2.0);

            // Displacement operator N_zeta: u(zeta) = u_mid + (t/2)*zeta*(theta x a3)
            rDisplacementOperator(0, index)     = Ni;
            rDisplacementOperator(1, index + 1) = Ni;
            rDisplacementOperator(2, index + 2) = Ni;

            const double coeff = Ni * zeta * (Thickness / 2.0);
            rDisplacementOperator(0, index + 4) =  coeff * normal_z;
            rDisplacementOperator(0, index + 5) = -coeff * normal_y;
            rDisplacementOperator(1, index + 3) = -coeff * normal_z;
            rDisplacementOperator(1, index + 5) =  coeff * normal_x;
            rDisplacementOperator(2, index + 3) =  coeff * normal_y;
            rDisplacementOperator(2, index + 4) = -coeff * normal_x;
        }
    }

    void CouplingNitsche6pCondition::CalculateAll(
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
        const SizeType number_of_nodes_total = number_of_nodes_master + number_of_nodes_slave;
        const SizeType mat_size_master = 6 * number_of_nodes_master;
        const SizeType mat_size = 6 * number_of_nodes_total;

        // Full local tangent matrix
        Matrix K = ZeroMatrix(mat_size, mat_size);

        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

        const Properties& r_props_master = GetProperties().GetSubProperties().front();
        const Properties& r_props_slave = GetProperties().GetSubProperties().back();
        const double thickness_master = r_props_master[THICKNESS];
        const double thickness_slave = r_props_slave[THICKNESS];

        Matrix D_local_master, D_local_slave;
        CalculateLocalConstitutiveMatrix(r_props_master, D_local_master);
        CalculateLocalConstitutiveMatrix(r_props_slave, D_local_slave);

        // 2-point Gauss quadrature through the thickness (matches Shell6pElement).
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

            // Reference-configuration surface arc-length Jacobian 
            array_1d<double, 3> local_tangent_master;
            GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent_master);
            const array_1d<double, 3> t_vec_master =
                local_tangent_master[0] * kinematics_master.BaseVector1 + local_tangent_master[1] * kinematics_master.BaseVector2;
            const double surface_jacobian = norm_2(t_vec_master);

            const double integration_weight = integration_points[point_number].Weight();

            // Master and slave patches are, in general, parametrized
            // independently and can have OPPOSITE director (a3) orientation
            // at the shared interface even when their tangent planes coincide
            const bool same_director_direction = inner_prod(kinematics_master.NormalVector, kinematics_slave.NormalVector) > 0.0;

            // Rotation-jump penalty 
            if (couple_rotation_x || couple_rotation_y || couple_rotation_z)
            {
                const Matrix& N_master_vals = r_geometry_master.ShapeFunctionsValues();
                const Matrix& N_slave_vals = r_geometry_slave.ShapeFunctionsValues();

                Matrix H_rot = ZeroMatrix(3, mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; ++i) {
                    const IndexType index = 6 * i;
                    if (couple_rotation_x) H_rot(0, index + 3) = N_master_vals(point_number, i);
                    if (couple_rotation_y) H_rot(1, index + 4) = N_master_vals(point_number, i);
                    if (couple_rotation_z) H_rot(2, index + 5) = N_master_vals(point_number, i);
                }
                for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
                    const IndexType index = 6 * (i + number_of_nodes_master);
                    if (couple_rotation_x) H_rot(0, index + 3) = -N_slave_vals(point_number, i);
                    if (couple_rotation_y) H_rot(1, index + 4) = -N_slave_vals(point_number, i);
                    if (couple_rotation_z) H_rot(2, index + 5) = -N_slave_vals(point_number, i);
                }

                noalias(K) += prod(trans(H_rot), H_rot) * (penalty_rotation * integration_weight * surface_jacobian);
            }

            // Through-thickness natural Nitsche flux (displacement + bending rotation).
            for (IndexType gauss_index = 0; gauss_index < 2; ++gauss_index)
            {
                const double zeta = gauss_zeta[gauss_index];
                // Same physical fiber, expressed in the slave's own (possibly
                // director-flipped) zeta convention -- see note above.
                const double zeta_slave = same_director_direction ? zeta : -zeta;

                Matrix jacobian_inv_master, jacobian_inv_slave;
                double jacobian_det_master, jacobian_det_slave;
                CalculateThicknessJacobian(zeta, thickness_master, kinematics_master, normal_derivatives_master, jacobian_inv_master, jacobian_det_master);
                CalculateThicknessJacobian(zeta_slave, thickness_slave, kinematics_slave, normal_derivatives_slave, jacobian_inv_slave, jacobian_det_slave);

                array_1d<double, 3> unit_conormal_master, unit_conormal_slave;
                double area_scale_master, area_scale_slave;
                CalculateLateralConormal(point_number, jacobian_inv_master, jacobian_det_master, PatchType::Master, unit_conormal_master, area_scale_master);
                CalculateLateralConormal(point_number, jacobian_inv_slave, jacobian_det_slave, PatchType::Slave, unit_conormal_slave, area_scale_slave);

                Matrix B_master, N_op_master;
                CalculateBAndDisplacementOperator(point_number, zeta, thickness_master, jacobian_inv_master, normal_derivatives_master, kinematics_master, PatchType::Master, B_master, N_op_master);
                Matrix B_slave, N_op_slave;
                CalculateBAndDisplacementOperator(point_number, zeta_slave, thickness_slave, jacobian_inv_slave, normal_derivatives_slave, kinematics_slave, PatchType::Slave, B_slave, N_op_slave);

                // Voigt traction matrices T_n (3x6), Voigt order [xx,yy,zz,xy,yz,xz].
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
                Matrix N_combined = ZeroMatrix(3, mat_size);
                for (SizeType c = 0; c < mat_size_master; ++c) {
                    column(F_combined, c) = column(F_master, c);
                    column(N_combined, c) = column(N_op_master, c);
                }
                for (SizeType c = 0; c < mat_size - mat_size_master; ++c) {
                    column(F_combined, c + mat_size_master) = -column(F_slave, c);
                    column(N_combined, c + mat_size_master) = -column(N_op_slave, c);
                }

                const double dW = integration_weight * surface_jacobian * area_scale_master * gauss_weight[gauss_index];

                const Matrix nitsche_block = prod(trans(F_combined), N_combined) + prod(trans(N_combined), F_combined);
                noalias(K) += nitsche_block * (-0.5 * dW);
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
            for (IndexType i = 0; i < number_of_nodes_master; ++i) {
                const array_1d<double, 3>& disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
                const array_1d<double, 3>& rot = r_geometry_master[i].FastGetSolutionStepValue(ROTATION);
                const IndexType index = 6 * i;
                d_current[index] = disp[0]; d_current[index + 1] = disp[1]; d_current[index + 2] = disp[2];
                d_current[index + 3] = rot[0]; d_current[index + 4] = rot[1]; d_current[index + 5] = rot[2];
            }
            for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
                const array_1d<double, 3>& disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
                const array_1d<double, 3>& rot = r_geometry_slave[i].FastGetSolutionStepValue(ROTATION);
                const IndexType index = 6 * (i + number_of_nodes_master);
                d_current[index] = disp[0]; d_current[index + 1] = disp[1]; d_current[index + 2] = disp[2];
                d_current[index + 3] = rot[0]; d_current[index + 4] = rot[1]; d_current[index + 5] = rot[2];
            }

            noalias(rRightHandSideVector) = -prod(K, d_current);
        }

        KRATOS_CATCH("")
    }

    int CouplingNitsche6pCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(NITSCHE_STABILIZATION_FACTOR))
            << "No NITSCHE_STABILIZATION_FACTOR defined in property of CouplingNitsche6pCondition" << std::endl;

        KRATOS_ERROR_IF(GetProperties().GetSubProperties().size() < 1)
            << "CouplingNitsche6pCondition needs master/slave material sub-properties (YOUNG_MODULUS, "
            << "POISSON_RATIO, THICKNESS) attached to its Properties" << std::endl;
        return 0;
    }

    void CouplingNitsche6pCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const IndexType number_of_nodes_master = r_geometry_master.size();
        const IndexType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 6 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(6 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 6;
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 6 * (i + number_of_nodes_master);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = r_node.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = r_node.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = r_node.GetDof(ROTATION_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CouplingNitsche6pCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const IndexType number_of_nodes_master = r_geometry_master.size();
        const IndexType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(6 * (number_of_nodes_master + number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_X));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Y));
            rElementalDofList.push_back(r_node.pGetDof(ROTATION_Z));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
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
