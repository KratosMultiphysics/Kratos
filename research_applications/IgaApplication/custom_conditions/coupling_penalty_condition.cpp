//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_conditions/coupling_penalty_condition.h"

namespace Kratos
{
    void CouplingPenaltyCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        // Size definitions
        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size = 3 * (number_of_nodes_master + number_of_nodes_slave);

        // Memory allocation
        if (CalculateStiffnessMatrixFlag) {
            if (rLeftHandSideMatrix.size1() != mat_size) {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        }
        if (CalculateResidualVectorFlag) {
            if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
        }

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        const bool integrate_conservative = GetProperties().Has(INTEGRATE_CONSERVATIVE)
            ? GetProperties()[INTEGRATE_CONSERVATIVE]
            : false;
        if (integrate_conservative) {
            DeterminantOfJacobianInitial(r_geometry_master, determinant_jacobian_vector);
        }
        else {
            r_geometry_master.DeterminantOfJacobian(determinant_jacobian_vector);
        }

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues();
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues();

            //FOR DISPLACEMENTS
            Matrix H = ZeroMatrix(3, mat_size);
            for (IndexType i = 0; i < number_of_nodes_master; ++i)
            {
                IndexType index = 3 * i;
                if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = N_master(point_number, i);
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = N_master(point_number, i);
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = N_master(point_number, i);
            }

            for (IndexType i = 0; i < number_of_nodes_slave; ++i)
            {
                IndexType index = 3 * (i + number_of_nodes_master);
                if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = -N_slave(point_number, i);
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = -N_slave(point_number, i);
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = -N_slave(point_number, i);
            }

            // Differential area
            const double penalty_integration = penalty * integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

            // Rotation coupling
            if (Is(IgaFlags::FIX_ROTATION_X))
            {
                Vector phi_r = ZeroVector(mat_size);
                Matrix phi_rs = ZeroMatrix(mat_size, mat_size);
                array_1d<double, 2> diff_phi;

                CalculateRotationalShapeFunctions(point_number, phi_r, phi_rs, diff_phi);

                if (CalculateStiffnessMatrixFlag) {
                    for (IndexType i = 0; i < mat_size; ++i)
                    {
                        for (IndexType j = 0; j < mat_size; ++j)
                        {
                            rLeftHandSideMatrix(i, j) = (phi_r(i) * phi_r(j) + diff_phi(0) * phi_rs(i, j)) * penalty_integration;
                        }
                    }
                }

                if (CalculateResidualVectorFlag) {
                    for (IndexType i = 0; i < mat_size; ++i)
                    {
                        rRightHandSideVector[i] = (diff_phi(0) * phi_r(i)) * penalty_integration;
                    }
                }
            }

            // Assembly
            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += prod(trans(H), H) * penalty_integration;
            }
            if (CalculateResidualVectorFlag) {

                Vector u(mat_size);
                for (IndexType i = 0; i < number_of_nodes_master; ++i)
                {
                    const array_1d<double, 3> disp = r_geometry_master[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * i;
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }
                for (IndexType i = 0; i < number_of_nodes_slave; ++i)
                {
                    const array_1d<double, 3> disp = r_geometry_slave[i].FastGetSolutionStepValue(DISPLACEMENT);
                    IndexType index = 3 * (i + number_of_nodes_master);
                    u[index]     = disp[0];
                    u[index + 1] = disp[1];
                    u[index + 2] = disp[2];
                }

                noalias(rRightHandSideVector) -= prod(prod(trans(H), H), u) * penalty_integration;
            }
        }

        KRATOS_CATCH("")
    }

    void CouplingPenaltyCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType nb_nodes = rGeometry.PointsNumber();

        Matrix J = ZeroMatrix(working_space_dimension, local_space_dimension);
        for (IndexType pnt = 0; pnt < nb_integration_points; pnt++)
        {
            const Matrix& r_DN_De = rGeometry.ShapeFunctionsLocalGradients()[pnt];
            J.clear();
            for (IndexType i = 0; i < nb_nodes; ++i) {
                const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        J(k, m) += value * r_DN_De(i, m);
                    }
                }
            }

            //Compute the tangent and  the normal to the boundary vector
            array_1d<double, 3> local_tangent;
            GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent);

            array_1d<double, 3> a_1 = column(J, 0);
            array_1d<double, 3> a_2 = column(J, 1);

            rDeterminantOfJacobian[pnt] = norm_2(a_1 * local_tangent[0] + a_2 * local_tangent[1]);
        }
    }

    void CouplingPenaltyCondition::CalculateRotationalShapeFunctions(
        IndexType IntegrationPointIndex,
        Vector &phi_r, 
        Matrix &phi_rs, 
        array_1d<double, 2> &diff_phi)
    {
        // compute rotation (master)
        array_1d<double, 3> local_tangent_master;
        GetGeometry().GetGeometryPart(0).Calculate(LOCAL_TANGENT, local_tangent_master);

        const IntegrationMethod integration_method_master = GetGeometry().GetGeometryPart(0).GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_master = GetGeometry().GetGeometryPart(0).ShapeFunctionsLocalGradients(integration_method_master);
        const Matrix& shape_functions_gradients_master = r_shape_functions_gradients_master(IntegrationPointIndex);

        const SizeType number_of_nodes_master = GetGeometry().GetGeometryPart(0).size();

        Vector phi_r_master = ZeroVector(number_of_nodes_master * 3);
        Matrix phi_rs_master = ZeroMatrix(number_of_nodes_master * 3, number_of_nodes_master * 3);
        array_1d<double, 2> phi_master;
        array_1d<double, 3> trim_tangents_master;

        CalculateRotation(IntegrationPointIndex, shape_functions_gradients_master, phi_r_master, phi_rs_master, phi_master, trim_tangents_master, local_tangent_master, true);

        // compute rotation (slave)
        array_1d<double, 3> local_tangent_slave;
        GetGeometry().GetGeometryPart(1).Calculate(LOCAL_TANGENT, local_tangent_slave);

        const IntegrationMethod integration_method_slave = GetGeometry().GetGeometryPart(1).GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_slave = GetGeometry().GetGeometryPart(1).ShapeFunctionsLocalGradients(integration_method_slave);
        const Matrix& shape_functions_gradients_slave = r_shape_functions_gradients_slave(IntegrationPointIndex);

        const SizeType number_of_nodes_slave = GetGeometry().GetGeometryPart(1).size();

        Vector phi_r_slave = ZeroVector(number_of_nodes_slave * 3);
        Matrix phi_rs_slave = ZeroMatrix(number_of_nodes_slave * 3, number_of_nodes_slave * 3);
        array_1d<double, 2> phi_slave;
        array_1d<double, 3> trim_tangents_slave;

        CalculateRotation(IntegrationPointIndex, shape_functions_gradients_slave, phi_r_slave, phi_rs_slave, phi_slave, trim_tangents_slave, local_tangent_slave, false);

        // compute phi_r, phi_rs and diff_phi
        bool opposite_direction_of_trims = true;
        if (inner_prod(trim_tangents_master, trim_tangents_slave) > 0) // tangents have the same direction (assumption for coordinates system changes)
        {
            opposite_direction_of_trims = false;
        }

        if (opposite_direction_of_trims)
        {
            diff_phi = phi_slave + phi_master;
        }
        else
        {
            diff_phi = phi_slave - phi_master;
        }
        
        for (IndexType i = 0; i < phi_r_master.size(); i++)
        {
            phi_r(i) = phi_r_master(i);
        }

        SizeType index = phi_r_master.size();
        for (IndexType i = 0; i < phi_r_slave.size(); i++)
        {
            if (opposite_direction_of_trims)
            {
                phi_r(i + index) = phi_r_slave(i);
            }
            else
            {
                phi_r(i + index) = -phi_r_slave(i);
            }
        }

        for (IndexType i = 0; i < phi_rs_master.size1(); i++)
        {
            for (IndexType j = 0; j < phi_rs_master.size2(); j++)
            {
                phi_rs(i, j) = phi_rs_master(i, j);
            }
        }

        SizeType index_1 = phi_rs_master.size1();
        SizeType index_2 = phi_rs_master.size2();
        for (IndexType i = 0; i < phi_rs_slave.size1(); i++)
        {
            for (IndexType j = 0; j < phi_rs_slave.size2(); j++)
            {
                if (opposite_direction_of_trims)
                    phi_rs(i + index_1, j + index_2) = phi_rs_slave(i, j);
                else
                    phi_rs(i + index_1, j + index_2) = -phi_rs_slave(i, j);
            }
        }
    } 

    void CouplingPenaltyCondition::CalculateRotation(
        IndexType IntegrationPointIndex,
        const Matrix &rShapeFunctionGradientValues,
        Vector &phi_r,
        Matrix &phi_rs,
        array_1d<double, 2> &phi,
        array_1d<double, 3> &trim_tangent,
        const Vector &local_tangent,
        const bool master)
    {
        KRATOS_TRY

        const SizeType number_of_points = rShapeFunctionGradientValues.size1();
        
        // compute the initialize base vectors of master or slave 
        Vector g10 = ZeroVector(3);
        Vector g20 = ZeroVector(3);
        Vector g30 = ZeroVector(3);
        if (master)
        {
            for (SizeType i = 0; i < GetGeometry().GetGeometryPart(0).size(); ++i){
                g10[0] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 0);
                g10[1] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 0);
                g10[2] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 0);

                g20[0] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 1);
                g20[1] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 1);
                g20[2] += (GetGeometry().GetGeometryPart(0).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 1);

                MathUtils<double>::CrossProduct(g30, g10, g20);
                g30 = g30 / norm_2(g30);
            }
        }
        else
        {
            for (SizeType i = 0; i < GetGeometry().GetGeometryPart(1).size(); ++i){
                g10[0] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 0);
                g10[1] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 0);
                g10[2] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 0);

                g20[0] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 1);
                g20[1] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 1);
                g20[2] += (GetGeometry().GetGeometryPart(1).GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 1);

                MathUtils<double>::CrossProduct(g30, g10, g20);
                g30 = g30 / norm_2(g30);
            }
        }

        // compute the actual base vectors of master or slave
        array_1d<double, 3> g1, g2, g3;
        Matrix J;

        if (master)
        {
            GetGeometry().GetGeometryPart(0).Jacobian(J, IntegrationPointIndex);
        }
        else 
        {
            GetGeometry().GetGeometryPart(1).Jacobian(J, IntegrationPointIndex);
        }

        g1 = column(J, 0);
        g2 = column(J, 1);

        MathUtils<double>::CrossProduct(g3, g1, g2);
        g3 = g3 / norm_2(g3);

        // compute the tangent (T2) and the normal (T1) to the boundary vector
        array_1d<double, 3> T1, T2;
        T2 = local_tangent[0] * g10 + local_tangent[1] * g20;
        trim_tangent = T2;
        MathUtils<double>::CrossProduct(T1, T2, g30);
        T2 = T2 / norm_2(T2);
        T1 = T1 / norm_2(T1);

        // compute the a3 displacement
        array_1d<double, 3> w = g3 - g30;
        array_1d<double, 3> sinus_omega_vector;
        MathUtils<double>::CrossProduct(sinus_omega_vector, g30, w);

        array_1d<double, 2> sinus_omega;
        sinus_omega(0) = inner_prod(sinus_omega_vector, T2);
        sinus_omega(1) = inner_prod(sinus_omega_vector, T1);

        array_1d<double, 3> omega;
        if (sinus_omega(0) > 1.0)
            sinus_omega(0) = 0.999999;
        if (sinus_omega(1) > 1.0)
            sinus_omega(1) = 0.999999;
        omega(0) = asin(sinus_omega(0));
        omega(1) = asin(sinus_omega(1));

        phi(0) = omega(0);
        phi(1) = omega(1);

        // compute variation of the a3 
        array_1d<double, 3> t3 = g3;
        array_1d<double, 3> tilde_t3; 
        MathUtils<double>::CrossProduct(tilde_t3, g1, g2);
        double length_t3 = norm_2(tilde_t3);

        std::vector<array_1d<double, 3>> t3_r(number_of_points * 3);
        std::vector<array_1d<double, 3>> tilde_3_r(number_of_points * 3);
        Vector line_t3_r = ZeroVector(number_of_points * 3);
        std::vector<array_1d<double, 3>> sinus_omega_r(number_of_points * 3);

        for (IndexType n = 0; n < number_of_points; n++)
        {
            for (IndexType i = 0; i < 3; i++)
            {
                int nb_dof = n * 3 + i;

                //variations of the basis vectors
                array_1d<double, 3> a1_r = ZeroVector(3);
                array_1d<double, 3> a2_r = ZeroVector(3);

                a1_r(i) = rShapeFunctionGradientValues(n, 0);
                a2_r(i) = rShapeFunctionGradientValues(n, 1);
                
                array_1d<double, 3> a1_r__g2, g1__a2_r = ZeroVector(3);
                MathUtils<double>::CrossProduct(a1_r__g2, a1_r, g2);
                MathUtils<double>::CrossProduct(g1__a2_r, g1, a2_r);

                //variation of the non normalized local vector
                tilde_3_r[nb_dof] = a1_r__g2 + g1__a2_r;
                line_t3_r[nb_dof] = inner_prod(t3, tilde_3_r[nb_dof]);
                t3_r[nb_dof] = tilde_3_r[nb_dof] / length_t3 - line_t3_r[nb_dof] * t3 / length_t3;

                MathUtils<double>::CrossProduct(sinus_omega_r[nb_dof], g30, t3_r[nb_dof]);
                phi_r(nb_dof) = 1.0 / sqrt(1.0 - pow(sinus_omega(0), 2))*inner_prod(sinus_omega_r[nb_dof], T2);
            }
        }

        for (IndexType n = 0; n < number_of_points; n++)
        {
            for (IndexType i = 0; i < 3; i++)
            {
                int nb_dof_n = n * 3 + i;
                
                //variations of the basis vectors
                array_1d<double, 3> a1_r_n = ZeroVector(3);
                array_1d<double, 3> a2_r_n = ZeroVector(3);

                a1_r_n(i) = rShapeFunctionGradientValues(n, 0);
                a2_r_n(i) = rShapeFunctionGradientValues(n, 1);

                for (IndexType m = 0; m < number_of_points; m++)
                {
                    for (IndexType j = 0; j < 3; j++)
                    {
                        int nb_dof_m = m * 3 + j;

                        //variations of the basis vectors
                        array_1d<double, 3> a1_r_m = ZeroVector(3);
                        array_1d<double, 3> a2_r_m = ZeroVector(3);

                        a1_r_m(j) = rShapeFunctionGradientValues(m, 0);
                        a2_r_m(j) = rShapeFunctionGradientValues(m, 1);

                        //variation of the non normalized local vector
                        array_1d<double, 3> a1_r_n__a2_r_m, a1_r_m__a2_r_n = ZeroVector(3);
                        MathUtils<double>::CrossProduct(a1_r_n__a2_r_m, a1_r_n, a2_r_m);
                        MathUtils<double>::CrossProduct(a1_r_m__a2_r_n, a1_r_m, a2_r_n);

                        array_1d<double, 3> tilde_t3_rs = a1_r_n__a2_r_m + a1_r_m__a2_r_n;
                        double line_t3_rs = inner_prod(t3_r[nb_dof_m], tilde_3_r[nb_dof_n]) + inner_prod(t3, tilde_t3_rs);

                        array_1d<double, 3> t3_rs = (tilde_t3_rs*length_t3 - line_t3_r[nb_dof_m] * tilde_3_r[nb_dof_n]) / pow(length_t3, 2)
                            - line_t3_rs * t3 / length_t3 - line_t3_r[nb_dof_n] * (t3_r[nb_dof_m] * length_t3 - line_t3_r[nb_dof_m] * t3) / pow(length_t3, 2);

                        array_1d<double, 3> sinus_omega_rs = ZeroVector(3);
                        MathUtils<double>::CrossProduct(sinus_omega_rs, g30, t3_rs);

                        phi_rs(n * 3 + i, m * 3 + j) = inner_prod(sinus_omega_rs, T2) / sqrt(1.0 - pow(sinus_omega(0), 2))
                            + inner_prod(sinus_omega_r[nb_dof_m], T2)*inner_prod(sinus_omega_r[nb_dof_n], T2)*sinus_omega(0) / pow(1.0
                                - pow(sinus_omega(0), 2), 1.5);
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }     

    int CouplingPenaltyCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyCondition" << std::endl;
        return 0;
    }

    void CouplingPenaltyCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto& r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto& r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        if (rResult.size() != 3 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(3 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry_master[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = 3 * (i + number_of_nodes_master);
            const auto& r_node = r_geometry_slave[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CouplingPenaltyCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const auto r_geometry_master = GetGeometry().GetGeometryPart(0);
        const auto r_geometry_slave = GetGeometry().GetGeometryPart(1);

        const SizeType number_of_nodes_master = r_geometry_master.size();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * (number_of_nodes_master + number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave.GetPoint(i);
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


