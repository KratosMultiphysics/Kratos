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
#include "custom_conditions/support_lagrange_condition.h"


namespace Kratos
{
    void SupportLagrangeCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (Is(IgaFlags::FIX_ROTATION_X))
        {
            const auto& r_geometry = GetGeometry();

            // Size definitions
            const SizeType number_of_nodes = r_geometry.size();

            // Integration
            const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

            // Prepare memory
            if (m_phi_r_vector.size() != r_number_of_integration_points)
                m_phi_r_vector.resize(r_number_of_integration_points);

            for (IndexType point_number = 0; point_number < r_number_of_integration_points; point_number++)
            {
                Vector phi_r = ZeroVector(3 * number_of_nodes);
                Matrix phi_rs = ZeroMatrix(3 * number_of_nodes, 3 * number_of_nodes);
                array_1d<double, 2> diff_phi;

                CalculateRotationalShapeFunctions(point_number, phi_r, phi_rs, diff_phi);

                m_phi_r_vector[point_number] = phi_r;
            }
        }

        KRATOS_CATCH("")
    }

    void SupportLagrangeCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();

        // Integration
        const typename GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // Determinant of jacobian 
        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        const bool integrate_conservative = GetProperties().Has(INTEGRATE_CONSERVATIVE)
            ? GetProperties()[INTEGRATE_CONSERVATIVE]
            : false;
        if (integrate_conservative) {
            DeterminantOfJacobianInitial(r_geometry, determinant_jacobian_vector);
        }
        else {
            r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);
        }

        if (!Is(IgaFlags::FIX_ROTATION_X))
        {
            const SizeType mat_size = number_of_non_zero_nodes * 6;

            Matrix LHS = ZeroMatrix(mat_size, mat_size);
            Vector u(mat_size);

            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            // non zero node counter for Lagrange Multipliers.
            IndexType counter_n = 0;
            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                // Differential area, being 1 for points.
                const double integration = (r_geometry.GetGeometryParent(0).GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Point)
                    ? 1
                    : integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

                // loop over Lagrange Multipliers
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    // non zero node counter for for displacements.
                    IndexType counter_m = 0;
                    if (r_N(point_number, i) > shape_function_tolerance) {
                        // loop over shape functions of displacements
                        for (IndexType j = 0; j < number_of_nodes; j++) {
                            if (r_N(point_number, j) > shape_function_tolerance) {
                                const double NN = r_N(point_number, j) * r_N(point_number, i) * integration;

                                // indices in local stiffness matrix.
                                const IndexType ibase = counter_n * 3 + 3 * (number_of_non_zero_nodes);
                                const IndexType jbase = counter_m * 3;

                                // Matrix in following shape:
                                // |0 H^T|
                                // |H 0  |

                                LHS(ibase, jbase) = NN;
                                LHS(ibase + 1, jbase + 1) = NN;
                                LHS(ibase + 2, jbase + 2) = NN;

                                LHS(jbase, ibase) = NN;
                                LHS(jbase + 1, ibase + 1) = NN;
                                LHS(jbase + 2, ibase + 2) = NN;

                                counter_m++;
                            }
                        }
                        counter_n++;
                    }
                }
                if (CalculateStiffnessMatrixFlag) {
                    noalias(rLeftHandSideMatrix) += LHS;
                }

                if (CalculateResidualVectorFlag) {
                    const array_1d<double, 3>& displacement = (Has(DISPLACEMENT))
                        ? this->GetValue(DISPLACEMENT)
                        : ZeroVector(3);

                    IndexType counter = 0;
                    for (IndexType i = 0; i < number_of_nodes; i++) {
                        for (IndexType n = 0; n < r_N.size1(); ++n) {
                            if (r_N(n, i) > shape_function_tolerance) {
                                const array_1d<double, 3>& r_disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

                                IndexType index = 3 * counter;
                                u[index] = (r_disp[0] - displacement[0]);
                                u[index + 1] = (r_disp[1] - displacement[1]);
                                u[index + 2] = (r_disp[2] - displacement[2]);
                                counter++;
                            }
                        }
                    }
                    for (IndexType i = 0; i < number_of_nodes; i++) {
                        for (IndexType n = 0; n < r_N.size1(); ++n) {
                            if (r_N(n, i) > shape_function_tolerance) {
                                const array_1d<double, 3>& r_l_m = r_geometry[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                                IndexType index = 3 * (counter);
                                u[index] = r_l_m[0];
                                u[index + 1] = r_l_m[1];
                                u[index + 2] = r_l_m[2];
                                counter++;
                            }
                        }
                    }

                    noalias(rRightHandSideVector) -= prod(LHS, u);
                }
            }
        }
        else if (Is(IgaFlags::FIX_ROTATION_X)) // Rotation coupling
        {
            const SizeType number_of_non_zero_DOFs = GetNumberOfNonZeroDOFsRotation();
            const SizeType mat_size = 3*number_of_nodes + 3*number_of_non_zero_nodes + number_of_non_zero_DOFs;

            Matrix LHS = ZeroMatrix(mat_size, mat_size);
            Vector u(mat_size);

            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            // non zero node counter for Lagrange Multipliers.
            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                // Differential area, being 1 for points.
                const double integration = (r_geometry.GetGeometryParent(0).GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Point)
                    ? 1
                    : integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

                // loop over Lagrange Multipliers
                IndexType counter_n = 0;
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    // non zero node counter for for displacements.
                    IndexType counter_m = 0;
                    if (r_N(point_number, i) > shape_function_tolerance) {
                        // loop over shape functions of displacements
                        for (IndexType j = 0; j < number_of_nodes; j++) {
                                const double NN = r_N(point_number, j) * r_N(point_number, i) * integration;

                                // indices in local stiffness matrix.
                                const IndexType ibase = counter_n * 3 + 3 * (number_of_non_zero_nodes);
                                const IndexType jbase = counter_m * 3;

                                // Matrix in following shape:
                                // |0 H^T|
                                // |H 0  |

                                LHS(ibase, jbase) = NN;
                                LHS(ibase + 1, jbase + 1) = NN;
                                LHS(ibase + 2, jbase + 2) = NN;

                                LHS(jbase, ibase) = NN;
                                LHS(jbase + 1, ibase + 1) = NN;
                                LHS(jbase + 2, ibase + 2) = NN;

                                counter_m++;
                        }
                        counter_n++;
                    }
                }

                // Rotation
                Vector H_Mu;
                H_Mu.resize(number_of_non_zero_DOFs);
                H_Mu = ZeroVector(number_of_non_zero_DOFs);

                IndexType index_mu = 0;
                for (IndexType i = 0; i < number_of_nodes; ++i)
                {
                    if (r_N(point_number, i) > shape_function_tolerance) {
                        H_Mu[index_mu] = r_N(point_number, i);
                        index_mu++;
                    }
                }

                Vector phi_r = ZeroVector(3 * number_of_nodes);
                phi_r = m_phi_r_vector[point_number];

                // loop over Lagrange Multipliers
                for (IndexType i = 0; i < number_of_non_zero_DOFs; ++i)
                {
                    // loop over shape functions of rotations
                    for (IndexType j = 0; j < 3 * number_of_nodes; ++j) 
                    {
                        const double Hphi = phi_r[j] * H_Mu[i] * integration;
                        const IndexType ibase = i + 3*number_of_nodes + 3*number_of_non_zero_nodes;
                            
                        // Matrix in following shape:
                        // |0 H^T M^T|
                        // |H 0   0  |
                        // |M 0   0  |

                        LHS(ibase,     j)     = Hphi;
                        LHS(j,     ibase)     = Hphi;
                    }
                }

                if (CalculateStiffnessMatrixFlag) {
                    noalias(rLeftHandSideMatrix) += LHS;
                }

                if (CalculateResidualVectorFlag) {
                    const array_1d<double, 3>& displacement = (Has(DISPLACEMENT))
                        ? this->GetValue(DISPLACEMENT)
                        : ZeroVector(3);

                    IndexType counter = 0;
                    for (IndexType i = 0; i < number_of_nodes; i++) {
                        const array_1d<double, 3>& r_disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

                        IndexType index = 3 * counter;
                        u[index] = (r_disp[0] - displacement[0]);
                        u[index + 1] = (r_disp[1] - displacement[1]);
                        u[index + 2] = (r_disp[2] - displacement[2]);
                        counter++;
                    }
                    for (IndexType i = 0; i < number_of_nodes; i++) {
                        for (IndexType n = 0; n < r_N.size1(); ++n) {
                            if (r_N(n, i) > shape_function_tolerance) {
                                const array_1d<double, 3>& r_l_m = r_geometry[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                                IndexType index = 3 * (counter);
                                u[index] = r_l_m[0];
                                u[index + 1] = r_l_m[1];
                                u[index + 2] = r_l_m[2];
                                counter++;
                            }
                        }
                    }

                    // rotation
                    IndexType index_rot = 3 * (counter);
                    for (IndexType m = 0; m < phi_r.size(); ++m) {
                        IndexType a = m / 3;
                        IndexType b = m % 3;

                        if ((std::abs(phi_r(m)) > shape_function_tolerance) && (r_N(0, a) > shape_function_tolerance)) {
                            if(b == 0){
                                u[index_rot] = r_geometry[a].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_X);
                                index_rot++;
                            }
                            else if(b == 1){
                                u[index_rot] = r_geometry[a].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_Y);
                                index_rot++;
                            }
                            else if(b == 2){
                                u[index_rot] = r_geometry[a].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_Z);
                                index_rot++;
                            }
                        }
                    }
                    noalias(rRightHandSideVector) -= prod(LHS, u);
                }
            }
        }
        KRATOS_CATCH("")
    }

    void SupportLagrangeCondition::DeterminantOfJacobianInitial(
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

    void SupportLagrangeCondition::CalculateRotationalShapeFunctions(
        IndexType IntegrationPointIndex,
        Vector &phi_r, 
        Matrix &phi_rs, 
        array_1d<double, 2> &diff_phi)
    {
        // compute rotation (support)
        array_1d<double, 3> local_tangent_support;
        GetGeometry().Calculate(LOCAL_TANGENT, local_tangent_support);

        const IntegrationMethod integration_method_support = GetGeometry().GetDefaultIntegrationMethod();
        const GeometryType::ShapeFunctionsGradientsType& r_shape_functions_gradients_support = GetGeometry().ShapeFunctionsLocalGradients(integration_method_support);
        const Matrix& shape_functions_gradients_support = r_shape_functions_gradients_support(IntegrationPointIndex);

        const SizeType number_of_nodes_support = GetGeometry().size();

        Vector phi_r_support = ZeroVector(number_of_nodes_support * 3);
        Matrix phi_rs_support = ZeroMatrix(number_of_nodes_support * 3, number_of_nodes_support * 3);
        array_1d<double, 2> phi_support;
        array_1d<double, 3> trim_tangents_support;

        CalculateRotation(IntegrationPointIndex, shape_functions_gradients_support, phi_r_support, phi_rs_support, phi_support, trim_tangents_support, local_tangent_support);

        //TO DO
        diff_phi = phi_support;
        
        for (IndexType i = 0; i < phi_r_support.size(); i++)
        {
            phi_r(i) = phi_r_support(i);
        }

        for (IndexType i = 0; i < phi_rs_support.size1(); i++)
        {
            for (IndexType j = 0; j < phi_rs_support.size2(); j++)
            {
                phi_rs(i, j) = phi_rs_support(i, j);
            }
        }
    } 

    void SupportLagrangeCondition::CalculateRotation(
        IndexType IntegrationPointIndex,
        const Matrix &rShapeFunctionGradientValues,
        Vector &phi_r,
        Matrix &phi_rs,
        array_1d<double, 2> &phi,
        array_1d<double, 3> &trim_tangent,
        const Vector &local_tangent)
    {
        KRATOS_TRY

        const SizeType number_of_points = rShapeFunctionGradientValues.size1();
        
        // compute the initialize base vectors of master or slave 
        Vector g10 = ZeroVector(3);
        Vector g20 = ZeroVector(3);
        Vector g30 = ZeroVector(3);

        for (SizeType i = 0; i < GetGeometry().size(); ++i){
            g10[0] += (GetGeometry().GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 0);
            g10[1] += (GetGeometry().GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 0);
            g10[2] += (GetGeometry().GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 0);

            g20[0] += (GetGeometry().GetPoint( i ).X0()) * rShapeFunctionGradientValues(i, 1);
            g20[1] += (GetGeometry().GetPoint( i ).Y0()) * rShapeFunctionGradientValues(i, 1);
            g20[2] += (GetGeometry().GetPoint( i ).Z0()) * rShapeFunctionGradientValues(i, 1);

            MathUtils<double>::CrossProduct(g30, g10, g20);
            g30 = g30 / norm_2(g30);
        }
    
        // compute the actual base vectors of master or slave
        array_1d<double, 3> g1, g2, g3;
        Matrix J;

        GetGeometry().Jacobian(J, IntegrationPointIndex);

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

        // KRATOS_WATCH(T2)
        // KRATOS_WATCH(g3)
        // KRATOS_WATCH(g30)

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

    void SupportLagrangeCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();

        if (!Is(IgaFlags::FIX_ROTATION_X))
        {
            if (rResult.size() != 6 * number_of_non_zero_nodes)
                rResult.resize(6 * number_of_non_zero_nodes, false);

            IndexType counter = 0;
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType n = 0; n < r_N.size1(); ++n) {
                    if (r_N(n, i) > shape_function_tolerance) {
                        const IndexType index = counter * 3;
                        const auto& r_node = r_geometry[i];
                        rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                        rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
                        counter++;
                    }
                }
            }

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType n = 0; n < r_N.size1(); ++n) {
                    if (r_N(n, i) > shape_function_tolerance) {
                        const IndexType index = 3 * (counter);
                        const auto& r_node = r_geometry[i];
                        rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
                        rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
                        rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
                        counter++;
                    }
                }
            }
        }
        else
        {
            if (rResult.size() != 3*number_of_nodes + 3*number_of_non_zero_nodes + GetNumberOfNonZeroDOFsRotation())
                rResult.resize(3*number_of_nodes + 3*number_of_non_zero_nodes + GetNumberOfNonZeroDOFsRotation(), false);

            IndexType counter = 0;
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const IndexType index = counter * 3;
                const auto& r_node = r_geometry[i];
                rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
                rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
                counter++;
            }

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType n = 0; n < r_N.size1(); ++n) {
                    if (r_N(n, i) > shape_function_tolerance) {
                        const IndexType index = 3 * (counter);
                        const auto& r_node = r_geometry[i];
                        rResult[index]     = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
                        rResult[index + 1] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
                        rResult[index + 2] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
                        counter++;
                    }
                }
            }

            const Vector phi_r = m_phi_r_vector[0];
            IndexType index_rot = 3 * (counter);
            for (IndexType m = 0; m < phi_r.size(); ++m) {
                IndexType a = m / 3;
                IndexType b = m % 3;

                if ((std::abs(phi_r(m)) > shape_function_tolerance) && (r_N(0, a) > shape_function_tolerance)) {
                    const auto& r_node = r_geometry[a];
                    if(b == 0){
                        rResult[index_rot] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_X).EquationId();
                        index_rot++;
                    }
                    else if(b == 1){
                        rResult[index_rot] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_Y).EquationId();
                        index_rot++;
                    }
                    else if(b == 2){
                        rResult[index_rot] = r_node.GetDof(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_Z).EquationId();
                        index_rot++;
                    }
                }
            }
        }
    }

    void SupportLagrangeCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        if (!Is(IgaFlags::FIX_ROTATION_X))
        {
            rElementalDofList.reserve(6 * GetNumberOfNonZeroNodes());

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType n = 0; n < r_N.size1(); ++n) {
                    if (r_N(n, i) > shape_function_tolerance) {
                        const auto& r_node = r_geometry[i];
                        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
                        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
                        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
                    }
                }
            }

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType n = 0; n < r_N.size1(); ++n) {
                    if (r_N(n, i) > shape_function_tolerance) {
                        const auto& r_node = r_geometry.GetPoint(i);
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
                    }
                }
            }
        }
        else
        {
            rElementalDofList.reserve(3*number_of_nodes + 3*GetNumberOfNonZeroNodes() + GetNumberOfNonZeroDOFsRotation());

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                const auto& r_node = r_geometry[i];
                rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
                rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
                rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            }

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                for (IndexType n = 0; n < r_N.size1(); ++n) {
                    if (r_N(n, i) > shape_function_tolerance) {
                        const auto& r_node = r_geometry.GetPoint(i);
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
                    }
                }
            }

            const Vector phi_r = m_phi_r_vector[0];
            for (IndexType m = 0; m < phi_r.size(); ++m) {
                IndexType a = m / 3;
                IndexType b = m % 3;

                if ((std::abs(phi_r(m)) > shape_function_tolerance) && (r_N(0, a) > shape_function_tolerance)) {
                    const auto& r_node = r_geometry[a];
                    if(b == 0){
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_X));
                    }
                    else if(b == 1){
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_Y));
                    }
                    else if(b == 2){
                        rElementalDofList.push_back(r_node.pGetDof(VECTOR_LAGRANGE_MULTIPLIER_ROTATION_Z));
                    }
                }
            }

        }
    }

    std::size_t SupportLagrangeCondition::GetNumberOfNonZeroNodes() const {
        const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

        SizeType counter = 0;
        for (IndexType n = 0; n < r_N.size1(); ++n) {
            for (IndexType m = 0; m < r_N.size2(); ++m) {
                if (r_N(n, m) > shape_function_tolerance) {
                    counter++;
                }
            }
        }
        return counter;
    }

    std::size_t SupportLagrangeCondition::GetNumberOfNonZeroDOFsRotation() const {
        const Vector phi_r = m_phi_r_vector[0];
        const Matrix N = GetGeometry().ShapeFunctionsValues();

        SizeType counter = 0;
        for (IndexType m = 0; m < phi_r.size(); ++m) {
            IndexType a = m / 3;

            if ((std::abs(phi_r(m)) > shape_function_tolerance) && (N(0, a) > shape_function_tolerance)) {
                counter++;
            }
        }    

        return counter;
    }

} // Namespace Kratos
