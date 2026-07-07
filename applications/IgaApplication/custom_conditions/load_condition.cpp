//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes

// External includes

// Project includes
#include "custom_conditions/load_condition.h"


namespace Kratos
{
    void LoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = r_geometry.WorkingSpaceDimension() * number_of_nodes;

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

        // Calculation of Force vector
        if (CalculateResidualVectorFlag) {
            Vector f = ZeroVector(mat_size);

            // Integration
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

            Vector determinat_jacobian_vector(integration_points.size());
            r_geometry.DeterminantOfJacobian(determinat_jacobian_vector);

            // initial determinant of jacobian for dead load
            Vector determinat_jacobian_vector_initial(integration_points.size());
            if (this->Has(DEAD_LOAD))
            {
                DeterminantOfJacobianInitial(r_geometry, determinat_jacobian_vector_initial);
            }

            // Shape function values for all integration points
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                // Differential area
                const double integration_weight = integration_points[point_number].Weight();
                const double d_weight = integration_weight * determinat_jacobian_vector[point_number];

                // Split only due to different existing variable names
                // No check included, which checks correctness of variable

                // Dead load. Dead load is dependent on the initial area.
                if (this->Has(DEAD_LOAD))
                {
                    const array_1d<double, 3>& dead_load = this->GetValue(DEAD_LOAD);

                    const double d0_weight = integration_weight * determinat_jacobian_vector_initial[point_number];

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += dead_load[0] * r_N(point_number, i) * d0_weight;
                        f[index + 1] += dead_load[1] * r_N(point_number, i) * d0_weight;
                        f[index + 2] += dead_load[2] * r_N(point_number, i) * d0_weight;
                    }
                }

                // Point loads
                if (this->Has(POINT_LOAD))
                {
                    const array_1d<double, 3>& point_load = this->GetValue(POINT_LOAD);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += point_load[0] * r_N(point_number, i);
                        f[index + 1] += point_load[1] * r_N(point_number, i);
                        f[index + 2] += point_load[2] * r_N(point_number, i);
                    }
                }

                // Line loads
                if (this->Has(LINE_LOAD))
                {
                    const array_1d<double, 3>& line_load = this->GetValue(LINE_LOAD);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += line_load[0] * r_N(point_number, i) * d_weight;
                        f[index + 1] += line_load[1] * r_N(point_number, i) * d_weight;
                        f[index + 2] += line_load[2] * r_N(point_number, i) * d_weight;
                    }
                }

                // Moment line loads
                if (this->Has(MOMENT_LINE_LOAD))
                {
                    const array_1d<double, 3>& moment_line_load = this->GetValue(MOMENT_LINE_LOAD);

                    Vector phi_r = ZeroVector(mat_size);
                    Matrix phi_rs = ZeroMatrix(mat_size, mat_size);
                    array_1d<double, 2> diff_phi;

                    CalculateRotationalShapeFunctions(point_number, phi_r, phi_rs, diff_phi);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += moment_line_load[0] * phi_r(index) * d_weight;
                        f[index + 1] += moment_line_load[1] * phi_r(index + 1) * d_weight;
                        f[index + 2] += moment_line_load[2] * phi_r(index + 2) * d_weight;
                    }
                }

                // Surface loads
                if (this->Has(SURFACE_LOAD))
                {
                    const array_1d<double, 3>& surface_load = this->GetValue(SURFACE_LOAD);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += surface_load[0] * r_N(point_number, i) * d_weight;
                        f[index + 1] += surface_load[1] * r_N(point_number, i) * d_weight;
                        f[index + 2] += surface_load[2] * r_N(point_number, i) * d_weight;
                    }
                }

                // Pressure loads
                if (this->Has(PRESSURE))
                {
                    const double pressure = this->GetValue(PRESSURE);

                    array_1d<double, 3> normal = r_geometry.Normal(point_number);
                    normal = normal / norm_2(normal);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += normal[0] * pressure * r_N(point_number, i) * d_weight;
                        f[index + 1] += normal[1] * pressure * r_N(point_number, i) * d_weight;
                        f[index + 2] += normal[2] * pressure * r_N(point_number, i) * d_weight;
                    }
                }

                // Pressure follower loads
                if (this->Has(PRESSURE_FOLLOWER_LOAD))
                {
                    const double pressure_follower_load = this->GetValue(PRESSURE_FOLLOWER_LOAD);
                    const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(point_number);

                    array_1d<double, 3> normal = r_geometry.Normal(point_number);
                    normal = normal / norm_2(normal);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += normal[0] * pressure_follower_load * r_N(point_number, i) * d_weight;
                        f[index + 1] += normal[1] * pressure_follower_load * r_N(point_number, i) * d_weight;
                        f[index + 2] += normal[2] * pressure_follower_load * r_N(point_number, i) * d_weight;
                    }

                    // compute the load stiffness matrix due to the follower load
                    // a. compute the basis functions and its first derivative
                    Matrix N = ZeroMatrix(3, mat_size);
                    Matrix DN_DXi = ZeroMatrix(3, mat_size);
                    Matrix DN_DEta = ZeroMatrix(3, mat_size);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;

                        N(0, index) = r_N(point_number, i);
                        N(1, index + 1) = r_N(point_number, i);
                        N(2, index + 2) = r_N(point_number, i);

                        DN_DXi(0, index) = r_DN_De(i, 0);
                        DN_DXi(1, index + 1) = r_DN_De(i, 0);
                        DN_DXi(2, index + 2) = r_DN_De(i, 0);

                        DN_DEta(0, index) = r_DN_De(i, 1);
                        DN_DEta(1, index + 1) = r_DN_De(i, 1);
                        DN_DEta(2, index + 2) = r_DN_De(i, 1);
                    }

                    // b. get the base vectors a1 and a2
                    Matrix J;
                    r_geometry.Jacobian(J, point_number);

                    array_1d<double, 3> a1 = column(J, 0);
                    array_1d<double, 3> a2 = column(J, 1);

                    // c. formulate the cross product (skew symmetric) matrices a1x and a2x
                    Matrix a1x = ZeroMatrix(3,3);
                    Matrix a2x = ZeroMatrix(3,3);

                    a1x(0,1) = -a1[2];
                    a1x(0,2) = a1[1];
                    a1x(1,2) = -a1[0];
                    a1x(1,0) = -a1x(0,1);
                    a1x(2,0) = -a1x(0,2);
                    a1x(2,1) = -a1x(1,2);

                    a2x(0,1) = a2[2];
                    a2x(0,2) = -a2[1];
                    a2x(1,2) = a2[0];
                    a2x(1,0) = a2x(0,1);
                    a2x(2,0) = a2x(0,2);
                    a2x(2,1) = a2x(1,2);

                    // d. compute the left hand side
                    noalias(rLeftHandSideMatrix) -= integration_weight * pressure_follower_load * (prod(trans(N), (prod(a2x, DN_DXi) + prod(a1x, DN_DEta))));
                }

                // Assembly
                noalias(rRightHandSideVector) += f;
            }
        }
    }

    void LoadCondition::DeterminantOfJacobianInitial(
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

            rDeterminantOfJacobian[pnt] = MathUtils<double>::GeneralizedDet(J);
        }
    }

    void LoadCondition::CalculateRotationalShapeFunctions(
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

    void LoadCondition::CalculateRotation(
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
                phi_r(nb_dof) = 1.0 / sqrt(1.0 - pow(sinus_omega(0), 2))*inner_prod(sinus_omega_r[nb_dof], T2); //tangent
                // phi_r(nb_dof) = 1.0 / sqrt(1.0 - pow(sinus_omega(1), 2))*inner_prod(sinus_omega_r[nb_dof], T1); //normal
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

    void LoadCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void LoadCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    };

} // Namespace Kratos
