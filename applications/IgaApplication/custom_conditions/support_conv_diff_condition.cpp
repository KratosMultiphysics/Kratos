//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_conv_diff_condition.h"

namespace Kratos
{
    void SupportConvDiffCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        const double penalty = GetProperties()[PENALTY_FACTOR];

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
        
        GeometryType::JacobiansType J0;
        const unsigned int dim = DN_De[0].size2();
        Matrix DN_DX(number_of_nodes,dim);
        Matrix InvJ0(dim,dim);
        
        r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());

        double DetJ0;

        Vector GP_parameter_coord(3); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]);

        //------------------------------------------------------
        // PRINT BOUNADRY GAUSS POINT TO FILE
        //------------------------------------------------------

        std::ofstream output_file("txt_files/support_boundary_GPs_3D.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z()  << std::endl;
            output_file.close();
        }

        typedef typename GeometryType::CoordinatesArrayType  CoordinatesArrayType;
        
        
        // Compute the normals
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 3> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        // r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        // double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // KRATOS_WATCH(tangent_parameter_space)
        // // NEW FOR GENERAL JACOBIAN
        // normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        // normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT
        // normal_parameter_space[2] = 0;
        
        r_geometry.Calculate(NORMAL, normal_parameter_space);

        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);

        // KRATOS_WATCH(normal_parameter_space)

        std::ofstream outputFile("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << GP_parameter_coord[0] << " " << GP_parameter_coord[1] << " " << GP_parameter_coord[2] <<"\n";
        outputFile.close();
        

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            const Matrix& N = r_geometry.ShapeFunctionsValues();

            Matrix Jacobian = ZeroMatrix(dim,dim);
            Jacobian(0,0) = J0[0](0,0);
            Jacobian(0,1) = J0[0](0,1);
            Jacobian(1,0) = J0[0](1,0);
            Jacobian(1,1) = J0[0](1,1);
            if (dim > 2) {
                Jacobian(0,2) = J0[0](0,2);
                Jacobian(1,2) = J0[0](1,2);
                Jacobian(2,0) = J0[0](2,0);
                Jacobian(2,1) = J0[0](2,1);
                Jacobian(2,2) = J0[0](2,2);
            }

            // Calculating inverse jacobian and jacobian determinant
            MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

            // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            noalias(DN_DX) = prod(DN_De[point_number],InvJ0);
            
            Matrix H = ZeroMatrix(1, number_of_nodes);
            Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
            Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
            Vector H_vector = ZeroVector(number_of_nodes);
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i)            = N(point_number, i);
                H_vector(i)        = N(point_number, i); 

                for (IndexType idim = 0; idim < dim; idim++) {
                    DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];           
                    DN_dot_n_vec(i)  += DN_DX(i, idim) * normal_physical_space[idim];
                } 
            }
            // Differential area
            double penalty_integration = penalty * integration_points[point_number].Weight() * std::abs(DetJ0);

            // Collins, Lozinsky & Scovazzi innovation
            double Guglielmo_innovation = 1.0;  // = 1 -> Penalty approach
                                                // = -1 -> Free-penalty approach
            if (penalty == -1.0) {
                penalty_integration = 0.0;
                Guglielmo_innovation = -1.0;
            }
            

            // Assembly
            noalias(rLeftHandSideMatrix) -= prod(trans(H), H) * penalty_integration;
            // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
            noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n)                        * integration_points[point_number].Weight() * std::abs(DetJ0) ;
            // Of the Dirichlet BCs -(GRAD_w* n,u) 
            // noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H) * integration_points[point_number].Weight() * std::abs(DetJ0) ;
            
            // Additional term Burman Penalty-free for convection-diffusion
            if (penalty == -1.0) {
                array_1d<double, 3> v;
                v[0]=this->GetValue(VELOCITY_X);
                v[1]=this->GetValue(VELOCITY_Y);

                if (inner_prod(normal_physical_space, v) < 0 && norm_2(v) > 300.0 ) {
                    noalias(rLeftHandSideMatrix) += prod(trans(H), H) * inner_prod(normal_physical_space, v)  * integration_points[point_number].Weight() * std::abs(DetJ0) ;
                }
                else{
                    // Penalty-free POISSON
                    noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H) * integration_points[point_number].Weight() * std::abs(DetJ0) ;
                }
            }
            

            if (CalculateResidualVectorFlag) {
                
                const double u_D_scalar = this->GetValue(TEMPERATURE);

                noalias(rRightHandSideVector) -=  H_vector * u_D_scalar * penalty_integration;
                // Of the Dirichlet BCs
                // noalias(rRightHandSideVector) -= Guglielmo_innovation * DN_dot_n_vec * u_D_scalar * integration_points[point_number].Weight() * std::abs(DetJ0);

                // Additional term Burman Penalty-free for convection-diffusion
                if (penalty == -1.0) {
                    array_1d<double, 3> v;
                    v[0]=this->GetValue(VELOCITY_X);
                    v[1]=this->GetValue(VELOCITY_Y);
                    if (inner_prod(normal_physical_space, v) < 0 && norm_2(v) > 300.0) {
                        noalias(rRightHandSideVector) += H_vector * u_D_scalar * inner_prod(normal_physical_space, v)  * integration_points[point_number].Weight() * std::abs(DetJ0) ;
                    }else{
                        // Penalty-free POISSON
                        noalias(rRightHandSideVector) -= Guglielmo_innovation * DN_dot_n_vec * u_D_scalar * integration_points[point_number].Weight() * std::abs(DetJ0);
                    }
                }

                Vector temp(number_of_nodes);
                // RHS = ExtForces - K*temp;
                for (unsigned int i = 0; i < number_of_nodes; i++) {
                    temp[i] = r_geometry[i].GetSolutionStepValue(TEMPERATURE);
                }
                // RHS -= K*temp
                noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

            }
        }
        // for (unsigned int i = 0; i < number_of_nodes; i++) {
        //     std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
        //     outputFile << r_geometry[i].GetId() << "  " << r_geometry[i].GetDof(TEMPERATURE).EquationId() <<"\n";
        //     outputFile.close();
        // }
        KRATOS_CATCH("")
    }

    int SupportConvDiffCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportConvDiffCondition" << std::endl;
        return 0;
    }

    void SupportConvDiffCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() !=  number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            // const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[i] = r_node.GetDof(TEMPERATURE).EquationId();
        }
    }

    void SupportConvDiffCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
        }
    };


} // Namespace Kratos
