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

        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);

        // Shape function derivatives (NEW) 
        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        const unsigned int dim = 2;
        Matrix DN_DX(number_of_nodes,3);
        Matrix InvJ0(dim,dim);

        // Compute the normals
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 2> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        Vector GP_parameter_coord(2); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]);
        
        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);

        // Stampa su file esterno le coordinate (projection[0],projection[1])
        std::ofstream outputFile("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << GP_parameter_coord[0] << " " << GP_parameter_coord[1]  <<"\n";
        outputFile.close();
        

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            const Matrix& N = r_geometry.ShapeFunctionsValues();

            Matrix Jacobian = ZeroMatrix(2,2);
            Jacobian(0,0) = J0[point_number](0,0);
            Jacobian(0,1) = J0[point_number](0,1);
            Jacobian(1,0) = J0[point_number](1,0);
            Jacobian(1,1) = J0[point_number](1,1);

            // Calculating inverse jacobian and jacobian determinant
            MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
            Matrix InvJ0_23 = ZeroMatrix(2,3);
            InvJ0_23(0,0) = InvJ0(0,0);
            InvJ0_23(0,1) = InvJ0(0,1);
            InvJ0_23(1,0) = InvJ0(1,0);
            InvJ0_23(1,1) = InvJ0(1,1);
            InvJ0_23(0,2) = 0;
            InvJ0_23(1,2) = 0;

            // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            noalias(DN_DX) = prod(DN_De[point_number],InvJ0_23);
            
            Matrix H = ZeroMatrix(1, number_of_nodes);
            Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
            Vector DN_dot_n_vector = ZeroVector(number_of_nodes);
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i)                = N(point_number, i);
                DN_dot_n(0, i)         = DN_DX(i, 0) * normal_physical_space[0] + DN_DX(i, 1) * normal_physical_space[1] ;
                DN_dot_n_vector(i)     = DN_DX(i, 0) * normal_physical_space[0] + DN_DX(i, 1) * normal_physical_space[1] ;
            }

            // Differential area
            double penalty_integration = penalty * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);

            // Guglielmo innovaction
            double Guglielmo_innovation = -1.0;  // = 1 -> Penalty approach
                                                 // = -1 -> Free-penalty approach
            if (Guglielmo_innovation == -1.0) {
                penalty_integration = 0.0;
            }
            // double tau = 0.0 ;
            // // laplacian N
            // Matrix laplacian_of_H = ZeroMatrix(1, number_of_nodes);
            // Matrix gradient_of_laplacian_of_H_dot_n = ZeroMatrix(1, number_of_nodes);
            // const Matrix DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, 0, r_geometry.GetDefaultIntegrationMethod());
            // const Matrix DDDN_DDDe = r_geometry.ShapeFunctionDerivatives(3, 0, r_geometry.GetDefaultIntegrationMethod());
            // for(unsigned int i=0; i<number_of_nodes; i++) {
            //         laplacian_of_H(0, i)                    = DDN_DDe(i,0) + DDN_DDe(i,2) ;
            //         gradient_of_laplacian_of_H_dot_n(0, i)  = (DDDN_DDDe(i,0) + DDDN_DDDe(i,2)) * normal_physical_space[0] + 
            //                                                   (DDDN_DDDe(i,1) + DDDN_DDDe(i,3)) * normal_physical_space[1] ;
            // }

            // Assembly
            noalias(rLeftHandSideMatrix) -= prod(trans(H), H) * penalty_integration;
            // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
            noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n)                        * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            // Of the Dirichlet BCs -(GRAD_w * n,u) 
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            // ________________________________________________________________________________________________________________________________________
            // // Convection Additional Terms
            // noalias(rLeftHandSideMatrix) -=  tau * prod(trans(H), DN_dot_n)                          * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            // noalias(rLeftHandSideMatrix) +=  tau * prod(trans(DN_dot_n), laplacian_of_H)             * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            // noalias(rLeftHandSideMatrix) -=  tau * prod(trans(H), gradient_of_laplacian_of_H_dot_n ) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;


            if (CalculateResidualVectorFlag) {
                
                // double& temperature = this->GetValue(TEMPERATURE) ;

                // double temperature = GP_parameter_coord[0]-GP_parameter_coord[1];
                double temperature = GP_parameter_coord[0]*GP_parameter_coord[0] + GP_parameter_coord[1]*GP_parameter_coord[1];
                // double temperature = sin(GP_parameter_coord[0]) * sinh(GP_parameter_coord[1]); // CAMBIA!!!!
                // double temperature = GP_parameter_coord[0]*GP_parameter_coord[0]*GP_parameter_coord[0] + GP_parameter_coord[1]*GP_parameter_coord[1]*GP_parameter_coord[1] ;
                // double temperature = GP_parameter_coord[0]*GP_parameter_coord[0]*GP_parameter_coord[0]*GP_parameter_coord[0]  + GP_parameter_coord[1]*GP_parameter_coord[1]*GP_parameter_coord[1] *GP_parameter_coord[1] ;
                
                Vector u_D(number_of_nodes);
                for (IndexType i = 0; i < number_of_nodes; ++i)
                {
                    
                    u_D[i] = - temperature;
                }
                noalias(rRightHandSideVector) += prod(prod(trans(H), H), u_D) * penalty_integration;
                // Of the Dirichlet BCs
                noalias(rRightHandSideVector) += Guglielmo_innovation * DN_dot_n_vector * (-temperature) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);
                // // Convection Additional Terms
                // noalias(rRightHandSideVector) -= tau *                  DN_dot_n_vector * (-temperature) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);
                // noalias(rRightHandSideVector) += tau *               prod(prod(trans(DN_dot_n), H), u_D) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);

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
        //     // KRATOS_WATCH(r_geometry[i].GetDof(TEMPERATURE).EquationId() )
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
