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
#include "custom_conditions/support_laplacian_condition.h"

namespace Kratos
{
    void SupportLaplacianCondition::CalculateAll(
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
        
        const unsigned int dim = DN_De[0].size2();
        Matrix DN_DX(number_of_nodes,dim);

        GeometryType::JacobiansType J0;
        Matrix InvJ0(dim,dim);
        r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
        Matrix Jacobian = ZeroMatrix(3,3);
        double DetJ0;
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);
        Jacobian(2,2) = 1.0;
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
        array_1d<double, 3> tangent_parameter_space;
        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        Vector add_factor = prod(Jacobian, tangent_parameter_space);
        add_factor[2] = 0.0;
        DetJ0 = norm_2(add_factor);

        //------------------------------------------------------
        // PRINT BOUNADRY GAUSS POINT TO FILE
        //------------------------------------------------------

        std::ofstream output_file("txt_files/support_boundary_GPs_3D.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z()  << std::endl;
            output_file.close();
        }
        //------------------------------------------------------

        typedef typename GeometryType::CoordinatesArrayType  CoordinatesArrayType;
        
        // Compute the normals
        array_1d<double, 3> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;
        r_geometry.Calculate(NORMAL, normal_parameter_space);

        normal_physical_space = normal_parameter_space; // prod(trans(J0[0]),normal_parameter_space);
        
        const Matrix& N = r_geometry.ShapeFunctionsValues();

        // Time-related variables
        const double theta = 1.0; 

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = DN_De[0];
        
        Matrix H = ZeroMatrix(1, number_of_nodes);
        Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
        Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
        Vector H_vector = ZeroVector(number_of_nodes);
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)            = N(0, i);
            H_vector(i)        = N(0, i); 

            for (IndexType idim = 0; idim < dim; idim++) {
                DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];           
                DN_dot_n_vec(i)  += DN_DX(i, idim) * normal_physical_space[idim];
            } 
            
        }
        // Differential area
        double penalty_integration = penalty * integration_points[0].Weight() * std::abs(DetJ0); //  * std::abs(DetJ0);

        // Collins, Lozinsky & Scovazzi innovation
        double Guglielmo_innovation = 1.0;  // = 1 -> Penalty approach
                                            // = -1 -> Free-penalty approach
        if (penalty == -1.0) {
            penalty_integration = 0.0;
            Guglielmo_innovation = -1.0;
        }
        

        // Assembly
        noalias(rLeftHandSideMatrix) -= theta * prod(trans(H), H) * penalty_integration;
        // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
        noalias(rLeftHandSideMatrix) -= theta * prod(trans(H), DN_dot_n)                        * integration_points[0].Weight()*std::abs(DetJ0); // * std::abs(DetJ0) ;
        // Of the Dirichlet BCs -(GRAD_w* n,u) 
        noalias(rLeftHandSideMatrix) -= theta * Guglielmo_innovation * prod(trans(DN_dot_n), H) * integration_points[0].Weight()*std::abs(DetJ0); // * std::abs(DetJ0) ;


        const double u_D_scalar = this->GetValue(TEMPERATURE);
        const double u_D_scalar_old = this->GetValue(TEMPERATURE_OLD_IT);
        
        noalias(rRightHandSideVector) -=  theta * H_vector * u_D_scalar * penalty_integration;
        noalias(rRightHandSideVector) -=  (1.0-theta) * H_vector * u_D_scalar_old * penalty_integration;
        // Of the Dirichlet BCs
        noalias(rRightHandSideVector) -= theta * Guglielmo_innovation * DN_dot_n_vec * u_D_scalar * integration_points[0].Weight()*std::abs(DetJ0); // * std::abs(DetJ0);
        noalias(rRightHandSideVector) -= (1.0-theta) * Guglielmo_innovation * DN_dot_n_vec * u_D_scalar_old * integration_points[0].Weight()*std::abs(DetJ0);

        // RHS contributions for (1 - theta)-weighted terms
        double temperature_previous = 0.0;
        double grad_u_previous_dot_normal = 0.0;
        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            temperature_previous += r_geometry[i].FastGetSolutionStepValue(TEMPERATURE, 1) * N(0, i);     // T^(n)
            for (unsigned int idim = 0; idim < dim; ++idim) {
                grad_u_previous_dot_normal += DN_DX(i, idim) * normal_physical_space[idim] * r_geometry[i].FastGetSolutionStepValue(TEMPERATURE, 1);
            }
        }
  
        // Residual for (1 - theta)-weighted terms
        noalias(rRightHandSideVector) += (1.0 - theta) * H_vector * temperature_previous * penalty_integration;
        noalias(rRightHandSideVector) += (1.0 - theta) * H_vector * grad_u_previous_dot_normal * integration_points[0].Weight()*std::abs(DetJ0); 
        noalias(rRightHandSideVector) += (1.0 - theta) * Guglielmo_innovation * DN_dot_n_vec * temperature_previous * integration_points[0].Weight()*std::abs(DetJ0);

        Vector temp(number_of_nodes);
        // RHS = ExtForces - K*temp;
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            temp[i] = r_geometry[i].GetSolutionStepValue(TEMPERATURE);
        }
        // RHS -= K*temp
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

        // for (unsigned int i = 0; i < number_of_nodes; i++) {
        //     std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
        //     outputFile << r_geometry[i].GetId() << "  " << r_geometry[i].GetDof(TEMPERATURE).EquationId() <<"\n";
        //     outputFile.close();
        // }
        KRATOS_CATCH("")
    }

    int SupportLaplacianCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportLaplacianCondition" << std::endl;
        return 0;
    }

    void SupportLaplacianCondition::EquationIdVector(
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

    void SupportLaplacianCondition::GetDofList(
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

    void SupportLaplacianCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        const double u_D_scalar = this->GetValue(TEMPERATURE);
        SetValue(TEMPERATURE_OLD_IT, u_D_scalar);
    }


} // Namespace Kratos