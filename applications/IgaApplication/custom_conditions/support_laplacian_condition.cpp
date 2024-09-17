//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
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
    void SupportLaplacianCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType mat_size = GetGeometry().size() * 1;

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);

        KRATOS_CATCH("")
    }
    
    
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
        
        // Compute the normals
        array_1d<double, 3> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
        normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);

        normal_physical_space = normal_parameter_space; // prod(trans(J0[0]),normal_parameter_space);
        
        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            const Matrix& N = r_geometry.ShapeFunctionsValues();

            // Calculating the cartesian derivatives
            noalias(DN_DX) = DN_De[point_number]; // prod(DN_De[point_number],InvJ0);
            
            Matrix H = ZeroMatrix(1, number_of_nodes);
            Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
            Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
            Vector H_vector = ZeroVector(number_of_nodes);
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i)            = N(point_number, i);
                H_vector(i)        = N(point_number, i); 

                for (unsigned int idim = 0; idim < dim; idim++) {
                    DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];           
                    DN_dot_n_vec(i)  += DN_DX(i, idim) * normal_physical_space[idim];
                } 
                
            }
            // Differential area
            double penalty_integration = penalty * integration_points[point_number].Weight(); //  * std::abs(DetJ0);

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
            noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n) * integration_points[point_number].Weight(); // * std::abs(DetJ0) ;

            // Assembly Dirichlet BCs -(GRAD_w* n,u) 
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H) * integration_points[point_number].Weight(); // * std::abs(DetJ0) ;


            if (CalculateResidualVectorFlag) {
                
                const double u_D_scalar = this->GetValue(TEMPERATURE);

                noalias(rRightHandSideVector) -=  H_vector * u_D_scalar * penalty_integration;
                noalias(rRightHandSideVector) -= Guglielmo_innovation * DN_dot_n_vec * u_D_scalar * integration_points[point_number].Weight(); // * std::abs(DetJ0);

                Vector temp(number_of_nodes);
                // RHS = ExtForces - K*temp;
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    temp[i] = r_geometry[i].GetSolutionStepValue(TEMPERATURE);
                }
                // RHS -= K*temp
                noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

            }
        }
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


} // Namespace Kratos
