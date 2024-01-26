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
#include "custom_conditions/sbm_laplacian_condition.h"
#include "containers/model.h"

namespace Kratos
{
    void SBMLaplacianCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
        Condition candidateClosestSkinSegment2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1];

        
        KRATOS_TRY
        const auto& r_geometry = this->GetGeometry();
        const SizeType number_of_nodes = r_geometry.PointsNumber();
        if (rRightHandSideVector.size() != number_of_nodes) {
            rRightHandSideVector.resize(number_of_nodes, false);
        }
        if (rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes) {
            rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
        }

        noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);

        double penalty = GetProperties()[PENALTY_FACTOR];

        // Read the refinements.iga.json
        const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");
        int insertions = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        double h = 2.0/(insertions+1) ;

        // Modify the penalty factor: penalty/h
        penalty = penalty/(h);

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        
        r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);  // = 1

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        const unsigned int dim = 2;
        Matrix DN_DX(number_of_nodes,dim);
        Matrix DN_DPSI(number_of_nodes,dim);
        Matrix InvJ0(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            // Note that integration_points.size() = 1   &   number_of_nodes = 9
            // Differential area
            double penalty_integration = penalty * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);

            // Get the parameter coordinates
            Vector GP_parameter_coord(2); 
            GP_parameter_coord = prod(r_geometry.Center(),J0[point_number]);


            // std::ifstream file("txt_files/true_points.txt");    // Read the true points from the mdpa
            // std::vector<double> x_true_boundary;
            // std::vector<double> y_true_boundary;
            // double x, y;
            // while (file >> x >> y) { x_true_boundary.push_back(x); y_true_boundary.push_back(y); }
            // file.close();
            // int index_min_distance = 0; // Initialization of the index of the true node closest to the Gauss point.
            // double min_distance_squared = 1e14;
            // for (int i = 1; i < x_true_boundary.size(); i++) {
            //     double current_distance_squared = (GP_parameter_coord[0]-x_true_boundary[i])*(GP_parameter_coord[0]-x_true_boundary[i]) +  
            //           (GP_parameter_coord[1]-y_true_boundary[i])*(GP_parameter_coord[1]-y_true_boundary[i]) ;
            //     if ( current_distance_squared  <  min_distance_squared ) {
            //             min_distance_squared = current_distance_squared;
            //             index_min_distance = i ;
            //           }
            // }
            // Vector projection(2);
            // projection[0] = x_true_boundary[index_min_distance] ;
            // projection[1] = y_true_boundary[index_min_distance] ;
            
            // Obtaining the projection from the closest skin segment
            Vector projection(2);
            projection[0] = candidateClosestSkinSegment1.GetGeometry()[0].X() ;
            projection[1] = candidateClosestSkinSegment1.GetGeometry()[0].Y() ;

            // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
            std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
            outputFile << projection[0] << " " << projection[1] << " "  << GP_parameter_coord[0] << " " << GP_parameter_coord[1] <<"\n";
            outputFile.close();

            Vector d(2);
            d[0] = projection[0] - GP_parameter_coord[0];
            d[1] = projection[1] - GP_parameter_coord[1];
            // d[0] = 0;
            // d[1] = 0;

            const Matrix& N = r_geometry.ShapeFunctionsValues();

            Matrix Jacobian = ZeroMatrix(2,2);
            Jacobian(0,0) = J0[point_number](0,0);
            Jacobian(0,1) = J0[point_number](0,1);
            Jacobian(1,0) = J0[point_number](1,0);
            Jacobian(1,1) = J0[point_number](1,1);

            // Calculating inverse jacobian and jacobian determinant
            MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

            // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
            noalias(DN_DX) = prod(DN_De[point_number],InvJ0);

            const Matrix& DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, point_number, this->GetIntegrationMethod());
            const Matrix& DDDN_DDDe = r_geometry.ShapeFunctionDerivatives(3, point_number, this->GetIntegrationMethod());
            const Matrix& DDDDN_DDDDe = r_geometry.ShapeFunctionDerivatives(4, point_number, this->GetIntegrationMethod());
            const Matrix& DDDDDN_DDDDDe = r_geometry.ShapeFunctionDerivatives(5, point_number, this->GetIntegrationMethod());

            // Calculating the PARAMETER SPACE derivatives
            Matrix Identity_Matrix = ZeroMatrix(2,2);
            Identity_Matrix(0,0) = 1.0;
            Identity_Matrix(1,1) = 1.0;
            noalias(DN_DPSI) = prod(DN_De[point_number],Identity_Matrix);

            // Compute the normals
            array_1d<double, 3> tangent_parameter_space;
            array_1d<double, 3> normal_parameter_space;
            array_1d<double, 2> normal_physical_space;

            r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
            double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
            
            // NEW FOR GENERAL JACOBIAN
            normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
            normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT

            normal_physical_space = prod(trans(J0[0]),normal_parameter_space);
            
            // Guglielmo innovaction
            double Guglielmo_innovation = -1.0;  // = 1 -> Penalty approach
                                                // = -1 -> Free-penalty approach
            if (Guglielmo_innovation < 0.0) {
                penalty_integration = 0.0;
            }

            
            Matrix H = ZeroMatrix(1, number_of_nodes);
            Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
            Matrix DN_dot_n_parameter = ZeroMatrix(1, number_of_nodes);
            // Need to modify H so that we introduce the Taylor expansion
            Matrix H_gradient_term = ZeroMatrix(1, number_of_nodes);
            Matrix H_hessian_term = ZeroMatrix(1, number_of_nodes);
            Matrix H_3rdTayor_term = ZeroMatrix(1, number_of_nodes);
            Matrix H_4thTayor_term = ZeroMatrix(1, number_of_nodes);
            Matrix H_5thTayor_term = ZeroMatrix(1, number_of_nodes);

            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i)               = N(point_number, i);
                DN_dot_n(0, i)        = DN_DX(i, 0) * normal_physical_space[0] + DN_DX(i, 1) * normal_physical_space[1] ;
                // Taylor expansion in the parameter space
                H_gradient_term(0, i) = DN_DPSI(i, 0) * d[0] + DN_DPSI(i, 1) * d[1] ; 
                H_hessian_term(0, i) = 1.0/2.0 *  d[0]*d[0]*DDN_DDe(i,0) + 1.0 * d[0]*d[1]*DDN_DDe(i,1) + 1.0/2.0 * d[1]*d[1]*DDN_DDe(i,2)  ;
                H_3rdTayor_term(0, i) = 1.0/6.0 * DDDN_DDDe(i,0) * d[0]*d[0]*d[0] + 1.0/2.0 * DDDN_DDDe(i,1)*d[0]*d[0]*d[1] + 1.0/2.0* DDDN_DDDe(i,2)*d[0]*d[1]*d[1] + 1.0/6.0 *DDDN_DDDe(i,3)*d[1]*d[1]*d[1] ;
                H_4thTayor_term(0, i) = 1.0/24.0 * DDDDN_DDDDe(i,0) * d[0]*d[0]*d[0]*d[0] + 1.0/6.0 * DDDDN_DDDDe(i,1) * d[0]*d[0]*d[0]*d[1] + 1.0/4.0 * DDDDN_DDDDe(i,2) * d[0]*d[0]*d[1]*d[1] + 1.0/6.0 * DDDDN_DDDDe(i,3) * d[0]*d[1]*d[1]*d[1] + 1.0/24.0 * DDDDN_DDDDe(i,4) * d[1]*d[1]*d[1]*d[1];
                H_5thTayor_term(0, i) = 1.0/120.0 * DDDDDN_DDDDDe(i,0) * d[0]*d[0]*d[0]*d[0]*d[0] + 1.0/24.0 * DDDDDN_DDDDDe(i,1) * d[0]*d[0]*d[0]*d[0]*d[1] + 1.0/12.0 * DDDDDN_DDDDDe(i,2) * d[0]*d[0]*d[0]*d[1]*d[1] + 1.0/12.0 * DDDDDN_DDDDDe(i,3) * d[0]*d[0]*d[1]*d[1]*d[1] + 1.0/24.0 * DDDDDN_DDDDDe(i,4) * d[0]*d[1]*d[1]*d[1]*d[1] + 1.0/120.0 * DDDDDN_DDDDDe(i,5) * d[1]*d[1]*d[1]*d[1]*d[1];
                // H_hessian_term(0, i) = 0;
                // H_3rdTayor_term(0, i) = 0;
                H_4thTayor_term(0, i) = 0;
                H_5thTayor_term(0, i) = 0;

            }
            // Assembly

            // Termine -(GRAD_w * n, u + GRAD_u * d + ...)
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H)                * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_gradient_term)  * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_hessian_term)   * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_3rdTayor_term)  * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_4thTayor_term)  * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_5thTayor_term)  * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;

            // Termine -(w,GRAD_u * n) from integration by parts -> Fundamental !! 
            noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n)                                      * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]) ;
            // SBM terms (Taylor Expansion) + alpha * (w + GRAD_w * d + ..., u + GRAD_u * d + ...)
            noalias(rLeftHandSideMatrix) += prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H              ) * penalty_integration ;
            noalias(rLeftHandSideMatrix) += prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H_gradient_term) * penalty_integration ;
            noalias(rLeftHandSideMatrix) += prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H_hessian_term ) * penalty_integration ;
            noalias(rLeftHandSideMatrix) += prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H_3rdTayor_term) * penalty_integration ;
            noalias(rLeftHandSideMatrix) += prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H_4thTayor_term) * penalty_integration ;
            noalias(rLeftHandSideMatrix) += prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H_5thTayor_term) * penalty_integration ;


            if (CalculateResidualVectorFlag) {
                
                // const double& temperature1 = Has(TEMPERATURE)
                //     ? this->GetValue(TEMPERATURE)
                //     : 0.0;
                // KRATOS_WATCH(temperature1)
                
                // const double temperature =  projection[0]*(projection[0]-2.0) * projection[1]*(projection[1]-2.0) ;// 1.0 ; // 1.0*(1.0-2.0)*1.0*(1.0-2.0) ;
                // const double temperature = projection[0]-projection[1];
                // double temperature = GP_parameter_coord[0]-GP_parameter_coord[1];
                // const double temperature = projection[0]*projection[0] + projection[1]*projection[1];
                // const double temperature = projection[0] * projection[1];
                const double temperature = sin(projection[0]) * sinh(projection[1]) ;
                // const double temperature = sin(GP_parameter_coord[0]) * sinh(GP_parameter_coord[1]) ;
                // const double temperature = projection[0]*projection[0]*projection[0] + projection[1]*projection[1]*projection[1] ;
                // const double temperature = projection[0]*projection[0]*projection[0]*projection[0] + projection[1]*projection[1]*projection[1]*projection[1] ;
                
                Vector u_D(number_of_nodes);

                for (IndexType i = 0; i < number_of_nodes; ++i)
                {
                    // const double temper = 0; // r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);      // Always = 0.  
                    // u_D[i] = (temper - temperature);
                    u_D[i] = -temperature ;
                }
                noalias(rRightHandSideVector) -= prod(prod(trans(H + H_gradient_term + H_hessian_term + H_3rdTayor_term + H_4thTayor_term + H_5thTayor_term), H), u_D) * penalty_integration;
                // Dirichlet BCs
                noalias(rRightHandSideVector) += Guglielmo_innovation * prod(prod(trans(DN_dot_n), H), u_D) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);



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
















    int SBMLaplacianCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SBMLaplacianCondition" << std::endl;
        return 0;
    }

    void SBMLaplacianCondition::EquationIdVector( // Essential
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        // std::cout << "EquationIdVector: " << this->Id() << std::endl;

        if (rResult.size() !=  number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rResult[i] = r_node.GetDof(TEMPERATURE).EquationId();
        }
    }

    void SBMLaplacianCondition::GetDofList(  // Essential
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(number_of_nodes);
        // What ".reserve" does?

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
        }
    }


    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters SBMLaplacianCondition::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        std::ifstream infile(rDataFileName);

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    };


} // Namespace Kratos
