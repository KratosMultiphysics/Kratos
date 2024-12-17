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
//                   Andrea Gorgi
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
        std::string boundaryConditionType = this->GetValue(BOUNDARY_CONDITION_TYPE);
        
        if (boundaryConditionType == "dirichlet")
        {
            this-> CalculateAllDirichlet(   rLeftHandSideMatrix,
                                            rRightHandSideVector,
                                            rCurrentProcessInfo,
                                            CalculateStiffnessMatrixFlag,
                                            CalculateResidualVectorFlag);
        }
        else if (boundaryConditionType == "neumann"){
            this-> CalculateAllNeumann(     rLeftHandSideMatrix,
                                            rRightHandSideVector,
                                            rCurrentProcessInfo,
                                            CalculateStiffnessMatrixFlag,
                                            CalculateResidualVectorFlag);
        }
        else{
            KRATOS_ERROR << "error in SBM_LAPLACIAN_CONDITION, no BOUNDARY_CONDITION_TYPE available" << std::endl;
        }
    }


    //_________________________________________________________________________________________________________________________________________
    // DIRICHLET CONDITION
    //_________________________________________________________________________________________________________________________________________

    void SBMLaplacianCondition::CalculateAllDirichlet(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;

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

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
        
        // Initialize DN_DX
        const unsigned int dim = DN_De[0].size2();
        Matrix DN_DX(number_of_nodes,dim);
        
        Vector meshSize_uv = this->GetValue(MARKER_MESHES);
        double h = std::min(meshSize_uv[0], meshSize_uv[1]);
        if (dim == 3) {h = std::min(h,  meshSize_uv[2]);}

        // Compute basis function order (Note: it is not allow to use different orders in different directions)
        if (dim == 3) {
            mbasisFunctionsOrder = std::cbrt(DN_De[0].size1()) - 1;
        } else {
            mbasisFunctionsOrder = std::sqrt(DN_De[0].size1()) - 1;
        }

        // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
        penalty = mbasisFunctionsOrder * mbasisFunctionsOrder * penalty / h;

        // Matrix InvJ0(dim,dim);
        // Vector determinant_jacobian_vector(integration_points.size());
        // r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);  // = 1
        // GeometryType::JacobiansType J0;
        // r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
        // double DetJ0;

        // Obtaining the projection from the closest skin segment

        // Find the closest node in condition
        auto GP_parameter_coord = r_geometry.Center();
        int closestNodeId;
        if (dim > 2) {
            double incumbent_dist = 1e16;
            for (int i = 0; i < dim; i++) {
                if (norm_2(candidateClosestSkinSegment1.GetGeometry()[i]-GP_parameter_coord) < incumbent_dist) {
                    incumbent_dist = norm_2(candidateClosestSkinSegment1.GetGeometry()[i]-GP_parameter_coord);
                    closestNodeId = i;
                }
            }
        }else{
            closestNodeId = 0;
        }
        Vector projection(3);
        projection[0] = candidateClosestSkinSegment1.GetGeometry()[closestNodeId].X() ;
        projection[1] = candidateClosestSkinSegment1.GetGeometry()[closestNodeId].Y() ;
        projection[2] = candidateClosestSkinSegment1.GetGeometry()[closestNodeId].Z() ;

        // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
        std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
        outputFile << projection[0] << " " << projection[1] << " " << projection[2] << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() <<"\n";
        outputFile.close();

        d.resize(3);
        d[0] = projection[0] - r_geometry.Center().X();
        d[1] = projection[1] - r_geometry.Center().Y();
        d[2] = projection[2] - r_geometry.Center().Z();
        // d[0] = 0;
        // d[1] = 0;
        // d[2] = 0;

        const Matrix& N = r_geometry.ShapeFunctionsValues();

        // Matrix Jacobian = ZeroMatrix(dim,dim);
        // Jacobian(0,0) = J0[0](0,0);
        // Jacobian(0,1) = J0[0](0,1);
        // Jacobian(1,0) = J0[0](1,0);
        // Jacobian(1,1) = J0[0](1,1);
        // if (dim > 2) {
        //     Jacobian(0,2) = J0[0](0,2);
        //     Jacobian(1,2) = J0[0](1,2);
        //     Jacobian(2,0) = J0[0](2,0);
        //     Jacobian(2,1) = J0[0](2,1);
        //     Jacobian(2,2) = J0[0](2,2);
        // }
        // // Calculating inverse jacobian and jacobian determinant
        // MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // Differential area
        double penalty_integration = penalty * integration_points[0].Weight() ; // * std::abs(DetJ0);

        // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = DN_De[0]; // prod(DN_De[0],InvJ0);

        // Compute the normals
        // array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 3> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        r_geometry.Calculate(NORMAL, normal_parameter_space);

        normal_physical_space = normal_parameter_space; // prod(trans(J0[0]),normal_parameter_space);
        
        // Collins, Lozinsky & Scovazzi innovation
        double Guglielmo_innovation = 1.0;  // = 1 -> Penalty approach
                                            // = -1 -> Free-penalty approach
        if (penalty == -1.0) {
            penalty_integration = 0.0;
            Guglielmo_innovation = -1.0;
        }

        Matrix H = ZeroMatrix(1, number_of_nodes);
        Vector H_vec = ZeroVector(number_of_nodes);
        Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
        Vector DN_dot_n_vec = ZeroVector(number_of_nodes);

        Vector H_sum_vec = ZeroVector(number_of_nodes);

        // Compute all the derivatives of the basis functions involved
        for (int n = 1; n <= mbasisFunctionsOrder; n++) {
            mShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
        }

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)               = N(0, i);
            H_vec(i)              = N(0, i);
            for (IndexType idim = 0; idim < dim; idim++) {
                    DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];
                    DN_dot_n_vec(i)  += DN_DX(i, idim) * normal_physical_space[idim];         
            } 
            // Reset for each node
            double H_taylor_term = 0.0; 

            if (dim == 2) {
                for (int n = 1; n <= mbasisFunctionsOrder; n++) {
                    // Retrieve the appropriate derivative for the term
                    Matrix& shapeFunctionDerivatives = mShapeFunctionDerivatives[n-1];
                    for (int k = 0; k <= n; k++) {
                        int n_k = n - k;
                        double derivative = shapeFunctionDerivatives(i,k); 
                        // Compute the Taylor term for this derivative
                        H_taylor_term += computeTaylorTerm(derivative, d[0], n_k, d[1], k);
                    }
                }
            } else {
                // 3D
                for (int n = 1; n <= mbasisFunctionsOrder; n++) {
                    Matrix& shapeFunctionDerivatives = mShapeFunctionDerivatives[n-1];
                    
                    int countDerivativeId = 0;
                    // Loop over blocks of derivatives in x
                    for (int k_x = n; k_x >= 0; k_x--) {
                        // Loop over the possible derivatives in y
                        for (int k_y = n - k_x; k_y >= 0; k_y--) {
                            
                            // derivatives in z
                            int k_z = n - k_x - k_y;
                            double derivative = shapeFunctionDerivatives(i,countDerivativeId); 

                            H_taylor_term += computeTaylorTerm3D(derivative, d[0], k_x, d[1], k_y, d[2], k_z);
                            countDerivativeId++;
                        }
                    }
                }
            }
            
            H_sum(0,i) = H_taylor_term + H(0,i);
            H_sum_vec(i) = H_taylor_term + H(0,i);
        }

        // Assembly
        // Termine -(GRAD_w * n, u + GRAD_u * d + ...)
        noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_sum)  * integration_points[0].Weight() ; // * std::abs(DetJ0) ;
        // Termine -(w,GRAD_u * n) from integration by parts -> Fundamental !! 
        noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n)                             * integration_points[0].Weight() ; // * std::abs(DetJ0) ;
        // SBM terms (Taylor Expansion) + alpha * (w + GRAD_w * d + ..., u + GRAD_u * d + ...)
        noalias(rLeftHandSideMatrix) += prod(trans(H_sum), H_sum) * penalty_integration ;


        if (CalculateResidualVectorFlag) {
            
            // double temperature = GP_parameter_coord[0]-GP_parameter_coord[1];
            // const double temperature = sin(projection[0]) * sinh(projection[1]) ;
            // const double temperature = sin(r_geometry.Center().X()) * sinh(r_geometry.Center().Y()) ;
            // const double temperature = projection[0]*projection[0]*projection[0] + projection[1]*projection[1]*projection[1] ;
            // const double temperature = projection[0]*projection[0]*projection[0]*projection[0] + projection[1]*projection[1]*projection[1]*projection[1] ;
            // const double u_D_scalar = sin(sqrt(2)*r_geometry.Center().X()) * sinh(r_geometry.Center().Y()) * cosh(r_geometry.Center().Z()) ;
            // const double u_D_scalar = r_geometry.Center().X() + r_geometry.Center().Y() + r_geometry.Center().Z() + 10;

            // FINAL ONE:
            const double u_D_scalar = candidateClosestSkinSegment1.GetGeometry()[closestNodeId].GetValue(TEMPERATURE);

            noalias(rRightHandSideVector) += H_sum_vec * u_D_scalar * penalty_integration;
            // Dirichlet BCs
            noalias(rRightHandSideVector) -= Guglielmo_innovation * DN_dot_n_vec * u_D_scalar * integration_points[0].Weight() ; // * std::abs(DetJ0) ;

            Vector temp(number_of_nodes);
            // RHS = ExtForces - K*temp;
            for (unsigned int i = 0; i < number_of_nodes; i++) {
                temp[i] = r_geometry[i].GetSolutionStepValue(TEMPERATURE);
            }
            // RHS -= K*temp
            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
            
        }

        // for (unsigned int i = 0; i < number_of_nodes; i++) {
        //     std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
        //     outputFile << r_geometry[i].GetId() << "  " << r_geometry[i].GetDof(TEMPERATURE).EquationId() <<"\n";
        //     outputFile.close();
        //     // KRATOS_WATCH(r_geometry[i].GetDof(TEMPERATURE).EquationId() )
        // }
        KRATOS_CATCH("")
    }

    //_________________________________________________________________________________________________________________________________________
    // NEUMANN CONDITION
    //_________________________________________________________________________________________________________________________________________

    void SBMLaplacianCondition::CalculateAllNeumann(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;

        // loopIdentifier is inner or outer
        std::string loopIdentifier = this->GetValue(IDENTIFIER);

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
        
        // Initialize DN_DX
        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        const unsigned int dim = DN_De[0].size2();
        Matrix DN_DX(number_of_nodes,dim);

        // Compute basis function order (Note: it is not allow to use different orders in different directions)
        if (dim == 3) {
            mbasisFunctionsOrder = std::cbrt(DN_De[0].size1()) - 1;
        } else {
            mbasisFunctionsOrder = std::sqrt(DN_De[0].size1()) - 1;
        }

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // Vector determinant_jacobian_vector(integration_points.size());
        // r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);  // = 1
        // // Initialize Jacobian
        // GeometryType::JacobiansType J0;
        // Matrix InvJ0(dim,dim);
        // r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        // double DetJ0;
        // Matrix Jacobian = ZeroMatrix(2,2);
        // Jacobian(0,0) = J0[point_number](0,0);
        // Jacobian(0,1) = J0[point_number](0,1);
        // Jacobian(1,0) = J0[point_number](1,0);
        // Jacobian(1,1) = J0[point_number](1,1);
        // MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            // Obtaining the projection from the closest skin segment
            Vector projection(3);
            projection[0] = candidateClosestSkinSegment1.GetGeometry()[0].X() ;
            projection[1] = candidateClosestSkinSegment1.GetGeometry()[0].Y() ;
            projection[2] = candidateClosestSkinSegment1.GetGeometry()[0].Z() ;

            Vector d(3);
            d[0] = projection[0] - r_geometry.Center().X();
            d[1] = projection[1] - r_geometry.Center().Y();
            d[2] = projection[2] - r_geometry.Center().Z();
            // d[0] = 0;
            // d[1] = 0;
            // d[2] = 0;

            const Matrix& N = r_geometry.ShapeFunctionsValues();

            noalias(DN_DX) = DN_De[0]; // prod(DN_De[point_number],InvJ0);

            // Compute the normal
            array_1d<double, 3> normal_parameter_space;
            r_geometry.Calculate(NORMAL, normal_parameter_space);

            // Neumann BC, the true normal is needed
            array_1d<double, 3> true_n;
            if (dim == 2) {
                // Need also the second closest condition in 2D
                Condition candidateClosestSkinSegment2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1] ;
                array_1d<double,3> vectorSkinSegment1 = candidateClosestSkinSegment1.GetGeometry()[1] - candidateClosestSkinSegment1.GetGeometry()[0];
                array_1d<double,3> vectorSkinSegment2 = candidateClosestSkinSegment2.GetGeometry()[1] - candidateClosestSkinSegment2.GetGeometry()[0];
                array_1d<double,3> vectorOutOfPlane = ZeroVector(3);
                vectorOutOfPlane[2] = 1.0;
                
                array_1d<double,3> crossProductSkinSegment1;
                array_1d<double,3> crossProductSkinSegment2; 
                MathUtils<double>::CrossProduct(crossProductSkinSegment1, vectorOutOfPlane, vectorSkinSegment1);
                MathUtils<double>::CrossProduct(crossProductSkinSegment2, vectorOutOfPlane, vectorSkinSegment2);
                
                true_n = crossProductSkinSegment1 / MathUtils<double>::Norm(crossProductSkinSegment1) + crossProductSkinSegment2 / MathUtils<double>::Norm(crossProductSkinSegment2);
                if (loopIdentifier == "inner") {
                    true_n = true_n / MathUtils<double>::Norm(true_n) ;
                } else { // outer
                    true_n = - true_n / MathUtils<double>::Norm(true_n) ;
                }
            }
            else {
                // 3D CASE
                array_1d<double,3> vectorSkinSegment1 = candidateClosestSkinSegment1.GetGeometry()[1] - candidateClosestSkinSegment1.GetGeometry()[0];
                array_1d<double,3> vectorSkinSegment2 = candidateClosestSkinSegment1.GetGeometry()[2] - candidateClosestSkinSegment1.GetGeometry()[1];
                MathUtils<double>::CrossProduct(true_n, vectorSkinSegment1, vectorSkinSegment2);

                if (loopIdentifier == "inner") {
                    true_n = true_n / MathUtils<double>::Norm(true_n) ;
                } else { // outer
                    true_n = - true_n / MathUtils<double>::Norm(true_n) ;
                }
            }

            // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
            std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
            outputFile << projection[0] << " " << projection[1] << " " << projection[2] << " "  << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() << " "
                                        << true_n[0] << " " << true_n[1]  << " " << true_n[2] << "\n";
            outputFile.close();

            // Compute all the derivatives of the basis functions involved
            std::vector<Matrix> nShapeFunctionDerivatives;
            for (int n = 1; n <= mbasisFunctionsOrder; n++) {
                nShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, point_number, this->GetIntegrationMethod()));
            }

            // Neumann (Taylor expansion of the gradient)
            Matrix H = ZeroMatrix(1, number_of_nodes);
            Matrix HgradX = ZeroMatrix(1, number_of_nodes);
            Matrix HgradY = ZeroMatrix(1, number_of_nodes);
            Matrix HgradZ = ZeroMatrix(1, number_of_nodes);

            Matrix DN_dot_n_tilde = ZeroMatrix(1, number_of_nodes);
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i) = N(point_number, i);
                // grad N cdot n_tilde
                for (IndexType idim = 0; idim < dim; idim++) {
                    DN_dot_n_tilde(0, i)  += DN_DX(i, idim) * normal_parameter_space[idim];
                } 
                // Reset for each control point
                double H_taylor_term_X = 0.0; 
                double H_taylor_term_Y = 0.0; 
                double H_taylor_term_Z = 0.0; 

                if (dim == 2) {
                    for (int n = 2; n <= mbasisFunctionsOrder; n++) {
                        // Retrieve the appropriate derivative for the term
                        Matrix& shapeFunctionDerivatives = nShapeFunctionDerivatives[n-1];
                        for (int k = 0; k <= n-1; k++) {
                            int n_k = n - 1 - k;
                            double derivative = shapeFunctionDerivatives(i,k); 
                            // Compute the Taylor term for this derivative
                            H_taylor_term_X += computeTaylorTerm(derivative, d[0], n_k, d[1], k);
                        }
                        for (int k = 0; k <= n-1; k++) {
                            int n_k = n - 1 - k;
                            double derivative = shapeFunctionDerivatives(i,k+1); 
                            // Compute the Taylor term for this derivative
                            H_taylor_term_Y += computeTaylorTerm(derivative, d[0], n_k, d[1], k);
                        }
                    }
                } else {
                    // 3D
                    for (int n = 2; n <= mbasisFunctionsOrder; n++) {
                        Matrix& shapeFunctionDerivatives = nShapeFunctionDerivatives[n-1];
                    
                        int countDerivativeId = 0;
                        // Loop over blocks of derivatives in x
                        for (int k_x = n; k_x >= 0; k_x--) {
                            // Loop over the possible derivatives in y
                            for (int k_y = n - k_x; k_y >= 0; k_y--) {
        
                                // derivatives in z
                                int k_z = n - k_x - k_y;
                                double derivative = shapeFunctionDerivatives(i,countDerivativeId); 
                                
                                if (k_x >= 1) {
                                    H_taylor_term_X += computeTaylorTerm3D(derivative, d[0], k_x-1, d[1], k_y, d[2], k_z);
                                }
                                if (k_y >= 1) {
                                    H_taylor_term_Y += computeTaylorTerm3D(derivative, d[0], k_x, d[1], k_y-1, d[2], k_z);
                                }
                                if (k_z >= 1) {
                                    H_taylor_term_Z += computeTaylorTerm3D(derivative, d[0], k_x, d[1], k_y, d[2], k_z-1);
                                }     
                                countDerivativeId++;
                            }
                        }
                        // H_taylor_term_X += shapeFunctionDerivatives(i,0)*d[0] + shapeFunctionDerivatives(i,1)*d[1] + shapeFunctionDerivatives(i,2)*d[2] ;
                        // H_taylor_term_Y += shapeFunctionDerivatives(i,1)*d[0] + shapeFunctionDerivatives(i,3)*d[1] + shapeFunctionDerivatives(i,4)*d[2] ;
                        // H_taylor_term_Z += shapeFunctionDerivatives(i,2)*d[0] + shapeFunctionDerivatives(i,4)*d[1] + shapeFunctionDerivatives(i,5)*d[2] ;
                    }
                }
                
                HgradX(0,i) = H_taylor_term_X + DN_DX(i, 0);
                HgradY(0,i) = H_taylor_term_Y + DN_DX(i, 1);
                HgradZ(0,i) = H_taylor_term_Z + DN_DX(i, 2);
            }    

            // dot product n cdot n_tilde
            double n_ntilde = true_n[0] * normal_parameter_space[0] + true_n[1] * normal_parameter_space[1] + true_n[2] * normal_parameter_space[2];

            // dot product grad cdot n
            Matrix HgradNdot_n = ZeroMatrix(1, number_of_nodes);
            HgradNdot_n = HgradX * true_n[0] + HgradY * true_n[1] + HgradZ * true_n[2];

            // compute Neumann contributions
            noalias(rLeftHandSideMatrix) += prod(trans(H), HgradNdot_n)  * n_ntilde   * integration_points[point_number].Weight(); // * std::abs(determinant_jacobian_vector[point_number]) ;
            noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n_tilde)            * integration_points[point_number].Weight() ; // * std::abs(DetJ0) ;

            if (CalculateResidualVectorFlag) {
                
                // FINAL ONE:
                // const double temperature = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(TEMPERATURE);
                
                Vector t_N(number_of_nodes);

                for (IndexType i = 0; i < number_of_nodes; ++i)
                {
                    // Fundamental for use Manufactured solution -> how to change this?
                    t_N[i] = cos(projection[0]) * sinh(projection[1]) * true_n[0] + sin(projection[0]) * cosh(projection[1])  * true_n[1] ;
                    
                    // 3D flux  
                    // true sol: sin(sqrt(2)*x)*sinh(y)*cosh(z)"
                    // t_N[i] = sqrt(2) * cos(sqrt(2)*projection[0]) * sinh(projection[1]) *  cosh(projection[2]) * true_n[0] + 
                    //                    sin(sqrt(2)*projection[0]) * cosh(projection[1]) *  cosh(projection[2]) * true_n[1] +
                    //                    sin(sqrt(2)*projection[0]) * sinh(projection[1]) *  sinh(projection[2]) * true_n[2];
                    // double x = r_geometry.Center().X(); double y = r_geometry.Center().Y(); double z = r_geometry.Center().Z();
                    // t_N[i] = sqrt(2) * cos(sqrt(2)*x) * sinh(y) *  cosh(z) * true_n[0] + 
                    //                    sin(sqrt(2)*x) * cosh(y) *  cosh(z) * true_n[1] +
                    //                    sin(sqrt(2)*x) * sinh(y) *  sinh(z) * true_n[2];

                    // true sol: x+y+z
                    // t_N[i] = true_n[0]+true_n[1]+true_n[2];

                    // t_N[i] = 2*projection[0] * true_n[0] + 2*projection[1]  * true_n[1] ;
                    // t_N[i] = 1 * true_n[0] -1  * true_n[1] ;

                    // double x = r_geometry.Center().X(); double y = r_geometry.Center().Y();
                    // t_N[i] = cos(x) * sinh(y) * true_n[0] + sin(x) * cosh(y)  * true_n[1] ;
                }
                // Neumann Contributions
                noalias(rRightHandSideVector) += prod(prod(trans(H), H), t_N) * n_ntilde * integration_points[point_number].Weight(); // * std::abs(determinant_jacobian_vector[point_number]);

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
    
    
    


    unsigned long long SBMLaplacianCondition::factorial(int n) 
    {
        if (n == 0) return 1;
        unsigned long long result = 1;
        for (int i = 2; i <= n; ++i) result *= i;
        return result;
    }

    // Function to compute a single term in the Taylor expansion
    double SBMLaplacianCondition::computeTaylorTerm(double derivative, double dx, int n_k, double dy, int k)
    {
        return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (factorial(k) * factorial(n_k));    
    }

    double SBMLaplacianCondition::computeTaylorTerm3D(double derivative, double dx, int k_x, double dy, int k_y, double dz, int k_z)
    {   
        return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (factorial(k_x) * factorial(k_y) * factorial(k_z));    
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


    void SBMLaplacianCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        double sol_on_surrogate = 0.0;

        for (IndexType i = 0; i < r_geometry.size(); ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value = r_geometry[i].GetSolutionStepValue(TEMPERATURE);
            sol_on_surrogate += r_N(0, i) * output_solution_step_value;
        }      

        double rOutput = 0;
        std::vector<Vector> integration_point_list_on_true_boundary  = this->GetValue(INTEGRATION_POINTS);
        std::vector<double> integration_weight_list_on_true_boundary = this->GetValue(INTEGRATION_WEIGHTS);
        
        std::ofstream output_file("txt_files/results_on_true_boundary.txt", std::ios::app);
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        // Loop over the number of skin boundary integration points
        for (int i_gauss = 0; i_gauss < integration_weight_list_on_true_boundary.size(); i_gauss++){
            // For each of them we compute the solution at the true boundary reversing the Taylor expansion
            Vector curr_gauss_point_on_true = integration_point_list_on_true_boundary[i_gauss];
            Vector d = curr_gauss_point_on_true - r_geometry.Center();

            double solution_on_true = 0.0; 

            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                double output_solution_step_value = r_geometry[i].GetSolutionStepValue(TEMPERATURE);

                double H_taylor_term = 0.0; // Reset for each node
                for (int n = 1; n <= mbasisFunctionsOrder; n++) {
                    // Retrieve the appropriate derivative for the term
                    Matrix& shapeFunctionDerivatives = mShapeFunctionDerivatives[n-1];
                    for (int k = 0; k <= n; k++) {
                        int n_k = n - k;
                        double derivative = shapeFunctionDerivatives(i,k); 
                        // Compute the Taylor term for this derivative
                        H_taylor_term += computeTaylorTerm(derivative, d[0], n_k, d[1], k);
                    }
                }

                solution_on_true += (r_N(0,i) + H_taylor_term) * output_solution_step_value;
            } 
            output_file << solution_on_true << " "<< integration_weight_list_on_true_boundary[i_gauss] << " " 
                        << curr_gauss_point_on_true[0] << " " << curr_gauss_point_on_true[1]  << " " << curr_gauss_point_on_true[2]  << " "
                        << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() 
                        << std::endl;
        }
        output_file.close();

        std::ofstream output_file2("txt_files/boundary_GPs_3D.txt", std::ios::app);
        if (output_file2.is_open()) {
            output_file2 << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file2 << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z()  << std::endl;
            output_file2.close();
        }

        // std::ofstream output_file("txt_files/boundary_results_GPs.txt", std::ios::app);
        // if (output_file.is_open()) {
        //     output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        //     output_file << rOutput << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << integration_points[0].Weight()*factor << " " << d[0] << " " << d[1] << std::endl;
        //     output_file.close();
        // }
    }


} // Namespace Kratos
