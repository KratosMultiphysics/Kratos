// //    |  /           |
// //    ' /   __| _` | __|  _ \   __|
// //    . \  |   (   | |   (   |\__ `
// //   _|\_\_|  \__,_|\__|\___/ ____/
// //                   Multi-Physics
// //
// //  License:         BSD License
// //                   Kratos default license: kratos/license.txt
// //
// //  Main authors:    Andea Gorgi
// //                  
// //

// // System includes

// // External includes

// // Project includes
// #include "custom_conditions/sbm_contact_2D_condition.h"
// #include "includes/global_pointer_variables.h"

// namespace Kratos
// {

//     void SbmContact2DCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
//     {
//         InitializeMaterial();


//         //****************************** */
//         // INITIALIZE BM PARAMETES
//         //******************************** */
//         const auto& r_geometry_master = GetMasterGeometry();
//         const auto& r_geometry_slave = GetSlaveGeometry();

//         NodeType projection_node_master;         NodeType projection_node_slave;
//         //
//         projection_node_master = r_geometry_master.GetValue(NEIGHBOUR_NODES)[0];
//         projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];

//         mDistanceMaster[0] = projection_node_master.X() - r_geometry_master.Center().X();   mDistanceSlave[0] = projection_node_slave.X() - r_geometry_slave.Center().X();             
//         mDistanceMaster[1] = projection_node_master.Y() - r_geometry_master.Center().Y();   mDistanceSlave[1] = projection_node_slave.Y() - r_geometry_slave.Center().Y();
        
//         Vector normal_true_master = projection_node_master.GetValue(NORMAL);
//         Vector normal_true_slave = projection_node_slave.GetValue(NORMAL);
        
//         this->SetValue(NORMAL_MASTER, normal_true_master);
//         this->SetValue(NORMAL_SLAVE, normal_true_slave);

//         Vector true_normal_mean = normal_true_master*0.5 - normal_true_slave*0.5;

//         this->SetValue(NORMAL, true_normal_mean);
        
//         // // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
//         // std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::app);
//         // outputFile << projection_node_master.X() << " " << projection_node_master.Y() << " " 
//         //            << projection_node_slave.X() << " " << projection_node_slave.Y() <<"\n";


//         std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::app);
//         outputFile << projection_node_master.X() << " " << projection_node_master.Y() << " " 
//                    << r_geometry_master.Center().X() << " " << r_geometry_master.Center().Y() <<"\n";
//         outputFile.close();
//     }


//     void SbmContact2DCondition::InitializeMaterial()
//     {
//         KRATOS_TRY
//         if ((*mpPropMaster)[CONSTITUTIVE_LAW] != nullptr ) {
//             const GeometryType& r_geometry_master = GetMasterGeometry();

//             const auto& N_values = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());

//             mpConstitutiveLawMaster = (*mpPropMaster)[CONSTITUTIVE_LAW]->Clone();
//             mpConstitutiveLawMaster->InitializeMaterial( *mpPropMaster, r_geometry_master, row(N_values , 0 ));

//         } else
//             KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;


//         if ( (*mpPropSlave)[CONSTITUTIVE_LAW] != nullptr ) {
//             const GeometryType& r_geometry_slave = GetSlaveGeometry();

//             const auto& N_values = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

//             mpConstitutiveLawSlave = (*mpPropSlave)[CONSTITUTIVE_LAW]->Clone();
//             mpConstitutiveLawSlave->InitializeMaterial( *mpPropSlave, r_geometry_slave, row(N_values , 0 ));

//         } else
//             KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

//         KRATOS_CATCH( "" );

//     }
//     void SbmContact2DCondition::CalculateAll(
//         MatrixType& rLeftHandSideMatrix,
//         VectorType& rRightHandSideVector,
//         const ProcessInfo& rCurrentProcessInfo,
//         const bool CalculateStiffnessMatrixFlag,
//         const bool CalculateResidualVectorFlag
//     )
//     {
//         KRATOS_TRY
//         int bc_type = (*mpPropMaster)[ACTIVATION_LEVEL];

//         // INITIALIZE AND RESIZE 

//         double penalty_master = 1e-4;// 50;// (*mpPropMaster)[PENALTY_FACTOR];
//         double penalty_slave =  1e-4;// 50;//(*mpPropSlave)[PENALTY_FACTOR];

//         Vector penalty_vector = ZeroVector(2);
//         penalty_vector[0] = penalty_master; penalty_vector[1] = penalty_slave;

//         this->SetValue(PENALTY, penalty_vector);
//         const SizeType dim = 2;

//         const auto& r_geometry_master = GetMasterGeometry();
//         const SizeType number_of_nodes_master = r_geometry_master.size();

//         const auto& r_geometry_slave = GetSlaveGeometry();
//         const SizeType number_of_nodes_slave = r_geometry_slave.size();

//         const SizeType mat_size_1 = (number_of_nodes_master+number_of_nodes_slave) * dim;
//         const SizeType mat_size_2 = (number_of_nodes_master+number_of_nodes_slave) * dim;
//         //resizing as needed the LHS
//         if(rLeftHandSideMatrix.size1() != mat_size_1 || rLeftHandSideMatrix.size2() != mat_size_2 )
//             rLeftHandSideMatrix.resize(mat_size_1,mat_size_2,false);

//         noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size_1,mat_size_2); //resetting LHS
        
//         // resizing as needed the RHS
//         if(rRightHandSideVector.size() != mat_size_1)
//             rRightHandSideVector.resize(mat_size_1,false);
//         noalias(rRightHandSideVector) = ZeroVector(mat_size_1); //resetting RHS

//         // Integration
//         const GeometryType::IntegrationPointsArrayType& integration_points_master = r_geometry_master.IntegrationPoints();

//         // Shape function derivatives 
//         // MASTER ------------------------------------------------
//         // Initialize Jacobian
//         GeometryType::JacobiansType J0_master;
//         // Initialize DN_DX
//         Matrix DN_DX_master(number_of_nodes_master,2);
//         Matrix InvJ0_master(dim,dim);

//         const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
//         r_geometry_master.Jacobian(J0_master,this->GetIntegrationMethod());

//         double DetJ0_master;

//         Matrix Jacobian = ZeroMatrix(2,2);
//         Jacobian(0,0) = J0_master[0](0,0); Jacobian(0,1) = J0_master[0](0,1);
//         Jacobian(1,0) = J0_master[0](1,0); Jacobian(1,1) = J0_master[0](1,1);

//         // // Calculating inverse jacobian and jacobian determinant
//         MathUtils<double>::InvertMatrix(Jacobian,InvJ0_master,DetJ0_master);

//         KRATOS_ERROR_IF(std::abs(DetJ0_master -1.0) > 1e-11) << ":::[SbmContact2DCondition]::: Error. Master: Not coincident physical and parameter spaces" 
//                         << "Jacobian:   \n" << Jacobian << std::endl;

//         const Matrix& N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());

//         // // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
//         noalias(DN_DX_master) = DN_De_master[0];

//         // // SLAVE ------------------------------------------------
//         // Initialize Jacobian
//         GeometryType::JacobiansType J0_slave;
//         // Initialize DN_DX
//         Matrix DN_DX_slave(number_of_nodes_slave,2);
//         Matrix InvJ0_slave(dim,dim);

//         const GeometryType::ShapeFunctionsGradientsType& DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
//         r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        
//         double DetJ0_slave;
//         Matrix Jacobian_slave = ZeroMatrix(2,2);
//         Jacobian_slave(0,0) = J0_slave[0](0,0); Jacobian_slave(0,1) = J0_slave[0](0,1);
//         Jacobian_slave(1,0) = J0_slave[0](1,0); Jacobian_slave(1,1) = J0_slave[0](1,1);

//         // Calculating inverse jacobian and jacobian determinant
//         MathUtils<double>::InvertMatrix(Jacobian_slave,InvJ0_slave,DetJ0_slave);  //DetJ0 is not the true one for NURBS geometries

//         KRATOS_ERROR_IF(std::abs(DetJ0_master -1.0) > 1e-11) << ":::[SbmContact2DCondition]::: Error. Slave: Not coincident physical and parameter spaces" 
//                         << "detJ" << DetJ0_slave << std::endl;

//         const Matrix& N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

//         // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
//         // noalias(DN_DX_slave) = prod(DN_De_slave[0],InvJ0_slave);
//         noalias(DN_DX_slave) = DN_De_slave[0];

//         r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());


//         // //------------------- REFERENCE WEIGHT

//         const double thickness_master = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;
//         const double integration_weight_master = r_geometry_master.IntegrationPoints()[0].Weight() * thickness_master;

//         const double thickness_slave = (*mpPropSlave).Has(THICKNESS) ? (*mpPropSlave)[THICKNESS] : 1.0;
//         double integration_weight_slave = r_geometry_slave.IntegrationPoints()[0].Weight() * thickness_slave;

//         const Vector normal_surrogate_master_3D = mNormalMaster; 

//         Vector normal_surrogate_master_2D(2); normal_surrogate_master_2D[0] = normal_surrogate_master_3D[0]; normal_surrogate_master_2D[1] = normal_surrogate_master_3D[1];

//         const Vector normal_surrogate_slave_3D = mNormalSlave;

//         Vector normal_surrogate_slave_2D(2); normal_surrogate_slave_2D[0] = normal_surrogate_slave_3D[0]; normal_surrogate_slave_2D[1] = normal_surrogate_slave_3D[1];

//         // //------------------- DB MATRICES 

//         // // MODIFIED
//         Matrix B_master = ZeroMatrix(3,number_of_nodes_master);

//         Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);

//         CalculateB(B_master, DN_DX_master, number_of_nodes_master);

//         CalculateB(B_slave, DN_DX_slave, number_of_nodes_slave);

//         // //---------- GET CONSTITUTIVE MATRICES  
//         const Matrix r_D_master = GetValue(CONSTITUTIVE_MATRIX_MASTER);

//         const Matrix r_D_slave  = GetValue(CONSTITUTIVE_MATRIX_SLAVE);

//         //  ----------------------------------------------------------------------------
//         //  ----------------------------------------------------------------------------
//         //  Shfited Boundary Method:: Retrieve infos
//         //  ----------------------------------------------------------------------------
//         //  ----------------------------------------------------------------------------
//         array_1d<double, 3> normal_true_master;  array_1d<double, 3> normal_true_slave;
//         array_1d<double, 3> tau_true_master;     array_1d<double, 3> tau_true_slave;
//         double n_ntilde_master;                  double n_ntilde_slave;

//         auto& projection_node_master = (r_geometry_master.GetValue(NEIGHBOUR_NODES)[0]);
//         auto& projection_node_slave  = (projection_node_master.GetValue(NEIGHBOUR_NODES)[0]);

//         //
//         // projection_node_master = r_geometry_master.GetValue(NEIGHBOUR_NODES)[0];
//         normal_true_master = projection_node_master.GetValue(NORMAL);
//         tau_true_master = projection_node_master.GetValue(LOCAL_TANGENT);

//         if (std::abs(inner_prod(normal_surrogate_master_3D, normal_true_master)) < 1e-8)
//         {
//             // normal_true_master = normal_surrogate_master_3D;
//             return;
//         }

//         this->SetValue(NORMAL_MASTER, normal_true_master);

//         // projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];
//         normal_true_slave = projection_node_slave.GetValue(NORMAL);
//         tau_true_slave = projection_node_slave.GetValue(LOCAL_TANGENT);

//         this->SetValue(NORMAL_SLAVE, normal_true_slave);
//         if (std::abs(inner_prod(normal_surrogate_slave_3D, normal_true_slave)) < 1e-8) 
//         {
//             // normal_true_slave = normal_surrogate_slave_3D;
//             return;
//         }

//         double curvature_true_master = projection_node_master.GetValue(CURVATURE);
//         double curvature_true_slave = projection_node_slave.GetValue(CURVATURE);

//         double curvature_sign_master = inner_prod(normal_true_master, mDistanceMaster) > 0 ? +1.0 : -1.0;
//         double curvature_sign_slave = inner_prod(normal_true_slave, mDistanceSlave) > 0 ? +1.0 : -1.0;

//         double corrective_factor_curvature_on_slave = 1.0;
//         double corrective_factor_curvature_on_master = 1.0;

//         n_ntilde_master = inner_prod(normal_surrogate_master_3D, normal_true_master);
//         n_ntilde_slave = inner_prod(normal_surrogate_slave_3D, normal_true_slave);

//         // KRATOS_WATCH(projection_node_master)
//         // KRATOS_WATCH(r_geometry_master.Center())
//         // KRATOS_WATCH(normal_surrogate_master_3D)
//         // KRATOS_WATCH(normal_true_master)
//         // KRATOS_WATCH("---------------")
        

//         KRATOS_ERROR_IF(n_ntilde_master <= 0) << ":::[SbmContact2DCondition]::: Error. Master: Negative n_ntilde_master" 
//                     << r_geometry_master.Center() << normal_surrogate_master_3D << normal_true_master <<  std::endl;
//         KRATOS_ERROR_IF(n_ntilde_slave <= 0) << ":::[SbmContact2DCondition]::: Error. Slave: Negative n_ntilde_slave" 
//                     << r_geometry_slave.Center() << normal_surrogate_slave_3D << normal_true_slave <<  std::endl;

//         Matrix H_master = ZeroMatrix(1, number_of_nodes_master);       Matrix H_slave = ZeroMatrix(1, number_of_nodes_slave);
//         Matrix H_sum_master = ZeroMatrix(1, number_of_nodes_master);   Matrix H_sum_slave = ZeroMatrix(1, number_of_nodes_slave);

//         // Compute all the derivatives of the basis functions involved
//         std::vector<Matrix> n_shape_function_derivatives_master; std::vector<Matrix> n_shape_function_derivatives_slave;
//         int basis_functions_order_master = std::sqrt(DN_De_master[0].size1()) - 1;
//         int basis_functions_order_slave = std::sqrt(DN_De_slave[0].size1()) - 1;

//         //TODO:
//         // double x = r_geometry_master.Center().X();
//         // double y = r_geometry_master.Center().Y();
//         // if ((x >= 0.3 && x <= 0.7))
//         //     basis_functions_order_master = 1;
//         //     basis_functions_order_slave = 1;


//         for (int n = 1; n <= basis_functions_order_master; n++) {
//             n_shape_function_derivatives_master.push_back(r_geometry_master.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
//         }
//         for (int n = 1; n <= basis_functions_order_slave; n++) {
//             n_shape_function_derivatives_slave.push_back(r_geometry_slave.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
//         }


//         // Neumann (Taylor expansion of the gradient)
//         Matrix H_grad_master = ZeroMatrix(number_of_nodes_master, 2);
//         Matrix H_extension_master = ZeroMatrix(number_of_nodes_master, 2);
//         for (IndexType i = 0; i < number_of_nodes_master; ++i)
//         {
//             H_master(0, i) = N_master(0, i);
//             double H_taylor_term_X = 0.0; // Reset for each node
//             double H_taylor_term_Y = 0.0; 
//             for (int n = 2; n <= basis_functions_order_master; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = n_shape_function_derivatives_master[n-1];
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceMaster[0], n_k, mDistanceMaster[1], k);
//                 }
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k+1); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceMaster[0], n_k, mDistanceMaster[1], k);
//                 }
//             }
//             H_grad_master(i,0) = DN_DX_master(i,0) + H_taylor_term_X;
//             H_grad_master(i,1) = DN_DX_master(i,1) + H_taylor_term_Y;

//             H_extension_master(i,0) = H_taylor_term_X;
//             H_extension_master(i,1) = H_taylor_term_Y;
//         }                                                     

//         Matrix B_sum_master = ZeroMatrix(3, number_of_nodes_master);
//         Matrix B_extension_master = ZeroMatrix(3, number_of_nodes_master);

//         CalculateB(B_sum_master, H_grad_master, number_of_nodes_master);
//         CalculateB(B_extension_master, H_extension_master, number_of_nodes_master);

//         // slave
//         Matrix H_grad_slave = ZeroMatrix(number_of_nodes_slave, 2);
//         Matrix H_extension_slave = ZeroMatrix(number_of_nodes_slave, 2);
//         for (IndexType i = 0; i < number_of_nodes_slave; ++i)
//         {
//             H_slave(0, i) = N_slave(0, i);
//             double H_taylor_term_X = 0.0; // Reset for each node
//             double H_taylor_term_Y = 0.0; 
//             for (int n = 2; n <= basis_functions_order_slave; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = n_shape_function_derivatives_slave[n-1];
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceSlave[0], n_k, mDistanceSlave[1], k);
//                 }
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k+1); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceSlave[0], n_k, mDistanceSlave[1], k);
//                 }
//             }
//             H_grad_slave(i,0) = DN_DX_slave(i,0) + H_taylor_term_X;
//             H_grad_slave(i,1) = DN_DX_slave(i,1) + H_taylor_term_Y;

//             H_extension_slave(i,0) = H_taylor_term_X;
//             H_extension_slave(i,1) = H_taylor_term_Y;
//         }                                                     

//         Matrix B_sum_slave = ZeroMatrix(3, number_of_nodes_slave);
//         Matrix B_extension_slave = ZeroMatrix(3, number_of_nodes_slave);

//         CalculateB(B_sum_slave, H_grad_slave, number_of_nodes_slave);
//         CalculateB(B_extension_slave, H_extension_slave, number_of_nodes_slave);

//         // COMPUTE H_sum
//         // master
//         for (IndexType i = 0; i < number_of_nodes_master; ++i)
//         {
//             double H_taylor_term = 0.0; // Reset for each node
//             for (int n = 1; n <= basis_functions_order_master; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = n_shape_function_derivatives_master[n-1];
//                 for (int k = 0; k <= n; k++) {
//                     int n_k = n - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term += ComputeTaylorTerm(derivative, mDistanceMaster[0], n_k, mDistanceMaster[1], k);
//                 }
//             }
//             H_sum_master(0,i) = H_taylor_term + N_master(0, i);
//         }
//         // slave
//         for (IndexType i = 0; i < number_of_nodes_slave; ++i)
//         {
//             double H_taylor_term = 0.0; // Reset for each node
//             for (int n = 1; n <= basis_functions_order_slave; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = n_shape_function_derivatives_slave[n-1];
//                 for (int k = 0; k <= n; k++) {
//                     int n_k = n - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term += ComputeTaylorTerm(derivative, mDistanceSlave[0], n_k, mDistanceSlave[1], k);
//                 }
//             }
//             H_sum_slave(0,i) = H_taylor_term + N_slave(0, i);
//         }
                                                                                                                
//         Vector meshSize_uv = this->GetValue(MARKER_MESHES);
//         double h = std::min(meshSize_uv[0], meshSize_uv[1]);

//         // Guglielmo innovaction
//         double theta = -1.0;//-1.0;  // = 1 -> Penalty approach
//                                                 // = -1 -> Free-penalty approach
//         // if (penalty == -1.0) {
//         //     penalty = 0.0;
//         //     theta = -1.0;
//         // }

//         // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
//         penalty_master = penalty_master * h/ basis_functions_order_master /basis_functions_order_master;
//         penalty_slave = penalty_slave * h/ basis_functions_order_slave /basis_functions_order_slave;

//         // //-----------
//         // // COMPUTE NORMAL GAP AND CHECK CONDITION
//         Vector initial_gap_vector = projection_node_master - projection_node_slave;
//         const double initial_normal_gap_true_master = inner_prod(projection_node_slave - projection_node_master, normal_true_master);
//         const double initial_normal_gap_true_slave = inner_prod(projection_node_master - projection_node_slave, normal_true_slave);
        
//         double curvature_sign_gap = inner_prod(normal_true_slave, initial_gap_vector) > 0 ? +1.0 : -1.0;
        
//         // // Assembly
//         Matrix DB_master = prod(r_D_master,B_master);
//         Matrix DB_slave = prod(r_D_slave,B_slave);
//         Matrix DB_sum_master = prod(r_D_master,B_sum_master);
//         Matrix DB_sum_slave = prod(r_D_slave,B_sum_slave);

//         Matrix DB_extension_master = prod(r_D_master,B_extension_master);
//         Matrix DB_extension_slave = prod(r_D_slave,B_extension_slave);

//         // //***************
//         // FLUX TERMS  */
//         // Flux Term on Master
        
//         // if (r_geometry_master.GetValue(NEIGHBOUR_NODES)(0)->GetValue(IS_ASSEMBLED))
//         // {
//             // KRATOS_WATCH(r_geometry_master.Center())
//             // KRATOS_WATCH(projection_node_master)
//             // KRATOS_WATCH(projection_node_slave)
//             // KRATOS_WATCH(r_geometry_slave.Center());
//         //     return;
//         // }

//         double jacobian_determinant = 1/n_ntilde_slave*n_ntilde_master * std::abs(inner_prod(normal_true_slave, normal_true_master));
//         integration_weight_slave = jacobian_determinant*integration_weight_master;
        
//         // n_ntilde_master = inner_prod(normal_surrogate_master_3D, normal_true_master)/(1+curvature_sign_master*curvature_true_master*norm_2(mDistanceMaster));
//         // n_ntilde_slave = inner_prod(normal_surrogate_slave_3D, normal_true_slave)/(1+curvature_sign_slave*curvature_true_slave*norm_2(mDistanceSlave));
//         corrective_factor_curvature_on_slave = 1/(1+curvature_sign_master*curvature_true_master*norm_2(mDistanceMaster))
//                                             * (1+curvature_sign_slave*curvature_true_slave*norm_2(mDistanceSlave));
//         corrective_factor_curvature_on_master = 1/(1+curvature_sign_slave*curvature_true_slave*norm_2(mDistanceSlave))
//                                             * (1+curvature_sign_master*curvature_true_master*norm_2(mDistanceMaster));

//         integration_weight_slave *= (1+curvature_sign_slave*curvature_true_slave*norm_2(mDistanceSlave));
//         integration_weight_slave /= (1+curvature_sign_gap*curvature_true_slave*norm_2(initial_gap_vector));
//         integration_weight_slave /= (1+curvature_sign_master*curvature_true_master*norm_2(mDistanceMaster));
//         this->SetValue(INTEGRATION_WEIGHT, integration_weight_master);
//         this->SetValue(GAMMA_CONTACT, integration_weight_slave);
//         // KRATOS_WATCH(integration_weight_slave)
//         // KRATOS_WATCH(integration_weight_master)

//         // -[sigma_1(u) \dot n_tilde] \dot * w_1)
//         // if (!(this->GetValue(ACTIVATION_LEVEL) != 1 && this->GetValue(ACTIVATION_LEVEL) != 3))
//         // if (this->GetValue(ACTIVATION_LEVEL) == 0)
//         // {
        
//         for (IndexType i = 0; i < number_of_nodes_master; i++)
//         {
//             for (IndexType j = 0; j < number_of_nodes_master; j++)
//             {
//                 for (IndexType idim = 0; idim < 2; idim++) {
//                     const int iglob = 2*i+idim;
//                     for (IndexType jdim = 0; jdim < 2; jdim++) {
//                         const int jglob = 2*j+jdim;
//                         if (bc_type == 0)
//                         {
//                             Vector sigma_u_n(2);
//                             sigma_u_n[0] = (DB_master(0, jglob)* normal_surrogate_master_2D[0] + DB_master(2, jglob)* normal_surrogate_master_2D[1]);
//                             sigma_u_n[1] = (DB_master(2, jglob)* normal_surrogate_master_2D[0] + DB_master(1, jglob)* normal_surrogate_master_2D[1]);

//                             rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_u_n[idim] * integration_weight_master;

//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* normal_true_master[0] + DB_sum_master(2, jglob)* normal_true_master[1]);
//                             extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* normal_true_master[0] + DB_sum_master(1, jglob)* normal_true_master[1]);

//                             rLeftHandSideMatrix(iglob, jglob) += H_master(0,i)*extension_sigma_u_n[idim] * n_ntilde_master * integration_weight_master;
//                         }
//                         else if (bc_type == 1)
//                         { // NEW VERSION 

//                             const double tau_n_tilde_master = inner_prod(tau_true_master, normal_surrogate_master_2D);
                        
//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_extension_master(0, jglob)* normal_true_master[0] + DB_extension_master(2, jglob)* normal_true_master[1]);
//                             extension_sigma_u_n[1] = (DB_extension_master(2, jglob)* normal_true_master[0] + DB_extension_master(1, jglob)* normal_true_master[1]);

//                             rLeftHandSideMatrix(iglob, jglob) += H_master(0,i)*extension_sigma_u_n[idim] * n_ntilde_master * integration_weight_master;




//                             // flux term
//                             Vector sigma_u_t(2);
//                             sigma_u_t[0] = (DB_master(0, jglob)* tau_true_master[0] + DB_master(2, jglob)* tau_true_master[1]);
//                             sigma_u_t[1] = (DB_master(2, jglob)* tau_true_master[0] + DB_master(1, jglob)* tau_true_master[1]);

//                             // rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_u_t[idim] * tau_n_tilde_master * integration_weight_master;

//                             const double sigma_u_t_t =  sigma_u_t[0]*tau_true_master[0] + sigma_u_t[1]*tau_true_master[1];

//                             // const double sigma_u_t_n =  sigma_u_t[0]*normal_true_master[0] + sigma_u_t[1]*normal_true_master[1];

//                             const double extension_sigma_u_n_t = extension_sigma_u_n[0]*tau_true_master[0] + extension_sigma_u_n[1]*tau_true_master[1];

//                             rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_u_t_t * tau_true_master[idim] * tau_n_tilde_master * integration_weight_master;

//                             //  rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_u_t_n * normal_true_master[idim] * tau_n_tilde_master * integration_weight_master;
//                             rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*extension_sigma_u_n_t * normal_true_master[idim] * tau_n_tilde_master *integration_weight_master;
//                         }
//                     }
//                 }
//             }
//         }
//         // // KRATOS_WATCH(r_geometry_master.Center())
//         // // KRATOS_WATCH(projection_node_master)
//         // // KRATOS_WATCH(projection_node_slave)
//         // // KRATOS_WATCH(normal_surrogate_master_3D)
//         // // KRATOS_WATCH(normal_true_master)
//         // // KRATOS_WATCH("-------------------")
//         // }
//         // // // // else
//         // // // // {KRATOS_WATCH(r_geometry_slave.Center());
//         // // // // exit(0);
//         // // // // }
        
//         // // // // Flux Term on Slave

//         // // // // FIXME:

        
//         // 
//         // if (this->GetValue(ACTIVATION_LEVEL) != 0)
//         for (IndexType i = 0; i < number_of_nodes_slave; i++)
//             for (IndexType idim = 0; idim < 2; idim++) 
//             {
//                 const int iglob = 2*i+idim + 2*number_of_nodes_master;
//                 for (IndexType j = 0; j < number_of_nodes_slave; j++)
//                 {
//                     for (IndexType jdim = 0; jdim < 2; jdim++) {
//                         const int jglob = 2*j+jdim + 2*number_of_nodes_master;
//                         const int j_index = 2*j+jdim;

//                         if (bc_type == 0)
//                         {
//                             Vector sigma_u_n(2);
//                             sigma_u_n[0] = (DB_slave(0, j_index)* normal_surrogate_slave_2D[0] + DB_slave(2, j_index)* normal_surrogate_slave_2D[1]);
//                             sigma_u_n[1] = (DB_slave(2, j_index)* normal_surrogate_slave_2D[0] + DB_slave(1, j_index)* normal_surrogate_slave_2D[1]);
                            
//                             double first_term = -H_slave(0,i)*sigma_u_n[idim] * integration_weight_slave;

//                             rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*sigma_u_n[idim] *integration_weight_slave;

//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* normal_true_slave[0] + DB_sum_slave(2, j_index)* normal_true_slave[1]);
//                             extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* normal_true_slave[0] + DB_sum_slave(1, j_index)* normal_true_slave[1]);

//                             rLeftHandSideMatrix(iglob, jglob) += H_slave(0,i)*extension_sigma_u_n[idim] * n_ntilde_slave * integration_weight_slave;

//                             double second_term = H_slave(0,i)*extension_sigma_u_n[idim] * n_ntilde_slave *integration_weight_slave;

//                         }
//                         else if (bc_type == 1)
//                         { // NEW VERSION 

//                             const double tau_n_tilde_slave = inner_prod(tau_true_slave, normal_surrogate_slave_2D);
                        
//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_extension_slave(0, j_index)* normal_true_slave[0] + DB_extension_slave(2, j_index)* normal_true_slave[1]);
//                             extension_sigma_u_n[1] = (DB_extension_slave(2, j_index)* normal_true_slave[0] + DB_extension_slave(1, j_index)* normal_true_slave[1]);

//                             rLeftHandSideMatrix(iglob, jglob) += H_slave(0,i)*extension_sigma_u_n[idim] * n_ntilde_slave * integration_weight_slave;




//                             // flux term
//                             Vector sigma_u_t(2);
//                             sigma_u_t[0] = (DB_slave(0, j_index)* tau_true_slave[0] + DB_slave(2, j_index)* tau_true_slave[1]);
//                             sigma_u_t[1] = (DB_slave(2, j_index)* tau_true_slave[0] + DB_slave(1, j_index)* tau_true_slave[1]);

//                             // rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*sigma_u_t[idim] * tau_n_tilde_slave * integration_weight_slave;

//                             const double sigma_u_t_t =  sigma_u_t[0]*tau_true_slave[0] + sigma_u_t[1]*tau_true_slave[1];

//                             // const double sigma_u_t_n =  sigma_u_t[0]*normal_true_master[0] + sigma_u_t[1]*normal_true_master[1];

//                             const double extension_sigma_u_n_t = extension_sigma_u_n[0]*tau_true_slave[0] + extension_sigma_u_n[1]*tau_true_slave[1];

//                             rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*sigma_u_t_t * tau_true_slave[idim] * tau_n_tilde_slave * integration_weight_slave;

//                             //  rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_u_t_n * normal_true_master[idim] * tau_n_tilde_master * integration_weight_slave;
//                             rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*extension_sigma_u_n_t * normal_true_slave[idim] * tau_n_tilde_slave *integration_weight_slave;
//                         }
//                     }
//                 }
                
//             }
//         n_ntilde_master *= 0;
//         n_ntilde_slave *= 2;
//         // exit(0);
//         // //***************
//         // NITSCHE CONSTANT TERMS  */

//         int active_constant_nitsche_terms = 1;
//         // Nitsche Constant Term on Master

//         // -0.5*[theta*gamma_1*(n_ntilde_1) (E(sigma_1(u)) \dot n_1) \dot n_1) *  (E(sigma_1(v)) \dot n_1) \dot n_1)] 
//         // if (this->GetValue(ACTIVATION_LEVEL) != 0)
//         for (IndexType i = 0; i < number_of_nodes_master; i++)
//             for (IndexType idim = 0; idim < 2; idim++) {
//                 const int iglob = 2*i+idim;
//                 Vector extension_sigma_w_n(2);
//                 extension_sigma_w_n[0] = (DB_sum_master(0, iglob)* normal_true_master[0] + DB_sum_master(2, iglob)* normal_true_master[1]);
//                 extension_sigma_w_n[1] = (DB_sum_master(2, iglob)* normal_true_master[0] + DB_sum_master(1, iglob)* normal_true_master[1]);
//                 double extension_sigma_w_n_n = extension_sigma_w_n[0]*normal_true_master[0] + extension_sigma_w_n[1]*normal_true_master[1];

//                 Vector sigma_w_n(2);
//                 sigma_w_n[0] = (DB_master(0, iglob)* normal_true_master[0] + DB_master(2, iglob)* normal_true_master[1]);
//                 sigma_w_n[1] = (DB_master(2, iglob)* normal_true_master[0] + DB_master(1, iglob)* normal_true_master[1]);
//                 double sigma_w_n_n = sigma_w_n[0]*normal_true_master[0] + sigma_w_n[1]*normal_true_master[1];

//                 for (IndexType j = 0; j < number_of_nodes_master; j++) {
//                     for (IndexType jdim = 0; jdim < 2; jdim++) {
//                         const int jglob = 2*j+jdim;
                        
//                         Vector extension_sigma_u_n(2);
//                         extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* normal_true_master[0] + DB_sum_master(2, jglob)* normal_true_master[1]);
//                         extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* normal_true_master[0] + DB_sum_master(1, jglob)* normal_true_master[1]);
//                         double extension_sigma_u_n_n = extension_sigma_u_n[0]*normal_true_master[0] + extension_sigma_u_n[1]*normal_true_master[1];

//                         rLeftHandSideMatrix(iglob, jglob) -= active_constant_nitsche_terms*0.5*theta*penalty_master * n_ntilde_master 
//                                                             * extension_sigma_u_n_n * sigma_w_n_n
//                                                             * integration_weight_master;
                        
//                     }
//                 }
//             }

//         // Nitsche Constant Term on Slave

//         // -0.5*[theta*gamma_2*(n_ntilde_2) (E(sigma_2(u)) \dot n_2) \dot n_2) *  (E(sigma_2(v)) \dot n_2) \dot n_2)] 
//         // if (!projection_node_slave.GetValue(NEIGHBOUR_GEOMETRY)->GetValue(IS_ASSEMBLED))
//         // if (this->GetValue(ACTIVATION_LEVEL) != 0)
//             for (IndexType i = 0; i < number_of_nodes_slave; i++)
//                 for (IndexType idim = 0; idim < 2; idim++) {
//                     const int iglob = 2*i+idim + 2*number_of_nodes_master;
//                     const int i_index = 2*i+idim;

//                     Vector extension_sigma_w_n(2);
//                     extension_sigma_w_n[0] = (DB_sum_slave(0, i_index)* normal_true_slave[0] + DB_sum_slave(2, i_index)* normal_true_slave[1]);
//                     extension_sigma_w_n[1] = (DB_sum_slave(2, i_index)* normal_true_slave[0] + DB_sum_slave(1, i_index)* normal_true_slave[1]);
//                     double extension_sigma_w_n_n = extension_sigma_w_n[0]*normal_true_slave[0] + extension_sigma_w_n[1]*normal_true_slave[1];
                    
//                     Vector sigma_w_n(2);
//                     sigma_w_n[1] = (DB_slave(2, i_index)* normal_true_slave[0] + DB_slave(1, i_index)* normal_true_slave[1]);
//                     sigma_w_n[0] = (DB_slave(0, i_index)* normal_true_slave[0] + DB_slave(2, i_index)* normal_true_slave[1]);
//                     for (IndexType j = 0; j < number_of_nodes_slave; j++) {
//                         double sigma_w_n_n = sigma_w_n[0]*normal_true_slave[0] + sigma_w_n[1]*normal_true_slave[1];
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim + 2*number_of_nodes_master;
//                             const int j_index = 2*j+jdim;

//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* normal_true_slave[0] + DB_sum_slave(2, j_index)* normal_true_slave[1]);
//                             extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* normal_true_slave[0] + DB_sum_slave(1, j_index)* normal_true_slave[1]);
//                             double extension_sigma_u_n_n = extension_sigma_u_n[0]*normal_true_slave[0] + extension_sigma_u_n[1]*normal_true_slave[1];

                            
//                             rLeftHandSideMatrix(iglob, jglob) += active_constant_nitsche_terms*0.5*theta*penalty_slave * n_ntilde_slave 
//                                                                 * extension_sigma_u_n_n * sigma_w_n_n
//                                                                 * integration_weight_slave;
                                                            
//                         }
//                     }
//                 }

//         //***************
//         // NITSCHE ACTIVATION TERMS  */
//         // Nitsche Activation Term on Master

//         // +0.5*(n_ntilde_1) * 
//         //      [ 1/gamma_1 * (E(u_1) - E(u_2)) \dot n_1 
//         //       + E(sigma_1(u_1))] 
//         //      * (w_1 + w_2|P - theta*gamma_1 * E(sigma_1(v)))
//         // 
//         if ((this->GetValue(ACTIVATION_LEVEL) == 1 || this->GetValue(ACTIVATION_LEVEL) == 3)) {
//             // v_1 related
//             for (IndexType i = 0; i < number_of_nodes_master; i++)
//             {
//                 for (IndexType idim = 0; idim < 2; idim++) {
//                     const int iglob = 2*i+idim;
//                     Vector extension_sigma_w_n(2);
//                     extension_sigma_w_n[0] = (DB_sum_master(0, iglob)* normal_true_master[0] + DB_sum_master(2, iglob)* normal_true_master[1]);
//                     extension_sigma_w_n[1] = (DB_sum_master(2, iglob)* normal_true_master[0] + DB_sum_master(1, iglob)* normal_true_master[1]);
//                     double extension_sigma_w_n_n = extension_sigma_w_n[0]*normal_true_master[0] + extension_sigma_w_n[1]*normal_true_master[1];

//                     Vector sigma_w_n(2);
//                     sigma_w_n[0] = (DB_master(0, iglob)* normal_true_master[0] + DB_master(2, iglob)* normal_true_master[1]);
//                     sigma_w_n[1] = (DB_master(2, iglob)* normal_true_master[0] + DB_master(1, iglob)* normal_true_master[1]);
//                     double sigma_w_n_n = sigma_w_n[0]*normal_true_master[0] + sigma_w_n[1]*normal_true_master[1];

//                     for (IndexType j = 0; j < number_of_nodes_master; j++){
//                         // u_1 related
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim;
                            
//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* normal_true_master[0] + DB_sum_master(2, jglob)* normal_true_master[1]);
//                             extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* normal_true_master[0] + DB_sum_master(1, jglob)* normal_true_master[1]);
//                             double extension_sigma_u_n_n = extension_sigma_u_n[0]*normal_true_master[0] + extension_sigma_u_n[1]*normal_true_master[1];

//                             // E(u_1)*v1_n
//                             rLeftHandSideMatrix(iglob, jglob) += 0.5*n_ntilde_master/penalty_master * H_sum_master(0,j) * normal_true_master[jdim]
//                                                                 * H_master(0, i) *normal_true_master[idim]
//                                                                 * integration_weight_master;
                            
//                             // -E(sigma(u_1))*v1_m
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_master * extension_sigma_u_n_n 
//                                                                 * H_master(0, i) *normal_true_master[idim]
//                                                                 * integration_weight_master;
                            
//                             // - theta*gamma_1 E(u_1)* E(sigma(v))
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_master*theta
//                                                                 * H_sum_master(0,j) * normal_true_master[jdim]
//                                                                 * sigma_w_n_n
//                                                                 * integration_weight_master;
                                                                
                            
//                             // + theta*gamma_1 E(sigma(u_1))* E(sigma(v))
//                             rLeftHandSideMatrix(iglob, jglob) += active_constant_nitsche_terms*0.5*n_ntilde_master*theta*penalty_master 
//                                                                 * extension_sigma_u_n_n
//                                                                 * sigma_w_n_n
//                                                                 * integration_weight_master;
                            
//                         }
//                     }
//                     // u_2 related
//                     for (IndexType j = 0; j < number_of_nodes_slave; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim + 2*number_of_nodes_master;
                            
//                             // -E(u_2)*v1_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_master/penalty_master 
//                                                                 * H_sum_slave(0,j) * normal_true_master[jdim]
//                                                                 * H_master(0, i) *normal_true_master[idim]
//                                                                 * integration_weight_master;
                            
                            
//                             // - theta*gamma_1 -E(u_2)* E(sigma(v))
//                             rLeftHandSideMatrix(iglob, jglob) += 0.5*n_ntilde_master*theta
//                                                                 * H_sum_slave(0,j) * normal_true_master[jdim]
//                                                                 * sigma_w_n_n
//                                                                 * integration_weight_master;

//                         }
//                     }

//                     if (CalculateResidualVectorFlag) {

//                         // g_n1 * v1_n
//                         rRightHandSideVector[iglob] += 0.5*n_ntilde_master/penalty_master 
//                                                         * initial_normal_gap_true_master
//                                                         * H_master(0, i) *normal_true_master[idim]
//                                                         * integration_weight_master;
                        
//                         // - theta*gamma_1 -E(u_2)* E(sigma(v))
//                         rRightHandSideVector[iglob] -= 0.5*n_ntilde_master*theta
//                                                         * initial_normal_gap_true_master
//                                                         * sigma_w_n_n
//                                                         * integration_weight_master;
//                     }
//                 }
//             }

//             // v_2 related
//             for (IndexType i = 0; i < number_of_nodes_slave; i++)
//             {
//                 for (IndexType idim = 0; idim < 2; idim++) {
//                     const int iglob = 2*i+idim +2*number_of_nodes_master;

//                     // u_1 related
//                     for (IndexType j = 0; j < number_of_nodes_master; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim;
                            
//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_master(0, jglob)* normal_true_master[0] + DB_sum_master(2, jglob)* normal_true_master[1]);
//                             extension_sigma_u_n[1] = (DB_sum_master(2, jglob)* normal_true_master[0] + DB_sum_master(1, jglob)* normal_true_master[1]);
//                             double extension_sigma_u_n_n = extension_sigma_u_n[0]*normal_true_master[0] + extension_sigma_u_n[1]*normal_true_master[1];

//                             // E(u_1)*v2_n
//                             rLeftHandSideMatrix(iglob, jglob) += 0.5*n_ntilde_master/penalty_master 
//                                                                 * H_sum_master(0,j) * normal_true_master[jdim]
//                                                                 * H_slave(0, i) *normal_true_slave[idim]
//                                                                 * integration_weight_master * corrective_factor_curvature_on_slave;
                            
//                             // -E(sigma(u_1))*v2_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_master * extension_sigma_u_n_n 
//                                                                 * H_slave(0, i) *normal_true_slave[idim]
//                                                                 * integration_weight_master * corrective_factor_curvature_on_slave;
                        
//                         }
//                     }
//                     // u_2 related
//                     for (IndexType j = 0; j < number_of_nodes_slave; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim+2*number_of_nodes_master;
                        
//                             // -E(u_2)*v2_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_master/penalty_master 
//                                                                 * H_sum_slave(0,j) * normal_true_master[jdim]
//                                                                 * H_slave(0, i) *normal_true_slave[idim]
//                                                                 * integration_weight_master * corrective_factor_curvature_on_slave;
//                         }
//                     }

//                     if (CalculateResidualVectorFlag) {

//                         // g_n1 * v2_n
//                         rRightHandSideVector[iglob] += 0.5*n_ntilde_master/penalty_master 
//                                                         * initial_normal_gap_true_master
//                                                         * H_slave(0, i) *normal_true_slave[idim]
//                                                         * integration_weight_master * corrective_factor_curvature_on_slave;  
//                     }
//                 }
//             }

//         }

//         // Nitsche Activation Term on Slave ----------------------------------------------------

//         // +0.5*(n_ntilde_2) * 
//         //      [ 1/gamma_2 * (E(u_2) - E(u_1)) \dot n_2 
//         //       + E(sigma_2(u_2))] 
//         //      * (w_2 + w_1|P - theta*gamma_2 * E(sigma_2(v)))
//         // 
     
//         if ((this->GetValue(ACTIVATION_LEVEL) == 2 || this->GetValue(ACTIVATION_LEVEL) == 3)) {
            
//             // v_2 related
//             for (IndexType i = 0; i < number_of_nodes_slave; i++)
//             {
//                 for (IndexType idim = 0; idim < 2; idim++) {
//                     const int iglob = 2*i+idim + 2*number_of_nodes_master;
//                     const int i_index = 2*i+idim;
//                     Vector extension_sigma_w_n(2);
//                     extension_sigma_w_n[0] = (DB_sum_slave(0, i_index)* normal_true_slave[0] + DB_sum_slave(2, i_index)* normal_true_slave[1]);
//                     extension_sigma_w_n[1] = (DB_sum_slave(2, i_index)* normal_true_slave[0] + DB_sum_slave(1, i_index)* normal_true_slave[1]);
//                     double extension_sigma_w_n_n = extension_sigma_w_n[0]*normal_true_slave[0] + extension_sigma_w_n[1]*normal_true_slave[1];


//                     Vector sigma_w_n(2);
//                     sigma_w_n[1] = (DB_slave(2, i_index)* normal_true_slave[0] + DB_slave(1, i_index)* normal_true_slave[1]);
//                     sigma_w_n[0] = (DB_slave(0, i_index)* normal_true_slave[0] + DB_slave(2, i_index)* normal_true_slave[1]);
//                     double sigma_w_n_n = sigma_w_n[0]*normal_true_slave[0] + sigma_w_n[1]*normal_true_slave[1];

//                     // u_2 related
//                     for (IndexType j = 0; j < number_of_nodes_slave; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim + 2*number_of_nodes_master;
//                             const int j_index = 2*j+jdim;
                            
//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* normal_true_slave[0] + DB_sum_slave(2, j_index)* normal_true_slave[1]);
//                             extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* normal_true_slave[0] + DB_sum_slave(1, j_index)* normal_true_slave[1]);
//                             double extension_sigma_u_n_n = extension_sigma_u_n[0]*normal_true_slave[0] + extension_sigma_u_n[1]*normal_true_slave[1];

//                             // E(u_2)*v2_n
//                             rLeftHandSideMatrix(iglob, jglob) += 0.5*n_ntilde_slave/penalty_slave 
//                                                                 * H_sum_slave(0,j) * normal_true_slave[jdim]
//                                                                 * H_slave(0, i) *normal_true_slave[idim]
//                                                                 * integration_weight_slave;
                                                                
                            
//                             // -E(sigma(u_2))*v2_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_slave * extension_sigma_u_n_n 
//                                                                 * H_slave(0, i) *normal_true_slave[idim]
//                                                                 * integration_weight_slave;
                            
//                             // - theta*gamma_2 E(u_2)* E(sigma(v))
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_slave*theta
//                                                                 * H_sum_slave(0,j) * normal_true_slave[jdim]
//                                                                 * sigma_w_n_n
//                                                                 * integration_weight_slave;
                            

//                             // + theta*gamma_1 E(sigma(u_2))* E(sigma(v))
//                             rLeftHandSideMatrix(iglob, jglob) += active_constant_nitsche_terms*0.5*n_ntilde_slave*theta*penalty_slave 
//                                                                 * extension_sigma_u_n_n
//                                                                 * sigma_w_n_n
//                                                                 * integration_weight_slave;

//                         }
//                     }
//                     // u_1 related
//                     for (IndexType j = 0; j < number_of_nodes_master; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim;
                            
//                             // -E(u_1)*v2_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_slave/penalty_slave 
//                                                                 * H_sum_master(0,j) * normal_true_slave[jdim]
//                                                                 * H_slave(0, i) *normal_true_slave[idim]
//                                                                 * integration_weight_slave;
                            
                            
//                             // - theta*gamma_2 -E(u_1)* E(sigma(v2))
//                             rLeftHandSideMatrix(iglob, jglob) += 0.5*n_ntilde_slave*theta
//                                                                 * H_sum_master(0,j) * normal_true_slave[jdim]
//                                                                 * sigma_w_n_n
//                                                                 * integration_weight_slave;
//                         }
//                     }

//                     if (CalculateResidualVectorFlag) {

//                         // g_n2 * v2_n
//                         rRightHandSideVector[iglob] += 0.5*n_ntilde_slave/penalty_slave 
//                                                         * initial_normal_gap_true_slave
//                                                         * H_slave(0, i) *normal_true_slave[idim]
//                                                         * integration_weight_slave;
                        
//                         // - theta*gamma_2 g_n2* E(sigma(v2))
//                         rRightHandSideVector[iglob] -= 0.5*n_ntilde_slave*theta
//                                                         * initial_normal_gap_true_slave
//                                                         * sigma_w_n_n
//                                                         * integration_weight_slave;
//                     }
//                 }
//             }

//             // v_1 related
//             for (IndexType i = 0; i < number_of_nodes_master; i++)
//             {
//                 for (IndexType idim = 0; idim < 2; idim++) {
//                     const int iglob = 2*i+idim;

//                     // u_2 related
//                     for (IndexType j = 0; j < number_of_nodes_slave; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim + 2*number_of_nodes_master;
//                             const int j_index = 2*j+jdim;
                            
//                             Vector extension_sigma_u_n(2);
//                             extension_sigma_u_n[0] = (DB_sum_slave(0, j_index)* normal_true_slave[0] + DB_sum_slave(2, j_index)* normal_true_slave[1]);
//                             extension_sigma_u_n[1] = (DB_sum_slave(2, j_index)* normal_true_slave[0] + DB_sum_slave(1, j_index)* normal_true_slave[1]);
//                             double extension_sigma_u_n_n = extension_sigma_u_n[0]*normal_true_slave[0] + extension_sigma_u_n[1]*normal_true_slave[1];
//                             // E(u_2)*v1_n
//                             rLeftHandSideMatrix(iglob, jglob) += 0.5*n_ntilde_slave/penalty_slave 
//                                                                 * H_sum_slave(0,j) * normal_true_slave[jdim]
//                                                                 * H_master(0, i) *normal_true_master[idim]
//                                                                 * integration_weight_slave * corrective_factor_curvature_on_master;
                            
//                             // -E(sigma(u_2))*v1_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_slave * extension_sigma_u_n_n 
//                                                                 * H_master(0, i) *normal_true_master[idim]
//                                                                 * integration_weight_slave * corrective_factor_curvature_on_master;
//                         }
//                     }
//                     // u_1 related
//                     for (IndexType j = 0; j < number_of_nodes_master; j++){
//                         for (IndexType jdim = 0; jdim < 2; jdim++) {
//                             const int jglob = 2*j+jdim;
                        
//                             // -E(u_1)*v1_n
//                             rLeftHandSideMatrix(iglob, jglob) -= 0.5*n_ntilde_slave/penalty_slave 
//                                                                 * H_sum_master(0,j) * normal_true_slave[jdim]
//                                                                 * H_master(0, i) *normal_true_master[idim]
//                                                                 * integration_weight_slave * corrective_factor_curvature_on_master;
//                         }
//                     }

//                     if (CalculateResidualVectorFlag) {

//                         // g_n2 * v1_n
//                         rRightHandSideVector[iglob] += 0.5*n_ntilde_slave/penalty_slave 
//                                                         * initial_normal_gap_true_slave
//                                                         * H_master(0, i) *normal_true_master[idim]
//                                                         * integration_weight_slave* corrective_factor_curvature_on_master;  
//                     }
//                 }
//             }
//         }
//         projection_node_slave.GetValue(NEIGHBOUR_GEOMETRY)->SetValue(IS_ASSEMBLED, true);

//         Vector temp = ZeroVector(mat_size_2);

//         GetValuesVector(temp, 2);

//         // RHS = ExtForces - K*temp;
//         noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

//         for (unsigned int i = 0; i < number_of_nodes_slave; i++) {

//             std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
//             outputFile << r_geometry_slave[i].GetId() << "  " <<r_geometry_slave[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
//             outputFile.close();
//         }

//         // for (unsigned int i = 0; i < number_of_nodes_master; i++) {

//         //     std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
//         //     outputFile << r_geometry_master[i].GetId() << "  " <<r_geometry_master[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
//         //     outputFile.close();
//         // }

//         // /////////////////////////////////////////////////////////////////////////////////////////
     
//         KRATOS_CATCH("")
//     }

//     int SbmContact2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
//     {
//         KRATOS_ERROR_IF_NOT((*mpPropMaster).Has(PENALTY_FACTOR))
//             << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
//         return 0;
//     }

//     void SbmContact2DCondition::EquationIdVector(
//         EquationIdVectorType& rResult,
//         const ProcessInfo& rCurrentProcessInfo
//     ) const
//     {
//         const auto& r_geometry_master = GetMasterGeometry();
//         const SizeType number_of_nodes_master = r_geometry_master.size();

//         const auto& r_geometry_slave = GetSlaveGeometry();
//         const SizeType number_of_nodes_slave = r_geometry_slave.size();

//         const SizeType number_of_nodes = (number_of_nodes_master + number_of_nodes_slave);

//         if (rResult.size() != 2 * number_of_nodes)
//             rResult.resize(2 * number_of_nodes, false);

//         for (IndexType i = 0; i < number_of_nodes_master; ++i) {
//             const IndexType index = i * 2;
//             const auto& r_node = r_geometry_master[i];
//             rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
//             rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
//         }

//         const int shift = number_of_nodes_master*2;
//         for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
//             const IndexType index = i * 2;
//             const auto& r_node = r_geometry_slave[i];
//             rResult[index+shift] = r_node.GetDof(DISPLACEMENT_X).EquationId();
//             rResult[index + 1+shift] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
//         }
//     }

//     void SbmContact2DCondition::GetDofList(
//         DofsVectorType& rElementalDofList,
//         const ProcessInfo& rCurrentProcessInfo
//     ) const
//     {
//         const auto& r_geometry_master = GetMasterGeometry();
//         const SizeType number_of_nodes_master = r_geometry_master.size();

//         const auto& r_geometry_slave = GetSlaveGeometry();
//         const SizeType number_of_nodes_slave = r_geometry_slave.size();

//         rElementalDofList.resize(0);
//         rElementalDofList.reserve(2 * (number_of_nodes_master+number_of_nodes_slave));

//         for (IndexType i = 0; i < number_of_nodes_master; ++i) {
//             const auto& r_node = r_geometry_master[i];
//             rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
//             rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
//         }

//         for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
//             const auto& r_node = r_geometry_slave[i];
//             rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
//             rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
//         }

//     };


//     void SbmContact2DCondition::GetValuesVector(
//         Vector& rValues, IndexType index) const
//     {
//         if (index == QuadraturePointCouplingGeometry2D<Point>::Master) {
            
//             SizeType number_of_control_points = GetMasterGeometry().size();
//             const SizeType mat_size = number_of_control_points * 2;

//             if (rValues.size() != mat_size)
//                 rValues.resize(mat_size, false);

//             for (IndexType i = 0; i < number_of_control_points; ++i)
//             {
//                 const array_1d<double, 3>& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
//                 IndexType index = i * 2;

//                 rValues[index] = displacement[0];
//                 rValues[index + 1] = displacement[1];
//             }
//         } else if (index == QuadraturePointCouplingGeometry2D<Point>::Slave) {
//             SizeType number_of_control_points = GetSlaveGeometry().size();
//             const SizeType mat_size = number_of_control_points * 2;

//             if (rValues.size() != mat_size)
//                 rValues.resize(mat_size, false);

//             for (IndexType i = 0; i < number_of_control_points; ++i)
//             {
//                 const array_1d<double, 3>& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
//                 IndexType index = i * 2;

//                 rValues[index] = displacement[0];
//                 rValues[index + 1] = displacement[1];
//             }
//         } else if (index == 2) {
            
//             SizeType number_of_control_points_master = GetMasterGeometry().size();
//             const SizeType mat_size_master = number_of_control_points_master * 2;
//             SizeType number_of_control_points_slave = GetSlaveGeometry().size();
//             const SizeType mat_size_slave = number_of_control_points_slave * 2;
//             const SizeType mat_size = mat_size_master + mat_size_slave;

//             if (rValues.size() != mat_size)
//                 rValues.resize(mat_size, false);

//             for (IndexType i = 0; i < number_of_control_points_master; ++i)
//             {
//                 const array_1d<double, 3>& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
//                 IndexType index = i * 2;

//                 rValues[index] = displacement[0];
//                 rValues[index + 1] = displacement[1];
//             }

//             for (IndexType i = 0; i < number_of_control_points_slave; ++i)
//             {
//                 const array_1d<double, 3>& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
//                 IndexType index = i * 2 + mat_size_master;

//                 rValues[index] = displacement[0];
//                 rValues[index + 1] = displacement[1];
//             }
//         }
//     }

//     void SbmContact2DCondition::GetStrainVector(
//         Vector& strainVector, IndexType index) const
//     {
//         const auto& r_geometry = GetGeometry().GetGeometryPart(index);

//         const SizeType dim = 2;//r_geometry.WorkingSpaceDimension();

//         const SizeType number_of_nodes = r_geometry.size();


//         // Shape function derivatives 
//         // Initialize Jacobian
//         GeometryType::JacobiansType J0;

//         // Initialize DN_DX
//         Matrix DN_DX(number_of_nodes,2);
//         Matrix InvJ0(dim,dim);

//         const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
//         r_geometry.Jacobian(J0,this->GetIntegrationMethod());

//         Matrix Jacobian = ZeroMatrix(2,2);
//         Jacobian(0,0) = J0[0](0,0);
//         Jacobian(0,1) = J0[0](0,1);
//         Jacobian(1,0) = J0[0](1,0);
//         Jacobian(1,1) = J0[0](1,1);

//         // Calculating inverse jacobian and jacobian determinant
//         double DetJ0;
//         MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);


//         // retrieve shape function values
//         Vector displacements(2*number_of_nodes);
//         GetValuesVector(displacements, index);

//         // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
//         noalias(DN_DX) = prod(DN_De[0],InvJ0);

//         // MODIFIED
//         Matrix B = ZeroMatrix(3,number_of_nodes);

//         CalculateB(B, DN_DX, number_of_nodes);

//         if (strainVector.size() != dim) strainVector.resize(dim);

//         strainVector = prod(B, displacements);
//     }

//     /**
//      * @brief 
//      * 
//      * @param StrainVector 
//      * @param index 
//      * @param rCurrentProcessInfo 
//      */
//     void SbmContact2DCondition::SetConstitutiveVariables(
//         Vector& StrainVector, IndexType index, const Kratos::ProcessInfo& rCurrentProcessInfo) 
//     {
//         const auto& r_geometry = GetGeometry().GetGeometryPart(index);
//         const SizeType dim = 2;//r_geometry.WorkingSpaceDimension();

//         ConstitutiveLaw::Pointer rpConstitutiveLaw = GetConstitutiveLaw(index);

//         PropertiesType r_prop = GetProperty(index);

//         ConstitutiveLaw::Parameters Values(r_geometry, r_prop, rCurrentProcessInfo);

//         const SizeType strain_size = rpConstitutiveLaw->GetStrainSize();
//         // Set constitutive law flags:
//         Flags& ConstitutiveLawOptions=Values.GetOptions();

//         ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
//         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
//         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

//         ConstitutiveVariables this_constitutive_variables(strain_size);
    
//         Values.SetStrainVector(StrainVector);
//         Values.SetStressVector(this_constitutive_variables.StressVector);

//         Values.SetConstitutiveMatrix(this_constitutive_variables.D);
//         rpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

//         if (index == 0) {
//             this->SetValue(CONSTITUTIVE_MATRIX_MASTER, Values.GetConstitutiveMatrix());
//             this->SetValue(STRESS_MASTER, Values.GetStressVector());
//         } 
//         else if (index == 1)
//         {
//             this->SetValue(CONSTITUTIVE_MATRIX_SLAVE, Values.GetConstitutiveMatrix());
//             this->SetValue(STRESS_SLAVE, Values.GetStressVector());
//         }
        
//     }

//     void SbmContact2DCondition::SetGap() 
//     {
//         NodeType projection_node_master;         NodeType projection_node_slave;
//         //
//         projection_node_master = GetMasterGeometry().GetValue(NEIGHBOUR_NODES)[0];
//         projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];

//         Vector normal_physical_space_master = this->GetValue(NORMAL_MASTER);

//         Vector skin_master_coord_deformed = projection_node_master + GetValue(DISPLACEMENT_MASTER);
//         Vector skin_slave_coord_deformed  = projection_node_slave + GetValue(DISPLACEMENT_SLAVE);
        
//         SetValue(SKIN_MASTER_COORDINATES, projection_node_master);
//         SetValue(SKIN_SLAVE_COORDINATES, projection_node_slave);

//         Vector gap_deformed = skin_slave_coord_deformed - skin_master_coord_deformed;

//         SetValue(GAP, gap_deformed);
//     }

//     void SbmContact2DCondition::CalculateB(
//         Matrix& rB, 
//         Matrix& r_DN_DX,
//         const SizeType number_of_control_points) const
//     {
        
//         const SizeType mat_size = number_of_control_points * 2;

//         if (rB.size1() != 3 || rB.size2() != mat_size)
//             rB.resize(3, mat_size);
//         noalias(rB) = ZeroMatrix(3, mat_size);

//         for (IndexType r = 0; r < mat_size; r++)
//         {
//             // local node number kr and dof direction dirr
//             IndexType kr = r / 2;
//             IndexType dirr = r % 2;

//             rB(0, r) = r_DN_DX(kr,0) * (1-dirr);
//             rB(1, r) = r_DN_DX(kr,1) * dirr;
//             rB(2, r) = r_DN_DX(kr,0) * (dirr) + r_DN_DX(kr,1) * (1-dirr);
//         }
//     }


//     void SbmContact2DCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
//     {
//         ConstitutiveLaw::Parameters constitutive_law_parameters_master(
//             GetMasterGeometry(), *mpPropMaster, rCurrentProcessInfo);
        
//         // Master
//         const auto& r_geometry_master = GetMasterGeometry();
//         const SizeType number_of_nodes_master = r_geometry_master.size();

//         // Slave
//         const auto& r_geometry_slave = GetSlaveGeometry();
//         const SizeType number_of_nodes_slave = r_geometry_slave.size();


//         NodeType projection_node_master;         NodeType projection_node_slave;
//         //
//         projection_node_master = r_geometry_master.GetValue(NEIGHBOUR_NODES)[0];
//         projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];
//         projection_node_slave.GetValue(NEIGHBOUR_GEOMETRY)->SetValue(IS_ASSEMBLED, false);

//         mpConstitutiveLawMaster->FinalizeMaterialResponse(constitutive_law_parameters_master, ConstitutiveLaw::StressMeasure_PK2);

//         ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
//             GetSlaveGeometry(), *mpPropSlave, rCurrentProcessInfo);

//         mpConstitutiveLawSlave->FinalizeMaterialResponse(constitutive_law_parameters_slave, ConstitutiveLaw::StressMeasure_PK2);

//         //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
//         Vector normal_true_master = GetValue(NORMAL_MASTER);
//         Vector stress_vector_master = GetValue(STRESS_MASTER);

//         Vector sigma_n(2);

//         sigma_n[0] = stress_vector_master[0]*normal_true_master[0] + stress_vector_master[2]*normal_true_master[1];
//         sigma_n[1] = stress_vector_master[2]*normal_true_master[0] + stress_vector_master[1]*normal_true_master[1];

//         SetValue(NORMAL_STRESS, sigma_n);
//     }

//     void SbmContact2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

//         // InitializeMaterial();
//         ConstitutiveLaw::Parameters constitutive_law_parameters_master(
//             GetMasterGeometry(), (*mpPropMaster), rCurrentProcessInfo);

//         mpConstitutiveLawMaster->InitializeMaterialResponse(constitutive_law_parameters_master, ConstitutiveLaw::StressMeasure_PK2);

//         // InitializeMaterial(); //slave
//         ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
//             GetSlaveGeometry(), (*mpPropSlave), rCurrentProcessInfo);

//         mpConstitutiveLawSlave->InitializeMaterialResponse(constitutive_law_parameters_slave, ConstitutiveLaw::StressMeasure_PK2);

//         this->InitializeNonLinearIteration(rCurrentProcessInfo);

//         // #############################################################
//         SetValue(ACTIVATION_LEVEL, 0);  
//         SetValue(YOUNG_MODULUS_MASTER, (*mpPropMaster)[YOUNG_MODULUS]);
//         SetValue(YOUNG_MODULUS_SLAVE, (*mpPropSlave)[YOUNG_MODULUS]);

//         const auto& r_geometry_master = GetMasterGeometry();
//         const SizeType number_of_nodes_master = r_geometry_master.size();

//         const auto& r_geometry_slave = GetSlaveGeometry();
//         const SizeType number_of_nodes_slave = r_geometry_slave.size();
//     }

//     /**
//      * @brief 
//      * 
//      * @param rCurrentProcessInfo 
//      */
//     void SbmContact2DCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo){

//         // Master
//         const auto& r_geometry_master = GetMasterGeometry();
//         const SizeType number_of_nodes_master = r_geometry_master.size();

//         // Slave
//         const auto& r_geometry_slave = GetSlaveGeometry();
//         const SizeType number_of_nodes_slave = r_geometry_slave.size();


//         NodeType projection_node_master;         
//         //
//         projection_node_master = r_geometry_master.GetValue(NEIGHBOUR_NODES)[0];
//         auto& projection_node_slave = projection_node_master.GetValue(NEIGHBOUR_NODES)[0];
//         projection_node_slave.GetValue(NEIGHBOUR_GEOMETRY)->SetValue(IS_ASSEMBLED, false);

//         auto& non_const_node = const_cast<NodeType&>(*(r_geometry_master.GetValue(NEIGHBOUR_NODES)(0)));
//         non_const_node.SetValue(IS_ASSEMBLED, false);

//         // master
//         const SizeType dim = 2; //r_geometry_master.WorkingSpaceDimension();

//         Matrix DN_DX_master(number_of_nodes_master,2); 
//         Matrix DN_DX_slave(number_of_nodes_slave,2); 
//         // //-----------------------------------------------------------------------------

//         // NORMALS
//         // master
//         array_1d<double, 3> tangent; 
//         array_1d<double, 3> normal; 

//         r_geometry_master.Calculate(LOCAL_TANGENT, tangent); // Gives the result in the parameter space !!
//         double magnitude = std::sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
        
//         // NEW FOR GENERAL JACOBIAN
//         normal[0] = + tangent[1] / magnitude;
//         normal[1] = - tangent[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT
//         normal[2] = 0.0;

//         tangent /= norm_2(tangent);

//         // ##
//         mNormalMaster = normal;
//         // ##
//         // slave
//         r_geometry_slave.Calculate(LOCAL_TANGENT, tangent); // Gives the result in the parameter space !!
//         magnitude = std::sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
        
//         // NEW FOR GENERAL JACOBIAN
//         normal[0] = + tangent[1] / magnitude;
//         normal[1] = - tangent[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT

//         tangent /= norm_2(tangent);
//         // ##

//         mNormalSlave = normal;
//         // ##

//         // SET REFERENCE WEIGHT

//         const double thickness = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;

//         const double IntToReferenceWeight = r_geometry_master.IntegrationPoints()[0].Weight() * thickness;

//         // ##
//         SetValue(INTEGRATION_WEIGHT, IntToReferenceWeight);
//         // ################################################################
//         // SBM PARAMETERS
//         // ################################################################
//         const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
//         int basis_functions_order_master = std::sqrt(DN_De_master[0].size1()) - 1;
//         std::vector<Matrix> shape_function_derivatives_master;

//         const GeometryType::ShapeFunctionsGradientsType& DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
//         int basis_functions_order_slave = std::sqrt(DN_De_slave[0].size1()) - 1;
//         std::vector<Matrix> shape_function_derivatives_slave;

//         // SET DISPLACEMENTS on the true boundary
//         // master
//         Vector coefficient_master(dim*number_of_nodes_master);
//         GetValuesVector(coefficient_master, QuadraturePointCouplingGeometry2D<Point>::Master);

//         const Matrix& N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());
//         Matrix H_sum = ZeroMatrix(1, number_of_nodes_master);
//         // Compute all the derivatives of the basis functions involved
//         for (int n = 1; n <= basis_functions_order_master; n++) {
//             shape_function_derivatives_master.push_back(r_geometry_master.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
//         }

//         for (IndexType i = 0; i < number_of_nodes_master; ++i)
//         {
//             double H_taylor_term = 0.0; // Reset for each node
//             for (int n = 1; n <= basis_functions_order_master; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = shape_function_derivatives_master[n-1];
//                 for (int k = 0; k <= n; k++) {
//                     int n_k = n - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term += ComputeTaylorTerm(derivative, mDistanceMaster[0], n_k, mDistanceMaster[1], k);
//                 }
//             }
//             H_sum(0,i) = H_taylor_term + N_master(0, i);
//         }
//         // reset matrix for matrix to vector product
//         Matrix H_master_true = ZeroMatrix(dim, dim*number_of_nodes_master);
//         for (IndexType i = 0; i < number_of_nodes_master; ++i)
//         {
//             for (IndexType i_dim = 0; i_dim < dim; i_dim++) {
//                 H_master_true(i_dim, dim*i+i_dim) = H_sum(0, i);
//             }
//         }

//         Vector displacement_master_true_sub = prod(H_master_true,coefficient_master);

//         Vector displacement_master_true = ZeroVector(3); displacement_master_true[0] = displacement_master_true_sub[0]; displacement_master_true[1] = displacement_master_true_sub[1];

//         // ##
//         this->SetValue(DISPLACEMENT_MASTER, displacement_master_true);
//         // ##
//         //  slave
//         Vector coefficient_slave(dim*number_of_nodes_slave);
//         GetValuesVector(coefficient_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);
//         const Matrix& N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());
//         H_sum = ZeroMatrix(1, number_of_nodes_slave);
//         // Compute all the derivatives of the basis functions involved
//         for (int n = 1; n <= basis_functions_order_slave; n++) {
//             shape_function_derivatives_slave.push_back(r_geometry_slave.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
//         }

//         for (IndexType i = 0; i < number_of_nodes_slave; ++i)
//         {
//             double H_taylor_term = 0.0; // Reset for each node
//             for (int n = 1; n <= basis_functions_order_slave; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = shape_function_derivatives_slave[n-1];
//                 for (int k = 0; k <= n; k++) {
//                     int n_k = n - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term += ComputeTaylorTerm(derivative, mDistanceSlave[0], n_k, mDistanceSlave[1], k);
//                 }
//             }
//             H_sum(0,i) = H_taylor_term + N_slave(0, i);
//         }
//         // reset matrix for matrix to vector product
//         Matrix H_slave_true = ZeroMatrix(dim, dim*number_of_nodes_slave);
//         for (IndexType i = 0; i < number_of_nodes_slave; ++i)
//         {
//             for (IndexType i_dim = 0; i_dim < dim; i_dim++) {
//                 H_slave_true(i_dim, dim*i+i_dim) = H_sum(0, i);
//             }
//         }
//         Vector displacement_slave_true_sub = prod(H_slave_true,coefficient_slave);

//         Vector displacement_slave_true = ZeroVector(3); displacement_slave_true[0] = displacement_slave_true_sub[0]; displacement_slave_true[1] = displacement_slave_true_sub[1];

//         // ##
//         this->SetValue(DISPLACEMENT_SLAVE, displacement_slave_true);
//         // ##

//         // SET STRAINS/STRESSES/CONSTITUTIVE MATRIX 
//         // master
//         // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
//         noalias(DN_DX_master) = DN_De_master[0];
//         Matrix H_grad_master = ZeroMatrix(number_of_nodes_master, 2); Matrix H_grad_slave = ZeroMatrix(number_of_nodes_slave, 2);
//         std::vector<Matrix> n_shape_function_derivatives_master; std::vector<Matrix> n_shape_function_derivatives_slave;

//         for (int n = 1; n <= basis_functions_order_master; n++) {
//             n_shape_function_derivatives_master.push_back(r_geometry_master.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
//         }
//         for (int n = 1; n <= basis_functions_order_slave; n++) {
//             n_shape_function_derivatives_slave.push_back(r_geometry_slave.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
//         }
//         for (IndexType i = 0; i < number_of_nodes_master; ++i)
//         {
//             double H_taylor_term_X = 0.0; // Reset for each node
//             double H_taylor_term_Y = 0.0; 
//             for (int n = 2; n <= basis_functions_order_master; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = n_shape_function_derivatives_master[n-1];
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceMaster[0], n_k, mDistanceMaster[1], k);
//                 }
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k+1); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceMaster[0], n_k, mDistanceMaster[1], k);
//                 }
//             }
//             H_grad_master(i,0) = H_taylor_term_X + DN_DX_master(i,0);
//             H_grad_master(i,1) = H_taylor_term_Y + DN_DX_master(i,1);
//         }                                                     

//         Matrix B_sum_master = ZeroMatrix(3, number_of_nodes_master);
//         CalculateB(B_sum_master, H_grad_master, number_of_nodes_master);

//         Vector strain_vector_master = prod(B_sum_master, coefficient_master);
        
        
//         // ##
//         this->SetValue(STRAIN_MASTER, strain_vector_master);
//         this->SetConstitutiveVariables(strain_vector_master, 0, rCurrentProcessInfo);
//         // ##
//         // SLAVE
//         noalias(DN_DX_slave) = DN_De_slave[0];
        
//         for (IndexType i = 0; i < number_of_nodes_slave; ++i)
//         {
//             double H_taylor_term_X = 0.0; // Reset for each node
//             double H_taylor_term_Y = 0.0; 
//             for (int n = 2; n <= basis_functions_order_slave; n++) {
//                 // Retrieve the appropriate derivative for the term
//                 Matrix& shapeFunctionDerivatives = n_shape_function_derivatives_slave[n-1];
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceSlave[0], n_k, mDistanceSlave[1], k);
//                 }
//                 for (int k = 0; k <= n-1; k++) {
//                     int n_k = n - 1 - k;
//                     double derivative = shapeFunctionDerivatives(i,k+1); 
//                     // Compute the Taylor term for this derivative
//                     H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceSlave[0], n_k, mDistanceSlave[1], k);
//                 }
//             }
//             H_grad_slave(i,0) = H_taylor_term_X + DN_DX_slave(i,0);
//             H_grad_slave(i,1) = H_taylor_term_Y + DN_DX_slave(i,1);
//         }                                                     

//         Matrix B_sum_slave = ZeroMatrix(3, number_of_nodes_slave);
//         CalculateB(B_sum_slave, H_grad_slave, number_of_nodes_slave);

//         Vector strain_vector_slave = prod(B_sum_slave, coefficient_slave);
        
//         // ##
//         this->SetValue(STRAIN_SLAVE, strain_vector_slave);
//         this->SetConstitutiveVariables(strain_vector_slave, 1, rCurrentProcessInfo);
//         // ##

//         SetGap();
//     }

//     /**
//      * @brief 
//      * 
//      * @param rCurrentProcessInfo 
//      */
//     void SbmContact2DCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo){

//         this->InitializeNonLinearIteration(rCurrentProcessInfo);
//     }


    


//     //----------------------------------------------------------------------------------
//     const Matrix SbmContact2DCondition::GetConstitutiveMatrix(IndexType index, Matrix& r_B, GeometryType r_geometry,
//                                                                    Vector& old_displacement, const Kratos::ProcessInfo& rCurrentProcessInfo,
//                                                                    Vector& stress_vector) {

//         ConstitutiveLaw::Pointer rpConstitutiveLaw = GetConstitutiveLaw(index);

//         PropertiesType r_prop = GetProperty(index);

//         ConstitutiveLaw::Parameters Values(r_geometry, r_prop, rCurrentProcessInfo);

//         const SizeType strain_size = rpConstitutiveLaw->GetStrainSize();
//         // Set constitutive law flags:
//         Flags& ConstitutiveLawOptions=Values.GetOptions();

//         ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
//         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
//         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
//         ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

//         ConstitutiveVariables this_constitutive_variables(strain_size);

//         Vector old_strain = prod(r_B,old_displacement);
    
//         // Values.SetStrainVector(this_constitutive_variables.StrainVector);
//         Values.SetStrainVector(old_strain);
//         Values.SetStressVector(this_constitutive_variables.StressVector);

//         Values.SetConstitutiveMatrix(this_constitutive_variables.D);
//         rpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        
//         stress_vector = Values.GetStressVector();

//         return Values.GetConstitutiveMatrix();
        
//     }



//     void SbmContact2DCondition::GetDeformed(const Matrix& N, Vector& reference_position, Vector& displacement, Vector& deformed_position) {
        
//         // compute the old displacement for the two points in the pair
//         double displacement_x = 0.0; double displacement_y = 0.0;
        
//         for (IndexType i = 0; i < GetMasterGeometry().size(); ++i)
//         {
//             // KRATOS_WATCH(r_geometry[i])
//             double output_solution_step_value_x = displacement[2*i];
//             double output_solution_step_value_y = displacement[2*i+1];
//             displacement_x += N(0, i) * output_solution_step_value_x;
//             displacement_y += N(0, i) * output_solution_step_value_y;
//         } 
//         deformed_position[0] = reference_position[0] + displacement_x; 
//         deformed_position[1] = reference_position[1] + displacement_y;
//     }


//     void SbmContact2DCondition::CalculateOnIntegrationPoints(
//         const Variable<double>& rVariable,
//         std::vector<double>& rValues,
//         const ProcessInfo& rCurrentProcessInfo
//         )
//     {
//         if (rValues.size() != 1)
//             rValues.resize(1);
//         if (rVariable == INTEGRATION_WEIGHT)
//             rValues[0] = this->GetValue(INTEGRATION_WEIGHT);
//         else if (rVariable == NORMAL_GAP)
//             rValues[0] = this->GetValue(NORMAL_GAP);
//     }

//     void SbmContact2DCondition::CalculateOnIntegrationPoints(
//         const Variable<Vector>& rVariable,
//         std::vector<Vector>& rValues,
//         const ProcessInfo& rCurrentProcessInfo
//         )
//     {
//         KRATOS_ERROR << "CalculateOnIntegrationPoints not ready for the SBM!!";
//         const SizeType dimension = 2;//GetMasterGeometry().WorkingSpaceDimension();

//         if (rValues.size() != dimension)
//             rValues.resize(dimension);
//         if (rVariable == NORMAL_STRESS)
//         {
//             Vector normal_physical_space = this->GetValue(NORMAL);
//             Vector sigma_voigt_master = GetValue(STRESS_MASTER);

//             Vector normal_stress_master(2);

//             normal_stress_master[0] = sigma_voigt_master[0]*normal_physical_space[0] + sigma_voigt_master[2]*normal_physical_space[1];
//             normal_stress_master[1] = sigma_voigt_master[2]*normal_physical_space[0] + sigma_voigt_master[1]*normal_physical_space[1];
            
//             this->SetValue(NORMAL_STRESS, normal_stress_master);

//             rValues[0] = normal_stress_master;
//         }
//         else if (rVariable == NORMAL_MASTER) 
//             rValues[0] = this->GetValue(NORMAL_MASTER);
//     }




//      // Function to compute a single term in the Taylor expansion
//     double SbmContact2DCondition::ComputeTaylorTerm(double derivative, double dx, int n_k, double dy, int k)
//     {
//         return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (Factorial(k) * Factorial(n_k));    
//     }

//     unsigned long long SbmContact2DCondition::Factorial(int n) 
//     {
//         if (n == 0) return 1;
//         unsigned long long result = 1;
//         for (int i = 2; i <= n; ++i) result *= i;
//         return result;
//     }

// } // Namespace Kratos