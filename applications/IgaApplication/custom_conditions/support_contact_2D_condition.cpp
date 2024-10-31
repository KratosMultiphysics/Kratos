//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andea Gorgi
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_contact_2D_condition.h"

namespace Kratos
{

    void SupportContact2DCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        InitializeMaterial();
    }


    void SupportContact2DCondition::InitializeMaterial()
    {
        KRATOS_TRY
        if ((*mpPropMaster)[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry_master = GetMasterGeometry();

            const auto& N_values = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());

            mpConstitutiveLawMaster = (*mpPropMaster)[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLawMaster->InitializeMaterial( *mpPropMaster, r_geometry_master, row(N_values , 0 ));

        } else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;


        if ( (*mpPropSlave)[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry_slave = GetSlaveGeometry();

            const auto& N_values = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

            mpConstitutiveLawSlave = (*mpPropSlave)[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLawSlave->InitializeMaterial( *mpPropSlave, r_geometry_slave, row(N_values , 0 ));

        } else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        KRATOS_CATCH( "" );

    }
    void SupportContact2DCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        // INITIALIZE AND RESIZE 

        const double penalty = (*mpPropMaster)[PENALTY_FACTOR];

        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size_1 = (number_of_nodes_master+number_of_nodes_slave) * 2;
        const SizeType mat_size_2 = (number_of_nodes_master+number_of_nodes_slave) * 2;
        //resizing as needed the LHS
        if(rLeftHandSideMatrix.size1() != mat_size_1 || rLeftHandSideMatrix.size2() != mat_size_2 )
            rLeftHandSideMatrix.resize(mat_size_1,mat_size_2,false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size_1,mat_size_2); //resetting LHS
        
        // resizing as needed the RHS
        if(rRightHandSideVector.size() != mat_size_1)
            rRightHandSideVector.resize(mat_size_1,false);
        noalias(rRightHandSideVector) = ZeroVector(mat_size_1); //resetting RHS

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points_master = r_geometry_master.IntegrationPoints();

        // Shape function derivatives 
        // MASTER ------------------------------------------------
        // Initialize Jacobian
        GeometryType::JacobiansType J0_master;
        // Initialize DN_DX
        const unsigned int dim = 2;
        Matrix DN_DX_master(number_of_nodes_master,2);
        Matrix InvJ0_master(dim,dim);

        // Compute the normals
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 3> tangent_physical_space;
        array_1d<double, 2> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        r_geometry_master.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT

        const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry_master.Jacobian(J0_master,this->GetIntegrationMethod());

        double DetJ0_master;

        Vector GP_parameter_coord(2); 
        GP_parameter_coord = r_geometry_master.Center();

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0_master[0](0,0);
        Jacobian(0,1) = J0_master[0](0,1);
        Jacobian(1,0) = J0_master[0](1,0);
        Jacobian(1,1) = J0_master[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0_master,DetJ0_master);


        // COMPUTE TRUE JACOBIAN DETERMINANT AND PHYSICAL NORMAL
        Vector add_factor = prod(Jacobian, tangent_parameter_space);

        DetJ0_master = norm_2(add_factor);

        normal_physical_space = prod(trans(InvJ0_master),normal_parameter_space);
        normal_physical_space /= norm_2(normal_physical_space);

        tangent_physical_space = prod(trans(InvJ0_master),tangent_parameter_space);
        tangent_physical_space /= norm_2(tangent_physical_space);

        // retrieve shape function values

        Vector old_displacement_master(2*number_of_nodes_master);
        GetValuesVector(old_displacement_master, QuadraturePointCouplingGeometry2D<Point>::Master);

        const Matrix& N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX_master) = prod(DN_De_master[0],InvJ0_master);

        // SLAVE ------------------------------------------------
        // Initialize Jacobian
        GeometryType::JacobiansType J0_slave;
        // Initialize DN_DX
        Matrix DN_DX_slave(number_of_nodes_slave,2);
        Matrix InvJ0_slave(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        
        double DetJ0_slave;
        Matrix Jacobian_slave = ZeroMatrix(2,2);
        Jacobian_slave(0,0) = J0_slave[0](0,0);
        Jacobian_slave(0,1) = J0_slave[0](0,1);
        Jacobian_slave(1,0) = J0_slave[0](1,0);
        Jacobian_slave(1,1) = J0_slave[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian_slave,InvJ0_slave,DetJ0_slave);  //DetJ0 is not the true one for NURBS geometries

         // retrieve shape function value
        Vector old_displacement_slave(2*number_of_nodes_slave);
        GetValuesVector(old_displacement_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);

        const Matrix& N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX_slave) = prod(DN_De_slave[0],InvJ0_slave);

        Vector GP_parameter_coord_on_slave(2); 
        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        GP_parameter_coord_on_slave = r_geometry_slave.Center();

        //------------------- REFERENCE WEIGHT

        const double thickness = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;

        const double IntToReferenceWeight = integration_points_master[0].Weight() * std::abs(DetJ0_master) * thickness;

        SetValue(INTEGRATION_WEIGHT, IntToReferenceWeight);

        //------------------- DB MATRICES 

        // MODIFIED
        Matrix B_master = ZeroMatrix(3,number_of_nodes_master);

        Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);

        CalculateB(B_master, DN_DX_master, number_of_nodes_master);


        CalculateB(B_slave, DN_DX_slave, number_of_nodes_slave);

        //---------- GET CONSTITUTIVE MATRICES  
        Vector stress_vector_master;
        const Matrix r_D_master = GetConstitutiveMatrix(0, B_master, r_geometry_master,
                                                         old_displacement_master, rCurrentProcessInfo, stress_vector_master); //master

        
        Vector stress_vector_slave;
        const Matrix r_D_slave  = GetConstitutiveMatrix(1, B_slave, r_geometry_slave,
                                                         old_displacement_slave, rCurrentProcessInfo, stress_vector_slave); // slave
    

        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------
        //  CHECK CONTACT CONDITION FOR THE PAIR
        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------

        // compute the normal gap
        Vector master_deformed_old(2); 
        GetDeformed(N_master, GP_parameter_coord, old_displacement_master, master_deformed_old);

        Vector slave_deformed_old(2); 
        GetDeformed(N_slave, GP_parameter_coord_on_slave, old_displacement_slave, slave_deformed_old);

        //-----------
        // COMPUTE NORMAL GAP AND CHECK CONDITION
        const double normal_gap = inner_prod(GP_parameter_coord_on_slave - GP_parameter_coord, normal_physical_space);

        // bool temp_contact_is_active = false;

        Matrix DB_master = prod(r_D_master,B_master);

        CheckCriteria(master_deformed_old, slave_deformed_old, old_displacement_master, normal_physical_space, DB_master, 
                      tangent_parameter_space, normal_physical_space, stress_vector_master); 
                
        Matrix H_master = ZeroMatrix(1, number_of_nodes_master);
        for (IndexType i = 0; i < number_of_nodes_master; ++i)
        {
            H_master(0, i)            = N_master(0, i);
        }

        Matrix H_slave = ZeroMatrix(1, number_of_nodes_slave);
        for (IndexType i = 0; i < number_of_nodes_slave; ++i)
        {
            H_slave(0, i)            = N_slave(0, i);
        }


        // Differential area
        double penalty_integration = penalty * IntToReferenceWeight;

        // Guglielmo innovaction
        double Guglielmo_innovation = 1.0;  // = 1 -> Penalty approach
                                                // = -1 -> Free-penalty approach
        if (penalty == -1.0) {
            penalty_integration = 0.0;
            Guglielmo_innovation = -1.0;
        }

        // Assembly

        if (true) {

            // MASTER
            Matrix DB_master = prod(r_D_master,B_master);
            for (IndexType i = 0; i < number_of_nodes_master; i++) {
                for (IndexType j = 0; j < number_of_nodes_master; j++) {
                    
                    for (IndexType idim = 0; idim < 2; idim++) {
                        // rLeftHandSideMatrix(2*i+idim, 2*j+idim) -= H_master(0,i)*H_master(0,j)* penalty_integration;
                        const int id1 = 2*idim;
                        const int iglob = 2*i+idim;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int id2 = (id1+2)%3;
                            const int jglob = 2*j+jdim;

                            const double sigma_n_u_kdim = (DB_master(id1, jglob)* normal_physical_space[0] + DB_master(id2, jglob)* normal_physical_space[1]);

                            const double sigma_n_u_dot_n = sigma_n_u_kdim*normal_physical_space[jdim];
                            const double sigma_n_u_dot_tau = sigma_n_u_kdim*tangent_physical_space[jdim]; // eventually for imposed forces on tangent direction

                            const double sigma_n_w_kdim = (DB_master(id1, 2*i+jdim)* normal_physical_space[0] + DB_master(id2, 2*i+jdim)* normal_physical_space[1]);

                            const double sigma_n_w_dot_n = sigma_n_w_kdim*normal_physical_space[jdim]; 
                            const double sigma_n_w_dot_tau = sigma_n_w_kdim*tangent_physical_space[jdim];//eventually for friction


                            rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*normal_physical_space[idim]
                                                                *H_master(0,j)*normal_physical_space[jdim]*penalty_integration;

                            rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_n_u_dot_n* normal_physical_space[jdim] * IntToReferenceWeight;

                            rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*H_master(0,j)* normal_physical_space[jdim]*sigma_n_w_dot_n* normal_physical_space[idim]* IntToReferenceWeight;
                        }

                    }
                }
            }

            // SLAVE
            Matrix DB_slave = prod(r_D_slave,B_slave);
            const int shift = number_of_nodes_master*2;
            for (IndexType i = 0; i < number_of_nodes_master; i++) {
                for (IndexType j = 0; j < number_of_nodes_slave; j++) {
                    
                    for (IndexType idim = 0; idim < 2; idim++) {
                        const int id1 = 2*idim;
                        const int iglob = 2*i+idim;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int id2 = (id1+2)%3;
                            const int jglob = 2*j+jdim + shift;

                            double sigma_n_w_dot_n = 0.0;
                            double sigma_n_w_dot_tau = 0.0;
                            for (IndexType kdim = 0; kdim < 2; kdim++) {
                                const double sigma_n_w_kdim = (DB_master(id1, 2*i+kdim)* normal_physical_space[0] + DB_master(id2, 2*i+kdim)* normal_physical_space[1]);

                                sigma_n_w_dot_n += sigma_n_w_kdim*normal_physical_space[kdim]; 
                                sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[kdim];//eventually for friction
                            }


                            rLeftHandSideMatrix(iglob, jglob) += H_master(0,i)*normal_physical_space[idim]
                                                                *H_slave(0,j)*normal_physical_space[jdim]*penalty_integration;

                            rLeftHandSideMatrix(iglob, jglob) += Guglielmo_innovation*H_slave(0,j)* normal_physical_space[jdim]*sigma_n_w_dot_n* normal_physical_space[idim]* IntToReferenceWeight;

                        }

                    }
                }
            }
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // INTEGRATION ON SLAVE
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // ONLY SLAVE (substitute /sigma_1 *n  =sigma_2*n)
            for (IndexType i = 0; i < number_of_nodes_slave; i++) {
                for (IndexType j = 0; j < number_of_nodes_master; j++) {
                    
                    for (IndexType idim = 0; idim < 2; idim++) {
                        const int id1 = 2*idim;
                        const int iglob = 2*i+idim + shift;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int id2 = (id1+2)%3;
                            const int jglob = 2*j+jdim;

                            
                            double sigma_n_u_dot_n = 0.0;
                            double sigma_n_u_dot_tau = 0.0;
                            for (IndexType kdim = 0; kdim < 2; kdim++) {
                                const int kglob = 2*j+kdim;
                                const double sigma_n_u_kdim = -(DB_master(id1, kglob)* normal_physical_space[0] + DB_master(id2, kglob)* normal_physical_space[1]);

                                sigma_n_u_dot_n += sigma_n_u_kdim*normal_physical_space[kdim];
                                sigma_n_u_dot_tau += sigma_n_u_kdim*tangent_physical_space[kdim]; // eventually for imposed forces on tangent direction
                            }

                            // i am imposing zero in the other direction

                            rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*sigma_n_u_dot_n* normal_physical_space[jdim] * IntToReferenceWeight;
                            
                        }

                    }
                }
            }

            //#################
            // MASTER
            // for (IndexType i = 0; i < number_of_nodes_slave; i++) {
            //     for (IndexType j = 0; j < number_of_nodes_master; j++) {
                    
            //         for (IndexType idim = 0; idim < 2; idim++) {
            //             // rLeftHandSideMatrix(2*i+idim, 2*j+idim) -= H_master(0,i)*H_master(0,j)* penalty_integration;
            //             const int id1 = 2*idim;
            //             const int iglob = 2*i+idim + shift;

            //             for (IndexType jdim = 0; jdim < 2; jdim++) {
            //                 const int id2 = (id1+2)%3;
            //                 const int jglob = 2*j+jdim;

            //                 const double sigma_n_w_kdim = (DB_slave(id1, 2*i+jdim)* normal_physical_space[0] + DB_slave(id2, 2*i+jdim)* normal_physical_space[1]);

            //                 const double sigma_n_w_dot_n = sigma_n_w_kdim*normal_physical_space[jdim]; 
            //                 const double sigma_n_w_dot_tau = sigma_n_w_kdim*tangent_physical_space[jdim];//eventually for friction


            //                 rLeftHandSideMatrix(iglob, jglob) += H_slave(0,i)*normal_physical_space[idim]
            //                                                     *H_master(0,j)*normal_physical_space[jdim]*penalty_integration;

            //                 rLeftHandSideMatrix(iglob, jglob) += Guglielmo_innovation*H_master(0,j)* normal_physical_space[jdim]*sigma_n_w_dot_n* normal_physical_space[idim]* IntToReferenceWeight;
            //             }

            //         }
            //     }
            // }

            // // SLAVE
            // for (IndexType i = 0; i < number_of_nodes_slave; i++) {
            //     for (IndexType j = 0; j < number_of_nodes_slave; j++) {
                    
            //         for (IndexType idim = 0; idim < 2; idim++) {
            //             const int id1 = 2*idim;
            //             const int iglob = 2*i+idim + shift;

            //             for (IndexType jdim = 0; jdim < 2; jdim++) {
            //                 const int id2 = (id1+2)%3;
            //                 const int jglob = 2*j+jdim + shift;

            //                 double sigma_n_w_dot_n = 0.0;
            //                 double sigma_n_w_dot_tau = 0.0;
            //                 for (IndexType kdim = 0; kdim < 2; kdim++) {
            //                     const double sigma_n_w_kdim = (DB_slave(id1, 2*i+kdim)* normal_physical_space[0] + DB_slave(id2, 2*i+kdim)* normal_physical_space[1]);

            //                     sigma_n_w_dot_n += sigma_n_w_kdim*normal_physical_space[kdim]; 
            //                     sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[kdim];//eventually for friction
            //                 }


            //                 rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*normal_physical_space[idim]
            //                                                     *H_slave(0,j)*normal_physical_space[jdim]*penalty_integration;

            //                 rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*H_slave(0,j)* normal_physical_space[jdim]*sigma_n_w_dot_n* normal_physical_space[idim]* IntToReferenceWeight;

            //             }

            //         }
                // }
            // }
            //#####################3


            if (CalculateResidualVectorFlag) {
                
                Vector gn = ZeroVector(2); //->GetValue(DISPLACEMENT);

                gn[0] = normal_gap*normal_physical_space[0];
                gn[1] = normal_gap*normal_physical_space[1];

                for (IndexType i = 0; i < number_of_nodes_master; i++) {

                    for (IndexType idim = 0; idim < 2; idim++) {

                        rRightHandSideVector[2*i+idim] -= H_master(0,i)*gn[idim]* penalty_integration;
                        const int id1 = idim*2;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int id2 = (id1+2)%3;

                            double sigma_n_w_dot_n = 0.0;
                            double sigma_n_w_dot_tau = 0.0;

                            const double sigma_n_w_kdim = (DB_master(id1, 2*i+jdim)* normal_physical_space[0] + DB_master(id2, 2*i+jdim)* normal_physical_space[1]);

                            sigma_n_w_dot_n += sigma_n_w_kdim*normal_physical_space[jdim]; 
                            sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[jdim];//eventually for friction


                            rRightHandSideVector(2*i+idim) -= Guglielmo_innovation*gn[jdim]*sigma_n_w_dot_n * normal_physical_space[idim] * IntToReferenceWeight;
                        }

                    }
                }
            }

            
            Vector temp = ZeroVector(mat_size_2);

            GetValuesVector(temp, 2);

            // RHS = ExtForces - K*temp;
            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
            // }
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // INTEGRATION ON SLAVE
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // else if (integrationDomain == 1) // integration on slave
            // {
            //     // ONLY SLAVE (substitute /sigma_1 *n  =sigma_2*n)
            //     Matrix DB_master = prod(r_D_master,B_master);
            //     Matrix DB_slave = prod(r_D_slave,B_slave);
            //     const int shift = number_of_nodes_master*2;
            //     for (IndexType i = 0; i < number_of_nodes_master; i++) {
            //         for (IndexType j = 0; j < number_of_nodes_slave; j++) {
                        
            //             for (IndexType idim = 0; idim < 2; idim++) {
            //                 const int id1 = 2*idim;
            //                 const int iglob = 2*i+idim;

            //                 for (IndexType jdim = 0; jdim < 2; jdim++) {
            //                     const int id2 = (id1+2)%3;
            //                     const int jglob = 2*j+jdim + shift;

                                
            //                     double sigma_n_u_dot_n = 0.0;
            //                     double sigma_n_u_dot_tau = 0.0;
            //                     for (IndexType kdim = 0; kdim < 2; kdim++) {
            //                         const int kglob = 2*j+kdim;
            //                         const double sigma_n_u_kdim = (DB_slave(id1, kglob)* normal_physical_space[0] + DB_slave(id2, kglob)* normal_physical_space[1]);

            //                         sigma_n_u_dot_n += sigma_n_u_kdim*normal_physical_space[kdim];
            //                         sigma_n_u_dot_tau += sigma_n_u_kdim*tangent_physical_space[kdim]; // eventually for imposed forces on tangent direction
            //                     }

            //                     // i am imposing zero in the other direction

            //                     rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_n_u_dot_n* normal_physical_space[jdim] * IntToReferenceWeight;
                                
            //                 }

            //             }
            //         }
            //     }


                // if (CalculateResidualVectorFlag) {
                    
                //     Vector fn = ZeroVector(2); 

                //     fn[0] = 0.0;
                //     fn[1] = 0.0;

                //     for (IndexType i = 0; i < number_of_nodes_master; i++) {

                //         for (IndexType idim = 0; idim < 2; idim++) {

                //             rRightHandSideVector(2*i+idim) -= H_master(0,i)*fn[idim]* IntToReferenceWeight;

                //         }
                //     }
                    
                //     Vector temp = ZeroVector(mat_size_2);

                //     GetValuesVector(temp, 2);

                //     // RHS = ExtForces - K*temp;
                //     noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
                // } 
            // }
            // KRATOS_WATCH("!!!!!!!! ACTIVE CONTACT")
        } else { // CONTACT NOT ACTIVE FOR THE PAIR

            // KRATOS_WATCH("NOT ACTIVE CONTACT")
        }
        KRATOS_CATCH("")
    }

    int SupportContact2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT((*mpPropMaster).Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SupportContact2DCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType number_of_nodes = (number_of_nodes_master + number_of_nodes_slave);

        if (rResult.size() != 2 * number_of_nodes)
            rResult.resize(2 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry_master[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }

        const int shift = number_of_nodes_master*2;
        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry_slave[i];
            rResult[index+shift] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1+shift] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void SupportContact2DCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * (number_of_nodes_master+number_of_nodes_slave));

        for (IndexType i = 0; i < number_of_nodes_master; ++i) {
            const auto& r_node = r_geometry_master[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }

        for (IndexType i = 0; i < number_of_nodes_slave; ++i) {
            const auto& r_node = r_geometry_slave[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }

    };


    void SupportContact2DCondition::GetValuesVector(
        Vector& rValues, IndexType index) const
    {
        if (index == QuadraturePointCouplingGeometry2D<Point>::Master) {
            
            SizeType number_of_control_points = GetMasterGeometry().size();
            const SizeType mat_size = number_of_control_points * 2;

            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points; ++i)
            {
                const array_1d<double, 2 >& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        } else if (index == QuadraturePointCouplingGeometry2D<Point>::Slave) {
            SizeType number_of_control_points = GetSlaveGeometry().size();
            const SizeType mat_size = number_of_control_points * 2;

            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points; ++i)
            {
                const array_1d<double, 2 >& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        } else if (index == 2) {
            
            SizeType number_of_control_points_master = GetMasterGeometry().size();
            const SizeType mat_size_master = number_of_control_points_master * 2;
            SizeType number_of_control_points_slave = GetSlaveGeometry().size();
            const SizeType mat_size_slave = number_of_control_points_slave * 2;
            const SizeType mat_size = mat_size_master + mat_size_slave;

            if (rValues.size() != mat_size)
                rValues.resize(mat_size, false);

            for (IndexType i = 0; i < number_of_control_points_master; ++i)
            {
                const array_1d<double, 2 >& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }

            for (IndexType i = 0; i < number_of_control_points_slave; ++i)
            {
                const array_1d<double, 2 >& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2 + mat_size_master;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        }
    }

    void SupportContact2DCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX,
        const SizeType number_of_control_points) const
    {
        
        const SizeType mat_size = number_of_control_points * 2;

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 2;
            IndexType dirr = r % 2;

            rB(0, r) = r_DN_DX(kr,0) * (1-dirr);
            rB(1, r) = r_DN_DX(kr,1) * dirr;
            rB(2, r) = r_DN_DX(kr,0) * (dirr) + r_DN_DX(kr,1) * (1-dirr);
        }
    }


    void SupportContact2DCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetMasterGeometry(), *mpPropMaster, rCurrentProcessInfo);

        mpConstitutiveLawMaster->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        
    /////////////////////////
    //     const auto& r_geometry = GetMasterGeometry();
    //     const SizeType nb_nodes = r_geometry.size();

    //     // Integration Points
    //     const GeometryType::IntegrationPointsArrayType& integration_points_master = r_geometry.IntegrationPoints();
    //     // Shape function values
    //     const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    //     GeometryType::JacobiansType J0;
    //     r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    //     // Get the parameter coordinates
    //     Vector GP_parameter_coord(2); 
    //     GP_parameter_coord = r_geometry.Center(); // Only one Integration Points 

    //     double x_coord_gauss_point = 0;
    //     double y_coord_gauss_point = 0;
    //     double rOutput_x = 0;
    //     double rOutput_y = 0;

    //     for (IndexType i = 0; i < nb_nodes; ++i)
    //     {
    //         // KRATOS_WATCH(r_geometry[i])
    //         double output_solution_step_value_x = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_X);
    //         double output_solution_step_value_y = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_Y);
    //         rOutput_x += r_N(0, i) * output_solution_step_value_x;
    //         rOutput_y += r_N(0, i) * output_solution_step_value_y;
    //         x_coord_gauss_point += r_N(0, i) * r_geometry[i].X0();
    //         y_coord_gauss_point += r_N(0, i) * r_geometry[i].Y0();
    //     }   

    // //////////////////////////
    // /////////////// SLAVE 
    // ////////////////////////
    //     const auto& r_geometry_slave = GetSlaveGeometry();
    //     const SizeType nb_nodes_slave = r_geometry_slave.size();

    //     // Integration Points
    //     const GeometryType::IntegrationPointsArrayType& integration_points_slave = r_geometry_slave.IntegrationPoints();
    //     // Shape function values
    //     const Matrix& r_N_slave = r_geometry_slave.ShapeFunctionsValues();

    //     GeometryType::JacobiansType J0_slave;
    //     r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
    //     // Get the parameter coordinates
    //     Vector GP_parameter_coord_slave(2); 
    //     GP_parameter_coord_slave = r_geometry_slave.Center(); // Only one Integration Points 

    //     double x_coord_gauss_point_slave = 0;
    //     double y_coord_gauss_point_slave = 0;
    //     double rOutput_slave_x = 0;
    //     double rOutput_slave_y = 0;

    //     for (IndexType i = 0; i < nb_nodes_slave; ++i)
    //     {
    //         // KRATOS_WATCH(r_geometry[i])
    //         double output_solution_step_value_x = r_geometry_slave[i].GetSolutionStepValue(DISPLACEMENT_X);
    //         double output_solution_step_value_y = r_geometry_slave[i].GetSolutionStepValue(DISPLACEMENT_Y);
    //         rOutput_slave_x += r_N_slave(0, i) * output_solution_step_value_x;
    //         rOutput_slave_y += r_N_slave(0, i) * output_solution_step_value_y;
    //         x_coord_gauss_point_slave += r_N_slave(0, i) * r_geometry_slave[i].X0();
    //         y_coord_gauss_point_slave += r_N_slave(0, i) * r_geometry_slave[i].Y0();
    //     }     

    //     if (m_contact_is_active) {
    //     // if (true) {
    //         if (x_coord_gauss_point <= 2.005 || GetGeometry().GetValue(ACTIVATION_LEVEL) == 0) {
    //             std::ofstream gap_file("txt_files/gap_output_master.txt", std::ios::app);
    //             if (gap_file.is_open()) {
    //                 gap_file << std::scientific << std::setprecision(14); // Set precision to 10^-14

    //                 const double x_final_pos_master = x_coord_gauss_point + rOutput_x;
    //                 const double x_final_pos_slave = x_coord_gauss_point_slave + rOutput_slave_x;

    //                 // gap_file << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point + rOutput_x << " " << y_coord_gauss_point + rOutput_y << " " <<integration_points_master[0].Weight() << std::endl;
    //                 gap_file << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points_master[0].Weight() << std::endl;
    //                 gap_file.close();

    //                 // KRATOS_WATCH(y_coord_gauss_point+rOutput_y)
    //                 // KRATOS_WATCH(y_coord_gauss_point_slave + rOutput_slave_y)
    //             }  

    //             std::ofstream gap_file2("txt_files/gap_output_slave.txt", std::ios::app);
    //             if (gap_file2.is_open()) {
    //                 gap_file2 << std::scientific << std::setprecision(14); // Set precision to 10^-14

    //                 const double x_final_pos_master = x_coord_gauss_point + rOutput_x;
    //                 const double x_final_pos_slave = x_coord_gauss_point_slave + rOutput_slave_x;

    //                 const double y_final_pos_slave = y_coord_gauss_point_slave + rOutput_slave_y;


    //                 // gap_file2 << x_final_pos_master-x_final_pos_slave << " " << x_final_pos_slave << " " << y_final_pos_slave << " " <<integration_points_master[0].Weight() << std::endl;
    //                 gap_file2 << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point_slave << " " << y_coord_gauss_point_slave << " " <<integration_points_master[0].Weight() << std::endl;
    //                 gap_file2.close();
    //             }  
    //         }
    //         else if (x_coord_gauss_point >2.005 || GetGeometry().GetValue(ACTIVATION_LEVEL) == 1){
    //             std::ofstream gap_file("txt_files/gap_output_slave.txt", std::ios::app);
    //             if (gap_file.is_open()) {
    //                 gap_file << std::scientific << std::setprecision(14); // Set precision to 10^-14

    //                 const double x_final_pos_master = x_coord_gauss_point + rOutput_x;
    //                 const double x_final_pos_slave = x_coord_gauss_point_slave + rOutput_slave_x;

    //                 // gap_file << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point + rOutput_x << " " << y_coord_gauss_point + rOutput_y << " " <<integration_points_master[0].Weight() << std::endl;
    //                 gap_file << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points_master[0].Weight() << std::endl;
    //                 gap_file.close();

    //                 // KRATOS_WATCH(y_coord_gauss_point+rOutput_y)
    //                 // KRATOS_WATCH(y_coord_gauss_point_slave + rOutput_slave_y)
    //             }  

    //             std::ofstream gap_file2("txt_files/gap_output_master.txt", std::ios::app);
    //             if (gap_file2.is_open()) {
    //                 gap_file2 << std::scientific << std::setprecision(14); // Set precision to 10^-14

    //                 const double x_final_pos_master = x_coord_gauss_point + rOutput_x;
    //                 const double x_final_pos_slave = x_coord_gauss_point_slave + rOutput_slave_x;

    //                 const double y_final_pos_slave = y_coord_gauss_point_slave + rOutput_slave_y;


    //                 // gap_file2 << x_final_pos_master-x_final_pos_slave << " " << x_final_pos_slave << " " << y_final_pos_slave << " " <<integration_points_master[0].Weight() << std::endl;
                    
    //                 gap_file2 << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point_slave << " " << y_coord_gauss_point_slave << " " <<integration_points_master[0].Weight() << std::endl;
                    
    //                 gap_file2.close();
    //             }  
    //         }

    //         std::ofstream disp_file("txt_files/contact_displacement.txt", std::ios::app);

    //         disp_file << x_coord_gauss_point << " " << x_coord_gauss_point + rOutput_x << " " << y_coord_gauss_point << " " << y_coord_gauss_point + rOutput_y << std::endl;

    //         disp_file << x_coord_gauss_point_slave << " " << x_coord_gauss_point_slave + rOutput_slave_x << " " << y_coord_gauss_point_slave << " " << y_coord_gauss_point_slave + rOutput_slave_y << std::endl;
        
    //     }

    //     if (x_coord_gauss_point <= 2.0) {
    //     std::ofstream output_file("txt_files/output_results_GPs_master.txt", std::ios::app);
    //         if (output_file.is_open()) {
    //             output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
    //             output_file << rOutput_x << " " << rOutput_y << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points_master[0].Weight() << std::endl;
    //             output_file.close();
    //         } 
    //     } else {
    //         std::ofstream output_file("txt_files/output_results_GPs_slave.txt", std::ios::app);
    //         if (output_file.is_open()) {
    //             output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
    //             output_file << rOutput_x << " " << rOutput_y << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points_master[0].Weight() << std::endl;
    //             output_file.close();
    //         } 
    //     }


    //     std::ofstream outputFile("txt_files/Gauss_Point_coordinates.txt", std::ios::app);
    //     if (!outputFile.is_open())
    //     {
    //         std::cerr << "Failed to open the file for writing." << std::endl;
    //         return;
    //     }
    //     outputFile << std::setprecision(14); // Set precision to 10^-14
    //     outputFile << x_coord_gauss_point << "  " << y_coord_gauss_point <<"\n";
    //     outputFile.close();


        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        const auto& r_geometry = GetMasterGeometry();
        const SizeType nb_nodes = r_geometry.size();
        const SizeType mat_size = nb_nodes * 2;

        // Shape function derivatives (NEW) 
        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        const unsigned int dim = 2;
        Matrix DN_DX(nb_nodes,2);
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
        // MODIFIED
        Vector old_displacement(mat_size);
        GetValuesVector(old_displacement, 0);
        
        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        // MODIFIED
        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX,nb_nodes);

        normal_physical_space = prod(trans(InvJ0),normal_parameter_space);

        normal_physical_space /= norm_2(normal_physical_space);

        // GET STRESS VECTOR
        Vector stress_vector_master;
        const Matrix r_D_master = GetConstitutiveMatrix(0, B, r_geometry,
                                                         old_displacement, rCurrentProcessInfo, stress_vector_master); //master

        Vector sigma_n(2);

        sigma_n[0] = stress_vector_master[0]*normal_physical_space[0] + stress_vector_master[2]*normal_physical_space[1];
        sigma_n[1] = stress_vector_master[2]*normal_physical_space[0] + stress_vector_master[1]*normal_physical_space[1];

        //-----------------------------------------
        // Vector sigma_n = ZeroVector(2);
        // Matrix DB = prod(r_D_master,B);
        // for (IndexType j = 0; j < r_geometry.size(); j++) {

        //     for (IndexType jdim = 0; jdim < 2; jdim++) {
        //         const int jglob = 2*j+jdim;

        //         sigma_n[0] += (DB(0, jglob)* normal_physical_space[0] + DB(2, jglob)* normal_physical_space[1])*old_displacement[jglob];
        //         sigma_n[1] += (DB(2, jglob)* normal_physical_space[0] + DB(1, jglob)* normal_physical_space[1])*old_displacement[jglob];
        //     }
        // }
        //2222222222222222222222222222222222222222222222222222222222222222222222222

        SetValue(NORMAL_STRESS, sigma_n);
        // //---------------------
    }

    void SupportContact2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

        // InitializeMaterial();
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetMasterGeometry(), (*mpPropMaster), rCurrentProcessInfo);

        mpConstitutiveLawMaster->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

    }

    //----------------------------------------------------------------------------------
    // void SupportContact2DCondition::GetDBMatrix(IndexType index, GeometryType::Pointer r_geometry, 
    //                                             Vector& old_displacement, Matrix& DB, const Kratos::ProcessInfo& rCurrentProcessInfo) {
        
    //     const int number_of_nodes = r_geometry->size();
    //     if (old_displacement.size() != number_of_nodes) old_displacement.resize(number_of_nodes);

    //     //  COMPUTE B MATRIX


    // }


    


    //----------------------------------------------------------------------------------
    const Matrix SupportContact2DCondition::GetConstitutiveMatrix(IndexType index, Matrix& r_B, GeometryType r_geometry,
                                                                   Vector& old_displacement, const Kratos::ProcessInfo& rCurrentProcessInfo,
                                                                   Vector& stress_vector) {

        ConstitutiveLaw::Pointer rpConstitutiveLaw = GetConstitutiveLaw(index);

        PropertiesType r_prop = GetProperty(index);

        ConstitutiveLaw::Parameters Values(r_geometry, r_prop, rCurrentProcessInfo);

        const SizeType strain_size = rpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(r_B,old_displacement);
    
        // Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);

        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        rpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        
        stress_vector = Values.GetStressVector();

        return Values.GetConstitutiveMatrix();
        
    }



    void SupportContact2DCondition::GetDeformed(const Matrix& N, Vector& reference_position, Vector& displacement, Vector& deformed_position) {
        
        // compute the old displacement for the two points in the pair
        double displacement_x = 0.0; double displacement_y = 0.0;
        
        for (IndexType i = 0; i < GetMasterGeometry().size(); ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value_x = displacement[2*i];
            double output_solution_step_value_y = displacement[2*i+1];
            displacement_x += N(0, i) * output_solution_step_value_x;
            displacement_y += N(0, i) * output_solution_step_value_y;
        } 
        deformed_position[0] = reference_position[0] + displacement_x; 
        deformed_position[1] = reference_position[1] + displacement_y;
    }

    //--------------- CHECK CRITERIA
    bool SupportContact2DCondition::CheckCriteria(const Vector& deformed_pos_master, 
                                                  const Vector& deformed_pos_slave, const Vector& displacement_master,
                                                  const Vector& normal, const Matrix& DB_master,
                                                  array_1d<double, 3>& local_tangent, array_1d<double, 2>& old_normal,
                                                  const Vector stress_vector_master){
        const double toll = 1e-12;

        const double toll_tangential_distance = 1e-1;

        double flag_normal_gap = inner_prod(deformed_pos_slave - deformed_pos_master, normal);

        double tangential_gap = norm_2((deformed_pos_slave - deformed_pos_master) - flag_normal_gap * normal);


        // COMPUTE THE UPDATED NORMAL
        //need du/dt = du/dx dx_curve/dt + du/dx dx_curve/dt (for each component)
        int nNodesMaster = GetMasterGeometry().size();
        Vector old_displacement_master(nNodesMaster*2);
        GetValuesVector(old_displacement_master, QuadraturePointCouplingGeometry2D<Point>::Master);

        Vector new_normal(3);
        new_normal = old_normal;
        if ( (norm_2(old_displacement_master) > 1e-12 && GetGeometry().GetValue(ACTIVATION_LEVEL) == 0)) {

            const int derivative_order = 1;
            Matrix shapeFunctionDerivatives_master = GetMasterGeometry().ShapeFunctionDerivatives(derivative_order, 0, this->GetIntegrationMethod());

            //compute du/dx 
            double dux_dx = 0.0; double dux_dy = 0.0;
            double duy_dx = 0.0; double duy_dy = 0.0;
            for (int i = 0; i < nNodesMaster; i++) {
                dux_dx += old_displacement_master[2*i]*shapeFunctionDerivatives_master(i, 0);
                dux_dy += old_displacement_master[2*i]*shapeFunctionDerivatives_master(i ,1);

                duy_dx += old_displacement_master[2*i+1]*shapeFunctionDerivatives_master(i, 0);
                duy_dy += old_displacement_master[2*i+1]*shapeFunctionDerivatives_master(i ,1);
            }

            Vector du_dt = ZeroVector(3);
            du_dt[0] = dux_dx * local_tangent[0] + dux_dy * local_tangent[1];
            du_dt[1] = duy_dx * local_tangent[0] + duy_dy * local_tangent[1];

            // remove and make it general !!!!!!!!!!! (ONLY WORKS FOR THE SQUARE/RECTANGULAR CASE)
            double magnitude = std::sqrt(local_tangent[0] * local_tangent[0] + local_tangent[1] * local_tangent[1]);
            //----
            Vector normal_displacement(3);

            Vector out_of_plane = ZeroVector(3); out_of_plane[2] = 1.0;

            // KRATOS_WATCH(deformed_pos_master)
            // KRATOS_WATCH(duy_dx)

            MathUtils<double>::CrossProduct(normal_displacement, du_dt, out_of_plane);

            // KRATOS_WATCH(normal_displacement)

            
            // new_normal = old_normal *magnitude  + normal_displacement;

            // new_normal = new_normal / norm_2(new_normal);

            new_normal = old_normal;

            // KRATOS_WATCH(new_normal)

            std::ofstream output_file("txt_files/new_normal.txt", std::ios::app);
            output_file << std::scientific << std::setprecision(14);

            output_file << deformed_pos_master[0] << " " << deformed_pos_master[1] << " " <<
                           new_normal[0] << " " << new_normal[1] << std::endl;

            // exit(0);
        }




        // COMPUTE TRUE NORMAL STRESS

        double true_normal_stress = 0.0;
        // Vector sigma_n = ZeroVector(2);
        // for (IndexType j = 0; j < GetMasterGeometry().size(); j++) {

        //     for (IndexType jdim = 0; jdim < 2; jdim++) {
        //         const int jglob = 2*j+jdim;

        //         sigma_n[0] += (DB_master(0, jglob)* new_normal[0] + DB_master(2, jglob)* new_normal[1])*displacement_master[jglob];
        //         sigma_n[1] += (DB_master(2, jglob)* new_normal[0] + DB_master(1, jglob)* new_normal[1])*displacement_master[jglob];
        //     }
        // }
        // true_normal_stress = sigma_n[0]*new_normal[0] + sigma_n[1]*new_normal[1];

        true_normal_stress = (stress_vector_master[0]* new_normal[0] + stress_vector_master[2]* new_normal[1])*new_normal[0] +
                                      (stress_vector_master[2]* new_normal[0] + stress_vector_master[1]* new_normal[1])*new_normal[1];


        // if (integrationDomain == 0) true_normal_stress = -true_normal_stress;
        if (-true_normal_stress -(*mpPropMaster)[YOUNG_MODULUS]*flag_normal_gap > toll) {
        // if (-(*mpPropMaster)[YOUNG_MODULUS]*flag_normal_gap > toll) {
            // if (tangential_gap > toll_tangential_distance) 
            // {
            //     m_contact_is_active = false;
            //     return false;
            // }

            // KRATOS_WATCH(displacement_master)

            m_contact_is_active = true;

            return true;
 
        } //else if (GetMasterGeometry().Center()[1] < 1.5001 && GetMasterGeometry().Center()[1] > 0.49999) {
        

        //     KRATOS_WATCH("nnnnnnnoooooooooooooooooooooo") 
        //     KRATOS_WATCH(true_normal_stress)
        // }

        m_contact_is_active = false;
        return false;
    }

} // Namespace Kratos
