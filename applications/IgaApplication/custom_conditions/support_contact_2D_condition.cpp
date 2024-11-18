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

        // this->Set(ACTIVE, true);

        SetValue(ACTIVATION_LEVEL, 0);
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

        // //----------------------------------------------
        // {
        //     double normal_gap = this->GetValue(NORMAL_GAP);
        //     Vector normal_stress = this->GetValue(STRESS_MASTER);
        //     Vector normal = this->GetValue(NORMAL_MASTER);
        //     double weight = this->GetValue(INTEGRATION_WEIGHT);
 
        //     // KRATOS_WATCH(i_cond->GetProperties()[YOUNG_MODULUS])
        //     double yound_modulus = 10.0;


        //     // double true_normal_stress = inner_prod(normal_stress, normal);

        //     double true_normal_stress = (normal_stress[0]* normal[0] + normal_stress[2]* normal[1])*normal[0] +
        //                               (normal_stress[2]* normal[0] + normal_stress[1]* normal[1])*normal[1];

        //     double check_value = -(true_normal_stress+yound_modulus*normal_gap);

        //     if (check_value > 0) SetValue(ACTIVATION_LEVEL, 1);
        //     else SetValue(ACTIVATION_LEVEL, 0);
        // }
        //----------------------------------------------------------------------

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

        const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry_master.Jacobian(J0_master,this->GetIntegrationMethod());

        double DetJ0_master;

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0_master[0](0,0); Jacobian(0,1) = J0_master[0](0,1);
        Jacobian(1,0) = J0_master[0](1,0); Jacobian(1,1) = J0_master[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0_master,DetJ0_master);

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
        Jacobian_slave(0,0) = J0_slave[0](0,0); Jacobian_slave(0,1) = J0_slave[0](0,1);
        Jacobian_slave(1,0) = J0_slave[0](1,0); Jacobian_slave(1,1) = J0_slave[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian_slave,InvJ0_slave,DetJ0_slave);  //DetJ0 is not the true one for NURBS geometries

        const Matrix& N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX_slave) = prod(DN_De_slave[0],InvJ0_slave);

        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        

        //------------------- REFERENCE WEIGHT

        const double thickness = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;

        const double IntToReferenceWeight = GetValue(INTEGRATION_WEIGHT);

        const Vector normal_physical_space_3D = GetValue(NORMAL_MASTER); SetValue(NORMAL, normal_physical_space_3D);

        Vector normal_physical_space(2); normal_physical_space[0] = normal_physical_space_3D[0]; normal_physical_space[1] = normal_physical_space_3D[1];

        //------------------- DB MATRICES 

        // MODIFIED
        Matrix B_master = ZeroMatrix(3,number_of_nodes_master);

        Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);

        CalculateB(B_master, DN_DX_master, number_of_nodes_master);


        CalculateB(B_slave, DN_DX_slave, number_of_nodes_slave);

        //---------- GET CONSTITUTIVE MATRICES  
        Vector stress_vector_master = GetValue(STRESS_MASTER);
        const Matrix r_D_master = GetValue(CONSTITUTIVE_MATRIX_MASTER);

        
        Vector stress_vector_slave = GetValue(STRESS_SLAVE);
        const Matrix r_D_slave  = GetValue(CONSTITUTIVE_MATRIX_SLAVE);
    

        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------
        //  CHECK CONTACT CONDITION FOR THE PAIR
        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------

        // compute the normal gap
        Vector master_deformed_old(3); 
        Vector GP_parameter_coord = r_geometry_master.Center();
        master_deformed_old = GP_parameter_coord + GetValue(DISPLACEMENT_MASTER);

        Vector slave_deformed_old(3); 
        Vector GP_parameter_coord_on_slave = r_geometry_slave.Center();
        slave_deformed_old = GP_parameter_coord_on_slave + GetValue(DISPLACEMENT_SLAVE);

        //-----------
        // COMPUTE NORMAL GAP AND CHECK CONDITION
        const double normal_gap = inner_prod(GP_parameter_coord_on_slave - GP_parameter_coord, normal_physical_space_3D);

        // bool temp_contact_is_active = false;

        Matrix DB_master = prod(r_D_master,B_master);
                
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
        // if (this->Is(ACTIVE)) {
        // KRATOS_WATCH(this->GetValue(ACTIVATION_LEVEL))
        // KRATOS_WATCH(GetValue(ACTIVATION_LEVEL))
        if (this->GetValue(ACTIVATION_LEVEL) == 1) {

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

                            // // PENALTY FREE g_n = 0
                            // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                            // //*********************************************** */
                            Vector sigma_w_n(2);
                            sigma_w_n[0] = (DB_master(0, iglob)* normal_physical_space[0] + DB_master(2, iglob)* normal_physical_space[1]);
                            sigma_w_n[1] = (DB_master(2, iglob)* normal_physical_space[0] + DB_master(1, iglob)* normal_physical_space[1]);

                            // double sigma_w_n_dot_n = sigma_w_n[idim] * normal_physical_space[idim];

                            double sigma_w_n_dot_n = inner_prod(sigma_w_n, normal_physical_space);

                            // u_1 \dot n = H_master[j] * n[jdim]
                            //----------------------------------
                            rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*sigma_w_n_dot_n * H_master(0,j) * normal_physical_space[jdim] * IntToReferenceWeight;
                            

                            // // PENALTY g_n = 0
                            // //*********************************************** */
                            rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*normal_physical_space[idim]
                                                                *H_master(0,j)*normal_physical_space[jdim]*penalty_integration;

                            // FLUX 
                            // [sigma_1(u) \dot n] \dot n * (-w_1 \dot n)
                            //*********************************************** */
                            Vector sigma_u_n(2);
                            sigma_u_n[0] = (DB_master(0, jglob)* normal_physical_space[0] + DB_master(2, jglob)* normal_physical_space[1]);
                            sigma_u_n[1] = (DB_master(2, jglob)* normal_physical_space[0] + DB_master(1, jglob)* normal_physical_space[1]);

                            double sigma_u_n_dot_n = inner_prod(sigma_u_n, normal_physical_space);

                            // w_1 \dot n = H_master[i] * n[idim]
                            //----------------------------------
                            rLeftHandSideMatrix(iglob, jglob) -= sigma_u_n_dot_n * H_master(0,i) * normal_physical_space[idim] * IntToReferenceWeight;


                //##########################################################################33
                            // const double sigma_n_u_kdim = (DB_master(id1, jglob)* normal_physical_space[0] + DB_master(id2, jglob)* normal_physical_space[1]);

                            // const double sigma_n_u_dot_n = sigma_n_u_kdim*normal_physical_space[jdim];
                            // // const double sigma_n_u_dot_tau = sigma_n_u_kdim*tangent_physical_space[jdim]; // eventually for imposed forces on tangent direction

                            // // sigma(w) \dot n
                            // const double sigma_n_w_kdim = (DB_master(id1, 2*i+jdim)* normal_physical_space[0] + DB_master(id2, 2*i+jdim)* normal_physical_space[1]);

                            // // (sigma(w) \dot n) \dot n 
                            // const double sigma_n_w_dot_n = sigma_n_w_kdim*normal_physical_space[jdim]; 

                            // // const double sigma_n_w_dot_tau = sigma_n_w_kdim*tangent_physical_space[jdim];//eventually for friction


                            // // rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*normal_physical_space[idim]
                            // //                                     *H_master(0,j)*normal_physical_space[jdim]*penalty_integration;

                            // // rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_n_u_dot_n* normal_physical_space[jdim] * IntToReferenceWeight;

                            // rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*H_master(0,j)* normal_physical_space[jdim]*sigma_n_w_dot_n* normal_physical_space[idim]* IntToReferenceWeight;
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

                            // PENALTY FREE g_n = 0
                            // [\sigma_1(w) \dot n] \dot n (+u_2 \dot n)
                            //*********************************************** */
                            Vector sigma_w_n(2);
                            sigma_w_n[0] = (DB_master(0, iglob)* normal_physical_space[0] + DB_master(2, iglob)* normal_physical_space[1]);
                            sigma_w_n[1] = (DB_master(2, iglob)* normal_physical_space[0] + DB_master(1, iglob)* normal_physical_space[1]);

                            // double sigma_w_n_dot_n = sigma_w_n[idim] * normal_physical_space[idim];; 
                            
                            double sigma_w_n_dot_n = inner_prod(sigma_w_n, normal_physical_space);

                            // u_1 \dot n = H_slave[j] * n[jdim]
                            //----------------------------------
                            rLeftHandSideMatrix(iglob, jglob) += Guglielmo_innovation*sigma_w_n_dot_n * H_slave(0,j) * normal_physical_space[jdim] * IntToReferenceWeight;
                            

                            // PENALTY g_n = 0
                            //*********************************************** */
                            rLeftHandSideMatrix(iglob, jglob) += H_master(0,i)*normal_physical_space[idim]
                                                                *H_slave(0,j)*normal_physical_space[jdim]*penalty_integration;





                            // double sigma_n_w_dot_n = 0.0;
                            // double sigma_n_w_dot_tau = 0.0;
                            // for (IndexType kdim = 0; kdim < 2; kdim++) {
                            //     const double sigma_n_w_kdim = (DB_master(id1, 2*i+kdim)* normal_physical_space[0] + DB_master(id2, 2*i+kdim)* normal_physical_space[1]);

                            //     sigma_n_w_dot_n += sigma_n_w_kdim*normal_physical_space[kdim]; 
                            //     // sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[kdim];//eventually for friction
                            // }


                            // // rLeftHandSideMatrix(iglob, jglob) += H_master(0,i)*normal_physical_space[idim]
                            // //                                     *H_slave(0,j)*normal_physical_space[jdim]*penalty_integration;

                            // rLeftHandSideMatrix(iglob, jglob) += Guglielmo_innovation*H_slave(0,j)* normal_physical_space[jdim]*sigma_n_w_dot_n* normal_physical_space[idim]* IntToReferenceWeight;

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

                            // // // FLUX 
                            // // [\sigma_1(u) \dot n] \dot n * (w_2 \dot n)
                            Vector sigma_u_n(2);
                            sigma_u_n[0] = (DB_master(0, jglob)* normal_physical_space[0] + DB_master(2, jglob)* normal_physical_space[1]);
                            sigma_u_n[1] = (DB_master(2, jglob)* normal_physical_space[0] + DB_master(1, jglob)* normal_physical_space[1]);

                            double sigma_u_n_dot_n = inner_prod(sigma_u_n, normal_physical_space);

                            // w_2 \dot n = H_slave[i] * n[idim]
                            //----------------------------------
                            rLeftHandSideMatrix(iglob, jglob) += sigma_u_n_dot_n * H_slave(0, i) * normal_physical_space[idim] * IntToReferenceWeight;



                            
                            // double sigma_n_u_dot_n = 0.0;
                            // double sigma_n_u_dot_tau = 0.0;
                            // for (IndexType kdim = 0; kdim < 2; kdim++) {
                            //     const int kglob = 2*j+kdim;
                            //     const double sigma_n_u_kdim = -(DB_master(id1, kglob)* normal_physical_space[0] + DB_master(id2, kglob)* normal_physical_space[1]);

                            //     sigma_n_u_dot_n += sigma_n_u_kdim*normal_physical_space[kdim];
                            //     // sigma_n_u_dot_tau += sigma_n_u_kdim*tangent_physical_space[kdim]; // eventually for imposed forces on tangent direction
                            // }

                            // // i am imposing zero in the other direction

                            // rLeftHandSideMatrix(iglob, jglob) -= H_slave(0,i)*sigma_n_u_dot_n* normal_physical_space[jdim] * IntToReferenceWeight;
                            
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
                
                Vector gn = ZeroVector(2); 

                gn[0] = normal_gap*normal_physical_space[0];
                gn[1] = normal_gap*normal_physical_space[1];

                for (IndexType i = 0; i < number_of_nodes_master; i++) {

                    for (IndexType idim = 0; idim < 2; idim++) {

                        rRightHandSideVector[2*i+idim] -= H_master(0,i)*gn[idim]* penalty_integration;
                        const int id1 = idim*2;
                        const int iglob = 2*i+idim;

                        // // PENALTY FREE g_n = 0
                        // // rhs -> [\sigma_1(w) \dot n] \dot n (-g_{n,0})
                        // //*********************************************** */
                        Vector sigma_w_n(2);
                        sigma_w_n[0] = (DB_master(0, iglob)* normal_physical_space[0] + DB_master(2, iglob)* normal_physical_space[1]);
                        sigma_w_n[1] = (DB_master(2, iglob)* normal_physical_space[0] + DB_master(1, iglob)* normal_physical_space[1]);

                        double sigma_w_n_dot_n = inner_prod(sigma_w_n, normal_physical_space);
                        // double sigma_w_n_dot_n = sigma_w_n[idim] * normal_physical_space[idim];

                        rRightHandSideVector(iglob) -= Guglielmo_innovation*sigma_w_n_dot_n * IntToReferenceWeight *normal_gap;
                                                        //  * gn[jdim];


                                            


                        // for (IndexType jdim = 0; jdim < 2; jdim++) {
                        //     const int id2 = (id1+2)%3;

                            // double sigma_n_w_dot_n = 0.0;
                            // double sigma_n_w_dot_tau = 0.0;

                            // // const double sigma_n_w_kdim = (DB_master(id1, 2*i+jdim)* normal_physical_space[0] + DB_master(id2, 2*i+jdim)* normal_physical_space[1]);
                            // // sigma_n_w_dot_n = sigma_n_w_kdim*normal_physical_space[jdim]; 

                            // Vector sigma_n_w(2);
                            // sigma_n_w[0] = (DB_master(0, 2*i+jdim)* normal_physical_space[0] + DB_master(2, 2*i+jdim)* normal_physical_space[1]);
                            // sigma_n_w[1] = (DB_master(2, 2*i+jdim)* normal_physical_space[0] + DB_master(1, 2*i+jdim)* normal_physical_space[1]);

                            // const double sigma_n_dot_n_w = inner_prod(sigma_n_w, normal_physical_space); 

                            // const double sigma_n_dot_n_n_w_jdim = sigma_n_dot_n_w * normal_physical_space[jdim];
                            
                            // // sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[jdim];//eventually for friction


                            // rRightHandSideVector(2*i+idim) -= Guglielmo_innovation*gn[jdim]*sigma_n_dot_n_n_w_jdim * normal_physical_space[idim] * IntToReferenceWeight;
                        // }

                    }
                }
            }

            
            Vector temp = ZeroVector(mat_size_2);

            GetValuesVector(temp, 2);

            // RHS = ExtForces - K*temp;
            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
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
                const array_1d<double, 3>& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
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
                const array_1d<double, 3>& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
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
                const array_1d<double, 3>& displacement = GetMasterGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }

            for (IndexType i = 0; i < number_of_control_points_slave; ++i)
            {
                const array_1d<double, 3>& displacement = GetSlaveGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                IndexType index = i * 2 + mat_size_master;

                rValues[index] = displacement[0];
                rValues[index + 1] = displacement[1];
            }
        }
    }

    void SupportContact2DCondition::GetStrainVector(
        Vector& strainVector, IndexType index) const
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(index);

        const SizeType dim = 2;//r_geometry.WorkingSpaceDimension();

        const SizeType number_of_nodes = r_geometry.size();


        // Shape function derivatives 
        // Initialize Jacobian
        GeometryType::JacobiansType J0;

        // Initialize DN_DX
        Matrix DN_DX(number_of_nodes,2);
        Matrix InvJ0(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        double DetJ0;
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);


        // retrieve shape function values
        Vector displacements(2*number_of_nodes);
        GetValuesVector(displacements, index);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        // MODIFIED
        Matrix B = ZeroMatrix(3,number_of_nodes);

        CalculateB(B, DN_DX, number_of_nodes);

        if (strainVector.size() != dim) strainVector.resize(dim);

        strainVector = prod(B, displacements);
    }

    /**
     * @brief 
     * 
     * @param StrainVector 
     * @param index 
     * @param rCurrentProcessInfo 
     */
    void SupportContact2DCondition::SetConstitutiveVariables(
        Vector& StrainVector, IndexType index, const Kratos::ProcessInfo& rCurrentProcessInfo) 
    {
        const auto& r_geometry = GetGeometry().GetGeometryPart(index);
        const SizeType dim = 2;//r_geometry.WorkingSpaceDimension();

        ConstitutiveLaw::Pointer rpConstitutiveLaw = GetConstitutiveLaw(index);

        PropertiesType r_prop = GetProperty(index);

        ConstitutiveLaw::Parameters Values(r_geometry, r_prop, rCurrentProcessInfo);

        const SizeType strain_size = rpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);
    
        Values.SetStrainVector(StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);

        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        rpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        if (index == 0) {
            this->SetValue(CONSTITUTIVE_MATRIX_MASTER, Values.GetConstitutiveMatrix());
            this->SetValue(STRESS_MASTER, Values.GetStressVector());
        } 
        else if (index == 1)
        {
            this->SetValue(CONSTITUTIVE_MATRIX_SLAVE, Values.GetConstitutiveMatrix());
            this->SetValue(STRESS_SLAVE, Values.GetStressVector());
        }
        
    }

    void SupportContact2DCondition::SetNormalGap() 
    {
        Vector normal_physical_space_master = this->GetValue(NORMAL_MASTER);

        Vector GP_master_coord_deformed = GetMasterGeometry().Center() + GetValue(DISPLACEMENT_MASTER);
        Vector GP_slave_coord_deformed  = GetSlaveGeometry().Center() + GetValue(DISPLACEMENT_SLAVE);

        double normal_gap = inner_prod(GP_slave_coord_deformed - GP_master_coord_deformed, normal_physical_space_master);

        SetValue(NORMAL_GAP, normal_gap);
        // double tangential_gap = norm_2((deformed_pos_slave - deformed_pos_master) - normal_gap * normal_physical_space);
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
        ConstitutiveLaw::Parameters constitutive_law_parameters_master(
            GetMasterGeometry(), *mpPropMaster, rCurrentProcessInfo);

        mpConstitutiveLawMaster->FinalizeMaterialResponse(constitutive_law_parameters_master, ConstitutiveLaw::StressMeasure_Cauchy);

        ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
            GetSlaveGeometry(), *mpPropSlave, rCurrentProcessInfo);

        mpConstitutiveLawSlave->FinalizeMaterialResponse(constitutive_law_parameters_slave, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        Vector normal_physical_space = GetValue(NORMAL);
        Vector stress_vector_master = GetValue(STRESS_MASTER);

        Vector sigma_n(2);

        sigma_n[0] = stress_vector_master[0]*normal_physical_space[0] + stress_vector_master[2]*normal_physical_space[1];
        sigma_n[1] = stress_vector_master[2]*normal_physical_space[0] + stress_vector_master[1]*normal_physical_space[1];

        SetValue(NORMAL_STRESS, sigma_n);
        // //---------------------
    }

    void SupportContact2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

        // InitializeMaterial();
        ConstitutiveLaw::Parameters constitutive_law_parameters_master(
            GetMasterGeometry(), (*mpPropMaster), rCurrentProcessInfo);

        mpConstitutiveLawMaster->InitializeMaterialResponse(constitutive_law_parameters_master, ConstitutiveLaw::StressMeasure_Cauchy);

        // InitializeMaterial(); //slave
        ConstitutiveLaw::Parameters constitutive_law_parameters_slave(
            GetSlaveGeometry(), (*mpPropSlave), rCurrentProcessInfo);

        mpConstitutiveLawSlave->InitializeMaterialResponse(constitutive_law_parameters_slave, ConstitutiveLaw::StressMeasure_Cauchy);


        this->InitializeNonLinearIteration(rCurrentProcessInfo);

    }

    /**
     * @brief 
     * 
     * @param rCurrentProcessInfo 
     */
    void SupportContact2DCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo){

        // Master
        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        // Slave
        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();


        // JACOBIANS
        // master
        GeometryType::JacobiansType J0_master; double DetJ0_master;
        const SizeType dim = 2; //r_geometry_master.WorkingSpaceDimension();

        Matrix DN_DX_master(number_of_nodes_master,2); Matrix InvJ0_master(3,3);
        r_geometry_master.Jacobian(J0_master,this->GetIntegrationMethod());

        Matrix Jacobian_master = ZeroMatrix(3,3);
        Jacobian_master(0,0) = J0_master[0](0,0); Jacobian_master(0,1) = J0_master[0](0,1);
        Jacobian_master(1,0) = J0_master[0](1,0); Jacobian_master(1,1) = J0_master[0](1,1);
        Jacobian_master(2,2) = 1.0;

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian_master,InvJ0_master,DetJ0_master);

        Matrix sub_inv_jacobian_master = ZeroMatrix(2,2);
        sub_inv_jacobian_master(0,0) = InvJ0_master(0,0);
        sub_inv_jacobian_master(1,0) = InvJ0_master(1,0);
        sub_inv_jacobian_master(0,1) = InvJ0_master(0,1);
        sub_inv_jacobian_master(1,1) = InvJ0_master(1,1);

        // slave
        GeometryType::JacobiansType J0_slave; double DetJ0_slave;

        Matrix DN_DX_slave(number_of_nodes_slave,2); Matrix InvJ0_slave(3,3);
        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());

        Matrix Jacobian_slave = ZeroMatrix(3,3);
        Jacobian_slave(0,0) = J0_slave[0](0,0); Jacobian_slave(0,1) = J0_slave[0](0,1);
        Jacobian_slave(1,0) = J0_slave[0](1,0); Jacobian_slave(1,1) = J0_slave[0](1,1);
        Jacobian_slave(2,2) = 1.0;

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian_slave,InvJ0_slave,DetJ0_slave);

        Matrix sub_inv_jacobian_slave = ZeroMatrix(2,2);
        sub_inv_jacobian_slave(0,0) = InvJ0_slave(0,0);
        sub_inv_jacobian_slave(1,0) = InvJ0_slave(1,0);
        sub_inv_jacobian_slave(0,1) = InvJ0_slave(0,1);
        sub_inv_jacobian_slave(1,1) = InvJ0_slave(1,1);
        //-----------------------------------------------------------------------------

        // NORMALS
        // master
        array_1d<double, 3> tangent_parameter_space; array_1d<double, 3> tangent_physical_space;
        array_1d<double, 3> normal_physical_space; array_1d<double, 3> normal_parameter_space;

        r_geometry_master.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT
        normal_parameter_space[2] = 0.0;

        normal_physical_space = prod(trans(InvJ0_master),normal_parameter_space);
        normal_physical_space[2] = 0.0;
        normal_physical_space /= norm_2(normal_physical_space);

        tangent_physical_space = prod(trans(InvJ0_master),tangent_parameter_space);
        tangent_physical_space[2] = 0.0;
        tangent_physical_space /= norm_2(tangent_physical_space);

        // COMPUTE TRUE JACOBIAN DETERMINANT AND PHYSICAL NORMAL
        Vector add_factor = prod(Jacobian_master, tangent_parameter_space);

        DetJ0_master = norm_2(add_factor);

        // ##
        this->SetValue(NORMAL_MASTER, normal_physical_space);
        // ##
        // slave
        r_geometry_slave.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT

        normal_physical_space = prod(trans(InvJ0_slave),normal_parameter_space);
        normal_physical_space[2] = 0.0;
        normal_physical_space /= norm_2(normal_physical_space);

        tangent_physical_space = prod(trans(InvJ0_slave),tangent_parameter_space);
        tangent_physical_space[2] = 0.0;
        tangent_physical_space /= norm_2(tangent_physical_space);
        // ##
        this->SetValue(NORMAL_SLAVE, normal_physical_space);
        // ##

        // SET REFERENCE WEIGHT

        const double thickness = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;

        const double IntToReferenceWeight = r_geometry_master.IntegrationPoints()[0].Weight() * std::abs(DetJ0_master) * thickness;

        // ##
        SetValue(INTEGRATION_WEIGHT, IntToReferenceWeight);
        // ##

        // SET DISPLACEMENTS
        // master
        Vector coefficient_master(dim*number_of_nodes_master);
        GetValuesVector(coefficient_master, QuadraturePointCouplingGeometry2D<Point>::Master);

        const Matrix& N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());
        Matrix H_master = ZeroMatrix(dim, dim*number_of_nodes_master);
        for (IndexType i = 0; i < number_of_nodes_master; ++i)
        {
            for (IndexType i_dim = 0; i_dim < dim; i_dim++) {
                H_master(i_dim, dim*i+i_dim) = N_master(0, i);
            }
        }

        Vector displacement_master_sub = prod(H_master,coefficient_master);

        Vector displacement_master = ZeroVector(3); displacement_master[0] = displacement_master_sub[0]; displacement_master[1] = displacement_master_sub[1];

        // ##
        this->SetValue(DISPLACEMENT_MASTER, displacement_master);
        // ##
        // slave
        Vector coefficient_slave(dim*number_of_nodes_slave);
        GetValuesVector(coefficient_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);

        const Matrix& N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());
        Matrix H_slave = ZeroMatrix(dim, dim*number_of_nodes_slave);
        for (IndexType i = 0; i < number_of_nodes_slave; ++i)
        {
            for (IndexType i_dim = 0; i_dim < dim; i_dim++) {
                H_slave(i_dim, dim*i+i_dim) = N_slave(0, i);
            }
        }

        Vector displacement_slave_sub = prod(H_slave,coefficient_slave);

        Vector displacement_slave = ZeroVector(3); displacement_slave[0] = displacement_slave_sub[0]; displacement_slave[1] = displacement_slave_sub[1];

        // ##
        this->SetValue(DISPLACEMENT_SLAVE, displacement_slave);
        // ##

        // SET STRAINS/STRESSES/CONSTITUTIVE MATRIX 
        // master
        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        const GeometryType::ShapeFunctionsGradientsType& DN_De_master = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        noalias(DN_DX_master) = prod(DN_De_master[0],sub_inv_jacobian_master);
        Matrix B_master = ZeroMatrix(3,number_of_nodes_master);
        CalculateB(B_master, DN_DX_master, number_of_nodes_master);

        Vector strain_vector_master = prod(B_master, coefficient_master);
        
        
        // ##
        this->SetValue(STRAIN_MASTER, strain_vector_master);
        this->SetConstitutiveVariables(strain_vector_master, 0, rCurrentProcessInfo);
        // ##


        const GeometryType::ShapeFunctionsGradientsType& DN_De_slave = r_geometry_slave.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        noalias(DN_DX_slave) = prod(DN_De_slave[0],sub_inv_jacobian_slave);
        Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);
        CalculateB(B_slave, DN_DX_slave, number_of_nodes_slave);
        Vector strain_vector_slave = prod(B_slave, coefficient_slave);
        

        // ##
        this->SetValue(STRAIN_SLAVE, strain_vector_slave);
        this->SetConstitutiveVariables(strain_vector_slave, 1, rCurrentProcessInfo);
        // ##

        SetNormalGap();
    
    }

    /**
     * @brief 
     * 
     * @param rCurrentProcessInfo 
     */
    void SupportContact2DCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo){

        this->InitializeNonLinearIteration(rCurrentProcessInfo);
    }


    


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

        SetValue(NORMAL_GAP, flag_normal_gap);

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

    void SupportContact2DCondition::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        if (rValues.size() != 1)
            rValues.resize(1);
        if (rVariable == INTEGRATION_WEIGHT)
            rValues[0] = this->GetValue(INTEGRATION_WEIGHT);
        else if (rVariable == NORMAL_GAP)
            rValues[0] = this->GetValue(NORMAL_GAP);
    }

    void SupportContact2DCondition::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const SizeType dimension = 2;//GetMasterGeometry().WorkingSpaceDimension();

        if (rValues.size() != dimension)
            rValues.resize(dimension);
        if (rVariable == NORMAL_STRESS)
        {
            Vector normal_physical_space = this->GetValue(NORMAL);
            Vector sigma_voigt_master = GetValue(STRESS_MASTER);

            Vector normal_stress_master(2);

            normal_stress_master[0] = sigma_voigt_master[0]*normal_physical_space[0] + sigma_voigt_master[2]*normal_physical_space[1];
            normal_stress_master[1] = sigma_voigt_master[2]*normal_physical_space[0] + sigma_voigt_master[1]*normal_physical_space[1];
            
            this->SetValue(NORMAL_STRESS, normal_stress_master);

            rValues[0] = normal_stress_master;
        }
        else if (rVariable == NORMAL_MASTER) 
            rValues[0] = this->GetValue(NORMAL_MASTER);
    }

} // Namespace Kratos