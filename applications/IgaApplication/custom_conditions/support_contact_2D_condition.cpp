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

        const int integrationDomain = GetGeometry().GetValue(ACTIVATION_LEVEL);

        const double penalty = (*mpPropMaster)[PENALTY_FACTOR];

        const auto& r_geometry_master = GetMasterGeometry();
        const SizeType number_of_nodes_master = r_geometry_master.size();

        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType number_of_nodes_slave = r_geometry_slave.size();

        const SizeType mat_size_1 = (number_of_nodes_master) * 2;
        const SizeType mat_size_2 = (number_of_nodes_master+number_of_nodes_slave) * 2;
        //resizing as needed the LHS
        if(rLeftHandSideMatrix.size1() != mat_size_1 && rLeftHandSideMatrix.size2() != mat_size_2 )
            rLeftHandSideMatrix.resize(mat_size_1,mat_size_2,false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size_1,mat_size_2); //resetting LHS
        
        // resizing as needed the RHS
        if(rRightHandSideVector.size() != mat_size_1)
            rRightHandSideVector.resize(mat_size_1,false);
        noalias(rRightHandSideVector) = ZeroVector(mat_size_1); //resetting RHS

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints();

        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        r_geometry_master.DeterminantOfJacobian(determinant_jacobian_vector);

        // Shape function derivatives (NEW) 
        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        GeometryType::JacobiansType J0_slave;
        // Initialize DN_DX
        const unsigned int dim = 2;
        Matrix DN_DX(number_of_nodes_master,2);
        Matrix InvJ0(dim,dim);

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

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry_master.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry_master.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        Vector GP_parameter_coord(2); 
        GP_parameter_coord = prod(r_geometry_master.Center(),J0[0]);
        
        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);
        tangent_physical_space = prod(trans(J0[0]),tangent_parameter_space);

        // MODIFIED
        Vector old_displacement_master(number_of_nodes_master);
        GetValuesVector(old_displacement_master, QuadraturePointCouplingGeometry2D<Point>::Master);

        Vector old_displacement_slave(number_of_nodes_slave);
        GetValuesVector(old_displacement_slave, QuadraturePointCouplingGeometry2D<Point>::Slave);
        
        /////
         // ---------------------MODIFIED 
        std::ofstream outputFile("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << GP_parameter_coord[0] << " " << GP_parameter_coord[1]  <<"\n";
        outputFile.close();
        

        const Matrix& N_master = r_geometry_master.ShapeFunctionsValues(this->GetIntegrationMethod());
        const Matrix& N_slave = r_geometry_slave.ShapeFunctionsValues(this->GetIntegrationMethod());

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        const double thickness = (*mpPropMaster).Has(THICKNESS) ? (*mpPropMaster)[THICKNESS] : 1.0;

        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

        // MODIFIED
        Matrix B_master = ZeroMatrix(3,number_of_nodes_master);

        Matrix B_slave = ZeroMatrix(3,number_of_nodes_slave);

        CalculateB(B_master, DN_DX, number_of_nodes_master);


        CalculateB(B_slave, DN_DX, number_of_nodes_slave);


     //---------- MODIFIED ----------------------------------------------------------------

        // MASTER
        ConstitutiveLaw::Parameters Values_master(r_geometry_master, *mpPropMaster, rCurrentProcessInfo);

        const SizeType strain_size_master = mpConstitutiveLawMaster->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptionsMaster=Values_master.GetOptions();

        ConstitutiveLawOptionsMaster.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptionsMaster.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptionsMaster.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptionsMaster.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables_master(strain_size_master);

        Vector old_strain_master = prod(B_master,old_displacement_master);
    
        // Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values_master.SetStrainVector(old_strain_master);

        Values_master.SetStressVector(this_constitutive_variables_master.StressVector);
        Values_master.SetConstitutiveMatrix(this_constitutive_variables_master.D);
        mpConstitutiveLawMaster->CalculateMaterialResponse(Values_master, ConstitutiveLaw::StressMeasure_PK2);

        const Matrix& r_D_master = Values_master.GetConstitutiveMatrix();

        // ------ SLAVE ---------------------------

        ConstitutiveLaw::Parameters Values_slave(r_geometry_slave, *mpPropSlave, rCurrentProcessInfo);

        const SizeType strain_size_slave = mpConstitutiveLawSlave->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptionsSlave=Values_slave.GetOptions();

        ConstitutiveLawOptionsSlave.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptionsSlave.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptionsSlave.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptionsSlave.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables_slave(strain_size_slave);

        Vector old_strain_slave = prod(B_slave,old_displacement_slave);
    
        // Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values_slave.SetStrainVector(old_strain_slave);
        Values_slave.SetStressVector(this_constitutive_variables_slave.StressVector);
        Values_slave.SetConstitutiveMatrix(this_constitutive_variables_slave.D);
        mpConstitutiveLawSlave->CalculateMaterialResponse(Values_slave, ConstitutiveLaw::StressMeasure_PK2);

        const Matrix& r_D_slave = Values_slave.GetConstitutiveMatrix();

        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------
        //  CHECK CONTACT CONDITION FOR THE PAIR
        //  ----------------------------------------------------------------------------
        //  ----------------------------------------------------------------------------

        // compute the normal gap
        Vector GP_parameter_coord_on_slave(2); 
        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        GP_parameter_coord_on_slave = prod(r_geometry_slave.Center(),J0_slave[0]);

        const double normal_gap = inner_prod(GP_parameter_coord_on_slave - GP_parameter_coord, normal_physical_space);

        const Vector& master_stress_vector = Values_master.GetStressVector();

        const double normal_master_stress = (master_stress_vector[0]*normal_physical_space[0]+master_stress_vector[1]*normal_physical_space[1])*normal_physical_space[0]+
                                            (master_stress_vector[1]*normal_physical_space[0]+master_stress_vector[1]*normal_physical_space[2])*normal_physical_space[1];

        bool contact_is_active = true;

        const double toll = 1e-12;
        if (normal_master_stress - (*mpPropMaster)[YOUNG_MODULUS]*normal_gap > toll) {
            contact_is_active = true;
        }
                
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

        if (contact_is_active) {

            if (integrationDomain == 0) // integration on master
            {
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

                                
                                double sigma_n_u_dot_n = 0.0;
                                double sigma_n_u_dot_tau = 0.0;

                                double sigma_n_w_dot_n = 0.0;
                                double sigma_n_w_dot_tau = 0.0;
                                for (IndexType kdim = 0; kdim < 2; kdim++) {
                                    const int kglob = 2*j+kdim;
                                    const double sigma_n_u_kdim = (DB_master(id1, kglob)* normal_physical_space[0] + DB_master(id2, kglob)* normal_physical_space[1]);

                                    sigma_n_u_dot_n += sigma_n_u_kdim*normal_physical_space[kdim];
                                    sigma_n_u_dot_tau += sigma_n_u_kdim*tangent_physical_space[kdim]; // eventually for imposed forces on tangent direction

                                    const double sigma_n_w_kdim = (DB_master(id1, 2*i+kdim)* normal_physical_space[0] + DB_master(id2, 2*i+kdim)* normal_physical_space[1]);

                                    sigma_n_w_dot_n += sigma_n_w_kdim*normal_physical_space[kdim]; 
                                    sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[kdim];//eventually for friction
                                }


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
                                for (IndexType kdim = 0; kdim < 2; kdim++) {
                                    const double sigma_n_w_kdim = (DB_master(id1, 2*i+kdim)* normal_physical_space[0] + DB_master(id2, 2*i+kdim)* normal_physical_space[1]);

                                    sigma_n_w_dot_n += sigma_n_w_kdim*normal_physical_space[kdim]; 
                                    sigma_n_w_dot_tau += sigma_n_w_kdim*tangent_physical_space[kdim];//eventually for friction
                                }


                                rRightHandSideVector(2*i+idim) -= Guglielmo_innovation*gn[jdim]*sigma_n_w_dot_n * normal_physical_space[idim] * IntToReferenceWeight;
                            }

                        }
                    }
                    
                    Vector temp = ZeroVector(mat_size_2);

                    GetValuesVector(temp, 2);

                    // RHS = ExtForces - K*temp;
                    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
                } 
            }
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // INTEGRATION ON SLAVE
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            else if (integrationDomain == 1) // integration on slave
            {
                // ONLY SLAVE (substitute /sigma_1 *n  =sigma_2*n)
                Matrix DB_master = prod(r_D_master,B_master);
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

                                
                                double sigma_n_u_dot_n = 0.0;
                                double sigma_n_u_dot_tau = 0.0;
                                for (IndexType kdim = 0; kdim < 2; kdim++) {
                                    const int kglob = 2*j+kdim;
                                    const double sigma_n_u_kdim = +(DB_slave(id1, kglob)* normal_physical_space[0] + DB_slave(id2, kglob)* normal_physical_space[1]);

                                    sigma_n_u_dot_n += sigma_n_u_kdim*normal_physical_space[kdim];
                                    sigma_n_u_dot_tau += sigma_n_u_kdim*tangent_physical_space[kdim]; // eventually for imposed forces on tangent direction
                                }

                                // i am imposing zero in the other direction

                                rLeftHandSideMatrix(iglob, jglob) -= H_master(0,i)*sigma_n_u_dot_n* normal_physical_space[jdim] * IntToReferenceWeight;

                                
                            }

                        }
                    }
                }

                if (CalculateResidualVectorFlag) {
                    
                    Vector fn = ZeroVector(2); 

                    fn[0] = 0.0;
                    fn[1] = 0.0;

                    for (IndexType i = 0; i < number_of_nodes_master; i++) {

                        for (IndexType idim = 0; idim < 2; idim++) {

                            rRightHandSideVector(2*i+idim) -= H_master(0,i)*fn[idim]* IntToReferenceWeight;

                        }
                    }
                    
                    Vector temp = ZeroVector(mat_size_2);

                    GetValuesVector(temp, 2);

                    // RHS = ExtForces - K*temp;
                    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
                } 
            }
        } else { // CONTACT NOT ACTIVE FOR THE PAIR

            KRATOS_WATCH("NOT ACTIVE CONTACT")
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
        
        const SizeType mat_size = number_of_control_points * 3;

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

        mpConstitutiveLawMaster->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

        
    /////////////////////////
        const auto& r_geometry = GetMasterGeometry();
        const SizeType nb_nodes = r_geometry.size();

        // Integration Points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        GeometryType::JacobiansType J0;
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        // Get the parameter coordinates
        Vector GP_parameter_coord(2); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]); // Only one Integration Points 

        double x_coord_gauss_point = 0;
        double y_coord_gauss_point = 0;
        double rOutput = 0;

        for (IndexType i = 0; i < nb_nodes; ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_X);
            rOutput += r_N(0, i) * output_solution_step_value;
            x_coord_gauss_point += r_N(0, i) * r_geometry[i].X0();
            y_coord_gauss_point += r_N(0, i) * r_geometry[i].Y0();
        }   

    //////////////////////////
    /////////////// SLAVE 
    ////////////////////////
        const auto& r_geometry_slave = GetSlaveGeometry();
        const SizeType nb_nodes_slave = r_geometry_slave.size();

        // Integration Points
        const GeometryType::IntegrationPointsArrayType& integration_points_slave = r_geometry_slave.IntegrationPoints();
        // Shape function values
        const Matrix& r_N_slave = r_geometry_slave.ShapeFunctionsValues();

        GeometryType::JacobiansType J0_slave;
        r_geometry_slave.Jacobian(J0_slave,this->GetIntegrationMethod());
        // Get the parameter coordinates
        Vector GP_parameter_coord_slave(2); 
        GP_parameter_coord_slave = prod(r_geometry_slave.Center(),J0_slave[0]); // Only one Integration Points 

        double x_coord_gauss_point_slave = 0;
        double y_coord_gauss_point_slave = 0;
        double rOutput_slave = 0;

        for (IndexType i = 0; i < nb_nodes_slave; ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value = r_geometry_slave[i].GetSolutionStepValue(DISPLACEMENT_X);
            rOutput_slave += r_N_slave(0, i) * output_solution_step_value;
            x_coord_gauss_point_slave += r_N_slave(0, i) * r_geometry_slave[i].X0();
            y_coord_gauss_point_slave += r_N_slave(0, i) * r_geometry_slave[i].Y0();
        }     


        std::ofstream gap_file("txt_files/gap_output.txt", std::ios::app);
        if (gap_file.is_open()) {
            gap_file << std::scientific << std::setprecision(14); // Set precision to 10^-14

            const double x_final_pos_master = x_coord_gauss_point + rOutput;
            const double x_final_pos_slave = x_coord_gauss_point_slave + rOutput_slave;

            gap_file << x_final_pos_master-x_final_pos_slave << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points[0].Weight() << std::endl;
            gap_file.close();
        }  


        std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << rOutput << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points[0].Weight() << std::endl;
            output_file.close();
        } 


        std::ofstream outputFile("txt_files/Gauss_Point_coordinates.txt", std::ios::app);
        if (!outputFile.is_open())
        {
            std::cerr << "Failed to open the file for writing." << std::endl;
            return;
        }
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << x_coord_gauss_point << "  " << y_coord_gauss_point <<"\n";
        outputFile.close();
    }

    void SupportContact2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetMasterGeometry(), (*mpPropMaster), rCurrentProcessInfo);

        mpConstitutiveLawMaster->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_PK2);

    }

} // Namespace Kratos
