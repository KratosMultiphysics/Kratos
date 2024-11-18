
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
#include "custom_conditions/support_solid_2D_condition.h"

namespace Kratos
{

    void SupportSolid2DCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        InitializeMaterial();
    }


    void SupportSolid2DCondition::InitializeMaterial()
    {
        KRATOS_TRY
        if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry = GetGeometry();
            const Properties& r_properties = GetProperties();
            const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

        } else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        KRATOS_CATCH( "" );

    }
    void SupportSolid2DCondition::CalculateAll(
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

        const SizeType mat_size = number_of_nodes * 2;
        //resizing as needed the LHS
        if(rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size,mat_size,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
        
        // resizing as needed the RHS
        if(rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size,false);
        noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

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
        Matrix DN_DX(number_of_nodes,2);
        Matrix InvJ0(3,3);

        // Compute the normals
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 3> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT
        normal_parameter_space[2] = 0.0;

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;
        
        // MODIFIED
        Vector old_displacement(mat_size);
        GetValuesVector(old_displacement);
        
        const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        Matrix Jacobian = ZeroMatrix(3,3);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);
        Jacobian(2,2) = 1.0;

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        Vector add_factor = prod(Jacobian, tangent_parameter_space);
        add_factor[2] = 0.0;

        DetJ0 = norm_2(add_factor);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix sub_inv_jacobian = ZeroMatrix(2,2);
        sub_inv_jacobian(0,0) = InvJ0(0,0);
        sub_inv_jacobian(1,0) = InvJ0(1,0);
        sub_inv_jacobian(0,1) = InvJ0(0,1);
        sub_inv_jacobian(1,1) = InvJ0(1,1);
        noalias(DN_DX) = prod(DN_De[0],sub_inv_jacobian);

        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

        SetValue(INTEGRATION_WEIGHT, IntToReferenceWeight);

        // MODIFIED
        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX);

        normal_physical_space = prod(trans(InvJ0),normal_parameter_space);
        normal_physical_space[2] = 0.0;

        normal_physical_space /= norm_2(normal_physical_space);

        SetValue(NORMAL, normal_physical_space);


    //---------- MODIFIED ----------------------------------------------------------------


        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B,old_displacement);
    
        // Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStrainVector(old_strain);

        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        const Matrix& r_D = Values.GetConstitutiveMatrix();
        //---------------------


            
        Matrix H = ZeroMatrix(1, number_of_nodes);
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)            = N(0, i);
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

        Matrix DB = prod(r_D,B);
        double integration_factor = IntToReferenceWeight;

        if (this->Has(DIRECTION)){
            // ASSIGN BC BY DIRECTION
            //--------------------------------------------------------------------------------------------
            Vector direction = this->GetValue(DIRECTION);

            for (IndexType i = 0; i < number_of_nodes; i++) {
                for (IndexType j = 0; j < number_of_nodes; j++) {
                    
                    for (IndexType idim = 0; idim < 2; idim++) {
                        const int iglob = 2*i+idim;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int jglob = 2*j+jdim;

                            // PENALTY TERM
                            rLeftHandSideMatrix(iglob, jglob) -= H(0,i)*H(0,j)* penalty_integration * direction[idim] * direction[jdim];

                             // FLUX 
                            // [sigma(u) \dot n] \dot n * (-w \dot n)
                            //*********************************************** */
                            Vector sigma_u_n = ZeroVector(3);
                            sigma_u_n[0] = DB(0, jglob)*normal_physical_space[0] + DB(2, jglob)*normal_physical_space[1];
                            sigma_u_n[1] = DB(2, jglob)*normal_physical_space[0] + DB(1, jglob)*normal_physical_space[1];

                            double sigma_u_n_dot_direction = inner_prod(sigma_u_n, direction);

                            rLeftHandSideMatrix(iglob, jglob) -= H(0,i) * sigma_u_n_dot_direction * direction[idim] * integration_factor;

                            // // PENALTY FREE g_n = 0
                            // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                            // //*********************************************** */
                            Vector sigma_w_n = ZeroVector(3);
                            sigma_w_n[0] = (DB(0, iglob)* normal_physical_space[0] + DB(2, iglob)* normal_physical_space[1]);
                            sigma_w_n[1] = (DB(2, iglob)* normal_physical_space[0] + DB(1, iglob)* normal_physical_space[1]);

                            double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                            rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*H(0,j) * sigma_w_n_dot_direction * direction[jdim] * integration_factor;
                        }

                    }
                }
            }

            if (CalculateResidualVectorFlag) {
                
                double displacement_module = this->GetValue(MODULE);
                
                for (IndexType i = 0; i < number_of_nodes; i++) {

                    for (IndexType idim = 0; idim < 2; idim++) {
                        const int iglob = 2*i+idim;

                        rRightHandSideVector(iglob) -= H(0,i) * direction[idim] * displacement_module * penalty_integration;

                        // // PENALTY FREE g_n = 0
                        // // rhs -> [\sigma_1(w) \dot n] \dot n (-g_{n,0})
                        // //*********************************************** */
                        Vector sigma_w_n = ZeroVector(3);
                        sigma_w_n[0] = (DB(0, iglob)* normal_physical_space[0] + DB(2, iglob)* normal_physical_space[1]);
                        sigma_w_n[1] = (DB(2, iglob)* normal_physical_space[0] + DB(1, iglob)* normal_physical_space[1]);

                        double sigma_w_n_dot_n = inner_prod(sigma_w_n, direction);

                        rRightHandSideVector(iglob) -= Guglielmo_innovation*sigma_w_n_dot_n * integration_factor *displacement_module;

                    }
                }
            }
        }
        else {
            // ASSIGN BC BY COMPONENTS 
            //--------------------------------------------------------------------------------------------
            for (IndexType i = 0; i < number_of_nodes; i++) {
                for (IndexType j = 0; j < number_of_nodes; j++) {
                    
                    for (IndexType idim = 0; idim < 2; idim++) {
                        const int id1 = 2*idim;
                        const int iglob = 2*i+idim;

                        rLeftHandSideMatrix(2*i+idim, 2*j+idim) -= H(0,i)*H(0,j)* penalty_integration;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int id2 = (id1+2)%3;
                            const int jglob = 2*j+jdim;
                            rLeftHandSideMatrix(iglob, jglob) -= H(0,i)*(DB(id1, jglob)* normal_physical_space[0] + DB(id2, jglob)* normal_physical_space[1]) * integration_factor;

                            rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*H(0,j)*(DB(id1, 2*i+jdim)* normal_physical_space[0] + DB(id2, 2*i+jdim)* normal_physical_space[1]) * integration_factor;
                        }

                    }
                }
            }
            
            if (CalculateResidualVectorFlag) {
                
                Vector u_D = ZeroVector(2); //->GetValue(DISPLACEMENT);
                
                // double x = GP_parameter_coord[0];
                // double y = GP_parameter_coord[1];

                u_D[0] = this->GetValue(DISPLACEMENT_X);
                u_D[1] = this->GetValue(DISPLACEMENT_Y);

                for (IndexType i = 0; i < number_of_nodes; i++) {

                    for (IndexType idim = 0; idim < 2; idim++) {

                        rRightHandSideVector[2*i+idim] -= H(0,i)*u_D[idim]* penalty_integration;
                        const int id1 = idim*2;

                        for (IndexType jdim = 0; jdim < 2; jdim++) {
                            const int id2 = (id1+2)%3;
                            rRightHandSideVector(2*i+idim) -= Guglielmo_innovation*u_D[jdim]*(DB(id1, 2*i+jdim)* normal_physical_space[0] + DB(id2, 2*i+jdim)* normal_physical_space[1]) * integration_factor;
                        }

                    }


                }
            }
        }
            
        Vector temp = ZeroVector(number_of_nodes);

        GetValuesVector(temp);

        // RHS = ExtForces - K*temp;
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
        
        // exit(0);
        KRATOS_CATCH("")
    }

    int SupportSolid2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SupportSolid2DCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 2 * number_of_nodes)
            rResult.resize(2 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void SupportSolid2DCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }
    };


    void SupportSolid2DCondition::GetValuesVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    void SupportSolid2DCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
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


    void SupportSolid2DCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        const auto& r_geometry = GetGeometry();
        const SizeType nb_nodes = r_geometry.size();
        const SizeType mat_size = nb_nodes * 2;

        // Shape function derivatives (NEW) 
        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        const unsigned int dim = 2;
        Matrix DN_DX(nb_nodes,2);
        Matrix InvJ0(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;
        // MODIFIED
        Vector old_displacement(mat_size);
        GetValuesVector(old_displacement);
        
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

        CalculateB(B, DN_DX);

        array_1d<double, 3> normal_physical_space = GetValue(NORMAL);

        // GET STRESS VECTOR
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B,old_displacement);
    
        // Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStrainVector(old_strain);

        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        const Vector sigma = Values.GetStressVector();
        Vector sigma_n(2);

        sigma_n[0] = sigma[0]*normal_physical_space[0] + sigma[2]*normal_physical_space[1];
        sigma_n[1] = sigma[2]*normal_physical_space[0] + sigma[1]*normal_physical_space[1];

        //-----------------------------------------
        // Vector sigma_n = ZeroVector(2);
        // const Matrix D = Values.GetConstitutiveMatrix();
        // Matrix DB = prod(D,B);
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

    void SupportSolid2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
    }

} // Namespace Kratos
