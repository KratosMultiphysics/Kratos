//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Polytimi Zis
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_solid_3D_condition.h"

namespace Kratos
{

    void SupportSolid3DCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
                
        // const auto& integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // //Constitutive Law initialisation
        // if ( mConstitutiveLawVector.size() != integration_points.size() )
        //     mConstitutiveLawVector.resize(integration_points.size());
            
        InitializeMaterial();
    }




void SupportSolid3DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    
    bool required = false;
        if (mpConstitutiveLaw->RequiresInitializeMaterialResponse()) {
            required = true;
        }
    if (required) {
        const auto& r_geom = GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        const Properties& r_properties = GetProperties();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geom,r_properties,rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points

        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);
        KRATOS_WATCH(N_values)
        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints();

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, 0, mThisIntegrationMethod);

        // Compute constitutive law variables
        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, integration_points);

        
        // Call the constitutive law to update material variables
        mpConstitutiveLaw->InitializeMaterialResponse(Values, GetStressMeasure());

        // TODO: Deprecated, remove this
        mpConstitutiveLaw->InitializeSolutionStep( r_properties, r_geom, row( N_values, 0 ), rCurrentProcessInfo);
        
    }
}    
    void SupportSolid3DCondition::InitializeMaterial()
    {

        KRATOS_TRY
        if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
            const GeometryType& r_geometry = GetGeometry();
            const Properties& r_properties = GetProperties();
            const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
            KRATOS_WATCH(N_values)

            mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

        } else
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        KRATOS_CATCH( "" );

    }

    /***********************************************************************************/
    /***********************************************************************************/

    ConstitutiveLaw::StressMeasure SupportSolid3DCondition::GetStressMeasure() const
    {
    return ConstitutiveLaw::StressMeasure_PK2;
    }


    void SupportSolid3DCondition::CalculateAll(
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

        const SizeType mat_size = number_of_nodes * 3;
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
        const unsigned int dim = 3;
        Matrix DN_DX(number_of_nodes,3);
        Matrix InvJ0(dim,dim);

        // Compute the normals
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 3> normal_physical_space;
        array_1d<double, 3> normal_parameter_space;

        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + 
                             tangent_parameter_space[1] * tangent_parameter_space[1] +
                             tangent_parameter_space[2] * tangent_parameter_space[2]);

        // // NEW FOR GENERAL JACOBIAN
        // normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        // normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT
        // normal_parameter_space[2] = 0.0; // Since this is the parameter space normal

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;
        
        r_geometry.Calculate(NORMAL, normal_parameter_space);
        //KRATOS_WATCH(normal_parameter_space)


        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);
        Vector GP_parameter_coord(3); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]);
        // KRATOS_WATCH(r_geometry.Center())
        // KRATOS_WATCH(J0[0])

        // MODIFIED
        Vector old_displacement(mat_size);
        GetValuesVector(old_displacement);

    
        const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        KRATOS_WATCH(N)
        
        Matrix Jacobian = ZeroMatrix(3,3);
        Jacobian(0, 0) = J0[0](0, 0);
        Jacobian(0, 1) = J0[0](0, 1);
        Jacobian(1, 0) = J0[0](1, 0);
        Jacobian(1, 1) = J0[0](1, 1);
        Jacobian(0,2) = J0[0](0,2);
        Jacobian(1,2) = J0[0](1,2);
        Jacobian(2,0) = J0[0](2,0);
        Jacobian(2,1) = J0[0](2,1);
        Jacobian(2, 2) = J0[0](2, 2);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);


        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

        // MODIFIED
        Matrix B = ZeroMatrix(6, mat_size);

        CalculateB(B, DN_DX);


    //---------- MODIFIED ----------------------------------------------------------------


        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B,old_displacement);
    
        // Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStrainVector(old_strain);

        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        const Matrix& r_D = Values.GetConstitutiveMatrix();
        //---------------------
       
        Matrix H = ZeroMatrix(1, number_of_nodes);
        // Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
        // Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
        // Vector H_vector = ZeroVector(number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)            = N(0, i);
            // H_vector(i)        = N(0, i); 
            // for (IndexType idim = 0; idim < dim; idim++) {
            //     DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];           
            //     DN_dot_n_vec(i)  += DN_DX(i, idim) * normal_physical_space[idim];
            // } 
                
        }
        // KRATOS_WATCH(DetJ0)
        // KRATOS_WATCH(penalty)
        // Differential area
        double penalty_integration = penalty * integration_points[0].Weight() * std::abs(DetJ0);

        // KRATOS_WATCH(penalty_integration)

        // Collins, Lozinsky & Scovazzi innovation
        double Guglielmo_innovation = 1.0;  // = 1 -> Penalty approach
                                                // = -1 -> Free-penalty approach
        if (penalty == -1.0) {
            penalty_integration = 0.0;
            Guglielmo_innovation = -1.0;
        }
            
        ConstitutiveLaw::StressVectorType &r_stress_vector = Values.GetStressVector();
        //Values.GetStressVector(rStressVector);
        // KRATOS_WATCH(r_stress_vector);


        // KRATOS_WATCH(r_stress_vector)
        Matrix DB = prod(r_D, B);
        // KRATOS_WATCH(DB)
        double integration_factor = IntToReferenceWeight;
        // KRATOS_WATCH(integration_factor)

        
        for (IndexType i = 0; i < number_of_nodes; i++) {
            for (IndexType j = 0; j < number_of_nodes; j++) {

                for (IndexType idim = 0; idim < 3; idim++) {  // Loop over 3 dimensions
                    //const int id1 = (3 * idim) % 6;
                    const int iglob = 3 * i + idim;
                    rLeftHandSideMatrix(iglob, 3 * j + idim) -= H(0, i) * H(0, j) * penalty_integration;
                    
                    const int id1 = (6-idim)%6;

                    for (IndexType jdim = 0; jdim < 3; jdim++) {  // Loop over 3 dimensions
                        //const int id1 = idim + idim * 3 + jdim;
                        //const int id2 = (id1 + 3) % 6;
                        const int jglob = 3 * j + jdim; //pass

                            if (idim == 0) {
                                rLeftHandSideMatrix(iglob, jglob) -= H(0, i) * (DB(0, jglob) * normal_parameter_space[0] +
                                    DB(3, jglob) * normal_parameter_space[1] +
                                    DB(4, jglob) * normal_parameter_space[2]) * integration_factor;

                    //             rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation * H(0, j) * (DB(0, 3 * i + jdim) * normal_parameter_space[0] +
                    //                 DB(3, 3 * i + jdim) * normal_parameter_space[1] +
                    //                 DB(4, 3 * i + jdim) * normal_parameter_space[2]) * integration_factor;
                            }
                            else if (idim == 1) {
                                rLeftHandSideMatrix(iglob, jglob) -= H(0, i) * (DB(3, jglob) * normal_parameter_space[0] +
                                    DB(1, jglob) * normal_parameter_space[1] +
                                    DB(5, jglob) * normal_parameter_space[2]) * integration_factor;

                    //             rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation * H(0, j) * (DB(3, 3 * i + jdim) * normal_parameter_space[0] +
                    //                 DB(1, 3 * i + jdim) * normal_parameter_space[1] +
                    //                 DB(5, 3 * i + jdim) * normal_parameter_space[2]) * integration_factor;
                            }
                            else if (idim == 2) {
                                rLeftHandSideMatrix(iglob, jglob) -= H(0, i) * (DB(4, jglob) * normal_parameter_space[0] +
                                    DB(5, jglob) * normal_parameter_space[1] +
                                    DB(2, jglob) * normal_parameter_space[2]) * integration_factor;

                    //             rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation * H(0, j) * (DB(4, 3 * i + jdim) * normal_parameter_space[0] +
                    //                 DB(5, 3 * i + jdim) * normal_parameter_space[1] +
                    //                 DB(2, 3 * i + jdim) * normal_parameter_space[2]) * integration_factor;
                            }           

                    //     //KRATOS_WATCH(rLeftHandSideMatrix(iglob, jglob))
                    }
                }
            }
        }


    // KRATOS_WATCH(rLeftHandSideMatrix)
        
        if (CalculateResidualVectorFlag) {
            
            // const double& temperature = Has(TEMPERATURE)
            //     ? this->GetValue(TEMPERATURE)
            //     : 0.0;
            
            Vector u_D = ZeroVector(3); //->GetValue(DISPLACEMENT);
            
            double x = GP_parameter_coord[0];
            double y = GP_parameter_coord[1];
            double z = GP_parameter_coord[2];


            // u_D[0] = cos(x)*sinh(y);
            // u_D[1] = sin(x)*cosh(y);

            // const auto& u_D[0] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
            // u_D[1] = this->GetValue(DISPLACEMENT_Y);
            // u_D[2] = this->GetValue(DISPLACEMENT_Z);

            //u_D = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);

            // u_D[0] = this->GetValue(DISPLACEMENT_X);
            // u_D[1] = this->GetValue(DISPLACEMENT_Y);
            // u_D[2] = this->GetValue(DISPLACEMENT_Z);

            
            u_D = this->GetValue(DISPLACEMENT); //will retrieve something already stored but not useful, may be trash
            // u_D = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
            
            // u_D = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1); // previous displacement
            // GetValuesVector(u_D);
            // KRATOS_WATCH(u_D)
            // KRATOS_WATCH(penalty)
            // KRATOS_WATCH(Guglielmo_innovation)
            // KRATOS_WATCH(integration_factor)

            const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(r_stress_vector);

            for (IndexType i = 0; i < number_of_nodes; i++) {

                for (IndexType idim = 0; idim < 3; idim++) {

                    rRightHandSideVector[3*i+idim] -= H(0, i) * u_D[idim] * penalty_integration;
                    
                    
                    // const int iglob = 3 * i + idim;


                    // const int id1 = (6-idim)%6;

                    // for (IndexType jdim = 0; jdim < 3; jdim++) {

                    //         //rRightHandSideVector[iglob] -= Guglielmo_innovation * u_D[jdim] * prod(stress_tensor, normal_parameter_space) * integration_factor;
                    //         if (idim == 0) {
                    //             rRightHandSideVector[iglob] -= Guglielmo_innovation * u_D[jdim] * (
                    //                 DB(0, 3*i+jdim)* normal_parameter_space[0] + 
                    //                 DB(3, 3*i+jdim) * normal_parameter_space[1] +
                    //                 DB(4, 3*i+jdim) * normal_parameter_space[2]) * integration_factor;

                    //         }
                    //         else if (idim == 1) {
                    //             rRightHandSideVector[iglob] -= Guglielmo_innovation * u_D[jdim] * (
                    //                 DB(3, 3*i+jdim)* normal_parameter_space[0] + 
                    //                 DB(1, 3*i+jdim) * normal_parameter_space[1] +
                    //                 DB(5, 3*i+jdim) * normal_parameter_space[2]) * integration_factor;
                    //         }
                    //         else if (idim == 2) {
                    //             rRightHandSideVector[iglob] -= Guglielmo_innovation * u_D[jdim] * (
                    //                 DB(4, 3*i+jdim)* normal_parameter_space[0] + 
                    //                 DB(5, 3*i+jdim) * normal_parameter_space[1] +
                    //                 DB(2, 3*i+jdim) * normal_parameter_space[2]) * integration_factor;
                    //         } 
                    // }
                }   
            }
            // KRATOS_WATCH(rRightHandSideVector)


            Vector temp = ZeroVector(number_of_nodes);

            GetValuesVector(temp);
            // KRATOS_WATCH(temp)
            // RHS = ExtForces - K*temp;
            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
            
            //KRATOS_WATCH(rRightHandSideVector)

            // exit(0);
        }
        KRATOS_CATCH("")
    }

    int SupportSolid3DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportSolid3DCondition" << std::endl;
        return 0;
    }

    void SupportSolid3DCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();

        }
    }

    void SupportSolid3DCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));

        }
    };


    void SupportSolid3DCondition::GetValuesVector(
        Vector& rValues) const  
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];

        }
        
    }
    void SupportSolid3DCondition::CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        )
        {

        }
    void SupportSolid3DCondition::CalculateB(
            Matrix& rB, 
            Matrix& r_DN_DX) const
        {
            const SizeType number_of_control_points = GetGeometry().size();
            const SizeType mat_size = number_of_control_points * 3;

            // KRATOS_WATCH(number_of_control_points)

            if (rB.size1() != 6 || rB.size2() != mat_size)
                rB.resize(6, mat_size);
            noalias(rB) = ZeroMatrix(6, mat_size);

            for (IndexType r = 0; r < mat_size; r++)
            {
                // local node number kr and dof direction dirr
                IndexType kr = r / 3;
                IndexType dirr = r % 3;

                // rB(0, r) = r_DN_DX(kr, 0) * (dirr == 0); // dN/dx
                // rB(1, r) = r_DN_DX(kr, 1) * (dirr == 1); // dN/dy
                // rB(2, r) = r_DN_DX(kr, 2) * (dirr == 2); // dN/dz

                // rB(3, r) = r_DN_DX(kr, 0) * (dirr == 1); // dN/dy
                // rB(3, r) = r_DN_DX(kr, 1) * (dirr == 0); // dN/dx

                // rB(4, r) = r_DN_DX(kr, 0) * (dirr == 2); // dN/dz
                // rB(4, r) = r_DN_DX(kr, 2) * (dirr == 0); // dN/dx

                // rB(5, r) = r_DN_DX(kr, 1) * (dirr == 2); // dN/dz
                // rB(5, r) = r_DN_DX(kr, 2) * (dirr == 1); // dN/dy
                // Filling the B matrix

                if (dirr == 0) {
                    rB(0, r)= r_DN_DX(kr, 0); // dN/dx
                    rB(3, r) = r_DN_DX(kr, 1); // dN/dy
                    rB(4, r) = r_DN_DX(kr, 2); // dN/dz
                }
                else if (dirr == 1) {
                    rB(1, r) = r_DN_DX(kr, 1); // dN/dy
                    rB(3, r) = r_DN_DX(kr, 0); // dN/dx
                    rB(5, r) = r_DN_DX(kr, 2); // dN/dz
                }
                else if (dirr == 2) {
                    rB(2, r) = r_DN_DX(kr, 2); // dN/dz
                    rB(4, r) = r_DN_DX(kr, 0); // dN/dx
                    rB(5, r) = r_DN_DX(kr, 1); // dN/dy
                }
            }
        //KRATOS_WATCH(rB)    
        }


/***********************************************************************************/
/***********************************************************************************/

void SupportSolid3DCondition::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}


bool SupportSolid3DCondition::UseElementProvidedStrain() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void SupportSolid3DCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    
    bool required = false;
        if (mpConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
            required = true;
        }
    if (required) {
        const auto& r_geom = GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        const Properties& r_properties = GetProperties();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geom,r_properties,rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points

        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);
        KRATOS_WATCH(N_values)
        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints();

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, 0, mThisIntegrationMethod);

        // Compute constitutive law variables
        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, integration_points);

        
        // Call the constitutive law to update material variables
        mpConstitutiveLaw->FinalizeMaterialResponse(Values, GetStressMeasure());

        // TODO: Deprecated, remove this
        mpConstitutiveLaw->FinalizeSolutionStep( r_properties, r_geom, row( N_values, 0 ), rCurrentProcessInfo);
        
    }
}    
} // Namespace Kratos
