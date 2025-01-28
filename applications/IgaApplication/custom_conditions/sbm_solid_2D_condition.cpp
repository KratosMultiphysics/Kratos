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
#include "custom_conditions/sbm_solid_2D_condition.h"

namespace Kratos
{
    void SBMSolid2DCondition:: Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        InitializeMaterial();
    }


    void SBMSolid2DCondition::InitializeMaterial()
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

    void SBMSolid2DCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
        // Condition candidateClosestSkinSegment2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1];
        const auto& r_geometry = this->GetGeometry();
        const SizeType number_of_nodes = r_geometry.PointsNumber();

        const SizeType mat_size = number_of_nodes * 2;
        //resizing as needed the LHS
        if(rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size,mat_size,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
        
        // resizing as needed the RHS
        if(rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size,false);
        noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

        double penalty = GetProperties()[PENALTY_FACTOR];

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
        const unsigned int dim = DN_De[0].size2();
        
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

        // Integration
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        
        r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);  // = 1

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX

        // TAYLOR EXPANSION TERM ----------------------------------------------------------------------

        Matrix H_sum = ZeroMatrix(1, number_of_nodes);

        // Obtaining the projection from the closest skin segment
        Vector projection(2);
        projection[0] = candidateClosestSkinSegment1.GetGeometry()[0].X() ;
        projection[1] = candidateClosestSkinSegment1.GetGeometry()[0].Y() ;

        // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
        std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
        outputFile << projection[0] << " " << projection[1] << " "  << r_geometry.Center().X() << " " << r_geometry.Center().Y() <<"\n";
        outputFile.close();

        Vector d(2);
        d[0] = projection[0] - r_geometry.Center().X();
        d[1] = projection[1] - r_geometry.Center().Y();
        // d[0] = 0;
        // d[1] = 0;
        const Matrix& N = r_geometry.ShapeFunctionsValues();

        // Compute all the derivatives of the basis functions involved
        for (int n = 1; n <= mbasisFunctionsOrder; n++) {
            mShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
        }

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
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
            H_sum(0,i) = H_taylor_term + N(0, i);;
        }
            
        Matrix DN_DX(number_of_nodes,3);
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

        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        Vector GP_parameter_coord(2); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]);
        
        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);

        // MODIFIED

        Vector old_displacement(mat_size);
        GetValuesVector(old_displacement);


        // Stampa su file esterno le coordinate (projection[0],projection[1])
        std::ofstream outputFile1("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile1 << std::setprecision(14); // Set precision to 10^-14
        outputFile1 << GP_parameter_coord[0] << " " << GP_parameter_coord[1]  <<"\n";
        outputFile1.close();
        
        // const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant 
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
        Matrix InvJ0_23 = ZeroMatrix(2,3);
        InvJ0_23(0,0) = InvJ0(0,0);
        InvJ0_23(0,1) = InvJ0(0,1);
        InvJ0_23(1,0) = InvJ0(1,0);
        InvJ0_23(1,1) = InvJ0(1,1);
        InvJ0_23(0,2) = 0;
        InvJ0_23(1,2) = 0;

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0_23);
        
        Matrix H = ZeroMatrix(1, number_of_nodes);
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)            = N(0, i);
        }


        // Differential area

        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;
        double penalty_integration = penalty * IntToReferenceWeight;

        // Guglielmo innovaction
        double Guglielmo_innovation = -1.0;  // = 1 -> Penalty approach
                                                // = -1 -> Free-penalty approach
        if (Guglielmo_innovation == -1.0) {
            penalty_integration = 0.0;
        }

        // COMPUTE THE EXTENSIONS OF THE BASIS FUNCTIONS FROM SURROGOGATE -> TRUE


        // Assembly

        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX);

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
        Values.SetConstitutiveMatrix(this_constitutive_variables.ConstitutiveMatrix);

        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2); 

        const Vector& r_stress_vector = Values.GetStressVector();
        const Matrix& r_D = Values.GetConstitutiveMatrix();
        //-----------------------------------------------------------------------------------

        Matrix DB = prod(r_D,B);
        double integration_factor = IntToReferenceWeight;
        for (IndexType i = 0; i < number_of_nodes; i++) {
            for (IndexType j = 0; j < number_of_nodes; j++) {
                
                for (IndexType idim = 0; idim < 2; idim++) {
                    rLeftHandSideMatrix(2*i+idim, 2*j+idim) -= H_sum(0,i)*H_sum(0,j)* penalty_integration;
                    const int id1 = 2*idim;
                    const int iglob = 2*i+idim;

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int id2 = (id1+2)%3;
                        const int jglob = 2*j+jdim;
                        rLeftHandSideMatrix(iglob, jglob) -= H(0,i)*(DB(id1, jglob)* normal_parameter_space[0] + DB(id2, jglob)* normal_parameter_space[1]) * integration_factor;

                        rLeftHandSideMatrix(iglob, jglob) -= Guglielmo_innovation*H_sum(0,j)*(DB(id1, 2*i+jdim)* normal_parameter_space[0] + DB(id2, 2*i+jdim)* normal_parameter_space[1]) * integration_factor;
                    }

                }
            }
        }
        
        // // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
        // // Of the Dirichlet BCs -(GRAD_w* n,u) 
        
        if (CalculateResidualVectorFlag) {
            
            // const double& temperature = Has(TEMPERATURE)
            //     ? this->GetValue(TEMPERATURE)
            //     : 0.0;
            
            Vector u_D = ZeroVector(2); //->GetValue(DISPLACEMENT);
            

            u_D[0] = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(DISPLACEMENT_X);
            u_D[1] = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(DISPLACEMENT_Y);

            // double x = GetGeometry().Center().X(); double y = GetGeometry().Center().Y();

            // u_D[0] = -cos(x)*sinh(y);
            // u_D[1] = sin(x)*cosh(y);
            


            for (IndexType i = 0; i < number_of_nodes; i++) {

                for (IndexType idim = 0; idim < 2; idim++) {

                    rRightHandSideVector[2*i+idim] -= H_sum(0,i)*u_D[idim]* penalty_integration;
                    const int id1 = idim*2;

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int id2 = (id1+2)%3;
                        rRightHandSideVector(2*i+idim) -= Guglielmo_innovation*u_D[jdim]*(DB(id1, 2*i+jdim)* normal_parameter_space[0] + DB(id2, 2*i+jdim)* normal_parameter_space[1]) * integration_factor;
                    }

                }


            }
            
            // noalias(rRightHandSideVector) += prod(prod(trans(H), H), u_D) * penalty_integration;
            // // Of the Dirichlet BCs
            // noalias(rRightHandSideVector) += Guglielmo_innovation * prod(prod(trans(DN_dot_n), H), u_D) * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);
            
            // Vector temp = ZeroVector(number_of_nodes);

            // GetValuesVector(temp);

            // RHS = ExtForces - K*temp;
            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,old_displacement);
            
            // exit(0);
        }
        KRATOS_CATCH("")
    }

    int SBMSolid2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SBMSolid2DCondition::EquationIdVector(
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

    void SBMSolid2DCondition::GetDofList(
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


    void SBMSolid2DCondition::GetValuesVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 2 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters SBMSolid2DCondition::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        std::ifstream infile(rDataFileName);

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    };


    void SBMSolid2DCondition::CalculateB(
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

    unsigned long long SBMSolid2DCondition::factorial(int n) 
    {
        if (n == 0) return 1;
        unsigned long long result = 1;
        for (int i = 2; i <= n; ++i) result *= i;
        return result;
    }

    // Function to compute a single term in the Taylor expansion
    double SBMSolid2DCondition::computeTaylorTerm(double derivative, double dx, int n_k, double dy, int k)
    {
        return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (factorial(k) * factorial(n_k));    
    }


    void SBMSolid2DCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        double sol_on_surrogate_x = 0.0;
        double sol_on_surrogate_y = 0.0;

        for (IndexType i = 0; i < r_geometry.size(); ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value_x = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_X);
            double output_solution_step_value_y = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_Y);
            sol_on_surrogate_x += r_N(0, i) * output_solution_step_value_x;
            sol_on_surrogate_y += r_N(0, i) * output_solution_step_value_y;
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

            double solution_on_true_x = 0.0; 
            double solution_on_true_y = 0.0;

            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                double output_solution_step_value_x = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_X);
                double output_solution_step_value_y = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_Y);

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

                solution_on_true_x += (r_N(0,i) + H_taylor_term) * output_solution_step_value_x;
                solution_on_true_y += (r_N(0,i) + H_taylor_term) * output_solution_step_value_y;
            } 
            output_file << solution_on_true_x << " " << solution_on_true_y << " "<< integration_weight_list_on_true_boundary[i_gauss] << " " 
                            << curr_gauss_point_on_true[0] << " " << curr_gauss_point_on_true[1]  
                            << " " << curr_gauss_point_on_true[2]   << std::endl;
        }
        output_file.close();

    }

} // Namespace Kratos
