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
#include "custom_conditions/sbm_load_solid_2D_condition.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

    void SBMLoadSolid2DCondition:: Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        InitializeMaterial();
    }


    void SBMLoadSolid2DCondition::InitializeMaterial()
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


    void SBMLoadSolid2DCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        // loopIdentifier is inner or outer
        std::string loopIdentifier = this->GetValue(IDENTIFIER);
        
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

        // compute surrogate normal and tangent
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 3> normal_parameter_space;
        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude; 
        normal_parameter_space[2] = 0; 

        // retrieve projections
        //------------------------------------------------------------------------------
        Vector projection(2);
        array_1d<double, 3> true_n;
        array_1d<double, 3> true_tau;
        double n_ntilde;
        Vector d(2);
        NodeType projection_node;

        // if (this->GetGeometry().GetValue(NEIGHBOUR_NODES).size() != 0) 
        // {
        //     projection_node = this->GetGeometry().GetValue(NEIGHBOUR_NODES)[0];

        //     true_n = projection_node.GetValue(NORMAL);
        //     true_tau = projection_node.GetValue(LOCAL_TANGENT);
        // }
        if (this->GetValue(NEIGHBOUR_NODES).size() != 0) 
        {
            projection_node = this->GetValue(NEIGHBOUR_NODES)[0];

            true_n = projection_node.GetValue(NORMAL);
            true_tau = projection_node.GetValue(LOCAL_TANGENT);

        } else if (this->GetValue(NEIGHBOUR_CONDITIONS).size() == 2)
        {
            Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
            Condition candidateClosestSkinSegment2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1];

            projection_node = candidateClosestSkinSegment1.GetGeometry()[0];

            // Neumann BC, the true normal is needed
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
            
            // compute true tau (tangential unit vector)
            MathUtils<double>::CrossProduct(true_tau, true_n, vectorOutOfPlane); 
        } else
            KRATOS_ERROR << ":::[SBMLoadSolid2DCondition]::: Error, no NEIGHBOUR NODE OR CONDITIONS SPECIFIED" << std::endl;
        

        projection[0] = projection_node.X();                     
        projection[1] = projection_node.Y() ;
        d[0] = projection[0] - r_geometry.Center().X();          
        d[1] = projection[1] - r_geometry.Center().Y();
        // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
        std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
        outputFile << projection[0] << " " << projection[1] << " "  << r_geometry.Center().X() << " " << r_geometry.Center().Y() <<"\n";
        outputFile.close();

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
        Matrix InvJ0(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        int basisFunctionsOrder;

        Vector meshSize_uv = this->GetValue(MARKER_MESHES);
        double h = std::min(meshSize_uv[0], meshSize_uv[1]);
        if (dim == 3) {h = std::min(h,  meshSize_uv[2]);}

        // Compute basis function order (Note: it is not allow to use different orders in different directions)
        basisFunctionsOrder = std::sqrt(DN_De[0].size1()) - 1;
        
        

        const Matrix& N = r_geometry.ShapeFunctionsValues();

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        Matrix H = ZeroMatrix(1, number_of_nodes);
        Matrix H_sum = ZeroMatrix(1, number_of_nodes);

        // Compute all the derivatives of the basis functions involved
        std::vector<Matrix> nShapeFunctionDerivatives;
        for (int n = 1; n <= basisFunctionsOrder; n++) {
            nShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
        }


        // Neumann (Taylor expansion of the gradient)
        Matrix Hgrad = ZeroMatrix(number_of_nodes, 2);
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i) = N(0, i);
            double H_taylor_term_X = 0.0; // Reset for each node
            double H_taylor_term_Y = 0.0; 
            for (int n = 2; n <= basisFunctionsOrder; n++) {
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
            Hgrad(i,0) = H_taylor_term_X + DN_DX(i,0);
            Hgrad(i,1) = H_taylor_term_Y + DN_DX(i,1);
        }    

        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

        Vector old_displacement(mat_size);
        GetValuesVector(old_displacement);

        //----------------------------------------------

        Matrix B_sum = ZeroMatrix(3,mat_size);

        CalculateB(B_sum, Hgrad);

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

        const Matrix& r_D = Values.GetConstitutiveMatrix();
        //---------------------
        

        //------------------------------
        Matrix DB_sum = prod(r_D, B_sum); //

        Matrix DB = prod(r_D, B); //


        // dot product n cdot n_tilde
        n_ntilde = true_n[0] * normal_parameter_space[0] + true_n[1] * normal_parameter_space[1];
        // double tau_ntilde = true_tau[0] * normal_parameter_space[0] + true_tau[1] * normal_parameter_space[1];


        double integration_factor = IntToReferenceWeight;
        for (IndexType i = 0; i < number_of_nodes; i++) {
            for (IndexType j = 0; j < number_of_nodes; j++) {
                
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int id1 = 2*idim;
                    const int iglob = 2*i+idim;

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int id2 = (id1+2)%3;
                        const int jglob = 2*j+jdim;
                        // -()
                        rLeftHandSideMatrix(iglob, jglob) -= H(0,i)*(DB(id1, jglob)* normal_parameter_space[0] + DB(id2, jglob)* normal_parameter_space[1]) * integration_factor;

                        rLeftHandSideMatrix(iglob, jglob) += H(0,i)*(DB_sum(id1, jglob)* true_n[0] + DB_sum(id2, jglob)* true_n[1]) * integration_factor * n_ntilde;
                    }

                }
            }
        }

        // // Assembly     
        if (CalculateResidualVectorFlag) {

            double nu = this->GetProperties().GetValue(POISSON_RATIO);
            double E = this->GetProperties().GetValue(YOUNG_MODULUS);
            Vector g_N = ZeroVector(2);

            // When "analysis_type" is "linear" temper = 0
            // const double x = projection[0];
            // const double y = projection[1];

            // g_N[0] = E/(1+nu)*(-sin(x)*sinh(y)) * true_n[0] + E/(1+nu)*(cos(x)*cosh(y)) * true_n[1]; 
            // g_N[1] = E/(1+nu)*(cos(x)*cosh(y)) * true_n[0] + E/(1+nu)*(sin(x)*sinh(y)) * true_n[1]; 

            // g_N[0] = E/(1-nu)*(sin(x)*sinh(y)) * true_n[0]; 
            // g_N[1] = E/(1-nu)*(sin(x)*sinh(y))  * true_n[1]; 


            // g_N[0] = this->GetValue(FORCE_X); 
            // g_N[1] = this->GetValue(FORCE_Y); 

            g_N[0] = projection_node.GetValue(FORCE_X);
            g_N[1] = projection_node.GetValue(FORCE_Y);
            
            for (IndexType i = 0; i < number_of_nodes; i++) {
                
                for (IndexType zdim = 0; zdim < 2; zdim++) {
                    
                    rRightHandSideVector[2*i+zdim] += H(0,i)*g_N[zdim] * n_ntilde  * IntToReferenceWeight;

                }
            }

            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,old_displacement);
            
        }

        for (unsigned int i = 0; i < number_of_nodes; i++) {

                std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
                outputFile << r_geometry[i].GetId() << "  " <<r_geometry[i].GetDof(DISPLACEMENT_X).EquationId() <<"\n";
                outputFile.close();
            }
        KRATOS_CATCH("")
    }

    int SBMLoadSolid2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    unsigned long long SBMLoadSolid2DCondition::factorial(int n) 
    {
        if (n == 0) return 1;
        unsigned long long result = 1;
        for (int i = 2; i <= n; ++i) result *= i;
        return result;
    }

    // Function to compute a single term in the Taylor expansion
    double SBMLoadSolid2DCondition::computeTaylorTerm(double derivative, double dx, int n_k, double dy, int k)
    {
        return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (factorial(k) * factorial(n_k));    
    }


    void SBMLoadSolid2DCondition::EquationIdVector(
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

    void SBMLoadSolid2DCondition::GetDofList(
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


    void SBMLoadSolid2DCondition::GetValuesVector(
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

    void SBMLoadSolid2DCondition::CalculateB(
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

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters SBMLoadSolid2DCondition::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        std::ifstream infile(rDataFileName);

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    };


} // Namespace Kratos
