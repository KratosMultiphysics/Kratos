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
#include "custom_conditions/sbm_plain_stress_condition.h"

namespace Kratos
{
    void SBMPlainStressCondition::CalculateAll(
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

        // Read the refinements.iga.json
        const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");
        int insertions = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        double h = 2.0/(insertions+1) ;

        int basisFunctionsOrder = refinements_parameters["refinements"][0]["parameters"]["increase_degree_u"].GetInt()+1;

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


        // TAYLOR EXPANSION TERM ----------------------------------------------------------------------

        Matrix H_sum = ZeroMatrix(1, number_of_nodes);

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
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
            std::vector<Matrix> nShapeFunctionDerivatives;
            for (int n = 1; n <= basisFunctionsOrder; n++) {
                nShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, point_number, this->GetIntegrationMethod()));
            }

            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                double H_taylor_term = 0.0; // Reset for each node
                for (int n = 1; n <= basisFunctionsOrder; n++) {
                    // Retrieve the appropriate derivative for the term
                    Matrix& shapeFunctionDerivatives = nShapeFunctionDerivatives[n-1];
                    for (int k = 0; k <= n; k++) {
                        int n_k = n - k;
                        double derivative = shapeFunctionDerivatives(i,k); 
                        // Compute the Taylor term for this derivative
                        H_taylor_term += computeTaylorTerm(derivative, d[0], n_k, d[1], k);
                    }
                }
                H_sum(0,i) = H_taylor_term + N(point_number, i);;
            }
        }
            
        
            // ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ



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

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        Vector GP_parameter_coord(2); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]);
        
        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);

        // MODIFIED
        double nu = this->GetProperties().GetValue(POISSON_RATIO);
        double E = this->GetProperties().GetValue(YOUNG_MODULUS);
        Matrix D = ZeroMatrix(3,3);
        D(0,0) = 1; 
        D(0,1) = nu;
        D(1,0) = nu;
        D(1,1) = 1;
        D(2,2) = (1-nu)/2;
        D *= E/(1-nu*nu);


        // Stampa su file esterno le coordinate (projection[0],projection[1])
        std::ofstream outputFile("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << GP_parameter_coord[0] << " " << GP_parameter_coord[1]  <<"\n";
        outputFile.close();
        
        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            Matrix Jacobian = ZeroMatrix(2,2);
            Jacobian(0,0) = J0[point_number](0,0);
            Jacobian(0,1) = J0[point_number](0,1);
            Jacobian(1,0) = J0[point_number](1,0);
            Jacobian(1,1) = J0[point_number](1,1);

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
            noalias(DN_DX) = prod(DN_De[point_number],InvJ0_23);
            
            Matrix H = ZeroMatrix(1, number_of_nodes);
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i)            = N(point_number, i);
            }


            // Differential area
            double penalty_integration = penalty * integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);

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

            Matrix DB = prod(D,B);
            double integration_factor = integration_points[point_number].Weight() * std::abs(determinant_jacobian_vector[point_number]);
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

                // const double temperature = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(TEMPERATURE);

                u_D[0] = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(DISPLACEMENT_X);
                u_D[1] = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(DISPLACEMENT_Y);
                


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
                
                Vector temp = ZeroVector(number_of_nodes);

                GetValuesVector(temp);

                // RHS = ExtForces - K*temp;
                noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
                
                // exit(0);
            }
        }
        KRATOS_CATCH("")
    }

    int SBMPlainStressCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SBMPlainStressCondition::EquationIdVector(
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

    void SBMPlainStressCondition::GetDofList(
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


    void SBMPlainStressCondition::GetValuesVector(
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
    Parameters SBMPlainStressCondition::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        std::ifstream infile(rDataFileName);

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    };


    void SBMPlainStressCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
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

    unsigned long long SBMPlainStressCondition::factorial(int n) 
    {
        if (n == 0) return 1;
        unsigned long long result = 1;
        for (int i = 2; i <= n; ++i) result *= i;
        return result;
    }

    // Function to compute a single term in the Taylor expansion
    double SBMPlainStressCondition::computeTaylorTerm(double derivative, double dx, int n_k, double dy, int k)
    {
        return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (factorial(k) * factorial(n_k));    
    }

} // Namespace Kratos
