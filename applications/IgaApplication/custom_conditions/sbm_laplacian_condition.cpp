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
//                   
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
        Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
        // Condition candidateClosestSkinSegment2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1];

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
        Matrix DN_DX(number_of_nodes,dim);
        Matrix DN_DPSI(number_of_nodes,dim);
        Matrix InvJ0(dim,dim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;


        // Differential area
        double penalty_integration = penalty * integration_points[0].Weight() * std::abs(determinant_jacobian_vector[0]);
        
        // Obtaining the projection from the closest skin segment
        Vector projection(2);
        projection[0] = candidateClosestSkinSegment1.GetGeometry()[0].X() ;
        projection[1] = candidateClosestSkinSegment1.GetGeometry()[0].Y() ;

        // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
        std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
        outputFile << projection[0] << " " << projection[1] << " "  << r_geometry.Center().X() << " " << r_geometry.Center().Y() <<"\n";
        outputFile.close();

        d.resize(2);
        d[0] = projection[0] - r_geometry.Center().X();
        d[1] = projection[1] - r_geometry.Center().Y();
        // d[0] = 0;
        // d[1] = 0;

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

        // Calculating the PARAMETER SPACE derivatives
        Matrix Identity_Matrix = ZeroMatrix(2,2);
        Identity_Matrix(0,0) = 1.0;
        Identity_Matrix(1,1) = 1.0;
        noalias(DN_DPSI) = prod(DN_De[0],Identity_Matrix);

        // Compute the normals
        array_1d<double, 3> tangent_parameter_space;
        array_1d<double, 2> normal_physical_space;

        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
        double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
        
        // NEW FOR GENERAL JACOBIAN
        normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
        normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude; 

        normal_physical_space = prod(trans(J0[0]),normal_parameter_space);
        
        // Collins, Lozinsky & Scovazzi innovation
        double Guglielmo_innovation = -1.0;  // = 1 -> Penalty approach
                                            // = -1 -> Free-penalty approach
        if (penalty == -1.0) {
            penalty_integration = 0.0;
            Guglielmo_innovation = -1.0;
        }

        Matrix H = ZeroMatrix(1, number_of_nodes);
        Vector H_vec = ZeroVector(number_of_nodes);
        Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
        Matrix DN_dot_n_parameter = ZeroMatrix(1, number_of_nodes);

        // Matrix H_sum = ZeroMatrix(1, number_of_nodes);
        Vector H_sum_vec = ZeroVector(number_of_nodes);

        // Compute all the derivatives of the basis functions involved
        std::vector<Matrix> nShapeFunctionDerivatives;
        for (int n = 1; n <= basisFunctionsOrder; n++) {
            nShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
        }

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)               = N(0, i);
            H_vec(i)              = N(0, i);
            DN_dot_n(0, i)        = DN_DX(i, 0) * normal_physical_space[0] + DN_DX(i, 1) * normal_physical_space[1] ;
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
            H_sum(0,i) = H_taylor_term + H(0,i);
            H_sum_vec(i) = H_taylor_term + H(0,i);
        }

        // Assembly
        // Termine -(GRAD_w * n, u + GRAD_u * d + ...)
        noalias(rLeftHandSideMatrix) -= Guglielmo_innovation * prod(trans(DN_dot_n), H_sum)  * integration_points[0].Weight() * std::abs(determinant_jacobian_vector[0]) ;
        // Termine -(w,GRAD_u * n) from integration by parts -> Fundamental !! 
        noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n)                             * integration_points[0].Weight() * std::abs(determinant_jacobian_vector[0]) ;
        // SBM terms (Taylor Expansion) + alpha * (w + GRAD_w * d + ..., u + GRAD_u * d + ...)
        noalias(rLeftHandSideMatrix) += prod(trans(H_sum), H_sum) * penalty_integration ;


        if (CalculateResidualVectorFlag) {
            
            // double temperature = GP_parameter_coord[0]-GP_parameter_coord[1];
            // const double temperature = sin(projection[0]) * sinh(projection[1]) ;

            

            // const double temperature = sin(r_geometry.Center().X()) * sinh(r_geometry.Center().Y()) ;
            // const double temperature = projection[0]*projection[0]*projection[0] + projection[1]*projection[1]*projection[1] ;
            // const double temperature = projection[0]*projection[0]*projection[0]*projection[0] + projection[1]*projection[1]*projection[1]*projection[1] ;
            
            // FINAL ONE:
            const double temperature = candidateClosestSkinSegment1.GetGeometry()[0].GetValue(TEMPERATURE);
            
            Vector u_D(number_of_nodes);

            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                // const double temper = 0; // r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);      // Always = 0.  
                // u_D[i] = (temper - temperature);
                u_D[i] = -temperature ;
            }
            // noalias(rRightHandSideVector) -= prod(prod(trans(H_sum), H), u_D) * penalty_integration;
            noalias(rRightHandSideVector) -= H_sum_vec * (-temperature) * penalty_integration;
            // Dirichlet BCs
            noalias(rRightHandSideVector) += Guglielmo_innovation * prod(prod(trans(DN_dot_n), H), u_D) * integration_points[0].Weight() * std::abs(determinant_jacobian_vector[0]);



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
        const SizeType nb_nodes = r_geometry.size();
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
        // Shape function values
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const ProcessInfo& r_process_info = rCurrentProcessInfo;
        const Variable<double>& r_unknown_var = TEMPERATURE;

        double rOutput = 0;
        for (IndexType i = 0; i < r_geometry.size(); ++i)
        {
            double output_solution_step_value = r_geometry[i].GetSolutionStepValue(r_unknown_var);
            rOutput = rOutput + H_sum(0,i)*output_solution_step_value;
        }
        
        // the true normal is needed
        Condition candidateClosestSkinSegment1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
        Condition candidateClosestSkinSegment2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1] ;
        array_1d<double,3> vectorSkinSegment1 = candidateClosestSkinSegment1.GetGeometry()[1] - candidateClosestSkinSegment1.GetGeometry()[0];
        array_1d<double,3> vectorSkinSegment2 = candidateClosestSkinSegment2.GetGeometry()[1] - candidateClosestSkinSegment2.GetGeometry()[0];
        array_1d<double,3> vectorOutOfPlane = ZeroVector(3);
        vectorOutOfPlane[2] = 1.0;
        array_1d<double,3> crossProductSkinSegment1;
        array_1d<double,3> crossProductSkinSegment2; 
        MathUtils<double>::CrossProduct(crossProductSkinSegment1, vectorOutOfPlane, vectorSkinSegment1);
        MathUtils<double>::CrossProduct(crossProductSkinSegment2, vectorOutOfPlane, vectorSkinSegment2);
        
        array_1d<double, 3> true_n = crossProductSkinSegment1 / MathUtils<double>::Norm(crossProductSkinSegment1) + crossProductSkinSegment2 / MathUtils<double>::Norm(crossProductSkinSegment2);
        std::string loopIdentifier = this->GetValue(IDENTIFIER);
        if (loopIdentifier == "inner") {
            true_n = true_n / MathUtils<double>::Norm(true_n) ;
        } else { // outer
            true_n = - true_n / MathUtils<double>::Norm(true_n) ;
        }
        // dot product n cdot n_tilde
        double n_ntilde = true_n[0] * normal_parameter_space[0] + true_n[1] * normal_parameter_space[1];
        // KRATOS_WATCH('\n')
        // KRATOS_WATCH(r_geometry.Center())
        // KRATOS_WATCH(true_n)
        // KRATOS_WATCH(normal_parameter_space)
        // // double n_x_ntilde = std::abs(true_n[0] * normal_parameter_space[1] - true_n[1] * normal_parameter_space[0]);
        // KRATOS_WATCH(n_ntilde)
        // KRATOS_WATCH(std::sqrt(1-n_ntilde*n_ntilde))
        // rOutput = rOutput*n_ntilde;

        double x1 = candidateClosestSkinSegment2.GetGeometry()[0].X();
        double y1 = candidateClosestSkinSegment2.GetGeometry()[0].Y();
        double x2 = candidateClosestSkinSegment2.GetGeometry()[1].X();
        double y2 = candidateClosestSkinSegment2.GetGeometry()[1].Y();
        double x3 = candidateClosestSkinSegment1.GetGeometry()[1].X();
        double y3 = candidateClosestSkinSegment1.GetGeometry()[1].Y();
        // double x1 = 1.25;
        // double y1 = 1.0;
        // double x2 = 1.0;
        // double y2 = 1.25;
        // double x3 = 0.75;
        // double y3 = 1.0;
        KRATOS_WATCH(candidateClosestSkinSegment2.GetGeometry()[0])
        KRATOS_WATCH(candidateClosestSkinSegment2.GetGeometry()[1])
        KRATOS_WATCH(candidateClosestSkinSegment1.GetGeometry()[1])

        double numA = (x1*x1 + y1*y1) * (y2 - y3) + (x2*x2 + y2*y2) * (y3 - y1) + (x3*x3 + y3*y3) * (y1 - y2);
        double numB = (x1*x1 + y1*y1) * (x3 - x2) + (x2*x2 + y2*y2) * (x1 - x3) + (x3*x3 + y3*y3) * (x2 - x1);
        double den = 2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2);
        
        KRATOS_WATCH(den)
        double curvature;
        if (std::abs(den) < 1e-15) {
            curvature = 0;
        }
        else {
            double x_c = numA / den;
            double y_c = numB / den;
            KRATOS_WATCH(x_c)
            KRATOS_WATCH(y_c)
            curvature = 1/(std::sqrt((x2-x_c)*(x2-x_c)+(y2-y_c)*(y2-y_c))) ;
        }

        double factor = n_ntilde / (1+norm_2(d)*curvature);
        
        KRATOS_WATCH(1/curvature)
        // exit(0);

        std::ofstream output_file("txt_files/boundary_results_GPs.txt", std::ios::app);
        if (output_file.is_open()) {
            output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
            output_file << rOutput << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << integration_points[0].Weight()*factor << " " << d[0] << " " << d[1] << std::endl;
            output_file.close();
        }
    }


} // Namespace Kratos
