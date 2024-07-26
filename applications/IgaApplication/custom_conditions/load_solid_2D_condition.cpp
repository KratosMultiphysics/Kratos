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
#include "custom_conditions/load_solid_2D_condition.h"

namespace Kratos
{
    void LoadSolid2DCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

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

        // Stampa su file esterno le coordinate (projection[0],projection[1])
        std::ofstream outputFile("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << GP_parameter_coord[0] << " " << GP_parameter_coord[1]  <<"\n";
        outputFile.close();

        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

        
        
        const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

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

        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0_23);
        
        Matrix H = ZeroMatrix(1, number_of_nodes);
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i)            = N(0, i);
        }


        // // Assembly     
        if (CalculateResidualVectorFlag) {

            double nu = this->GetProperties().GetValue(POISSON_RATIO);
            double E = this->GetProperties().GetValue(YOUNG_MODULUS);
            Vector g_N = ZeroVector(2);

            // When "analysis_type" is "linear" temper = 0
            const double x = GP_parameter_coord[0];
            const double y = GP_parameter_coord[1];

            // g_N[0] = E/(1+nu)*(-sin(x)*sinh(y)) * normal_parameter_space[0] + E/(1+nu)*(cos(x)*cosh(y)) * normal_parameter_space[1]; 
            // g_N[1] = E/(1+nu)*(cos(x)*cosh(y)) * normal_parameter_space[0] + E/(1+nu)*(sin(x)*sinh(y)) * normal_parameter_space[1]; 

            g_N[0] = this->GetValue(FORCE_X); 
            g_N[1] = this->GetValue(FORCE_Y); 
            
            for (IndexType i = 0; i < number_of_nodes; i++) {
                
                for (IndexType zdim = 0; zdim < 2; zdim++) {
                    
                    rRightHandSideVector[2*i+zdim] += H(0,i)*g_N[zdim] * IntToReferenceWeight;

                }
            }
            
            Vector temp = ZeroVector(number_of_nodes);


            GetValuesVector(temp);

            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
            
        }
        KRATOS_CATCH("")
    }

    int LoadSolid2DCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void LoadSolid2DCondition::EquationIdVector(
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

    void LoadSolid2DCondition::GetDofList(
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

    void LoadSolid2DCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
       
    /////////////////////////
        const auto& r_geometry = GetGeometry();
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

        // std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
        // if (output_file.is_open()) {
        //     output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        //     output_file << rOutput << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points[0].Weight() << std::endl;
        //     output_file.close();
        // } 


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


    void LoadSolid2DCondition::GetValuesVector(
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

} // Namespace Kratos
