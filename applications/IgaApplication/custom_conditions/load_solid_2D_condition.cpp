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
     void LoadSolid2DCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        InitializeMaterial();
    }


    void LoadSolid2DCondition::InitializeMaterial()
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


    void LoadSolid2DCondition::CalculateB(
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

    void LoadSolid2DCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){

        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
    }


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
        Matrix DN_DX(number_of_nodes,2);
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
        GP_parameter_coord = r_geometry.Center();
        
        // Stampa su file esterno le coordinate (projection[0],projection[1])
        std::ofstream outputFile("txt_files/boundary_GPs.txt", std::ios::app);
        outputFile << std::setprecision(14); // Set precision to 10^-14
        outputFile << GP_parameter_coord[0] << " " << GP_parameter_coord[1]  <<"\n";
        outputFile.close();

        const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        Vector add_factor = prod(Jacobian, tangent_parameter_space);

        DetJ0 = norm_2(add_factor);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)

        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

        const double IntToReferenceWeight = integration_points[0].Weight() * std::abs(DetJ0) * thickness;

        SetValue(INTEGRATION_WEIGHT, IntToReferenceWeight);

        normal_physical_space = prod(trans(InvJ0),normal_parameter_space);

        normal_physical_space /= norm_2(normal_physical_space);
        
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
            // const double x = GP_parameter_coord[0];
            // const double y = GP_parameter_coord[1];

            // g_N[0] = E/(1+nu)*(-sin(x)*sinh(y)) * normal_parameter_space[0] + E/(1+nu)*(cos(x)*cosh(y)) * normal_parameter_space[1]; 
            // g_N[1] = E/(1+nu)*(cos(x)*cosh(y)) * normal_parameter_space[0] + E/(1+nu)*(sin(x)*sinh(y)) * normal_parameter_space[1]; 

            // g_N[0] = E/(1-nu)*(sin(x)*sinh(y)) * normal_physical_space[0]; 
            // g_N[1] = E/(1-nu)*(sin(x)*sinh(y))  * normal_physical_space[1]; 

            g_N[0] = this->GetValue(FORCE_X); 
            g_N[1] = this->GetValue(FORCE_Y); 

            // SetValue(NORMAL_STRESS, g_N);
            
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
       ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        
    // /////////////////////////
    //     const auto& r_geometry = GetGeometry();
    //     const SizeType nb_nodes = r_geometry.size();

    //     // Integration Points
    //     const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    //     // Shape function values
    //     const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    //     GeometryType::JacobiansType J0;
    //     r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    //     // Get the parameter coordinates
    //     Vector GP_parameter_coord(2); 
    //     GP_parameter_coord = r_geometry.Center(); // Only one Integration Points 

    //     double x_coord_gauss_point = 0;
    //     double y_coord_gauss_point = 0;
    //     double rOutput = 0;

    //     for (IndexType i = 0; i < nb_nodes; ++i)
    //     {
    //         // KRATOS_WATCH(r_geometry[i])
    //         double output_solution_step_value = r_geometry[i].GetSolutionStepValue(DISPLACEMENT_X);
    //         rOutput += r_N(0, i) * output_solution_step_value;
    //         x_coord_gauss_point += r_N(0, i) * r_geometry[i].X0();
    //         y_coord_gauss_point += r_N(0, i) * r_geometry[i].Y0();
    //     }        

    //     // ########################################################################333
    //     double x = GP_parameter_coord[0];
    //     double y = GP_parameter_coord[1];
    //     // double true_sol = -cos(x)*sinh(y);
    //     Vector approx_lhs = ZeroVector(9);
    //     Vector approx_rhs = ZeroVector(9);

    //     for (IndexType i = 0; i < 9; i++) {
    //         for (IndexType j = 0; j < 9; j++) {
                
    //             for (IndexType idim = 0; idim < 1; idim++) {
    //                 const int id1 = 2*idim;
    //                 const int iglob = 2*i+idim;

    //                 approx_lhs[i] -= r_N(0,i)*r_N(0,j)* r_geometry[j].GetSolutionStepValue(DISPLACEMENT_X);
    //             }
    //         }
    //     }

    //     Vector u_D(2);
    //     u_D[0] = this->GetValue(DISPLACEMENT_X);
    //     u_D[1] = this->GetValue(DISPLACEMENT_Y);

    //     for (IndexType i = 0; i < 9; i++) {

    //         for (IndexType idim = 0; idim < 1; idim++) {

    //             approx_rhs[i] -= r_N(0,i)*u_D[idim];
    //         }
    //     }

    //     // KRATOS_WATCH(approx_lhs - approx_rhs)

    //     // ########################################################################333

    //     // KRATOS_WATCH(err)

    //     // std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
    //     // if (output_file.is_open()) {
    //     //     output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
    //     //     output_file << rOutput << " " << x_coord_gauss_point << " " << y_coord_gauss_point << " " <<integration_points[0].Weight() << std::endl;
    //     //     output_file.close();
    //     // } 


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

        normal_physical_space = prod(trans(InvJ0),normal_parameter_space);

        normal_physical_space /= norm_2(normal_physical_space);

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
