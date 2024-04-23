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

// System includes

// External includes

// Project includes
#include "custom_conditions/sbm_support_lagrange_condition.h"


namespace Kratos
{
    void SBMSupportLagrangeCondition::CalculateAll(
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
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();
        const SizeType mat_size = number_of_non_zero_nodes * 2;

        Matrix LHS = ZeroMatrix(mat_size, mat_size);
        Vector u(mat_size);

        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        // Integration
        const typename GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

        // Determinant of jacobian
        // Determine the integration: conservative -> initial; non-conservative -> current
        Vector determinant_jacobian_vector(integration_points.size());
        r_geometry.DeterminantOfJacobian(determinant_jacobian_vector);

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        Vector GP_parameter_coord(2); 
        GP_parameter_coord = prod(r_geometry.Center(),J0[0]);


        //__________________________________________________________________________________________________________________________________

        // GENERAL COMPUTATION OF THE PROJECTION
        std::ifstream file("txt_files/true_points.txt");    // Read the true points from the mdpa
        std::vector<double> x_true_boundary;
        std::vector<double> y_true_boundary;
        double x, y;
        while (file >> x >> y) { x_true_boundary.push_back(x); y_true_boundary.push_back(y); }
        file.close();
        int index_min_distance = 0; // Initialization of the index of the true node closest to the Gauss point.
        double min_distance_squared = 1e14;
        for (int i = 1; i < x_true_boundary.size(); i++) {
            double current_distance_squared = (GP_parameter_coord[0]-x_true_boundary[i])*(GP_parameter_coord[0]-x_true_boundary[i]) +  
                    (GP_parameter_coord[1]-y_true_boundary[i])*(GP_parameter_coord[1]-y_true_boundary[i]) ;
            if ( current_distance_squared  <  min_distance_squared ) {
                    min_distance_squared = current_distance_squared;
                    index_min_distance = i ;
                    }
        }
        Vector projection(2);
        projection[0] = x_true_boundary[index_min_distance] ;
        projection[1] = y_true_boundary[index_min_distance] ;
        
        // Stampa su file esterno le coordinate (projection[0],projection[1])
        std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
        outputFile << projection[0] << " " << projection[1] << " "  << GP_parameter_coord[0] << " " << GP_parameter_coord[1] <<"\n";
        outputFile.close();

        Vector d(2);
        d[0] = projection[0] - GP_parameter_coord[0];
        d[1] = projection[1] - GP_parameter_coord[1];

        const unsigned int dim = 2;
        Matrix InvJ0(dim,dim);
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double DetJ0;

        const Matrix& N = r_geometry.ShapeFunctionsValues();

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        const Matrix& DN_De = r_geometry.ShapeFunctionDerivatives(1, 0, this->GetIntegrationMethod());
        const Matrix& DDN_DDe = r_geometry.ShapeFunctionDerivatives(2, 0, this->GetIntegrationMethod());
        const Matrix& DDDN_DDDe = r_geometry.ShapeFunctionDerivatives(3, 0, this->GetIntegrationMethod());
        //__________________________________________________________________________________________________________________________________
        
        
        // non zero node counter for Lagrange Multipliers.
        IndexType counter_n = 0;
        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            // Differential area, being 1 for points.
            const double integration = integration_points[point_number].Weight() * determinant_jacobian_vector[point_number];

            //__________________________________________________________________________________________________________________________________
            
            Matrix H = ZeroMatrix(1, number_of_nodes);
            Matrix H_gradient_term = ZeroMatrix(1, number_of_nodes);
            Matrix H_hessian_term = ZeroMatrix(1, number_of_nodes);
            Matrix H_3rdTayor_term = ZeroMatrix(1, number_of_nodes);
            
            for (IndexType i = 0; i < number_of_nodes; ++i)
            {
                H(0, i)               = N(point_number, i);
                // Taylor expansion in the parameter space
                H_gradient_term(0, i) = DN_De(i,0) * d[0] + DN_De(i,1) * d[1] ; 
                H_hessian_term(0, i) = 0.5 * ( d[0]*d[0]*DDN_DDe(i,0) + 2.0 * d[0]*d[1]*DDN_DDe(i,1) + d[1]*d[1]*DDN_DDe(i,2) ) ;
                H_3rdTayor_term(0, i) = 1.0/6.0 * DDDN_DDDe(i,0) * d[0]*d[0]*d[0] + 3.0/6.0 * DDDN_DDDe(i,1)*d[0]*d[0]*d[1] + 3.0/6.0* DDDN_DDDe(i,2)*d[0]*d[1]*d[1] + 1.0/6.0 *DDDN_DDDe(i,3)*d[1]*d[1]*d[1] ;
            }
            //__________________________________________________________________________________________________________________________________
            // loop over Lagrange Multipliers
            for (IndexType i = 0; i < number_of_nodes; i++) {
                // non zero node counter for for displacements.
                IndexType counter_m = 0;
                if (r_N(point_number, i) > shape_function_tolerance) { // || H_gradient_term(0, i) > shape_function_tolerance || H_hessian_term(0, i) > shape_function_tolerance || H_3rdTayor_term(0, i) > shape_function_tolerance) { 
                    // loop over shape functions of displacements
                    for (IndexType j = 0; j < number_of_nodes; j++) {
                        if (r_N(point_number, j) > shape_function_tolerance ) { // || H_gradient_term(0, j) > shape_function_tolerance || H_hessian_term(0, j) > shape_function_tolerance || H_3rdTayor_term(0, j) > shape_function_tolerance) {
                            // const double NN = r_N(point_number, j) * r_N(point_number, i) * integration;
                            const double NN = (r_N(point_number, j) + H_gradient_term(0, j) + H_hessian_term(0, j) + H_3rdTayor_term(0, j)) * (r_N(point_number, i) + H_gradient_term(0, i) + H_hessian_term(0, i) + H_3rdTayor_term(0, i)) * integration;

                            // indices in local stiffness matrix.
                            const IndexType ibase = counter_n * 1 + 1 * (number_of_non_zero_nodes);
                            const IndexType jbase = counter_m * 1;

                            // Matrix in following shape:
                            // |0 H^T|
                            // |H 0  |

                            LHS(ibase, jbase) = NN;

                            LHS(jbase, ibase) = NN;

                            counter_m++;
                        }
                    }
                    counter_n++;
                }
            }

            if (CalculateStiffnessMatrixFlag) {
                noalias(rLeftHandSideMatrix) += LHS;
            }

            if (CalculateResidualVectorFlag) {
                // const double& displacement = (Has(TEMPERATURE))
                //     ? this->GetValue(TEMPERATURE)
                //     : 0.0;

                // const double displacement = sin(GP_parameter_coord[0]) * sinh(GP_parameter_coord[1]);
                const double displacement = sin(projection[0]) * sinh(projection[1]) ;
                // const double displacement = GP_parameter_coord[0] - GP_parameter_coord[1];
                // const double displacement = projection[0] - projection[1] ;
                // KRATOS_WATCH(displacement)
                // KRATOS_WATCH(sin(GP_parameter_coord[0]) * sinh(GP_parameter_coord[1]))

                IndexType counter = 0;
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    for (IndexType n = 0; n < r_N.size1(); ++n) {
                        if (r_N(n, i) > shape_function_tolerance) {
                            const double r_disp = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
                            
                            u[counter] = (r_disp - displacement);
                            counter++;
                        }
                    }
                }
                for (IndexType i = 0; i < number_of_nodes; i++) {
                    for (IndexType n = 0; n < r_N.size1(); ++n) {
                        if (r_N(n, i) > shape_function_tolerance) {
                            const double r_l_m = r_geometry[i].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
                            u[counter] = r_l_m;
                            counter++;
                        }
                    }
                }

                noalias(rRightHandSideVector) -= prod(LHS, u);
            }
        }
        KRATOS_CATCH("")
    }

    void SBMSupportLagrangeCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType number_of_non_zero_nodes = GetNumberOfNonZeroNodes();

        if (rResult.size() != 2 * number_of_non_zero_nodes)
            rResult.resize(2 * number_of_non_zero_nodes, false);

        IndexType counter = 0;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rResult[counter]     = r_node.GetDof(TEMPERATURE).EquationId();
                    counter++;
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rResult[counter]     = r_node.GetDof(SCALAR_LAGRANGE_MULTIPLIER).EquationId();
                    counter++;
                }
            }
        }
    }

    void SBMSupportLagrangeCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * GetNumberOfNonZeroNodes());

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry[i];
                    rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
                }
            }
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType n = 0; n < r_N.size1(); ++n) {
                if (r_N(n, i) > shape_function_tolerance) {
                    const auto& r_node = r_geometry.GetPoint(i);
                    rElementalDofList.push_back(r_node.pGetDof(SCALAR_LAGRANGE_MULTIPLIER));
                }
            }
        }
    }

    std::size_t SBMSupportLagrangeCondition::GetNumberOfNonZeroNodes() const {
        const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

        SizeType counter = 0;
        for (IndexType n = 0; n < r_N.size1(); ++n) {
            for (IndexType m = 0; m < r_N.size2(); ++m) {
                if (r_N(n, m) > shape_function_tolerance) {
                    counter++;
                }
            }
        }
        return counter;
    }

} // Namespace Kratos
