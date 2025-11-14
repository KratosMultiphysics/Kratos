//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
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
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_conditions/support_laplacian_condition.h"


namespace Kratos
{
void SupportLaplacianCondition::CalculateLocalSystem(
MatrixType& rLeftHandSideMatrix,
VectorType& rRightHandSideVector,
const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetGeometry().size() * 1;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SupportLaplacianCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType mat_size = GetGeometry().size() * 1;

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    
    const double penalty = GetProperties()[PENALTY_FACTOR];

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Integration
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    const unsigned int dim = r_DN_De[0].size2();
    Matrix DN_DX(number_of_nodes,dim);
    
    // Compute the normals
    array_1d<double, 3> normal_physical_space;
    array_1d<double, 3> normal_parameter_space;

    normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
    if(dim == 3) {
        r_geometry.Calculate(NORMAL, normal_parameter_space);
    }
    normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);

    normal_physical_space = normal_parameter_space; // prod(trans(J0[0]),normal_parameter_space);
    
    const Matrix& N = r_geometry.ShapeFunctionsValues();

    // Calculating the cartesian derivatives
    noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[point_number],InvJ0);

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    // Jacobian matrix cause J0 is  3x2 and we need 3x3
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0; // 2D case

    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector determinant_factor = prod(Jacobian, tangent_parameter_space);
    determinant_factor[2] = 0.0; // 2D case
    double det_J0 = norm_2(determinant_factor);

    if (dim == 3) {
        Matrix tangent_matrix;
        r_geometry.Calculate(LOCAL_TANGENT_MATRIX, tangent_matrix);  // 3x2
        
        array_1d<double,3> t1, t2;
        for (std::size_t i = 0; i < 3; ++i) {
            t1[i] = tangent_matrix(i, 0);
            t2[i] = tangent_matrix(i, 1);
        }
        // Cross product of the two tangents
        array_1d<double, 3> det_vector = MathUtils<double>::CrossProduct(t1, t2);
        // Norm gives the surface integration factor
        det_J0 = norm_2(det_vector);
    }
    
    Matrix H = ZeroMatrix(1, number_of_nodes);
    Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        H(0, i)            = N(0, i);
        for (unsigned int idim = 0; idim < dim; idim++) {
            DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];    
        } 
    }
    // Differential area
    double penalty_integration = penalty * integration_points[0].Weight();

    // Collins, Lozinsky & Scovazzi innovation
    double nitsche_penalty = 1.0;  // = 1 -> Penalty approach
                                    // = -1 -> Free-penalty approach
    if (penalty == -1.0) {
        penalty_integration = 0.0;
        nitsche_penalty = -1.0;
    }
    
    // Assembly
    noalias(rLeftHandSideMatrix) -= prod(trans(H), H) * penalty_integration * std::abs(det_J0);
    // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
    noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n) * integration_points[0].Weight() * std::abs(det_J0);

    // Assembly Dirichlet BCs -(GRAD_w* n,u) 
    noalias(rLeftHandSideMatrix) -= nitsche_penalty * prod(trans(DN_dot_n), H) * integration_points[0].Weight() * std::abs(det_J0);
        
}

void SupportLaplacianCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType mat_size = GetGeometry().size() * 1;

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();
    
    const double penalty = GetProperties()[PENALTY_FACTOR];

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Integration
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    const unsigned int dim = r_DN_De[0].size2();
    Matrix DN_DX(number_of_nodes,dim);
    
    // Compute the normals
    array_1d<double, 3> normal_physical_space;
    array_1d<double, 3> normal_parameter_space;

    normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
    if(dim == 3) {
        r_geometry.Calculate(NORMAL, normal_parameter_space);
    }
    normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);
    normal_physical_space = normal_parameter_space;
    
    const Matrix& N = r_geometry.ShapeFunctionsValues();

    // Calculating the cartesian derivatives
    noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[point_number],InvJ0);

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    // Jacobian matrix cause J0 is  3x2 and we need 3x3
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0; // 2D case

    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector determinant_factor = prod(Jacobian, tangent_parameter_space);
    determinant_factor[2] = 0.0; // 2D case
    double det_J0 = norm_2(determinant_factor);

    if (dim == 3) {
        Matrix tangent_matrix;
        r_geometry.Calculate(LOCAL_TANGENT_MATRIX, tangent_matrix);  // 3x2

        array_1d<double,3> t1, t2;
        for (std::size_t i = 0; i < 3; ++i) {
            t1[i] = tangent_matrix(i, 0);
            t2[i] = tangent_matrix(i, 1);
        }
        // Cross product of the two tangents
        array_1d<double, 3> det_vector = MathUtils<double>::CrossProduct(t1, t2);
        // Norm gives the surface integration factor
        det_J0 = norm_2(det_vector);
    }
    
    Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
    Vector H_vector = ZeroVector(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        H_vector(i)        = N(0, i); 
        for (unsigned int idim = 0; idim < dim; idim++) {
            DN_dot_n_vec(i)  += DN_DX(i, idim) * normal_physical_space[idim];
        } 
        
    }
    // Differential area
    double penalty_integration = penalty * integration_points[0].Weight();

    // Collins, Lozinsky & Scovazzi innovation
    double nitsche_penalty = 1.0;  // = 1 -> Penalty approach
                                    // = -1 -> Free-penalty approach
    if (penalty == -1.0) {
        penalty_integration = 0.0;
        nitsche_penalty = -1.0;
    }
        
    const double u_D_scalar = this->GetValue(r_unknown_var);

    noalias(rRightHandSideVector) -=  H_vector * u_D_scalar * penalty_integration * std::abs(det_J0);
    noalias(rRightHandSideVector) -= nitsche_penalty * DN_dot_n_vec * u_D_scalar * integration_points[0].Weight() * std::abs(det_J0);

    Vector temperature_old_iteration(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; i++) {
        temperature_old_iteration[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
    }
    // Corresponding RHS
    noalias(rRightHandSideVector) += H_vector * inner_prod(H_vector, temperature_old_iteration) * penalty_integration * std::abs(det_J0);
    // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
    noalias(rRightHandSideVector) += H_vector * inner_prod(DN_dot_n_vec, temperature_old_iteration) * integration_points[0].Weight() * std::abs(det_J0) ;
    // Assembly Dirichlet BCs -(GRAD_w* n,u) 
    noalias(rRightHandSideVector) += nitsche_penalty * DN_dot_n_vec * inner_prod(H_vector,temperature_old_iteration) * integration_points[0].Weight() * std::abs(det_J0) ;


    // for (unsigned int i = 0; i < number_of_nodes; i++) {
    //     std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
    //     outputFile << r_geometry[i].GetId() << "  " << r_geometry[i].GetDof(TEMPERATURE).EquationId() <<"\n";
    //     outputFile.close();
    // }
}

int SupportLaplacianCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << "No penalty factor (PENALTY_FACTOR) defined in property of SupportLaplacianCondition" << std::endl;
    return 0;
}

void SupportLaplacianCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    if (rResult.size() !=  number_of_nodes)
        rResult.resize(number_of_nodes, false);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto& r_node = r_geometry[i];
        rResult[i] = r_node.GetDof(r_unknown_var).EquationId();
    }
}

void SupportLaplacianCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto& r_node = r_geometry[i];
        rElementalDofList.push_back(r_node.pGetDof(r_unknown_var));
    }
};

void SupportLaplacianCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // array_1d<double, 3> normal_parameter_space = - GetGeometry().Normal(0, GetIntegrationMethod());
    // const GeometryType::ShapeFunctionsGradientsType& r_DN_De = GetGeometry().ShapeFunctionsLocalGradients(GetGeometry().GetDefaultIntegrationMethod());
    // const unsigned int dim = r_DN_De[0].size2();
    // if(dim == 3) {
    //     GetGeometry().Calculate(NORMAL, normal_parameter_space);
    // }
    // SetValue(NORMAL_TO_WALL, normal_parameter_space);
}


} // Namespace Kratos
