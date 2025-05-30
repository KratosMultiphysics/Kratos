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
//                   Andrea Gorgi
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "custom_conditions/sbm_laplacian_condition_dirichlet.h"


namespace Kratos
{

void SbmLaplacianConditionDirichlet::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}

void SbmLaplacianConditionDirichlet::CalculateLocalSystem(
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

void SbmLaplacianConditionDirichlet::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();
    
    Vector mesh_size_uv = this->GetValue(MARKER_MESHES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    double penalty = GetProperties()[PENALTY_FACTOR];
    // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
    mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * penalty / h;

    // https://doi.org/10.1016/j.cma.2023.116301 (A penalty-free Shifted Boundary Method of arbitrary order)
    mNitschePenalty = 1.0;   // = 1.0 -> Penalty approach
                                    // = -1.0 -> Free-penalty approach
    if (penalty == -1.0) {
        mPenalty = 0.0;
        mNitschePenalty = -1.0;
    }
    // Compute the normals
    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;
}

void SbmLaplacianConditionDirichlet::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    // Retrieve projection
    Condition candidate_closest_skin_segment_1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
    // Find the closest node in condition
    int closestNodeId;
    if (mDim > 2) {
        double incumbent_dist = 1e16;
        // Loop over the three nodes of the closest skin element
        for (unsigned int i = 0; i < 3; i++) {
            if (norm_2(candidate_closest_skin_segment_1.GetGeometry()[i]-r_geometry.Center()) < incumbent_dist) {
                incumbent_dist = norm_2(candidate_closest_skin_segment_1.GetGeometry()[i]-r_geometry.Center());
                closestNodeId = i;
            }
        }
    } else {
        closestNodeId = 0;
    }
    mpProjectionNode = &candidate_closest_skin_segment_1.GetGeometry()[closestNodeId] ;

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = mpProjectionNode->Coordinates() - r_geometry.Center().Coordinates();
}

void SbmLaplacianConditionDirichlet::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();

    if (rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes) {
        rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Integration
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    Matrix DN_DX(number_of_nodes,mDim);

    // Differential area
    double penalty_integration = mPenalty * r_integration_points[0].Weight() ;

    // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[0],InvJ0);

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    // Jacobian matrix cause J0 is  3x2 and we need 3x3
    Matrix jacobian_matrix = ZeroMatrix(3,3);
    jacobian_matrix(0,0) = J0[0](0,0);
    jacobian_matrix(0,1) = J0[0](0,1);
    jacobian_matrix(1,0) = J0[0](1,0);
    jacobian_matrix(1,1) = J0[0](1,1);
    jacobian_matrix(2,2) = 1.0; // 2D case

    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector determinant_factor = prod(jacobian_matrix, tangent_parameter_space);
    determinant_factor[2] = 0.0; // 2D case
    const double det_J0 = norm_2(determinant_factor);

    Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);
    Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
    
    const Matrix& H = r_geometry.ShapeFunctionsValues();
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        for (IndexType idim = 0; idim < mDim; idim++) {
                DN_dot_n_vec(i)  += DN_DX(i, idim) * mNormalPhysicalSpace[idim];         
        } 
    }
    noalias(row(DN_dot_n, 0)) = DN_dot_n_vec;

    // compute Taylor expansion contribution: H_sum_vec
    Vector H_sum_vec = ZeroVector(number_of_nodes);
    ComputeTaylorExpansionContribution (H_sum_vec);
    
    Matrix H_sum = ZeroMatrix(1, number_of_nodes);
    noalias(row(H_sum, 0)) = H_sum_vec;

    // Assembly
    // -(GRAD_w * n, u + GRAD_u * d + ...)
    noalias(rLeftHandSideMatrix) -= mNitschePenalty * prod(trans(DN_dot_n), H_sum)  * r_integration_points[0].Weight() * std::abs(det_J0) ;
    // -(w,GRAD_u * n) from integration by parts -> Fundamental !! 
    noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n)                        * r_integration_points[0].Weight() * std::abs(det_J0) ;
    // SBM terms (Taylor Expansion) + alpha * (w + GRAD_w * d + ..., u + GRAD_u * d + ...)
    noalias(rLeftHandSideMatrix) += prod(trans(H_sum), H_sum) * penalty_integration * std::abs(det_J0);
}


void SbmLaplacianConditionDirichlet::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{ 
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != number_of_nodes) {
        rRightHandSideVector.resize(number_of_nodes, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);

    // Integration
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    Matrix DN_DX(number_of_nodes,mDim);

    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    // Jacobian matrix cause J0 is  3x2 and we need 3x3
    Matrix jacobian_matrix = ZeroMatrix(3,3);
    jacobian_matrix(0,0) = J0[0](0,0);
    jacobian_matrix(0,1) = J0[0](0,1);
    jacobian_matrix(1,0) = J0[0](1,0);
    jacobian_matrix(1,1) = J0[0](1,1);
    jacobian_matrix(2,2) = 1.0; // 2D case

    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector determinant_factor = prod(jacobian_matrix, tangent_parameter_space);
    determinant_factor[2] = 0.0; // 2D case
    const double det_J0 = norm_2(determinant_factor);

    // Differential area
    double penalty_integration = mPenalty * r_integration_points[0].Weight();

    // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[0],InvJ0);

    const Matrix& H = r_geometry.ShapeFunctionsValues();
    Vector DN_dot_n_vec = ZeroVector(number_of_nodes);
    Vector H_vec = ZeroVector(number_of_nodes);
    
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        H_vec(i) = H(0, i);
        for (IndexType idim = 0; idim < mDim; idim++) {
                DN_dot_n_vec(i)  += DN_DX(i, idim) * mNormalPhysicalSpace[idim];         
        } 
    }

    // compute Taylor expansion contribution: H_sum_vec
    Vector H_sum_vec = ZeroVector(number_of_nodes);
    ComputeTaylorExpansionContribution (H_sum_vec);
    
    Matrix H_sum = ZeroMatrix(1, number_of_nodes);
    noalias(row(H_sum, 0)) = H_sum_vec;

    // Assembly
    const double u_D_scalar = mpProjectionNode->GetValue(r_unknown_var);

    noalias(rRightHandSideVector) += H_sum_vec * u_D_scalar * penalty_integration * std::abs(det_J0);
    // Dirichlet BCs
    noalias(rRightHandSideVector) -= mNitschePenalty * DN_dot_n_vec * u_D_scalar * r_integration_points[0].Weight() * std::abs(det_J0) ;


    Vector temperature_old_iteration(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; i++) {
        temperature_old_iteration[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
    }
    // Corresponding RHS
    noalias(rRightHandSideVector) += mNitschePenalty * DN_dot_n_vec * inner_prod(H_sum_vec,temperature_old_iteration)  * r_integration_points[0].Weight() * std::abs(det_J0);
    // Assembly of the integration by parts term -(w,GRAD_u * n) -> Fundamental !!
    noalias(rRightHandSideVector) += H_vec * inner_prod(DN_dot_n_vec,temperature_old_iteration) * r_integration_points[0].Weight() * std::abs(det_J0) ;
    noalias(rRightHandSideVector) -= H_sum_vec * inner_prod(H_sum_vec,temperature_old_iteration) * penalty_integration * std::abs(det_J0) ;

}

void SbmLaplacianConditionDirichlet::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const Matrix& N = r_geometry.ShapeFunctionsValues();

    if (H_sum_vec.size() != number_of_nodes)
    {
        H_sum_vec = ZeroVector(number_of_nodes);
    }

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        // Reset for each node
        double H_taylor_term = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    double derivative = r_shape_function_derivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        
                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i,countDerivativeId); 

                        H_taylor_term += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + N(0,i);
    }
}

// Function to compute a single term in the Taylor expansion
double SbmLaplacianConditionDirichlet::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double SbmLaplacianConditionDirichlet::ComputeTaylorTerm3D(
    const double derivative, 
    const double dx, 
    const IndexType k_x, 
    const double dy, 
    const IndexType k_y, 
    const double dz, 
    const IndexType k_z)
{   
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));    
}


int SbmLaplacianConditionDirichlet::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << "No penalty factor (PENALTY_FACTOR) defined in property of SbmLaplacianConditionDirichlet" << std::endl;
    return 0;
}

void SbmLaplacianConditionDirichlet::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
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

void SbmLaplacianConditionDirichlet::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
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
}

} // Namespace Kratos
