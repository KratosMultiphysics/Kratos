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
#include "custom_conditions/sbm_laplacian_condition_neumann.h"


namespace Kratos
{

void SbmLaplacianConditionNeumann::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}

void SbmLaplacianConditionNeumann::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

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

    // Integration
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    Matrix DN_DX(number_of_nodes,mDim);

    // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[0],InvJ0);
    
    const Matrix& H = r_geometry.ShapeFunctionsValues();

    // compute gradient Taylor expansion contribution: grad_H_sum
    Matrix grad_H_sum = ZeroMatrix(3, number_of_nodes);
    ComputeGradientTaylorExpansionContribution(grad_H_sum);
    
    // dot product grad cdot n
    Matrix H_grad_N_dot_n = ZeroMatrix(1, number_of_nodes);
    noalias(row(H_grad_N_dot_n, 0)) = row(grad_H_sum, 0) * mTrueNormal[0] + row(grad_H_sum, 1) * mTrueNormal[1] + row(grad_H_sum, 2) * mTrueNormal[2];

    Matrix DN_dot_n_tilde = ZeroMatrix(1, number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        // grad N cdot n_tilde
        for (IndexType idim = 0; idim < mDim; idim++) {
            DN_dot_n_tilde(0, i)  += DN_DX(i, idim) * mNormalParameterSpace[idim];
        } 
    }
    
    // Assembly
    // compute Neumann contributions
    noalias(rLeftHandSideMatrix) += prod(trans(H), H_grad_N_dot_n) * mTrueDotSurrogateNormal * r_integration_points[0].Weight(); // * std::abs(determinant_jacobian_vector[point_number]) ;
    noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n_tilde)                           * r_integration_points[0].Weight() ; // * std::abs(DetJ0) ;
    
    Vector t_N(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        t_N[i] = mpProjectionNode->GetValue(HEAT_FLUX);
    }
    // Neumann Contributions
    noalias(rRightHandSideVector) += prod(prod(trans(H), H), t_N) * mTrueDotSurrogateNormal * r_integration_points[0].Weight(); // * std::abs(determinant_jacobian_vector[point_number]);

    Vector temp(number_of_nodes);
    // RHS = ExtForces - K*temp;
    for (IndexType i = 0; i < number_of_nodes; i++) {
        temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
    }
    // RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("")

}

void SbmLaplacianConditionNeumann::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    // Compute the normals
    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace; // prod(trans(J0[0]),mNormalParameterSpace);
}

void SbmLaplacianConditionNeumann::InitializeSbmMemberVariables()
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

    // loopIdentifier is inner or outer
    std::string loopIdentifier = this->GetValue(IDENTIFIER);
    if (mDim == 2) {
        // Need also the second closest condition in 2D
        Condition candidate_closest_skin_segment_2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1] ;
        array_1d<double,3> vector_skin_segment_1 = candidate_closest_skin_segment_1.GetGeometry()[1] - candidate_closest_skin_segment_1.GetGeometry()[0];
        array_1d<double,3> vector_skin_segment_2 = candidate_closest_skin_segment_2.GetGeometry()[1] - candidate_closest_skin_segment_2.GetGeometry()[0];
        array_1d<double,3> vector_out_of_plane = ZeroVector(3);
        vector_out_of_plane[2] = 1.0;
        
        array_1d<double,3> crossProductSkinSegment1;
        array_1d<double,3> crossProductSkinSegment2; 
        MathUtils<double>::CrossProduct(crossProductSkinSegment1, vector_out_of_plane, vector_skin_segment_1);
        MathUtils<double>::CrossProduct(crossProductSkinSegment2, vector_out_of_plane, vector_skin_segment_2);
        
        mTrueNormal = crossProductSkinSegment1 / MathUtils<double>::Norm(crossProductSkinSegment1) + crossProductSkinSegment2 / MathUtils<double>::Norm(crossProductSkinSegment2);
        if (loopIdentifier == "inner") {
            mTrueNormal = mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        } else { // outer
            mTrueNormal = - mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        }
    } else {
        // 3D CASE
        array_1d<double,3> vector_skin_segment_1 = candidate_closest_skin_segment_1.GetGeometry()[1] - candidate_closest_skin_segment_1.GetGeometry()[0];
        array_1d<double,3> vector_skin_segment_2 = candidate_closest_skin_segment_1.GetGeometry()[2] - candidate_closest_skin_segment_1.GetGeometry()[1];
        MathUtils<double>::CrossProduct(mTrueNormal, vector_skin_segment_1, vector_skin_segment_2);

        if (loopIdentifier == "inner") {
            mTrueNormal = mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        } else { // outer
            mTrueNormal = - mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        }
    }

    // dot product n dot n_tilde
    mTrueDotSurrogateNormal = inner_prod(mNormalParameterSpace, mTrueNormal);
}

void SbmLaplacianConditionNeumann::CalculateLeftHandSide(
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

    // Calculating the PHYSICAL SPACE derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[0],InvJ0);
    
    const Matrix& H = r_geometry.ShapeFunctionsValues();

    // compute gradient Taylor expansion contribution: grad_H_sum
    Matrix grad_H_sum = ZeroMatrix(3, number_of_nodes);
    ComputeGradientTaylorExpansionContribution(grad_H_sum);
    
    // dot product grad cdot n
    Matrix H_grad_N_dot_n = ZeroMatrix(1, number_of_nodes);
    noalias(row(H_grad_N_dot_n, 0)) = row(grad_H_sum, 0) * mTrueNormal[0] + row(grad_H_sum, 1) * mTrueNormal[1] + row(grad_H_sum, 2) * mTrueNormal[2];

    Matrix DN_dot_n_tilde = ZeroMatrix(1, number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        // grad N cdot n_tilde
        for (IndexType idim = 0; idim < mDim; idim++) {
            DN_dot_n_tilde(0, i)  += DN_DX(i, idim) * mNormalParameterSpace[idim];
        } 
    }
    
    // Assembly
    // compute Neumann contributions
    noalias(rLeftHandSideMatrix) += prod(trans(H), H_grad_N_dot_n) * mTrueDotSurrogateNormal * r_integration_points[0].Weight(); // * std::abs(determinant_jacobian_vector[point_number]) ;
    noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n_tilde)                           * r_integration_points[0].Weight() ; // * std::abs(DetJ0) ;
}


void SbmLaplacianConditionNeumann::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != number_of_nodes) {
        rRightHandSideVector.resize(number_of_nodes, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);

    // Integration
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    
    const Matrix& H = r_geometry.ShapeFunctionsValues();
    
    // Assembly
    Vector t_N(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        t_N[i] = mpProjectionNode->GetValue(HEAT_FLUX);
    }
    // Neumann Contributions
    noalias(rRightHandSideVector) += prod(prod(trans(H), H), t_N) * mTrueDotSurrogateNormal * r_integration_points[0].Weight(); // * std::abs(determinant_jacobian_vector[point_number]);
}

void SbmLaplacianConditionNeumann::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_nodes)
    {
        grad_H_sum.resize(3, number_of_nodes);
    }

    // Neumann (Taylor expansion of the gradient)
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        // Reset for each control point
        double H_taylor_term_X = 0.0; 
        double H_taylor_term_Y = 0.0; 
        double H_taylor_term_Z = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (int n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
            
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (int k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (int k_y = n - k_x; k_y >= 0; k_y--) {

                        // derivatives in z
                        int k_z = n - k_x - k_y;
                        double derivative = shapeFunctionDerivatives(i,countDerivativeId); 
                        
                        if (k_x >= 1) {
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x-1, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y-1, mDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z-1);
                        }     
                        countDerivativeId++;
                    }
                }
            }
        }
        grad_H_sum(0,i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1,i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2,i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else 
            grad_H_sum(2,i) = 0;
    }    
}

// Function to compute a single term in the Taylor expansion
double SbmLaplacianConditionNeumann::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double SbmLaplacianConditionNeumann::ComputeTaylorTerm3D(
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


int SbmLaplacianConditionNeumann::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    return 0;
}

void SbmLaplacianConditionNeumann::EquationIdVector(
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

void SbmLaplacianConditionNeumann::GetDofList(
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
