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

// Application includes
#include "custom_conditions/sbm_fluid_condition.h"


namespace Kratos
{
void SBMFluidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    const unsigned int dim = DN_De[0].size2();
    const SizeType mat_size = number_of_nodes * (dim+1);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
        
    std::string boundaryConditionTypeStr = this->GetValue(BOUNDARY_CONDITION_TYPE);
    BoundaryConditionType boundaryConditionType = GetBoundaryConditionType(boundaryConditionTypeStr);

    if (boundaryConditionType == BoundaryConditionType::Dirichlet)
    {
        this-> CalculateAllDirichlet(rLeftHandSideMatrix,
                                     rRightHandSideVector,
                                     rCurrentProcessInfo);
    } else if (boundaryConditionType == BoundaryConditionType::Neumann) 
    {
        this-> CalculateAllNeumann(rLeftHandSideMatrix,
                                   rRightHandSideVector,
                                   rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "error in SBM_LAPLACIAN_CONDITION, no BOUNDARY_CONDITION_TYPE available" << std::endl;
    }
}

SBMFluidCondition::BoundaryConditionType SBMFluidCondition::GetBoundaryConditionType(const std::string& rType)
{
    KRATOS_ERROR_IF_NOT(rType == "dirichlet" || rType == "neumann") << "Invalid boundary condition type."  << std::endl;
    if (rType == "dirichlet") {
        return BoundaryConditionType::Dirichlet;
    } else {    
        return BoundaryConditionType::Neumann;
    }
}

void SBMFluidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    const unsigned int dim = DN_De[0].size2();
    const SizeType mat_size = number_of_nodes * (dim+1);

    VectorType right_hand_side_vector;

    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2())
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLocalSystem(rLeftHandSideMatrix, right_hand_side_vector,
        rCurrentProcessInfo);
}

void SBMFluidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    const unsigned int dim = DN_De[0].size2();
    const SizeType mat_size = number_of_nodes * (dim+1);

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    MatrixType left_hand_side_matrix = ZeroMatrix(mat_size, mat_size);

    CalculateLocalSystem(left_hand_side_matrix, rRightHandSideVector,
        rCurrentProcessInfo);
}


//_________________________________________________________________________________________________________________________________________
// DIRICHLET CONDITION
//_________________________________________________________________________________________________________________________________________

void SBMFluidCondition::CalculateAllDirichlet(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    Condition candidate_closest_skin_segment_1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;

    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();

    double penalty = GetProperties()[PENALTY_FACTOR];

    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    // Initialize DN_DX
    const unsigned int dim = r_DN_De[0].size2();
    Matrix DN_DX(number_of_nodes,dim);
    noalias(DN_DX) = r_DN_De[0]; 

    Vector meshSize_uv = this->GetValue(MARKER_MESHES);
    double h = std::min(meshSize_uv[0], meshSize_uv[1]);
    if (dim == 3) {h = std::min(h,  meshSize_uv[2]);}

    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (dim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // constitutive law
    Matrix B = ZeroMatrix(3,number_of_nodes*dim);
    CalculateB(B, DN_DX);

    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    ConstitutiveVariables this_constitutive_variables(strain_size);
    Vector old_displacement(number_of_nodes*dim);
    GetValuesVector(old_displacement);

    Vector old_strain = prod(B,old_displacement);
    // Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStrainVector(old_strain);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    const Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix sigmaVoigt = Matrix(prod(r_D, B));

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
    penalty = mBasisFunctionsOrder * mBasisFunctionsOrder * penalty / h;

    // Find the closest node in condition
    int closestNodeId;
    if (dim > 2) {
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
    Vector projection(3);
    projection = candidate_closest_skin_segment_1.GetGeometry()[closestNodeId].Coordinates() ;

    // Print on external file the projection coordinates (projection[0],projection[1]) -> For PostProcess
    std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
    outputFile << projection[0] << " " << projection[1] << " " << projection[2] << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() <<"\n";
    outputFile.close();

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = projection - r_geometry.Center().Coordinates();

    const Matrix& N = r_geometry.ShapeFunctionsValues();

    // Differential area
    double penalty_integration = penalty * r_integration_points[0].Weight() ; // * std::abs(DetJ0);

    // Compute the normals
    array_1d<double, 3> normal_physical_space;

    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    normal_physical_space = mNormalParameterSpace; // prod(trans(J0[0]),mNormalParameterSpace);

    Matrix H = ZeroMatrix(1, number_of_nodes);
    Matrix DN_dot_n = ZeroMatrix(1, number_of_nodes);

    // Compute all the derivatives of the basis functions involved
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        mShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
    }

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        H(0, i) = N(0, i);
        for (IndexType idim = 0; idim < dim; idim++) {
                DN_dot_n(0, i)   += DN_DX(i, idim) * normal_physical_space[idim];
        } 
        // Reset for each node
        double H_taylor_term = 0.0; 

        if (dim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& r_shape_function_derivatives = mShapeFunctionDerivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    double derivative = r_shape_function_derivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term += computeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = mShapeFunctionDerivatives[n-1];
                
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        
                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i,countDerivativeId); 

                        H_taylor_term += computeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        
        mHsum(0,i) = H_taylor_term + H(0,i);
    }

    Vector n_tensor(2);
    n_tensor(0) = mNormalParameterSpace(0); // Component in x direction
    n_tensor(1) = mNormalParameterSpace(1); // Component in y direction
    // Assembly
    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType j = 0; j < number_of_nodes; j++) {
            for (IndexType idim = 0; idim < 2; idim++) {
                
                // Penalty term for the velocity
                rLeftHandSideMatrix(3*i+idim, 3*j+idim) += mHsum(0,i)*mHsum(0,j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {

                    // Extract the 2x2 block for the control point i from the sigma matrix.
                    Matrix sigma_block = ZeroMatrix(2, 2);
                    sigma_block(0, 0) = sigmaVoigt(0, 2*j+jdim);      // sigma(4 * j + 2*jdim, 0);     // sigma_xx for control point i.
                    sigma_block(0, 1) = sigmaVoigt(2, 2*j+jdim);      // sigma(4 * j + 2*jdim, 1);     // sigma_xy for control point i.
                    sigma_block(1, 0) = sigmaVoigt(2, 2*j+jdim);      // sigma(4 * j + 2*jdim + 1, 0); // sigma_yx (symmetric) for control point i.
                    sigma_block(1, 1) = sigmaVoigt(1, 2*j+jdim);      // sigma(4 * j + 2*jdim + 1, 1); // sigma_yy for control point i.
                    
                    // Compute the traction vector: sigma * n.
                    Vector traction = prod(sigma_block, n_tensor); // This results in a 2x1 vector.

                    // integration by parts velocity --> With Constitutive law
                    rLeftHandSideMatrix(3*i+idim, 3*j+jdim) -= H(0, i) * traction(idim) * r_integration_points[0].Weight();

                    // Nitsche term --> With Constitutive law
                    // rLeftHandSideMatrix(3*i+idim, 3*j+jdim) -= mHsum(0, j) * traction(jdim) * r_integration_points[0].Weight();
                }
                
                // integration by parts PRESSURE
                rLeftHandSideMatrix(3*i+idim, 3*j+2) += H(0,j)* ( H(0,i) * mNormalParameterSpace[idim] )
                        * r_integration_points[0].Weight();
            }
        }
    }

    const Vector u_D = candidate_closest_skin_segment_1.GetGeometry()[closestNodeId].GetValue(VELOCITY);

    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {       

            // Penalty term for the velocity
            rRightHandSideVector[3*i+idim] += mHsum(0,i)*u_D[idim]* penalty_integration;

            // Extract the 2x2 block for the control point i from the sigma matrix.
            Matrix sigma_block = ZeroMatrix(2, 2);
            sigma_block(0, 0) = sigmaVoigt(0, 2*i+i);         // sigma(4 * j + 2*jdim, 0);     // sigma_xx for control point i.
            sigma_block(0, 1) = sigmaVoigt(2, 2*i+idim);      // sigma(4 * j + 2*jdim, 1);     // sigma_xy for control point i.
            sigma_block(1, 0) = sigmaVoigt(2, 2*i+idim);      // sigma(4 * j + 2*jdim + 1, 0); // sigma_yx (symmetric) for control point i.
            sigma_block(1, 1) = sigmaVoigt(1, 2*i+idim);      // sigma(4 * j + 2*jdim + 1, 1); // sigma_yy for control point i.
            // Compute the traction vector: sigma * n.
            Vector traction = prod(sigma_block, n_tensor); // This results in a 2x1 vector.
            // Nitsche term --> With Constitutive law
            // rRightHandSideVector[3*i+idim] -= u_D[idim] * traction(idim) * r_integration_points[0].Weight();
        }
    }

    Vector temp = ZeroVector(number_of_nodes);
    GetValuesVector(temp);
    // RHS = ExtForces - K*temp;
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
    
    KRATOS_CATCH("")
}

//_________________________________________________________________________________________________________________________________________
// NEUMANN CONDITION
//_________________________________________________________________________________________________________________________________________

void SBMFluidCondition::CalculateAllNeumann(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    Condition candidate_closest_skin_segment_1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;

    // loopIdentifier is inner or outer
    std::string loopIdentifier = this->GetValue(IDENTIFIER);

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
    
    // Initialize DN_DX
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const unsigned int dim = r_DN_De[0].size2();
    Matrix DN_DX(number_of_nodes,dim);

    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (dim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();

    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number)
    {           
        // Find the closest node in condition
        int closestNodeId;
        if (dim > 2) {
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

        // Obtaining the projection from the closest skin segment
        Vector projection(3);
        projection = candidate_closest_skin_segment_1.GetGeometry()[closestNodeId].Coordinates() ;

        Vector mDistanceVector(3);
        noalias(mDistanceVector) = projection - r_geometry.Center().Coordinates();

        const Matrix& N = r_geometry.ShapeFunctionsValues();

        noalias(DN_DX) = r_DN_De[0]; // prod(r_DN_De[point_number],InvJ0);

        // Compute the normal
        mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
        mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);

        // Neumann BC, the true normal is needed
        array_1d<double, 3> true_n;
        if (dim == 2) {
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
            
            true_n = crossProductSkinSegment1 / MathUtils<double>::Norm(crossProductSkinSegment1) + crossProductSkinSegment2 / MathUtils<double>::Norm(crossProductSkinSegment2);
            if (loopIdentifier == "inner") {
                true_n = true_n / MathUtils<double>::Norm(true_n) ;
            } else { // outer
                true_n = - true_n / MathUtils<double>::Norm(true_n) ;
            }
        } else {
            // 3D CASE
            array_1d<double,3> vector_skin_segment_1 = candidate_closest_skin_segment_1.GetGeometry()[1] - candidate_closest_skin_segment_1.GetGeometry()[0];
            array_1d<double,3> vector_skin_segment_2 = candidate_closest_skin_segment_1.GetGeometry()[2] - candidate_closest_skin_segment_1.GetGeometry()[1];
            MathUtils<double>::CrossProduct(true_n, vector_skin_segment_1, vector_skin_segment_2);

            if (loopIdentifier == "inner") {
                true_n = true_n / MathUtils<double>::Norm(true_n) ;
            } else { // outer
                true_n = - true_n / MathUtils<double>::Norm(true_n) ;
            }
        }

        // Compute all the derivatives of the basis functions involved
        std::vector<Matrix> nShapeFunctionDerivatives;
        for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
            nShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, point_number, this->GetIntegrationMethod()));
        }

        // Neumann (Taylor expansion of the gradient)
        Matrix H = ZeroMatrix(1, number_of_nodes);
        Matrix HgradX = ZeroMatrix(1, number_of_nodes);
        Matrix HgradY = ZeroMatrix(1, number_of_nodes);
        Matrix HgradZ = ZeroMatrix(1, number_of_nodes);

        Matrix DN_dot_n_tilde = ZeroMatrix(1, number_of_nodes);
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            H(0, i) = N(point_number, i);
            // grad N cdot n_tilde
            for (unsigned int idim = 0; idim < dim; idim++) {
                DN_dot_n_tilde(0, i)  += DN_DX(i, idim) * mNormalParameterSpace[idim];
            } 
            // Reset for each control point
            double H_taylor_term_X = 0.0; 
            double H_taylor_term_Y = 0.0; 
            double H_taylor_term_Z = 0.0; 

            if (dim == 2) {
                for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                    // Retrieve the appropriate derivative for the term
                    Matrix& r_shape_function_derivatives = nShapeFunctionDerivatives[n-1];
                    for (IndexType k = 0; k <= n-1; k++) {
                        IndexType n_k = n - 1 - k;
                        double derivative = r_shape_function_derivatives(i,k); 
                        // Compute the Taylor term for this derivative
                        H_taylor_term_X += computeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                    }
                    for (IndexType k = 0; k <= n-1; k++) {
                        IndexType n_k = n - 1 - k;
                        double derivative = r_shape_function_derivatives(i,k+1); 
                        // Compute the Taylor term for this derivative
                        H_taylor_term_Y += computeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                    }
                }
            } else {
                // 3D
                for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                    Matrix& r_shape_function_derivatives = nShapeFunctionDerivatives[n-1];
                
                    int countDerivativeId = 0;
                    // Loop over blocks of derivatives in x
                    for (IndexType k_x = n; k_x >= 0; k_x--) {
                        // Loop over the possible derivatives in y
                        for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
    
                            // derivatives in z
                            IndexType k_z = n - k_x - k_y;
                            double derivative = r_shape_function_derivatives(i,countDerivativeId); 
                            
                            if (k_x >= 1) {
                                H_taylor_term_X += computeTaylorTerm3D(derivative, mDistanceVector[0], k_x-1, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                            }
                            if (k_y >= 1) {
                                H_taylor_term_Y += computeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y-1, mDistanceVector[2], k_z);
                            }
                            if (k_z >= 1) {
                                H_taylor_term_Z += computeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z-1);
                            }     
                            countDerivativeId++;
                        }
                    }
                }
            }
            
            HgradX(0,i) = H_taylor_term_X + DN_DX(i, 0);
            HgradY(0,i) = H_taylor_term_Y + DN_DX(i, 1);
            if (dim == 3) {
                HgradZ(0,i) = H_taylor_term_Z + DN_DX(i, 2);
            }
        }    

        // dot product n cdot n_tilde
        double n_ntilde = true_n[0] * mNormalParameterSpace[0] + true_n[1] * mNormalParameterSpace[1] + true_n[2] * mNormalParameterSpace[2];

        // dot product grad cdot n
        Matrix HgradNdot_n = ZeroMatrix(1, number_of_nodes);
        HgradNdot_n = HgradX * true_n[0] + HgradY * true_n[1] + HgradZ * true_n[2];

        // compute Neumann contributions
        noalias(rLeftHandSideMatrix) += prod(trans(H), HgradNdot_n)  * n_ntilde   * r_integration_points[point_number].Weight(); // * std::abs(determinant_jacobian_vector[point_number]) ;
        noalias(rLeftHandSideMatrix) -= prod(trans(H), DN_dot_n_tilde)            * r_integration_points[point_number].Weight() ; // * std::abs(DetJ0) ;
                            
        Vector t_N(number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            t_N[i] = candidate_closest_skin_segment_1.GetGeometry()[closestNodeId].GetValue(HEAT_FLUX);
        }
        // Neumann Contributions
        noalias(rRightHandSideVector) += prod(prod(trans(H), H), t_N) * n_ntilde * r_integration_points[point_number].Weight(); // * std::abs(determinant_jacobian_vector[point_number]);

        // Vector temp(number_of_nodes);
        // // RHS = ExtForces - K*temp;
        // for (IndexType i = 0; i < number_of_nodes; i++) {
        //     temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
        // }
        // RHS -= K*temp
        // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
    }

    KRATOS_CATCH("")
}

void SBMFluidCondition::CalculateB(
        Matrix& rB, 
        const ShapeDerivativesType& r_DN_DX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (rB.size1() != 3 || rB.size2() != mat_size)
        rB.resize(3, mat_size);

    noalias(rB) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // x-derivatives of shape functions -> relates to strain component ε_11 (xx component)
        rB(0, 2 * i)     = r_DN_DX(i, 0); // ∂N_i / ∂x
        // y-derivatives of shape functions -> relates to strain component ε_22 (yy component)
        rB(1, 2 * i + 1) = r_DN_DX(i, 1); // ∂N_i / ∂y
        // Symmetric shear strain component ε_12 (xy component)
        rB(2, 2 * i)     = r_DN_DX(i, 1); // ∂N_i / ∂y
        rB(2, 2 * i + 1) = r_DN_DX(i, 0); // ∂N_i / ∂x
    }
}

void SBMFluidCondition:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
}

void SBMFluidCondition::InitializeMaterial()
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




unsigned long long SBMFluidCondition::factorial(IndexType n) 
{
    if (n == 0) return 1;
    unsigned long long result = 1;
    for (IndexType i = 2; i <= n; ++i) result *= i;
    return result;
}

// Function to compute a single term in the Taylor expansion
double SBMFluidCondition::computeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (factorial(k) * factorial(n_k));    
}

double SBMFluidCondition::computeTaylorTerm3D(double derivative, double dx, IndexType k_x, double dy, IndexType k_y, double dz, IndexType k_z)
{   
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (factorial(k_x) * factorial(k_y) * factorial(k_z));    
}



int SBMFluidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << "No penalty factor (PENALTY_FACTOR) defined in property of SBMFluidCondition" << std::endl;
    return 0;
}

void SBMFluidCondition::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const int dim = 2;
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (dim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void SBMFluidCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(3 * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SBMFluidCondition::GetValuesVector(
    Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        IndexType index = i * 3;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X);
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y);
        rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue(PRESSURE);
    }
}

} // Namespace Kratos
