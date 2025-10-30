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
#include "custom_conditions/sbm_fluid_condition_dirichlet.h"

namespace Kratos
{

void SbmFluidConditionDirichlet::Initialize(const ProcessInfo& rCurrentProcessInfo)
{ 
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    InitializeMaterial();
}

void SbmFluidConditionDirichlet::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const std::size_t number_of_nodes = r_geometry.size();

    const std::size_t mat_size = number_of_nodes * (mDim+1);
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    Matrix DN_DX(number_of_nodes,mDim);
    noalias(DN_DX) = DN_De[0];

    // Compute the B matrix
    Matrix B = ZeroMatrix(3,number_of_nodes*mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(3);
    ApplyConstitutiveLaw(B, Values, constitutive_variables);
    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix DB_voigt = Matrix(prod(r_D, B));

    const Matrix& H = r_geometry.ShapeFunctionsValues();
    
    GeometryType::JacobiansType J0;
    Matrix InvJ0(mDim,mDim);
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    Matrix jacobian = ZeroMatrix(3,3);
    double det_J0;
    jacobian(0,0) = J0[0](0,0);
    jacobian(0,1) = J0[0](0,1);
    jacobian(1,0) = J0[0](1,0);
    jacobian(1,1) = J0[0](1,1);
    jacobian(2,2) = 1.0;
    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(jacobian,InvJ0,det_J0);
    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space);
    Vector add_factor = prod(jacobian, tangent_parameter_space);
    add_factor[2] = 0.0;
    det_J0 = norm_2(add_factor);

    const double penalty_integration = mPenalty * r_integration_points[0].Weight() * std::abs(det_J0);
    const double integration_weight = r_integration_points[0].Weight() * std::abs(det_J0);

    // Compute the pressure & velocity at the previous iteration
    double pressure_old_iteration = 0.0;
    Vector velocity_old_iteration = ZeroVector(2);
    for(unsigned int j = 0; j < number_of_nodes; ++j) {
        pressure_old_iteration    += r_geometry[j].GetSolutionStepValue(PRESSURE) * H(0,j);
        velocity_old_iteration[0] += r_geometry[j].GetSolutionStepValue(VELOCITY_X) * mHsum(0,j);
        velocity_old_iteration[1] += r_geometry[j].GetSolutionStepValue(VELOCITY_Y) * mHsum(0,j);
    }

    Vector n_tensor(2);
    n_tensor(0) = mNormalParameterSpace(0);
    n_tensor(1) = mNormalParameterSpace(1);

    // Compute the traction vector: sigma * n using r_stress_vector
    Matrix stress_old = ZeroMatrix(2, 2);
    stress_old(0, 0) = r_stress_vector[0];      stress_old(0, 1) = r_stress_vector[2];      
    stress_old(1, 0) = r_stress_vector[2];      stress_old(1, 1) = r_stress_vector[1];         
    Vector traction_old_iteration = prod(stress_old, n_tensor); // This results in a 2x1 vector.

    Matrix DB_contribution_w = ZeroMatrix(2, 2);
    Matrix DB_contribution = ZeroMatrix(2, 2);

    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            DB_contribution_w(0, 0) = DB_voigt(0, 2*i+idim);
            DB_contribution_w(0, 1) = DB_voigt(2, 2*i+idim);
            DB_contribution_w(1, 0) = DB_voigt(2, 2*i+idim);
            DB_contribution_w(1, 1) = DB_voigt(1, 2*i+idim);

            for (IndexType j = 0; j < number_of_nodes; j++) {
                // Compute the traction vector: sigma * n.
                Vector traction_nitsche_w = prod(DB_contribution_w, n_tensor);

                // Penalty term
                rLeftHandSideMatrix(3*i+idim, 3*j+idim) += mHsum(0,i)*mHsum(0,j)* penalty_integration;
                
                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    // Extract the 2x2 block for the control point i from the sigma matrix.
                    DB_contribution(0, 0) = DB_voigt(0, 2*j+jdim);
                    DB_contribution(0, 1) = DB_voigt(2, 2*j+jdim);
                    DB_contribution(1, 0) = DB_voigt(2, 2*j+jdim);
                    DB_contribution(1, 1) = DB_voigt(1, 2*j+jdim);
                    // Compute the traction vector: sigma * n.
                    Vector traction = prod(DB_contribution, n_tensor);

                    // integration by parts velocity < v cdot (DB cdot n) >
                    rLeftHandSideMatrix(3*i+idim, 3*j+jdim) -= H(0, i) * traction(idim) * integration_weight;
                    
                    // skew-symmetric Nitsche term
                    rLeftHandSideMatrix(3*i+idim, 3*j+jdim) += mHsum(0, j) * traction_nitsche_w(jdim) * integration_weight;
                }

                // integration by parts PRESSURE
                rLeftHandSideMatrix(3*i+idim, 3*j+2) += H(0,j)* ( H(0,i) * mNormalParameterSpace[idim] )
                        * integration_weight;
                
                // Nitsche term --> q term
                rLeftHandSideMatrix(3*j+mDim, 3*i+idim) -= H(0,j)* ( mHsum(0,i) * mNormalParameterSpace[idim] )
                        * integration_weight;
            }
        }

        // --- RHS corresponding terms ---
        for (IndexType idim = 0; idim < 2; idim++) {
            // Penalty term for the velocity
            rRightHandSideVector(3*i+idim) -= mHsum(0,i)* velocity_old_iteration[idim] * penalty_integration;
            // integration by parts velocity
            rRightHandSideVector(3*i+idim) += H(0,i) * traction_old_iteration(idim) * integration_weight;
            // integration by parts PRESSURE
            rRightHandSideVector(3*i+idim) -= pressure_old_iteration * ( H(0,i) * mNormalParameterSpace[idim] ) * integration_weight;
        
            // skew-symmetric Nitsche term
            DB_contribution(0, 0) = DB_voigt(0, 2*i+idim); 
            DB_contribution(0, 1) = DB_voigt(2, 2*i+idim); 
            DB_contribution(1, 0) = DB_voigt(2, 2*i+idim);
            DB_contribution(1, 1) = DB_voigt(1, 2*i+idim); 
            // Compute the traction vector: sigma * n.
            Vector traction = prod(DB_contribution, n_tensor);
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector(3*i+idim) -= velocity_old_iteration[jdim] * traction(jdim) * integration_weight;
            }
            // Nitsche term --> q term
            rRightHandSideVector(3*i+mDim) += velocity_old_iteration[idim] * ( H(0,i) * mNormalParameterSpace[idim] )
                        * integration_weight;
        }
    }
    
    // Get the projection node velocity
    const Vector u_D = mpProjectionNode->GetValue(VELOCITY);

    for (IndexType i = 0; i < number_of_nodes; i++) {

        for (IndexType idim = 0; idim < 2; idim++) {
            
            // Penalty term for the velocity
            rRightHandSideVector[3*i+idim] += mHsum(0,i) * u_D[idim] * penalty_integration;

            // Extract the 2x2 block for the control point i from the sigma matrix.
            Matrix sigma_block = ZeroMatrix(2, 2);

            sigma_block(0, 0) = DB_voigt(0, 2*i+idim);
            sigma_block(0, 1) = DB_voigt(2, 2*i+idim);
            sigma_block(1, 0) = DB_voigt(2, 2*i+idim);
            sigma_block(1, 1) = DB_voigt(1, 2*i+idim);
            // Compute the traction vector: sigma * n.
            Vector traction = prod(sigma_block, n_tensor);
            // skew-symmetric Nitsche term
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector[3*i+idim] += u_D[jdim] * traction(jdim) * integration_weight;
            }
            // Nitsche term --> q term
            rRightHandSideVector[3*i+mDim] -= u_D[idim] * H(0,i)*mNormalParameterSpace[idim] * integration_weight;

        }
    }
    KRATOS_CATCH("")
}

void SbmFluidConditionDirichlet::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim > 2)
        << "Stokes problem is not available in 3D" << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    mPenalty = GetProperties()[PENALTY_FACTOR];
    KRATOS_ERROR_IF(mPenalty == -1.0)
        << "Penalty-free formulation is not available for the Stokes problem" << std::endl;

    // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
    mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * mPenalty / h;

    // Compute the normals
    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;
}

void SbmFluidConditionDirichlet::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const std::size_t number_of_nodes = r_geometry.size();

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

    // Compute all the derivatives of the basis functions involved
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        mShapeFunctionDerivatives.push_back(r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod()));
    }
    const Matrix& H = r_geometry.ShapeFunctionsValues();

    // Compute the Hsum matrix
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        // Reset for each node
        double H_taylor_term = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& r_shape_function_derivatives = mShapeFunctionDerivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    double derivative = r_shape_function_derivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        }
        mHsum(0,i) = H_taylor_term + H(0,i);
    }
}

void SbmFluidConditionDirichlet::CalculateB(
        Matrix& rB, 
        const ShapeDerivativesType& r_DN_DX) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

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

void SbmFluidConditionDirichlet::ApplyConstitutiveLaw(
        const Matrix& rB, 
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables) const
{
    const std::size_t number_of_nodes = GetGeometry().size();

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    Vector old_displacement(number_of_nodes*mDim);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(rB,old_displacement);
    rValues.SetStrainVector(old_strain);
    rValues.SetStressVector(rConstitutiveVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void SbmFluidConditionDirichlet::ApplyConstitutiveLawTrue(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
{
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    
    rValues.SetStrainVector(rStrain);
    rValues.SetStressVector(rConstitutiVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
}

void SbmFluidConditionDirichlet::InitializeMaterial()
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

int SbmFluidConditionDirichlet::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
        << "No penalty factor (PENALTY_FACTOR) defined in property of SbmFluidConditionDirichlet" << std::endl;
    return 0;
}

void SbmFluidConditionDirichlet::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const std::size_t number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (mDim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (mDim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void SbmFluidConditionDirichlet::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const std::size_t number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve((mDim + 1) * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SbmFluidConditionDirichlet::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3>& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * mDim;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
    }
}

// Function to compute a single term in the Taylor expansion
double SbmFluidConditionDirichlet::ComputeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}


void SbmFluidConditionDirichlet::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
#pragma omp critical
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = number_of_nodes * mDim;
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    const Matrix integration_point_list_on_true_boundary  = this->GetValue(INTEGRATION_POINTS);
    const Vector integration_weight_list_on_true_boundary = this->GetValue(INTEGRATION_WEIGHTS);
    const Matrix normals_on_true = this->GetValue(INTEGRATION_NORMALS);
    
    // Create a matrix with the same schema used in other IGA fluid conditions
    // [0] w, [1] fx_tot, [2] fy_tot, [3] fx_visc, [4] fy_visc, [5] fx_pres, [6] fy_pres,
    // [7] nx, [8] ny, [9] x_gp, [10] y_gp
    const std::size_t num_results = integration_weight_list_on_true_boundary.size();
    Matrix integration_results(num_results, 11, 0.0);

    double pressure_max_min = 0.0;

    // if (integration_weight_list_on_true_boundary.size() == 0) {
    //     std::cout << "No integration points on true boundary found for condition with ID " << this->Id() << std::endl;
    // }

    for (int i_gauss = 0; i_gauss < integration_weight_list_on_true_boundary.size(); ++i_gauss) 
    {
        const Vector& gp = row(integration_point_list_on_true_boundary, i_gauss);
        const double weight = integration_weight_list_on_true_boundary[i_gauss];
        Vector d = gp - r_geometry.Center();

        // normals
        array_1d<double,3> true_normal;
        true_normal[0] = normals_on_true(i_gauss,0);
        true_normal[1] = normals_on_true(i_gauss,1);
        true_normal[2] = normals_on_true(i_gauss,2);

        // === 1. Pressure at true boundary ===
        double p_true = 0.0;
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            double p_i = r_geometry[i].GetSolutionStepValue(PRESSURE);

            double H_taylor = 0.0;
            for (int n = 1; n <= mBasisFunctionsOrder; ++n) {
                Matrix& p_derivatives = mShapeFunctionDerivatives[n - 1];
                for (int k = 0; k <= n; ++k) {
                    int n_k = n - k;
                    double deriv = p_derivatives(i, k);
                    H_taylor += ComputeTaylorTerm(deriv, d[0], n_k, d[1], k);
                }
            }
            p_true += (r_N(0, i) + H_taylor) * p_i;
        }

        // Andrea's version of the code
        Matrix grad_H_sum_transposed = ZeroMatrix(3, number_of_nodes);
        ComputeGradientTaylorExpansionContribution(grad_H_sum_transposed);
        Matrix grad_H_sum = trans(grad_H_sum_transposed);
        Matrix B_sum = ZeroMatrix(mDim,mat_size);
        CalculateB(B_sum, grad_H_sum);
        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        // obtain the old stress vector on the true boundary (on the projection node)
        ConstitutiveLaw::Parameters values_true(r_geometry, GetProperties(), rCurrentProcessInfo);
        Vector old_strain_on_true = prod(B_sum,old_displacement_coefficient_vector);

        const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLawTrue(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);
        const Vector& r_stress_vector_on_true = values_true.GetStressVector();

            // // Use the normal at the projection node if available (preferred over segment-based average)
            // array_1d<double, 3> true_normal = ZeroVector(3);
            // std::string loop_identifier = this->GetValue(IDENTIFIER);
            // if (mpProjectionNode != nullptr) {
            //     true_normal = mpProjectionNode->GetValue(NORMAL);
            //     const double nrm = MathUtils<double>::Norm(true_normal);
            //     // Enforce unit normal at projection node
            //     KRATOS_ERROR_IF(std::abs(nrm - 1.0) > 1e-12)
            //         << "Projection node NORMAL is not unit-length. Norm = " << nrm
            //         << ". Node Id: " << mpProjectionNode->Id() << std::endl;
            //     // Flip if inner
            //     if (loop_identifier == "outer") {
            //         true_normal = -true_normal;
            //     }
            // } else {
            //     KRATOS_ERROR << "Projection node not set for condition with ID " << this->Id() << std::endl;
            // }



            // std::string loop_identifier = this->GetValue(IDENTIFIER);
            // // Need also the second closest condition in 2D
            // Condition candidate_closest_skin_segment_1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
            // Condition candidate_closest_skin_segment_2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1] ;
            // array_1d<double,3> vector_skin_segment_1 = candidate_closest_skin_segment_1.GetGeometry()[1] - candidate_closest_skin_segment_1.GetGeometry()[0];
            // array_1d<double,3> vector_skin_segment_2 = candidate_closest_skin_segment_2.GetGeometry()[1] - candidate_closest_skin_segment_2.GetGeometry()[0];
            // array_1d<double,3> vector_out_of_plane = ZeroVector(3);
            // vector_out_of_plane[2] = 1.0;
            
            // array_1d<double,3> crossProductSkinSegment1;
            // array_1d<double,3> crossProductSkinSegment2; 
            // MathUtils<double>::CrossProduct(crossProductSkinSegment1, vector_out_of_plane, vector_skin_segment_1);
            // MathUtils<double>::CrossProduct(crossProductSkinSegment2, vector_out_of_plane, vector_skin_segment_2);
            
            // array_1d<double, 3> true_normal = crossProductSkinSegment1 / MathUtils<double>::Norm(crossProductSkinSegment1) + crossProductSkinSegment2 / MathUtils<double>::Norm(crossProductSkinSegment2);
            // if (loop_identifier == "inner") {
            //     true_normal = -true_normal / MathUtils<double>::Norm(true_normal) ; // TODO: check this
            // } else { // outer
            //     true_normal = true_normal / MathUtils<double>::Norm(true_normal) ; // TODO: check this
            // }
        
        
        // === 4. Compute traction = σ · n ===
        Vector normal_stress_true_old = ZeroVector(3);
        normal_stress_true_old[0] = (r_stress_vector_on_true[0] * true_normal[0] + r_stress_vector_on_true[2] * true_normal[1]);
        normal_stress_true_old[1] = (r_stress_vector_on_true[2] * true_normal[0] + r_stress_vector_on_true[1] * true_normal[1]);
        
        // === 5. Compute traction contributions ===
        array_1d<double,2> t_visc = ZeroVector(2);
        t_visc[0] = normal_stress_true_old[0];
        t_visc[1] = normal_stress_true_old[1];

        array_1d<double,2> t_pres = ZeroVector(2);
        t_pres[0] = -p_true * true_normal[0];
        t_pres[1] = -p_true * true_normal[1];

        array_1d<double,2> t_tot = t_visc + t_pres;

        // === 6. Output in standard 11-column layout ===
        integration_results(i_gauss, 0)  = weight;
        integration_results(i_gauss, 1)  = t_tot[0];
        integration_results(i_gauss, 2)  = t_tot[1];
        integration_results(i_gauss, 3)  = t_visc[0];
        integration_results(i_gauss, 4)  = t_visc[1];
        integration_results(i_gauss, 5)  = t_pres[0];
        integration_results(i_gauss, 6)  = t_pres[1];
        integration_results(i_gauss, 7)  = true_normal[0];
        integration_results(i_gauss, 8)  = true_normal[1];
        integration_results(i_gauss, 9)  = gp[0];
        integration_results(i_gauss, 10) = gp[1];

        if (gp[0]>0.2499968 || gp[0]<0.15001)
        {
            // KRATOS_WATCH("First or Last point on true boundary");
            pressure_max_min = p_true;
        }

    }

    // Save to condition as a matrix
    this->SetValue(RESULTS_ON_TRUE_BOUNDARY, integration_results);
    this->SetValue(PRESSURE, pressure_max_min);
}  
}

    



void SbmFluidConditionDirichlet::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_control_points = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points)
    {
        grad_H_sum.resize(3, number_of_control_points);
    }

    // Neumann (Taylor expansion of the gradient)
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each control point
        double H_taylor_term_X = 0.0; 
        double H_taylor_term_Y = 0.0; 
        double H_taylor_term_Z = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = mShapeFunctionDerivatives[n-1];
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
        }
        grad_H_sum(0,i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1,i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2,i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else 
            grad_H_sum(2,i) = 0;
    }    
}



} // Namespace Kratos
