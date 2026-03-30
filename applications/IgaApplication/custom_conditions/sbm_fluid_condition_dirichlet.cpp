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

    const std::size_t block_size = mDim + 1;
    const std::size_t mat_size = number_of_nodes * block_size;
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
    const SizeType strain_size = (mDim == 3) ? 6 : 3;

    // Compute the B matrix
    Matrix B = ZeroMatrix(strain_size, number_of_nodes * mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(strain_size);
    ApplyConstitutiveLaw(B, Values, constitutive_variables);
    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix DB_voigt = Matrix(prod(r_D, B));

    const Matrix& H = r_geometry.ShapeFunctionsValues();
    
    // Use geometry Jacobian determinant to get the correct measure (length/area).
    const double det_J0 = std::abs(r_geometry.DeterminantOfJacobian(0, r_geometry.GetDefaultIntegrationMethod()));

    const double penalty_integration = mPenalty * r_integration_points[0].Weight() * std::abs(det_J0);
    const double integration_weight = r_integration_points[0].Weight() * std::abs(det_J0);

    // Compute the pressure & velocity at the previous iteration
    double pressure_old_iteration = 0.0;
    Vector velocity_old_iteration = ZeroVector(mDim);
    for(unsigned int j = 0; j < number_of_nodes; ++j) {
        pressure_old_iteration    += r_geometry[j].GetSolutionStepValue(PRESSURE) * H(0,j);
        const auto& r_velocity = r_geometry[j].GetSolutionStepValue(VELOCITY);
        for (IndexType d = 0; d < mDim; ++d) {
            velocity_old_iteration[d] += r_velocity[d] * mHsum(0,j);
        }
    }

    Vector n_tensor(mDim);
    for (IndexType d = 0; d < mDim; ++d) {
        n_tensor(d) = mNormalParameterSpace(d);
    }

    // Compute the traction vector: sigma * n using r_stress_vector
    Matrix stress_old = ZeroMatrix(mDim, mDim);
    if (mDim == 2) {
        stress_old(0, 0) = r_stress_vector[0];
        stress_old(1, 1) = r_stress_vector[1];
        stress_old(0, 1) = r_stress_vector[2];
        stress_old(1, 0) = r_stress_vector[2];
    } else {
        // 3D Voigt order: xx, yy, zz, xy, yz, xz.
        stress_old(0, 0) = r_stress_vector[0];
        stress_old(1, 1) = r_stress_vector[1];
        stress_old(2, 2) = r_stress_vector[2];
        stress_old(0, 1) = r_stress_vector[3];
        stress_old(1, 0) = r_stress_vector[3];
        stress_old(1, 2) = r_stress_vector[4];
        stress_old(2, 1) = r_stress_vector[4];
        stress_old(0, 2) = r_stress_vector[5];
        stress_old(2, 0) = r_stress_vector[5];
    }
    Vector traction_old_iteration = prod(stress_old, n_tensor); // This results in a 2x1 vector.

    Matrix DB_contribution_w = ZeroMatrix(mDim, mDim);
    Matrix DB_contribution = ZeroMatrix(mDim, mDim);

    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType idim = 0; idim < mDim; idim++) {
            const IndexType col_w = i * mDim + idim;
            BuildStressFromVoigtColumn(DB_contribution_w, DB_voigt, col_w);

            for (IndexType j = 0; j < number_of_nodes; j++) {
                // Compute the traction vector: sigma * n.
                Vector traction_nitsche_w = prod(DB_contribution_w, n_tensor);

                // Penalty term
                rLeftHandSideMatrix(i * block_size + idim, j * block_size + idim) +=
                    mHsum(0,i) * mHsum(0,j) * penalty_integration;
                
                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const IndexType col = j * mDim + jdim;
                    BuildStressFromVoigtColumn(DB_contribution, DB_voigt, col);
                    Vector traction = prod(DB_contribution, n_tensor);

                    // integration by parts velocity < v cdot (DB cdot n) >
                    rLeftHandSideMatrix(i * block_size + idim, j * block_size + jdim) -=
                        H(0, i) * traction(idim) * integration_weight;
                    
                    // skew-symmetric Nitsche term
                    rLeftHandSideMatrix(i * block_size + idim, j * block_size + jdim) +=
                        mHsum(0, j) * traction_nitsche_w(jdim) * integration_weight;
                }

                // integration by parts PRESSURE
                rLeftHandSideMatrix(i * block_size + idim, j * block_size + mDim) +=
                    H(0,j) * ( H(0,i) * mNormalParameterSpace[idim] )
                        * integration_weight;
                
                // Nitsche term --> q term
                rLeftHandSideMatrix(j * block_size + mDim, i * block_size + idim) -=
                    H(0,j) * ( mHsum(0,i) * mNormalParameterSpace[idim] )
                        * integration_weight;
            }
        }

        // --- RHS corresponding terms ---
        for (IndexType idim = 0; idim < mDim; idim++) {
            // Penalty term for the velocity
            rRightHandSideVector(i * block_size + idim) -=
                mHsum(0,i) * velocity_old_iteration[idim] * penalty_integration;
            // integration by parts velocity
            rRightHandSideVector(i * block_size + idim) +=
                H(0,i) * traction_old_iteration(idim) * integration_weight;
            // integration by parts PRESSURE
            rRightHandSideVector(i * block_size + idim) -=
                pressure_old_iteration * ( H(0,i) * mNormalParameterSpace[idim] ) * integration_weight;
        
            // skew-symmetric Nitsche term
            BuildStressFromVoigtColumn(DB_contribution, DB_voigt, i * mDim + idim);
            Vector traction = prod(DB_contribution, n_tensor);
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector(i * block_size + idim) -=
                    velocity_old_iteration[jdim] * traction(jdim) * integration_weight;
            }
            // Nitsche term --> q term
            rRightHandSideVector(i * block_size + mDim) +=
                velocity_old_iteration[idim] * ( H(0,i) * mNormalParameterSpace[idim] )
                        * integration_weight;
        }
    }
    
    // Get the projection node velocity
    const Vector u_D = mpProjectionNode->GetValue(VELOCITY);

    for (IndexType i = 0; i < number_of_nodes; i++) {

        for (IndexType idim = 0; idim < mDim; idim++) {
            
            // Penalty term for the velocity
            rRightHandSideVector[i * block_size + idim] +=
                mHsum(0,i) * u_D[idim] * penalty_integration;

            // Extract the 2x2 block for the control point i from the sigma matrix.
            Matrix sigma_block = ZeroMatrix(mDim, mDim);
            BuildStressFromVoigtColumn(sigma_block, DB_voigt, i * mDim + idim);
            Vector traction = prod(sigma_block, n_tensor);
            // skew-symmetric Nitsche term
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector[i * block_size + idim] +=
                    u_D[jdim] * traction(jdim) * integration_weight;
            }
            // Nitsche term --> q term
            rRightHandSideVector[i * block_size + mDim] -=
                u_D[idim] * H(0,i) * mNormalParameterSpace[idim] * integration_weight;

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
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    if (mBasisFunctionsOrder != 2)
        mBasisFunctionsOrder *= 2;

    mPenalty = GetProperties()[PENALTY_FACTOR];

    // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
    mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * mPenalty / h;

    // Compute the normals
    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
    if (mDim == 3) {
        r_geometry.Calculate(NORMAL, mNormalParameterSpace);
    }
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

    mDistanceVector.resize(mDim);
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
            mHsum(0,i) = H_taylor_term + H(0,i);
        } else {
            // 3D Taylor expansion for the velocity dofs
            for (int n = 1; n <= int(mBasisFunctionsOrder); n++) {
                Matrix& r_shape_function_derivatives = mShapeFunctionDerivatives[n-1];

                int countDerivativeId = 0;
                for (int k_x = n; k_x >= 0; k_x--) {
                    for (int k_y = n - k_x; k_y >= 0; k_y--) {
                        int k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i, countDerivativeId);

                        H_taylor_term += ComputeTaylorTerm3D(
                            derivative,
                            mDistanceVector[0], k_x,
                            mDistanceVector[1], k_y,
                            mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
            mHsum(0,i) = H_taylor_term + H(0,i);
        }
    }
}

void SbmFluidConditionDirichlet::CalculateB(
        Matrix& rB, 
        const ShapeDerivativesType& r_DN_DX) const
{
    const std::size_t number_of_control_points = GetGeometry().size();
    const std::size_t mat_size = number_of_control_points * mDim;
    const std::size_t strain_size = (mDim == 3) ? 6 : 3;

    // Resize B matrix to Voigt strain size and appropriate number of columns.
    if (rB.size1() != strain_size || rB.size2() != mat_size)
        rB.resize(strain_size, mat_size);

    noalias(rB) = ZeroMatrix(strain_size, mat_size);

    if (mDim == 2) {
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            rB(0, 2 * i)     = r_DN_DX(i, 0);
            rB(1, 2 * i + 1) = r_DN_DX(i, 1);
            rB(2, 2 * i)     = r_DN_DX(i, 1);
            rB(2, 2 * i + 1) = r_DN_DX(i, 0);
        }
    } else {
        // 3D small-strain Voigt order: xx, yy, zz, xy, yz, xz.
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            rB(0, 3 * i)     = r_DN_DX(i, 0);
            rB(1, 3 * i + 1) = r_DN_DX(i, 1);
            rB(2, 3 * i + 2) = r_DN_DX(i, 2);
            rB(3, 3 * i)     = r_DN_DX(i, 1);
            rB(3, 3 * i + 1) = r_DN_DX(i, 0);
            rB(4, 3 * i + 1) = r_DN_DX(i, 2);
            rB(4, 3 * i + 2) = r_DN_DX(i, 1);
            rB(5, 3 * i)     = r_DN_DX(i, 2);
            rB(5, 3 * i + 2) = r_DN_DX(i, 0);
        }
    }
}

void SbmFluidConditionDirichlet::BuildStressFromVoigtColumn(
    Matrix& rSigma,
    const Matrix& rDBVoigt,
    const IndexType Column) const
{
    noalias(rSigma) = ZeroMatrix(mDim, mDim);

    if (mDim == 2) {
        rSigma(0, 0) = rDBVoigt(0, Column);
        rSigma(1, 1) = rDBVoigt(1, Column);
        rSigma(0, 1) = rDBVoigt(2, Column);
        rSigma(1, 0) = rDBVoigt(2, Column);
    } else {
        // 3D Voigt order: xx, yy, zz, xy, yz, xz.
        rSigma(0, 0) = rDBVoigt(0, Column);
        rSigma(1, 1) = rDBVoigt(1, Column);
        rSigma(2, 2) = rDBVoigt(2, Column);
        rSigma(0, 1) = rDBVoigt(3, Column);
        rSigma(1, 0) = rDBVoigt(3, Column);
        rSigma(1, 2) = rDBVoigt(4, Column);
        rSigma(2, 1) = rDBVoigt(4, Column);
        rSigma(0, 2) = rDBVoigt(5, Column);
        rSigma(2, 0) = rDBVoigt(5, Column);
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
        if (mDim > 2) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        }
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

        for (IndexType d = 0; d < mDim; ++d) {
            rValues[index + d] = velocity[d];
        }
    }
}

void SbmFluidConditionDirichlet::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    #pragma omp critical
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();

        double pressure_on_surrogate = 0.0;
        array_1d<double, 3> velocity_on_surrogate = ZeroVector(3);
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const double N_i = r_N(0, i);
            pressure_on_surrogate += N_i * r_geometry[i].GetSolutionStepValue(PRESSURE);

            const array_1d<double, 3>& r_velocity_i = r_geometry[i].GetSolutionStepValue(VELOCITY);
            for (IndexType d = 0; d < mDim; ++d) {
                velocity_on_surrogate[d] += N_i * r_velocity_i[d];
            }
        }

        this->SetValue(VELOCITY_X, velocity_on_surrogate[0]);
        this->SetValue(VELOCITY_Y, velocity_on_surrogate[1]);
        this->SetValue(VELOCITY_Z, velocity_on_surrogate[2]);
        this->SetValue(PRESSURE, pressure_on_surrogate);
    }

    KRATOS_CATCH("")
}

// Function to compute a single term in the Taylor expansion
double SbmFluidConditionDirichlet::ComputeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double SbmFluidConditionDirichlet::ComputeTaylorTerm3D(
    double derivative,
    double dx, int k_x,
    double dy, int k_y,
    double dz, int k_z)
{
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) /
        (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
}


} // Namespace Kratos
