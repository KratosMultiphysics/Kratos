
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
#include "custom_conditions/sbm_load_solid_condition.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

void SbmLoadSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}


void SbmLoadSolidCondition::InitializeMaterial()
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

void SbmLoadSolidCondition::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "SbmLoadSolidCondition momentarily only supports 2D and 3D conditions, but the current dimension is" << mDim << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    // Compute the normals
    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());

    if (mDim == 3) {
        r_geometry.Calculate(NORMAL, mNormalParameterSpace);
    }

    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);

    mNormalPhysicalSpace = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpace);

    // calculate the integration weight
    // reading integration point
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    if (mDim == 2) {
        // Compute the local tangent
        array_1d<double, 3> tangent_parameter_space;
        r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space

        Matrix InvJ0(mDim,mDim);
        double detJ0;
        Matrix Jacobian;
        CalculateInitialJacobian(r_geometry, Jacobian);
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

        Vector add_factor = ZeroVector(3);
        add_factor[0] = Jacobian(0,0) * tangent_parameter_space[0] + Jacobian(0,1) * tangent_parameter_space[1];
        add_factor[1] = Jacobian(1,0) * tangent_parameter_space[0] + Jacobian(1,1) * tangent_parameter_space[1];
        detJ0 = norm_2(add_factor);
        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
        const double int_to_reference_weight = r_integration_points[0].Weight() * std::abs(detJ0) * thickness;
        SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);
    }
    else
    {
        const double det_J0 = std::abs(r_geometry.DeterminantOfJacobian(0, r_geometry.GetDefaultIntegrationMethod()));
        const double int_to_reference_weight = r_integration_points[0].Weight() * det_J0;

        SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);
    }
}

void SbmLoadSolidCondition::InitializeSbmMemberVariables()
{
    auto& r_geometry = this->GetGeometry();
    std::string loopIdentifier = this->GetValue(IDENTIFIER);

    // NURBS case
    if (this->GetValue(NEIGHBOUR_NODES).size() != 0) 
    {
        mpProjectionNode = &r_geometry.GetValue(NEIGHBOUR_NODES)[0];

        mTrueNormal = mpProjectionNode->GetValue(NORMAL);

        if (loopIdentifier == "inner")
            mTrueNormal = -mTrueNormal;
            
        mDistanceVector.resize(3);
        noalias(mDistanceVector) = mpProjectionNode->Coordinates() - r_geometry.Center().Coordinates();

        // dot product n dot n_tilde
        mTrueDotSurrogateNormal = inner_prod(mNormalParameterSpace, mTrueNormal);
        return;
    }

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

    // calculate the integration weight
    // loopIdentifier is inner or outer
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

void SbmLoadSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetGeometry().size() * mDim;

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

void SbmLoadSolidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const SizeType mat_size = number_of_control_points * mDim;
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points, mDim);

    // Calculating the physical derivatives (it is avoided storing them to minimize storage)
    if (mDim == 2) {
        Matrix InvJ0(mDim, mDim);
        double detJ0;
        Matrix Jacobian;
        CalculateInitialJacobian(r_geometry, Jacobian);
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

        Matrix sub_inv_jacobian = ZeroMatrix(2,2);
        sub_inv_jacobian(0,0) = InvJ0(0,0);
        sub_inv_jacobian(1,0) = InvJ0(1,0);
        sub_inv_jacobian(0,1) = InvJ0(0,1);
        sub_inv_jacobian(1,1) = InvJ0(1,1);
        noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);
    } else {
        noalias(DN_DX) = r_DN_De[0];
    }

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B = ZeroMatrix(strain_size,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);

    const Matrix& r_D = Values.GetConstitutiveMatrix();

    const Matrix DB = prod(r_D,B);

    // compute Taylor expansion contribution: H_sum_vec
    Matrix grad_H_sum_transposed = ZeroMatrix(strain_size, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(mDistanceVector, grad_H_sum_transposed);

    Matrix grad_H_sum = trans(grad_H_sum_transposed);

    Matrix B_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_sum, grad_H_sum);
    Matrix DB_sum = prod(r_D, B_sum); //

    for (IndexType i = 0; i < number_of_control_points; i++) {
        for (IndexType j = 0; j < number_of_control_points; j++) {
            
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim;

                    // FLUX STANDARD TERM
                    Vector sigma_u_n;
                    Vector stress_column_u = column(DB, jglob);
                    CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n);

                    Vector sigma_u_n_sum;
                    Vector stress_column_u_sum = column(DB_sum, jglob);
                    CalculateTraction(stress_column_u_sum, mTrueNormal, sigma_u_n_sum);

                    rLeftHandSideMatrix(iglob, jglob) -= r_N(0,i)*sigma_u_n[idim] * int_to_reference_weight;
                
                    // SBM TERM
                    rLeftHandSideMatrix(iglob, jglob) += r_N(0,i)*sigma_u_n_sum[idim] 
                                                        * int_to_reference_weight * mTrueDotSurrogateNormal;
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void SbmLoadSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const SizeType mat_size = number_of_control_points * mDim;
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points, mDim);

    // Calculating the physical derivatives (it is avoided storing them to minimize storage)
    if (mDim == 2) {
        Matrix InvJ0(mDim, mDim);
        double detJ0;
        Matrix Jacobian;
        CalculateInitialJacobian(r_geometry, Jacobian);
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

        Matrix sub_inv_jacobian = ZeroMatrix(2,2);
        sub_inv_jacobian(0,0) = InvJ0(0,0);
        sub_inv_jacobian(1,0) = InvJ0(1,0);
        sub_inv_jacobian(0,1) = InvJ0(0,1);
        sub_inv_jacobian(1,1) = InvJ0(1,1);
        noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);
    } else {
        noalias(DN_DX) = r_DN_De[0];
    }

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B = ZeroMatrix(strain_size,mat_size);
    CalculateB(B, DN_DX);

    // compute Taylor expansion contribution: H_sum_vec
    Matrix grad_H_sum_transposed = ZeroMatrix(strain_size, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(mDistanceVector, grad_H_sum_transposed);

    Matrix grad_H_sum = trans(grad_H_sum_transposed);

    Matrix B_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_sum, grad_H_sum);

    // Obtain the tangent costitutive law matrix

    ConstitutiveLaw::Parameters values_surrogate(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    const SizeType strain_size_surrogate = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_surrogate(strain_size_surrogate);
    ApplyConstitutiveLaw(mat_size, old_strain, values_surrogate, this_constitutive_variables_surrogate);

    const Matrix& r_D = values_surrogate.GetConstitutiveMatrix();
    const Vector& r_stress_vector = values_surrogate.GetStressVector();
    
    // obtain the old stress vector on the true boundary (on the projection node)
    ConstitutiveLaw::Parameters values_true(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_strain_on_true = prod(B_sum,old_displacement_coefficient_vector);

    const SizeType strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Vector& r_stress_vector_on_true = values_true.GetStressVector();

    // Assembly
    const Matrix DB = prod(r_D,B);
    const Matrix DB_sum = prod(r_D, B_sum); //

    Vector g_N = ZeroVector(3);
    g_N = mpProjectionNode->GetValue(FORCE);

    //FIXME:
    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);


    //2D
    // const double x = mpProjectionNode->X();
    // const double y = mpProjectionNode->Y();

    // // // // cosinusoidal
    // g_N[0] = E/(1-nu)*(sin(x)*sinh(y)) * mTrueNormal[0]; 
    // g_N[1] = E/(1-nu)*(sin(x)*sinh(y)) * mTrueNormal[1]; 


    //3D 
    const double x = mpProjectionNode->X();
    const double y = mpProjectionNode->Y();
    const double z = mpProjectionNode->Z();

    const double c_vol = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    const double c_shear = E / (2.0 * (1.0 + nu));

    const double div_u =
    -std::sin(x) * std::sinh(y) * std::cosh(z)
    -std::sin(y) * std::sinh(z) * std::cosh(x)
    -std::sin(z) * std::sinh(x) * std::cosh(y);

    const double sigma_xx =
        c_vol * div_u
        - 2.0 * c_shear * std::sin(x) * std::sinh(y) * std::cosh(z);

    const double sigma_yy =
        c_vol * div_u
        - 2.0 * c_shear * std::sin(y) * std::sinh(z) * std::cosh(x);

    const double sigma_zz =
        c_vol * div_u
        - 2.0 * c_shear * std::sin(z) * std::sinh(x) * std::cosh(y);

    const double sigma_xy =
        c_shear * (
            std::cos(x) * std::cosh(y) * std::cosh(z)
        + std::cos(y) * std::sinh(z) * std::sinh(x)
        );

    const double sigma_yz =
        c_shear * (
            std::cos(y) * std::cosh(z) * std::cosh(x)
        + std::cos(z) * std::sinh(x) * std::sinh(y)
        );

    const double sigma_xz =
        c_shear * (
            std::cos(x) * std::sinh(y) * std::sinh(z)
        + std::cos(z) * std::cosh(x) * std::cosh(y)
        );

    array_1d<double, 3> traction;

    traction[0] =
        sigma_xx * mTrueNormal[0]
        + sigma_xy * mTrueNormal[1]
        + sigma_xz * mTrueNormal[2];
    
    traction[1] =
        sigma_xy * mTrueNormal[0]
        + sigma_yy * mTrueNormal[1]
        + sigma_yz * mTrueNormal[2];
    
    traction[2] =
        sigma_xz * mTrueNormal[0]
        + sigma_yz * mTrueNormal[1]
        + sigma_zz * mTrueNormal[2];

    g_N = traction; //FIXME::


    Vector normal_stress_old = ZeroVector(3);
    CalculateTraction(r_stress_vector, mNormalPhysicalSpace, normal_stress_old);

    Vector normal_stress_true_old = ZeroVector(3);
    CalculateTraction(r_stress_vector_on_true, mTrueNormal, normal_stress_true_old);
   
    for (IndexType i = 0; i < number_of_control_points; i++) {
        
        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iglob = mDim*i+idim;

            Vector sigma_w_n;
            Vector stress_column_w = column(DB, iglob);
            CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

            // External load term
            rRightHandSideVector[iglob] += r_N(0,i) * g_N[idim] * mTrueDotSurrogateNormal * int_to_reference_weight;

            // Residual terms
            // FLUX STANDARD TERM
            rRightHandSideVector(iglob) += r_N(0,i) * normal_stress_old[idim] * int_to_reference_weight;
        
            // // SBM TERM
            rRightHandSideVector(iglob) -= r_N(0,i)*normal_stress_true_old[idim] * int_to_reference_weight * mTrueDotSurrogateNormal;
        }
    }
    KRATOS_CATCH("")
}

    int SbmLoadSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SbmLoadSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != mDim * number_of_control_points)
            rResult.resize(mDim * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * mDim;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            if (mDim == 3) {
                rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
            }
        }
    }

    void SbmLoadSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(mDim * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            if (mDim == 3) {
                rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
            }
        }
    };


    void SbmLoadSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * mDim;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * mDim;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            if (mDim == 3) {
                rValues[index + 2] = displacement[2];
            }
        }
    }

    void SbmLoadSolidCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType dim = mDim;
        const SizeType strain_size = (dim == 2) ? 3 : 6;
        const SizeType mat_size = number_of_control_points * dim;

        if (rB.size1() != strain_size || rB.size2() != mat_size) {
            rB.resize(strain_size, mat_size, false);
        }

        noalias(rB) = ZeroMatrix(strain_size, mat_size);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * dim;

            if (dim == 2) {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy

                rB(2, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(2, index + 1) = r_DN_DX(i, 0);
            } else {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy
                rB(2, index + 2) = r_DN_DX(i, 2); // ezz

                rB(3, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(3, index + 1) = r_DN_DX(i, 0);

                rB(4, index + 1) = r_DN_DX(i, 2); // gamma_yz
                rB(4, index + 2) = r_DN_DX(i, 1);

                rB(5, index + 0) = r_DN_DX(i, 2); // gamma_xz
                rB(5, index + 2) = r_DN_DX(i, 0);
            }
        }
    }

    void SbmLoadSolidCondition::CalculateTraction(
        const Vector& rStressVector,
        const array_1d<double, 3>& rNormal,
        Vector& rTraction) const
    {
        if (rTraction.size() != 3) {
            rTraction.resize(3, false);
        }

        noalias(rTraction) = ZeroVector(3);

        if (mDim == 2) {
            rTraction[0] = rStressVector[0] * rNormal[0]
                         + rStressVector[2] * rNormal[1];

            rTraction[1] = rStressVector[2] * rNormal[0]
                         + rStressVector[1] * rNormal[1];
        } else {
            rTraction[0] = rStressVector[0] * rNormal[0]
                         + rStressVector[3] * rNormal[1]
                         + rStressVector[5] * rNormal[2];

            rTraction[1] = rStressVector[3] * rNormal[0]
                         + rStressVector[1] * rNormal[1]
                         + rStressVector[4] * rNormal[2];

            rTraction[2] = rStressVector[5] * rNormal[0]
                         + rStressVector[4] * rNormal[1]
                         + rStressVector[2] * rNormal[2];
        }
    }

    void SbmLoadSolidCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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


    void SbmLoadSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
        
        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * mDim;

        // Initialize DN_DX
        Matrix DN_DX(number_of_control_points, mDim);

        const GeometryType::ShapeFunctionsGradientsType& r_DN_De =
            r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        Vector displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(displacement_coefficient_vector);

        // Compute physical derivatives
        if (mDim == 2) {
            Matrix jacobian;
            CalculateInitialJacobian(r_geometry, jacobian);

            Matrix inv_jacobian(mDim, mDim);
            double det_jacobian = 0.0;
            MathUtils<double>::InvertMatrix(jacobian, inv_jacobian, det_jacobian);

            noalias(DN_DX) = prod(r_DN_De[0], inv_jacobian);
        } else {
            noalias(DN_DX) = r_DN_De[0];
        }

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix B = ZeroMatrix(strain_size,mat_size);
        CalculateB(B, DN_DX);

        // GET STRESS VECTOR
        ConstitutiveLaw::Parameters values_surrogate(r_geometry, GetProperties(), rCurrentProcessInfo);
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=values_surrogate.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B,displacement_coefficient_vector);
    
        values_surrogate.SetStrainVector(old_strain);
        values_surrogate.SetStressVector(this_constitutive_variables.StressVector);
        values_surrogate.SetConstitutiveMatrix(this_constitutive_variables.D);
        mpConstitutiveLaw->CalculateMaterialResponse(values_surrogate, ConstitutiveLaw::StressMeasure_Cauchy);

        const Vector sigma = values_surrogate.GetStressVector();
        Vector sigma_n;
        CalculateTraction(sigma, mNormalPhysicalSpace, sigma_n);

        SetValue(NORMAL_STRESS, sigma_n);

        if (mDim == 2) {
            SetValue(CAUCHY_STRESS_XX, sigma[0]);
            SetValue(CAUCHY_STRESS_YY, sigma[1]);
            SetValue(CAUCHY_STRESS_XY, sigma[2]);
        } else {
            SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(sigma));
        }
    }

void SbmLoadSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

}

void SbmLoadSolidCondition::ComputeGradientTaylorExpansionContribution(const Vector& rDistanceVector, Matrix& grad_H_sum)
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
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
            
                IndexType countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (int k_x = static_cast<int>(n); k_x >= 0; --k_x) {
                    // Loop over the possible derivatives in y
                    for (int k_y = static_cast<int>(n) - k_x; k_y >= 0; --k_y) {

                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
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
double SbmLoadSolidCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double SbmLoadSolidCondition::ComputeTaylorTerm3D(
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

} // Namespace Kratos
