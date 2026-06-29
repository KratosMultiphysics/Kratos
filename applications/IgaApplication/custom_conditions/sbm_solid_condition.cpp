
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
#include "custom_conditions/sbm_solid_condition.h"

namespace Kratos
{

void SbmSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}


void SbmSolidCondition::InitializeMaterial()
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

void SbmSolidCondition::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "SbmSolidCondition momentarily only supports 2D and 3D conditions, but the current dimension is" << mDim << std::endl;

    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}

    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder*=3; //use complete Talor expansion
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder*=2; //use complete Talor expansion
    }

    double penalty = GetProperties()[PENALTY_FACTOR];

    // https://doi.org/10.1016/j.cma.2023.116301 (A penalty-free Shifted Boundary Method of arbitrary order)
    mNitschePenalty = 1.0;   // = 1.0 -> Penalty approach
                                    // = -1.0 -> Free-penalty approach
    if (penalty == -1.0) {
        mPenalty = 0.0;
        mNitschePenalty = -1.0;
    }
    else
    {
        // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
        mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * penalty / h;
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

void SbmSolidCondition::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    // Retrieve projection
    Condition candidate_closest_skin_segment_1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
    // Find the closest node in condition
    int closestNodeId = 0;
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

    this->SetValue(PROJECTION_NODE_COORDINATES, mpProjectionNode->Coordinates());

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = mpProjectionNode->Coordinates() - r_geometry.Center().Coordinates();
}

void SbmSolidCondition::CalculateLocalSystem(
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

void SbmSolidCondition::CalculateLeftHandSide(
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

    // Differential area
    double penalty_integration = mPenalty * int_to_reference_weight;


    const Matrix DB = prod(r_D,B);

    // compute Taylor expansion contribution: H_sum_vec
    Vector H_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution (H_sum_vec);

    // Assembly
    if (mpProjectionNode->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = mpProjectionNode->GetValue(DIRECTION);

        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType j = 0; j < number_of_control_points; j++) {

                for (IndexType idim = 0; idim < mDim; idim++) {
                    const int iglob = mDim*i+idim;

                    for (IndexType jdim = 0; jdim < mDim; jdim++) {
                        const int jglob = mDim*j+jdim;

                        // PENALTY TERM
                        rLeftHandSideMatrix(iglob, jglob) += r_N(0,i)*H_sum_vec(j)* penalty_integration * direction[idim] * direction[jdim];

                        // FLUX
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        Vector sigma_u_n;
                        Vector stress_column_u = column(DB, jglob);
                        CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n);

                        double sigma_u_n_dot_direction = inner_prod(sigma_u_n, direction);

                        rLeftHandSideMatrix(iglob, jglob) -= r_N(0,i) * sigma_u_n_dot_direction * direction[idim] * int_to_reference_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        Vector sigma_w_n;
                        Vector stress_column_w = column(DB, iglob);
                        CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

                        double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                        rLeftHandSideMatrix(iglob, jglob) -= mNitschePenalty*H_sum_vec(j) * sigma_w_n_dot_direction * direction[jdim] * int_to_reference_weight;
                    }

                }
            }
        }
    }
    else {
        // ASSIGN BC BY COMPONENTS
        //--------------------------------------------------------------------------------------------
        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType j = 0; j < number_of_control_points; j++) {

                for (IndexType idim = 0; idim < mDim; idim++) {
                    const int iglob = mDim*i+idim;

                    // PENALTY TERM
                    rLeftHandSideMatrix(iglob, mDim*j+idim) += r_N(0,i)*H_sum_vec(j)* penalty_integration;

                    Vector sigma_w_n;
                    Vector stress_column_w = column(DB, iglob);
                    CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

                    for (IndexType jdim = 0; jdim < mDim; jdim++) {
                        const int jglob = mDim*j+jdim;

                        Vector sigma_u_n;
                        Vector stress_u_column = column(DB, jglob);
                        CalculateTraction(stress_u_column, mNormalPhysicalSpace, sigma_u_n);

                        // FLUX
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= r_N(0,i)*sigma_u_n[idim] * int_to_reference_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= mNitschePenalty*H_sum_vec(j)*sigma_w_n[jdim] * int_to_reference_weight;
                    }

                }
            }
        }
    }
    KRATOS_CATCH("")
}

void SbmSolidCondition::CalculateRightHandSide(
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
    Matrix DN_DX(number_of_control_points,mDim);

    // Calculating the physical derivatives (it is avoided storing them to minimize storage)
    if (mDim == 2) {
        Matrix InvJ0(mDim,mDim);
        double detJ0;
        Matrix Jacobian;
        CalculateInitialJacobian(r_geometry, Jacobian);
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

        Matrix sub_inv_jacobian = ZeroMatrix(mDim,mDim);
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                sub_inv_jacobian(i,j) = InvJ0(i,j);
            }
        }
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
    const Vector& r_stress_vector = Values.GetStressVector();

    // Differential area
    double penalty_integration = mPenalty * int_to_reference_weight;


    // Assembly

    const Matrix DB = prod(r_D,B);
    // compute Taylor expansion contribution: H_sum_vec
    Vector H_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution (H_sum_vec);

    Vector old_displacement = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        old_displacement[0] += H_sum_vec(i) * old_displacement_coefficient_vector[mDim*i];
        old_displacement[1] += H_sum_vec(i) * old_displacement_coefficient_vector[mDim*i + 1];

        if (mDim == 3) {
            old_displacement[2] += H_sum_vec(i) * old_displacement_coefficient_vector[mDim*i + 2];
        }
    }

    if (mpProjectionNode->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = mpProjectionNode->GetValue(DIRECTION);
        direction /= norm_2(direction);

        const Vector& r_displacement = mpProjectionNode->GetValue(DISPLACEMENT);

        const double prescribed_displacement_direction =
            inner_prod(r_displacement, direction);

        Vector projected_displacement = prescribed_displacement_direction * direction;

        KRATOS_ERROR_IF(norm_2(r_displacement - projected_displacement) > 1e-12)
            << "::[SbmSolidCondition]:: Error: the projection node displacement is not aligned "
            << "with the prescribed direction." << std::endl;

        const double displacement_module = norm_2(r_displacement);

        const double old_displacement_direction = inner_prod(old_displacement, direction);

        for (IndexType i = 0; i < number_of_control_points; i++) {

            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                rRightHandSideVector(iglob) += r_N(0,i) * direction[idim] * (displacement_module-old_displacement_direction) * penalty_integration;

                // // PENALTY FREE g_n = 0
                // // rhs -> [\sigma_1(w) \dot n] \dot n (-g_{n,0})
                // //*********************************************** */
                Vector sigma_w_n;
                Vector stress_column_w = column(DB, iglob);
                CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

                double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                //PENALTY FREE
                // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                rRightHandSideVector(iglob) -= mNitschePenalty*sigma_w_n_dot_direction * int_to_reference_weight *(displacement_module - old_displacement_direction);

                // residual terms

                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                CalculateTraction(r_stress_vector, mNormalPhysicalSpace, old_stress_normal);

                double old_stress_normal_dot_direction = inner_prod(old_stress_normal, direction);
                rRightHandSideVector(iglob) += r_N(0,i) * old_stress_normal_dot_direction * direction[idim] * int_to_reference_weight;
            }
        }
    }
    else {
        // ASSIGN BC BY COMPONENTS
        //--------------------------------------------------------------------------------------------

        Vector u_D = mpProjectionNode->GetValue(DISPLACEMENT);

        for (IndexType i = 0; i < number_of_control_points; i++) {

            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                rRightHandSideVector[iglob] += r_N(0,i)*(u_D-old_displacement)[idim]* penalty_integration;

                Vector sigma_w_n;
                Vector stress_column_w = column(DB, iglob);
                CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);


                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    rRightHandSideVector(iglob) -= mNitschePenalty*(u_D[jdim]-old_displacement[jdim])*sigma_w_n[jdim] * int_to_reference_weight;
                }

                // residual terms
                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                CalculateTraction(r_stress_vector, mNormalPhysicalSpace, old_stress_normal);

                rRightHandSideVector(iglob) += r_N(0,i) * old_stress_normal[idim] * int_to_reference_weight;
            }
        }
    }

    KRATOS_CATCH("")
}

    int SbmSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SbmSolidCondition::EquationIdVector(
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

    void SbmSolidCondition::GetDofList(
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


    void SbmSolidCondition::GetSolutionCoefficientVector(
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

    void SbmSolidCondition::CalculateB(
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

    void SbmSolidCondition::CalculateTraction(
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

    void SbmSolidCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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


    void SbmSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
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
        Matrix DN_DX(number_of_control_points,mDim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        // MODIFIED
        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // Calculating the physical derivatives (it is avoided storing them to minimize storage)
        if (mDim == 2) {
            Matrix InvJ0(mDim,mDim);
            double detJ0;
            Matrix Jacobian;
            CalculateInitialJacobian(r_geometry, Jacobian);
            MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);
            noalias(DN_DX) = prod(DN_De[0],InvJ0);
        } else {
            noalias(DN_DX) = DN_De[0];
        }

        // MODIFIED
        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix B = ZeroMatrix(strain_size,mat_size);

        CalculateB(B, DN_DX);

        // GET STRESS VECTOR
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B,old_displacement);

        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        const Vector sigma = Values.GetStressVector();
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
        // //---------------------
    }

void SbmSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void SbmSolidCondition::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_control_points = r_geometry.PointsNumber();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    if (H_sum_vec.size() != number_of_control_points)
    {
        H_sum_vec = ZeroVector(number_of_control_points);
    }

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    for (IndexType i = 0; i < number_of_control_points; ++i)
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
                for (int k_x = static_cast<int>(n); k_x >= 0; --k_x) {
                    // Loop over the possible derivatives in y
                    for (int k_y = static_cast<int>(n) - k_x; k_y >= 0; --k_y) {

                        // derivatives in z
                        const IndexType k_z = n - static_cast<IndexType>(k_x) - static_cast<IndexType>(k_y);
                        double derivative = r_shape_function_derivatives(i,countDerivativeId);

                        H_taylor_term += ComputeTaylorTerm3D(derivative, mDistanceVector[0], static_cast<IndexType>(k_x), mDistanceVector[1], static_cast<IndexType>(k_y), mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0,i);
    }
}

// Function to compute a single term in the Taylor expansion
double SbmSolidCondition::ComputeTaylorTerm(
    const double derivative,
    const double dx,
    const IndexType n_k,
    const double dy,
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

double SbmSolidCondition::ComputeTaylorTerm3D(
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
