
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

    KRATOS_ERROR_IF(mDim != 2) << "SbmSolidCondition momentarily only supports 2D conditions, but the current dimension is" << mDim << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
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
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpace);

    // calculate the integration weight
    // reading integration point
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Initialize Jacobian
    Matrix InvJ0(3,3);
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // Compute the local tangents
    array_1d<double, 3> tangent_parameter_space;

    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space

    // compute complete jacobian transformation including parameter->physical space transformation
    double detJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

    Vector add_factor = prod(Jacobian, tangent_parameter_space); //additional factor to the determinant of the jacobian for the parameter->physical space transformation
    add_factor[2] = 0.0; 
    detJ0 = norm_2(add_factor);

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double int_to_reference_weight = r_integration_points[0].Weight() * std::abs(detJ0) * thickness;

    SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);
}

void SbmSolidCondition::InitializeSbmMemberVariables()
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

void SbmSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetGeometry().size() * 2;

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
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(3,3);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // compute complete jacobian transformation including parameter->physical space transformation
    double detJ0;
    Matrix jacobian = ZeroMatrix(3,3);
    jacobian(0,0) = J0[0](0,0);
    jacobian(0,1) = J0[0](0,1);
    jacobian(1,0) = J0[0](1,0);
    jacobian(1,1) = J0[0](1,1);
    jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(jacobian,InvJ0,detJ0);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    Matrix sub_inv_jacobian = ZeroMatrix(2,2);
    sub_inv_jacobian(0,0) = InvJ0(0,0);
    sub_inv_jacobian(1,0) = InvJ0(1,0);
    sub_inv_jacobian(0,1) = InvJ0(0,1);
    sub_inv_jacobian(1,1) = InvJ0(1,1);
    noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);

    Matrix B = ZeroMatrix(mDim,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);

    const Matrix& r_D = Values.GetConstitutiveMatrix();

    // Differential area
    double penalty_integration = mPenalty * int_to_reference_weight;

    // Collins, Lozinsky & Scovazzi innovation
    double nitsche_penalty = -1.0;  // = 1 -> Penalty approach
                                        // = -1 -> Free-penalty approach
    if (mPenalty == -1.0) {
        penalty_integration = 0.0;
        nitsche_penalty = -1.0;
    }

    const Matrix DB = prod(r_D,B);

    // compute Taylor expansion contribution: H_sum_vec
    Vector H_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution (H_sum_vec);

    // Assembly
    if (this->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = this->GetValue(DIRECTION);

        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType j = 0; j < number_of_control_points; j++) {
                
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int iglob = 2*i+idim;

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int jglob = 2*j+jdim;

                        // PENALTY TERM
                        rLeftHandSideMatrix(iglob, jglob) += r_N(0,i)*H_sum_vec(j)* penalty_integration * direction[idim] * direction[jdim];

                        // FLUX 
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        Vector sigma_u_n = ZeroVector(3);
                        sigma_u_n[0] = DB(0, jglob)*mNormalPhysicalSpace[0] + DB(2, jglob)*mNormalPhysicalSpace[1];
                        sigma_u_n[1] = DB(2, jglob)*mNormalPhysicalSpace[0] + DB(1, jglob)*mNormalPhysicalSpace[1];

                        double sigma_u_n_dot_direction = inner_prod(sigma_u_n, direction);

                        rLeftHandSideMatrix(iglob, jglob) -= r_N(0,i) * sigma_u_n_dot_direction * direction[idim] * int_to_reference_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        Vector sigma_w_n = ZeroVector(3);
                        sigma_w_n[0] = (DB(0, iglob)* mNormalPhysicalSpace[0] + DB(2, iglob)* mNormalPhysicalSpace[1]);
                        sigma_w_n[1] = (DB(2, iglob)* mNormalPhysicalSpace[0] + DB(1, iglob)* mNormalPhysicalSpace[1]);

                        double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                        rLeftHandSideMatrix(iglob, jglob) -= nitsche_penalty*H_sum_vec(j) * sigma_w_n_dot_direction * direction[jdim] * int_to_reference_weight;
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
                
                for (IndexType idim = 0; idim < 2; idim++) {
                    const int id1 = 2*idim;
                    const int iglob = 2*i+idim;

                    // PENALTY TERM
                    rLeftHandSideMatrix(iglob, 2*j+idim) += r_N(0,i)*H_sum_vec(j)* penalty_integration;

                    Vector sigma_w_n = ZeroVector(3);
                    sigma_w_n[0] = (DB(0, iglob)* mNormalPhysicalSpace[0] + DB(2, iglob)* mNormalPhysicalSpace[1]);
                    sigma_w_n[1] = (DB(2, iglob)* mNormalPhysicalSpace[0] + DB(1, iglob)* mNormalPhysicalSpace[1]);

                    for (IndexType jdim = 0; jdim < 2; jdim++) {
                        const int id2 = (id1+2)%3;
                        const int jglob = 2*j+jdim;

                        // FLUX 
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= r_N(0,i)*(DB(id1, jglob)* mNormalPhysicalSpace[0] + DB(id2, jglob)* mNormalPhysicalSpace[1]) * int_to_reference_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= nitsche_penalty*H_sum_vec(j)*sigma_w_n[jdim] * int_to_reference_weight;
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
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(3,3);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // compute complete jacobian transformation including parameter->physical space transformation
    double detJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    Matrix sub_inv_jacobian = ZeroMatrix(2,2);
    sub_inv_jacobian(0,0) = InvJ0(0,0);
    sub_inv_jacobian(1,0) = InvJ0(1,0);
    sub_inv_jacobian(0,1) = InvJ0(0,1);
    sub_inv_jacobian(1,1) = InvJ0(1,1);
    noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);

    Matrix B = ZeroMatrix(mDim,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, Values, this_constitutive_variables);

    const Matrix& r_D = Values.GetConstitutiveMatrix();
    const Vector& r_stress_vector = Values.GetStressVector();

    // Differential area
    double penalty_integration = mPenalty * int_to_reference_weight;

    // Collins, Lozinsky & Scovazzi innovation
    double nitsche_penalty = -1.0;  // = 1 -> Penalty approach
                                        // = -1 -> Free-penalty approach
    if (mPenalty == -1.0) {
        penalty_integration = 0.0;
        nitsche_penalty = -1.0;
    }

    // Assembly

    const Matrix DB = prod(r_D,B);
    // compute Taylor expansion contribution: H_sum_vec
    Vector H_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution (H_sum_vec);

    Vector old_displacement = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        old_displacement[0] += H_sum_vec(i) * old_displacement_coefficient_vector[2*i];
        old_displacement[1] += H_sum_vec(i) * old_displacement_coefficient_vector[2*i + 1];
    }

    if (this->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = this->GetValue(DIRECTION);
        const Vector displacement = this->GetValue(DISPLACEMENT); //already direction
        const double displacement_module = inner_prod(displacement, direction);

        const double old_displacement_direction = inner_prod(old_displacement, direction);
            
        for (IndexType i = 0; i < number_of_control_points; i++) {

            for (IndexType idim = 0; idim < 2; idim++) {
                const int iglob = 2*i+idim;

                rRightHandSideVector(iglob) += r_N(0,i) * direction[idim] * (displacement_module-old_displacement_direction) * penalty_integration;

                // // PENALTY FREE g_n = 0
                // // rhs -> [\sigma_1(w) \dot n] \dot n (-g_{n,0})
                // //*********************************************** */
                Vector sigma_w_n = ZeroVector(3);
                sigma_w_n[0] = (DB(0, iglob)* mNormalPhysicalSpace[0] + DB(2, iglob)* mNormalPhysicalSpace[1]);
                sigma_w_n[1] = (DB(2, iglob)* mNormalPhysicalSpace[0] + DB(1, iglob)* mNormalPhysicalSpace[1]);

                double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                //PENALTY FREE
                // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                rRightHandSideVector(iglob) -= nitsche_penalty*sigma_w_n_dot_direction * int_to_reference_weight *(displacement_module - old_displacement_direction);

                // residual terms

                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                old_stress_normal[0] = r_stress_vector[0]*mNormalPhysicalSpace[0] + r_stress_vector[2]*mNormalPhysicalSpace[1];
                old_stress_normal[1] = r_stress_vector[2]*mNormalPhysicalSpace[0] + r_stress_vector[1]*mNormalPhysicalSpace[1];

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

            for (IndexType idim = 0; idim < 2; idim++) {
                const int iglob = 2*i+idim;

                rRightHandSideVector[iglob] += r_N(0,i)*(u_D-old_displacement)[idim]* penalty_integration;

                Vector sigma_w_n = ZeroVector(3);
                sigma_w_n[0] = (DB(0, iglob)* mNormalPhysicalSpace[0] + DB(2, iglob)* mNormalPhysicalSpace[1]);
                sigma_w_n[1] = (DB(2, iglob)* mNormalPhysicalSpace[0] + DB(1, iglob)* mNormalPhysicalSpace[1]);

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    rRightHandSideVector(iglob) -= nitsche_penalty*(u_D[jdim]-old_displacement[jdim])*sigma_w_n[jdim] * int_to_reference_weight;
                }

                // residual terms
                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                old_stress_normal[0] = r_stress_vector[0]*mNormalPhysicalSpace[0] + r_stress_vector[2]*mNormalPhysicalSpace[1];
                old_stress_normal[1] = r_stress_vector[2]*mNormalPhysicalSpace[0] + r_stress_vector[1]*mNormalPhysicalSpace[1];

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

        if (rResult.size() != 2 * number_of_control_points)
            rResult.resize(2 * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
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
        rElementalDofList.reserve(2 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }
    };


    void SbmSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    void SbmSolidCondition::CalculateB(
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
        const SizeType mat_size = number_of_control_points * 2;

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        const unsigned int mDim = 2;
        Matrix DN_DX(number_of_control_points,2);
        Matrix InvJ0(mDim,mDim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        r_geometry.Jacobian(J0,this->GetIntegrationMethod());
        double detJ0;
        // MODIFIED
        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);
        
        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        // MODIFIED
        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX);

        // GET STRESS VECTOR
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
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
        Vector sigma_n(2);

        sigma_n[0] = sigma[0]*mNormalPhysicalSpace[0] + sigma[2]*mNormalPhysicalSpace[1];
        sigma_n[1] = sigma[2]*mNormalPhysicalSpace[0] + sigma[1]*mNormalPhysicalSpace[1];

        SetValue(NORMAL_STRESS, sigma_n);

        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
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
