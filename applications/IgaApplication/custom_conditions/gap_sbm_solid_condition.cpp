
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

// System includes

// External includes

// Project includes
#include "custom_conditions/gap_sbm_solid_condition.h"

namespace Kratos
{

void GapSbmSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeMaterial();
    InitializeSbmMemberVariables();
}


void GapSbmSolidCondition::InitializeMaterial()
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

void GapSbmSolidCondition::InitializeMemberVariables()
{
    // // Compute class memeber variables
    const auto& r_geometry = GetGeometry();

    const auto& r_projected_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_DN_De = r_projected_geometry.ShapeFunctionsLocalGradients(r_projected_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "GapSbmSolidCondition momentarily only supports 2D and 3D conditions, but the current dimension is" << mDim << std::endl;
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder *= 3; 
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
        mBasisFunctionsOrder *= 2; 
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
        // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH) (/4 because of the local increase *2 for the enhanced shift operator)
        mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder /4* penalty / h;
        if (mDim == 3)
        {
            mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder /9* penalty / h;
        }
    }
    // Compute the normals
    mNormalParameterSpace = r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpace);

    // calculate the integration weight
    // reading integration point
    if (mDim == 2)
    {
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
        const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
        const double integration_weight = r_integration_points[0].Weight()*thickness;
        SetValue(INTEGRATION_WEIGHT, integration_weight);
    }
    else //mDim == 3
    {
        const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
        const double integration_weight = r_integration_points[0].Weight();
        SetValue(INTEGRATION_WEIGHT, integration_weight);        
    }
    
}

void GapSbmSolidCondition::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();
    
    Vector center = ZeroVector(3);
    center[0] = 1.0; center[1] = 1.0; center[2] = 1.0;
    Vector true_normal = r_geometry.Center().Coordinates() - center;
    true_normal /= MathUtils<double>::Norm(true_normal);
    
    // if (norm_2(mNormalPhysicalSpace - true_normal) > 5e-1)
    // {
    //     SetValue(NORMAL, mNormalPhysicalSpace);
    //     KRATOS_WATCH(mNormalPhysicalSpace)
    //     KRATOS_WATCH(true_normal)
    // //     KRATOS_WATCH("error with the normal")
    // //     exit(0);
    // // }
    // if (inner_prod(mNormalPhysicalSpace, mDistanceVector) < 0.0) {
    //     // KRATOS_WATCH("ERROR WITH NORMALS!!!!!!!!!!!!!!!!!!!!!!!")
    //     // KRATOS_WATCH(true_normal)
    //     // KRATOS_WATCH(mNormalPhysicalSpace)
    //     // KRATOS_WATCH("--------------------")

    //     mNormalPhysicalSpace *= -1;
    //     SetValue(NORMAL, mNormalPhysicalSpace);
    //     exit(0);
    //     this->SetValue(ACTIVATION_LEVEL, 4.0);
    // }

}

void GapSbmSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const std::size_t mat_size = GetValue(NEIGHBOUR_GEOMETRIES)[0]->size() * mDim;

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

void GapSbmSolidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_true_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const std::size_t mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(N_sum_vec);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true position
    ConstitutiveLaw::Parameters values_true(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);
    const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Matrix& r_D_on_true = values_true.GetConstitutiveMatrix();

    Matrix DB_sum = prod(r_D_on_true, B_sum);

    // Differential area
    double penalty_integration = mPenalty * integration_weight;

    // ASSEMBLE
    //-----------------------------------------------------
    if (this->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = this->GetValue(DIRECTION);

        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType j = 0; j < number_of_control_points; j++) {
                
                for (IndexType idim = 0; idim < mDim; idim++) {
                    const int iglob = mDim*i+idim;

                    for (IndexType jdim = 0; jdim < mDim; jdim++) {
                        const int jglob = mDim*j+jdim;
    
                        // PENALTY TERM
                        rLeftHandSideMatrix(iglob, jglob) += N_sum_vec(i)*N_sum_vec(j)* penalty_integration * direction[idim] * direction[jdim];

                        // FLUX 
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        Vector sigma_u_n;
                        Vector stress_column_u = column(DB_sum, jglob);
                        CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n);

                        double sigma_u_n_dot_direction = inner_prod(sigma_u_n, direction);

                        rLeftHandSideMatrix(iglob, jglob) -= N_sum_vec(i) * sigma_u_n_dot_direction * direction[idim] * integration_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        Vector sigma_w_n;
                        Vector stress_column_w = column(DB_sum, iglob);
                        CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

                        double sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);
                        rLeftHandSideMatrix(iglob, jglob) -= mNitschePenalty*N_sum_vec(j) * sigma_w_n_dot_direction * direction[jdim] * integration_weight;
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
                    rLeftHandSideMatrix(iglob, mDim*j+idim) += N_sum_vec(i)*N_sum_vec(j)* penalty_integration;

                    Vector sigma_w_n = ZeroVector(3);
                    Vector stress_column_w = column(DB_sum, iglob);
                    CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);
                    
                    for (IndexType jdim = 0; jdim < mDim; jdim++) {
                        const int jglob = mDim*j+jdim;

                        Vector sigma_u_n;
                        Vector stress_u_column = column(DB_sum, jglob);
                        CalculateTraction(stress_u_column, mNormalPhysicalSpace, sigma_u_n);

                        // FLUX 
                        // [sigma(u) \dot n] \dot n * (-w \dot n)
                        //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= N_sum_vec(i)*sigma_u_n[idim] * integration_weight;

                        // // PENALTY FREE g_n = 0
                        // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                        // //*********************************************** */
                        rLeftHandSideMatrix(iglob, jglob) -= mNitschePenalty*N_sum_vec(j)*sigma_w_n[jdim] * integration_weight;
                    }

                }
            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const auto& r_true_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();

    // reading integration points and local gradients
    const std::size_t mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(N_sum_vec);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_sum = ZeroMatrix(strain_size,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true positio
    ConstitutiveLaw::Parameters values_true(r_true_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

    const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Matrix& r_D_on_true = values_true.GetConstitutiveMatrix();
    const Vector& r_stress_vector_on_true = values_true.GetStressVector();

    // compute the old_displacement solution on the true boundary
    Vector old_displacement = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        old_displacement[0] += N_sum_vec(i) * old_displacement_coefficient_vector[mDim*i];
        old_displacement[1] += N_sum_vec(i) * old_displacement_coefficient_vector[mDim*i + 1];

        if (mDim == 3) {
            old_displacement[2] += N_sum_vec(i) * old_displacement_coefficient_vector[mDim*i + 2];
        }
    }

    Matrix DB_sum = prod(r_D_on_true, B_sum);

    // Differential area
    double penalty_integration = mPenalty * integration_weight;

    // ASSEMBLE
    //-----------------------------------------------------
    Vector u_D = this->GetValue(DISPLACEMENT);

    //FIXME:
    u_D = r_true_geometry.GetValue(PROJECTION_NODE)->GetValue(DISPLACEMENT);

    if (this->Has(DIRECTION)){
        // ASSIGN BC BY DIRECTION
        //--------------------------------------------------------------------------------------------
        Vector direction = this->GetValue(DIRECTION);
        direction /= norm_2(direction);

        const Vector displacement = this->GetValue(DISPLACEMENT);

        const double prescribed_displacement_direction =
            inner_prod(displacement, direction);

        Vector projected_displacement = prescribed_displacement_direction * direction;

        KRATOS_ERROR_IF(norm_2(displacement - projected_displacement) > 1e-12)
            << "::[GapSbmSolidCondition]:: Error: the prescribed displacement is not aligned "
            << "with the prescribed direction." << std::endl;

        const double displacement_module = norm_2(displacement);

        const double old_displacement_direction = inner_prod(old_displacement, direction);
            
        for (IndexType i = 0; i < number_of_control_points; i++) {
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                rRightHandSideVector(iglob) += N_sum_vec(i) * direction[idim] * (displacement_module-old_displacement_direction) * penalty_integration;

                // // PENALTY FREE g_n = 0
                // // rhs -> [\sigma_1(w) \dot n] \dot n (-g_{n,0})
                // //*********************************************** */
                Vector sigma_w_n;
                Vector stress_column_w = column(DB_sum, iglob);
                CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

                double cut_sigma_w_n_dot_direction = inner_prod(sigma_w_n, direction);

                //PENALTY FREE
                // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                rRightHandSideVector(iglob) -= mNitschePenalty*cut_sigma_w_n_dot_direction * integration_weight *(displacement_module - old_displacement_direction);

                // residual terms

                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                CalculateTraction(r_stress_vector_on_true, mNormalPhysicalSpace, old_stress_normal);

                double old_stress_normal_dot_direction = inner_prod(old_stress_normal, direction);
                rRightHandSideVector(iglob) += N_sum_vec(i) * old_stress_normal_dot_direction * direction[idim] * integration_weight;
            }
        }
    }
    else {
        // ASSIGN BC BY COMPONENTS 
            //--------------------------------------------------------------------------------------------
        for (IndexType i = 0; i < number_of_control_points; i++) {
                
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                // PENALTY TERM
                rRightHandSideVector[iglob] += N_sum_vec(i)*(u_D - old_displacement)[idim]* penalty_integration;

                Vector sigma_w_n;
                Vector stress_column_w = column(DB_sum, iglob);
                CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n);

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    rRightHandSideVector(iglob) -= mNitschePenalty*(u_D[jdim]-old_displacement[jdim])*sigma_w_n[jdim] * integration_weight;
                }

                // residual terms
                // FLUX
                Vector old_stress_normal = ZeroVector(3);
                CalculateTraction(r_stress_vector_on_true, mNormalPhysicalSpace, old_stress_normal);
                
                rRightHandSideVector(iglob) += N_sum_vec(i) * old_stress_normal[idim] * integration_weight;

            }
        }
    }


    // // DEBUG PURPOSE ONLY
    // MatrixType lhs;
    // CalculateLeftHandSide(lhs, rCurrentProcessInfo);

    // const std::size_t size = lhs.size1();
    // KRATOS_ERROR_IF(size == 0) << "SolidCouplingCondition::CalculateRightHandSide: The left hand side matrix has zero size." << std::endl;

    // Vector plus_values;
    // GetSolutionCoefficientVector(plus_values);

    // Vector values(plus_values.size());
    // for (IndexType i = 0; i < plus_values.size(); ++i) {
    //     values[i] = plus_values[i];
    // }

    // if (rRightHandSideVector.size() != size) {
    //     rRightHandSideVector.resize(size, false);
    // }
    // noalias(rRightHandSideVector) = -prod(lhs, values);

    // return;

    for (unsigned int i = 0; i < this->GetValue(NEIGHBOUR_GEOMETRIES)[0]->size(); i++) {
        auto& r_node = (*this->GetValue(NEIGHBOUR_GEOMETRIES)[0])[i];
        // if (r_geometry[i].GetId() == 420) 
        // {
        //     KRATOS_WATCH(r_geometry[i].Coordinates())
        //     KRATOS_WATCH(DN_DX(i,0))
        //     KRATOS_WATCH(DN_DX(i,1))
        // }
    
        std::ofstream outputFile("txt_files/Id_active_control_points_condition.txt", std::ios::app);
        outputFile << r_node.GetId() << "  " << r_node.GetDof(DISPLACEMENT_Z).EquationId() <<"\n";
        outputFile.close();
    }

    KRATOS_CATCH("")
}

    int GapSbmSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void GapSbmSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const std::size_t number_of_control_points = r_geometry.size();

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

    void GapSbmSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const std::size_t number_of_control_points = r_geometry.size();

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


    void GapSbmSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const std::size_t number_of_control_points = r_geometry.size();
        const std::size_t mat_size = number_of_control_points * mDim;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = r_geometry[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * mDim;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            if (mDim == 3) {
                rValues[index + 2] = displacement[2];
            }
        }
    }

    void GapSbmSolidCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const SizeType number_of_control_points = r_surrogate_geometry.size();
        const SizeType mat_size = number_of_control_points * mDim;
        const SizeType strain_size = (mDim == 2) ? 3 : 6;
    
        if (rB.size1() != strain_size || rB.size2() != mat_size)
            rB.resize(strain_size, mat_size, false);
    
        noalias(rB) = ZeroMatrix(strain_size, mat_size);
    
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * mDim;
    
            if (mDim == 2) {
                rB(0, index + 0) = r_DN_DX(i, 0); // exx
                rB(1, index + 1) = r_DN_DX(i, 1); // eyy
    
                rB(2, index + 0) = r_DN_DX(i, 1); // gamma_xy
                rB(2, index + 1) = r_DN_DX(i, 0);
            }
            else if (mDim == 3) {
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

    void GapSbmSolidCondition::CalculateTraction(
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

    void GapSbmSolidCondition::ApplyConstitutiveLaw(std::size_t matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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


    void GapSbmSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        ConstitutiveLaw::Parameters constitutive_law_parameters(
            GetGeometry(), GetProperties(), rCurrentProcessInfo);

        mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * mDim;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
        Matrix grad_N_sum = trans(grad_N_sum_transposed);
        
        const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
        Matrix B_sum = ZeroMatrix(strain_size,mat_size);
        CalculateB(B_sum, grad_N_sum);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_true(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

        const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

        const Vector sigma = values_true.GetStressVector();
        Vector sigma_n = ZeroVector(mDim);

        if (mDim == 2) {
            sigma_n[0] = sigma[0] * mNormalPhysicalSpace[0]
                    + sigma[2] * mNormalPhysicalSpace[1];

            sigma_n[1] = sigma[2] * mNormalPhysicalSpace[0]
                    + sigma[1] * mNormalPhysicalSpace[1];
        } else {
            sigma_n[0] = sigma[0] * mNormalPhysicalSpace[0]
                    + sigma[3] * mNormalPhysicalSpace[1]
                    + sigma[5] * mNormalPhysicalSpace[2];

            sigma_n[1] = sigma[3] * mNormalPhysicalSpace[0]
                    + sigma[1] * mNormalPhysicalSpace[1]
                    + sigma[4] * mNormalPhysicalSpace[2];

            sigma_n[2] = sigma[5] * mNormalPhysicalSpace[0]
                    + sigma[4] * mNormalPhysicalSpace[1]
                    + sigma[2] * mNormalPhysicalSpace[2];
        }

        SetValue(NORMAL_STRESS, sigma_n);

        if (mDim == 2) {
            SetValue(CAUCHY_STRESS_XX, sigma[0]);
            SetValue(CAUCHY_STRESS_YY, sigma[1]);
            SetValue(CAUCHY_STRESS_XY, sigma[2]);
        } else {
            SetValue(CAUCHY_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(sigma));
        }
        // //---------------------

        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);

        array_1d<double, 3> current_displacement = ZeroVector(3);
        for (IndexType i = 0; i < number_of_control_points; ++i) {
            current_displacement[0] += N_sum_vec[i] * old_displacement_coefficient_vector[mDim*i];
            current_displacement[1] += N_sum_vec[i] * old_displacement_coefficient_vector[mDim*i + 1];
            if (mDim == 3) {
                current_displacement[2] += N_sum_vec[i] * old_displacement_coefficient_vector[mDim*i + 2];
            }
        }
        // SetValue(DISPLACEMENT, current_displacement);
    }

void GapSbmSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void GapSbmSolidCondition::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
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

void GapSbmSolidCondition::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
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
double GapSbmSolidCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double GapSbmSolidCondition::ComputeTaylorTerm3D(
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
