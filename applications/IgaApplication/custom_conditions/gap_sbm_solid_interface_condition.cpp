
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
#include "custom_conditions/gap_sbm_solid_interface_condition.h"
#include "geometries/coupling_geometry.h"

namespace Kratos
{

void GapSbmSolidInterfaceCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables(GetGeometry());
    InitializeSbmMemberVariables(GetGeometry());
}


void GapSbmSolidInterfaceCondition::InitializeMaterial()
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

void GapSbmSolidInterfaceCondition::InitializeMemberVariables(
    const GeometryType& rTrueGeometry)
{
    // // Compute class memeber variables
    const auto& r_geometry = rTrueGeometry;

    const auto& r_surrogate_geometry = GetGeometryPlus();
    const auto& r_DN_De = r_surrogate_geometry.ShapeFunctionsLocalGradients(r_surrogate_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2 && mDim != 3) << "GapSbmSolidInterfaceCondition momentarily only supports 2D and 3D conditions, but the current dimension is" << mDim << std::endl;
    
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
        // Modify the penalty factor: p^2 * penalty / h (NITSCHE APPROACH)
        mPenalty = mBasisFunctionsOrder * mBasisFunctionsOrder * penalty / h;
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

void GapSbmSolidInterfaceCondition::InitializeSbmMemberVariables(
    const GeometryType& rTrueGeometry)
{
    //TODO:
    const auto& r_geometry = rTrueGeometry;
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();

    mDistanceVectorPlus.resize(3);
    noalias(mDistanceVectorPlus) = r_geometry.Center().Coordinates() - r_surrogate_geometry_plus.Center().Coordinates();

    mDistanceVectorMinus.resize(3);
    noalias(mDistanceVectorMinus) = r_geometry.Center().Coordinates() - r_surrogate_geometry_minus.Center().Coordinates();

}

void GapSbmSolidInterfaceCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const std::size_t number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;

    const std::size_t mat_size = number_of_control_points * mDim;

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

void GapSbmSolidInterfaceCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();
    const auto& r_true_geometry = GetGeometry();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const std::size_t number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;

    const std::size_t mat_size_plus = number_of_control_points_plus * mDim;
    const std::size_t mat_size_minus = number_of_control_points_minus * mDim;
    const std::size_t mat_size = number_of_control_points * mDim;

    // reading integration points and local gradients
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec_plus = ZeroVector(number_of_control_points_plus);
    ComputeTaylorExpansionContribution(r_surrogate_geometry_plus, mDistanceVectorPlus, N_sum_vec_plus);

    Vector N_sum_vec_minus = ZeroVector(number_of_control_points_minus);
    ComputeTaylorExpansionContribution(r_surrogate_geometry_minus, mDistanceVectorMinus, N_sum_vec_minus);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed_plus = ZeroMatrix(3, number_of_control_points_plus);
    ComputeGradientTaylorExpansionContribution(r_surrogate_geometry_plus, mDistanceVectorPlus, grad_N_sum_transposed_plus);
    Matrix grad_N_sum_plus = trans(grad_N_sum_transposed_plus);

    Matrix grad_N_sum_transposed_minus = ZeroMatrix(3, number_of_control_points_minus);
    ComputeGradientTaylorExpansionContribution(r_surrogate_geometry_minus, mDistanceVectorMinus, grad_N_sum_transposed_minus);
    Matrix grad_N_sum_minus = trans(grad_N_sum_transposed_minus);


    // compute the B matrix
    const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_sum_plus = ZeroMatrix(strain_size,mat_size_plus);
    CalculateB(r_surrogate_geometry_plus, B_sum_plus, grad_N_sum_plus);

    Matrix B_sum_minus = ZeroMatrix(strain_size,mat_size_minus);
    CalculateB(r_surrogate_geometry_minus, B_sum_minus, grad_N_sum_minus);

    // obtain the tangent constitutive matrix at the true position for the plus side
    ConstitutiveLaw::Parameters values_true_plus(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector_plus(mat_size_plus);
    GetSolutionCoefficientVectorPlus(old_displacement_coefficient_vector_plus);
    Vector old_strain_on_true_plus = prod(B_sum_plus, old_displacement_coefficient_vector_plus);
    const std::size_t strain_size_true_plus = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true_plus(strain_size_true_plus);
    ApplyConstitutiveLaw(mat_size_plus, old_strain_on_true_plus, values_true_plus, this_constitutive_variables_true_plus);

    const Matrix& r_D_on_true_plus = values_true_plus.GetConstitutiveMatrix();

    // obtain the tangent constitutive matrix at the true position for the minus side
    ConstitutiveLaw::Parameters values_true_minus(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector_minus(mat_size_minus);
    GetSolutionCoefficientVectorMinus(old_displacement_coefficient_vector_minus);
    Vector old_strain_on_true_minus = prod(B_sum_minus, old_displacement_coefficient_vector_minus);
    const std::size_t strain_size_true_minus = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true_minus(strain_size_true_minus);
    ApplyConstitutiveLaw(mat_size_minus, old_strain_on_true_minus, values_true_minus, this_constitutive_variables_true_minus);

    const Matrix& r_D_on_true_minus = values_true_minus.GetConstitutiveMatrix();

    // compute the DB product
    Matrix DB_sum_plus = prod(r_D_on_true_plus, B_sum_plus);
    Matrix DB_sum_minus = prod(r_D_on_true_minus, B_sum_minus);

    // ASSEMBLE
    //-----------------------------------------------------
    const std::size_t shift_dof = mat_size_plus;
    // -w_plus * (sigma_plus + sigma_minus) /2
    for (IndexType i = 0; i < number_of_control_points_plus; i++) {
        // FIRST TERM: -w_plus * sigma_plus /2
        for (IndexType j = 0; j < number_of_control_points_plus; j++) {
            
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim;

                    Vector sigma_u_n_plus = ZeroVector(3);
                    Vector stress_column_u = column(DB_sum_plus, jglob);
                    CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n_plus);                    

                    rLeftHandSideMatrix(iglob, jglob) -= 0.5 * N_sum_vec_plus(i)
                                                    * sigma_u_n_plus[idim] 
                                                    * integration_weight;

                }
            }
        }

        // SECOND TERM: -w_plus * sigma_minus /2
        for (IndexType j = 0; j < number_of_control_points_minus; j++) {
            
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jloc = mDim*j+jdim;
                    const int jglob = mDim*j+jdim + shift_dof;

                    Vector sigma_u_n_minus = ZeroVector(3);
                    Vector stress_column_u = column(DB_sum_minus, jloc);
                    CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n_minus); 

                    rLeftHandSideMatrix(iglob, jglob) -= 0.5 * N_sum_vec_plus(i)
                                                    * sigma_u_n_minus[idim] 
                                                    * integration_weight;
                }
            }
        }
    }

    // +w_minus * (sigma_plus + sigma_minus) /2
    for (IndexType i = 0; i < number_of_control_points_minus; i++) {
        // FIRST TERM: +w_minus * sigma_plus /2
        for (IndexType j = 0; j < number_of_control_points_plus; j++) {
            
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim + shift_dof;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim;

                    Vector sigma_u_n_plus = ZeroVector(3);
                    Vector stress_column_u = column(DB_sum_plus, jglob);
                    CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n_plus); 

                    rLeftHandSideMatrix(iglob, jglob) += 0.5 * N_sum_vec_minus(i)
                                                    * sigma_u_n_plus[idim] 
                                                    * integration_weight;

                }
            }
        }

        // SECOND TERM: +w_minus * sigma_minus /2
        for (IndexType j = 0; j < number_of_control_points_minus; j++) {
            
            for (IndexType idim = 0; idim < mDim; idim++) {
                const int iglob = mDim*i+idim + shift_dof;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jloc = mDim*j+jdim;
                    const int jglob = mDim*j+jdim + shift_dof;

                    Vector sigma_u_n_minus = ZeroVector(3);
                    Vector stress_column_u = column(DB_sum_minus, jloc);
                    CalculateTraction(stress_column_u, mNormalPhysicalSpace, sigma_u_n_minus); 

                    rLeftHandSideMatrix(iglob, jglob) += 0.5 * N_sum_vec_minus(i)
                                                    * sigma_u_n_minus[idim] 
                                                    * integration_weight;
                }
            }
        }
    }


    // CONTINUITY TERMS
    // Differential area
    double penalty_integration = mPenalty * integration_weight;

    // ASSEMBLE
    //-----------------------------------------------------
    for (IndexType i = 0; i < number_of_control_points_plus; i++) {
        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iglob = mDim*i+idim;

            Vector sigma_w_n_plus = ZeroVector(3);
            Vector stress_column_w = column(DB_sum_plus, iglob);
            CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n_plus);

            for (IndexType j = 0; j < number_of_control_points_plus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, mDim*j+idim) += N_sum_vec_plus(i)*N_sum_vec_plus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*mNitschePenalty*N_sum_vec_plus(j)*sigma_w_n_plus[jdim] * integration_weight;
                }

            }

            // minus side
            for (IndexType j = 0; j < number_of_control_points_minus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, mDim*j+shift_dof+idim) -= N_sum_vec_plus(i)*N_sum_vec_minus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim + shift_dof;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += 0.5*mNitschePenalty*N_sum_vec_minus(j)*sigma_w_n_plus[jdim] * integration_weight;
                }
            }
        }
    }

    // minus side
    for (IndexType i = 0; i < number_of_control_points_minus; i++) {
        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iloc = mDim*i+idim;
            const int iglob = mDim*i+idim + shift_dof;

            Vector sigma_w_n_minus = ZeroVector(3);
            Vector stress_column_w = column(DB_sum_minus, iloc);
            CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n_minus);

            for (IndexType j = 0; j < number_of_control_points_plus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, mDim*j+idim) -= N_sum_vec_minus(i)*N_sum_vec_plus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*mNitschePenalty*N_sum_vec_plus(j)*sigma_w_n_minus[jdim] * integration_weight;
                }

            }

            // minus side
            for (IndexType j = 0; j < number_of_control_points_minus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, mDim*j+shift_dof+idim) += N_sum_vec_minus(i)*N_sum_vec_minus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < mDim; jdim++) {
                    const int jglob = mDim*j+jdim + shift_dof;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += 0.5*mNitschePenalty*N_sum_vec_minus(j)*sigma_w_n_minus[jdim] * integration_weight;
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmSolidInterfaceCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();
    const auto& r_true_geometry = GetGeometry();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const std::size_t number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;

    const std::size_t mat_size_plus = number_of_control_points_plus * mDim;
    const std::size_t mat_size_minus = number_of_control_points_minus * mDim;
    const std::size_t mat_size = number_of_control_points * mDim;

    // reading integration points and local gradients
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec_plus = ZeroVector(number_of_control_points_plus);
    ComputeTaylorExpansionContribution(r_surrogate_geometry_plus, mDistanceVectorPlus, N_sum_vec_plus);

    Vector N_sum_vec_minus = ZeroVector(number_of_control_points_minus);
    ComputeTaylorExpansionContribution(r_surrogate_geometry_minus, mDistanceVectorMinus, N_sum_vec_minus);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed_plus = ZeroMatrix(3, number_of_control_points_plus);
    ComputeGradientTaylorExpansionContribution(r_surrogate_geometry_plus, mDistanceVectorPlus, grad_N_sum_transposed_plus);
    Matrix grad_N_sum_plus = trans(grad_N_sum_transposed_plus);

    Matrix grad_N_sum_transposed_minus = ZeroMatrix(3, number_of_control_points_minus);
    ComputeGradientTaylorExpansionContribution(r_surrogate_geometry_minus, mDistanceVectorMinus, grad_N_sum_transposed_minus);
    Matrix grad_N_sum_minus = trans(grad_N_sum_transposed_minus);


    // compute the B matrix
    const std::size_t strain_size = mpConstitutiveLaw->GetStrainSize();
    Matrix B_sum_plus = ZeroMatrix(strain_size,mat_size_plus);
    CalculateB(r_surrogate_geometry_plus, B_sum_plus, grad_N_sum_plus);

    Matrix B_sum_minus = ZeroMatrix(strain_size,mat_size_minus);
    CalculateB(r_surrogate_geometry_minus, B_sum_minus, grad_N_sum_minus);

    // obtain the tangent constitutive matrix at the true position for the plus side
    ConstitutiveLaw::Parameters values_true_plus(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector_plus(mat_size_plus);
    GetSolutionCoefficientVectorPlus(old_displacement_coefficient_vector_plus);
    Vector old_strain_on_true_plus = prod(B_sum_plus, old_displacement_coefficient_vector_plus);
    const std::size_t strain_size_true_plus = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true_plus(strain_size_true_plus);
    ApplyConstitutiveLaw(mat_size_plus, old_strain_on_true_plus, values_true_plus, this_constitutive_variables_true_plus);

    const Vector& r_stress_on_true_plus = values_true_plus.GetStressVector();
    const Matrix& r_D_on_true_plus = values_true_plus.GetConstitutiveMatrix();
    Vector old_stress_plus = ZeroVector(3);
    CalculateTraction(r_stress_on_true_plus, mNormalPhysicalSpace, old_stress_plus);

    // obtain the tangent constitutive matrix at the true position for the minus side
    ConstitutiveLaw::Parameters values_true_minus(r_true_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector_minus(mat_size_minus);
    GetSolutionCoefficientVectorMinus(old_displacement_coefficient_vector_minus);
    Vector old_strain_on_true_minus = prod(B_sum_minus, old_displacement_coefficient_vector_minus);
    const std::size_t strain_size_true_minus = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true_minus(strain_size_true_minus);
    ApplyConstitutiveLaw(mat_size_minus, old_strain_on_true_minus, values_true_minus, this_constitutive_variables_true_minus);

    const Vector& r_stress_on_true_minus = values_true_minus.GetStressVector();
    const Matrix& r_D_on_true_minus = values_true_minus.GetConstitutiveMatrix();
    Vector old_stress_minus = ZeroVector(3);
    CalculateTraction(r_stress_on_true_minus, mNormalPhysicalSpace, old_stress_minus);

    // compute the old_displacement solution on the true boundary
    Vector old_displacement_plus = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points_plus; ++i) {
        old_displacement_plus[0] += N_sum_vec_plus(i) * old_displacement_coefficient_vector_plus[mDim*i];
        old_displacement_plus[1] += N_sum_vec_plus(i) * old_displacement_coefficient_vector_plus[mDim*i + 1];
        if (mDim == 3) {
            old_displacement_plus[2] += N_sum_vec_plus(i) * old_displacement_coefficient_vector_plus[mDim*i + 2];
        }
    }

    Vector old_displacement_minus = ZeroVector(3);
    for (IndexType i = 0; i < number_of_control_points_minus; ++i) {
        old_displacement_minus[0] += N_sum_vec_minus(i) * old_displacement_coefficient_vector_minus[mDim*i];
        old_displacement_minus[1] += N_sum_vec_minus(i) * old_displacement_coefficient_vector_minus[mDim*i + 1];
        if (mDim == 3) {
            old_displacement_minus[2] += N_sum_vec_minus(i) * old_displacement_coefficient_vector_minus[mDim*i + 2];
        }
    }

    // ASSEMBLE
    //-----------------------------------------------------
    const std::size_t shift_dof = mat_size_plus;
    // -w_plus * (sigma_plus + sigma_minus) /2
    for (IndexType i = 0; i < number_of_control_points_plus; i++) {

        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iglob = mDim*i+idim;

            rRightHandSideVector[iglob] += N_sum_vec_plus(i) * 0.5*(old_stress_plus + old_stress_minus)[idim] * integration_weight;

        }
    }

    // +w_minus * (sigma_plus + sigma_minus) /2
    for (IndexType i = 0; i < number_of_control_points_minus; i++) {

        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iglob = mDim*i+idim + shift_dof;

            rRightHandSideVector[iglob] -= N_sum_vec_minus(i) * 0.5*(old_stress_plus + old_stress_minus)[idim] * integration_weight;
        }
    }

    // continuity terms
    // Differential area
    double penalty_integration = mPenalty * integration_weight;

    // compute the DB product
    Matrix DB_sum_plus = prod(r_D_on_true_plus, B_sum_plus);
    Matrix DB_sum_minus = prod(r_D_on_true_minus, B_sum_minus);
    // ASSEMBLE
    //-----------------------------------------------------
    for (IndexType i = 0; i < number_of_control_points_plus; i++) {
        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iglob = mDim*i+idim;

            Vector sigma_w_n_plus = ZeroVector(3);
            Vector stress_column_w = column(DB_sum_plus, iglob);
            CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n_plus);

            // PENALTY RESIDUAL TERM
            rRightHandSideVector[iglob] -= N_sum_vec_plus(i) * (old_displacement_plus-old_displacement_minus)[idim] * penalty_integration;

            // PENALTY FREE RESIDUAL TERM
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector[iglob] += 0.5*mNitschePenalty * (old_displacement_plus-old_displacement_minus)[jdim]  
                                            * sigma_w_n_plus[jdim] * integration_weight;
            }
        }
    }


    // minus side
    for (IndexType i = 0; i < number_of_control_points_minus; i++) {
        for (IndexType idim = 0; idim < mDim; idim++) {
            const int iloc = mDim*i+idim;
            const int iglob = mDim*i+idim + shift_dof;

            Vector sigma_w_n_minus = ZeroVector(3);
            Vector stress_column_w = column(DB_sum_minus, iloc);
            CalculateTraction(stress_column_w, mNormalPhysicalSpace, sigma_w_n_minus);

            // PENALTY RESIDUAL TERM
            rRightHandSideVector[iglob] += N_sum_vec_minus(i) * (old_displacement_plus-old_displacement_minus)[idim] * penalty_integration;

            // PENALTY FREE RESIDUAL TERM
            for (IndexType jdim = 0; jdim < mDim; jdim++) {
                rRightHandSideVector[iglob] += 0.5*mNitschePenalty * (old_displacement_plus-old_displacement_minus)[jdim]  
                                            * sigma_w_n_minus[jdim] * integration_weight;
            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmSolidInterfaceCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
) const
{   
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const std::size_t number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;

    const std::size_t shift_dof = number_of_control_points_plus * mDim;

    if (rResult.size() != mDim * number_of_control_points)
        rResult.resize(mDim * number_of_control_points, false);
    
    // first the plus geometry
    for (IndexType i = 0; i < number_of_control_points_plus; ++i) {
        const IndexType index = i * mDim;
        const auto& r_node = r_surrogate_geometry_plus[i];
        rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        if (mDim == 3) {
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    // then the minus geometry
    for (IndexType i = 0; i < number_of_control_points_minus; ++i) {
        const IndexType index = i * mDim + shift_dof;
        const auto& r_node = r_surrogate_geometry_minus[i];
        rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        if (mDim == 3) {
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }
}

void GapSbmSolidInterfaceCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const std::size_t number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;

    rElementalDofList.resize(0);
    rElementalDofList.reserve(mDim * number_of_control_points);

    // first the plus geometry
    for (IndexType i = 0; i < number_of_control_points_plus; ++i) {
        const auto& r_node = r_surrogate_geometry_plus[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        if (mDim == 3) {
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    }

    // then the minus geometry
    for (IndexType i = 0; i < number_of_control_points_minus; ++i) {
        const auto& r_node = r_surrogate_geometry_minus[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        if (mDim == 3) {
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    }
};


void GapSbmSolidInterfaceCondition::GetSolutionCoefficientVectorPlus(
    Vector& rValues) const
{
    const auto& r_surrogate_geometry_plus = GetGeometryPlus();

    const std::size_t number_of_control_points_plus = r_surrogate_geometry_plus.size();

    const std::size_t mat_size = number_of_control_points_plus * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points_plus; ++i)
    {
        const array_1d<double, 3 >& displacement = r_surrogate_geometry_plus[i].GetSolutionStepValue(DISPLACEMENT);
        IndexType index = i * mDim;

        for (IndexType idim = 0; idim < mDim; ++idim) {
            rValues[index + idim] = displacement[idim];
        }
    }
}

void GapSbmSolidInterfaceCondition::GetSolutionCoefficientVectorMinus(
    Vector& rValues) const
{
    const auto& r_surrogate_geometry_minus = GetGeometryMinus();

    const std::size_t number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const std::size_t mat_size = number_of_control_points_minus * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points_minus; ++i)
    {
        const array_1d<double, 3 >& displacement = r_surrogate_geometry_minus[i].GetSolutionStepValue(DISPLACEMENT);
        IndexType index = i * mDim;

        for (IndexType idim = 0; idim < mDim; ++idim) {
            rValues[index + idim] = displacement[idim];
        }
    }
}

void GapSbmSolidInterfaceCondition::CalculateB(
    const GeometryType& rGeometry,
    Matrix& rB, 
    Matrix& r_DN_DX) const
{
    const SizeType number_of_control_points = rGeometry.size();
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

void GapSbmSolidInterfaceCondition::CalculateTraction(
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

void GapSbmSolidInterfaceCondition::ApplyConstitutiveLaw(std::size_t matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
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


void GapSbmSolidInterfaceCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void GapSbmSolidInterfaceCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    //--------------------------------------------------------------------------------------------
    // calculate the constitutive law response
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void GapSbmSolidInterfaceCondition::ComputeTaylorExpansionContribution(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector, 
    Vector& H_sum_vec)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();

    if (H_sum_vec.size() != number_of_control_points)
    {
        H_sum_vec = ZeroVector(number_of_control_points);
    }

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
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
                    H_taylor_term += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
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

                        H_taylor_term += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0,i);
    }
}

void GapSbmSolidInterfaceCondition::ComputeGradientTaylorExpansionContribution(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector, 
    Matrix& grad_H_sum)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(rGeometry.GetDefaultIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
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
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
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
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x-1, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y-1, rDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z-1);
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
double GapSbmSolidInterfaceCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double GapSbmSolidInterfaceCondition::ComputeTaylorTerm3D(
    const double derivative, 
    const double dx, 
    const IndexType k_x, 
    const double dy, 
    const IndexType k_y, 
    const double dz, 
    const IndexType k_z)
{   
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) *
        std::pow(dz, k_z) /
        (MathUtils<double>::Factorial(k_x) *
         MathUtils<double>::Factorial(k_y) *
         MathUtils<double>::Factorial(k_z));
}

GapSbmSolidInterfaceCondition::TaylorDerivativeData
GapSbmSolidInterfaceCondition::CalculateTaylorDerivativeData(
    const GeometryType& rGeometry) const
{
    TaylorDerivativeData derivative_data(mBasisFunctionsOrder);
    for (IndexType order = 1; order <= mBasisFunctionsOrder; ++order) {
        derivative_data[order - 1] = rGeometry.ShapeFunctionDerivatives(
            order,
            0,
            this->GetIntegrationMethod());
    }
    return derivative_data;
}

void GapSbmSolidInterfaceCondition::EvaluateTaylorExpansion(
    const GeometryType& rGeometry,
    const Vector& rDistanceVector,
    const TaylorDerivativeData& rDerivativeData,
    Vector* pShapeFunctions,
    Matrix* pShapeFunctionGradients) const
{
    const std::size_t number_of_control_points =
        rGeometry.PointsNumber();
    KRATOS_ERROR_IF(rDerivativeData.size() != mBasisFunctionsOrder)
        << "Invalid Taylor derivative data in interface condition #"
        << Id() << ".\n";

    if (pShapeFunctions != nullptr) {
        if (pShapeFunctions->size() != number_of_control_points) {
            pShapeFunctions->resize(number_of_control_points, false);
        }
        noalias(*pShapeFunctions) = row(
            rGeometry.ShapeFunctionsValues(),
            0);
    }

    if (pShapeFunctionGradients != nullptr) {
        if (pShapeFunctionGradients->size1() != 3 ||
            pShapeFunctionGradients->size2() !=
                number_of_control_points) {
            pShapeFunctionGradients->resize(
                3,
                number_of_control_points,
                false);
        }
        noalias(*pShapeFunctionGradients) =
            ZeroMatrix(3, number_of_control_points);
        const auto& r_local_gradients =
            rGeometry.ShapeFunctionsLocalGradients(
                rGeometry.GetDefaultIntegrationMethod())[0];
        for (std::size_t i = 0; i < number_of_control_points; ++i) {
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                (*pShapeFunctionGradients)(component, i) =
                    r_local_gradients(i, component);
            }
        }
    }

    std::array<std::vector<double>, 3> scaled_powers;
    for (std::size_t component = 0; component < 3; ++component) {
        scaled_powers[component].resize(mBasisFunctionsOrder + 1);
        scaled_powers[component][0] = 1.0;
        for (IndexType order = 1;
             order <= mBasisFunctionsOrder;
             ++order) {
            scaled_powers[component][order] =
                scaled_powers[component][order - 1] *
                rDistanceVector[component] /
                static_cast<double>(order);
        }
    }

    for (IndexType order = 1;
         order <= mBasisFunctionsOrder;
         ++order) {
        const Matrix& r_derivatives = rDerivativeData[order - 1];
        if (mDim == 2) {
            for (IndexType y_order = 0;
                 y_order <= order;
                 ++y_order) {
                const IndexType x_order = order - y_order;
                const double value_coefficient =
                    scaled_powers[0][x_order] *
                    scaled_powers[1][y_order];
                const double gradient_x_coefficient = x_order > 0
                    ? scaled_powers[0][x_order - 1] *
                        scaled_powers[1][y_order]
                    : 0.0;
                const double gradient_y_coefficient = y_order > 0
                    ? scaled_powers[0][x_order] *
                        scaled_powers[1][y_order - 1]
                    : 0.0;
                for (std::size_t i = 0;
                     i < number_of_control_points;
                     ++i) {
                    const double derivative =
                        r_derivatives(i, y_order);
                    if (pShapeFunctions != nullptr) {
                        (*pShapeFunctions)[i] +=
                            derivative * value_coefficient;
                    }
                    if (pShapeFunctionGradients != nullptr && order > 1) {
                        (*pShapeFunctionGradients)(0, i) +=
                            derivative * gradient_x_coefficient;
                        (*pShapeFunctionGradients)(1, i) +=
                            derivative * gradient_y_coefficient;
                    }
                }
            }
        } else {
            std::size_t derivative_index = 0;
            for (int x_order = static_cast<int>(order);
                 x_order >= 0;
                 --x_order) {
                for (int y_order =
                         static_cast<int>(order) - x_order;
                     y_order >= 0;
                     --y_order) {
                    const IndexType z_order = order -
                        static_cast<IndexType>(x_order) -
                        static_cast<IndexType>(y_order);
                    const double value_coefficient =
                        scaled_powers[0][x_order] *
                        scaled_powers[1][y_order] *
                        scaled_powers[2][z_order];
                    const double gradient_x_coefficient = x_order > 0
                        ? scaled_powers[0][x_order - 1] *
                            scaled_powers[1][y_order] *
                            scaled_powers[2][z_order]
                        : 0.0;
                    const double gradient_y_coefficient = y_order > 0
                        ? scaled_powers[0][x_order] *
                            scaled_powers[1][y_order - 1] *
                            scaled_powers[2][z_order]
                        : 0.0;
                    const double gradient_z_coefficient = z_order > 0
                        ? scaled_powers[0][x_order] *
                            scaled_powers[1][y_order] *
                            scaled_powers[2][z_order - 1]
                        : 0.0;
                    for (std::size_t i = 0;
                         i < number_of_control_points;
                         ++i) {
                        const double derivative =
                            r_derivatives(i, derivative_index);
                        if (pShapeFunctions != nullptr) {
                            (*pShapeFunctions)[i] +=
                                derivative * value_coefficient;
                        }
                        if (pShapeFunctionGradients != nullptr && order > 1) {
                            (*pShapeFunctionGradients)(0, i) +=
                                derivative * gradient_x_coefficient;
                            (*pShapeFunctionGradients)(1, i) +=
                                derivative * gradient_y_coefficient;
                            (*pShapeFunctionGradients)(2, i) +=
                                derivative * gradient_z_coefficient;
                        }
                    }
                    ++derivative_index;
                }
            }
        }
    }
}

GapSbmSolidInterfaceConditionBatched::
    GapSbmSolidInterfaceConditionBatched(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : GapSbmSolidInterfaceCondition(NewId, pGeometry)
{
    CompactQuadratureGeometries();
}

GapSbmSolidInterfaceConditionBatched::
    GapSbmSolidInterfaceConditionBatched(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : GapSbmSolidInterfaceCondition(
        NewId,
        pGeometry,
        pProperties)
{
    CompactQuadratureGeometries();
}

Condition::Pointer GapSbmSolidInterfaceConditionBatched::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<
        GapSbmSolidInterfaceConditionBatched>(
            NewId,
            pGeom,
            pProperties);
}

Condition::Pointer GapSbmSolidInterfaceConditionBatched::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << "GapSbmSolidInterfaceConditionBatched #" << Id()
        << " has no quadrature-point geometry parts.\n";

    std::vector<GeometryType::Pointer> quadrature_point_geometries;
    quadrature_point_geometries.reserve(NumberOfQuadraturePoints());
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        quadrature_point_geometries.push_back(
            GetGeometry().pGetGeometryPart(point_index));
    }

    auto p_coupling_geometry =
        Kratos::make_shared<CouplingGeometry<Node>>(
            std::move(quadrature_point_geometries));

    return Kratos::make_intrusive<
        GapSbmSolidInterfaceConditionBatched>(
            NewId,
            p_coupling_geometry,
            pProperties);
}

GeometryData::IntegrationMethod
GapSbmSolidInterfaceConditionBatched::GetIntegrationMethod() const
{
    if (NumberOfQuadraturePoints() > 0) {
        return GetQuadraturePointGeometry(0)
            .GetDefaultIntegrationMethod();
    }
    return BaseType::GetIntegrationMethod();
}

void GapSbmSolidInterfaceConditionBatched::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(NumberOfQuadraturePoints() == 0)
        << "GapSbmSolidInterfaceConditionBatched #" << Id()
        << " has no quadrature points.\n";
    KRATOS_ERROR_IF_NOT(Has(NEIGHBOUR_GEOMETRIES))
        << "GapSbmSolidInterfaceConditionBatched #" << Id()
        << " has no NEIGHBOUR_GEOMETRIES.\n";
    const auto& r_neighbour_geometries =
        GetValue(NEIGHBOUR_GEOMETRIES);
    KRATOS_ERROR_IF(r_neighbour_geometries.size() != 2)
        << "GapSbmSolidInterfaceConditionBatched #" << Id()
        << " requires exactly two ordered NEIGHBOUR_GEOMETRIES, got "
        << r_neighbour_geometries.size() << ".\n";
    KRATOS_ERROR_IF_NOT(r_neighbour_geometries[0])
        << "Null plus-side NEIGHBOUR_GEOMETRIES entry in condition #"
        << Id() << ".\n";
    KRATOS_ERROR_IF_NOT(r_neighbour_geometries[1])
        << "Null minus-side NEIGHBOUR_GEOMETRIES entry in condition #"
        << Id() << ".\n";
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CONSTITUTIVE_LAW))
        << "A constitutive law is required by condition #"
        << Id() << ".\n";
    KRATOS_ERROR_IF(GetProperties()[CONSTITUTIVE_LAW] == nullptr)
        << "The constitutive law is null in condition #"
        << Id() << ".\n";

    InitializeMemberVariables(GetQuadraturePointGeometry(0));

    const std::size_t number_of_points =
        NumberOfQuadraturePoints();

    KRATOS_ERROR_IF(
        mQuadraturePointReferenceWeights.size() != number_of_points ||
        mQuadraturePointCoordinates.size1() != number_of_points ||
        mQuadraturePointCoordinates.size2() != 3 ||
        mQuadraturePointNormals.size1() != number_of_points ||
        mQuadraturePointNormals.size2() != 3)
        << "Missing compact quadrature data in "
        << "GapSbmSolidInterfaceConditionBatched #" << Id() << ".\n";

    mConstitutiveLawVector.resize(number_of_points);
    mQuadraturePointWeights.resize(number_of_points, false);

    double total_integration_weight = 0.0;
    for (std::size_t point_index = 0;
         point_index < number_of_points;
         ++point_index) {
        const auto& r_quadrature_geometry =
            GetQuadraturePointGeometry(point_index);
        double integration_weight =
            mQuadraturePointReferenceWeights[point_index];
        if (mDim == 2) {
            integration_weight *= GetProperties().Has(THICKNESS)
                ? GetProperties()[THICKNESS]
                : 1.0;
        }
        KRATOS_ERROR_IF(integration_weight < 0.0)
            << "Negative fitted interface integration weight "
            << integration_weight << " at point " << point_index
            << " of condition #" << Id()
            << ". Negative weights are not corrected.\n";
        KRATOS_ERROR_IF(
            !std::isfinite(integration_weight) ||
            integration_weight <= 0.0)
            << "Invalid integration weight " << integration_weight
            << " at point " << point_index << " of condition #"
            << Id() << ".\n";
        mQuadraturePointWeights[point_index] = integration_weight;
        total_integration_weight += integration_weight;

        mConstitutiveLawVector[point_index] =
            GetProperties()[CONSTITUTIVE_LAW]->Clone();
        const auto& r_shape_function_values =
            r_quadrature_geometry.ShapeFunctionsValues(
                GetIntegrationMethod());
        mConstitutiveLawVector[point_index]->InitializeMaterial(
            GetProperties(),
            r_quadrature_geometry,
            row(r_shape_function_values, 0));
    }

    SetValue(INTEGRATION_WEIGHT, total_integration_weight);
    mpConstitutiveLaw = mConstitutiveLawVector.front();

    KRATOS_CATCH("")
}

void GapSbmSolidInterfaceConditionBatched::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        &rLeftHandSideMatrix,
        &rRightHandSideVector,
        rCurrentProcessInfo);
}

void GapSbmSolidInterfaceConditionBatched::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        &rLeftHandSideMatrix,
        nullptr,
        rCurrentProcessInfo);
}

void GapSbmSolidInterfaceConditionBatched::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateAllContributions(
        nullptr,
        &rRightHandSideVector,
        rCurrentProcessInfo);
}

void GapSbmSolidInterfaceConditionBatched::CalculateAllContributions(
    MatrixType* pLeftHandSideMatrix,
    VectorType* pRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_neighbour_geometries =
        GetValue(NEIGHBOUR_GEOMETRIES);
    KRATOS_ERROR_IF(r_neighbour_geometries.size() != 2)
        << "GapSbmSolidInterfaceConditionBatched #" << Id()
        << " requires exactly two ordered neighbour geometries.\n";
    const auto& r_geometry_plus = *r_neighbour_geometries[0];
    const auto& r_geometry_minus = *r_neighbour_geometries[1];
    const std::size_t number_of_control_points_plus =
        r_geometry_plus.PointsNumber();
    const std::size_t number_of_control_points_minus =
        r_geometry_minus.PointsNumber();
    const std::size_t matrix_size_plus =
        number_of_control_points_plus * mDim;
    const std::size_t matrix_size_minus =
        number_of_control_points_minus * mDim;
    const std::size_t matrix_size =
        matrix_size_plus + matrix_size_minus;

    if (pLeftHandSideMatrix != nullptr) {
        if (pLeftHandSideMatrix->size1() != matrix_size ||
            pLeftHandSideMatrix->size2() != matrix_size) {
            pLeftHandSideMatrix->resize(
                matrix_size,
                matrix_size,
                false);
        }
        noalias(*pLeftHandSideMatrix) =
            ZeroMatrix(matrix_size, matrix_size);
    }
    if (pRightHandSideVector != nullptr) {
        if (pRightHandSideVector->size() != matrix_size) {
            pRightHandSideVector->resize(matrix_size, false);
        }
        noalias(*pRightHandSideVector) = ZeroVector(matrix_size);
    }

    KRATOS_ERROR_IF(
        mConstitutiveLawVector.size() != NumberOfQuadraturePoints() ||
        mQuadraturePointWeights.size() != NumberOfQuadraturePoints() ||
        mQuadraturePointNormals.size1() != NumberOfQuadraturePoints() ||
        mQuadraturePointNormals.size2() != 3)
        << "Inconsistent quadrature cache in "
        << "GapSbmSolidInterfaceConditionBatched #" << Id()
        << ".\n";

    Vector displacement_plus(matrix_size_plus);
    Vector displacement_minus(matrix_size_minus);
    GetSolutionCoefficientVectorPlus(displacement_plus);
    GetSolutionCoefficientVectorMinus(displacement_minus);

    KRATOS_ERROR_IF(mConstitutiveLawVector.empty())
        << "No constitutive laws in condition #" << Id() << ".\n";
    const std::size_t strain_size =
        mConstitutiveLawVector.front()->GetStrainSize();
    Matrix B_plus(strain_size, matrix_size_plus);
    Matrix B_minus(strain_size, matrix_size_minus);
    Matrix traction_tangent(mDim, matrix_size);
    Vector signed_shape_function(matrix_size);
    Vector stress_column(strain_size);
    Vector traction_column(3);
    Vector traction_plus(3);
    Vector traction_minus(3);
    Vector displacement_jump(3);
    Vector shape_functions_plus(number_of_control_points_plus);
    Vector shape_functions_minus(number_of_control_points_minus);
    Matrix shape_function_gradients_plus(
        3,
        number_of_control_points_plus);
    Matrix shape_function_gradients_minus(
        3,
        number_of_control_points_minus);
    const TaylorDerivativeData taylor_derivatives_plus =
        CalculateTaylorDerivativeData(r_geometry_plus);
    const TaylorDerivativeData taylor_derivatives_minus =
        CalculateTaylorDerivativeData(r_geometry_minus);

    ConstitutiveVariables variables_plus(strain_size);
    ConstitutiveVariables variables_minus(strain_size);

    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        auto& p_constitutive_law =
            mConstitutiveLawVector[point_index];
        KRATOS_ERROR_IF_NOT(p_constitutive_law)
            << "Null constitutive law at point " << point_index
            << " of condition #" << Id() << ".\n";
        KRATOS_ERROR_IF(
            p_constitutive_law->GetStrainSize() != strain_size)
            << "Inconsistent strain size at point " << point_index
            << " of condition #" << Id() << ".\n";

        SetQuadraturePointDistances(point_index);

        EvaluateTaylorExpansion(
            r_geometry_plus,
            mDistanceVectorPlus,
            taylor_derivatives_plus,
            &shape_functions_plus,
            &shape_function_gradients_plus);
        EvaluateTaylorExpansion(
            r_geometry_minus,
            mDistanceVectorMinus,
            taylor_derivatives_minus,
            &shape_functions_minus,
            &shape_function_gradients_minus);

        FillB(
            B_plus,
            shape_function_gradients_plus,
            number_of_control_points_plus);
        FillB(
            B_minus,
            shape_function_gradients_minus,
            number_of_control_points_minus);

        for (std::size_t strain_index = 0;
             strain_index < strain_size;
             ++strain_index) {
            double strain_plus = 0.0;
            for (std::size_t column = 0;
                 column < matrix_size_plus;
                 ++column) {
                strain_plus +=
                    B_plus(strain_index, column) *
                    displacement_plus[column];
            }
            variables_plus.StrainVector[strain_index] = strain_plus;
            variables_plus.StressVector[strain_index] = 0.0;

            double strain_minus = 0.0;
            for (std::size_t column = 0;
                 column < matrix_size_minus;
                 ++column) {
                strain_minus +=
                    B_minus(strain_index, column) *
                    displacement_minus[column];
            }
            variables_minus.StrainVector[strain_index] = strain_minus;
            variables_minus.StressVector[strain_index] = 0.0;
        }
        noalias(variables_plus.D) =
            ZeroMatrix(strain_size, strain_size);
        noalias(variables_minus.D) =
            ZeroMatrix(strain_size, strain_size);

        const auto& r_quadrature_geometry =
            GetQuadraturePointGeometry(point_index);
        ConstitutiveLaw::Parameters parameters_plus(
            r_quadrature_geometry,
            GetProperties(),
            rCurrentProcessInfo);
        Flags& r_options_plus = parameters_plus.GetOptions();
        r_options_plus.Set(
            ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN,
            true);
        r_options_plus.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options_plus.Set(
            ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
            true);
        parameters_plus.SetStrainVector(
            variables_plus.StrainVector);
        parameters_plus.SetStressVector(
            variables_plus.StressVector);
        parameters_plus.SetConstitutiveMatrix(variables_plus.D);
        p_constitutive_law->CalculateMaterialResponse(
            parameters_plus,
            ConstitutiveLaw::StressMeasure_Cauchy);

        ConstitutiveLaw::Parameters parameters_minus(
            r_quadrature_geometry,
            GetProperties(),
            rCurrentProcessInfo);
        Flags& r_options_minus = parameters_minus.GetOptions();
        r_options_minus.Set(
            ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN,
            true);
        r_options_minus.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options_minus.Set(
            ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
            true);
        parameters_minus.SetStrainVector(
            variables_minus.StrainVector);
        parameters_minus.SetStressVector(
            variables_minus.StressVector);
        parameters_minus.SetConstitutiveMatrix(variables_minus.D);
        p_constitutive_law->CalculateMaterialResponse(
            parameters_minus,
            ConstitutiveLaw::StressMeasure_Cauchy);

        array_1d<double, 3> normal = ZeroVector(3);
        for (std::size_t component = 0; component < 3; ++component) {
            normal[component] =
                mQuadraturePointNormals(point_index, component);
        }
        CalculateTraction(
            variables_plus.StressVector,
            normal,
            traction_plus);
        CalculateTraction(
            variables_minus.StressVector,
            normal,
            traction_minus);

        for (std::size_t column = 0;
             column < matrix_size_plus;
             ++column) {
            for (std::size_t row = 0; row < strain_size; ++row) {
                double value = 0.0;
                for (std::size_t inner = 0;
                     inner < strain_size;
                     ++inner) {
                    value += variables_plus.D(row, inner) *
                        B_plus(inner, column);
                }
                stress_column[row] = value;
            }
            CalculateTraction(stress_column, normal, traction_column);
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                traction_tangent(component, column) =
                    traction_column[component];
            }
        }
        for (std::size_t column = 0;
             column < matrix_size_minus;
             ++column) {
            for (std::size_t row = 0; row < strain_size; ++row) {
                double value = 0.0;
                for (std::size_t inner = 0;
                     inner < strain_size;
                     ++inner) {
                    value += variables_minus.D(row, inner) *
                        B_minus(inner, column);
                }
                stress_column[row] = value;
            }
            CalculateTraction(stress_column, normal, traction_column);
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                traction_tangent(
                    component,
                    matrix_size_plus + column) =
                    traction_column[component];
            }
        }

        noalias(displacement_jump) = ZeroVector(3);
        for (std::size_t i = 0;
            i < number_of_control_points_plus;
             ++i) {
            const double shape_function =
                shape_functions_plus[i];
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                const std::size_t index = i * mDim + component;
                signed_shape_function[index] = shape_function;
                displacement_jump[component] +=
                    shape_function * displacement_plus[index];
            }
        }
        for (std::size_t i = 0;
            i < number_of_control_points_minus;
             ++i) {
            const double shape_function =
                shape_functions_minus[i];
            for (std::size_t component = 0;
                 component < mDim;
                 ++component) {
                const std::size_t local_index =
                    i * mDim + component;
                const std::size_t global_index =
                    matrix_size_plus + local_index;
                signed_shape_function[global_index] =
                    -shape_function;
                displacement_jump[component] -=
                    shape_function * displacement_minus[local_index];
            }
        }

        const double integration_weight =
            mQuadraturePointWeights[point_index];

        if (pLeftHandSideMatrix != nullptr) {
            for (std::size_t row = 0; row < matrix_size; ++row) {
                const std::size_t row_component = row % mDim;
                const double signed_shape_row =
                    signed_shape_function[row];

                for (std::size_t column = 0;
                     column < matrix_size;
                     ++column) {
                    const std::size_t column_component =
                        column % mDim;
                    double contribution =
                        -0.5 * signed_shape_row *
                            traction_tangent(
                                row_component,
                                column)
                        -0.5 * mNitschePenalty *
                            traction_tangent(
                                column_component,
                                row) *
                            signed_shape_function[column];
                    if (row_component == column_component) {
                        contribution +=
                            mPenalty * signed_shape_row *
                            signed_shape_function[column];
                    }
                    (*pLeftHandSideMatrix)(row, column) +=
                        integration_weight * contribution;
                }
            }
        }

        if (pRightHandSideVector != nullptr) {
            for (std::size_t row = 0; row < matrix_size; ++row) {
                const std::size_t row_component = row % mDim;
                double contribution =
                    signed_shape_function[row] *
                    (0.5 *
                        (traction_plus[row_component] +
                         traction_minus[row_component])
                     - mPenalty *
                        displacement_jump[row_component]);
                for (std::size_t component = 0;
                     component < mDim;
                     ++component) {
                    contribution +=
                        0.5 * mNitschePenalty *
                        traction_tangent(component, row) *
                        displacement_jump[component];
                }
                (*pRightHandSideVector)[row] +=
                    integration_weight * contribution;
            }
        }
    }

    KRATOS_CATCH("")
}

void GapSbmSolidInterfaceConditionBatched::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        ConstitutiveLaw::Parameters constitutive_parameters(
            GetQuadraturePointGeometry(point_index),
            GetProperties(),
            rCurrentProcessInfo);
        mConstitutiveLawVector[point_index]
            ->InitializeMaterialResponse(
                constitutive_parameters,
                ConstitutiveLaw::StressMeasure_Cauchy);
    }
}

void GapSbmSolidInterfaceConditionBatched::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    for (std::size_t point_index = 0;
         point_index < NumberOfQuadraturePoints();
         ++point_index) {
        ConstitutiveLaw::Parameters constitutive_parameters(
            GetQuadraturePointGeometry(point_index),
            GetProperties(),
            rCurrentProcessInfo);
        mConstitutiveLawVector[point_index]
            ->FinalizeMaterialResponse(
                constitutive_parameters,
                ConstitutiveLaw::StressMeasure_Cauchy);
    }
}

std::size_t
GapSbmSolidInterfaceConditionBatched::NumberOfQuadraturePoints() const
{
    return GetGeometry().NumberOfGeometryParts();
}

const GapSbmSolidInterfaceConditionBatched::GeometryType&
GapSbmSolidInterfaceConditionBatched::GetQuadraturePointGeometry(
    const std::size_t PointIndex) const
{
    KRATOS_ERROR_IF(PointIndex >= NumberOfQuadraturePoints())
        << "Quadrature point index " << PointIndex
        << " is out of range [0, " << NumberOfQuadraturePoints()
        << ") in GapSbmSolidInterfaceConditionBatched #"
        << Id() << ".\n";
    return GetGeometry().GetGeometryPart(PointIndex);
}

void GapSbmSolidInterfaceConditionBatched::
    CompactQuadratureGeometries()
{
    const std::size_t number_of_points =
        GetGeometry().NumberOfGeometryParts();
    if (number_of_points == 0 ||
        mQuadraturePointCoordinates.size1() == number_of_points) {
        return;
    }

    mQuadraturePointReferenceWeights.resize(number_of_points, false);
    mQuadraturePointCoordinates.resize(number_of_points, 3, false);
    mQuadraturePointNormals.resize(number_of_points, 3, false);
    auto p_representative_geometry =
        GetGeometry().pGetGeometryPart(0);

    for (std::size_t point_index = 0;
         point_index < number_of_points;
         ++point_index) {
        auto p_quadrature_geometry =
            GetGeometry().pGetGeometryPart(point_index);
        const auto integration_method =
            p_quadrature_geometry->GetDefaultIntegrationMethod();
        const auto& r_integration_points =
            p_quadrature_geometry->IntegrationPoints(integration_method);
        KRATOS_ERROR_IF(r_integration_points.size() != 1)
            << "Each geometry part of "
            << "GapSbmSolidInterfaceConditionBatched #" << Id()
            << " must contain exactly one integration point.\n";
        mQuadraturePointReferenceWeights[point_index] =
            p_quadrature_geometry->Has(INTEGRATION_WEIGHT)
                ? p_quadrature_geometry->GetValue(INTEGRATION_WEIGHT)
                : r_integration_points.front().Weight();
        const auto center = p_quadrature_geometry->Center();
        for (std::size_t component = 0; component < 3; ++component) {
            mQuadraturePointCoordinates(point_index, component) =
                center[component];
        }
        auto normal = p_quadrature_geometry->Normal(
            0,
            integration_method);
        const double normal_norm = MathUtils<double>::Norm(normal);
        KRATOS_ERROR_IF(normal_norm <= 0.0)
            << "Zero interface normal at point " << point_index
            << " of condition #" << Id() << ".\n";
        normal /= normal_norm;
        for (std::size_t component = 0; component < 3; ++component) {
            mQuadraturePointNormals(point_index, component) =
                normal[component];
        }
        if (point_index > 0) {
            GetGeometry().SetGeometryPart(
                point_index,
                p_representative_geometry);
        }
    }
}

void GapSbmSolidInterfaceConditionBatched::
    SetQuadraturePointDistances(const std::size_t PointIndex)
{
    const auto& r_neighbour_geometries =
        GetValue(NEIGHBOUR_GEOMETRIES);
    const auto center_plus = r_neighbour_geometries[0]->Center();
    const auto center_minus = r_neighbour_geometries[1]->Center();
    if (mDistanceVectorPlus.size() != 3) {
        mDistanceVectorPlus.resize(3, false);
        mDistanceVectorMinus.resize(3, false);
    }
    for (std::size_t component = 0; component < 3; ++component) {
        mDistanceVectorPlus[component] =
            mQuadraturePointCoordinates(PointIndex, component) -
            center_plus[component];
        mDistanceVectorMinus[component] =
            mQuadraturePointCoordinates(PointIndex, component) -
            center_minus[component];
    }
}

void GapSbmSolidInterfaceConditionBatched::FillB(
    Matrix& rB,
    const Matrix& rShapeFunctionGradients,
    const std::size_t NumberOfControlPoints) const
{
    const std::size_t matrix_size = NumberOfControlPoints * mDim;
    const std::size_t strain_size = mDim == 2 ? 3 : 6;
    KRATOS_ERROR_IF(
        rShapeFunctionGradients.size1() != 3 ||
        rShapeFunctionGradients.size2() != NumberOfControlPoints)
        << "Invalid Taylor gradients in condition #" << Id()
        << ".\n";

    if (rB.size1() != strain_size ||
        rB.size2() != matrix_size) {
        rB.resize(strain_size, matrix_size, false);
    }
    noalias(rB) = ZeroMatrix(strain_size, matrix_size);

    for (std::size_t i = 0; i < NumberOfControlPoints; ++i) {
        const std::size_t index = i * mDim;
        const double derivative_x =
            rShapeFunctionGradients(0, i);
        const double derivative_y =
            rShapeFunctionGradients(1, i);
        if (mDim == 2) {
            rB(0, index) = derivative_x;
            rB(1, index + 1) = derivative_y;
            rB(2, index) = derivative_y;
            rB(2, index + 1) = derivative_x;
        } else {
            const double derivative_z =
                rShapeFunctionGradients(2, i);
            rB(0, index) = derivative_x;
            rB(1, index + 1) = derivative_y;
            rB(2, index + 2) = derivative_z;
            rB(3, index) = derivative_y;
            rB(3, index + 1) = derivative_x;
            rB(4, index + 1) = derivative_z;
            rB(4, index + 2) = derivative_y;
            rB(5, index) = derivative_z;
            rB(5, index + 2) = derivative_x;
        }
    }
}

} // Namespace Kratos
