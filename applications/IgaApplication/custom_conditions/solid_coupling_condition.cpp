//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi

// System includes
#include <limits>

// External includes

// Project includes
#include "custom_conditions/solid_coupling_condition.h"
#include "includes/constitutive_law.h"
#include "includes/variables.h"
#include "iga_application_variables.h"

namespace Kratos
{

void SolidCouplingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    
    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;
    const double integration_weight = r_integration_points[0].Weight()*thickness;

    SetValue(INTEGRATION_WEIGHT, integration_weight);

    if (mIsGapSbmCoupling) {
        const auto& r_true = GetGeometry();
        const auto& r_surrogate_B = GetGeometryMirror();

        mDistanceVectorB.resize(3);
        noalias(mDistanceVectorB) = r_true.Center().Coordinates() - r_surrogate_B.Center().Coordinates();
    }
}

void SolidCouplingCondition::InitializeMaterial()
{
    KRATOS_TRY
    if (GetProperties().Has(CONSTITUTIVE_LAW) && GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial(r_properties, r_geometry, row(N_values, 0));
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the condition with ID " << this->Id() << std::endl;
    }
    KRATOS_CATCH("")
}

void SolidCouplingCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(
        constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void SolidCouplingCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->FinalizeMaterialResponse(
        constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}

void SolidCouplingCondition::InitializeMemberVariables()
{
    const auto& r_geometry_patchA = GetGeometry();
    const auto& r_geometry_patchB = GetGeometryMirror();

    // Basis function order of the Taylor expansion is decided by the r_geometry_patchB
    const auto& r_DN_De_B = r_geometry_patchB.ShapeFunctionsLocalGradients(r_geometry_patchB.GetDefaultIntegrationMethod());
    KRATOS_ERROR_IF(r_DN_De_B.empty()) << "SolidCouplingCondition requires gradients on the surrogate patch." << std::endl;
    mDim = r_DN_De_B[0].size2();
    KRATOS_ERROR_IF(mDim != 2) << "SolidCouplingCondition momentarily only supports 2D conditions, but the current dimension is " << mDim << std::endl;
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De_B[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De_B[0].size1()) - 1;
    }
    mBasisFunctionsOrder *= 2;

    const auto& r_DN_De_A = r_geometry_patchA.ShapeFunctionsLocalGradients(r_geometry_patchA.GetDefaultIntegrationMethod());
    KRATOS_ERROR_IF(r_DN_De_A.empty()) << "SolidCouplingCondition requires gradients on the true patch." << std::endl;
    IndexType basis_order_A = 0;
    if (mDim == 3) {
        basis_order_A = std::cbrt(r_DN_De_A[0].size1()) - 1;
    } else {
        basis_order_A = std::sqrt(r_DN_De_A[0].size1()) - 1;
    }
    basis_order_A *= 2;
    mBasisFunctionsOrder = std::min(mBasisFunctionsOrder, basis_order_A);

    // if (norm_2(r_geometry_patchA.Center()-r_geometry_patchB.Center()) > 1e-12)
        mIsGapSbmCoupling = true;

    mNormalParameterSpaceA = -r_geometry_patchA.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceA /= MathUtils<double>::Norm(mNormalParameterSpaceA);
    mNormalPhysicalSpaceA = mNormalParameterSpaceA;

    mNormalParameterSpaceB = -r_geometry_patchB.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceB /= MathUtils<double>::Norm(mNormalParameterSpaceB);
    mNormalPhysicalSpaceB = mNormalParameterSpaceB;

    KRATOS_ERROR_IF(std::abs(inner_prod(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB)+1) > 1e-12 && !mIsGapSbmCoupling)
        << "SolidCouplingCondition: normals are not opposite." << std::endl;
    
    SetValue(NORMAL, mNormalPhysicalSpaceA);
    SetValue(PROJECTION_NODE_COORDINATES, r_geometry_patchB.Center());
}

void SolidCouplingCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{   
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const std::size_t n_patch_A = r_patch_A.PointsNumber();
    const std::size_t n_patch_B = r_patch_B.PointsNumber();
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const std::size_t dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const std::size_t dofs_per_node = dim;
    const std::size_t total_size = dofs_per_node * (n_patch_A + n_patch_B);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != total_size)
        rLeftHandSideMatrix.resize(total_size,total_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size,total_size); //resetting LHS

    if (rRightHandSideVector.size() != total_size) {
        rRightHandSideVector.resize(total_size, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(total_size);

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void SolidCouplingCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const std::size_t n_patch_A = r_patch_A.PointsNumber();
    const std::size_t n_patch_B = r_patch_B.PointsNumber();
    const auto& rDN_De_A_dim = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const std::size_t dim = rDN_De_A_dim.empty() ? 2 : rDN_De_A_dim[0].size2();
    const std::size_t dofs_per_node = dim;

    KRATOS_ERROR_IF(dim != 2) << "SolidCouplingCondition currently only supports 2D." << std::endl;

    const std::size_t mat_size_A = dofs_per_node * n_patch_A;
    const std::size_t mat_size_B = dofs_per_node * n_patch_B;

    const std::size_t total_size = dofs_per_node * (n_patch_A + n_patch_B);

    KRATOS_ERROR_IF(total_size == 0) << "SolidCouplingCondition found empty geometry." << std::endl;

    if (rLeftHandSideMatrix.size1() != total_size || rLeftHandSideMatrix.size2() != total_size) {
        rLeftHandSideMatrix.resize(total_size, total_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);

    const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    const auto& integration_points_patch_B = r_patch_B.IntegrationPoints();

    KRATOS_ERROR_IF(integration_points_patch_A.empty() || integration_points_patch_B.empty())
        << "SolidCouplingCondition expects at least one integration point." << std::endl;

    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_B =
        r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());

    KRATOS_ERROR_IF(rDN_De_A.size() == 0 || rDN_De_B.size() == 0)
        << "Shape function gradients not available for SolidCouplingCondition. Configure the modeler with derivative order >= 2." << std::endl;

    KRATOS_ERROR_IF(rDN_De_A[0].size1() != n_patch_A)
        << "SolidCouplingCondition: gradient matrix rows mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << rDN_De_A[0].size1() << "." << std::endl;
    KRATOS_ERROR_IF(rDN_De_B[0].size1() != n_patch_B)
        << "SolidCouplingCondition: gradient matrix rows mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << rDN_De_B[0].size1() << "." << std::endl;

    const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    Matrix N_patch_B = r_patch_B.ShapeFunctionsValues();

    if (mIsGapSbmCoupling) {
        Vector N_sum_B_vector(n_patch_B);
        ComputeTaylorExpansionContribution(r_patch_B, mDistanceVectorB, N_sum_B_vector);
        for (IndexType j = 0; j < n_patch_B; ++j) {
            N_patch_B(0, j) = N_sum_B_vector[j];
        }
    }

    KRATOS_ERROR_IF(N_patch_A.size2() != n_patch_A)
        << "SolidCouplingCondition: shape function values columns mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << N_patch_A.size2() << "." << std::endl;
    KRATOS_ERROR_IF(N_patch_B.size2() != n_patch_B)
        << "SolidCouplingCondition: shape function values columns mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << N_patch_B.size2() << "." << std::endl;

    const unsigned int dim_patch_A = rDN_De_A[0].size2();
    const unsigned int dim_patch_B = rDN_De_B[0].size2();

    KRATOS_ERROR_IF(dim_patch_A < 2 || dim_patch_B < 2)
        << "SolidCouplingCondition requires at least 2D parameter space." << std::endl;

    // Jacobians to compute the interface measure
    GeometryType::JacobiansType J_patch_A;
    r_patch_A.Jacobian(J_patch_A, r_patch_A.GetDefaultIntegrationMethod());
    Matrix jacobian_patch_A = ZeroMatrix(3, 3);
    jacobian_patch_A(0, 0) = J_patch_A[0](0, 0);
    jacobian_patch_A(0, 1) = J_patch_A[0](0, 1);
    jacobian_patch_A(1, 0) = J_patch_A[0](1, 0);
    jacobian_patch_A(1, 1) = J_patch_A[0](1, 1);
    jacobian_patch_A(2, 2) = 1.0;
    array_1d<double, 3> tangent_patch_A;
    r_patch_A.Calculate(LOCAL_TANGENT, tangent_patch_A);
    Vector determinant_factor_patch_A = prod(jacobian_patch_A, tangent_patch_A);
    determinant_factor_patch_A[2] = 0.0;
    const double detJ_patch_A = norm_2(determinant_factor_patch_A);
    KRATOS_ERROR_IF(detJ_patch_A <= std::numeric_limits<double>::epsilon())
        << "SolidCouplingCondition: degenerate Jacobian for patch A (determinant ~ 0)." << std::endl;

    const double integration_weight = GetValue(INTEGRATION_WEIGHT)* detJ_patch_A; //FIXME:

    KRATOS_ERROR_IF(std::abs(integration_weight) <= std::numeric_limits<double>::epsilon())
        << "SolidCouplingCondition: zero integration integration_weight on patch A." << std::endl;

    double penalty_factor = 0;
    double nitsche_penalty_free = 1.0;
    if (GetProperties().Has(PENALTY_FACTOR)) {
        penalty_factor = GetProperties()[PENALTY_FACTOR];
        // penalty_factor = 1000; 
        if (penalty_factor <= 0.0) {
            penalty_factor = 0.0;
            nitsche_penalty_free = -1.0;
        }
    } else {
        KRATOS_ERROR << "SolidCouplingCondition requires PENALTY_FACTOR to be defined in the Properties." << std::endl;
    }

    double characteristic_length;
    if (Has(KNOT_SPAN_SIZES)) {
        const Vector& knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
        if (!knot_span_sizes.empty()) {
            characteristic_length = knot_span_sizes[0];
            if (knot_span_sizes.size() > 1) {
                characteristic_length = std::min(characteristic_length, knot_span_sizes[1]);
            }
        } else {
            KRATOS_ERROR << "KNOT_SPAN_SIZES is empty in SolidCouplingCondition" << std::endl;
        }
    } else {
        KRATOS_ERROR << "KNOT_SPAN_SIZES is not defined in SolidCouplingCondition" << std::endl;
    }

    penalty_factor = penalty_factor / characteristic_length;

    Matrix DN_patch_A = rDN_De_A[0];
    Matrix DN_patch_B = rDN_De_B[0];

    // compute Taylor expansion contribution: H_sum_vec (only for gap-sbm)
    Vector N_vec_A = row(N_patch_A, 0);
    Vector N_sum_vec_B = row(N_patch_B, 0);
    if (mIsGapSbmCoupling) {
        N_sum_vec_B = ZeroVector(n_patch_B);
        ComputeTaylorExpansionContribution(r_patch_B, mDistanceVectorB, N_sum_vec_B);
    }

    // compute Taylor expansion contribution: grad_H_sum (only for gap-sbm)
    Matrix grad_N_A = DN_patch_A;
    Matrix grad_N_sum_B = DN_patch_B;
    if (mIsGapSbmCoupling) {
        Matrix grad_N_sum_transposed_B = ZeroMatrix(3, n_patch_B);
        ComputeGradientTaylorExpansionContribution(r_patch_B, mDistanceVectorB, grad_N_sum_transposed_B);
        grad_N_sum_B = trans(grad_N_sum_transposed_B);
    }


    // compute the B matrix
    Matrix B_A = ZeroMatrix(mDim, mat_size_A);
    CalculateB(r_patch_A, B_A, grad_N_A);

    Matrix B_sum_B = ZeroMatrix(mDim, mat_size_B);
    CalculateB(r_patch_B, B_sum_B, grad_N_sum_B);

    // obtain the tangent constitutive matrix for the body fitted side
    ConstitutiveLaw::Parameters values_A(r_patch_A, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector_A(mat_size_A);
    GetSolutionCoefficientVectorA(old_displacement_coefficient_vector_A);
    Vector old_strain_A = prod(B_A, old_displacement_coefficient_vector_A);
    const std::size_t strain_size_A = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_A(strain_size_A);
    ApplyConstitutiveLaw(mat_size_A, old_strain_A, values_A, this_constitutive_variables_A);

    const Matrix& r_D_A = values_A.GetConstitutiveMatrix();

    // obtain the tangent constitutive matrix at the true position for the minus side
    ConstitutiveLaw::Parameters values_true_B(r_patch_B, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector_B(mat_size_B);
    GetSolutionCoefficientVectorB(old_displacement_coefficient_vector_B);
    Vector old_strain_on_true_B = prod(B_sum_B, old_displacement_coefficient_vector_B);
    const std::size_t strain_size_true_B = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true_B(strain_size_true_B);
    ApplyConstitutiveLaw(mat_size_B, old_strain_on_true_B, values_true_B, this_constitutive_variables_true_B);

    const Matrix& r_D_on_true_B = values_true_B.GetConstitutiveMatrix();

    // compute the DB product
    Matrix DB_A = prod(r_D_A, B_A);
    Matrix DB_sum_B = prod(r_D_on_true_B, B_sum_B);

    // // ASSEMBLE
    // //-----------------------------------------------------
    const std::size_t shift_dof = mat_size_A;
    // -w_A * (sigma_A + sigma_B) 
    for (IndexType i = 0; i < n_patch_A; i++) {
        // FIRST TERM: -w_A * sigma_A 
        for (IndexType j = 0; j < n_patch_A; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 2*i+idim;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*N_vec_A(i)
                                                    * (DB_A(id1, jglob)* mNormalPhysicalSpaceA[0] + DB_A(id2, jglob)* mNormalPhysicalSpaceA[1]) 
                                                    * integration_weight;
                }
            }
        }

        // SECOND TERM: -w_A * sigma_B 
        for (IndexType j = 0; j < n_patch_B; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 2*i+idim;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jloc = 2*j+jdim;
                    const int jglob = 2*j+jdim + shift_dof;

                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*N_vec_A(i)
                                                    * (DB_sum_B(id1, jloc)* mNormalPhysicalSpaceA[0] + DB_sum_B(id2, jloc)* mNormalPhysicalSpaceA[1]) 
                                                    * integration_weight;
                }
            }
        }
    }

    // +w_B * (sigma_A + sigma_B) 
    for (IndexType i = 0; i < n_patch_B; i++) {
        // FIRST TERM: +w_B * sigma_A 
        for (IndexType j = 0; j < n_patch_A; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 2*i+idim + shift_dof;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;

                    rLeftHandSideMatrix(iglob, jglob) += 0.5*N_sum_vec_B(i)
                                                    * (DB_A(id1, jglob)* mNormalPhysicalSpaceA[0] + DB_A(id2, jglob)* mNormalPhysicalSpaceA[1]) 
                                                    * integration_weight;

                }
            }
        }

        // SECOND TERM: +w_B * sigma_B 
        for (IndexType j = 0; j < n_patch_B; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 2*i+idim + shift_dof;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jloc = 2*j+jdim;
                    const int jglob = 2*j+jdim + shift_dof;

                    rLeftHandSideMatrix(iglob, jglob) += 0.5*N_sum_vec_B(i)
                                                    * (DB_sum_B(id1, jloc)* mNormalPhysicalSpaceA[0] + DB_sum_B(id2, jloc)* mNormalPhysicalSpaceA[1]) 
                                                    * integration_weight;
                }
            }
        }
    }

    // // // CONTINUITY TERMS
    // // // Differential area
    double penalty_integration = penalty_factor * integration_weight;

    // // ASSEMBLE
    // //-----------------------------------------------------
    for (IndexType i = 0; i < n_patch_A; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            const int id1 = 2*idim;
            const int iglob = 2*i+idim;

            Vector Cut_sigma_w_n_A = ZeroVector(3);
            Cut_sigma_w_n_A[0] = (DB_A(0, iglob)* mNormalPhysicalSpaceA[0] + DB_A(2, iglob)* mNormalPhysicalSpaceA[1]);
            Cut_sigma_w_n_A[1] = (DB_A(2, iglob)* mNormalPhysicalSpaceA[0] + DB_A(1, iglob)* mNormalPhysicalSpaceA[1]);

            for (IndexType j = 0; j < n_patch_A; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 2*j+idim) += N_vec_A(i)*N_vec_A(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*nitsche_penalty_free*N_vec_A(j)*Cut_sigma_w_n_A[jdim] * integration_weight;
                }

            }

            // minus side
            for (IndexType j = 0; j < n_patch_B; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 2*j+shift_dof+idim) -= N_vec_A(i)*N_sum_vec_B(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim + shift_dof;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += 0.5*nitsche_penalty_free*N_sum_vec_B(j)*Cut_sigma_w_n_A[jdim] * integration_weight;
                }
            }
        }
    }

    // minus side
    for (IndexType i = 0; i < n_patch_B; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            const int iloc = 2*i+idim;
            const int iglob = 2*i+idim + shift_dof;

            Vector Cut_sigma_w_n_B = ZeroVector(3);
            Cut_sigma_w_n_B[0] = (DB_sum_B(0, iloc)* mNormalPhysicalSpaceA[0] + DB_sum_B(2, iloc)* mNormalPhysicalSpaceA[1]);
            Cut_sigma_w_n_B[1] = (DB_sum_B(2, iloc)* mNormalPhysicalSpaceA[0] + DB_sum_B(1, iloc)* mNormalPhysicalSpaceA[1]);

            for (IndexType j = 0; j < n_patch_A; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 2*j+idim) -= N_sum_vec_B(i)*N_vec_A(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int jglob = 2*j+jdim;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*nitsche_penalty_free*N_vec_A(j)*Cut_sigma_w_n_B[jdim] * integration_weight;
                }

            }

            // minus side
            for (IndexType j = 0; j < n_patch_B; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 2*j+shift_dof+idim) += N_sum_vec_B(i)*N_sum_vec_B(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int jglob = 2*j+jdim + shift_dof;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += 0.5*nitsche_penalty_free*N_sum_vec_B(j)*Cut_sigma_w_n_B[jdim] * integration_weight;
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void SolidCouplingCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType lhs;
    CalculateLeftHandSide(lhs, rCurrentProcessInfo);

    const std::size_t size = lhs.size1();
    KRATOS_ERROR_IF(size == 0) << "SolidCouplingCondition::CalculateRightHandSide: The left hand side matrix has zero size." << std::endl;

    Vector plus_values;
    GetSolutionCoefficientVectorA(plus_values);
    Vector minus_values;
    GetSolutionCoefficientVectorB(minus_values);

    Vector values(plus_values.size() + minus_values.size());
    for (IndexType i = 0; i < plus_values.size(); ++i) {
        values[i] = plus_values[i];
    }
    for (IndexType i = 0; i < minus_values.size(); ++i) {
        values[plus_values.size() + i] = minus_values[i];
    }

    if (rRightHandSideVector.size() != size) {
        rRightHandSideVector.resize(size, false);
    }
    noalias(rRightHandSideVector) = -prod(lhs, values);
}

void SolidCouplingCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometryA = GetGeometry();
    const auto& r_geometryB = GetGeometryMirror();

    const std::size_t number_of_control_points_A = r_geometryA.size();
    const std::size_t number_of_control_points_B = r_geometryB.size();

    const std::size_t number_of_control_points = number_of_control_points_A + number_of_control_points_B;

    const std::size_t shift_dof = number_of_control_points_A * 2;

    if (rResult.size() != 2 * number_of_control_points)
        rResult.resize(2 * number_of_control_points, false);
    
    // first the plus geometry
    for (IndexType i = 0; i < number_of_control_points_A; ++i) {
        const IndexType index = i * 2;
        const auto& r_node = r_geometryA[i];
        rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
    }

    // then the minus geometry
    for (IndexType i = 0; i < number_of_control_points_B; ++i) {
        const IndexType index = i * 2 + shift_dof;
        const auto& r_node = r_geometryB[i];
        rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
    }
}

void SolidCouplingCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();

    const std::size_t number_of_control_points_A = r_patch_A.size();
    const std::size_t number_of_control_points_B = r_patch_B.size();

    const std::size_t number_of_control_points = number_of_control_points_A + number_of_control_points_B;

    rElementalDofList.resize(0);
    rElementalDofList.reserve(2 * number_of_control_points);

    // first the plus geometry
    for (IndexType i = 0; i < number_of_control_points_A; ++i) {
        const auto& r_node = r_patch_A[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
    }

    // then the minus geometry
    for (IndexType i = 0; i < number_of_control_points_B; ++i) {
        const auto& r_node = r_patch_B[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
    }
}

int SolidCouplingCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    return Condition::Check(rCurrentProcessInfo);
}

void SolidCouplingCondition::GetSolutionCoefficientVectorA(
    Vector& rValues) const
{
    const auto& r_patch_A = GetGeometry();
    const std::size_t n_patch_A = r_patch_A.size();
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const std::size_t dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const std::size_t nvals = n_patch_A * dim;
    if (rValues.size() != nvals) {
        rValues.resize(nvals, false);
    }
    IndexType k = 0;
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const auto& u = r_patch_A[i].GetSolutionStepValue(DISPLACEMENT);
        rValues[k++] = u[0];
        rValues[k++] = u[1];
        if (dim == 3) rValues[k++] = u[2];
    }
}

void SolidCouplingCondition::GetSolutionCoefficientVectorB(
    Vector& rValues) const
{
    const auto& r_patch_B = GetGeometryMirror();
    const std::size_t n_patch_B = r_patch_B.size();
    const auto& rDN_De_B = r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());
    const std::size_t dim = rDN_De_B.empty() ? 2 : rDN_De_B[0].size2();
    const std::size_t nvals = n_patch_B * dim;
    if (rValues.size() != nvals) {
        rValues.resize(nvals, false);
    }
    IndexType k = 0;
    for (IndexType i = 0; i < n_patch_B; ++i) {
        const auto& u = r_patch_B[i].GetSolutionStepValue(DISPLACEMENT);
        rValues[k++] = u[0];
        rValues[k++] = u[1];
        if (dim == 3) rValues[k++] = u[2];
    }
}

void SolidCouplingCondition::ComputeTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Vector& H_sum_vec)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    if (H_sum_vec.size() != number_of_control_points) H_sum_vec.resize(number_of_control_points);
    // derivatives up to order mBasisFunctionsOrder
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term = 0.0;
        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    const double derivative = r_shape_function_derivatives(i, k);
                    H_taylor_term += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n - 1];
                int countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        const double derivative = r_shape_function_derivatives(i, countDerivativeId);
                        H_taylor_term += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0, i);
    }
}

void SolidCouplingCondition::ComputeGradientTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Matrix& grad_H_sum)
{
    const std::size_t number_of_control_points = rGeometry.PointsNumber();
    const auto& r_DN_De = rGeometry.ShapeFunctionsLocalGradients(rGeometry.GetDefaultIntegrationMethod());
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n - 1] = rGeometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points) grad_H_sum.resize(3, number_of_control_points);
    for (IndexType i = 0; i < number_of_control_points; ++i) {
        double H_taylor_term_X = 0.0, H_taylor_term_Y = 0.0, H_taylor_term_Z = 0.0;
        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n - 1];
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k);
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n - 1; k++) {
                    IndexType n_k = n - 1 - k;
                    const double derivative = shapeFunctionDerivatives(i, k + 1);
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
            }
        } else {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n - 1];
                IndexType countDerivativeId = 0;
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        IndexType k_z = n - k_x - k_y;
                        const double derivative = shapeFunctionDerivatives(i, countDerivativeId);
                        if (k_x >= 1) H_taylor_term_X += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x - 1, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        if (k_y >= 1) H_taylor_term_Y += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y - 1, rDistanceVector[2], k_z);
                        if (k_z >= 1) H_taylor_term_Z += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z - 1);
                        countDerivativeId++;
                    }
                }
            }
        }
        grad_H_sum(0, i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1, i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2, i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else
            grad_H_sum(2, i) = 0.0;
    }
}

double SolidCouplingCondition::ComputeTaylorTerm(const double derivative, const double dx, const IndexType n_k, const double dy, const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));
}

double SolidCouplingCondition::ComputeTaylorTerm3D(const double derivative, const double dx, const IndexType k_x, const double dy, const IndexType k_y, const double dz, const IndexType k_z)
{
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));
}

void SolidCouplingCondition::CalculateB(
    const GeometryType& rGeometry,
    Matrix& rB,
    Matrix& r_DN_DX) const
{
    const std::size_t number_of_control_points = rGeometry.size();
    const std::size_t mat_size = number_of_control_points * 2;

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

void SolidCouplingCondition::ApplyConstitutiveLaw(
    std::size_t matSize,
    Vector& rStrain,
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveVariables& rConstitutiVariables)
{
    Flags& ConstitutiveLawOptions = rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    rValues.SetStrainVector(rStrain);
    rValues.SetStressVector(rConstitutiVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy);
}

} // namespace Kratos
