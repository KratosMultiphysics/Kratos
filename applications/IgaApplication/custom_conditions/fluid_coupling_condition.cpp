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

// System includes

// External includes

// Project includes
#include "custom_conditions/fluid_coupling_condition.h"
#include "includes/variables.h"
#include "iga_application_variables.h"

namespace Kratos
{
typedef std::size_t SizeType;
typedef std::size_t IndexType;


void FluidCouplingCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeMaterial();
    const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    const double integration_weight = r_integration_points[0].Weight();
    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void FluidCouplingCondition::InitializeMemberVariables()
{
    const auto& r_geometry_patchA = GetGeometry();
    const auto& r_geometry_patchB = GetGeometryMirror();

    KRATOS_ERROR_IF(norm_2(r_geometry_patchA.Center()-r_geometry_patchB.Center()) > 1e-12)
        << "FluidCouplingCondition found non matching geometries." << std::endl;

    mNormalParameterSpaceA = r_geometry_patchA.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceA /= MathUtils<double>::Norm(mNormalParameterSpaceA);
    mNormalPhysicalSpaceA = mNormalParameterSpaceA;

    mNormalParameterSpaceB = r_geometry_patchB.Normal(0, GetIntegrationMethod());
    mNormalParameterSpaceB /= MathUtils<double>::Norm(mNormalParameterSpaceB);
    mNormalPhysicalSpaceB = mNormalParameterSpaceB;

    KRATOS_ERROR_IF(std::abs(inner_prod(mNormalPhysicalSpaceA, mNormalPhysicalSpaceB)+1) > 1e-12)
        << "FluidCouplingCondition: normals are not opposite." << std::endl;

    // Estimate basis function order from number of control points and dimension (as in SupportFluidCondition)
    const auto& rDN_De_A = r_geometry_patchA.ShapeFunctionsLocalGradients(r_geometry_patchA.GetDefaultIntegrationMethod());
    const unsigned int dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    // store dimension for later use (e.g., in ApplyConstitutiveLaw)
    mDim = dim;
    const SizeType n_ctrl = rDN_De_A.empty() ? r_geometry_patchA.PointsNumber() : rDN_De_A[0].size1();
    if (dim == 3) {
        mBasisFunctionsOrder = static_cast<IndexType>(std::cbrt(static_cast<double>(n_ctrl)) - 1.0);
    } else {
        mBasisFunctionsOrder = static_cast<IndexType>(std::sqrt(static_cast<double>(n_ctrl)) - 1.0);
    }
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);

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
}

void FluidCouplingCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Pre-size and initialize LHS/RHS consistent with this condition DOFs
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    KRATOS_ERROR_IF(n_patch_A + n_patch_B == 0)
        << "FluidCouplingCondition::CalculateLocalSystem: empty geometry." << std::endl;

    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const unsigned int dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const SizeType dofs_per_node = dim + 1; // velocities (dim) + pressure
    const SizeType total_size = dofs_per_node * (n_patch_A + n_patch_B);

    if (rLeftHandSideMatrix.size1() != total_size || rLeftHandSideMatrix.size2() != total_size) {
        rLeftHandSideMatrix.resize(total_size, total_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);

    if (rRightHandSideVector.size() != total_size) {
        rRightHandSideVector.resize(total_size, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(total_size);

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void FluidCouplingCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    // Each node: mDim velocity components + 1 pressure
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A_all =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const unsigned int dim = rDN_De_A_all[0].size2();
    const SizeType dofs_per_node = dim + 1;
    const SizeType total_size = dofs_per_node * (n_patch_A + n_patch_B);

    KRATOS_ERROR_IF(total_size == 0) << "FluidCouplingCondition found empty geometry." << std::endl;

    if (rLeftHandSideMatrix.size1() != total_size || rLeftHandSideMatrix.size2() != total_size) {
        rLeftHandSideMatrix.resize(total_size, total_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);

    const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    const auto& integration_points_patch_B = r_patch_B.IntegrationPoints();

    KRATOS_ERROR_IF(integration_points_patch_A.empty() || integration_points_patch_B.empty())
        << "FluidCouplingCondition expects at least one integration point." << std::endl;

    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_B =
        r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());

    KRATOS_ERROR_IF(rDN_De_A.size() == 0 || rDN_De_B.size() == 0)
        << "Shape function gradients not available for FluidCouplingCondition. Configure the modeler with derivative order >= 2." << std::endl;

    KRATOS_ERROR_IF(rDN_De_A[0].size1() != n_patch_A)
        << "FluidCouplingCondition: gradient matrix rows mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << rDN_De_A[0].size1() << "." << std::endl;
    KRATOS_ERROR_IF(rDN_De_B[0].size1() != n_patch_B)
        << "FluidCouplingCondition: gradient matrix rows mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << rDN_De_B[0].size1() << "." << std::endl;

    const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    const Matrix& N_patch_B = r_patch_B.ShapeFunctionsValues();

    KRATOS_ERROR_IF(N_patch_A.size2() != n_patch_A)
        << "FluidCouplingCondition: shape function values columns mismatch for patch A. Expected "
        << n_patch_A << ", obtained " << N_patch_A.size2() << "." << std::endl;
    KRATOS_ERROR_IF(N_patch_B.size2() != n_patch_B)
        << "FluidCouplingCondition: shape function values columns mismatch for patch B. Expected "
        << n_patch_B << ", obtained " << N_patch_B.size2() << "." << std::endl;

    const unsigned int dim_patch_A = rDN_De_A[0].size2();
    const unsigned int dim_patch_B = rDN_De_B[0].size2();

    KRATOS_ERROR_IF(dim_patch_A < 2 || dim_patch_B < 2)
        << "FluidCouplingCondition requires at least 2D parameter space." << std::endl;

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
        << "FluidCouplingCondition: degenerate Jacobian for patch A (determinant ~ 0)." << std::endl;
    const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;
    KRATOS_ERROR_IF(std::abs(weight) <= std::numeric_limits<double>::epsilon())
        << "FluidCouplingCondition: zero integration weight on patch A." << std::endl;

    double penalty_factor;
    double nitsche_stabilization_factor = 0.0;//-1.0;
    if (GetProperties().Has(PENALTY_FACTOR)) {
        // penalty_factor = -1.0; // GetProperties()[PENALTY_FACTOR];
        penalty_factor = 0.0;//1e3; 
    } else {
        KRATOS_ERROR << "FluidCouplingCondition requires PENALTY_FACTOR to be defined in the Properties." << std::endl;
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
            KRATOS_ERROR << "KNOT_SPAN_SIZES is empty in FluidCouplingCondition" << std::endl;
        }
    }
    if (characteristic_length <= 0.0) {
        KRATOS_ERROR << "KNOT_SPAN_SIZES is < 0 in FluidCouplingCondition" << std::endl;
    }
    const double penalty_over_h = penalty_factor * mBasisFunctionsOrder * mBasisFunctionsOrder / characteristic_length;

    // Helper to compute index in block with (node, component, side)
    auto idxA = [dofs_per_node](IndexType node, IndexType idim){ return node * dofs_per_node + idim; };
    auto idxB = [dofs_per_node, n_patch_A](IndexType node, IndexType idim){ return n_patch_A * dofs_per_node + node * dofs_per_node + idim; };


    // --- Precompute traction operators on each side: Tn = (D * B) applied to normal ---
    Matrix DN_DX_A = rDN_De_A[0];
    Matrix DN_DX_B = rDN_De_B[0];

    Matrix B_A(3, n_patch_A*dim, 0.0), B_B(3, n_patch_B*dim, 0.0);
    CalculateB(B_A, DN_DX_A); // same shape as in SupportFluidCondition
    CalculateB(B_B, DN_DX_B);

    // constitutive law A
    ConstitutiveLaw::Parameters Values_A(r_patch_A, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables_A(3);
    ApplyConstitutiveLaw(r_patch_A, B_A, Values_A, constitutive_variables_A);
    const Matrix& r_D_A = Values_A.GetConstitutiveMatrix();
    Matrix DB_A = Matrix(prod(r_D_A, B_A));

    // constitutive law B
    ConstitutiveLaw::Parameters Values_B(r_patch_B, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables_B(3);
    ApplyConstitutiveLaw(r_patch_B, B_B, Values_B, constitutive_variables_B);
    const Matrix& r_D_B = Values_B.GetConstitutiveMatrix();
    Matrix DB_B = Matrix(prod(r_D_B, B_B));

    array_1d<double,2> normal_A; 
    normal_A[0]=mNormalPhysicalSpaceA[0]; 
    normal_A[1]=mNormalPhysicalSpaceA[1];
    array_1d<double,2> normal_B; 
    normal_B[0]=mNormalPhysicalSpaceB[0]; 
    normal_B[1]=mNormalPhysicalSpaceB[1];












    // ANDREA's CODE:
    const auto& r_surrogate_geometry_plus = GetGeometry();
    const auto& r_surrogate_geometry_minus = GetGeometryMirror();
    const auto& r_true_geometry = GetGeometry();

    const SizeType number_of_control_points_plus = r_surrogate_geometry_plus.size();
    const SizeType number_of_control_points_minus = r_surrogate_geometry_minus.size();

    const SizeType number_of_control_points = number_of_control_points_plus + number_of_control_points_minus;

    const SizeType mat_size_plus = number_of_control_points_plus * mDim;
    const SizeType mat_size_minus = number_of_control_points_minus * mDim;

    // reading integration points and local gradients
    const double integration_weight = weight;

    // // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec_plus = row(N_patch_A, 0); 
    Vector N_sum_vec_minus = row(N_patch_B, 0);

    Matrix grad_N_sum_plus = rDN_De_A[0]; 
    Matrix grad_N_sum_minus = rDN_De_B[0];


    // compute the B matrix
    Matrix B_sum_plus = ZeroMatrix(mDim,mat_size_plus);
    CalculateB(B_sum_plus, grad_N_sum_plus);

    Matrix B_sum_minus = ZeroMatrix(mDim,mat_size_minus);
    CalculateB(B_sum_minus, grad_N_sum_minus);

    // compute the DB product
    Matrix DB_sum_plus = prod(r_D_A, B_sum_plus);
    Matrix DB_sum_minus = prod(r_D_B, B_sum_minus);

    // ASSEMBLE
    //-----------------------------------------------------
    const SizeType shift_dof = mat_size_plus + number_of_control_points_plus; // shift for minus side DOFs
    // -w_plus * (sigma_plus + sigma_minus) /2
    for (IndexType i = 0; i < number_of_control_points_plus; i++) {
        // FIRST TERM: -w_plus * sigma_plus /2
        for (IndexType j = 0; j < number_of_control_points_plus; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 3*i+idim;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob_DB = 2*j+jdim;
                    const int jglob = 3*j+jdim;

                    rLeftHandSideMatrix(iglob, jglob) -= 0.5 * N_sum_vec_plus(i)
                                                    * (DB_sum_plus(id1, jglob_DB)* normal_A[0] + DB_sum_plus(id2, jglob_DB)* normal_A[1]) 
                                                    * integration_weight;
                }
            }
        }

        // SECOND TERM: -w_plus * sigma_minus /2
        for (IndexType j = 0; j < number_of_control_points_minus; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 3*i+idim;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jloc = 2*j+jdim;
                    const int jglob = 3*j+jdim + shift_dof;

                    rLeftHandSideMatrix(iglob, jglob) -= 0.5 * N_sum_vec_plus(i)
                                                    * (DB_sum_minus(id1, jloc)* normal_A[0] + DB_sum_minus(id2, jloc)* normal_A[1]) 
                                                    * integration_weight;
                }
            }
        }
    }

    // +w_minus * (sigma_plus + sigma_minus) /2
    for (IndexType i = 0; i < number_of_control_points_minus; i++) {
        // FIRST TERM: +w_minus * sigma_plus /2
        for (IndexType j = 0; j < number_of_control_points_plus; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 3*i+idim + shift_dof;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob_DB = 2*j+jdim;
                    const int jglob = 3*j+jdim;

                    rLeftHandSideMatrix(iglob, jglob) += 0.5 * N_sum_vec_minus(i)
                                                    * (DB_sum_plus(id1, jglob_DB)* normal_A[0] + DB_sum_plus(id2, jglob_DB)* normal_A[1]) 
                                                    * integration_weight;
                }
            }
        }

        // SECOND TERM: +w_minus * sigma_minus /2
        for (IndexType j = 0; j < number_of_control_points_minus; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 3*i+idim + shift_dof;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jloc = 2*j+jdim;
                    const int jglob = 3*j+jdim + shift_dof;

                    rLeftHandSideMatrix(iglob, jglob) += 0.5 * N_sum_vec_minus(i)
                                                    * (DB_sum_minus(id1, jloc)* normal_A[0] + DB_sum_minus(id2, jloc)* normal_A[1]) 
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
        for (IndexType idim = 0; idim < 2; idim++) {
            const int iglob = 3*i+idim;
            const int iglob_DB = 2*i+idim;

            Vector tau_w_n_A = ZeroVector(3);
            tau_w_n_A[0] = (DB_sum_plus(0, iglob_DB)* normal_A[0] + DB_sum_plus(2, iglob_DB)* normal_A[1]);
            tau_w_n_A[1] = (DB_sum_plus(2, iglob_DB)* normal_A[0] + DB_sum_plus(1, iglob_DB)* normal_A[1]);

            for (IndexType j = 0; j < number_of_control_points_plus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 3*j+idim) += N_sum_vec_plus(i)*N_sum_vec_plus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int jglob = 3*j+jdim;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*mNitschePenalty*N_sum_vec_plus(j)*tau_w_n_A[jdim] * integration_weight;
                }

            }

            // minus side
            for (IndexType j = 0; j < number_of_control_points_minus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 3*j+shift_dof+idim) -= N_sum_vec_plus(i)*N_sum_vec_minus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int jglob = 3*j+jdim + shift_dof;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += 0.5*mNitschePenalty*N_sum_vec_minus(j)*tau_w_n_A[jdim] * integration_weight;
                }
            }
        }
    }

    // minus side
    for (IndexType i = 0; i < number_of_control_points_minus; i++) {
        for (IndexType idim = 0; idim < 2; idim++) {
            const int iloc = 2*i+idim;
            const int iglob = 3*i+idim + shift_dof;

            Vector tau_w_n_B = ZeroVector(3);
            tau_w_n_B[0] = (DB_sum_minus(0, iloc)* normal_A[0] + DB_sum_minus(2, iloc)* normal_A[1]);
            tau_w_n_B[1] = (DB_sum_minus(2, iloc)* normal_A[0] + DB_sum_minus(1, iloc)* normal_A[1]);

            for (IndexType j = 0; j < number_of_control_points_plus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 3*j+idim) -= N_sum_vec_minus(i)*N_sum_vec_plus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int jglob = 3*j+jdim;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) -= 0.5*mNitschePenalty*N_sum_vec_plus(j)*tau_w_n_B[jdim] * integration_weight;
                }

            }

            // minus side
            for (IndexType j = 0; j < number_of_control_points_minus; j++) {
                // PENALTY TERM
                rLeftHandSideMatrix(iglob, 3*j+shift_dof+idim) += N_sum_vec_minus(i)*N_sum_vec_minus(j)* penalty_integration;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int jglob = 3*j+jdim + shift_dof;

                    // // PENALTY FREE g_n = 0
                    // // [\sigma_1(w) \dot n] \dot n (-u_1 \dot n)
                    // //*********************************************** */
                    rLeftHandSideMatrix(iglob, jglob) += 0.5*mNitschePenalty*N_sum_vec_minus(j)*tau_w_n_B[jdim] * integration_weight;
                }
            }
        }
    }



    // // --- first block: -[[w]], {{sigma n}} 
    // for (IndexType idim = 0; idim < dim; ++idim) {
    //     // rows: test from A  -> contributes - w_A * 0.5*(tau_A n_A + tau_B n_B)
    //     for (IndexType i = 0; i < n_patch_A; ++i) {
    //         const double vA = N_patch_A(0, i);

    //         // columns from A -> - w_A * 0.5*(tau_A n_A)
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution = ZeroMatrix(2, 2);
    //                 DB_contribution(0, 0) = DB_A(0, 2*j+jdim);
    //                 DB_contribution(0, 1) = DB_A(2, 2*j+jdim);
    //                 DB_contribution(1, 0) = DB_A(2, 2*j+jdim);
    //                 DB_contribution(1, 1) = DB_A(1, 2*j+jdim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction = prod(DB_contribution, normal_A);
    //                 rLeftHandSideMatrix(idxA(i, idim), idxA(j, jdim)) -= 0.5 * vA * traction(idim) * weight;
    //             }
    //         }
    //         // columns from B -> - w_A * 0.5*(tau_B n_B)
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution = ZeroMatrix(2, 2);
    //                 DB_contribution(0, 0) = DB_B(0, 2*j+jdim);
    //                 DB_contribution(0, 1) = DB_B(2, 2*j+jdim);
    //                 DB_contribution(1, 0) = DB_B(2, 2*j+jdim);
    //                 DB_contribution(1, 1) = DB_B(1, 2*j+jdim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction = prod(DB_contribution, normal_B);
    //                 rLeftHandSideMatrix(idxA(i, idim), idxB(j, jdim)) -= 0.5 * vA * traction(idim) * weight;
    //             }
    //         }

    //         // // pressure contribution to traction: -[[w]] · 0.5*(-p_A n_A - p_B n_B)
    //         // // -> for A rows: - w_A · 0.5*(-p_A n_A - p_B n_B) = +0.5 w_A · (p_A n_A + p_B n_B)
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     rLeftHandSideMatrix(idxA(i, idim), idxA(j, dim)) += weight * 0.5 * vA * normal_A[idim] * N_patch_A(0, j);
    //         // }
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     rLeftHandSideMatrix(idxA(i, idim), idxB(j, dim)) += weight * 0.5 * vA * normal_B[idim] * N_patch_B(0, j);
    //         // }
    //     }



    //     // rows: test from B  -> contributes + w_B * 0.5*(tau_A n_A + tau_B n_B)
    //     for (IndexType i = 0; i < n_patch_B; ++i) {
    //         const double vB = N_patch_B(0, i);

    //         // columns from A -> + w_B * 0.5*(tau_A n_A)
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution = ZeroMatrix(2, 2);
    //                 DB_contribution(0, 0) = DB_A(0, 2*j+jdim);
    //                 DB_contribution(0, 1) = DB_A(2, 2*j+jdim);
    //                 DB_contribution(1, 0) = DB_A(2, 2*j+jdim);
    //                 DB_contribution(1, 1) = DB_A(1, 2*j+jdim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction = prod(DB_contribution, normal_A);
    //                 rLeftHandSideMatrix(idxB(i, idim), idxA(j, jdim)) += 0.5 * vB * traction(idim) * weight;
    //             }
    //         }
    //         // columns from B -> + w_B * 0.5*(tau_B n_B)
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution = ZeroMatrix(2, 2);  
    //                 DB_contribution(0, 0) = DB_B(0, 2*j+jdim);
    //                 DB_contribution(0, 1) = DB_B(2, 2*j+jdim);
    //                 DB_contribution(1, 0) = DB_B(2, 2*j+jdim);
    //                 DB_contribution(1, 1) = DB_B(1, 2*j+jdim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction = prod(DB_contribution, normal_B);
    //                 rLeftHandSideMatrix(idxB(i, idim), idxB(j, jdim)) += 0.5 * vB * traction(idim) * weight;
    //             }
    //         }

    //         // // pressure contribution to traction for B rows: + w_B · 0.5*(-p_A n_A - p_B n_B)
    //         // // -> gives -0.5 w_B · (p_A n_A + p_B n_B) in the matrix
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     rLeftHandSideMatrix(idxB(i, idim), idxA(j, dim)) -= weight * 0.5 * vB * normal_A[idim] * N_patch_A(0, j);
    //         // }
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     rLeftHandSideMatrix(idxB(i, idim), idxB(j, dim)) -= weight * 0.5 * vB * normal_B[idim] * N_patch_B(0, j);
    //         // }
    //     }
    // }



    // // Nitsche: + {{tau(w)}}, [[u]] 
    // double factor = 0.5*nitsche_stabilization_factor;

    // // (tau(w_A)+tau(w_B)) * 0.5 * (u_A - u_B)
    // for (IndexType idim = 0; idim < dim; ++idim) {
    //     // tau(w_A)  -> contributes + tau(w_A)*n_A *0.5 (u_A - u_B)
    //     for (IndexType i = 0; i < n_patch_A; ++i) {
    //         // columns from A -> + tau(w_A)*n_A *0.5 * u_A
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             const double uA = N_patch_A(0, j);
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution_w = ZeroMatrix(2, 2);
    //                 DB_contribution_w(0, 0) = DB_A(0, 2*i+idim);
    //                 DB_contribution_w(0, 1) = DB_A(2, 2*i+idim);
    //                 DB_contribution_w(1, 0) = DB_A(2, 2*i+idim);
    //                 DB_contribution_w(1, 1) = DB_A(1, 2*i+idim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction_nitsche_w = prod(DB_contribution_w, normal_A);
    //                 rLeftHandSideMatrix(idxA(i, idim), idxA(j, jdim)) += factor * uA * traction_nitsche_w(jdim) * weight;
    //             }
    //         }
    //         // columns from B -> - tau(w_A)*n_A *0.5 * u_B
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             const double uB = N_patch_B(0, j);
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution_w = ZeroMatrix(2, 2);
    //                 DB_contribution_w(0, 0) = DB_A(0, 2*i+idim);
    //                 DB_contribution_w(0, 1) = DB_A(2, 2*i+idim);
    //                 DB_contribution_w(1, 0) = DB_A(2, 2*i+idim);
    //                 DB_contribution_w(1, 1) = DB_A(1, 2*i+idim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction_nitsche_w = prod(DB_contribution_w, normal_A);
    //                 rLeftHandSideMatrix(idxA(i, idim), idxB(j, jdim)) -= factor * uB * traction_nitsche_w(jdim) * weight;
    //             }
    //         }

    //         // // pressure contribution to traction: + {{-q n}}, [[u]] 
    //         // // -> for A rows: - q_A*n_A · 0.5*(u_A)
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     const double uA = N_patch_A(0, j);
    //         //     rLeftHandSideMatrix(idxA(i, dim), idxA(j, idim)) -= weight * factor * normal_A[idim] * N_patch_A(0, i) * uA ;
    //         // } 
    //         // // -> for A rows: + q_A*n_A · 0.5*(u_B) 
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     const double uB = N_patch_B(0, j);
    //         //     rLeftHandSideMatrix(idxA(i, dim), idxB(j, idim)) += weight * factor * normal_A[idim] * N_patch_A(0, i) * uB ;
    //         // }
    //     }



    //     // tau(w_B)  -> contributes + tau(w_B)*n_B *0.5 (u_A - u_B)
    //     for (IndexType i = 0; i < n_patch_B; ++i) {
    //         // columns from A -> + tau(w_B)*n_B *0.5 * u_A
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             const double uA = N_patch_A(0, j);
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution_w_B = ZeroMatrix(2, 2);
    //                 DB_contribution_w_B(0, 0) = DB_B(0, 2*i+idim);
    //                 DB_contribution_w_B(0, 1) = DB_B(2, 2*i+idim);
    //                 DB_contribution_w_B(1, 0) = DB_B(2, 2*i+idim);
    //                 DB_contribution_w_B(1, 1) = DB_B(1, 2*i+idim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction_nitsche_w = prod(DB_contribution_w_B, normal_B);
    //                 rLeftHandSideMatrix(idxB(i, idim), idxA(j, jdim)) += factor * uA * traction_nitsche_w(jdim) * weight;
    //             }
    //         }
    //         // columns from B -> - tau(w_B)*n_B *0.5 * u_B
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             const double uB = N_patch_B(0, j);
    //             for (IndexType jdim = 0; jdim < 2; jdim++) {
    //                 Matrix DB_contribution_w_B = ZeroMatrix(2, 2);
    //                 DB_contribution_w_B(0, 0) = DB_B(0, 2*i+idim);
    //                 DB_contribution_w_B(0, 1) = DB_B(2, 2*i+idim);
    //                 DB_contribution_w_B(1, 0) = DB_B(2, 2*i+idim);
    //                 DB_contribution_w_B(1, 1) = DB_B(1, 2*i+idim);
    //                 // Compute the traction vector: sigma * n.
    //                 Vector traction_nitsche_w = prod(DB_contribution_w_B, normal_B);
    //                 rLeftHandSideMatrix(idxB(i, idim), idxB(j, jdim)) -= factor * uB * traction_nitsche_w(jdim) * weight;
    //             }
    //         }

    //         // // pressure contribution to traction: + {{-q n}}, [[u]] 
    //         // // -> for A rows: - q_B*n_B · 0.5*(u_A)
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     const double uA = N_patch_A(0, j);
    //         //     rLeftHandSideMatrix(idxB(i, dim), idxA(j, idim)) -= weight * factor * normal_B[idim] * N_patch_B(0, i) * uA ;
    //         // } 
    //         // // -> for A rows: + q_B*n_B · 0.5*(u_B) 
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     const double uB = N_patch_B(0, j);
    //         //     rLeftHandSideMatrix(idxB(i, dim), idxB(j, idim)) += weight * factor * normal_B[idim] * N_patch_B(0, i) * uB ;
    //         // }
    //     }
    // }

    
    // // Penalty term
    // if (penalty_over_h != 0.0) {

    //     // // new new
    //     // const double g = penalty_over_h * weight;
    //     // // unit normals already built
    //     // array_1d<double,2> nA = normal_A;                  // n_B = -n_A holds
    //     // array_1d<double,2> nB = normal_B;
    //     // // projectors
    //     // double Pn[2][2] = { { nA[0]*nA[0], nA[0]*nA[1] },
    //     //                     { nA[1]*nA[0], nA[1]*nA[1] } };
    //     // double Pt[2][2] = { { 1.0 - Pn[0][0],   -Pn[0][1] },
    //     //                     {   -Pn[1][0],    1.0 - Pn[1][1] } };
    //     // // choose weights (set both to 1.0 to enforce full continuity)
    //     // const double beta_n = 1.0;
    //     // const double beta_t = 1.0;
    //     // // P = beta_n * nn^T + beta_t * (I - nn^T)
    //     // double P[2][2] = {
    //     // { beta_n*Pn[0][0] + beta_t*Pt[0][0], beta_n*Pn[0][1] + beta_t*Pt[0][1] },
    //     // { beta_n*Pn[1][0] + beta_t*Pt[1][0], beta_n*Pn[1][1] + beta_t*Pt[1][1] }
    //     // };
    //     // for (IndexType i=0;i<n_patch_A;++i){
    //     //     for (IndexType j=0;j<n_patch_A;++j){
    //     //         const double Nij = N_patch_A(0,i)*N_patch_A(0,j);
    //     //         for (unsigned p=0;p<2;++p) for (unsigned q=0;q<2;++q)
    //     //         rLeftHandSideMatrix(idxA(i,p), idxA(j,q)) += g * P[p][q] * Nij;        // AA
    //     //     }
    //     //     for (IndexType j=0;j<n_patch_B;++j){
    //     //         const double Nij = N_patch_A(0,i)*N_patch_B(0,j);
    //     //         for (unsigned p=0;p<2;++p) for (unsigned q=0;q<2;++q)
    //     //         rLeftHandSideMatrix(idxA(i,p), idxB(j,q)) -= g * P[p][q] * Nij;        // AB (minus)
    //     //     }
    //     //     }
    //     //     for (IndexType i=0;i<n_patch_B;++i){
    //     //     for (IndexType j=0;j<n_patch_A;++j){
    //     //         const double Nij = N_patch_B(0,i)*N_patch_A(0,j);
    //     //         for (unsigned p=0;p<2;++p) for (unsigned q=0;q<2;++q)
    //     //         rLeftHandSideMatrix(idxB(i,p), idxA(j,q)) -= g * P[p][q] * Nij;        // BA (minus)
    //     //     }
    //     //     for (IndexType j=0;j<n_patch_B;++j){
    //     //         const double Nij = N_patch_B(0,i)*N_patch_B(0,j);
    //     //         for (unsigned p=0;p<2;++p) for (unsigned q=0;q<2;++q)
    //     //         rLeftHandSideMatrix(idxB(i,p), idxB(j,q)) += g * P[p][q] * Nij;        // BB
    //     //     }
    //     // }

    //     // penalty only continuity
    //     const double g = penalty_over_h * weight;
    //     // enforce continuity component-wise: AA += g N_i N_j, AB -= g N_i M_j, ...
    //     for (IndexType idim = 0; idim < dim; ++idim) {
    //         // A–A
    //         for (IndexType i=0;i<n_patch_A;++i){
    //             const double Ni = N_patch_A(0,i);
    //             const IndexType row = idxA(i,idim);
    //             for (IndexType j=0;j<n_patch_A;++j){
    //                 const double Nj = N_patch_A(0,j);
    //                 rLeftHandSideMatrix(row, idxA(j,idim)) += g * Ni * Nj;
    //             }
    //             // A–B
    //             for (IndexType j=0;j<n_patch_B;++j){
    //                 const double Mj = N_patch_B(0,j);
    //                 rLeftHandSideMatrix(row, idxB(j,idim)) -= g * Ni * Mj;
    //             }
    //         }
    //         // B–A and B–B
    //         for (IndexType i=0;i<n_patch_B;++i){
    //             const double Mi = N_patch_B(0,i);
    //             const IndexType row = idxB(i,idim);
    //             for (IndexType j=0;j<n_patch_A;++j){
    //                 const double Nj = N_patch_A(0,j);
    //                 rLeftHandSideMatrix(row, idxA(j,idim)) -= g * Mi * Nj;
    //             }
    //             for (IndexType j=0;j<n_patch_B;++j){
    //                 const double Mj = N_patch_B(0,j);
    //                 rLeftHandSideMatrix(row, idxB(j,idim)) += g * Mi * Mj;
    //             }
    //         }
    //     }


    // }

    // // // Zero out all pressure (q) rows and columns to avoid spurious coupling
    // // {
    // //     const SizeType n_total = rLeftHandSideMatrix.size1();
    // //     // Zero rows/cols for pressures on patch A
    // //     for (IndexType i = 0; i < n_patch_A; ++i) {
    // //         const IndexType p_index = idxA(i, dim);
    // //         // zero row
    // //         for (IndexType j = 0; j < n_total; ++j) {
    // //             rLeftHandSideMatrix(p_index, j) = 0.0;
    // //         }
    // //         // zero column
    // //         for (IndexType j = 0; j < n_total; ++j) {
    // //             rLeftHandSideMatrix(j, p_index) = 0.0;
    // //         }
    // //     }
    // //     // Zero rows/cols for pressures on patch B
    // //     for (IndexType i = 0; i < n_patch_B; ++i) {
    // //         const IndexType p_index = idxB(i, dim);
    // //         // zero row
    // //         for (IndexType j = 0; j < n_total; ++j) {
    // //             rLeftHandSideMatrix(p_index, j) = 0.0;
    // //         }
    // //         // zero column
    // //         for (IndexType j = 0; j < n_total; ++j) {
    // //             rLeftHandSideMatrix(j, p_index) = 0.0;
    // //         }
    // //     }
    // // }    
    

    // // Symmetry check of the resulting LHS (lightweight, early-exit)
    // {
    //     const SizeType n = rLeftHandSideMatrix.size1();
    //     double max_asym = 0.0;
    //     const double tol = 1e-12; // tight tolerance for assembly-level symmetry
    //     for (SizeType i = 0; i < n; ++i) {
    //         for (SizeType j = i + 1; j < n; ++j) {
    //             const double diff = std::abs(rLeftHandSideMatrix(i, j) - rLeftHandSideMatrix(j, i));
    //             if (diff > max_asym) {
    //                 max_asym = diff;
    //                 if (max_asym > tol) break; // early exit on first violation
    //             }
    //         }
    //         if (max_asym > tol) break;
    //     }
    //     KRATOS_INFO_IF("maiooo", max_asym > tol)
    //         << "Assembled LHS is not symmetric. Max asymmetry = " << max_asym << std::endl;
    // }
    

    KRATOS_CATCH("")
}

void FluidCouplingCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and initialize RHS vector (dim velocities + pressure per node on both patches)
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_A =
        r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& rDN_De_B =
        r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());
    const unsigned int dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const SizeType dofs_per_node = dim + 1; // velocities (dim) + pressure
    const SizeType total_size = dofs_per_node * (n_patch_A + n_patch_B);

    if (rRightHandSideVector.size() != total_size) {
        rRightHandSideVector.resize(total_size, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(total_size);


    MatrixType lhs;
    CalculateLeftHandSide(lhs, rCurrentProcessInfo);
    const SizeType size = lhs.size1();
    KRATOS_ERROR_IF(size == 0) << "FluidCouplingCondition::CalculateRightHandSide: The left hand side matrix has zero size." << std::endl;

    auto idxA = [dofs_per_node](IndexType node, IndexType idim){ return node * dofs_per_node + idim; };
    auto idxB = [dofs_per_node, n_patch_A](IndexType node, IndexType idim){ return n_patch_A * dofs_per_node + node * dofs_per_node + idim; };
    Vector valsA, valsB;
    GetVelocityCoefficientVectorA(valsA);
    GetVelocityCoefficientVectorB(valsB);
    Vector values_full(size, false);
    noalias(values_full) = ZeroVector(size);

    // Fill velocity and pressure entries from current iteration values
    for (IndexType i = 0; i < n_patch_A; ++i) {
        for (IndexType c = 0; c < dim; ++c) {
            values_full[idxA(i, c)] = valsA[i * dim + c];
        }
        // pressure A
        values_full[idxA(i, dim)] = r_patch_A[i].FastGetSolutionStepValue(PRESSURE);
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        for (IndexType c = 0; c < dim; ++c) {
            values_full[idxB(i, c)] = valsB[i * dim + c];
        }
        // pressure B
        values_full[idxB(i, dim)] = r_patch_B[i].FastGetSolutionStepValue(PRESSURE);
    }

    if (rRightHandSideVector.size() != size) {
        rRightHandSideVector.resize(size, false);
    }
    noalias(rRightHandSideVector) = -prod(lhs, values_full);





    // Proper RHS

    // const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    // const auto& integration_points_patch_B = r_patch_B.IntegrationPoints();
    // KRATOS_ERROR_IF(integration_points_patch_A.empty() || integration_points_patch_B.empty())
    //     << "FluidCouplingCondition expects at least one integration point." << std::endl;
    // const Matrix& N_patch_A = r_patch_A.ShapeFunctionsValues();
    // const Matrix& N_patch_B = r_patch_B.ShapeFunctionsValues();
    // KRATOS_ERROR_IF(N_patch_A.size2() != n_patch_A)
    //     << "FluidCouplingCondition: shape function values columns mismatch for patch A. Expected "
    //     << n_patch_A << ", obtained " << N_patch_A.size2() << "." << std::endl;
    // KRATOS_ERROR_IF(N_patch_B.size2() != n_patch_B)
    //     << "FluidCouplingCondition: shape function values columns mismatch for patch B. Expected "
    //     << n_patch_B << ", obtained " << N_patch_B.size2() << "." << std::endl;

    // const unsigned int dim_patch_A = rDN_De_A[0].size2();
    // const unsigned int dim_patch_B = rDN_De_B[0].size2();

    // KRATOS_ERROR_IF(dim_patch_A < 2 || dim_patch_B < 2)
    //     << "FluidCouplingCondition requires at least 2D parameter space." << std::endl;

    // // Jacobians to compute the interface measure
    // GeometryType::JacobiansType J_patch_A;
    // r_patch_A.Jacobian(J_patch_A, r_patch_A.GetDefaultIntegrationMethod());
    // Matrix jacobian_patch_A = ZeroMatrix(3, 3);
    // jacobian_patch_A(0, 0) = J_patch_A[0](0, 0);
    // jacobian_patch_A(0, 1) = J_patch_A[0](0, 1);
    // jacobian_patch_A(1, 0) = J_patch_A[0](1, 0);
    // jacobian_patch_A(1, 1) = J_patch_A[0](1, 1);
    // jacobian_patch_A(2, 2) = 1.0;
    // array_1d<double, 3> tangent_patch_A;
    // r_patch_A.Calculate(LOCAL_TANGENT, tangent_patch_A);
    // Vector determinant_factor_patch_A = prod(jacobian_patch_A, tangent_patch_A);
    // determinant_factor_patch_A[2] = 0.0;
    // const double detJ_patch_A = norm_2(determinant_factor_patch_A);
    // KRATOS_ERROR_IF(detJ_patch_A <= std::numeric_limits<double>::epsilon())
    //     << "FluidCouplingCondition: degenerate Jacobian for patch A (determinant ~ 0)." << std::endl;
    // const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;
    // KRATOS_ERROR_IF(std::abs(weight) <= std::numeric_limits<double>::epsilon())
    //     << "FluidCouplingCondition: zero integration weight on patch A." << std::endl;

    // double penalty_factor;
    // double nitsche_stabilization_factor = -1.0;
    // if (GetProperties().Has(PENALTY_FACTOR)) {
    //     // penalty_factor = -1.0; // GetProperties()[PENALTY_FACTOR];
    //     penalty_factor = 1e3; 
    // } else {
    //     KRATOS_ERROR << "FluidCouplingCondition requires PENALTY_FACTOR to be defined in the Properties." << std::endl;
    // }



    // double characteristic_length;
    // if (Has(KNOT_SPAN_SIZES)) {
    //     const Vector& knot_span_sizes = GetValue(KNOT_SPAN_SIZES);
    //     if (!knot_span_sizes.empty()) {
    //         characteristic_length = knot_span_sizes[0];
    //         if (knot_span_sizes.size() > 1) {
    //             characteristic_length = std::min(characteristic_length, knot_span_sizes[1]);
    //         }
    //     } else {
    //         KRATOS_ERROR << "KNOT_SPAN_SIZES is empty in FluidCouplingCondition" << std::endl;
    //     }
    // }
    // if (characteristic_length <= 0.0) {
    //     KRATOS_ERROR << "KNOT_SPAN_SIZES is < 0 in FluidCouplingCondition" << std::endl;
    // }
    // const double penalty_over_h = penalty_factor * mBasisFunctionsOrder * mBasisFunctionsOrder / characteristic_length;

    // // // Helper to compute index in block with (node, component, side)
    // auto idxA = [dofs_per_node](IndexType node, IndexType idim){ return node * dofs_per_node + idim; };
    // auto idxB = [dofs_per_node, n_patch_A](IndexType node, IndexType idim){ return n_patch_A * dofs_per_node + node * dofs_per_node + idim; };


    // // --- Precompute traction operators on each side: Tn = (D * B) applied to normal ---
    // Matrix DN_DX_A = rDN_De_A[0];
    // Matrix DN_DX_B = rDN_De_B[0];

    // Matrix B_A(3, n_patch_A*dim, 0.0), B_B(3, n_patch_B*dim, 0.0);
    // CalculateB(B_A, DN_DX_A); // same shape as in SupportFluidCondition
    // CalculateB(B_B, DN_DX_B);

    // // constitutive law A
    // ConstitutiveLaw::Parameters Values_A(r_patch_A, GetProperties(), rCurrentProcessInfo);
    // ConstitutiveVariables constitutive_variables_A(3);
    // ApplyConstitutiveLaw(r_patch_A, B_A, Values_A, constitutive_variables_A);
    // const Matrix& r_D_A = Values_A.GetConstitutiveMatrix();
    // Matrix DB_A = Matrix(prod(r_D_A, B_A));

    // // constitutive law B
    // ConstitutiveLaw::Parameters Values_B(r_patch_B, GetProperties(), rCurrentProcessInfo);
    // ConstitutiveVariables constitutive_variables_B(3);
    // ApplyConstitutiveLaw(r_patch_B, B_B, Values_B, constitutive_variables_B);
    // const Matrix& r_D_B = Values_B.GetConstitutiveMatrix();
    // Matrix DB_B = Matrix(prod(r_D_B, B_B));

    // // 2D traction mapping from Voigt to vector: t = [ [1,0,0],[0,1,0],[0,0,1] ]* ... then dot n
    // array_1d<double,2> normal_A; normal_A[0]=mNormalPhysicalSpaceA[0]; normal_A[1]=mNormalPhysicalSpaceA[1];
    // array_1d<double,2> normal_B; normal_B[0]=mNormalPhysicalSpaceB[0]; normal_B[1]=mNormalPhysicalSpaceB[1];

    // // Compute the traction vector: sigma * n using r_stress_vector
    // Vector& r_stress_vector_A = Values_A.GetStressVector();
    // Matrix stress_old_A = ZeroMatrix(2, 2);
    // stress_old_A(0, 0) = r_stress_vector_A[0];      stress_old_A(0, 1) = r_stress_vector_A[2];      
    // stress_old_A(1, 0) = r_stress_vector_A[2];      stress_old_A(1, 1) = r_stress_vector_A[1];         
    // Vector traction_current_iteration_A = prod(stress_old_A, normal_A); 

    // Vector& r_stress_vector_B = Values_B.GetStressVector();
    // Matrix stress_old_B = ZeroMatrix(2, 2);
    // stress_old_B(0, 0) = r_stress_vector_B[0];      stress_old_B(0, 1) = r_stress_vector_B[2];      
    // stress_old_B(1, 0) = r_stress_vector_B[2];      stress_old_B(1, 1) = r_stress_vector_B[1];         
    // Vector traction_current_iteration_B = prod(stress_old_B, normal_B); 
    


    // // --- first block: -[[w]], {{sigma n}} 
    // for (IndexType idim = 0; idim < dim; ++idim) {
    //     // rows: test from A  -> contributes - w_A * 0.5*(tau_A n_A + tau_B n_B)
    //     for (IndexType i = 0; i < n_patch_A; ++i) {
    //         const double vA = N_patch_A(0, i);
    //         // // columns from A -> - w_A * 0.5*(tau_A n_A)
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     for (IndexType jdim = 0; jdim < 2; jdim++) {
    //         //         Matrix DB_contribution = ZeroMatrix(2, 2);
    //         //         DB_contribution(0, 0) = DB_A(0, 2*j+jdim);
    //         //         DB_contribution(0, 1) = DB_A(2, 2*j+jdim);
    //         //         DB_contribution(1, 0) = DB_A(2, 2*j+jdim);
    //         //         DB_contribution(1, 1) = DB_A(1, 2*j+jdim);
    //         //         // Compute the traction vector: sigma * n.
    //         //         Vector traction = prod(DB_contribution, normal_A);
    //         //         rLeftHandSideMatrix(idxA(i, idim), idxA(j, jdim)) -= 0.5 * vA * traction(idim) * weight;
    //         //     }
    //         // }
    //         rRightHandSideVector(idxA(i, idim)) += 0.5 * vA * traction_current_iteration_A(idim) * weight;
    //         // // columns from B -> - w_A * 0.5*(tau_B n_B)
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     for (IndexType jdim = 0; jdim < 2; jdim++) {
    //         //         Matrix DB_contribution = ZeroMatrix(2, 2);
    //         //         DB_contribution(0, 0) = DB_B(0, 2*j+jdim);
    //         //         DB_contribution(0, 1) = DB_B(2, 2*j+jdim);
    //         //         DB_contribution(1, 0) = DB_B(2, 2*j+jdim);
    //         //         DB_contribution(1, 1) = DB_B(1, 2*j+jdim);
    //         //         // Compute the traction vector: sigma * n.
    //         //         Vector traction = prod(DB_contribution, normal_B);
    //         //         rLeftHandSideMatrix(idxA(i, idim), idxB(j, jdim)) -= 0.5 * vA * traction(idim) * weight;
    //         //     }
    //         // }
    //         rRightHandSideVector(idxA(i, idim)) += 0.5 * vA * traction_current_iteration_B(idim) * weight;

    //         // // pressure contribution to traction: -[[w]] · 0.5*(-p_A n_A - p_B n_B)
    //         // // -> for A rows: - w_A · 0.5*(-p_A n_A - p_B n_B) = +0.5 w_A · (p_A n_A + p_B n_B)
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     rLeftHandSideMatrix(idxA(i, idim), idxA(j, dim)) += weight * 0.5 * vA * normal_A[idim] * N_patch_A(0, j);
    //         // }
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     rLeftHandSideMatrix(idxA(i, idim), idxB(j, dim)) += weight * 0.5 * vA * normal_B[idim] * N_patch_B(0, j);
    //         // }
    //     }

    //     // rows: test from B  -> contributes + w_B * 0.5*(tau_A n_A + tau_B n_B)
    //     for (IndexType i = 0; i < n_patch_B; ++i) {
    //         const double vB = N_patch_B(0, i);

    //         // // columns from A -> + w_B * 0.5*(tau_A n_A)
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     for (IndexType jdim = 0; jdim < 2; jdim++) {
    //         //         Matrix DB_contribution = ZeroMatrix(2, 2);
    //         //         DB_contribution(0, 0) = DB_A(0, 2*j+jdim);
    //         //         DB_contribution(0, 1) = DB_A(2, 2*j+jdim);
    //         //         DB_contribution(1, 0) = DB_A(2, 2*j+jdim);
    //         //         DB_contribution(1, 1) = DB_A(1, 2*j+jdim);
    //         //         // Compute the traction vector: sigma * n.
    //         //         Vector traction = prod(DB_contribution, normal_A);
    //         //         rLeftHandSideMatrix(idxB(i, idim), idxA(j, jdim)) += 0.5 * vB * traction(idim) * weight;
    //         //     }
    //         // }
    //         rRightHandSideVector(idxB(i, idim)) -= 0.5 * vB * traction_current_iteration_A(idim) * weight;

    //         // columns from B -> + w_B * 0.5*(tau_B n_B)
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     for (IndexType jdim = 0; jdim < 2; jdim++) {
    //         //         Matrix DB_contribution = ZeroMatrix(2, 2);  
    //         //         DB_contribution(0, 0) = DB_B(0, 2*j+jdim);
    //         //         DB_contribution(0, 1) = DB_B(2, 2*j+jdim);
    //         //         DB_contribution(1, 0) = DB_B(2, 2*j+jdim);
    //         //         DB_contribution(1, 1) = DB_B(1, 2*j+jdim);
    //         //         // Compute the traction vector: sigma * n.
    //         //         Vector traction = prod(DB_contribution, normal_B);
    //         //         rLeftHandSideMatrix(idxB(i, idim), idxB(j, jdim)) += 0.5 * vB * traction(idim) * weight;
    //         //     }
    //         // }
    //         rRightHandSideVector(idxB(i, idim)) -= 0.5 * vB * traction_current_iteration_B(idim) * weight;

    //         // // pressure contribution to traction for B rows: + w_B · 0.5*(-p_A n_A - p_B n_B)
    //         // // -> gives -0.5 w_B · (p_A n_A + p_B n_B) in the matrix
    //         // for (IndexType j = 0; j < n_patch_A; ++j) {
    //         //     rLeftHandSideMatrix(idxB(i, idim), idxA(j, dim)) -= weight * 0.5 * vB * normal_A[idim] * N_patch_A(0, j);
    //         // }
    //         // for (IndexType j = 0; j < n_patch_B; ++j) {
    //         //     rLeftHandSideMatrix(idxB(i, idim), idxB(j, dim)) -= weight * 0.5 * vB * normal_B[idim] * N_patch_B(0, j);
    //         // }
    //     }
    // }




    






    
    // noalias(rRightHandSideVector) = ZeroVector(size);
    // // Add analytical flux contribution on both sides (manufactured solution forcing)
    // // Using p_D = sin(x)*cos(y) and viscous traction t = 2*sym(grad u) * n (nu=1)
    // {
    //     GeometryType::JacobiansType J_patch_A;
    //     r_patch_A.Jacobian(J_patch_A, r_patch_A.GetDefaultIntegrationMethod());
    //     Matrix jacobian_patch_A = ZeroMatrix(3, 3);
    //     jacobian_patch_A(0, 0) = J_patch_A[0](0, 0);
    //     jacobian_patch_A(0, 1) = J_patch_A[0](0, 1);
    //     jacobian_patch_A(1, 0) = J_patch_A[0](1, 0);
    //     jacobian_patch_A(1, 1) = J_patch_A[0](1, 1);
    //     jacobian_patch_A(2, 2) = 1.0;
    //     array_1d<double, 3> tangent_patch_A;
    //     r_patch_A.Calculate(LOCAL_TANGENT, tangent_patch_A);
    //     Vector determinant_factor_patch_A = prod(jacobian_patch_A, tangent_patch_A);
    //     determinant_factor_patch_A[2] = 0.0;
    //     const double detJ_patch_A = norm_2(determinant_factor_patch_A);
    //     KRATOS_ERROR_IF(detJ_patch_A <= std::numeric_limits<double>::epsilon())
    //         << "FluidCouplingCondition: degenerate Jacobian for patch A (determinant ~ 0)." << std::endl;
    //     const auto& integration_points_patch_A = r_patch_A.IntegrationPoints();
    //     const double weight = integration_points_patch_A[0].Weight() * detJ_patch_A;
        


    //     const auto center = r_patch_A.Center();
    //     const double x = center.X();
    //     const double y = center.Y();
    //     const double p_D = 0; //std::sin(x) * std::cos(y);

    //     // grad u and symmetric part (2D)
    //     Matrix grad_u(dim, dim, 0.0);
    //     if (dim >= 2) {
    //         // grad_u(0,0) = std::sinh(x) * std::sinh(y);
    //         // grad_u(0,1) = std::cosh(x) * std::cosh(y);
    //         // grad_u(1,0) = -std::cosh(x) * std::cosh(y);
    //         // grad_u(1,1) = -std::sinh(x) * std::sinh(y);
    //         grad_u(0,0) = 0;
    //         grad_u(0,1) = 1;
    //         grad_u(1,0) = 1;
    //         grad_u(1,1) = 0;
    //     }
    //     // Side A
    //     {
    //         Matrix sym_grad_u(dim, dim, 0.0);
    //         for (IndexType i = 0; i < dim; ++i)
    //             for (IndexType j = 0; j < dim; ++j)
    //                 sym_grad_u(i,j) = 0.5 * (grad_u(i,j) + grad_u(j,i));

    //         array_1d<double,3> nA = mNormalPhysicalSpaceA;
    //         array_1d<double,2> normal_A; normal_A[0]=nA[0]; normal_A[1]=nA[1];
    //         Vector tN(dim, 0.0);
    //         for (IndexType i = 0; i < dim; ++i)
    //             for (IndexType j = 0; j < dim; ++j)
    //                 tN[i] += 2.0 * sym_grad_u(i,j) * normal_A[j];

    //         const Matrix& N_A = r_patch_A.ShapeFunctionsValues();
    //         for (IndexType j = 0; j < n_patch_A; ++j) {
    //             for (IndexType idim = 0; idim < dim; ++idim) {
    //                 // pressure term: - < v·n, p_D >
    //                 rRightHandSideVector(idxA(j, idim)) -= p_D * (N_A(0,j) * normal_A[idim]) * weight;
    //                 // viscous traction term: + < v, 2*sym(grad u) n >
    //                 rRightHandSideVector(idxA(j, idim)) += N_A(0,j) * tN[idim] * weight;
    //             }
    //         }
    //     }

    //     // Side B
    //     {
    //         Matrix sym_grad_u(dim, dim, 0.0);
    //         for (IndexType i = 0; i < dim; ++i)
    //             for (IndexType j = 0; j < dim; ++j)
    //                 sym_grad_u(i,j) = 0.5 * (grad_u(i,j) + grad_u(j,i));

    //         array_1d<double,3> nB = mNormalPhysicalSpaceB;
    //         array_1d<double,2> normal_B; normal_B[0]=nB[0]; normal_B[1]=nB[1];
    //         Vector tN(dim, 0.0);
    //         for (IndexType i = 0; i < dim; ++i)
    //             for (IndexType j = 0; j < dim; ++j)
    //                 tN[i] += 2.0 * sym_grad_u(i,j) * normal_B[j];

    //         const Matrix& N_B = r_patch_B.ShapeFunctionsValues();
    //         for (IndexType j = 0; j < n_patch_B; ++j) {
    //             for (IndexType idim = 0; idim < dim; ++idim) {
    //                 rRightHandSideVector(idxB(j, idim)) -= p_D * (N_B(0,j) * normal_B[idim]) * weight;
    //                 rRightHandSideVector(idxB(j, idim)) += N_B(0,j) * tN[idim] * weight;
    //             }
    //         }
    //     }
    // }


}

void FluidCouplingCondition::InitializeMaterial()
{
    KRATOS_TRY
    if (GetProperties().Has(CONSTITUTIVE_LAW) && GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial(r_properties, r_geometry, row(N_values, 0));
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }
    KRATOS_CATCH("")
}

void FluidCouplingCondition::ApplyConstitutiveLaw(
        const GeometryType& rGeometry,
        const Matrix& rB, 
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables) const
{
    const SizeType number_of_nodes = rGeometry.size();

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    Vector old_displacement(number_of_nodes * mDim, 0.0);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double,3>& v = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
        const IndexType base = i * mDim;
        old_displacement[base + 0] = v[0];
        if (mDim > 1) old_displacement[base + 1] = v[1];
        if (mDim > 2) old_displacement[base + 2] = v[2];
    }
    Vector old_strain = prod(rB,old_displacement);
    rValues.SetStrainVector(old_strain);
    rValues.SetStressVector(rConstitutiveVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveVariables.D);

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void FluidCouplingCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometryA = GetGeometry();
    const auto& r_geometryB = GetGeometryMirror();
    const SizeType n_patch_A = r_geometryA.PointsNumber();
    const SizeType n_patch_B = r_geometryB.PointsNumber();

    // Determine velocity dimension from shape gradients (consistent with LHS)
    const auto& rDN_De_A = r_geometryA.ShapeFunctionsLocalGradients(r_geometryA.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();

    const SizeType dofs_per_node = dim + 1; // dim velocities + pressure
    const SizeType total_size = dofs_per_node * (n_patch_A + n_patch_B);

    if (rResult.size() != total_size) {
        rResult.resize(total_size, false);
    }

    SizeType k = 0;
    // Validate DOFs on A side
    for (IndexType i = 0; i < n_patch_A; ++i) {
        KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(VELOCITY_X)) << "Node missing VELOCITY_X DOF in A, node id: " << r_geometryA[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(VELOCITY_Y)) << "Node missing VELOCITY_Y DOF in A, node id: " << r_geometryA[i].Id() << std::endl;
        if (dim == 3) KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(VELOCITY_Z)) << "Node missing VELOCITY_Z DOF in A, node id: " << r_geometryA[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryA[i].HasDofFor(PRESSURE))   << "Node missing PRESSURE DOF in A, node id: "   << r_geometryA[i].Id() << std::endl;
    }
    // Validate DOFs on B side
    for (IndexType i = 0; i < n_patch_B; ++i) {
        KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(VELOCITY_X)) << "Node missing VELOCITY_X DOF in B, node id: " << r_geometryB[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(VELOCITY_Y)) << "Node missing VELOCITY_Y DOF in B, node id: " << r_geometryB[i].Id() << std::endl;
        if (dim == 3) KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(VELOCITY_Z)) << "Node missing VELOCITY_Z DOF in B, node id: " << r_geometryB[i].Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometryB[i].HasDofFor(PRESSURE))   << "Node missing PRESSURE DOF in B, node id: "   << r_geometryB[i].Id() << std::endl;
    }

    // Fill EIDs and track min/max for debugging
    int min_eid = std::numeric_limits<int>::max();
    int max_eid = std::numeric_limits<int>::lowest();
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const int ex = r_geometryA[i].GetDof(VELOCITY_X).EquationId();   min_eid = std::min(min_eid, ex); max_eid = std::max(max_eid, ex); rResult[k++] = ex;
        const int ey = r_geometryA[i].GetDof(VELOCITY_Y).EquationId();   min_eid = std::min(min_eid, ey); max_eid = std::max(max_eid, ey); rResult[k++] = ey;
        if (dim == 3) {
            const int ez = r_geometryA[i].GetDof(VELOCITY_Z).EquationId(); min_eid = std::min(min_eid, ez); max_eid = std::max(max_eid, ez); rResult[k++] = ez;
        }
        const int ep = r_geometryA[i].GetDof(PRESSURE).EquationId();     min_eid = std::min(min_eid, ep); max_eid = std::max(max_eid, ep); rResult[k++] = ep;
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        const int ex = r_geometryB[i].GetDof(VELOCITY_X).EquationId();   min_eid = std::min(min_eid, ex); max_eid = std::max(max_eid, ex); rResult[k++] = ex;
        const int ey = r_geometryB[i].GetDof(VELOCITY_Y).EquationId();   min_eid = std::min(min_eid, ey); max_eid = std::max(max_eid, ey); rResult[k++] = ey;
        if (dim == 3) {
            const int ez = r_geometryB[i].GetDof(VELOCITY_Z).EquationId(); min_eid = std::min(min_eid, ez); max_eid = std::max(max_eid, ez); rResult[k++] = ez;
        }
        const int ep = r_geometryB[i].GetDof(PRESSURE).EquationId();     min_eid = std::min(min_eid, ep); max_eid = std::max(max_eid, ep); rResult[k++] = ep;
    }

}

void FluidCouplingCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_patch_A = GetGeometry();
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const SizeType n_patch_B = r_patch_B.PointsNumber();

    // Determine velocity dimension from shape gradients (consistent with LHS)
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();

    rElementalDofList.resize(0);
    rElementalDofList.reserve((dim + 1) * (n_patch_A + n_patch_B));

    for (IndexType i = 0; i < n_patch_A; ++i) {
        rElementalDofList.push_back(r_patch_A[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(r_patch_A[i].pGetDof(VELOCITY_Y));
        if (dim == 3) rElementalDofList.push_back(r_patch_A[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(r_patch_A[i].pGetDof(PRESSURE));
    }
    for (IndexType i = 0; i < n_patch_B; ++i) {
        rElementalDofList.push_back(r_patch_B[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(r_patch_B[i].pGetDof(VELOCITY_Y));
        if (dim == 3) rElementalDofList.push_back(r_patch_B[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(r_patch_B[i].pGetDof(PRESSURE));
    }
}

int FluidCouplingCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    return Condition::Check(rCurrentProcessInfo);
}

void FluidCouplingCondition::GetVelocityCoefficientVectorA(
    Vector& rValues) const
{
    const auto& r_patch_A = GetGeometry();
    const SizeType n_patch_A = r_patch_A.PointsNumber();
    const auto& rDN_De_A = r_patch_A.ShapeFunctionsLocalGradients(r_patch_A.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_A.empty() ? 2 : rDN_De_A[0].size2();
    const SizeType nvals = n_patch_A * dim;
    if (rValues.size() != nvals) rValues.resize(nvals, false);
    IndexType k = 0;
    for (IndexType i = 0; i < n_patch_A; ++i) {
        const auto& v = r_patch_A[i].FastGetSolutionStepValue(VELOCITY);
        rValues[k++] = v[0];
        rValues[k++] = v[1];
        if (dim == 3) rValues[k++] = v[2];
    }
}

void FluidCouplingCondition::GetVelocityCoefficientVectorB(
    Vector& rValues) const
{
    const auto& r_patch_B = GetGeometryMirror();
    const SizeType n_patch_B = r_patch_B.PointsNumber();
    const auto& rDN_De_B = r_patch_B.ShapeFunctionsLocalGradients(r_patch_B.GetDefaultIntegrationMethod());
    const SizeType dim = rDN_De_B.empty() ? 2 : rDN_De_B[0].size2();
    const SizeType nvals = n_patch_B * dim;
    if (rValues.size() != nvals) rValues.resize(nvals, false);
    IndexType k = 0;
    for (IndexType i = 0; i < n_patch_B; ++i) {
        const auto& v = r_patch_B[i].FastGetSolutionStepValue(VELOCITY);
        rValues[k++] = v[0];
        rValues[k++] = v[1];
        if (dim == 3) rValues[k++] = v[2];
    }
}

// Builds the 3x(2*n) Voigt B-matrix from shape derivatives r_DN_DX (n x 2)
void FluidCouplingCondition::CalculateB(Matrix& rB, const ShapeDerivativesType& r_DN_DX) const
{
    const SizeType number_of_control_points = r_DN_DX.size1();
    const SizeType dim = r_DN_DX.size2();
    KRATOS_ERROR_IF(dim < 2) << "FluidCouplingCondition::CalculateB expects 2D derivatives (size2>=2)." << std::endl;

    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    if (rB.size1() != 3 || rB.size2() != mat_size)
        rB.resize(3, mat_size, false);

    noalias(rB) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        // epsilon_xx
        rB(0, 2 * i)     = r_DN_DX(i, 0); // dNi/dx
        // epsilon_yy
        rB(1, 2 * i + 1) = r_DN_DX(i, 1); // dNi/dy
        // // epsilon_xy (symmetrized)
        rB(2, 2 * i)     = r_DN_DX(i, 1); // dNi/dy
        rB(2, 2 * i + 1) = r_DN_DX(i, 0); // dNi/dx



    }
}

} // namespace Kratos
