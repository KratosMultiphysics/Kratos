// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "utilities/atomic_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"
#include "includes/cfd_variables.h" //FIXME: For the OSS_SWITCH. Move to a more suitable place

// Application includes
#include "custom_elements/small_displacement_mixed_volumetric_strain_oss_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedVolumetricStrainOssElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssElement>(
        NewId,
        GetGeometry().Create(ThisNodes),
        pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainOssElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssElement>(
        NewId,
        pGeom,
        pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainOssElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedVolumetricStrainOssElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssElement>(
        NewId,
        GetGeometry().Create(rThisNodes),
        pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // const auto& r_geometry = GetGeometry();
    // const SizeType n_nodes = r_geometry.PointsNumber();
    // const SizeType dim = r_geometry.WorkingSpaceDimension();
    // const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // // Check if the Orthogonal SubScales (OSS) are active
    // const bool oss_switch = rCurrentProcessInfo.Has(OSS_SWITCH) ? rCurrentProcessInfo[OSS_SWITCH] : false;

    // // Get dynamic data if required
    // const double density = mIsDynamic ? GetProperties()[DENSITY] : 0.0;
    // const double dt = mIsDynamic ? rCurrentProcessInfo[DELTA_TIME] : 0.0;

    // // Create the kinematics container and fill the nodal data
    // KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    // for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //     const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //     for (IndexType d = 0; d < dim; ++d) {
    //         kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //     }
    //     kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    // }

    // // Set the constitutive law values
    // ConstitutiveVariables constitutive_variables(strain_size);
    // ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    // auto& r_cons_law_options = cons_law_values.GetOptions();
    // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // // Calculate and update the Gauss points dynamic subscale values
    // Vector m_T;
    // Vector voigt_identity;
    // if (mIsDynamic) {
    //     voigt_identity =  ZeroVector(strain_size);
    //     for (IndexType d = 0; d < dim; ++d) {
    //         voigt_identity[d] = 1.0;
    //     }
    //     m_T = prod(voigt_identity, mAnisotropyTensor);
    // }

    // Vector aux_u_s;
    // array_1d<double,3> u_proj_gauss;           
    // for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
    //     // Recompute the kinematics
    //     CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //     // Set the constitutive variables
    //     CalculateConstitutiveVariables(
    //         kinematic_variables,
    //         constitutive_variables,
    //         cons_law_values,
    //         i_gauss,r_integration_points,
    //         ConstitutiveLaw::StressMeasure_Cauchy);

    //     // Call the constitutive law to update material variables
    //     mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);

    //     // Update the dynamic subscale values
    //     if (mIsDynamic) {
    //         // Get current Gauss point data
    //         const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
    //         const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
    //         const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rCurrentProcessInfo);
    //         const Vector grad_eps = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);

    //         // Calculate the n+1 subscale with current database
    //         aux_u_s = (density / std::pow(dt,2)) * (2.0*mDisplacementSubscale1[i_gauss] - mDisplacementSubscale2[i_gauss]);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             aux_u_s[d] += body_force[d];
    //             aux_u_s[d] += bulk_modulus * grad_eps[d];
    //         }

    //         // Substract the projection in case the OSS is active
    //         if (oss_switch) {
    //             // Calculate the displacement residual projection at current Gauss point
    //             noalias(u_proj_gauss) = ZeroVector(3);
    //             for (IndexType i = 0; i < n_nodes; ++i) {
    //                 const double N_i = kinematic_variables.N[i];
    //                 // const auto& r_u_proj_i = r_geometry[i].GetValue(DISPLACEMENT_PROJECTION);
    //                 const auto& r_u_proj_i = GetNodeDisplacementProjectionValue(i);
    //                 noalias(u_proj_gauss) += N_i * r_u_proj_i;
    //             }
    //             // Substract the nodal projection to get the orthogonal subscale value
    //             for (IndexType d = 0; d < dim; ++d) {
    //                 aux_u_s[d] -= u_proj_gauss[d];
    //             }
    //         }

    //         // Multiply by the stabilization parameter
    //         aux_u_s *= tau_1;

    //         // Update subscale vector values
    //         mDisplacementSubscale2[i_gauss] = mDisplacementSubscale1[i_gauss]; // n-1 subscale value update
    //         mDisplacementSubscale1[i_gauss] = aux_u_s; // n subscale value update
    //     }
    // }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // const auto& r_geometry = GetGeometry();
    // const SizeType dim = r_geometry.WorkingSpaceDimension();
    // const SizeType n_nodes = r_geometry.PointsNumber();
    // const SizeType block_size = dim + 1;
    // const SizeType matrix_size = block_size * n_nodes;
    // const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // // Check if the Orthogonal SubScales (OSS) are active
    // const bool oss_switch = rCurrentProcessInfo.Has(OSS_SWITCH) ? rCurrentProcessInfo[OSS_SWITCH] : false;

    // // Get dynamic data if required
    // const double density = mIsDynamic ? GetProperties()[DENSITY] : 0.0;
    // const double dt = mIsDynamic ? rCurrentProcessInfo[DELTA_TIME] : 0.0;

    // // Check RHS size
    // if (rRightHandSideVector.size() != matrix_size) {
    //     rRightHandSideVector.resize(matrix_size, false);
    // }

    // // Check LHS size
    // if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
    //     rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    // }

    // // Create the kinematics container and fill the nodal data
    // KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    // for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //     const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //     for (IndexType d = 0; d < dim; ++d) {
    //         kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //     }
    //     kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    // }

    // // Create the constitutive variables and values containers
    // ConstitutiveVariables constitutive_variables(strain_size);
    // ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    // auto& r_cons_law_options = cons_law_values.GetOptions();
    // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // // Set auxiliary arrays
    // Vector voigt_identity =  ZeroVector(strain_size);
    // for (IndexType d = 0; d < dim; ++d) {
    //     voigt_identity[d] = 1.0;
    // }
    // double rhs_mass_i;
    // Vector rhs_mom_i(dim);
    // Matrix lhs_uu_ij(dim,dim);
    // Vector lhs_ue_ij(dim);
    // Vector lhs_eu_ij(dim);
    // double lhs_ee_ij;
    // Vector G_i(dim);
    // Vector G_j(dim);
    // Vector psi_i(dim);
    // Vector psi_j(dim);
    // Vector transBi_C_invT_m(dim);
    // Matrix B_i(strain_size, dim);
    // Matrix B_j(strain_size, dim);

    // // Calculate the anisotropy tensor products
    // const Vector m_T = prod(voigt_identity, mAnisotropyTensor);
    // const Vector invT_m = prod(mInverseAnisotropyTensor, voigt_identity);

    // // Calculate the RHS and LHS contributions
    // rLeftHandSideMatrix.clear();
    // rRightHandSideVector.clear();

    // const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    // const double thickness = (dim == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    // for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
    //     // Calculate kinematics
    //     CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
    //     const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

    //     // Calculate the constitutive response
    //     CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

    //     // Calculate tau_1 stabilization constant
    //     const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rCurrentProcessInfo);

    //     // Calculate tau_2 stabilization constant
    //     const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
    //     const double shear_modulus = CalculateShearModulus(constitutive_variables.D);
    //     const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

    //     // Calculate and add the LHS contributions
    //     const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
    //     const Vector C_invT_m_voigt = prod(cons_law_values.GetConstitutiveMatrix(), invT_m);
    //     const Vector grad_eps = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);
    //     const Vector tau_1_body = tau_1 * body_force;

    //     // Calculate the dynamic subscale terms for current Gauss point
    //     Vector aux_u_s_acc;
    //     if (mIsDynamic) {
    //         aux_u_s_acc = (density / std::pow(dt, 2)) *(2.0 * mDisplacementSubscale1[i_gauss] - mDisplacementSubscale2[i_gauss]);
    //     }

    //     for (IndexType i = 0; i < n_nodes; ++i) {
    //         const double N_i = kinematic_variables.N[i];
    //         noalias(G_i) = row(kinematic_variables.DN_DX, i);
    //         for (IndexType k = 0;  k < strain_size; ++k) {
    //             for (IndexType l = 0; l < dim; ++l) {
    //                 B_i(k,l) = kinematic_variables.B(k, i * dim + l);
    //             }
    //         }
    //         noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));
    //         noalias(transBi_C_invT_m) = prod(trans(B_i), C_invT_m_voigt);

    //         // Add momentum body force RHS contribution
    //         noalias(rhs_mom_i) = w_gauss * N_i * body_force;
    //         // Add momentum internal force RHS contribution
    //         // Note that this includes both the deviatoric and volumetric internal force RHS contributions
    //         noalias(rhs_mom_i) -= w_gauss * prod(trans(B_i), cons_law_values.GetStressVector());
    //         // Add the divergence mass stabilization term (grad(eps_vol) x grad(eps_vol)) to the RHS
    //         rhs_mass_i = w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, grad_eps);
    //         // Add the divergence mass stabilization term (grad(eps_vol) x body_force) to the RHS
    //         // rhs_mass_i -= w_gauss * bulk_modulus * inner_prod(G_i, tau_1_body);
    //         rhs_mass_i += w_gauss * bulk_modulus * inner_prod(G_i, tau_1_body); //TODO: Cross-check this sign with Ramon and Ricc
            
    //         // Add the dynamic subscale terms to the RHS
    //         if (mIsDynamic) {
    //             noalias(rhs_mom_i) += w_gauss * (1.0 - tau_1 * density / std::pow(dt, 2)) * N_i * aux_u_s_acc;
    //             noalias(rhs_mom_i) -= w_gauss * (density / std::pow(dt, 2)) * N_i * tau_1_body;
    //             noalias(rhs_mom_i) -= w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * grad_eps;
    //             rhs_mass_i += w_gauss * bulk_modulus * tau_1 * inner_prod(G_i, aux_u_s_acc);
    //         }

    //         // Add the OSS projection terms
    //         if (oss_switch) {
    //             // Get OSS projection values at current Gauss point
    //             double eps_proj_gauss = 0.0;
    //             Vector u_proj_gauss = ZeroVector(dim);
    //             for (IndexType j = 0; j < n_nodes; ++j) {
    //                 const double N_j = kinematic_variables.N[j];
    //                 // const auto& r_u_s_proj = r_geometry[j].GetValue(DISPLACEMENT_PROJECTION);
    //                 // eps_proj_gauss += N_j * r_geometry[j].GetValue(VOLUMETRIC_STRAIN_PROJECTION);
    //                 const auto& r_u_s_proj = GetNodeDisplacementProjectionValue(j);
    //                 eps_proj_gauss += N_j * GetNodeVolumetricStrainProjectionValue(j);
    //                 for (IndexType d = 0; d < dim; ++d) {
    //                     u_proj_gauss[d] += N_j * r_u_s_proj[d];
    //                 }
    //             }
    //             // Substract the projection values to current residual
    //             rhs_mass_i -= w_gauss * bulk_modulus * tau_1 * inner_prod(G_i, u_proj_gauss);
    //             rhs_mass_i -= w_gauss * N_i * bulk_modulus * tau_2 * eps_proj_gauss;
    //             noalias(rhs_mom_i) -= w_gauss * bulk_modulus * tau_2 * psi_i * eps_proj_gauss;
    //             if (mIsDynamic) {
    //                 noalias(rhs_mom_i) += w_gauss * tau_1 * (density / std::pow(dt, 2)) * N_i * u_proj_gauss;
    //             }
    //         }

    //         for (IndexType j = 0; j < n_nodes; ++j) {
    //             const double N_j = kinematic_variables.N[j];
    //             noalias(G_j) = row(kinematic_variables.DN_DX, j);
    //             for (IndexType k = 0;  k < strain_size; ++k) {
    //                 for (IndexType l = 0; l < dim; ++l) {
    //                     B_j(k,l) = kinematic_variables.B(k, j * dim + l);
    //                 }
    //             }
    //             noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_j)));

    //             // Add the extra terms in the RHS momentum equation
    //             Vector disp_j(dim);
    //             for (unsigned int d = 0; d < dim; ++d) {
    //                 disp_j[d] = kinematic_variables.Displacements[j * dim + d];
    //             }
    //             const double aux_vol_j = (N_j * kinematic_variables.VolumetricNodalStrains[j] - inner_prod(psi_j, disp_j));
    //             noalias(rhs_mom_i) += (w_gauss / dim) * transBi_C_invT_m * aux_vol_j;
    //             noalias(rhs_mom_i) -= w_gauss * (1 - tau_2) * bulk_modulus * psi_i * aux_vol_j;

    //             rhs_mass_i += w_gauss * (1 - tau_2) * bulk_modulus * N_i * aux_vol_j;

    //             // Add momentum internal force LHS contribution
    //             noalias(lhs_uu_ij) = w_gauss * prod(trans(B_i), Matrix(prod(cons_law_values.GetConstitutiveMatrix(), B_j)));
    //             noalias(lhs_uu_ij) -= w_gauss * (1 - tau_2) * bulk_modulus * outer_prod(psi_i, trans(psi_j));
    //             // Add momentum volumetric strain LHS contribution
    //             noalias(lhs_ue_ij) = w_gauss * (1 - tau_2) * bulk_modulus * psi_i * N_j;
    //             // Add mass conservation displacement divergence LHS contribution
    //             noalias(lhs_eu_ij) = w_gauss * (1 - tau_2) * bulk_modulus * N_i * trans(psi_j);
    //             // Add mass conservation volumetric strain LHS contribution
    //             lhs_ee_ij = - w_gauss * (1 - tau_2) * bulk_modulus * N_i * N_j;
    //             lhs_ee_ij -= w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, G_j);
    //             // Add the dynamic subscale terms to the LHS
    //             if (mIsDynamic) {
    //                 noalias(lhs_ue_ij) += w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * G_j;
    //             }

    //             // Assemble the LHS contributions
    //             rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim) += lhs_ee_ij;
    //             for (IndexType d = 0; d < dim; ++d) {
    //                 rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += lhs_ue_ij[d];
    //                 rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) += lhs_eu_ij[d];
    //                 for (IndexType d2 = 0; d2 < dim; ++d2) {
    //                     rLeftHandSideMatrix(i * block_size + d, j * block_size + d2) += lhs_uu_ij(d, d2);
    //                 }
    //             }
    //         }

    //         // Assemble the RHS contributions
    //         rRightHandSideVector[i * block_size + dim] += rhs_mass_i;
    //         for (IndexType d = 0; d < dim; ++d) {
    //             rRightHandSideVector[i * block_size + d] += rhs_mom_i[d];
    //         }
    //     }
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // const auto& r_geometry = GetGeometry();
    // const SizeType dim = r_geometry.WorkingSpaceDimension();
    // const SizeType n_nodes = r_geometry.PointsNumber();
    // const SizeType block_size = dim + 1;
    // const SizeType matrix_size = block_size * n_nodes;
    // const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // // Check if the Orthogonal SubScales (OSS) are active
    // const bool oss_switch = rCurrentProcessInfo.Has(OSS_SWITCH) ? rCurrentProcessInfo[OSS_SWITCH] : false;

    // // Get dynamic data if required
    // const double density = mIsDynamic ? GetProperties()[DENSITY] : 0.0;
    // const double dt = mIsDynamic ? rCurrentProcessInfo[DELTA_TIME] : 0.0;

    // // Check RHS size
    // if (rRightHandSideVector.size() != matrix_size) {
    //     rRightHandSideVector.resize(matrix_size, false);
    // }

    // // Create the kinematics container and fill the nodal data
    // KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    // for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //     const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //     for (IndexType d = 0; d < dim; ++d) {
    //         kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //     }
    //     kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    // }

    // // Create the constitutive variables and values containers
    // ConstitutiveVariables constitutive_variables(strain_size);
    // ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    // auto& r_cons_law_options = cons_law_values.GetOptions();
    // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // // Set auxiliary arrays
    // Vector voigt_identity =  ZeroVector(strain_size);
    // for (IndexType d = 0; d < dim; ++d) {
    //     voigt_identity[d] = 1.0;
    // }
    // double rhs_mass_i;
    // Vector rhs_mom_i(dim);
    // Vector G_i(dim);
    // Vector G_j(dim);
    // Vector psi_i(dim);
    // Vector psi_j(dim);
    // Vector transBi_C_invT_m(dim);
    // Matrix B_i(strain_size, dim);
    // Matrix B_j(strain_size, dim);

    // // Calculate the anisotropy tensor products
    // const Vector m_T = prod(voigt_identity, mAnisotropyTensor);
    // const Vector invT_m = prod(mInverseAnisotropyTensor, voigt_identity);

    // // Calculate the RHS contributions
    // rRightHandSideVector.clear();

    // const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    // const double thickness = (dim == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    // for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
    //     // Calculate kinematics
    //     CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
    //     const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

    //     // Calculate the constitutive response
    //     CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

    //     // Calculate tau_1 stabilization constant
    //     const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rCurrentProcessInfo);

    //     // Calculate tau_2 stabilization constant
    //     const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
    //     const double shear_modulus = CalculateShearModulus(constitutive_variables.D);
    //     const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

    //     // Calculate and add the LHS contributions
    //     const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
    //     const Vector C_invT_m_voigt = prod(cons_law_values.GetConstitutiveMatrix(), invT_m);
    //     const Vector grad_eps = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);
    //     const Vector tau_1_body = tau_1 * body_force;

    //     // Calculate the dynamic subscale terms for current Gauss point
    //     Vector aux_u_s_acc;
    //     if (mIsDynamic) {
    //         aux_u_s_acc = (density / std::pow(dt, 2)) *(2.0 * mDisplacementSubscale1[i_gauss] - mDisplacementSubscale2[i_gauss]);
    //     }

    //     for (IndexType i = 0; i < n_nodes; ++i) {
    //         const double N_i = kinematic_variables.N[i];
    //         noalias(G_i) = row(kinematic_variables.DN_DX, i);
    //         for (IndexType k = 0;  k < strain_size; ++k) {
    //             for (IndexType l = 0; l < dim; ++l) {
    //                 B_i(k,l) = kinematic_variables.B(k, i * dim + l);
    //             }
    //         }
    //         noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));
    //         noalias(transBi_C_invT_m) = prod(trans(B_i), C_invT_m_voigt);

    //         // Add momentum body force RHS contribution
    //         noalias(rhs_mom_i) = w_gauss * N_i * body_force;
    //         // Add momentum internal force RHS contribution
    //         // Note that this includes both the deviatoric and volumetric internal force RHS contributions
    //         noalias(rhs_mom_i) -= w_gauss * prod(trans(B_i), cons_law_values.GetStressVector());
    //         // Add the divergence mass stabilization term (grad(eps_vol) x grad(eps_vol)) to the RHS
    //         rhs_mass_i = w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, grad_eps);
    //         // Add the divergence mass stabilization term (grad(eps_vol) x body_force) to the RHS
    //         // rhs_mass_i -= w_gauss * bulk_modulus * inner_prod(G_i, tau_1_body);
    //         rhs_mass_i += w_gauss * bulk_modulus * inner_prod(G_i, tau_1_body); // TODO: Cross-check this sign with Ramon and Ricc
            
    //         // Add the dynamic subscale terms to the RHS
    //         if (mIsDynamic) {
    //             if (!oss_switch) {
    //                 noalias(rhs_mom_i) += w_gauss * (1.0 - tau_1 * density / std::pow(dt, 2)) * N_i * aux_u_s_acc;
    //                 noalias(rhs_mom_i) -= w_gauss * (density / std::pow(dt, 2)) * N_i * tau_1_body;
    //                 noalias(rhs_mom_i) -= w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * grad_eps;
    //             }
    //             rhs_mass_i += w_gauss * bulk_modulus * tau_1 * inner_prod(G_i, aux_u_s_acc);
    //         }

    //         // Add the OSS projection terms
    //         if (oss_switch) {
    //             // Get OSS projection values at current Gauss point
    //             double eps_proj_gauss = 0.0;
    //             Vector u_proj_gauss = ZeroVector(dim);
    //             for (IndexType j = 0; j < n_nodes; ++j) {
    //                 const double N_j = kinematic_variables.N[j];
    //                 // const auto& r_u_s_proj = r_geometry[j].GetValue(DISPLACEMENT_PROJECTION);
    //                 // eps_proj_gauss += N_j * r_geometry[j].GetValue(VOLUMETRIC_STRAIN_PROJECTION);
    //                 const auto& r_u_s_proj = GetNodeDisplacementProjectionValue(j);
    //                 eps_proj_gauss += N_j * GetNodeVolumetricStrainProjectionValue(j);
    //                 for (IndexType d = 0; d < dim; ++d) {
    //                     u_proj_gauss[d] += N_j * r_u_s_proj[d];
    //                 }
    //             }
    //             // Substract the projection values to current residual
    //             rhs_mass_i -= w_gauss * bulk_modulus * tau_1 * inner_prod(G_i, u_proj_gauss);
    //             rhs_mass_i -= w_gauss * N_i * bulk_modulus * tau_2 * eps_proj_gauss;
    //             noalias(rhs_mom_i) -= w_gauss * bulk_modulus * tau_2 * psi_i * eps_proj_gauss;
    //         }

    //         for (IndexType j = 0; j < n_nodes; ++j) {
    //             const double N_j = kinematic_variables.N[j];
    //             for (IndexType k = 0;  k < strain_size; ++k) {
    //                 for (IndexType l = 0; l < dim; ++l) {
    //                     B_j(k,l) = kinematic_variables.B(k, j * dim + l);
    //                 }
    //             }
    //             noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_j)));

    //             // Add the extra terms in the RHS momentum equation
    //             Vector disp_j(dim);
    //             for (unsigned int d = 0; d < dim; ++d) {
    //                 disp_j[d] = kinematic_variables.Displacements[j * dim + d];
    //             }
    //             const double aux_vol_j = (N_j * kinematic_variables.VolumetricNodalStrains[j] - inner_prod(psi_j, disp_j));
    //             noalias(rhs_mom_i) += (w_gauss / dim) * transBi_C_invT_m * aux_vol_j;
    //             noalias(rhs_mom_i) -= w_gauss * (1 - tau_2) * bulk_modulus * psi_i * aux_vol_j;

    //             rhs_mass_i += w_gauss * (1 - tau_2) * bulk_modulus * N_i * aux_vol_j;
    //         }

    //         // Assemble the RHS contributions
    //         rRightHandSideVector[i * block_size + dim] += rhs_mass_i;
    //         for (IndexType d = 0; d < dim; ++d) {
    //             rRightHandSideVector[i * block_size + d] += rhs_mom_i[d];
    //         }
    //     }
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    // if (rVariable == VOLUMETRIC_STRAIN_PROJECTION) {
    //     // Get geometry data
    //     auto& r_geometry = GetGeometry();
    //     const auto& r_prop = GetProperties();
    //     const SizeType n_nodes = r_geometry.PointsNumber();
    //     const SizeType dim = r_geometry.WorkingSpaceDimension();
    //     const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    //     const SizeType n_gauss = r_integration_points.size();
    //     const SizeType strain_size = r_prop.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    //     // Calculate the mass matrix (lumped or consistent)
    //     const double density = r_prop[DENSITY]; //TODO: Leave this for the weighted projection
    //     const double thickness = (dim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    //     // Create the kinematics container and fill the nodal data
    //     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //         }
    //         kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     }

    //     // Calculate the kinematic constraint residual projection
    //     Vector voigt_identity = ZeroVector(strain_size);
    //     for (IndexType d = 0; d < dim; ++d) {
    //         voigt_identity[d] = 1.0;
    //     }
    //     Vector psi_j(dim);
    //     Vector aux_proj = ZeroVector(n_nodes);
    //     Matrix B_j(strain_size, dim);
    //     for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
    //         // Calculate kinematics
    //         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //         // Calculate current integration point weight
    //         // Note that this already includes the thickness
    //         const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

    //         // Calculate current Gauss point displacement subscale projection contribution
    //         double aux_res;
    //         for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //             const double N_i = kinematic_variables.N[i_node];
    //             for (IndexType j_node = 0; j_node < n_nodes; ++j_node) {
    //                 for (IndexType k = 0;  k < strain_size; ++k) {
    //                     for (IndexType l = 0; l < dim; ++l) {
    //                         B_j(k,l) = kinematic_variables.B(k, j_node * dim + l);
    //                     }
    //                 }
    //                 noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_j)));

    //                 aux_res = 0.0;
    //                 const double N_j = kinematic_variables.N[j_node];
    //                 for (IndexType d = 0; d < dim; ++d) {
    //                     aux_res += psi_j[d] * kinematic_variables.Displacements[j_node*dim + d];
    //                 }
    //                 aux_res -= N_j * kinematic_variables.VolumetricNodalStrains[j_node];
    //             }
    //             aux_proj[i_node] += w_gauss * N_i * aux_res;
    //         }
    //     }

    //     // Nodal assembly of the projection contributions
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         // AtomicAdd(r_geometry[i_node].GetValue(VOLUMETRIC_STRAIN_PROJECTION), aux_proj[i_node]);
    //         AtomicAdd(GetNodeVolumetricStrainProjectionValue(i_node), aux_proj[i_node]);
    //     }

    // } else {
    //     KRATOS_ERROR << "Variable not implemented." << std::endl;
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    // if (rVariable == DISPLACEMENT_PROJECTION) {
    //     // Get geometry data
    //     auto& r_geometry = GetGeometry();
    //     const auto& r_prop = GetProperties();
    //     const SizeType n_nodes = r_geometry.PointsNumber();
    //     const SizeType dim = r_geometry.WorkingSpaceDimension();
    //     const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    //     const SizeType n_gauss = r_integration_points.size();
    //     const SizeType strain_size = r_prop.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    //     // Calculate the mass matrix (lumped or consistent)
    //     const double density = r_prop[DENSITY]; //TODO: Leave this for the weighted projection
    //     const double thickness = (dim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    //     // Create the kinematics container and fill the nodal data
    //     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //         }
    //         kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     }

    //     // Set te constitutive law values
    //     ConstitutiveVariables constitutive_variables(strain_size);
    //     ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_prop, rCurrentProcessInfo);
    //     auto &r_cons_law_options = cons_law_values.GetOptions();
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    //     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    //     // Calculate the momentum equation projection
    //     Vector G_j(dim);
    //     Vector grad_eps(dim);
    //     Vector aux_proj = ZeroVector(n_nodes*dim);
    //     for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
    //         // Calculate kinematics
    //         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //         // Calculate the constitutive response
    //         CalculateConstitutiveVariables(
    //             kinematic_variables,
    //             constitutive_variables,
    //             cons_law_values,
    //             i_gauss,
    //             r_integration_points,
    //             ConstitutiveLaw::StressMeasure_Cauchy);

    //         // Calculate bulk modulus
    //         const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);

    //         // Calculate current integration point weight
    //         // Note that this already includes the thickness
    //         const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

    //         // Get current Gauss point body force value
    //         const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);

    //         // Calculate current Gauss point displacement subscale projection contribution
    //         for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //             noalias(grad_eps) = ZeroVector(dim);
    //             const double N_i = kinematic_variables.N[i_node];
    //             for (IndexType j_node = 0; j_node < n_nodes; ++j_node) {
    //                 noalias(G_j) = row(kinematic_variables.DN_DX, j_node);
    //                 const double eps_j = kinematic_variables.VolumetricNodalStrains[j_node];
    //                 noalias(grad_eps) += G_j * eps_j;
    //             }

    //             for (IndexType d = 0; d < dim; ++d) {
    //                 aux_proj[i_node*dim + d] += w_gauss * N_i * (body_force[d] + bulk_modulus * grad_eps[d]);
    //             }
    //         }
    //     }

    //     // Nodal assembly of the projection contributions
    //     array_1d<double,3> i_aux_proj;
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         i_aux_proj = ZeroVector(3);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             i_aux_proj[d] = aux_proj[i_node*dim + d];
    //         }
    //         // AtomicAdd(r_geometry[i_node].GetValue(DISPLACEMENT_PROJECTION), i_aux_proj);
    //         AtomicAdd(GetNodeDisplacementProjectionValue(i_node), i_aux_proj);
    //     }

    // } else {
    //     KRATOS_ERROR << "Variable not implemented." << std::endl;
    // }
}

// /***********************************************************************************/
// /***********************************************************************************/

// void SmallDisplacementMixedVolumetricStrainOssElement::CalculateOrthogonalSubScalesOperator(
//     MatrixType &rOrthogonalSubScalesOperator,
//     const ProcessInfo& rProcessInfo) const
// {
//     const auto &r_geometry = GetGeometry();
//     const SizeType dim = r_geometry.WorkingSpaceDimension();
//     const SizeType n_nodes = r_geometry.PointsNumber();
//     const SizeType block_size = dim + 1;
//     const SizeType matrix_size = block_size * n_nodes;
//     const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

//     // Check LHS size
//     if (rOrthogonalSubScalesOperator.size1() != matrix_size || rOrthogonalSubScalesOperator.size2() != matrix_size) {
//         rOrthogonalSubScalesOperator.resize(matrix_size, matrix_size, false);
//     }

//     // Initialize the OSS operator with zeros
//     rOrthogonalSubScalesOperator.clear();

//     // Create the kinematics container
//     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

//     // Create the constitutive variables and values containers
//     ConstitutiveVariables constitutive_variables(strain_size);
//     ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rProcessInfo);
//     auto& r_cons_law_options = cons_law_values.GetOptions();
//     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
//     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
//     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

//     // Set auxiliary arrays
//     Vector voigt_identity =  ZeroVector(strain_size);
//     for (IndexType d = 0; d < dim; ++d) {
//         voigt_identity[d] = 1.0;
//     }
//     Vector G_i(dim);
//     Vector psi_i(dim);
//     Matrix B_i(strain_size, dim);

//     // Calculate the anisotropy tensor products
//     const Vector m_T = prod(voigt_identity, mAnisotropyTensor);

//     const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
//     const double thickness = (dim == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
//     const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
//     for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
//         // Calculate kinematics
//         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
//         const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

//         // Calculate the constitutive response
//         CalculateConstitutiveVariables(
//             kinematic_variables,
//             constitutive_variables,
//             cons_law_values,
//             i_gauss,
//             r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
//             ConstitutiveLaw::StressMeasure_Cauchy);

//         // Calculate tau_1 stabilization constant
//         const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rProcessInfo);

//         // Calculate tau_2 stabilization constant
//         const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
//         const double shear_modulus = CalculateShearModulus(constitutive_variables.D);
//         const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

//         const double aux_w_kappa_tau_1 = w_gauss * bulk_modulus * tau_1;
//         const double aux_w_kappa_tau_2 = w_gauss * bulk_modulus * tau_2;

//         for (IndexType i = 0; i < n_nodes; ++i) {
//             const double N_i = kinematic_variables.N[i];
//             noalias(G_i) = row(kinematic_variables.DN_DX, i);
//             for (IndexType k = 0;  k < strain_size; ++k) {
//                 for (IndexType l = 0; l < dim; ++l) {
//                     B_i(k,l) = kinematic_variables.B(k, i * dim + l);
//                 }
//             }
//             noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));

//             for (IndexType j = 0; j < n_nodes; ++j) {
//                 const double N_j = kinematic_variables.N[j];
//                 rOrthogonalSubScalesOperator(i*block_size + dim, j*block_size + dim) += aux_w_kappa_tau_2 * N_i * N_j;
//                 for (IndexType d = 0; d < dim; ++d) {
//                     rOrthogonalSubScalesOperator(i*block_size + dim, j*block_size + d) += aux_w_kappa_tau_1 * G_i[d] * N_j;
//                     rOrthogonalSubScalesOperator(i*block_size + d, j*block_size + dim) += aux_w_kappa_tau_2 * psi_i[d] * N_j;
//                 }
//             }
//         }
//     }
// }

// /***********************************************************************************/
// /***********************************************************************************/

// void SmallDisplacementMixedVolumetricStrainOssElement::CalculateOrthogonalSubScalesLumpedProjectionOperator(
//         MatrixType& rOrthogonalSubScalesLumpedProjectionOperator,
//         const ProcessInfo& rProcessInfo) const
// {
//     const auto &r_geometry = GetGeometry();
//     const SizeType dim = r_geometry.WorkingSpaceDimension();
//     const SizeType n_nodes = r_geometry.PointsNumber();
//     const SizeType block_size = dim + 1;
//     const SizeType matrix_size = block_size * n_nodes;
//     const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

//     // Check LHS size
//     if (rOrthogonalSubScalesLumpedProjectionOperator.size1() != matrix_size || rOrthogonalSubScalesLumpedProjectionOperator.size2() != matrix_size) {
//         rOrthogonalSubScalesLumpedProjectionOperator.resize(matrix_size, matrix_size, false);
//     }

//     // Initialize the OSS operator with zeros
//     rOrthogonalSubScalesLumpedProjectionOperator.clear();

//     // Create the kinematics container
//     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

//     // Create the constitutive variables and values containers
//     ConstitutiveVariables constitutive_variables(strain_size);
//     ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rProcessInfo);
//     auto& r_cons_law_options = cons_law_values.GetOptions();
//     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
//     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
//     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

//     // Set auxiliary arrays
//     Vector voigt_identity =  ZeroVector(strain_size);
//     for (IndexType d = 0; d < dim; ++d) {
//         voigt_identity[d] = 1.0;
//     }

//     // Calculate the anisotropy tensor products
//     const Vector m_T = prod(voigt_identity, mAnisotropyTensor);

//     const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
//     const double thickness = (dim == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
//     const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
//     for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
//         // Calculate kinematics
//         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
//         const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

//         // Calculate the constitutive response
//         CalculateConstitutiveVariables(
//             kinematic_variables,
//             constitutive_variables,
//             cons_law_values,
//             i_gauss,
//             r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
//             ConstitutiveLaw::StressMeasure_Cauchy);

//         // Calculate tau_1 stabilization constant
//         const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rProcessInfo);

//         // Calculate tau_2 stabilization constant
//         const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
//         const double shear_modulus = CalculateShearModulus(constitutive_variables.D);
//         const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

//         double sum_N_j;
//         const double aux_w_kappa_tau_1 = w_gauss * bulk_modulus * tau_1;
//         const double aux_w_kappa_tau_2 = w_gauss * bulk_modulus * tau_2;
//         for (IndexType i = 0; i < n_nodes; ++i) {
//             sum_N_j = 0.0;
//             const double N_i = kinematic_variables.N[i];
//             for (IndexType j = 0; j < n_nodes; ++j) {
//                 sum_N_j += kinematic_variables.N[j];
//             }
//             for (IndexType d = 0; d < dim; ++d) {
//                 rOrthogonalSubScalesLumpedProjectionOperator(i*block_size + d, i*block_size + d) -= aux_w_kappa_tau_1 * N_i * sum_N_j;
//             }
//             rOrthogonalSubScalesLumpedProjectionOperator(i*block_size + dim, i*block_size + dim) -= aux_w_kappa_tau_2 * N_i * sum_N_j;
//         }
//     }
// }

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedVolumetricStrainOssElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Call the base element check
    int check = BaseType::Check(rCurrentProcessInfo);

    // Check that the element also includes the OSS projection variables
    const auto& r_geometry = this->GetGeometry();
    for (const auto& r_node : r_geometry) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT_PROJECTION,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUMETRIC_STRAIN_PROJECTION,r_node)
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    // const auto& r_geometry = GetGeometry();
    // const SizeType n_nodes = r_geometry.PointsNumber();
    // const SizeType dim = r_geometry.WorkingSpaceDimension();
    // const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // // Check if the Orthogonal SubScales (OSS) are active
    // const bool oss_switch = rCurrentProcessInfo.Has(OSS_SWITCH) ? rCurrentProcessInfo[OSS_SWITCH] : false;

    // // Resize output container
    // const SizeType n_gauss = r_integration_points.size();
    // if (rOutput.size() != n_gauss) {
    //     rOutput.resize(n_gauss);
    // }

    // if (mConstitutiveLawVector[0]->Has(rVariable)) {
    //     GetValueOnConstitutiveLaw(rVariable, rOutput);
    // } else if (rVariable == VOLUMETRIC_STRAIN_SUBSCALE) {
    //     // Create the kinematics container and fill the nodal data
    //     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //         }
    //         kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     }

    //     // Set the constitutive law values
    //     ConstitutiveVariables constitutive_variables(strain_size);
    //     ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    //     auto& r_cons_law_options = cons_law_values.GetOptions();
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    //     // Calculate the Gauss points volumetric strain subscale values
    //     Vector voigt_identity = ZeroVector(strain_size);
    //     for (IndexType d = 0; d < dim; ++d) {
    //         voigt_identity[d] = 1.0;
    //     }

    //     Vector psi_i(dim);
    //     Vector disp_i(dim);
    //     Matrix B_i(strain_size, dim);
    //     for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
    //         // Initialize current value
    //         rOutput[i_gauss] = 0.0;

    //         // Recompute the kinematics
    //         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //         // Set the constitutive variables
    //         CalculateConstitutiveVariables(
    //             kinematic_variables,
    //             constitutive_variables,
    //             cons_law_values,
    //             i_gauss,r_integration_points,
    //             ConstitutiveLaw::StressMeasure_Cauchy);

    //         // Calculate tau_2 stabilization constant
    //         const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
    //         const double shear_modulus = CalculateShearModulus(constitutive_variables.D);
    //         const double tau_2 = std::min(1.0e-2, 4.0 * shear_modulus / bulk_modulus);

    //         // Calculate current gauss point subscale value
    //         for (IndexType i = 0; i < n_nodes; ++i) {
    //             for (IndexType d = 0; d < dim; ++d) {
    //                 for (IndexType k = 0;  k < strain_size; ++k) {
    //                     B_i(k,d) = kinematic_variables.B(k, i * dim + d);
    //                 }
    //                 disp_i[d] = kinematic_variables.Displacements(i*n_nodes + d);
    //             }
    //             noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));
    //             const double eps_vol_i = kinematic_variables.N[i] * kinematic_variables.VolumetricNodalStrains[i];
    //             rOutput[i_gauss] += tau_2 * (inner_prod(psi_i, disp_i) - eps_vol_i);
    //             if (oss_switch) {
    //                 // const double eps_vol_proj_i = kinematic_variables.N[i] * r_geometry[i].GetValue(VOLUMETRIC_STRAIN_PROJECTION);
    //                 const double eps_vol_proj_i = kinematic_variables.N[i] * GetNodeVolumetricStrainProjectionValue(i);
    //                 rOutput[i_gauss] -= tau_2 * eps_vol_proj_i;
    //             }
    //         }
    //     }
    // } else {
    //     CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    // const auto& r_geometry = GetGeometry();
    // const SizeType n_nodes = r_geometry.PointsNumber();
    // const SizeType dim = r_geometry.WorkingSpaceDimension();
    // const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // // Check if the Orthogonal SubScales (OSS) are active
    // const bool oss_switch = rCurrentProcessInfo.Has(OSS_SWITCH) ? rCurrentProcessInfo[OSS_SWITCH] : false;

    // // Get dynamic data if required
    // const double density = mIsDynamic ? GetProperties()[DENSITY] : 0.0;
    // const double dt = mIsDynamic ? rCurrentProcessInfo[DELTA_TIME] : 0.0;

    // // Resize output container
    // const SizeType n_gauss = r_integration_points.size();
    // if (rOutput.size() != n_gauss) {
    //     rOutput.resize(n_gauss);
    // }

    // if (mConstitutiveLawVector[0]->Has(rVariable)) {
    //     GetValueOnConstitutiveLaw(rVariable, rOutput);
    // } else if (rVariable == DISPLACEMENT_SUBSCALE) {
    //     // Create the kinematics container and fill the nodal data
    //     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //         }
    //         kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     }

    //     // Set the constitutive law values
    //     ConstitutiveVariables constitutive_variables(strain_size);
    //     ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    //     auto& r_cons_law_options = cons_law_values.GetOptions();
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    //     // Calculate the Gauss points displacement subscale values
    //     Vector m_T;
    //     Vector voigt_identity;
    //     if (mIsDynamic) {
    //         voigt_identity =  ZeroVector(strain_size);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             voigt_identity[d] = 1.0;
    //         }
    //         m_T = prod(voigt_identity, mAnisotropyTensor);
    //     }

    //     Vector aux_u_s(dim);
    //     Vector grad_eps(dim);
    //     array_1d<double, 3> body_force;
    //     array_1d<double, 3> u_proj_gauss;
    //     for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
    //         // Recompute the kinematics
    //         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //         // Set the constitutive variables
    //         CalculateConstitutiveVariables(
    //             kinematic_variables,
    //             constitutive_variables,
    //             cons_law_values,
    //             i_gauss,r_integration_points,
    //             ConstitutiveLaw::StressMeasure_Cauchy);

    //         // Calculate tau_1 stabilization constant
    //         const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rCurrentProcessInfo);

    //         // Get current Gauss point data
    //         const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
    //         noalias(body_force) = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
    //         noalias(grad_eps) = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);

    //         // If dynamic subscales are active, add the dynamic subscale component
    //         aux_u_s = bulk_modulus * grad_eps;
    //         for (IndexType d = 0; d < dim; ++d) {
    //             aux_u_s[d] += body_force[d];
    //         }
    //         if (mIsDynamic) {
    //             aux_u_s += (density / std::pow(dt,2)) * (2.0*mDisplacementSubscale1[i_gauss] - mDisplacementSubscale2[i_gauss]);
    //         }

    //         // Substract the projection in case the OSS is active
    //         if (oss_switch) {
    //             // Calculate the displacement residual projection at current Gauss point
    //             noalias(u_proj_gauss) = ZeroVector(3);
    //             for (IndexType i = 0; i < n_nodes; ++i) {
    //                 const double N_i = kinematic_variables.N[i];
    //                 // const auto& r_u_proj_i = r_geometry[i].GetValue(DISPLACEMENT_PROJECTION);
    //                 const auto& r_u_proj_i = GetNodeDisplacementProjectionValue(i);
    //                 noalias(u_proj_gauss) += N_i * r_u_proj_i;
    //             }
    //             // Substract the nodal projection to get the orthogonal subscale value
    //             for (IndexType d = 0; d < dim; ++d) {
    //                 aux_u_s[d] -= u_proj_gauss[d];
    //             }
    //         }

    //         // Multiply by the stabilization parameter
    //         aux_u_s *= tau_1;

    //         // Save current value in the output container
    //         for (IndexType d = 0; d < dim; ++d) {
    //             rOutput[i_gauss][d] = aux_u_s[d];
    //         }
    //     }
    // } else {
    //     CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedVolumetricStrainOssElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","dynamic"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["CAUCHY_STRESS_VECTOR"],
            "nodal_historical"       : ["DISPLACEMENT","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION"],
            "nodal_non_historical"   : [],
            "entity"                 : ["DISPLACEMENT_SUBSCALE", "VOLUMETRIC_STRAIN_SUBSCALE"]
        },
        "required_variables"         : ["DISPLACEMENT","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION"],
        "required_dofs"              : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Quadrilateral2D4", "Tetrahedra3D4","Hexahedra3D8"],
        "element_integrates_in_time" : false,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","PlaneStress","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements a mixed displacement - volumetric strain formulation with Variational MultiScales (VMS) stabilization. This formulation is capable to deal with materials in the incompressible limit as well as with anisotropy."
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","VOLUMETRIC_STRAIN"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

} // Namespace Kratos
