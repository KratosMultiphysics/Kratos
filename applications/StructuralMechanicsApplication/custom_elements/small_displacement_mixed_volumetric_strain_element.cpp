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
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedVolumetricStrainElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void SmallDisplacementMixedVolumetricStrainElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType dof_size = n_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType eps_vol_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN);

    IndexType aux_index = 0;
    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
    } else {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType dof_size  = n_nodes*(dim+1);

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    } else if (dim == 3) {
        for(IndexType i = 0; i < n_nodes; ++i){
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * (dim + 1) + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::Initialize()
{
    KRATOS_TRY

    // Integration method initialization
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
    const auto& r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Constitutive Law Vector initialisation
    if (mConstitutiveLawVector.size() != r_integration_points.size()) {
        mConstitutiveLawVector.resize(r_integration_points.size());
    }

    // Initialize material
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Set te constitutive law values
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector strain(strain_size);
    Vector stress(strain_size);
    Matrix cons_matrix(strain_size, strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    cons_law_values.SetStrainVector(strain);
    cons_law_values.SetStressVector(stress);
    cons_law_values.SetConstitutiveMatrix(cons_matrix);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Set te constitutive law values
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector strain(strain_size);
    Vector stress(strain_size);
    Matrix cons_matrix(strain_size, strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    cons_law_values.SetStrainVector(strain);
    cons_law_values.SetStressVector(stress);
    cons_law_values.SetConstitutiveMatrix(cons_matrix);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set auxiliary arrays
    Vector voigt_identity =  ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the linearised bulk modulus
        const double bulk_modulus = CalculateLinearisedBulkModulus(kinematic_variables, constitutive_variables);

        // Calculate stabilization constants
        const double h = CalculateElementSize(kinematic_variables);
        const double shear_modulus = CalculateLinearisedShearModulus(kinematic_variables, constitutive_variables);
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

        // Calculate the volumetric residual
        double div_u_gauss = 0.0;
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            for (IndexType d = 0; d < dim; ++d) {
                div_u_gauss += kinematic_variables.DN_DX(i_node, d) * kinematic_variables.Displacements(i_node * dim + d);
            }
        }
        const double eps_vol_gauss = inner_prod(kinematic_variables.N, kinematic_variables.VolumetricNodalStrains);
        const double vol_residual =  eps_vol_gauss - div_u_gauss;

        // Get force vectors
        const auto internal_force = prod(trans(kinematic_variables.B), cons_law_values.GetStressVector());
        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);

        // Calculate and add the LHS contributions
        const double w_alpha = w_gauss / dim;
        const double w_tau_1_k = w_gauss * tau_1 * bulk_modulus;
        const double w_1_tau_2_k =  w_gauss * (1 - tau_2) * bulk_modulus;

        const Matrix transB_C = prod(trans(kinematic_variables.B), cons_law_values.GetConstitutiveMatrix());
        const Matrix transB_C_B = prod(transB_C, kinematic_variables.B);

        const Vector C_m_voigt = prod(cons_law_values.GetConstitutiveMatrix(), voigt_identity);
        const Matrix C_m = MathUtils<double>::StressVectorToTensor(C_m_voigt); // Expansion from Voigt to tensor notation
        const Vector transB_C_m = prod(trans(kinematic_variables.B), C_m_voigt);
        const Vector grad_eps = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);
        const Vector C_m_grad_eps = prod(C_m, grad_eps);

        for (IndexType i = 0; i < n_nodes; ++i) {
            // Mass conservation volumetric residual RHS term
            // (already includes mass conservation volumetric strain contribution (eps_vol x eps_vol))
            // (already includes mass conservation divergence contribution (eps_vol x div(u)))
            double& r_rhs_mass_row = rRightHandSideVector[i * block_size + dim];
            r_rhs_mass_row += w_1_tau_2_k * kinematic_variables.N[i] * vol_residual;

            const Vector G_I = row(kinematic_variables.DN_DX, i);
            for (IndexType d = 0; d < dim; ++d) {
                double& r_rhs_mom_row = rRightHandSideVector[i * block_size + d];
                // Add momentum body force RHS contribution
                r_rhs_mom_row += w_gauss * kinematic_variables.N[i] * body_force[d];
                // Add momentum internal force RHS contribution
                // Note that this includes both the deviatoric and volumetric internal force RHS contributions
                r_rhs_mom_row -= w_gauss * internal_force(i * dim + d);
                // Add the extra terms in the RHS momentum equation
                r_rhs_mom_row += w_alpha * transB_C_m(i * dim + d) * vol_residual;
                r_rhs_mom_row -= w_1_tau_2_k * G_I(d) * vol_residual;
                // Add the divergence mass stabilization term (grad(eps_vol) x grad(eps_vol)) to the RHS
                r_rhs_mass_row += w_tau_1_k * kinematic_variables.DN_DX(i, d) * C_m_grad_eps(d);
                // Add the divergence mass stabilization term (grad(eps_vol) x body_force) to the RHS
                r_rhs_mass_row -= w_tau_1_k * kinematic_variables.DN_DX(i, d) * body_force[d];

                for (IndexType j = 0; j < n_nodes; ++j) {
                    const Vector G_J = row(kinematic_variables.DN_DX, j);
                    // Add mass conservation displacement divergence LHS contribution
                    rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) += w_1_tau_2_k * kinematic_variables.N[i] * G_J(d);
                    // Add momentum volumetric strain LHS contribution
                    rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += w_1_tau_2_k * G_I(d) * kinematic_variables.N(j);
                    // Add momentum internal force LHS contribution
                    for (IndexType d2 = 0; d2 < dim; ++d2) {
                        double& r_LHS_mom_mom = rLeftHandSideMatrix(i * block_size + d, j * block_size + d2);
                        r_LHS_mom_mom += w_gauss * transB_C_B(i * dim + d, j * dim + d2);
                        r_LHS_mom_mom -= w_1_tau_2_k * G_I(d) * G_J(d2);
                    }
                }
            }

            // Add mass conservation volumetric strain LHS contribution
            const Vector G_I_C_m = prod(G_I, C_m);
            for (IndexType j = 0; j < n_nodes; ++j) {
                const Vector G_J = row(kinematic_variables.DN_DX, j);
                double& r_LHS_mass_mass = rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim);
                r_LHS_mass_mass -= w_1_tau_2_k * kinematic_variables.N[i] * kinematic_variables.N[j];
                r_LHS_mass_mass -= w_tau_1_k * inner_prod(G_I_C_m, G_J);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set auxiliary arrays
    Vector voigt_identity =  ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }

    // Calculate the LHS contributions
    rLeftHandSideMatrix.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the linearised bulk modulus
        const double bulk_modulus = CalculateLinearisedBulkModulus(kinematic_variables, constitutive_variables);

        // Calculate stabilization constants
        const double h = CalculateElementSize(kinematic_variables);
        const double shear_modulus = CalculateLinearisedShearModulus(kinematic_variables, constitutive_variables);
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

        // Calculate and add the LHS contributions
        const double w_tau_1_k = w_gauss * tau_1 * bulk_modulus;
        const double w_1_tau_2_k =  w_gauss * (1 - tau_2) * bulk_modulus;

        const Matrix transB_C = prod(trans(kinematic_variables.B), cons_law_values.GetConstitutiveMatrix());
        const Matrix transB_C_B = prod(transB_C, kinematic_variables.B);

        const Vector C_m_voigt = prod(cons_law_values.GetConstitutiveMatrix(), voigt_identity);
        const Matrix C_m = MathUtils<double>::StressVectorToTensor(C_m_voigt);

        for (IndexType i = 0; i < n_nodes; ++i) {
            const Vector G_I = row(kinematic_variables.DN_DX, i);
            const Vector G_I_C_m = prod(G_I, C_m);
            for (IndexType j = 0; j < n_nodes; ++j) {
                const Vector G_J = row(kinematic_variables.DN_DX, j);
                const Matrix G_I_transG_J = outer_prod(G_I, G_J);
                // Add mass conservation volumetric strain contribution
                double& r_LHS_mass_mass = rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim);
                r_LHS_mass_mass -= w_1_tau_2_k * kinematic_variables.N[i] * kinematic_variables.N[j];
                r_LHS_mass_mass -= w_tau_1_k * inner_prod(G_I_C_m, G_J);
                for (IndexType d = 0; d < dim; ++d) {
                    // Add mass conservation displacement divergence contribution
                    rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) += w_1_tau_2_k * kinematic_variables.N[i] * G_J(d);
                    // Add momentum internal force contribution
                    for (IndexType d2 = 0; d2 < dim; ++d2) {
                        double& r_LHS_mom_mom = rLeftHandSideMatrix(i * block_size + d, j * block_size + d2);
                        r_LHS_mom_mom += w_gauss * transB_C_B(i * dim + d, j * dim + d2);
                        r_LHS_mom_mom -= w_1_tau_2_k * G_I_transG_J(d, d2);
                    }
                    // Add momentum volumetric strain contribution
                    rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += w_1_tau_2_k * G_I(d) * kinematic_variables.N(j);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set auxiliary arrays
    Vector voigt_identity =  ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }

    // Calculate and add the RHS contributions
    rRightHandSideVector.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the linearised bulk modulus
        const double bulk_modulus = CalculateLinearisedBulkModulus(kinematic_variables, constitutive_variables);

        // Calculate stabilization constants
        const double h = CalculateElementSize(kinematic_variables);
        const double shear_modulus = CalculateLinearisedShearModulus(kinematic_variables, constitutive_variables);
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = std::min(1.0e-2, 4.0*shear_modulus/bulk_modulus);

        // Calculate the volumetric residual
        double div_u_gauss = 0.0;
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            for (IndexType d = 0; d < dim; ++d) {
                div_u_gauss += kinematic_variables.DN_DX(i_node, d) * kinematic_variables.Displacements(i_node * dim + d);
            }
        }
        const double eps_vol_gauss = inner_prod(kinematic_variables.N, kinematic_variables.VolumetricNodalStrains);
        const double vol_residual =  eps_vol_gauss - div_u_gauss;

        // Get force vectors
        const auto internal_force = prod(trans(kinematic_variables.B), cons_law_values.GetStressVector());
        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);

        // Add the RHS contributions
        const double w_alpha = w_gauss / dim;
        const double w_tau_1_k = w_gauss * tau_1 * bulk_modulus;
        const double w_1_tau_2_k =  w_gauss * (1 - tau_2) * bulk_modulus;

        const Vector C_m_voigt = prod(cons_law_values.GetConstitutiveMatrix(), voigt_identity);
        const Matrix C_m = MathUtils<double>::StressVectorToTensor(C_m_voigt); // Expansion from Voigt to tensor notation
        const Matrix transB_C = prod(trans(kinematic_variables.B), cons_law_values.GetConstitutiveMatrix());
        const Vector transB_C_m = prod(trans(kinematic_variables.B), C_m_voigt);
        const Vector grad_eps = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);
        const Vector C_m_grad_eps = prod(C_m, grad_eps);

        for (IndexType i = 0; i < n_nodes; ++i) {
            // Mass conservation equation terms (more efficient)
            // (already includes mass conservation volumetric strain contribution (eps_vol x eps_vol))
            // (already includes mass conservation divergence contribution (eps_vol x div(u)))
            double& r_rhs_mass_row = rRightHandSideVector[i * block_size + dim];
            r_rhs_mass_row += w_1_tau_2_k * kinematic_variables.N[i] * vol_residual;

            const Vector G_I = row(kinematic_variables.DN_DX, i);
            // Momentum equation terms
            for (IndexType d = 0; d < dim; ++d) {
                double& r_rhs_mom_row = rRightHandSideVector[i * block_size + d];
                // Add momentum body force RHS contribution
                r_rhs_mom_row += w_gauss * kinematic_variables.N[i] * body_force[d];
                // Add momentum internal force RHS contribution
                // Note that this includes both the deviatoric and volumetric internal force RHS contributions
                r_rhs_mom_row -= w_gauss * internal_force(i * dim + d);
                // Add the extra terms in the RHS momentum equation
                r_rhs_mom_row += w_alpha * transB_C_m(i * dim + d) * vol_residual;
                r_rhs_mom_row -= w_1_tau_2_k * G_I(d) * vol_residual;
                // Add the divergence mass stabilization term (grad(eps_vol) x grad(eps_vol)) to the RHS
                r_rhs_mass_row += w_tau_1_k * kinematic_variables.DN_DX(i, d) * C_m_grad_eps(d);
                // Add the divergence mass stabilization term (grad(eps_vol) x body_force) to the RHS
                r_rhs_mass_row -= w_tau_1_k * kinematic_variables.DN_DX(i, d) * body_force[d];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::InitializeMaterial()
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry = GetGeometry();
        const auto& r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        IndexType aux = 0;
        for (auto &it_gauss_pt : mConstitutiveLawVector) {
            it_gauss_pt = (r_properties[CONSTITUTIVE_LAW])->Clone();
            (it_gauss_pt)->InitializeMaterial(r_properties, r_geometry, row(r_N_values, aux));
            aux++;
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementMixedVolumetricStrainElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetStrainVector(rThisKinematicVariables.EquivalentStrain); // equivalent total strain
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // assuming that det(F) is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); // assuming that F is computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> SmallDisplacementMixedVolumetricStrainElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    array_1d<double, 3> body_force;
    for (IndexType i = 0; i < 3; ++i) {
        body_force[i] = 0.0;
    }

    const auto& r_properties = GetProperties();
    const double density = r_properties.Has(DENSITY) ? r_properties[DENSITY] : 0.0;

    if (r_properties.Has(VOLUME_ACCELERATION)) {
        noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];
    }

    const auto& r_geometry = GetGeometry();
    if(r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION)) {
        Vector N;
        N = r_geometry.ShapeFunctionsValues(N, rIntegrationPoints[PointNumber].Coordinates());
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        }
    }

    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // Calculate the inverse Jacobian
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geometry,
        r_integration_points[PointNumber],
        rThisKinematicVariables.J0);
    MathUtils<double>::InvertMatrix(
        rThisKinematicVariables.J0,
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.detJ0);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
        << "Element ID: " << this->Id() << " is inverted. det(J0) = " << rThisKinematicVariables.detJ0 << std::endl;

    // Calculate the shape functions gradients
    GeometryUtils::ShapeFunctionsGradients(
        r_geometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[PointNumber],
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.DN_DX);

    // Compute B
    CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.DN_DX);

    // Calculate the equivalent total strain
    CalculateEquivalentStrain(rThisKinematicVariables);

    // Compute equivalent F
    ComputeEquivalentF(rThisKinematicVariables.F, rThisKinematicVariables.EquivalentStrain);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*2    ) = rDN_DX(i, 0);
            rB(1, i*2 + 1) = rDN_DX(i, 1);
            rB(2, i*2    ) = rDN_DX(i, 1);
            rB(2, i*2 + 1) = rDN_DX(i, 0);
        }
    } else if(dimension == 3) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*3    ) = rDN_DX(i, 0);
            rB(1, i*3 + 1) = rDN_DX(i, 1);
            rB(2, i*3 + 2) = rDN_DX(i, 2);
            rB(3, i*3    ) = rDN_DX(i, 1);
            rB(3, i*3 + 1) = rDN_DX(i, 0);
            rB(4, i*3 + 1) = rDN_DX(i, 2);
            rB(4, i*3 + 2) = rDN_DX(i, 1);
            rB(5, i*3    ) = rDN_DX(i, 2);
            rB(5, i*3 + 2) = rDN_DX(i, 0);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    const auto& r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    // Add the deviatoric contribution to the equivalent strain
    const Vector total_strain = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    rThisKinematicVariables.EquivalentStrain = prod(rThisKinematicVariables.DevStrainOp, total_strain);

    // Interpolate and add the nodal volumetric strain contribution
    double gauss_vol_strain = 0.0;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        gauss_vol_strain += rThisKinematicVariables.N[i_node] * rThisKinematicVariables.VolumetricNodalStrains[i_node];
    }
    for (IndexType d = 0; d < dim; ++d) {
        rThisKinematicVariables.EquivalentStrain[d] += (1.0/dim) * gauss_vol_strain;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::ComputeEquivalentF(
    Matrix &rF,
    const Vector &rStrainTensor) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();

    if(dim == 2) {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(2);
        rF(1,0) = 0.5*rStrainTensor(2);
        rF(1,1) = 1.0+rStrainTensor(1);
    } else {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(3);
        rF(0,2) = 0.5*rStrainTensor(5);
        rF(1,0) = 0.5*rStrainTensor(3);
        rF(1,1) = 1.0+rStrainTensor(1);
        rF(1,2) = 0.5*rStrainTensor(4);
        rF(2,0) = 0.5*rStrainTensor(5);
        rF(2,1) = 0.5*rStrainTensor(4);
        rF(2,2) = 1.0+rStrainTensor(2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateElementSize(const KinematicVariables &rThisKinematicVariables) const
{
    const auto& r_geometry = GetGeometry();
    switch (r_geometry.GetGeometryType())
    {
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
        {
            const BoundedMatrix<double,3,2> aux_grads(rThisKinematicVariables.DN_DX);
            return ElementSizeCalculator<2,3>::GradientsElementSize(aux_grads);
        }
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
        {
            const BoundedMatrix<double,4,3> aux_grads(rThisKinematicVariables.DN_DX);
            return ElementSizeCalculator<3,4>::GradientsElementSize(aux_grads);
        }
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
        {
            return ElementSizeCalculator<2,4>::AverageElementSize(r_geometry);
        }
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:
        {
            return ElementSizeCalculator<3,8>::AverageElementSize(r_geometry);
        }
        default:
        {
            KRATOS_ERROR << "Asking for a non-implemented geometry element size calculation." << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateApproximatedBulkModulus(
    const ProcessInfo& rCurrentProcessInfo,
    const IndexType PointNumber,
    const Vector& rN
    ) const
{
    const auto& r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    // Calculate the bulk modulus with a fake volumetric strain field
    SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector aux_vol_stress_vect(strain_size);
    Vector aux_vol_strain_vect = ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        aux_vol_strain_vect[d] = 1.0;
    }

    // Call the constitutive law to get the material response of the fake volumetric strain field
    Matrix deformation_gradient(dim, dim);
    ConstitutiveLaw::Parameters aux_cons_law_values(r_geom, GetProperties(), rCurrentProcessInfo);
    auto& r_aux_cons_law_options = aux_cons_law_values.GetOptions();
    r_aux_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_aux_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_aux_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    aux_cons_law_values.SetShapeFunctionsValues(rN);
    aux_cons_law_values.SetStrainVector(aux_vol_strain_vect);
    aux_cons_law_values.SetStressVector(aux_vol_stress_vect);
    ComputeEquivalentF(deformation_gradient, aux_vol_strain_vect);
    aux_cons_law_values.SetDeformationGradientF(deformation_gradient);
    aux_cons_law_values.SetDeterminantF(MathUtils<double>::Det(deformation_gradient));
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(aux_cons_law_values);

    double aux_vol_strain = 0.0;
    double aux_vol_stress = 0.0;
    const double alpha = r_geom.WorkingSpaceDimension();
    for (IndexType d = 0; d < dim; ++d) {
        aux_vol_strain += aux_vol_strain_vect[d];
        aux_vol_stress += aux_cons_law_values.GetStressVector()[d];
    }
    aux_vol_stress /= alpha;

    return aux_vol_stress * aux_vol_strain / std::pow(aux_vol_strain, 2);
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateLinearisedBulkModulus(
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables) const
{
    const auto& r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    double bulk_modulus = 0.0;
    for (IndexType i = 0; i < dim; ++i) {
        for (IndexType j = 0; j < dim; ++j) {
            bulk_modulus += rThisConstitutiveVariables.D(i,j);
        }
    }

    return bulk_modulus / std::pow(dim, 2);
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateLinearisedShearModulus(
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables) const
{
    const auto& r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Get the volumetric constitutive tensor as the volumetric projection of the complete tensor
    // Note that this constitutive tensor has been previously computed using the complete effective strain
    const Matrix dev_cons_tensor = prod(rThisKinematicVariables.DevStrainOp, rThisConstitutiveVariables.D);

    // Apply a deviatoric strain perturbation and compute the stress using the deviatoric constitutive tensor
    const double eps_dev_strain = 1.0e-5;
    Vector eps_dev_strain_vect = ZeroVector(strain_size);
    for (IndexType i = dim; i < strain_size; ++i) {
        eps_dev_strain_vect[i] = eps_dev_strain;
    }
    const Vector eps_dev_stress_vect = prod(dev_cons_tensor, eps_dev_strain_vect);

    // Calculate and return the linearised shear modulus
    return inner_prod(eps_dev_stress_vect, eps_dev_strain_vect) / inner_prod(eps_dev_strain_vect, eps_dev_strain_vect);
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedVolumetricStrainElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int base_element_check = SmallDisplacementMixedVolumetricStrainElement::BaseType::Check(rCurrentProcessInfo);

    return base_element_check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
        // Create and initialize element variables:
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
            }
            kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Create the constitutive variables and values containers
        ConstitutiveVariables constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);


        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR) {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            } else {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
            }

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != strain_size) {
                rOutput[i_gauss].resize(strain_size, false);
            }
            rOutput[i_gauss] = constitutive_variables.StressVector;
        }
    } else {
        KRATOS_ERROR << "Integration point variable " + rVariable.Name() + " not implemented yet." << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("mConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

} // Namespace Kratos
