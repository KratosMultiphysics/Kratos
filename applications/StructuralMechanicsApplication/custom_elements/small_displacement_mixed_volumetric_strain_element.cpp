// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "utilities/atomic_utilities.h"
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

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
    const ProcessInfo& rCurrentProcessInfo) const
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
    const ProcessInfo& rCurrentProcessInfo) const
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

void SmallDisplacementMixedVolumetricStrainElement::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Integration method initialization
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // Constitutive Law Vector initialisation
        const SizeType n_gauss = r_integration_points.size();
        if (mConstitutiveLawVector.size() != n_gauss) {
            mConstitutiveLawVector.resize(n_gauss);
        }

        // Initialize material
        InitializeMaterial();

        // Calculate the and save the anisotropy transformation tensor
        CalculateAnisotropyTensor(rCurrentProcessInfo);
        CalculateInverseAnisotropyTensor();

        // Check if the problem is dynamic
        // This is required in order to add (or not) the non-acceleration dynamic terms
        unsigned int aux = 0;
        const auto& r_geom = this->GetGeometry();
        for (const auto& r_node : r_geom) {
            if (r_node.SolutionStepsDataHas(ACCELERATION)) {
                aux++;
            }
        }
        mIsDynamic = aux == r_geom.PointsNumber() ? true : false;

        // Initialize the displacement subscale containers
        if (mIsDynamic) {
            // Resize the arrays to the number of Gauss point (one subscale value per Gauss point)
            if (mDisplacementSubscale1.size() != n_gauss) {
                mDisplacementSubscale1.resize(n_gauss);
            }
            if (mDisplacementSubscale2.size() != n_gauss) {
                mDisplacementSubscale2.resize(n_gauss);
            }
            // Initialize to zero the subscale history
            const Vector aux_zero = ZeroVector(r_geom.WorkingSpaceDimension());
            for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                mDisplacementSubscale1[i_gauss] = aux_zero;
                mDisplacementSubscale2[i_gauss] = aux_zero;
            }
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set the constitutive law values
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Finalize the material response and update the Gauss points dynamic subscale values
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_integration_points,
            ConstitutiveLaw::StressMeasure_Cauchy);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);

        // Update the dynamic subscale values
        if (mIsDynamic) {
            UpdateGaussPointDisplacementSubscaleHistory(kinematic_variables, constitutive_variables, rCurrentProcessInfo, i_gauss);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
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

    // Set the auxiliary Gauss point variable container
    GaussPointAuxiliaryVariables gauss_point_auxiliary_variables(this, dim, strain_size);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the auxiliary Gauss point variables
        CalculateGaussPointAuxiliaryVariables(
            gauss_point_auxiliary_variables,
            kinematic_variables,
            constitutive_variables,
            rCurrentProcessInfo,
            i_gauss);

        // Assemble the current Gauss point contribution
        CalculateLocalSystemGaussPointContribution(
            rRightHandSideVector,
            rLeftHandSideMatrix,
            kinematic_variables,
            constitutive_variables,
            gauss_point_auxiliary_variables);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
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

    // Set the auxiliary Gauss point variable container
    GaussPointAuxiliaryVariables gauss_point_auxiliary_variables(this, dim, strain_size);

    // Calculate the LHS contributions
    rLeftHandSideMatrix.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the auxiliary Gauss point variables
        CalculateGaussPointAuxiliaryVariables(
            gauss_point_auxiliary_variables,
            kinematic_variables,
            constitutive_variables,
            rCurrentProcessInfo,
            i_gauss);

        // Assemble the current Gauss point contribution
        CalculateLeftHandSideGaussPointContribution(
            rLeftHandSideMatrix,
            kinematic_variables,
            constitutive_variables,
            gauss_point_auxiliary_variables);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
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

    // Set the auxiliary Gauss point variable container
    GaussPointAuxiliaryVariables gauss_point_auxiliary_variables(this, dim, strain_size);

    // Calculate the RHS contributions
    rRightHandSideVector.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the auxiliary Gauss point variables
        CalculateGaussPointAuxiliaryVariables(
            gauss_point_auxiliary_variables,
            kinematic_variables,
            constitutive_variables,
            rCurrentProcessInfo,
            i_gauss);

        // Assemble the current Gauss point contribution
        CalculateRightHandSideGaussPointContribution(
            rRightHandSideVector,
            kinematic_variables,
            constitutive_variables,
            gauss_point_auxiliary_variables);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLocalSystemGaussPointContribution(
    VectorType &rRightHandSideVector,
    MatrixType &rLeftHandSideMatrix,
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables,
    const GaussPointAuxiliaryVariables &rThisGaussPointAuxiliaryVariables) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    double rhs_mass_i;
    Vector rhs_mom_i(dim);
    Matrix lhs_uu_ij(dim,dim);
    Vector lhs_ue_ij(dim);
    Vector lhs_eu_ij(dim);
    double lhs_ee_ij;
    Vector G_i(dim);
    Vector G_j(dim);
    Vector psi_i(dim);
    Vector psi_j(dim);
    Vector transBi_C_invT_m(dim);
    Matrix B_i(strain_size, dim);
    Matrix B_j(strain_size, dim);

    const double dt = rThisGaussPointAuxiliaryVariables.Dt;
    const double density = rThisGaussPointAuxiliaryVariables.Rho;
    const double tau_1 = rThisGaussPointAuxiliaryVariables.Tau1;
    const double tau_2 = rThisGaussPointAuxiliaryVariables.Tau2;
    const double w_gauss = rThisGaussPointAuxiliaryVariables.Weight;
    const double bulk_modulus = rThisGaussPointAuxiliaryVariables.BulkModulus;
    const auto& voigt_identity = rThisGaussPointAuxiliaryVariables.m;
    const auto& body_force = rThisGaussPointAuxiliaryVariables.BodyForce;
    const auto& grad_eps = rThisGaussPointAuxiliaryVariables.VolumetricStrainGradient;
    const Vector C_invT_m_voigt = prod(rThisConstitutiveVariables.D, rThisGaussPointAuxiliaryVariables.invT_m);

    for (IndexType i = 0; i < n_nodes; ++i) {
        const double N_i = rThisKinematicVariables.N[i];
        noalias(G_i) = row(rThisKinematicVariables.DN_DX, i);
        for (IndexType k = 0;  k < strain_size; ++k) {
            for (IndexType l = 0; l < dim; ++l) {
                B_i(k,l) = rThisKinematicVariables.B(k, i * dim + l);
            }
        }
        noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));
        noalias(transBi_C_invT_m) = prod(trans(B_i), C_invT_m_voigt);

        // Add momentum body force RHS contribution
        noalias(rhs_mom_i) = w_gauss * N_i * body_force;
        // Add momentum internal force RHS contribution
        // Note that this includes both the deviatoric and volumetric internal force RHS contributions
        noalias(rhs_mom_i) -= w_gauss * prod(trans(B_i), rThisConstitutiveVariables.StressVector);
        // Add the divergence mass stabilization term (grad(eps_vol) x grad(eps_vol)) to the RHS
        rhs_mass_i = w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, grad_eps);
        // Add the divergence mass stabilization term (grad(eps_vol) x body_force) to the RHS
        // rhs_mass_i -= w_gauss * bulk_modulus * inner_prod(G_i, tau_1*body_force);
        rhs_mass_i += w_gauss * bulk_modulus * inner_prod(G_i, tau_1*body_force); //TODO: Cross-check this sign with Ramon and Ricc

        // Add the dynamic subscale terms to the RHS
        if (mIsDynamic) {
            const auto &aux_u_s_acc = rThisGaussPointAuxiliaryVariables.DisplacementSubscaleAcceleration;
            noalias(rhs_mom_i) += w_gauss * (1.0 - tau_1 * density / std::pow(dt, 2)) * N_i * aux_u_s_acc;
            noalias(rhs_mom_i) -= w_gauss * (density / std::pow(dt, 2)) * N_i * tau_1*body_force;
            noalias(rhs_mom_i) -= w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * grad_eps;
            rhs_mass_i += w_gauss * bulk_modulus * tau_1 * inner_prod(G_i, aux_u_s_acc);
        }

        for (IndexType j = 0; j < n_nodes; ++j) {
            const double N_j = rThisKinematicVariables.N[j];
            noalias(G_j) = row(rThisKinematicVariables.DN_DX, j);
            for (IndexType k = 0;  k < strain_size; ++k) {
                for (IndexType l = 0; l < dim; ++l) {
                    B_j(k,l) = rThisKinematicVariables.B(k, j * dim + l);
                }
            }
            noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_j)));

            // Add the extra terms in the RHS momentum equation
            Vector disp_j(dim);
            for (unsigned int d = 0; d < dim; ++d) {
                disp_j[d] = rThisKinematicVariables.Displacements[j * dim + d];
            }
            const double aux_vol_j = (N_j * rThisKinematicVariables.VolumetricNodalStrains[j] - inner_prod(psi_j, disp_j));
            noalias(rhs_mom_i) += (w_gauss / dim) * transBi_C_invT_m * aux_vol_j;
            noalias(rhs_mom_i) -= w_gauss * (1 - tau_2) * bulk_modulus * psi_i * aux_vol_j;

            rhs_mass_i += w_gauss * (1 - tau_2) * bulk_modulus * N_i * aux_vol_j;

            // Add momentum internal force LHS contribution
            noalias(lhs_uu_ij) = w_gauss * prod(trans(B_i), Matrix(prod(rThisConstitutiveVariables.D, B_j)));
            noalias(lhs_uu_ij) -= w_gauss * (1 - tau_2) * bulk_modulus * outer_prod(psi_i, trans(psi_j));
            // Add momentum volumetric strain LHS contribution
            noalias(lhs_ue_ij) = w_gauss * (1 - tau_2) * bulk_modulus * psi_i * N_j;
            // Add mass conservation displacement divergence LHS contribution
            noalias(lhs_eu_ij) = w_gauss * (1 - tau_2) * bulk_modulus * N_i * trans(psi_j);
            // Add mass conservation volumetric strain LHS contribution
            lhs_ee_ij = - w_gauss * (1 - tau_2) * bulk_modulus * N_i * N_j;
            lhs_ee_ij -= w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, G_j);
            // Add the dynamic subscale terms to the LHS
            if (mIsDynamic) {
                noalias(lhs_ue_ij) += w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * G_j;
            }

            // Assemble the LHS contributions
            rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim) += lhs_ee_ij;
            for (IndexType d = 0; d < dim; ++d) {
                rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += lhs_ue_ij[d];
                rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) += lhs_eu_ij[d];
                for (IndexType d2 = 0; d2 < dim; ++d2) {
                    rLeftHandSideMatrix(i * block_size + d, j * block_size + d2) += lhs_uu_ij(d, d2);
                }
            }
        }

        // Assemble the RHS contributions
        rRightHandSideVector[i * block_size + dim] += rhs_mass_i;
        for (IndexType d = 0; d < dim; ++d) {
            rRightHandSideVector[i * block_size + d] += rhs_mom_i[d];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateRightHandSideGaussPointContribution(
    VectorType &rRightHandSideVector,
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables,
    const GaussPointAuxiliaryVariables& rThisGaussPointAuxiliaryVariables) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    double rhs_mass_i;
    Vector rhs_mom_i(dim);
    Vector G_i(dim);
    Vector psi_i(dim);
    Vector psi_j(dim);
    Vector transBi_C_invT_m(dim);
    Matrix B_i(strain_size, dim);
    Matrix B_j(strain_size, dim);

    const double dt = rThisGaussPointAuxiliaryVariables.Dt;
    const double density = rThisGaussPointAuxiliaryVariables.Rho;
    const double tau_1 = rThisGaussPointAuxiliaryVariables.Tau1;
    const double tau_2 = rThisGaussPointAuxiliaryVariables.Tau2;
    const double w_gauss = rThisGaussPointAuxiliaryVariables.Weight;
    const double bulk_modulus = rThisGaussPointAuxiliaryVariables.BulkModulus;
    const auto& voigt_identity = rThisGaussPointAuxiliaryVariables.m;
    const auto& body_force = rThisGaussPointAuxiliaryVariables.BodyForce;
    const auto& grad_eps = rThisGaussPointAuxiliaryVariables.VolumetricStrainGradient;
    const Vector C_invT_m_voigt = prod(rThisConstitutiveVariables.D, rThisGaussPointAuxiliaryVariables.invT_m);

    for (IndexType i = 0; i < n_nodes; ++i) {
        const double N_i = rThisKinematicVariables.N[i];
        noalias(G_i) = row(rThisKinematicVariables.DN_DX, i);
        for (IndexType k = 0;  k < strain_size; ++k) {
            for (IndexType l = 0; l < dim; ++l) {
                B_i(k, l) = rThisKinematicVariables.B(k, i * dim + l);
            }
        }
        noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));
        noalias(transBi_C_invT_m) = prod(trans(B_i), C_invT_m_voigt);

        // Add momentum body force RHS contribution
        noalias(rhs_mom_i) = w_gauss * N_i * body_force;
        // Add momentum internal force RHS contribution
        // Note that this includes both the deviatoric and volumetric internal force RHS contributions
        noalias(rhs_mom_i) -= w_gauss * prod(trans(B_i), rThisConstitutiveVariables.StressVector);
        // Add the divergence mass stabilization term (grad(eps_vol) x grad(eps_vol)) to the RHS
        rhs_mass_i = w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, grad_eps);
        // Add the divergence mass stabilization term (grad(eps_vol) x body_force) to the RHS
        // rhs_mass_i -= w_gauss * bulk_modulus * inner_prod(G_i, tau_1*body_force);
        rhs_mass_i += w_gauss * bulk_modulus * inner_prod(G_i, tau_1*body_force); // TODO: Cross-check this sign with Ramon and Ricc

        // Add the dynamic subscale terms to the RHS
        if (mIsDynamic) {
            const auto& aux_u_s_acc = rThisGaussPointAuxiliaryVariables.DisplacementSubscaleAcceleration;
            noalias(rhs_mom_i) += w_gauss * (1.0 - tau_1 * density / std::pow(dt, 2)) * N_i * aux_u_s_acc;
            noalias(rhs_mom_i) -= w_gauss * (density / std::pow(dt, 2)) * N_i * tau_1*body_force;
            noalias(rhs_mom_i) -= w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * grad_eps;
            rhs_mass_i += w_gauss * bulk_modulus * tau_1 * inner_prod(G_i, aux_u_s_acc);
        }

        for (IndexType j = 0; j < n_nodes; ++j) {
            const double N_j = rThisKinematicVariables.N[j];
            for (IndexType k = 0;  k < strain_size; ++k) {
                for (IndexType l = 0; l < dim; ++l) {
                    B_j(k,l) = rThisKinematicVariables.B(k, j * dim + l);
                }
            }
            noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_j)));

            // Add the extra terms in the RHS momentum equation
            Vector disp_j(dim);
            for (unsigned int d = 0; d < dim; ++d) {
                disp_j[d] = rThisKinematicVariables.Displacements[j * dim + d];
            }
            const double aux_vol_j = (N_j * rThisKinematicVariables.VolumetricNodalStrains[j] - inner_prod(psi_j, disp_j));
            noalias(rhs_mom_i) += (w_gauss / dim) * transBi_C_invT_m * aux_vol_j;
            noalias(rhs_mom_i) -= w_gauss * (1 - tau_2) * bulk_modulus * psi_i * aux_vol_j;

            rhs_mass_i += w_gauss * (1 - tau_2) * bulk_modulus * N_i * aux_vol_j;
        }

        // Assemble the RHS contributions
        rRightHandSideVector[i * block_size + dim] += rhs_mass_i;
        for (IndexType d = 0; d < dim; ++d) {
            rRightHandSideVector[i * block_size + d] += rhs_mom_i[d];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLeftHandSideGaussPointContribution(
    MatrixType &rLeftHandSideMatrix,
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables,
    const GaussPointAuxiliaryVariables &rThisGaussPointAuxiliaryVariables) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    Matrix lhs_uu_ij(dim, dim);
    Vector lhs_ue_ij(dim);
    Vector lhs_eu_ij(dim);
    double lhs_ee_ij;
    Vector G_i(dim);
    Vector G_j(dim);
    Vector psi_i(dim);
    Vector psi_j(dim);
    Vector transBi_C_invT_m(dim);
    Matrix B_i(strain_size, dim);
    Matrix B_j(strain_size, dim);

    const double dt = rThisGaussPointAuxiliaryVariables.Dt;
    const double density = rThisGaussPointAuxiliaryVariables.Rho;
    const double tau_1 = rThisGaussPointAuxiliaryVariables.Tau1;
    const double tau_2 = rThisGaussPointAuxiliaryVariables.Tau2;
    const double w_gauss = rThisGaussPointAuxiliaryVariables.Weight;
    const double bulk_modulus = rThisGaussPointAuxiliaryVariables.BulkModulus;
    const auto &voigt_identity = rThisGaussPointAuxiliaryVariables.m;
    const Vector C_invT_m_voigt = prod(rThisConstitutiveVariables.D, rThisGaussPointAuxiliaryVariables.invT_m);

    for (IndexType i = 0; i < n_nodes; ++i) {
        const double N_i = rThisKinematicVariables.N[i];
        noalias(G_i) = row(rThisKinematicVariables.DN_DX, i);
        for (IndexType k = 0; k < strain_size; ++k) {
            for (IndexType l = 0; l < dim; ++l) {
                B_i(k, l) = rThisKinematicVariables.B(k, i * dim + l);
            }
        }
        noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));

        for (IndexType j = 0; j < n_nodes; ++j) {
            const double N_j = rThisKinematicVariables.N[j];
            noalias(G_j) = row(rThisKinematicVariables.DN_DX, j);
            for (IndexType k = 0; k < strain_size; ++k) {
                for (IndexType l = 0; l < dim; ++l) {
                    B_j(k, l) = rThisKinematicVariables.B(k, j * dim + l);
                }
            }
            noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_j)));

            // Add momentum internal force LHS contribution
            noalias(lhs_uu_ij) = w_gauss * prod(trans(B_i), Matrix(prod(rThisConstitutiveVariables.D, B_j)));
            noalias(lhs_uu_ij) -= w_gauss * (1 - tau_2) * bulk_modulus * outer_prod(psi_i, trans(psi_j));
            // Add momentum volumetric strain LHS contribution
            noalias(lhs_ue_ij) = w_gauss * (1 - tau_2) * bulk_modulus * psi_i * N_j;
            // Add mass conservation displacement divergence LHS contribution
            noalias(lhs_eu_ij) = w_gauss * (1 - tau_2) * bulk_modulus * N_i * trans(psi_j);
            // Add mass conservation volumetric strain LHS contribution
            lhs_ee_ij = -w_gauss * (1 - tau_2) * bulk_modulus * N_i * N_j;
            lhs_ee_ij -= w_gauss * std::pow(bulk_modulus, 2) * tau_1 * inner_prod(G_i, G_j);
            // Add the dynamic subscale terms to the LHS
            if (mIsDynamic) {
                noalias(lhs_ue_ij) += w_gauss * bulk_modulus * tau_1 * (density / std::pow(dt, 2)) * N_i * G_j;
            }

            // Assemble the LHS contributions
            rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim) += lhs_ee_ij;
            for (IndexType d = 0; d < dim; ++d) {
                rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += lhs_ue_ij[d];
                rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) += lhs_eu_ij[d];
                for (IndexType d2 = 0; d2 < dim; ++d2) {
                    rLeftHandSideMatrix(i * block_size + d, j * block_size + d2) += lhs_uu_ij(d, d2);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::UpdateGaussPointDisplacementSubscaleHistory(
    const KinematicVariables& rThisKinematicVariables,
    const ConstitutiveVariables& rThisConstitutiveVariables,
    const ProcessInfo& rProcessInfo,
    const IndexType PointIndex)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector m_T(strain_size);
    Vector voigt_identity = ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }
    m_T = prod(voigt_identity, mAnisotropyTensor);

    // Get current Gauss point data
    const double bulk_modulus = CalculateBulkModulus(rThisConstitutiveVariables.D);
    const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), PointIndex);
    const double tau_1 = CalculateTau1(m_T, rThisKinematicVariables, rThisConstitutiveVariables, rProcessInfo);
    const Vector grad_eps = prod(trans(rThisKinematicVariables.DN_DX), rThisKinematicVariables.VolumetricNodalStrains);

    // Calculate the n+1 subscale with current database
    const double dt = rProcessInfo[DELTA_TIME];
    const double density = GetProperties()[DENSITY];
    Vector aux_u_s = (density / std::pow(dt,2)) * (2.0*mDisplacementSubscale1[PointIndex] - mDisplacementSubscale2[PointIndex]);
    for (IndexType d = 0; d < dim; ++d) {
        aux_u_s[d] += body_force[d];
        aux_u_s[d] += bulk_modulus * grad_eps[d];
    }

    // Multiply by the stabilization parameter
    aux_u_s *= tau_1;

    // Update subscale vector values
    mDisplacementSubscale2[PointIndex] = mDisplacementSubscale1[PointIndex]; // n-1 subscale value update
    mDisplacementSubscale1[PointIndex] = aux_u_s; // n subscale value update
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType dim = r_geom.WorkingSpaceDimension();
    const SizeType n_nodes = r_geom.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Clear matrix
    if (rMassMatrix.size1() != matrix_size || rMassMatrix.size2() != matrix_size) {
        rMassMatrix.resize(matrix_size, matrix_size, false);
    }
    rMassMatrix.clear();

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geom[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geom[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geom, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set the auxiliary Gauss point variable container
    GaussPointAuxiliaryVariables gauss_point_auxiliary_variables(this, dim, strain_size);

    // Calculate the consistent mass matrix
    const double density = r_prop[DENSITY];
    const double thickness = (dim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;
    const bool compute_lumped_mass_matrix = r_prop.Has(COMPUTE_LUMPED_MASS_MATRIX) ? r_prop[COMPUTE_LUMPED_MASS_MATRIX] : false;

    if (compute_lumped_mass_matrix) {
        KRATOS_ERROR << "Mass matrix is not symmetric and cannot be lumped." << std::endl;
    } else {
        Matrix J0;
        Vector G_i(dim);
        const SizeType n_gauss = r_geom.IntegrationPointsNumber(GetIntegrationMethod());
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Calculate the constitutive response
            // Note that this is required in order to get the effective bulk modulus
            CalculateConstitutiveVariables(
                kinematic_variables,
                constitutive_variables,
                cons_law_values,
                i_gauss,
                r_geom.IntegrationPoints(this->GetIntegrationMethod()),
                ConstitutiveLaw::StressMeasure_Cauchy);

            // Calculate the auxiliary Gauss point variables
            CalculateGaussPointAuxiliaryVariables(
                gauss_point_auxiliary_variables,
                kinematic_variables,
                constitutive_variables,
                rCurrentProcessInfo,
                i_gauss);

            // Get required data
            const double tau_1 = gauss_point_auxiliary_variables.Tau1;
            const double w_gauss = thickness * gauss_point_auxiliary_variables.Weight;
            const double bulk_modulus = gauss_point_auxiliary_variables.BulkModulus;

            // Check if the problem is dynamic to add the dynamic subscale component
            // Note that this check avoids this contribution to be added to the eigenvalue problem mass matrix
            double mass_factor = 1.0;
            if (mIsDynamic) {
                const double dt = gauss_point_auxiliary_variables.Dt;
                mass_factor -= tau_1 / std::pow(dt,2);
            }

            // Assemble nodal contributions
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                const double N_i = kinematic_variables.N[i_node];
                noalias(G_i) = row(kinematic_variables.DN_DX, i_node);
                for (IndexType j_node = 0; j_node < n_nodes; ++j_node) {
                    const double N_j = kinematic_variables.N[j_node];
                    const double aux_1 = density * mass_factor * w_gauss * N_i * N_j;
                    const double aux_2 = bulk_modulus * tau_1 * w_gauss * N_j;
                    for (IndexType d = 0; d < dim; ++d) {
                        rMassMatrix(i_node*block_size + d, j_node*block_size + d) += aux_1;
                        rMassMatrix(i_node*block_size + dim, j_node*block_size + d) += aux_2 * G_i[d]; // Consistent (but not symmetric) mass matrix contribution
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetSecondDerivativesVector(
    Vector &rValues,
    int Step) const
{
    const auto &r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();
    const SizeType n_nodes = r_geom.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;

    if (rValues.size() != matrix_size) {
        rValues.resize(matrix_size, false);
    }

    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const SizeType index = i_node * block_size;
        const auto& r_acc = r_geom[i_node].FastGetSolutionStepValue(ACCELERATION, Step);
        for(IndexType d = 0; d < dim; ++d) {
            rValues[index + d] = r_acc[d];
        }
        rValues[index + dim] = 0.0;
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
    ) const
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
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

Vector SmallDisplacementMixedVolumetricStrainElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const auto aux_body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
    Vector body_force(dim);
    for (IndexType d = 0; d < dim; ++d) {
        body_force(d) = aux_body_force(d);
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

void SmallDisplacementMixedVolumetricStrainElement::CalculateGaussPointAuxiliaryVariables(
    GaussPointAuxiliaryVariables& rThisGaussPointAuxiliaryVariables,
    const KinematicVariables& rThisKinematicVariables,
    const ConstitutiveVariables& rThisConstitutiveVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const IndexType PointNumber) const
{
    // Get auxiliary references
    double& dt = rThisGaussPointAuxiliaryVariables.Dt;
    double& rho = rThisGaussPointAuxiliaryVariables.Rho;
    double& tau_1 = rThisGaussPointAuxiliaryVariables.Tau1;
    double& tau_2 = rThisGaussPointAuxiliaryVariables.Tau2;
    double& w_gauss = rThisGaussPointAuxiliaryVariables.Weight;
    double& bulk_modulus = rThisGaussPointAuxiliaryVariables.BulkModulus;
    auto& body_force = rThisGaussPointAuxiliaryVariables.BodyForce;
    auto& grad_eps = rThisGaussPointAuxiliaryVariables.VolumetricStrainGradient;
    auto& us_acc = rThisGaussPointAuxiliaryVariables.DisplacementSubscaleAcceleration;
    auto& m_T = rThisGaussPointAuxiliaryVariables.m_T;

    // Calculate Gauss point weight
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    const double thickness = (r_geometry.WorkingSpaceDimension() == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
    w_gauss = thickness * rThisKinematicVariables.detJ0 * r_integration_points[PointNumber].Weight();

    // Calculate tau_1 stabilization constant
    tau_1 = CalculateTau1(m_T, rThisKinematicVariables, rThisConstitutiveVariables, rCurrentProcessInfo);

    // Calculate tau_2 stabilization constant
    tau_2 = CalculateTau2(rThisConstitutiveVariables);

    // Calculate the bulk modulus
    bulk_modulus = CalculateBulkModulus(rThisConstitutiveVariables.D);

    // Calculate and add the LHS contributions
    body_force = GetBodyForce(r_integration_points, PointNumber);
    grad_eps = prod(trans(rThisKinematicVariables.DN_DX), rThisKinematicVariables.VolumetricNodalStrains);

    // Get dynamic data if required
    if (mIsDynamic) {
        rho = GetProperties()[DENSITY];
        dt = rCurrentProcessInfo[DELTA_TIME];
        us_acc = (rho / std::pow(dt, 2)) * (2.0 * mDisplacementSubscale1[PointNumber] - mDisplacementSubscale2[PointNumber]);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX
    ) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    const auto& r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    const SizeType dim = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Set the identity in Voigt notation
    Vector voigt_identity = ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }

    // Calculate the total strain
    Vector eq_strain = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);

    // Calculate the volumetric contribution
    const Vector T_B_U = prod(mAnisotropyTensor, eq_strain);
    double gauss_vol_strain = inner_prod(voigt_identity, T_B_U);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        gauss_vol_strain -= rThisKinematicVariables.N[i_node] * rThisKinematicVariables.VolumetricNodalStrains[i_node];
    }
    gauss_vol_strain /= dim;

    // Substract the volumetric contribution to the total strain
    eq_strain -= gauss_vol_strain * prod(mInverseAnisotropyTensor, voigt_identity);

    // Save the obtained value in the kinematics container
    rThisKinematicVariables.EquivalentStrain = eq_strain;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateAnisotropyTensor(const ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Set a null equivalent strain field to calculate the initial C tensor
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    kinematic_variables.EquivalentStrain = ZeroVector(strain_size);

    // Calculate a fake material response and calculate the initial C tensor
    // Note that we take the first Gauss point as the T tensor must be identical for all the Gauss pts
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, 0, r_integration_points);
    mConstitutiveLawVector[0]->InitializeMaterialResponseCauchy(cons_law_values);
    mConstitutiveLawVector[0]->CalculateMaterialResponse(cons_law_values, ConstitutiveLaw::StressMeasure_Cauchy);

    // Calculate the closest isotropic tensor
    // The k_iso and mu_iso coefficients are obtained from the Frobenius norm minimization of the closest isotropic tensor to C
    // Such minimization is constrained such that the compressibility of the new C_iso is equivalent to the input anisotropic C one
    Vector u = ZeroVector(strain_size);
    for (IndexType i = 0; i < dim; ++i) {
        u(i) = 1.0 / std::sqrt(dim);
    }
    const Matrix J = outer_prod(u, u);
    Matrix I = ZeroMatrix(strain_size);
    for (IndexType i = 0; i < strain_size; ++i){
        I(i, i) = i < dim ? 1.0 : 0.5;
    }
    const Matrix K = I - J;
    const Matrix &rC = constitutive_variables.D;
    const double k_iso = CalculateBulkModulus(rC);
    const double mu_iso = CalculateShearModulus(rC);
    const Matrix C_iso = 3.0 * (dim * k_iso / 3.0) * J + 2.0 * mu_iso * K;

    // Calculate the square root of the C and closest isotropic C tensors
    Matrix a;
    Matrix b;
    const double tolerance = 1.0e-12;
    const unsigned int max_iterations = 100;
    const bool is_converged_a = MathUtils<double>::MatrixSquareRoot(rC, a, tolerance, max_iterations);
    KRATOS_WARNING_IF("SmallDisplacementMixedVolumetricStrainElement", !is_converged_a) << "Element " << Id() << " anisotropic tensor square root did not converge.";
    const bool is_converged_b = MathUtils<double>::MatrixSquareRoot(C_iso, b, tolerance, max_iterations);
    KRATOS_WARNING_IF("SmallDisplacementMixedVolumetricStrainElement", !is_converged_b) << "Element " << Id() << " isotropic tensor square root did not converge.";

    // Calculate the anisotropy tensor T as inv(b)*a
    Matrix inv_b;
    double det_b;
    MathUtils<double>::InvertMatrix(b, inv_b, det_b);
    mAnisotropyTensor = prod(inv_b, a);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateInverseAnisotropyTensor()
{
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    double aux_det;
    mInverseAnisotropyTensor = ZeroMatrix(strain_size, strain_size);
    MathUtils<double>::InvertMatrix(mAnisotropyTensor, mInverseAnisotropyTensor, aux_det);
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateBulkModulus(const Matrix &rC) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    double bulk_modulus = 0.0;
    for (SizeType i = 0; i < dim; ++i) {
        for (SizeType j = 0; j < dim; ++j) {
            bulk_modulus += rC(i,j);
        }
    }
    return bulk_modulus / std::pow(dim, 2);
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateShearModulus(const Matrix &rC) const
{
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    return (strain_size == 3) ?
        0.2 * (rC(0,0) - 2.0*rC(0,1) + rC(1,1) + rC(2,2)) :
        (4.0 / 33.0)*(rC(0,0) - rC(0,1) - rC(0,2) + rC(1,1) - rC(1,2) + rC(2,2) + (3.0/4.0)*(rC(3,3) + rC(4,4) + rC(5,5)));
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateTau1(
    const Vector &rVoigtIdAnysotropyMatrixProd,
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables,
    const ProcessInfo &rCurrentProcessInfo) const
{
    const auto& r_prop = GetProperties();
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = r_prop.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    Matrix B_i(strain_size, dim);
    Matrix aux = ZeroMatrix(dim,dim);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        for (IndexType k = 0;  k < strain_size; ++k) {
            for (IndexType l = 0; l < dim; ++l) {
                B_i(k,l) = rThisKinematicVariables.B(k, i_node * dim + l);
            }
        }
        noalias(aux) += prod(trans(B_i), Matrix(prod(rThisConstitutiveVariables.D, B_i)));
    }
    double det;
    const double c_1 = 1.0;
    Matrix tau_1_mat(dim, dim);
    MathUtils<double>::InvertMatrix(aux, tau_1_mat, det);
    const Vector tau_1_vect = MathUtils<double>::SymmetricTensorToVector(tau_1_mat);
    const double m_T_tau_1 = c_1*inner_prod(rVoigtIdAnysotropyMatrixProd, tau_1_vect);

    // Uncomment this to use the stabilisation parameters defined in the paper
    // const double c_1 = 1.0;
    // const double G = CalculateShearModulus(rThisConstitutiveVariables.D);
    // const double h = n_nodes == 3 ? ElementSizeCalculator<2, 3>::AverageElementSize(this->GetGeometry()) : ElementSizeCalculator<2, 4>::AverageElementSize(this->GetGeometry());
    // const double m_T_tau_1 = c_1 * std::pow(h,2) / G;

    double tau_1;
    if (mIsDynamic) {
        const double density = r_prop[DENSITY];
        const double dt = rCurrentProcessInfo[DELTA_TIME];
        tau_1 = 1.0 / (density / std::pow(dt, 2) + 1.0 / m_T_tau_1);
    } else {
        tau_1 = m_T_tau_1;
    }

    return tau_1;
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateTau2(const ConstitutiveVariables& rThisConstitutiveVariables) const
{
    const double c_2 = 4.0;
    const double max_tau_2 = 1.0e-2;
    const double bulk_modulus = CalculateBulkModulus(rThisConstitutiveVariables.D);
    const double shear_modulus = CalculateShearModulus(rThisConstitutiveVariables.D);
    return std::min(max_tau_2, c_2 * shear_modulus / bulk_modulus);
    // return c_2 * shear_modulus / (shear_modulus + bulk_modulus); // Uncomment this to use the stabilisation parameters defined in the paper
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    ) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedVolumetricStrainElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = SmallDisplacementMixedVolumetricStrainElement::BaseType::Check(rCurrentProcessInfo);

    // Base check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto& r_geometry = this->GetGeometry();
    for ( IndexType i = 0; i < r_geometry.size(); i++ ) {
        const NodeType& r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUMETRIC_STRAIN,r_node)
        KRATOS_CHECK_DOF_IN_NODE(VOLUMETRIC_STRAIN, r_node)
    }

    // Check if the density is provided for the calculation of the intertial terms and dynamic subscales
    if (mIsDynamic) {
        const auto& r_prop = GetProperties();
        KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY needs to be provided for the calculation of dynamic subscales!" << std::endl;
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Resize output container
    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else if (rVariable == VOLUMETRIC_STRAIN_SUBSCALE) {
        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
            }
            kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Set the constitutive law values
        ConstitutiveVariables constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        // Calculate the Gauss points volumetric strain subscale values
        Vector voigt_identity = ZeroVector(strain_size);
        for (IndexType d = 0; d < dim; ++d) {
            voigt_identity[d] = 1.0;
        }

        Vector psi_i(dim);
        Vector disp_i(dim);
        Matrix B_i(strain_size, dim);
        for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
            // Initialize current value
            rOutput[i_gauss] = 0.0;

            // Recompute the kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Set the constitutive variables
            CalculateConstitutiveVariables(
                kinematic_variables,
                constitutive_variables,
                cons_law_values,
                i_gauss,r_integration_points,
                ConstitutiveLaw::StressMeasure_Cauchy);

            // Calculate tau_2 stabilization constant
            const double tau_2 = CalculateTau2(constitutive_variables);

            // Calculate current gauss point subscale value
            for (IndexType i = 0; i < n_nodes; ++i) {
                for (IndexType d = 0; d < dim; ++d) {
                    for (IndexType k = 0;  k < strain_size; ++k) {
                        B_i(k,d) = kinematic_variables.B(k, i * dim + d);
                    }
                    disp_i[d] = kinematic_variables.Displacements(i*n_nodes + d);
                }
                noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(mAnisotropyTensor, B_i)));
                const double eps_vol_i = kinematic_variables.N[i] * kinematic_variables.VolumetricNodalStrains[i];
                rOutput[i_gauss] += tau_2 * (inner_prod(psi_i, disp_i) - eps_vol_i);
            }
        }
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
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
                // Compute material response
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            } else {
                // Compute material response
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
            }

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != strain_size) {
                rOutput[i_gauss].resize(strain_size, false);
            }
            rOutput[i_gauss] = constitutive_variables.StressVector;
        }
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Get dynamic data if required
    const double density = mIsDynamic ? GetProperties()[DENSITY] : 0.0;
    const double dt = mIsDynamic ? rCurrentProcessInfo[DELTA_TIME] : 0.0;

    // Resize output container
    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else if (rVariable == DISPLACEMENT_SUBSCALE) {
        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
            }
            kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Set the constitutive law values
        ConstitutiveVariables constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        // Calculate the Gauss points displacement subscale values
        Vector m_T;
        Vector voigt_identity;
        if (mIsDynamic) {
            voigt_identity =  ZeroVector(strain_size);
            for (IndexType d = 0; d < dim; ++d) {
                voigt_identity[d] = 1.0;
            }
            m_T = prod(voigt_identity, mAnisotropyTensor);
        }

        Vector aux_u_s(dim);
        Vector grad_eps(dim);
        array_1d<double, 3> body_force;
        array_1d<double, 3> u_proj_gauss;
        for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
            // Recompute the kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Set the constitutive variables
            CalculateConstitutiveVariables(
                kinematic_variables,
                constitutive_variables,
                cons_law_values,
                i_gauss,r_integration_points,
                ConstitutiveLaw::StressMeasure_Cauchy);

            // Calculate tau_1 stabilization constant
            const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rCurrentProcessInfo);

            // Get current Gauss point data
            const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
            noalias(body_force) = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
            noalias(grad_eps) = prod(trans(kinematic_variables.DN_DX), kinematic_variables.VolumetricNodalStrains);

            // If dynamic subscales are active, add the dynamic subscale component
            aux_u_s = bulk_modulus * grad_eps;
            for (IndexType d = 0; d < dim; ++d) {
                aux_u_s[d] += body_force[d];
            }
            if (mIsDynamic) {
                aux_u_s += (density / std::pow(dt,2)) * (2.0*mDisplacementSubscale1[i_gauss] - mDisplacementSubscale2[i_gauss]);
            }

            // Multiply by the stabilization parameter
            aux_u_s *= tau_1;

            // Save current value in the output container
            for (IndexType d = 0; d < dim; ++d) {
                rOutput[i_gauss][d] = aux_u_s[d];
            }
        }
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedVolumetricStrainElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","dynamic"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["CAUCHY_STRESS_VECTOR"],
            "nodal_historical"       : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
            "nodal_non_historical"   : [],
            "entity"                 : ["DISPLACEMENT_SUBSCALE", "VOLUMETRIC_STRAIN_SUBSCALE"]
        },
        "required_variables"         : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
        "required_dofs"              : [],
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
            "This element implements a mixed displacement - volumetric strain formulation with Variational MultiScales (VMS) stabilization. This formulation is capable to deal with materials in the incompressible limit as well as with anisotropy. Note that the resulting mass matrix is not symmetric."
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

void SmallDisplacementMixedVolumetricStrainElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.save("AnisotropyTensor", mAnisotropyTensor);
    rSerializer.save("InverseAnisotropyTensor", mInverseAnisotropyTensor);
    // rSerializer.save("DisplacementSubscale1", mDisplacementSubscale1);
    // rSerializer.save("DisplacementSubscale1", mDisplacementSubscale2);
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
    rSerializer.load("AnisotropyTensor", mAnisotropyTensor);
    rSerializer.load("InverseAnisotropyTensor", mInverseAnisotropyTensor);
    // rSerializer.load("DisplacementSubscale1", mDisplacementSubscale1);
    // rSerializer.load("DisplacementSubscale1", mDisplacementSubscale2);
}

} // Namespace Kratos
