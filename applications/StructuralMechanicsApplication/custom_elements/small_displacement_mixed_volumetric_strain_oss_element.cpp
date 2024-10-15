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

    auto p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssElement>(
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

void SmallDisplacementMixedVolumetricStrainOssElement::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    // Deactivate the dynamic subscales as these are not implemented in the OSS elements yet
    mIsDynamic = 0;
    mDisplacementSubscale1.resize(0);
    mDisplacementSubscale2.resize(0);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateRightHandSide(
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

    // Allocate a temporary matrix for the OSS projection operator
    MatrixType oss_proj_op = ZeroMatrix(matrix_size, matrix_size);

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

        // Assemble the current Gauss point OSS projection operator
        // Note that this calculates the OSS stabilization operator to be applied to the RHS
        CalculateOssStabilizationOperatorGaussPointContribution(
            oss_proj_op,
            kinematic_variables,
            gauss_point_auxiliary_variables);
    }

    // Substract the OSS projections to the current residual
    VectorType proj_vect(matrix_size);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_us_proj = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
        for (IndexType d = 0; d < dim; ++d) {
            proj_vect(i_node*block_size + d) = r_us_proj[d];
        }
        proj_vect(i_node*block_size + dim) = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
    }
    rRightHandSideVector += prod(oss_proj_op, proj_vect);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateLocalSystem(
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

    // Allocate a temporary matrix for the OSS projection operator
    MatrixType oss_proj_op = ZeroMatrix(matrix_size, matrix_size);

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

        // Assemble the current Gauss point OSS projection operator
        // Note that this calculates the OSS stabilization operator to be applied to the RHS
        CalculateOssStabilizationOperatorGaussPointContribution(
            oss_proj_op,
            kinematic_variables,
            gauss_point_auxiliary_variables);
    }

    // Substract the OSS projections to the current residual
    VectorType proj_vect(matrix_size);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_us_proj = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
        for (IndexType d = 0; d < dim; ++d) {
            proj_vect(i_node*block_size + d) = r_us_proj[d];
        }
        proj_vect(i_node*block_size + dim) = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
    }
    rRightHandSideVector += prod(oss_proj_op, proj_vect);
}


/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_prop = GetProperties();
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;

    // Check LHS size and initialize
    if (rMassMatrix.size1() != matrix_size || rMassMatrix.size2() != matrix_size) {
        rMassMatrix.resize(matrix_size, matrix_size, false);
    }
    rMassMatrix.clear();

    // Calculate the consistent mass matrix
    // Note that we cannot use the base element implementation as it includes some extra terms that are 0 in the OSS case
    const double density = r_prop[DENSITY];
    const double thickness = (dim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;
    const double aux_rho_thickness = density * thickness;
    const bool compute_lumped_mass_matrix = r_prop.Has(COMPUTE_LUMPED_MASS_MATRIX) ? r_prop[COMPUTE_LUMPED_MASS_MATRIX] : false;

    if (compute_lumped_mass_matrix) {
        KRATOS_ERROR << "Lumped mass matrix is not implemented." << std::endl;
    } else {
        Matrix J0;
        const auto& r_integration_method = this->GetIntegrationMethod();
        const auto& r_N_container = r_geometry.ShapeFunctionsValues(r_integration_method);
        const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
        const SizeType n_gauss = r_integration_points.size();
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate current integration point weight
            // Note that this includes the density and the thickness
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[i_gauss], J0);
            const double detJ0 = MathUtils<double>::Det(J0);
            const double w_gauss = aux_rho_thickness * detJ0 * r_integration_points[i_gauss].Weight();

            // Assemble nodal contributions
            const auto &r_N = row(r_N_container, i_gauss);
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                for (IndexType j_node = 0; j_node < n_nodes; ++j_node) {
                    const double aux = w_gauss * r_N[i_node] * r_N[j_node];
                    for (IndexType d = 0; d < dim; ++d) {
                        rMassMatrix(i_node*block_size + d, j_node*block_size + d) += aux;
                    }
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VOLUMETRIC_STRAIN_PROJECTION) {
        // Get geometry data
        auto& r_geometry = GetGeometry();
        const auto& r_prop = GetProperties();
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
        const SizeType n_gauss = r_integration_points.size();
        const SizeType strain_size = r_prop.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

        // Calculate the mass matrix (lumped or consistent)
        // const double density = r_prop[DENSITY]; //TODO: Leave this for the weighted projection
        const double thickness = (dim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
            }
            kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Calculate the kinematic constraint residual projection
        Vector voigt_identity = ZeroVector(strain_size);
        for (IndexType d = 0; d < dim; ++d) {
            voigt_identity[d] = 1.0;
        }
        Vector psi_j(dim);
        Vector aux_proj = ZeroVector(n_nodes);
        Matrix B_j(strain_size, dim);
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Calculate current integration point weight
            // Note that this already includes the thickness
            const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

            // Calculate current Gauss point volumetric strain subscale projection contribution
            double aux_res = 0.0;
            for (IndexType j_node = 0; j_node < n_nodes; ++j_node) {
                for (IndexType k = 0;  k < strain_size; ++k) {
                    for (IndexType l = 0; l < dim; ++l) {
                        B_j(k,l) = kinematic_variables.B(k, j_node * dim + l);
                    }
                }
                const double N_j = kinematic_variables.N[j_node];
                noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(GetAnisotropyTensor(), B_j)));
                for (IndexType d = 0; d < dim; ++d) {
                    aux_res += psi_j[d] * kinematic_variables.Displacements[j_node*dim + d];
                }
                aux_res -= N_j * kinematic_variables.VolumetricNodalStrains[j_node];
            }

            // Nodal assembly of the projection contributions
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                const double N_i = kinematic_variables.N[i_node];
                AtomicAdd(r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION), w_gauss*N_i*aux_res);
            }
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == DISPLACEMENT_PROJECTION) {
        // Get geometry data
        auto& r_geometry = GetGeometry();
        const auto& r_prop = GetProperties();
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
        const SizeType n_gauss = r_integration_points.size();
        const SizeType strain_size = r_prop.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

        // Calculate the mass matrix (lumped or consistent)
        // const double density = r_prop[DENSITY]; //TODO: Leave this for the weighted projection

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
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_prop, rCurrentProcessInfo);
        auto &r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        // Set the auxiliary Gauss point variable container
        GaussPointAuxiliaryVariables gauss_point_auxiliary_variables(this, dim, strain_size);

        // Calculate the momentum equation projection
        array_1d<double,3> aux_proj;
        array_1d<double,3> w_Ni_aux_proj;
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Calculate the constitutive response
            CalculateConstitutiveVariables(
                kinematic_variables,
                constitutive_variables,
                cons_law_values,
                i_gauss,
                r_integration_points,
                ConstitutiveLaw::StressMeasure_Cauchy);

            // Calculate the auxiliary Gauss point variables
            CalculateGaussPointAuxiliaryVariables(
                gauss_point_auxiliary_variables,
                kinematic_variables,
                constitutive_variables,
                rCurrentProcessInfo,
                i_gauss);

            // Calculate current Gauss point displacement subscale projection contribution
            aux_proj = ZeroVector(3);
            const auto& r_body_force = gauss_point_auxiliary_variables.BodyForce;
            const auto& r_grad_eps = gauss_point_auxiliary_variables.VolumetricStrainGradient;
            for (IndexType d = 0; d < dim; ++d) {
                aux_proj[d] = r_body_force[d] + gauss_point_auxiliary_variables.BulkModulus * r_grad_eps[d];
            }

            // Nodal assembly of the projection contributions
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                w_Ni_aux_proj = gauss_point_auxiliary_variables.Weight * kinematic_variables.N[i_node] * aux_proj;
                AtomicAdd(r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT_PROJECTION), w_Ni_aux_proj);
            }
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::UpdateGaussPointDisplacementSubscaleHistory(
    const KinematicVariables &rThisKinematicVariables,
    const ConstitutiveVariables &rThisConstitutiveVariables,
    const ProcessInfo &rProcessInfo,
    const IndexType PointIndex)
{
    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector m_T(strain_size);
    Vector voigt_identity = ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d){
        voigt_identity[d] = 1.0;
    }
    m_T = prod(voigt_identity, GetAnisotropyTensor());

    // Get current Gauss point data
    const double bulk_modulus = CalculateBulkModulus(rThisConstitutiveVariables.D);
    const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), PointIndex);
    const double tau_1 = CalculateTau1(m_T, rThisKinematicVariables, rThisConstitutiveVariables, rProcessInfo);
    const Vector grad_eps = prod(trans(rThisKinematicVariables.DN_DX), rThisKinematicVariables.VolumetricNodalStrains);

    // Calculate the n+1 subscale with current database
    const double dt = rProcessInfo[DELTA_TIME];
    const double density = GetProperties()[DENSITY];
    Vector aux_u_s = (density / std::pow(dt, 2)) * (2.0 * mDisplacementSubscale1[PointIndex] - mDisplacementSubscale2[PointIndex]);
    for (IndexType d = 0; d < dim; ++d) {
        aux_u_s[d] += body_force[d];
        aux_u_s[d] += bulk_modulus * grad_eps[d];
    }

    // Calculate the displacement residual projection at current Gauss point
    array_1d<double, 3> us_proj_gauss = ZeroVector(3);
    for (IndexType i = 0; i < n_nodes; ++i) {
        const double N_i = rThisKinematicVariables.N[i];
        const auto &r_u_proj_i = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
        noalias(us_proj_gauss) += N_i * r_u_proj_i;
    }

    // Substract the nodal projection to get the orthogonal subscale value
    for (IndexType d = 0; d < dim; ++d) {
        aux_u_s[d] -= us_proj_gauss[d];
    }

    // Multiply by the stabilization parameter
    aux_u_s *= tau_1;

    // Update subscale vector values
    mDisplacementSubscale2[PointIndex] = mDisplacementSubscale1[PointIndex]; // n-1 subscale value update
    mDisplacementSubscale1[PointIndex] = aux_u_s; // n subscale value update
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateOssStabilizationOperatorGaussPointContribution(
    MatrixType& rOrthogonalSubScalesOperator,
    const KinematicVariables& rThisKinematicVariables,
    const GaussPointAuxiliaryVariables& rThisGaussPointAuxiliaryVariables) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    const double w_kappa = rThisGaussPointAuxiliaryVariables.Weight * rThisGaussPointAuxiliaryVariables.BulkModulus;
    const double aux_w_kappa_tau_1 = w_kappa * rThisGaussPointAuxiliaryVariables.Tau1;
    const double aux_w_kappa_tau_2 = w_kappa * rThisGaussPointAuxiliaryVariables.Tau2;

    Vector G_i(dim);
    Vector psi_i(dim);
    Matrix B_i(strain_size, dim);
    for (IndexType i = 0; i < n_nodes; ++i) {
        const double N_i = rThisKinematicVariables.N[i];
        noalias(G_i) = row(rThisKinematicVariables.DN_DX, i);
        for (IndexType k = 0;  k < strain_size; ++k) {
            for (IndexType l = 0; l < dim; ++l) {
                B_i(k, l) = rThisKinematicVariables.B(k, i * dim + l);
            }
        }
        noalias(psi_i) = prod(trans(rThisGaussPointAuxiliaryVariables.m), Matrix(prod(GetAnisotropyTensor(), B_i)));

        for (IndexType j = 0; j < n_nodes; ++j) {
            const double N_j = rThisKinematicVariables.N[j];
            rOrthogonalSubScalesOperator(i*block_size + dim, j*block_size + dim) -= aux_w_kappa_tau_2 * N_i * N_j;
            for (IndexType d = 0; d < dim; ++d) {
                rOrthogonalSubScalesOperator(i*block_size + dim, j*block_size + d) -= aux_w_kappa_tau_1 * G_i[d] * N_j;
                rOrthogonalSubScalesOperator(i*block_size + d, j*block_size + dim) += aux_w_kappa_tau_2 * psi_i[d] * N_j;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int SmallDisplacementMixedVolumetricStrainOssElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

    if (rVariable == VOLUMETRIC_STRAIN_SUBSCALE) {
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
                noalias(psi_i) = prod(trans(voigt_identity), Matrix(prod(GetAnisotropyTensor(), B_i)));
                const double eps_vol_i = kinematic_variables.N[i] * kinematic_variables.VolumetricNodalStrains[i];
                const double eps_vol_proj_i = kinematic_variables.N[i] * r_geometry[i].FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
                rOutput[i_gauss] += tau_2 * (inner_prod(psi_i, disp_i) - eps_vol_i);
                rOutput[i_gauss] -= tau_2 * eps_vol_proj_i;
            }
        }
    } else {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
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

    else if (rVariable == DISPLACEMENT_SUBSCALE) {
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
        if (IsDynamic()) {
            voigt_identity =  ZeroVector(strain_size);
            for (IndexType d = 0; d < dim; ++d) {
                voigt_identity[d] = 1.0;
            }
            m_T = prod(voigt_identity, GetAnisotropyTensor());
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
            if (IsDynamic()) {
                const double density = GetProperties()[DENSITY];
                const double dt = rCurrentProcessInfo[DELTA_TIME];
                aux_u_s += (density / std::pow(dt,2)) * (2.0*mDisplacementSubscale1[i_gauss] - mDisplacementSubscale2[i_gauss]);
            }

            // Calculate the displacement residual projection at current Gauss point
            noalias(u_proj_gauss) = ZeroVector(3);
            for (IndexType i = 0; i < n_nodes; ++i) {
                const double N_i = kinematic_variables.N[i];
                const auto& r_u_proj_i = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
                noalias(u_proj_gauss) += N_i * r_u_proj_i;
            }

            // Substract the nodal projection to get the orthogonal subscale value
            for (IndexType d = 0; d < dim; ++d) {
                aux_u_s[d] -= u_proj_gauss[d];
            }

            // Multiply by the stabilization parameter
            aux_u_s *= tau_1;

            // Save current value in the output container
            for (IndexType d = 0; d < dim; ++d) {
                rOutput[i_gauss][d] = aux_u_s[d];
            }
        }
    } else {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
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
            "This element implements a mixed displacement - volumetric strain formulation with Variational MultiScales (VMS) stabilization.This formulation is capable to deal with materials in the incompressible limit as well as with anisotropy. As a difference to the base element, this element uses an OSS stabilization approach with linearised projections (note that this requires the corresponding scheme to calculate the projections at the end of the iteration). Thanks to the OSS approach, the resulting mass matrix is symmetric. However, this element cannot be used for the eigenvalue analysis (use the non-linear OSS instead)."
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
