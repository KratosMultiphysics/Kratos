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
#include "custom_elements/small_displacement_mixed_volumetric_strain_oss_non_linear_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedVolumetricStrainOssNonLinearElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssNonLinearElement>(
        NewId,
        GetGeometry().Create(ThisNodes),
        pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainOssNonLinearElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssNonLinearElement>(
        NewId,
        pGeom,
        pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainOssNonLinearElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedVolumetricStrainOssNonLinearElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainOssNonLinearElement>(
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

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType dof_size = 2.0*n_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType eps_vol_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN);
    const IndexType disp_proj_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_PROJECTION_X);
    const IndexType eps_vol_proj_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN_PROJECTION);

    IndexType aux_index = 0;
    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_PROJECTION_X, disp_proj_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_PROJECTION_Y, disp_proj_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN_PROJECTION, eps_vol_proj_pos).EquationId();
        }
    } else {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_PROJECTION_X, disp_proj_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_PROJECTION_Y, disp_proj_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_PROJECTION_Z, disp_proj_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN_PROJECTION, eps_vol_proj_pos).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType block_size = dim + 1;
    const IndexType block_size_full = 2 * block_size;
    const IndexType matrix_size = block_size * n_nodes;
    const IndexType matrix_size_full = n_nodes * block_size_full;

    if (rElementalDofList.size() != matrix_size_full){
        rElementalDofList.resize(matrix_size_full);
    }

    if (dim == 2) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * block_size] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * block_size + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
            rElementalDofList[matrix_size + i * block_size] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_X);
            rElementalDofList[matrix_size + i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_Y);
            rElementalDofList[matrix_size + i * block_size + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN_PROJECTION);
        }
    } else if (dim == 3) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * block_size] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * block_size + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * block_size + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
            rElementalDofList[matrix_size + i * block_size] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_X);
            rElementalDofList[matrix_size + i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_Y);
            rElementalDofList[matrix_size + i * block_size + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_Z);
            rElementalDofList[matrix_size + i * block_size + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN_PROJECTION);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType block_size_full = 2.0 * block_size;
    const SizeType matrix_size_full = block_size_full * n_nodes;

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size_full) {
        rRightHandSideVector.resize(matrix_size_full, false);
    }

    // Check LHS size and initialize
    if (rLeftHandSideMatrix.size1() != matrix_size_full || rLeftHandSideMatrix.size2() != matrix_size_full) {
        rLeftHandSideMatrix.resize(matrix_size_full, matrix_size_full, false);
    }
    rLeftHandSideMatrix.clear();

    // Calculate the LHS matrix
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // Calculate the RHS vector
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType block_size_full = 2.0 * block_size;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType matrix_size_full = block_size_full * n_nodes;

    // Check LHS size and initialize
    if (rLeftHandSideMatrix.size1() != matrix_size_full || rLeftHandSideMatrix.size2() != matrix_size_full) {
        rLeftHandSideMatrix.resize(matrix_size_full, matrix_size_full, false);
    }
    rLeftHandSideMatrix.clear();

    // Call the base element to get the stiffness matrix
    MatrixType aux_stiffnes(matrix_size, matrix_size);
    BaseType::CalculateLeftHandSide(
        aux_stiffnes,
        rCurrentProcessInfo);

    // Call the base element OSS operator
    // Note that this calculates the LHS of the OSS stabilization operator
    MatrixType aux_oss_stab_operator(matrix_size, matrix_size);
    CalculateOrthogonalSubScalesStabilizationOperator(
        aux_oss_stab_operator,
        rCurrentProcessInfo);

    // Call the element OSS operator
    // Note that this calculates the OSS operator to be applied to the RHS (transpose to the stabilization operator)
    MatrixType aux_oss_operator(matrix_size, matrix_size);
    CalculateOrthogonalSubScalesOperator(
        aux_oss_operator,
        rCurrentProcessInfo);

    // Call the element OSS lumped projection operator
    // Note that this calculates the lumped projection operator to be applied to the RHS
    MatrixType aux_lumped_mass_operator(matrix_size, matrix_size);
    CalculateOrthogonalSubScalesLumpedProjectionOperator(
        aux_lumped_mass_operator,
        rCurrentProcessInfo);

    // Assemble the extended LHS
    for (IndexType i = 0; i < matrix_size; ++i) {
        for (IndexType j = 0; j < matrix_size; ++j) {
            rLeftHandSideMatrix(i,j) = aux_stiffnes(i,j);
            rLeftHandSideMatrix(i, matrix_size + j) = -aux_oss_stab_operator(i,j);
            rLeftHandSideMatrix(matrix_size + i, j) = -aux_oss_operator(i,j);
            // rLeftHandSideMatrix(matrix_size + j, i) = aux_oss_stab_operator(i,j);
            rLeftHandSideMatrix(matrix_size + i, matrix_size + j) = -aux_lumped_mass_operator(i,j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType block_size_full = 2.0 * block_size;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType matrix_size_full = block_size_full * n_nodes;

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size_full) {
        rRightHandSideVector.resize(matrix_size_full, false);
    }
    rRightHandSideVector.clear();

    // Call the base element to get the standard residual vector
    // Note that this already includes the substraction of the OSS projections
    VectorType aux_rhs(block_size * n_nodes);
    BaseType::CalculateRightHandSide(
        aux_rhs,
        rCurrentProcessInfo);

    // Call the element OSS operator
    // Note that this calculates the OSS operator to be applied to the RHS (transpose to the stabilization operator)
    MatrixType aux_oss_operator(matrix_size, matrix_size);
    CalculateOrthogonalSubScalesOperator(
        aux_oss_operator,
        rCurrentProcessInfo);

    // Call the element OSS lumped projection operator
    // Note that this calculates the lumped projection operator to be applied to the RHS
    MatrixType aux_lumped_mass_operator(matrix_size, matrix_size);
    CalculateOrthogonalSubScalesLumpedProjectionOperator(
        aux_lumped_mass_operator,
        rCurrentProcessInfo);

    // Get nodal values auxiliary arrays
    VectorType aux_unk(matrix_size);
    VectorType aux_projection(matrix_size);
    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_node = r_geometry[i_node];
            const auto& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            const auto& r_disp_proj = r_node.FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
            aux_unk[i_node * block_size] = r_disp[0];
            aux_unk[i_node * block_size + 1] = r_disp[1];
            aux_unk[i_node * block_size + 2] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
            aux_projection[i_node * block_size] = r_disp_proj[0];
            aux_projection[i_node * block_size + 1] = r_disp_proj[1];
            aux_projection[i_node * block_size + 2] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
        }
    } else {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_node = r_geometry[i_node];
            const auto& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            const auto& r_disp_proj = r_node.FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
            aux_unk[i_node * block_size] = r_disp[0];
            aux_unk[i_node * block_size + 1] = r_disp[1];
            aux_unk[i_node * block_size + 2] = r_disp[2];
            aux_unk[i_node * block_size + 3] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
            aux_projection[i_node * block_size] = r_disp_proj[0];
            aux_projection[i_node * block_size + 1] = r_disp_proj[1];
            aux_projection[i_node * block_size + 2] = r_disp_proj[2];
            aux_projection[i_node * block_size + 3] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
        }
    }

    // Calculate auxiliary arrays
    const VectorType oss_stab_times_unk = prod(aux_oss_operator, aux_unk);
    const VectorType lumped_mass_times_proj = prod(aux_lumped_mass_operator, aux_projection);

    // Assemble the extended modal analysis RHS
    for (IndexType i = 0; i < matrix_size; ++i) {
        rRightHandSideVector(i) = aux_rhs(i); // Note that this already includes the projections substraction
        rRightHandSideVector(matrix_size + i) = lumped_mass_times_proj(i) + oss_stab_times_unk(i);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_prop = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size_full = 2.0 * block_size * n_nodes;

    // Check LHS size and initialize
    if (rMassMatrix.size1() != matrix_size_full || rMassMatrix.size2() != matrix_size_full) {
        rMassMatrix.resize(matrix_size_full, matrix_size_full, false);
    }
    rMassMatrix.clear();

    // Calculate the consistent mass matrix
    // Note that we cannot use the base element implementation as this needs to include the zeros in the projections rows
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

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::GetSecondDerivativesVector(
    Vector &rValues,
    int Step) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType block_size_full = 2 * block_size;
    const SizeType matrix_size_full = block_size_full * n_nodes;

    if (rValues.size() != matrix_size_full) {
        rValues.resize(matrix_size_full, false);
    }
    rValues.clear();

    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const SizeType index = i_node * block_size;
        const auto& r_acc = r_geometry[i_node].FastGetSolutionStepValue(ACCELERATION, Step);
        for(IndexType d = 0; d < dim; ++d) {
            rValues[index + d] = r_acc[d];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int SmallDisplacementMixedVolumetricStrainOssNonLinearElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto& r_geometry = this->GetGeometry();
    for ( IndexType i = 0; i < r_geometry.size(); i++ ) {
        const NodeType& r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUMETRIC_STRAIN_PROJECTION, r_node)
        KRATOS_CHECK_DOF_IN_NODE(VOLUMETRIC_STRAIN_PROJECTION, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT_PROJECTION, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_PROJECTION_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_PROJECTION_Y, r_node)
        if (rCurrentProcessInfo[DOMAIN_SIZE] == 3) {
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_PROJECTION_Z, r_node)
        }
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedVolumetricStrainOssNonLinearElement::GetSpecifications() const
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
        "required_variables"         : ["DISPLACEMENT","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION","VOLUMETRIC_STRAIN_PROJECTION"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Quadrilateral2D4", "Tetrahedra3D4","Hexahedra3D8"],
        "element_integrates_in_time" : true,
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
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_PROJECTION_X", "DISPLACEMENT_PROJECTION_Y", "VOLUMETRIC_STRAIN", "VOLUMETRIC_STRAIN_PROJECTION"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "DISPLACEMENT_PROJECTION_X", "DISPLACEMENT_PROJECTION_Y", "DISPLACEMENT_PROJECTION_Z", "VOLUMETRIC_STRAIN", "VOLUMETRIC_STRAIN_PROJECTION"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateOrthogonalSubScalesOperator(
    MatrixType &rOrthogonalSubScalesOperator,
    const ProcessInfo &rProcessInfo) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rOrthogonalSubScalesOperator.size1() != matrix_size || rOrthogonalSubScalesOperator.size2() != matrix_size) {
        rOrthogonalSubScalesOperator.resize(matrix_size, matrix_size, false);
    }

    // Initialize the OSS operator with zeros
    rOrthogonalSubScalesOperator.clear();

    // Create the kinematics container
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set auxiliary arrays
    Vector voigt_identity =  ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }
    Vector G_j(dim);
    Vector psi_j(dim);
    Matrix B_j(strain_size, dim);

    // Calculate the anisotropy tensor products
    const Vector m_T = prod(voigt_identity, GetAnisotropyTensor());

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const double thickness = (dim == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate tau_1 stabilization constant
        const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rProcessInfo);

        // Calculate tau_2 stabilization constant
        const double tau_2 = CalculateTau2(constitutive_variables);

        const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
        const double aux_w_kappa_tau_1 = w_gauss * bulk_modulus * tau_1;
        const double aux_w_kappa_tau_2 = w_gauss * bulk_modulus * tau_2;

        for (IndexType j = 0; j < n_nodes; ++j) {
            const double N_j = kinematic_variables.N[j];
            noalias(G_j) = row(kinematic_variables.DN_DX, j);
            for (IndexType k = 0;  k < strain_size; ++k) {
                for (IndexType l = 0; l < dim; ++l) {
                    B_j(k,l) = kinematic_variables.B(k, j * dim + l);
                }
            }
            noalias(psi_j) = prod(trans(voigt_identity), Matrix(prod(GetAnisotropyTensor(), B_j)));

            for (IndexType i = 0; i < n_nodes; ++i) {
                const double N_i = kinematic_variables.N[i];
                rOrthogonalSubScalesOperator(i * block_size + dim, j * block_size + dim) -= aux_w_kappa_tau_2 * N_i * N_j; // Note that we multiply by kappa*tau_2 to symmetrize
                for (IndexType d = 0; d < dim; ++d) {
                    rOrthogonalSubScalesOperator(i * block_size + d, j * block_size + dim) -= aux_w_kappa_tau_1 * N_i * G_j[d]; // Note that we multiply by tau_1 to symmetrize
                    rOrthogonalSubScalesOperator(i * block_size + dim, j * block_size + d) += aux_w_kappa_tau_2 * N_i * psi_j[d]; // Note that we multiply by kappa*tau_2 to symmetrize
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateOrthogonalSubScalesStabilizationOperator(
    MatrixType &rOrthogonalSubScalesStabilizationOperator,
    const ProcessInfo &rProcessInfo) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rOrthogonalSubScalesStabilizationOperator.size1() != matrix_size || rOrthogonalSubScalesStabilizationOperator.size2() != matrix_size) {
        rOrthogonalSubScalesStabilizationOperator.resize(matrix_size, matrix_size, false);
    }

    // Initialize the OSS operator with zeros
    rOrthogonalSubScalesStabilizationOperator.clear();

    // Create the kinematics container
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set the auxiliary Gauss point variable container
    GaussPointAuxiliaryVariables gauss_point_auxiliary_variables(this, dim, strain_size);

    // Calculate and assemble the integration points contribution
    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the auxiliary Gauss point variables
        CalculateGaussPointAuxiliaryVariables(
            gauss_point_auxiliary_variables,
            kinematic_variables,
            constitutive_variables,
            rProcessInfo,
            i_gauss);

        // Assemble the current Gauss point OSS projection operator
        CalculateOssStabilizationOperatorGaussPointContribution(
            rOrthogonalSubScalesStabilizationOperator,
            kinematic_variables,
            gauss_point_auxiliary_variables);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::CalculateOrthogonalSubScalesLumpedProjectionOperator(
    MatrixType& rOrthogonalSubScalesLumpedProjectionOperator,
    const ProcessInfo& rProcessInfo) const
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rOrthogonalSubScalesLumpedProjectionOperator.size1() != matrix_size || rOrthogonalSubScalesLumpedProjectionOperator.size2() != matrix_size) {
        rOrthogonalSubScalesLumpedProjectionOperator.resize(matrix_size, matrix_size, false);
    }

    // Initialize the OSS operator with zeros
    rOrthogonalSubScalesLumpedProjectionOperator.clear();

    // Create the kinematics container
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set auxiliary arrays
    Vector voigt_identity =  ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d) {
        voigt_identity[d] = 1.0;
    }

    // Calculate the anisotropy tensor products
    const Vector m_T = prod(voigt_identity, GetAnisotropyTensor());

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const double thickness = (dim == 2 && GetProperties().Has(THICKNESS)) ? GetProperties()[THICKNESS] : 1.0;
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = thickness * kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate tau_1 stabilization constant
        const double tau_1 = CalculateTau1(m_T, kinematic_variables, constitutive_variables, rProcessInfo);

        // Calculate tau_2 stabilization constant
        const double tau_2 = CalculateTau2(constitutive_variables);

        double sum_N_j;
        const double bulk_modulus = CalculateBulkModulus(constitutive_variables.D);
        // const double aux_w_kappa_tau_1 = w_gauss * bulk_modulus * tau_1;
        const double aux_w_tau_1 = w_gauss * tau_1;
        const double aux_w_kappa_tau_2 = w_gauss * bulk_modulus * tau_2;
        for (IndexType i = 0; i < n_nodes; ++i) {
            sum_N_j = 0.0;
            const double N_i = kinematic_variables.N[i];
            for (IndexType j = 0; j < n_nodes; ++j) {
                sum_N_j += kinematic_variables.N[j];
            }
            for (IndexType d = 0; d < dim; ++d) {
                rOrthogonalSubScalesLumpedProjectionOperator(i * block_size + d, i * block_size + d) += aux_w_tau_1 * N_i * sum_N_j; // Note that we multiply by tau_1 to scale the projection with the symmetrization factors (see CalculateOrthogonalSubScalesOperator)
            }
            rOrthogonalSubScalesLumpedProjectionOperator(i * block_size + dim, i * block_size + dim) -= aux_w_kappa_tau_2 * N_i * sum_N_j; // Note that we multiply by kappa*tau_2 to scale the projection with the symmetrization factors (see CalculateOrthogonalSubScalesOperator)
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainOssNonLinearElement::BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainOssNonLinearElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainOssNonLinearElement::BaseType);
}

} // Namespace Kratos
