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
#include "custom_elements/small_displacement_mixed_volumetric_strain_modal_analysis_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedVolumetricStrainModalAnalysisElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainModalAnalysisElement>(
        NewId,
        GetGeometry().Create(ThisNodes),
        pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainModalAnalysisElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainModalAnalysisElement>(
        NewId,
        pGeom,
        pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainModalAnalysisElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedVolumetricStrainModalAnalysisElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainModalAnalysisElement>(
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

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::EquationIdVector(
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

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType block_size = 2.0 * (dim + 1);
    const IndexType dof_size = n_nodes * block_size;

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * block_size] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * block_size + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
            rElementalDofList[i * block_size + 3] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_X);
            rElementalDofList[i * block_size + 4] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_Y);
            rElementalDofList[i * block_size + 5] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN_PROJECTION);
        }
    } else if (dim == 3) {
        for(IndexType i = 0; i < n_nodes; ++i){
            rElementalDofList[i * block_size] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * block_size + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * block_size + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
            rElementalDofList[i * block_size + 4] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_X);
            rElementalDofList[i * block_size + 5] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_Y);
            rElementalDofList[i * block_size + 6] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_PROJECTION_Z);
            rElementalDofList[i * block_size + 7] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN_PROJECTION);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
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

    // Check LHS size and initialize
    if (rLeftHandSideMatrix.size1() != matrix_size_full || rLeftHandSideMatrix.size2() != matrix_size_full) {
        rLeftHandSideMatrix.resize(matrix_size_full, matrix_size_full, false);
    }
    rLeftHandSideMatrix.clear();

    // Calculate the LHS matrix
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // Calculate the RHS vector
    // IndexType count = 0;
    // VectorType aux_data(matrix_size_full);
    // for (const auto& r_node : r_geometry) {
    //     const auto& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
    //     const auto& r_disp_proj = r_node.FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
    //     aux_data[count++] = r_disp[0];
    //     aux_data[count++] = r_disp[1];
    //     if (dim == 3) {
    //         aux_data[count++] = r_disp[2];
    //     }
    //     aux_data[count++] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     aux_data[count++] = r_disp_proj[0];
    //     aux_data[count++] = r_disp_proj[1];
    //     if (dim == 3) {
    //         aux_data[count++] = r_disp_proj[2];
    //     }
    //     aux_data[count++] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
    // }

    VectorType aux_data(matrix_size_full);
    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            IndexType i_row_full = i_node * block_size_full;
            const auto& r_node = r_geometry[i_node];
            const auto& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            const auto& r_disp_proj = r_node.FastGetSolutionStepValue(DISPLACEMENT_PROJECTION);
            aux_data[i_row_full] = r_disp[0];
            aux_data[i_row_full + 1] = r_disp[1];
            aux_data[i_row_full + 2] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
            aux_data[i_row_full + 3] = r_disp_proj[0];
            aux_data[i_row_full + 4] = r_disp_proj[1];
            aux_data[i_row_full + 5] = r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN_PROJECTION);
        }
    } else {
        KRATOS_ERROR << "No 3D implementation yet." << std::endl;
    }

    rRightHandSideVector = prod(rLeftHandSideMatrix, aux_data);
    rRightHandSideVector *= -1.0;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::CalculateLeftHandSide(
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
    MatrixType aux_oss_operator(matrix_size, matrix_size);
    BaseType::CalculateOrthogonalSubScalesOperator(
        aux_oss_operator,
        rCurrentProcessInfo);

    // Call the base element OSS lumped projection operator
    MatrixType aux_lumped_mass_operator(matrix_size, matrix_size);
    BaseType::CalculateOrthogonalSubScalesLumpedProjectionOperator(
        aux_lumped_mass_operator,
        rCurrentProcessInfo);

    // Assemble the extended modal analysis LHS
    for (IndexType i = 0; i < n_nodes; ++i) {
        IndexType i_row = i*block_size;
        IndexType i_row_full = i*block_size_full;
        for (IndexType j = 0; j < n_nodes; ++j) {
            IndexType j_col = j*block_size;
            IndexType j_col_full = j*block_size_full;
            for (IndexType l = 0; l < block_size; ++l) {
                for (IndexType m = 0; m < block_size; ++m) {
                    rLeftHandSideMatrix(i_row_full+l, j_col_full+m) = aux_stiffnes(i_row+l, j_col+m);
                    rLeftHandSideMatrix(i_row_full+l, j_col_full+block_size+m) = aux_oss_operator(i_row+l, j_col+m);
                    rLeftHandSideMatrix(i_row_full+block_size+l, j_col_full+m) = aux_oss_operator(j_col+m, i_row+l);
                    // rLeftHandSideMatrix(i*block_size_full + block_size + l, j*block_size_full + block_size + m) = aux_lumped_mass_operator(i*block_size + l, j*block_size + m);
                }
                rLeftHandSideMatrix(i_row_full+block_size+l, i_row_full+block_size+l) = aux_lumped_mass_operator(i_row+l, i_row+l);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::CalculateRightHandSide(
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

    MatrixType aux_LHS(matrix_size_full, matrix_size_full);
    CalculateLocalSystem(aux_LHS, rRightHandSideVector, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType matrix_size_full = 2.0 * block_size * n_nodes;

    // Check LHS size and initialize
    if (rMassMatrix.size1() != matrix_size_full || rMassMatrix.size2() != matrix_size_full) {
        rMassMatrix.resize(matrix_size_full, matrix_size_full, false);
    }
    rMassMatrix.clear();

    // Call the base element to get the stiffness matrix
    MatrixType aux_mass(matrix_size, matrix_size);
    BaseType::CalculateMassMatrix(
        aux_mass,
        rCurrentProcessInfo);

    // Assemble the extended mass matrix
    for (IndexType i = 0; i < matrix_size; ++i) {
        for (IndexType j = 0; j < matrix_size; ++j) {
            rMassMatrix(i,j) = aux_mass(i,j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::GetSecondDerivativesVector(
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
        const SizeType index = i_node * block_size_full;
        const auto& r_acc = r_geometry[i_node].FastGetSolutionStepValue(ACCELERATION, Step);
        for(IndexType d = 0; d < dim; ++d) {
            rValues[index + d] = r_acc[d];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int SmallDisplacementMixedVolumetricStrainModalAnalysisElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_PROJECTION_Z, r_node)
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedVolumetricStrainModalAnalysisElement::GetSpecifications() const
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
        "required_dofs"              : ["DISPLACEMENT","VOLUMETRIC_STRAIN","DISPLACEMENT_PROJECTION","VOLUMETRIC_STRAIN_PROJECTION"],
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

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainModalAnalysisElement::BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainModalAnalysisElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainModalAnalysisElement::BaseType);
}

} // Namespace Kratos
