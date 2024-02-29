// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "displacement_shifted_boundary_condition.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    :Condition(NewId, ThisNodes)
{
}

DisplacementShiftedBoundaryCondition::~DisplacementShiftedBoundaryCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer DisplacementShiftedBoundaryCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DisplacementShiftedBoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer DisplacementShiftedBoundaryCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DisplacementShiftedBoundaryCondition>(NewId, pGeom, pProperties);
}

void DisplacementShiftedBoundaryCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the equation ids. vector
    const auto &r_geometry = this->GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
    const SizeType local_size = n_dim * n_nodes;
    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    const SizeType disp_x_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    if (n_dim == 2) {
        for (std::size_t i = 0; i < n_nodes; ++i){
            rResult[i*n_dim] = r_geometry[i].GetDof(DISPLACEMENT_X, disp_x_pos).EquationId();
            rResult[i*n_dim +1] = r_geometry[i].GetDof(DISPLACEMENT_Y, disp_x_pos + 1).EquationId();
        }
    } else {
        for (std::size_t i = 0; i < n_nodes; ++i){
            rResult[i*n_dim] = r_geometry[i].GetDof(DISPLACEMENT_X, disp_x_pos).EquationId();
            rResult[i*n_dim + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y, disp_x_pos + 1).EquationId();
            rResult[i*n_dim + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z, disp_x_pos + 2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

void DisplacementShiftedBoundaryCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the DOFs vector
    const auto& r_geometry = this->GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
    const SizeType local_size = n_dim * n_nodes;
    if (rConditionalDofList.size() != local_size){
        rConditionalDofList.resize(local_size);
    }

    // Fill the DOFs vector from the condition nodes
    if (n_dim == 2) {
        for (std::size_t i = 0; i < n_nodes; ++i) {
            rConditionalDofList[i*n_dim] = r_geometry[i].pGetDof(DISPLACEMENT_X);
            rConditionalDofList[i*n_dim + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
        }
    } else {
        for (std::size_t i = 0; i < n_nodes; ++i) {
            rConditionalDofList[i*n_dim] = r_geometry[i].pGetDof(DISPLACEMENT_X);
            rConditionalDofList[i*n_dim + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
            rConditionalDofList[i*n_dim + 2] = r_geometry[i].pGetDof(DISPLACEMENT_Z);
        }
    }

        KRATOS_CATCH("")
    }

void DisplacementShiftedBoundaryCondition::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    // Get the parent constitutive law from the properties and clone them to be current properties
    // Note that we are assuming that the parent properties were assigned as properties when creating the condition
    const auto& r_prop = GetProperties();
    const auto p_cons_law = r_prop[CONSTITUTIVE_LAW];
    const auto& r_N = this->GetValue(SHAPE_FUNCTIONS_VECTOR);
    if (p_cons_law != nullptr) {
        this->SetValue(CONSTITUTIVE_LAW, p_cons_law->Clone());
        this->GetValue(CONSTITUTIVE_LAW)->InitializeMaterial(r_prop, GetGeometry(), r_N);
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified in condition " << this->Id() << "." << std::endl;
    }
}

void DisplacementShiftedBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get problem dimensions
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Check (and resize) LHS and RHS matrix
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t local_size = n_nodes * n_dim;
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size, false);
    }
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(local_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

    // Get meshless geometry data
    const double w = GetValue(INTEGRATION_WEIGHT);
    const auto& r_N = GetValue(SHAPE_FUNCTIONS_VECTOR);
    const auto& r_DN_DX = GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);
    array_1d<double,3> normal = GetValue(NORMAL);
    normal /= norm_2(normal);

    // Get unknown values
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (std::size_t d = 0; d < n_dim; ++d) {
            unknown_values(i_node*n_dim + d) = r_disp[d];
        }
    }

    // Calculate the material response to get the constitutive matrix
    const auto p_cons_law = this->GetValue(CONSTITUTIVE_LAW);
    const SizeType strain_size = p_cons_law->GetStrainSize();

    Vector stress_vect(strain_size);
    Matrix B_mat(strain_size, local_size);
    Matrix C_mat(strain_size, strain_size);
    StructuralMechanicsElementUtilities::CalculateB(*this, r_DN_DX, B_mat);
    Vector strain_vect = prod(B_mat, unknown_values);

    ConstitutiveLaw::Parameters cons_law_values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    cons_law_values.SetShapeFunctionsValues(r_N);
    cons_law_values.SetStressVector(stress_vect);
    cons_law_values.SetStrainVector(strain_vect);
    cons_law_values.SetConstitutiveMatrix(C_mat);

    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

    // Get Dirichlet BC imposition data
    const double h = GetValue(ELEMENT_H);
    const auto& r_bc_val = GetValue(DISPLACEMENT);
    const double gamma = rCurrentProcessInfo[PENALTY_COEFFICIENT];

    // Calculate the Nitsche BC imposition contribution
    // 1. Add Nitsche penalty term
    double aux_1;
    double aux_2;
    const double rho_C = norm_frobenius(C_mat); //TODO: GS uses the spectral radius in here
    const double aux_weight = w * gamma * rho_C / h;
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        aux_1 = aux_weight * r_N[i_node];
        for (std::size_t d = 0; d < n_dim; ++d) {
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                aux_2 = aux_1 * r_N[j_node];
                rLeftHandSideMatrix(i_node*n_dim + d, j_node*n_dim + d) += aux_2;
                rRightHandSideVector(i_node*n_dim + d) -= aux_2 * unknown_values(j_node*n_dim + d);
            }
            rRightHandSideVector(i_node*n_dim + d) += aux_1 * r_bc_val[d];
        }
    }

    // 2. Add Nitsche stabilization term
    Matrix aux_N(n_dim, local_size);
    Matrix transB_C_proj(local_size, n_dim);
    CalculateAuxShapeFunctionsMatrix(r_N, aux_N);
    CalculateBtransCProjectionLinearisation(C_mat, B_mat, normal, transB_C_proj);
    const Matrix stab_lhs = prod(transB_C_proj, aux_N);
    Vector stab_rhs_bc = ZeroVector(local_size);
    for (IndexType i = 0; i < local_size; ++i) {
        for (IndexType d = 0; d < n_dim; ++d) {
            stab_rhs_bc[i] += transB_C_proj(i,d) * r_bc_val[d];
        }
    }
    const Vector stab_rhs_unk = prod(stab_lhs, unknown_values);
    rLeftHandSideMatrix += w * stab_lhs;
    rRightHandSideVector += w * stab_rhs_bc;
    rRightHandSideVector -= w * stab_rhs_unk;

    KRATOS_CATCH("")
}

void DisplacementShiftedBoundaryCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get problem dimensions
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Check (and resize) LHS matrix
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t local_size = n_nodes * n_dim;
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_nodes,n_nodes);

    // Get meshless geometry data
    const double w = GetValue(INTEGRATION_WEIGHT);
    const auto& r_N = GetValue(SHAPE_FUNCTIONS_VECTOR);
    const auto& r_DN_DX = GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);
    array_1d<double,3> normal = GetValue(NORMAL);
    normal /= norm_2(normal);

    // Get unknown values
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (std::size_t d = 0; d < n_dim; ++d) {
            unknown_values(i_node*n_dim + d) = r_disp[d];
        }
    }

    // Calculate the material response to get the constitutive matrix
    const auto p_cons_law = this->GetValue(CONSTITUTIVE_LAW);
    const SizeType strain_size = p_cons_law->GetStrainSize();

    Vector stress_vect(strain_size);
    Matrix B_mat(strain_size, local_size);
    Matrix C_mat(strain_size, strain_size);
    StructuralMechanicsElementUtilities::CalculateB(*this, r_DN_DX, B_mat);
    Vector strain_vect = prod(B_mat, unknown_values);

    ConstitutiveLaw::Parameters cons_law_values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    cons_law_values.SetShapeFunctionsValues(r_N);
    cons_law_values.SetStressVector(stress_vect);
    cons_law_values.SetStrainVector(strain_vect);
    cons_law_values.SetConstitutiveMatrix(C_mat);

    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

    // Get Dirichlet BC imposition data
    const double h = GetValue(ELEMENT_H);
    const double gamma = rCurrentProcessInfo[PENALTY_COEFFICIENT];

    // Calculate the Nitsche BC imposition contribution
    // 1. Add Nitsche penalty term
    double aux_1;
    double aux_2;
    const double rho_C = norm_frobenius(C_mat); //TODO: GS uses the spectral radius in here
    const double aux_weight = w * gamma * rho_C / h;
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        aux_1 = aux_weight * r_N[i_node];
        for (std::size_t d = 0; d < n_dim; ++d) {
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                aux_2 = aux_1 * r_N[j_node];
                rLeftHandSideMatrix(i_node*n_dim + d, j_node*n_dim + d) += aux_2;
            }
        }
    }

    // 2. Add Nitsche stabilization term
    Matrix aux_N(n_dim, local_size);
    Matrix transB_C_proj(local_size, n_dim);
    CalculateAuxShapeFunctionsMatrix(r_N, aux_N);
    CalculateBtransCProjectionLinearisation(C_mat, B_mat, normal, transB_C_proj);
    const Matrix stab_lhs = prod(transB_C_proj, aux_N);
    rLeftHandSideMatrix += w * stab_lhs;

    KRATOS_CATCH("")
}

void DisplacementShiftedBoundaryCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get problem dimensions
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Check (and resize) LHS and RHS matrix
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t local_size = n_nodes * n_dim;
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(local_size);

    // Get meshless geometry data
    const double w = GetValue(INTEGRATION_WEIGHT);
    const auto& r_N = GetValue(SHAPE_FUNCTIONS_VECTOR);
    const auto& r_DN_DX = GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);
    array_1d<double,3> normal = GetValue(NORMAL);
    normal /= norm_2(normal);

    // Get unknown values
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (std::size_t d = 0; d < n_dim; ++d) {
            unknown_values(i_node*n_dim + d) = r_disp[d];
        }
    }

    // Calculate the material response to get the constitutive matrix
    const auto p_cons_law = this->GetValue(CONSTITUTIVE_LAW);
    const SizeType strain_size = p_cons_law->GetStrainSize();

    Vector stress_vect(strain_size);
    Matrix B_mat(strain_size, local_size);
    Matrix C_mat(strain_size, strain_size);
    StructuralMechanicsElementUtilities::CalculateB(*this, r_DN_DX, B_mat);
    Vector strain_vect = prod(B_mat, unknown_values);

    ConstitutiveLaw::Parameters cons_law_values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    cons_law_values.SetShapeFunctionsValues(r_N);
    cons_law_values.SetStressVector(stress_vect);
    cons_law_values.SetStrainVector(strain_vect);
    cons_law_values.SetConstitutiveMatrix(C_mat);

    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

    // Get Dirichlet BC imposition data
    const double h = GetValue(ELEMENT_H);
    const auto& r_bc_val = GetValue(DISPLACEMENT);
    const double gamma = rCurrentProcessInfo[PENALTY_COEFFICIENT];

    // Calculate the Nitsche BC imposition contribution
    // 1. Add Nitsche penalty term
    double aux_1;
    double aux_2;
    const double rho_C = norm_frobenius(C_mat); //TODO: GS uses the spectral radius in here
    const double aux_weight = w * gamma * rho_C / h;
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        aux_1 = aux_weight * r_N[i_node];
        for (std::size_t d = 0; d < n_dim; ++d) {
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                aux_2 = aux_1 * r_N[j_node];
                rRightHandSideVector(i_node*n_dim + d) -= aux_2 * unknown_values(j_node*n_dim + d);
            }
            rRightHandSideVector(i_node*n_dim + d) += aux_1 * r_bc_val[d];
        }
    }

    // 2. Add Nitsche stabilization term
    Matrix aux_N(n_dim, local_size);
    Matrix transB_C_proj(local_size, n_dim);
    CalculateAuxShapeFunctionsMatrix(r_N, aux_N);
    CalculateBtransCProjectionLinearisation(C_mat, B_mat, normal, transB_C_proj);
    const Matrix stab_lhs = prod(transB_C_proj, aux_N);
    Vector stab_rhs_bc = ZeroVector(local_size);
    for (IndexType i = 0; i < local_size; ++i) {
        for (IndexType d = 0; d < n_dim; ++d) {
            stab_rhs_bc[i] += transB_C_proj(i,d) * r_bc_val[d];
        }
    }
    const Vector stab_rhs_unk = prod(stab_lhs, unknown_values);
    rRightHandSideVector += w * stab_rhs_bc;
    rRightHandSideVector -= w * stab_rhs_unk;
    KRATOS_CATCH("")
}

int DisplacementShiftedBoundaryCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("")
}

void DisplacementShiftedBoundaryCondition::CalculateBtransCProjectionLinearisation(
    const Matrix& rC,
    const Matrix& rB,
    const array_1d<double, 3>& rUnitNormal,
    Matrix& rAuxMat)
{
    const SizeType local_size = rB.size2();
    const SizeType strain_size = rB.size1();
    const Matrix aux_transBC = prod(trans(rB), rC);
    if (strain_size == 3) {
        for (std::size_t j = 0; j < local_size; ++j) {
            rAuxMat(j,0) = rUnitNormal[0]*aux_transBC(j,0) + rUnitNormal[1]*aux_transBC(j,2);
        }
        for (std::size_t j = 0; j < local_size; ++j) {
            rAuxMat(j,1) = rUnitNormal[0]*aux_transBC(j,2) + rUnitNormal[1]*aux_transBC(j,1);
        }
    } else {
        for (std::size_t j = 0; j < local_size; ++j) {
            rAuxMat(j,0) = rUnitNormal[0]*aux_transBC(j,0) + rUnitNormal[1]*aux_transBC(j,3) + rUnitNormal[2]*aux_transBC(j,5);
        }
        for (std::size_t j = 0; j < local_size; ++j) {
            rAuxMat(j,1) = rUnitNormal[0]*aux_transBC(j,3) + rUnitNormal[1]*aux_transBC(j,1) + rUnitNormal[2]*aux_transBC(j,4);
        }
        for (std::size_t j = 0; j < local_size; ++j) {
            rAuxMat(j,2) = rUnitNormal[0]*aux_transBC(j,5) + rUnitNormal[1]*aux_transBC(j,4) + rUnitNormal[2]*aux_transBC(j,2);
        }
    }
}

void DisplacementShiftedBoundaryCondition::CalculateAuxShapeFunctionsMatrix(
    const Vector& rN,
    Matrix& rAuxMat)
{
    const SizeType n_nodes = rN.size();
    const SizeType n_dim = rAuxMat.size1();
    const SizeType local_size = rAuxMat.size2();
    rAuxMat = ZeroMatrix(n_dim, local_size);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        for (IndexType d = 0; d < n_dim; ++d) {
            rAuxMat(d, i_node*n_dim + d) = rN[i_node];
        }
    }
}

// Input and Output ///////////////////////////////////////////////////////////

std::string DisplacementShiftedBoundaryCondition::Info() const
{
    std::stringstream buffer;
    buffer << "DisplacementShiftedBoundaryCondition #" << Id();
    return buffer.str();
}

void DisplacementShiftedBoundaryCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DisplacementShiftedBoundaryCondition #" << Id();
}

void DisplacementShiftedBoundaryCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "DisplacementShiftedBoundaryCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition():
    Condition()
{
}

void DisplacementShiftedBoundaryCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void DisplacementShiftedBoundaryCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
