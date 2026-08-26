// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Ricky Aristio
//

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "utilities/integration_utilities.h"
#include "utilities/math_utils.h"

// Application includes
#include "displacement_shell_shifted_boundary_condition.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

DisplacementShellShiftedBoundaryCondition::DisplacementShellShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

DisplacementShellShiftedBoundaryCondition::DisplacementShellShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

DisplacementShellShiftedBoundaryCondition::DisplacementShellShiftedBoundaryCondition(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    :Condition(NewId, ThisNodes)
{
}

DisplacementShellShiftedBoundaryCondition::~DisplacementShellShiftedBoundaryCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer DisplacementShellShiftedBoundaryCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DisplacementShellShiftedBoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer DisplacementShellShiftedBoundaryCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DisplacementShellShiftedBoundaryCondition>(NewId, pGeom, pProperties);
}

void DisplacementShellShiftedBoundaryCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the equation ids. vector
    const auto &r_geometry = this->GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType n_dim = 6;
    const SizeType local_size = n_dim * n_nodes;

    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    const SizeType disp_x_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    for (std::size_t i = 0; i < n_nodes; ++i){
        rResult[i*n_dim] = r_geometry[i].GetDof(DISPLACEMENT_X, disp_x_pos).EquationId();
        rResult[i*n_dim +1] = r_geometry[i].GetDof(DISPLACEMENT_Y, disp_x_pos + 1).EquationId();
        rResult[i*n_dim +2] = r_geometry[i].GetDof(DISPLACEMENT_Z, disp_x_pos + 2).EquationId();
        rResult[i*n_dim +3] = r_geometry[i].GetDof(ROTATION_X, disp_x_pos + 3).EquationId();
        rResult[i*n_dim +4] = r_geometry[i].GetDof(ROTATION_Y, disp_x_pos + 4).EquationId();
        rResult[i*n_dim +5] = r_geometry[i].GetDof(ROTATION_Z, disp_x_pos + 5).EquationId();
    }

    KRATOS_CATCH("")
}

void DisplacementShellShiftedBoundaryCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the DOFs vector
    const auto& r_geometry = this->GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType n_dim = 6;
    const SizeType local_size = n_dim * n_nodes;
    if (rConditionalDofList.size() != local_size){
        rConditionalDofList.resize(local_size);
    }

    // Fill the DOFs vector from the condition nodes
    for (std::size_t i = 0; i < n_nodes; ++i) {
        rConditionalDofList[i*n_dim] = r_geometry[i].pGetDof(DISPLACEMENT_X);
        rConditionalDofList[i*n_dim + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
        rConditionalDofList[i*n_dim + 2] = r_geometry[i].pGetDof(DISPLACEMENT_Z);
        rConditionalDofList[i*n_dim + 3] = r_geometry[i].pGetDof(ROTATION_X);
        rConditionalDofList[i*n_dim + 4] = r_geometry[i].pGetDof(ROTATION_Y);
        rConditionalDofList[i*n_dim + 5] = r_geometry[i].pGetDof(ROTATION_Z);
    }

    KRATOS_CATCH("")
}

void DisplacementShellShiftedBoundaryCondition::Initialize(const ProcessInfo &rCurrentProcessInfo)
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

void DisplacementShellShiftedBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get problem dimensions
    const SizeType n_dim = 6;

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
    Vector unknown_values = ZeroVector(local_size);

    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_rot = r_geometry[i_node].FastGetSolutionStepValue(ROTATION);
        for (std::size_t d = 0; d < 3; ++d) {
            unknown_values(i_node*6 + d) = r_disp[d];
            unknown_values(i_node*6 + d + 3) = r_rot[d];
        }
    }

    // Calculate the material response to get the constitutive matrix
    const auto p_cons_law = this->GetValue(CONSTITUTIVE_LAW);
    const SizeType strain_size = p_cons_law->GetStrainSize();

    Vector stress_vect(strain_size);
    Matrix B_mat(strain_size, 2*n_nodes);
    Matrix C_mat(strain_size, strain_size);

    // ConstitutiveLaw::VoigtSizeMatrixType D;
    Matrix D(8, 8);
    Matrix B_mat_shell(8, 6*n_nodes);
    noalias(D) = ZeroMatrix(8, 8);
    noalias(B_mat_shell) = ZeroMatrix(8, 6*n_nodes);

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

    double thickness = this->GetProperties()[THICKNESS];

    for(std::size_t i = 0; i < C_mat.size1(); ++i)
    {
        for(std::size_t j = 0; j < C_mat.size1(); ++j)
        {
            D(i,j) = C_mat(i,j) * thickness;
            D(i + 3,j +3) = C_mat(i,j) / 12.0 * thickness * thickness * thickness;
        }
    }

    //determine longest side length (eqn 22) - for shear correction
    Vector P12 = Vector(r_geometry[0] - r_geometry[1]);
    Vector P13 = Vector(r_geometry[0] - r_geometry[2]);
    Vector P23 = Vector(r_geometry[1] - r_geometry[2]);

    double h_e = std::sqrt(inner_prod(P12, P12));
    double edge_length = std::sqrt(inner_prod(P13, P13));
    if (edge_length > h_e) {
        h_e = edge_length;
    }
    edge_length = std::sqrt(inner_prod(P23, P23));
    if (edge_length > h_e) {
        h_e = edge_length;
    }

    double alpha = 0.1;
    double hMean = thickness;

    // Write Stenberg shear stabilisation coefficient
    double shearStabilisation = (hMean*hMean)
                              / (hMean*hMean + alpha*h_e*h_e);
    
    D(6,6) = D(2,2) * 5.0 / 6.0 * shearStabilisation; //plane strain!
    D(7,7) = D(2,2) * 5.0 / 6.0 * shearStabilisation; //stenberg is missing

    // Local in-plane orthonormal frame of the parent background element 
    const array_1d<double,3>& local_e1 = GetValue(LOCAL_AXIS_1);
    const array_1d<double,3>& local_e2 = GetValue(LOCAL_AXIS_2);
    array_1d<double,3> local_e3;
    MathUtils<double>::CrossProduct(local_e3, local_e1, local_e2);

    for ( IndexType i = 0; i < n_nodes; ++i ) {
        const IndexType initial_index = i*6;
        for (std::size_t k = 0; k < 3; ++k) {
            B_mat_shell(0, initial_index + k) = r_DN_DX(i, 0) * local_e1[k];
            B_mat_shell(1, initial_index + k) = r_DN_DX(i, 1) * local_e2[k];
            B_mat_shell(2, initial_index + k) = r_DN_DX(i, 1) * local_e1[k] + r_DN_DX(i, 0) * local_e2[k];

            B_mat_shell(3, initial_index + 3 + k) = r_DN_DX(i, 0) * local_e2[k];
            B_mat_shell(4, initial_index + 3 + k) = -r_DN_DX(i, 1) * local_e1[k];
            B_mat_shell(5, initial_index + 3 + k) = -r_DN_DX(i, 0) * local_e1[k] + r_DN_DX(i, 1) * local_e2[k];
        }
    }

    // Shear rows: match whichever shear formulation the interior element actually uses
    // (SHELL_SHEAR_FORMULATION, default DSG=0).
    const int shell_shear_formulation = GetProperties().Has(SHELL_SHEAR_FORMULATION)
        ? GetProperties()[SHELL_SHEAR_FORMULATION]
        : 0;
    if (shell_shear_formulation == 0 && Has(SBM_DSG_SHEAR_OPERATOR)) {
        const auto& r_dsg_shear = GetValue(SBM_DSG_SHEAR_OPERATOR);
        for (std::size_t j = 0; j < 6*n_nodes; ++j) {
            B_mat_shell(6, j) = r_dsg_shear(0, j);
            B_mat_shell(7, j) = r_dsg_shear(1, j);
        }
    } else {
        // Naive (frame-correct) formula: used whenever SHELL_SHEAR_FORMULATION selects
        // naive/CST shear instead of DSG
        for ( IndexType i = 0; i < n_nodes; ++i ) {
            const IndexType initial_index = i*6;
            for (std::size_t k = 0; k < 3; ++k) {
                B_mat_shell(6, initial_index + k)     = r_DN_DX(i, 0) * local_e3[k];
                B_mat_shell(6, initial_index + 3 + k) = r_N[i] * local_e2[k];
                B_mat_shell(7, initial_index + k)     = r_DN_DX(i, 1) * local_e3[k];
                B_mat_shell(7, initial_index + 3 + k) = -r_N[i] * local_e1[k];
            }
        }
    }

    // Get Dirichlet BC imposition data
    const double h = GetValue(ELEMENT_H);
    const auto& r_bc_val = GetValue(DISPLACEMENT);
    const auto& r_bc_val_rot = GetValue(ROTATION);
    Vector bc_values = ZeroVector(6);
    for (std::size_t d = 0; d < 3; ++d) {
        bc_values[d] = r_bc_val[d];
        bc_values[d + 3] = r_bc_val_rot[d];
    }
    const double gamma = GetProperties().Has(PENALTY_COEFFICIENT)
        ? GetProperties()[PENALTY_COEFFICIENT]
        : 0.0;

    // Get the per-DOF constrained mask (ux,uy,uz,rx,ry,rz)
    array_1d<double, 6> constrained_dofs;
    if (GetProperties().Has(SBM_CONSTRAINED_DOFS)) {
        const auto& r_constrained_dofs = GetProperties()[SBM_CONSTRAINED_DOFS];
        KRATOS_ERROR_IF(r_constrained_dofs.size() != 6)
            << "SBM_CONSTRAINED_DOFS must have size 6 (ux,uy,uz,rx,ry,rz) in Properties " << GetProperties().Id() << std::endl;
        for (std::size_t d = 0; d < 6; ++d) {
            constrained_dofs[d] = r_constrained_dofs[d];
        }
    } else {
        for (std::size_t d = 0; d < 6; ++d) {
            constrained_dofs[d] = 1.0;
        }
    }

    // Surrogate-boundary flux formulation choice 
    const int sbm_formulation_type = GetProperties().Has(SBM_FORMULATION_TYPE)
        ? GetProperties()[SBM_FORMULATION_TYPE]
        : 0;
    const double mls_only_flag = (sbm_formulation_type == 0) ? 1.0 : 0.0;

    // Calculate the Nitsche BC imposition contribution
    // 1. Add Nitsche penalty term (translation + rotation)
    double aux_1;
    double aux_2;

    // Per-DOF-type penalty scaling.
    array_1d<double, 6> rho_dof;
    for (std::size_t d = 0; d < 6; ++d) {
        double sum_val = 0.0;
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            const IndexType j = i_node * n_dim + d;
            double val = 0.0;
            for (std::size_t p = 0; p < 8; ++p) {
                for (std::size_t q = 0; q < 8; ++q) {
                    val += B_mat_shell(p, j) * D(p, q) * B_mat_shell(q, j);
                }
            }
            sum_val += val;
        }
        rho_dof[d] = std::sqrt(std::abs(sum_val) / static_cast<double>(n_nodes));
    }

    array_1d<double, 6> aux_weight_dof;
    for (std::size_t d = 0; d < 6; ++d) {
        aux_weight_dof[d] = w * gamma * rho_dof[d] / h * mls_only_flag;
    }

    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        for (std::size_t d = 0; d < 6; ++d) {
            if (constrained_dofs[d] == 0.0) {
                continue;
            }
            aux_1 = aux_weight_dof[d] * r_N[i_node];
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                aux_2 = aux_1 * r_N[j_node];
                rLeftHandSideMatrix(i_node*n_dim + d, j_node*n_dim + d) += aux_2;
                rRightHandSideVector(i_node*n_dim + d) -= aux_2 * unknown_values(j_node*n_dim + d);
            }
            rRightHandSideVector(i_node*n_dim + d) += aux_1 * bc_values[d];
        }
    }

    // 2. Add Nitsche stabilization term (adjoint consistency)
    Matrix aux_N(6, local_size);
    Matrix transB_C_proj = ZeroMatrix(local_size, 6); //first variation traction
    CalculateAuxShapeFunctionsMatrix(r_N, constrained_dofs, aux_N);
    CalculateBtransCProjectionLinearisation(D, B_mat_shell, normal, local_e1, local_e2, local_e3, transB_C_proj);
    const Matrix stab_lhs = prod(transB_C_proj, aux_N);

    Vector stab_rhs_bc = ZeroVector(local_size);
    for (IndexType i = 0; i < local_size; ++i) {
        for (IndexType d = 0; d < n_dim; ++d) {
            stab_rhs_bc[i] += transB_C_proj(i,d) * bc_values[d];
        }
    }

    // Residual correction (current-state piece)
    const Vector stab_rhs_unk = prod(stab_lhs, unknown_values);
    rLeftHandSideMatrix -= w * stab_lhs * mls_only_flag;
    rRightHandSideVector += w * stab_rhs_bc * mls_only_flag;
    rRightHandSideVector += w * stab_rhs_unk * mls_only_flag;

    KRATOS_CATCH("")
}

void DisplacementShellShiftedBoundaryCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Delegate to CalculateLocalSystem so this can never drift out of sync with it.
    VectorType dummy_rhs;
    CalculateLocalSystem(rLeftHandSideMatrix, dummy_rhs, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void DisplacementShellShiftedBoundaryCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Delegate to CalculateLocalSystem so this can never drift out of sync with it.
    MatrixType dummy_lhs;
    CalculateLocalSystem(dummy_lhs, rRightHandSideVector, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

int DisplacementShellShiftedBoundaryCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("")
}

void DisplacementShellShiftedBoundaryCondition::CalculateBtransCProjectionLinearisation(
    const Matrix& rC,
    const Matrix& rB,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rLocalE1,
    const array_1d<double, 3>& rLocalE2,
    const array_1d<double, 3>& rLocalE3,
    Matrix& rAuxMat)
{
    const SizeType local_size = rB.size2();
    const Matrix aux_transBC = prod(trans(rB), rC);

    // rUnitNormal is the (global-frame, curvature-aware) in-plane conormal
    const double normal_e1 = inner_prod(rUnitNormal, rLocalE1);
    const double normal_e2 = inner_prod(rUnitNormal, rLocalE2);

    for (std::size_t j = 0; j < local_size; ++j) {
        const double t1 = normal_e1*aux_transBC(j,0) + normal_e2*aux_transBC(j,2); // membrane traction along e1
        const double t2 = normal_e1*aux_transBC(j,2) + normal_e2*aux_transBC(j,1); // membrane traction along e2
        const double m1 = normal_e1*aux_transBC(j,3) + normal_e2*aux_transBC(j,5); // moment, e2-rotation-conjugate
        const double m2 = normal_e1*aux_transBC(j,5) + normal_e2*aux_transBC(j,4); // moment, e1-rotation-conjugate
        const double q  = normal_e1*aux_transBC(j,6) + normal_e2*aux_transBC(j,7); // shear traction along e3

        for (std::size_t k = 0; k < 3; ++k) {
            rAuxMat(j, k)     = t1*rLocalE1[k] + t2*rLocalE2[k] + q*rLocalE3[k];
            rAuxMat(j, 3 + k) = m2*rLocalE1[k] + m1*rLocalE2[k];
        }
    }
}

void DisplacementShellShiftedBoundaryCondition::CalculateAuxShapeFunctionsMatrix(
    const Vector& rN,
    const array_1d<double, 6>& rConstrainedDofs,
    Matrix& rAuxMat)
{
    const SizeType n_nodes = rN.size();
    const SizeType n_dim = rAuxMat.size1();
    const SizeType local_size = rAuxMat.size2();
    rAuxMat = ZeroMatrix(n_dim, local_size);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        for (IndexType d = 0; d < 6; ++d) {
            if (rConstrainedDofs[d] == 0.0) {
                continue;
            }
            rAuxMat(d, i_node*n_dim + d) = rN[i_node];
        }
    }
}

// Input and Output ///////////////////////////////////////////////////////////

std::string DisplacementShellShiftedBoundaryCondition::Info() const
{
    std::stringstream buffer;
    buffer << "DisplacementShellShiftedBoundaryCondition #" << Id();
    return buffer.str();
}

void DisplacementShellShiftedBoundaryCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DisplacementShellShiftedBoundaryCondition #" << Id();
}

void DisplacementShellShiftedBoundaryCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "DisplacementShellShiftedBoundaryCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

DisplacementShellShiftedBoundaryCondition::DisplacementShellShiftedBoundaryCondition():
    Condition()
{
}

void DisplacementShellShiftedBoundaryCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void DisplacementShellShiftedBoundaryCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
