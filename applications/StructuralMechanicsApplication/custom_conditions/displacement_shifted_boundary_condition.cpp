// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//


#include "displacement_shifted_boundary_condition.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry)
    : Condition(
        NewId,
        pGeometry)
{
}

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(
        NewId,
        pGeometry,
        pProperties)
{
}

DisplacementShiftedBoundaryCondition::DisplacementShiftedBoundaryCondition(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    :Condition(
        NewId,
        ThisNodes)
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
    const auto &r_geometry = this->GetGeometry();
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

void DisplacementShiftedBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get problem dimension
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Check (and resize) LHS and RHS matrix
    const auto &r_geometry = this->GetGeometry();
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

    // Interpolate conductivity and get unknown values
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (std::size_t d = 0; d < n_dim; ++d) {
            unknown_values(i_node*n_dim + d) = r_disp[d];
        }
    }

    // Get Dirichlet BC imposition data
    const double h = GetValue(ELEMENT_H);
    const auto& r_bc_val = GetValue(DISPLACEMENT);
    const double gamma = rCurrentProcessInfo[PENALTY_DIRICHLET];

    // Calculate the Nitsche BC imposition contribution
    double aux_1;
    double aux_2;
    const double aux_weight = w * gamma; //TODO: Temporary penalty constant (this is not unit-consisten)
    // const double aux_weight_stab = w * k;
    DenseVector<double> i_node_grad(rCurrentProcessInfo[DOMAIN_SIZE]);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        //TODO: So far we are only adding the penalty part of the Nitsche
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

    KRATOS_CATCH("")
}

void DisplacementShiftedBoundaryCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_nodes,n_nodes);

    KRATOS_CATCH("")
}

void DisplacementShiftedBoundaryCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes,false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(n_nodes);

    KRATOS_CATCH("")
}

int DisplacementShiftedBoundaryCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("")
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
