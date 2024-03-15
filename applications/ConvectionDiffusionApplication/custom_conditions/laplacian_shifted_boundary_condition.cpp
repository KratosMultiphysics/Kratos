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


#include "laplacian_shifted_boundary_condition.h"
#include "utilities/integration_utilities.h"
#include "includes/convection_diffusion_settings.h"
#include "convection_diffusion_application_variables.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

LaplacianShiftedBoundaryCondition::LaplacianShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

LaplacianShiftedBoundaryCondition::LaplacianShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

LaplacianShiftedBoundaryCondition::LaplacianShiftedBoundaryCondition(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    :Condition(NewId, ThisNodes)
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer LaplacianShiftedBoundaryCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer LaplacianShiftedBoundaryCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryCondition>(NewId, pGeom, pProperties);
}

void LaplacianShiftedBoundaryCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the equation ids. vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    if (rResult.size() != n_nodes) {
        rResult.resize(n_nodes, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    for (std::size_t i = 0; i < n_nodes; ++i){
        rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the DOFs vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    if (rConditionalDofList.size() != n_nodes){
        rConditionalDofList.resize(n_nodes);
    }

    // Fill the DOFs vector from the condition nodes
    for (std::size_t i = 0; i < n_nodes; ++i){
        rConditionalDofList[i] = r_geometry[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const auto& r_flux_var = r_conv_diff_settings.GetSurfaceSourceVariable();
    const auto& r_diffusivity_var = r_conv_diff_settings.GetDiffusionVariable();

    // Check (and resize) LHS and RHS matrix
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes,false);
    }
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(n_nodes);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_nodes,n_nodes);

    // Get meshless geometry data
    const double w = GetValue(INTEGRATION_WEIGHT);
    const auto& r_N = GetValue(SHAPE_FUNCTIONS_VECTOR);
    const auto& r_DN_DX = GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);
    array_1d<double,3> normal = GetValue(NORMAL);
    normal /= norm_2(normal);

    // Interpolate conductivity and get unknown values
    double k = 0.0;
    Vector unknown_values(n_nodes);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        k += r_N[i_node] * r_geometry[i_node].FastGetSolutionStepValue(r_diffusivity_var);
        unknown_values(i_node) = r_geometry[i_node].FastGetSolutionStepValue(r_unknown_var);
    }

    // Check user-defined BCs
    const bool is_neumann = Has(r_flux_var);
    const bool is_dirichlet = Has(r_unknown_var);
    KRATOS_ERROR_IF(is_dirichlet && is_neumann) << "Inconsistent BCs. Condition " << Id() << " has both unknown and flux values to be imposed." << std::endl;

    // Calculate boundary integration point contribution
    if (is_dirichlet) {
        // Get Dirichlet BC imposition data
        const double h = GetValue(ELEMENT_H);
        const double& r_bc_val = GetValue(r_unknown_var);
        const double gamma = rCurrentProcessInfo[PENALTY_COEFFICIENT];

        // Calculate the Nitsche BC imposition contribution
        double aux_1;
        double aux_stab;
        const double aux_weight = w * k * gamma / h;
        const double aux_weight_stab = w * k;
        DenseVector<double> i_node_grad(rCurrentProcessInfo[DOMAIN_SIZE]);
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            aux_1 = aux_weight * r_N(i_node);
            i_node_grad = row(r_DN_DX, i_node);
            aux_stab = 0.0;
            for (std::size_t d = 0; d < i_node_grad.size(); ++d) {
                aux_stab += i_node_grad(d) * normal(d);
            }
            aux_stab *= aux_weight_stab;
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                rLeftHandSideMatrix(i_node, j_node) += (aux_1 - aux_stab) * r_N[j_node];
                rRightHandSideVector(i_node) -= (aux_1 - aux_stab) * r_N[j_node] * unknown_values(j_node);
            }
            rRightHandSideVector(i_node) += aux_1 * r_bc_val;
            rRightHandSideVector(i_node) -= aux_stab * r_bc_val;
        }
    } else if (is_neumann) {
        // Get Neumann BC imposition data
        const double& r_bc_flux_val = GetValue(r_flux_var);

        // Add the Neumann BC imposition (RUBEN)
        double aux_1;
        double aux_2;
        const std::size_t dim = rCurrentProcessInfo[DOMAIN_SIZE];
        DenseVector<double> j_node_grad(dim);
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            aux_1 = w * r_N(i_node);
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                noalias(j_node_grad) = row(r_DN_DX, j_node);
                for (std::size_t d = 0; d < dim; ++d) {
                    aux_2 = aux_1 * k * j_node_grad(d) * normal(d);
                    rRightHandSideVector(i_node) -= aux_2 * unknown_values(j_node);
                    rLeftHandSideMatrix(i_node, j_node) += aux_2;
                }
            }
            rRightHandSideVector(i_node) += aux_1 * r_bc_flux_val;
        }
    } else {
        KRATOS_ERROR << "No BCs are specified in condition " << Id() << "." << std::endl;
    }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::CalculateLeftHandSide(
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

void LaplacianShiftedBoundaryCondition::CalculateRightHandSide(
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

int LaplacianShiftedBoundaryCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("")
}

// Input and Output ///////////////////////////////////////////////////////////

std::string LaplacianShiftedBoundaryCondition::Info() const
{
    std::stringstream buffer;
    buffer << "LaplacianShiftedBoundaryCondition #" << Id();
    return buffer.str();
}

void LaplacianShiftedBoundaryCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "LaplacianShiftedBoundaryCondition #" << Id();
}

void LaplacianShiftedBoundaryCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "LaplacianShiftedBoundaryCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

void LaplacianShiftedBoundaryCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void LaplacianShiftedBoundaryCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
