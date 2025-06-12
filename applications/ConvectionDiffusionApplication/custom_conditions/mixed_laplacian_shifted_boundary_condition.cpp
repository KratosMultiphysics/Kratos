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


#include "mixed_laplacian_shifted_boundary_condition.h"
#include "utilities/integration_utilities.h"
#include "includes/convection_diffusion_settings.h"
#include "convection_diffusion_application_variables.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

MixedLaplacianShiftedBoundaryCondition::MixedLaplacianShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

MixedLaplacianShiftedBoundaryCondition::MixedLaplacianShiftedBoundaryCondition(
    IndexType NewId,
    Geometry<Node>::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

MixedLaplacianShiftedBoundaryCondition::MixedLaplacianShiftedBoundaryCondition(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    :Condition(NewId, ThisNodes)
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer MixedLaplacianShiftedBoundaryCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedLaplacianShiftedBoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer MixedLaplacianShiftedBoundaryCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<MixedLaplacianShiftedBoundaryCondition>(NewId, pGeom, pProperties);
}

void MixedLaplacianShiftedBoundaryCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_conv_diff_settings.GetGradientVariable();

    // Resize the equation ids. vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t dim = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const std::size_t block_size = dim + 1;
    const std::size_t local_size = block_size * n_nodes;
    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    std::size_t i = 0;
    const auto& r_gradient_var_x = KratosComponents<Variable<double>>::Get(r_gradient_var.Name() + "_X");
    const auto& r_gradient_var_y = KratosComponents<Variable<double>>::Get(r_gradient_var.Name() + "_Y");
    const auto& r_gradient_var_z = KratosComponents<Variable<double>>::Get(r_gradient_var.Name() + "_Z");
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node){
        rResult[i++] = r_geometry[i_node].GetDof(r_unknown_var).EquationId();
        rResult[i++] = r_geometry[i_node].GetDof(r_gradient_var_x).EquationId();
        rResult[i++] = r_geometry[i_node].GetDof(r_gradient_var_y).EquationId();
        if (dim == 3) {
            rResult[i++] = r_geometry[i_node].GetDof(r_gradient_var_z).EquationId();
        }
    }

    KRATOS_CATCH("")
}

void MixedLaplacianShiftedBoundaryCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_conv_diff_settings.GetGradientVariable();

    // Resize the equation ids. vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t dim = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const std::size_t block_size = dim + 1;
    const std::size_t local_size = block_size * n_nodes;
    if (rConditionalDofList.size() != local_size) {
        rConditionalDofList.resize(local_size);
    }

    // Fill the equation ids. vector from the condition DOFs
    std::size_t i = 0;
    const auto& r_gradient_var_x = KratosComponents<Variable<double>>::Get(r_gradient_var.Name() + "_X");
    const auto& r_gradient_var_y = KratosComponents<Variable<double>>::Get(r_gradient_var.Name() + "_Y");
    const auto& r_gradient_var_z = KratosComponents<Variable<double>>::Get(r_gradient_var.Name() + "_Z");
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node){
        rConditionalDofList[i++] = r_geometry[i_node].pGetDof(r_unknown_var);
        rConditionalDofList[i++] = r_geometry[i_node].pGetDof(r_gradient_var_x);
        rConditionalDofList[i++] = r_geometry[i_node].pGetDof(r_gradient_var_y);
        if (dim == 3) {
            rConditionalDofList[i++] = r_geometry[i_node].pGetDof(r_gradient_var_z);
        }
    }

    KRATOS_CATCH("")
}

void MixedLaplacianShiftedBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const auto& r_gradient_var = r_conv_diff_settings.GetGradientVariable();
    const auto& r_flux_var = r_conv_diff_settings.GetSurfaceSourceVariable();
    const auto& r_diffusivity_var = r_conv_diff_settings.GetDiffusionVariable();

    // Check (and resize) LHS and RHS matrix
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t dim = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const std::size_t block_size = dim + 1;
    const std::size_t local_size = block_size * n_nodes;
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
    double k = 0.0;
    Vector unknown_values(n_nodes);
    std::vector<array_1d<double,3>> gradient_values(n_nodes);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        k += r_N[i_node] * r_geometry[i_node].FastGetSolutionStepValue(r_diffusivity_var);
        unknown_values(i_node) = r_geometry[i_node].FastGetSolutionStepValue(r_unknown_var);
        gradient_values[i_node] = r_geometry[i_node].FastGetSolutionStepValue(r_gradient_var);
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
                rLeftHandSideMatrix(i_node*block_size, j_node*block_size) += (aux_1 - aux_stab) * r_N[j_node];
                rRightHandSideVector(i_node*block_size) -= (aux_1 - aux_stab) * r_N[j_node] * unknown_values(j_node);
            }
            rRightHandSideVector(i_node*block_size) += aux_1 * r_bc_val;
            rRightHandSideVector(i_node*block_size) -= aux_stab * r_bc_val;
        }
    } else if (is_neumann) {
        // Get Neumann BC imposition data
        const double& r_bc_flux_val = GetValue(r_flux_var);

        // Add the Neumann BC imposition
        double aux_1;
        double aux_2;
        array_1d<double,3> grad_j_node;
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            aux_1 = w * r_N(i_node);
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                aux_2 = aux_1 * k * r_N(j_node);
                noalias(grad_j_node) = gradient_values[j_node];
                for (std::size_t d = 0; d < dim; ++d) {
                    rRightHandSideVector(i_node*block_size) -= aux_2 * grad_j_node(d) * normal(d);
                    rLeftHandSideMatrix(i_node*block_size, j_node*block_size + 1 + d) += aux_2 * normal(d);
                }
            }
            rRightHandSideVector(i_node*block_size) += aux_1 * r_bc_flux_val;
        }
    } else {
        KRATOS_ERROR << "No BCs are specified in condition " << Id() << "." << std::endl;
    }

    KRATOS_CATCH("")
}

void MixedLaplacianShiftedBoundaryCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t dim = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const std::size_t block_size = dim + 1;
    const std::size_t local_size = block_size * n_nodes;
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Initialize LHS matrix
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

    KRATOS_CATCH("")
}

void MixedLaplacianShiftedBoundaryCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t dim = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const std::size_t block_size = dim + 1;
    const std::size_t local_size = block_size * n_nodes;
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size, false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(local_size);

    KRATOS_CATCH("")
}

int MixedLaplacianShiftedBoundaryCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("")
}

// Input and Output ///////////////////////////////////////////////////////////

std::string MixedLaplacianShiftedBoundaryCondition::Info() const
{
    std::stringstream buffer;
    buffer << "MixedLaplacianShiftedBoundaryCondition #" << Id();
    return buffer.str();
}

void MixedLaplacianShiftedBoundaryCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "MixedLaplacianShiftedBoundaryCondition #" << Id();
}

void MixedLaplacianShiftedBoundaryCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "MixedLaplacianShiftedBoundaryCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

void MixedLaplacianShiftedBoundaryCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void MixedLaplacianShiftedBoundaryCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
