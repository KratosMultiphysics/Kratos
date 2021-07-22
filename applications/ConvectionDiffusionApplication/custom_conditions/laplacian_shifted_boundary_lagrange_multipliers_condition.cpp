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


#include "laplacian_shifted_boundary_lagrange_multipliers_condition.h"
#include "utilities/integration_utilities.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

LaplacianShiftedBoundaryLagrangeMultipliersCondition::LaplacianShiftedBoundaryLagrangeMultipliersCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry)
    : Condition(
        NewId,
        pGeometry)
{
}

LaplacianShiftedBoundaryLagrangeMultipliersCondition::LaplacianShiftedBoundaryLagrangeMultipliersCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(
        NewId,
        pGeometry,
        pProperties)
{
}

LaplacianShiftedBoundaryLagrangeMultipliersCondition::LaplacianShiftedBoundaryLagrangeMultipliersCondition(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    :Condition(
        NewId,
        ThisNodes)
{
}

LaplacianShiftedBoundaryLagrangeMultipliersCondition::~LaplacianShiftedBoundaryLagrangeMultipliersCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer LaplacianShiftedBoundaryLagrangeMultipliersCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryLagrangeMultipliersCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer LaplacianShiftedBoundaryLagrangeMultipliersCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryLagrangeMultipliersCondition>(NewId, pGeom, pProperties);
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::EquationIdVector(
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
    if (rResult.size() != 2*n_nodes) {
        rResult.resize(2*n_nodes, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    for (std::size_t i = 0; i < n_nodes; ++i){
        rResult[2*i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
        rResult[2*i+1] = r_geometry[i].GetDof(SCALAR_LAGRANGE_MULTIPLIER).EquationId();
    }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::GetDofList(
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
    if (rConditionalDofList.size() != 2*n_nodes){
        rConditionalDofList.resize(2*n_nodes);
    }

    // Fill the DOFs vector from the condition nodes
    for (std::size_t i = 0; i < n_nodes; ++i){
        rConditionalDofList[2*i] = r_geometry[i].pGetDof(r_unknown_var);
        rConditionalDofList[2*i+1] = r_geometry[i].pGetDof(SCALAR_LAGRANGE_MULTIPLIER);
    }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::CalculateLocalSystem(
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
    const std::size_t local_size = 2*n_nodes;
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size,false);
    }
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(local_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size,local_size);

    // Get meshless geometry data
    const double w = GetValue(INTEGRATION_WEIGHT);
    //FIXME: Find variables for these
    const auto& r_N = GetValue(BDF_COEFFICIENTS);
    const auto& r_DN_DX = GetValue(LOCAL_AXES_MATRIX);
    array_1d<double,3>& r_normal = GetValue(NORMAL);
    r_normal /= norm_2(r_normal);

    // Interpolate conductivity and get unknown values
    double k = 0.0;
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        k += r_N[i_node] * r_geometry[i_node].FastGetSolutionStepValue(r_diffusivity_var);
        unknown_values(2*i_node) = r_geometry[i_node].FastGetSolutionStepValue(r_unknown_var);
        unknown_values(2*i_node+1) = r_geometry[i_node].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
    }

    // Check user-defined BCs
    const bool is_neumann = Has(r_flux_var);
    const bool is_dirichlet = Has(r_unknown_var);
    KRATOS_ERROR_IF(is_dirichlet && is_neumann) << "Inconsistent BCs. Condition " << Id() << " has both unknown and flux values to be imposed." << std::endl;

    // Calculate boundary integration point contribution
    if (is_dirichlet) {
        // Calculate the Lagrange multipliers contribution
        double aux_1;
        double aux_2;
        double i_node_grad_proj;
        double j_node_grad_proj;
        const double h = GetValue(ELEMENT_H);
        const double& r_bc_val = GetValue(r_unknown_var);
        const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
        DenseVector<double> i_node_grad(n_dim);
        DenseVector<double> j_node_grad(n_dim);
        const double stab_1 = 1.0e-6 * std::pow(h,2);
        const double stab_2 = 1.0e-6 * std::pow(h,2);
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            aux_1 = w * r_N(i_node);
            i_node_grad_proj = 0.0;
            i_node_grad = row(r_DN_DX, i_node);
            for (std::size_t d = 0; d < n_dim; ++d) {
                i_node_grad_proj += i_node_grad(d) * r_normal(d);
            }
            i_node_grad_proj *= w;
            for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
                aux_2 = aux_1 * r_N[j_node];
                j_node_grad_proj = 0.0;
                j_node_grad = row(r_DN_DX, j_node);
                for (std::size_t d = 0; d < n_dim; ++d) {
                    j_node_grad_proj += j_node_grad(d) * r_normal(d);
                }
                // u - lambda contribution
                rLeftHandSideMatrix(2*i_node, 2*j_node+1) -= aux_2;
                rRightHandSideVector(2*i_node) += aux_2 * unknown_values(2*j_node+1);
                // lambda - u contribution
                rLeftHandSideMatrix(2*i_node+1, 2*j_node) -= aux_2;
                rRightHandSideVector(2*i_node+1) += aux_2 * unknown_values(2*j_node);

                // stabilization - term 1
                rLeftHandSideMatrix(2*i_node, 2*j_node) -= stab_1 * i_node_grad_proj * j_node_grad_proj;
                rRightHandSideVector(2*i_node) += stab_1 * i_node_grad_proj * j_node_grad_proj * unknown_values(2*j_node);

                rLeftHandSideMatrix(2*i_node, 2*j_node + 1) -= stab_1 * i_node_grad_proj * r_N[j_node];
                rRightHandSideVector(2*i_node) += stab_1 * i_node_grad_proj * r_N[j_node] * unknown_values(2*j_node + 1);

                rLeftHandSideMatrix(2*i_node + 1, 2*j_node) -= stab_1 * aux_1 * j_node_grad_proj;
                rRightHandSideVector(2*i_node + 1) += stab_1 * aux_1 * j_node_grad_proj * unknown_values(2*j_node);

                rLeftHandSideMatrix(2*i_node + 1, 2*j_node + 1) -= stab_1 * aux_1 * r_N[j_node];
                rRightHandSideVector(2*i_node + 1) += stab_1 * aux_1 * r_N[j_node] * unknown_values(2*j_node + 1);

                // stabilization - term 2
                rLeftHandSideMatrix(2*i_node, 2*j_node) += stab_2 * aux_2;
                rRightHandSideVector(2*i_node) -= stab_2 * aux_2 * unknown_values(2*j_node);
            }
            // lambda - u contribution (BC value)
            rRightHandSideVector(2*i_node+1) -= aux_1 * r_bc_val;

            // stabilization - term 2 (BC value)
            rRightHandSideVector(2*i_node) += stab_2 * aux_1 * r_bc_val;
        }
    } else if (is_neumann) {
        // Get Neumann BC imposition data
        const double& r_bc_flux_val = GetValue(r_flux_var);

        // // // Add the Nitsche Neumann weak BC imposition (RICCARDO)
        // // for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        // //     rRightHandSideVector(i_node) += w * r_N[i_node] * r_bc_flux_val;
        // // }

        // // Add the Nitsche Neumann BC imposition (RUBEN)
        // double aux_1;
        // double aux_2;
        // const std::size_t dim = rCurrentProcessInfo[DOMAIN_SIZE];
        // DenseVector<double> j_node_grad(dim);
        // for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        //     aux_1 = w * r_N[i_node];
        //     for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
        //         aux_2 = 0.0;
        //         j_node_grad = row(r_DN_DX, j_node);
        //         for (std::size_t d = 0; d < dim; ++d) {
        //             aux_2 += j_node_grad(d) * r_normal(d);
        //         }
        //         aux_2 *= k;
        //         rLeftHandSideMatrix(i_node, j_node) += aux_1 * aux_2;
        //         rRightHandSideVector(i_node) -= aux_1 * aux_2 * unknown_values(j_node);
        //     }
        //     rRightHandSideVector(i_node) += aux_1 * r_bc_flux_val;
        // }

        // // // Add the Nitsche Neumann BC imposition
        // // double aux_1;
        // // double aux_2;
        // // const double h = GetValue(ELEMENT_H);
        // // const double gamma = rCurrentProcessInfo[INITIAL_PENALTY];
        // // DenseVector<double> i_node_grad(rCurrentProcessInfo[DOMAIN_SIZE]);
        // // DenseVector<double> j_node_grad(rCurrentProcessInfo[DOMAIN_SIZE]);
        // // for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        // //     aux_1 = r_N[i_node];
        // //     i_node_grad = row(r_DN_DX, i_node);
        // //     for (std::size_t d = 0; d < i_node_grad.size(); ++d) {
        // //         // aux_1 -= gamma * h * i_node_grad(d) * r_normal(d);
        // //         aux_1 += gamma * h * i_node_grad(d) * r_normal(d);
        // //     }
        // //     aux_1 *= w;
        // //     for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
        // //         aux_2 = 0.0;
        // //         j_node_grad = row(r_DN_DX, j_node);
        // //         for (std::size_t d = 0; d < j_node_grad.size(); ++d) {
        // //             aux_2 += j_node_grad(d) * r_normal(d);
        // //         }
        // //         aux_2 *= k;
        // //         rLeftHandSideMatrix(i_node, j_node) += aux_1 * aux_2;
        // //         rRightHandSideVector(i_node) -= aux_1 * aux_2 * unknown_values(j_node);
        // //     }
        // //     rRightHandSideVector(i_node) += aux_1 * r_bc_flux_val;
        // // }
    } else {
        KRATOS_ERROR << "No BCs are specified in condition " << Id() << "." << std::endl;
    }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const std::size_t local_size = 2*n_nodes;
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size,local_size);

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const std::size_t local_size = 2*n_nodes;
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size,false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(local_size);

    KRATOS_CATCH("")
}

int LaplacianShiftedBoundaryLagrangeMultipliersCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("")
}

// Input and Output ///////////////////////////////////////////////////////////

std::string LaplacianShiftedBoundaryLagrangeMultipliersCondition::Info() const
{
    std::stringstream buffer;
    buffer << "LaplacianShiftedBoundaryLagrangeMultipliersCondition #" << Id();
    return buffer.str();
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "LaplacianShiftedBoundaryLagrangeMultipliersCondition #" << Id();
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "LaplacianShiftedBoundaryLagrangeMultipliersCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

LaplacianShiftedBoundaryLagrangeMultipliersCondition::LaplacianShiftedBoundaryLagrangeMultipliersCondition():
    Condition()
{
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void LaplacianShiftedBoundaryLagrangeMultipliersCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
