//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "shifted_boundary_wall_condition.h"


namespace Kratos
{

// Public Operations //////////////////////////////////////////////////////////

void ShiftedBoundaryWallCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    /*auto& r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const auto& r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const auto& r_flux_var = r_conv_diff_settings.GetSurfaceSourceVariable();
    const auto& r_diffusivity_var = r_conv_diff_settings.GetDiffusionVariable();*/

    // Get problem dimensions
    const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Check (and resize) LHS and RHS matrix
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t local_size = n_nodes * (n_dim+1);
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size, false);
    }
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(local_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

    /*
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
    */

    KRATOS_CATCH("")
}

void ShiftedBoundaryWallCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType rRightHandSideVector;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

void ShiftedBoundaryWallCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType rLeftHandSideMatrix;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

void ShiftedBoundaryWallCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the equation ids. vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
    const  std::size_t local_size = (n_dim+1) * n_nodes;
    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    // TODO check Navier Stokes!!!
    // Fill the equation ids. vector from the condition DOFs
    const  std::size_t u_x_pos = r_geometry[0].GetDofPosition(VELOCITY_X);
    if (n_dim == 2) {
        for (std::size_t i = 0; i < n_nodes; ++i){
            rResult[i*n_dim] = r_geometry[i].GetDof(VELOCITY_X, u_x_pos).EquationId();
            rResult[i*n_dim +1] = r_geometry[i].GetDof(VELOCITY_Y, u_x_pos+1).EquationId();
            rResult[i*n_dim +2] = r_geometry[i].GetDof(PRESSURE, u_x_pos+2).EquationId();
        }
    } else {
        for (std::size_t i = 0; i < n_nodes; ++i){
            rResult[i*n_dim] = r_geometry[i].GetDof(VELOCITY_X, u_x_pos).EquationId();
            rResult[i*n_dim +1] = r_geometry[i].GetDof(VELOCITY_Y, u_x_pos+1).EquationId();
            rResult[i*n_dim +2] = r_geometry[i].GetDof(VELOCITY_Z, u_x_pos+2).EquationId();
            rResult[i*n_dim +3] = r_geometry[i].GetDof(PRESSURE, u_x_pos+3).EquationId();
        }
    }

    KRATOS_CATCH("")
}

void ShiftedBoundaryWallCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the DOFs vector
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
    const  std::size_t local_size = (n_dim+1) * n_nodes;
    if (rConditionalDofList.size() != local_size){
        rConditionalDofList.resize(local_size);
    }

    // TODO check Navier Stokes!!!
    // Fill the DOFs vector from the condition nodes
    if (n_dim == 2) {
        for (std::size_t i = 0; i < n_nodes; ++i) {
            rConditionalDofList[i*n_dim] = r_geometry[i].pGetDof(VELOCITY_X);
            rConditionalDofList[i*n_dim +1] = r_geometry[i].pGetDof(VELOCITY_Y);
            rConditionalDofList[i*n_dim +2] = r_geometry[i].pGetDof(PRESSURE);
        }
    } else {
        for (std::size_t i = 0; i < n_nodes; ++i) {
            rConditionalDofList[i*n_dim] = r_geometry[i].pGetDof(VELOCITY_X);
            rConditionalDofList[i*n_dim +1] = r_geometry[i].pGetDof(VELOCITY_Y);
            rConditionalDofList[i*n_dim +2] = r_geometry[i].pGetDof(VELOCITY_Z);
            rConditionalDofList[i*n_dim +3] = r_geometry[i].pGetDof(PRESSURE);
        }
    }

    KRATOS_CATCH("")
}

int ShiftedBoundaryWallCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // TODO check certain variables as in NavierStokesWallCondition? like viscosity?

    //Note that the base check is intentionally not called to avoid calling the base geometry check
    /*int check = BaseType::Check(rCurrentProcessInfo);
    if (check != 0) {
        return check;
    } else {
        // Checks on nodes
        const auto& r_geom = BaseType::GetGeometry();
        for (const auto& r_node : r_geom) {
            if(r_node.SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false){
                KRATOS_ERROR << "missing DYNAMIC_VISCOSITY variable on solution step data for node " << r_node.Id();
            }
        }

        return check;
    }*/
    return 0;

    KRATOS_CATCH("");
}

// Protected Operations //////////////////////////////////////////////////////////

// template<unsigned int TDim, unsigned int TNumNodes, class TWallModel>
// void ShiftedBoundaryWallCondition<TDim,TNumNodes,TWallModel>::ComputeGaussPointNavierSlipRHSContribution(
//     array_1d<double,LocalSize>& rRightHandSideVector,
//     const ConditionDataStruct& rDataStruct )
// {
//     KRATOS_TRY

//     const auto& r_geom = this->GetGeometry();
//     const double zero_tol = 1.0e-12;
//     array_1d<double,3> nodal_normal;
//     BoundedMatrix<double, TNumNodes, TNumNodes> nodal_projection_matrix;
//     for (unsigned int i_node = 0; i_node < TNumNodes; i_node++){
//         // Finding the nodal tangential projection matrix (I - n x n)
//         nodal_normal = r_geom[i_node].FastGetSolutionStepValue(NORMAL);
//         const double norm = norm_2(nodal_normal);
//         if (norm > zero_tol) {
//             nodal_normal /= norm;
//         }
//         FluidElementUtilities<3>::SetTangentialProjectionMatrix(nodal_normal, nodal_projection_matrix);

//         // Finding the coefficent to relate velocity to drag
//         const double viscosity = r_geom[i_node].GetSolutionStepValue(DYNAMIC_VISCOSITY);
//         const double navier_slip_length = r_geom[i_node].GetValue(SLIP_LENGTH);
//         KRATOS_ERROR_IF_NOT( navier_slip_length > 0.0 ) << "Negative or zero slip length was defined" << std::endl;
//         const double nodal_beta = viscosity / navier_slip_length;

//         const array_1d<double, TNumNodes> N = rDataStruct.N;
//         const double wGauss = rDataStruct.wGauss;
//         Vector interpolated_velocity = ZeroVector(TNumNodes);
//         for( unsigned int comp = 0; comp < TNumNodes; comp++){
//             for (unsigned int i = 0; i < TNumNodes; i++){
//                 // necessary because VELOCITY with 3 entries even in 2D case
//                 interpolated_velocity[i] -= N[comp] * r_geom[comp].FastGetSolutionStepValue(VELOCITY)[i];
//             }
//         }
//         // Application of the nodal projection matrix
//         const array_1d<double,TNumNodes> nodal_entry_rhs = prod( nodal_projection_matrix, (wGauss * N[i_node] * nodal_beta * interpolated_velocity) );
//         for (unsigned int entry = 0; entry < TNumNodes; entry++){
//             rRightHandSideVector( i_node*(TNumNodes+1) + entry ) += nodal_entry_rhs[entry];
//         }
//     }

//     KRATOS_CATCH("")
// }


// template<unsigned int TDim, unsigned int TNumNodes, class TWallModel>
// void ShiftedBoundaryWallCondition<TDim,TNumNodes,TWallModel>::ComputeGaussPointNavierSlipLHSContribution(
//     BoundedMatrix<double,LocalSize,LocalSize>& rLeftHandSideMatrix,
//     const ConditionDataStruct& rDataStruct )
// {
//     KRATOS_TRY

//     const GeometryType& r_geom = this->GetGeometry();

//     array_1d<double, TNumNodes> N = rDataStruct.N;
//     const double wGauss = rDataStruct.wGauss;

//     for(unsigned int inode = 0; inode < TNumNodes; inode++){

//         // finding the nodal projection matrix nodal_projection_matrix = ( [I] - (na)(na) )
//         BoundedMatrix<double, TNumNodes, TNumNodes> nodal_projection_matrix;
//         array_1d<double,3> nodal_normal = r_geom[inode].FastGetSolutionStepValue(NORMAL);
//         double sum_of_squares = 0.0;
//         for (unsigned int j = 0; j < 3; j++){
//             sum_of_squares += nodal_normal[j] * nodal_normal[j];
//         }
//         nodal_normal /= sqrt(sum_of_squares);
//         FluidElementUtilities<3>::SetTangentialProjectionMatrix( nodal_normal, nodal_projection_matrix );

//         // finding the coefficent to relate velocity to drag
//         const double viscosity = r_geom[inode].GetSolutionStepValue(DYNAMIC_VISCOSITY);
//         const double navier_slip_length = r_geom[inode].GetValue(SLIP_LENGTH);
//         KRATOS_ERROR_IF_NOT( navier_slip_length > 0.0 ) << "Negative or zero slip length was defined" << std::endl;
//         const double nodal_beta = viscosity / navier_slip_length;

//         for(unsigned int jnode = 0; jnode < TNumNodes; jnode++){

//             const BoundedMatrix<double, TNumNodes, TNumNodes> nodal_lhs_contribution = wGauss * nodal_beta * N[inode] * N[jnode] * nodal_projection_matrix;

//             for( unsigned int i = 0; i < TNumNodes; i++){
//                 for( unsigned int j = 0; j < TNumNodes; j++){

//                     const unsigned int istart = inode * (TNumNodes+1);
//                     const unsigned int jstart = jnode * (TNumNodes+1);
//                     rLeftHandSideMatrix(istart + i, jstart + j) += nodal_lhs_contribution(i,j);
//                 }
//             }
//         }
//     }

//     KRATOS_CATCH("")
// }

} // namespace Kratos
