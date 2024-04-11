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
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "shifted_boundary_wall_condition.h"
#include <boost/numeric/ublas/vector.hpp>
#include <cstddef>
#include <string>


namespace Kratos
{

// Public Operations //////////////////////////////////////////////////////////

template<std::size_t TDim>
void ShiftedBoundaryWallCondition<TDim>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the equation IDs vector
    const auto &r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
    const std::size_t local_size = (n_dim+1) * n_nodes;
    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    // Get position of DOF for first velocity component and pressure
    const std::size_t ux_pos = r_geometry[0].GetDofPosition(VELOCITY_X);
    const std::size_t p_pos = r_geometry[0].GetDofPosition(PRESSURE);

    // Fill the equation IDs vector from the condition DOFs
    std::size_t i_local = 0;
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node)
    {
        rResult[i_local++] = r_geometry[i_node].GetDof(VELOCITY_X, ux_pos).EquationId();
        rResult[i_local++] = r_geometry[i_node].GetDof(VELOCITY_Y, ux_pos+1).EquationId();
        if (n_dim == 3) rResult[i_local++] = r_geometry[i_node].GetDof(VELOCITY_Z, ux_pos+2).EquationId();
        rResult[i_local++] = r_geometry[i_node].GetDof(PRESSURE, p_pos).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void ShiftedBoundaryWallCondition<TDim>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Resize the DOFs vector
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];
    const std::size_t local_size = (n_dim+1) * n_nodes;
    if (rConditionalDofList.size() != local_size){
        rConditionalDofList.resize(local_size);
    }

    // Get position of DOF for first velocity component and pressure
    const std::size_t ux_pos = r_geometry[0].GetDofPosition(VELOCITY_X);
    const std::size_t p_pos = r_geometry[0].GetDofPosition(PRESSURE);

    // Fill the DOFs vector from the condition nodes
    std::size_t i_local = 0;
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node)
    {
        rConditionalDofList[i_local++] = r_geometry[i_node].pGetDof(VELOCITY_X, ux_pos);
        rConditionalDofList[i_local++] = r_geometry[i_node].pGetDof(VELOCITY_Y, ux_pos+1);
        if (n_dim == 3) rConditionalDofList[i_local++] = r_geometry[i_node].pGetDof(VELOCITY_Z, ux_pos+2);
        rConditionalDofList[i_local++] = r_geometry[i_node].pGetDof(PRESSURE, p_pos);
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void ShiftedBoundaryWallCondition<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    /*

    // Get problem dimensions
    const std::size_t n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Check (and resize) LHS and RHS matrix
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();  // Number of cloud nodes which contribute to the values at the integration point
    const std::size_t block_size = n_dim + 1;
    const std::size_t local_size = n_nodes * block_size;    // All unknowns of every cloud node
    const std::size_t strain_size = (n_dim-1) * 3;          // see WeaklyCompressibleNavierStokes
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size, false);
    }
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(local_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

    // Get meshless geometry data (for the integration point)
    const double parent_size = GetValue(ELEMENT_H);                     // parent element size
    const double w_int = GetValue(INTEGRATION_WEIGHT);                  // integration weight for the integration point
    const auto& r_N = GetValue(SHAPE_FUNCTIONS_VECTOR);                 // shape function values for all cloud points
    const auto& r_DN_DX = GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);    // shape function spacial derivatives for all cloud points
    array_1d<double,3> normal = GetValue(NORMAL);                       // integration weight ate the integration point
    normal /= norm_2(normal);

    // Obtain the previous iteration velocity and pressure solution (and subtract the embedded nodal velocity for FM-ALE) for all cloud nodes
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_velocity= r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
        //const auto& r_embedded_velocity = r_geometry[i_node].GetValue(EMBEDDED_VELOCITY);
        for (std::size_t d = 0; d < n_dim; ++d) {
            unknown_values[i_node*block_size + d] = r_velocity[d];  // - r_embedded_velocity[d];
        }
        unknown_values[i_node*block_size + n_dim] = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE);
    }

    // Set whether the shear stress term is adjoint consistent (1.0) or adjoint inconsistent (-1.0) for Nitsche imposition
    const double adjoint_consistency = -1.0;

    // Get process data
    const double slip_length = rCurrentProcessInfo.GetValue(SLIP_LENGTH);
    const double penalty = 1.0 / rCurrentProcessInfo.GetValue(PENALTY_COEFFICIENT);
    const double delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);

    // Get material data  //TODO ??? produced by the constitutive law; dynamic viscosity from Properties?
    const double effective_viscosity = 0.0;

    // Using Nitsche imposition of Navier-Slip boundary condition (Winter, 2018)
    // NOTE: Implement the SKEW-SYMMETRIC adjoint if needed in the future?

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the Nitsche slip normal imposition penalty coefficient
    const double pen_coeff = this->ComputeSlipNormalPenaltyCoefficient(r_N, delta_time, penalty, parent_size);

    // Add slip normal penalty contribution of the integration point
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
            for (std::size_t di = 0; di < n_dim; ++di) {
                const std::size_t row = i_node * block_size + di;
                for (std::size_t dj = 0; dj < n_dim; ++dj){
                    const std::size_t col = j_node * block_size + dj;
                    double lhs_ij = pen_coeff * w_int * r_N(i_node) * normal(di) * normal(dj) * r_N(j_node);
                    rLeftHandSideMatrix(row, col) += lhs_ij;
                    rRightHandSideVector(row) -= lhs_ij * unknown_values(col);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    // Compute slip normal symmetric counterpart penalty contribution
    Matrix aux_LHS = ZeroMatrix(local_size, local_size);

    //TODO CHECK
    // Fill the pressure to Voigt notation operator normal projected matrix
    Matrix trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(local_size, n_dim);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
        for (unsigned int di = 0; di < n_dim; ++di){
            trans_pres_to_voigt_matrix_normal_op(i_node*block_size + n_dim, di) = r_N(i_node)*normal(di);
        }
    }

    // Set the shape functions auxiliary matrix
    Matrix N_mat = ZeroMatrix(n_dim, local_size);
    for (unsigned int i = 0; i < n_nodes; ++i){
        for (unsigned int comp = 0; comp < n_dim; ++comp){
            N_mat(comp, i*block_size + comp) = r_N(i);
        }
    }

    // Set the current Gauss pt. strain matrix
    Matrix B_matrix = ZeroMatrix(strain_size, local_size);
    FluidElementUtilities<n_nodes>::GetStrainMatrix(r_DN_DX, B_matrix);

    // Set the normal projection matrix (n x n)
    Matrix normal_proj_matrix;
    FluidElementUtilities<3>::SetNormalProjectionMatrix(normal, normal_proj_matrix);

    // Get the normal projection matrix in Voigt notation
    Matrix voigt_normal_proj_matrix = ZeroMatrix(n_dim, strain_size);
    FluidElementUtilities<3>::VoigtTransformForProduct(normal, voigt_normal_proj_matrix);

    // Compute some Gauss pt. auxiliary matrices
    const Matrix aux_matrix_BC = prod(trans(B_matrix), trans(rData.C));
    const Matrix aux_matrix_APnorm = prod(trans(voigt_normal_proj_matrix), normal_proj_matrix);
    const Matrix aux_matrix_BCAPnorm = prod(aux_matrix_BC, aux_matrix_APnorm);

    // Contribution coming from the shear stress operator
    noalias(aux_LHS) -= adjoint_consistency * w_int * prod(aux_matrix_BCAPnorm, N_mat);

    // Contribution coming from the pressure terms
    const Matrix aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, normal_proj_matrix);
    noalias(aux_LHS) -= w_int * prod(aux_matrix_VPnorm, N_mat);

    // Add slip normal symmetric counterpart contribution of the integration point
    noalias(rLeftHandSideMatrix) += aux_LHS;
    noalias(rRightHandSideVector) -= prod(aux_LHS, unknown_values);

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the Nitsche tangential imposition penalty coefficients
    std::pair<const double, const double> pen_coeffs = this->ComputeSlipTangentialPenaltyCoefficients(slip_length, penalty, parent_size);

    // Declare auxiliar matrices
    Matrix aux_LHS_1 = ZeroMatrix(local_size, local_size); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    Matrix aux_LHS_2 = ZeroMatrix(local_size, local_size); // Adds the contribution generated by the viscous shear force generated by the velocity

    //TODO CHECK
    // Set the shape functions auxiliar matrices
    BoundedMatrix<double, n_dim, local_size> N_mat = ZeroMatrix(n_dim, local_size);
    for (unsigned int i = 0; i < n_nodes; ++i){
        for (unsigned int comp = 0; comp < n_dim; ++comp){
            N_mat(comp, i*block_size + comp) = aux_N(i);
        }
    }
    BoundedMatrix<double, local_size, n_dim> N_mat_trans = trans(N_mat);

    // Set the tangential projection matrix (I - n x n)
    BoundedMatrix<double, n_dim, n_dim> tang_proj_matrix;
    FluidElementUtilities<3>::SetTangentialProjectionMatrix(normal, tang_proj_matrix);

    // Set the current Gauss pt. strain matrix
    BoundedMatrix<double, strain_size, local_size> B_matrix = ZeroMatrix(strain_size, local_size);
    FluidElementUtilities<local_size>::GetStrainMatrix(r_DN_DX, B_matrix);

    // Get the normal projection matrix in Voigt notation
    BoundedMatrix<double, n_dim, strain_size> voigt_normal_proj_matrix = ZeroMatrix(n_dim, strain_size);
    FluidElementUtilities<3>::VoigtTransformForProduct(normal, voigt_normal_proj_matrix);

    // Compute some Gauss pt. auxiliary matrices
    const BoundedMatrix<double, strain_size, local_size> aux_matrix_CB = prod(rData.C, B_matrix);
    const BoundedMatrix<double, strain_size, n_dim> aux_matrix_PtangA = prod(tang_proj_matrix, voigt_normal_proj_matrix);
    const BoundedMatrix<double, local_size, n_dim> aux_matrix_PtangACB = prod(aux_matrix_PtangA, aux_matrix_CB);

    // Contribution coming from the traction vector tangential component
    noalias(aux_LHS_1) += pen_coeffs.first*weight*prod(N_mat_trans, aux_matrix_PtangACB);

    // Contribution coming from the shear force generated by the velocity jump
    const BoundedMatrix<double, local_size, n_dim> aux_matrix_N_trans_tang = prod(N_mat_trans, tang_proj_matrix);
    noalias(aux_LHS_2) += pen_coeffs.second*weight*prod(aux_matrix_N_trans_tang, N_mat);

    // Add slip tangential penalty contribution of the integration point
    // NOTE for mesh movement (FM-ALE) the level set velocity contribution needs to be added here as well!
    noalias(rLeftHandSideMatrix) += aux_LHS_1;
    noalias(rLeftHandSideMatrix) += aux_LHS_2;
    noalias(rRightHandSideVector) -= prod(aux_LHS_1, unknown_values);
    noalias(rRightHandSideVector) -= prod(aux_LHS_2, unknown_values);

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the coefficients
    std::pair<const double, const double> nitsche_coeffs = this->ComputeSlipTangentialNitscheCoefficients(slip_length, penalty, parent_size);

    // Declare auxiliary matrices
    Matrix aux_LHS_1 = ZeroMatrix(local_size, local_size); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    Matrix aux_LHS_2 = ZeroMatrix(local_size, local_size); // Adds the contribution generated by the viscous shear force generated by the velocity

    //TODO CHECK
    // Set the shape functions auxiliary matrices
    BoundedMatrix<double, n_dim, local_size> N_mat = ZeroMatrix(n_dim, local_size);
    for (unsigned int i = 0; i < n_nodes; ++i){
        for (unsigned int comp = 0; comp < n_dim; ++comp){
            N_mat(comp, i*BlockSize + comp) = r_N(i);
        }
    }

    // Set the current Gauss pt. strain matrix
    BoundedMatrix<double, strain_size, local_size> B_matrix = ZeroMatrix(strain_size, local_size);
    FluidElementUtilities<3>::GetStrainMatrix(aux_DN_DX, B_matrix);

    // Set the tangential projection matrix (I - n x n)
    BoundedMatrix<double, n_dim, n_dim> tang_proj_matrix;
    FluidElementUtilities<3>::SetTangentialProjectionMatrix(normal, tang_proj_matrix);

    // Get the normal projection matrix in Voigt notation
    BoundedMatrix<double, n_dim, strain_size> voigt_normal_proj_matrix = ZeroMatrix(n_dim, strain_size);
    FluidElementUtilities<3>::VoigtTransformForProduct(normal, voigt_normal_proj_matrix);

    // Compute some Gauss pt. auxiliar matrices
    const BoundedMatrix<double, local_size, n_dim> aux_matrix_BtransAtrans = prod(trans(B_matrix), trans(voigt_normal_proj_matrix));
    const BoundedMatrix<double, local_size, n_dim> aux_matrix_BtransAtransPtan = prod(aux_matrix_BtransAtrans, tang_proj_matrix);
    const BoundedMatrix<double, strain_size, local_size> aux_matrix_CB = prod(rData.C, B_matrix);
    const BoundedMatrix<double, n_dim, local_size> aux_matrix_ACB = prod(voigt_normal_proj_matrix, aux_matrix_CB);
    const BoundedMatrix<double, local_size, local_size> aux_matrix_BtransAtransPtanACB = prod(aux_matrix_BtransAtransPtan, aux_matrix_ACB);

    // Contribution coming from the traction vector tangencial component
    noalias(aux_LHS_1) -= adjoint_consistency * nitsche_coeffs.first * w_int * aux_matrix_BtransAtransPtanACB;

    // Contribution coming from the shear force generated by the velocity jump
    noalias(aux_LHS_2) -= adjoint_consistency * nitsche_coeffs.second * w_int * prod(aux_matrix_BtransAtransPtan, N_mat);

    // Add slip tangential symmetric counterpart contribution of the integration point
    // NOTE for mesh movement (FM-ALE) the level set velocity contribution needs to be added here as well!
    noalias(rLeftHandSideMatrix) += aux_LHS_1;
    noalias(rLeftHandSideMatrix) += aux_LHS_2;
    noalias(rRightHandSideVector) -= prod(aux_LHS_1, unknown_values);
    noalias(rRightHandSideVector) -= prod(aux_LHS_2, unknown_values);

    */

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void ShiftedBoundaryWallCondition<TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType rRightHandSideVector;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template<std::size_t TDim>
void ShiftedBoundaryWallCondition<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType rLeftHandSideMatrix;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template<std::size_t TDim>
int ShiftedBoundaryWallCondition<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    //TODO check certain variables as in NavierStokesWallCondition? like viscosity?
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

    //NOTE that the base check is intentionally not called to avoid calling the base geometry check
    return 0;

    KRATOS_CATCH("");
}

// Protected Operations //////////////////////////////////////////////////////////
template<std::size_t TDim>
double ShiftedBoundaryWallCondition<TDim>::ComputeSlipNormalPenaltyCoefficient(
    const Vector& rN,
    const double DeltaTime,
    const double Penalty,
    const double ParentSize) const
{
    // Get the velocity and density for the integration point
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    double int_pt_rho = rN(0) * r_geometry[0].FastGetSolutionStepValue(DENSITY);
    array_1d<double,TDim> int_pt_v = rN(0) * r_geometry[0].FastGetSolutionStepValue(VELOCITY);
    for (unsigned int i_node = 1;  i_node < n_nodes; ++i_node) {
        int_pt_rho += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(DENSITY);
        noalias(int_pt_v) += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
    }
    const double int_pt_v_norm = norm_2(int_pt_v);

    const double eff_mu = 0.0;  //TODO

    // Compute the Nitsche coefficient (including the Winter stabilization term)
    const double coeff = (eff_mu + eff_mu + int_pt_rho*int_pt_v_norm*ParentSize + int_pt_rho*ParentSize*ParentSize/DeltaTime) / (ParentSize*Penalty);

    return coeff;
}

template<std::size_t TDim>
std::pair<const double, const double> ShiftedBoundaryWallCondition<TDim>::ComputeSlipTangentialPenaltyCoefficients(
    const double SlipLength,
    const double Penalty,
    const double ParentSize) const
{
    const double eff_mu = 0.0;  //TODO

    const double coeff_1 = SlipLength / (SlipLength + Penalty*ParentSize);
    const double coeff_2 = eff_mu / (SlipLength + Penalty*ParentSize);

    std::pair<const double, const double> coeffs(coeff_1, coeff_2);
    return coeffs;
}

template<std::size_t TDim>
std::pair<const double, const double> ShiftedBoundaryWallCondition<TDim>::ComputeSlipTangentialNitscheCoefficients(
    const double SlipLength,
    const double Penalty,
    const double ParentSize) const
{
    const double eff_mu = 0.0;  //TODO

    const double coeff_1 = SlipLength*Penalty*ParentSize / (SlipLength + Penalty*ParentSize);
    const double coeff_2 = eff_mu*Penalty*ParentSize / (SlipLength + Penalty*ParentSize);

    std::pair<const double, const double> coeffs(coeff_1, coeff_2);
    return coeffs;
}

template class ShiftedBoundaryWallCondition<2>;
template class ShiftedBoundaryWallCondition<3>;

} // namespace Kratos
