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
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "input_output/logger.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "shifted_boundary_wall_condition.h"
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
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
    const std::size_t local_size = BlockSize * n_nodes;
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
        if (TDim == 3) rResult[i_local++] = r_geometry[i_node].GetDof(VELOCITY_Z, ux_pos+2).EquationId();
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
    const std::size_t local_size = BlockSize * n_nodes;
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
        if (TDim == 3) rConditionalDofList[i_local++] = r_geometry[i_node].pGetDof(VELOCITY_Z, ux_pos+2);
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

    // Check (and resize) LHS and RHS matrix
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();  // Number of cloud nodes which contribute to the values at the integration point
    const std::size_t local_size = n_nodes * BlockSize;     // All unknowns of every cloud node
    if (rRightHandSideVector.size() != local_size) {
        rRightHandSideVector.resize(local_size, false);
    }
    if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    // Set LHS and RHS to zero
    noalias(rRightHandSideVector) = ZeroVector(local_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

    AddNitscheImposition(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

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
    // SLIP_LENGTH and PENALTY_COEFFICIENT and EMBEDDED_CHARACT_LENGTH of ProcessInfo
    // constitutive law
    // EMBEDDED_VELOCITY as variable at nodes?
    // check whether condition has ELEMENT_H, INTEGRATION_WEIGHT, SHAPE_FUNCTIONS_VECTOR, SHAPE_FUNCTIONS_GRADIENT_MATRIX and NORMAL

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
void ShiftedBoundaryWallCondition<TDim>::AddNitscheImposition(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = this->GetGeometry();
    const std::size_t local_size = rLHS.size1();
    const std::size_t n_nodes = local_size / BlockSize;

    // Get meshless geometry data (for the integration point)
    const double parent_size = this->GetValue(ELEMENT_H);               // parent element size
    const double weight = GetValue(INTEGRATION_WEIGHT);                 // integration weight for the integration point
    const auto& r_N = GetValue(SHAPE_FUNCTIONS_VECTOR);                 // shape function values for all cloud points (VectorType)
    const auto& r_DN_DX = GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);    // shape function spacial derivatives for all cloud points (MatrixType)
    array_1d<double,3> normal = GetValue(NORMAL);                       // boundary normal at the integration point
    normal /= norm_2(normal);

    // Set whether the shear stress term is adjoint consistent (1.0) or adjoint inconsistent (-1.0) for Nitsche imposition
    // NOTE that adjoint inconsistent formulation enjoys improved inf-sup stability for any value of the penalty constant (gamma) according to Winter et al. (2018)
    const double adjoint_consistency = -1.0;

    // Get process data
    const double slip_length = rCurrentProcessInfo.GetValue(SLIP_LENGTH);
    const double gamma_penalty = 1.0 / rCurrentProcessInfo.GetValue(PENALTY_COEFFICIENT);
    const double gamma_penalty_shear = 1.0 / rCurrentProcessInfo.GetValue(PENALTY_COEFFICIENT_TANGENTIAL);
    const double delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const double charact_length = rCurrentProcessInfo.GetValue(EMBEDDED_CHARACT_LENGTH);

    // Obtain the previous iteration velocity and pressure solution (and subtract the embedded nodal velocity EMBEDDED_VELOCITY for FM-ALE) for all cloud nodes
    Vector unknown_values(local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_velocity= r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
        //const auto& r_embedded_velocity = r_geometry[i_node].GetValue(EMBEDDED_VELOCITY);
        for (std::size_t d = 0; d < TDim; ++d) {
            unknown_values[i_node*BlockSize + d] = r_velocity[d];  // - r_embedded_velocity[d];  TODO
        }
        unknown_values[i_node*BlockSize + TDim] = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE);
    }

    // Get material data, calculated by the constitutive law
    double effective_viscosity = 0.0;

    /////////////////////////////////////////////////////////////////////////////////
    // Get B matrix, C matrix, Voigt normal projection matrix and some auxilary matrices

    // Set the current integration point's strain matrix (SF derivatives at all cloud points)
    Matrix B_matrix = ZeroMatrix(VoigtSize, local_size);
    CalculateStrainMatrix(r_DN_DX, n_nodes, B_matrix);  // see FluidElementUtilities

    // Get C Matrix from constitutive law, see FluidElementData and FluidElement::CalculateMaterialResponse
    const auto p_constitutive_law =  this->GetProperties()[CONSTITUTIVE_LAW];
    ConstitutiveLaw::Parameters constitutive_law_values(r_geometry, this->GetProperties(), rCurrentProcessInfo);
    Vector strain_rate = prod(B_matrix, unknown_values);
    Vector shear_stress(VoigtSize);
    Matrix C_matrix = ZeroMatrix(VoigtSize, VoigtSize);

    // Set constitutive law flags
    Flags& r_cl_options = constitutive_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);  // not being used (?)
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set constitutive law values
    constitutive_law_values.SetShapeFunctionsValues(r_N);
    constitutive_law_values.SetShapeFunctionsDerivatives(r_DN_DX);
    constitutive_law_values.SetStrainVector(strain_rate);
    constitutive_law_values.SetStressVector(shear_stress);    //this is an ouput parameter
    constitutive_law_values.SetConstitutiveMatrix(C_matrix);  //this is an ouput parameter

    // Calculate material response and effective viscosity
    //NOTE that here we assume that only one constitutive law is employed for all of the integration points.
    //this is ok under the hypothesis that no history dependent behaviour is employed
    p_constitutive_law->CalculateMaterialResponseCauchy(constitutive_law_values);
    p_constitutive_law->CalculateValue(constitutive_law_values, EFFECTIVE_VISCOSITY, effective_viscosity);

    // Auxiliary matrix C*B
    const Matrix aux_matrix_CB = prod(C_matrix, B_matrix);

    // Set the shape functions auxiliary matrix
    Matrix N_matrix = ZeroMatrix(TDim, local_size);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node){
        for (std::size_t d = 0; d < TDim; ++d){
            N_matrix(d, i_node*BlockSize + d) = r_N(i_node);
        }
    }
    const Matrix N_matrix_trans = trans(N_matrix);

    // Get the normal projection matrix in Voigt notation (A_matrix)
    BoundedMatrix<double, TDim, VoigtSize> voigt_normal_proj_matrix;
    FluidElementUtilities<3>::VoigtTransformForProduct(normal, voigt_normal_proj_matrix);

    // Auxilary matrix for adding contributions
    MatrixType aux_LHS = ZeroMatrix(local_size, local_size);

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the Nitsche normal imposition penalty coefficient
    const double pen_coeff_normal = this->ComputeSlipNormalPenaltyCoefficient(r_N, delta_time, gamma_penalty, parent_size, effective_viscosity);

    // Add normal velocity penalty contribution of the integration point
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
        for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
            for (std::size_t di = 0; di < TDim; ++di) {
                const std::size_t row = i_node * BlockSize + di;
                for (std::size_t dj = 0; dj < TDim; ++dj){
                    const std::size_t col = j_node * BlockSize + dj;
                    aux_LHS(row, col) += pen_coeff_normal * weight * r_N(i_node) * normal(di) * normal(dj) * r_N(j_node);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the Nitsche normal symmetric counterpart stabilization

    // Fill the pressure to Voigt notation operator normal projected matrix
    Matrix trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(local_size, TDim);
    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node){
        for (std::size_t d = 0; d < TDim; ++d){
            trans_pres_to_voigt_matrix_normal_op(i_node*BlockSize + TDim, d) = r_N(i_node)*normal(d);
        }
    }

    // Set the normal projection matrix (n x n)
    BoundedMatrix<double, TDim, TDim> normal_proj_matrix;
    FluidElementUtilities<3>::SetNormalProjectionMatrix(normal, normal_proj_matrix);

    // Compute some integration point auxiliary matrices
    const Matrix aux_matrix_APnorm = prod(trans(voigt_normal_proj_matrix), normal_proj_matrix);
    const Matrix aux_matrix_BCAPnorm = prod(trans(aux_matrix_CB), aux_matrix_APnorm);

    // Contribution coming from the shear stress operator (traction vector normal component)
    noalias(aux_LHS) -= adjoint_consistency * weight * prod(aux_matrix_BCAPnorm, N_matrix);

    // Contribution coming from the pressure terms
    const Matrix aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, normal_proj_matrix);
    noalias(aux_LHS) -= weight * prod(aux_matrix_VPnorm, N_matrix);

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the Nitsche slip tangential penalty coefficients
    std::pair<const double, const double> pen_coeffs_tang = this->ComputeSlipTangentialPenaltyCoefficients(slip_length, gamma_penalty, gamma_penalty_shear, parent_size, effective_viscosity);

    // Set the tangential projection matrix (I - n x n)
    BoundedMatrix<double, TDim, TDim> tang_proj_matrix;
    FluidElementUtilities<3>::SetTangentialProjectionMatrix(normal, tang_proj_matrix);

    // Compute some integration point auxiliary matrices
    const Matrix aux_matrix_N_trans_tang = prod(N_matrix_trans, tang_proj_matrix);
    const Matrix aux_matrix_ACB = prod(voigt_normal_proj_matrix, aux_matrix_CB);

    // Contribution coming from the traction vector tangential component (no contribution for no slip)
    noalias(aux_LHS) += pen_coeffs_tang.first  * weight * prod(aux_matrix_N_trans_tang, aux_matrix_ACB);

    // Contribution coming from the tangential velocity (towards zero for slip)
    noalias(aux_LHS) += pen_coeffs_tang.second * weight * prod(aux_matrix_N_trans_tang, N_matrix);

    /////////////////////////////////////////////////////////////////////////////////
    // Compute the Nitsche slip tangential symmetric counterpart stabilization
    std::pair<const double, const double> nitsche_coeffs_tang = this->ComputeSlipTangentialNitscheCoefficients(slip_length, gamma_penalty_shear, charact_length, effective_viscosity);

    // Compute some integration point auxiliary matrices
    const Matrix aux_matrix_BtransAtrans = prod(trans(B_matrix), trans(voigt_normal_proj_matrix));
    const Matrix aux_matrix_BtransAtransPtan = prod(aux_matrix_BtransAtrans, tang_proj_matrix);

    // Contribution coming from the traction vector tangential component (no contribution for no slip)
    noalias(aux_LHS) -= adjoint_consistency * nitsche_coeffs_tang.first  * weight * 2 * prod(aux_matrix_BtransAtransPtan, aux_matrix_ACB);

    // Contribution coming from the tangential velocity (towards zero for slip)
    noalias(aux_LHS) -= adjoint_consistency * nitsche_coeffs_tang.second * weight * 2 * prod(aux_matrix_BtransAtransPtan, N_matrix);

    /////////////////////////////////////////////////////////////////////////////////
    //KRATOS_WATCH(aux_LHS);

    // Add Nitsche Navier-slip contributions of the integration point to the local system (residual-based formulation with RHS=f_gamma-LHS*previous_solution)
    // NOTE for mesh movement (FM-ALE) the level set velocity contribution needs to be added here as well! TODO
    noalias(rLHS) += aux_LHS;
    noalias(rRHS) -= prod(aux_LHS, unknown_values);
}

template<std::size_t TDim>
void ShiftedBoundaryWallCondition<TDim>::CalculateStrainMatrix(
    const Matrix& rDN_DX,
    const std::size_t NumNodes,
    Matrix& rB)
{
    rB.clear();

    if(TDim == 2) {
        for (std::size_t i_node = 0; i_node < NumNodes; ++i_node ) {
            const std::size_t col = BlockSize * i_node;
            rB(0, col  ) = rDN_DX(i_node, 0);
            rB(1, col+1) = rDN_DX(i_node, 1);
            rB(2, col  ) = rDN_DX(i_node, 1);
            rB(2, col+1) = rDN_DX(i_node, 0);
        }
    } else if(TDim == 3) {
        for (std::size_t i_node = 0; i_node < NumNodes; ++i_node ) {
            const std::size_t col = BlockSize * i_node;
            rB(0, col  ) = rDN_DX(i_node, 0);
            rB(1, col+1) = rDN_DX(i_node, 1);
            rB(2, col+2) = rDN_DX(i_node, 2);
            rB(3, col  ) = rDN_DX(i_node, 1);
            rB(3, col+1) = rDN_DX(i_node, 0);
            rB(4, col+1) = rDN_DX(i_node, 2);
            rB(4, col+2) = rDN_DX(i_node, 1);
            rB(5, col  ) = rDN_DX(i_node, 2);
            rB(5, col+2) = rDN_DX(i_node, 0);
        }
    }
}

template<std::size_t TDim>
double ShiftedBoundaryWallCondition<TDim>::ComputeSlipNormalPenaltyCoefficient(
    const Vector& rN,
    const double DeltaTime,
    const double Gamma,
    const double ParentSize,
    const double EffectiveViscosity) const
{
    // Get the velocity and density for the integration point
    const auto& r_geometry = this->GetGeometry();
    const std::size_t n_nodes = r_geometry.PointsNumber();
    double int_pt_rho = rN(0) * r_geometry[0].FastGetSolutionStepValue(DENSITY);
    array_1d<double,3> int_pt_v = rN(0) * r_geometry[0].FastGetSolutionStepValue(VELOCITY);
    for (std::size_t i_node = 1;  i_node < n_nodes; ++i_node) {
        int_pt_rho += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(DENSITY);
        int_pt_v += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
    }
    const double int_pt_v_norm = norm_2(int_pt_v);

    const double stab_constant_u = EffectiveViscosity + int_pt_rho*int_pt_v_norm*ParentSize / 6.0 + int_pt_rho*ParentSize*ParentSize/DeltaTime / 12.0;

    // Compute the Nitsche coefficient (including the Schott et al. (doi: 10.1002/fld.4218) stabilization term)
    const double coeff = (EffectiveViscosity + stab_constant_u) / (Gamma * ParentSize);

    return coeff;
}

template<std::size_t TDim>
std::pair<const double, const double> ShiftedBoundaryWallCondition<TDim>::ComputeSlipTangentialPenaltyCoefficients(
    const double SlipLength,
    const double Gamma,
    const double GammaShear,
    const double ParentSize,
    const double EffectiveViscosity) const
{
    // // Get the velocity and density for the integration point
    // const auto& r_geometry = this->GetGeometry();
    // const std::size_t n_nodes = r_geometry.PointsNumber();
    // double int_pt_rho = rN(0) * r_geometry[0].FastGetSolutionStepValue(DENSITY);
    // array_1d<double,3> int_pt_v = rN(0) * r_geometry[0].FastGetSolutionStepValue(VELOCITY);
    // for (std::size_t i_node = 1;  i_node < n_nodes; ++i_node) {
    //     int_pt_rho += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(DENSITY);
    //     int_pt_v += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
    // }
    // const double int_pt_v_norm = norm_2(int_pt_v);

    // const double stab_constant_u = EffectiveViscosity + int_pt_rho*int_pt_v_norm*ParentSize / 6.0 + int_pt_rho*ParentSize*ParentSize/DeltaTime / 12.0;

    // const double penalty_coeff = 1.0 / (SlipLength + Gamma*ParentSize);

    // ShearStab (a)
    //const double coeff_1 = 0.0;
    // Winter et al. (2018): * SlipLength;
    // const double coeff_1 =  penalty_coeff * SlipLength;
    // ShearStab (b)
    // const double coeff_1 = penalty_coeff * SlipLength * 1.0/GammaShear * CharactLength / ParentSize;
    // ShearStab (c)
    // const double coeff_1 = penalty_coeff * SlipLength * 1.0/GammaShear * int_pt_rho*std::pow(int_pt_v_norm,2)*DeltaTime/EffectiveViscosity;
    // ShearStab (d)
    // const double coeff_1 = penalty_coeff * SlipLength * 1.0/GammaShear * int_pt_rho*int_pt_v_norm*ParentSize/EffectiveViscosity;
    // ShearStab (e)
    // const double coeff_1 = penalty_coeff * SlipLength * 1.0/GammaShear;
    // ShearStab (f)
    // const double coeff_1 = penalty_coeff * SlipLength * 1.0/GammaShear * int_pt_v_norm*DeltaTime/ParentSize;

    // Development version 6.2e
    const double penalty_coeff = 1.0 / ( GammaShear * (SlipLength + Gamma*ParentSize/GammaShear) );
    const double coeff_1 = penalty_coeff * SlipLength;
    const double coeff_2 = penalty_coeff * EffectiveViscosity;  // + stab_constant_u);  // Winter et al. (2018): * EffectiveViscosity;

    std::pair<const double, const double> coefficients(coeff_1, coeff_2);
    return coefficients;
}

template<std::size_t TDim>
std::pair<const double, const double> ShiftedBoundaryWallCondition<TDim>::ComputeSlipTangentialNitscheCoefficients(
    const double SlipLength,
    const double GammaShear,
    const double CharactLength,
    const double EffectiveViscosity) const
{
    // // Get the velocity and density for the integration point
    // const auto& r_geometry = this->GetGeometry();
    // const std::size_t n_nodes = r_geometry.PointsNumber();
    // double int_pt_rho = rN(0) * r_geometry[0].FastGetSolutionStepValue(DENSITY);
    // array_1d<double,3> int_pt_v = rN(0) * r_geometry[0].FastGetSolutionStepValue(VELOCITY);
    // for (std::size_t i_node = 1;  i_node < n_nodes; ++i_node) {
    //     int_pt_rho += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(DENSITY);
    //     int_pt_v += rN(i_node) * r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
    // }
    // const double int_pt_v_norm = norm_2(int_pt_v);

    // const double b_ref = 1.0;

    // const double stab_coeff = 1 / (SlipLength + 1);  // Winter et al. (2018): GammaShear * ParentSize / (SlipLength + GammaShear*ParentSize);

    // Winter et al. (2018): * SlipLength
    // const double coeff_1 = stab_coeff * SlipLength * GammaShear * ParentSize;
    // ShearStab (1) - dev2.3
    // const double coeff_1 = stab_coeff * SlipLength * ParentSize/GammaShear;
    // ShearStab (2.2)
    // const double coeff_1 = stab_coeff * SlipLength * CharactLength*b_ref/ParentSize * 1/GammaShear;
    // ShearStab (3.2)
    // const double coeff_1 = stab_coeff * SlipLength * int_pt_rho*std::pow(int_pt_v_norm,2)*ParentSize*DeltaTime/EffectiveViscosity * 1/GammaShear;
    // ShearStab (4.2)
    // const double coeff_1 = stab_coeff * SlipLength * int_pt_rho*int_pt_v_norm*std::pow(ParentSize,2)/EffectiveViscosity * 1/GammaShear;
    // ShearStab (5.2)
    // const double coeff_1 = stab_coeff * SlipLength * int_pt_v_norm*DeltaTime * 1/GammaShear;
    // ShearStab (6.2)
    // const double coeff_1 = stab_coeff * SlipLength * CharactLength * 1/GammaShear;

    // std::string grad_coeff = "Stabilization coefficient: " + std::to_string(coeff_1);
    // KRATOS_WATCH(grad_coeff);

    // Development version 6.2e
    const double stab_coeff = CharactLength / ( GammaShear * (SlipLength + CharactLength/GammaShear) );
    const double coeff_1 = stab_coeff * SlipLength;
    const double coeff_2 = stab_coeff * EffectiveViscosity;

    std::pair<const double, const double> coefficients(coeff_1, coeff_2);
    return coefficients;
}

template class ShiftedBoundaryWallCondition<2>;
template class ShiftedBoundaryWallCondition<3>;

} // namespace Kratos
