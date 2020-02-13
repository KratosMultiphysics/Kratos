//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla (based on Elisa Magliozzi previous work)
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_elements/compressible_navier_stokes_explicit.h"

namespace Kratos {

template <>
void CompressibleNavierStokesExplicit<2>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 4;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (ElementalDofList.size() != dof_size) {
        ElementalDofList.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<3>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 5;
    unsigned int dof_size = n_nodes * block_size;

    if (ElementalDofList.size() != dof_size) {
        ElementalDofList.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto &r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Z, mom_pos + 2);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TBlockSize, unsigned int TNumNodes>
int CompressibleNavierStokesExplicit<TDim, TBlockSize, TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if (ErrorCode != 0) {
        return ErrorCode;
    }

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(MOMENTUM);
    KRATOS_CHECK_VARIABLE_KEY(TOTAL_ENERGY);
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(CONDUCTIVITY);
    KRATOS_CHECK_VARIABLE_KEY(SPECIFIC_HEAT);
    KRATOS_CHECK_VARIABLE_KEY(HEAT_CAPACITY_RATIO);
    KRATOS_CHECK_VARIABLE_KEY(BODY_FORCE);
    KRATOS_CHECK_VARIABLE_KEY(EXTERNAL_PRESSURE);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY)) << "Missing DENSITY variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM)) << "Missing MOMENTUM variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE)) << "Missing BODY_FORCE variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_PRESSURE)) << "Missing EXTERNAL_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();

        // Activate as soon as we start using the explicit DOF based strategy
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing DENSITY DOF in node ", this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) || this->GetGeometry()[i].HasDofFor(MOMENTUM_Y)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
        if (TDim == 3) {
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_Z)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
        }
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing TOTAL_ENERGY DOF in node ", this->GetGeometry()[i].Id();
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TBlockSize, unsigned int TNumNodes>
void CompressibleNavierStokesExplicit<TDim, TBlockSize, TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Getting data for the given geometry
    const auto& r_geometry = GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, rData.DN_DX, rData.N, rData.volume);

    // Compute element size
    rData.h = ComputeH(rData.DN_DX);

    // Database access to all of the variables needed
    Properties &r_properties = this->GetProperties();
    rData.nu = r_properties.GetValue(KINEMATIC_VISCOSITY);
    rData.mu = r_properties.GetValue(DYNAMIC_VISCOSITY);
    rData.lambda = r_properties.GetValue(CONDUCTIVITY);
    rData.c_v = r_properties.GetValue(SPECIFIC_HEAT);
    rData.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        const array_1d<double, 3> &r_momentum = r_geometry[i].FastGetSolutionStepValue(MOMENTUM);
        const array_1d<double, 3> &r_body_force = r_geometry[i].FastGetSolutionStepValue(BODY_FORCE);

        for (unsigned int k = 0; k < TDim; ++k) {
            rData.U(i, k + 1) = r_momentum[k];
            rData.f_ext(i, k) = r_body_force[k];
        }
        rData.U(i, 0) = r_geometry[i].FastGetSolutionStepValue(DENSITY);
        rData.U(i, TDim + 1) = r_geometry[i].FastGetSolutionStepValue(TOTAL_ENERGY);
        rData.r(i) = r_geometry[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);
    }

    // Get shock capturing viscosity and heat conductivity
    CalculateShockCapturingValues(rData);
}

template <unsigned int TDim, unsigned int TBlockSize, unsigned int TNumNodes>
double CompressibleNavierStokesExplicit<TDim, TBlockSize, TNumNodes>::ComputeH(BoundedMatrix<double,TNumNodes, TDim>& rDN_DX)
{
    double h = 0.0;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        double h_inv = 0.0;
        for (unsigned int k = 0; k < TDim; ++k) {
            h_inv += rDN_DX(i,k) * rDN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    h = sqrt(h) / static_cast<double>(TNumNodes);
    return h;
}

template <>
void CompressibleNavierStokesExplicit<2>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 4;
    constexpr unsigned int matrix_size = n_nodes * block_size;

    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false); //false says not to preserve existing storage!!
    }

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes> &r = data.r;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double v_sc = data.nu_sc;
    const double k_sc = data.lambda_sc;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Solution vector values from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);

    // Hardcoded shape functions for linear triangular element
    // This is explicitly done to minimize the allocation and matrix acceses
    // The notation N_i_j means shape function for node j in Gauss pt. i
    const double one_sixt = 1.0/6.0;
    const double two_third = 2.0/3.0;
    const double N_0_0 = one_sixt;
    const double N_0_1 = one_sixt;
    const double N_0_2 = two_third;
    const double N_1_0 = one_sixt;
    const double N_1_1 = two_third;
    const double N_1_2 = one_sixt;
    const double N_2_0 = two_third;
    const double N_2_1 = one_sixt;
    const double N_2_2 = one_sixt;

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    //substitute_rhs_2D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template<>
void CompressibleNavierStokesExplicit<3>::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 5;
    constexpr unsigned int matrix_size = n_nodes * block_size;

    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false); //false says not to preserve existing storage!!
    }

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const array_1d<double, n_nodes> &r = data.r;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double mu = data.mu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double v_sc = data.nu_sc;
    const double k_sc = data.lambda_sc;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Solution vector values from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_0_4 = data.U(0, 4);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_1_4 = data.U(1, 4);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_2_4 = data.U(2, 4);
    const double &U_3_0 = data.U(3, 0);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);
    const double &U_3_4 = data.U(3, 4);

    // Hardcoded shape functions for linear tetrahedra element
    // This is explicitly done to minimize the alocation and matrix acceses
    // The notation N_i_j means shape function for node j in Gauss pt. i
    const double N_0_0 = 0.58541020;
    const double N_0_1 = 0.13819660;
    const double N_0_2 = 0.13819660;
    const double N_0_3 = 0.13819660;
    const double N_1_0 = 0.13819660;
    const double N_1_1 = 0.58541020;
    const double N_1_2 = 0.13819660;
    const double N_1_3 = 0.13819660;
    const double N_2_0 = 0.13819660;
    const double N_2_1 = 0.13819660;
    const double N_2_2 = 0.58541020;
    const double N_2_3 = 0.13819660;
    const double N_3_0 = 0.13819660;
    const double N_3_1 = 0.13819660;
    const double N_3_2 = 0.13819660;
    const double N_3_3 = 0.58541020;

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    //substitute_rhs_3D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 2;
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 4;

    // Calculate the explicit residual vector
    VectorType rhs;
    CalculateRightHandSide(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[i_node * block_size];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[i_node * block_size + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[i_node * block_size + 3];
    }
}

template <>
void CompressibleNavierStokesExplicit<3>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 3;
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 5;

    // Calculate the explicit residual vector
    VectorType rhs;
    CalculateRightHandSide(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[i_node * block_size];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[i_node * block_size + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[i_node * block_size + 4];
    }
}

// TODO: We still require to decide the shock capturing technique
template <>
void CompressibleNavierStokesExplicit<2>::CalculateShockCapturingValues(ElementDataStruct &rData) const
{
    rData.nu_sc = 0.0;
    rData.lambda_sc = 0.0;
}

// TODO: We still require to decide the shock capturing technique
template <>
void CompressibleNavierStokesExplicit<3>::CalculateShockCapturingValues(ElementDataStruct &rData) const
{
    rData.nu_sc = 0.0;
    rData.lambda_sc = 0.0;
}

// template<>
// double CompressibleNavierStokesExplicit<2>::ShockCapturingViscosity(const ElementDataStruct& rData) const
// {
//     const int dim = 2;
//     const int n_nodes = 3;
//     const int block_size = dim + 2;

//     const double h = rData.h;                                // Characteristic element size
//     const double alpha = 0.8;                               // Algorithm constant
//     const double tol = 0.001;

//     const BoundedMatrix<double, n_nodes, block_size>& U = rData.U;
//     const BoundedMatrix<double, n_nodes, dim>& f_ext = rData.f_ext;
//     const double gamma = rData.gamma;
//     double v_sc = 0.0;                                      //Shock capturing viscosity
//     BoundedMatrix<double, dim, 1> res_m;
//     res_m(0,0) =0; res_m(1,0) =0;

//     // Solution vector values from nodal data
//     // This is intentionally done in this way to limit the matrix acceses
//     // The notation U_i_j DOF j value in node i
//     const double &U_0_0 = rData.U(0, 0);
//     const double &U_0_1 = rData.U(0, 1);
//     const double &U_0_2 = rData.U(0, 2);
//     const double &U_0_3 = rData.U(0, 3);
//     const double &U_1_0 = rData.U(1, 0);
//     const double &U_1_1 = rData.U(1, 1);
//     const double &U_1_2 = rData.U(1, 2);
//     const double &U_1_3 = rData.U(1, 3);
//     const double &U_2_0 = rData.U(2, 0);
//     const double &U_2_1 = rData.U(2, 1);
//     const double &U_2_2 = rData.U(2, 2);
//     const double &U_2_3 = rData.U(2, 3);

//     // Hardcoded shape functions for linear triangular element
//     // This is explicitly done to minimize the allocation and matrix acceses
//     // The notation N_i_j means shape function for node j in Gauss pt. i
//     const double one_sixt = 1.0 / 6.0;
//     const double two_third = 2.0 / 3.0;
//     const double N_0_0 = one_sixt;
//     const double N_0_1 = one_sixt;
//     const double N_0_2 = two_third;
//     const double N_1_0 = one_sixt;
//     const double N_1_1 = two_third;
//     const double N_1_2 = one_sixt;
//     const double N_2_0 = two_third;
//     const double N_2_1 = one_sixt;
//     const double N_2_2 = one_sixt;

//     // Hardcoded shape functions gradients for linear triangular element
//     // This is explicitly done to minimize the matrix acceses
//     // The notation DN_i_j means shape function for node i in dimension j
//     const double &DN_DX_0_0 = rData.DN_DX(0, 0);
//     const double &DN_DX_0_1 = rData.DN_DX(0, 1);
//     const double &DN_DX_1_0 = rData.DN_DX(1, 0);
//     const double &DN_DX_1_1 = rData.DN_DX(1, 1);
//     const double &DN_DX_2_0 = rData.DN_DX(2, 0);
//     const double &DN_DX_2_1 = rData.DN_DX(2, 1);

//     // Get shape function values
//     const array_1d<double, n_nodes>& N = rData.N;
//     const BoundedMatrix<double, n_nodes, dim>& DN = rData.DN_DX;

//     // Auxiliary variables used in the calculation of the RHS
//     const array_1d<double, dim> f_gauss = prod(trans(f_ext), N);
//     const array_1d<double, block_size> U_gauss = prod(trans(U), N);
//     const BoundedMatrix<double,block_size,dim> grad_U = prod(trans(U), DN);     // Dfi/Dxj

//     //substitute_res_m_2D

//     double norm_res_m;
//     norm_res_m = sqrt(res_m(0,0)*res_m(0,0)+res_m(1,0)*res_m(1,0));

//     double norm_gradm = 0.0;                                    // Frobenius norm of momentum gradient
//     for (unsigned int i = 1; i < dim + 1; ++i){
//         for (unsigned int j = 0; j < dim; ++j) {
//             norm_gradm += grad_U(i,j)*grad_U(i,j);
//         }
//     }
//     norm_gradm = sqrt(norm_gradm);

//     if (norm_gradm>tol) {
//         v_sc = 0.5*h*alpha*(norm_res_m/norm_gradm);
//     }

//     return v_sc;
// }

// template<>
// double CompressibleNavierStokesExplicit<3>::ShockCapturingViscosity(const ElementDataStruct& rData) const
// {
//     const int dim = 3;
//     const int n_nodes = 4;
//     const int block_size = dim + 2;

//     const double h = rData.h;                                // Characteristic element size
//     const double alpha = 0.8;                               // Algorithm constant
//     const double tol = 0.001;

//     const BoundedMatrix<double, n_nodes, block_size>& U = rData.U;
//     const BoundedMatrix<double, n_nodes, dim>& f_ext = rData.f_ext;
//     const double gamma = rData.gamma;
//     double v_sc = 0.0;                                      //Shock capturing viscosity
//     BoundedMatrix<double, dim, 1> res_m;
//     res_m(0,0)= 0; res_m(1,0)= 0; res_m(2,0)= 0;

//     // Solution vector values from nodal data
//     // This is intentionally done in this way to limit the matrix acceses
//     // The notation U_i_j DOF j value in node i
//     const double &U_0_0 = rData.U(0, 0);
//     const double &U_0_1 = rData.U(0, 1);
//     const double &U_0_2 = rData.U(0, 2);
//     const double &U_0_3 = rData.U(0, 3);
//     const double &U_0_4 = rData.U(0, 4);
//     const double &U_1_0 = rData.U(1, 0);
//     const double &U_1_1 = rData.U(1, 1);
//     const double &U_1_2 = rData.U(1, 2);
//     const double &U_1_3 = rData.U(1, 3);
//     const double &U_1_4 = rData.U(1, 4);
//     const double &U_2_0 = rData.U(2, 0);
//     const double &U_2_1 = rData.U(2, 1);
//     const double &U_2_2 = rData.U(2, 2);
//     const double &U_2_3 = rData.U(2, 3);
//     const double &U_2_4 = rData.U(2, 4);
//     const double &U_3_0 = rData.U(3, 0);
//     const double &U_3_1 = rData.U(3, 1);
//     const double &U_3_2 = rData.U(3, 2);
//     const double &U_3_3 = rData.U(3, 3);
//     const double &U_3_4 = rData.U(3, 4);

//     // Hardcoded shape functions for linear tetrahedra element
//     // This is explicitly done to minimize the alocation and matrix acceses
//     // The notation N_i_j means shape function for node j in Gauss pt. i
//     const double N_0_0 = 0.58541020;
//     const double N_0_1 = 0.13819660;
//     const double N_0_2 = 0.13819660;
//     const double N_0_3 = 0.13819660;
//     const double N_1_0 = 0.13819660;
//     const double N_1_1 = 0.58541020;
//     const double N_1_2 = 0.13819660;
//     const double N_1_3 = 0.13819660;
//     const double N_2_0 = 0.13819660;
//     const double N_2_1 = 0.13819660;
//     const double N_2_2 = 0.58541020;
//     const double N_2_3 = 0.13819660;
//     const double N_3_0 = 0.13819660;
//     const double N_3_1 = 0.13819660;
//     const double N_3_2 = 0.13819660;
//     const double N_3_3 = 0.58541020;

//     // Hardcoded shape functions gradients for linear tetrahedra element
//     // This is explicitly done to minimize the matrix acceses
//     // The notation DN_i_j means shape function for node i in dimension j
//     const double &DN_DX_0_0 = rData.DN_DX(0, 0);
//     const double &DN_DX_0_1 = rData.DN_DX(0, 1);
//     const double &DN_DX_0_2 = rData.DN_DX(0, 2);
//     const double &DN_DX_1_0 = rData.DN_DX(1, 0);
//     const double &DN_DX_1_1 = rData.DN_DX(1, 1);
//     const double &DN_DX_1_2 = rData.DN_DX(1, 2);
//     const double &DN_DX_2_0 = rData.DN_DX(2, 0);
//     const double &DN_DX_2_1 = rData.DN_DX(2, 1);
//     const double &DN_DX_2_2 = rData.DN_DX(2, 2);
//     const double &DN_DX_3_0 = rData.DN_DX(3, 0);
//     const double &DN_DX_3_1 = rData.DN_DX(3, 1);
//     const double &DN_DX_3_2 = rData.DN_DX(3, 2);

//     // Get shape function values
//     const array_1d<double, n_nodes>& N = rData.N;
//     const BoundedMatrix<double, n_nodes, dim>& DN = rData.DN_DX;

//     // Auxiliary variables used in the calculation of the RHS
//     const array_1d<double, dim> f_gauss = prod(trans(f_ext), N);
//     const array_1d<double, block_size> U_gauss = prod(trans(U), N);
//     const BoundedMatrix<double, block_size, dim> grad_U = prod(trans(U), DN);     // Dfi/Dxj

//     //substitute_res_m_3D

//     double norm_res_m;
//     norm_res_m = sqrt(res_m(0,0)*res_m(0,0)+res_m(1,0)*res_m(1,0)+res_m(2,0)*res_m(2,0));

//     double norm_gradm = 0.0;                                    // Frobenius norm of momentum gradient
//     for (unsigned int i=1; i<dim+1; i++){
//         for (unsigned int j=0; j<dim; j++) {
//             norm_gradm += grad_U(i,j)*grad_U(i,j);
//         }
//     }
//     norm_gradm = sqrt(norm_gradm);

//     if (norm_gradm>tol) {
//         v_sc = 0.5*h*alpha*(norm_res_m/norm_gradm);
//     }

//     return v_sc;
// }

// template<>
// double CompressibleNavierStokesExplicit<2>::ShockCapturingConductivity(const ElementDataStruct& rData) const
// {
//     const int dim = 2;
//     const int n_nodes = 3;
//     const int block_size = dim + 2;

//     const double h = rData.h;                                // Characteristic element size
//     const double alpha = 0.8;                               // Algorithm constant
//     const double tol = 0.001;

//     const BoundedMatrix<double, n_nodes, block_size>& U = rData.U;
//     const BoundedMatrix<double, n_nodes, dim>& f_ext = rData.f_ext;
//     const array_1d<double, n_nodes>& r = rData.r;
//     const double gamma = rData.gamma;
//     double k_sc = 0.0;          // Shock Capturing Conductivity
//     BoundedMatrix<double, dim, 1> res_e;
//     res_e(0,0) = 0;

//     // Solution vector values from nodal data
//     // This is intentionally done in this way to limit the matrix acceses
//     // The notation U_i_j DOF j value in node i
//     const double &U_0_0 = rData.U(0, 0);
//     const double &U_0_1 = rData.U(0, 1);
//     const double &U_0_2 = rData.U(0, 2);
//     const double &U_0_3 = rData.U(0, 3);
//     const double &U_1_0 = rData.U(1, 0);
//     const double &U_1_1 = rData.U(1, 1);
//     const double &U_1_2 = rData.U(1, 2);
//     const double &U_1_3 = rData.U(1, 3);
//     const double &U_2_0 = rData.U(2, 0);
//     const double &U_2_1 = rData.U(2, 1);
//     const double &U_2_2 = rData.U(2, 2);
//     const double &U_2_3 = rData.U(2, 3);

//     // Hardcoded shape functions for linear triangular element
//     // This is explicitly done to minimize the allocation and matrix acceses
//     // The notation N_i_j means shape function for node j in Gauss pt. i
//     const double one_sixt = 1.0 / 6.0;
//     const double two_third = 2.0 / 3.0;
//     const double N_0_0 = one_sixt;
//     const double N_0_1 = one_sixt;
//     const double N_0_2 = two_third;
//     const double N_1_0 = one_sixt;
//     const double N_1_1 = two_third;
//     const double N_1_2 = one_sixt;
//     const double N_2_0 = two_third;
//     const double N_2_1 = one_sixt;
//     const double N_2_2 = one_sixt;

//     // Hardcoded shape functions gradients for linear triangular element
//     // This is explicitly done to minimize the matrix acceses
//     // The notation DN_i_j means shape function for node i in dimension j
//     const double &DN_DX_0_0 = rData.DN_DX(0, 0);
//     const double &DN_DX_0_1 = rData.DN_DX(0, 1);
//     const double &DN_DX_1_0 = rData.DN_DX(1, 0);
//     const double &DN_DX_1_1 = rData.DN_DX(1, 1);
//     const double &DN_DX_2_0 = rData.DN_DX(2, 0);
//     const double &DN_DX_2_1 = rData.DN_DX(2, 1);

//     // Get shape function values
//     const array_1d<double, n_nodes>& N = rData.N;
//     const BoundedMatrix<double, n_nodes, dim>& DN = rData.DN_DX;

//     // Auxiliary variables used in the calculation of the RHS
//     const array_1d<double, dim> f_gauss = prod(trans(f_ext), N);
//     const array_1d<double, block_size> U_gauss = prod(trans(U), N);
//     const BoundedMatrix<double, block_size, dim> grad_U = prod(trans(U), DN);     // Dfi/Dxj

//     //substitute_res_e_2D

//     double norm_res_e;
//     norm_res_e = sqrt(res_e(0,0)*res_e(0,0));

//     double norm_grade = 0.0;              // Frobenius norm of total energy gradient
//     for (unsigned int i=0; i<dim; i++) {
//         norm_grade += grad_U(dim+1,i)*grad_U(dim+1,i);
//     }
//     norm_grade = sqrt(norm_grade);

//     if (norm_grade > tol) {
//         k_sc = 0.5*h*alpha*(norm_res_e/norm_grade);
//     }

//     return k_sc;
// }

// template<>
// double CompressibleNavierStokesExplicit<3>::ShockCapturingConductivity(const ElementDataStruct& rData) const
// {
//     const int dim = 3;
//     const int n_nodes = 4;
//     const int block_size = dim + 2;

//     const double h = rData.h;                                // Characteristic element size
//     const double alpha = 0.8;                               // Algorithm constant
//     const double tol = 0.001;

//     const BoundedMatrix<double, n_nodes, block_size>& U = rData.U;
//     const BoundedMatrix<double, n_nodes, dim>& f_ext = rData.f_ext;
//     const array_1d<double, n_nodes>& r = rData.r;
//     const double gamma = rData.gamma;
//     double k_sc = 0.0;          // Shock Capturing Conductivity
//     BoundedMatrix<double, dim, 1> res_e;
//     res_e(0,0) = 0;

//     // Solution vector values from nodal data
//     // This is intentionally done in this way to limit the matrix acceses
//     // The notation U_i_j DOF j value in node i
//     const double &U_0_0 = rData.U(0, 0);
//     const double &U_0_1 = rData.U(0, 1);
//     const double &U_0_2 = rData.U(0, 2);
//     const double &U_0_3 = rData.U(0, 3);
//     const double &U_0_4 = rData.U(0, 4);
//     const double &U_1_0 = rData.U(1, 0);
//     const double &U_1_1 = rData.U(1, 1);
//     const double &U_1_2 = rData.U(1, 2);
//     const double &U_1_3 = rData.U(1, 3);
//     const double &U_1_4 = rData.U(1, 4);
//     const double &U_2_0 = rData.U(2, 0);
//     const double &U_2_1 = rData.U(2, 1);
//     const double &U_2_2 = rData.U(2, 2);
//     const double &U_2_3 = rData.U(2, 3);
//     const double &U_2_4 = rData.U(2, 4);
//     const double &U_3_0 = rData.U(3, 0);
//     const double &U_3_1 = rData.U(3, 1);
//     const double &U_3_2 = rData.U(3, 2);
//     const double &U_3_3 = rData.U(3, 3);
//     const double &U_3_4 = rData.U(3, 4);

//     // Hardcoded shape functions for linear tetrahedra element
//     // This is explicitly done to minimize the alocation and matrix acceses
//     // The notation N_i_j means shape function for node j in Gauss pt. i
//     const double N_0_0 = 0.58541020;
//     const double N_0_1 = 0.13819660;
//     const double N_0_2 = 0.13819660;
//     const double N_0_3 = 0.13819660;
//     const double N_1_0 = 0.13819660;
//     const double N_1_1 = 0.58541020;
//     const double N_1_2 = 0.13819660;
//     const double N_1_3 = 0.13819660;
//     const double N_2_0 = 0.13819660;
//     const double N_2_1 = 0.13819660;
//     const double N_2_2 = 0.58541020;
//     const double N_2_3 = 0.13819660;
//     const double N_3_0 = 0.13819660;
//     const double N_3_1 = 0.13819660;
//     const double N_3_2 = 0.13819660;
//     const double N_3_3 = 0.58541020;

//     // Hardcoded shape functions gradients for linear tetrahedra element
//     // This is explicitly done to minimize the matrix acceses
//     // The notation DN_i_j means shape function for node i in dimension j
//     const double &DN_DX_0_0 = rData.DN_DX(0, 0);
//     const double &DN_DX_0_1 = rData.DN_DX(0, 1);
//     const double &DN_DX_0_2 = rData.DN_DX(0, 2);
//     const double &DN_DX_1_0 = rData.DN_DX(1, 0);
//     const double &DN_DX_1_1 = rData.DN_DX(1, 1);
//     const double &DN_DX_1_2 = rData.DN_DX(1, 2);
//     const double &DN_DX_2_0 = rData.DN_DX(2, 0);
//     const double &DN_DX_2_1 = rData.DN_DX(2, 1);
//     const double &DN_DX_2_2 = rData.DN_DX(2, 2);
//     const double &DN_DX_3_0 = rData.DN_DX(3, 0);
//     const double &DN_DX_3_1 = rData.DN_DX(3, 1);
//     const double &DN_DX_3_2 = rData.DN_DX(3, 2);

//     // Get shape function values
//     const array_1d<double, n_nodes>& N = rData.N;
//     const BoundedMatrix<double, n_nodes, dim>& DN = rData.DN_DX;

//     // Auxiliary variables used in the calculation of the RHS
//     const array_1d<double, dim> f_gauss = prod(trans(f_ext), N);
//     const array_1d<double, block_size> U_gauss = prod(trans(U), N);
//     const BoundedMatrix<double, block_size, dim> grad_U = prod(trans(U), DN);     // Dfi/Dxj

//     //substitute_res_e_3D

//     double norm_res_e;
//     norm_res_e = sqrt(res_e(0,0)*res_e(0,0));

//     double norm_grade = 0.0;              // Frobenius norm of total energy gradient
//     for (unsigned int i=0; i<dim; i++) {
//         norm_grade += grad_U(dim+1,i)*grad_U(dim+1,i);
//     }
//     norm_grade = sqrt(norm_grade);

//     if (norm_grade > tol) {
//         k_sc = 0.5*h*alpha*(norm_res_e/norm_grade);
//     }

//     return k_sc;
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNavierStokesExplicit<2>;
template class CompressibleNavierStokesExplicit<3>;

}
