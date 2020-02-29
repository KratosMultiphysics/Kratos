//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
void CompressibleNavierStokesExplicit<2>::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 4;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (rResult.size() != dof_size) {
        rResult.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNavierStokesExplicit<3>::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 5;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (rResult.size() != dof_size) {
        rResult.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Z, mom_pos + 2).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
int CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::Check(const ProcessInfo &rCurrentProcessInfo)
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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Getting data for the given geometry
    const auto& r_geometry = GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, rData.DN_DX, rData.N, rData.volume);

    // Compute element size
    rData.h = CalculateElementSize(rData.DN_DX);

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

template <unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
double CompressibleNavierStokesExplicit<TDim, TNumNodes, TBlockSize>::CalculateElementSize(const BoundedMatrix<double,TNumNodes, TDim>& rDN_DX)
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
void CompressibleNavierStokesExplicit<2>::CalculateRightHandSideInternal(
    BoundedVector<double, 12> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

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
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template<>
void CompressibleNavierStokesExplicit<3>::CalculateRightHandSideInternal(
    BoundedVector<double, 20> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

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
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("")
}

template <>
void CompressibleNavierStokesExplicit<2>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 2;
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 4;

    // Calculate the explicit residual vector
    BoundedVector<double, 12> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + 3];
    }
}

template <>
void CompressibleNavierStokesExplicit<3>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 3;
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 5;

    // Calculate the explicit residual vector
    BoundedVector<double, 20> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + 4];
    }
}

template <>
void CompressibleNavierStokesExplicit<2>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 4;

    // Initialize and fill the mass matrix values
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six; rMassMatrix(0, 4) = one_twelve; rMassMatrix(0, 8) = one_twelve;
    rMassMatrix(1, 1) = one_six; rMassMatrix(1, 5) = one_twelve; rMassMatrix(1, 9) = one_twelve;
    rMassMatrix(2, 2) = one_six; rMassMatrix(2, 6) = one_twelve; rMassMatrix(2, 10) = one_twelve;
    rMassMatrix(3, 3) = one_six; rMassMatrix(3, 7) = one_twelve; rMassMatrix(3, 11) = one_twelve;
    rMassMatrix(4, 0) = one_twelve; rMassMatrix(4, 4) = one_six; rMassMatrix(4, 8) = one_twelve;
    rMassMatrix(5, 1) = one_twelve; rMassMatrix(5, 5) = one_six; rMassMatrix(5, 9) = one_twelve;
    rMassMatrix(6, 2) = one_twelve; rMassMatrix(6, 6) = one_six; rMassMatrix(6, 10) = one_twelve;
    rMassMatrix(7, 3) = one_twelve; rMassMatrix(7, 7) = one_six; rMassMatrix(7, 11) = one_twelve;
    rMassMatrix(8, 0) = one_twelve; rMassMatrix(8, 4) = one_twelve; rMassMatrix(8, 8) = one_six;
    rMassMatrix(9, 1) = one_twelve; rMassMatrix(9, 5) = one_twelve; rMassMatrix(9, 9) = one_six;
    rMassMatrix(10, 2) = one_twelve; rMassMatrix(10, 6) = one_twelve; rMassMatrix(10, 10) = one_six;
    rMassMatrix(11, 3) = one_twelve; rMassMatrix(11, 7) = one_twelve; rMassMatrix(11, 11) = one_six;

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

template <>
void CompressibleNavierStokesExplicit<3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 5;

    // Initialize and fill the mass matrix values
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_ten; rMassMatrix(0, 5) = one_twenty; rMassMatrix(0, 10) = one_twenty; rMassMatrix(0,15) = one_twenty;
    rMassMatrix(1, 1) = one_ten; rMassMatrix(1, 6) = one_twenty; rMassMatrix(1, 11) = one_twenty; rMassMatrix(1,16) = one_twenty;
    rMassMatrix(2, 2) = one_ten; rMassMatrix(2, 7) = one_twenty; rMassMatrix(2, 12) = one_twenty; rMassMatrix(2,17) = one_twenty;
    rMassMatrix(3, 3) = one_ten; rMassMatrix(3, 8) = one_twenty; rMassMatrix(3, 13) = one_twenty; rMassMatrix(3,18) = one_twenty;
    rMassMatrix(4, 4) = one_ten; rMassMatrix(4, 9) = one_twenty; rMassMatrix(4, 14) = one_twenty; rMassMatrix(4,19) = one_twenty;
    rMassMatrix(5, 0) = one_twenty; rMassMatrix(5, 5) = one_ten; rMassMatrix(5, 10) = one_twenty; rMassMatrix(5,15) = one_twenty;
    rMassMatrix(6, 1) = one_twenty; rMassMatrix(6, 6) = one_ten; rMassMatrix(6, 11) = one_twenty; rMassMatrix(6,16) = one_twenty;
    rMassMatrix(7, 2) = one_twenty; rMassMatrix(7, 7) = one_ten; rMassMatrix(7, 12) = one_twenty; rMassMatrix(7,17) = one_twenty;
    rMassMatrix(8, 3) = one_twenty; rMassMatrix(8, 8) = one_ten; rMassMatrix(8, 13) = one_twenty; rMassMatrix(8,18) = one_twenty;
    rMassMatrix(9, 4) = one_twenty; rMassMatrix(9, 9) = one_ten; rMassMatrix(9, 14) = one_twenty; rMassMatrix(9,19) = one_twenty;
    rMassMatrix(10, 0) = one_twenty; rMassMatrix(10, 5) = one_twenty; rMassMatrix(10, 10) = one_ten; rMassMatrix(10,15) = one_twenty;
    rMassMatrix(11, 1) = one_twenty; rMassMatrix(11, 6) = one_twenty; rMassMatrix(11, 11) = one_ten; rMassMatrix(11,16) = one_twenty;
    rMassMatrix(12, 2) = one_twenty; rMassMatrix(12, 7) = one_twenty; rMassMatrix(12, 12) = one_ten; rMassMatrix(12,17) = one_twenty;
    rMassMatrix(13, 3) = one_twenty; rMassMatrix(13, 8) = one_twenty; rMassMatrix(13, 13) = one_ten; rMassMatrix(13,18) = one_twenty;
    rMassMatrix(14, 4) = one_twenty; rMassMatrix(14, 9) = one_twenty; rMassMatrix(14, 14) = one_ten; rMassMatrix(14,19) = one_twenty;
    rMassMatrix(15, 0) = one_twenty; rMassMatrix(15, 5) = one_twenty; rMassMatrix(15, 10) = one_twenty; rMassMatrix(15,15) = one_ten;
    rMassMatrix(16, 1) = one_twenty; rMassMatrix(16, 6) = one_twenty; rMassMatrix(16, 11) = one_twenty; rMassMatrix(16,16) = one_ten;
    rMassMatrix(17, 2) = one_twenty; rMassMatrix(17, 7) = one_twenty; rMassMatrix(17, 12) = one_twenty; rMassMatrix(17,17) = one_ten;
    rMassMatrix(18, 3) = one_twenty; rMassMatrix(18, 8) = one_twenty; rMassMatrix(18, 13) = one_twenty; rMassMatrix(18,18) = one_ten;
    rMassMatrix(19, 4) = one_twenty; rMassMatrix(19, 9) = one_twenty; rMassMatrix(19, 14) = one_twenty; rMassMatrix(19,19) = one_ten;

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Volume();
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNavierStokesExplicit<2>;
template class CompressibleNavierStokesExplicit<3>;

}
