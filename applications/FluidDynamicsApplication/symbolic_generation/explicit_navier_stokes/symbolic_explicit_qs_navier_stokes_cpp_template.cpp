//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

// Application includes
#include "symbolic_explicit_qs_navier_stokes.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
SymbolicExplicitQSNavierStokes<TDim,TNumNodes>::SymbolicExplicitQSNavierStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
SymbolicExplicitQSNavierStokes<TDim,TNumNodes>::SymbolicExplicitQSNavierStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void SymbolicExplicitQSNavierStokes<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const unsigned int dof_size = TNumNodes * (TDim+1);
    if (rResult.size() != dof_size)
        rResult.resize(dof_size, false);

    const GeometryType& r_geometry = this->GetGeometry();
    unsigned int local_index = 0;

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        if (TDim == 3) rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE,ppos).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void SymbolicExplicitQSNavierStokes<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const unsigned int dof_size = TNumNodes * (TDim+1);
    if (rElementalDofList.size() != dof_size)
        rElementalDofList.resize(dof_size);

    const GeometryType& r_geometry = this->GetGeometry();
    unsigned int local_index = 0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y);
        if (TDim == 3) rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void SymbolicExplicitQSNavierStokes<2>::CalculateRightHandSideInternal(
    BoundedVector<double, 9> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &h = data.h;
    const auto &f = data.forcing;
    const auto &nu = data.nu;
    const auto &vconv = data.velocity_convective;
    const auto &v = data.velocity;
    const auto &p = data.pressure;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

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

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void SymbolicExplicitQSNavierStokes<3>::CalculateRightHandSideInternal(
    BoundedVector<double, 16> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &h = data.h;
    const auto &f = data.forcing;
    const auto &nu = data.nu;
    const auto &vconv = data.velocity_convective;
    const auto &v = data.velocity;
    const auto &p = data.pressure;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Hardcoded shape functions gradients for linear triangular element
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

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void SymbolicExplicitQSNavierStokes<2>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr IndexType dim = 2;
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 3;

    // Calculate the explicit residual vector
    BoundedVector<double, 9> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + d];
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void SymbolicExplicitQSNavierStokes<3>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr IndexType dim = 3;
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 4;

    // Calculate the explicit residual vector
    BoundedVector<double, 16> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + d];
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void SymbolicExplicitQSNavierStokes<2>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 3;

    // Initialize and fill the mass matrix values
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six;    rMassMatrix(0, 3) = one_twelve; rMassMatrix(0, 6) = one_twelve;
    rMassMatrix(1, 1) = one_six;    rMassMatrix(1, 4) = one_twelve; rMassMatrix(1, 7) = one_twelve;
    rMassMatrix(2, 2) = one_six;    rMassMatrix(2, 5) = one_twelve; rMassMatrix(2, 8) = one_twelve;
    rMassMatrix(3, 0) = one_twelve; rMassMatrix(3, 3) = one_six;    rMassMatrix(3, 6) = one_twelve;
    rMassMatrix(4, 1) = one_twelve; rMassMatrix(4, 4) = one_six;    rMassMatrix(4, 7) = one_twelve;
    rMassMatrix(5, 2) = one_twelve; rMassMatrix(5, 5) = one_six;    rMassMatrix(5, 8) = one_twelve;
    rMassMatrix(6, 0) = one_twelve; rMassMatrix(6, 3) = one_twelve; rMassMatrix(6, 6) = one_six;
    rMassMatrix(7, 1) = one_twelve; rMassMatrix(7, 4) = one_twelve; rMassMatrix(7, 7) = one_six;
    rMassMatrix(8, 2) = one_twelve; rMassMatrix(8, 5) = one_twelve; rMassMatrix(8, 8) = one_six;

    // Here we assume that all the Gauss points have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

/***********************************************************************************/

template<>
void SymbolicExplicitQSNavierStokes<3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 4;

    // Initialize and fill the mass matrix values
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_ten;     rMassMatrix(0, 4) = one_twenty;  rMassMatrix(0, 8) = one_twenty;   rMassMatrix(0,12) = one_twenty;
    rMassMatrix(1, 1) = one_ten;     rMassMatrix(1, 5) = one_twenty;  rMassMatrix(1, 9) = one_twenty;   rMassMatrix(1,13) = one_twenty;
    rMassMatrix(2, 2) = one_ten;     rMassMatrix(2, 6) = one_twenty;  rMassMatrix(2, 10) = one_twenty;  rMassMatrix(2,14) = one_twenty;
    rMassMatrix(3, 3) = one_ten;     rMassMatrix(3, 7) = one_twenty;  rMassMatrix(3, 11) = one_twenty;  rMassMatrix(3,15) = one_twenty;
    rMassMatrix(4, 0) = one_twenty;  rMassMatrix(4, 4) = one_ten;     rMassMatrix(4, 8) = one_twenty;   rMassMatrix(4,12) = one_twenty;
    rMassMatrix(5, 1) = one_twenty;  rMassMatrix(5, 5) = one_ten;     rMassMatrix(5, 9) = one_twenty;   rMassMatrix(5,13) = one_twenty;
    rMassMatrix(6, 2) = one_twenty;  rMassMatrix(6, 6) = one_ten;     rMassMatrix(6, 10) = one_twenty;  rMassMatrix(6,14) = one_twenty;
    rMassMatrix(7, 3) = one_twenty;  rMassMatrix(7, 7) = one_ten;     rMassMatrix(7, 11) = one_twenty;  rMassMatrix(7,15) = one_twenty;
    rMassMatrix(8, 0) = one_twenty;  rMassMatrix(8, 4) = one_twenty;  rMassMatrix(8, 8) = one_ten;      rMassMatrix(8,12) = one_twenty;
    rMassMatrix(9, 1) = one_twenty;  rMassMatrix(9, 5) = one_twenty;  rMassMatrix(9, 9) = one_ten;      rMassMatrix(9,13) = one_twenty;
    rMassMatrix(10, 2) = one_twenty; rMassMatrix(10, 6) = one_twenty; rMassMatrix(10, 10) = one_ten;    rMassMatrix(10,14) = one_twenty;
    rMassMatrix(11, 3) = one_twenty; rMassMatrix(11, 7) = one_twenty; rMassMatrix(11, 11) = one_ten;    rMassMatrix(11,15) = one_twenty;
    rMassMatrix(12, 0) = one_twenty; rMassMatrix(12, 4) = one_twenty; rMassMatrix(12, 8) = one_twenty;  rMassMatrix(12,12) = one_ten;
    rMassMatrix(13, 1) = one_twenty; rMassMatrix(13, 5) = one_twenty; rMassMatrix(13, 9) = one_twenty;  rMassMatrix(13,13) = one_ten;
    rMassMatrix(14, 2) = one_twenty; rMassMatrix(14, 6) = one_twenty; rMassMatrix(14, 10) = one_twenty; rMassMatrix(14,14) = one_ten;
    rMassMatrix(15, 3) = one_twenty; rMassMatrix(15, 7) = one_twenty; rMassMatrix(15, 11) = one_twenty; rMassMatrix(15,15) = one_ten;

    // Here we assume that all the Gauss points have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
int SymbolicExplicitQSNavierStokes<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;
    return 0;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void SymbolicExplicitQSNavierStokes<TDim, TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
double SymbolicExplicitQSNavierStokes<TDim, TNumNodes>::CalculateElementSize(
    const BoundedMatrix<double,TNumNodes, TDim>& rDN_DX)
{
    KRATOS_TRY;

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

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class SymbolicExplicitQSNavierStokes<2>;
template class SymbolicExplicitQSNavierStokes<3>;

/***********************************************************************************/
/***********************************************************************************/

}