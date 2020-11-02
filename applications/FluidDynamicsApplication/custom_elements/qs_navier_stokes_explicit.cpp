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
#include "qs_navier_stokes_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
QSNavierStokesExplicit<TDim,TNumNodes>::QSNavierStokesExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
QSNavierStokesExplicit<TDim,TNumNodes>::QSNavierStokesExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesExplicit<TDim,TNumNodes>::EquationIdVector(
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
void QSNavierStokesExplicit<TDim,TNumNodes>::GetDofList(
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
void QSNavierStokesExplicit<2>::CalculateRightHandSideInternal(
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

    const double crRightHandSideBoundedVector0 =             0.25*f(1,0);
const double crRightHandSideBoundedVector1 =             0.166666666666667*vconv(0,0);
const double crRightHandSideBoundedVector2 =             0.166666666666667*vconv(2,0);
const double crRightHandSideBoundedVector3 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + 0.666666666666667*vconv(1,0);
const double crRightHandSideBoundedVector4 =             DN_DX_0_0*v(0,0) + DN_DX_1_0*v(1,0) + DN_DX_2_0*v(2,0);
const double crRightHandSideBoundedVector5 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector6 =             -0.166666666666667*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector7 =             0.166666666666667*vconv(0,1);
const double crRightHandSideBoundedVector8 =             0.166666666666667*vconv(2,1);
const double crRightHandSideBoundedVector9 =             crRightHandSideBoundedVector7 + crRightHandSideBoundedVector8 + 0.666666666666667*vconv(1,1);
const double crRightHandSideBoundedVector10 =             DN_DX_0_1*v(0,0) + DN_DX_1_1*v(1,0) + DN_DX_2_1*v(2,0);
const double crRightHandSideBoundedVector11 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector12 =             -0.166666666666667*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector13 =             0.25*f(2,0);
const double crRightHandSideBoundedVector14 =             0.166666666666667*p[0];
const double crRightHandSideBoundedVector15 =             0.166666666666667*p[1];
const double crRightHandSideBoundedVector16 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector15 + 0.666666666666667*p[2];
const double crRightHandSideBoundedVector17 =             0.166666666666667*p[2];
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector17 + 0.666666666666667*p[1];
const double crRightHandSideBoundedVector19 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector17 + 0.666666666666667*p[0];
const double crRightHandSideBoundedVector20 =             3*crRightHandSideBoundedVector4*nu;
const double crRightHandSideBoundedVector21 =             0.166666666666667*vconv(1,0);
const double crRightHandSideBoundedVector22 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector21 + 0.666666666666667*vconv(2,0);
const double crRightHandSideBoundedVector23 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector24 =             -0.166666666666667*crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector21 + 0.666666666666667*vconv(0,0);
const double crRightHandSideBoundedVector26 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector27 =             0.166666666666667*vconv(1,1);
const double crRightHandSideBoundedVector28 =             crRightHandSideBoundedVector27 + crRightHandSideBoundedVector7 + 0.666666666666667*vconv(2,1);
const double crRightHandSideBoundedVector29 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector30 =             -0.166666666666667*crRightHandSideBoundedVector29;
const double crRightHandSideBoundedVector31 =             crRightHandSideBoundedVector27 + crRightHandSideBoundedVector8 + 0.666666666666667*vconv(0,1);
const double crRightHandSideBoundedVector32 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector31;
const double crRightHandSideBoundedVector33 =             DN_DX_0_0*v(0,1) + DN_DX_1_0*v(1,1) + DN_DX_2_0*v(2,1);
const double crRightHandSideBoundedVector34 =             3*nu*(crRightHandSideBoundedVector10 + crRightHandSideBoundedVector33);
const double crRightHandSideBoundedVector35 =             0.25*f(1,1);
const double crRightHandSideBoundedVector36 =             crRightHandSideBoundedVector3*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector37 =             -0.166666666666667*crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector38 =             DN_DX_0_1*v(0,1) + DN_DX_1_1*v(1,1) + DN_DX_2_1*v(2,1);
const double crRightHandSideBoundedVector39 =             crRightHandSideBoundedVector38*crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector40 =             -0.166666666666667*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector41 =             0.25*f(2,1);
const double crRightHandSideBoundedVector42 =             3*crRightHandSideBoundedVector38*nu;
const double crRightHandSideBoundedVector43 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector44 =             -0.166666666666667*crRightHandSideBoundedVector43;
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector28*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector47 =             -0.166666666666667*crRightHandSideBoundedVector46;
const double crRightHandSideBoundedVector48 =             crRightHandSideBoundedVector31*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector38 + crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector50 =             -0.166666666666667*crRightHandSideBoundedVector26 - 0.166666666666667*crRightHandSideBoundedVector32 + 0.25*f(0,0);
const double crRightHandSideBoundedVector51 =             -0.166666666666667*crRightHandSideBoundedVector45 - 0.166666666666667*crRightHandSideBoundedVector48 + 0.25*f(0,1);
const double crRightHandSideBoundedVector52 =             -1.0*crRightHandSideBoundedVector49;
            rRightHandSideBoundedVector[0]=DN_DX_0_0*crRightHandSideBoundedVector16 + DN_DX_0_0*crRightHandSideBoundedVector18 + DN_DX_0_0*crRightHandSideBoundedVector19 - DN_DX_0_0*crRightHandSideBoundedVector20 - DN_DX_0_1*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13 + crRightHandSideBoundedVector24 - 0.666666666666667*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector30 - 0.666666666666667*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector6 + 0.5*f(0,0);
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector34 + DN_DX_0_1*crRightHandSideBoundedVector16 + DN_DX_0_1*crRightHandSideBoundedVector18 + DN_DX_0_1*crRightHandSideBoundedVector19 - DN_DX_0_1*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector35 + crRightHandSideBoundedVector37 + crRightHandSideBoundedVector40 + crRightHandSideBoundedVector41 + crRightHandSideBoundedVector44 - 0.666666666666667*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector47 - 0.666666666666667*crRightHandSideBoundedVector48 + 0.5*f(0,1);
            rRightHandSideBoundedVector[2]=-1.0*crRightHandSideBoundedVector49;
            rRightHandSideBoundedVector[3]=DN_DX_1_0*crRightHandSideBoundedVector16 + DN_DX_1_0*crRightHandSideBoundedVector18 + DN_DX_1_0*crRightHandSideBoundedVector19 - DN_DX_1_0*crRightHandSideBoundedVector20 - DN_DX_1_1*crRightHandSideBoundedVector34 - 0.666666666666667*crRightHandSideBoundedVector11 + crRightHandSideBoundedVector13 + crRightHandSideBoundedVector24 + crRightHandSideBoundedVector30 - 0.666666666666667*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector50 + 0.5*f(1,0);
            rRightHandSideBoundedVector[4]=-DN_DX_1_0*crRightHandSideBoundedVector34 + DN_DX_1_1*crRightHandSideBoundedVector16 + DN_DX_1_1*crRightHandSideBoundedVector18 + DN_DX_1_1*crRightHandSideBoundedVector19 - DN_DX_1_1*crRightHandSideBoundedVector42 - 0.666666666666667*crRightHandSideBoundedVector36 - 0.666666666666667*crRightHandSideBoundedVector39 + crRightHandSideBoundedVector41 + crRightHandSideBoundedVector44 + crRightHandSideBoundedVector47 + crRightHandSideBoundedVector51 + 0.5*f(1,1);
            rRightHandSideBoundedVector[5]=crRightHandSideBoundedVector52;
            rRightHandSideBoundedVector[6]=DN_DX_2_0*crRightHandSideBoundedVector16 + DN_DX_2_0*crRightHandSideBoundedVector18 + DN_DX_2_0*crRightHandSideBoundedVector19 - DN_DX_2_0*crRightHandSideBoundedVector20 - DN_DX_2_1*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector12 - 0.666666666666667*crRightHandSideBoundedVector23 - 0.666666666666667*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector50 + crRightHandSideBoundedVector6 + 0.5*f(2,0);
            rRightHandSideBoundedVector[7]=-DN_DX_2_0*crRightHandSideBoundedVector34 + DN_DX_2_1*crRightHandSideBoundedVector16 + DN_DX_2_1*crRightHandSideBoundedVector18 + DN_DX_2_1*crRightHandSideBoundedVector19 - DN_DX_2_1*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector35 + crRightHandSideBoundedVector37 + crRightHandSideBoundedVector40 - 0.666666666666667*crRightHandSideBoundedVector43 - 0.666666666666667*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector51 + 0.5*f(2,1);
            rRightHandSideBoundedVector[8]=crRightHandSideBoundedVector52;


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3>::CalculateRightHandSideInternal(
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

    const double crRightHandSideBoundedVector0 =             0.19999999899376*f(1,0);
const double crRightHandSideBoundedVector1 =             0.1381966*vconv(0,0);
const double crRightHandSideBoundedVector2 =             0.1381966*vconv(2,0);
const double crRightHandSideBoundedVector3 =             0.1381966*vconv(3,0);
const double crRightHandSideBoundedVector4 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector5 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector4 + 0.5854102*vconv(1,0);
const double crRightHandSideBoundedVector6 =             DN_DX_0_0*v(0,0);
const double crRightHandSideBoundedVector7 =             DN_DX_1_0*v(1,0);
const double crRightHandSideBoundedVector8 =             DN_DX_2_0*v(2,0);
const double crRightHandSideBoundedVector9 =             DN_DX_3_0*v(3,0);
const double crRightHandSideBoundedVector10 =             crRightHandSideBoundedVector6 + crRightHandSideBoundedVector7 + crRightHandSideBoundedVector8 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector12 =             -0.1381966*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector13 =             0.1381966*vconv(0,1);
const double crRightHandSideBoundedVector14 =             0.1381966*vconv(2,1);
const double crRightHandSideBoundedVector15 =             0.1381966*vconv(3,1);
const double crRightHandSideBoundedVector16 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector17 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector16 + 0.5854102*vconv(1,1);
const double crRightHandSideBoundedVector18 =             DN_DX_0_1*v(0,0) + DN_DX_1_1*v(1,0) + DN_DX_2_1*v(2,0) + DN_DX_3_1*v(3,0);
const double crRightHandSideBoundedVector19 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector20 =             -0.1381966*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector21 =             0.1381966*vconv(0,2);
const double crRightHandSideBoundedVector22 =             0.1381966*vconv(2,2);
const double crRightHandSideBoundedVector23 =             0.1381966*vconv(3,2);
const double crRightHandSideBoundedVector24 =             crRightHandSideBoundedVector22 + crRightHandSideBoundedVector23;
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector21 + crRightHandSideBoundedVector24 + 0.5854102*vconv(1,2);
const double crRightHandSideBoundedVector26 =             DN_DX_0_2*v(0,0) + DN_DX_1_2*v(1,0) + DN_DX_2_2*v(2,0) + DN_DX_3_2*v(3,0);
const double crRightHandSideBoundedVector27 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector28 =             -0.1381966*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector29 =             0.19999999899376*f(2,0);
const double crRightHandSideBoundedVector30 =             0.19999999899376*f(3,0);
const double crRightHandSideBoundedVector31 =             0.1381966*p[0];
const double crRightHandSideBoundedVector32 =             0.1381966*p[1];
const double crRightHandSideBoundedVector33 =             crRightHandSideBoundedVector31 + crRightHandSideBoundedVector32;
const double crRightHandSideBoundedVector34 =             0.1381966*p[2];
const double crRightHandSideBoundedVector35 =             crRightHandSideBoundedVector33 + crRightHandSideBoundedVector34 + 0.5854102*p[3];
const double crRightHandSideBoundedVector36 =             0.1381966*p[3];
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector33 + crRightHandSideBoundedVector36 + 0.5854102*p[2];
const double crRightHandSideBoundedVector38 =             crRightHandSideBoundedVector34 + crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector39 =             crRightHandSideBoundedVector31 + crRightHandSideBoundedVector38 + 0.5854102*p[1];
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector32 + crRightHandSideBoundedVector38 + 0.5854102*p[0];
const double crRightHandSideBoundedVector41 =             4*crRightHandSideBoundedVector10*nu;
const double crRightHandSideBoundedVector42 =             0.1381966*vconv(1,0);
const double crRightHandSideBoundedVector43 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector44 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector43 + 0.5854102*vconv(3,0);
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector44;
const double crRightHandSideBoundedVector46 =             -0.1381966*crRightHandSideBoundedVector45;
const double crRightHandSideBoundedVector47 =             crRightHandSideBoundedVector3 + crRightHandSideBoundedVector43 + 0.5854102*vconv(2,0);
const double crRightHandSideBoundedVector48 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector49 =             -0.1381966*crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector4 + crRightHandSideBoundedVector42 + 0.5854102*vconv(0,0);
const double crRightHandSideBoundedVector51 =             crRightHandSideBoundedVector10*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector52 =             0.1381966*vconv(1,1);
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector52;
const double crRightHandSideBoundedVector54 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector53 + 0.5854102*vconv(3,1);
const double crRightHandSideBoundedVector55 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector56 =             -0.1381966*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector57 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector53 + 0.5854102*vconv(2,1);
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector59 =             -0.1381966*crRightHandSideBoundedVector58;
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector16 + crRightHandSideBoundedVector52 + 0.5854102*vconv(0,1);
const double crRightHandSideBoundedVector61 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector62 =             0.1381966*vconv(1,2);
const double crRightHandSideBoundedVector63 =             crRightHandSideBoundedVector21 + crRightHandSideBoundedVector62;
const double crRightHandSideBoundedVector64 =             crRightHandSideBoundedVector22 + crRightHandSideBoundedVector63 + 0.5854102*vconv(3,2);
const double crRightHandSideBoundedVector65 =             crRightHandSideBoundedVector26*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector66 =             -0.1381966*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector67 =             crRightHandSideBoundedVector23 + crRightHandSideBoundedVector63 + 0.5854102*vconv(2,2);
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector26*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector69 =             -0.1381966*crRightHandSideBoundedVector68;
const double crRightHandSideBoundedVector70 =             crRightHandSideBoundedVector24 + crRightHandSideBoundedVector62 + 0.5854102*vconv(0,2);
const double crRightHandSideBoundedVector71 =             crRightHandSideBoundedVector26*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector72 =             4*DN_DX_0_1*nu;
const double crRightHandSideBoundedVector73 =             DN_DX_0_0*v(0,1) + DN_DX_1_0*v(1,1) + DN_DX_2_0*v(2,1) + DN_DX_3_0*v(3,1);
const double crRightHandSideBoundedVector74 =             crRightHandSideBoundedVector18 + crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector75 =             4*DN_DX_0_2*nu;
const double crRightHandSideBoundedVector76 =             DN_DX_0_0*v(0,2) + DN_DX_1_0*v(1,2) + DN_DX_2_0*v(2,2) + DN_DX_3_0*v(3,2);
const double crRightHandSideBoundedVector77 =             crRightHandSideBoundedVector26 + crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector78 =             0.19999999899376*f(1,1);
const double crRightHandSideBoundedVector79 =             crRightHandSideBoundedVector5*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector80 =             -0.1381966*crRightHandSideBoundedVector79;
const double crRightHandSideBoundedVector81 =             DN_DX_0_1*v(0,1);
const double crRightHandSideBoundedVector82 =             DN_DX_1_1*v(1,1);
const double crRightHandSideBoundedVector83 =             DN_DX_2_1*v(2,1);
const double crRightHandSideBoundedVector84 =             DN_DX_3_1*v(3,1);
const double crRightHandSideBoundedVector85 =             crRightHandSideBoundedVector81 + crRightHandSideBoundedVector82 + crRightHandSideBoundedVector83 + crRightHandSideBoundedVector84;
const double crRightHandSideBoundedVector86 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector87 =             -0.1381966*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector88 =             DN_DX_0_2*v(0,1) + DN_DX_1_2*v(1,1) + DN_DX_2_2*v(2,1) + DN_DX_3_2*v(3,1);
const double crRightHandSideBoundedVector89 =             crRightHandSideBoundedVector25*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             -0.1381966*crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector91 =             0.19999999899376*f(2,1);
const double crRightHandSideBoundedVector92 =             0.19999999899376*f(3,1);
const double crRightHandSideBoundedVector93 =             4*crRightHandSideBoundedVector85*nu;
const double crRightHandSideBoundedVector94 =             crRightHandSideBoundedVector44*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector95 =             -0.1381966*crRightHandSideBoundedVector94;
const double crRightHandSideBoundedVector96 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector97 =             -0.1381966*crRightHandSideBoundedVector96;
const double crRightHandSideBoundedVector98 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector99 =             crRightHandSideBoundedVector54*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector100 =             -0.1381966*crRightHandSideBoundedVector99;
const double crRightHandSideBoundedVector101 =             crRightHandSideBoundedVector57*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector102 =             -0.1381966*crRightHandSideBoundedVector101;
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector60*crRightHandSideBoundedVector85;
const double crRightHandSideBoundedVector104 =             crRightHandSideBoundedVector64*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector105 =             -0.1381966*crRightHandSideBoundedVector104;
const double crRightHandSideBoundedVector106 =             crRightHandSideBoundedVector67*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector107 =             -0.1381966*crRightHandSideBoundedVector106;
const double crRightHandSideBoundedVector108 =             crRightHandSideBoundedVector70*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector109 =             4*DN_DX_0_0*nu;
const double crRightHandSideBoundedVector110 =             DN_DX_0_1*v(0,2) + DN_DX_1_1*v(1,2) + DN_DX_2_1*v(2,2) + DN_DX_3_1*v(3,2);
const double crRightHandSideBoundedVector111 =             crRightHandSideBoundedVector110 + crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector112 =             0.19999999899376*f(1,2);
const double crRightHandSideBoundedVector113 =             crRightHandSideBoundedVector5*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector114 =             -0.1381966*crRightHandSideBoundedVector113;
const double crRightHandSideBoundedVector115 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector116 =             -0.1381966*crRightHandSideBoundedVector115;
const double crRightHandSideBoundedVector117 =             DN_DX_0_2*v(0,2) + DN_DX_1_2*v(1,2) + DN_DX_2_2*v(2,2) + DN_DX_3_2*v(3,2);
const double crRightHandSideBoundedVector118 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector119 =             -0.1381966*crRightHandSideBoundedVector118;
const double crRightHandSideBoundedVector120 =             0.19999999899376*f(2,2);
const double crRightHandSideBoundedVector121 =             0.19999999899376*f(3,2);
const double crRightHandSideBoundedVector122 =             4*crRightHandSideBoundedVector117*nu;
const double crRightHandSideBoundedVector123 =             crRightHandSideBoundedVector44*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector124 =             -0.1381966*crRightHandSideBoundedVector123;
const double crRightHandSideBoundedVector125 =             crRightHandSideBoundedVector47*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector126 =             -0.1381966*crRightHandSideBoundedVector125;
const double crRightHandSideBoundedVector127 =             crRightHandSideBoundedVector50*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector128 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector54;
const double crRightHandSideBoundedVector129 =             -0.1381966*crRightHandSideBoundedVector128;
const double crRightHandSideBoundedVector130 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector57;
const double crRightHandSideBoundedVector131 =             -0.1381966*crRightHandSideBoundedVector130;
const double crRightHandSideBoundedVector132 =             crRightHandSideBoundedVector110*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector133 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector134 =             -0.1381966*crRightHandSideBoundedVector133;
const double crRightHandSideBoundedVector135 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector67;
const double crRightHandSideBoundedVector136 =             -0.1381966*crRightHandSideBoundedVector135;
const double crRightHandSideBoundedVector137 =             crRightHandSideBoundedVector117*crRightHandSideBoundedVector70;
const double crRightHandSideBoundedVector138 =             -1.0*crRightHandSideBoundedVector117 - 1.0*crRightHandSideBoundedVector6 - 1.0*crRightHandSideBoundedVector7 - 1.0*crRightHandSideBoundedVector8 - 1.0*crRightHandSideBoundedVector81 - 1.0*crRightHandSideBoundedVector82 - 1.0*crRightHandSideBoundedVector83 - 1.0*crRightHandSideBoundedVector84 - 1.0*crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector139 =             0.19999999899376*f(0,0);
const double crRightHandSideBoundedVector140 =             -0.1381966*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector141 =             -0.1381966*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector142 =             -0.1381966*crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector143 =             4*DN_DX_1_1*nu;
const double crRightHandSideBoundedVector144 =             4*DN_DX_1_2*nu;
const double crRightHandSideBoundedVector145 =             0.19999999899376*f(0,1);
const double crRightHandSideBoundedVector146 =             -0.1381966*crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector147 =             -0.1381966*crRightHandSideBoundedVector103;
const double crRightHandSideBoundedVector148 =             -0.1381966*crRightHandSideBoundedVector108;
const double crRightHandSideBoundedVector149 =             4*DN_DX_1_0*nu;
const double crRightHandSideBoundedVector150 =             0.19999999899376*f(0,2);
const double crRightHandSideBoundedVector151 =             -0.1381966*crRightHandSideBoundedVector127;
const double crRightHandSideBoundedVector152 =             -0.1381966*crRightHandSideBoundedVector132;
const double crRightHandSideBoundedVector153 =             -0.1381966*crRightHandSideBoundedVector137;
const double crRightHandSideBoundedVector154 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector139 + crRightHandSideBoundedVector140 + crRightHandSideBoundedVector141 + crRightHandSideBoundedVector142 + crRightHandSideBoundedVector20 + crRightHandSideBoundedVector28;
const double crRightHandSideBoundedVector155 =             4*DN_DX_2_1*nu;
const double crRightHandSideBoundedVector156 =             4*DN_DX_2_2*nu;
const double crRightHandSideBoundedVector157 =             crRightHandSideBoundedVector145 + crRightHandSideBoundedVector146 + crRightHandSideBoundedVector147 + crRightHandSideBoundedVector148 + crRightHandSideBoundedVector78 + crRightHandSideBoundedVector80 + crRightHandSideBoundedVector87 + crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector158 =             4*DN_DX_2_0*nu;
const double crRightHandSideBoundedVector159 =             crRightHandSideBoundedVector112 + crRightHandSideBoundedVector114 + crRightHandSideBoundedVector116 + crRightHandSideBoundedVector119 + crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151 + crRightHandSideBoundedVector152 + crRightHandSideBoundedVector153;
const double crRightHandSideBoundedVector160 =             4*DN_DX_3_1*nu;
const double crRightHandSideBoundedVector161 =             4*DN_DX_3_2*nu;
const double crRightHandSideBoundedVector162 =             4*DN_DX_3_0*nu;
            rRightHandSideBoundedVector[0]=DN_DX_0_0*crRightHandSideBoundedVector35 + DN_DX_0_0*crRightHandSideBoundedVector37 + DN_DX_0_0*crRightHandSideBoundedVector39 + DN_DX_0_0*crRightHandSideBoundedVector40 - DN_DX_0_0*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector20 + crRightHandSideBoundedVector28 + crRightHandSideBoundedVector29 + crRightHandSideBoundedVector30 + crRightHandSideBoundedVector46 + crRightHandSideBoundedVector49 - 0.5854102*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector56 + crRightHandSideBoundedVector59 - 0.5854102*crRightHandSideBoundedVector61 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector69 - 0.5854102*crRightHandSideBoundedVector71 - crRightHandSideBoundedVector72*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector75*crRightHandSideBoundedVector77 + 0.40000000301872*f(0,0);
            rRightHandSideBoundedVector[1]=DN_DX_0_1*crRightHandSideBoundedVector35 + DN_DX_0_1*crRightHandSideBoundedVector37 + DN_DX_0_1*crRightHandSideBoundedVector39 + DN_DX_0_1*crRightHandSideBoundedVector40 - DN_DX_0_1*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector100 + crRightHandSideBoundedVector102 - 0.5854102*crRightHandSideBoundedVector103 + crRightHandSideBoundedVector105 + crRightHandSideBoundedVector107 - 0.5854102*crRightHandSideBoundedVector108 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector75 + crRightHandSideBoundedVector78 + crRightHandSideBoundedVector80 + crRightHandSideBoundedVector87 + crRightHandSideBoundedVector90 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92 + crRightHandSideBoundedVector95 + crRightHandSideBoundedVector97 - 0.5854102*crRightHandSideBoundedVector98 + 0.40000000301872*f(0,1);
            rRightHandSideBoundedVector[2]=-DN_DX_0_2*crRightHandSideBoundedVector122 + DN_DX_0_2*crRightHandSideBoundedVector35 + DN_DX_0_2*crRightHandSideBoundedVector37 + DN_DX_0_2*crRightHandSideBoundedVector39 + DN_DX_0_2*crRightHandSideBoundedVector40 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector77 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector112 + crRightHandSideBoundedVector114 + crRightHandSideBoundedVector116 + crRightHandSideBoundedVector119 + crRightHandSideBoundedVector120 + crRightHandSideBoundedVector121 + crRightHandSideBoundedVector124 + crRightHandSideBoundedVector126 - 0.5854102*crRightHandSideBoundedVector127 + crRightHandSideBoundedVector129 + crRightHandSideBoundedVector131 - 0.5854102*crRightHandSideBoundedVector132 + crRightHandSideBoundedVector134 + crRightHandSideBoundedVector136 - 0.5854102*crRightHandSideBoundedVector137 + 0.40000000301872*f(0,2);
            rRightHandSideBoundedVector[3]=crRightHandSideBoundedVector138;
            rRightHandSideBoundedVector[4]=DN_DX_1_0*crRightHandSideBoundedVector35 + DN_DX_1_0*crRightHandSideBoundedVector37 + DN_DX_1_0*crRightHandSideBoundedVector39 + DN_DX_1_0*crRightHandSideBoundedVector40 - DN_DX_1_0*crRightHandSideBoundedVector41 - 0.5854102*crRightHandSideBoundedVector11 + crRightHandSideBoundedVector139 + crRightHandSideBoundedVector140 + crRightHandSideBoundedVector141 + crRightHandSideBoundedVector142 - crRightHandSideBoundedVector143*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector144*crRightHandSideBoundedVector77 - 0.5854102*crRightHandSideBoundedVector19 - 0.5854102*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector29 + crRightHandSideBoundedVector30 + crRightHandSideBoundedVector46 + crRightHandSideBoundedVector49 + crRightHandSideBoundedVector56 + crRightHandSideBoundedVector59 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector69 + 0.40000000301872*f(1,0);
            rRightHandSideBoundedVector[5]=DN_DX_1_1*crRightHandSideBoundedVector35 + DN_DX_1_1*crRightHandSideBoundedVector37 + DN_DX_1_1*crRightHandSideBoundedVector39 + DN_DX_1_1*crRightHandSideBoundedVector40 - DN_DX_1_1*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector100 + crRightHandSideBoundedVector102 + crRightHandSideBoundedVector105 + crRightHandSideBoundedVector107 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector144 + crRightHandSideBoundedVector145 + crRightHandSideBoundedVector146 + crRightHandSideBoundedVector147 + crRightHandSideBoundedVector148 - crRightHandSideBoundedVector149*crRightHandSideBoundedVector74 - 0.5854102*crRightHandSideBoundedVector79 - 0.5854102*crRightHandSideBoundedVector86 - 0.5854102*crRightHandSideBoundedVector89 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92 + crRightHandSideBoundedVector95 + crRightHandSideBoundedVector97 + 0.40000000301872*f(1,1);
            rRightHandSideBoundedVector[6]=-DN_DX_1_2*crRightHandSideBoundedVector122 + DN_DX_1_2*crRightHandSideBoundedVector35 + DN_DX_1_2*crRightHandSideBoundedVector37 + DN_DX_1_2*crRightHandSideBoundedVector39 + DN_DX_1_2*crRightHandSideBoundedVector40 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector143 - 0.5854102*crRightHandSideBoundedVector113 - 0.5854102*crRightHandSideBoundedVector115 - 0.5854102*crRightHandSideBoundedVector118 + crRightHandSideBoundedVector120 + crRightHandSideBoundedVector121 + crRightHandSideBoundedVector124 + crRightHandSideBoundedVector126 + crRightHandSideBoundedVector129 + crRightHandSideBoundedVector131 + crRightHandSideBoundedVector134 + crRightHandSideBoundedVector136 - crRightHandSideBoundedVector149*crRightHandSideBoundedVector77 + crRightHandSideBoundedVector150 + crRightHandSideBoundedVector151 + crRightHandSideBoundedVector152 + crRightHandSideBoundedVector153 + 0.40000000301872*f(1,2);
            rRightHandSideBoundedVector[7]=crRightHandSideBoundedVector138;
            rRightHandSideBoundedVector[8]=DN_DX_2_0*crRightHandSideBoundedVector35 + DN_DX_2_0*crRightHandSideBoundedVector37 + DN_DX_2_0*crRightHandSideBoundedVector39 + DN_DX_2_0*crRightHandSideBoundedVector40 - DN_DX_2_0*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector154 - crRightHandSideBoundedVector155*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector156*crRightHandSideBoundedVector77 + crRightHandSideBoundedVector46 - 0.5854102*crRightHandSideBoundedVector48 + crRightHandSideBoundedVector56 - 0.5854102*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector66 - 0.5854102*crRightHandSideBoundedVector68 + 0.40000000301872*f(2,0) + 0.19999999899376*f(3,0);
            rRightHandSideBoundedVector[9]=DN_DX_2_1*crRightHandSideBoundedVector35 + DN_DX_2_1*crRightHandSideBoundedVector37 + DN_DX_2_1*crRightHandSideBoundedVector39 + DN_DX_2_1*crRightHandSideBoundedVector40 - DN_DX_2_1*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector100 - 0.5854102*crRightHandSideBoundedVector101 + crRightHandSideBoundedVector105 - 0.5854102*crRightHandSideBoundedVector106 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector156 + crRightHandSideBoundedVector157 - crRightHandSideBoundedVector158*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector95 - 0.5854102*crRightHandSideBoundedVector96 + 0.40000000301872*f(2,1) + 0.19999999899376*f(3,1);
            rRightHandSideBoundedVector[10]=-DN_DX_2_2*crRightHandSideBoundedVector122 + DN_DX_2_2*crRightHandSideBoundedVector35 + DN_DX_2_2*crRightHandSideBoundedVector37 + DN_DX_2_2*crRightHandSideBoundedVector39 + DN_DX_2_2*crRightHandSideBoundedVector40 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector155 + crRightHandSideBoundedVector124 - 0.5854102*crRightHandSideBoundedVector125 + crRightHandSideBoundedVector129 - 0.5854102*crRightHandSideBoundedVector130 + crRightHandSideBoundedVector134 - 0.5854102*crRightHandSideBoundedVector135 - crRightHandSideBoundedVector158*crRightHandSideBoundedVector77 + crRightHandSideBoundedVector159 + 0.40000000301872*f(2,2) + 0.19999999899376*f(3,2);
            rRightHandSideBoundedVector[11]=crRightHandSideBoundedVector138;
            rRightHandSideBoundedVector[12]=DN_DX_3_0*crRightHandSideBoundedVector35 + DN_DX_3_0*crRightHandSideBoundedVector37 + DN_DX_3_0*crRightHandSideBoundedVector39 + DN_DX_3_0*crRightHandSideBoundedVector40 - DN_DX_3_0*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector154 - crRightHandSideBoundedVector160*crRightHandSideBoundedVector74 - crRightHandSideBoundedVector161*crRightHandSideBoundedVector77 - 0.5854102*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector49 - 0.5854102*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector59 - 0.5854102*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector69 + 0.19999999899376*f(2,0) + 0.40000000301872*f(3,0);
            rRightHandSideBoundedVector[13]=DN_DX_3_1*crRightHandSideBoundedVector35 + DN_DX_3_1*crRightHandSideBoundedVector37 + DN_DX_3_1*crRightHandSideBoundedVector39 + DN_DX_3_1*crRightHandSideBoundedVector40 - DN_DX_3_1*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector102 - 0.5854102*crRightHandSideBoundedVector104 + crRightHandSideBoundedVector107 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector161 + crRightHandSideBoundedVector157 - crRightHandSideBoundedVector162*crRightHandSideBoundedVector74 - 0.5854102*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector97 - 0.5854102*crRightHandSideBoundedVector99 + 0.19999999899376*f(2,1) + 0.40000000301872*f(3,1);
            rRightHandSideBoundedVector[14]=-DN_DX_3_2*crRightHandSideBoundedVector122 + DN_DX_3_2*crRightHandSideBoundedVector35 + DN_DX_3_2*crRightHandSideBoundedVector37 + DN_DX_3_2*crRightHandSideBoundedVector39 + DN_DX_3_2*crRightHandSideBoundedVector40 - crRightHandSideBoundedVector111*crRightHandSideBoundedVector160 - 0.5854102*crRightHandSideBoundedVector123 + crRightHandSideBoundedVector126 - 0.5854102*crRightHandSideBoundedVector128 + crRightHandSideBoundedVector131 - 0.5854102*crRightHandSideBoundedVector133 + crRightHandSideBoundedVector136 + crRightHandSideBoundedVector159 - crRightHandSideBoundedVector162*crRightHandSideBoundedVector77 + 0.19999999899376*f(2,2) + 0.40000000301872*f(3,2);
            rRightHandSideBoundedVector[15]=crRightHandSideBoundedVector138;


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2>::AddExplicitContribution(
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
void QSNavierStokesExplicit<3>::AddExplicitContribution(
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
void QSNavierStokesExplicit<2>::CalculateMassMatrix(
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
void QSNavierStokesExplicit<3>::CalculateMassMatrix(
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
int QSNavierStokesExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
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
void QSNavierStokesExplicit<TDim, TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
double QSNavierStokesExplicit<TDim, TNumNodes>::CalculateElementSize(
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

template class QSNavierStokesExplicit<2>;
template class QSNavierStokesExplicit<3>;

/***********************************************************************************/
/***********************************************************************************/

}