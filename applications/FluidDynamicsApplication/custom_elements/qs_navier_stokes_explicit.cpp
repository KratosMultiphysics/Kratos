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

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            this->FractionalVelocityEquationIdVector(rResult,rCurrentProcessInfo);
            break;
        }
        case 3:
        {
            this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
            break;
        }
        case 4:
        {
            this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
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

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
        {
        case 1:
        {
            this->GetFractionalVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        case 5:
        {
            this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        case 6:
        {
            this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesExplicit<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // TODO: we should check and correct the shapes.

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 3:
        {
            this->CalculateLocalPressureSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
            break;
        }
        case 4:
        {
            this->CalculateLocalEndOfStepSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::CalculateLocalFractionalVelocitySystem(
    BoundedVector<double, 6> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_velocity_convective;
    const auto &fracvn = data.fractional_velocity_old;
    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

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
const double crRightHandSideBoundedVector1 =             0.166666666666667*fracvconv(0,0);
const double crRightHandSideBoundedVector2 =             0.166666666666667*fracvconv(2,0);
const double crRightHandSideBoundedVector3 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector2 + 0.666666666666667*fracvconv(1,0);
const double crRightHandSideBoundedVector4 =             DN_DX_0_0*fracv(0,0) + DN_DX_1_0*fracv(1,0) + DN_DX_2_0*fracv(2,0);
const double crRightHandSideBoundedVector5 =             0.166666666666667*fracvconv(0,1);
const double crRightHandSideBoundedVector6 =             0.166666666666667*fracvconv(2,1);
const double crRightHandSideBoundedVector7 =             crRightHandSideBoundedVector5 + crRightHandSideBoundedVector6 + 0.666666666666667*fracvconv(1,1);
const double crRightHandSideBoundedVector8 =             DN_DX_0_1*fracv(0,0) + DN_DX_1_1*fracv(1,0) + DN_DX_2_1*fracv(2,0);
const double crRightHandSideBoundedVector9 =             rho*(crRightHandSideBoundedVector3*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector7*crRightHandSideBoundedVector8);
const double crRightHandSideBoundedVector10 =             -0.166666666666667*crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             0.25*f(2,0);
const double crRightHandSideBoundedVector12 =             DN_DX_0_0*gamma;
const double crRightHandSideBoundedVector13 =             0.166666666666667*pn[0];
const double crRightHandSideBoundedVector14 =             0.166666666666667*pn[1];
const double crRightHandSideBoundedVector15 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector14 + 0.666666666666667*pn[2];
const double crRightHandSideBoundedVector16 =             0.166666666666667*pn[2];
const double crRightHandSideBoundedVector17 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector16 + 0.666666666666667*pn[1];
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector16 + 0.666666666666667*pn[0];
const double crRightHandSideBoundedVector19 =             3*crRightHandSideBoundedVector4*nu;
const double crRightHandSideBoundedVector20 =             DN_DX_0_0*fracv(0,1) + DN_DX_1_0*fracv(1,1) + DN_DX_2_0*fracv(2,1);
const double crRightHandSideBoundedVector21 =             3*nu*(crRightHandSideBoundedVector20 + crRightHandSideBoundedVector8);
const double crRightHandSideBoundedVector22 =             0.166666666666667*fracvconv(1,0);
const double crRightHandSideBoundedVector23 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector22 + 0.666666666666667*fracvconv(2,0);
const double crRightHandSideBoundedVector24 =             0.166666666666667*fracvconv(1,1);
const double crRightHandSideBoundedVector25 =             crRightHandSideBoundedVector24 + crRightHandSideBoundedVector5 + 0.666666666666667*fracvconv(2,1);
const double crRightHandSideBoundedVector26 =             rho*(crRightHandSideBoundedVector23*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector25*crRightHandSideBoundedVector8);
const double crRightHandSideBoundedVector27 =             -0.166666666666667*crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector28 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector22 + 0.666666666666667*fracvconv(0,0);
const double crRightHandSideBoundedVector29 =             crRightHandSideBoundedVector24 + crRightHandSideBoundedVector6 + 0.666666666666667*fracvconv(0,1);
const double crRightHandSideBoundedVector30 =             rho*(crRightHandSideBoundedVector28*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector8);
const double crRightHandSideBoundedVector31 =             0.25*f(1,1);
const double crRightHandSideBoundedVector32 =             DN_DX_0_1*fracv(0,1) + DN_DX_1_1*fracv(1,1) + DN_DX_2_1*fracv(2,1);
const double crRightHandSideBoundedVector33 =             rho*(crRightHandSideBoundedVector20*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector32*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector34 =             -0.166666666666667*crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector35 =             0.25*f(2,1);
const double crRightHandSideBoundedVector36 =             DN_DX_0_1*gamma;
const double crRightHandSideBoundedVector37 =             3*crRightHandSideBoundedVector32*nu;
const double crRightHandSideBoundedVector38 =             rho*(crRightHandSideBoundedVector20*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector25*crRightHandSideBoundedVector32);
const double crRightHandSideBoundedVector39 =             -0.166666666666667*crRightHandSideBoundedVector38;
const double crRightHandSideBoundedVector40 =             rho*(crRightHandSideBoundedVector20*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector32);
const double crRightHandSideBoundedVector41 =             -0.166666666666667*crRightHandSideBoundedVector30 + 0.25*f(0,0);
const double crRightHandSideBoundedVector42 =             DN_DX_1_0*gamma;
const double crRightHandSideBoundedVector43 =             -0.166666666666667*crRightHandSideBoundedVector40 + 0.25*f(0,1);
const double crRightHandSideBoundedVector44 =             DN_DX_1_1*gamma;
const double crRightHandSideBoundedVector45 =             DN_DX_2_0*gamma;
const double crRightHandSideBoundedVector46 =             DN_DX_2_1*gamma;
            rRightHandSideBoundedVector[0]=-DN_DX_0_0*crRightHandSideBoundedVector19 - DN_DX_0_1*crRightHandSideBoundedVector21 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector10 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector15 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector17 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector18 + crRightHandSideBoundedVector27 - 0.666666666666667*crRightHandSideBoundedVector30 + 0.5*f(0,0);
            rRightHandSideBoundedVector[1]=-DN_DX_0_0*crRightHandSideBoundedVector21 - DN_DX_0_1*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector36 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector36 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector36 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector34 + crRightHandSideBoundedVector35 + crRightHandSideBoundedVector39 - 0.666666666666667*crRightHandSideBoundedVector40 + 0.5*f(0,1);
            rRightHandSideBoundedVector[2]=-DN_DX_1_0*crRightHandSideBoundedVector19 - DN_DX_1_1*crRightHandSideBoundedVector21 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector42 + crRightHandSideBoundedVector27 + crRightHandSideBoundedVector41 - 0.666666666666667*crRightHandSideBoundedVector9 + 0.5*f(1,0);
            rRightHandSideBoundedVector[3]=-DN_DX_1_0*crRightHandSideBoundedVector21 - DN_DX_1_1*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector44 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector44 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector44 - 0.666666666666667*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector35 + crRightHandSideBoundedVector39 + crRightHandSideBoundedVector43 + 0.5*f(1,1);
            rRightHandSideBoundedVector[4]=-DN_DX_2_0*crRightHandSideBoundedVector19 - DN_DX_2_1*crRightHandSideBoundedVector21 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector10 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector45 - 0.666666666666667*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector41 + 0.5*f(2,0);
            rRightHandSideBoundedVector[5]=-DN_DX_2_0*crRightHandSideBoundedVector21 - DN_DX_2_1*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector15*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector31 + crRightHandSideBoundedVector34 - 0.666666666666667*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector43 + 0.5*f(2,1);


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::CalculateLocalFractionalVelocitySystem(
    BoundedVector<double, 12> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_velocity_convective;
    const auto &fracvn = data.fractional_velocity_old;
    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

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
const double crRightHandSideBoundedVector1 =             0.1381966*fracvconv(0,0);
const double crRightHandSideBoundedVector2 =             0.1381966*fracvconv(2,0);
const double crRightHandSideBoundedVector3 =             0.1381966*fracvconv(3,0);
const double crRightHandSideBoundedVector4 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3;
const double crRightHandSideBoundedVector5 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector4 + 0.5854102*fracvconv(1,0);
const double crRightHandSideBoundedVector6 =             DN_DX_0_0*fracv(0,0) + DN_DX_1_0*fracv(1,0) + DN_DX_2_0*fracv(2,0) + DN_DX_3_0*fracv(3,0);
const double crRightHandSideBoundedVector7 =             0.1381966*fracvconv(0,1);
const double crRightHandSideBoundedVector8 =             0.1381966*fracvconv(2,1);
const double crRightHandSideBoundedVector9 =             0.1381966*fracvconv(3,1);
const double crRightHandSideBoundedVector10 =             crRightHandSideBoundedVector8 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             crRightHandSideBoundedVector10 + crRightHandSideBoundedVector7 + 0.5854102*fracvconv(1,1);
const double crRightHandSideBoundedVector12 =             DN_DX_0_1*fracv(0,0) + DN_DX_1_1*fracv(1,0) + DN_DX_2_1*fracv(2,0) + DN_DX_3_1*fracv(3,0);
const double crRightHandSideBoundedVector13 =             0.1381966*fracvconv(0,2);
const double crRightHandSideBoundedVector14 =             0.1381966*fracvconv(2,2);
const double crRightHandSideBoundedVector15 =             0.1381966*fracvconv(3,2);
const double crRightHandSideBoundedVector16 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector17 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector16 + 0.5854102*fracvconv(1,2);
const double crRightHandSideBoundedVector18 =             DN_DX_0_2*fracv(0,0) + DN_DX_1_2*fracv(1,0) + DN_DX_2_2*fracv(2,0) + DN_DX_3_2*fracv(3,0);
const double crRightHandSideBoundedVector19 =             rho*(crRightHandSideBoundedVector11*crRightHandSideBoundedVector12 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector18 + crRightHandSideBoundedVector5*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector20 =             -0.1381966*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector21 =             0.19999999899376*f(2,0);
const double crRightHandSideBoundedVector22 =             0.19999999899376*f(3,0);
const double crRightHandSideBoundedVector23 =             DN_DX_0_0*gamma;
const double crRightHandSideBoundedVector24 =             0.1381966*pn[0];
const double crRightHandSideBoundedVector25 =             0.1381966*pn[1];
const double crRightHandSideBoundedVector26 =             crRightHandSideBoundedVector24 + crRightHandSideBoundedVector25;
const double crRightHandSideBoundedVector27 =             0.1381966*pn[2];
const double crRightHandSideBoundedVector28 =             crRightHandSideBoundedVector26 + crRightHandSideBoundedVector27 + 0.5854102*pn[3];
const double crRightHandSideBoundedVector29 =             0.1381966*pn[3];
const double crRightHandSideBoundedVector30 =             crRightHandSideBoundedVector26 + crRightHandSideBoundedVector29 + 0.5854102*pn[2];
const double crRightHandSideBoundedVector31 =             crRightHandSideBoundedVector27 + crRightHandSideBoundedVector29;
const double crRightHandSideBoundedVector32 =             crRightHandSideBoundedVector24 + crRightHandSideBoundedVector31 + 0.5854102*pn[1];
const double crRightHandSideBoundedVector33 =             crRightHandSideBoundedVector25 + crRightHandSideBoundedVector31 + 0.5854102*pn[0];
const double crRightHandSideBoundedVector34 =             4*crRightHandSideBoundedVector6*nu;
const double crRightHandSideBoundedVector35 =             4*DN_DX_0_1*nu;
const double crRightHandSideBoundedVector36 =             DN_DX_0_0*fracv(0,1) + DN_DX_1_0*fracv(1,1) + DN_DX_2_0*fracv(2,1) + DN_DX_3_0*fracv(3,1);
const double crRightHandSideBoundedVector37 =             crRightHandSideBoundedVector12 + crRightHandSideBoundedVector36;
const double crRightHandSideBoundedVector38 =             4*DN_DX_0_2*nu;
const double crRightHandSideBoundedVector39 =             DN_DX_0_0*fracv(0,2) + DN_DX_1_0*fracv(1,2) + DN_DX_2_0*fracv(2,2) + DN_DX_3_0*fracv(3,2);
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector18 + crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector41 =             0.1381966*fracvconv(1,0);
const double crRightHandSideBoundedVector42 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector41;
const double crRightHandSideBoundedVector43 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector42 + 0.5854102*fracvconv(3,0);
const double crRightHandSideBoundedVector44 =             0.1381966*fracvconv(1,1);
const double crRightHandSideBoundedVector45 =             crRightHandSideBoundedVector44 + crRightHandSideBoundedVector7;
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector45 + crRightHandSideBoundedVector8 + 0.5854102*fracvconv(3,1);
const double crRightHandSideBoundedVector47 =             0.1381966*fracvconv(1,2);
const double crRightHandSideBoundedVector48 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector47;
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector48 + 0.5854102*fracvconv(3,2);
const double crRightHandSideBoundedVector50 =             rho*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector43*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector51 =             -0.1381966*crRightHandSideBoundedVector50;
const double crRightHandSideBoundedVector52 =             crRightHandSideBoundedVector3 + crRightHandSideBoundedVector42 + 0.5854102*fracvconv(2,0);
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector45 + crRightHandSideBoundedVector9 + 0.5854102*fracvconv(2,1);
const double crRightHandSideBoundedVector54 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector48 + 0.5854102*fracvconv(2,2);
const double crRightHandSideBoundedVector55 =             rho*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector54 + crRightHandSideBoundedVector52*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector56 =             -0.1381966*crRightHandSideBoundedVector55;
const double crRightHandSideBoundedVector57 =             crRightHandSideBoundedVector4 + crRightHandSideBoundedVector41 + 0.5854102*fracvconv(0,0);
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector10 + crRightHandSideBoundedVector44 + 0.5854102*fracvconv(0,1);
const double crRightHandSideBoundedVector59 =             crRightHandSideBoundedVector16 + crRightHandSideBoundedVector47 + 0.5854102*fracvconv(0,2);
const double crRightHandSideBoundedVector60 =             rho*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector59 + crRightHandSideBoundedVector57*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector61 =             0.19999999899376*f(1,1);
const double crRightHandSideBoundedVector62 =             DN_DX_0_1*fracv(0,1) + DN_DX_1_1*fracv(1,1) + DN_DX_2_1*fracv(2,1) + DN_DX_3_1*fracv(3,1);
const double crRightHandSideBoundedVector63 =             DN_DX_0_2*fracv(0,1) + DN_DX_1_2*fracv(1,1) + DN_DX_2_2*fracv(2,1) + DN_DX_3_2*fracv(3,1);
const double crRightHandSideBoundedVector64 =             rho*(crRightHandSideBoundedVector11*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector36*crRightHandSideBoundedVector5);
const double crRightHandSideBoundedVector65 =             -0.1381966*crRightHandSideBoundedVector64;
const double crRightHandSideBoundedVector66 =             0.19999999899376*f(2,1);
const double crRightHandSideBoundedVector67 =             0.19999999899376*f(3,1);
const double crRightHandSideBoundedVector68 =             DN_DX_0_1*gamma;
const double crRightHandSideBoundedVector69 =             4*crRightHandSideBoundedVector62*nu;
const double crRightHandSideBoundedVector70 =             4*DN_DX_0_0*nu;
const double crRightHandSideBoundedVector71 =             DN_DX_0_1*fracv(0,2) + DN_DX_1_1*fracv(1,2) + DN_DX_2_1*fracv(2,2) + DN_DX_3_1*fracv(3,2);
const double crRightHandSideBoundedVector72 =             crRightHandSideBoundedVector63 + crRightHandSideBoundedVector71;
const double crRightHandSideBoundedVector73 =             rho*(crRightHandSideBoundedVector36*crRightHandSideBoundedVector43 + crRightHandSideBoundedVector46*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector49*crRightHandSideBoundedVector63);
const double crRightHandSideBoundedVector74 =             -0.1381966*crRightHandSideBoundedVector73;
const double crRightHandSideBoundedVector75 =             rho*(crRightHandSideBoundedVector36*crRightHandSideBoundedVector52 + crRightHandSideBoundedVector53*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector54*crRightHandSideBoundedVector63);
const double crRightHandSideBoundedVector76 =             -0.1381966*crRightHandSideBoundedVector75;
const double crRightHandSideBoundedVector77 =             rho*(crRightHandSideBoundedVector36*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector58*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector63);
const double crRightHandSideBoundedVector78 =             0.19999999899376*f(1,2);
const double crRightHandSideBoundedVector79 =             DN_DX_0_2*fracv(0,2) + DN_DX_1_2*fracv(1,2) + DN_DX_2_2*fracv(2,2) + DN_DX_3_2*fracv(3,2);
const double crRightHandSideBoundedVector80 =             rho*(crRightHandSideBoundedVector11*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector17*crRightHandSideBoundedVector79 + crRightHandSideBoundedVector39*crRightHandSideBoundedVector5);
const double crRightHandSideBoundedVector81 =             -0.1381966*crRightHandSideBoundedVector80;
const double crRightHandSideBoundedVector82 =             0.19999999899376*f(2,2);
const double crRightHandSideBoundedVector83 =             0.19999999899376*f(3,2);
const double crRightHandSideBoundedVector84 =             DN_DX_0_2*gamma;
const double crRightHandSideBoundedVector85 =             4*crRightHandSideBoundedVector79*nu;
const double crRightHandSideBoundedVector86 =             rho*(crRightHandSideBoundedVector39*crRightHandSideBoundedVector43 + crRightHandSideBoundedVector46*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector49*crRightHandSideBoundedVector79);
const double crRightHandSideBoundedVector87 =             -0.1381966*crRightHandSideBoundedVector86;
const double crRightHandSideBoundedVector88 =             rho*(crRightHandSideBoundedVector39*crRightHandSideBoundedVector52 + crRightHandSideBoundedVector53*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector54*crRightHandSideBoundedVector79);
const double crRightHandSideBoundedVector89 =             -0.1381966*crRightHandSideBoundedVector88;
const double crRightHandSideBoundedVector90 =             rho*(crRightHandSideBoundedVector39*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector58*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector79);
const double crRightHandSideBoundedVector91 =             0.19999999899376*f(0,0);
const double crRightHandSideBoundedVector92 =             -0.1381966*crRightHandSideBoundedVector60;
const double crRightHandSideBoundedVector93 =             DN_DX_1_0*gamma;
const double crRightHandSideBoundedVector94 =             4*DN_DX_1_1*nu;
const double crRightHandSideBoundedVector95 =             4*DN_DX_1_2*nu;
const double crRightHandSideBoundedVector96 =             0.19999999899376*f(0,1);
const double crRightHandSideBoundedVector97 =             -0.1381966*crRightHandSideBoundedVector77;
const double crRightHandSideBoundedVector98 =             DN_DX_1_1*gamma;
const double crRightHandSideBoundedVector99 =             4*DN_DX_1_0*nu;
const double crRightHandSideBoundedVector100 =             0.19999999899376*f(0,2);
const double crRightHandSideBoundedVector101 =             -0.1381966*crRightHandSideBoundedVector90;
const double crRightHandSideBoundedVector102 =             DN_DX_1_2*gamma;
const double crRightHandSideBoundedVector103 =             crRightHandSideBoundedVector0 + crRightHandSideBoundedVector20 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92;
const double crRightHandSideBoundedVector104 =             DN_DX_2_0*gamma;
const double crRightHandSideBoundedVector105 =             4*DN_DX_2_1*nu;
const double crRightHandSideBoundedVector106 =             4*DN_DX_2_2*nu;
const double crRightHandSideBoundedVector107 =             crRightHandSideBoundedVector61 + crRightHandSideBoundedVector65 + crRightHandSideBoundedVector96 + crRightHandSideBoundedVector97;
const double crRightHandSideBoundedVector108 =             DN_DX_2_1*gamma;
const double crRightHandSideBoundedVector109 =             4*DN_DX_2_0*nu;
const double crRightHandSideBoundedVector110 =             crRightHandSideBoundedVector100 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector78 + crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector111 =             DN_DX_2_2*gamma;
const double crRightHandSideBoundedVector112 =             DN_DX_3_0*gamma;
const double crRightHandSideBoundedVector113 =             4*DN_DX_3_1*nu;
const double crRightHandSideBoundedVector114 =             4*DN_DX_3_2*nu;
const double crRightHandSideBoundedVector115 =             DN_DX_3_1*gamma;
const double crRightHandSideBoundedVector116 =             4*DN_DX_3_0*nu;
const double crRightHandSideBoundedVector117 =             DN_DX_3_2*gamma;
            rRightHandSideBoundedVector[0]=-DN_DX_0_0*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector0 + crRightHandSideBoundedVector20 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector23*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector35*crRightHandSideBoundedVector37 - crRightHandSideBoundedVector38*crRightHandSideBoundedVector40 + crRightHandSideBoundedVector51 + crRightHandSideBoundedVector56 - 0.5854102*crRightHandSideBoundedVector60 + 0.40000000301872*f(0,0);
            rRightHandSideBoundedVector[1]=-DN_DX_0_1*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector28*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector30*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector32*crRightHandSideBoundedVector68 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector68 - crRightHandSideBoundedVector37*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector38*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector61 + crRightHandSideBoundedVector65 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector67 + crRightHandSideBoundedVector74 + crRightHandSideBoundedVector76 - 0.5854102*crRightHandSideBoundedVector77 + 0.40000000301872*f(0,1);
            rRightHandSideBoundedVector[2]=-DN_DX_0_2*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector28*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector30*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector32*crRightHandSideBoundedVector84 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector84 - crRightHandSideBoundedVector35*crRightHandSideBoundedVector72 - crRightHandSideBoundedVector40*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector78 + crRightHandSideBoundedVector81 + crRightHandSideBoundedVector82 + crRightHandSideBoundedVector83 + crRightHandSideBoundedVector87 + crRightHandSideBoundedVector89 - 0.5854102*crRightHandSideBoundedVector90 + 0.40000000301872*f(0,2);
            rRightHandSideBoundedVector[3]=-DN_DX_1_0*crRightHandSideBoundedVector34 - 0.5854102*crRightHandSideBoundedVector19 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector28*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector30*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector32*crRightHandSideBoundedVector93 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector93 - crRightHandSideBoundedVector37*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector40*crRightHandSideBoundedVector95 + crRightHandSideBoundedVector51 + crRightHandSideBoundedVector56 + crRightHandSideBoundedVector91 + crRightHandSideBoundedVector92 + 0.40000000301872*f(1,0);
            rRightHandSideBoundedVector[4]=-DN_DX_1_1*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector28*crRightHandSideBoundedVector98 + crRightHandSideBoundedVector30*crRightHandSideBoundedVector98 + crRightHandSideBoundedVector32*crRightHandSideBoundedVector98 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector98 - crRightHandSideBoundedVector37*crRightHandSideBoundedVector99 - 0.5854102*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector67 - crRightHandSideBoundedVector72*crRightHandSideBoundedVector95 + crRightHandSideBoundedVector74 + crRightHandSideBoundedVector76 + crRightHandSideBoundedVector96 + crRightHandSideBoundedVector97 + 0.40000000301872*f(1,1);
            rRightHandSideBoundedVector[5]=-DN_DX_1_2*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector100 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector102*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector102*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector102*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector102*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector40*crRightHandSideBoundedVector99 - crRightHandSideBoundedVector72*crRightHandSideBoundedVector94 - 0.5854102*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector82 + crRightHandSideBoundedVector83 + crRightHandSideBoundedVector87 + crRightHandSideBoundedVector89 + 0.40000000301872*f(1,2);
            rRightHandSideBoundedVector[6]=-DN_DX_2_0*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector103 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector104*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector105*crRightHandSideBoundedVector37 - crRightHandSideBoundedVector106*crRightHandSideBoundedVector40 + crRightHandSideBoundedVector51 - 0.5854102*crRightHandSideBoundedVector55 + 0.40000000301872*f(2,0) + 0.19999999899376*f(3,0);
            rRightHandSideBoundedVector[7]=-DN_DX_2_1*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector106*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector107 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector108*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector74 - 0.5854102*crRightHandSideBoundedVector75 + 0.40000000301872*f(2,1) + 0.19999999899376*f(3,1);
            rRightHandSideBoundedVector[8]=-DN_DX_2_2*crRightHandSideBoundedVector85 - crRightHandSideBoundedVector105*crRightHandSideBoundedVector72 - crRightHandSideBoundedVector109*crRightHandSideBoundedVector40 + crRightHandSideBoundedVector110 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector111*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector87 - 0.5854102*crRightHandSideBoundedVector88 + 0.40000000301872*f(2,2) + 0.19999999899376*f(3,2);
            rRightHandSideBoundedVector[9]=-DN_DX_3_0*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector103 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector113*crRightHandSideBoundedVector37 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector40 - 0.5854102*crRightHandSideBoundedVector50 + crRightHandSideBoundedVector56 + 0.19999999899376*f(2,0) + 0.40000000301872*f(3,0);
            rRightHandSideBoundedVector[10]=-DN_DX_3_1*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector107 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector115*crRightHandSideBoundedVector33 - crRightHandSideBoundedVector116*crRightHandSideBoundedVector37 - 0.5854102*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector76 + 0.19999999899376*f(2,1) + 0.40000000301872*f(3,1);
            rRightHandSideBoundedVector[11]=-DN_DX_3_2*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector110 - crRightHandSideBoundedVector113*crRightHandSideBoundedVector72 - crRightHandSideBoundedVector116*crRightHandSideBoundedVector40 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector28 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector117*crRightHandSideBoundedVector33 - 0.5854102*crRightHandSideBoundedVector86 + crRightHandSideBoundedVector89 + 0.19999999899376*f(2,2) + 0.40000000301872*f(3,2);


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_velocity_convective;
    const auto &fracvn = data.fractional_velocity_old;
    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

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

    const double crRightHandSideVector0 =             DN_DX_0_0*fracv(0,0);
const double crRightHandSideVector1 =             DN_DX_0_1*fracv(0,1);
const double crRightHandSideVector2 =             DN_DX_1_0*fracv(1,0);
const double crRightHandSideVector3 =             DN_DX_1_1*fracv(1,1);
const double crRightHandSideVector4 =             DN_DX_2_0*fracv(2,0);
const double crRightHandSideVector5 =             DN_DX_2_1*fracv(2,1);
const double crRightHandSideVector6 =             3*dt/rho;
const double crRightHandSideVector7 =             gamma*pn[0] - p[0];
const double crRightHandSideVector8 =             gamma*pn[1] - p[1];
const double crRightHandSideVector9 =             gamma*pn[2] - p[2];
const double crRightHandSideVector10 =             DN_DX_0_0*crRightHandSideVector7 + DN_DX_1_0*crRightHandSideVector8 + DN_DX_2_0*crRightHandSideVector9;
const double crRightHandSideVector11 =             DN_DX_0_1*crRightHandSideVector7 + DN_DX_1_1*crRightHandSideVector8 + DN_DX_2_1*crRightHandSideVector9;
const double crRightHandSideVector12 =             1.0*crRightHandSideVector0 + 1.0*crRightHandSideVector1 + 1.0*crRightHandSideVector2 + 1.0*crRightHandSideVector3 + 1.0*crRightHandSideVector4 + 1.0*crRightHandSideVector5;
            rRightHandSideVector[0]=1.0*crRightHandSideVector0 + 1.0*crRightHandSideVector1 + 1.0*crRightHandSideVector2 + 1.0*crRightHandSideVector3 + 1.0*crRightHandSideVector4 + 1.0*crRightHandSideVector5 - crRightHandSideVector6*(DN_DX_0_0*crRightHandSideVector10 + DN_DX_0_1*crRightHandSideVector11);
            rRightHandSideVector[1]=crRightHandSideVector12 - crRightHandSideVector6*(DN_DX_1_0*crRightHandSideVector10 + DN_DX_1_1*crRightHandSideVector11);
            rRightHandSideVector[2]=crRightHandSideVector12 - crRightHandSideVector6*(DN_DX_2_0*crRightHandSideVector10 + DN_DX_2_1*crRightHandSideVector11);

    const double crLeftHandSideMatrix0 =             3*dt/rho;
const double crLeftHandSideMatrix1 =             -crLeftHandSideMatrix0*(DN_DX_0_0*DN_DX_1_0 + DN_DX_0_1*DN_DX_1_1);
const double crLeftHandSideMatrix2 =             -crLeftHandSideMatrix0*(DN_DX_0_0*DN_DX_2_0 + DN_DX_0_1*DN_DX_2_1);
const double crLeftHandSideMatrix3 =             -crLeftHandSideMatrix0*(DN_DX_1_0*DN_DX_2_0 + DN_DX_1_1*DN_DX_2_1);
            rLeftHandSideMatrix(0,0)=-crLeftHandSideMatrix0*(pow(DN_DX_0_0, 2) + pow(DN_DX_0_1, 2));
            rLeftHandSideMatrix(0,1)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(0,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(1,0)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(1,1)=-crLeftHandSideMatrix0*(pow(DN_DX_1_0, 2) + pow(DN_DX_1_1, 2));
            rLeftHandSideMatrix(1,2)=crLeftHandSideMatrix3;
            rLeftHandSideMatrix(2,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,1)=crLeftHandSideMatrix3;
            rLeftHandSideMatrix(2,2)=-crLeftHandSideMatrix0*(pow(DN_DX_2_0, 2) + pow(DN_DX_2_1, 2));


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_velocity_convective;
    const auto &fracvn = data.fractional_velocity_old;
    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

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

    const double crRightHandSideVector0 =             1.0*DN_DX_0_0*fracv(0,0) + 1.0*DN_DX_0_1*fracv(0,1) + 1.0*DN_DX_0_2*fracv(0,2) + 1.0*DN_DX_1_0*fracv(1,0) + 1.0*DN_DX_1_1*fracv(1,1) + 1.0*DN_DX_1_2*fracv(1,2) + 1.0*DN_DX_2_0*fracv(2,0) + 1.0*DN_DX_2_1*fracv(2,1) + 1.0*DN_DX_2_2*fracv(2,2) + 1.0*DN_DX_3_0*fracv(3,0) + 1.0*DN_DX_3_1*fracv(3,1) + 1.0*DN_DX_3_2*fracv(3,2);
const double crRightHandSideVector1 =             4*dt/rho;
const double crRightHandSideVector2 =             gamma*pn[0] - p[0];
const double crRightHandSideVector3 =             gamma*pn[1] - p[1];
const double crRightHandSideVector4 =             gamma*pn[2] - p[2];
const double crRightHandSideVector5 =             gamma*pn[3] - p[3];
const double crRightHandSideVector6 =             DN_DX_0_0*crRightHandSideVector2 + DN_DX_1_0*crRightHandSideVector3 + DN_DX_2_0*crRightHandSideVector4 + DN_DX_3_0*crRightHandSideVector5;
const double crRightHandSideVector7 =             DN_DX_0_1*crRightHandSideVector2 + DN_DX_1_1*crRightHandSideVector3 + DN_DX_2_1*crRightHandSideVector4 + DN_DX_3_1*crRightHandSideVector5;
const double crRightHandSideVector8 =             DN_DX_0_2*crRightHandSideVector2 + DN_DX_1_2*crRightHandSideVector3 + DN_DX_2_2*crRightHandSideVector4 + DN_DX_3_2*crRightHandSideVector5;
            rRightHandSideVector[0]=crRightHandSideVector0 - crRightHandSideVector1*(DN_DX_0_0*crRightHandSideVector6 + DN_DX_0_1*crRightHandSideVector7 + DN_DX_0_2*crRightHandSideVector8);
            rRightHandSideVector[1]=crRightHandSideVector0 - crRightHandSideVector1*(DN_DX_1_0*crRightHandSideVector6 + DN_DX_1_1*crRightHandSideVector7 + DN_DX_1_2*crRightHandSideVector8);
            rRightHandSideVector[2]=crRightHandSideVector0 - crRightHandSideVector1*(DN_DX_2_0*crRightHandSideVector6 + DN_DX_2_1*crRightHandSideVector7 + DN_DX_2_2*crRightHandSideVector8);
            rRightHandSideVector[3]=crRightHandSideVector0 - crRightHandSideVector1*(DN_DX_3_0*crRightHandSideVector6 + DN_DX_3_1*crRightHandSideVector7 + DN_DX_3_2*crRightHandSideVector8);

    const double crLeftHandSideMatrix0 =             4*dt/rho;
const double crLeftHandSideMatrix1 =             -crLeftHandSideMatrix0*(DN_DX_0_0*DN_DX_1_0 + DN_DX_0_1*DN_DX_1_1 + DN_DX_0_2*DN_DX_1_2);
const double crLeftHandSideMatrix2 =             -crLeftHandSideMatrix0*(DN_DX_0_0*DN_DX_2_0 + DN_DX_0_1*DN_DX_2_1 + DN_DX_0_2*DN_DX_2_2);
const double crLeftHandSideMatrix3 =             -crLeftHandSideMatrix0*(DN_DX_0_0*DN_DX_3_0 + DN_DX_0_1*DN_DX_3_1 + DN_DX_0_2*DN_DX_3_2);
const double crLeftHandSideMatrix4 =             -crLeftHandSideMatrix0*(DN_DX_1_0*DN_DX_2_0 + DN_DX_1_1*DN_DX_2_1 + DN_DX_1_2*DN_DX_2_2);
const double crLeftHandSideMatrix5 =             -crLeftHandSideMatrix0*(DN_DX_1_0*DN_DX_3_0 + DN_DX_1_1*DN_DX_3_1 + DN_DX_1_2*DN_DX_3_2);
const double crLeftHandSideMatrix6 =             -crLeftHandSideMatrix0*(DN_DX_2_0*DN_DX_3_0 + DN_DX_2_1*DN_DX_3_1 + DN_DX_2_2*DN_DX_3_2);
            rLeftHandSideMatrix(0,0)=-crLeftHandSideMatrix0*(pow(DN_DX_0_0, 2) + pow(DN_DX_0_1, 2) + pow(DN_DX_0_2, 2));
            rLeftHandSideMatrix(0,1)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(0,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(0,3)=crLeftHandSideMatrix3;
            rLeftHandSideMatrix(1,0)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(1,1)=-crLeftHandSideMatrix0*(pow(DN_DX_1_0, 2) + pow(DN_DX_1_1, 2) + pow(DN_DX_1_2, 2));
            rLeftHandSideMatrix(1,2)=crLeftHandSideMatrix4;
            rLeftHandSideMatrix(1,3)=crLeftHandSideMatrix5;
            rLeftHandSideMatrix(2,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,1)=crLeftHandSideMatrix4;
            rLeftHandSideMatrix(2,2)=-crLeftHandSideMatrix0*(pow(DN_DX_2_0, 2) + pow(DN_DX_2_1, 2) + pow(DN_DX_2_2, 2));
            rLeftHandSideMatrix(2,3)=crLeftHandSideMatrix6;
            rLeftHandSideMatrix(3,0)=crLeftHandSideMatrix3;
            rLeftHandSideMatrix(3,1)=crLeftHandSideMatrix5;
            rLeftHandSideMatrix(3,2)=crLeftHandSideMatrix6;
            rLeftHandSideMatrix(3,3)=-crLeftHandSideMatrix0*(pow(DN_DX_3_0, 2) + pow(DN_DX_3_1, 2) + pow(DN_DX_3_2, 2));


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::CalculateLocalEndOfStepSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_velocity_convective;
    const auto &fracvn = data.fractional_velocity_old;
    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

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

    const double crRightHandSideVector0 =             0.166666666666667*p[0];
const double crRightHandSideVector1 =             0.166666666666667*p[1];
const double crRightHandSideVector2 =             0.166666666666667*pn[0];
const double crRightHandSideVector3 =             0.166666666666667*pn[1];
const double crRightHandSideVector4 =             crRightHandSideVector0 + crRightHandSideVector1 - gamma*(crRightHandSideVector2 + crRightHandSideVector3 + 0.666666666666667*pn[2]) + 0.666666666666667*p[2];
const double crRightHandSideVector5 =             0.166666666666667*p[2];
const double crRightHandSideVector6 =             0.166666666666667*pn[2];
const double crRightHandSideVector7 =             crRightHandSideVector0 + crRightHandSideVector5 - gamma*(crRightHandSideVector2 + crRightHandSideVector6 + 0.666666666666667*pn[1]) + 0.666666666666667*p[1];
const double crRightHandSideVector8 =             crRightHandSideVector1 + crRightHandSideVector5 - gamma*(crRightHandSideVector3 + crRightHandSideVector6 + 0.666666666666667*pn[0]) + 0.666666666666667*p[0];
const double crRightHandSideVector9 =             1.0/dt;
const double crRightHandSideVector10 =             0.166666666666667*fracv(1,0);
const double crRightHandSideVector11 =             -0.166666666666667*v(1,0);
const double crRightHandSideVector12 =             0.166666666666667*fracv(0,0);
const double crRightHandSideVector13 =             -0.166666666666667*v(0,0);
const double crRightHandSideVector14 =             crRightHandSideVector9*rho*(crRightHandSideVector10 + crRightHandSideVector11 + crRightHandSideVector12 + crRightHandSideVector13 + 0.666666666666667*fracv(2,0) - 0.666666666666667*v(2,0));
const double crRightHandSideVector15 =             0.166666666666667*crRightHandSideVector14;
const double crRightHandSideVector16 =             0.166666666666667*fracv(2,0) - 0.166666666666667*v(2,0);
const double crRightHandSideVector17 =             crRightHandSideVector9*rho*(crRightHandSideVector12 + crRightHandSideVector13 + crRightHandSideVector16 + 0.666666666666667*fracv(1,0) - 0.666666666666667*v(1,0));
const double crRightHandSideVector18 =             0.166666666666667*crRightHandSideVector17;
const double crRightHandSideVector19 =             crRightHandSideVector9*rho*(crRightHandSideVector10 + crRightHandSideVector11 + crRightHandSideVector16 + 0.666666666666667*fracv(0,0) - 0.666666666666667*v(0,0));
const double crRightHandSideVector20 =             0.166666666666667*fracv(1,1);
const double crRightHandSideVector21 =             -0.166666666666667*v(1,1);
const double crRightHandSideVector22 =             0.166666666666667*fracv(0,1);
const double crRightHandSideVector23 =             -0.166666666666667*v(0,1);
const double crRightHandSideVector24 =             crRightHandSideVector9*rho*(crRightHandSideVector20 + crRightHandSideVector21 + crRightHandSideVector22 + crRightHandSideVector23 + 0.666666666666667*fracv(2,1) - 0.666666666666667*v(2,1));
const double crRightHandSideVector25 =             0.166666666666667*crRightHandSideVector24;
const double crRightHandSideVector26 =             0.166666666666667*fracv(2,1) - 0.166666666666667*v(2,1);
const double crRightHandSideVector27 =             crRightHandSideVector9*rho*(crRightHandSideVector22 + crRightHandSideVector23 + crRightHandSideVector26 + 0.666666666666667*fracv(1,1) - 0.666666666666667*v(1,1));
const double crRightHandSideVector28 =             0.166666666666667*crRightHandSideVector27;
const double crRightHandSideVector29 =             crRightHandSideVector9*rho*(crRightHandSideVector20 + crRightHandSideVector21 + crRightHandSideVector26 + 0.666666666666667*fracv(0,1) - 0.666666666666667*v(0,1));
const double crRightHandSideVector30 =             0.166666666666667*crRightHandSideVector19;
const double crRightHandSideVector31 =             0.166666666666667*crRightHandSideVector29;
            rRightHandSideVector[0]=DN_DX_0_0*crRightHandSideVector4 + DN_DX_0_0*crRightHandSideVector7 + DN_DX_0_0*crRightHandSideVector8 + crRightHandSideVector15 + crRightHandSideVector18 + 0.666666666666667*crRightHandSideVector19;
            rRightHandSideVector[1]=DN_DX_0_1*crRightHandSideVector4 + DN_DX_0_1*crRightHandSideVector7 + DN_DX_0_1*crRightHandSideVector8 + crRightHandSideVector25 + crRightHandSideVector28 + 0.666666666666667*crRightHandSideVector29;
            rRightHandSideVector[2]=DN_DX_1_0*crRightHandSideVector4 + DN_DX_1_0*crRightHandSideVector7 + DN_DX_1_0*crRightHandSideVector8 + crRightHandSideVector15 + 0.666666666666667*crRightHandSideVector17 + crRightHandSideVector30;
            rRightHandSideVector[3]=DN_DX_1_1*crRightHandSideVector4 + DN_DX_1_1*crRightHandSideVector7 + DN_DX_1_1*crRightHandSideVector8 + crRightHandSideVector25 + 0.666666666666667*crRightHandSideVector27 + crRightHandSideVector31;
            rRightHandSideVector[4]=DN_DX_2_0*crRightHandSideVector4 + DN_DX_2_0*crRightHandSideVector7 + DN_DX_2_0*crRightHandSideVector8 + 0.666666666666667*crRightHandSideVector14 + crRightHandSideVector18 + crRightHandSideVector30;
            rRightHandSideVector[5]=DN_DX_2_1*crRightHandSideVector4 + DN_DX_2_1*crRightHandSideVector7 + DN_DX_2_1*crRightHandSideVector8 + 0.666666666666667*crRightHandSideVector24 + crRightHandSideVector28 + crRightHandSideVector31;

    const double crLeftHandSideMatrix0 =             rho/dt;
const double crLeftHandSideMatrix1 =             0.5*crLeftHandSideMatrix0;
const double crLeftHandSideMatrix2 =             0.25*crLeftHandSideMatrix0;
            rLeftHandSideMatrix(0,0)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(0,1)=0;
            rLeftHandSideMatrix(0,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(0,3)=0;
            rLeftHandSideMatrix(0,4)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(0,5)=0;
            rLeftHandSideMatrix(1,0)=0;
            rLeftHandSideMatrix(1,1)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(1,2)=0;
            rLeftHandSideMatrix(1,3)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(1,4)=0;
            rLeftHandSideMatrix(1,5)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,1)=0;
            rLeftHandSideMatrix(2,2)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(2,3)=0;
            rLeftHandSideMatrix(2,4)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,5)=0;
            rLeftHandSideMatrix(3,0)=0;
            rLeftHandSideMatrix(3,1)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(3,2)=0;
            rLeftHandSideMatrix(3,3)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(3,4)=0;
            rLeftHandSideMatrix(3,5)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(4,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(4,1)=0;
            rLeftHandSideMatrix(4,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(4,3)=0;
            rLeftHandSideMatrix(4,4)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(4,5)=0;
            rLeftHandSideMatrix(5,0)=0;
            rLeftHandSideMatrix(5,1)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(5,2)=0;
            rLeftHandSideMatrix(5,3)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(5,4)=0;
            rLeftHandSideMatrix(5,5)=crLeftHandSideMatrix1;


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::CalculateLocalEndOfStepSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_velocity_convective;
    const auto &fracvn = data.fractional_velocity_old;
    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

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

    const double crRightHandSideVector0 =             1.0/dt;
const double crRightHandSideVector1 =             0.1381966*fracv(2,0);
const double crRightHandSideVector2 =             -0.1381966*v(2,0);
const double crRightHandSideVector3 =             0.1381966*fracv(0,0);
const double crRightHandSideVector4 =             0.1381966*fracv(1,0);
const double crRightHandSideVector5 =             -0.1381966*v(0,0);
const double crRightHandSideVector6 =             -0.1381966*v(1,0);
const double crRightHandSideVector7 =             crRightHandSideVector0*rho*(crRightHandSideVector1 + crRightHandSideVector2 + crRightHandSideVector3 + crRightHandSideVector4 + crRightHandSideVector5 + crRightHandSideVector6 + 0.5854102*fracv(3,0) - 0.5854102*v(3,0));
const double crRightHandSideVector8 =             0.1381966*crRightHandSideVector7;
const double crRightHandSideVector9 =             0.1381966*fracv(3,0);
const double crRightHandSideVector10 =             -0.1381966*v(3,0);
const double crRightHandSideVector11 =             crRightHandSideVector0*rho*(crRightHandSideVector10 + crRightHandSideVector3 + crRightHandSideVector4 + crRightHandSideVector5 + crRightHandSideVector6 + crRightHandSideVector9 + 0.5854102*fracv(2,0) - 0.5854102*v(2,0));
const double crRightHandSideVector12 =             0.1381966*crRightHandSideVector11;
const double crRightHandSideVector13 =             crRightHandSideVector12 + crRightHandSideVector8;
const double crRightHandSideVector14 =             0.1381966*p[0];
const double crRightHandSideVector15 =             0.1381966*p[1];
const double crRightHandSideVector16 =             crRightHandSideVector14 + crRightHandSideVector15;
const double crRightHandSideVector17 =             0.1381966*p[2];
const double crRightHandSideVector18 =             0.1381966*pn[0];
const double crRightHandSideVector19 =             0.1381966*pn[1];
const double crRightHandSideVector20 =             crRightHandSideVector18 + crRightHandSideVector19;
const double crRightHandSideVector21 =             0.1381966*pn[2];
const double crRightHandSideVector22 =             crRightHandSideVector16 + crRightHandSideVector17 - gamma*(crRightHandSideVector20 + crRightHandSideVector21 + 0.5854102*pn[3]) + 0.5854102*p[3];
const double crRightHandSideVector23 =             0.1381966*p[3];
const double crRightHandSideVector24 =             0.1381966*pn[3];
const double crRightHandSideVector25 =             crRightHandSideVector16 + crRightHandSideVector23 - gamma*(crRightHandSideVector20 + crRightHandSideVector24 + 0.5854102*pn[2]) + 0.5854102*p[2];
const double crRightHandSideVector26 =             crRightHandSideVector17 + crRightHandSideVector23;
const double crRightHandSideVector27 =             crRightHandSideVector21 + crRightHandSideVector24;
const double crRightHandSideVector28 =             crRightHandSideVector14 + crRightHandSideVector26 - gamma*(crRightHandSideVector18 + crRightHandSideVector27 + 0.5854102*pn[1]) + 0.5854102*p[1];
const double crRightHandSideVector29 =             crRightHandSideVector15 + crRightHandSideVector26 - gamma*(crRightHandSideVector19 + crRightHandSideVector27 + 0.5854102*pn[0]) + 0.5854102*p[0];
const double crRightHandSideVector30 =             crRightHandSideVector1 + crRightHandSideVector10 + crRightHandSideVector2 + crRightHandSideVector9;
const double crRightHandSideVector31 =             crRightHandSideVector0*rho*(crRightHandSideVector3 + crRightHandSideVector30 + crRightHandSideVector5 + 0.5854102*fracv(1,0) - 0.5854102*v(1,0));
const double crRightHandSideVector32 =             0.1381966*crRightHandSideVector31;
const double crRightHandSideVector33 =             crRightHandSideVector0*rho*(crRightHandSideVector30 + crRightHandSideVector4 + crRightHandSideVector6 + 0.5854102*fracv(0,0) - 0.5854102*v(0,0));
const double crRightHandSideVector34 =             0.1381966*fracv(2,1);
const double crRightHandSideVector35 =             -0.1381966*v(2,1);
const double crRightHandSideVector36 =             0.1381966*fracv(0,1);
const double crRightHandSideVector37 =             0.1381966*fracv(1,1);
const double crRightHandSideVector38 =             -0.1381966*v(0,1);
const double crRightHandSideVector39 =             -0.1381966*v(1,1);
const double crRightHandSideVector40 =             crRightHandSideVector0*rho*(crRightHandSideVector34 + crRightHandSideVector35 + crRightHandSideVector36 + crRightHandSideVector37 + crRightHandSideVector38 + crRightHandSideVector39 + 0.5854102*fracv(3,1) - 0.5854102*v(3,1));
const double crRightHandSideVector41 =             0.1381966*crRightHandSideVector40;
const double crRightHandSideVector42 =             0.1381966*fracv(3,1);
const double crRightHandSideVector43 =             -0.1381966*v(3,1);
const double crRightHandSideVector44 =             crRightHandSideVector0*rho*(crRightHandSideVector36 + crRightHandSideVector37 + crRightHandSideVector38 + crRightHandSideVector39 + crRightHandSideVector42 + crRightHandSideVector43 + 0.5854102*fracv(2,1) - 0.5854102*v(2,1));
const double crRightHandSideVector45 =             0.1381966*crRightHandSideVector44;
const double crRightHandSideVector46 =             crRightHandSideVector41 + crRightHandSideVector45;
const double crRightHandSideVector47 =             crRightHandSideVector34 + crRightHandSideVector35 + crRightHandSideVector42 + crRightHandSideVector43;
const double crRightHandSideVector48 =             crRightHandSideVector0*rho*(crRightHandSideVector36 + crRightHandSideVector38 + crRightHandSideVector47 + 0.5854102*fracv(1,1) - 0.5854102*v(1,1));
const double crRightHandSideVector49 =             0.1381966*crRightHandSideVector48;
const double crRightHandSideVector50 =             crRightHandSideVector0*rho*(crRightHandSideVector37 + crRightHandSideVector39 + crRightHandSideVector47 + 0.5854102*fracv(0,1) - 0.5854102*v(0,1));
const double crRightHandSideVector51 =             0.1381966*fracv(2,2);
const double crRightHandSideVector52 =             -0.1381966*v(2,2);
const double crRightHandSideVector53 =             0.1381966*fracv(0,2);
const double crRightHandSideVector54 =             0.1381966*fracv(1,2);
const double crRightHandSideVector55 =             -0.1381966*v(0,2);
const double crRightHandSideVector56 =             -0.1381966*v(1,2);
const double crRightHandSideVector57 =             crRightHandSideVector0*rho*(crRightHandSideVector51 + crRightHandSideVector52 + crRightHandSideVector53 + crRightHandSideVector54 + crRightHandSideVector55 + crRightHandSideVector56 + 0.5854102*fracv(3,2) - 0.5854102*v(3,2));
const double crRightHandSideVector58 =             0.1381966*crRightHandSideVector57;
const double crRightHandSideVector59 =             0.1381966*fracv(3,2);
const double crRightHandSideVector60 =             -0.1381966*v(3,2);
const double crRightHandSideVector61 =             crRightHandSideVector0*rho*(crRightHandSideVector53 + crRightHandSideVector54 + crRightHandSideVector55 + crRightHandSideVector56 + crRightHandSideVector59 + crRightHandSideVector60 + 0.5854102*fracv(2,2) - 0.5854102*v(2,2));
const double crRightHandSideVector62 =             0.1381966*crRightHandSideVector61;
const double crRightHandSideVector63 =             crRightHandSideVector58 + crRightHandSideVector62;
const double crRightHandSideVector64 =             crRightHandSideVector51 + crRightHandSideVector52 + crRightHandSideVector59 + crRightHandSideVector60;
const double crRightHandSideVector65 =             crRightHandSideVector0*rho*(crRightHandSideVector53 + crRightHandSideVector55 + crRightHandSideVector64 + 0.5854102*fracv(1,2) - 0.5854102*v(1,2));
const double crRightHandSideVector66 =             0.1381966*crRightHandSideVector65;
const double crRightHandSideVector67 =             crRightHandSideVector0*rho*(crRightHandSideVector54 + crRightHandSideVector56 + crRightHandSideVector64 + 0.5854102*fracv(0,2) - 0.5854102*v(0,2));
const double crRightHandSideVector68 =             0.1381966*crRightHandSideVector33;
const double crRightHandSideVector69 =             0.1381966*crRightHandSideVector50;
const double crRightHandSideVector70 =             0.1381966*crRightHandSideVector67;
const double crRightHandSideVector71 =             crRightHandSideVector32 + crRightHandSideVector68;
const double crRightHandSideVector72 =             crRightHandSideVector49 + crRightHandSideVector69;
const double crRightHandSideVector73 =             crRightHandSideVector66 + crRightHandSideVector70;
            rRightHandSideVector[0]=DN_DX_0_0*crRightHandSideVector22 + DN_DX_0_0*crRightHandSideVector25 + DN_DX_0_0*crRightHandSideVector28 + DN_DX_0_0*crRightHandSideVector29 + crRightHandSideVector13 + crRightHandSideVector32 + 0.5854102*crRightHandSideVector33;
            rRightHandSideVector[1]=DN_DX_0_1*crRightHandSideVector22 + DN_DX_0_1*crRightHandSideVector25 + DN_DX_0_1*crRightHandSideVector28 + DN_DX_0_1*crRightHandSideVector29 + crRightHandSideVector46 + crRightHandSideVector49 + 0.5854102*crRightHandSideVector50;
            rRightHandSideVector[2]=DN_DX_0_2*crRightHandSideVector22 + DN_DX_0_2*crRightHandSideVector25 + DN_DX_0_2*crRightHandSideVector28 + DN_DX_0_2*crRightHandSideVector29 + crRightHandSideVector63 + crRightHandSideVector66 + 0.5854102*crRightHandSideVector67;
            rRightHandSideVector[3]=DN_DX_1_0*crRightHandSideVector22 + DN_DX_1_0*crRightHandSideVector25 + DN_DX_1_0*crRightHandSideVector28 + DN_DX_1_0*crRightHandSideVector29 + crRightHandSideVector13 + 0.5854102*crRightHandSideVector31 + crRightHandSideVector68;
            rRightHandSideVector[4]=DN_DX_1_1*crRightHandSideVector22 + DN_DX_1_1*crRightHandSideVector25 + DN_DX_1_1*crRightHandSideVector28 + DN_DX_1_1*crRightHandSideVector29 + crRightHandSideVector46 + 0.5854102*crRightHandSideVector48 + crRightHandSideVector69;
            rRightHandSideVector[5]=DN_DX_1_2*crRightHandSideVector22 + DN_DX_1_2*crRightHandSideVector25 + DN_DX_1_2*crRightHandSideVector28 + DN_DX_1_2*crRightHandSideVector29 + crRightHandSideVector63 + 0.5854102*crRightHandSideVector65 + crRightHandSideVector70;
            rRightHandSideVector[6]=DN_DX_2_0*crRightHandSideVector22 + DN_DX_2_0*crRightHandSideVector25 + DN_DX_2_0*crRightHandSideVector28 + DN_DX_2_0*crRightHandSideVector29 + 0.5854102*crRightHandSideVector11 + crRightHandSideVector71 + crRightHandSideVector8;
            rRightHandSideVector[7]=DN_DX_2_1*crRightHandSideVector22 + DN_DX_2_1*crRightHandSideVector25 + DN_DX_2_1*crRightHandSideVector28 + DN_DX_2_1*crRightHandSideVector29 + crRightHandSideVector41 + 0.5854102*crRightHandSideVector44 + crRightHandSideVector72;
            rRightHandSideVector[8]=DN_DX_2_2*crRightHandSideVector22 + DN_DX_2_2*crRightHandSideVector25 + DN_DX_2_2*crRightHandSideVector28 + DN_DX_2_2*crRightHandSideVector29 + crRightHandSideVector58 + 0.5854102*crRightHandSideVector61 + crRightHandSideVector73;
            rRightHandSideVector[9]=DN_DX_3_0*crRightHandSideVector22 + DN_DX_3_0*crRightHandSideVector25 + DN_DX_3_0*crRightHandSideVector28 + DN_DX_3_0*crRightHandSideVector29 + crRightHandSideVector12 + 0.5854102*crRightHandSideVector7 + crRightHandSideVector71;
            rRightHandSideVector[10]=DN_DX_3_1*crRightHandSideVector22 + DN_DX_3_1*crRightHandSideVector25 + DN_DX_3_1*crRightHandSideVector28 + DN_DX_3_1*crRightHandSideVector29 + 0.5854102*crRightHandSideVector40 + crRightHandSideVector45 + crRightHandSideVector72;
            rRightHandSideVector[11]=DN_DX_3_2*crRightHandSideVector22 + DN_DX_3_2*crRightHandSideVector25 + DN_DX_3_2*crRightHandSideVector28 + DN_DX_3_2*crRightHandSideVector29 + 0.5854102*crRightHandSideVector57 + crRightHandSideVector62 + crRightHandSideVector73;

    const double crLeftHandSideMatrix0 =             rho/dt;
const double crLeftHandSideMatrix1 =             0.40000000301872*crLeftHandSideMatrix0;
const double crLeftHandSideMatrix2 =             0.19999999899376*crLeftHandSideMatrix0;
            rLeftHandSideMatrix(0,0)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(0,1)=0;
            rLeftHandSideMatrix(0,2)=0;
            rLeftHandSideMatrix(0,3)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(0,4)=0;
            rLeftHandSideMatrix(0,5)=0;
            rLeftHandSideMatrix(0,6)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(0,7)=0;
            rLeftHandSideMatrix(0,8)=0;
            rLeftHandSideMatrix(0,9)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(0,10)=0;
            rLeftHandSideMatrix(0,11)=0;
            rLeftHandSideMatrix(1,0)=0;
            rLeftHandSideMatrix(1,1)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(1,2)=0;
            rLeftHandSideMatrix(1,3)=0;
            rLeftHandSideMatrix(1,4)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(1,5)=0;
            rLeftHandSideMatrix(1,6)=0;
            rLeftHandSideMatrix(1,7)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(1,8)=0;
            rLeftHandSideMatrix(1,9)=0;
            rLeftHandSideMatrix(1,10)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(1,11)=0;
            rLeftHandSideMatrix(2,0)=0;
            rLeftHandSideMatrix(2,1)=0;
            rLeftHandSideMatrix(2,2)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(2,3)=0;
            rLeftHandSideMatrix(2,4)=0;
            rLeftHandSideMatrix(2,5)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,6)=0;
            rLeftHandSideMatrix(2,7)=0;
            rLeftHandSideMatrix(2,8)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,9)=0;
            rLeftHandSideMatrix(2,10)=0;
            rLeftHandSideMatrix(2,11)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(3,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(3,1)=0;
            rLeftHandSideMatrix(3,2)=0;
            rLeftHandSideMatrix(3,3)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(3,4)=0;
            rLeftHandSideMatrix(3,5)=0;
            rLeftHandSideMatrix(3,6)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(3,7)=0;
            rLeftHandSideMatrix(3,8)=0;
            rLeftHandSideMatrix(3,9)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(3,10)=0;
            rLeftHandSideMatrix(3,11)=0;
            rLeftHandSideMatrix(4,0)=0;
            rLeftHandSideMatrix(4,1)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(4,2)=0;
            rLeftHandSideMatrix(4,3)=0;
            rLeftHandSideMatrix(4,4)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(4,5)=0;
            rLeftHandSideMatrix(4,6)=0;
            rLeftHandSideMatrix(4,7)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(4,8)=0;
            rLeftHandSideMatrix(4,9)=0;
            rLeftHandSideMatrix(4,10)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(4,11)=0;
            rLeftHandSideMatrix(5,0)=0;
            rLeftHandSideMatrix(5,1)=0;
            rLeftHandSideMatrix(5,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(5,3)=0;
            rLeftHandSideMatrix(5,4)=0;
            rLeftHandSideMatrix(5,5)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(5,6)=0;
            rLeftHandSideMatrix(5,7)=0;
            rLeftHandSideMatrix(5,8)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(5,9)=0;
            rLeftHandSideMatrix(5,10)=0;
            rLeftHandSideMatrix(5,11)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(6,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(6,1)=0;
            rLeftHandSideMatrix(6,2)=0;
            rLeftHandSideMatrix(6,3)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(6,4)=0;
            rLeftHandSideMatrix(6,5)=0;
            rLeftHandSideMatrix(6,6)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(6,7)=0;
            rLeftHandSideMatrix(6,8)=0;
            rLeftHandSideMatrix(6,9)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(6,10)=0;
            rLeftHandSideMatrix(6,11)=0;
            rLeftHandSideMatrix(7,0)=0;
            rLeftHandSideMatrix(7,1)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(7,2)=0;
            rLeftHandSideMatrix(7,3)=0;
            rLeftHandSideMatrix(7,4)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(7,5)=0;
            rLeftHandSideMatrix(7,6)=0;
            rLeftHandSideMatrix(7,7)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(7,8)=0;
            rLeftHandSideMatrix(7,9)=0;
            rLeftHandSideMatrix(7,10)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(7,11)=0;
            rLeftHandSideMatrix(8,0)=0;
            rLeftHandSideMatrix(8,1)=0;
            rLeftHandSideMatrix(8,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(8,3)=0;
            rLeftHandSideMatrix(8,4)=0;
            rLeftHandSideMatrix(8,5)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(8,6)=0;
            rLeftHandSideMatrix(8,7)=0;
            rLeftHandSideMatrix(8,8)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(8,9)=0;
            rLeftHandSideMatrix(8,10)=0;
            rLeftHandSideMatrix(8,11)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(9,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(9,1)=0;
            rLeftHandSideMatrix(9,2)=0;
            rLeftHandSideMatrix(9,3)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(9,4)=0;
            rLeftHandSideMatrix(9,5)=0;
            rLeftHandSideMatrix(9,6)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(9,7)=0;
            rLeftHandSideMatrix(9,8)=0;
            rLeftHandSideMatrix(9,9)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(9,10)=0;
            rLeftHandSideMatrix(9,11)=0;
            rLeftHandSideMatrix(10,0)=0;
            rLeftHandSideMatrix(10,1)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(10,2)=0;
            rLeftHandSideMatrix(10,3)=0;
            rLeftHandSideMatrix(10,4)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(10,5)=0;
            rLeftHandSideMatrix(10,6)=0;
            rLeftHandSideMatrix(10,7)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(10,8)=0;
            rLeftHandSideMatrix(10,9)=0;
            rLeftHandSideMatrix(10,10)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(10,11)=0;
            rLeftHandSideMatrix(11,0)=0;
            rLeftHandSideMatrix(11,1)=0;
            rLeftHandSideMatrix(11,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(11,3)=0;
            rLeftHandSideMatrix(11,4)=0;
            rLeftHandSideMatrix(11,5)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(11,6)=0;
            rLeftHandSideMatrix(11,7)=0;
            rLeftHandSideMatrix(11,8)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(11,9)=0;
            rLeftHandSideMatrix(11,10)=0;
            rLeftHandSideMatrix(11,11)=crLeftHandSideMatrix1;


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            constexpr IndexType dim = 2;
            constexpr IndexType n_nodes = 3;
            constexpr IndexType block_size = 2;

            // Calculate the explicit residual vector
            BoundedVector<double, 6> rhs;
            this->CalculateLocalFractionalVelocitySystem(rhs, rCurrentProcessInfo);

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
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            constexpr IndexType dim = 3;
            constexpr IndexType n_nodes = 4;
            constexpr IndexType block_size = 3;

            // Calculate the explicit residual vector
            BoundedVector<double, 12> rhs;
            this->CalculateLocalFractionalVelocitySystem(rhs, rCurrentProcessInfo);

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
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }


    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::CalculateMassMatrix(
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
void QSNavierStokesExplicit<3,4>::CalculateMassMatrix(
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
void QSNavierStokesExplicit<TDim,TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
double QSNavierStokesExplicit<TDim,TNumNodes>::CalculateElementSize(
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

template<>
void QSNavierStokesExplicit<2,3>::FractionalVelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACT_VEL_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y,xpos+1).EquationId();
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::FractionalVelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACT_VEL_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Z,xpos+2).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::VelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::VelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesExplicit<TDim,TNumNodes>::PressureEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes);

    const unsigned int pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (SizeType i = 0; i < NumNodes; ++i)
        rResult[i] = rGeom[i].GetDof(PRESSURE,pos).EquationId();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::GetFractionalVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Y);
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::GetFractionalVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Z);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesExplicit<2,3>::GetVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesExplicit<3,4>::GetVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesExplicit<TDim,TNumNodes>::GetPressureDofList(DofsVectorType& rElementalDofList,
                                              const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class QSNavierStokesExplicit<2,3>;
template class QSNavierStokesExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}