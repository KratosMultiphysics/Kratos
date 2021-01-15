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
#include "qs_navier_stokes_semiexplicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
QSNavierStokesSemiExplicit<TDim,TNumNodes>::QSNavierStokesSemiExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
QSNavierStokesSemiExplicit<TDim,TNumNodes>::QSNavierStokesSemiExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            // this->FractionalVelocityEquationIdVector(rResult,rCurrentProcessInfo);
            this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
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
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
        {
        case 1:
        {
            // this->GetFractionalVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        case 3:
        {
            this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        case 4:
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
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
void QSNavierStokesSemiExplicit<2,3>::CalculateLocalFractionalVelocitySystem(
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
    const auto &fracvconv = data.fractional_convective_velocity;

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

    const double crRightHandSideBoundedVector0 =             1.0/rho;
const double crRightHandSideBoundedVector1 =             0.25*f(1,0);
const double crRightHandSideBoundedVector2 =             0.166666666666667*vconv(0,0);
const double crRightHandSideBoundedVector3 =             0.166666666666667*vconv(2,0);
const double crRightHandSideBoundedVector4 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector3 + 0.666666666666667*vconv(1,0);
const double crRightHandSideBoundedVector5 =             DN_DX_0_0*v(0,0) + DN_DX_1_0*v(1,0) + DN_DX_2_0*v(2,0);
const double crRightHandSideBoundedVector6 =             0.166666666666667*vconv(0,1);
const double crRightHandSideBoundedVector7 =             0.166666666666667*vconv(2,1);
const double crRightHandSideBoundedVector8 =             crRightHandSideBoundedVector6 + crRightHandSideBoundedVector7 + 0.666666666666667*vconv(1,1);
const double crRightHandSideBoundedVector9 =             DN_DX_0_1*v(0,0) + DN_DX_1_1*v(1,0) + DN_DX_2_1*v(2,0);
const double crRightHandSideBoundedVector10 =             rho*(crRightHandSideBoundedVector4*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector8*crRightHandSideBoundedVector9);
const double crRightHandSideBoundedVector11 =             -0.166666666666667*crRightHandSideBoundedVector10;
const double crRightHandSideBoundedVector12 =             0.25*f(2,0);
const double crRightHandSideBoundedVector13 =             DN_DX_0_0*gamma;
const double crRightHandSideBoundedVector14 =             0.166666666666667*pn[0];
const double crRightHandSideBoundedVector15 =             0.166666666666667*pn[1];
const double crRightHandSideBoundedVector16 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector15 + 0.666666666666667*pn[2];
const double crRightHandSideBoundedVector17 =             0.166666666666667*pn[2];
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector17 + 0.666666666666667*pn[1];
const double crRightHandSideBoundedVector19 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector17 + 0.666666666666667*pn[0];
const double crRightHandSideBoundedVector20 =             3*crRightHandSideBoundedVector5*nu;
const double crRightHandSideBoundedVector21 =             DN_DX_0_0*v(0,1) + DN_DX_1_0*v(1,1) + DN_DX_2_0*v(2,1);
const double crRightHandSideBoundedVector22 =             3*nu*(crRightHandSideBoundedVector21 + crRightHandSideBoundedVector9);
const double crRightHandSideBoundedVector23 =             0.166666666666667*vconv(1,0);
const double crRightHandSideBoundedVector24 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector23 + 0.666666666666667*vconv(2,0);
const double crRightHandSideBoundedVector25 =             0.166666666666667*vconv(1,1);
const double crRightHandSideBoundedVector26 =             crRightHandSideBoundedVector25 + crRightHandSideBoundedVector6 + 0.666666666666667*vconv(2,1);
const double crRightHandSideBoundedVector27 =             rho*(crRightHandSideBoundedVector24*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector26*crRightHandSideBoundedVector9);
const double crRightHandSideBoundedVector28 =             -0.166666666666667*crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector29 =             crRightHandSideBoundedVector23 + crRightHandSideBoundedVector3 + 0.666666666666667*vconv(0,0);
const double crRightHandSideBoundedVector30 =             crRightHandSideBoundedVector25 + crRightHandSideBoundedVector7 + 0.666666666666667*vconv(0,1);
const double crRightHandSideBoundedVector31 =             rho*(crRightHandSideBoundedVector29*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector30*crRightHandSideBoundedVector9);
const double crRightHandSideBoundedVector32 =             0.25*f(1,1);
const double crRightHandSideBoundedVector33 =             DN_DX_0_1*v(0,1) + DN_DX_1_1*v(1,1) + DN_DX_2_1*v(2,1);
const double crRightHandSideBoundedVector34 =             rho*(crRightHandSideBoundedVector21*crRightHandSideBoundedVector4 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector8);
const double crRightHandSideBoundedVector35 =             -0.166666666666667*crRightHandSideBoundedVector34;
const double crRightHandSideBoundedVector36 =             0.25*f(2,1);
const double crRightHandSideBoundedVector37 =             DN_DX_0_1*gamma;
const double crRightHandSideBoundedVector38 =             3*crRightHandSideBoundedVector33*nu;
const double crRightHandSideBoundedVector39 =             rho*(crRightHandSideBoundedVector21*crRightHandSideBoundedVector24 + crRightHandSideBoundedVector26*crRightHandSideBoundedVector33);
const double crRightHandSideBoundedVector40 =             -0.166666666666667*crRightHandSideBoundedVector39;
const double crRightHandSideBoundedVector41 =             rho*(crRightHandSideBoundedVector21*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector30*crRightHandSideBoundedVector33);
const double crRightHandSideBoundedVector42 =             -0.166666666666667*crRightHandSideBoundedVector31 + 0.25*f(0,0);
const double crRightHandSideBoundedVector43 =             DN_DX_1_0*gamma;
const double crRightHandSideBoundedVector44 =             -0.166666666666667*crRightHandSideBoundedVector41 + 0.25*f(0,1);
const double crRightHandSideBoundedVector45 =             DN_DX_1_1*gamma;
const double crRightHandSideBoundedVector46 =             DN_DX_2_0*gamma;
const double crRightHandSideBoundedVector47 =             DN_DX_2_1*gamma;
            rRightHandSideBoundedVector[0]=crRightHandSideBoundedVector0*(-DN_DX_0_0*crRightHandSideBoundedVector20 - DN_DX_0_1*crRightHandSideBoundedVector22 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector16 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector18 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector19 + crRightHandSideBoundedVector28 - 0.666666666666667*crRightHandSideBoundedVector31 + 0.5*f(0,0));
            rRightHandSideBoundedVector[1]=crRightHandSideBoundedVector0*(-DN_DX_0_0*crRightHandSideBoundedVector22 - DN_DX_0_1*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector37 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector35 + crRightHandSideBoundedVector36 + crRightHandSideBoundedVector40 - 0.666666666666667*crRightHandSideBoundedVector41 + 0.5*f(0,1));
            rRightHandSideBoundedVector[2]=crRightHandSideBoundedVector0*(-DN_DX_1_0*crRightHandSideBoundedVector20 - DN_DX_1_1*crRightHandSideBoundedVector22 - 0.666666666666667*crRightHandSideBoundedVector10 + crRightHandSideBoundedVector12 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector43 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector43 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector43 + crRightHandSideBoundedVector28 + crRightHandSideBoundedVector42 + 0.5*f(1,0));
            rRightHandSideBoundedVector[3]=crRightHandSideBoundedVector0*(-DN_DX_1_0*crRightHandSideBoundedVector22 - DN_DX_1_1*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector45 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector45 - 0.666666666666667*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector36 + crRightHandSideBoundedVector40 + crRightHandSideBoundedVector44 + 0.5*f(1,1));
            rRightHandSideBoundedVector[4]=crRightHandSideBoundedVector0*(-DN_DX_2_0*crRightHandSideBoundedVector20 - DN_DX_2_1*crRightHandSideBoundedVector22 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector11 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector46 - 0.666666666666667*crRightHandSideBoundedVector27 + crRightHandSideBoundedVector42 + 0.5*f(2,0));
            rRightHandSideBoundedVector[5]=crRightHandSideBoundedVector0*(-DN_DX_2_0*crRightHandSideBoundedVector22 - DN_DX_2_1*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector32 + crRightHandSideBoundedVector35 - 0.666666666666667*crRightHandSideBoundedVector39 + crRightHandSideBoundedVector44 + 0.5*f(2,1));


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateLocalFractionalVelocitySystem(
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
    const auto &fracvconv = data.fractional_convective_velocity;

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

    const double crRightHandSideBoundedVector0 =             1.0/rho;
const double crRightHandSideBoundedVector1 =             0.19999999899376*f(1,0);
const double crRightHandSideBoundedVector2 =             0.1381966*vconv(0,0);
const double crRightHandSideBoundedVector3 =             0.1381966*vconv(2,0);
const double crRightHandSideBoundedVector4 =             0.1381966*vconv(3,0);
const double crRightHandSideBoundedVector5 =             crRightHandSideBoundedVector3 + crRightHandSideBoundedVector4;
const double crRightHandSideBoundedVector6 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector5 + 0.5854102*vconv(1,0);
const double crRightHandSideBoundedVector7 =             DN_DX_0_0*v(0,0) + DN_DX_1_0*v(1,0) + DN_DX_2_0*v(2,0) + DN_DX_3_0*v(3,0);
const double crRightHandSideBoundedVector8 =             0.1381966*vconv(0,1);
const double crRightHandSideBoundedVector9 =             0.1381966*vconv(2,1);
const double crRightHandSideBoundedVector10 =             0.1381966*vconv(3,1);
const double crRightHandSideBoundedVector11 =             crRightHandSideBoundedVector10 + crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector12 =             crRightHandSideBoundedVector11 + crRightHandSideBoundedVector8 + 0.5854102*vconv(1,1);
const double crRightHandSideBoundedVector13 =             DN_DX_0_1*v(0,0) + DN_DX_1_1*v(1,0) + DN_DX_2_1*v(2,0) + DN_DX_3_1*v(3,0);
const double crRightHandSideBoundedVector14 =             0.1381966*vconv(0,2);
const double crRightHandSideBoundedVector15 =             0.1381966*vconv(2,2);
const double crRightHandSideBoundedVector16 =             0.1381966*vconv(3,2);
const double crRightHandSideBoundedVector17 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector18 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector17 + 0.5854102*vconv(1,2);
const double crRightHandSideBoundedVector19 =             DN_DX_0_2*v(0,0) + DN_DX_1_2*v(1,0) + DN_DX_2_2*v(2,0) + DN_DX_3_2*v(3,0);
const double crRightHandSideBoundedVector20 =             rho*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector13 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector19 + crRightHandSideBoundedVector6*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector21 =             -0.1381966*crRightHandSideBoundedVector20;
const double crRightHandSideBoundedVector22 =             0.19999999899376*f(2,0);
const double crRightHandSideBoundedVector23 =             0.19999999899376*f(3,0);
const double crRightHandSideBoundedVector24 =             DN_DX_0_0*gamma;
const double crRightHandSideBoundedVector25 =             0.1381966*pn[0];
const double crRightHandSideBoundedVector26 =             0.1381966*pn[1];
const double crRightHandSideBoundedVector27 =             crRightHandSideBoundedVector25 + crRightHandSideBoundedVector26;
const double crRightHandSideBoundedVector28 =             0.1381966*pn[2];
const double crRightHandSideBoundedVector29 =             crRightHandSideBoundedVector27 + crRightHandSideBoundedVector28 + 0.5854102*pn[3];
const double crRightHandSideBoundedVector30 =             0.1381966*pn[3];
const double crRightHandSideBoundedVector31 =             crRightHandSideBoundedVector27 + crRightHandSideBoundedVector30 + 0.5854102*pn[2];
const double crRightHandSideBoundedVector32 =             crRightHandSideBoundedVector28 + crRightHandSideBoundedVector30;
const double crRightHandSideBoundedVector33 =             crRightHandSideBoundedVector25 + crRightHandSideBoundedVector32 + 0.5854102*pn[1];
const double crRightHandSideBoundedVector34 =             crRightHandSideBoundedVector26 + crRightHandSideBoundedVector32 + 0.5854102*pn[0];
const double crRightHandSideBoundedVector35 =             4*crRightHandSideBoundedVector7*nu;
const double crRightHandSideBoundedVector36 =             4*DN_DX_0_1*nu;
const double crRightHandSideBoundedVector37 =             DN_DX_0_0*v(0,1) + DN_DX_1_0*v(1,1) + DN_DX_2_0*v(2,1) + DN_DX_3_0*v(3,1);
const double crRightHandSideBoundedVector38 =             crRightHandSideBoundedVector13 + crRightHandSideBoundedVector37;
const double crRightHandSideBoundedVector39 =             4*DN_DX_0_2*nu;
const double crRightHandSideBoundedVector40 =             DN_DX_0_0*v(0,2) + DN_DX_1_0*v(1,2) + DN_DX_2_0*v(2,2) + DN_DX_3_0*v(3,2);
const double crRightHandSideBoundedVector41 =             crRightHandSideBoundedVector19 + crRightHandSideBoundedVector40;
const double crRightHandSideBoundedVector42 =             0.1381966*vconv(1,0);
const double crRightHandSideBoundedVector43 =             crRightHandSideBoundedVector2 + crRightHandSideBoundedVector42;
const double crRightHandSideBoundedVector44 =             crRightHandSideBoundedVector3 + crRightHandSideBoundedVector43 + 0.5854102*vconv(3,0);
const double crRightHandSideBoundedVector45 =             0.1381966*vconv(1,1);
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector45 + crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector47 =             crRightHandSideBoundedVector46 + crRightHandSideBoundedVector9 + 0.5854102*vconv(3,1);
const double crRightHandSideBoundedVector48 =             0.1381966*vconv(1,2);
const double crRightHandSideBoundedVector49 =             crRightHandSideBoundedVector14 + crRightHandSideBoundedVector48;
const double crRightHandSideBoundedVector50 =             crRightHandSideBoundedVector15 + crRightHandSideBoundedVector49 + 0.5854102*vconv(3,2);
const double crRightHandSideBoundedVector51 =             rho*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector47 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector50 + crRightHandSideBoundedVector44*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector52 =             -0.1381966*crRightHandSideBoundedVector51;
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector4 + crRightHandSideBoundedVector43 + 0.5854102*vconv(2,0);
const double crRightHandSideBoundedVector54 =             crRightHandSideBoundedVector10 + crRightHandSideBoundedVector46 + 0.5854102*vconv(2,1);
const double crRightHandSideBoundedVector55 =             crRightHandSideBoundedVector16 + crRightHandSideBoundedVector49 + 0.5854102*vconv(2,2);
const double crRightHandSideBoundedVector56 =             rho*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector54 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector55 + crRightHandSideBoundedVector53*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector57 =             -0.1381966*crRightHandSideBoundedVector56;
const double crRightHandSideBoundedVector58 =             crRightHandSideBoundedVector42 + crRightHandSideBoundedVector5 + 0.5854102*vconv(0,0);
const double crRightHandSideBoundedVector59 =             crRightHandSideBoundedVector11 + crRightHandSideBoundedVector45 + 0.5854102*vconv(0,1);
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector17 + crRightHandSideBoundedVector48 + 0.5854102*vconv(0,2);
const double crRightHandSideBoundedVector61 =             rho*(crRightHandSideBoundedVector13*crRightHandSideBoundedVector59 + crRightHandSideBoundedVector19*crRightHandSideBoundedVector60 + crRightHandSideBoundedVector58*crRightHandSideBoundedVector7);
const double crRightHandSideBoundedVector62 =             0.19999999899376*f(1,1);
const double crRightHandSideBoundedVector63 =             DN_DX_0_1*v(0,1) + DN_DX_1_1*v(1,1) + DN_DX_2_1*v(2,1) + DN_DX_3_1*v(3,1);
const double crRightHandSideBoundedVector64 =             DN_DX_0_2*v(0,1) + DN_DX_1_2*v(1,1) + DN_DX_2_2*v(2,1) + DN_DX_3_2*v(3,1);
const double crRightHandSideBoundedVector65 =             rho*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector64 + crRightHandSideBoundedVector37*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector66 =             -0.1381966*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector67 =             0.19999999899376*f(2,1);
const double crRightHandSideBoundedVector68 =             0.19999999899376*f(3,1);
const double crRightHandSideBoundedVector69 =             DN_DX_0_1*gamma;
const double crRightHandSideBoundedVector70 =             4*crRightHandSideBoundedVector63*nu;
const double crRightHandSideBoundedVector71 =             4*DN_DX_0_0*nu;
const double crRightHandSideBoundedVector72 =             DN_DX_0_1*v(0,2) + DN_DX_1_1*v(1,2) + DN_DX_2_1*v(2,2) + DN_DX_3_1*v(3,2);
const double crRightHandSideBoundedVector73 =             crRightHandSideBoundedVector64 + crRightHandSideBoundedVector72;
const double crRightHandSideBoundedVector74 =             rho*(crRightHandSideBoundedVector37*crRightHandSideBoundedVector44 + crRightHandSideBoundedVector47*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector50*crRightHandSideBoundedVector64);
const double crRightHandSideBoundedVector75 =             -0.1381966*crRightHandSideBoundedVector74;
const double crRightHandSideBoundedVector76 =             rho*(crRightHandSideBoundedVector37*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector54*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector55*crRightHandSideBoundedVector64);
const double crRightHandSideBoundedVector77 =             -0.1381966*crRightHandSideBoundedVector76;
const double crRightHandSideBoundedVector78 =             rho*(crRightHandSideBoundedVector37*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector64);
const double crRightHandSideBoundedVector79 =             0.19999999899376*f(1,2);
const double crRightHandSideBoundedVector80 =             DN_DX_0_2*v(0,2) + DN_DX_1_2*v(1,2) + DN_DX_2_2*v(2,2) + DN_DX_3_2*v(3,2);
const double crRightHandSideBoundedVector81 =             rho*(crRightHandSideBoundedVector12*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector18*crRightHandSideBoundedVector80 + crRightHandSideBoundedVector40*crRightHandSideBoundedVector6);
const double crRightHandSideBoundedVector82 =             -0.1381966*crRightHandSideBoundedVector81;
const double crRightHandSideBoundedVector83 =             0.19999999899376*f(2,2);
const double crRightHandSideBoundedVector84 =             0.19999999899376*f(3,2);
const double crRightHandSideBoundedVector85 =             DN_DX_0_2*gamma;
const double crRightHandSideBoundedVector86 =             4*crRightHandSideBoundedVector80*nu;
const double crRightHandSideBoundedVector87 =             rho*(crRightHandSideBoundedVector40*crRightHandSideBoundedVector44 + crRightHandSideBoundedVector47*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector50*crRightHandSideBoundedVector80);
const double crRightHandSideBoundedVector88 =             -0.1381966*crRightHandSideBoundedVector87;
const double crRightHandSideBoundedVector89 =             rho*(crRightHandSideBoundedVector40*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector54*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector55*crRightHandSideBoundedVector80);
const double crRightHandSideBoundedVector90 =             -0.1381966*crRightHandSideBoundedVector89;
const double crRightHandSideBoundedVector91 =             rho*(crRightHandSideBoundedVector40*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector59*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector80);
const double crRightHandSideBoundedVector92 =             0.19999999899376*f(0,0);
const double crRightHandSideBoundedVector93 =             -0.1381966*crRightHandSideBoundedVector61;
const double crRightHandSideBoundedVector94 =             DN_DX_1_0*gamma;
const double crRightHandSideBoundedVector95 =             4*DN_DX_1_1*nu;
const double crRightHandSideBoundedVector96 =             4*DN_DX_1_2*nu;
const double crRightHandSideBoundedVector97 =             0.19999999899376*f(0,1);
const double crRightHandSideBoundedVector98 =             -0.1381966*crRightHandSideBoundedVector78;
const double crRightHandSideBoundedVector99 =             DN_DX_1_1*gamma;
const double crRightHandSideBoundedVector100 =             4*DN_DX_1_0*nu;
const double crRightHandSideBoundedVector101 =             0.19999999899376*f(0,2);
const double crRightHandSideBoundedVector102 =             -0.1381966*crRightHandSideBoundedVector91;
const double crRightHandSideBoundedVector103 =             DN_DX_1_2*gamma;
const double crRightHandSideBoundedVector104 =             crRightHandSideBoundedVector1 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector92 + crRightHandSideBoundedVector93;
const double crRightHandSideBoundedVector105 =             DN_DX_2_0*gamma;
const double crRightHandSideBoundedVector106 =             4*DN_DX_2_1*nu;
const double crRightHandSideBoundedVector107 =             4*DN_DX_2_2*nu;
const double crRightHandSideBoundedVector108 =             crRightHandSideBoundedVector62 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector97 + crRightHandSideBoundedVector98;
const double crRightHandSideBoundedVector109 =             DN_DX_2_1*gamma;
const double crRightHandSideBoundedVector110 =             4*DN_DX_2_0*nu;
const double crRightHandSideBoundedVector111 =             crRightHandSideBoundedVector101 + crRightHandSideBoundedVector102 + crRightHandSideBoundedVector79 + crRightHandSideBoundedVector82;
const double crRightHandSideBoundedVector112 =             DN_DX_2_2*gamma;
const double crRightHandSideBoundedVector113 =             DN_DX_3_0*gamma;
const double crRightHandSideBoundedVector114 =             4*DN_DX_3_1*nu;
const double crRightHandSideBoundedVector115 =             4*DN_DX_3_2*nu;
const double crRightHandSideBoundedVector116 =             DN_DX_3_1*gamma;
const double crRightHandSideBoundedVector117 =             4*DN_DX_3_0*nu;
const double crRightHandSideBoundedVector118 =             DN_DX_3_2*gamma;
            rRightHandSideBoundedVector[0]=crRightHandSideBoundedVector0*(-DN_DX_0_0*crRightHandSideBoundedVector35 + crRightHandSideBoundedVector1 + crRightHandSideBoundedVector21 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector23 + crRightHandSideBoundedVector24*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector24*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector24*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector24*crRightHandSideBoundedVector34 - crRightHandSideBoundedVector36*crRightHandSideBoundedVector38 - crRightHandSideBoundedVector39*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector52 + crRightHandSideBoundedVector57 - 0.5854102*crRightHandSideBoundedVector61 + 0.40000000301872*f(0,0));
            rRightHandSideBoundedVector[1]=crRightHandSideBoundedVector0*(-DN_DX_0_1*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector31*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector34*crRightHandSideBoundedVector69 - crRightHandSideBoundedVector38*crRightHandSideBoundedVector71 - crRightHandSideBoundedVector39*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector62 + crRightHandSideBoundedVector66 + crRightHandSideBoundedVector67 + crRightHandSideBoundedVector68 + crRightHandSideBoundedVector75 + crRightHandSideBoundedVector77 - 0.5854102*crRightHandSideBoundedVector78 + 0.40000000301872*f(0,1));
            rRightHandSideBoundedVector[2]=crRightHandSideBoundedVector0*(-DN_DX_0_2*crRightHandSideBoundedVector86 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector31*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector85 + crRightHandSideBoundedVector34*crRightHandSideBoundedVector85 - crRightHandSideBoundedVector36*crRightHandSideBoundedVector73 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector79 + crRightHandSideBoundedVector82 + crRightHandSideBoundedVector83 + crRightHandSideBoundedVector84 + crRightHandSideBoundedVector88 + crRightHandSideBoundedVector90 - 0.5854102*crRightHandSideBoundedVector91 + 0.40000000301872*f(0,2));
            rRightHandSideBoundedVector[3]=crRightHandSideBoundedVector0*(-DN_DX_1_0*crRightHandSideBoundedVector35 - 0.5854102*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector22 + crRightHandSideBoundedVector23 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector31*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector94 + crRightHandSideBoundedVector34*crRightHandSideBoundedVector94 - crRightHandSideBoundedVector38*crRightHandSideBoundedVector95 - crRightHandSideBoundedVector41*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector52 + crRightHandSideBoundedVector57 + crRightHandSideBoundedVector92 + crRightHandSideBoundedVector93 + 0.40000000301872*f(1,0));
            rRightHandSideBoundedVector[4]=crRightHandSideBoundedVector0*(-DN_DX_1_1*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector100*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector29*crRightHandSideBoundedVector99 + crRightHandSideBoundedVector31*crRightHandSideBoundedVector99 + crRightHandSideBoundedVector33*crRightHandSideBoundedVector99 + crRightHandSideBoundedVector34*crRightHandSideBoundedVector99 - 0.5854102*crRightHandSideBoundedVector65 + crRightHandSideBoundedVector67 + crRightHandSideBoundedVector68 - crRightHandSideBoundedVector73*crRightHandSideBoundedVector96 + crRightHandSideBoundedVector75 + crRightHandSideBoundedVector77 + crRightHandSideBoundedVector97 + crRightHandSideBoundedVector98 + 0.40000000301872*f(1,1));
            rRightHandSideBoundedVector[5]=crRightHandSideBoundedVector0*(-DN_DX_1_2*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector100*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector101 + crRightHandSideBoundedVector102 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector103*crRightHandSideBoundedVector34 - crRightHandSideBoundedVector73*crRightHandSideBoundedVector95 - 0.5854102*crRightHandSideBoundedVector81 + crRightHandSideBoundedVector83 + crRightHandSideBoundedVector84 + crRightHandSideBoundedVector88 + crRightHandSideBoundedVector90 + 0.40000000301872*f(1,2));
            rRightHandSideBoundedVector[6]=crRightHandSideBoundedVector0*(-DN_DX_2_0*crRightHandSideBoundedVector35 + crRightHandSideBoundedVector104 + crRightHandSideBoundedVector105*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector105*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector105*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector105*crRightHandSideBoundedVector34 - crRightHandSideBoundedVector106*crRightHandSideBoundedVector38 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector23 + crRightHandSideBoundedVector52 - 0.5854102*crRightHandSideBoundedVector56 + 0.40000000301872*f(2,0));
            rRightHandSideBoundedVector[7]=crRightHandSideBoundedVector0*(-DN_DX_2_1*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector107*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector108 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector109*crRightHandSideBoundedVector34 - crRightHandSideBoundedVector110*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector68 + crRightHandSideBoundedVector75 - 0.5854102*crRightHandSideBoundedVector76 + 0.40000000301872*f(2,1));
            rRightHandSideBoundedVector[8]=crRightHandSideBoundedVector0*(-DN_DX_2_2*crRightHandSideBoundedVector86 - crRightHandSideBoundedVector106*crRightHandSideBoundedVector73 - crRightHandSideBoundedVector110*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector111 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector112*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector84 + crRightHandSideBoundedVector88 - 0.5854102*crRightHandSideBoundedVector89 + 0.40000000301872*f(2,2));
            rRightHandSideBoundedVector[9]=crRightHandSideBoundedVector0*(-DN_DX_3_0*crRightHandSideBoundedVector35 + crRightHandSideBoundedVector104 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector113*crRightHandSideBoundedVector34 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector38 - crRightHandSideBoundedVector115*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector22 - 0.5854102*crRightHandSideBoundedVector51 + crRightHandSideBoundedVector57 + 0.40000000301872*f(3,0));
            rRightHandSideBoundedVector[10]=crRightHandSideBoundedVector0*(-DN_DX_3_1*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector108 - crRightHandSideBoundedVector115*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector116*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector116*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector116*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector116*crRightHandSideBoundedVector34 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector67 - 0.5854102*crRightHandSideBoundedVector74 + crRightHandSideBoundedVector77 + 0.40000000301872*f(3,1));
            rRightHandSideBoundedVector[11]=crRightHandSideBoundedVector0*(-DN_DX_3_2*crRightHandSideBoundedVector86 + crRightHandSideBoundedVector111 - crRightHandSideBoundedVector114*crRightHandSideBoundedVector73 - crRightHandSideBoundedVector117*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector29 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector33 + crRightHandSideBoundedVector118*crRightHandSideBoundedVector34 + crRightHandSideBoundedVector83 - 0.5854102*crRightHandSideBoundedVector87 + crRightHandSideBoundedVector90 + 0.40000000301872*f(3,2));


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != n_nodes )
        rLeftHandSideMatrix.resize(n_nodes,n_nodes);
    rLeftHandSideMatrix = ZeroMatrix(n_nodes,n_nodes);
    if( rRightHandSideVector.size() != n_nodes )
        rRightHandSideVector.resize(n_nodes);
    rRightHandSideVector = ZeroVector(n_nodes);

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_convective_velocity;

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

    const double crRightHandSideVector0 =             3*DN_DX_0_0*p[0] + 3*DN_DX_1_0*p[1] + 3*DN_DX_2_0*p[2];
const double crRightHandSideVector1 =             3*DN_DX_0_1*p[0] + 3*DN_DX_1_1*p[1] + 3*DN_DX_2_1*p[2];
const double crRightHandSideVector2 =             rho*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1))/dt;
const double crRightHandSideVector3 =             3*gamma;
const double crRightHandSideVector4 =             DN_DX_0_0*pn[0] + DN_DX_1_0*pn[1] + DN_DX_2_0*pn[2];
const double crRightHandSideVector5 =             DN_DX_0_1*pn[0] + DN_DX_1_1*pn[1] + DN_DX_2_1*pn[2];
const double crRightHandSideVector6 =             -1.0*crRightHandSideVector2;
            rRightHandSideVector[0]=-DN_DX_0_0*crRightHandSideVector0 - DN_DX_0_1*crRightHandSideVector1 - 1.0*crRightHandSideVector2 + crRightHandSideVector3*(DN_DX_0_0*crRightHandSideVector4 + DN_DX_0_1*crRightHandSideVector5);
            rRightHandSideVector[1]=-DN_DX_1_0*crRightHandSideVector0 - DN_DX_1_1*crRightHandSideVector1 + crRightHandSideVector3*(DN_DX_1_0*crRightHandSideVector4 + DN_DX_1_1*crRightHandSideVector5) + crRightHandSideVector6;
            rRightHandSideVector[2]=-DN_DX_2_0*crRightHandSideVector0 - DN_DX_2_1*crRightHandSideVector1 + crRightHandSideVector3*(DN_DX_2_0*crRightHandSideVector4 + DN_DX_2_1*crRightHandSideVector5) + crRightHandSideVector6;

    const double crLeftHandSideMatrix0 =             3*(DN_DX_0_0*DN_DX_1_0 + DN_DX_0_1*DN_DX_1_1);
const double crLeftHandSideMatrix1 =             3*(DN_DX_0_0*DN_DX_2_0 + DN_DX_0_1*DN_DX_2_1);
const double crLeftHandSideMatrix2 =             3*(DN_DX_1_0*DN_DX_2_0 + DN_DX_1_1*DN_DX_2_1);
            rLeftHandSideMatrix(0,0)=3*(pow(DN_DX_0_0, 2) + pow(DN_DX_0_1, 2));
            rLeftHandSideMatrix(0,1)=crLeftHandSideMatrix0;
            rLeftHandSideMatrix(0,2)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(1,0)=crLeftHandSideMatrix0;
            rLeftHandSideMatrix(1,1)=3*(pow(DN_DX_1_0, 2) + pow(DN_DX_1_1, 2));
            rLeftHandSideMatrix(1,2)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,0)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(2,1)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(2,2)=3*(pow(DN_DX_2_0, 2) + pow(DN_DX_2_1, 2));


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != n_nodes )
        rLeftHandSideMatrix.resize(n_nodes,n_nodes);
    rLeftHandSideMatrix = ZeroMatrix(n_nodes,n_nodes);
    if( rRightHandSideVector.size() != n_nodes )
        rRightHandSideVector.resize(n_nodes);
    rRightHandSideVector = ZeroVector(n_nodes);

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_convective_velocity;

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

    const double crRightHandSideVector0 =             4*DN_DX_0_0*p[0] + 4*DN_DX_1_0*p[1] + 4*DN_DX_2_0*p[2] + 4*DN_DX_3_0*p[3];
const double crRightHandSideVector1 =             4*DN_DX_0_1*p[0] + 4*DN_DX_1_1*p[1] + 4*DN_DX_2_1*p[2] + 4*DN_DX_3_1*p[3];
const double crRightHandSideVector2 =             4*DN_DX_0_2*p[0] + 4*DN_DX_1_2*p[1] + 4*DN_DX_2_2*p[2] + 4*DN_DX_3_2*p[3];
const double crRightHandSideVector3 =             -1.0*rho*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2))/dt;
const double crRightHandSideVector4 =             4*gamma;
const double crRightHandSideVector5 =             DN_DX_0_0*pn[0] + DN_DX_1_0*pn[1] + DN_DX_2_0*pn[2] + DN_DX_3_0*pn[3];
const double crRightHandSideVector6 =             DN_DX_0_1*pn[0] + DN_DX_1_1*pn[1] + DN_DX_2_1*pn[2] + DN_DX_3_1*pn[3];
const double crRightHandSideVector7 =             DN_DX_0_2*pn[0] + DN_DX_1_2*pn[1] + DN_DX_2_2*pn[2] + DN_DX_3_2*pn[3];
            rRightHandSideVector[0]=-DN_DX_0_0*crRightHandSideVector0 - DN_DX_0_1*crRightHandSideVector1 - DN_DX_0_2*crRightHandSideVector2 + crRightHandSideVector3 + crRightHandSideVector4*(DN_DX_0_0*crRightHandSideVector5 + DN_DX_0_1*crRightHandSideVector6 + DN_DX_0_2*crRightHandSideVector7);
            rRightHandSideVector[1]=-DN_DX_1_0*crRightHandSideVector0 - DN_DX_1_1*crRightHandSideVector1 - DN_DX_1_2*crRightHandSideVector2 + crRightHandSideVector3 + crRightHandSideVector4*(DN_DX_1_0*crRightHandSideVector5 + DN_DX_1_1*crRightHandSideVector6 + DN_DX_1_2*crRightHandSideVector7);
            rRightHandSideVector[2]=-DN_DX_2_0*crRightHandSideVector0 - DN_DX_2_1*crRightHandSideVector1 - DN_DX_2_2*crRightHandSideVector2 + crRightHandSideVector3 + crRightHandSideVector4*(DN_DX_2_0*crRightHandSideVector5 + DN_DX_2_1*crRightHandSideVector6 + DN_DX_2_2*crRightHandSideVector7);
            rRightHandSideVector[3]=-DN_DX_3_0*crRightHandSideVector0 - DN_DX_3_1*crRightHandSideVector1 - DN_DX_3_2*crRightHandSideVector2 + crRightHandSideVector3 + crRightHandSideVector4*(DN_DX_3_0*crRightHandSideVector5 + DN_DX_3_1*crRightHandSideVector6 + DN_DX_3_2*crRightHandSideVector7);

    const double crLeftHandSideMatrix0 =             4*(DN_DX_0_0*DN_DX_1_0 + DN_DX_0_1*DN_DX_1_1 + DN_DX_0_2*DN_DX_1_2);
const double crLeftHandSideMatrix1 =             4*(DN_DX_0_0*DN_DX_2_0 + DN_DX_0_1*DN_DX_2_1 + DN_DX_0_2*DN_DX_2_2);
const double crLeftHandSideMatrix2 =             4*(DN_DX_0_0*DN_DX_3_0 + DN_DX_0_1*DN_DX_3_1 + DN_DX_0_2*DN_DX_3_2);
const double crLeftHandSideMatrix3 =             4*(DN_DX_1_0*DN_DX_2_0 + DN_DX_1_1*DN_DX_2_1 + DN_DX_1_2*DN_DX_2_2);
const double crLeftHandSideMatrix4 =             4*(DN_DX_1_0*DN_DX_3_0 + DN_DX_1_1*DN_DX_3_1 + DN_DX_1_2*DN_DX_3_2);
const double crLeftHandSideMatrix5 =             4*(DN_DX_2_0*DN_DX_3_0 + DN_DX_2_1*DN_DX_3_1 + DN_DX_2_2*DN_DX_3_2);
            rLeftHandSideMatrix(0,0)=4*(pow(DN_DX_0_0, 2) + pow(DN_DX_0_1, 2) + pow(DN_DX_0_2, 2));
            rLeftHandSideMatrix(0,1)=crLeftHandSideMatrix0;
            rLeftHandSideMatrix(0,2)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(0,3)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(1,0)=crLeftHandSideMatrix0;
            rLeftHandSideMatrix(1,1)=4*(pow(DN_DX_1_0, 2) + pow(DN_DX_1_1, 2) + pow(DN_DX_1_2, 2));
            rLeftHandSideMatrix(1,2)=crLeftHandSideMatrix3;
            rLeftHandSideMatrix(1,3)=crLeftHandSideMatrix4;
            rLeftHandSideMatrix(2,0)=crLeftHandSideMatrix1;
            rLeftHandSideMatrix(2,1)=crLeftHandSideMatrix3;
            rLeftHandSideMatrix(2,2)=4*(pow(DN_DX_2_0, 2) + pow(DN_DX_2_1, 2) + pow(DN_DX_2_2, 2));
            rLeftHandSideMatrix(2,3)=crLeftHandSideMatrix5;
            rLeftHandSideMatrix(3,0)=crLeftHandSideMatrix2;
            rLeftHandSideMatrix(3,1)=crLeftHandSideMatrix4;
            rLeftHandSideMatrix(3,2)=crLeftHandSideMatrix5;
            rLeftHandSideMatrix(3,3)=4*(pow(DN_DX_3_0, 2) + pow(DN_DX_3_1, 2) + pow(DN_DX_3_2, 2));


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::CalculateLocalEndOfStepSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "Calling the CalculateLocalEndOfStepSystem() method for the semiexplicit Navier-Stokes element. You should call the Calculate method instead.";

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::AddExplicitContribution(
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
void QSNavierStokesSemiExplicit<3,4>::AddExplicitContribution(
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
void QSNavierStokesSemiExplicit<2,3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 2;

    // Initialize and fill the mass matrix values
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six;    rMassMatrix(0, 2) = one_twelve; rMassMatrix(0, 4) = one_twelve;
    rMassMatrix(1, 1) = one_six;    rMassMatrix(1, 3) = one_twelve; rMassMatrix(1, 5) = one_twelve;
    rMassMatrix(2, 0) = one_twelve; rMassMatrix(2, 2) = one_six;    rMassMatrix(2, 4) = one_twelve;
    rMassMatrix(3, 1) = one_twelve; rMassMatrix(3, 3) = one_six;    rMassMatrix(3, 5) = one_twelve;
    rMassMatrix(4, 0) = one_twelve; rMassMatrix(4, 2) = one_twelve; rMassMatrix(4, 4) = one_six;
    rMassMatrix(5, 1) = one_twelve; rMassMatrix(5, 3) = one_twelve; rMassMatrix(5, 5) = one_six;

    // Here we assume that all the Gauss points have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 3;

    // Initialize and fill the mass matrix values
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_ten;     rMassMatrix(0, 3) = one_twenty;  rMassMatrix(0, 6) = one_twenty;  rMassMatrix(0,9) = one_twenty;
    rMassMatrix(1, 1) = one_ten;     rMassMatrix(1, 4) = one_twenty;  rMassMatrix(1, 7) = one_twenty;  rMassMatrix(1,10) = one_twenty;
    rMassMatrix(2, 2) = one_ten;     rMassMatrix(2, 5) = one_twenty;  rMassMatrix(2, 8) = one_twenty;  rMassMatrix(2,11) = one_twenty;
    rMassMatrix(3, 0) = one_twenty;  rMassMatrix(3, 3) = one_ten;     rMassMatrix(3, 6) = one_twenty;  rMassMatrix(3,9) = one_twenty;
    rMassMatrix(4, 1) = one_twenty;  rMassMatrix(4, 4) = one_ten;     rMassMatrix(4, 7) = one_twenty;  rMassMatrix(4,10) = one_twenty;
    rMassMatrix(5, 2) = one_twenty;  rMassMatrix(5, 5) = one_ten;     rMassMatrix(5, 8) = one_twenty;  rMassMatrix(5,11) = one_twenty;
    rMassMatrix(6, 0) = one_twenty;  rMassMatrix(6, 3) = one_twenty;  rMassMatrix(6, 6) = one_ten;     rMassMatrix(6,9) = one_twenty;
    rMassMatrix(7, 1) = one_twenty;  rMassMatrix(7, 4) = one_twenty;  rMassMatrix(7, 7) = one_ten;     rMassMatrix(7,10) = one_twenty;
    rMassMatrix(8, 2) = one_twenty;  rMassMatrix(8, 5) = one_twenty;  rMassMatrix(8, 8) = one_ten;     rMassMatrix(8,11) = one_twenty;
    rMassMatrix(9, 0) = one_twenty;  rMassMatrix(9, 3) = one_twenty;  rMassMatrix(9, 6) = one_twenty;  rMassMatrix(9,9) = one_ten;
    rMassMatrix(10, 1) = one_twenty; rMassMatrix(10, 4) = one_twenty; rMassMatrix(10, 7) = one_twenty; rMassMatrix(10,10) = one_ten;
    rMassMatrix(11, 2) = one_twenty; rMassMatrix(11, 5) = one_twenty; rMassMatrix(11, 8) = one_twenty; rMassMatrix(11,11) = one_ten;

    // Here we assume that all the Gauss points have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Define local variables
    const unsigned int local_size = 6; // two components velocity * three nodes
    const double one_third_area = GetGeometry().Area() / 3.0;
    // Initialize and calculate elemental lumped mass vector
    if (rLumpedMassVector.size() != local_size) {
        rLumpedMassVector.resize(local_size, false);
    }
    for (IndexType i = 0; i < local_size; ++i) {
        rLumpedMassVector(i) = one_third_area;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Define local variables
    const unsigned int local_size = 12; // three components velocity * four nodes
    const double one_fourth_volume = GetGeometry().Volume() / 4.0;
    // Initialize and calculate elemental lumped mass vector
    if (rLumpedMassVector.size() != local_size) {
        rLumpedMassVector.resize(local_size, false);
    }
    for (IndexType i = 0; i < local_size; ++i) {
        rLumpedMassVector(i) = one_fourth_volume;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::Calculate(
    const Variable<array_1d<double,3> > &rVariable,
    array_1d<double,3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == VELOCITY)
    {
        constexpr unsigned int nodes_per_dimension = 6;

        // Check sizes and initialize
        VectorType rRightHandSideVector;
        if( rRightHandSideVector.size() != nodes_per_dimension )
            rRightHandSideVector.resize(nodes_per_dimension);
        rRightHandSideVector = ZeroVector(nodes_per_dimension);

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Substitute the formulation symbols by the data structure values
        const auto &dt = data.dt;
        const auto &f = data.forcing;
        const auto &fracv = data.fractional_velocity;
        const auto &fracvconv = data.fractional_convective_velocity;

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

        const double crRightHandSideVector0 =             0.166666666666667*pn[0];
const double crRightHandSideVector1 =             0.166666666666667*pn[1];
const double crRightHandSideVector2 =             0.166666666666667*pn[2];
const double crRightHandSideVector3 =             dt*(-gamma*(crRightHandSideVector0 + crRightHandSideVector1 + 0.666666666666667*pn[2]) - gamma*(crRightHandSideVector0 + crRightHandSideVector2 + 0.666666666666667*pn[1]) - gamma*(crRightHandSideVector1 + crRightHandSideVector2 + 0.666666666666667*pn[0]) + 1.0*p[0] + 1.0*p[1] + 1.0*p[2])/rho;
            rRightHandSideVector[0]=DN_DX_0_0*crRightHandSideVector3;
            rRightHandSideVector[1]=DN_DX_0_1*crRightHandSideVector3;
            rRightHandSideVector[2]=DN_DX_1_0*crRightHandSideVector3;
            rRightHandSideVector[3]=DN_DX_1_1*crRightHandSideVector3;
            rRightHandSideVector[4]=DN_DX_2_0*crRightHandSideVector3;
            rRightHandSideVector[5]=DN_DX_2_1*crRightHandSideVector3;


        SizeType Index = 0;
        GeometryType& rGeom = this->GetGeometry();

        for (SizeType i = 0; i < 3; ++i) // loop over nodes
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < 2; ++d) // loop over dimensions
            {
                rTemp[d] += rRightHandSideVector[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::Calculate(
    const Variable<array_1d<double,3> > &rVariable,
    array_1d<double,3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == VELOCITY)
    {

        constexpr unsigned int nodes_per_dimension = 12;

        // Check sizes and initialize
        VectorType rRightHandSideVector;
        if( rRightHandSideVector.size() != nodes_per_dimension )
            rRightHandSideVector.resize(nodes_per_dimension);
        rRightHandSideVector = ZeroVector(nodes_per_dimension);

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Substitute the formulation symbols by the data structure values
        const auto &dt = data.dt;
        const auto &f = data.forcing;
        const auto &fracv = data.fractional_velocity;
        const auto &fracvconv = data.fractional_convective_velocity;

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

        const double crRightHandSideVector0 =             0.1381966*pn[0];
const double crRightHandSideVector1 =             0.1381966*pn[1];
const double crRightHandSideVector2 =             crRightHandSideVector0 + crRightHandSideVector1;
const double crRightHandSideVector3 =             0.1381966*pn[2];
const double crRightHandSideVector4 =             0.1381966*pn[3];
const double crRightHandSideVector5 =             crRightHandSideVector3 + crRightHandSideVector4;
const double crRightHandSideVector6 =             dt*(-gamma*(crRightHandSideVector0 + crRightHandSideVector5 + 0.5854102*pn[1]) - gamma*(crRightHandSideVector1 + crRightHandSideVector5 + 0.5854102*pn[0]) - gamma*(crRightHandSideVector2 + crRightHandSideVector3 + 0.5854102*pn[3]) - gamma*(crRightHandSideVector2 + crRightHandSideVector4 + 0.5854102*pn[2]) + 1.0*p[0] + 1.0*p[1] + 1.0*p[2] + 1.0*p[3])/rho;
            rRightHandSideVector[0]=DN_DX_0_0*crRightHandSideVector6;
            rRightHandSideVector[1]=DN_DX_0_1*crRightHandSideVector6;
            rRightHandSideVector[2]=DN_DX_0_2*crRightHandSideVector6;
            rRightHandSideVector[3]=DN_DX_1_0*crRightHandSideVector6;
            rRightHandSideVector[4]=DN_DX_1_1*crRightHandSideVector6;
            rRightHandSideVector[5]=DN_DX_1_2*crRightHandSideVector6;
            rRightHandSideVector[6]=DN_DX_2_0*crRightHandSideVector6;
            rRightHandSideVector[7]=DN_DX_2_1*crRightHandSideVector6;
            rRightHandSideVector[8]=DN_DX_2_2*crRightHandSideVector6;
            rRightHandSideVector[9]=DN_DX_3_0*crRightHandSideVector6;
            rRightHandSideVector[10]=DN_DX_3_1*crRightHandSideVector6;
            rRightHandSideVector[11]=DN_DX_3_2*crRightHandSideVector6;


        SizeType Index = 0;
        GeometryType& rGeom = this->GetGeometry();

        for (SizeType i = 0; i < 4; ++i) // loopover nodes
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < 3; ++d) // loop over dimensions
            {
                rTemp[d] += rRightHandSideVector[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
int QSNavierStokesSemiExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
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
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Getting data for the given geometry and integration method
    const auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    array_1d<double,TNumNodes> N_aux;
    GeometryUtils::CalculateGeometryData(r_geometry,rData.DN_DX,N_aux,rData.volume);
    rData.N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize some scalar data
    rData.lumping_factor = 1.00 / double(TNumNodes);
    rData.mu = 0.0;
    rData.dt = rCurrentProcessInfo[DELTA_TIME];
    // rData.dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    rData.rho = 0.0;
    // Commented code should become useful when adding stabilization and to replace theta computation
    double theta = rCurrentProcessInfo[RUNGE_KUTTA_STEP];
    if (theta == 1) {
        theta = 0.0;
    }
    else if (theta == 2 or theta == 3) {
        theta = 0.5;
    }
    else {
        theta = 1;
    }
    // double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    // if (theta == 0.0) {
    //     rData.explicit_step_coefficient = 0.0;
    // }
    // else {
    //     rData.explicit_step_coefficient = 1.0/((theta)*r_process_info[DELTA_TIME]);
    // }

    for(unsigned int node_element = 0; node_element<local_size; node_element++) {
        // Assign variables to local variables
        const auto& r_body_force = r_geometry[node_element].FastGetSolutionStepValue(BODY_FORCE);
        const auto& r_body_force_old = r_geometry[node_element].FastGetSolutionStepValue(BODY_FORCE,1);
        const auto& r_velocity = r_geometry[node_element].FastGetSolutionStepValue(VELOCITY);
        const auto& r_velocity_old = r_geometry[node_element].FastGetSolutionStepValue(VELOCITY,1);
        // Observations
        // * unknown acceleration approximated as (u-u_old)*explicit_step_coefficient = (u-u_old)/((theta)*dt)
        //   observe that for theta = 0.0, u = u_old and explicit_step_coefficient = 0

        for (unsigned int k = 0; k < TDim; ++k) {
            // forcing term: interpolation exploiting theta
            rData.forcing(node_element,k) = (1-theta) * r_body_force_old[k] + theta * r_body_force[k];
            // velocity previous time step
            rData.velocity_old(node_element,k) = r_velocity_old[k];
            // fractional velocity current time step
            // observe that VELOCITY and not FRACT_VEL is being used as DoF, to avoid saving additional DoFs
            rData.fractional_velocity(node_element,k) = r_velocity[k];
            // convective fractional velocity is fractional velocity of current explicit step
            // observe that VELOCITY and not FRACT_VEL is being used as dof, to avoid saving additional dofs
            rData.fractional_convective_velocity(node_element,k) = r_velocity[k];
        }
        // pressure current and previous time step
        rData.pressure[node_element] = r_geometry[node_element].FastGetSolutionStepValue(PRESSURE);
        rData.pressure_old[node_element] = r_geometry[node_element].FastGetSolutionStepValue(PRESSURE,1);
        // rData.oss_projection[node_element] = r_geometry[node_element].FastGetSolutionStepValue();
        rData.mu += r_geometry[node_element].FastGetSolutionStepValue(VISCOSITY);
        rData.rho += r_geometry[node_element].FastGetSolutionStepValue(DENSITY);
    }
    // divide by number of nodes scalar data
    rData.mu *= rData.lumping_factor;
    rData.rho *= rData.lumping_factor;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
double QSNavierStokesSemiExplicit<TDim,TNumNodes>::CalculateElementSize(
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
void QSNavierStokesSemiExplicit<2,3>::FractionalVelocityEquationIdVector(
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
void QSNavierStokesSemiExplicit<3,4>::FractionalVelocityEquationIdVector(
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
void QSNavierStokesSemiExplicit<2,3>::VelocityEquationIdVector(
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
void QSNavierStokesSemiExplicit<3,4>::VelocityEquationIdVector(
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
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::PressureEquationIdVector(
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
void QSNavierStokesSemiExplicit<2,3>::GetFractionalVelocityDofList(
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
void QSNavierStokesSemiExplicit<3,4>::GetFractionalVelocityDofList(
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
void QSNavierStokesSemiExplicit<2,3>::GetVelocityDofList(
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
void QSNavierStokesSemiExplicit<3,4>::GetVelocityDofList(
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
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::GetPressureDofList(DofsVectorType& rElementalDofList,
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

template class QSNavierStokesSemiExplicit<2,3>;
template class QSNavierStokesSemiExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}