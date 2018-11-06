#if !defined(KRATOS_CELLULAR_FLOW_PARTIAL_DERIVATIVES_H)
#define KRATOS_CELLULAR_FLOW_PARTIAL_DERIVATIVES_H

/* System includes */
#include "cellular_flow_partial_derivatives.h"

namespace Kratos
{
void CellularFlowPartialDerivatives::UpdateCoordinates(const double time, const array_1d<double, 3>& coor)
{
    mSinOmegaT[i_thread] = std::sin(mOmegaUOverL * time);
    mCosOmegaT[i_thread] = std::cos(mOmegaUOverL * time);
    mSinPiX0[i_thread]   = std::sin(mOneOverL * coor[0]);
    mCosPiX0[i_thread]   = std::cos(mOneOverL * coor[0]);
    mSinPiX1[i_thread]   = std::sin(mOneOverL * coor[1]);
    mCosPiX1[i_thread]   = std::cos(mOneOverL * coor[1]);
}

// Values

double CellularFlowPartialDerivatives::U0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowPartialDerivatives::U1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowPartialDerivatives::U2(const int i){return 0.0;}

// First-order derivatives

double CellularFlowPartialDerivatives::U0DT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return   mU * mK * mOmegaUOverL * mCosOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U0D0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowPartialDerivatives::U0D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowPartialDerivatives::U0D2(const int i){return 0.0;}

double CellularFlowPartialDerivatives::U1DT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mCosPiX0[i] * mSinPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U1D0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowPartialDerivatives::U1D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowPartialDerivatives::U1D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2DT(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D0(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D1(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D2(const int i){return 0.0;}

// Second-order derivatives

double CellularFlowPartialDerivatives::U0DTDT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mOmegaUOverL * mSinOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U0DTD0(const int i)
{
    if (mOmega == 0.0){
            return 0.0;
    }
    else {
        return   mU * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U0DTD1(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U0DTD2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U0D0D0(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowPartialDerivatives::U0D0D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowPartialDerivatives::U0D0D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U0D1D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowPartialDerivatives::U0D1D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U0D2D2(const int i){return 0.0;}

double CellularFlowPartialDerivatives::U1DTDT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mOmegaUOverL * mSinOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U1DTD0(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return   mU * mK * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U1DTD1(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowPartialDerivatives::U1DTD2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U1D0D0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowPartialDerivatives::U1D0D1(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowPartialDerivatives::U1D0D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U1D1D1(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowPartialDerivatives::U1D1D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U1D2D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2DTDT(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2DTD0(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2DTD1(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2DTD2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D0D0(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D0D1(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D0D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D1D1(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D1D2(const int i){return 0.0;}
double CellularFlowPartialDerivatives::U2D2D2(const int i){return 0.0;}

//***************************************************************************************************************
//***************************************************************************************************************
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_CELLULAR_FLOW_PARTIAL_DERIVATIVES_H defined
