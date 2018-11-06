#include "cellular_flow_field.h"

namespace Kratos
{

void CellularFlowField::ResizeVectorsForParallelism(const int n_threads)
{
    mSinOmegaT.resize(n_threads);
    mCosOmegaT.resize(n_threads);
    mSinPiX0.resize(n_threads);
    mCosPiX0.resize(n_threads);
    mSinPiX1.resize(n_threads);
    mCosPiX1.resize(n_threads);
    mCoordinatesAreUpToDate.resize(n_threads);
    for (int i = 0; i < n_threads; i++){
        mCoordinatesAreUpToDate[i] = false;
    }
}
void CellularFlowField::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mSinOmegaT[i_thread] = std::sin(mOmegaUOverL * time);
        mCosOmegaT[i_thread] = std::cos(mOmegaUOverL * time);
        mSinPiX0[i_thread]   = std::sin(mOneOverL * coor[0]);
        mCosPiX0[i_thread]   = std::cos(mOneOverL * coor[0]);
        mSinPiX1[i_thread]   = std::sin(mOneOverL * coor[1]);
        mCosPiX1[i_thread]   = std::cos(mOneOverL * coor[1]);
    }
}

void CellularFlowField::UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mSinOmegaT[i_thread] = std::sin(mOmegaUOverL * time);
        mCosOmegaT[i_thread] = std::cos(mOmegaUOverL * time);
        mSinPiX0[i_thread]   = std::sin(mOneOverL * coor[0]);
        mCosPiX0[i_thread]   = std::cos(mOneOverL * coor[0]);
        mSinPiX1[i_thread]   = std::sin(mOneOverL * coor[1]);
        mCosPiX1[i_thread]   = std::cos(mOneOverL * coor[1]);
    }
}

void CellularFlowField::LockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 1;
}

void CellularFlowField::UnlockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 0;
}
// Values

double CellularFlowField::U0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U2(const int i){return 0.0;}

// First-order derivatives

double CellularFlowField::U0DT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return   mU * mK * mOmegaUOverL * mCosOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowField::U0D0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U0D2(const int i){return 0.0;}

double CellularFlowField::U1DT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mCosPiX0[i] * mSinPiX1[i];
    }
}
double CellularFlowField::U1D0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1D2(const int i){return 0.0;}
double CellularFlowField::U2DT(const int i){return 0.0;}
double CellularFlowField::U2D0(const int i){return 0.0;}
double CellularFlowField::U2D1(const int i){return 0.0;}
double CellularFlowField::U2D2(const int i){return 0.0;}

// Second-order derivatives

double CellularFlowField::U0DTDT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mOmegaUOverL * mSinOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowField::U0DTD0(const int i)
{
    if (mOmega == 0.0){
            return 0.0;
    }
    else {
        return   mU * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowField::U0DTD1(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
    }
}
double CellularFlowField::U0DTD2(const int i){return 0.0;}
double CellularFlowField::U0D0D0(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D0D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U0D0D2(const int i){return 0.0;}
double CellularFlowField::U0D1D1(const int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D1D2(const int i){return 0.0;}
double CellularFlowField::U0D2D2(const int i){return 0.0;}

double CellularFlowField::U1DTDT(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mOmegaUOverL * mSinOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowField::U1DTD0(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return   mU * mK * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mSinPiX0[i] * mSinPiX1[i];
    }
}
double CellularFlowField::U1DTD1(const int i)
{   if (mOmega == 0.0){
        return 0.0;
    }
    else {
        return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mOneOverL * mCosPiX0[i] * mCosPiX1[i];
    }
}
double CellularFlowField::U1DTD2(const int i){return 0.0;}
double CellularFlowField::U1D0D0(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D0D1(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1D0D2(const int i){return 0.0;}
double CellularFlowField::U1D1D1(const int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mOneOverL * mOneOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D1D2(const int i){return 0.0;}
double CellularFlowField::U1D2D2(const int i){return 0.0;}
double CellularFlowField::U2DTDT(const int i){return 0.0;}
double CellularFlowField::U2DTD0(const int i){return 0.0;}
double CellularFlowField::U2DTD1(const int i){return 0.0;}
double CellularFlowField::U2DTD2(const int i){return 0.0;}
double CellularFlowField::U2D0D0(const int i){return 0.0;}
double CellularFlowField::U2D0D1(const int i){return 0.0;}
double CellularFlowField::U2D0D2(const int i){return 0.0;}
double CellularFlowField::U2D1D1(const int i){return 0.0;}
double CellularFlowField::U2D1D2(const int i){return 0.0;}
double CellularFlowField::U2D2D2(const int i){return 0.0;}

} // namespace Kratos.




