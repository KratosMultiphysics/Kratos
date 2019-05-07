#include "ethier_flow_field.h"

namespace Kratos
{

void EthierVelocityField::ResizeVectorsForParallelism(const int n_threads)
{
    mExpD2T.resize(n_threads);
    mExpAX.resize(n_threads);
    mExpAZ.resize(n_threads);
    mExpAY.resize(n_threads);
    mSinAXDY.resize(n_threads);
    mCosAXDY.resize(n_threads);
    mSinAYDZ.resize(n_threads);
    mCosAYDZ.resize(n_threads);
    mSinAZDX.resize(n_threads);
    mCosAZDX.resize(n_threads);
    mCoordinatesAreUpToDate.resize(n_threads);
    for (int i = 0; i < n_threads; i++){
        mCoordinatesAreUpToDate[i] = false;
    }
}

void EthierVelocityField::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mExpD2T[i_thread]  = std::exp(- mD * mD * time);
        mExpAX[i_thread]   = std::exp(mA * coor[0]);
        mExpAY[i_thread]   = std::exp(mA * coor[1]);
        mExpAZ[i_thread]   = std::exp(mA * coor[2]);
        mSinAXDY[i_thread] = std::sin(mA * coor[0] + mD * coor[1]);
        mCosAXDY[i_thread] = std::cos(mA * coor[0] + mD * coor[1]);
        mSinAYDZ[i_thread] = std::sin(mA * coor[1] + mD * coor[2]);
        mCosAYDZ[i_thread] = std::cos(mA * coor[1] + mD * coor[2]);
        mSinAZDX[i_thread] = std::sin(mA * coor[2] + mD * coor[0]);
        mCosAZDX[i_thread] = std::cos(mA * coor[2] + mD * coor[0]);
    }
}

void EthierVelocityField::UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mExpD2T[i_thread]  = std::exp(- mD * mD * time);
        mExpAX[i_thread]   = std::exp(mA * coor[0]);
        mExpAY[i_thread]   = std::exp(mA * coor[1]);
        mExpAZ[i_thread]   = std::exp(mA * coor[2]);
        mSinAXDY[i_thread] = std::sin(mA * coor[0] + mD * coor[1]);
        mCosAXDY[i_thread] = std::cos(mA * coor[0] + mD * coor[1]);
        mSinAYDZ[i_thread] = std::sin(mA * coor[1] + mD * coor[2]);
        mCosAYDZ[i_thread] = std::cos(mA * coor[1] + mD * coor[2]);
        mSinAZDX[i_thread] = std::sin(mA * coor[2] + mD * coor[0]);
        mCosAZDX[i_thread] = std::cos(mA * coor[2] + mD * coor[0]);
    }
}

void EthierVelocityField::LockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 1;
}

void EthierVelocityField::UnlockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 0;
}

// Values

double EthierVelocityField::U0(const int i)
{
    return - mA * (mExpAX[i] * mSinAYDZ[i] + mExpAZ[i] * mCosAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U1(const int i)
{
    return - mA * (mExpAY[i] * mSinAZDX[i] + mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U2(const int i)
{
    return - mA * (mExpAZ[i] * mSinAXDY[i] + mExpAY[i] * mCosAZDX[i]) * mExpD2T[i];
}

// First-order derivatives

double EthierVelocityField::U0DT(const int i)
{
    return - mD * mD * U0(i);
}
double EthierVelocityField::U0D0(const int i)
{
    return - mA * (mA * mExpAX[i] * mSinAYDZ[i] - mA * mExpAZ[i] * mSinAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D1(const int i)
{
    return - mA * (mA * mExpAX[i] * mCosAYDZ[i] - mD * mExpAZ[i] * mSinAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D2(const int i)
{
    return - mA * (mD * mExpAX[i] * mCosAYDZ[i] + mA * mExpAZ[i] * mCosAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U1DT(const int i)
{
    return - mD * mD * U1(i);
}
double EthierVelocityField::U1D0(const int i)
{
    return - mA * (mD * mExpAY[i] * mCosAZDX[i] + mA * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D1(const int i)
{
    return - mA * (mA * mExpAY[i] * mSinAZDX[i] - mA * mExpAX[i] * mSinAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D2(const int i)
{
    return - mA * (mA * mExpAY[i] * mCosAZDX[i] - mD * mExpAX[i] * mSinAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U2DT(const int i)
{
    return - mD * mD * U2(i);
}
double EthierVelocityField::U2D0(const int i)
{
    return - mA * (mA * mExpAZ[i] * mCosAXDY[i] - mD * mExpAY[i] * mSinAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D1(const int i)
{
    return - mA * (mD * mExpAZ[i] * mCosAXDY[i] + mA * mExpAY[i] * mCosAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D2(const int i)
{
    return - mA * (mA * mExpAZ[i] * mSinAXDY[i] - mA * mExpAY[i] * mSinAZDX[i]) * mExpD2T[i];
}

// Second-order derivatives

double EthierVelocityField::U0DTDT(const int i)
{
    return - mD * mD * U0DT(i);
}
double EthierVelocityField::U0DTD0(const int i)
{
    return 000000;
}
double EthierVelocityField::U0DTD1(const int i)
{
    return 000000;
}
double EthierVelocityField::U0DTD2(const int i)
{
    return 000000;
}
double EthierVelocityField::U0D0D0(const int i)
{
    return - mA * (mA * mA * mExpAX[i] * mSinAYDZ[i] - mA * mA * mExpAZ[i] * mCosAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D0D1(const int i)
{
    return - mA * (mA * mA * mExpAX[i] * mCosAYDZ[i] - mA * mD * mExpAZ[i] * mSinAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D0D2(const int i)
{
    return - mA * (mA * mD * mExpAX[i] * mCosAYDZ[i] - mA * mA * mExpAZ[i] * mSinAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D1D1(const int i)
{
    return - mA * (- mA * mA * mExpAX[i] * mSinAYDZ[i] - mD * mD * mExpAZ[i] * mCosAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D1D2(const int i)
{
    return - mA * (- mA * mD * mExpAX[i] * mSinAYDZ[i] - mA * mD * mExpAZ[i] * mSinAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U0D2D2(const int i)
{
    return - mA * (- mD * mD * mExpAX[i] * mSinAYDZ[i] + mA * mA * mExpAZ[i] * mCosAXDY[i]) * mExpD2T[i];
}
double EthierVelocityField::U1DTDT(const int i)
{
    return - mD * mD * U1DT(i);
}
double EthierVelocityField::U1DTD0(const int i)
{
    return 000000;
}
double EthierVelocityField::U1DTD1(const int i)
{
    return 000000;
}
double EthierVelocityField::U1DTD2(const int i)
{
    return 000000;
}
double EthierVelocityField::U1D0D0(const int i)
{
    return - mA * (- mD * mD * mExpAY[i] * mSinAZDX[i] + mA * mA * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D0D1(const int i)
{
    return - mA * (mA * mA * mExpAY[i] * mSinAZDX[i] - mA * mA * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D0D2(const int i)
{
    return - mA * (mA * mA * mExpAY[i] * mCosAZDX[i] - mA * mD * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D1D1(const int i)
{
    return - mA * (mA * mA * mExpAY[i] * mSinAZDX[i] - mA * mA * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D1D2(const int i)
{
    return - mA * (mA * mA * mExpAY[i] * mCosAZDX[i] - mA * mD * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U1D2D2(const int i)
{
    return - mA * (- mA * mA * mExpAY[i] * mSinAZDX[i] - mD * mD * mExpAX[i] * mCosAYDZ[i]) * mExpD2T[i];
}
double EthierVelocityField::U2DTDT(const int i)
{
    return - mD * mD * U2DT(i);
}
double EthierVelocityField::U2DTD0(const int i)
{
    return 000000;
}
double EthierVelocityField::U2DTD1(const int i)
{
    return 000000;
}
double EthierVelocityField::U2DTD2(const int i)
{
    return 000000;
}
double EthierVelocityField::U2D0D0(const int i)
{
    return - mA * (- mA * mA * mExpAZ[i] * mSinAXDY[i] - mD * mD * mExpAY[i] * mCosAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D0D1(const int i)
{
    return - mA * (mD * mD * mExpAZ[i] * mSinAXDY[i] - mA * mA * mExpAY[i] * mSinAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D0D2(const int i)
{
    return - mA * (mA * mA * mExpAZ[i] * mSinAXDY[i] - mA * mA * mExpAY[i] * mCosAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D1D1(const int i)
{
    return - mA * (- mD * mD * mExpAZ[i] * mSinAXDY[i] + mA * mA * mExpAY[i] * mCosAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D1D2(const int i)
{
    return - mA * (mA * mD * mExpAZ[i] * mCosAXDY[i] - mA * mA * mExpAY[i] * mSinAZDX[i]) * mExpD2T[i];
}
double EthierVelocityField::U2D2D2(const int i)
{
    return - mA * (mA * mA * mExpAZ[i] * mSinAXDY[i] - mA * mA * mExpAY[i] * mCosAZDX[i]) * mExpD2T[i];
}

void EthierPressureField::ResizeVectorsForParallelism(const int n_threads)
{
    mExpD2T.resize(n_threads);
    mExpAX.resize(n_threads);
    mExpAZ.resize(n_threads);
    mExpAY.resize(n_threads);
    mSinAXDY.resize(n_threads);
    mCosAXDY.resize(n_threads);
    mSinAYDZ.resize(n_threads);
    mCosAYDZ.resize(n_threads);
    mSinAZDX.resize(n_threads);
    mCosAZDX.resize(n_threads);
    mCoordinatesAreUpToDate.resize(n_threads);
    for (int i = 0; i < n_threads; i++){
        mCoordinatesAreUpToDate[i] = false;
    }
}

double EthierPressureField::Evaluate(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    return this->U(i_thread);
}

//***************************************************************************************************************
//***************************************************************************************************************

double EthierPressureField::CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    return 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

void EthierPressureField::CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& gradient, const int i_thread)
{
    gradient[0] = 0.0;
    gradient[1] = 0.0;
    gradient[2] = 0.0;
}

//***************************************************************************************************************
//***************************************************************************************************************

void EthierPressureField::CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& laplacian, const int i_thread)
{
    laplacian[0] = 0.0;
    laplacian[1] = 0.0;
    laplacian[2] = 0.0;
}

double EthierPressureField::U(const int i)
{
    return - 0.5 * std::pow(mA, 2) * (std::pow(mExpAX[i], 2)
                                    + std::pow(mExpAY[i], 2)
                                    + std::pow(mExpAZ[i], 2)
                                    + 2 * mSinAXDY[i] * mCosAZDX[i] * mExpAY[i] * mExpAZ[i]
                                    + 2 * mSinAYDZ[i] * mCosAXDY[i] * mExpAZ[i] * mExpAX[i]
                                    + 2 * mSinAZDX[i] * mCosAYDZ[i] * mExpAX[i] * mExpAY[i]) * std::pow(mExpD2T[i], 2);
}

// First-order derivatives
double EthierPressureField::UDT(const int i){return 0.0;}
double EthierPressureField::UD1(const int i){return 0.0;}
double EthierPressureField::UD2(const int i){return 0.0;}
double EthierPressureField::UD3(const int i){return 0.0;}

// Second-order derivatives
double EthierPressureField::UDTDT(const int i){return 0.0;}
double EthierPressureField::UDTD0(const int i){return 0.0;}
double EthierPressureField::UDTD1(const int i){return 0.0;}
double EthierPressureField::UDTD2(const int i){return 0.0;}
double EthierPressureField::UD0D0(const int i){return 0.0;}
double EthierPressureField::UD0D1(const int i){return 0.0;}
double EthierPressureField::UD0D2(const int i){return 0.0;}
double EthierPressureField::UD1D1(const int i){return 0.0;}
double EthierPressureField::UD1D2(const int i){return 0.0;}
double EthierPressureField::UD2D2(const int i){return 0.0;}

} // namespace Kratos.




