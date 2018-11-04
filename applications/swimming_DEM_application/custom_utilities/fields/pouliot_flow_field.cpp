#include "pouliot_flow_field.h"

namespace Kratos
{

void PouliotFlowField::ResizeVectorsForParallelism(const int n_threads)
{
    mExpX.resize(n_threads);
    mExpY.resize(n_threads);
    mCoordinatesAreUpToDate.resize(n_threads);
    for (int i = 0; i < n_threads; i++){
        mCoordinatesAreUpToDate[i] = false;
    }
}
void PouliotFlowField::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mExpX[i_thread] = std::exp(- 25 * coor[0]);
        mExpY[i_thread] = std::exp(- 25 * coor[1]);
    }
}

void PouliotFlowField::UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mExpX[i_thread] = std::exp(- 25 * coor[0]);
        mExpY[i_thread] = std::exp(- 25 * coor[1]);
    }
}

void PouliotFlowField::LockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 1;
}

void PouliotFlowField::UnlockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 0;
}

// Values

double PouliotFlowField::U0(const int i)
{
    return mExpX[i] + mExpY[i];
}
double PouliotFlowField::U1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2(const int i)
{
    return 0.0;
}

// First-order derivatives

double PouliotFlowField::U0DT(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0D0(const int i)
{
    return - 25 * mExpX[i];
}
double PouliotFlowField::U0D1(const int i)
{
    return - 25 * mExpY[i];
}
double PouliotFlowField::U0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1DT(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2DT(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D2(const int i)
{
    return 0.0;
}

// Second-order derivatives

double PouliotFlowField::U0DTDT(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0DTD0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0DTD1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0DTD2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0D0D0(const int i)
{
    return 625 * mExpX[i];
}
double PouliotFlowField::U0D0D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0D0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0D1D1(const int i)
{
    return 625 * mExpY[i];
}
double PouliotFlowField::U0D1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U0D2D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1DTDT(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1DTD0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1DTD1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1DTD2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D0D0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D0D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D1D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U1D2D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2DTDT(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2DTD0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2DTD1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2DTD2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D0D0(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D0D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D1D1(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField::U2D2D2(const int i)
{
    return 0.0;
}

} // namespace Kratos.




