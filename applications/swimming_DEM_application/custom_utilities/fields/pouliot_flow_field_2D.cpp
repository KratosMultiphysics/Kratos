#include "pouliot_flow_field_2D.h"

namespace Kratos
{

void PouliotFlowField2D::ResizeVectorsForParallelism(const int n_threads)
{
    mX.resize(n_threads);
    mY.resize(n_threads);
    mCoordinatesAreUpToDate.resize(n_threads);
    for (int i = 0; i < n_threads; i++){
        mCoordinatesAreUpToDate[i] = false;
    }
}
void PouliotFlowField2D::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mX[i_thread] = coor[0];
        mY[i_thread] = coor[1];
    }
}

void PouliotFlowField2D::UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread)
{
    if (!mCoordinatesAreUpToDate[i_thread]){
        mX[i_thread] = coor[0];
        mY[i_thread] = coor[1];
    }
}

void PouliotFlowField2D::LockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 1;
}

void PouliotFlowField2D::UnlockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 0;
}

// Values

double PouliotFlowField2D::U0(const int i)
{
    const double x = mX[i];
    const double y = mY[i];
    return x*x*x + y*y*y;
}
double PouliotFlowField2D::U1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2(const int i)
{
    return 0.0;
}

// First-order derivatives

double PouliotFlowField2D::U0DT(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0D0(const int i)
{
    const double x = mX[i];
    return 3*x*x;
}
double PouliotFlowField2D::U0D1(const int i)
{
    const double y = mY[i];
    return 3*y*y;
}
double PouliotFlowField2D::U0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1DT(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2DT(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D2(const int i)
{
    return 0.0;
}

// Second-order derivatives

double PouliotFlowField2D::U0DTDT(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0DTD0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0DTD1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0DTD2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0D0D0(const int i)
{
    const double x = mX[i];
    return 6*x;
}
double PouliotFlowField2D::U0D0D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0D0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0D1D1(const int i)
{
    const double y = mY[i];
    return 6*y;
}
double PouliotFlowField2D::U0D1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U0D2D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1DTDT(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1DTD0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1DTD1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1DTD2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D0D0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D0D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D1D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U1D2D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2DTDT(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2DTD0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2DTD1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2DTD2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D0D0(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D0D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D0D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D1D1(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D1D2(const int i)
{
    return 0.0;
}
double PouliotFlowField2D::U2D2D2(const int i)
{
    return 0.0;
}

} // namespace Kratos.




