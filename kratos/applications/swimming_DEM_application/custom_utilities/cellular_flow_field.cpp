#include "cellular_flow_field.h"

namespace Kratos
{

void CellularFlowField::ResizeVectorsForParallelism(const unsigned int n_threads)
{
    mSinOmegaT.resize(n_threads);
    mCosOmegaT.resize(n_threads);
    mSinPiX0.resize(n_threads);
    mCosPiX0.resize(n_threads);
    mSinPiX1.resize(n_threads);
    mCosPiX1.resize(n_threads);
}

void CellularFlowField::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const unsigned int n_threads, const unsigned int i_thread)
{
    mSinOmegaT[i_thread] = std::sin(mOmegaUOverL * time);
    mCosOmegaT[i_thread] = std::cos(mOmegaUOverL * time);
    mSinPiX0[i_thread]   = std::sin(mPiOverL * coor[0]);
    mCosPiX0[i_thread]   = std::cos(mPiOverL * coor[0]);
    mSinPiX1[i_thread]   = std::sin(mPiOverL * coor[1]);
    mCosPiX1[i_thread]   = std::cos(mPiOverL * coor[1]);
}

// Values

double CellularFlowField::U0(unsigned int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1(unsigned int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U2(unsigned int i){return 0.0;}

// First-order derivatives

double CellularFlowField::U0DT(unsigned int i)
{
    return   mU * mK * mOmegaUOverL * mCosOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D0(unsigned int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D1(unsigned int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U0D2(unsigned int i){return 0.0;}

double CellularFlowField::U1DT(unsigned int i)
{
    return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D0(unsigned int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D1(unsigned int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1D2(unsigned int i){return 0.0;}
double CellularFlowField::U2DT(unsigned int i){return 0.0;}
double CellularFlowField::U2D0(unsigned int i){return 0.0;}
double CellularFlowField::U2D1(unsigned int i){return 0.0;}
double CellularFlowField::U2D2(unsigned int i){return 0.0;}

// Second-order derivatives

double CellularFlowField::U0DTDT(unsigned int i)
{
    return - mU * mK * mOmegaUOverL * mOmegaUOverL * mSinOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0DTD0(unsigned int i)
{
    return   mU * mOmegaUOverL * mCosOmegaT[i] * mPiOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0DTD1(unsigned int i)
{
    return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mPiOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U0DTD2(unsigned int i){return 0.0;}
double CellularFlowField::U0D0D0(unsigned int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mPiOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D0D1(unsigned int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mPiOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U0D0D2(unsigned int i){return 0.0;}
double CellularFlowField::U0D1D1(unsigned int i)
{
    return - mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mPiOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U0D1D2(unsigned int i){return 0.0;}
double CellularFlowField::U0D2D2(unsigned int i){return 0.0;}

double CellularFlowField::U1DTDT(unsigned int i)
{
    return - mU * mK * mOmegaUOverL * mOmegaUOverL * mSinOmegaT[i] * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1DTD0(unsigned int i)
{
    return   mU * mK * mOmegaUOverL * mCosOmegaT[i] * mPiOverL * mSinPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1DTD1(unsigned int i)
{
    return - mU * mK * mOmegaUOverL * mCosOmegaT[i] * mPiOverL * mCosPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1DTD2(unsigned int i){return 0.0;}
double CellularFlowField::U1D0D0(unsigned int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mPiOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D0D1(unsigned int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mPiOverL * mSinPiX0[i] * mCosPiX1[i];
}
double CellularFlowField::U1D0D2(unsigned int i){return 0.0;}
double CellularFlowField::U1D1D1(unsigned int i)
{
    return   mU * (1.0 + mK * mSinOmegaT[i]) * mPiOverL * mPiOverL * mCosPiX0[i] * mSinPiX1[i];
}
double CellularFlowField::U1D1D2(unsigned int i){return 0.0;}
double CellularFlowField::U1D2D2(unsigned int i){return 0.0;}
double CellularFlowField::U2DTDT(unsigned int i){return 0.0;}
double CellularFlowField::U2DTD0(unsigned int i){return 0.0;}
double CellularFlowField::U2DTD1(unsigned int i){return 0.0;}
double CellularFlowField::U2DTD2(unsigned int i){return 0.0;}
double CellularFlowField::U2D0D0(unsigned int i){return 0.0;}
double CellularFlowField::U2D0D1(unsigned int i){return 0.0;}
double CellularFlowField::U2D0D2(unsigned int i){return 0.0;}
double CellularFlowField::U2D1D1(unsigned int i){return 0.0;}
double CellularFlowField::U2D1D2(unsigned int i){return 0.0;}
double CellularFlowField::U2D2D2(unsigned int i){return 0.0;}

} // namespace Kratos.




