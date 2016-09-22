#include "cellular_flow_field.h"

namespace Kratos
{

void CellularFlowField::UpdateCoordinates(const double time, const array_1d<double, 3>& coor)
{
    mSinOmegaT = std::sin(mOmegaUOverL);
    mCosOmegaT = std::cos(mOmegaUOverL);
    mSinPiX0   = std::sin(mPiOverL * coor[0]);
    mCosPiX0   = std::cos(mPiOverL * coor[0]);
    mSinPiX1   = std::sin(mPiOverL * coor[1]);
    mCosPiX1   = std::cos(mPiOverL * coor[1]);
}

// Values

double CellularFlowField::U0()
{
    return   mU * (1.0 + mK * mSinOmegaT) * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U1()
{
    return - mU * (1.0 + mK * mSinOmegaT) * mCosPiX0 * mSinPiX1;
}
double CellularFlowField::U2(){return 0.0;}

// First-order derivatives

double CellularFlowField::U0DT()
{
    return   mU * mOmegaUOverL * mCosOmegaT * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U0D0()
{
    return   mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mCosPiX0 * mCosPiX1;
}
double CellularFlowField::U0D1()
{
    return - mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mSinPiX0 * mSinPiX1;
}
double CellularFlowField::U0D2(){return 0.0;}

double CellularFlowField::U1DT()
{
    return - mU * mOmegaUOverL * mCosOmegaT * mCosPiX0 * mSinPiX1;
}
double CellularFlowField::U1D0()
{
    return   mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mSinPiX0 * mSinPiX1;
}
double CellularFlowField::U1D1()
{
    return - mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mCosPiX0 * mCosPiX1;
}
double CellularFlowField::U1D2(){return 0.0;}
double CellularFlowField::U2DT(){return 0.0;}
double CellularFlowField::U2D0(){return 0.0;}
double CellularFlowField::U2D1(){return 0.0;}
double CellularFlowField::U2D2(){return 0.0;}

// Second-order derivatives

double CellularFlowField::U0DTDT()
{
    return - mU * mOmegaUOverL * mOmegaUOverL * mSinOmegaT * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U0DTD0()
{
    return   mU * mOmegaUOverL * mCosOmegaT * mPiOverL * mCosPiX0 * mCosPiX1;
}
double CellularFlowField::U0DTD1()
{
    return - mU * mOmegaUOverL * mCosOmegaT * mPiOverL * mSinPiX0 * mSinPiX1;
}
double CellularFlowField::U0DTD2(){return 0.0;}
double CellularFlowField::U0D0D0()
{
    return - mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mPiOverL * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U0D0D1()
{
    return - mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mPiOverL * mCosPiX0 * mSinPiX1;
}
double CellularFlowField::U0D0D2(){return 0.0;}
double CellularFlowField::U0D1D1()
{
    return - mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mPiOverL * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U0D1D2(){return 0.0;}
double CellularFlowField::U0D2D2(){return 0.0;}

double CellularFlowField::U1DTDT()
{
    return - mU * mOmegaUOverL * mOmegaUOverL * mSinOmegaT * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U1DTD0()
{
    return   mU * mOmegaUOverL * mCosOmegaT * mPiOverL * mSinPiX0 * mSinPiX1;
}
double CellularFlowField::U1DTD1()
{
    return - mU * mOmegaUOverL * mCosOmegaT * mPiOverL * mCosPiX0 * mCosPiX1;
}
double CellularFlowField::U1DTD2(){return 0.0;}
double CellularFlowField::U1D0D0()
{
    return   mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mPiOverL * mCosPiX0 * mSinPiX1;
}
double CellularFlowField::U1D0D1()
{
    return   mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mPiOverL * mSinPiX0 * mCosPiX1;
}
double CellularFlowField::U1D0D2(){return 0.0;}
double CellularFlowField::U1D1D1()
{
    return   mU * (1.0 + mK * mSinOmegaT) * mPiOverL * mPiOverL * mCosPiX0 * mSinPiX1;
}
double CellularFlowField::U1D1D2(){return 0.0;}
double CellularFlowField::U1D2D2(){return 0.0;}

double CellularFlowField::U2DTDT(){return 0.0;}
double CellularFlowField::U2DTD0(){return 0.0;}
double CellularFlowField::U2DTD1(){return 0.0;}
double CellularFlowField::U2DTD2(){return 0.0;}
double CellularFlowField::U2D0D0(){return 0.0;}
double CellularFlowField::U2D0D1(){return 0.0;}
double CellularFlowField::U2D0D2(){return 0.0;}
double CellularFlowField::U2D1D1(){return 0.0;}
double CellularFlowField::U2D1D2(){return 0.0;}
double CellularFlowField::U2D2D2(){return 0.0;}

} // namespace Kratos.




