#include "quadratic_steady_real_field.h"

namespace Kratos
{
    void QuadraticSteadyRealFlowField::ResizeVectorsForParallelism(const int n_threads)
    {
        // Position vector
        for (unsigned i = 0; i < mVars.size(); i++)
            mVars[i].resize(n_threads);

        mCoordinatesAreUpToDate.resize(n_threads);
        for (int i = 0; i < n_threads; i++)
        {
            mCoordinatesAreUpToDate[i] = false;
        }
    }
    void QuadraticSteadyRealFlowField::UpdateCoordinates(const double time, const array_1d<double, 3> &coor, const int i_thread)
    {
        if (!mCoordinatesAreUpToDate[i_thread])
        {
            mVars[0][i_thread] = (coor[0] - mX0[0]) * (coor[0] - mX0[0]); // x * x
            mVars[1][i_thread] = (coor[1] - mX0[1]) * (coor[1] - mX0[1]); // y * y
            mVars[2][i_thread] = (coor[2] - mX0[2]) * (coor[2] - mX0[2]); // z * z
            mVars[3][i_thread] = (coor[0] - mX0[0]) * (coor[1] - mX0[1]); // x * y
            mVars[4][i_thread] = (coor[0] - mX0[0]) * (coor[2] - mX0[2]); // x * z
            mVars[5][i_thread] = (coor[1] - mX0[1]) * (coor[2] - mX0[2]); // y * z
            mVars[6][i_thread] = coor[0] - mX0[0];                        // x
            mVars[7][i_thread] = coor[1] - mX0[1];                        // y
            mVars[8][i_thread] = coor[2] - mX0[2];                        // z
        }
    }

    void QuadraticSteadyRealFlowField::UpdateCoordinates(const double time, const DenseVector<double> &coor, const int i_thread)
    {
        if (!mCoordinatesAreUpToDate[i_thread])
        {
            mVars[0][i_thread] = (coor[0] - mX0[0]) * (coor[0] - mX0[0]); // x * x
            mVars[1][i_thread] = (coor[1] - mX0[1]) * (coor[1] - mX0[1]); // y * y
            mVars[2][i_thread] = (coor[2] - mX0[2]) * (coor[2] - mX0[2]); // z * z
            mVars[3][i_thread] = (coor[0] - mX0[0]) * (coor[1] - mX0[1]); // x * y
            mVars[4][i_thread] = (coor[0] - mX0[0]) * (coor[2] - mX0[2]); // x * z
            mVars[5][i_thread] = (coor[1] - mX0[1]) * (coor[2] - mX0[2]); // y * z
            mVars[6][i_thread] = coor[0] - mX0[0];                        // x
            mVars[7][i_thread] = coor[1] - mX0[1];                        // y
            mVars[8][i_thread] = coor[2] - mX0[2];                        // z
        }
    }

    void QuadraticSteadyRealFlowField::LockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 1;
    }

    void QuadraticSteadyRealFlowField::UnlockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 0;
    }

    // Values

    double QuadraticSteadyRealFlowField::U0(const int i)
    {
        double prod = 0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUxCoefficients[n];
        return prod;
    }
    double QuadraticSteadyRealFlowField::U1(const int i)
    {
        double prod = 0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUyCoefficients[n];
        return prod;
    }
    double QuadraticSteadyRealFlowField::U2(const int i)
    {
        double prod = 0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUzCoefficients[n];
        return prod;
    }

    // First-order derivatives

    double QuadraticSteadyRealFlowField::U0DT(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0D0(const int i)
    {
        return 2. * mUxCoefficients[0] * mVars[6][i] + mUxCoefficients[3] * mVars[7][i] + mUxCoefficients[4] * mVars[8][i] + mUxCoefficients[6];
    }
    double QuadraticSteadyRealFlowField::U0D1(const int i)
    {
        return 2. * mUxCoefficients[1] * mVars[7][i] + mUxCoefficients[3] * mVars[6][i] + mUxCoefficients[5] * mVars[8][i] + mUxCoefficients[7];
    }
    double QuadraticSteadyRealFlowField::U0D2(const int i)
    {
        return 2. * mUxCoefficients[2] * mVars[8][i] + mUxCoefficients[4] * mVars[6][i] + mUxCoefficients[5] * mVars[7][i] + mUxCoefficients[8];
    }
    double QuadraticSteadyRealFlowField::U1DT(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D0(const int i)
    {
        return 2. * mUyCoefficients[0] * mVars[6][i] + mUyCoefficients[3] * mVars[7][i] + mUyCoefficients[4] * mVars[8][i] + mUyCoefficients[6];
    }
    double QuadraticSteadyRealFlowField::U1D1(const int i)
    {
        return 2. * mUyCoefficients[1] * mVars[7][i] + mUyCoefficients[3] * mVars[6][i] + mUyCoefficients[5] * mVars[8][i] + mUyCoefficients[7];
    }
    double QuadraticSteadyRealFlowField::U1D2(const int i)
    {
        return 2. * mUyCoefficients[2] * mVars[8][i] + mUyCoefficients[4] * mVars[6][i] + mUyCoefficients[5] * mVars[7][i] + mUyCoefficients[8];
    }
    double QuadraticSteadyRealFlowField::U2DT(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D0(const int i)
    {
        return 2. * mUzCoefficients[0] * mVars[6][i] + mUzCoefficients[3] * mVars[7][i] + mUzCoefficients[4] * mVars[8][i] + mUzCoefficients[6];
    }
    double QuadraticSteadyRealFlowField::U2D1(const int i)
    {
        return 2. * mUzCoefficients[1] * mVars[7][i] + mUzCoefficients[3] * mVars[6][i] + mUzCoefficients[5] * mVars[8][i] + mUzCoefficients[7];
    }
    double QuadraticSteadyRealFlowField::U2D2(const int i)
    {
        return 2. * mUzCoefficients[2] * mVars[8][i] + mUzCoefficients[4] * mVars[6][i] + mUzCoefficients[5] * mVars[7][i] + mUzCoefficients[8];
    }

    // Second-order derivatives

    double QuadraticSteadyRealFlowField::U0DTDT(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0DTD0(const int i)
    {
        return 000000;
    }
    double QuadraticSteadyRealFlowField::U0DTD1(const int i)
    {
        return 000000;
    }
    double QuadraticSteadyRealFlowField::U0DTD2(const int i)
    {
        return 000000;
    }
    double QuadraticSteadyRealFlowField::U0D0D0(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0D0D1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0D0D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0D1D1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0D1D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U0D2D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1DTDT(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1DTD0(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1DTD1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1DTD2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D0D0(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D0D1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D0D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D1D1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D1D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U1D2D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2DTDT(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2DTD0(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2DTD1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2DTD2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D0D0(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D0D1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D0D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D1D1(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D1D2(const int i)
    {
        return 0.;
    }
    double QuadraticSteadyRealFlowField::U2D2D2(const int i)
    {
        return 0.;
    }

}