#include "quadratic_linear_time_real_field.h"

namespace Kratos
{
    void QuadraticLinearTimeRealFlowField::ResizeVectorsForParallelism(const int n_threads)
    {
        // Position vector
        for (unsigned i = 0; i < mVars.size(); i++)
            mVars[i].resize(n_threads);

        mCoordinatesAreUpToDate.resize(n_threads);
        mCurrentTime.resize(n_threads);
        for (int i = 0; i < n_threads; i++)
        {
            mCoordinatesAreUpToDate[i] = false;
            mCurrentTime[i] = 0.0;
        }
    }

    void QuadraticLinearTimeRealFlowField::UpdateCoordinates(const double time, const array_1d<double, 3> &coor, const int i_thread)
    {
        mCurrentTime[i_thread] = time;
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

    void QuadraticLinearTimeRealFlowField::UpdateCoordinates(const double time, const DenseVector<double> &coor, const int i_thread)
    {
        mCurrentTime[i_thread] = time;
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

    void QuadraticLinearTimeRealFlowField::LockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 1;
    }

    void QuadraticLinearTimeRealFlowField::UnlockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 0;
    }

    // Values

    double QuadraticLinearTimeRealFlowField::U0(const int i)
    {
        const double time = mCurrentTime[i];
        double prod = 0.0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * (mUxCoefficients0[n] + time * mUxCoefficients1[n]);
        return prod;
    }
    double QuadraticLinearTimeRealFlowField::U1(const int i)
    {
        const double time = mCurrentTime[i];
        double prod = 0.0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * (mUyCoefficients0[n] + time * mUyCoefficients1[n]);
        return prod;
    }
    double QuadraticLinearTimeRealFlowField::U2(const int i)
    {
        const double time = mCurrentTime[i];
        double prod = 0.0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * (mUzCoefficients0[n] + time * mUzCoefficients1[n]);
        return prod;
    }

    // First-order derivatives

    double QuadraticLinearTimeRealFlowField::U0DT(const int i)
    {
        double prod = 0.0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUxCoefficients1[n];
        return prod;
    }
    double QuadraticLinearTimeRealFlowField::U0D0(const int i)
    {
        const double time = mCurrentTime[i];
        const double c0 = mUxCoefficients0[0] + time * mUxCoefficients1[0];
        const double c3 = mUxCoefficients0[3] + time * mUxCoefficients1[3];
        const double c4 = mUxCoefficients0[4] + time * mUxCoefficients1[4];
        const double c6 = mUxCoefficients0[6] + time * mUxCoefficients1[6];
        return 2. * c0 * mVars[6][i] + c3 * mVars[7][i] + c4 * mVars[8][i] + c6;
    }
    double QuadraticLinearTimeRealFlowField::U0D1(const int i)
    {
        const double time = mCurrentTime[i];
        const double c1 = mUxCoefficients0[1] + time * mUxCoefficients1[1];
        const double c3 = mUxCoefficients0[3] + time * mUxCoefficients1[3];
        const double c5 = mUxCoefficients0[5] + time * mUxCoefficients1[5];
        const double c7 = mUxCoefficients0[7] + time * mUxCoefficients1[7];
        return 2. * c1 * mVars[7][i] + c3 * mVars[6][i] + c5 * mVars[8][i] + c7;
    }
    double QuadraticLinearTimeRealFlowField::U0D2(const int i)
    {
        const double time = mCurrentTime[i];
        const double c2 = mUxCoefficients0[2] + time * mUxCoefficients1[2];
        const double c4 = mUxCoefficients0[4] + time * mUxCoefficients1[4];
        const double c5 = mUxCoefficients0[5] + time * mUxCoefficients1[5];
        const double c8 = mUxCoefficients0[8] + time * mUxCoefficients1[8];
        return 2. * c2 * mVars[8][i] + c4 * mVars[6][i] + c5 * mVars[7][i] + c8;
    }
    double QuadraticLinearTimeRealFlowField::U1DT(const int i)
    {
        double prod = 0.0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUyCoefficients1[n];
        return prod;
    }
    double QuadraticLinearTimeRealFlowField::U1D0(const int i)
    {
        const double time = mCurrentTime[i];
        const double c0 = mUyCoefficients0[0] + time * mUyCoefficients1[0];
        const double c3 = mUyCoefficients0[3] + time * mUyCoefficients1[3];
        const double c4 = mUyCoefficients0[4] + time * mUyCoefficients1[4];
        const double c6 = mUyCoefficients0[6] + time * mUyCoefficients1[6];
        return 2. * c0 * mVars[6][i] + c3 * mVars[7][i] + c4 * mVars[8][i] + c6;
    }
    double QuadraticLinearTimeRealFlowField::U1D1(const int i)
    {
        const double time = mCurrentTime[i];
        const double c1 = mUyCoefficients0[1] + time * mUyCoefficients1[1];
        const double c3 = mUyCoefficients0[3] + time * mUyCoefficients1[3];
        const double c5 = mUyCoefficients0[5] + time * mUyCoefficients1[5];
        const double c7 = mUyCoefficients0[7] + time * mUyCoefficients1[7];
        return 2. * c1 * mVars[7][i] + c3 * mVars[6][i] + c5 * mVars[8][i] + c7;
    }
    double QuadraticLinearTimeRealFlowField::U1D2(const int i)
    {
        const double time = mCurrentTime[i];
        const double c2 = mUyCoefficients0[2] + time * mUyCoefficients1[2];
        const double c4 = mUyCoefficients0[4] + time * mUyCoefficients1[4];
        const double c5 = mUyCoefficients0[5] + time * mUyCoefficients1[5];
        const double c8 = mUyCoefficients0[8] + time * mUyCoefficients1[8];
        return 2. * c2 * mVars[8][i] + c4 * mVars[6][i] + c5 * mVars[7][i] + c8;
    }
    double QuadraticLinearTimeRealFlowField::U2DT(const int i)
    {
        double prod = 0.0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUzCoefficients1[n];
        return prod;
    }
    double QuadraticLinearTimeRealFlowField::U2D0(const int i)
    {
        const double time = mCurrentTime[i];
        const double c0 = mUzCoefficients0[0] + time * mUzCoefficients1[0];
        const double c3 = mUzCoefficients0[3] + time * mUzCoefficients1[3];
        const double c4 = mUzCoefficients0[4] + time * mUzCoefficients1[4];
        const double c6 = mUzCoefficients0[6] + time * mUzCoefficients1[6];
        return 2. * c0 * mVars[6][i] + c3 * mVars[7][i] + c4 * mVars[8][i] + c6;
    }
    double QuadraticLinearTimeRealFlowField::U2D1(const int i)
    {
        const double time = mCurrentTime[i];
        const double c1 = mUzCoefficients0[1] + time * mUzCoefficients1[1];
        const double c3 = mUzCoefficients0[3] + time * mUzCoefficients1[3];
        const double c5 = mUzCoefficients0[5] + time * mUzCoefficients1[5];
        const double c7 = mUzCoefficients0[7] + time * mUzCoefficients1[7];
        return 2. * c1 * mVars[7][i] + c3 * mVars[6][i] + c5 * mVars[8][i] + c7;
    }
    double QuadraticLinearTimeRealFlowField::U2D2(const int i)
    {
        const double time = mCurrentTime[i];
        const double c2 = mUzCoefficients0[2] + time * mUzCoefficients1[2];
        const double c4 = mUzCoefficients0[4] + time * mUzCoefficients1[4];
        const double c5 = mUzCoefficients0[5] + time * mUzCoefficients1[5];
        const double c8 = mUzCoefficients0[8] + time * mUzCoefficients1[8];
        return 2. * c2 * mVars[8][i] + c4 * mVars[6][i] + c5 * mVars[7][i] + c8;
    }

    // Second-order derivatives

    double QuadraticLinearTimeRealFlowField::U0DTDT(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0DTD0(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0DTD1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0DTD2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0D0D0(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0D0D1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0D0D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0D1D1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0D1D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U0D2D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1DTDT(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1DTD0(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1DTD1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1DTD2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1D0D0(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1D0D1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1D0D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1D1D1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1D1D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U1D2D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2DTDT(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2DTD0(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2DTD1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2DTD2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2D0D0(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2D0D1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2D0D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2D1D1(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2D1D2(const int i)
    {
        return 0.0;
    }
    double QuadraticLinearTimeRealFlowField::U2D2D2(const int i)
    {
        return 0.0;
    }

}
