#include "cubic_steady_real_field.h"

namespace Kratos
{
    void CubicSteadyRealFlowField::ResizeVectorsForParallelism(const int n_threads)
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
    void CubicSteadyRealFlowField::UpdateCoordinates(const double time, const array_1d<double, 3> &coor, const int i_thread)
    {
        if (!mCoordinatesAreUpToDate[i_thread])
        {
            double x = coor[0] - mX0[0];
            double y = coor[1] - mX0[1];
            double z = coor[2] - mX0[2];

            // Order 3
            mVars[0][i_thread] = x * x * x;  // x * x * x
            mVars[1][i_thread] = y * y * y;  // y * y * y
            mVars[2][i_thread] = z * z * z;  // z * z * z
            mVars[3][i_thread] = x * x * y;  // x * x * y
            mVars[4][i_thread] = x * x * z;  // x * x * z
            mVars[5][i_thread] = x * y * y;  // x * y * y
            mVars[6][i_thread] = y * y * z;  // y * y * z
            mVars[7][i_thread] = x * z * z;  // x * z * z
            mVars[8][i_thread] = y * z * z;  // y * z * z
            mVars[9][i_thread] = x * y * z;  // x * y * z

            // Order 2
            mVars[10][i_thread] = x * x;  // x * x
            mVars[11][i_thread] = y * y;  // y * y
            mVars[12][i_thread] = z * z;  // z * z
            mVars[13][i_thread] = x * y;  // x * y
            mVars[14][i_thread] = x * z;  // x * z
            mVars[15][i_thread] = y * z;  // y * z

            // Order 1
            mVars[16][i_thread] = x;  // x
            mVars[17][i_thread] = y;  // y
            mVars[18][i_thread] = z;  // z
        }
    }

    void CubicSteadyRealFlowField::UpdateCoordinates(const double time, const DenseVector<double> &coor, const int i_thread)
    {
        if (!mCoordinatesAreUpToDate[i_thread])
        {
            double x = coor[0] - mX0[0];
            double y = coor[1] - mX0[1];
            double z = coor[2] - mX0[2];

            // Order 3
            mVars[0][i_thread] = x * x * x;
            mVars[1][i_thread] = y * y * y;
            mVars[2][i_thread] = z * z * z;
            mVars[3][i_thread] = x * x * y;
            mVars[4][i_thread] = x * x * z;
            mVars[5][i_thread] = x * y * y;
            mVars[6][i_thread] = y * y * z;
            mVars[7][i_thread] = x * z * z;
            mVars[8][i_thread] = y * z * z;
            mVars[9][i_thread] = x * y * z;

            // Order 2
            mVars[10][i_thread] = x * x;
            mVars[11][i_thread] = y * y;
            mVars[12][i_thread] = z * z;
            mVars[13][i_thread] = x * y;
            mVars[14][i_thread] = x * z;
            mVars[15][i_thread] = y * z;

            // Order 1
            mVars[16][i_thread] = x;
            mVars[17][i_thread] = y;
            mVars[18][i_thread] = z;
        }
    }

    void CubicSteadyRealFlowField::LockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 1;
    }

    void CubicSteadyRealFlowField::UnlockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 0;
    }

    // Values

    double CubicSteadyRealFlowField::U0(const int i)
    {
        double prod = 0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUxCoefficients[n];
        return prod;
    }
    double CubicSteadyRealFlowField::U1(const int i)
    {
        double prod = 0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUyCoefficients[n];
        return prod;
    }
    double CubicSteadyRealFlowField::U2(const int i)
    {
        double prod = 0;
        for (unsigned n = 0; n < mVars.size(); n++)
            prod += mVars[n][i] * mUzCoefficients[n];
        return prod;
    }

    // First-order derivatives

    double CubicSteadyRealFlowField::U0DT(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0D0(const int i)
    {
        double u0d0 = 0.0;

        // Terms of order 3
        u0d0 += 3. * mUxCoefficients[0] * mVars[10][i];
        u0d0 += 2. * mUxCoefficients[3] * mVars[16][i] * mVars[17][i] + 2. * mUxCoefficients[4] * mVars[16][i] * mVars[18][i];
        u0d0 += mUxCoefficients[5] * mVars[11][i] + mUxCoefficients[7] * mVars[12][i] + mUxCoefficients[9] * mVars[17][i] * mVars[18][i];

        // Terms of order 2
        u0d0 += 2. * mUxCoefficients[10] * mVars[16][i] + mUxCoefficients[13] * mVars[17][i] + mUxCoefficients[14] * mVars[18][i];

        // Terms of order 1
        u0d0 += mUxCoefficients[16];

        return u0d0;
    }
    double CubicSteadyRealFlowField::U0D1(const int i)
    {
        double u0d1 = 0.0;

        // Terms of order 3
        u0d1 += 3. * mUxCoefficients[1] * mVars[11][i];
        u0d1 += 2. * mUxCoefficients[5] * mVars[16][i] * mVars[17][i] + 2. * mUxCoefficients[6] * mVars[17][i] * mVars[18][i];
        u0d1 += mUxCoefficients[3] * mVars[10][i] + mUxCoefficients[8] * mVars[12][i] + mUxCoefficients[9] * mVars[16][i] * mVars[18][i];

        // Terms of order 2
        u0d1 += 2. * mUxCoefficients[11] * mVars[17][i] + mUxCoefficients[13] * mVars[16][i] + mUxCoefficients[15] * mVars[18][i];

        // Terms of order 1
        u0d1 += mUxCoefficients[17];

        return u0d1;
    }
    double CubicSteadyRealFlowField::U0D2(const int i)
    {
        double u0d2 = 0.0;

        // Terms of order 3
        u0d2 += 3. * mUxCoefficients[2] * mVars[12][i];
        u0d2 += 2. * mUxCoefficients[7] * mVars[16][i] * mVars[18][i] + 2. * mUxCoefficients[8] * mVars[17][i] * mVars[18][i];
        u0d2 += mUxCoefficients[4] * mVars[10][i] + mUxCoefficients[6] * mVars[11][i] + mUxCoefficients[9] * mVars[16][i] * mVars[17][i];

        // Terms of order 2
        u0d2 += 2. * mUxCoefficients[12] * mVars[18][i] + mUxCoefficients[14] * mVars[16][i] + mUxCoefficients[15] * mVars[17][i];

        // Terms of order 1
        u0d2 += mUxCoefficients[18];

        return u0d2;
    }
    double CubicSteadyRealFlowField::U1DT(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D0(const int i)
    {
        double u1d0 = 0.0;

        // Terms of order 3
        u1d0 += 3. * mUyCoefficients[0] * mVars[10][i];
        u1d0 += 2. * mUyCoefficients[3] * mVars[16][i] * mVars[17][i] + 2. * mUyCoefficients[4] * mVars[16][i] * mVars[18][i];
        u1d0 += mUyCoefficients[5] * mVars[11][i] + mUyCoefficients[7] * mVars[12][i] + mUyCoefficients[9] * mVars[17][i] * mVars[18][i];

        // Terms of order 2
        u1d0 += 2. * mUyCoefficients[10] * mVars[16][i] + mUyCoefficients[13] * mVars[17][i] + mUyCoefficients[14] * mVars[18][i];

        // Terms of order 1
        u1d0 += mUyCoefficients[16];

        return u1d0;
    }
    double CubicSteadyRealFlowField::U1D1(const int i)
    {
        double u1d1 = 0.0;

        // Terms of order 3
        u1d1 += 3. * mUyCoefficients[1] * mVars[11][i];
        u1d1 += 2. * mUyCoefficients[5] * mVars[16][i] * mVars[17][i] + 2. * mUyCoefficients[6] * mVars[17][i] * mVars[18][i];
        u1d1 += mUyCoefficients[3] * mVars[10][i] + mUyCoefficients[8] * mVars[12][i] + mUyCoefficients[9] * mVars[16][i] * mVars[18][i];

        // Terms of order 2
        u1d1 += 2. * mUyCoefficients[11] * mVars[17][i] + mUyCoefficients[13] * mVars[16][i] + mUyCoefficients[15] * mVars[18][i];

        // Terms of order 1
        u1d1 += mUyCoefficients[17];

        return u1d1;
    }
    double CubicSteadyRealFlowField::U1D2(const int i)
    {
        double u1d2 = 0.0;

        // Terms of order 3
        u1d2 += 3. * mUyCoefficients[2] * mVars[12][i];
        u1d2 += 2. * mUyCoefficients[7] * mVars[16][i] * mVars[18][i] + 2. * mUyCoefficients[8] * mVars[17][i] * mVars[18][i];
        u1d2 += mUyCoefficients[4] * mVars[10][i] + mUyCoefficients[6] * mVars[11][i] + mUyCoefficients[9] * mVars[16][i] * mVars[17][i];

        // Terms of order 2
        u1d2 += 2. * mUyCoefficients[12] * mVars[18][i] + mUyCoefficients[14] * mVars[16][i] + mUyCoefficients[15] * mVars[17][i];

        // Terms of order 1
        u1d2 += mUyCoefficients[18];

        return u1d2;
    }
    double CubicSteadyRealFlowField::U2DT(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D0(const int i)
    {
        double u2d0 = 0.0;

        // Terms of order 3
        u2d0 += 3. * mUzCoefficients[0] * mVars[10][i];
        u2d0 += 2. * mUzCoefficients[3] * mVars[16][i] * mVars[17][i] + 2. * mUzCoefficients[4] * mVars[16][i] * mVars[18][i];
        u2d0 += mUzCoefficients[5] * mVars[11][i] + mUzCoefficients[7] * mVars[12][i] + mUzCoefficients[9] * mVars[17][i] * mVars[18][i];

        // Terms of order 2
        u2d0 += 2. * mUzCoefficients[10] * mVars[16][i] + mUzCoefficients[13] * mVars[17][i] + mUzCoefficients[14] * mVars[18][i];

        // Terms of order 1
        u2d0 += mUzCoefficients[16];

        return u2d0;
    }
    double CubicSteadyRealFlowField::U2D1(const int i)
    {
        double u1d1 = 0.0;

        // Terms of order 3
        u1d1 += 3. * mUzCoefficients[1] * mVars[11][i];
        u1d1 += 2. * mUzCoefficients[5] * mVars[16][i] * mVars[17][i] + 2. * mUzCoefficients[6] * mVars[17][i] * mVars[18][i];
        u1d1 += mUzCoefficients[3] * mVars[10][i] + mUzCoefficients[8] * mVars[12][i] + mUzCoefficients[9] * mVars[16][i] * mVars[18][i];

        // Terms of order 2
        u1d1 += 2. * mUzCoefficients[11] * mVars[17][i] + mUzCoefficients[13] * mVars[16][i] + mUzCoefficients[15] * mVars[18][i];

        // Terms of order 1
        u1d1 += mUzCoefficients[17];

        return u1d1;
    }
    double CubicSteadyRealFlowField::U2D2(const int i)
    {
        double u2d2 = 0.0;

        // Terms of order 3
        u2d2 += 3. * mUzCoefficients[2] * mVars[12][i];
        u2d2 += 2. * mUzCoefficients[7] * mVars[16][i] * mVars[18][i] + 2. * mUzCoefficients[8] * mVars[17][i] * mVars[18][i];
        u2d2 += mUzCoefficients[4] * mVars[10][i] + mUzCoefficients[6] * mVars[11][i] + mUzCoefficients[9] * mVars[16][i] * mVars[17][i];

        // Terms of order 2
        u2d2 += 2. * mUzCoefficients[12] * mVars[18][i] + mUzCoefficients[14] * mVars[16][i] + mUzCoefficients[15] * mVars[17][i];

        // Terms of order 1
        u2d2 += mUzCoefficients[18];

        return u2d2;
    }

    // Second-order derivatives

    double CubicSteadyRealFlowField::U0DTDT(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0DTD0(const int i)
    {
        return 000000;
    }
    double CubicSteadyRealFlowField::U0DTD1(const int i)
    {
        return 000000;
    }
    double CubicSteadyRealFlowField::U0DTD2(const int i)
    {
        return 000000;
    }
    double CubicSteadyRealFlowField::U0D0D0(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0D0D1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0D0D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0D1D1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0D1D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U0D2D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1DTDT(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1DTD0(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1DTD1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1DTD2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D0D0(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D0D1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D0D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D1D1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D1D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U1D2D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2DTDT(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2DTD0(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2DTD1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2DTD2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D0D0(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D0D1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D0D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D1D1(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D1D2(const int i)
    {
        return 0.;
    }
    double CubicSteadyRealFlowField::U2D2D2(const int i)
    {
        return 0.;
    }

}