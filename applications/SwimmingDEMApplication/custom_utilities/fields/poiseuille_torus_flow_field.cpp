#include "poiseuille_torus_flow_field.h"

namespace Kratos
{
    void PoiseuilleTorusFlowField::ResizeVectorsForParallelism(const int n_threads)
    {
        std::cout << "Resizing with n_threads = " << n_threads << std::endl;
        mXYDistance.resize(n_threads);
        mRho.resize(n_threads);
        mSin.resize(n_threads);
        mCos.resize(n_threads);
        mZ.resize(n_threads);
        mCommonTerm.resize(n_threads);

        mCoordinatesAreUpToDate.resize(n_threads);
        for (int i = 0; i < n_threads; i++)
        {
            mCoordinatesAreUpToDate[i] = false;
        }
    }

    void PoiseuilleTorusFlowField::UpdateCoordinates(const double time, const array_1d<double, 3> &coor, const int i_thread)
    {
        if (!mCoordinatesAreUpToDate[i_thread])
        {
            double xy_distance = std::sqrt(coor[0] * coor[0] + coor[1] * coor[1]);
            double local_xy_distance = xy_distance - mMajorRadius;

            mXYDistance[i_thread] = xy_distance;
            mRho[i_thread] = std::sqrt(local_xy_distance * local_xy_distance + coor[2] * coor[2]) / mMinorRadius;
            mSin[i_thread] = coor[0] / xy_distance;
            mCos[i_thread] = coor[1] / xy_distance;
            mZ[i_thread] = coor[2];
            mCommonTerm[i_thread] = (4. / mU0) * (mXYDistance[i_thread] - mMajorRadius) / (mMinorRadius * mMinorRadius) * mRho[i_thread];

            std::cout << "Updating coords for i_thread = " << i_thread << " for coor = " << coor << ": (mCoordinatesAreUpToDate = " << mCoordinatesAreUpToDate[i_thread] << ")" << std::endl;
            std::cout << "xy_distance = " << xy_distance << std::endl;
            std::cout << "local_xy    = " << local_xy_distance << std::endl;
            std::cout << "mXyDist     = " << mXYDistance[i_thread] << std::endl;
            std::cout << "rho         = " << mRho[i_thread] << std::endl;
            std::cout << "sin         = " << mSin[i_thread] << std::endl;
            std::cout << "cos         = " << mCos[i_thread] << std::endl;
            std::cout << "z           = " << mZ[i_thread] << std::endl;
            std::cout << "common_term = " << mCommonTerm[i_thread] << std::endl;
        }
    }

    void PoiseuilleTorusFlowField::UpdateCoordinates(const double time, const DenseVector<double> &coor, const int i_thread)
    {
        if (!mCoordinatesAreUpToDate[i_thread])
        {
            double xy_distance = std::sqrt(coor[0] * coor[0] + coor[1] * coor[1]);
            double local_xy_distance = xy_distance - mMajorRadius;

            mXYDistance[i_thread] = xy_distance;
            mRho[i_thread] = std::sqrt(local_xy_distance * local_xy_distance + coor[2] * coor[2]) / mMinorRadius;
            mSin[i_thread] = coor[0] / xy_distance;
            mCos[i_thread] = coor[1] / xy_distance;
            mZ[i_thread] = coor[2];
            mCommonTerm[i_thread] = (4. / mU0) * (mXYDistance[i_thread] - mMajorRadius) / (mMinorRadius * mMinorRadius) * mRho[i_thread];
        }
    }

    void PoiseuilleTorusFlowField::LockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 1;
    }

    void PoiseuilleTorusFlowField::UnlockCoordinates(const int i_thread)
    {
        mCoordinatesAreUpToDate[i_thread] = 0;
    }

    // Values

    double PoiseuilleTorusFlowField::U0(const int i_thread)
    {
        // double velocity_module = getVelocityModule(i_thread);
        return mU0 * mCos[i_thread];
    }

    double PoiseuilleTorusFlowField::U1(const int i_thread)
    {
        return -1. * mU0 * mSin[i_thread];
    }

    double PoiseuilleTorusFlowField::U2(const int i_thread)
    {
        return 0.0;
    }

    // First-order derivatives

    double PoiseuilleTorusFlowField::U0DT(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D0(const int i)
    {
        return mCos[i] * U1(i) / mXYDistance[i] - mCommonTerm[i] * mCos[i] * mSin[i];
    }
    double PoiseuilleTorusFlowField::U0D1(const int i)
    {
        return -mSin[i] * U1(i) / mXYDistance[i] - mCommonTerm[i] * mCos[i] * mCos[i];
    }
    double PoiseuilleTorusFlowField::U0D2(const int i)
    {
        return -4. * mU0 / (mMinorRadius * mMinorRadius) * mZ[i] * mCos[i] * mRho[i];
    }
    double PoiseuilleTorusFlowField::U1DT(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D0(const int i)
    {
        return -mCos[i] * U0(i) / mXYDistance[i] + mCommonTerm[i] * mSin[i] * mSin[i];
    }
    double PoiseuilleTorusFlowField::U1D1(const int i)
    {
        return mSin[i] * U0(i) / mXYDistance[i] + mCommonTerm[i] * mCos[i] * mSin[i];
    }
    double PoiseuilleTorusFlowField::U1D2(const int i)
    {
        return 4. * mU0 / (mMinorRadius * mMinorRadius) * mZ[i] * mSin[i] * mRho[i];
    }
    double PoiseuilleTorusFlowField::U2DT(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D2(const int i)
    {
        return 0.0;
    }

    // Second-order derivatives

    double PoiseuilleTorusFlowField::U0DTDT(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0DTD0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0DTD1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0DTD2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D0D0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D0D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D0D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D1D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D1D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U0D2D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1DTDT(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1DTD0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1DTD1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1DTD2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D0D0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D0D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D0D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D1D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D1D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U1D2D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2DTDT(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2DTD0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2DTD1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2DTD2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D0D0(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D0D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D0D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D1D1(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D1D2(const int i)
    {
        return 0.0;
    }
    double PoiseuilleTorusFlowField::U2D2D2(const int i)
    {
        return 0.0;
    }
}