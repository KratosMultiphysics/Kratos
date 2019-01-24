#include "product_of_sines_field.h"

namespace Kratos
{

void ProductOfSines::ResizeVectorsForParallelism(const int n_threads)
{
    mSin0.resize(n_threads);
    mCos0.resize(n_threads);
    mSin1.resize(n_threads);
    mCos1.resize(n_threads);
    mSin2.resize(n_threads);
    mCos2.resize(n_threads);

    mCoordinatesAreUpToDate.resize(n_threads);
    for (int i = 0; i < n_threads; i++){
        mCoordinatesAreUpToDate[i] = false;
    }
}
void ProductOfSines::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    (void)time;

    if (!mCoordinatesAreUpToDate[i_thread]){
        const double alpha_0 = mOmega * coor[0];
        const double alpha_1 = mOmega * coor[1];
        const double alpha_2 = mOmega * coor[2];
        mSin0[i_thread] = std::sin(alpha_0);
        mCos0[i_thread] = std::cos(alpha_0);
        mSin1[i_thread] = std::sin(alpha_1);
        mCos1[i_thread] = std::cos(alpha_1);
        mSin2[i_thread] = std::sin(alpha_2);
        mCos2[i_thread] = std::cos(alpha_2);
    }
}

void ProductOfSines::UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread)
{
    (void)time;

    if (!mCoordinatesAreUpToDate[i_thread]){
        const double alpha_0 = mOmega * coor[0];
        const double alpha_1 = mOmega * coor[1];
        const double alpha_2 = mOmega * coor[2];
        mSin0[i_thread] = std::sin(alpha_0);
        mCos0[i_thread] = std::cos(alpha_0);
        mSin1[i_thread] = std::sin(alpha_1);
        mCos1[i_thread] = std::cos(alpha_1);
        mSin2[i_thread] = std::sin(alpha_2);
        mCos2[i_thread] = std::cos(alpha_2);
    }
}

void ProductOfSines::LockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 1;
}

void ProductOfSines::UnlockCoordinates(const int i_thread)
{
    mCoordinatesAreUpToDate[i_thread] = 0;
}

// Values

double ProductOfSines::U0(const int i)
{
    return mSin0[i]  * mSin1[i]  * mSin2[i];
}

// First-order derivatives

double ProductOfSines::U0D0(const int i)
{
    return mOmega * mCos0[i]  * mSin1[i]  * mSin2[i];
}
double ProductOfSines::U0D1(const int i)
{
    return mOmega * mSin0[i]  * mCos1[i]  * mSin2[i];
}
double ProductOfSines::U0D2(const int i)
{
    return mOmega * mSin0[i]  * mSin1[i]  * mCos2[i];
}

// Second-order derivatives
double ProductOfSines::U0D0D0(const int i)
{
    return - mOmega * mOmega * mSin0[i]  * mSin1[i]  * mSin2[i];
}
double ProductOfSines::U0D0D1(const int i)
{
    return mOmega * mOmega * mCos0[i]  * mCos1[i]  * mSin2[i];
}
double ProductOfSines::U0D0D2(const int i)
{
    return mOmega * mOmega * mCos0[i]  * mSin1[i]  * mCos2[i];
}
double ProductOfSines::U0D1D1(const int i)
{
    return - mOmega * mOmega * mSin0[i]  * mSin1[i]  * mSin2[i];
}
double ProductOfSines::U0D1D2(const int i)
{
    return mOmega * mOmega * mSin0[i]  * mCos1[i]  * mCos2[i];
}
double ProductOfSines::U0D2D2(const int i)
{
    return - mOmega * mOmega * mSin0[i]  * mSin1[i]  * mSin2[i];
}

} // namespace Kratos.




