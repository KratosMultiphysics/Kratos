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
        mSin0[i_thread]  = std::sin(Globals::Pi * coor[0]);
        mCos0[i_thread]  = std::cos(Globals::Pi * coor[0]);
        mSin1[i_thread]  = std::sin(Globals::Pi * coor[1]);
        mCos1[i_thread]  = std::cos(Globals::Pi * coor[1]);
        mSin2[i_thread]  = std::sin(Globals::Pi * coor[2]);
        mCos2[i_thread]  = std::cos(Globals::Pi * coor[2]);
    }
}

void ProductOfSines::UpdateCoordinates(const double time, const vector<double>& coor, const int i_thread)
{
    (void)time;

    if (!mCoordinatesAreUpToDate[i_thread]){
        mSin0[i_thread]  = std::sin(Globals::Pi * coor[0]);
        mCos0[i_thread]  = std::cos(Globals::Pi * coor[0]);
        mSin1[i_thread]  = std::sin(Globals::Pi * coor[1]);
        mCos1[i_thread]  = std::cos(Globals::Pi * coor[1]);
        mSin2[i_thread]  = std::sin(Globals::Pi * coor[2]);
        mCos2[i_thread]  = std::cos(Globals::Pi * coor[2]);
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
    return Globals::Pi * mCos0[i]  * mSin1[i]  * mSin2[i];
}
double ProductOfSines::U0D1(const int i)
{
    return Globals::Pi * mSin0[i]  * mCos1[i]  * mSin2[i];
}
double ProductOfSines::U0D2(const int i)
{
    return Globals::Pi * mSin0[i]  * mSin1[i]  * mCos2[i];
}

// Second-order derivatives
double ProductOfSines::U0D0D0(const int i)
{
    return - Globals::Pi * Globals::Pi * mSin0[i]  * mSin1[i]  * mSin2[i];
}
double ProductOfSines::U0D0D1(const int i)
{
    return Globals::Pi * Globals::Pi * mCos0[i]  * mCos1[i]  * mSin2[i];
}
double ProductOfSines::U0D0D2(const int i)
{
    return Globals::Pi * Globals::Pi * mCos0[i]  * mSin1[i]  * mCos2[i];
}
double ProductOfSines::U0D1D1(const int i)
{
    return - Globals::Pi * Globals::Pi * mSin0[i]  * mSin1[i]  * mSin2[i];
}
double ProductOfSines::U0D1D2(const int i)
{
    return Globals::Pi * Globals::Pi * mSin0[i]  * mCos1[i]  * mCos2[i];
}
double ProductOfSines::U0D2D2(const int i)
{
    return - Globals::Pi * Globals::Pi * mSin0[i]  * mSin1[i]  * mSin2[i];
}

} // namespace Kratos.




