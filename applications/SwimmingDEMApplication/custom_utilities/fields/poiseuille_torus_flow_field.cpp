#include "custom_utilities/fields/poiseuille_torus_flow_field.h"

namespace Kratos
{
    double PoiseuilleTorus::getDistanceToCenter(const array_1d<double, 3>& coor)
    {
        double major_dist = std::sqrt(coor[0] * coor[0] + coor[1] * coor[1]);
        double minor_dist = major_dist - mMajorRadius;
        double distance = std::sqrt(minor_dist * minor_dist + coor[2] * coor[2]);

        return distance;
    }

    double PoiseuilleTorus::getThetaAngle(const array_1d<double, 3>& coor)
    {
        return std::atan2(coor[0], coor[1]);
    }

    double PoiseuilleTorus::getVelocityModule()
    {
        return mCenterVelocity * (1 - mRho * mRho);
    }

    void PoiseuilleTorus::UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread = 0)
    {
        mTheta = getThetaAngle(coor);
        mRho = getDistanceToCenter(coor) / mMinorRadius;

        mX = coor[0];
        mY = coor[1];
        mZ = coor[2];
    }

    double PoiseuilleTorus::U0(const int i_thread = 0)
    {
        double velocity_module = getVelocityModule();
        double velocity_x = 1.0 * velocity_module * std::cos(mTheta);
        return velocity_x;
    }

    double PoiseuilleTorus::U1(const int i_thread = 0)
    {
        double velocity_module = getVelocityModule();
        double velocity_x = -1.0 * velocity_module * std::sin(mTheta);
        return velocity_x;
    }

    double PoiseuilleTorus::U2(const int i_thread = 0)
    {
        return 0.0;
    }

    void PoiseuilleTorus::CalculateMaterialAcceleration(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& accel, const int i_thread = 0)
    {
        UpdateCoordinates(time, coor, i_thread);
        double velocity_module = getVelocityModule();
        double major_dist_2 = coor[0] * coor[0] + coor[1] * coor[1];

        accel[0] = -(coor[0] / major_dist_2) * velocity_module * velocity_module;
        accel[1] = -(coor[1] / major_dist_2) * velocity_module * velocity_module;
        accel[2] = 0.0;
    }
}