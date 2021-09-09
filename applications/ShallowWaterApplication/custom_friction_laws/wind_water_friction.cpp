//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "wind_water_friction.h"
#include "custom_utilities/shallow_water_utilities.h"


namespace Kratos
{

void WindWaterFriction::Initialize(const GeometryType& rGeometry, const ProcessInfo& rProcessInfo)
{
    mAirDensity = rProcessInfo[DENSITY_AIR];
    mWaterDensity = rProcessInfo[DENSITY];
}

double WindWaterFriction::CalculateLHS(
    const array_1d<double,3>& rInnerVelocity,
    const array_1d<double,3>& rOuterVelocity)
{
    const auto wind = rOuterVelocity - rInnerVelocity;
    const double abs_wind = norm_2(wind);
    double coefficient;
    if (abs_wind < 1.0)
    {
        coefficient = 0.5e-3 * std::pow(abs_wind, 1./5.);
    }
    else if (abs_wind < 15.0)
    {
        coefficient = 0.5e-3 * std::pow(abs_wind, 1./2.);
    }
    else
    {
        coefficient = 2.6e-3;
    }
    return mAirDensity / mWaterDensity * coefficient * abs_wind;
}

array_1d<double,3> WindWaterFriction::CalculateRHS(
    const array_1d<double,3>& rInnerVelocity,
    const array_1d<double,3>& rOuterVelocity)
{
    const auto wind = rOuterVelocity - rInnerVelocity;
    return wind * CalculateLHS(rInnerVelocity, rOuterVelocity);
}

}  // namespace Kratos
