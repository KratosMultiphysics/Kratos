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
#include "manning_law.h"
#include "custom_utilities/shallow_water_utilities.h"


namespace Kratos
{

void ManningLaw::Initialize(const GeometryType& rGeometry, const ProcessInfo& rProcessInfo)
{
    double manning = 0.0;
    for (auto& r_node : rGeometry)
    {
        manning += r_node.FastGetSolutionStepValue(MANNING);
    }
    manning /= rGeometry.size();
    mManning2 = std::pow(manning, 2);

    mEpsilon = rGeometry.Length() * rProcessInfo[RELATIVE_DRY_HEIGHT];
}

double ManningLaw::CalculateLHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    const double inv_height = ShallowWaterUtilities().InverseHeight(rHeight, mEpsilon);
    return mManning2 * norm_2(rVelocity) * std::pow(inv_height, 4.0 / 3.0);
}

array_1d<double,3> ManningLaw::CalculateRHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    return rVelocity * CalculateLHS(rHeight, rVelocity);
}

}  // namespace Kratos
