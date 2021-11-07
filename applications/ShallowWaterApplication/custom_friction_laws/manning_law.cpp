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
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

ManningLaw::ManningLaw(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    this->Initialize(rGeometry, rProperty, rProcessInfo);
}

void ManningLaw::Initialize(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    const double manning = rProperty.GetValue(MANNING);
    mManning2 = std::pow(manning, 2);
    mEpsilon = rGeometry.Length() * rProcessInfo[RELATIVE_DRY_HEIGHT];
}

double ManningLaw::CalculateLHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    const double inv_height = PhaseFunction::InverseHeight(rHeight, mEpsilon);
    return mManning2 * norm_2(rVelocity) * std::pow(inv_height, 4.0 / 3.0);
}

array_1d<double,3> ManningLaw::CalculateRHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    return rVelocity * CalculateLHS(rHeight, rVelocity);
}

}  // namespace Kratos
