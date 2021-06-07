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
#include "chezy_law.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"


namespace Kratos
{

void ChezyLaw::Initialize(const GeometryType& rGeometry, const ProcessInfo& rProcessInfo)
{
    double chezy = 0.0;
    for (auto& r_node : rGeometry)
    {
        chezy += r_node.FastGetSolutionStepValue(CHEZY);
    }
    chezy /= rGeometry.size();
    mCoefficient = 1. / std::pow(chezy, 2);

    mEpsilon = rGeometry.Length() * rProcessInfo[RELATIVE_DRY_HEIGHT];
}

double ChezyLaw::CalculateLHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    const double inv_height = ShallowWaterUtilities().InverseHeight(rHeight, mEpsilon);
    return mCoefficient * norm_2(rVelocity) * inv_height;
}

array_1d<double,3> ChezyLaw::CalculateRHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    return rVelocity * CalculateLHS(rHeight, rVelocity);
}

}  // namespace Kratos
