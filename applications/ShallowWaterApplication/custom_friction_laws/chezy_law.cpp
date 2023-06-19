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
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

ChezyLaw::ChezyLaw(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    this->Initialize(rGeometry, rProperty, rProcessInfo);
}

void ChezyLaw::Initialize(
    const GeometryType& rGeometry,
    const Properties& rProperty,
    const ProcessInfo& rProcessInfo)
{
    const double chezy = rProperty.GetValue(CHEZY);
    mCoefficient = 1. / std::pow(chezy, 2);
    mEpsilon = rGeometry.Length() * rProcessInfo[RELATIVE_DRY_HEIGHT];
}

double ChezyLaw::CalculateLHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    const double inv_height = PhaseFunction::InverseHeight(rHeight, mEpsilon);
    return mCoefficient * norm_2(rVelocity) * inv_height;
}

array_1d<double,3> ChezyLaw::CalculateRHS(const double& rHeight, const array_1d<double,3>& rVelocity)
{
    return rVelocity * CalculateLHS(rHeight, rVelocity);
}

}  // namespace Kratos
