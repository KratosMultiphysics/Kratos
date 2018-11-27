//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "includes/checks.h"
#include "time_discretization.h"

namespace Kratos {
namespace TimeDiscretization {

void BDF1::ComputeBDFCoefficients(const double DeltaTime, std::array<double, 2>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    rCoefficients[0] =  1.0/DeltaTime;
    rCoefficients[1] = -1.0/DeltaTime;
}

void BDF2::ComputeBDFCoefficients(const double DeltaTime, std::array<double, 3>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 1.0/(2.0*DeltaTime);

    rCoefficients[0] =  3.0 / denom;
    rCoefficients[1] = -4.0 / denom;
    rCoefficients[2] =  2.0 / denom;
}

void BDF2::ComputeBDFCoefficients(const double DeltaTime, const double PreviousDeltaTime, std::array<double, 3>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;
    KRATOS_ERROR_IF(PreviousDeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects PreviousDeltaTime > 0!" << std::endl;

    const double rho = PreviousDeltaTime / DeltaTime;
    double time_coeff = 1.0 / (DeltaTime * rho * rho + DeltaTime * rho);

    rCoefficients[0] =  time_coeff * (rho * rho + 2.0 * rho); // coefficient for step n+1 (3/2Dt if Dt is constant)
    rCoefficients[1] = -time_coeff * (rho * rho + 2.0 * rho + 1.0); // coefficient for step n (-4/2Dt if Dt is constant)
    rCoefficients[2] =  time_coeff; // coefficient for step n-1 (1/2Dt if Dt is constant)
}

void BDF3::ComputeBDFCoefficients(const double DeltaTime, std::array<double, 4>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 1.0/(6.0*DeltaTime);

    rCoefficients[0] =  11.0 / denom;
    rCoefficients[1] = -18.0 / denom;
    rCoefficients[2] =   9.0 / denom;
    rCoefficients[3] =  -2.0 / denom;
}

void BDF4::ComputeBDFCoefficients(const double DeltaTime, std::array<double, 5>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 1.0/(12.0*DeltaTime);

    rCoefficients[0] =  25.0 / denom;
    rCoefficients[1] = -48.0 / denom;
    rCoefficients[2] =  36.0 / denom;
    rCoefficients[3] = -16.0 / denom;
    rCoefficients[4] =   3.0 / denom;
}

void BDF5::ComputeBDFCoefficients(const double DeltaTime, std::array<double, 6>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 1.0/(60.0*DeltaTime);

    rCoefficients[0] =  137.0 / denom;
    rCoefficients[1] = -300.0 / denom;
    rCoefficients[2] =  300.0 / denom;
    rCoefficients[3] = -200.0 / denom;
    rCoefficients[4] =   75.0 / denom;
    rCoefficients[5] =  -12.0 / denom;
}

void BDF6::ComputeBDFCoefficients(const double DeltaTime, std::array<double, 7>& rCoefficients) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 1.0/(60.0*DeltaTime);

    rCoefficients[0] =  147.0 / denom;
    rCoefficients[1] = -360.0 / denom;
    rCoefficients[2] =  450.0 / denom;
    rCoefficients[3] = -400.0 / denom;
    rCoefficients[4] =  225.0 / denom;
    rCoefficients[5] =  -72.0 / denom;
    rCoefficients[6] =   10.0 / denom;
}

} // namespace TimeDiscretization.
}  // namespace Kratos.


