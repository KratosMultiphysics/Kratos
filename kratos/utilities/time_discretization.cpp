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

std::array<double, 2> BDF1::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    std::array<double, 2> coefficients;

    coefficients[0] =  1.0/DeltaTime;
    coefficients[1] = -1.0/DeltaTime;

    return coefficients;
}

std::array<double, 3> BDF2::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 2.0*DeltaTime;

    std::array<double, 3> coefficients;

    coefficients[0] =  3.0 / denom;
    coefficients[1] = -4.0 / denom;
    coefficients[2] =  2.0 / denom;

    return coefficients;
}

std::array<double, 3> BDF2::ComputeBDFCoefficients(const double DeltaTime, const double PreviousDeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;
    KRATOS_ERROR_IF(PreviousDeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects PreviousDeltaTime > 0!" << std::endl;

    const double rho = PreviousDeltaTime / DeltaTime;
    double time_coeff = 1.0 / (DeltaTime * rho * rho + DeltaTime * rho);

    std::array<double, 3> coefficients;

    coefficients[0] =  time_coeff * (rho * rho + 2.0 * rho); // coefficient for step n+1 (3/2Dt if Dt is constant)
    coefficients[1] = -time_coeff * (rho * rho + 2.0 * rho + 1.0); // coefficient for step n (-4/2Dt if Dt is constant)
    coefficients[2] =  time_coeff; // coefficient for step n-1 (1/2Dt if Dt is constant)

    return coefficients;
}

std::array<double, 4> BDF3::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 6.0*DeltaTime;

    std::array<double, 4> coefficients;

    coefficients[0] =  11.0 / denom;
    coefficients[1] = -18.0 / denom;
    coefficients[2] =   9.0 / denom;
    coefficients[3] =  -2.0 / denom;

    return coefficients;
}

std::array<double, 5> BDF4::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 12.0*DeltaTime;

    std::array<double, 5> coefficients;

    coefficients[0] =  25.0 / denom;
    coefficients[1] = -48.0 / denom;
    coefficients[2] =  36.0 / denom;
    coefficients[3] = -16.0 / denom;
    coefficients[4] =   3.0 / denom;

    return coefficients;
}

std::array<double, 6> BDF5::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 60.0*DeltaTime;

    std::array<double, 6> coefficients;

    coefficients[0] =  137.0 / denom;
    coefficients[1] = -300.0 / denom;
    coefficients[2] =  300.0 / denom;
    coefficients[3] = -200.0 / denom;
    coefficients[4] =   75.0 / denom;
    coefficients[5] =  -12.0 / denom;

    return coefficients;
}

std::array<double, 7> BDF6::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 60.0*DeltaTime;

    std::array<double, 7> coefficients;

    coefficients[0] =  147.0 / denom;
    coefficients[1] = -360.0 / denom;
    coefficients[2] =  450.0 / denom;
    coefficients[3] = -400.0 / denom;
    coefficients[4] =  225.0 / denom;
    coefficients[5] =  -72.0 / denom;
    coefficients[6] =   10.0 / denom;

    return coefficients;
}

} // namespace TimeDiscretization.
}  // namespace Kratos.


