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

// External includes

// Project includes
#include "time_discretization.h"

namespace Kratos {
namespace TimeDiscretization {

void BDF1::ComputeBDFCoefficients(const double DeltaTime, std::vector<double>& rCoefficients)
{

    if (rCoefficients.size() != 2) rCoefficients.resize(2);

    rCoefficients[0] =  1.0/DeltaTime;
    rCoefficients[1] = -1.0/DeltaTime;
}

void BDF2::ComputeBDFCoefficients(const double DeltaTime, std::vector<double>& rCoefficients)
{
    if (rCoefficients.size() != 3) rCoefficients.resize(3);

    const double denom = 1.0/(2.0*DeltaTime);

    rCoefficients[0] =  3.0 / denom;
    rCoefficients[1] = -4.0 / denom;
    rCoefficients[2] =  2.0 / denom;

}

void BDF3::ComputeBDFCoefficients(const double DeltaTime, std::vector<double>& rCoefficients)
{
    if (rCoefficients.size() != 4) rCoefficients.resize(4);

    const double denom = 1.0/(6.0*DeltaTime);

    rCoefficients[0] =  11.0 / denom;
    rCoefficients[1] = -18.0 / denom;
    rCoefficients[2] =   9.0 / denom;
    rCoefficients[3] =  -2.0 / denom;
}

void BDF4::ComputeBDFCoefficients(const double DeltaTime, std::vector<double>& rCoefficients)
{
    if (rCoefficients.size() != 5) rCoefficients.resize(5);

    const double denom = 1.0/(12.0*DeltaTime);

    rCoefficients[0] =  25.0 / denom;
    rCoefficients[1] = -48.0 / denom;
    rCoefficients[2] =  36.0 / denom;
    rCoefficients[3] = -16.0 / denom;
    rCoefficients[4] =   3.0 / denom;
}

void BDF5::ComputeBDFCoefficients(const double DeltaTime, std::vector<double>& rCoefficients)
{
    if (rCoefficients.size() != 6) rCoefficients.resize(6);

    const double denom = 1.0/(60.0*DeltaTime);

    rCoefficients[0] =  137.0 / denom;
    rCoefficients[1] = -300.0 / denom;
    rCoefficients[2] =  300.0 / denom;
    rCoefficients[3] = -200.0 / denom;
    rCoefficients[4] =   75.0 / denom;
    rCoefficients[5] =  -12.0 / denom;
}

void BDF6::ComputeBDFCoefficients(const double DeltaTime, std::vector<double>& rCoefficients)
{
    if (rCoefficients.size() != 7) rCoefficients.resize(7);

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


