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
#include "includes/process_info.h"
#include "includes/variables.h"
#include "input_output/logger.h"
#include "time_discretization.h"

namespace Kratos {
namespace TimeDiscretization {

void BDF::SetAuxBDFPointer(
    const unsigned int TimeOrder,
    unique_ptr<BDF> &rpAuxBDF)
{
    switch (TimeOrder) {
        case 1: {
            rpAuxBDF = Kratos::make_unique<BDF1>();
            break;
        } case 2: {
            rpAuxBDF = Kratos::make_unique<BDF2>();
            break;
        } case 3: {
            rpAuxBDF = Kratos::make_unique<BDF3>();
            break;
        } case 4: {
            rpAuxBDF = Kratos::make_unique<BDF4>();
            break;
        } case 5: {
            rpAuxBDF = Kratos::make_unique<BDF5>();
            break;
        } case 6: {
            rpAuxBDF = Kratos::make_unique<BDF6>();
            break;
        } default: {
            KRATOS_ERROR << "Asked for time order " << TimeOrder << ". Maximum time order is 6.";
        }
    }
}

Vector BDF::ComputeBDFCoefficients(double DeltaTime) const
{
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;
    return mpAuxBDF->ComputeBDFCoefficients(DeltaTime);
}

Vector BDF::ComputeBDFCoefficients(double DeltaTime, double PreviousDeltaTime) const
{
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;
    return mpAuxBDF->ComputeBDFCoefficients(DeltaTime, PreviousDeltaTime);
}

Vector BDF::ComputeBDFCoefficients(const ProcessInfo &rProcessInfo) const
{
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;
    return mpAuxBDF->ComputeBDFCoefficients(rProcessInfo);
}

const unsigned int BDF::GetTimeOrder() const
{
    return mTimeOrder;
}

Vector BDF1::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    Vector coefficients(2);

    coefficients[0] =  1.0/DeltaTime;
    coefficients[1] = -1.0/DeltaTime;

    return coefficients;
}

Vector BDF1::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

Vector BDF2::ComputeBDFCoefficients(const double DeltaTime, const double PreviousDeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double rho = PreviousDeltaTime / DeltaTime;
    double time_coeff = 1.0 / (DeltaTime * rho * rho + DeltaTime * rho);

    Vector coefficients(3);

    if (PreviousDeltaTime < std::numeric_limits<double>::epsilon()) {
        KRATOS_DETAIL("ComputeBDFCoefficients") << "previous delta-time is zero, using "
            << "constant time-step for computation of coefficients" << std::endl;
        coefficients[0] =  1.5 / DeltaTime;
        coefficients[1] = -2.0 / DeltaTime;
        coefficients[2] =  0.5 / DeltaTime;
    }
    else {
        coefficients[0] =  time_coeff * (rho * rho + 2.0 * rho); // coefficient for step n+1 (3/2Dt if Dt is constant)
        coefficients[1] = -time_coeff * (rho * rho + 2.0 * rho + 1.0); // coefficient for step n (-4/2Dt if Dt is constant)
        coefficients[2] =  time_coeff; // coefficient for step n-1 (1/2Dt if Dt is constant)
    }

    return coefficients;
}

Vector BDF2::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME],
                                  rProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME]);
}

Vector BDF3::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 6.0*DeltaTime;

    Vector coefficients(4);

    coefficients[0] =  11.0 / denom;
    coefficients[1] = -18.0 / denom;
    coefficients[2] =   9.0 / denom;
    coefficients[3] =  -2.0 / denom;

    return coefficients;
}

Vector BDF3::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

Vector BDF4::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 12.0*DeltaTime;

    Vector coefficients(5);

    coefficients[0] =  25.0 / denom;
    coefficients[1] = -48.0 / denom;
    coefficients[2] =  36.0 / denom;
    coefficients[3] = -16.0 / denom;
    coefficients[4] =   3.0 / denom;

    return coefficients;
}

Vector BDF4::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

Vector BDF5::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 60.0*DeltaTime;

    Vector coefficients(6);

    coefficients[0] =  137.0 / denom;
    coefficients[1] = -300.0 / denom;
    coefficients[2] =  300.0 / denom;
    coefficients[3] = -200.0 / denom;
    coefficients[4] =   75.0 / denom;
    coefficients[5] =  -12.0 / denom;

    return coefficients;
}

Vector BDF5::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

Vector BDF6::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 60.0*DeltaTime;

    Vector coefficients(7);

    coefficients[0] =  147.0 / denom;
    coefficients[1] = -360.0 / denom;
    coefficients[2] =  450.0 / denom;
    coefficients[3] = -400.0 / denom;
    coefficients[4] =  225.0 / denom;
    coefficients[5] =  -72.0 / denom;
    coefficients[6] =   10.0 / denom;

    return coefficients;
}

Vector BDF6::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

} // namespace TimeDiscretization.
}  // namespace Kratos.
