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
//                   Ruben Zorrilla
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
    const std::size_t TimeOrder,
    Kratos::unique_ptr<BDF> &rpAuxBDF)
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

void BDF::ComputeAndSaveBDFCoefficients(ProcessInfo &rProcessInfo) const
{
    // Check if the auxiliary BDF util pointer is set
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;

    // Compute the BDF coefficients
    const auto bdf_coeffs = this->ComputeBDFCoefficients(rProcessInfo);

    // Check ProcessInfo BDF coefficients vector size
    const unsigned int n_coefs = bdf_coeffs.size();
    auto &r_proc_inf_bdf_coeffs = rProcessInfo[BDF_COEFFICIENTS];
    if (r_proc_inf_bdf_coeffs.size() != n_coefs) {
        r_proc_inf_bdf_coeffs.resize(n_coefs);
    }

    // Save the computed BDF coefficients in the model part ProcessInfo
    for (std::size_t i_coeff = 0; i_coeff < n_coefs; ++i_coeff) {
        r_proc_inf_bdf_coeffs[i_coeff] = bdf_coeffs[i_coeff];
    }
}

std::vector<double> BDF::ComputeBDFCoefficients(double DeltaTime) const
{
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;
    return mpAuxBDF->ComputeBDFCoefficients(DeltaTime);
}

std::vector<double> BDF::ComputeBDFCoefficients(double DeltaTime, double PreviousDeltaTime) const
{
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;
    return mpAuxBDF->ComputeBDFCoefficients(DeltaTime, PreviousDeltaTime);
}

std::vector<double> BDF::ComputeBDFCoefficients(const ProcessInfo &rProcessInfo) const
{
    KRATOS_ERROR_IF(!mpAuxBDF)
        << "Pointer to auxiliary BDF class implementing the desired order is null" << std::endl;
    return mpAuxBDF->ComputeBDFCoefficients(rProcessInfo);
}

std::size_t BDF::GetTimeOrder() const
{
    return mTimeOrder;
}

std::vector<double> BDF1::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    std::vector<double> coefficients(2);

    coefficients[0] =  1.0/DeltaTime;
    coefficients[1] = -1.0/DeltaTime;

    return coefficients;
}

std::vector<double> BDF1::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

std::vector<double> BDF2::ComputeBDFCoefficients(const double DeltaTime, const double PreviousDeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double rho = PreviousDeltaTime / DeltaTime;
    double time_coeff = 1.0 / (DeltaTime * rho * rho + DeltaTime * rho);

    std::vector<double> coefficients(3);

    if (PreviousDeltaTime < std::numeric_limits<double>::epsilon()) {
        KRATOS_DETAIL("ComputeBDFCoefficients") << "previous delta-time is zero, using "
            << "BDF1 Coefficients \n";
        coefficients[0] =  1.0 / DeltaTime;
        coefficients[1] = -1.0 / DeltaTime;
        coefficients[2] =  0.0;
    }
    else {
        coefficients[0] =  time_coeff * (rho * rho + 2.0 * rho); // coefficient for step n+1 (3/2Dt if Dt is constant)
        coefficients[1] = -time_coeff * (rho * rho + 2.0 * rho + 1.0); // coefficient for step n (-4/2Dt if Dt is constant)
        coefficients[2] =  time_coeff; // coefficient for step n-1 (1/2Dt if Dt is constant)
    }

    return coefficients;
}

std::vector<double> BDF2::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME],
                                  rProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME]);
}

std::vector<double> BDF3::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 6.0*DeltaTime;

    std::vector<double> coefficients(4);

    coefficients[0] =  11.0 / denom;
    coefficients[1] = -18.0 / denom;
    coefficients[2] =   9.0 / denom;
    coefficients[3] =  -2.0 / denom;

    return coefficients;
}

std::vector<double> BDF3::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

std::vector<double> BDF4::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 12.0*DeltaTime;

    std::vector<double> coefficients(5);

    coefficients[0] =  25.0 / denom;
    coefficients[1] = -48.0 / denom;
    coefficients[2] =  36.0 / denom;
    coefficients[3] = -16.0 / denom;
    coefficients[4] =   3.0 / denom;

    return coefficients;
}

std::vector<double> BDF4::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

std::vector<double> BDF5::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 60.0*DeltaTime;

    std::vector<double> coefficients(6);

    coefficients[0] =  137.0 / denom;
    coefficients[1] = -300.0 / denom;
    coefficients[2] =  300.0 / denom;
    coefficients[3] = -200.0 / denom;
    coefficients[4] =   75.0 / denom;
    coefficients[5] =  -12.0 / denom;

    return coefficients;
}

std::vector<double> BDF5::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

std::vector<double> BDF6::ComputeBDFCoefficients(const double DeltaTime) const
{
    KRATOS_ERROR_IF(DeltaTime < std::numeric_limits<double>::epsilon())
        << "Expects DeltaTime > 0!" << std::endl;

    const double denom = 60.0*DeltaTime;

    std::vector<double> coefficients(7);

    coefficients[0] =  147.0 / denom;
    coefficients[1] = -360.0 / denom;
    coefficients[2] =  450.0 / denom;
    coefficients[3] = -400.0 / denom;
    coefficients[4] =  225.0 / denom;
    coefficients[5] =  -72.0 / denom;
    coefficients[6] =   10.0 / denom;

    return coefficients;
}

std::vector<double> BDF6::ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DELTA_TIME)) << "No DELTA_TIME "
        << "defined in the ProcessInfo!" << std::endl;
    return ComputeBDFCoefficients(rProcessInfo[DELTA_TIME]);
}

} // namespace TimeDiscretization.
}  // namespace Kratos.
