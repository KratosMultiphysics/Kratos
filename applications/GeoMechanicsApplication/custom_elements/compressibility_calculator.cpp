// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "compressibility_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

CompressibilityCalculator::CompressibilityCalculator(InputProvider rInputProvider)
    : mInputProvider(std::move(rInputProvider))
{
}

Matrix CompressibilityCalculator::LHSContribution() { return {}; }

Vector CompressibilityCalculator::RHSContribution() { return {}; }

std::pair<Matrix, Vector> CompressibilityCalculator::CalculateLeftAndRightHandSide()
{
    auto compressibility_matrix = CalculateCompressibilityMatrix();
    return {compressibility_matrix * mInputProvider.GetDtPressureCoefficient(), -prod(compressibility_matrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE))};
}

Matrix CompressibilityCalculator::CalculateCompressibilityMatrix() const
{
    const auto& r_N_container              = mInputProvider.GetNContainer();
    const auto& r_integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    auto        result = Matrix{ZeroMatrix{r_N_container.size2(), r_N_container.size2()}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < r_integration_coefficients.size(); ++integration_point_index) {
        const auto   N = Vector{row(r_N_container, integration_point_index)};
        const double BiotModulusInverse =
            CalculateBiotModulusInverse(mInputProvider.GetRetentionLaws()[integration_point_index]);
        result += GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            N, BiotModulusInverse, r_integration_coefficients[integration_point_index]);
    }
    return result;
}

double CompressibilityCalculator::CalculateBiotModulusInverse(const RetentionLaw::Pointer& rRetentionLaw) const
{
    const auto&  r_properties     = mInputProvider.GetElementProperties();
    const double biot_coefficient = r_properties[BIOT_COEFFICIENT];

    double bulk_fluid = TINY;
    if (!r_properties[IGNORE_UNDRAINED]) {
        bulk_fluid = r_properties[BULK_MODULUS_FLUID];
    }
    double result = (biot_coefficient - r_properties[POROSITY]) / r_properties[BULK_MODULUS_SOLID] +
                    r_properties[POROSITY] / bulk_fluid;

    RetentionLaw::Parameters retention_parameters(mInputProvider.GetElementProperties());
    const double degree_of_saturation = rRetentionLaw->CalculateSaturation(retention_parameters);
    const double derivative_of_saturation = rRetentionLaw->CalculateDerivativeOfSaturation(retention_parameters);

    result *= degree_of_saturation;
    result -= derivative_of_saturation * r_properties[POROSITY];
    return result;
}

} // namespace Kratos
