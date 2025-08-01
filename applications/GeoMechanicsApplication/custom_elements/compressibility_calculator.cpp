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

std::optional<Matrix> CompressibilityCalculator::LHSContribution()
{
    return std::make_optional(LHSContribution(CalculateCompressibilityMatrix()));
}

Vector CompressibilityCalculator::RHSContribution()
{
    return RHSContribution(CalculateCompressibilityMatrix());
}

Vector CompressibilityCalculator::RHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return -prod(rCompressibilityMatrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE));
}

Matrix CompressibilityCalculator::LHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return mInputProvider.GetMatrixScalarFactor() * rCompressibilityMatrix;
}

std::pair<std::optional<Matrix>, Vector> CompressibilityCalculator::LocalSystemContribution()
{
    const auto compressibility_matrix = CalculateCompressibilityMatrix();
    return {std::make_optional(LHSContribution(compressibility_matrix)), RHSContribution(compressibility_matrix)};
}

Matrix CompressibilityCalculator::CalculateCompressibilityMatrix() const
{
    const auto& r_N_container            = mInputProvider.GetNContainer();
    const auto& integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    const auto& fluid_pressures          = mInputProvider.GetFluidPressures();
    Matrix      result(r_N_container.size2(), r_N_container.size2(), 0.0);
    Vector      N(r_N_container.size2());

    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        noalias(N) = row(r_N_container, integration_point_index);
        noalias(result) += GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            N,
            CalculateBiotModulusInverse(mInputProvider.GetRetentionLaws()[integration_point_index],
                                        fluid_pressures[integration_point_index]),
            integration_coefficients[integration_point_index]);
    }
    return result;
}

double CompressibilityCalculator::CalculateBiotModulusInverse(const RetentionLaw::Pointer& rRetentionLaw,
                                                              double FluidPressure) const
{
    const auto&  r_properties     = mInputProvider.GetElementProperties();
    const double biot_coefficient = r_properties[BIOT_COEFFICIENT];

    double bulk_fluid = TINY;
    if (!r_properties[IGNORE_UNDRAINED]) {
        bulk_fluid = r_properties[BULK_MODULUS_FLUID];
    }
    double result = (biot_coefficient - r_properties[POROSITY]) / r_properties[BULK_MODULUS_SOLID] +
                    r_properties[POROSITY] / bulk_fluid;

    RetentionLaw::Parameters retention_parameters(r_properties);
    retention_parameters.SetFluidPressure(FluidPressure);
    result *= rRetentionLaw->CalculateSaturation(retention_parameters);
    result -= rRetentionLaw->CalculateDerivativeOfSaturation(retention_parameters) * r_properties[POROSITY];
    return result;
}

} // namespace Kratos
