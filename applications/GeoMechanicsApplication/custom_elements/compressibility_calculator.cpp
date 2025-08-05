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
template <unsigned int TNumNodes>
CompressibilityCalculator<TNumNodes>::CompressibilityCalculator(InputProvider rInputProvider)
    : mInputProvider(std::move(rInputProvider))
{
}

template <unsigned int TNumNodes>
std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> CompressibilityCalculator<TNumNodes>::LHSContribution()
{
    return std::make_optional(LHSContribution(CalculateCompressibilityMatrix()));
}

template <unsigned int TNumNodes>
BoundedVector<double, TNumNodes> CompressibilityCalculator<TNumNodes>::RHSContribution()
{
    return RHSContribution(CalculateCompressibilityMatrix());
}

template <unsigned int TNumNodes>
BoundedVector<double, TNumNodes> CompressibilityCalculator<TNumNodes>::RHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return -prod(rCompressibilityMatrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE));
}

template <unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> CompressibilityCalculator<TNumNodes>::LHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return mInputProvider.GetMatrixScalarFactor() * rCompressibilityMatrix;
}

template <unsigned int TNumNodes>
std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> CompressibilityCalculator<TNumNodes>::LocalSystemContribution()
{
    const auto compressibility_matrix = CalculateCompressibilityMatrix();
    return {std::make_optional(LHSContribution(compressibility_matrix)), RHSContribution(compressibility_matrix)};
}

template <unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> CompressibilityCalculator<TNumNodes>::CalculateCompressibilityMatrix() const
{
    const auto& r_N_container            = mInputProvider.GetNContainer();
    const auto& integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    const auto& fluid_pressures          = mInputProvider.GetFluidPressures();
    BoundedMatrix<double, TNumNodes, TNumNodes> result;
    BoundedVector<double, TNumNodes>            N;

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

template <unsigned int TNumNodes>
double CompressibilityCalculator<TNumNodes>::CalculateBiotModulusInverse(const RetentionLaw::Pointer& rRetentionLaw,
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
