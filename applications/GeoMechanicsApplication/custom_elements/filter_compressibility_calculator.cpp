// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Richard Faasse
//
#include "filter_compressibility_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{
template <unsigned int TNumNodes>
FilterCompressibilityCalculator<TNumNodes>::FilterCompressibilityCalculator(InputProvider AnInputProvider)
    : mInputProvider(std::move(AnInputProvider))
{
}

template <unsigned int TNumNodes>
std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> FilterCompressibilityCalculator<TNumNodes>::LHSContribution()
{
    return std::make_optional(LHSContribution(CalculateCompressibilityMatrix()));
}

template <unsigned int TNumNodes>
BoundedVector<double, TNumNodes> FilterCompressibilityCalculator<TNumNodes>::RHSContribution()
{
    return RHSContribution(CalculateCompressibilityMatrix());
}

template <unsigned int TNumNodes>
BoundedVector<double, TNumNodes> FilterCompressibilityCalculator<TNumNodes>::RHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return -prod(rCompressibilityMatrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE));
}

template <unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> FilterCompressibilityCalculator<TNumNodes>::LHSContribution(
    const Matrix& rCompressibilityMatrix) const
{
    return mInputProvider.GetMatrixScalarFactor() * rCompressibilityMatrix;
}

template <unsigned int TNumNodes>
std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> FilterCompressibilityCalculator<TNumNodes>::LocalSystemContribution()
{
    const auto compressibility_matrix = CalculateCompressibilityMatrix();
    return {std::make_optional(LHSContribution(compressibility_matrix)), RHSContribution(compressibility_matrix)};
}

template <unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> FilterCompressibilityCalculator<TNumNodes>::CalculateCompressibilityMatrix() const
{
    const auto& r_N_container            = mInputProvider.GetNContainer();
    const auto& integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    const auto& projected_gravity_on_integration_points =
        mInputProvider.GetProjectedGravityAtIntegrationPoints();
    auto result = Matrix{ZeroMatrix{r_N_container.size2(), r_N_container.size2()}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        const auto N = Vector{row(r_N_container, integration_point_index)};
        result += GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            N, CalculateElasticCapacity(projected_gravity_on_integration_points[integration_point_index][0]),
            integration_coefficients[integration_point_index]);
    }
    return result;
}

template <unsigned int TNumNodes>
double FilterCompressibilityCalculator<TNumNodes>::CalculateElasticCapacity(double ProjectedGravity) const
{
    const auto& r_properties = mInputProvider.GetElementProperties();
    return 1.0 / (r_properties[DENSITY_WATER] * ProjectedGravity * r_properties[FILTER_LENGTH]) +
           1.0 / r_properties[BULK_MODULUS_FLUID];
}

} // namespace Kratos
