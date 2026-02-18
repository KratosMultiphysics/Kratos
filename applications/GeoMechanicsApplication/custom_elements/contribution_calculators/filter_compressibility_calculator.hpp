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

#pragma once

#include "contribution_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_aliases.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

#include <utility>
#include <vector>

namespace Kratos
{
template <unsigned int TNumNodes>
class FilterCompressibilityCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
        InputProvider(Geo::PropertiesGetter                GetElementProperties,
                      std::function<const Matrix&()>       GetNContainer,
                      Geo::IntegrationCoefficientsGetter   GetIntegrationCoefficients,
                      std::function<std::vector<Vector>()> GetProjectedGravityAtIntegrationPoints,
                      std::function<double()>              GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf)
            : GetElementProperties(std::move(GetElementProperties)),
              GetNContainer(std::move(GetNContainer)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetProjectedGravityAtIntegrationPoints(std::move(GetProjectedGravityAtIntegrationPoints)),
              GetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              GetNodalValues(std::move(GetNodalValuesOf))
        {
        }

        Geo::PropertiesGetter                          GetElementProperties;
        std::function<const Matrix&()>                 GetNContainer;
        Geo::IntegrationCoefficientsGetter             GetIntegrationCoefficients;
        std::function<std::vector<Vector>()>           GetProjectedGravityAtIntegrationPoints;
        std::function<double()>                        GetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)> GetNodalValues;
    };

    explicit FilterCompressibilityCalculator(InputProvider AnInputProvider)
        : mInputProvider(std::move(AnInputProvider))
    {
    }

    std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>> LHSContribution() override
    {
        return std::make_optional(LHSContribution(CalculateCompressibilityMatrix()));
    }

    BoundedVector<double, TNumNodes> RHSContribution() override
    {
        return RHSContribution(CalculateCompressibilityMatrix());
    }

    std::pair<std::optional<BoundedMatrix<double, TNumNodes, TNumNodes>>, BoundedVector<double, TNumNodes>> LocalSystemContribution() override
    {
        const auto compressibility_matrix = CalculateCompressibilityMatrix();
        return {std::make_optional(LHSContribution(compressibility_matrix)),
                RHSContribution(compressibility_matrix)};
    }

private:
    [[nodiscard]] BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCompressibilityMatrix() const
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

    [[nodiscard]] double CalculateElasticCapacity(double ProjectedGravity) const
    {
        const auto& r_properties = mInputProvider.GetElementProperties();
        return 1.0 / (r_properties[DENSITY_WATER] * ProjectedGravity * r_properties[FILTER_LENGTH]) +
               1.0 / r_properties[BULK_MODULUS_FLUID];
    }

    [[nodiscard]] BoundedVector<double, TNumNodes> RHSContribution(
        const BoundedMatrix<double, TNumNodes, TNumNodes>& rCompressibilityMatrix) const
    {
        return -prod(rCompressibilityMatrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE));
    }

    [[nodiscard]] BoundedMatrix<double, TNumNodes, TNumNodes> LHSContribution(const Matrix& rCompressibilityMatrix) const
    {
        return mInputProvider.GetMatrixScalarFactor() * rCompressibilityMatrix;
    }

    InputProvider mInputProvider;
};

} // namespace Kratos
