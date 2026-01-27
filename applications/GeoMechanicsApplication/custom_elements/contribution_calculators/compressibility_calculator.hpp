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

#pragma once

#include "contribution_calculator.h"
#include "custom_retention/retention_law.h"
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
class CompressibilityCalculator : public ContributionCalculator<TNumNodes>
{
public:
    struct InputProvider {
        InputProvider(Geo::PropertiesGetter                          GetElementProperties,
                      Geo::RetentionLawsGetter                       GetRetentionLaws,
                      std::function<const Matrix&()>                 GetNContainer,
                      Geo::IntegrationCoefficientsGetter             GetIntegrationCoefficients,
                      std::function<double()>                        GetMatrixScalarFactor,
                      std::function<Vector(const Variable<double>&)> GetNodalValuesOf,
                      std::function<std::vector<double>()>           GetFluidPressures)
            : GetElementProperties(std::move(GetElementProperties)),
              GetRetentionLaws(std::move(GetRetentionLaws)),
              GetNContainer(std::move(GetNContainer)),
              GetIntegrationCoefficients(std::move(GetIntegrationCoefficients)),
              GetMatrixScalarFactor(std::move(GetMatrixScalarFactor)),
              GetNodalValues(std::move(GetNodalValuesOf)),
              GetFluidPressures(std::move(GetFluidPressures))
        {
        }

        Geo::PropertiesGetter                          GetElementProperties;
        Geo::RetentionLawsGetter                       GetRetentionLaws;
        std::function<const Matrix&()>                 GetNContainer;
        Geo::IntegrationCoefficientsGetter             GetIntegrationCoefficients;
        std::function<double()>                        GetMatrixScalarFactor;
        std::function<Vector(const Variable<double>&)> GetNodalValues;
        std::function<std::vector<double>()>           GetFluidPressures;
    };

    explicit CompressibilityCalculator(InputProvider rInputProvider)
        : mInputProvider(std::move(rInputProvider))
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
        const auto& fluid_pressures          = mInputProvider.GetFluidPressures();
        BoundedMatrix<double, TNumNodes, TNumNodes> result = ZeroMatrix(TNumNodes, TNumNodes);
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

    [[nodiscard]] double CalculateBiotModulusInverse(const RetentionLaw::Pointer& rRetentionLaw,
                                                     double FluidPresssure) const
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
        retention_parameters.SetFluidPressure(FluidPresssure);
        result *= rRetentionLaw->CalculateSaturation(retention_parameters);
        result -= rRetentionLaw->CalculateDerivativeOfSaturation(retention_parameters) * r_properties[POROSITY];
        return result;
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
