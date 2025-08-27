// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//                   Anne van de Graaf
//
#pragma once

#include "geo_aliases.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class Element;

class IntegrationCoefficientModifier
{
public:
    virtual ~IntegrationCoefficientModifier()                                           = default;
    [[nodiscard]] virtual std::unique_ptr<IntegrationCoefficientModifier> Clone() const = 0;

    virtual double operator()(double                           IntegrationCoefficient,
                              const Geo::IntegrationPointType& rIntegrationPoint,
                              const Element&                   rElement) const = 0;
};

class KRATOS_API(GEO_MECHANICS_APPLICATION) IntegrationCoefficientsCalculator
{
public:
    explicit IntegrationCoefficientsCalculator(std::unique_ptr<IntegrationCoefficientModifier> = nullptr);

    [[nodiscard]] std::unique_ptr<IntegrationCoefficientModifier> CloneModifier() const;

    template <typename OutputContainer = std::vector<double>>
    OutputContainer Run(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                        const Vector&                          rDetJs,
                        const Element*                         pElement = nullptr) const
    {
        auto result                            = OutputContainer(rIntegrationPoints.size());
        auto calculate_integration_coefficient = [](const auto& rIntegrationPoint, auto DetJ) {
            return rIntegrationPoint.Weight() * DetJ;
        };
        std::ranges::transform(rIntegrationPoints, rDetJs, result.begin(), calculate_integration_coefficient);

        if (mCoefficientModifier && pElement) {
            auto apply_modifier = [&element = *pElement, &modifier = *mCoefficientModifier](
                                      auto IntegrationCoefficient, const auto& rIntegrationPoint) {
                return modifier(IntegrationCoefficient, rIntegrationPoint, element);
            };
            std::ranges::transform(result, rIntegrationPoints, result.begin(), apply_modifier);
        }

        return result;
    }

private:
    std::unique_ptr<IntegrationCoefficientModifier> mCoefficientModifier;

    friend class Serializer;
    void save(const Serializer& rSerializer) const;
    void load(const Serializer& rSerializer) const;
};

} // namespace Kratos