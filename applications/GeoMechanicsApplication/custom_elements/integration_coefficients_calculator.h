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
//
#pragma once

#include "geo_aliases.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class Element;
class Serializer;

class KRATOS_API(GEO_MECHANICS_APPLICATION) IntegrationCoefficientsCalculator
{
public:
    virtual ~IntegrationCoefficientsCalculator() = default;

    [[nodiscard]] virtual Vector CalculateIntegrationCoefficients(const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
                                                                  const Vector& rDetJs,
                                                                  double        CrossArea,
                                                                  std::size_t LocalDimension = 0) const
    {
        KRATOS_ERROR
            << "IntegrationCoefficientsCalculator::CalculateIntegrationCoefficients is called."
            << std::endl;
    }

    [[nodiscard]] virtual std::vector<double> CalculateIntegrationCoefficients(
        const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
        const Vector&                                     rDetJs,
        const Geometry<Node>&                             rGeometry) const;

    [[nodiscard]] virtual std::unique_ptr<IntegrationCoefficientsCalculator> Clone() const = 0;

private:
    [[nodiscard]] virtual double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                 double DetJ,
                                                                 const Geometry<Node>& rGeometry) const
    {
        KRATOS_ERROR
            << "IntegrationCoefficientsCalculator::CalculateIntegrationCoefficient is called." << std::endl;
    }

    friend class Serializer;
    virtual void save(Serializer& rSerializer) const = 0;
    virtual void load(Serializer& rSerializer)       = 0;
};

// Prototype for a simpler design
class IntegrationCoefficientModifier
{
public:
    virtual ~IntegrationCoefficientModifier()                                           = default;
    [[nodiscard]] virtual std::unique_ptr<IntegrationCoefficientModifier> Clone() const = 0;

    virtual double operator()(double                           IntegrationCoefficient,
                              const Geo::IntegrationPointType& rIntegrationPoint,
                              const Element&                   rElement) const = 0;
};

class KRATOS_API(GEO_MECHANICS_APPLICATION) CalculateIntegrationCoefficients0
{
public:
    explicit CalculateIntegrationCoefficients0(std::unique_ptr<IntegrationCoefficientModifier> = nullptr);

    [[nodiscard]] std::unique_ptr<IntegrationCoefficientModifier> CloneModifier() const;

    template <typename OutputContainer = std::vector<double>>
    OutputContainer Run(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                        const Vector&                          rDetJs,
                        const Element*                         pElement = nullptr) const
    {
        auto result = OutputContainer(rIntegrationPoints.size());
        auto calculate_integration_coefficient = [](const auto& rIntegrationPoint, auto DetJ) {
            return rIntegrationPoint.Weight() * DetJ;
        };
        std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), rDetJs.begin(),
                       result.begin(), calculate_integration_coefficient);

        if (mCoefficientModifier && pElement) {
            auto apply_modifier = [&element = *pElement, &modifier = *mCoefficientModifier](
                                      auto IntegrationCoefficient, const auto& rIntegrationPoint) {
                return modifier(IntegrationCoefficient, rIntegrationPoint, element);
            };
            std::transform(result.cbegin(), result.cend(), rIntegrationPoints.begin(),
                           result.begin(), apply_modifier);
        }

        return result;
    }

private:
    std::unique_ptr<IntegrationCoefficientModifier> mCoefficientModifier;
};

} // namespace Kratos