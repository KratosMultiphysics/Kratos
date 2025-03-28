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

#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class Serializer;

class KRATOS_API(GEO_MECHANICS_APPLICATION) IntegrationCoefficientsCalculator
{
public:
    virtual ~IntegrationCoefficientsCalculator() = default;

    typedef Node               NodeType;
    typedef Geometry<NodeType> GeometryType;

    [[nodiscard]] virtual Vector CalculateIntegrationCoefficients(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
                                                                  const Vector& rDetJs,
                                                                  std::size_t   LocalDimension,
                                                                  double        CrossArea) const
    {
        KRATOS_ERROR
            << "IntegrationCoefficientsCalculator::CalculateIntegrationCoefficients is called."
            << std::endl;
    }

    [[nodiscard]] virtual std::vector<double> CalculateIntegrationCoefficients(
        const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
        const Vector&                                   rDetJs,
        const Geometry<Node>&                           rGeometry) const;

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

} // namespace Kratos