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
//                   Marjan Fathian
//
#pragma once

#include "geo_mechanics_application_constants.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class StressStatePolicy
{
public:
    virtual ~StressStatePolicy() = default;

    [[nodiscard]] virtual Matrix CalculateBMatrix(const Matrix&         rGradNpT,
                                                  const Vector&         rNp,
                                                  const Geometry<Node>& rGeometry) const = 0;
    [[nodiscard]] virtual double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                 double DetJ,
                                                                 const Geometry<Node>& rGeometry) const = 0;
    [[nodiscard]] virtual Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const = 0;
    [[nodiscard]] virtual std::unique_ptr<StressStatePolicy> Clone() const = 0;
};

} // namespace Kratos