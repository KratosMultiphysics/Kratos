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
    [[nodiscard]] virtual Matrix CalculateBMatrix(const Matrix&         GradNpT,
                                                  const Vector&         Np,
                                                  const Geometry<Node>& rGeometry) const = 0;
    [[nodiscard]] virtual double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                 double detJ,
                                                                 const Geometry<Node>& rGeometry) const = 0;
    [[nodiscard]] virtual Vector CalculateGreenLagrangeStrain(const Matrix& rTotalDeformationGradient) const = 0;
    [[nodiscard]] virtual std::unique_ptr<StressStatePolicy> Clone() const = 0;

    virtual ~StressStatePolicy() = default;
};

} // namespace Kratos