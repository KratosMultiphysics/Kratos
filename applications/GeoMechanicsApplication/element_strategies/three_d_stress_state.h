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

#include "stress_state_strategy.h"

namespace Kratos
{

class ThreeDStressState : public StressStateStrategy
{
public:
    void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector&, const Geometry<Node>& rGeometry) override;
    double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointsArrayType& IntegrationPoints,
                                           unsigned int PointNumber,
                                           double detJ,
                                           const Geometry<Node>& rGeometry) override;
    std::unique_ptr<StressStateStrategy> Create() override;
};

} // namespace Kratos
