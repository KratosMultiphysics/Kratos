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

#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class StressStateStrategy
{
public:
    virtual void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np, const Geometry<Node>& rGeometry) = 0;
    virtual ~StressStateStrategy() = default;
};

} // namespace Kratos
