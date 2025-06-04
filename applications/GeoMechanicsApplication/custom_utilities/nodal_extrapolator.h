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

class NodalExtrapolator
{
public:
    using GeometryType = Geometry<Node>;

    [[nodiscard]] virtual Matrix CalculateElementExtrapolationMatrix(
        const GeometryType& rGeometry, const GeometryData::IntegrationMethod& rIntegrationMethod) const = 0;
    virtual ~NodalExtrapolator() = default;
};

} // namespace Kratos
