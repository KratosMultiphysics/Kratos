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
#include "includes/define.h"
#include "includes/node.h"
#include <cstddef>

namespace Kratos
{

class NodalExtrapolator
{
public:
    using NodeType     = Node;
    using GeometryType = Geometry<NodeType>;
    using SizeType     = std::size_t;
    using IndexType    = std::size_t;

    Matrix CalculateElementExtrapolationMatrix(GeometryType& rGeometry,
                                               GeometryData::IntegrationMethod IntegrationMethod) const;
    void   CheckIfGeometryIsSupported(const GeometryType& r_this_geometry) const;
    std::unique_ptr<NodalExtrapolator::GeometryType> CreateLowerOrderGeometry(GeometryType& rGeometry) const;
};

} // namespace Kratos
