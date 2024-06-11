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

    struct TLSType {
        Vector vector_J;
        Vector N;
    };

    Matrix CalculateElementExtrapolationMatrix(GeometryType& r_this_geometry,
                                               GeometryData::IntegrationMethod this_integration_method) const;
};

} // namespace Kratos
