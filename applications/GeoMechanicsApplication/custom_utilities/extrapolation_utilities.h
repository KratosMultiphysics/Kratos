// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Gennady Markelov
//

#pragma once

#include "geometries/geometry.h"
#include "includes/kratos_export_api.h"

namespace Kratos
{
class Node;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ExtrapolationUtilities
{
public:
    [[nodiscard]] static Matrix CalculateExtrapolationMatrix(const Geometry<Node>& rGeometry,
                                                             const GeometryData::IntegrationMethod& rIntegrationMethod);

    [[nodiscard]] static std::vector<std::optional<Vector>> CalculateNodalStresses(
        const std::vector<std::size_t>&        node_ids,
        const Geometry<Node>&                  rGeometry,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        const std::vector<Vector>&             rIntegrationPointStresses);
};

} // namespace Kratos