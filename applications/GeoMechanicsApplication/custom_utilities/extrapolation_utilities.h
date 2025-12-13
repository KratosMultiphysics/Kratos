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
class NodalExtrapolator;
class Node;
class Element;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ExtrapolationUtilities
{
public:
    [[nodiscard]] static Matrix CalculateExtrapolationMatrix(const Element& rElement);

    [[nodiscard]] static std::vector<std::optional<Vector>> CalculateNodalVectors(
        const std::vector<std::size_t>& rNodeIds, const Element& rElement, const std::vector<Vector>& rVectorsAtIntegrationPoints);
};

} // namespace Kratos