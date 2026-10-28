//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#pragma once

#include "includes/model_part.h"

namespace Kratos {
/** @brief Physical gap quadrature carried by the existing parent FEM nodes.
 * Store the same rule for material initialization (GI_GAUSS_1) and fluid
 * assembly (GI_GAUSS_2). No interpolation nodes or unknowns are introduced.
 */
class KRATOS_API(KRATOS_CORE) GapQuadratureUtility
{
public:
    static Geometry<Node>::Pointer Create(
        Geometry<Node>& rParent,
        const Geometry<Node>::IntegrationPointsArrayType& rPhysicalPoints);
};
} // namespace Kratos
