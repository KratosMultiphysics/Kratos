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
//                   Anne van de Graaf
//

#pragma once

#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{

class LineInterfaceGeometry : public Geometry<Node>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineInterfaceGeometry);

    [[nodiscard]] Geometry<Node>::Pointer Create(const Geometry<Node>::PointsArrayType& rThisPoints) const override
    {
        return make_shared<LineInterfaceGeometry>();
    }
};

} // namespace Kratos
