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

    LineInterfaceGeometry() = default;

    explicit LineInterfaceGeometry(const Geometry<Node>::PointsArrayType& rThisPoints)
        : Geometry<Node>(rThisPoints)
    {
    }

    LineInterfaceGeometry(const IndexType NewGeometryId, const Geometry<Node>::PointsArrayType& rThisPoints)
        : Geometry<Node>(NewGeometryId, rThisPoints)
    {
    }

    [[nodiscard]] Geometry<Node>::Pointer Create(const Geometry<Node>::PointsArrayType& rThisPoints) const override
    {
        const auto id = IndexType{0};
        return Create(id, rThisPoints);
    }

    [[nodiscard]] Geometry<Node>::Pointer Create(const IndexType        NewGeometryId,
                                                 const PointsArrayType& rThisPoints) const override
    {
        return std::make_shared<LineInterfaceGeometry>(NewGeometryId, rThisPoints);
    }

    double ShapeFunctionValue(IndexType ShapeFunctionIndex, const CoordinatesArrayType& rCoordinates) const override
    {
        Line2D2<Node> line(this->Points());
        return line.ShapeFunctionValue(ShapeFunctionIndex, rCoordinates);
    }
};

} // namespace Kratos
