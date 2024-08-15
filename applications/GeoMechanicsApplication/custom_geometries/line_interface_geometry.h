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
        : LineInterfaceGeometry(0, rThisPoints)
    {
    }

    LineInterfaceGeometry(const IndexType NewGeometryId, const Geometry<Node>::PointsArrayType& rThisPoints)
        : Geometry<Node>(NewGeometryId, rThisPoints)
    {
        KRATOS_ERROR_IF_NOT((rThisPoints.size() == 4) || (rThisPoints.size() == 6)) << "Number of nodes must be four or six\n";
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
        PointsArrayType points = this->Points();
        points.resize(points.size() / 2);

        switch (points.size()) {
        case 2: {
            const Line2D2<Node> line(points);
            return line.ShapeFunctionValue(ShapeFunctionIndex % 2, rCoordinates);
        }
        case 3: {
            const Line2D3<Node> line(points);
            return line.ShapeFunctionValue(ShapeFunctionIndex % 3, rCoordinates);
        }
        }
    }
};

} // namespace Kratos
