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
        KRATOS_ERROR_IF_NOT((rThisPoints.size() == 4) || (rThisPoints.size() == 6))
            << "Number of nodes must be four or six\n";

        auto points = this->Points();
        points.resize(points.size() / 2);
        if (points.size() == 2) {
            mLineGeometry = std::make_unique<Line2D2<Node>>(points);
        } else {
            mLineGeometry = std::make_unique<Line2D3<Node>>(points);
        }
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

    [[nodiscard]] double ShapeFunctionValue(IndexType                   ShapeFunctionIndex,
                                            const CoordinatesArrayType& rCoordinates) const override
    {
        return mLineGeometry->ShapeFunctionValue(ShapeFunctionIndex, rCoordinates);
    }

    Vector& ShapeFunctionsValues(Vector& rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        rResult.resize(mLineGeometry->PointsNumber());
        for (auto i = IndexType{0}; i < mLineGeometry->PointsNumber(); ++i) {
            rResult[i] = this->ShapeFunctionValue(i, rCoordinates);
        }
        return rResult;
    }

    Matrix& ShapeFunctionsLocalGradients(Matrix& rResult, const CoordinatesArrayType& rPoint) const override
    {
        return mLineGeometry->ShapeFunctionsLocalGradients(rResult, rPoint);
    }

private:
    std::unique_ptr<Geometry<Node>> mLineGeometry;
};

} // namespace Kratos
