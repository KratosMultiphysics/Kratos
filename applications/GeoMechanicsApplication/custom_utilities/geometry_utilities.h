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

#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "includes/kratos_export_api.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeometryUtilities
{
public:
    static Matrix Calculate2DRotationMatrixForLineGeometry(const Geometry<Node>& rGeometry,
                                                           const array_1d<double, 3>& rLocalCoordinate);
    static Matrix Calculate3DRotationMatrixForPlaneGeometry(const Geometry<Node>& rGeometry,
                                                            const array_1d<double, 3>& rLocalCoordinate);
    static std::vector<std::size_t> GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry);

    [[nodiscard]] static auto GetReversedNodesForGeometryFamily(const GeometryData::KratosGeometryFamily& rGeometryFamily,
                                                                const auto& rOrderedNodes)
    {
        auto result = rOrderedNodes;

        // For line geometries we want to reverse all 'corner points', while for surfaces we don't
        // change the starting node, but only reverse the order of the rest of the corner points.
        auto begin_of_corner_points = rGeometryFamily == GeometryData::KratosGeometryFamily::Kratos_Linear
                                          ? result.begin()
                                          : result.begin() + 1;
        auto end_of_corner_points = result.begin() + GetNumberOfCornerPoints(rGeometryFamily);

        std::reverse(begin_of_corner_points, end_of_corner_points);
        std::reverse(end_of_corner_points, result.end());

        return result;
    }

private:
    [[nodiscard]] static std::size_t GetNumberOfCornerPoints(const GeometryData::KratosGeometryFamily& rGeometryFamily);
};

} // namespace Kratos
