//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti

// System includes
#include <algorithm>
#include <cmath>
#include <limits>

// Project includes
#include "custom_utilities/mapping_triangulation_utilities.h"

namespace Kratos
{
namespace MappingTriangulationUtilities
{
namespace
{

using PolygonType = std::vector<CoordinatesArrayType>;

double SignedPolygonArea(const PolygonType& rPolygon)
{
    double twice_area = 0.0;
    for (std::size_t i = 0; i < rPolygon.size(); ++i) {
        const auto& r_point = rPolygon[i];
        const auto& r_next_point = rPolygon[(i + 1) % rPolygon.size()];
        twice_area += r_point[0] * r_next_point[1] -
                      r_next_point[0] * r_point[1];
    }
    return 0.5 * twice_area;
}

double TriangleArea(const TriangleType& rTriangle)
{
    return 0.5 * std::abs(
        (rTriangle[1][0] - rTriangle[0][0]) *
            (rTriangle[2][1] - rTriangle[0][1]) -
        (rTriangle[2][0] - rTriangle[0][0]) *
            (rTriangle[1][1] - rTriangle[0][1]));
}

void RemoveDegenerateVertices(PolygonType& rPolygon, const double Tolerance)
{
    PolygonType unique_points;
    unique_points.reserve(rPolygon.size());
    for (const auto& r_point : rPolygon) {
        if (unique_points.empty() ||
            norm_2(r_point - unique_points.back()) > Tolerance) {
            unique_points.push_back(r_point);
        }
    }
    if (unique_points.size() > 1 &&
        norm_2(unique_points.front() - unique_points.back()) <= Tolerance) {
        unique_points.pop_back();
    }

    bool removed_vertex = true;
    while (removed_vertex && unique_points.size() >= 3) {
        removed_vertex = false;
        for (std::size_t i = 0; i < unique_points.size(); ++i) {
            const auto& r_previous =
                unique_points[(i + unique_points.size() - 1) % unique_points.size()];
            const auto& r_current = unique_points[i];
            const auto& r_next = unique_points[(i + 1) % unique_points.size()];

            const double cross_product =
                (r_current[0] - r_previous[0]) * (r_next[1] - r_current[1]) -
                (r_current[1] - r_previous[1]) * (r_next[0] - r_current[0]);
            const double edge_scale = std::max(
                norm_2(r_current - r_previous) * norm_2(r_next - r_current),
                Tolerance);

            if (std::abs(cross_product) <= Tolerance * edge_scale) {
                unique_points.erase(unique_points.begin() + i);
                removed_vertex = true;
                break;
            }
        }
    }
    rPolygon = std::move(unique_points);
}

PolygonType ClipAgainstBoundary(
    const PolygonType& rPolygon,
    const std::size_t Direction,
    const double Boundary,
    const bool KeepGreater,
    const double Tolerance)
{
    PolygonType clipped_polygon;
    if (rPolygon.empty()) {
        return clipped_polygon;
    }
    clipped_polygon.reserve(rPolygon.size() + 1);

    const auto is_inside = [=](const CoordinatesArrayType& rPoint) {
        return KeepGreater
            ? rPoint[Direction] >= Boundary - Tolerance
            : rPoint[Direction] <= Boundary + Tolerance;
    };

    const auto intersection = [=](
        const CoordinatesArrayType& rFirst,
        const CoordinatesArrayType& rSecond)
    {
        const double denominator = rSecond[Direction] - rFirst[Direction];
        KRATOS_ERROR_IF(std::abs(denominator) <=
            std::numeric_limits<double>::epsilon())
            << "Cannot intersect an edge parallel to a clipping boundary.\n";

        const double position = (Boundary - rFirst[Direction]) / denominator;
        CoordinatesArrayType point = rFirst + position * (rSecond - rFirst);
        point[Direction] = Boundary;
        point[2] = 0.0;
        return point;
    };

    CoordinatesArrayType previous_point = rPolygon.back();
    bool previous_is_inside = is_inside(previous_point);
    for (const auto& r_current_point : rPolygon) {
        const bool current_is_inside = is_inside(r_current_point);
        if (current_is_inside) {
            if (!previous_is_inside) {
                clipped_polygon.push_back(intersection(previous_point, r_current_point));
            }
            clipped_polygon.push_back(r_current_point);
        } else if (previous_is_inside) {
            clipped_polygon.push_back(intersection(previous_point, r_current_point));
        }
        previous_point = r_current_point;
        previous_is_inside = current_is_inside;
    }
    return clipped_polygon;
}

PolygonType ClipTriangleWithRectangle(
    const TriangleType& rTriangle,
    const double UMin,
    const double UMax,
    const double VMin,
    const double VMax,
    const double Tolerance)
{
    PolygonType polygon{rTriangle.begin(), rTriangle.end()};
    polygon = ClipAgainstBoundary(polygon, 0, UMin, true, Tolerance);
    polygon = ClipAgainstBoundary(polygon, 0, UMax, false, Tolerance);
    polygon = ClipAgainstBoundary(polygon, 1, VMin, true, Tolerance);
    polygon = ClipAgainstBoundary(polygon, 1, VMax, false, Tolerance);
    RemoveDegenerateVertices(polygon, Tolerance);
    return polygon;
}

void TriangulatePolygon(
    PolygonType& rPolygon,
    std::vector<TriangleType>& rTriangles,
    const double AreaTolerance)
{
    if (rPolygon.size() < 3) {
        return;
    }
    if (SignedPolygonArea(rPolygon) < 0.0) {
        std::reverse(rPolygon.begin(), rPolygon.end());
    }

    for (std::size_t i = 1; i + 1 < rPolygon.size(); ++i) {
        TriangleType triangle{{rPolygon[0], rPolygon[i], rPolygon[i + 1]}};
        if (TriangleArea(triangle) > AreaTolerance) {
            rTriangles.push_back(std::move(triangle));
        }
    }
}

} // unnamed namespace

void Triangulation(
    const TriangleType& rOriginalTriangleCoordinates,
    GeometryPointerType pMasterGeometry,
    std::vector<TriangleType>& rNewTriangles)
{
    KRATOS_ERROR_IF_NOT(pMasterGeometry)
        << "A valid master geometry is required for triangle subdivision.\n";

    auto p_background_geometry = pMasterGeometry->pGetGeometryPart(
        GeometryType::BACKGROUND_GEOMETRY_INDEX);
    KRATOS_ERROR_IF_NOT(p_background_geometry)
        << "The master geometry has no background geometry.\n";

    std::vector<double> knot_spans_u;
    std::vector<double> knot_spans_v;
    p_background_geometry->SpansLocalSpace(knot_spans_u, 0);
    p_background_geometry->SpansLocalSpace(knot_spans_v, 1);
    KRATOS_ERROR_IF(knot_spans_u.size() < 2 || knot_spans_v.size() < 2)
        << "At least two span boundaries are required in each direction.\n";

    const auto [triangle_u_min_it, triangle_u_max_it] = std::minmax_element(
        rOriginalTriangleCoordinates.begin(), rOriginalTriangleCoordinates.end(),
        [](const CoordinatesArrayType& rA, const CoordinatesArrayType& rB) {
            return rA[0] < rB[0];
        });
    const auto [triangle_v_min_it, triangle_v_max_it] = std::minmax_element(
        rOriginalTriangleCoordinates.begin(), rOriginalTriangleCoordinates.end(),
        [](const CoordinatesArrayType& rA, const CoordinatesArrayType& rB) {
            return rA[1] < rB[1];
        });

    const double parameter_scale = std::max({
        1.0,
        std::abs(knot_spans_u.back() - knot_spans_u.front()),
        std::abs(knot_spans_v.back() - knot_spans_v.front())});
    const double coordinate_tolerance = 1e-12 * parameter_scale;
    const double original_area = TriangleArea(rOriginalTriangleCoordinates);
    const double area_tolerance = 1e-14 * parameter_scale * parameter_scale;
    KRATOS_ERROR_IF(original_area <= area_tolerance)
        << "Cannot subdivide a degenerate triangle with area " << original_area << ".\n";

    rNewTriangles.clear();
    for (std::size_t i = 0; i + 1 < knot_spans_u.size(); ++i) {
        const double u_min = knot_spans_u[i];
        const double u_max = knot_spans_u[i + 1];
        if (u_max < (*triangle_u_min_it)[0] - coordinate_tolerance ||
            u_min > (*triangle_u_max_it)[0] + coordinate_tolerance) {
            continue;
        }

        for (std::size_t j = 0; j + 1 < knot_spans_v.size(); ++j) {
            const double v_min = knot_spans_v[j];
            const double v_max = knot_spans_v[j + 1];
            if (v_max < (*triangle_v_min_it)[1] - coordinate_tolerance ||
                v_min > (*triangle_v_max_it)[1] + coordinate_tolerance) {
                continue;
            }

            auto clipped_polygon = ClipTriangleWithRectangle(
                rOriginalTriangleCoordinates,
                u_min, u_max, v_min, v_max, coordinate_tolerance);
            TriangulatePolygon(clipped_polygon, rNewTriangles, area_tolerance);
        }
    }

    double triangulated_area = 0.0;
    for (const auto& r_triangle : rNewTriangles) {
        triangulated_area += TriangleArea(r_triangle);
    }
    constexpr double relative_area_tolerance = 1e-6;
    KRATOS_ERROR_IF(
        std::abs(triangulated_area - original_area) / original_area >
        relative_area_tolerance)
        << "Knot-span subdivision does not conserve triangle area. Original area: "
        << original_area << ", triangulated area: " << triangulated_area << ".\n";
}

} // namespace MappingTriangulationUtilities
} // namespace Kratos
