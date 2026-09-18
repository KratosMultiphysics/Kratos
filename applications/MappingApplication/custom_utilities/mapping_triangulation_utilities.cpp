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
#include "custom_utilities/brep_clipper_utilities.h"

namespace Kratos
{
namespace MappingTriangulationUtilities
{
namespace
{

/*
 * Computes the signed area of a polygon in the local parameter plane.
 * rPolygon: Ordered polygon vertices in parametric coordinates.
 * Returns a positive area for counter-clockwise ordering and a negative area for clockwise ordering.
 */
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

/*
 * Computes the absolute area of a triangle in the local parameter plane.
 * rTriangle: Three triangle vertices in parametric coordinates.
 * Returns the non-negative triangle area.
 */
double TriangleArea(const TriangleType& rTriangle)
{
    return 0.5 * std::abs(
        (rTriangle[1][0] - rTriangle[0][0]) *
            (rTriangle[2][1] - rTriangle[0][1]) -
        (rTriangle[2][0] - rTriangle[0][0]) *
            (rTriangle[1][1] - rTriangle[0][1]));
}

/*
 * Removes consecutive duplicate and collinear vertices from a polygon.
 * rPolygon: Polygon modified in place.
 * Tolerance: Distance and collinearity tolerance.
 */
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

/*
 * Clips a polygon against one axis-aligned half-plane.
 * rPolygon: Polygon to clip.
 * Direction: Coordinate index: zero for u and one for v.
 * Boundary: Coordinate value of the clipping line.
 * KeepGreater: Whether values greater than the boundary are retained.
 * Tolerance: Tolerance used by the inside test.
 * Returns the polygon retained by the specified half-plane.
 */
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

/*
 * Clips a triangle to an axis-aligned knot-span rectangle.
 * rTriangle: Triangle to clip in parametric coordinates.
 * UMin: Lower u boundary.
 * UMax: Upper u boundary.
 * VMin: Lower v boundary.
 * VMax: Upper v boundary.
 * Tolerance: Clipping and degenerate-vertex tolerance.
 * Returns the polygon representing the part of the triangle inside the rectangle.
 */
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

/*
 * Triangulates a simple ordered polygon using ear clipping while discarding
 * degenerate triangles.
 * rPolygon: Polygon reordered in place when its orientation is clockwise.
 * rTriangles: Generated triangles appended to this container.
 * AreaTolerance: Minimum area required to retain a generated triangle.
 */
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

    PolygonType remaining_vertices = rPolygon;
    while (remaining_vertices.size() > 3) {
        bool ear_was_found = false;
        for (std::size_t i = 0; i < remaining_vertices.size(); ++i) {
            const std::size_t previous_index =
                (i + remaining_vertices.size() - 1) % remaining_vertices.size();
            const std::size_t next_index = (i + 1) % remaining_vertices.size();
            const auto& r_previous = remaining_vertices[previous_index];
            const auto& r_current = remaining_vertices[i];
            const auto& r_next = remaining_vertices[next_index];

            const double cross_product =
                (r_current[0] - r_previous[0]) *
                    (r_next[1] - r_current[1]) -
                (r_current[1] - r_previous[1]) *
                    (r_next[0] - r_current[0]);
            if (cross_product <= 2.0 * AreaTolerance) {
                continue;
            }

            const auto is_inside_ear = [&r_previous, &r_current, &r_next,
                AreaTolerance](const CoordinatesArrayType& rPoint) {
                const auto edge_cross = [&rPoint](
                    const CoordinatesArrayType& rFirst,
                    const CoordinatesArrayType& rSecond) {
                    return (rSecond[0] - rFirst[0]) *
                               (rPoint[1] - rFirst[1]) -
                           (rSecond[1] - rFirst[1]) *
                               (rPoint[0] - rFirst[0]);
                };
                return edge_cross(r_previous, r_current) >= -AreaTolerance &&
                       edge_cross(r_current, r_next) >= -AreaTolerance &&
                       edge_cross(r_next, r_previous) >= -AreaTolerance;
            };

            bool contains_other_vertex = false;
            for (std::size_t j = 0; j < remaining_vertices.size(); ++j) {
                if (j != previous_index && j != i && j != next_index &&
                    is_inside_ear(remaining_vertices[j])) {
                    contains_other_vertex = true;
                    break;
                }
            }
            if (contains_other_vertex) {
                continue;
            }

            rTriangles.push_back({{r_previous, r_current, r_next}});
            remaining_vertices.erase(remaining_vertices.begin() + i);
            ear_was_found = true;
            break;
        }

        KRATOS_ERROR_IF_NOT(ear_was_found)
            << "Failed to find an ear while triangulating a simple polygon.\n";
    }

    TriangleType last_triangle{{
        remaining_vertices[0], remaining_vertices[1], remaining_vertices[2]}};
    if (TriangleArea(last_triangle) > AreaTolerance) {
        rTriangles.push_back(std::move(last_triangle));
    }
}

/*
 * Orders parametric polygon vertices counter-clockwise about their centroid.
 * rVertices: Polygon vertices reordered in place.
 */
void SortVerticesCounterClockwise(PolygonType& rVertices)
{
    CoordinatesArrayType centroid = ZeroVector(3);
    for (const auto& r_vertex : rVertices) {
        centroid += r_vertex;
    }
    centroid /= static_cast<double>(rVertices.size());

    std::sort(
        rVertices.begin(), rVertices.end(),
        [&centroid](const CoordinatesArrayType& rA, const CoordinatesArrayType& rB) {
            return std::atan2(rA[1] - centroid[1], rA[0] - centroid[0]) <
                   std::atan2(rB[1] - centroid[1], rB[0] - centroid[0]);
        });
}

} // unnamed namespace

/*
 * Sorts a convex parametric polygon and triangulates it using a fan.
 * rPolygon: Polygon sorted counter-clockwise in place.
 * rTriangles: Generated three-vertex polygons appended to this container.
 */
void TriangulatePolygonFan(
    PolygonType& rPolygon,
    std::vector<PolygonType>& rTriangles)
{
    if (rPolygon.size() < 3) {
        return;
    }

    SortVerticesCounterClockwise(rPolygon);
    double coordinate_scale = 1.0;
    for (const auto& r_point : rPolygon) {
        coordinate_scale = std::max({
            coordinate_scale,
            std::abs(r_point[0]),
            std::abs(r_point[1])});
    }
    constexpr double relative_area_tolerance = 1e-14;
    if (std::abs(SignedPolygonArea(rPolygon)) <=
        relative_area_tolerance * coordinate_scale * coordinate_scale) {
        return;
    }

    for (std::size_t i = 1; i + 1 < rPolygon.size(); ++i) {
        rTriangles.push_back({rPolygon[0], rPolygon[i], rPolygon[i + 1]});
    }
}

/*
 * Decomposes and triangulates a clipped region that can contain holes.
 * rTrimmedRegion: Retained region represented by integer Clipper paths.
 * Factor: Conversion factor between parametric and integer coordinates.
 * rTriangles: Generated parametric triangles appended to this container.
 */
void TriangulateTrimmedRegion(
    const Clipper2Lib::Paths64& rTrimmedRegion,
    const double Factor,
    std::vector<Matrix>& rTriangles)
{
    if (rTrimmedRegion.empty()) {
        return;
    }

    std::vector<std::int64_t> strip_boundaries;
    std::int64_t minimum_y = std::numeric_limits<std::int64_t>::max();
    std::int64_t maximum_y = std::numeric_limits<std::int64_t>::lowest();
    for (const auto& r_path : rTrimmedRegion) {
        for (const auto& r_point : r_path) {
            strip_boundaries.push_back(r_point.x);
            minimum_y = std::min(minimum_y, r_point.y);
            maximum_y = std::max(maximum_y, r_point.y);
        }
    }

    std::sort(strip_boundaries.begin(), strip_boundaries.end());
    strip_boundaries.erase(
        std::unique(strip_boundaries.begin(), strip_boundaries.end()),
        strip_boundaries.end());

    if (strip_boundaries.size() < 2 || minimum_y >= maximum_y) {
        return;
    }

    for (std::size_t i = 0; i + 1 < strip_boundaries.size(); ++i) {
        const std::int64_t left = strip_boundaries[i];
        const std::int64_t right = strip_boundaries[i + 1];
        if (left == right) {
            continue;
        }

        const Clipper2Lib::Rect64 strip_rectangle(
            left, minimum_y, right, maximum_y);
        Clipper2Lib::Clipper64 strip_clipper;
        strip_clipper.AddSubject(rTrimmedRegion);
        strip_clipper.AddClip(
            Clipper2Lib::Paths64{strip_rectangle.AsPath()});

        Clipper2Lib::Paths64 strip_regions;
        strip_clipper.Execute(
            Clipper2Lib::ClipType::Intersection,
            Clipper2Lib::FillRule::NonZero,
            strip_regions);

        for (const auto& r_strip_region : strip_regions) {
            const double region_area = std::abs(
                Clipper2Lib::Area(r_strip_region));
            if (r_strip_region.size() >= 3 && region_area > 0.0) {
                PolygonType strip_polygon;
                strip_polygon.reserve(r_strip_region.size());
                for (const auto& r_point : r_strip_region) {
                    const auto point =
                        BrepTrimmingUtilities<false>::IntPointToDoublePoint(
                            r_point, Factor);
                    strip_polygon.push_back(
                        CoordinatesArrayType{point[0], point[1], 0.0});
                }

                const double coordinate_tolerance =
                    std::numeric_limits<double>::epsilon() *
                    std::max(1.0, Factor * std::sqrt(region_area));
                RemoveDegenerateVertices(
                    strip_polygon, coordinate_tolerance);

                std::vector<TriangleType> strip_triangles;
                TriangulatePolygon(
                    strip_polygon,
                    strip_triangles,
                    coordinate_tolerance * coordinate_tolerance);
                for (const auto& r_triangle : strip_triangles) {
                    Matrix triangle(3, 2);
                    for (std::size_t j = 0; j < 3; ++j) {
                        triangle(j, 0) = r_triangle[j][0];
                        triangle(j, 1) = r_triangle[j][1];
                    }
                    rTriangles.push_back(std::move(triangle));
                }
            }
        }
    }
}

/*
 * Clips candidate triangles against outer trimming loops and inner holes.
 * rCandidateTriangles: Parametric triangles to clip.
 * rBrepSurface: BRep surface providing the trimming loops.
 * Factor: Conversion factor between parametric and integer coordinates.
 * Returns triangles covering only the retained trimmed region.
 */
std::vector<PolygonType> ClipTrianglesWithTrimmingLoops(
    const std::vector<PolygonType>& rCandidateTriangles,
    const BrepSurfaceType& rBrepSurface,
    const double Factor)
{
    const auto all_loops = BrepClipperUtilities::CreateAllLoops(
        rBrepSurface, Factor);
    Clipper2Lib::Paths64 outer_paths;
    Clipper2Lib::Paths64 inner_paths;
    BrepClipperUtilities::SplitOuterAndInnerPaths(
        all_loops, outer_paths, inner_paths);
    if (outer_paths.empty()) {
        return rCandidateTriangles;
    }

    std::vector<PolygonType> clipped_triangles;
    for (std::size_t candidate_index = 0;
        candidate_index < rCandidateTriangles.size();
        ++candidate_index) {
        const auto& r_candidate_triangle =
            rCandidateTriangles[candidate_index];
        KRATOS_ERROR_IF(r_candidate_triangle.size() != 3)
            << "Expected a candidate triangle with 3 vertices, got "
            << r_candidate_triangle.size() << ".\n";

        Clipper2Lib::Paths64 triangle_paths(1);
        for (const auto& r_point : r_candidate_triangle) {
            triangle_paths[0].push_back(
                BrepTrimmingUtilities<false>::ToIntPoint(
                    r_point[0], r_point[1], Factor));
        }

        const double triangle_area = std::abs(Clipper2Lib::Area(triangle_paths));
        KRATOS_ERROR_IF(
            triangle_area <= std::numeric_limits<double>::epsilon())
            << "Cannot process a degenerate triangle with area "
            << triangle_area << ".\n";

        const auto trimmed_region =
            BrepClipperUtilities::ClipPathsWithTrimmedDomain(
                triangle_paths, outer_paths, inner_paths);
        if (trimmed_region.empty()) {
            continue;
        }

        const double trimmed_area = std::abs(Clipper2Lib::Area(trimmed_region));
        constexpr double relative_area_tolerance = 1e-6;
        if (std::abs(triangle_area - trimmed_area) / triangle_area <
            relative_area_tolerance) {
            clipped_triangles.push_back(r_candidate_triangle);
            continue;
        }

        std::vector<Matrix> triangulated_region;
        TriangulateTrimmedRegion(trimmed_region, Factor, triangulated_region);
        for (const auto& r_triangle : triangulated_region) {
            PolygonType triangle_vertices(3);
            for (std::size_t i = 0; i < 3; ++i) {
                triangle_vertices[i][0] = r_triangle(i, 0);
                triangle_vertices[i][1] = r_triangle(i, 1);
                triangle_vertices[i][2] = 0.0;
            }
            SortVerticesCounterClockwise(triangle_vertices);
            clipped_triangles.push_back(std::move(triangle_vertices));
        }
    }

    return clipped_triangles;
}

/*
 * Subdivides a parametric triangle at every intersected knot span.
 * rOriginalTriangleCoordinates: Triangle to subdivide.
 * pMasterGeometry: Geometry whose background surface provides the knot spans.
 * rNewTriangles: Subdivided triangles replacing the previous contents.
 */
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
