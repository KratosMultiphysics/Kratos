//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "geometries/line_3d_2.h"
#include "utilities/intersection_utilities.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
GeometricalProjectionUtilities::DistanceComputed GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius(
    double& rDistance,
    const Geometry<Node>& rSegment,
    const Point& rPoint,
    const double Radius,
    const double Tolerance
    )
{
    // If radius is zero, we keep the distance as it is
    if (Radius < std::numeric_limits<double>::epsilon()) {
        return DistanceComputed::NO_RADIUS;
    } else {
        Line3D2<Node> line(rSegment.Points()); // NOTE: TO ENSURE THAT IT IS 3D IN CASE IS DECLARED AS 2D
        Point projected_point;
        const double projected_distance = GeometricalProjectionUtilities::FastProjectOnLine(line, rPoint, projected_point);
        typename Geometry<Node>::CoordinatesArrayType projected_local;
        // If projection is inside, just remove the radius
        if (line.IsInside(projected_point.Coordinates(), projected_local, Tolerance)) {
            rDistance = projected_distance - Radius;
            return DistanceComputed::RADIUS_PROJECTED;
        } else { // Othwerise we compute the distance to the closest node and compute the difference with the "radius cylinder"
            // Distances to the nodes
            const Point::Pointer point = Kratos::make_shared<Point>(rPoint.Coordinates());
            const Point::Pointer point_a = Kratos::make_shared<Point>(line[0].Coordinates());
            const Point::Pointer point_b = Kratos::make_shared<Point>(line[1].Coordinates());
            const double distance_a = rPoint.Distance(*point_a);
            const double distance_b = rPoint.Distance(*point_b);

            // Positive distance. Remove distance to parallel line in Radius
            if (projected_distance > Radius) {
                array_1d<double, 3> vector_line = line[1].Coordinates() - line[0].Coordinates();
                vector_line /= norm_2(vector_line);
                const double N_line = Radius/projected_distance;
                const Point::Pointer point_distance_r_projection = Kratos::make_shared<Point>((1.0 - N_line) * projected_point.Coordinates() + N_line * rPoint.Coordinates());
                const Point::Pointer aux_point_parallel = Kratos::make_shared<Point>(point_distance_r_projection->Coordinates() + vector_line);

                // Parallel line to the segment
                Geometry<Point>::PointsArrayType points_array_parallel_line;
                points_array_parallel_line.reserve(2);
                points_array_parallel_line.push_back(point_distance_r_projection);
                points_array_parallel_line.push_back(aux_point_parallel);
                Line3D2<Point> parallel_line(points_array_parallel_line);

                // Line from the point to the segment
                Geometry<Point>::PointsArrayType points_array_distance;
                points_array_distance.reserve(2);
                points_array_distance.push_back(point);
                if (distance_a < distance_b) {
                    points_array_distance.push_back(point_a);
                } else {
                    points_array_distance.push_back(point_b);
                }
                Line3D2<Point> distance_line(points_array_distance);

                // Compute intersection
                const auto line_intersection = Line3D2<Point>(IntersectionUtilities::ComputeShortestLineBetweenTwoLines(parallel_line, distance_line));
                const Point intersection_point(line_intersection.Center().Coordinates());

                // Compute distance
                if (distance_a < distance_b) {
                    rDistance = distance_a - point_a->Distance(intersection_point);
                } else {
                    rDistance = distance_b - point_b->Distance(intersection_point);
                }
                return DistanceComputed::RADIUS_NOT_PROJECTED_OUTSIDE;
            } else { // Negative distance
                array_1d<double, 3> projection_vector = rPoint.Coordinates() - projected_point.Coordinates();
                projection_vector /= norm_2(projection_vector);

                // Parallel projected point
                const array_1d<double, 3> aux_parallel_projection = (distance_a < distance_b) ? point_a->Coordinates() + projection_vector * Radius : point_b->Coordinates() + projection_vector * Radius;
                Point parallel_projection_point(aux_parallel_projection);

                // Compute distance
                rDistance = - point->Distance(parallel_projection_point);
                return DistanceComputed::RADIUS_NOT_PROJECTED_INSIDE;
            }
        }
    }
    return DistanceComputed::PROJECTION_ERROR;
}

} // namespace Kratos