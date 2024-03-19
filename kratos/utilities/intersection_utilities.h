//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/point.h"
#include "containers/pointer_vector.h"
#include "utilities/math_utils.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class IntersectionUtilities
 * @ingroup KratosCore
 * @brief Utilities to compute intersections between different geometries
 * @details This class provides static methods to check if there is
 * intersections between different entities, and if there is, give back
 * the intersection points.
 * @author Ruben Zorrilla
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) IntersectionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IntersectionUtilities
    KRATOS_CLASS_POINTER_DEFINITION( IntersectionUtilities );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    IntersectionUtilities(){}

    /// Destructor
    virtual ~IntersectionUtilities(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * Find the 3D intersection of a line (bounded) with a triangle (bounded)
     * @param rTriangleGeometry Is the triangle to intersect
     * @param rLinePoint1 Coordinates of the first point of the intersecting line
     * @param rLinePoint2 Coordinates of the second point of the intersecting line
     * @return rIntersectionPoint The intersection point coordinates
     * @return The intersection type index:
     * -1 (the triangle is degenerate)
     * 0 (disjoint - no intersection)
     * 1 (intersect in a unique point)
     * 2 (are in the same plane)
     */
    template <class TGeometryType>
    static int ComputeTriangleLineIntersection(
        const TGeometryType& rTriangleGeometry,
        const array_1d<double,3>& rLinePoint1,
        const array_1d<double,3>& rLinePoint2,
        array_1d<double,3>& rIntersectionPoint,
        const double epsilon = 1e-12) {

        // This is the adaption of the implementation provided in:
        // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
        // Based on Tomas Möller & Ben Trumbore (1997) Fast, Minimum Storage Ray-Triangle Intersection, Journal of Graphics Tools, 2:1, 21-28, DOI: 10.1080/10867651.1997.10487468 

        // Get triangle edge vectors and plane normal
        const array_1d<double,3> u = rTriangleGeometry[1] - rTriangleGeometry[0];
        const array_1d<double,3> v = rTriangleGeometry[2] - rTriangleGeometry[0];
        array_1d<double,3> n;
        MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(n,u,v);

        // Check if the triangle is degenerate (do not deal with this case)
        if (MathUtils<double>::Norm3(n) < epsilon){
            return -1;
        }

        const array_1d<double,3> dir = rLinePoint2 - rLinePoint1; // Edge direction vector
        const array_1d<double,3> w_0 = rLinePoint1 - rTriangleGeometry[0];
        const double a = -inner_prod(n,w_0);
        const double b = inner_prod(n,dir);

        // Check if the ray is parallel to the triangle plane
        if (std::abs(b) < epsilon){
            if (a == 0.0){
                return 2;    // Edge lies in the triangle plane
            } else {
                return 0;    // Edge does not lie in the triangle plane
            }
        }

        // If the edge is not parallel, compute the intersection point
        const double r = a / b;
        if (r < 0.0){
            return 0;    // Edge goes away from triangle
        } else if (r > 1.0) {
            return 0;    // Edge goes away from triangle
        }

        rIntersectionPoint = rLinePoint1 + r*dir;

        // Check if the intersection point is inside the triangle
        if (PointInTriangle(rTriangleGeometry[0], rTriangleGeometry[1], rTriangleGeometry[2], rIntersectionPoint)) {
            return 1;
        }
        return 0;
    }

    /**
     * Find the 2D intersection of a line (bounded) with a triangle (bounded)
     * @param rTriangle Is the triangle to intersect
     * @param rPoint0 Coordinates of the first point of the intersecting line
     * @param rPoint1 Coordinates of the second point of the intersecting line
     * @return If there is intersection
     */
    template <class TGeometryType>
    static bool TriangleLineIntersection2D(
        const TGeometryType& rTriangle,
        const array_1d<double,3>& rPoint0,
        const array_1d<double,3>& rPoint1)
    {
        return TriangleLineIntersection2D(rTriangle[0], rTriangle[1], rTriangle[2], rPoint0, rPoint1);
    }

    /**
     * Find the 2D intersection of a line (bounded) with a triangle (bounded)
     * @param rVert1 The first vertex of the triangle to intersect
     * @param rVert2 The second vertex of the triangle to intersect
     * @param rVert3 The third vertex of the triangle to intersect
     * @param rPoint1 Coordinates of the first point of the intersecting line
     * @param rPoint2 Coordinates of the second point of the intersecting line
     * @return If there is intersection
     */
    static bool TriangleLineIntersection2D(
        const array_1d<double,3>& rVert0,
        const array_1d<double,3>& rVert1,
        const array_1d<double,3>& rVert2,
        const array_1d<double,3>& rPoint0,
        const array_1d<double,3>& rPoint1)
    {
        // Check the intersection of each edge against the intersecting object
        array_1d<double,3> int_point;
        if (ComputeLineLineIntersection(rVert0, rVert1, rPoint0, rPoint1, int_point)) return true;
        if (ComputeLineLineIntersection(rVert1, rVert2, rPoint0, rPoint1, int_point)) return true;
        if (ComputeLineLineIntersection(rVert2, rVert0, rPoint0, rPoint1, int_point)) return true;

        // Let check second geometry is inside the first one.
        if (PointInTriangle(rVert0, rVert1, rVert2, rPoint0)) return true;

        return false;
    }

    /**
     * Check if a point is inside a 2D triangle
     * @brief This uses the Cramer's rule for solving a linear system and obtain the barycentric coordinates
     * @details 
     * @param rVert0 The first vertex of the triangle to intersect
     * @param rVert1 The second vertex of the triangle to intersect
     * @param rVert2 The third vertex of the triangle to intersect
     * @param rPoint Coordinates of the point
     * @return The intersection type index:
     * -1 (the triangle is degenerate)
     *  0 (disjoint - no intersection)
     *  1 (intersect)
     */
    static bool PointInTriangle(
        const array_1d<double,3>& rVert0,
        const array_1d<double,3>& rVert1,
        const array_1d<double,3>& rVert2,
        const array_1d<double,3>& rPoint,
        const double Tolerance = std::numeric_limits<double>::epsilon())
    {
        const array_1d<double,3> u = rVert1 - rVert0;
        const array_1d<double,3> v = rVert2 - rVert0;
        const array_1d<double,3> w = rPoint - rVert0;
        
        const double uu = inner_prod(u, u);
        const double uv = inner_prod(u, v);
        const double vv = inner_prod(v, v);
        const double wu = inner_prod(w, u);
        const double wv = inner_prod(w, v);
        const double denom = uv * uv - uu * vv;

        const double xi  = (uv * wv - vv * wu) / denom;
        const double eta = (uv * wu - uu * wv) / denom;

        if (xi < -Tolerance) return false;
        if (eta < -Tolerance) return false;
        if (xi + eta > 1.0 + Tolerance) return false;
        return true;
    }

    /**
     * @brief Find the 3D intersection of a line (bounded) with a triangle (bounded) in the same plane
     * @tparam TGeometryType The geometry type
     * @tparam TCoordinatesType The type of coordinates
     * @tparam TConsiderInsidePoints If considering inside points or just the intersections in the faces
     * @param [in] rTriangleGeometry Is the tetrahedra to intersect
     * @param [in] rLinePoint1 Coordinates of the first point of the intersecting line
     * @param [in] rLinePoint2 Coordinates of the second point of the intersecting line
     * @param [out] rIntersectionPoint1 The first intersection point coordinates
     * @param [out] rIntersectionPoint2 The second intersection point coordinates
     * @param [out] rSolution The intersection type index:
     *         NO_INTERSECTION (disjoint - no intersection)
     *         TWO_POINTS_INTERSECTION (intersect in two points)
     *         ONE_POINT_INTERSECTION (intersect in one point)
     * @param Epsilon The tolerance
     */
    template <class TGeometryType, class TCoordinatesType>
    static int ComputeTriangleLineIntersectionInTheSamePlane( //static IntersectionUtilitiesTetrahedraLineIntersectionStatus ComputeTriangleLineIntersectionInTheSamePlane(
        const TGeometryType& rTriangleGeometry,
        const TCoordinatesType& rLinePoint1,
        const TCoordinatesType& rLinePoint2,
        TCoordinatesType& rIntersectionPoint1,
        TCoordinatesType& rIntersectionPoint2,
        int& rSolution,//IntersectionUtilitiesTetrahedraLineIntersectionStatus& rSolution,
        const double Epsilon = 1e-12
        ) 
    {
        // Lambda function to check the inside of a line corrected
        auto is_inside_projected = [&Epsilon] (auto& rGeometry, const TCoordinatesType& rPoint) -> bool {
            // We compute the distance, if it is not in the plane we project
            const Point point_to_project(rPoint);
            Point point_projected;
            const double distance = GeometricalProjectionUtilities::FastProjectOnLine(rGeometry, point_to_project, point_projected);

            // We check if we are on the plane
            if (std::abs(distance) > Epsilon * rGeometry.Length()) {
                return false;
            }
            array_1d<double, 3> local_coordinates;
            return rGeometry.IsInside(point_projected, local_coordinates);
        };

        // Compute intersection with the edges
        for (auto& r_edge : rTriangleGeometry.GenerateEdges()) {
            const auto& r_edge_point_1 = r_edge[0].Coordinates();
            const auto& r_edge_point_2 = r_edge[1].Coordinates();
            array_1d<double, 3> intersection_point_1, intersection_point_2;
            const auto check_1 = ComputeLineLineIntersection(rLinePoint1, rLinePoint2, r_edge_point_1, r_edge_point_2, intersection_point_1, Epsilon);
            const auto check_2 = ComputeLineLineIntersection(r_edge_point_1, r_edge_point_2, rLinePoint1, rLinePoint2, intersection_point_2, Epsilon);
            if (check_1 == 0 && check_2 == 0) continue; // No intersection
            array_1d<double, 3> intersection_point = check_1 != 0 ? intersection_point_1 : intersection_point_2;
            if (check_1 == 2 || check_2 == 2) { // Aligned
                // Check actually are aligned
                array_1d<double, 3> vector_line = r_edge_point_2 - r_edge_point_1;
                vector_line /= norm_2(vector_line);
                array_1d<double, 3> diff_coor_1 = rLinePoint1 - r_edge_point_1;
                const double diff_coor_1_norm = norm_2(diff_coor_1);
                if (diff_coor_1_norm > std::numeric_limits<double>::epsilon()) {
                    diff_coor_1 /= diff_coor_1_norm;
                } else {
                    diff_coor_1 = rLinePoint1 - r_edge_point_2;
                    diff_coor_1 /= norm_2(diff_coor_1);
                }
                array_1d<double, 3> diff_coor_2 = rLinePoint2 - r_edge_point_1;
                const double diff_coor_2_norm = norm_2(diff_coor_2);
                if (diff_coor_2_norm > std::numeric_limits<double>::epsilon()) {
                    diff_coor_2 /= diff_coor_2_norm;
                } else {
                    diff_coor_2 = rLinePoint2 - r_edge_point_2;
                    diff_coor_2 /= norm_2(diff_coor_2);
                }
                const double diff1m = norm_2(diff_coor_1 - vector_line);
                const double diff1p = norm_2(diff_coor_1 + vector_line);
                const double diff2m = norm_2(diff_coor_2 - vector_line);
                const double diff2p = norm_2(diff_coor_2 + vector_line);

                // Now we compute the intersection
                if ((diff1m < Epsilon || diff1p < Epsilon) && (diff2m < Epsilon || diff2p < Epsilon)) {
                    // First point
                    if (is_inside_projected(r_edge, rLinePoint1)) { // Is inside the line
                        if (rSolution == 0) {//if (rSolution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION) {
                            noalias(rIntersectionPoint1) = rLinePoint1;
                            rSolution = 2;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION;
                        } else {
                            if (norm_2(rIntersectionPoint1 - rLinePoint1) > Epsilon) { // Must be different from the first one
                                noalias(rIntersectionPoint2) = rLinePoint1;
                                rSolution = 1;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION;
                                break;
                            }
                        }
                    } else { // Is in the border of the line
                        if (rSolution == 0) {//if (rSolution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION) {
                            noalias(rIntersectionPoint1) = norm_2(r_edge_point_1 - rLinePoint1) <  norm_2(r_edge_point_2 - rLinePoint1) ? r_edge_point_1 : r_edge_point_2;
                            rSolution = 2;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION;
                        } else {
                            noalias(intersection_point) = norm_2(r_edge_point_1 - rLinePoint1) <  norm_2(r_edge_point_2 - rLinePoint1) ? r_edge_point_1 : r_edge_point_2;
                            if (norm_2(rIntersectionPoint1 - intersection_point) > Epsilon) { // Must be different from the first one
                                noalias(rIntersectionPoint2) = intersection_point;
                                rSolution = 1;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION;
                                break;
                            }
                        }
                    }
                    // Second point
                    if (rSolution == 2) {//if (rSolution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION) {
                        if (is_inside_projected(r_edge, rLinePoint2)) { // Is inside the line
                            if (norm_2(rIntersectionPoint1 - rLinePoint2) > Epsilon) { // Must be different from the first one
                                noalias(rIntersectionPoint2) = rLinePoint2;
                                rSolution = 1;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION;
                                break;
                            }
                        } else { // Is in the border of the line
                            noalias(intersection_point) = norm_2(r_edge_point_1 - rLinePoint2) <  norm_2(r_edge_point_2 - rLinePoint2) ? r_edge_point_1 : r_edge_point_2;
                            if (norm_2(rIntersectionPoint1 - intersection_point) > Epsilon) { // Must be different from the first one
                                noalias(rIntersectionPoint2) = intersection_point;
                                rSolution = 1;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION;
                                break;
                            }
                        }
                    } else { // We are done
                        break;
                    }
                }
            } else { // Direct intersection
                if (rSolution == 0) {//if (rSolution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION) {
                    noalias(rIntersectionPoint1) = intersection_point;
                    rSolution = 2;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION;
                } else {
                    if (norm_2(rIntersectionPoint1 - intersection_point) > Epsilon) { // Must be different from the first one
                        noalias(rIntersectionPoint2) = intersection_point;
                        rSolution = 1;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION;
                        break;
                    }
                }
            }
        }

        return rSolution;
    }

    /**
     * @brief Find the 3D intersection of a line (bounded) with a tetrahedra (bounded)
     * @tparam TGeometryType The geometry type
     * @tparam TCoordinatesType The type of coordinates
     * @tparam TConsiderInsidePoints If considering inside points or just the intersections in the faces
     * @param [in] rTetrahedraGeometry Is the tetrahedra to intersect
     * @param [in] rLinePoint1 Coordinates of the first point of the intersecting line
     * @param [in] rLinePoint2 Coordinates of the second point of the intersecting line
     * @param [out] rIntersectionPoint1 The first intersection point coordinates
     * @param [out] rIntersectionPoint2 The second intersection point coordinates
     * @return The intersection type index:
     *         NO_INTERSECTION (disjoint - no intersection)
     *         TWO_POINTS_INTERSECTION (intersect in two points)
     *         ONE_POINT_INTERSECTION (intersect in one point)
     *         TWO_POINTS_INTERSECTION_BOTH_INSIDE (intersect in two points inside the tetrahedra)
     *         TWO_POINTS_INTERSECTION_ONE_INSIDE (intersect in two points, one inside the tetrahedra)
     *         FIRST_CORNER (intersect in the first corner of the tetrahedra)
     *         SECOND_CORNER (intersect in the second corner of the tetrahedra)
     *         THIRD_CORNER (intersect in the thid corner of the tetrahedra)
     *         FOURTH_CORNER (intersect in the fourth corner of the tetrahedra)
     * Equivalent enum:
     *   enum class IntersectionUtilitiesTetrahedraLineIntersectionStatus
     *   {
     *       NO_INTERSECTION = 0,                     // (disjoint - no intersection)
     *       TWO_POINTS_INTERSECTION = 1,             // (intersect in two points)
     *       ONE_POINT_INTERSECTION = 2,              // (intersect in one point)
     *       TWO_POINTS_INTERSECTION_BOTH_INSIDE = 3, // (intersect in two points inside the tetrahedra)
     *       TWO_POINTS_INTERSECTION_ONE_INSIDE = 4,  // (intersect in two points, one inside the tetrahedra)
     *       FIRST_CORNER = 5,                        // (intersect in the first corner of the tetrahedra)
     *       SECOND_CORNER = 6,                       // (intersect in the second corner of the tetrahedra)
     *       THIRD_CORNER = 7,                        // (intersect in the thid corner of the tetrahedra)
     *       FOURTH_CORNER = 8                        // (intersect in the fourth corner of the tetrahedra)
     *   };
     * @param Epsilon The tolerance
     */
    template <class TGeometryType, class TCoordinatesType, bool TConsiderInsidePoints = true>
    static int ComputeTetrahedraLineIntersection(//static IntersectionUtilitiesTetrahedraLineIntersectionStatus ComputeTetrahedraLineIntersection(
        const TGeometryType& rTetrahedraGeometry,
        const TCoordinatesType& rLinePoint1,
        const TCoordinatesType& rLinePoint2,
        TCoordinatesType& rIntersectionPoint1,
        TCoordinatesType& rIntersectionPoint2,
        const double Epsilon = 1e-12
        ) 
    {
        int solution = 0;// IntersectionUtilitiesTetrahedraLineIntersectionStatus solution = IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION;
        for (auto& r_face : rTetrahedraGeometry.GenerateFaces()) {
            array_1d<double,3> intersection_point;
            const int face_solution = ComputeTriangleLineIntersection(r_face, rLinePoint1, rLinePoint2, intersection_point, Epsilon);
            if (face_solution == 1) { // The line intersects the face
                if (solution == 0) {// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION) {
                    noalias(rIntersectionPoint1) = intersection_point;
                    solution = 2;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION;
                } else {
                    if (norm_2(rIntersectionPoint1 - intersection_point) > Epsilon) { // Must be different from the first one
                        noalias(rIntersectionPoint2) = intersection_point;
                        solution = 1;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION;
                        break;
                    }
                }
            } else if (face_solution == 2) { // The line is coincident with the face
                ComputeTriangleLineIntersectionInTheSamePlane(r_face, rLinePoint1, rLinePoint2, rIntersectionPoint1, rIntersectionPoint2, solution, Epsilon);
                if (solution == 1) break;// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION) break;
            }
        }

        // Check if points are inside 
        if constexpr (TConsiderInsidePoints) {
            if (solution == 0) {// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION) {
                array_1d<double,3> local_coordinates;
                if (rTetrahedraGeometry.IsInside(rLinePoint1, local_coordinates)) {
                    noalias(rIntersectionPoint1) = rLinePoint1;
                    solution = 4;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_ONE_INSIDE;
                }
                if (rTetrahedraGeometry.IsInside(rLinePoint2, local_coordinates)) {
                    if (solution == 0) {// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::NO_INTERSECTION) {
                        noalias(rIntersectionPoint1) = rLinePoint2;
                        solution = 4;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_ONE_INSIDE;
                    } else {
                        noalias(rIntersectionPoint2) = rLinePoint2;
                        solution = 3;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_BOTH_INSIDE;
                    }
                }
            } else if (solution == 2) {// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION) {
                array_1d<double,3> local_coordinates;
                if (rTetrahedraGeometry.IsInside(rLinePoint1, local_coordinates)) {
                    if (norm_2(rIntersectionPoint1 - rLinePoint1) > Epsilon) { // Must be different from the first one
                        noalias(rIntersectionPoint2) = rLinePoint1;
                        solution = 4;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_ONE_INSIDE;
                    }
                } 
                if (solution == 2) {// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION) {
                    if (rTetrahedraGeometry.IsInside(rLinePoint2, local_coordinates)) {
                        if (norm_2(rIntersectionPoint1 - rLinePoint2) > Epsilon) {  // Must be different from the first one
                            noalias(rIntersectionPoint2) = rLinePoint2;
                            solution = 4;// IntersectionUtilitiesTetrahedraLineIntersectionStatus::TWO_POINTS_INTERSECTION_ONE_INSIDE;
                        }
                    }
                }
            }
        }

        // Checking if any node of the tetrahedra
        if (solution == 2) {// if (solution == IntersectionUtilitiesTetrahedraLineIntersectionStatus::ONE_POINT_INTERSECTION) {
            // Detect the node of the tetrahedra and directly assign
            int index_node = -1;
            for (int i_node = 0; i_node < 4; ++i_node) {
                if (norm_2(rTetrahedraGeometry[i_node].Coordinates() - rIntersectionPoint1) < Epsilon) {
                    index_node = i_node;
                    break;
                }
            }
            // Return index
            if (index_node > -1) {
                return index_node + 5;
                //return static_cast<IntersectionUtilitiesTetrahedraLineIntersectionStatus>(index_node + 5);
            }
        }

        return solution;
    }

    /**
     * @brief Calculates the line to line intersection (shortest line). If line is length 0, it is considered a point and therefore there is intersection (3D version)
     * @details Calculate the line segment PaPb that is the shortest route between two lines P1P2 and P3P4. Calculate also the values of mua and mub where
     *    Pa = P1 + mua (P2 - P1)
     *    Pb = P3 + mub (P4 - P3)
     *    http://paulbourke.net/geometry/pointlineplane/
     * @param rSegment1 The first segment
     * @param rSegment2 The second segment
     * @tparam TGeometryType The geometry type
     * @return Return empty points array if no solution exists. Otherwise returns the line intersection points array
     */
    template<class TGeometryType>
    static PointerVector<Point> ComputeShortestLineBetweenTwoLines(
        const TGeometryType& rSegment1,
        const TGeometryType& rSegment2
        )  
    {
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // Check geometry type
        KRATOS_ERROR_IF_NOT((rSegment1.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear && rSegment1.PointsNumber() == 2)) << "The first geometry type is not correct, it is suppossed to be a linear line" << std::endl;
        KRATOS_ERROR_IF_NOT((rSegment2.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear && rSegment2.PointsNumber() == 2)) << "The second geometry type is not correct, it is suppossed to be a linear line" << std::endl;

        // Resulting line segment
        auto resulting_line = PointerVector<Point>();

        // Variable definitions
        array_1d<double, 3> p13,p43,p21;
        double d1343,d4321,d1321,d4343,d2121;
        double mua, mub;
        double numer,denom;

        // Points segments
        const Point& p1 = rSegment1[0];
        const Point& p2 = rSegment1[1];
        const Point& p3 = rSegment2[0];
        const Point& p4 = rSegment2[1];

        p13[0] = p1.X() - p3.X();
        p13[1] = p1.Y() - p3.Y();
        p13[2] = p1.Z() - p3.Z();

        p43[0] = p4.X() - p3.X();
        p43[1] = p4.Y() - p3.Y();
        p43[2] = p4.Z() - p3.Z();
        if (std::abs(p43[0]) < zero_tolerance && std::abs(p43[1]) < zero_tolerance && std::abs(p43[2]) < zero_tolerance)
            return resulting_line;

        p21[0] = p2.X() - p1.X();
        p21[1] = p2.Y() - p1.Y();
        p21[2] = p2.Z() - p1.Z();
        if (std::abs(p21[0]) < zero_tolerance && std::abs(p21[1]) < zero_tolerance && std::abs(p21[2]) < zero_tolerance)
            return resulting_line;

        d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
        d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
        d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
        d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
        d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

        denom = d2121 * d4343 - d4321 * d4321;
        auto pa = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        auto pb = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
        if (std::abs(denom) < zero_tolerance) { // Parallel lines, infinite solutions. Projecting points and getting one perpendicular line
            // Projection auxiliary variables
            Point projected_point;
            array_1d<double,3> local_coords;
            // Projecting first segment
            GeometricalProjectionUtilities::FastProjectOnLine(rSegment2, rSegment1[0], projected_point);
            if (rSegment2.IsInside(projected_point, local_coords)) {
                pa->Coordinates() = rSegment1[0].Coordinates();
                pb->Coordinates() = projected_point;
            } else {
                GeometricalProjectionUtilities::FastProjectOnLine(rSegment2, rSegment1[1], projected_point);
                if (rSegment2.IsInside(projected_point, local_coords)) {
                    pa->Coordinates() = rSegment1[1].Coordinates();
                    pb->Coordinates() = projected_point;
                } else { // Trying to project second segment
                    GeometricalProjectionUtilities::FastProjectOnLine(rSegment1, rSegment2[0], projected_point);
                    if (rSegment1.IsInside(projected_point, local_coords)) {
                        pa->Coordinates() = rSegment2[0].Coordinates();
                        pb->Coordinates() = projected_point;
                    } else {
                        GeometricalProjectionUtilities::FastProjectOnLine(rSegment1, rSegment2[1], projected_point);
                        if (rSegment1.IsInside(projected_point, local_coords)) {
                            pa->Coordinates() = rSegment2[1].Coordinates();
                            pb->Coordinates() = projected_point;
                        } else { // Parallel and not projection possible
                            return resulting_line;
                        }
                    }
                }
            }
        } else {
            numer = d1343 * d4321 - d1321 * d4343;

            mua = numer / denom;
            mub = (d1343 + d4321 * mua) / d4343;

            pa->X() = p1.X() + mua * p21[0];
            pa->Y() = p1.Y() + mua * p21[1];
            pa->Z() = p1.Z() + mua * p21[2];
            pb->X() = p3.X() + mub * p43[0];
            pb->Y() = p3.Y() + mub * p43[1];
            pb->Z() = p3.Z() + mub * p43[2];
        }

        resulting_line.push_back(pa);
        resulting_line.push_back(pb);
        return resulting_line;
    }

    /**
     * Find the 2D intersection of two lines (both bounded)
     * @param rLineGeometry Is the line to intersect
     * @param rLinePoint1 Coordinates of the first point of the intersecting line
     * @param rLinePoint2 Coordinates of the second point of the intersecting line
     * @return rIntersectionPoint The intersection point coordinates
     * @return The intersection type index:
     * 0 (disjoint - no intersection)
     * 1 (intersect in a unique point)
     * 2 (overlap)
     * 3 (intersect in one endpoint)
     */
    template <class TGeometryType>
    static int ComputeLineLineIntersection(
        const TGeometryType& rLineGeometry,
        const array_1d<double,3>& rLinePoint0,
        const array_1d<double,3>& rLinePoint1,
        array_1d<double,3>& rIntersectionPoint,
        const double epsilon = 1e-12)
    {
        return ComputeLineLineIntersection(
            rLineGeometry[0], rLineGeometry[1], rLinePoint0, rLinePoint1, rIntersectionPoint, epsilon);
    }

    /**
     * Find the 2D intersection of two lines (both bounded)
     * @param rLine1Point0 Coordinates of the first point of the first line
     * @param rLine1Point1 Coordinates of the second point of the first line
     * @param rLine2Point0 Coordinates of the first point of the second line
     * @param rLine2Point1 Coordinates of the second point of the second line
     * @return rIntersectionPoint The intersection point coordinates
     * @return The intersection type index:
     * 0 (disjoint - no intersection)
     * 1 (intersect in a unique point)
     * 2 (overlap)
     * 3 (intersect in one endpoint)
     */
    static int ComputeLineLineIntersection(
        const array_1d<double,3>& rLine1Point0,
        const array_1d<double,3>& rLine1Point1,
        const array_1d<double,3>& rLine2Point0,
        const array_1d<double,3>& rLine2Point1,
        array_1d<double,3>& rIntersectionPoint,
        const double epsilon = 1e-12)
    {
        const array_1d<double,3> r = rLine1Point1 - rLine1Point0;
        const array_1d<double,3> s = rLine2Point1 - rLine2Point0;
        const array_1d<double,3> q_p = rLine2Point0 - rLine1Point0;        // q - p

        const double aux_1 = CrossProd2D(r,s);
        const double aux_2 = CrossProd2D(q_p,r);
        const double aux_3 = CrossProd2D(q_p,s);

        // Check intersection cases
        // Check that the lines are not collinear
        if (std::abs(aux_1) < epsilon && std::abs(aux_2) < epsilon){
            const double aux_4 = inner_prod(r,r);
            const double aux_5 = inner_prod(s,r);
            const double t_0 = inner_prod(q_p,r)/aux_4;
            const double t_1 = t_0 + aux_5/aux_4;
            if (aux_5 < 0.0){
                if (t_1 >= 0.0 && t_0 <= 1.0){
                    return 2;    // The lines are collinear and overlapping
                }
            } else {
                if (t_0 >= 0.0 && t_1 <= 1.0){
                    return 2;    // The lines are collinear and overlapping
                }
            }
        // Check that the lines are parallel and non-intersecting
        } else if (std::abs(aux_1) < epsilon && std::abs(aux_2) > epsilon){
            return 0;
        // Check that the lines intersect
        } else if (std::abs(aux_1) > epsilon){
            const double u = aux_2/aux_1;
            const double t = aux_3/aux_1;
            if (((u >= 0.0) && (u <= 1.0)) && ((t >= 0.0) && (t <= 1.0))){
                rIntersectionPoint = rLine2Point0 + u*s;
                // Check if the intersection occurs in one of the end points
                if (u < epsilon || (1.0 - u) < epsilon) {
                    return 3;
                } else {
                    return 1;
                }
            }
        }
        // Otherwise, the lines are non-parallel but do not intersect
        return 0;
    }

    /**
     * @brief Find the 3D intersection of a plane (infinite) with a segment (bounded)
     * @param rPlaneBasePoint Base point of the plane to intersect with
     * @param rPlaneNormal Normal vector of the plane to intersect with
     * @param rLinePoint1 Coordinates of the first point of the segment
     * @param rLinePoint2 Coordinates of the second point of the segment
     * @param rIntersectionPoint The intersection point coordinates
     * @return The intersection type index:
     * 0 (parallel or out of bounds - no intersection)
     * 1 (unique intersection point)
     * 2 (edge and plane coincide - no intersection)
     */
    static int ComputePlaneLineIntersection(
        const array_1d<double,3>& rPlaneBasePoint,
        const array_1d<double,3>& rPlaneNormal,
        const array_1d<double,3>& rLinePoint1,
        const array_1d<double,3>& rLinePoint2,
        array_1d<double,3>& rIntersectionPoint,
        const double epsilon = 1e-12)
    {
        // This is the adaption of the implementation provided in:
        // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
        // (Intersection of a Segment with a Plane)

        // Get direction vector of edge
        const array_1d<double,3> line_dir = rLinePoint2 - rLinePoint1;

        // Check if the segment is parallel to the plane or even coincides with it
        const double a = inner_prod(rPlaneNormal,( rPlaneBasePoint - rLinePoint1 ));
        const double b = inner_prod(rPlaneNormal,line_dir);
        if (std::abs(b) < epsilon){
            if (std::abs(a) < epsilon){
                return 2;    // Segment lies in the plane
            } else {
                return 0;    // Segment does not lie in the plane, but is parallel to it
            }
        }

        // Compute the intersection point and check if it is inside the bounds of the segment
        const double r = a / b;
        if (r < 0.0 - epsilon){
            return 0;    // Intersection point lies outside the bounds of the segment
        } else if (r > 1.0 + epsilon) {
            return 0;    // Intersection point lies outside the bounds of the segment
        }
        rIntersectionPoint = rLinePoint1 + r * line_dir;

        return 1;
    }

    /**
     * @brief Compute a segment box intersection
     * Provided the minimum and maximum points of a box cell, this method checks if
     * the segment intersects it. If it does intersect, it returns the rIntersectionPointpoint as well.
     * Note that the cell box is assumed to be aligned to the cartesian axes.
     * Adapted from: https://www.3dkingdoms.com/weekly/weekly.php?a=3
     * @param rBoxPoint0 Minimum point of the box cell
     * @param rBoxPoint1 Maximum point of the box cell
     * @param rLinePoint0 Segment origin point
     * @param rLinePoint1 Segment end point
     * @return int Returns 0 if there is no intersection and 1 otherwise
     */
    static int ComputeLineBoxIntersection(
        const array_1d<double,3> &rBoxPoint0,
        const array_1d<double,3> &rBoxPoint1,
        const array_1d<double,3> &rLinePoint0,
        const array_1d<double,3> &rLinePoint1)
    {
        array_1d<double,3> intersection_point = ZeroVector(3);

        if (rLinePoint1[0] < rBoxPoint0[0] && rLinePoint0[0] < rBoxPoint0[0]) return false;
        if (rLinePoint1[0] > rBoxPoint1[0] && rLinePoint0[0] > rBoxPoint1[0]) return false;
        if (rLinePoint1[1] < rBoxPoint0[1] && rLinePoint0[1] < rBoxPoint0[1]) return false;
        if (rLinePoint1[1] > rBoxPoint1[1] && rLinePoint0[1] > rBoxPoint1[1]) return false;
        if (rLinePoint1[2] < rBoxPoint0[2] && rLinePoint0[2] < rBoxPoint0[2]) return false;
        if (rLinePoint1[2] > rBoxPoint1[2] && rLinePoint0[2] > rBoxPoint1[2]) return false;
        if (rLinePoint0[0] > rBoxPoint0[0] && rLinePoint0[0] < rBoxPoint1[0] &&
            rLinePoint0[1] > rBoxPoint0[1] && rLinePoint0[1] < rBoxPoint1[1] &&
            rLinePoint0[2] > rBoxPoint0[2] && rLinePoint0[2] < rBoxPoint1[2]) {
            return true;
        }
        if ((GetLineBoxIntersection(rLinePoint0[0]-rBoxPoint0[0], rLinePoint1[0]-rBoxPoint0[0], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 1 )) ||
            (GetLineBoxIntersection(rLinePoint0[1]-rBoxPoint0[1], rLinePoint1[1]-rBoxPoint0[1], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 2 )) ||
            (GetLineBoxIntersection(rLinePoint0[2]-rBoxPoint0[2], rLinePoint1[2]-rBoxPoint0[2], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 3 )) ||
            (GetLineBoxIntersection(rLinePoint0[0]-rBoxPoint1[0], rLinePoint1[0]-rBoxPoint1[0], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 1 )) ||
            (GetLineBoxIntersection(rLinePoint0[1]-rBoxPoint1[1], rLinePoint1[1]-rBoxPoint1[1], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 2 )) ||
            (GetLineBoxIntersection(rLinePoint0[2]-rBoxPoint1[2], rLinePoint1[2]-rBoxPoint1[2], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 3 ))){
            return true;
        }

        return false;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This inline function computes the 2D cross product between two arrays
     * @param a First vector
     * @param b Second vector
     * @return The 2D cross product value
     */

    static inline double CrossProd2D(const array_1d<double,3> &a, const array_1d<double,3> &b){
        return (a(0)*b(1) - a(1)*b(0));
    }

    static inline int GetLineBoxIntersection(
        const double Dist1,
        const double Dist2,
        const array_1d<double,3> &rPoint1,
        const array_1d<double,3> &rPoint2,
        array_1d<double,3> &rIntersectionPoint)
    {
        if ((Dist1 * Dist2) >= 0.0){
            return 0;
        }
        // if ( Dist1 == Dist2) return 0;
        if (std::abs(Dist1-Dist2) < 1e-12){
            return 0;
        }
        rIntersectionPoint = rPoint1 + (rPoint2-rPoint1)*(-Dist1/(Dist2-Dist1));
        return 1;
    }

    static inline int InBox(
        const array_1d<double,3> &rIntersectionPoint,
        const array_1d<double,3> &rBoxPoint0,
        const array_1d<double,3> &rBoxPoint1,
        const unsigned int Axis)
    {
        if ( Axis==1 && rIntersectionPoint[2] > rBoxPoint0[2] && rIntersectionPoint[2] < rBoxPoint1[2] && rIntersectionPoint[1] > rBoxPoint0[1] && rIntersectionPoint[1] < rBoxPoint1[1]) return 1;
        if ( Axis==2 && rIntersectionPoint[2] > rBoxPoint0[2] && rIntersectionPoint[2] < rBoxPoint1[2] && rIntersectionPoint[0] > rBoxPoint0[0] && rIntersectionPoint[0] < rBoxPoint1[0]) return 1;
        if ( Axis==3 && rIntersectionPoint[0] > rBoxPoint0[0] && rIntersectionPoint[0] < rBoxPoint1[0] && rIntersectionPoint[1] > rBoxPoint0[1] && rIntersectionPoint[1] < rBoxPoint1[1]) return 1;
        return 0;
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class IntersectionUtilities */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/