//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// Project includes
#include "queso/embedding/geometry_query.h"

namespace queso {

    bool GeometryQuery::IsWithinBoundingBox(const PointType& rPoint) const {
        return mTree.IsWithinBoundingBox(rPoint);
    }

    std::pair<bool, bool> GeometryQuery::IsInside( const Ray_AABB_primitive& rRay ) const {
        if( mMeshIsClosed ){
            return IsInsideClosed(rRay);
        } else {
            return IsInsideOpen(rRay);
        }
    }

    bool GeometryQuery::DoIntersect(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const {
        AABB_primitive aabb(rLowerBound, rUpperBound);
        auto result = mTree.Query(aabb);

        const double snap_tolerance = RelativeSnapTolerance(rLowerBound, rUpperBound, Tolerance);
        for( auto r : result){
            const auto& p1 = mTriangleMesh.P1(r);
            const auto& p2 = mTriangleMesh.P2(r);
            const auto& p3 = mTriangleMesh.P3(r);
            if( aabb.intersect(p1, p2, p3, snap_tolerance) ){
                return true;
            }
        }
        return false;
    }

    Unique<std::vector<IndexType>> GeometryQuery::GetIntersectedTriangleIds(
            const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const {

        // Perform fast search based on aabb tree. Conservative search.
        AABB_primitive aabb(rLowerBound, rUpperBound);
        auto potential_intersections = mTree.Query(aabb);

        auto intersected_triangle_ids = MakeUnique<std::vector<IndexType>>();
        intersected_triangle_ids->reserve(potential_intersections.size());
        for( auto triangle_id : potential_intersections){
            const auto& p1 = mTriangleMesh.P1(triangle_id);
            const auto& p2 = mTriangleMesh.P2(triangle_id);
            const auto& p3 = mTriangleMesh.P3(triangle_id);

            // Perform actual intersection test.
            if( aabb.intersect(p1, p2, p3, Tolerance) ){
                intersected_triangle_ids->push_back(triangle_id);
            }
        }

        return intersected_triangle_ids;
    }

    std::pair<bool, bool> GeometryQuery::IsInsideOpen( const Ray_AABB_primitive& rRay ) const {
        double min_distance = MAXD;
        bool is_inside = false;
        auto potential_intersections = mTree.Query(rRay);
        for( auto r : potential_intersections){
            const auto& p1 = mTriangleMesh.P1(r);
            const auto& p2 = mTriangleMesh.P2(r);
            const auto& p3 = mTriangleMesh.P3(r);
            double t, u, v;
            bool back_facing, parallel;
            if( rRay.intersect(p1, p2, p3, t, u, v, back_facing, parallel) ) {
                if( !parallel ){
                    double sum_u_v = u+v;
                    if( t < ZEROTOL ){ // origin lies on boundary
                        return std::make_pair(false, true);
                    }
                    // Ray shoots through boundary.
                    if( u < 0.0+ZEROTOL || v < 0.0+ZEROTOL || sum_u_v > 1.0-ZEROTOL ){
                        return std::make_pair(false, false);
                    }
                    if( t < min_distance ){
                        is_inside = back_facing;
                        min_distance = t;
                    }
                }
            }
        }

        return std::make_pair(is_inside, true);
    }

    std::pair<bool, bool> GeometryQuery::IsInsideClosed( const Ray_AABB_primitive& rRay ) const {
        // Get potential ray intersections from AABB tree.
        auto potential_intersections = mTree.Query(rRay);
        IndexType intersection_count = 0;
        for( auto r : potential_intersections){
            const auto& p1 = mTriangleMesh.P1(r);
            const auto& p2 = mTriangleMesh.P2(r);
            const auto& p3 = mTriangleMesh.P3(r);
            double t, u, v;
            bool back_facing, parallel;
            // Test for actual intersection, if area is area is not zero.
            if( mTriangleMesh.Area(r) > 100.0*ZEROTOL
                    && rRay.intersect(p1, p2, p3, t, u, v, back_facing, parallel) ) {
                intersection_count++;
                const double sum_u_v = u+v;
                if( t < ZEROTOL ){ // Origin lies on boundary
                    return std::make_pair(false, true);
                }
                // Ray shoots through boundary of triangle.
                if( u < 0.0+ZEROTOL || v < 0.0+ZEROTOL || sum_u_v > 1.0-ZEROTOL ){
                    return std::make_pair(false, false);
                }
                if( parallel ){ // Triangle is parallel to ray.
                    return std::make_pair(false, false);
                }
            }
        }
        if( intersection_count % 2 == 1){ // If intersection count is odd, return true.
            return std::make_pair(true, true);
        } else {
            return std::make_pair(false, true);
        }
    }

} // End queso namespace