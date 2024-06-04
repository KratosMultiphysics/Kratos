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

#ifndef GEOMETRY_QUERY_INCLUDE_H
#define GEOMETRY_QUERY_INCLUDE_H

/// Project includes
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/ray_aabb_primitive.h"
#include "queso/embedding/aabb_primitive.h"
#include "queso/embedding/aabb_tree.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Polygon
 * @author Manuel Messmer
 * @brief Provides functions to query geometry tests. Uses AABB tree to accelerate this process.
*/
class GeometryQuery {

public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    GeometryQuery(const TriangleMeshInterface& rTriangleMesh, bool MeshIsClosed = true)
        : mTriangleMesh(rTriangleMesh), mTree(rTriangleMesh), mMeshIsClosed(MeshIsClosed)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns true, if rPoint is within the bounding box of the given triangle mesh.
    /// @param rPoint
    /// @return bool.
    bool IsWithinBoundingBox(const PointType& rPoint) const;

    /// @brief Ray tracing to check, if a Point is inside or outside of the given triangle mesh.
    /// @details Calls: IsInsideOpen or IsInsideClosed depending on mMeshIsClosed.
    /// @param rRay Ray.
    /// @return std::pair<bool, bool> first-is_inside second-test_successful.
    std::pair<bool, bool> IsInside( const Ray_AABB_primitive& rRay ) const;

    /// @brief Returns true, if the AABB intersects with the given triangle mesh.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param Tolerance Reduces size of AABB.
    /// @return bool
    bool DoIntersect(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const;

    /// @brief Returns a vector of ids of all triangles that intersect with the AABB.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param Tolerance Reduces size of AABB.
    /// @return Unique<std::vector<IndexType>>
    Unique<std::vector<IndexType>> GetIntersectedTriangleIds(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const;

private:
    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Perfoms classical ray tracing. IsInside=true, if the number of intersected triangles is odd.
    /// @param rRay Ray.
    /// @return std::pair<bool, bool> first-is_inside second-test_successful
    std::pair<bool, bool> IsInsideClosed( const Ray_AABB_primitive& rRay ) const;

    /// @brief IsInside=true, if the nearest triangle is back facing.
    /// @param rRay Ray.
    /// @return std::pair<bool, bool> first-is_inside second-test_successful
    std::pair<bool, bool> IsInsideOpen( const Ray_AABB_primitive& rRay ) const;

    ///@}
    ///@name Private Members
    ///@{

    const TriangleMeshInterface& mTriangleMesh;
    AABB_tree mTree;
    bool mMeshIsClosed;
    ///@}

}; // End GeometryQuery class
///@} End QuESo classes

} // End namespace queso

#endif // GEOMETRY_QUERY_INCLUDE_H