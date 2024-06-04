
//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef RAY_AABB_PRIMITIVE_INCLUDE_H
#define RAY_AABB_PRIMITIVE_INCLUDE_H

//// Project includes
#include "queso/embedding/aabb_primitive.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Ray_AABB_primitive
 * @author Manuel Messmer
 * @brief  Ray to be used in AABB_tree.
 *         Provides functions to check for ray-triangle and ray-aabb intersections.
*/
class Ray_AABB_primitive : public AABB_primitive_base {

public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life cycle
    ///@{

    ///Constructor.
    ///@param Origin
    ///@param Direction Must be always positive: (x>0, y>0, z>0).
    Ray_AABB_primitive(const Vector3d& Origin, const Vector3d& Direction) :
        mOrigin(Origin), mDirection(Direction)
    {
        mInvDirection[0] = 1.0 / mDirection[0];
        mInvDirection[1] = 1.0 / mDirection[1];
        mInvDirection[2] = 1.0 / mDirection[2];

        mPositiveDir = true;
        if( Direction[0] <= 0 || Direction[1] <= 0 || Direction[2] <= 0 ){
            mPositiveDir = false;
            mSign[0] = (mInvDirection[0] < 0);
            mSign[1] = (mInvDirection[1] < 0);
            mSign[2] = (mInvDirection[2] < 0);
        }
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if ray intersects aabb.
    ///@param aabb
    ///@see intersect_general and intersect_positive_ray.
    ///@return bool
    bool intersect(const AABB_primitive &aabb) const override;

    /// @brief Returns true if ray intersect aabb.
    /// @param aabb
    /// @return bool
    bool intersect_general(const AABB_primitive &aabb) const;

    /// @brief Returns true if ray intersect aabb. Direction of ray must be positive (x>0,y>0,z>0).
    ///        This implementation is slightly faster than intersect_general.
    /// @param aabb
    /// @return bool
    bool intersect_positive_ray(const AABB_primitive &aabb) const;

    ///@brief Returns true if ray intersects triangle (only checks intersections in positive direction of ray).
    ///@param v0 Triangle Vertex 1
    ///@param v1 Triangle Vertex 2
    ///@param v2 Triangle Vertex 3
    ///@param t Distance to intersection.
    ///@param u Parametric coordinate 1.
    ///@param v Parametric coordinate 2.
    ///@param BackFacing True is triangle is back facing.
    ///@param Parallel True if triangle is parallel.
    ///@todo Put Backfacing and Parallel into one enum. Also add OriginIsOnBoundary.
    ///@return bool
    bool intersect( const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                    double &t, double &u, double &v, bool& BackFacing, bool& Parallel) const;

    ///@brief Returns true if triangle is parallel to ray.
    ///@param v0 Triangle Vertex 1
    ///@param v1 Triangle Vertex 2
    ///@param v2 Triangle Vertex 3
    ///@param Tolerance Default: ZEROTOL.
    ///@return bool
    bool is_parallel(const Vector3d &v0, const Vector3d &v1, const Vector3d &v2, double Tolerance = ZEROTOL ) const;

    ///@}
private:

    ///@name Private Members
    ///@{
    Vector3d mOrigin;
    Vector3d mDirection;
    Vector3d mInvDirection;
    Vector3i mSign;

    bool mPositiveDir;
    ///@}
}; // End Ray_AABB_primitive
///@} // End QuESo classes

} // End namespace queso

#endif