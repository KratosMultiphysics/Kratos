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

//// Project includes
#include "queso/embedding/ray_aabb_primitive.h"
#include "queso/utilities/math_utilities.hpp"

namespace queso {

// bool Ray_AABB_primitive::inside(const AABB_primitive &aabb) const {
//     if(    mOrigin[0] >= aabb.lowerBound[0]
//         || mOrigin[1] >= aabb.lowerBound[1]
//         || mOrigin[2] >= aabb.lowerBound[2]
//         || mOrigin[0] <= aabb.upperBound[0]
//         || mOrigin[1] <= aabb.upperBound[1]
//         || mOrigin[2] <= aabb.upperBound[2] ) {
//             return true;
//         }

//     return true;
// }

bool Ray_AABB_primitive::intersect(const AABB_primitive &aabb) const
{
    if( mPositiveDir ){
        // More efficient, but expects ray direction to be in positive direction: x>0, y>0, z>0.
        return intersect_positive_ray(aabb);
    }
    else {
        /// Works with all Ray directions.
        return intersect_general(aabb);
    }
}

bool Ray_AABB_primitive::intersect_positive_ray(const AABB_primitive &aabb) const {
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    double lower_0 = aabb.lowerBound[0];
    double lower_1 = aabb.lowerBound[1];

    double upper_0 = aabb.upperBound[0];
    double upper_1 = aabb.upperBound[1];

    double origin_0 = mOrigin[0];
    double origin_1 = mOrigin[1];

    double inv_direction_0 = mInvDirection[0];
    double inv_direction_1 = mInvDirection[1];

    double lower_2 = aabb.lowerBound[2];
    double upper_2 = aabb.upperBound[2];
    double origin_2 = mOrigin[2];

    // Check if origin of lies inside aabb.
    if(    origin_0 >= lower_0
        && origin_1 >= lower_1
        && origin_2 >= lower_2
        && origin_0 <= upper_0
        && origin_1 <= upper_1
        && origin_2 <= upper_2 ) {
            return true;
        }

    tmin = (lower_0 - origin_0) * inv_direction_0;
    tymax = (upper_1 - origin_1) * inv_direction_1;
    if(tmin > tymax) {
        return false;
    }

    tmax = (upper_0 - origin_0) * inv_direction_0;
    tymin = (lower_1 - origin_1) * inv_direction_1;


    if( tymin > tmax )
        return false;

    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    double inv_direction_2 = mInvDirection[2];
    tzmin = (lower_2 - origin_2) * inv_direction_2;
    if((tzmin > tmax))
        return false;

    tzmax = (upper_2 - origin_2) * inv_direction_2;
    if( tmin > tzmax )
        return false;

    if (tzmin > tmin)
        tmin = tzmin;

    if (tmin < 0) {
        return false;
    }

    return true;
}

bool Ray_AABB_primitive::intersect_general(const AABB_primitive &aabb) const {
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    std::array<Vector3d, 2> bounds = {aabb.lowerBound, aabb.upperBound};
    double lower_0 = bounds[mSign[0]][0];
    double lower_1 = bounds[mSign[1]][1];

    double upper_0 = bounds[1-mSign[0]][0];
    double upper_1 = bounds[1-mSign[1]][1];

    double origin_0 = mOrigin[0];
    double origin_1 = mOrigin[1];

    double inv_direction_0 = mInvDirection[0];
    double inv_direction_1 = mInvDirection[1];

    double lower_2 = bounds[mSign[2]][2];
    double upper_2 = bounds[1-mSign[2]][2];
    double origin_2 = mOrigin[2];

    //Check if origin of lies inside aabb.
    if(    origin_0 >= aabb.lowerBound[0]
        && origin_1 >= aabb.lowerBound[1]
        && origin_2 >= aabb.lowerBound[2]
        && origin_0 <= aabb.upperBound[0]
        && origin_1 <= aabb.upperBound[1]
        && origin_2 <= aabb.upperBound[2] ) {
            return true;
        }

    tmin = (lower_0 - origin_0) * inv_direction_0;
    tymax = (upper_1 - origin_1) * inv_direction_1;
    if(tmin > tymax) {
        return false;
    }

    tmax = (upper_0 - origin_0) * inv_direction_0;
    tymin = (lower_1 - origin_1) * inv_direction_1;


    if( tymin > tmax )
        return false;

    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    double inv_direction_2 = mInvDirection[2];
    tzmin = (lower_2 - origin_2) * inv_direction_2;
    if((tzmin > tmax))
        return false;

    tzmax = (upper_2 - origin_2) * inv_direction_2;
    if( tmin > tzmax )
        return false;

    if (tzmin > tmin)
        tmin = tzmin;

    if (tmin < 0) {
        return false;
    }

    return true;
}

bool Ray_AABB_primitive::intersect( const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                double &t, double &u, double &v, bool& BackFacing, bool& Parallel) const {

    // Substraction: v1-v0 and v2-v0
    Vector3d v0v1 = Math::Subtract( v1, v0 );
    Vector3d v0v2 = Math::Subtract( v2, v0 );

    // Cross product: mDirection x v0v2
    Vector3d pvec = Math::Cross(mDirection, v0v2);

    // Dot product: v0v1 * pvec
    double det = Math::Dot(v0v1, pvec);

    // Set default flag values
    BackFacing = false;
    Parallel = false;

    // If det is zero, triangle is parallel to ray.
    if (std::abs(det) < 10.0*ZEROTOL){
        Parallel = true;
        return true;
    }

    // If det is smaller than zero triangle is back facing.
    if( det < ZEROTOL )
        BackFacing = true;

    // Get inverse of determinant.
    double invDet = 1 / det;

    // Substraction: mOrigin - v0
    Vector3d tvec = Math::Subtract( mOrigin, v0 );

    // Dot product x invDet: (tvec * pvec) * invDet
    u = Math::Dot(tvec, pvec) * invDet;

    if (u < -ZEROTOL || u > 1+ZEROTOL)
        return false;

    // Cross product: tvec x v0v1
    Vector3d qvec = Math::Cross(tvec, v0v1);

    // Dot product x invDet: (mDirection * qvec) * invDet
    v = Math::Dot(mDirection, qvec) * invDet;

    if (v < -ZEROTOL || u + v > 1+ZEROTOL)
        return false;

    // Dot product x invDet: (v0v2 * qvec) * invDet
    t = Math::Dot(v0v2, qvec) * invDet;

    // Return false if ray intersects in negative direction.
    if( t < -ZEROTOL )
        return false;

    return true;
}

bool Ray_AABB_primitive::is_parallel( const Vector3d &v0, const Vector3d &v1, const Vector3d &v2, double Tolerance) const {
    // Substraction: v1-v0 and v2-v0
    Vector3d v0v1 = Math::Subtract( v1, v0 );
    Vector3d v0v2 = Math::Subtract( v2, v0 );

    // Cross product: mDirection x v0v2
    Vector3d pvec = Math::Cross(mDirection, v0v2);

    // Dot product: v0v1 * pvec
    double det = Math::Dot(v0v1, pvec);

    // If det is zero, triangle is parallel to ray.
    if (std::abs(det) < Tolerance){
        return true;
    }
    return false;
}

} // End namespace queso