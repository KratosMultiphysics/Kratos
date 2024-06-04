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
#include "queso/embedding/aabb_primitive.h"

namespace queso {

bool AABB_primitive::intersect(const AABB_primitive &aabb) const  {
    for (unsigned int i = 0; i < 3; ++i) {
        if (aabb.upperBound[i] < lowerBound[i] || aabb.lowerBound[i] > upperBound[i] ) {
            return false;
        }
    }
    return true;
}

bool AABB_primitive::intersect(const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                               double tolerance) const {

    /// Get extent of aabb.
    Vector3d extent{(upperBound[0] - lowerBound[0] )/2.0 - tolerance,
                    (upperBound[1] - lowerBound[1] )/2.0 - tolerance,
                    (upperBound[2] - lowerBound[2] )/2.0 - tolerance};

    /// Translate triangle to origin.
    Vector3d v0_orig{v0[0] - centre[0], v0[1] - centre[1], v0[2] - centre[2]};
    Vector3d v1_orig{v1[0] - centre[0], v1[1] - centre[1], v1[2] - centre[2]};
    Vector3d v2_orig{v2[0] - centre[0], v2[1] - centre[1], v2[2] - centre[2]};

    // Compute the edge vectors of the triangle  (ABC). Line between vertices.
    Vector3d f0{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    Vector3d f1{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
    Vector3d f2{v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]};

    // Compute the face normals of the AABB, because the AABB
    // AABB is axis algined by definition.
    Vector3d u0{1.0, 0.0, 0.0};
    Vector3d u1{0.0, 1.0, 0.0};
    Vector3d u2{0.0, 0.0, 1.0};

    // There are a total of 13 axis to test!

    // First test (u0, u1, u2) vs. (f0, f1, f2). 9 tests in total.
    // u0 vs f0.
    Vector3d axis_u0_f0{u0[1]*f0[2] - u0[2]*f0[1], u0[2]*f0[0] - u0[0]*f0[2], u0[0]*f0[1] - u0[1]*f0[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f0) ){
        return false;
    }

    // u0 vs f1.
    Vector3d axis_u0_f1{u0[1]*f1[2] - u0[2]*f1[1], u0[2]*f1[0] - u0[0]*f1[2], u0[0]*f1[1] - u0[1]*f1[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f1) ){
        return false;
    }

    // u0 vs f2.
    Vector3d axis_u0_f2{u0[1]*f2[2] - u0[2]*f2[1], u0[2]*f2[0] - u0[0]*f2[2], u0[0]*f2[1] - u0[1]*f2[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f2) ){
        return false;
    }

    // u1 vs f0.
    Vector3d axis_u1_f0{u1[1]*f0[2] - u1[2]*f0[1], u1[2]*f0[0] - u1[0]*f0[2], u1[0]*f0[1] - u1[1]*f0[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f0) ){
        return false;
    }

    // u1 vs f1.
    Vector3d axis_u1_f1{u1[1]*f1[2] - u1[2]*f1[1], u1[2]*f1[0] - u1[0]*f1[2], u1[0]*f1[1] - u1[1]*f1[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f1) ){
        return false;
    }

    // u1 vs f2.
    Vector3d axis_u1_f2{u1[1]*f2[2] - u1[2]*f2[1], u1[2]*f2[0] - u1[0]*f2[2], u1[0]*f2[1] - u1[1]*f2[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f2) ){
        return false;
    }

    // u2 vs f0.
    Vector3d axis_u2_f0{u2[1]*f0[2] - u2[2]*f0[1], u2[2]*f0[0] - u2[0]*f0[2], u2[0]*f0[1] - u2[1]*f0[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f0) ){
        return false;
    }

    // u2 vs f1.
    Vector3d axis_u2_f1{u2[1]*f1[2] - u2[2]*f1[1], u2[2]*f1[0] - u2[0]*f1[2], u2[0]*f1[1] - u2[1]*f1[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f1) ){
        return false;
    }

    // u2 vs f2.
    Vector3d axis_u2_f2{u2[1]*f2[2] - u2[2]*f2[1], u2[2]*f2[0] - u2[0]*f2[2], u2[0]*f2[1] - u2[1]*f2[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f2) ){
        return false;
    }

    // Test face normals of aabb. 3 Tests.
    // axis1: (1, 0, 0)
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u0) ){
        return false;
    }
    // axis2: (0, 1, 0)
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u1) ){
        return false;
    }

    // axis3 (0, 0, 1)
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u2) ){
        return false;
    }

    // Test face normal of triangle
    Vector3d triangle_normal{f0[1]*f1[2] - f0[2]*f1[1], f0[2]*f1[0] - f0[0]*f1[2], f0[0]*f1[1] - f0[1]*f1[0]};
    if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, triangle_normal) ){
        return false;
    }

    // Return true of all checkes returned true.
    return true;
}


bool AABB_primitive::check_axis( const Vector3d &u0, const Vector3d &u1, const Vector3d &u2,
                                 const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                                 const Vector3d &extent, const Vector3d& test_axis ) const {

    // Project all 3 vertices of the triangle onto the test_axis.
    double pv0 = v0[0]*test_axis[0] + v0[1]*test_axis[1] + v0[2]*test_axis[2];
    double pv1 = v1[0]*test_axis[0] + v1[1]*test_axis[1] + v1[2]*test_axis[2];
    double pv2 = v2[0]*test_axis[0] + v2[1]*test_axis[1] + v2[2]*test_axis[2];

    // Project normals of aabb onto the test_axis.
    double pu0 = u0[0]*test_axis[0] + u0[1]*test_axis[1] + u0[2]*test_axis[2];
    double pu1 = u1[0]*test_axis[0] + u1[1]*test_axis[1] + u1[2]*test_axis[2];
    double pu2 = u2[0]*test_axis[0] + u2[1]*test_axis[1] + u2[2]*test_axis[2];

    // Project the AABB onto the seperating axis
    double r = extent[0] * std::abs(pu0) + extent[1] * std::abs(pu1) + extent[2] * std::abs(pu2);

    // Check if most extreme of the triangle points intersect r.
    if( std::max( -std::max({pv0, pv1, pv2}), std::min({pv0, pv1, pv2}) ) > r ){
        return false;
    }

    return true;
}

} // End namespace queso