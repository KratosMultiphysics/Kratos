//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// System includes
#include <limits>

// External includes

// Project includes
#include "geometries/triangle_3d_3.h"
#include "includes/kratos_application.h"
#include "modeler/internals/cartesian_ray.h"
#include "testing/testing.h"

namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(CartesianRayAddIntersectionSmallTri, KratosCoreFastSuite)
{
    Point point_1(2.22931774067e-05, 3.51532295186e-05, -3.080646269e-05);
    Point point_2(2.22931774067e-05, 3.51532295186e-05, 0.00033179965779);

    Point::Pointer p_point_tri_0 = Kratos::make_shared<Point>(2.3029014e-05,3.5707878e-05,-1.943157e-08);
    Point::Pointer p_point_tri_1 = Kratos::make_shared<Point>(2.0441657e-05,3.3120339e-05,-1.943157e-08);
    Point::Pointer p_point_tri_2 = Kratos::make_shared<Point>(2.0441657e-05,3.5707878e-05,-1.943157e-08);

    Triangle3D3<Point> triangle(p_point_tri_0, p_point_tri_1, p_point_tri_2);

    Internals::CartesianRay<Triangle3D3<Point>> ray(1, point_1, point_2);
    const double tol = 1e-12;
    ray.AddIntersection(triangle, tol);

    KRATOS_EXPECT_EQ(ray.GetIntersections().size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(CartesianRayAddIntersectionParallel, KratosCoreFastSuite)
{
    Point point_1(1, 1, 2);
    Point point_2(2, 3, 2);

    Point::Pointer p_point_tri_0 = Kratos::make_shared<Point>(0,0,0);
    Point::Pointer p_point_tri_1 = Kratos::make_shared<Point>(1,1,0);
    Point::Pointer p_point_tri_2 = Kratos::make_shared<Point>(0,1,0);

    Triangle3D3<Point> triangle(p_point_tri_0, p_point_tri_1, p_point_tri_2);

    Internals::CartesianRay<Triangle3D3<Point>> ray(1, point_1, point_2);
    const double tol = 1e-12;
    ray.AddIntersection(triangle, tol);

    KRATOS_EXPECT_EQ(ray.GetIntersections().size(), 0);
}

} // namespace Kratos::Testing
