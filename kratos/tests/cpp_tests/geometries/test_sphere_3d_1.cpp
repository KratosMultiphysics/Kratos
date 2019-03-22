//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/sphere_3d_1.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Factory functions

    /** Generates a point type sample Sphere3D1
    * @return Pointer to a Sphere3D1
    */
    Sphere3D1<Point>::Pointer GenerateSphere3D1CenterOrigin() {
        return Kratos::make_shared<Sphere3D1<Point>>(Kratos::make_shared<Point>(0.0, 0.0, 0.0));
    }

    Sphere3D1<Point>::Pointer GenerateSphere3D1CenterUnitXDirection() {
        return Kratos::make_shared<Sphere3D1<Point>>(Kratos::make_shared<Point>(1.0, 0.0, 0.0));
    }

    /** Generates a point type sample Line3Sphere3D1D2N.
    * @return Pointer to a Sphere3D1
    */
    Sphere3D1<Point>::Pointer GenerateSphere3D1WithPoint(Point::Pointer rPoint ) {
        return Kratos::make_shared<Sphere3D1<Point>>(rPoint);
    }

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Sphere3dEdgesNumber, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateSphere3D1CenterOrigin();
        KRATOS_CHECK_EQUAL(p_geom->EdgesNumber(), 0);
    }

    /** Checks if the number of faces is correct.
    * Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Sphere3dFacesNumber, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateSphere3D1CenterOrigin();
        KRATOS_CHECK_EQUAL(p_geom->FacesNumber(), 1);
    }

    /** Checks if the length of the line is calculated correctly.
    * Checks if the length of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Sphere3dSetRadius, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateSphere3D1CenterOrigin();
        double radius = 1.0;
        p_geom->SetRadius(radius);
        KRATOS_CHECK_EQUAL(p_geom->GetRadius(), 1.0);
        radius = 2.0;
        KRATOS_CHECK_EQUAL(p_geom->GetRadius(), 2.0);
    }

    /** Checks that the copy constructor works as expected.
     *  The new sphere should not share the same radius with the original one
    * Checks if the length of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(SphereCopyConstructor, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateSphere3D1CenterOrigin();
        double radius = 1.0;
        p_geom->SetRadius(radius);
        // use copy-constructor to create a new sphere
        auto p_cloned_geom = Kratos::make_shared<Sphere3D1<Point>>(*p_geom);
        double big_radius = 2.0;
        p_cloned_geom->SetRadius(big_radius);
        const double volume = 4.0/3 * Globals::Pi;
        const double big_volume = 32.0/3 * Globals::Pi;
        KRATOS_CHECK_DOUBLE_EQUAL(p_geom->Volume(), volume);
        KRATOS_CHECK_DOUBLE_EQUAL(p_cloned_geom->Volume(), big_volume);
    }

    /** Checks that the HasIntersection method is well implemented for
    * sphere-bounding box intersections.
    */
    KRATOS_TEST_CASE_IN_SUITE(Sphere3dTestIntersectionWithBox, KratosCoreGeometriesFastSuite) {
        const auto low_point = Point(-1.0, -1.0, -1.0);
        const auto high_point = Point(1.0, 1.0, 1.0);

        // test sphere that is fully inside:
        double radius = 0.5;
        auto p_geom_inside = GenerateSphere3D1CenterOrigin();
        p_geom_inside->SetRadius(radius);
        const bool has_intersection_centered_fully_inside = p_geom_inside->HasIntersection(low_point, high_point);
        KRATOS_CHECK(has_intersection_centered_fully_inside);

        // test sphere that is inside, but bigger than the box:
        radius = 2.0;
        const bool has_intersection_centered_partially_inside = p_geom_inside->HasIntersection(low_point, high_point);
        KRATOS_CHECK(has_intersection_centered_partially_inside);

        // test sphere that intersects with one face only:
        auto p_geom_one_face = GenerateSphere3D1WithPoint(Kratos::make_shared<Point>(2.0, 0.0, 0.0));
        radius = 1.1;
        p_geom_one_face->SetRadius(radius);
        const bool has_intersection_with_one_face = p_geom_one_face->HasIntersection(low_point, high_point);
        KRATOS_CHECK(has_intersection_with_one_face);

        // test sphere that is too small to intersect with one face:
        radius = 0.9;
        const bool does_not_intersect_with_one_face = !p_geom_one_face->HasIntersection(low_point, high_point);
        KRATOS_CHECK(does_not_intersect_with_one_face);

        // test sphere that intersects with one edge only:
        auto p_geom_one_edge = GenerateSphere3D1WithPoint(Kratos::make_shared<Point>(-2.0, 2.0, 0.0));
        radius = 1.5;
        p_geom_one_edge->SetRadius(radius);
        const bool has_intersection_with_one_edge = p_geom_one_edge->HasIntersection(low_point, high_point);
        KRATOS_CHECK(has_intersection_with_one_edge);

        // test sphere that is too small to intersect with one edge:
        radius = 1.0;
        const bool does_not_intersect_with_one_edge = !p_geom_one_edge->HasIntersection(low_point, high_point);
        KRATOS_CHECK(does_not_intersect_with_one_edge);

        // test sphere that intersects with one corner only:
        auto p_geom_one_corner = GenerateSphere3D1WithPoint(Kratos::make_shared<Point>(2.0, 2.0, 2.0));
        radius = 1.8;
        p_geom_one_corner->SetRadius(radius);
        const bool has_intersection_with_one_corner = p_geom_one_corner->HasIntersection(low_point, high_point);
        KRATOS_CHECK(has_intersection_with_one_corner);

        // test sphere that is too small to intersect with one corner:
        radius = 1.5;
        const bool does_not_intersect_with_one_corner = !p_geom_one_corner->HasIntersection(low_point, high_point);
        KRATOS_CHECK(does_not_intersect_with_one_corner);
    }

    /** Checks that an exception is properly thrown either when
     * doing 'GetRadius' before initializing the radius or when
     * doing 'SetRadius' for a second time.
    */
    KRATOS_TEST_CASE_IN_SUITE(SphereRadiusGetSetExceptions, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateSphere3D1CenterOrigin();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(p_geom->GetRadius(), "The Sphere3D1's radius address has not been properly set before calling the method 'GetRadius'. Use 'SetRadius' for that.");
        double radius = 1.0;
        p_geom->SetRadius(radius);
        KRATOS_CHECK_EXCEPTION_IS_THROWN(p_geom->SetRadius(radius), "Trying to set the sphere's radius for a second time to a value other than NULL. This operation must be done only once.");
    }

} // namespace Testing.
} // namespace Kratos.
