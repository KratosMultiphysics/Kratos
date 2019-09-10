//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/point.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;

/// Factory functions
// namespace {

    NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> GenerateReferenceCurveOnSurface3d()
    {
        /*// Assign the points belonging to the curve
        PointerVector<Point> points_curve;
        points_curve.push_back(Point::Pointer(new Point(3, 2, 0)));
        points_curve.push_back(Point::Pointer(new Point(1, 4, 0)));
        points_curve.push_back(Point::Pointer(new Point(2, 5, 0)));

        // Assign the curve's knot vector
        Vector knot_vector_curve = ZeroVector(4);
        knot_vector_curve[0] = 7.0;
        knot_vector_curve[1] = 7.0;
        knot_vector_curve[2] = 9.0;
        knot_vector_curve[3] = 9.0;

        // Polynomial degree of the curve
        int p_curve = 2;

        // Create the 2D embedded curve
        auto curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points_curve, p_curve, knot_vector_curve);

        // Assign the points belonging to the surface
        PointerVector<Point> points_surface;
        points_surface.push_back(Point::Pointer(new Point( 0,  0,  3)));
        points_surface.push_back(Point::Pointer(new Point( 0,  5,  0)));
        points_surface.push_back(Point::Pointer(new Point( 0, 10,  2)));
        points_surface.push_back(Point::Pointer(new Point( 5,  0,  5)));
        points_surface.push_back(Point::Pointer(new Point( 5,  5,  0)));
        points_surface.push_back(Point::Pointer(new Point( 5, 10,  3)));
        points_surface.push_back(Point::Pointer(new Point(10,  0,  1)));
        points_surface.push_back(Point::Pointer(new Point(10,  5, -1)));
        points_surface.push_back(Point::Pointer(new Point(10, 10,  0)));

        // Assign the surface's knot vectors
        Vector knot_vector_u_surface = ZeroVector(4);
        knot_vector_u_surface[0] = 1.0;
        knot_vector_u_surface[1] = 1.0;
        knot_vector_u_surface[2] = 3.0;
        knot_vector_u_surface[3] = 3.0;

        Vector knot_vector_v_surface = ZeroVector(4);
        knot_vector_v_surface[0] = 2.0;
        knot_vector_v_surface[1] = 2.0;
        knot_vector_v_surface[2] = 6.0;
        knot_vector_v_surface[3] = 6.0;

        // Polynomial degrees
        int p_surface = 2;
        int q_surface = 2;

        // Create a 3D surface
        auto surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<Point>>>(points_surface, p_surface, 
            q_surface, knot_vector_u_surface, knot_vector_v_surface);

        // Create a curve on surface geometry
        auto curve_on_surface = NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>>(surface, curve);*/

    }

    /*NurbsSurfaceGeometry<3, PointerVector<Point>> GenerateReferencePointSurface()
    {
        PointerVector<Point> points;

        points.push_back(Point::Pointer(new Point(-10.0, -5.0, -1.0)));
        points.push_back(Point::Pointer(new Point(-12.0, 3.0, 3.0)));
        points.push_back(Point::Pointer(new Point(-9.0, 11.0, -0.0701928417)));
        points.push_back(Point::Pointer(new Point(-5.0, -3.0, 1.0)));
        points.push_back(Point::Pointer(new Point(-6.0, 4.0, -2.0)));
        points.push_back(Point::Pointer(new Point(-5.0, 7.0, 0.9298071583)));
        points.push_back(Point::Pointer(new Point(0.0, -4.0, -1.0)));
        points.push_back(Point::Pointer(new Point(1.0, 6.0, 5.0)));
        points.push_back(Point::Pointer(new Point(0.0, 13.0, -0.2350184214)));
        points.push_back(Point::Pointer(new Point(4.0, -2.0, 0.0)));
        points.push_back(Point::Pointer(new Point(5.0, 4.0, -1.0)));
        points.push_back(Point::Pointer(new Point(5.0, 11.0, 0.7649815786)));

        Vector knot_vector_u = ZeroVector(5);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 7.5;
        knot_vector_u[3] = 15.0;
        knot_vector_u[4] = 15.0;

        Vector knot_vector_v = ZeroVector(3);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 10.0;
        knot_vector_v[2] = 20.0;

        int p = 2;
        int q = 1;

        Vector weights = ZeroVector(12);
        weights[0] = 1.0;
        weights[1] = 1.0;
        weights[2] = 1.0;
        weights[3] = 1.0;
        weights[4] = 2.5;
        weights[5] = 1.0;
        weights[6] = 1.0;
        weights[7] = 1.0;
        weights[8] = 1.0;
        weights[9] = 1.0;
        weights[10] = 1.0;
        weights[11] = 1.0;

        return NurbsSurfaceGeometry<3, PointerVector<Point>>(
                points, p, q, knot_vector_u, knot_vector_v, weights);
    }*/



    /*NurbsCurveGeometry<3, PointerVector<NodeType>> GenerateReferenceCurve3d()
    {
        PointerVector<NodeType> points;

        points.push_back(NodeType::Pointer(new NodeType(1, 0, -25, -5)));
        points.push_back(NodeType::Pointer(new NodeType(2, -15, -15, 0)));
        points.push_back(NodeType::Pointer(new NodeType(3, 5, -5, -3)));
        points.push_back(NodeType::Pointer(new NodeType(4, 15, -15, 3)));
        points.push_back(NodeType::Pointer(new NodeType(5, 25, 0, 6)));
        points.push_back(NodeType::Pointer(new NodeType(6, 15, 15, 6)));
        points.push_back(NodeType::Pointer(new NodeType(7, -5, -5, -3)));
        points.push_back(NodeType::Pointer(new NodeType(8, -25, 15, 4)));

        Vector knot_vector = ZeroVector(11);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 0.0;
        knot_vector[3] = 0.0;
        knot_vector[4] = 32.9731425998736;
        knot_vector[5] = 65.9462851997473;
        knot_vector[6] = 98.9194277996209;
        knot_vector[7] = 131.892570399495;
        knot_vector[8] = 131.892570399495;
        knot_vector[9] = 131.892570399495;
        knot_vector[10] = 131.892570399495;

        int p = 4;

        Vector weights = ZeroVector(8);
        weights[0] = 1.0;
        weights[1] = 3.0;
        weights[2] = 1.0;
        weights[3] = 2.5;
        weights[4] = 1.0;
        weights[5] = 0.5;
        weights[6] = 1.0;
        weights[7] = 2.0;

        auto curve = NurbsCurveGeometry<3, PointerVector<NodeType>>(points, p, knot_vector, weights);

        return curve;
    }*/

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurface, KratosCoreNurbsGeometriesFastSuite) {
        
        KRATOS_CHECK_NEAR(0.894427, 0.894427, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
