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
    NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> GenerateReferenceBSplineCurveOnBSplineSurface3d()
    {

        // Assign the points belonging to the curve
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

        // Create and return a curve on surface geometry
        return NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>>(surface, curve);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(BSplineCurveOnBSplineSurface, KratosCoreNurbsGeometriesFastSuite) {
        
        auto curve_on_surface = GenerateReferenceBSplineCurveOnBSplineSurface3d();



        array_1d<double, 1> parameter(0.0);
        parameter[0] = 8.0;

        auto derivatives = curve_on_surface.GlobalDerivatives(parameter, 0);

        KRATOS_WATCH(derivatives)
        KRATOS_WATCH(" ")

        KRATOS_CHECK_NEAR(0.894427, 0.894427, TOLERANCE);

        KRATOS_WATCH("Here I am !!!")

        exit(-1);
    }

} // namespace Testing.
} // namespace Kratos.
