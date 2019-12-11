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
//                   Tobias Teschemacher
//                   Andreas Apostolatos
//
//  Ported partly from the ANurbs library (https://github.com/oberbichler/ANurbs)
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
    // Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
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
        points_surface.push_back(Point::Pointer(new Point( 5,  0,  5)));
        points_surface.push_back(Point::Pointer(new Point(10,  0,  1)));
        points_surface.push_back(Point::Pointer(new Point( 0,  5,  0)));
        points_surface.push_back(Point::Pointer(new Point( 5,  5,  0)));
        points_surface.push_back(Point::Pointer(new Point(10,  5, -1)));
        points_surface.push_back(Point::Pointer(new Point( 0, 10,  2)));
        points_surface.push_back(Point::Pointer(new Point( 5, 10,  3)));
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

    NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> GenerateReferenceNurbsCurveOnNurbsSurface3d()
    {
        // Assign the points belonging to the curve
        PointerVector<Point> points_curve;
        points_curve.push_back(Point::Pointer(new Point(-9, -2, 0)));
        points_curve.push_back(Point::Pointer(new Point(-2, -4, 0)));
        points_curve.push_back(Point::Pointer(new Point(2, 4, 0)));
        points_curve.push_back(Point::Pointer(new Point(9, -2, 0)));

        // Assign the curve's knot vector
        Vector knot_vector_curve = ZeroVector(5);
        knot_vector_curve[0] = -1.0;
        knot_vector_curve[1] = -1.0;
        knot_vector_curve[2] = 4.5;
        knot_vector_curve[3] = 9.0;
        knot_vector_curve[4] = 9.0;

        // Assign the curve's weights
        Vector weights_curve = ZeroVector(4);
        weights_curve[0] = 1.0;
        weights_curve[1] = 3.9;
        weights_curve[2] = 6.1;
        weights_curve[3] = 1.0;

        // Polynomial degree of the curve
        int p_curve = 2;

        // Create the 2D embedded curve
        auto curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points_curve, p_curve, knot_vector_curve, weights_curve);

        // Assign the points belonging to the surface
        PointerVector<Point> points_surface;
        points_surface.push_back(Point::Pointer(new Point(0, -0.075, 0)));
        points_surface.push_back(Point::Pointer(new Point(0.075, -0.075, 0)));
        points_surface.push_back(Point::Pointer(new Point(0.075, 0, 0)));

        points_surface.push_back(Point::Pointer(new Point(0, -0.075, 0.075)));
        points_surface.push_back(Point::Pointer(new Point(0.075, -0.075, 0.075)));
        points_surface.push_back(Point::Pointer(new Point(0.075, 0, 0)));

        points_surface.push_back(Point::Pointer(new Point(0, 0, 0.075)));
        points_surface.push_back(Point::Pointer(new Point(0.075, 0, 0.075)));
        points_surface.push_back(Point::Pointer(new Point(0.075, 0, 0)));

        // Assign the weights belonging to the surface
        Vector weights_surface = ZeroVector(9);
        weights_surface[0] = 1;
        weights_surface[1] = 7.071067811865476e-01;
        weights_surface[2] = 1;

        weights_surface[3] = 7.071067811865476e-01;
        weights_surface[4] = 0.5;
        weights_surface[5] = 7.071067811865476e-01;

        weights_surface[6] = 1;
        weights_surface[7] = 7.071067811865476e-01;
        weights_surface[8] = 1;

        // Assign the surface's knot vectors
        Vector knot_vector_u_surface = ZeroVector(4);
        knot_vector_u_surface[0] = -10.0;
        knot_vector_u_surface[1] = -10.0;
        knot_vector_u_surface[2] = 10.0;
        knot_vector_u_surface[3] = 10.0;

        Vector knot_vector_v_surface = ZeroVector(4);
        knot_vector_v_surface[0] = -10.0;
        knot_vector_v_surface[1] = -10.0;
        knot_vector_v_surface[2] = 10.0;
        knot_vector_v_surface[3] = 10.0;

        // Polynomial degrees
        int p_surface = 2;
        int q_surface = 2;

        // Create a 3D surface
        auto surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<Point>>>(points_surface, p_surface, 
            q_surface, knot_vector_u_surface, knot_vector_v_surface, weights_surface);

        // Create and return a curve on surface geometry
        return NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>>(surface, curve);
    }

    ///// Tests
    // Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
    KRATOS_TEST_CASE_IN_SUITE(BSplineCurveOnBSplineSurface, KratosCoreNurbsGeometriesFastSuite) 
    {    
        // Create a B-Spline curve on a B-Spline surface
        auto curve_on_surface = GenerateReferenceBSplineCurveOnBSplineSurface3d();

        // Choose a parametric coordinate on the curve
        array_1d<double, 3> parameter = ZeroVector(3);
        parameter[0] = 8.0;

        // Compute the derivatives of the gradient of the curve on surface up to 2nd order
        std::vector<array_1d<double, 3>> derivatives;
        curve_on_surface.GlobalSpaceDerivatives(derivatives, parameter, 2);

        // Compare the position vectors and the gradients up to 2nd order
        std::vector<double> positionVct = {3.75, 4.375, 1.5063476563};
        std::vector<double> gradient1 = {-2.5, 3.75, -0.658203125};
        std::vector<double> gradient2 = {7.5, -1.25, 1.1621094};

        KRATOS_CHECK_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[2], gradient2, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnNurbsSurface, KratosCoreNurbsGeometriesFastSuite) 
    {
        // Create a Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateReferenceNurbsCurveOnNurbsSurface3d();

        // Evaluate the global coordinates of the curve at u = 4.0
        {
            array_1d<double, 3> parameter = ZeroVector(3);
            parameter[0] = 4.0;

            // Compute the global coordinates of the curve on surface
            std::vector<array_1d<double, 3>> global_coordinates;
            curve_on_surface.GlobalSpaceDerivatives(global_coordinates, parameter, 1);

            // Compare the position vectors
            std::vector<double> positionVct = {0.054174511426802, -0.034958783455029, 0.038314563432352};

            KRATOS_CHECK_VECTOR_NEAR(global_coordinates[0], positionVct, TOLERANCE);
        }

        // Evaluate the gradients of the curve at u = 0.0
        {
            array_1d<double, 3> parameter = ZeroVector(3);
            parameter[0] = 0.0;

            // Compute the derivatives of the gradient of the curve on surface up to 2nd order
            std::vector<array_1d<double, 3>> derivatives;
            curve_on_surface.GlobalSpaceDerivatives(derivatives, parameter, 1);

            // Compare the position vectors and the gradients up to 2nd order
            std::vector<double> positionVct = {0.032428916974017, -0.057744163816062, 0.035199103526597};
            std::vector<double> gradient1 = {0.013539362404523, 0.005729880482752, -0.003073933458934};

            KRATOS_CHECK_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
            KRATOS_CHECK_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
        }

        // Evaluate the gradients of the curve at u = 7.2109
        {
            array_1d<double, 3> parameter = ZeroVector(3);
            parameter[0] = 7.2109;

            // Compute the derivatives of the gradient of the curve on surface up to 2nd order
            std::vector<array_1d<double, 3>> derivatives;
            curve_on_surface.GlobalSpaceDerivatives(derivatives, parameter, 1);

            // Compare the position vectors and the gradients up to 2nd order
            std::vector<double> positionVct = {0.062263623548203, -0.021647529053958, 0.035771855815790};
            std::vector<double> gradient1 = {0.002960292593061, 0.001993598876310, -0.003946176422519};

            KRATOS_CHECK_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
            KRATOS_CHECK_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
        }
    }

} // namespace Testing.
} // namespace Kratos.
