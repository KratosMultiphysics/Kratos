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

typedef Node NodeType;

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

    NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> GenerateLineNurbsCurveOnSurface3d()
    {
        PointerVector<Point> points_curve(2);
        points_curve(0) = Kratos::make_shared<Point>(5.58207546402288, 3.9936436484419833, 0);
        points_curve(1) = Kratos::make_shared<Point>(3.1513351032106063, 1.5872612941918325, 0);

        Vector knot_vector_curve = ZeroVector(2);
        knot_vector_curve[0] = 0.0;
        knot_vector_curve[1] = 5.0;

        int p_curve = 1;

        // Create the 2D embedded curve
        auto curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(
            points_curve, p_curve, knot_vector_curve);

        // Assign the points belonging to the surface
        PointerVector<Point> points_surface(4);
        points_surface(0) = Kratos::make_shared<Point>(-6, 5, 0);
        points_surface(1) = Kratos::make_shared<Point>(0, -2, 0);
        points_surface(2) = Kratos::make_shared<Point>(39, 30, 0);
        points_surface(3) = Kratos::make_shared<Point>(9, 1, 0);

        // Assign the surface's knot vectors
        Vector knot_vector_u_surface = ZeroVector(2);
        knot_vector_u_surface[0] = 0.0;
        knot_vector_u_surface[1] = 10.0;

        Vector knot_vector_v_surface = ZeroVector(2);
        knot_vector_v_surface[0] = 0.0;
        knot_vector_v_surface[1] = 13.0;

        // Polynomial degrees
        int p_surface = 1;
        int q_surface = 1;

        // Create a 3D surface
        auto surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<Point>>>(points_surface, p_surface,
            q_surface, knot_vector_u_surface, knot_vector_v_surface);

        // Create and return a curve on surface geometry
        return NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>>(surface, curve);
    }

    NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> GenerateReferenceNurbsCOS3dforKnotIntersections(bool uniform_knot_spans = true)
    {
        // Assign the points belonging to the curve
        PointerVector<Point> points_curve;
        if (uniform_knot_spans) {
            points_curve.push_back(Kratos::make_shared<Point>(6, 7, 0));
            points_curve.push_back(Kratos::make_shared<Point>(6, 1, 0));
            points_curve.push_back(Kratos::make_shared<Point>(14, 9, 0));
            points_curve.push_back(Kratos::make_shared<Point>(14, 3, 0));
        }
        else {
            points_curve.push_back(Kratos::make_shared<Point>(0.3, 0.7, 0));
            points_curve.push_back(Kratos::make_shared<Point>(0.3, 0.1, 0));
            points_curve.push_back(Kratos::make_shared<Point>(0.7, 0.9, 0));
            points_curve.push_back(Kratos::make_shared<Point>(0.7, 0.3, 0));
        }

        // Assign the curve's knot vector
        Vector knot_vector_curve = ZeroVector(6);
        knot_vector_curve[0] = 0.0;
        knot_vector_curve[1] = 0.0;
        knot_vector_curve[2] = 0.0;
        if (uniform_knot_spans) {
            knot_vector_curve[3] = 23.313708498984759;
            knot_vector_curve[4] = 23.313708498984759;
            knot_vector_curve[5] = 23.313708498984759;
        }
        else {
            knot_vector_curve[3] = 2.0;
            knot_vector_curve[4] = 2.0;
            knot_vector_curve[5] = 2.0;
        }

        // Polynomial degree of the curve
        int p_curve = 3;

        // Create the 2D embedded curve
        auto curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points_curve, p_curve, knot_vector_curve);

        // Assign the points belonging to the surface
        PointerVector<Point> points_surface;
        points_surface.push_back(Point::Pointer(new Point(0, 0, 0)));
        points_surface.push_back(Point::Pointer(new Point(3.33333333333333, 0, 0)));
        points_surface.push_back(Point::Pointer(new Point(10, 0, 0)));
        points_surface.push_back(Point::Pointer(new Point(16.6666666666666, 0, 0)));
        points_surface.push_back(Point::Pointer(new Point(20, 0, 0)));

        points_surface.push_back(Point::Pointer(new Point(0, 5, 0)));
        points_surface.push_back(Point::Pointer(new Point(3.33333333333333, 5, 0)));
        points_surface.push_back(Point::Pointer(new Point(10, 5.0, 0)));
        points_surface.push_back(Point::Pointer(new Point(16.6666666666666, 5.0, 0)));
        points_surface.push_back(Point::Pointer(new Point(20, 5, 0)));

        points_surface.push_back(Point::Pointer(new Point(0, 10, 0)));
        points_surface.push_back(Point::Pointer(new Point(3.3333333333333326, 10, 0)));
        points_surface.push_back(Point::Pointer(new Point(10, 10.0, 0)));
        points_surface.push_back(Point::Pointer(new Point(16.6666666666666, 10.0, 0)));
        points_surface.push_back(Point::Pointer(new Point(20, 10, 0)));

        // Assign the surface's knot vectors
        Vector knot_vector_u_surface = ZeroVector(7);
        knot_vector_u_surface[0] = 0.0;
        knot_vector_u_surface[1] = 0.0;
        knot_vector_u_surface[2] = 0.0;
        knot_vector_u_surface[3] = 10.0;
        knot_vector_u_surface[4] = 20.0;
        knot_vector_u_surface[5] = 20.0;
        knot_vector_u_surface[6] = 20.0;

        Vector knot_vector_v_surface = ZeroVector(3);
        knot_vector_v_surface[0] = 0.0;
        knot_vector_v_surface[1] = 5.0;
        knot_vector_v_surface[2] = 10.0;


        if (!uniform_knot_spans) {
            knot_vector_u_surface *= 0.05;
            knot_vector_v_surface *= 0.1;
        }

        // Polynomial degrees
        int p_surface = 3;
        int q_surface = 1;

        // Create a 3D surface
        auto surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<Point>>>(points_surface, p_surface,
            q_surface, knot_vector_u_surface, knot_vector_v_surface);

        // Create and return a curve on surface geometry
        return NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>>(surface, curve);
    }

    NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> GenerateReferenceNurbs2dforKnotIntersections(std::vector<Vector>& brep_coordinates)
    {
        Vector knot_vector_u = ZeroVector(4); 
        Vector knot_vector_v = ZeroVector(4); 
        const SizeType p = 1;
        const SizeType q = 1;
        
        // lower
        PointerVector<Point> points_surface;

        knot_vector_u[0] = 0.0; knot_vector_u[1] = 1.0; knot_vector_u[2] = 2.0; knot_vector_u[3] = 3.0;  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = 1.0; knot_vector_v[2] = 2.0; knot_vector_v[3] = 3.0; 

        points_surface.push_back(Point::Pointer(new Point(0.0, 0.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(1.0, 0.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(2.0, 0.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(3.0, 0.0, 0.0)));
        //
        points_surface.push_back(Point::Pointer(new Point(0.0, 0.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(1.0, 1.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(2.0, 1.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(3.0, 1.0, 0.0)));
        //
        points_surface.push_back(Point::Pointer(new Point(0.0, 2.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(1.0, 2.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(2.0, 2.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(3.0, 2.0, 0.0)));
        //
        points_surface.push_back(Point::Pointer(new Point(0.0, 3.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(1.0, 3.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(2.0, 3.0, 0.0)));
        points_surface.push_back(Point::Pointer(new Point(3.0, 3.0, 0.0)));

        auto p_surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<Point>>>(points_surface, p, q, knot_vector_u, knot_vector_v);

        Vector knot_vector_curve = ZeroVector(2);
        knot_vector_curve[0] = 0.0;
        knot_vector_curve[1] = 3.0;

        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType segment;
        segment.push_back(Point::Pointer(new Point(brep_coordinates[0][0], brep_coordinates[0][1])));
        segment.push_back(Point::Pointer(new Point(brep_coordinates[1][0], brep_coordinates[1][1])));

        auto p_curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(segment, p, knot_vector_curve);

        // Create and return a curve on surface geometry
        return NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>>(p_surface, p_curve);
    }

    ///// Tests
    // Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
    KRATOS_TEST_CASE_IN_SUITE(BSplineCurveOnSurfaceBSpline, KratosCoreNurbsGeometriesFastSuite)
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

        KRATOS_EXPECT_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
        KRATOS_EXPECT_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
        KRATOS_EXPECT_VECTOR_NEAR(derivatives[2], gradient2, TOLERANCE);

        // Check kratos geometry family
        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);

    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceNurbs, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateReferenceNurbsCurveOnNurbsSurface3d();

        // Check kratos geometry family
        {
            const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
            const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
            KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
            KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);
        }

        // Evaluate the global coordinates of the curve at u = 4.0
        {
            array_1d<double, 3> parameter = ZeroVector(3);
            parameter[0] = 4.0;

            // Compute the global coordinates of the curve on surface
            std::vector<array_1d<double, 3>> global_coordinates;
            curve_on_surface.GlobalSpaceDerivatives(global_coordinates, parameter, 1);

            // Compare the position vectors
            std::vector<double> positionVct = {0.054174511426802, -0.034958783455029, 0.038314563432352};

            KRATOS_EXPECT_VECTOR_NEAR(global_coordinates[0], positionVct, TOLERANCE);
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

            KRATOS_EXPECT_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
            KRATOS_EXPECT_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
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

            KRATOS_EXPECT_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
            KRATOS_EXPECT_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(CurveOnSurfaceLineLength, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateLineNurbsCurveOnSurface3d();
        KRATOS_EXPECT_NEAR(curve_on_surface.Length(), 5, 1e-1);
    }

    // test intersection with background surface
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceIntersectionSpans, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateReferenceNurbsCOS3dforKnotIntersections();

        std::vector<double> spans;

        curve_on_surface.SpansLocalSpace(spans);

        // Test size
        KRATOS_EXPECT_EQ(spans.size(), 5);

        // Compare each value
        KRATOS_EXPECT_NEAR(spans[0], 0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[1], 4.02565, 1e-4);
        KRATOS_EXPECT_NEAR(spans[2], 11.6569, 1e-4);
        KRATOS_EXPECT_NEAR(spans[3], 19.2881, 1e-4);
        KRATOS_EXPECT_NEAR(spans[4], 23.3137, 1e-4);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);
    }

    // test intersection with background surface with not coinciding knot vectors
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceNonUniformKnotVectorsIntersectionSpans, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateReferenceNurbsCOS3dforKnotIntersections(false);

        std::vector<double> spans;

        curve_on_surface.SpansLocalSpace(spans);

        // Test size
        KRATOS_EXPECT_EQ(spans.size(), 5);


        const double scaling_factor = 23.313708498984759 / 2;
        // Compare each value
        KRATOS_EXPECT_NEAR(spans[0], 0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[1], 4.02565 / scaling_factor, 1e-4);
        KRATOS_EXPECT_NEAR(spans[2], 11.6569 / scaling_factor, 1e-4);
        KRATOS_EXPECT_NEAR(spans[3], 19.2881 / scaling_factor, 1e-4);
        KRATOS_EXPECT_NEAR(spans[4], 23.3137 / scaling_factor, 1e-4);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);
    }

    // test integration of curve on surface
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceCreateIntegrationPoints, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateReferenceNurbsCOS3dforKnotIntersections();

        // Check general information, input to ouput
        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = curve_on_surface.GetDefaultIntegrationInfo();
        curve_on_surface.CreateIntegrationPoints(integration_points, integration_info);

        KRATOS_EXPECT_EQ(integration_points.size(), 20);
        double length = 0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            length += integration_points[i].Weight();
        }
        KRATOS_EXPECT_NEAR(length, 23.313708498984759, TOLERANCE);
    }

    // test quadrature points of curve on surface
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceCreateQuadraturePoints, KratosCoreNurbsGeometriesFastSuite)
    {
        // Nurbs curve on a Nurbs surface
        auto curve_on_surface = GenerateReferenceNurbsCOS3dforKnotIntersections();

        // Check general information, input to ouput
        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = curve_on_surface.GetDefaultIntegrationInfo();
        curve_on_surface.CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Point>::GeometriesArrayType quadrature_points;
        curve_on_surface.CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        KRATOS_EXPECT_EQ(quadrature_points.size(), 20);
        double length_parameter_space = 0;
        double length_global_space = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                length_parameter_space += quadrature_points[i].IntegrationPoints()[j].Weight();
                length_global_space += quadrature_points[i].IntegrationPoints()[j].Weight()
                    * quadrature_points[i].DeterminantOfJacobian(j, quadrature_points[i].GetDefaultIntegrationMethod());
            }
        }
        KRATOS_EXPECT_NEAR(length_parameter_space, 23.313708498984759, TOLERANCE);
        KRATOS_EXPECT_NEAR(length_global_space, curve_on_surface.Length(), TOLERANCE);

        array_1d<double, 3> global_coords;
        array_1d<double, 3> local_coords;
        local_coords[0] = integration_points[2][0];
        local_coords[1] = integration_points[2][1];
        curve_on_surface.GlobalCoordinates(global_coords, local_coords);

        KRATOS_EXPECT_VECTOR_NEAR(quadrature_points[2].Center(), global_coords, TOLERANCE);

        // Check local tangent access
        array_1d<double, 3> local_tangent;
        std::vector<double> local_reference = {0.162409, -0.487862, 0.0};
        quadrature_points[2].Calculate(LOCAL_TANGENT, local_tangent);
        KRATOS_EXPECT_VECTOR_NEAR(local_tangent, local_reference, TOLERANCE);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry;

        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            KRATOS_EXPECT_EQ(quadrature_points[i].GetGeometryFamily(), geometry_family);
            KRATOS_EXPECT_EQ(quadrature_points[i].GetGeometryType(), geometry_type);
        }
    }

    // test intersection with background surface in the SBM scenario
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceIntersectionSpansExternalBodyFitted, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        std::vector<Vector> brep_coordinates(2);
        brep_coordinates[0].resize(2); 
        brep_coordinates[1].resize(2); 
        
        brep_coordinates[0][0] = 0.0;
        brep_coordinates[0][1] = 0.0;
        brep_coordinates[1][0] = 0.0;
        brep_coordinates[1][1] = 3.0;
        
        NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> curve_on_surface = GenerateReferenceNurbs2dforKnotIntersections(brep_coordinates);

        auto p_surface = curve_on_surface.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        std::vector<double> spans;

        curve_on_surface.SpansLocalSpaceSBM(spans);

        // Test size
        KRATOS_EXPECT_EQ(spans.size(), 4);

        // Compare each value
        KRATOS_EXPECT_NEAR(spans[0], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[1], 1.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[2], 2.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[3], 3.0, TOLERANCE);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);
    }

    // test intersection with background surface in the SBM scenario
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceIntersectionSpansExternalSBM, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        std::vector<Vector> brep_coordinates(2);
        brep_coordinates[0].resize(2); 
        brep_coordinates[1].resize(2); 
        
        brep_coordinates[0][0] = 1.0;
        brep_coordinates[0][1] = 0.0;
        brep_coordinates[1][0] = 2.0;
        brep_coordinates[1][1] = 0.0;
        
        NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> curve_on_surface = GenerateReferenceNurbs2dforKnotIntersections(brep_coordinates);

        auto p_surface = curve_on_surface.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        std::vector<double> spans;

        curve_on_surface.SpansLocalSpaceSBM(spans);

        // Test size
        KRATOS_EXPECT_EQ(spans.size(), 2);

        // Compare each value
        KRATOS_EXPECT_NEAR(spans[0], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[1], 3.0, TOLERANCE);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);
    }

    // test intersection with background surface in the SBM scenario
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceIntersectionSpansInternalSBM, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        std::vector<Vector> brep_coordinates(2);
        brep_coordinates[0].resize(2); 
        brep_coordinates[1].resize(2); 
        
        brep_coordinates[0][0] = 1.0;
        brep_coordinates[0][1] = 1.0;
        brep_coordinates[1][0] = 2.0;
        brep_coordinates[1][1] = 1.0;
        
        NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> curve_on_surface = GenerateReferenceNurbs2dforKnotIntersections(brep_coordinates);
        auto p_surface = curve_on_surface.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        std::vector<double> spans;

        curve_on_surface.SpansLocalSpaceSBM(spans);

        // Test size
        KRATOS_EXPECT_EQ(spans.size(), 2);

        // Compare each value
        KRATOS_EXPECT_NEAR(spans[0], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(spans[1], 3.0, TOLERANCE);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Nurbs;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface;
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(curve_on_surface.GetGeometryType(), geometry_type);
    }

    // test quadrature points of curve on surface
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceCreateQuadraturePointsSBMInternal, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        std::vector<Vector> brep_coordinates(2);
        brep_coordinates[0].resize(2); 
        brep_coordinates[1].resize(2); 
        
        brep_coordinates[0][0] = 1.0;
        brep_coordinates[0][1] = 1.0;
        brep_coordinates[1][0] = 2.0;
        brep_coordinates[1][1] = 1.0;
        
        NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> curve_on_surface = GenerateReferenceNurbs2dforKnotIntersections(brep_coordinates);
        auto p_surface = curve_on_surface.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        // Check general information, input to ouput
        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = curve_on_surface.GetDefaultIntegrationInfo();
        curve_on_surface.CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Point>::GeometriesArrayType quadrature_points;
        curve_on_surface.CreateQuadraturePointGeometriesSBM(quadrature_points, 3, integration_points, integration_info);

        KRATOS_EXPECT_EQ(quadrature_points.size(), 3);

        array_1d<double, 3> global_coords;
        array_1d<double, 3> local_coords;
        local_coords[0] = integration_points[2][0];
        local_coords[1] = integration_points[2][1];
        curve_on_surface.GlobalCoordinates(global_coords, local_coords);

        KRATOS_EXPECT_VECTOR_NEAR(quadrature_points[2].Center(), global_coords, TOLERANCE);

        std::vector<Vector> expected_cps(4);
        
        expected_cps[0].resize(2); expected_cps[1].resize(2); expected_cps[2].resize(2); expected_cps[3].resize(2);

        expected_cps[0][0] = 1.0; expected_cps[0][1] = 1.0; 
        expected_cps[1][0] = 2.0; expected_cps[1][1] = 1.0; 
        expected_cps[2][0] = 1.0; expected_cps[2][1] = 2.0; 
        expected_cps[3][0] = 2.0; expected_cps[3][1] = 2.0; 

        // Check the correct SpanU++
        KRATOS_EXPECT_EQ(quadrature_points[0].size(), 4);
        for (IndexType i = 0; i < quadrature_points[0].size(); i++){
            
            Vector cp_coordinates(2);            
            cp_coordinates[0] = quadrature_points[0][i].X();
            cp_coordinates[1] = quadrature_points[0][i].Y();

            KRATOS_EXPECT_NEAR(norm_2(cp_coordinates-expected_cps[i]), 0.0, TOLERANCE);
        
        }        
        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry;

        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            KRATOS_EXPECT_EQ(quadrature_points[i].GetGeometryFamily(), geometry_family);
            KRATOS_EXPECT_EQ(quadrature_points[i].GetGeometryType(), geometry_type);
        }
    }

    // test quadrature points of curve on surface
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveOnSurfaceCreateQuadraturePointsSBMExternal, KratosCoreNurbsGeometriesFastSuite)
    {
        // Create a Nurbs curve on a Nurbs surface
        std::vector<Vector> brep_coordinates(2);
        brep_coordinates[0].resize(2); 
        brep_coordinates[1].resize(2); 
        
        brep_coordinates[0][0] = 1.0;
        brep_coordinates[0][1] = 0.0;
        brep_coordinates[1][0] = 2.0;
        brep_coordinates[1][1] = 0.0;
        
        NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<Point>> curve_on_surface = GenerateReferenceNurbs2dforKnotIntersections(brep_coordinates);
        auto p_surface = curve_on_surface.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        // Check general information, input to ouput
        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = curve_on_surface.GetDefaultIntegrationInfo();
        curve_on_surface.CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Point>::GeometriesArrayType quadrature_points;
        curve_on_surface.CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        KRATOS_EXPECT_EQ(quadrature_points.size(), 3);

        array_1d<double, 3> global_coords;
        array_1d<double, 3> local_coords;
        local_coords[0] = integration_points[2][0];
        local_coords[1] = integration_points[2][1];
        curve_on_surface.GlobalCoordinates(global_coords, local_coords);

        KRATOS_EXPECT_VECTOR_NEAR(quadrature_points[2].Center(), global_coords, TOLERANCE);

        std::vector<Vector> expected_cps(4);
        
        expected_cps[0].resize(2); expected_cps[1].resize(2); expected_cps[2].resize(2); expected_cps[3].resize(2);

        expected_cps[0][0] = 1.0; expected_cps[0][1] = 0.0; 
        expected_cps[1][0] = 2.0; expected_cps[1][1] = 0.0; 
        expected_cps[2][0] = 1.0; expected_cps[2][1] = 1.0; 
        expected_cps[3][0] = 2.0; expected_cps[3][1] = 1.0; 

        // Check the correct SpanU++
        KRATOS_EXPECT_EQ(quadrature_points[0].size(), 4);
        for (IndexType i = 0; i < quadrature_points[0].size(); i++){
            
            Vector cp_coordinates(2);            
            cp_coordinates[0] = quadrature_points[0][i].X();
            cp_coordinates[1] = quadrature_points[0][i].Y();
            
            KRATOS_EXPECT_NEAR(norm_2(cp_coordinates-expected_cps[i]), 0.0, TOLERANCE);
        
        }        
        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry;

        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            KRATOS_EXPECT_EQ(quadrature_points[i].GetGeometryFamily(), geometry_family);
            KRATOS_EXPECT_EQ(quadrature_points[i].GetGeometryType(), geometry_type);
        }
    }


    
} // namespace Testing.
} // namespace Kratos.
