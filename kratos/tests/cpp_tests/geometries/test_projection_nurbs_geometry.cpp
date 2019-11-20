//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Apostolatos
//                   Tobias Teschemacher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;

    /// Factory functions
    NurbsCurveGeometry<2, PointerVector<Point>> GenerateReferenceCurveForProjection2d()
    {
        PointerVector<Point> points;

        points.push_back(Point::Pointer(new Point(-9, -2, 0)));
        points.push_back(Point::Pointer(new Point(-5, -3, 0)));
        points.push_back(Point::Pointer(new Point(-3, -5, 0)));
        points.push_back(Point::Pointer(new Point(-2, -4, 0)));
        points.push_back(Point::Pointer(new Point(2, 4, 0)));
        points.push_back(Point::Pointer(new Point(5, 3, 0)));
        points.push_back(Point::Pointer(new Point(9, -2, 0)));

        Vector knot_vector = ZeroVector(9);
        knot_vector[0] = -1.0;
        knot_vector[1] = -1.0;
        knot_vector[2] = -1.0;
        knot_vector[3] = -0.5;
        knot_vector[4] = 0.0;
        knot_vector[5] = 0.75;
        knot_vector[6] = 1.0;
        knot_vector[7] = 1.0;
        knot_vector[8] = 1.0;

        int p = 3;

        Vector weights = ZeroVector(7);
        weights[0] = 1.0;
        weights[1] = 2.0;
        weights[2] = 3.4;
        weights[3] = 1.0;
        weights[4] = 5.7;
        weights[5] = 4.3;
        weights[6] = 1.0;

        auto curve = NurbsCurveGeometry<2, PointerVector<Point>>(points, p, knot_vector, weights);

        return curve;
    }

    NurbsSurfaceGeometry<3, PointerVector<Point>> GenerateReferenceUnrefinedQuarterSphereGeometry()
    {
        PointerVector<Point> points;

        points.push_back(Point::Pointer(new Point(0.0, -0.075, 0.0)));
        points.push_back(Point::Pointer(new Point(0.075, -0.075, 0.0)));
        points.push_back(Point::Pointer(new Point(0.075, 0.0, 0.0)));

        points.push_back(Point::Pointer(new Point(0.0, -0.075, 0.075)));
        points.push_back(Point::Pointer(new Point(0.075, -0.075, 0.075)));
        points.push_back(Point::Pointer(new Point(0.075, 0.0, 0.0)));

        points.push_back(Point::Pointer(new Point(0.0, 0.0, 0.075)));
        points.push_back(Point::Pointer(new Point(0.075, 0.0, 0.075)));
        points.push_back(Point::Pointer(new Point(0.075, 0.0, 0.0)));

        Vector knot_vector_u = ZeroVector(4);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 1.0;
        knot_vector_u[3] = 1.0;

        Vector knot_vector_v = ZeroVector(4);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 1.0;
        knot_vector_v[3] = 1.0;

        int p = 2;
        int q = 2;

        Vector weights = ZeroVector(9);
        weights[0] = 1.0;
        weights[1] = 7.071067811865476e-01;
        weights[2] = 1.0;

        weights[3] = 7.071067811865476e-01;
        weights[4] = 5.000000000000001e-01;
        weights[5] = 7.071067811865476e-01;

        weights[6] = 1.0;
        weights[7] = 7.071067811865476e-01;
        weights[8] = 1.0;

        return NurbsSurfaceGeometry<3, PointerVector<Point>>(
            points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    /// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveGeometryProjection2d, KratosCoreNurbsGeometriesFastSuite) {
        auto curve = GenerateReferenceCurveForProjection2d();

        array_1d<double, 3> point;
        point[0] = -2.0;
        point[1] = 1.0;
        point[2] = 0.0;

        // Initialize projected point
        array_1d<double, 3> projected_point = ZeroVector(3);

        // Try projection with initial guess at u = 0.0
        array_1d<double, 3> parameter = ZeroVector(3);

        bool isConverged1 = ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(parameter, point, projected_point, curve);

        KRATOS_CHECK_EQUAL(isConverged1, true);

        KRATOS_CHECK_NEAR(parameter[0], 0.099395977882318, TOLERANCE);

        std::vector<double> projected_point1 = {-0.129744540301921, -0.044240249340891, 0.0};

        KRATOS_CHECK_VECTOR_NEAR(projected_point, projected_point1, TOLERANCE);

        // Try projection with initial guess at u = -1.0
        parameter[0] = curve.DomainInterval().MinParameter();

        bool isConverged2 = ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(parameter, point, projected_point, curve);

        KRATOS_CHECK_EQUAL(isConverged2, true);

        KRATOS_CHECK_NEAR(parameter[0], -0.788227217287371, TOLERANCE);

        std::vector<double> projected_point2 = {-4.694701201131293, -3.571229085898834, 0.0};

        KRATOS_CHECK_VECTOR_NEAR(projected_point, projected_point2, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceGeometryProjection3d, KratosCoreNurbsGeometriesFastSuite) {
        auto surface = GenerateReferenceUnrefinedQuarterSphereGeometry();

        // Initialize projected point and projected parameter
        array_1d<double, 3> projected_point = ZeroVector(3);
        array_1d<double, 3> parameter = ZeroVector(3);

        // Initialize the point to be projected onto the Nurbs surface
        array_1d<double, 3> point;

        // Initialize expected projected points
        std::vector<double> projected_point_expected1 = {0.043301282246217, -0.043301264160723, 0.043301264160723};
        std::vector<double> projected_point_expected2 = {0.074813167561162, -0.003740986465195, 0.003740986465195};

        // Try projection of point (0.06, -0.06, 0.06) far from the singularity
        point[0] = 0.06;
        point[1] = -0.06;
        point[2] = 0.06;

        // Try projection with initial guess at u = 0.5, v = 0.5
        parameter[0] = 0.5;
        parameter[1] = 0.5;

        // Try projection
        bool isConverged1 = ProjectionNurbsGeometryUtilities::NewtonRaphsonSurface(parameter, point, projected_point, surface);
        
        // Check the solution
        KRATOS_CHECK_EQUAL(isConverged1, true);
        KRATOS_CHECK_NEAR(parameter[0], 0.397197796315686, TOLERANCE);
        KRATOS_CHECK_NEAR(parameter[1], 0.5, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(projected_point, projected_point_expected1, TOLERANCE);

        // Try projection of point (0.1, -0.005, 0.005) close to the singularity
        point[0] = 0.1;
        point[1] = -0.005;
        point[2] = 0.005;

        // Try projection with initial guess at u = 0.5, v = 0.5
        parameter[0] = 0.5;
        parameter[1] = 0.5;

        // Try projection
        bool isConverged2 = ProjectionNurbsGeometryUtilities::NewtonRaphsonSurface(parameter, point, projected_point, surface, 100, 1e-12);

        // Check the solution
        KRATOS_CHECK_EQUAL(isConverged2, true);
        KRATOS_CHECK_NEAR(parameter[0], 0.950782, TOLERANCE);
        KRATOS_CHECK_NEAR(parameter[1], 0.500000, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(projected_point, projected_point_expected2, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.