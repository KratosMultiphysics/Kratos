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

typedef Node NodeType;

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

    NurbsSurfaceGeometry<3, PointerVector<Node>> GenerateReferenceC0HemisphereGeometry(ModelPart& rModelPart)
    {
        PointerVector<Node> points;

        points.push_back(rModelPart.CreateNewNode(1, 0.0, 0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(2, 0.0, 0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(3, 0.0, 0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(4, 0.0, 0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(5, 0.0, 0.25, 0.0));

        points.push_back(rModelPart.CreateNewNode(6, -0.25, 0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(7, -0.25, 0.25, 0.25));
        points.push_back(rModelPart.CreateNewNode(8, -1.5308084989341915e-17, 0.25, 0.25));
        points.push_back(rModelPart.CreateNewNode(9, 0.24999999999999994, 0.25, 0.25));
        points.push_back(rModelPart.CreateNewNode(10, 0.25, 0.25, 3.061616997868383e-17));

        points.push_back(rModelPart.CreateNewNode(11, -0.25, 1.5308084989341915e-17, 0.0));
        points.push_back(rModelPart.CreateNewNode(12, -0.25, 1.5308084989341915e-17, 0.25));
        points.push_back(rModelPart.CreateNewNode(13, -1.5308084989341915e-17, 1.5308084989341915e-17, 0.25));
        points.push_back(rModelPart.CreateNewNode(14, 0.24999999999999994, 1.5308084989341915e-17, 0.25));
        points.push_back(rModelPart.CreateNewNode(15, 0.25, 1.5308084989341915e-17, 3.061616997868383e-17));

        points.push_back(rModelPart.CreateNewNode(16, -0.25, -0.24999999999999994, 0.0));
        points.push_back(rModelPart.CreateNewNode(17, -0.25, -0.24999999999999994, 0.25));
        points.push_back(rModelPart.CreateNewNode(18, -1.5308084989341915e-17, -0.24999999999999994, 0.25));
        points.push_back(rModelPart.CreateNewNode(19, 0.24999999999999994, -0.24999999999999994, 0.25));
        points.push_back(rModelPart.CreateNewNode(20, 0.25, -0.24999999999999994, 3.061616997868383e-17));

        points.push_back(rModelPart.CreateNewNode(21, 0.0, -0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(22, 0.0, -0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(23, 0.0, -0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(24, 0.0, -0.25, 0.0));
        points.push_back(rModelPart.CreateNewNode(25, 0.0, -0.25, 0.0));

        Vector knot_vector_u = ZeroVector(8);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 0.39269908169872414;
        knot_vector_u[4] = 0.39269908169872414;
        knot_vector_u[5] = 0.78539816339744828;
        knot_vector_u[6] = 0.78539816339744828;
        knot_vector_u[7] = 0.78539816339744828;

        Vector knot_vector_v = knot_vector_u;

        const int p = 2;
        const int q = 2;

        Vector weights = ZeroVector(25);
        weights[0] = 1.0;
        weights[1] = 0.70710678118654757;
        weights[2] = 1.0;
        weights[3] = 0.70710678118654757;
        weights[4] = 1.0;

        weights[5] = 0.70710678118654757;
        weights[6] = 0.50000000000000011;
        weights[7] = 0.70710678118654757;
        weights[8] = 0.50000000000000011;
        weights[9] = 0.70710678118654757;

        weights[10] = 1.0;
        weights[11] = 0.70710678118654757;
        weights[12] = 1.0;
        weights[13] = 0.70710678118654757;
        weights[14] = 1.0;

        weights[15] = 0.70710678118654757;
        weights[16] = 0.50000000000000011;
        weights[17] = 0.70710678118654757;
        weights[18] = 0.50000000000000011;
        weights[19] = 0.70710678118654757;

        weights[20] = 1.0;
        weights[21] = 0.70710678118654757;
        weights[22] = 1.0;
        weights[23] = 0.70710678118654757;
        weights[24] = 1.0;

        return NurbsSurfaceGeometry<3, PointerVector<Node>>(
            points,
            p,
            q,
            knot_vector_u,
            knot_vector_v,
            weights);
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

        KRATOS_EXPECT_EQ(isConverged1, true);

        KRATOS_EXPECT_NEAR(parameter[0], 0.099395977882318, TOLERANCE);

        std::vector<double> projected_point1 = {-0.129744540301921, -0.044240249340891, 0.0};

        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point1, TOLERANCE);

        // Try projection with initial guess at u = -1.0
        parameter[0] = curve.DomainInterval().MinParameter();

        bool isConverged2 = ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(parameter, point, projected_point, curve);

        KRATOS_EXPECT_EQ(isConverged2, true);

        KRATOS_EXPECT_NEAR(parameter[0], -0.788227217287371, TOLERANCE);

        std::vector<double> projected_point2 = {-4.694701201131293, -3.571229085898834, 0.0};

        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point2, TOLERANCE);

        array_1d<double, 3> curve_extreme_point{-9.0, -2.0, 0.0};

        // Initialize projected point
        array_1d<double, 3> projected_point_extreme = ZeroVector(3);

        // Try projection with initial guess at u = 0.0
        array_1d<double, 3> parameter_extreme = ZeroVector(3);

        const bool is_converged_extreme = ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(parameter_extreme, curve_extreme_point, projected_point_extreme, curve);

        KRATOS_EXPECT_EQ(is_converged_extreme, true);

        KRATOS_EXPECT_NEAR(parameter_extreme[0], -1.0, TOLERANCE);

        array_1d<double, 3> projected_point_extreme1{-9.0, -2.0, 0.0};

        KRATOS_EXPECT_VECTOR_NEAR(projected_point_extreme, projected_point_extreme1, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveGeometryLevenbergMarquardtProjection2d, KratosCoreNurbsGeometriesFastSuite)
    {
        auto curve = GenerateReferenceCurveForProjection2d();

        array_1d<double, 3> point;
        point[0] = -2.0;
        point[1] = 1.0;
        point[2] = 0.0;

        array_1d<double, 3> projected_point = ZeroVector(3);

        // Test 1: Initial guess at the center of the parameter domain
        array_1d<double, 3> parameter = ZeroVector(3);

        bool is_converged_1 = ProjectionNurbsGeometryUtilities::LevenbergMarquardtCurve(parameter, point, projected_point, curve, 50, 1e-7);

        KRATOS_EXPECT_EQ(is_converged_1, true);

        KRATOS_EXPECT_NEAR(parameter[0],0.099395977882318, TOLERANCE);

        std::vector<double> projected_point_expected_1 = {
            -0.129744540301921,
            -0.044240249340891,
            0.0
        };

        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point_expected_1, TOLERANCE);
        

        // Test 2: Initial guess at the left boundary
        parameter = ZeroVector(3);
        parameter[0] = curve.DomainInterval().MinParameter();

        bool is_converged_2 = ProjectionNurbsGeometryUtilities::LevenbergMarquardtCurve(parameter, point, projected_point, curve, 50, 1e-7);

        KRATOS_EXPECT_EQ(is_converged_2, true);

        KRATOS_EXPECT_NEAR(parameter[0], -0.788227217287371, TOLERANCE);

        std::vector<double> projected_point_expected_2 = {
            -4.694701201131293,
            -3.571229085898834,
            0.0
        };

        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point_expected_2, TOLERANCE);

        // Test 3: Projection of a point already located on the curve end
        array_1d<double, 3> curve_extreme_point{
            -9.0,
            -2.0,
            0.0
        };

        array_1d<double, 3> projected_point_extreme = ZeroVector(3);
        array_1d<double, 3> parameter_extreme = ZeroVector(3);

        bool is_converged_3 = ProjectionNurbsGeometryUtilities::LevenbergMarquardtCurve(parameter_extreme, curve_extreme_point, projected_point_extreme, curve, 50, 1e-7);

        KRATOS_EXPECT_EQ(is_converged_3, true);

        KRATOS_EXPECT_NEAR(parameter_extreme[0], -1.0, TOLERANCE);

        array_1d<double, 3> projected_point_expected_3{
            -9.0,
            -2.0,
            0.0
        };

        KRATOS_EXPECT_VECTOR_NEAR(projected_point_extreme, projected_point_expected_3, TOLERANCE);
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
        KRATOS_EXPECT_EQ(isConverged1, true);
        KRATOS_EXPECT_NEAR(parameter[0], 0.397197796315686, TOLERANCE);
        KRATOS_EXPECT_NEAR(parameter[1], 0.5, TOLERANCE);
        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point_expected1, TOLERANCE);

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
        KRATOS_EXPECT_EQ(isConverged2, true);
        KRATOS_EXPECT_NEAR(parameter[0], 0.950782, TOLERANCE);
        KRATOS_EXPECT_NEAR(parameter[1], 0.500000, TOLERANCE);
        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point_expected2, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceGeometryLevenbergMarquardtProjection3d, KratosCoreNurbsGeometriesFastSuite)
    {
        auto surface = GenerateReferenceUnrefinedQuarterSphereGeometry();

        array_1d<double, 3> projected_point = ZeroVector(3);
        array_1d<double, 3> parameter = ZeroVector(3);
        array_1d<double, 3> point;

        const double lm_tolerance = 1e-5;

        std::vector<double> projected_point_expected1 = {
            0.043301282246217,
            -0.043301264160723,
            0.043301264160723
        };

        std::vector<double> projected_point_expected2 = {
            0.074813167561162,
            -0.003740986465195,
            0.003740986465195
        };

        // Projection far from the singularity
        point[0] = 0.06;
        point[1] = -0.06;
        point[2] = 0.06;

        parameter[0] = 0.5;
        parameter[1] = 0.5;

        bool is_converged_1 = ProjectionNurbsGeometryUtilities::LevenbergMarquardtSurface(parameter, point, projected_point, surface, 100, 1e-6);

        KRATOS_EXPECT_EQ(is_converged_1, true);
        KRATOS_EXPECT_NEAR(parameter[0], 0.397197796315686, lm_tolerance);
        KRATOS_EXPECT_NEAR(parameter[1], 0.5, lm_tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(projected_point, projected_point_expected1, lm_tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceGeometryC0HemisphereProjectionNewtonVsLevenbergMarquardt, KratosCoreNurbsGeometriesFastSuite)
    {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("TestModelPart");

        auto surface = GenerateReferenceC0HemisphereGeometry(model_part);

        array_1d<double, 3> point;
        point[0] = -0.25029993;
        point[1] =  0.00097306;
        point[2] =  0.10123734;

        array_1d<double, 3> projected_point_nr = ZeroVector(3);
        array_1d<double, 3> projected_point_lm = ZeroVector(3);

        array_1d<double, 3> parameter_nr = ZeroVector(3);
        array_1d<double, 3> parameter_lm = ZeroVector(3);

        parameter_nr[0] = 0.39269908169872414;
        parameter_nr[1] = 0.39269908169872414;

        parameter_lm[0] = 0.39269908169872414;
        parameter_lm[1] = 0.39269908169872414;
        
        // Newton-Raphson fails here
        const bool is_converged_nr = ProjectionNurbsGeometryUtilities::NewtonRaphsonSurface(parameter_nr, point, projected_point_nr, surface, 20, 1e-6);
        // LevenbergMarquardt works here
        const bool is_converged_lm = ProjectionNurbsGeometryUtilities::LevenbergMarquardtSurface(parameter_lm, point, projected_point_lm, surface, 50, 1e-6);

        KRATOS_EXPECT_EQ(is_converged_nr, false);
        KRATOS_EXPECT_EQ(is_converged_lm, true);

        KRATOS_EXPECT_NEAR(parameter_lm[0], 0.1, 1e-5);
        KRATOS_EXPECT_NEAR(parameter_lm[1], 0.39169908169872414, 1e-5);
    }

} // namespace Testing.
} // namespace Kratos.
