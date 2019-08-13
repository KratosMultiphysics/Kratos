//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_surface_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    /// Factory functions
    NurbsSurfaceGeometry<3, Point> GenerateReferencePointSurface()
    {
        NurbsSurfaceGeometry<3, Point>::PointsArrayType points;

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

        return NurbsSurfaceGeometry<3, Point>(
                points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    NurbsSurfaceGeometry<3, Point> GenerateReferencePieceOfCylinderNurbsSurface()
    {
        NurbsSurfaceGeometry<3, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(0, 10, 0)));
        points.push_back(Point::Pointer(new Point(6.6817863791929888, 10, 0)));
        points.push_back(Point::Pointer(new Point(9.2387953251128678, 3.8268343236508979, 0)));
        points.push_back(Point::Pointer(new Point(11.795804271032745, -2.3463313526982033, 0)));
        points.push_back(Point::Pointer(new Point(7.0710678118654755, -7.0710678118654755, 0)));
        points.push_back(Point::Pointer(new Point(0, 10, 10)));
        points.push_back(Point::Pointer(new Point(6.6817863791929888, 10, 10)));
        points.push_back(Point::Pointer(new Point(9.2387953251128678, 3.8268343236508979, 10)));
        points.push_back(Point::Pointer(new Point(11.795804271032745, -2.3463313526982033, 10)));
        points.push_back(Point::Pointer(new Point(7.0710678118654755, -7.0710678118654755)));

        Vector knot_vector_u = ZeroVector(6);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 11.780972450961723;
        knot_vector_u[3] = 11.780972450961723;
        knot_vector_u[4] = 23.561944901923447;
        knot_vector_u[5] = 23.561944901923447;

        Vector knot_vector_v = ZeroVector(2);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 10.0;

        int p = 2;
        int q = 1;

        Vector weights = ZeroVector(10);
        weights[0] = 1.0;
        weights[1] = 0.83146961230254524;
        weights[2] = 1.0;
        weights[3] = 0.83146961230254524;
        weights[4] = 1.0;
        weights[5] = 1.0;
        weights[6] = 0.83146961230254524;
        weights[7] = 1.0;
        weights[8] = 0.83146961230254524;
        weights[9] = 1.0;

        return NurbsSurfaceGeometry<3, Point>(
            points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    NurbsSurfaceGeometry<3, NodeType> GenerateReferenceNodeSurface() {
        Geometry<NodeType>::PointsArrayType points;

        points.push_back(NodeType::Pointer(new NodeType(1, 0, 5, 0)));
        points.push_back(NodeType::Pointer(new NodeType(2, 5, 5, 0)));
        points.push_back(NodeType::Pointer(new NodeType(3, 10, 5, -4)));
        points.push_back(NodeType::Pointer(new NodeType(4, 0, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(5, 5, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(6, 10, 0, -4)));

        Vector knot_u = ZeroVector(4);
        knot_u[0] = 0.0;
        knot_u[1] = 0.0;
        knot_u[2] = 10.0;
        knot_u[3] = 10.0;
        Vector knot_v = ZeroVector(2); 
        knot_v[0] = 0.0;
        knot_v[1] = 5.0;

        int p = 2;
        int q = 1;

        return NurbsSurfaceGeometry<3, NodeType>(
            points, p, q, knot_u, knot_v);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfacePoint, KratosCoreGeometriesFastSuite) {
        auto surface = GenerateReferencePointSurface();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), true);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 2);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 1);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 5);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 3);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 4);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 3);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 12);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.0;
        parameter[1] = 0.0;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCylinderSurface, KratosCoreGeometriesFastSuite) {
        auto surface = GenerateReferencePieceOfCylinderNurbsSurface();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), true);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 2);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 1);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 6);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 2);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 5);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 2);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 10);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 10.0;
        parameter[1] = 3.5;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
        double length = sqrt(result[0] * result[0] + result[1] * result[1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], parameter[1], TOLERANCE);

        auto derivatives = surface.GlobalDerivatives(parameter, 3);
        array_1d<double, 3> cross(0.0);
        array_1d<double, 3> colinear_vector(0.0);
        derivatives[0][2] = 0.0;
        MathUtils<double>::CrossProduct(cross, derivatives[1], derivatives[2]);
        MathUtils<double>::CrossProduct(colinear_vector, cross, derivatives[0]);
        KRATOS_CHECK_NEAR(norm_2(colinear_vector), 0, TOLERANCE);

        parameter[0] = 6.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        length = sqrt(result[0] * result[0] + result[1] * result[1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], parameter[1], TOLERANCE);

        parameter[0] = 0.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        length = sqrt(result[0] * result[0] + result[1] * result[1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], parameter[1], TOLERANCE);

        parameter[0] = 0;
        parameter[1] = 1.0;
        auto derivatives_2 = surface.GlobalDerivatives(parameter, 3);
        length = sqrt(derivatives_2[0][0] * derivatives_2[0][0]
            + derivatives_2[0][1] * derivatives_2[0][1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[0][2], parameter[1], TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][0]/norm_2(derivatives_2[1]), 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][2], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][2], 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceNode, KratosCoreGeometriesFastSuite) {
        auto surface = GenerateReferenceNodeSurface();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), false);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 2);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 1);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 4);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 2);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 3);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 2);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 6);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 10.0;
        parameter[1] = 3.5;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 1.5, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], -4.0, TOLERANCE);

        parameter[0] = 6.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 6.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 4.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], - 1.44, TOLERANCE);

        parameter[0] = 0.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 4.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 0.0, TOLERANCE);
    }
} // namespace Testing.
} // namespace Kratos.
