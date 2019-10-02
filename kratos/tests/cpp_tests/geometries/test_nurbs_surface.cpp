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
//                   Andreas Apostolatos
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_surface_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    /// Factory functions
    NurbsSurfaceGeometry<3, PointerVector<Point>> GenerateReferencePointSurface()
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
    }

    NurbsSurfaceGeometry<3, PointerVector<Point>> GenerateReferencePieceOfCylinderNurbsSurface()
    {
        PointerVector<Point> points;

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

        return NurbsSurfaceGeometry<3, PointerVector<Point>>(
            points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    NurbsSurfaceGeometry<3, PointerVector<NodeType>> GenerateReferenceNodeSurface() {
        PointerVector<NodeType> points;

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

        return NurbsSurfaceGeometry<3, PointerVector<NodeType>>(
            points, p, q, knot_u, knot_v);
    }

    NurbsSurfaceGeometry<3, PointerVector<Point>> GenerateReferenceQuarterSphereGeometry()
    {
        NurbsSurfaceGeometry<3, PointerVector<Point>>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(000000000000000, -7.500000000000000e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(1.260557702650798e-02, -7.500000000000000e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(3.897300454513554e-02, -6.844963719002463e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(6.844963719002463e-02, -3.897300454513554e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -1.260557702650798e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -7.500000000000000e-02, 1.430576773729096e-02)));
        points.push_back(Point::Pointer(new Point(1.260557702650798e-02, -7.500000000000000e-02, 1.430576773729096e-02)));
        points.push_back(Point::Pointer(new Point(3.897300454513554e-02, -6.844963719002463e-02, 1.305632815123102e-02)));
        points.push_back(Point::Pointer(new Point(6.844963719002464e-02, -3.897300454513555e-02, 7.433850013961258e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -1.260557702650798e-02, 2.404432761810054e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -6.746212024587489e-02, 4.484848098349960e-02)));
        points.push_back(Point::Pointer(new Point(1.260557702650798e-02, -6.746212024587490e-02, 4.484848098349960e-02)));
        points.push_back(Point::Pointer(new Point(3.897300454513555e-02, -6.157010206533272e-02, 4.093149669125690e-02)));
        points.push_back(Point::Pointer(new Point(6.844963719002464e-02, -3.505602025289285e-02, 2.330506737616474e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -1.133865270841226e-02, 7.537879754125100e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -4.484848098349960e-02, 6.746212024587489e-02)));
        points.push_back(Point::Pointer(new Point(1.260557702650798e-02, -4.484848098349960e-02, 6.746212024587490e-02)));
        points.push_back(Point::Pointer(new Point(3.897300454513555e-02, -4.093149669125690e-02, 6.157010206533272e-02)));
        points.push_back(Point::Pointer(new Point(6.844963719002464e-02, -2.330506737616474e-02, 3.505602025289285e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -7.537879754125100e-03, 1.133865270841226e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -1.430576773729096e-02, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(1.260557702650798e-02, -1.430576773729096e-02, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(3.897300454513554e-02, -1.305632815123102e-02, 6.844963719002463e-02)));
        points.push_back(Point::Pointer(new Point(6.844963719002464e-02, -7.433850013961258e-03, 3.897300454513555e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -2.404432761810054e-03, 1.260557702650798e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, 000000000000000, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(1.260557702650798e-02, 000000000000000, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(3.897300454513554e-02, 000000000000000, 6.844963719002463e-02)));
        points.push_back(Point::Pointer(new Point(6.844963719002463e-02, 000000000000000, 3.897300454513554e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 1.260557702650798e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        Vector knot_vector_u = ZeroVector(8);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 0.0;
        knot_vector_u[3] = 3.333333333333333e-01;
        knot_vector_u[4] = 6.666666666666666e-01;
        knot_vector_u[5] = 1.0;
        knot_vector_u[6] = 1.0;
        knot_vector_u[7] = 1.0;

        Vector knot_vector_v = ZeroVector(9);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 0.0;
        knot_vector_v[2] = 0.0;
        knot_vector_v[3] = 0.0;
        knot_vector_v[4] = 0.5;
        knot_vector_v[5] = 1.0;
        knot_vector_v[6] = 1.0;
        knot_vector_v[7] = 1.0;
        knot_vector_v[8] = 1.0;

        int p = 3;
        int q = 4;

        Vector weights = ZeroVector(36);
        weights[0] = 1.0;
        weights[1] = 9.349126180414551e-01;
        weights[2] = 8.481294420967285e-01;
        weights[3] = 8.481294420967285e-01;
        weights[4] = 9.349126180414551e-01;
        weights[5] = 1.0;

        weights[6] = 9.267766952966369e-01;
        weights[7] = 8.664552265395866e-01;
        weights[8] = 7.860266015301863e-01;
        weights[9] = 7.860266015301862e-01;
        weights[10] = 8.664552265395866e-01;
        weights[11] = 9.267766952966369e-01;

        weights[12] = 8.291456223588194e-01;
        weights[13] = 7.751787045370954e-01;
        weights[14] = 7.032228141081300e-01;
        weights[15] = 7.032228141081300e-01;
        weights[16] = 7.751787045370954e-01;
        weights[17] = 8.291456223588194e-01;

        weights[18] = 8.291456223588194e-01;
        weights[19] = 7.751787045370954e-01;
        weights[20] = 7.032228141081300e-01;
        weights[21] = 7.032228141081300e-01;
        weights[22] = 7.751787045370954e-01;
        weights[23] = 8.291456223588194e-01;

        weights[24] = 9.267766952966369e-01;
        weights[25] = 8.664552265395866e-01;
        weights[26] = 7.860266015301863e-01;
        weights[27] = 7.860266015301862e-01;
        weights[28] = 8.664552265395866e-01;
        weights[29] = 9.267766952966369e-01;

        weights[30] = 1.0;
        weights[31] = 9.349126180414551e-01;
        weights[32] = 8.481294420967285e-01;
        weights[33] = 8.481294420967285e-01;
        weights[34] = 9.349126180414551e-01;
        weights[35] = 1.0;

        return NurbsSurfaceGeometry<3, PointerVector<Point>>(
            points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfacePoint, KratosCoreNurbsGeometriesFastSuite) {
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

    KRATOS_TEST_CASE_IN_SUITE(NurbsCylinderSurface, KratosCoreNurbsGeometriesFastSuite) {
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

        std::vector<array_1d<double, 3>> derivatives;
        surface.GlobalSpaceDerivatives(derivatives, parameter, 3);
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
        std::vector<array_1d<double, 3>> derivatives_2;
        surface.GlobalSpaceDerivatives(derivatives_2, parameter, 3);
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

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceNode, KratosCoreNurbsGeometriesFastSuite) {
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

    KRATOS_TEST_CASE_IN_SUITE(NurbsQuarterSphereSurface, KratosCoreNurbsGeometriesFastSuite) {
        auto surface = GenerateReferenceQuarterSphereGeometry();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), true);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 3);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 4);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 8);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 9);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 6);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 6);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 36);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = .25;
        parameter[1] = .75;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 0.027607103217140, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], -0.025668761597520, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 0.064837971359442, TOLERANCE);

        std::vector<array_1d<double, 3>> derivatives;
        surface.GlobalSpaceDerivatives(derivatives, parameter, 6);

        // Compare the position vectors and the gradients up to 6th order
        std::vector<double> positionVct = {0.027607103217140, -0.025668761597520, 0.064837971359442};
        std::vector<double> gradient1 = {0.110787255345493, 0.016144510322666, -0.040780202579556};
        std::vector<double> gradient2 = {-0.000000000000000, 0.103008693927056, 0.040780202579556};
        std::vector<double> gradient3 = {-0.033227653950941, 0.070099959163148, -0.177068890809183};
        std::vector<double> gradient4 = {-0.000000000000000, -0.064787890764092, -0.025648935146511};
        std::vector<double> gradient5 = {-0.000000000000000, 0.030894683915335, -0.177068890809183};
        std::vector<double> gradient6 = {-0.470230208134118, 0.005450671128399, -0.013768114880420};
        std::vector<double> gradient7 = {-0.000000000000003, -0.281311009504776, -0.111368463360849};
        std::vector<double> gradient8 = {-0.000000000000001, -0.019431383220095, 0.111368463360851};
        std::vector<double> gradient9 = {-0.000000000000002, -0.437214546329255, 0.013768114880415};
        std::vector<double> gradient10 = {-0.356492597666529, -0.546390143805925, 1.380153396204004};
        std::vector<double> gradient11 = {-0.000000000000008, -0.021873533393088, -0.008659532403435};
        std::vector<double> gradient12 = {0.000000000000013, -0.084371662130872, 0.483565284894256};
        std::vector<double> gradient13 = {-0.000000000000012, 0.274988519785633, -0.008659532403469};
        std::vector<double> gradient14 = {-0.000000000000040, 0.331462646725662, 1.380153396204038};
        std::vector<double> gradient15 = {5.602352380836633, -0.970638457490299, 2.451782812662428};
        std::vector<double> gradient16 = {0.000000000000008, 2.192662660184668, 0.868055152065908};
        std::vector<double> gradient17 = {0.000000000000021, -0.006560377328662, 0.037599955385659};
        std::vector<double> gradient18 = {-0.000000000000028, 1.194008590043404, -0.037599955385641};
        std::vector<double> gradient19 = {-0.000000000000007, -0.208475274559328, -0.868055152065867};
        std::vector<double> gradient20 = {0.000000000000007, 5.209001702131119, -2.451782812662477};
        std::vector<double> gradient21 = {18.097994437108714, 8.870576986489658, -22.406621153367627};
        std::vector<double> gradient22 = {-0.000000000000121, 3.895170376708266, 1.542062431706231};
        std::vector<double> gradient23 = {0.000000000000011, 0.657630120689917, -3.769122103737656};
        std::vector<double> gradient24 = {-0.000000000000085, 0.092840969189012, -0.002923610707697};
        std::vector<double> gradient25 = {-0.000000000000014, -0.905206038526934, -3.769122103737637};
        std::vector<double> gradient26 = {-0.000000000000085, -3.276230582115929, 1.542062431706356};
        std::vector<double> gradient27 = {-0.000000000000027, -16.827303500317022, -22.406621153367290};

        KRATOS_CHECK_VECTOR_NEAR(derivatives[0], positionVct, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[1], gradient1, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[2], gradient2, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[3], gradient3, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[4], gradient4, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[5], gradient5, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[6], gradient6, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[7], gradient7, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[8], gradient8, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[9], gradient9, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[10], gradient10, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[11], gradient11, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[12], gradient12, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[13], gradient13, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[14], gradient14, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[15], gradient15, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[16], gradient16, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[17], gradient17, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[18], gradient18, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[19], gradient19, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[20], gradient20, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[21], gradient21, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[22], gradient22, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[23], gradient23, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[24], gradient24, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[25], gradient25, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[26], gradient26, TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(derivatives[27], gradient27, TOLERANCE);
    }
} // namespace Testing.
} // namespace Kratos.
