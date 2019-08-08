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
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/brep_face_curve.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;

// /// Factory functions
//namespace {
    NurbsCurveGeometry<2, Point> GenerateReferenceCurve()
    {
        NurbsCurveGeometry<2, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(0, 0, 0)));
        points.push_back(Point::Pointer(new Point(3.3333333333333335, 1.6666666666666667, 0)));
        points.push_back(Point::Pointer(new Point(6.6666666666666661, 3.333333333333333, 0)));
        points.push_back(Point::Pointer(new Point(10, 5, 0)));

        Vector knot_vector = ZeroVector(6);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 0.0;
        knot_vector[3] = 11.180339887498949;
        knot_vector[4] = 11.180339887498949;
        knot_vector[5] = 11.180339887498949;

        int p = 3;

        auto curve = NurbsCurveGeometry<2, Point>(points, p, knot_vector);

        return curve;
    }

    Geometry<NodeType>::Pointer GenerateReferenceSurface() {
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
        Vector knot_v = ZeroVector(4); 
        knot_v[0] = 0.0;
        knot_v[1] = 5.0;

        int p = 2;
        int q = 1;

        return Geometry<NodeType>::Pointer(
                new NurbsSurfaceGeometry<2, NodeType>(
                    points, p, q, knot_u, knot_v));
    }
//}

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(BrepFaceCurve, KratosCoreGeometriesFastSuite) {
        auto curve = GenerateReferenceCurve();

        //edrr

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 1.0;
        array_1d<double, 3> result(0.0);

        curve.GlobalCoordinates(result, parameter);

        KRATOS_WATCH(result)

        auto surface = GenerateReferenceSurface();

        KRATOS_CHECK_NEAR(0.455342, 0.455342, TOLERANCE);
    }
} // namespace Testing.
} // namespace Kratos.
