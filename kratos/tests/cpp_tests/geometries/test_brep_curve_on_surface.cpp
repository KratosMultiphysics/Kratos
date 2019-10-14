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
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_face_curve.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    NurbsSurfaceGeometry<3, NodeType>::Pointer GenerateReferenceNodeSurfacePointer() {
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

        return Kratos::make_shared < NurbsSurfaceGeometry<3, NodeType >>(points, p, q, knot_u, knot_v);
    }

    NurbsCurveGeometry<2, Point>::Pointer GenerateReferenceCurve2dPointer()
    {
        NurbsCurveGeometry<2, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(0, 0)));
        points.push_back(Point::Pointer(new Point(3.3333333333333335, 1.6666666666666667)));
        points.push_back(Point::Pointer(new Point(6.6666666666666661, 3.333333333333333)));
        points.push_back(Point::Pointer(new Point(10, 5)));

        Vector knot_vector = ZeroVector(6);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 0.0;
        knot_vector[3] = 11.180339887498949;
        knot_vector[4] = 11.180339887498949;
        knot_vector[5] = 11.180339887498949;

        int p = 3;

        return Kratos::make_shared<NurbsCurveGeometry<2, Point>>(points, p, knot_vector);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsBrepFaceCurve, KratosCoreGeometriesFastSuite) {
        auto p_surface = GenerateReferenceNodeSurfacePointer();
        auto p_curve = GenerateReferenceCurve2dPointer();

        auto brep_face_curve = BrepFaceCurve<NodeType, Point>(
            p_surface, p_curve);

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(brep_face_curve.Dimension(), 1);
        KRATOS_CHECK_EQUAL(brep_face_curve.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(brep_face_curve.LocalSpaceDimension(), 2);


        array_1d<double, 3> coords(3, 0.0);
        coords[0] = 1.0;
        Vector N;

        //brep_face_curve.ShapeFunctionsValues(N, coords);
        //KRATOS_WATCH(N)
        //Matrix DN_De;
        //brep_face_curve.ShapeFunctionsLocalGradients(DN_De, coords);
        //KRATOS_WATCH(DN_De)

        //array_1d<double, 3> result(3, 0.0);
        //brep_face_curve.GlobalCoordinates(result, coords);
        //KRATOS_WATCH(result)

        //auto results = p_brep_face_curve.GlobalDerivatives(coords, 3);
        //KRATOS_WATCH(results[0])
    }
} // namespace Testing.
} // namespace Kratos.
