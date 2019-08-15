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
#include "geometries/brep_face.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    NurbsSurfaceGeometry<3, NodeType>::Pointer GenerateReferenceNodeSurfaceHalfCirclePointer() {
        Geometry<NodeType>::PointsArrayType points;

        points.push_back(NodeType::Pointer(new NodeType(1, 10, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(2, 10, 0, 10)));
        points.push_back(NodeType::Pointer(new NodeType(3, 0, 0, 10)));
        points.push_back(NodeType::Pointer(new NodeType(4, -10, 0, 10)));
        points.push_back(NodeType::Pointer(new NodeType(5, -10, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(6, 10, 10, 0)));
        points.push_back(NodeType::Pointer(new NodeType(7, 10, 10, 10)));
        points.push_back(NodeType::Pointer(new NodeType(8, 0, 10, 10)));
        points.push_back(NodeType::Pointer(new NodeType(9, -10, 10, 10)));
        points.push_back(NodeType::Pointer(new NodeType(10, -10, 10, 0)));

        Vector knot_u = ZeroVector(6);
        knot_u[0] = 0.0;
        knot_u[1] = 0.0;
        knot_u[2] = 15.707963267948966;
        knot_u[3] = 15.707963267948966;
        knot_u[4] = 31.415926535897931;
        knot_u[5] = 31.415926535897931;
        Vector knot_v = ZeroVector(2); 
        knot_v[0] = 0.0;
        knot_v[1] = 10.0;

        int p = 2;
        int q = 1;

        Vector weights = ZeroVector(10);
        weights[0] = 1.0;
        weights[1] = 0.70710678118654757;
        weights[2] = 1.0;
        weights[3] = 0.70710678118654757;
        weights[4] = 1.0;
        weights[5] = 1.0;
        weights[6] = 0.70710678118654757;
        weights[7] = 1.0;
        weights[8] = 0.70710678118654757;
        weights[9] = 1.0;

        return Kratos::make_shared < NurbsSurfaceGeometry<3, NodeType >>(points, p, q, knot_u, knot_v, weights);
    }

    NurbsCurveGeometry<2, Point>::Pointer GenerateReference1Curve2dPointer()
    {
        NurbsCurveGeometry<2, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(0, 0)));
        points.push_back(Point::Pointer(new Point(31.415926535897931, 0)));

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 31.415926535897931;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, Point>>(points, p, knot_vector);
    }

    NurbsCurveGeometry<2, Point>::Pointer GenerateReference2Curve2dPointer()
    {
        NurbsCurveGeometry<2, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(31.415926535897931, 0)));
        points.push_back(Point::Pointer(new Point(31.415926535897931, 10)));

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 10.0;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, Point>>(points, p, knot_vector);
    }

    NurbsCurveGeometry<2, Point>::Pointer GenerateReference3Curve2dPointer()
    {
        NurbsCurveGeometry<2, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(31.415926535897931, 10)));
        points.push_back(Point::Pointer(new Point(28.500619950938574, 10)));
        points.push_back(Point::Pointer(new Point(26.033141879879658, 9.48077148945451)));
        points.push_back(Point::Pointer(new Point(23.323541783172974, 8.4455316972070165)));
        points.push_back(Point::Pointer(new Point(20.805481541235128, 7.4834723273026622)));
        points.push_back(Point::Pointer(new Point(18.157204808069945, 6.1025460603903392)));
        points.push_back(Point::Pointer(new Point(15.707963267948966, 5)));
        points.push_back(Point::Pointer(new Point(13.258721727828151, 3.8974539396097354)));
        points.push_back(Point::Pointer(new Point(10.610444994662638, 2.5165276726972694)));
        points.push_back(Point::Pointer(new Point(8.09238475272496,   1.5544683027929818)));
        points.push_back(Point::Pointer(new Point(5.382784656018293,  0.51922851054549679)));
        points.push_back(Point::Pointer(new Point(2.91530658495934,   0)));
        points.push_back(Point::Pointer(new Point(0, 0)));

        Vector knot_vector = ZeroVector(15);
        knot_vector[0] =  -34.16313447529194;
        knot_vector[1] =  -34.16313447529194;
        knot_vector[2] =  -34.16313447529194;
        knot_vector[3] =  -25.309419191442988;
        knot_vector[4] =  -25.309419191442988;
        knot_vector[5] =  -25.309419191442988;
        knot_vector[6] =  -17.081567237645665;
        knot_vector[7] =  -17.081567237645665;
        knot_vector[8] =  -17.081567237645665;
        knot_vector[9] =  -8.8537152838488975;
        knot_vector[10] = -8.8537152838488975;
        knot_vector[11] = -8.8537152838488975;
        knot_vector[12] = 0;
        knot_vector[13] = 0;
        knot_vector[14] = 0;

        int p = 3;

        return Kratos::make_shared<NurbsCurveGeometry<2, Point>>(points, p, knot_vector);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsBrepFace, KratosCoreGeometriesFastSuite) {
        auto p_surface = GenerateReferenceNodeSurfaceHalfCirclePointer();
        auto p_curve_1 = GenerateReference1Curve2dPointer();
        auto p_curve_2 = GenerateReference2Curve2dPointer();
        auto p_curve_3 = GenerateReference3Curve2dPointer();

        auto p_brep_face_curve_1 = Kratos::make_shared<BrepFaceCurve<NodeType, Point>>(
            p_surface, p_curve_1);
        auto p_brep_face_curve_2 = Kratos::make_shared<BrepFaceCurve<NodeType, Point>>(
            p_surface, p_curve_2);
        auto p_brep_face_curve_3 = Kratos::make_shared<BrepFaceCurve<NodeType, Point>>(
            p_surface, p_curve_3);

        BrepFace<NodeType, Point>::BrepFaceCurveLoopType outer_loop(3);
        outer_loop[0] = p_brep_face_curve_1;
        outer_loop[1] = p_brep_face_curve_2;
        outer_loop[2] = p_brep_face_curve_3;

        BrepFace<NodeType, Point>::BrepFaceCurveLoopArrayType outer_loops(1);
        outer_loops[0] = outer_loop;

        BrepFace<NodeType, Point>::BrepFaceCurveLoopArrayType inner_loops(0);
        //std::vector<std::vector<BrepFaceCurve<NodeType, Point>::Pointer>> inner_loops(0);

        auto brep_face = BrepFace<NodeType, Point>(
            p_surface, outer_loops, inner_loops);

        //// Check general information, input to ouput
        KRATOS_CHECK_EQUAL(brep_face.Dimension(), 2);
        KRATOS_CHECK_EQUAL(brep_face.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(brep_face.LocalSpaceDimension(), 2);

        KRATOS_WATCH(brep_face.Dimension())
            KRATOS_WATCH(brep_face.WorkingSpaceDimension())
            KRATOS_WATCH(brep_face.LocalSpaceDimension())

        //array_1d<double, 3> coords(3, 0.0);
        //coords[0] = 1.0;
        //Vector N;

        //p_brep_face_curve.ShapeFunctionsValues(N, coords);
        //KRATOS_WATCH(N)
        //Matrix DN_De;
        //p_brep_face_curve.ShapeFunctionsLocalGradients(DN_De, coords);
        //KRATOS_WATCH(DN_De)

        //array_1d<double, 3> result(3, 0.0);
        //p_brep_face_curve.GlobalCoordinates(result, coords);
        //KRATOS_WATCH(result)

        //auto results = p_brep_face_curve.GlobalDerivatives(coords, 3);
        //KRATOS_WATCH(results[0])
    }
} // namespace Testing.
} // namespace Kratos.
