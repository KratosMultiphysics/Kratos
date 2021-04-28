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
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_surface.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    typename NurbsSurfaceGeometry<3, PointerVector<NodeType>>::Pointer GenerateNurbsSurface() {
        Geometry<NodeType>::PointsArrayType points(4);

        points(0) = Kratos::make_shared<NodeType>(1, 0, 0, 0);
        points(1) = Kratos::make_shared<NodeType>(2, 10, 0, 10);
        points(2) = Kratos::make_shared<NodeType>(3, 10, 5, 10);
        points(3) = Kratos::make_shared<NodeType>(4, 0, 5, 10);

        Vector knot_u = ZeroVector(2);
        knot_u[0] = 0.0;
        knot_u[1] = 1.0;
        Vector knot_v = ZeroVector(2); 
        knot_v[0] = 0.0;
        knot_v[1] = 1.0;

        const int p = 1;
        const int q = 1;

        return Kratos::make_shared < NurbsSurfaceGeometry<3, PointerVector<NodeType> >>(points, p, q, knot_u, knot_v);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateNurbsCurve1()
    {
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType points(2);

        points(0) = Kratos::make_shared<NodeType>(1, 0, 0, 0);
        points(1) = Kratos::make_shared<NodeType>(2, 1, 0, 10);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateNurbsCurve2()
    {
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType points(2);

        points(0) = Kratos::make_shared<NodeType>(1, 1, 0, 0);
        points(1) = Kratos::make_shared<NodeType>(2, 1, 1, 10);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateNurbsCurve3()
    {
        NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType points(2);

        points(0) = Kratos::make_shared<NodeType>(1, 1, 1, 0);
        points(1) = Kratos::make_shared<NodeType>(2, 0, 0, 10);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename BrepSurface<PointerVector<NodeType>, PointerVector<Point>>::Pointer GenerateTrimmedBrepSurface()
    {
        auto p_surface = GenerateNurbsSurface();
        auto p_curve_1 = GenerateNurbsCurve1();
        auto p_curve_2 = GenerateNurbsCurve2();
        auto p_curve_3 = GenerateNurbsCurve3();

        auto p_brep_curve_on_surface_1 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, PointerVector<Point>>>(
            p_surface, p_curve_1);
        auto p_brep_curve_on_surface_2 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, PointerVector<Point>>>(
            p_surface, p_curve_2);
        auto p_brep_curve_on_surface_3 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, PointerVector<Point>>>(
            p_surface, p_curve_3);

        BrepSurface<PointerVector<NodeType>, PointerVector<Point>>::BrepCurveOnSurfaceLoopType outer_loop(3);
        outer_loop[0] = p_brep_curve_on_surface_1;
        outer_loop[1] = p_brep_curve_on_surface_2;
        outer_loop[2] = p_brep_curve_on_surface_3;

        BrepSurface<PointerVector<NodeType>, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType outer_loops(1);
        outer_loops[0] = outer_loop;
        BrepSurface<PointerVector<NodeType>, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType inner_loops(0);

        return Kratos::make_shared<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(
            p_surface, outer_loops, inner_loops);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsTrimmedBrepSurfaceIntegration, KratosCoreGeometriesFastSuite) {
        auto p_brep_surface = GenerateTrimmedBrepSurface();

        //// Check general information, input to ouput
        KRATOS_CHECK_EQUAL(p_brep_surface->Dimension(), 2);
        KRATOS_CHECK_EQUAL(p_brep_surface->WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(p_brep_surface->LocalSpaceDimension(), 2);

        typename Geometry<Node<3>>::IntegrationPointsArrayType integration_points;
        p_brep_surface->CreateIntegrationPoints(integration_points);

        typename Geometry<Node<3>>::GeometriesArrayType quadrature_points;
        p_brep_surface->CreateQuadraturePointGeometries(quadrature_points, 3, integration_points);

        double area = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                area += quadrature_points[i].IntegrationPoints()[j].Weight();
            }
        }
        KRATOS_CHECK_NEAR(area, 50, TOLERANCE);
    }
} // namespace Testing.
} // namespace Kratos.
