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
#include "containers/pointer_vector.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/brep_surface.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node NodeType;

    typename NurbsSurfaceGeometry<3, PointerVector<NodeType>>::Pointer GenerateNurbsSurface() {
        PointerVector<NodeType> points(4);

        points(0) = Kratos::make_intrusive<NodeType>(1, 0, 0, 0);
        points(1) = Kratos::make_intrusive<NodeType>(2, 0, 5, 0);
        points(2) = Kratos::make_intrusive<NodeType>(3, 10, 0, 0);
        points(3) = Kratos::make_intrusive<NodeType>(4, 10, 5, 0);

        Vector knot_u = ZeroVector(2);
        knot_u[0] = 0.0;
        knot_u[1] = 1.0;
        Vector knot_v = ZeroVector(2); 
        knot_v[0] = 0.0;
        knot_v[1] = 1.0;

        const int p = 1;
        const int q = 1;

        return Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType> >>(points, p, q, knot_u, knot_v);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateNurbsCurve1()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(0, 0, 0);
        points(1) = Kratos::make_shared<Point>(1, 0, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateNurbsCurve2()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(1, 0, 0);
        points(1) = Kratos::make_shared<Point>(1, 1, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateNurbsCurve3()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(1, 1, 0);
        points(1) = Kratos::make_shared<Point>(0, 0, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateTrimmedBrepSurface()
    {
        auto p_surface = GenerateNurbsSurface();
        auto p_curve_1 = GenerateNurbsCurve1();
        auto p_curve_2 = GenerateNurbsCurve2();
        auto p_curve_3 = GenerateNurbsCurve3();

        auto p_brep_curve_on_surface_1 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_1);
        auto p_brep_curve_on_surface_2 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_2);
        auto p_brep_curve_on_surface_3 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_3);

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopType outer_loop(3);
        outer_loop[0] = p_brep_curve_on_surface_1;
        outer_loop[1] = p_brep_curve_on_surface_2;
        outer_loop[2] = p_brep_curve_on_surface_3;

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType outer_loops(1);
        outer_loops[0] = outer_loop;
        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType inner_loops(0);

        return Kratos::make_shared<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, outer_loops, inner_loops);
    }

    typename NurbsSurfaceGeometry<3, PointerVector<NodeType>>::Pointer GenerateReferenceNurbsSurface() {
        PointerVector<NodeType> points(9);

        points(0) = Kratos::make_intrusive<NodeType>(1, 0, 0, 0);
        points(1) = Kratos::make_intrusive<NodeType>(2, 2.5, 0, 0);
        points(2) = Kratos::make_intrusive<NodeType>(3, 5, 0, 0);
        points(3) = Kratos::make_intrusive<NodeType>(4, 0, 0.5, 0);
        points(4) = Kratos::make_intrusive<NodeType>(5, 2.5, 0.5, 0);
        points(5) = Kratos::make_intrusive<NodeType>(6, 5, 0.5, 0);
        points(6) = Kratos::make_intrusive<NodeType>(7, 0, 1, 0);
        points(7) = Kratos::make_intrusive<NodeType>(8, 2.5, 1, 0);
        points(8) = Kratos::make_intrusive<NodeType>(9, 5, 1, 0);

        Vector knot_u = ZeroVector(4);
        knot_u[0] = 0.0;
        knot_u[1] = 0.0;
        knot_u[2] = 5.0;
        knot_u[3] = 5.0;
        Vector knot_v = ZeroVector(4); 
        knot_v[0] = 0.0;
        knot_v[1] = 0.0;
        knot_v[2] = 1.0;
        knot_v[3] = 1.0;

        const int p = 2;
        const int q = 2;

        return Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType> >>(points, p, q, knot_u, knot_v);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve1()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(4.333333333333333, 1, 0);
        points(1) = Kratos::make_shared<Point>(0, 1, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.666666666666666;
        knot_vector[1] = 5;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve2()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(0, 1, 0);
        points(1) = Kratos::make_shared<Point>(0, 0, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve3()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(0, 0, 0);
        points(1) = Kratos::make_shared<Point>(4.666666666666666, 0, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 4.666666666666666;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve4()
    {
        PointerVector<Point> points(4);

        points(0) = Kratos::make_shared<Point>(4.666666666666666, 0, 0);
        points(1) = Kratos::make_shared<Point>(4.555555555555555, 0.3333333333333333, 0);
        points(2) = Kratos::make_shared<Point>(4.444444444444444, 0.6666666666666666, 0);
        points(3) = Kratos::make_shared<Point>(4.333333333333333, 1, 0);

        Vector knot_vector = ZeroVector(6);
        knot_vector[0] = 1.0540925533894603;
        knot_vector[1] = 1.0540925533894603;
        knot_vector[2] = 1.0540925533894603;
        knot_vector[3] = 2.1081851067789192;
        knot_vector[4] = 2.1081851067789192;
        knot_vector[5] = 2.1081851067789192;

        int p = 3;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve5()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(0, 0, 0);
        points(1) = Kratos::make_shared<Point>(5, 0, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 5;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }
    
    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve6()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(5, 0, 0);
        points(1) = Kratos::make_shared<Point>(5, 1, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve7()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(5, 1, 0);
        points(1) = Kratos::make_shared<Point>(0, 1, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 5;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve8()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(0, 1, 0);
        points(1) = Kratos::make_shared<Point>(0, 0, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve9()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(3.0, 0.25, 0);
        points(1) = Kratos::make_shared<Point>(2.0, 0.25, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }
    
    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve10()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(2, 0.25, 0);
        points(1) = Kratos::make_shared<Point>(2, 0.75, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.5;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve11()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(2, 0.75, 0);
        points(1) = Kratos::make_shared<Point>(3, 0.75, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 1;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReferenceNurbsCurve12()
    {
        PointerVector<Point> points(2);

        points(0) = Kratos::make_shared<Point>(3, 0.75, 0);
        points(1) = Kratos::make_shared<Point>(3, 0.25, 0);

        Vector knot_vector = ZeroVector(2);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.5;

        int p = 1;

        return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    }

    typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateOuterInclinedTrimmedBrepSurface()
    {
        auto p_surface = GenerateReferenceNurbsSurface();
        auto p_curve_1 = GenerateReferenceNurbsCurve1();
        auto p_curve_2 = GenerateReferenceNurbsCurve2();
        auto p_curve_3 = GenerateReferenceNurbsCurve3();
        auto p_curve_4 = GenerateReferenceNurbsCurve4();

        auto p_brep_curve_on_surface_1 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_1);
        auto p_brep_curve_on_surface_2 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_2);
        auto p_brep_curve_on_surface_3 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_3);
        auto p_brep_curve_on_surface_4 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_4);

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopType outer_loop(4);
        outer_loop[0] = p_brep_curve_on_surface_1;
        outer_loop[1] = p_brep_curve_on_surface_2;
        outer_loop[2] = p_brep_curve_on_surface_3;
        outer_loop[3] = p_brep_curve_on_surface_4;

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType outer_loops(1);
        outer_loops[0] = outer_loop;
        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType inner_loops(0);

        return Kratos::make_shared<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, outer_loops, inner_loops);
    }

    typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateOuterInclinedRefinedTrimmedBrepSurface()
    {
        auto p_surface = GenerateReferenceNurbsSurface();

        //refinement
        std::vector<double> spans_local_space_u;
        p_surface->SpansLocalSpace(spans_local_space_u, 0);
        std::vector<double>  spans_local_space_v;
        p_surface->SpansLocalSpace(spans_local_space_v, 1);
        int nb_per_span = 6;

        std::vector<double> knots_to_insert_u;
        for (IndexType i = 0; i < spans_local_space_u.size() - 1; ++i) {
            const double delta_u = (spans_local_space_u[i + 1] - spans_local_space_u[i]) / (nb_per_span + 1);
            for (IndexType j = 1; j < nb_per_span + 1; ++j) {
                knots_to_insert_u.push_back(spans_local_space_u[i] + delta_u * j);
            }
        }
        std::vector<double>  knots_to_insert_v;
        for (IndexType i = 0; i < spans_local_space_v.size() - 1; ++i) {
            const double delta_v = (spans_local_space_v[i + 1] - spans_local_space_v[i]) / (nb_per_span + 1);
            for (IndexType j = 1; j < nb_per_span + 1; ++j) {
                knots_to_insert_v.push_back(spans_local_space_v[i] + delta_v * j);
            }
        }

        PointerVector<NodeType> PointsRefined;
        Vector KnotsURefined;
        Vector KnotsVRefined;
        Vector WeightsRefined;  

        NurbsSurfaceRefinementUtilities::KnotRefinementU(*p_surface, knots_to_insert_u,
            PointsRefined, KnotsURefined, WeightsRefined);
        p_surface->SetInternals(PointsRefined,
            p_surface->PolynomialDegreeU(), p_surface->PolynomialDegreeV(),
            KnotsURefined, p_surface->KnotsV(),
            WeightsRefined);

        NurbsSurfaceRefinementUtilities::KnotRefinementV(*p_surface, knots_to_insert_v,
            PointsRefined, KnotsVRefined, WeightsRefined);
        p_surface->SetInternals(PointsRefined,
            p_surface->PolynomialDegreeU(), p_surface->PolynomialDegreeV(),
            p_surface->KnotsU(), KnotsVRefined,
            WeightsRefined);

        auto p_curve_1 = GenerateReferenceNurbsCurve1();
        auto p_curve_2 = GenerateReferenceNurbsCurve2();
        auto p_curve_3 = GenerateReferenceNurbsCurve3();
        auto p_curve_4 = GenerateReferenceNurbsCurve4();

        auto p_brep_curve_on_surface_1 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_1);
        auto p_brep_curve_on_surface_2 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_2);
        auto p_brep_curve_on_surface_3 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_3);
        auto p_brep_curve_on_surface_4 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_4);

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopType outer_loop(4);
        outer_loop[0] = p_brep_curve_on_surface_1;
        outer_loop[1] = p_brep_curve_on_surface_2;
        outer_loop[2] = p_brep_curve_on_surface_3;
        outer_loop[3] = p_brep_curve_on_surface_4;

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType outer_loops(1);
        outer_loops[0] = outer_loop;
        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType inner_loops(0);

        return Kratos::make_shared<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, outer_loops, inner_loops);
    }

    typename BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer GenerateInnerRectRefinedTrimmedBrepSurface()
    {
        auto p_surface = GenerateReferenceNurbsSurface();

        //refinement
        std::vector<double> knots_to_insert_u;
        knots_to_insert_u.push_back(2.5);
        std::vector<double> knots_to_insert_v;
        knots_to_insert_v.push_back(0.5);

        PointerVector<NodeType> PointsRefined;
        Vector KnotsURefined;
        Vector KnotsVRefined;
        Vector WeightsRefined;  

        NurbsSurfaceRefinementUtilities::KnotRefinementU(*p_surface, knots_to_insert_u,
            PointsRefined, KnotsURefined, WeightsRefined);
        p_surface->SetInternals(PointsRefined,
            p_surface->PolynomialDegreeU(), p_surface->PolynomialDegreeV(),
            KnotsURefined, p_surface->KnotsV(),
            WeightsRefined);

        NurbsSurfaceRefinementUtilities::KnotRefinementV(*p_surface, knots_to_insert_v,
            PointsRefined, KnotsVRefined, WeightsRefined);
        p_surface->SetInternals(PointsRefined,
            p_surface->PolynomialDegreeU(), p_surface->PolynomialDegreeV(),
            p_surface->KnotsU(), KnotsVRefined,
            WeightsRefined);

        auto p_curve_5 = GenerateReferenceNurbsCurve5();
        auto p_curve_6 = GenerateReferenceNurbsCurve6();
        auto p_curve_7 = GenerateReferenceNurbsCurve7();
        auto p_curve_8 = GenerateReferenceNurbsCurve8();
        auto p_curve_9 = GenerateReferenceNurbsCurve9();
        auto p_curve_10 = GenerateReferenceNurbsCurve10();
        auto p_curve_11 = GenerateReferenceNurbsCurve11();
        auto p_curve_12 = GenerateReferenceNurbsCurve12();

        auto p_brep_curve_on_surface_5 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_5);
        auto p_brep_curve_on_surface_6 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_6);
        auto p_brep_curve_on_surface_7 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_7);
        auto p_brep_curve_on_surface_8 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_8);
        auto p_brep_curve_on_surface_9 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_9);
        auto p_brep_curve_on_surface_10 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_10);
        auto p_brep_curve_on_surface_11 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_11);
        auto p_brep_curve_on_surface_12 = Kratos::make_shared<BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, p_curve_12);

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopType outer_loop(4);
        outer_loop[0] = p_brep_curve_on_surface_5;
        outer_loop[1] = p_brep_curve_on_surface_6;
        outer_loop[2] = p_brep_curve_on_surface_7;
        outer_loop[3] = p_brep_curve_on_surface_8;

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopType inner_loop(4);
        inner_loop[0] = p_brep_curve_on_surface_9;
        inner_loop[1] = p_brep_curve_on_surface_10;
        inner_loop[2] = p_brep_curve_on_surface_11;
        inner_loop[3] = p_brep_curve_on_surface_12;

        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType outer_loops(1);
        outer_loops[0] = outer_loop;
        BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>::BrepCurveOnSurfaceLoopArrayType inner_loops(1);
        inner_loops[0] = inner_loop;

        return Kratos::make_shared<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(
            p_surface, outer_loops, inner_loops);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsTrimmedBrepSurfaceIntegration, KratosCoreGeometriesFastSuite) {
        auto p_brep_surface = GenerateTrimmedBrepSurface();

        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = p_brep_surface->GetDefaultIntegrationInfo();
        p_brep_surface->CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Node>::GeometriesArrayType quadrature_points;
        p_brep_surface->CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        double area = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                area += quadrature_points[i].IntegrationPoints()[j].Weight() * quadrature_points[i].DeterminantOfJacobian(0);
            }
        }
        KRATOS_EXPECT_NEAR(area, 25, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsOuterInclinedTrimmedBrepSurfaceIntegration, KratosCoreGeometriesFastSuite) {
        auto p_brep_surface = GenerateOuterInclinedTrimmedBrepSurface();

        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = p_brep_surface->GetDefaultIntegrationInfo();
        p_brep_surface->CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Node>::GeometriesArrayType quadrature_points;
        p_brep_surface->CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        double area = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                area += quadrature_points[i].IntegrationPoints()[j].Weight() * quadrature_points[i].DeterminantOfJacobian(0);
            }
        }
        KRATOS_EXPECT_NEAR(area, 4.5, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsOuterInclinedRefinedTrimmedBrepSurfaceIntegration, KratosCoreGeometriesFastSuite) {
        auto p_brep_surface = GenerateOuterInclinedRefinedTrimmedBrepSurface();

        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = p_brep_surface->GetDefaultIntegrationInfo();
        p_brep_surface->CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Node>::GeometriesArrayType quadrature_points;
        p_brep_surface->CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        double area = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                area += quadrature_points[i].IntegrationPoints()[j].Weight() * quadrature_points[i].DeterminantOfJacobian(0);
            }
        }
        KRATOS_EXPECT_NEAR(area, 4.5, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsInnerRectRefinedTrimmedBrepSurfaceIntegration, KratosCoreGeometriesFastSuite) {
        auto p_brep_surface = GenerateInnerRectRefinedTrimmedBrepSurface();

        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = p_brep_surface->GetDefaultIntegrationInfo();
        p_brep_surface->CreateIntegrationPoints(integration_points, integration_info);

        typename Geometry<Node>::GeometriesArrayType quadrature_points;
        p_brep_surface->CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);

        double area = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                area += quadrature_points[i].IntegrationPoints()[j].Weight() * quadrature_points[i].DeterminantOfJacobian(0);
            }
        }
        KRATOS_EXPECT_NEAR(area, 4.5, TOLERANCE);
    }
} // namespace Testing.
} // namespace Kratos.
