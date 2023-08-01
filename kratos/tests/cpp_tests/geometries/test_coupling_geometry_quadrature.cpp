//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License          BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/coupling_geometry.h"

namespace Kratos {
    namespace Testing {

        typedef Node NodeType;
        typedef Geometry<NodeType> GeometryType;

        /// Generates a coupling geometry containing two nurbs curve geometries.
        GeometryType::Pointer GenerateNurbsCurveCouplingGeometry() {
            /// Curve 1
            PointerVector<NodeType> points_1(2);
            points_1(0) = Kratos::make_intrusive<NodeType>(1, 0, 0, 0);
            points_1(1) = Kratos::make_intrusive<NodeType>(2, 1.0, 0, 0);

            Vector knot_vector_1 = ZeroVector(2);
            knot_vector_1[0] = 0.0;
            knot_vector_1[1] = 1.0;

            auto p_line_1 = Kratos::make_shared<
                NurbsCurveGeometry<2, PointerVector<NodeType>>>(
                    points_1, 1, knot_vector_1);

            /// Curve 2
            PointerVector<NodeType> points_2(3);
            points_2(0) = Kratos::make_intrusive<NodeType>(1, 0.2, 0, 0);
            points_2(1) = Kratos::make_intrusive<NodeType>(2, 0.8, 0, 0);
            points_2(2) = Kratos::make_intrusive<NodeType>(3, 1.1, 0, 0);

            Vector knot_vector_2 = ZeroVector(3);
            knot_vector_2[0] = 0.0;
            knot_vector_2[1] = 0.6;
            knot_vector_2[2] = 0.9;

            auto p_line_2 = Kratos::make_shared<
                NurbsCurveGeometry<2, PointerVector<NodeType>>>(
                    points_2, 1, knot_vector_2);

            /// Coupling Geometry
            return Kratos::make_shared<CouplingGeometry<NodeType>>(
                p_line_1, p_line_2);
        }

        /// Generates a sample coupling geometry containing two line geometries.
        GeometryType::Pointer GenerateNurbsCurveOnSurfaceCouplingGeometry() {
            /// NurbsSurface 1:
            PointerVector<NodeType> points_surface_1(4);
            points_surface_1(0) = Kratos::make_intrusive<NodeType>(1, 0, 0, 0);
            points_surface_1(1) = Kratos::make_intrusive<NodeType>(2, 2, 0, 0);
            points_surface_1(2) = Kratos::make_intrusive<NodeType>(3, 0, 2, 0);
            points_surface_1(3) = Kratos::make_intrusive<NodeType>(4, 2, 2, 0);
            Vector knot_vector_surface_u_1 = ZeroVector(2);
            knot_vector_surface_u_1[0] = 0.0;
            knot_vector_surface_u_1[1] = 2.0;
            Vector knot_vector_surface_v_1 = ZeroVector(2);
            knot_vector_surface_v_1[0] = 0.0;
            knot_vector_surface_v_1[1] = 1.0;
            auto p_surface_1 = Kratos::make_shared<
                NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(
                    points_surface_1, 1, 1, knot_vector_surface_u_1, knot_vector_surface_v_1);
            /// NurbsCurve 1:
            PointerVector<Point> points_curve_1(3);
            points_curve_1(0) = Kratos::make_shared<Point>(1, 0.5, 0);
            points_curve_1(1) = Kratos::make_shared<Point>(1, 1, 0);
            points_curve_1(2) = Kratos::make_shared<Point>(2, 1, 0);
            Vector knot_vector_curve_1 = ZeroVector(4);
            knot_vector_curve_1[0] = 0.0;
            knot_vector_curve_1[1] = 0.0;
            knot_vector_curve_1[2] = 1.570796;
            knot_vector_curve_1[3] = 1.570796;
            Vector weights_curve_1 = ZeroVector(3);
            weights_curve_1[0] = 1;
            weights_curve_1[1] = 0.70710678;
            weights_curve_1[2] = 1;
            auto p_curve_1 = Kratos::make_shared<
                NurbsCurveGeometry<2, PointerVector<Point>>>(
                    points_curve_1, 2, knot_vector_curve_1, weights_curve_1);

            /// NurbsCurveOnSurface 1:
            auto p_curve_on_surface_1 = Kratos::make_shared<
                NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<NodeType>>>(
                    p_surface_1, p_curve_1);

            /// NurbsSurface 2:
            PointerVector<NodeType> points_surface_2(10);
            points_surface_2(0) = Kratos::make_intrusive<NodeType>(1, 1, 0, 0);
            points_surface_2(1) = Kratos::make_intrusive<NodeType>(2, 3, 0, 0);
            points_surface_2(2) = Kratos::make_intrusive<NodeType>(3, 1, 0.33333333, 0);
            points_surface_2(3) = Kratos::make_intrusive<NodeType>(4, 3, 0.33333333, 0);
            points_surface_2(4) = Kratos::make_intrusive<NodeType>(5, 1, 1, 0);
            points_surface_2(5) = Kratos::make_intrusive<NodeType>(6, 3, 1, 0);
            points_surface_2(6) = Kratos::make_intrusive<NodeType>(7, 1, 1.66666666, 0);
            points_surface_2(7) = Kratos::make_intrusive<NodeType>(8, 3, 1.66666666, 0);
            points_surface_2(8) = Kratos::make_intrusive<NodeType>(9, 1, 2, 0);
            points_surface_2(9) = Kratos::make_intrusive<NodeType>(10, 3, 2, 0);
            Vector knot_vector_surface_u_2 = ZeroVector(2);
            knot_vector_surface_u_2[0] = 0.0;
            knot_vector_surface_u_2[1] = 2.0;
            Vector knot_vector_surface_v_2 = ZeroVector(6);
            knot_vector_surface_v_2[0] = 0.0;
            knot_vector_surface_v_2[1] = 0.0;
            knot_vector_surface_v_2[2] = 0.66666666;
            knot_vector_surface_v_2[3] = 1.33333333;
            knot_vector_surface_v_2[4] = 2.0;
            knot_vector_surface_v_2[5] = 2.0;
            auto p_surface_2 = Kratos::make_shared<
                NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(
                    points_surface_2, 1, 2, knot_vector_surface_u_2, knot_vector_surface_v_2);
            /// NurbsCurve 2:
            PointerVector<Point> points_curve_2(4);
            points_curve_2(0) = Kratos::make_shared<Point>(1, 2, 0);
            points_curve_2(1) = Kratos::make_shared<Point>(0.44771525, 2, 0);
            points_curve_2(2) = Kratos::make_shared<Point>(0, 1.55228474983, 0);
            points_curve_2(3) = Kratos::make_shared<Point>(0, 1, 0);
            Vector knot_vector_curve_2 = ZeroVector(6);
            knot_vector_curve_2[0] = -1.570796;
            knot_vector_curve_2[1] = -1.570796;
            knot_vector_curve_2[2] = -1.570796;
            knot_vector_curve_2[3] = 0.0;
            knot_vector_curve_2[4] = 0.0;
            knot_vector_curve_2[5] = 0.0;
            auto p_curve_2 = Kratos::make_shared<
                NurbsCurveGeometry<2, PointerVector<Point>>>(
                    points_curve_2, 3, knot_vector_curve_2);

            /// NurbsCurveOnSurface 2:
            auto p_curve_on_surface_2 = Kratos::make_shared<
                NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<NodeType>>>(
                    p_surface_2, p_curve_2);

            /// CouplingGeometry
            return Kratos::make_shared<CouplingGeometry<NodeType>>(
                p_curve_on_surface_1, p_curve_on_surface_2);
        }

        /// Create integration points on nurbs curve.
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryCurvesCreateIntegrationPoints, KratosCoreCouplingGeometriesFastSuite) {
            auto p_coupling_geometry = GenerateNurbsCurveCouplingGeometry();

            typename GeometryType::IntegrationPointsArrayType integration_points;
            IntegrationInfo integration_info = p_coupling_geometry->GetDefaultIntegrationInfo();
            p_coupling_geometry->CreateIntegrationPoints(integration_points, integration_info);

            // 2 span intersections and each time (p=1) + 1 = 2 integration points.
            KRATOS_CHECK_EQUAL(integration_points.size(), 4);

            double length = 0;
            for (IndexType i = 0; i < integration_points.size(); ++i) {
                length += integration_points[i].Weight();
            }
            KRATOS_CHECK_NEAR(length, 0.8, TOLERANCE);
        }

        /// Create quadrature point geometries on nurbs curve.
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryCurvesCreateQuadraturePoints, KratosCoreCouplingGeometriesFastSuite) {
            auto p_coupling_geometry = GenerateNurbsCurveCouplingGeometry();

            typename GeometryType::GeometriesArrayType quadrature_points;
            IntegrationInfo integration_info = p_coupling_geometry->GetDefaultIntegrationInfo();
            p_coupling_geometry->CreateQuadraturePointGeometries(quadrature_points, 2, integration_info);

            // 2 span intersections and each time (p=1) + 1 integration point.
            KRATOS_CHECK_EQUAL(quadrature_points.size(), 4);

            double length = 0;
            for (IndexType i = 0; i < quadrature_points.size(); ++i) {
                for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                    length += quadrature_points[i].IntegrationPoints()[j].Weight();
                }
            }
            KRATOS_CHECK_NEAR(length, 0.8, TOLERANCE);
        }

        /// Create integration point geometries on nurbs curve on surface.
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryCurvesOnSurfaceCreateIntegrationPoints, KratosCoreCouplingGeometriesFastSuite) {
            auto p_coupling_geometry = GenerateNurbsCurveOnSurfaceCouplingGeometry();

            std::vector<double> spans;
            p_coupling_geometry->SpansLocalSpace(spans, 0);
            KRATOS_CHECK_EQUAL(spans.size(), 3);
            KRATOS_CHECK_NEAR(spans[0], 0.0, TOLERANCE);
            KRATOS_CHECK_NEAR(spans[1], 0.355753, TOLERANCE);
            KRATOS_CHECK_NEAR(spans[2], 1.570796, TOLERANCE);

            typename GeometryType::IntegrationPointsArrayType integration_points;
            IntegrationInfo integration_info = p_coupling_geometry->GetDefaultIntegrationInfo();
            p_coupling_geometry->CreateIntegrationPoints(integration_points, integration_info);

            // 2 span intersections and each time (p=1) + (q=2) + 1 integration point.
            KRATOS_CHECK_EQUAL(integration_info.GetNumberOfIntegrationPointsPerSpan(0), 4);
            KRATOS_CHECK_EQUAL(integration_points.size(), 8);

            double length = 0;
            for (IndexType i = 0; i < integration_points.size(); ++i) {
                length += integration_points[i].Weight();
            }
            KRATOS_CHECK_NEAR(length, 1.5707963, TOLERANCE);

            // Check with modified number of integration points per span.
            integration_info.SetNumberOfIntegrationPointsPerSpan(0, 8);
            p_coupling_geometry->CreateIntegrationPoints(integration_points, integration_info);

            KRATOS_CHECK_EQUAL(integration_info.GetNumberOfIntegrationPointsPerSpan(0), 8);
            KRATOS_CHECK_EQUAL(integration_points.size(), 16);

            double length2 = 0;
            for (IndexType i = 0; i < integration_points.size(); ++i) {
                length2 += integration_points[i].Weight();
            }
            KRATOS_CHECK_NEAR(length2, 1.5707963, TOLERANCE);
        }

        /// Create quadrature points on nurbs curve on surface.
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryCurvesOnSurfaceCreateQuadraturePoints, KratosCoreCouplingGeometriesFastSuite) {
            auto p_coupling_geometry = GenerateNurbsCurveOnSurfaceCouplingGeometry();

            typename GeometryType::GeometriesArrayType quadrature_points;
            IntegrationInfo integration_info = p_coupling_geometry->GetDefaultIntegrationInfo();
            p_coupling_geometry->CreateQuadraturePointGeometries(quadrature_points, 2, integration_info);

            // 2 span intersections and each time (p=1) + (q=2) + 1 integration point.
            KRATOS_CHECK_EQUAL(quadrature_points.size(), 8);

            double length = 0;
            for (IndexType i = 0; i < quadrature_points.size(); ++i) {
                for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                    length += quadrature_points[i].IntegrationPoints()[j].Weight();
                }
            }
            KRATOS_CHECK_NEAR(length, 1.5707963, TOLERANCE);
        }
    } // namespace Testing.
} // namespace Kratos.
