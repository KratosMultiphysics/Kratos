//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Wataru Fukuda
//

// System includes
#include <limits>
#include <algorithm>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/thb_surface_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node NodeType;
    using THBSurface3D = THBSurfaceGeometry<3, PointerVector<NodeType>>;

    /// Factory functions

    THBSurface3D::Pointer GenerateReferenceNodeTHBSurface()
    {
        PointerVector<NodeType> points;

        points.push_back(NodeType::Pointer(new NodeType(1,  0,  5,  0)));
        points.push_back(NodeType::Pointer(new NodeType(2,  5,  5,  0)));
        points.push_back(NodeType::Pointer(new NodeType(3, 10,  5, -4)));
        points.push_back(NodeType::Pointer(new NodeType(4,  0,  0,  0)));
        points.push_back(NodeType::Pointer(new NodeType(5,  5,  0,  0)));
        points.push_back(NodeType::Pointer(new NodeType(6, 10,  0, -4)));

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

        return Kratos::make_shared<THBSurface3D>(points, p, q, knot_u, knot_v);
    }

    THBSurface3D::Pointer GenerateRefinedNodeTHBSurface()
    {
        auto p_thb = GenerateReferenceNodeTHBSurface();
        p_thb->AddLevelByHRefinement(1);
        p_thb->AddRefinementDomain(1, 0.0, 10.0, 0.0, 5.0);
        return p_thb;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceNodeInfo, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateReferenceNodeTHBSurface();

        KRATOS_EXPECT_EQ(p_thb->WorkingSpaceDimension(), 3);
        KRATOS_EXPECT_EQ(p_thb->LocalSpaceDimension(), 2);

        KRATOS_EXPECT_EQ(p_thb->PolynomialDegree(0), 2);
        KRATOS_EXPECT_EQ(p_thb->PolynomialDegree(1), 1);

        KRATOS_EXPECT_EQ(p_thb->NumberOfLevels(), 1);
        KRATOS_EXPECT_EQ(p_thb->PointsNumber(), 6);
        KRATOS_EXPECT_EQ(p_thb->NumberOfActiveFunctions(0), 6);
    }

    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceNodeGlobalCoordinates, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateReferenceNodeTHBSurface();

        array_1d<double, 3> parameter(0.0);
        array_1d<double, 3> result(0.0);

        parameter[0] = 10.0;
        parameter[1] = 3.5;
        p_thb->GlobalCoordinates(result, parameter);
        KRATOS_EXPECT_NEAR(result[0], 10.0,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[1],  1.5,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[2], -4.0,  TOLERANCE);

        parameter[0] = 6.0;
        parameter[1] = 1.0;
        p_thb->GlobalCoordinates(result, parameter);
        KRATOS_EXPECT_NEAR(result[0],  6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(result[1],  4.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(result[2], -1.44,  TOLERANCE);

        parameter[0] = 0.0;
        parameter[1] = 1.0;
        p_thb->GlobalCoordinates(result, parameter);
        KRATOS_EXPECT_NEAR(result[0], 0.0,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[1], 4.0,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[2], 0.0,  TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceNodeCreateIntegrationPoints, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateReferenceNodeTHBSurface();

        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = p_thb->GetDefaultIntegrationInfo();
        p_thb->CreateIntegrationPoints(integration_points, integration_info);

        KRATOS_EXPECT_EQ(integration_points.size(), 6);
        double area = 0.0;
        for (IndexType i = 0; i < integration_points.size(); ++i)
            area += integration_points[i].Weight();
        KRATOS_EXPECT_NEAR(area, 50.0, TOLERANCE);
    }

    /// Checks that AddLevelByHRefinement grows the level count and CP pool.
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceAddHRefinementLevel, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateReferenceNodeTHBSurface();
        p_thb->AddLevelByHRefinement(1);

        KRATOS_EXPECT_EQ(p_thb->NumberOfLevels(), 2);

        // 6 coarse + 12 fine = 18 total control points
        KRATOS_EXPECT_EQ(p_thb->PointsNumber(), 18);

        // Without a refinement domain, no coarse CP is deactivated
        KRATOS_EXPECT_EQ(p_thb->NumberOfActiveFunctions(0), 6);
        KRATOS_EXPECT_EQ(p_thb->NumberOfActiveFunctions(1), 12);
    }

    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceActiveFunctionsAfterFullRefinement, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();

        KRATOS_EXPECT_EQ(p_thb->NumberOfLevels(), 2);
        KRATOS_EXPECT_EQ(p_thb->NumberOfActiveFunctions(0),  0);  // all deactivated
        KRATOS_EXPECT_EQ(p_thb->NumberOfActiveFunctions(1), 12);  // all fine CPs active
    }

    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceNodeGlobalCoordinatesAfterRefinement, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();

        array_1d<double, 3> parameter(0.0);
        array_1d<double, 3> result(0.0);

        parameter[0] = 10.0;
        parameter[1] = 3.5;
        p_thb->GlobalCoordinates(result, parameter);
        KRATOS_EXPECT_NEAR(result[0], 10.0,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[1],  1.5,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[2], -4.0,  TOLERANCE);

        parameter[0] = 6.0;
        parameter[1] = 1.0;
        p_thb->GlobalCoordinates(result, parameter);
        KRATOS_EXPECT_NEAR(result[0],  6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(result[1],  4.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(result[2], -1.44,  TOLERANCE);

        parameter[0] = 0.0;
        parameter[1] = 1.0;
        p_thb->GlobalCoordinates(result, parameter);
        KRATOS_EXPECT_NEAR(result[0], 0.0,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[1], 4.0,  TOLERANCE);
        KRATOS_EXPECT_NEAR(result[2], 0.0,  TOLERANCE);
    }

    /// ActiveLevelAtPoint returns the finest level whose refinement domain contains (u, v).
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceActiveLevelAtPoint, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();

        // Interior points -> level 1
        KRATOS_EXPECT_EQ(p_thb->ActiveLevelAtPoint(5.0,  2.5), 1);
        KRATOS_EXPECT_EQ(p_thb->ActiveLevelAtPoint(1.0,  1.0), 1);
        KRATOS_EXPECT_EQ(p_thb->ActiveLevelAtPoint(9.0,  4.0), 1);

        // Boundary points also satisfy ">= MinU - eps" / "<= MaxU + eps" -> level 1
        KRATOS_EXPECT_EQ(p_thb->ActiveLevelAtPoint(0.0,  0.0), 1);
        KRATOS_EXPECT_EQ(p_thb->ActiveLevelAtPoint(10.0, 5.0), 1);
    }

    /// SpansLocalSpace
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceSpansLocalSpace, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();

        std::vector<double> spans_u, spans_v;
        p_thb->SpansLocalSpace(spans_u, 0);
        p_thb->SpansLocalSpace(spans_v, 1);

        const std::vector<double> expected_u = {0.0, 5.0, 10.0};
        const std::vector<double> expected_v = {0.0, 2.5, 5.0};

        KRATOS_EXPECT_EQ(spans_u.size(), expected_u.size());
        KRATOS_EXPECT_EQ(spans_v.size(), expected_v.size());
        for (std::size_t i = 0; i < expected_u.size(); ++i)
            KRATOS_EXPECT_NEAR(spans_u[i], expected_u[i], TOLERANCE);
        for (std::size_t i = 0; i < expected_v.size(); ++i)
            KRATOS_EXPECT_NEAR(spans_v[i], expected_v[i], TOLERANCE);
    }

    /// GetActiveCells returns one entry per active knot-span cell.
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceGetActiveCells, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();
        const auto cells = p_thb->GetActiveCells();

        KRATOS_EXPECT_EQ(cells.size(), 4);
    }

    /// Integration after full h-refinement: 4 fine cells x 3 x 2 = 24 Gauss points.
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceCreateIntegrationPointsRefined, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();

        typename Geometry<Node>::IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = p_thb->GetDefaultIntegrationInfo();
        p_thb->CreateIntegrationPoints(integration_points, integration_info);

        KRATOS_EXPECT_EQ(integration_points.size(), 24);
        double area = 0.0;
        for (IndexType i = 0; i < integration_points.size(); ++i)
            area += integration_points[i].Weight();
        KRATOS_EXPECT_NEAR(area, 50.0, TOLERANCE);
    }

    /// partition of unity
    KRATOS_TEST_CASE_IN_SUITE(THBSurfacePartitionOfUnity, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateRefinedNodeTHBSurface();

        // Test representative points covering both fine cells and boundaries
        const std::vector<std::pair<double,double>> test_points = {
            {1.0,  0.5},
            {4.0,  2.0},
            {6.0,  1.0},
            {9.0,  4.0},
            {5.0,  2.5},
        };

        const std::size_t n_pts = test_points.size();
        Matrix uv(n_pts, 2);
        for (std::size_t i = 0; i < n_pts; ++i) {
            uv(i, 0) = test_points[i].first;
            uv(i, 1) = test_points[i].second;
        }

        const Matrix N = p_thb->EvaluateShapeFunctions(uv);

        for (std::size_t i = 0; i < n_pts; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < N.size2(); ++j)
                sum += N(i, j);
            KRATOS_EXPECT_NEAR(sum, 1.0, 1e-10);
        }
    }

    /// ShapeFunctionsValuesAndCPIndices returns the same N values as the
    /// QuadraturePointGeometry created by CreateQuadraturePointGeometries,
    /// and also first and second derivatives when requested.
    KRATOS_TEST_CASE_IN_SUITE(THBSurfaceShapeFunctionsValuesAndCPIndices, KratosCoreNurbsGeometriesFastSuite)
    {
        auto p_thb = GenerateReferenceNodeTHBSurface();

        // Pick an interior parametric point (same as NurbsSurfaceNode test)
        array_1d<double, 3> local_coords(0.0);
        local_coords[0] = 6.0;
        local_coords[1] = 1.0;

        // Build a quadrature-point geometry at this location for reference
        typename THBSurface3D::IntegrationPointsArrayType integration_points;
        integration_points.push_back(typename THBSurface3D::IntegrationPointType(
            local_coords[0], local_coords[1], 1.0));

        IntegrationInfo dummy_info = p_thb->GetDefaultIntegrationInfo();

        typename Geometry<Node>::GeometriesArrayType qps;
        p_thb->CreateQuadraturePointGeometries(qps, 3, integration_points, dummy_info);

        KRATOS_EXPECT_EQ(qps.size(), 1);
        const auto& r_qp_geom = qps[0];

        // Call the THB-specific ShapeFunctionsValuesAndCPIndices
        std::vector<IndexType> cp_indices;
        Vector shape_functions;
        DenseVector<Matrix> derivatives;

        p_thb->ShapeFunctionsValuesAndCPIndices(
            local_coords,
            cp_indices,
            shape_functions,
            2,
            &derivatives);

        // Basic sanity
        KRATOS_EXPECT_GT(cp_indices.size(), 0);
        KRATOS_EXPECT_EQ(shape_functions.size(), cp_indices.size());
        KRATOS_EXPECT_EQ(derivatives.size(), 2);

        // For single-level THB on a single knot span, all 6 CPs are nonzero
        KRATOS_EXPECT_EQ(cp_indices.size(), 6);

        // All 6 node IDs (1-6) must appear in cp_indices
        std::vector<IndexType> sorted_ids = cp_indices;
        std::sort(sorted_ids.begin(), sorted_ids.end());
        for (IndexType id = 1; id <= 6; ++id)
            KRATOS_EXPECT_EQ(sorted_ids[id - 1], id);

        // Compare N values against the QP geometry row
        const Matrix& rN_ref = r_qp_geom.ShapeFunctionsValues();
        KRATOS_EXPECT_EQ(rN_ref.size1(), 1);
        KRATOS_EXPECT_EQ(rN_ref.size2(), shape_functions.size());
        for (IndexType j = 0; j < shape_functions.size(); ++j)
            KRATOS_EXPECT_NEAR(shape_functions[j], rN_ref(0, j), TOLERANCE);

        // Partition of unity for the returned values
        double sum = 0.0;
        for (IndexType j = 0; j < shape_functions.size(); ++j)
            sum += shape_functions[j];
        KRATOS_EXPECT_NEAR(sum, 1.0, TOLERANCE);

        // Compare 1st derivatives against the QP geometry gradient
        const Matrix& rDN_ref = r_qp_geom.ShapeFunctionLocalGradient(0);
        KRATOS_EXPECT_EQ(derivatives[0].size1(), shape_functions.size());
        KRATOS_EXPECT_EQ(derivatives[0].size2(), 2);
        for (IndexType j = 0; j < shape_functions.size(); ++j) {
            KRATOS_EXPECT_NEAR(derivatives[0](j, 0), rDN_ref(j, 0), TOLERANCE);
            KRATOS_EXPECT_NEAR(derivatives[0](j, 1), rDN_ref(j, 1), TOLERANCE);
        }

        // Compare 2nd derivatives against the QP geometry second derivatives
        const Matrix& rD2N_ref = r_qp_geom.ShapeFunctionDerivatives(2, 0);
        KRATOS_EXPECT_EQ(derivatives[1].size1(), shape_functions.size());
        KRATOS_EXPECT_EQ(derivatives[1].size2(), 3);
        for (IndexType j = 0; j < shape_functions.size(); ++j) {
            KRATOS_EXPECT_NEAR(derivatives[1](j, 0), rD2N_ref(j, 0), TOLERANCE);
            KRATOS_EXPECT_NEAR(derivatives[1](j, 1), rD2N_ref(j, 1), TOLERANCE);
            KRATOS_EXPECT_NEAR(derivatives[1](j, 2), rD2N_ref(j, 2), TOLERANCE);
        }
    }

} // namespace Testing.
} // namespace Kratos.
