#include "testing/testing.h"
#include "custom_utilities/iga_mapping_intersection_utilities.h"
#include "geometries/nurbs_surface_geometry.h"

#include <vector>
#include <cmath>

namespace Kratos::Testing {

namespace Utils = Kratos::IgaMappingIntersectionUtilities;
using CoordinatesArrayType = Utils::CoordinatesArrayType;

namespace {

static BrepSurface<PointerVector<Node>, false, PointerVector<Point>>::Pointer
CreateTestBrepSurface_NoTrim(
    ModelPart& rModelPart,
    const IndexType patch_id,
    const double x0,
    const double y0,
    const double z0,
    const double L,
    const IndexType node_id_offset)   // <<< NEW
{
    using NodeType = Node;
    using BrepSurfaceType = BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>;

    // ------------------------------------------------------------------
    // Create a very simple 2D degree-2 NURBS surface (3x3 control points)
    // ------------------------------------------------------------------
    Geometry<NodeType>::PointsArrayType ctrl_pts;

    IndexType id = node_id_offset;

    // 3x3 grid (u-direction, v-direction)
    // coordinates span [x0, x0+L] x [y0, y0+L]
    for (IndexType j = 0; j < 3; ++j) {
        for (IndexType i = 0; i < 3; ++i) {
            const double x = x0 + L * (static_cast<double>(i) / 2.0);
            const double y = y0 + L * (static_cast<double>(j) / 2.0);
            const double z = z0;

            ctrl_pts.push_back(rModelPart.CreateNewNode(id++, x, y, z));
        }
    }

    // Degree
    const int p = 2;
    const int q = 2;

    // Open uniform knot vectors for degree 2 with 3 control points:
    // [0 0 0 1 1 1]
    Vector knot_u(6);
    knot_u[0]=0.0; knot_u[1]=0.0; knot_u[2]=0.0; knot_u[3]=1.0; knot_u[4]=1.0; knot_u[5]=1.0;

    Vector knot_v(6);
    knot_v[0]=0.0; knot_v[1]=0.0; knot_v[2]=0.0; knot_v[3]=1.0; knot_v[4]=1.0; knot_v[5]=1.0;

    // Non-rational
    Vector weights(9);
    for (IndexType k = 0; k < 9; ++k) weights[k] = 1.0;

    auto p_surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(
        ctrl_pts, p, q, knot_u, knot_v, weights);

    // No trimming loops
    typename BrepSurfaceType::BrepCurveOnSurfaceLoopArrayType outer_loops(0);
    typename BrepSurfaceType::BrepCurveOnSurfaceLoopArrayType inner_loops(0);

    auto p_brep = Kratos::make_shared<BrepSurfaceType>(p_surface, outer_loops, inner_loops);

    // IMPORTANT: geometry id must be patch_id so BuildPatchCaches can find it
    p_brep->SetId(patch_id);

    // Register inside modelpart geometries
    rModelPart.AddGeometry(p_brep);

    return p_brep;
}

double SignedArea2D(const std::vector<CoordinatesArrayType>& r_vertices)
{
    const std::size_t n = r_vertices.size();
    double A = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        A += r_vertices[i][0] * r_vertices[j][1]
           - r_vertices[j][0] * r_vertices[i][1];
    }
    return 0.5 * A;
}

void CheckCCW(std::vector<CoordinatesArrayType>& r_vertices)
{
    Utils::SortVerticesCounterClockwise(r_vertices);
    KRATOS_EXPECT_GT(SignedArea2D(r_vertices), 0.0);
}

} // unnamed namespace

KRATOS_TEST_CASE_IN_SUITE(IgaMappingIntersection_SortVerticesCounterClockwise_Triangle, KratosMappingApplicationSerialTestSuite)
{
    std::vector<CoordinatesArrayType> tri(3);
    for (auto& v : tri) v = ZeroVector(3);

    tri[0][0] = 1.0; tri[0][1] = 0.0;
    tri[1][0] = 0.0; tri[1][1] = 0.0;
    tri[2][0] = 0.0; tri[2][1] = 1.0;

    CheckCCW(tri);
}

KRATOS_TEST_CASE_IN_SUITE(IgaMappingIntersection_SortVerticesCounterClockwise_Quad, KratosMappingApplicationSerialTestSuite)
{
    std::vector<CoordinatesArrayType> quad(4);
    for (auto& v : quad) v = ZeroVector(3);

    // square shuffled
    quad[0][0] = 1.0; quad[0][1] = 0.0;
    quad[1][0] = 0.0; quad[1][1] = 1.0;
    quad[2][0] = 0.0; quad[2][1] = 0.0;
    quad[3][0] = 1.0; quad[3][1] = 1.0;

    CheckCCW(quad);
}

KRATOS_TEST_CASE_IN_SUITE(IgaMappingIntersection_BuildPatchCaches_Basic, KratosMappingApplicationSerialTestSuite)
{
    using namespace Kratos;
    using namespace Kratos::IgaMappingIntersectionUtilities;

    Model model;
    ModelPart& r_mp = model.CreateModelPart("iga_domain");
    r_mp.CreateNewProperties(0);

    // --- Create one BrepSurface patch inside the modelpart ---
    const IndexType patch_id = 1;

    // surface origin + size
    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;
    const double L  = 1.0;

    // unique node id range for this patch
    const IndexType node_id_offset = 1;

    auto p_brep = CreateTestBrepSurface_NoTrim(
        r_mp,
        patch_id,
        x0, y0, z0,
        L,
        node_id_offset);

    // BuildPatchCaches needs patch ids
    std::vector<IndexType> patches_id{patch_id};

    const IndexType n_div = 10; // keep test fast
    auto cache = BuildPatchCaches(patches_id, r_mp, n_div);

    KRATOS_EXPECT_EQ(cache.size(), 1u);

    auto it = cache.find(patch_id);
    KRATOS_EXPECT_TRUE(it != cache.end());

    const auto& pc = it->second;

    KRATOS_EXPECT_TRUE(pc.p_bins != nullptr);
    KRATOS_EXPECT_EQ(pc.number_pts, n_div + 1);

    // (n_div+1)^2 points
    KRATOS_EXPECT_EQ(pc.points.size(), (n_div + 1) * (n_div + 1));
}

KRATOS_TEST_CASE_IN_SUITE(
    IgaMappingIntersection_GetPatchesWithProbableProjection_Basic,
    KratosMappingApplicationSerialTestSuite)
{
    using namespace Kratos;
    using namespace Kratos::IgaMappingIntersectionUtilities;

    Model model;
    ModelPart& r_model_part = model.CreateModelPart("iga_model_part");
    r_model_part.CreateNewProperties(0);

    const IndexType patch_id_1 = 1;
    const IndexType patch_id_2 = 2;

    // No gap: patch2 starts at x = 1.0
    const double patch2_offset_x = 1.0;

    IndexType next_node_id = 1;

    auto p_patch_1 = CreateTestBrepSurface_NoTrim(
        r_model_part, patch_id_1, 0.0, 0.0, 0.0, 1.0, next_node_id);
    next_node_id += 9;

    auto p_patch_2 = CreateTestBrepSurface_NoTrim(
        r_model_part, patch_id_2, 1.0, 0.0, 0.0, 1.0, next_node_id);

    std::vector<IndexType> patches_id = {patch_id_1, patch_id_2};

    const IndexType n_div = 20; // a bit finer helps at the shared boundary
    auto patch_cache = BuildPatchCaches(patches_id, r_model_part, n_div);

    KRATOS_WATCH(patch_cache.size())

    KRATOS_EXPECT_EQ(patch_cache.size(), 2);

    // ------------------------------------------------------------
    // Query 1: clearly inside patch 1
    CoordinatesArrayType element_center_1 = ZeroVector(3);
    element_center_1[0] = 0.25;
    element_center_1[1] = 0.50;
    element_center_1[2] = 0.0;

    const double search_radius_close = 0.20;

    auto out_1 = GetPatchesWithProbableProjection(
        patches_id, patch_cache, element_center_1, search_radius_close);

    KRATOS_EXPECT_EQ(out_1.size(), 1);
    KRATOS_EXPECT_EQ(out_1[0], patch_id_1);

    // ------------------------------------------------------------
    // Query 2: clearly inside patch 2
    CoordinatesArrayType element_center_2 = ZeroVector(3);
    element_center_2[0] = 1.75;
    element_center_2[1] = 0.50;
    element_center_2[2] = 0.0;

    auto out_2 = GetPatchesWithProbableProjection(
        patches_id, patch_cache, element_center_2, search_radius_close);

    KRATOS_EXPECT_EQ(out_2.size(), 1);
    KRATOS_EXPECT_EQ(out_2[0], patch_id_2);

    // ------------------------------------------------------------
    // Query 3: exactly on the shared boundary -> should hit BOTH
    CoordinatesArrayType element_center_mid = ZeroVector(3);
    element_center_mid[0] = 1.0;
    element_center_mid[1] = 0.50;
    element_center_mid[2] = 0.0;

    const double search_radius_mid = 0.05;

    auto out_mid = GetPatchesWithProbableProjection(
        patches_id, patch_cache, element_center_mid, search_radius_mid);

    KRATOS_EXPECT_EQ(out_mid.size(), 2);
    KRATOS_EXPECT_EQ(out_mid[0], patch_id_1);
    KRATOS_EXPECT_EQ(out_mid[1], patch_id_2);

    // ------------------------------------------------------------
    // Query 4: far from both -> fallback to full list
    CoordinatesArrayType element_center_far = ZeroVector(3);
    element_center_far[0] = 100.0;
    element_center_far[1] = 100.0;
    element_center_far[2] = 0.0;

    auto out_far = GetPatchesWithProbableProjection(
        patches_id, patch_cache, element_center_far, search_radius_close);

    KRATOS_EXPECT_EQ(out_far.size(), patches_id.size());
    KRATOS_EXPECT_EQ(out_far[0], patches_id[0]);
    KRATOS_EXPECT_EQ(out_far[1], patches_id[1]);
}

KRATOS_TEST_CASE_IN_SUITE(IgaMappingIntersection_FindInitialGuessNewtonRaphsonProjection_Basic, KratosMappingApplicationSerialTestSuite)
{
    using namespace Kratos;
    using namespace Kratos::IgaMappingIntersectionUtilities;

    Model model;
    ModelPart& r_mp = model.CreateModelPart("iga_domain");
    r_mp.CreateNewProperties(0);

    // --- Create BrepSurface patch ---
    const IndexType patch_id = 1;

    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;
    const double L  = 1.0;

    const IndexType node_id_offset = 1;

    auto p_brep = CreateTestBrepSurface_NoTrim(
        r_mp,
        patch_id,
        x0, y0, z0,
        L,
        node_id_offset);

    std::vector<IndexType> patches_id{patch_id};

    // --- Build cache ---
    const IndexType n_div = 10; // small but enough
    auto patch_cache = BuildPatchCaches(patches_id, r_mp, n_div);

    KRATOS_EXPECT_EQ(patch_cache.size(), 1u);

    // --- CASE 1: point close to the patch -> should return true ---
    CoordinatesArrayType slave_xyz = ZeroVector(3);
    slave_xyz[0] = 0.5;
    slave_xyz[1] = 0.5;
    slave_xyz[2] = 0.0;

    CoordinatesArrayType initial_guess = ZeroVector(3);

    const double search_radius = 0.25;

    const bool found = FindInitialGuessNewtonRaphsonProjection(
        slave_xyz,
        p_brep,          // GeometryPointerType
        patch_cache,
        initial_guess,
        search_radius);

    KRATOS_EXPECT_TRUE(found);

    // Check that initial guess is "reasonable"
    KRATOS_EXPECT_GE(initial_guess[0], 0.0);
    KRATOS_EXPECT_LE(initial_guess[0], 1.0);
    KRATOS_EXPECT_GE(initial_guess[1], 0.0);
    KRATOS_EXPECT_LE(initial_guess[1], 1.0);

    // --- CASE 2: far point -> should return false ---
    CoordinatesArrayType slave_xyz_far = ZeroVector(3);
    slave_xyz_far[0] = 100.0;
    slave_xyz_far[1] = 100.0;
    slave_xyz_far[2] = 0.0;

    CoordinatesArrayType initial_guess_far = ZeroVector(3);

    const bool found_far = FindInitialGuessNewtonRaphsonProjection(
        slave_xyz_far,
        p_brep,
        patch_cache,
        initial_guess_far,
        search_radius);

    KRATOS_EXPECT_FALSE(found_far);
}

KRATOS_TEST_CASE_IN_SUITE(IgaMappingIntersection_AreProjectionsOnParameterSpaceBoundary_Basic, KratosMappingApplicationSerialTestSuite)
{
    using namespace Kratos;
    using namespace Kratos::IgaMappingIntersectionUtilities;

    Model model;
    ModelPart& r_mp = model.CreateModelPart("iga_domain");
    r_mp.CreateNewProperties(0);

    // --- Create BrepSurface patch ---
    const IndexType patch_id = 1;

    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;
    const double L  = 1.0;

    const IndexType node_id_offset = 1;

    auto p_brep = CreateTestBrepSurface_NoTrim(
        r_mp,
        patch_id,
        x0, y0, z0,
        L,
        node_id_offset);

    // Get BACKGROUND nurbs surface from BrepSurface
    auto p_background = p_brep->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

    auto p_nurbs = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(
        p_background);

    KRATOS_EXPECT_TRUE(p_nurbs != nullptr);

    const auto& r_nurbs = *p_nurbs;

    // Read intervals to construct boundary points robustly
    const double u_min = r_nurbs.DomainIntervalU().MinParameter();
    const double u_max = r_nurbs.DomainIntervalU().MaxParameter();
    const double v_min = r_nurbs.DomainIntervalV().MinParameter();
    const double v_max = r_nurbs.DomainIntervalV().MaxParameter();

    // ============================================================
    // CASE 0: empty -> false
    // ============================================================
    {
        std::vector<CoordinatesArrayType> points;
        const bool is_boundary = AreProjectionsOnParameterSpaceBoundary(points, r_nurbs);
        KRATOS_EXPECT_FALSE(is_boundary);
    }

    // ============================================================
    // CASE 1: one point strictly inside -> false
    // ============================================================
    {
        std::vector<CoordinatesArrayType> points(1);
        points[0] = ZeroVector(3);
        points[0][0] = 0.5 * (u_min + u_max);
        points[0][1] = 0.5 * (v_min + v_max);

        const bool is_boundary = AreProjectionsOnParameterSpaceBoundary(points, r_nurbs);
        KRATOS_EXPECT_FALSE(is_boundary);
    }

    // ============================================================
    // CASE 2: one point on boundary (u=u_min) -> true
    // ============================================================
    {
        std::vector<CoordinatesArrayType> points(1);
        points[0] = ZeroVector(3);
        points[0][0] = u_min;                      // boundary
        points[0][1] = 0.5 * (v_min + v_max);

        const bool is_boundary = AreProjectionsOnParameterSpaceBoundary(points, r_nurbs);
        KRATOS_EXPECT_TRUE(is_boundary);
    }

    // ============================================================
    // CASE 3: two points on same boundary edge -> true
    // (both on v=v_max)
    // ============================================================
    {
        std::vector<CoordinatesArrayType> points(2);

        points[0] = ZeroVector(3);
        points[1] = ZeroVector(3);

        points[0][0] = 0.25 * (u_min + 3.0*u_max);
        points[1][0] = 0.75 * (u_min + 1.0*u_max);

        points[0][1] = v_max;
        points[1][1] = v_max;

        const bool is_boundary = AreProjectionsOnParameterSpaceBoundary(points, r_nurbs);
        KRATOS_EXPECT_TRUE(is_boundary);
    }

    // ============================================================
    // CASE 4: two points on different boundaries -> false
    // (one on u_min, one on u_max) => not "same edge"
    // ============================================================
    {
        std::vector<CoordinatesArrayType> points(2);

        points[0] = ZeroVector(3);
        points[1] = ZeroVector(3);

        points[0][0] = u_min;
        points[0][1] = 0.25 * (v_min + 3.0*v_max);

        points[1][0] = u_max;
        points[1][1] = 0.75 * (v_min + 1.0*v_max);

        const bool is_boundary = AreProjectionsOnParameterSpaceBoundary(points, r_nurbs);
        KRATOS_EXPECT_FALSE(is_boundary);
    }

    // ============================================================
    // CASE: p0 on u_min, p1 on v_min -> false (different edges)
    // ============================================================
    {
        std::vector<CoordinatesArrayType> points(2);

        points[0] = ZeroVector(3);
        points[1] = ZeroVector(3);

        // p0 on u = u_min
        points[0][0] = u_min;
        points[0][1] = 0.75 * (v_min + 1.0*v_max); // not on v boundary

        // p1 on v = v_min
        points[1][0] = 0.75 * (u_min + 1.0*u_max); // not on u boundary
        points[1][1] = v_min;

        const bool is_boundary = AreProjectionsOnParameterSpaceBoundary(points, r_nurbs);
        KRATOS_EXPECT_FALSE(is_boundary);
    }

}

KRATOS_TEST_CASE_IN_SUITE(
    IgaMappingIntersection_FindTriangleSegmentSurfaceIntersectionWithBisection_Basic,
    KratosMappingApplicationSerialTestSuite)
{
    using namespace Kratos;
    using namespace Kratos::IgaMappingIntersectionUtilities;

    Model model;
    ModelPart& r_mp = model.CreateModelPart("iga_domain");
    r_mp.CreateNewProperties(0);

    // --- Create one BrepSurface patch and register it ---
    const IndexType patch_id = 1;

    // IMPORTANT: use the same helper you already have
    auto p_brep = CreateTestBrepSurface_NoTrim(r_mp, patch_id, 0.0, 0.0, 0.0, 1.0, 1);
    r_mp.AddGeometry(p_brep);

    GeometryPointerType p_geom_master = p_brep;

    // --- Get parametric bounds (u_min,u_max,v_min,v_max) ---
    auto p_background = p_brep->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

    std::vector<double> knot_u, knot_v;
    p_background->SpansLocalSpace(knot_u, 0);
    p_background->SpansLocalSpace(knot_v, 1);

    KRATOS_ERROR_IF(knot_u.empty() || knot_v.empty())
        << "Empty knot vectors in test surface.\n";

    const double u_min = knot_u.front();
    const double u_max = knot_u.back();
    const double v_min = knot_v.front();
    const double v_max = knot_v.back();

    const double u_mid = 0.5 * (u_min + u_max);
    const double v_mid = 0.5 * (v_min + v_max);

    // --- Define inside point using GlobalCoordinates at (u_mid, v_mid) ---
    CoordinatesArrayType local_inside = ZeroVector(3);
    local_inside[0] = u_mid;
    local_inside[1] = v_mid;

    CoordinatesArrayType point_inside_xyz = ZeroVector(3);
    p_brep->GlobalCoordinates(point_inside_xyz, local_inside);

    // --- Define outside point by moving in +X direction (outside patch domain) ---
    CoordinatesArrayType point_outside_xyz = point_inside_xyz;
    point_outside_xyz[0] += 2.0; // move far outside in x

    // --- Provide an initial guess in local space ---
    CoordinatesArrayType initial_guess = ZeroVector(3);
    initial_guess[0] = u_mid;
    initial_guess[1] = v_mid;

    // --- Output of function ---
    CoordinatesArrayType intersection_local = ZeroVector(3);

    const bool found =
        FindTriangleSegmentSurfaceIntersectionWithBisection(
            p_geom_master,
            point_inside_xyz,
            point_outside_xyz,
            initial_guess,
            intersection_local);

    KRATOS_EXPECT_TRUE(found);

    // Expect intersection to be close to u_max boundary
    KRATOS_EXPECT_NEAR(intersection_local[0], u_max, 1e-4);

    // and keep v close to v_mid (since we moved only in X direction)
    KRATOS_EXPECT_NEAR(intersection_local[1], v_mid, 1e-3);

    // Extra verification: physical intersection is on the boundary 
    CoordinatesArrayType intersection_xyz = ZeroVector(3);
    p_brep->GlobalCoordinates(intersection_xyz, intersection_local);

    // The boundary physical point should have x ≈ max x on surface
    // (since u=u_max means right edge for your flat surface)
    // We also check it's near z-plane of the surface.
    KRATOS_EXPECT_NEAR(intersection_xyz[2], point_inside_xyz[2], 1e-8);
}

} // namespace Kratos::Testing
