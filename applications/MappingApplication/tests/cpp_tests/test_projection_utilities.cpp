//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/brep_surface.h"
#include "testing/testing.h"
#include "custom_utilities/projection_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos::Testing {

typedef Node NodeType;
typedef Geometry<NodeType> GeometryType;
typedef Kratos::QuadraturePointGeometry<NodeType, 3, 2> QuadraturePointGeometryType;
typedef QuadraturePointGeometryType::Pointer QuadraturePointGeometryPointerType;

namespace {

template<std::size_t TSize>
void SetEqIdsOnNodes(GeometryType& rGeometry,
                     const std::array<int, TSize>& rExpEquationIds)
{
    KRATOS_EXPECT_EQ(TSize, rGeometry.PointsNumber());

    std::size_t node_idx = 0;
    for (auto& r_node : rGeometry.Points()) {
        r_node.SetValue(INTERFACE_EQUATION_ID, rExpEquationIds[node_idx++]);
    }
}

template<std::size_t TSize>
void TestComputeProjection(const GeometryType& rGeometry,
                           const Point& rPointToProject,
                           const double LocalCoordTol,
                           const std::array<double, TSize>& rExpSFValues,
                           const std::array<int, TSize>& rExpEquationIds,
                           const double ExpProjectionDistance,
                           const ProjectionUtilities::PairingIndex ExpPairingIndex,
                           const bool ComputeApproximation,
                           const bool FullProjection)
{
    Vector sf_values;
    std::vector<int> eq_ids;
    double proj_dist;
    ProjectionUtilities::PairingIndex pairing_index;

    const bool is_full_projection = ProjectionUtilities::ComputeProjection(rGeometry, rPointToProject, LocalCoordTol, sf_values, eq_ids, proj_dist, pairing_index, ComputeApproximation);

    KRATOS_EXPECT_EQ(FullProjection, is_full_projection);
    KRATOS_EXPECT_EQ(ExpPairingIndex, pairing_index);

    if (pairing_index != ProjectionUtilities::PairingIndex::Unspecified) {
        KRATOS_EXPECT_NEAR(ExpProjectionDistance, proj_dist, 1e-6);

        KRATOS_EXPECT_EQ(rExpSFValues.size(), sf_values.size());
        KRATOS_EXPECT_EQ(rExpEquationIds.size(), rExpEquationIds.size());

        for (std::size_t i=0; i<TSize; ++i) {
            KRATOS_EXPECT_NEAR(rExpSFValues[i], sf_values[i], 1E-6);
            KRATOS_EXPECT_EQ(rExpEquationIds[i], eq_ids[i]);
        }
    }
}

}

GeometryType::Pointer GenerateBSplineCurveOnBSplineSurface2d()
    {

    // Assign the points belonging to the curve
    PointerVector<NodeType> control_points_curve;
    control_points_curve.push_back(Kratos::make_intrusive<NodeType>(1, 1.4142135623730951, 3, 0));
    control_points_curve.push_back(Kratos::make_intrusive<NodeType>(2, 0, 3, 0));

    // Assign the curve's knot vector
    Vector knot_vector_curve = ZeroVector(2);
    knot_vector_curve[0] = 0.0;
    knot_vector_curve[1] = 1.4142135623730951;

    // Polynomial degree of the curve
    int p_curve = 1;

    // Create the 2D embedded curve
    auto curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<NodeType>>>(control_points_curve, p_curve, knot_vector_curve);

    // Assign the points belonging to the surface
    PointerVector<NodeType> control_points_surface;
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(1, 0, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(2, 0, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(3, 0, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(4, 0, 1, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(5, 0.5, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(6, 0.5625, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(7, 0.6875, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(8, 0.75, 1, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(9, 1.5, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(10, 1.6875, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(11, 2.0625, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(12, 2.25, 1, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(13, 2, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(14, 2.25, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(15, 2.75, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(16, 3, 1, 0));


    // Assign the surface's knot vectors
    Vector knot_vector_u_surface = ZeroVector(5);
    knot_vector_u_surface[0] = 0.0;
    knot_vector_u_surface[1] = 0.0;
    knot_vector_u_surface[2] = 0.70710678118654757;
    knot_vector_u_surface[3] = 1.4142135623730951;
    knot_vector_u_surface[4] = 1.4142135623730951;

    Vector knot_vector_v_surface = ZeroVector(5);
    knot_vector_v_surface[0] = 0.0;
    knot_vector_v_surface[1] = 0.0;
    knot_vector_v_surface[2] = 1.5;
    knot_vector_v_surface[3] = 3.0;
    knot_vector_v_surface[4] = 3.0;

    // Polynomial degrees
    int p_surface = 2;
    int q_surface = 2;

    // Create a 2D surface
    auto surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(control_points_surface, p_surface,
            q_surface, knot_vector_u_surface, knot_vector_v_surface);

    // Create and return a curve on surface geometry
    GeometryType::Pointer nurbs_curve_on_surface = Kratos::make_shared<NurbsCurveOnSurfaceGeometry<3, PointerVector<NodeType>, PointerVector<NodeType>>>(surface, curve);
    return nurbs_curve_on_surface;
}

GeometryType::Pointer GenerateBrepSplineSurface2d()
{
    // Assign the points belonging to the surface
    PointerVector<NodeType> control_points_surface;
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(1, 0, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(2, 0, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(3, 0, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(4, 0, 1, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(5, 0.5, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(6, 0.5625, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(7, 0.6875, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(8, 0.75, 1, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(9, 1.5, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(10, 1.6875, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(11, 2.0625, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(12, 2.25, 1, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(13, 2, 0, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(14, 2.25, 0.25, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(15, 2.75, 0.75, 0));
    control_points_surface.push_back(Kratos::make_intrusive<NodeType>(16, 3, 1, 0));



    // Assign the surface's knot vectors
    Vector knot_vector_u_surface = ZeroVector(5);
    knot_vector_u_surface[0] = 0.0;
    knot_vector_u_surface[1] = 0.0;
    knot_vector_u_surface[2] = 0.70710678118654757;
    knot_vector_u_surface[3] = 1.4142135623730951;
    knot_vector_u_surface[4] = 1.4142135623730951;

    Vector knot_vector_v_surface = ZeroVector(5);
    knot_vector_v_surface[0] = 0.0;
    knot_vector_v_surface[1] = 0.0;
    knot_vector_v_surface[2] = 1.5;
    knot_vector_v_surface[3] = 3.0;
    knot_vector_v_surface[4] = 3.0;

    // Polynomial degrees
    int p_surface = 2;
    int q_surface = 2;

    // Create a 2D surface
    NurbsSurfaceGeometry<3, PointerVector<NodeType>>::Pointer nurbs_surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(control_points_surface, p_surface,
            q_surface, knot_vector_u_surface, knot_vector_v_surface);

    auto brep_surface = Kratos::make_shared<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(nurbs_surface);

    return brep_surface;
}

GeometryType::GeometriesArrayType CreateQuadraturePointsGeometries(
    GeometryType::Pointer pGeometry)
{
    GeometryType::GeometriesArrayType result_geometries; 

    // Get integration points and info
    Geometry<Node>::IntegrationPointsArrayType integration_points;
    IntegrationInfo integration_info = pGeometry->GetDefaultIntegrationInfo();
    pGeometry->CreateIntegrationPoints(integration_points, integration_info);

    // Resize the result vector to match integration points
    result_geometries.resize(integration_points.size());

    // Call the method from your implementation
    const std::size_t num_shape_derivatives = 6; // usually N and dN/dξ
    pGeometry->CreateQuadraturePointGeometries(
        result_geometries, num_shape_derivatives, integration_points, integration_info);

    return result_geometries;
}



KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Line_Inside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    const GeometryType::Pointer p_geom(Kratos::make_shared<Line2D2<NodeType>>(node_1, node_2));

    double proj_dist = 0.2;
    const Point point_to_project(0.25, proj_dist, 0.0);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Line_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<double, 2> exp_sf_values {0.75, 0.25};
    const std::array<int, 2> exp_eq_ids {35, 18};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Line_Outside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    const GeometryType::Pointer p_geom(Kratos::make_shared<Line2D2<NodeType>>(node_1, node_2));

    double proj_dist = 0.2;
    const Point point_to_project(-0.1, proj_dist, 0.0);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Line_Outside;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 2> exp_sf_values {1.1, -0.1};
    const std::array<int, 2> exp_eq_ids {35, 18};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Line_Closest_Point, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    const GeometryType::Pointer p_geom(Kratos::make_shared<Line2D2<NodeType>>(node_1, node_2));

    const Point point_to_project(-0.35, 0.2, 0.0);
    const double local_coord_tol = 0.0;
    double proj_dist = std::sqrt(std::pow(point_to_project[0],2) + std::pow(point_to_project[1],2) + std::pow(point_to_project[2],2)); // closest point is (0,0,0)
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Closest_Point;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 1> exp_sf_values {1.0};
    const std::array<int, 1> exp_eq_ids {35};
    const std::array<int, 2> all_eq_ids {35, 18};

    SetEqIdsOnNodes(*p_geom, all_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Line_Unspecified, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    const GeometryType::Pointer p_geom(Kratos::make_shared<Line2D2<NodeType>>(node_1, node_2));

    double proj_dist = 0.2;
    const Point point_to_project(-0.25, proj_dist, 0.0);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Unspecified;
    const bool compute_approximation = false;
    const bool full_projection = false;

    const std::array<double, 2> exp_sf_values {};
    const std::array<int, 2> exp_eq_ids {};

    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Triangle_Inside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    double proj_dist = 0.35;
    const Point point_to_project(0.5, 0.3, proj_dist);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Surface_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<double, 3> exp_sf_values {0.5, 0.2, 0.3};
    const std::array<int, 3> exp_eq_ids {35, 18, 108};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Triangle_Outside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    double proj_dist = 0.35;
    const Point point_to_project(1.1, 0.1, proj_dist);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Surface_Outside;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 3> exp_sf_values {-0.1, 1.0, 0.1};
    const std::array<int, 3> exp_eq_ids {35, 18, 108};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Triangle_Line, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    const Point point_to_project(1.1, 0.1, 0.35);
    const double local_coord_tol = 0.0;
    double proj_dist = std::sqrt(std::pow(0.1, 2) + std::pow(0.35,2));
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Line_Inside;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 2> exp_sf_values {0.9, 0.1};
    const std::array<int, 2> exp_eq_ids {18, 108};
    const std::array<int, 3> all_eq_ids {35, 18, 108};

    SetEqIdsOnNodes(*p_geom, all_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Triangle_Closest_Point, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    const Point point_to_project(1.1, -1.2, 0.0);
    const double local_coord_tol = 0.0;
    double proj_dist = std::sqrt(std::pow(0.1, 2) + std::pow(point_to_project[1], 2));
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Closest_Point;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 1> exp_sf_values {1.0};
    const std::array<int, 1> exp_eq_ids {18};
    const std::array<int, 3> all_eq_ids {35, 18, 108};

    SetEqIdsOnNodes(*p_geom, all_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Triangle_Unspecified, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Triangle3D3<NodeType>>(node_1, node_2, node_3));

    const Point point_to_project(1.1, -0.1, 0.0);
    const double local_coord_tol = 0.0;
    double proj_dist = std::sqrt(std::pow(0.1, 2) + std::pow(point_to_project[1], 2));
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Unspecified;
    const bool compute_approximation = false;
    const bool full_projection = false;

    const std::array<double, 1> exp_sf_values {};
    const std::array<int, 1> exp_eq_ids {};

    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Quadrilateral_Inside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Quadrilateral3D4<NodeType>>(node_1, node_2, node_3, node_4));

    double proj_dist = 0.35;
    const Point point_to_project(0.5, 0.3, proj_dist);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Surface_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<double, 4> exp_sf_values {0.35, 0.35, 0.15, 0.15};
    const std::array<int, 4> exp_eq_ids {35, 18, 108, 95};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Quadrilateral_Outside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, 1.0, 0.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Quadrilateral3D4<NodeType>>(node_1, node_2, node_3, node_4));

    double proj_dist = 0.35;
    const Point point_to_project(-0.1, -0.1, proj_dist);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Surface_Outside;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 4> exp_sf_values {1.21, -0.11, 0.01, -0.11};
    const std::array<int, 4> exp_eq_ids {35, 18, 108, 95};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Tetrahedra_Inside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.5, 1.0, 1.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Tetrahedra3D4<NodeType>>(node_1, node_2, node_3, node_4));

    const Point point_to_project(0.5, 0.3, 0.2);
    const double local_coord_tol = 0.2;
    double proj_dist = 1.4465476141489435058;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Volume_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<double, 4> exp_sf_values {0.4, 0.3, 0.1, 0.2};
    const std::array<int, 4> exp_eq_ids {35, 18, 108, 95};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Hexahedra_Inside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, 1.0, 0.0));
    auto node_5(Kratos::make_intrusive<NodeType>(5, 0.0, 0.0, 1.0));
    auto node_6(Kratos::make_intrusive<NodeType>(6, 1.0, 0.0, 1.0));
    auto node_7(Kratos::make_intrusive<NodeType>(7, 1.0, 1.0, 1.0));
    auto node_8(Kratos::make_intrusive<NodeType>(8, 0.0, 1.0, 1.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Hexahedra3D8<NodeType>>(node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8));

    const Point point_to_project(0.5, 0.3, 0.2);
    const double local_coord_tol = 0.2;
    double proj_dist = 0.36055512754639901241;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Volume_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<double, 8> exp_sf_values {0.28, 0.28, 0.12, 0.12, 0.07, 0.07, 0.03, 0.03};
    const std::array<int, 8> exp_eq_ids {35, 18, 108, 95, 12, 14, 19, 22};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Hexahedra_Outside, KratosMappingApplicationSerialTestSuite)
{
    auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
    auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
    auto node_3(Kratos::make_intrusive<NodeType>(3, 1.0, 1.0, 0.0));
    auto node_4(Kratos::make_intrusive<NodeType>(4, 0.0, 1.0, 0.0));
    auto node_5(Kratos::make_intrusive<NodeType>(5, 0.0, 0.0, 1.0));
    auto node_6(Kratos::make_intrusive<NodeType>(6, 1.0, 0.0, 1.0));
    auto node_7(Kratos::make_intrusive<NodeType>(7, 1.0, 1.0, 1.0));
    auto node_8(Kratos::make_intrusive<NodeType>(8, 0.0, 1.0, 1.0));

    const GeometryType::Pointer p_geom(Kratos::make_shared<Hexahedra3D8<NodeType>>(node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8));

    const Point point_to_project(0.5, 0.5, -0.1);
    const double local_coord_tol = 0.2;
    double proj_dist = 0.6;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Volume_Outside;
    const bool compute_approximation = true;
    const bool full_projection = false;

    const std::array<double, 8> exp_sf_values {0.275, 0.275, 0.275, 0.275, -0.025, -0.025, -0.025, -0.025};
    const std::array<int, 8> exp_eq_ids {35, 18, 108, 95, 12, 14, 19, 22};

    SetEqIdsOnNodes(*p_geom, exp_eq_ids);
    TestComputeProjection(*p_geom, point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Nurbs_Curve_On_Surface_Inside, KratosMappingApplicationSerialTestSuite)
{
    const auto nurbs_curve_on_surface = GenerateBSplineCurveOnBSplineSurface2d();
    GeometryType::GeometriesArrayType qp_geometries = CreateQuadraturePointsGeometries(nurbs_curve_on_surface);

    const Point point_to_project(2.0, 0.0, 0.0);
    const double local_coord_tol = 1e-6;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Line_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<int, 16> set_eq_ids {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    const std::array<double, 9> exp_sf_values {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    const std::array<int, 9> exp_eq_ids {4, 5, 6, 8, 9, 10, 12, 13, 14};

    GeometryType::Pointer surface = nurbs_curve_on_surface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
    SetEqIdsOnNodes(*surface, set_eq_ids);
    
    Vector sf_values;
    sf_values.clear();
    std::vector<int> eq_ids;
    double proj_dist = 0.0;

    qp_geometries(0)->SetGeometryParent(nurbs_curve_on_surface.get());
    TestComputeProjection(*qp_geometries(0), point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}


KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Nurbs_Curve_On_Surface_Outside, KratosMappingApplicationSerialTestSuite)
{
    const auto nurbs_curve_on_surface = GenerateBSplineCurveOnBSplineSurface2d();
    GeometryType::GeometriesArrayType qp_geometries = CreateQuadraturePointsGeometries(nurbs_curve_on_surface);

    const Point point_to_project(2.1, -0.1, 0.0);
    const double local_coord_tol = 0.2;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Line_Inside;
    const bool compute_approximation = true;
    const bool full_projection = true;

    const std::array<int, 16> set_eq_ids {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    const std::array<double, 9> exp_sf_values {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    const std::array<int, 9> exp_eq_ids {4, 5, 6, 8, 9, 10, 12, 13, 14};

    GeometryType::Pointer surface = nurbs_curve_on_surface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
    SetEqIdsOnNodes(*surface, set_eq_ids);
    
    Vector sf_values;
    sf_values.clear();
    std::vector<int> eq_ids;
    double proj_dist = std::pow(0.02, 0.5);

    qp_geometries(0)->SetGeometryParent(nurbs_curve_on_surface.get());
    TestComputeProjection(*qp_geometries(0), point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

KRATOS_TEST_CASE_IN_SUITE(ProjectionUtils_Nurbs_Surface_Inside, KratosMappingApplicationSerialTestSuite)
{
    const auto brep_surface = GenerateBrepSplineSurface2d();
    GeometryType::GeometriesArrayType qp_geometries = CreateQuadraturePointsGeometries(brep_surface);

    const Point point_to_project(1.0, 0.5, 0.0);
    const double local_coord_tol = 1e-6;
    ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Surface_Inside;
    const bool compute_approximation = false;
    const bool full_projection = true;

    const std::array<int, 16> set_eq_ids {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    const std::array<double, 9> exp_sf_values {0.02, 0.02, 0, 0.32, 0.32, 0, 0.16, 0.16, 0};
    const std::array<int, 9> exp_eq_ids {1, 2, 3, 5, 6, 7, 9, 10, 11};

    GeometryType::Pointer surface = brep_surface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
    SetEqIdsOnNodes(*surface, set_eq_ids);
    
    Vector sf_values;
    sf_values.clear();
    std::vector<int> eq_ids;
    double proj_dist = 0.0;

    qp_geometries(0)->SetGeometryParent(brep_surface.get());
    TestComputeProjection(*qp_geometries(0), point_to_project, local_coord_tol, exp_sf_values, exp_eq_ids, proj_dist, pairing_index, compute_approximation, full_projection);
}

}  // namespace Kratos::Testing