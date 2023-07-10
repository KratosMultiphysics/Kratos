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
#include "testing/testing.h"
#include "custom_utilities/projection_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos::Testing {

typedef Node NodeType;
typedef Geometry<NodeType> GeometryType;

namespace {

template<std::size_t TSize>
void SetEqIdsOnNodes(GeometryType& rGeometry,
                     const std::array<int, TSize>& rExpEquationIds)
{
    KRATOS_CHECK_EQUAL(TSize, rGeometry.PointsNumber());

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

    KRATOS_CHECK_EQUAL(FullProjection, is_full_projection);
    KRATOS_CHECK_EQUAL(ExpPairingIndex, pairing_index);

    if (pairing_index != ProjectionUtilities::PairingIndex::Unspecified) {
        KRATOS_CHECK_DOUBLE_EQUAL(ExpProjectionDistance, proj_dist);

        KRATOS_CHECK_EQUAL(rExpSFValues.size(), sf_values.size());
        KRATOS_CHECK_EQUAL(rExpEquationIds.size(), rExpEquationIds.size());

        for (std::size_t i=0; i<TSize; ++i) {
            KRATOS_CHECK_NEAR(rExpSFValues[i], sf_values[i], 1E-13);
            KRATOS_CHECK_EQUAL(rExpEquationIds[i], eq_ids[i]);
        }
    }
}

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

}  // namespace Kratos::Testing