// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Anne van de Graaf
//

#include "custom_geometries/interface_geometry.h"
#include "geometries/geometry_data.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

auto CreateTwoPlusTwoNoded2DLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 7.0, 0.2, 0.0));
    return InterfaceGeometry<Line2D2<Node>>{1, nodes};
}

auto CreateThreePlusThreeNoded2DLineInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 2.5, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, -1.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 7.0, 0.2, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 3.5, 0.4, 0.0));
    return InterfaceGeometry<Line2D3<Node>>{1, nodes};
}

auto CreateThreePlusThreeNoded3DSurfaceInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 1.0, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 1.0, 1.0, 0.5));
    return InterfaceGeometry<Triangle3D3<Node>>{1, nodes};
}

PointerVector<Node> CreateNodesForSixPlusSixNoded3DSurfaceInterface()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(5, 0.5, 0.5, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(6, 0.0, 0.5, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(7, 0.0, 0.0, 0.5));
    result.push_back(Kratos::make_intrusive<Node>(8, 1.0, 0.0, 0.5));
    result.push_back(Kratos::make_intrusive<Node>(9, 0.0, 1.0, 0.5));
    result.push_back(Kratos::make_intrusive<Node>(10, 0.5, 0.0, 0.5));
    result.push_back(Kratos::make_intrusive<Node>(11, 0.5, 0.5, 0.5));
    result.push_back(Kratos::make_intrusive<Node>(12, 0.0, 0.5, 0.5));
    return result;
}

auto CreateSixPlusSixNoded3DSurfaceInterfaceGeometry()
{
    return InterfaceGeometry<Triangle3D6<Node>>{1, CreateNodesForSixPlusSixNoded3DSurfaceInterface()};
}

auto CreateFourPlusFourNoded3DSurfaceInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 0.0, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 1.0, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(7, 1.0, 1.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(8, 0.0, 1.0, 0.5));
    return InterfaceGeometry<Quadrilateral3D4<Node>>{1, nodes};
}

auto CreateEightPlusEightNoded3DSurfaceInterfaceGeometry()
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(5, 0.5, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(6, 1.0, 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(7, 0.5, 1.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(8, 0.0, 0.5, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(9, 0.0, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(10, 1.0, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(11, 1.0, 1.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(12, 0.0, 1.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(13, 0.5, 0.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(14, 1.0, 0.5, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(15, 0.5, 1.0, 0.5));
    nodes.push_back(Kratos::make_intrusive<Node>(16, 0.0, 0.5, 0.5));
    return InterfaceGeometry<Quadrilateral3D8<Node>>{1, nodes};
}

void AssertNodeIdsOfGeometry(const Geometry<Node>& rGeometry, const std::vector<std::size_t>& rExpectedNodeIds)
{
    auto node_ids = std::vector<std::size_t>{};
    node_ids.reserve(rGeometry.size());
    std::transform(rGeometry.begin(), rGeometry.end(), std::back_inserter(node_ids),
                   [](const auto& rNode) { return rNode.Id(); });
    KRATOS_EXPECT_VECTOR_EQ(node_ids, rExpectedNodeIds)
}

} // namespace

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryIsAGeometry, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry        = InterfaceGeometry<Line2D3<Node>>();
    const auto p_base_geometry = dynamic_cast<const Geometry<Node>*>(&geometry);

    KRATOS_EXPECT_NE(p_base_geometry, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryCanBeConstructedGivenASetOfNullPointersToNodes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // The following constructor input data resembles what is done at element registration time
    const auto six_null_pointers_to_nodes = Geometry<Node>::PointsArrayType{6};

    const auto geometry = InterfaceGeometry<Line2D3<Node>>{six_null_pointers_to_nodes};

    KRATOS_EXPECT_EQ(geometry.PointsNumber(), 6);
    KRATOS_EXPECT_EQ(geometry.LocalSpaceDimension(), 1);
    KRATOS_EXPECT_EQ(geometry.WorkingSpaceDimension(), 2);
    KRATOS_EXPECT_EQ(geometry.GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(geometry.GetGeometryOrderType(), GeometryData::KratosGeometryOrderType::Kratos_Quadratic_Order);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryCanBeConstructedGivenASetOfNullPointersToNodesAndASurfaceGeometryType,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // The following constructor input data resembles what is done at element registration time
    const auto six_null_pointers_to_nodes = Geometry<Node>::PointsArrayType{6};

    const auto geometry = InterfaceGeometry<Triangle3D3<Node>>{six_null_pointers_to_nodes};

    KRATOS_EXPECT_EQ(geometry.PointsNumber(), 6);
    KRATOS_EXPECT_EQ(geometry.LocalSpaceDimension(), 2);
    KRATOS_EXPECT_EQ(geometry.WorkingSpaceDimension(), 3);
    KRATOS_EXPECT_EQ(geometry.GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Triangle);
    KRATOS_EXPECT_EQ(geometry.GetGeometryOrderType(), GeometryData::KratosGeometryOrderType::Kratos_Linear_Order);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Create_CreatesNewInstanceOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = InterfaceGeometry<Line2D2<Node>>();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));

    const auto p_new_geometry = geometry.Create(nodes);

    KRATOS_EXPECT_NE(p_new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const InterfaceGeometry<Line2D2<Node>>*>(p_new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(p_new_geometry->PointsNumber(), 4);
    KRATOS_EXPECT_EQ(p_new_geometry->Id(), 0);
    KRATOS_EXPECT_EQ(p_new_geometry->LocalSpaceDimension(), 1);
    KRATOS_EXPECT_EQ(p_new_geometry->WorkingSpaceDimension(), 2);
    KRATOS_EXPECT_EQ(p_new_geometry->GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(p_new_geometry->GetGeometryOrderType(),
                     GeometryData::KratosGeometryOrderType::Kratos_Linear_Order);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_CreateWithId_CreatesNewInstanceOfCorrectTypeAndIdForLinearLineMidGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          geometry = InterfaceGeometry<Line2D2<Node>>();
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(4, 0.0, 0.0, 0.0));

    constexpr auto new_geometry_id = 1;
    const auto     p_new_geometry  = geometry.Create(new_geometry_id, nodes);

    KRATOS_EXPECT_NE(p_new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const InterfaceGeometry<Line2D2<Node>>*>(p_new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(p_new_geometry->PointsNumber(), 4);
    KRATOS_EXPECT_EQ(p_new_geometry->Id(), new_geometry_id);
    KRATOS_EXPECT_EQ(p_new_geometry->LocalSpaceDimension(), 1);
    KRATOS_EXPECT_EQ(p_new_geometry->WorkingSpaceDimension(), 2);
    KRATOS_EXPECT_EQ(p_new_geometry->GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(p_new_geometry->GetGeometryOrderType(),
                     GeometryData::KratosGeometryOrderType::Kratos_Linear_Order);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_CreateWithId_CreatesNewInstanceOfCorrectTypeAndIdForQuadraticTriangularMidGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = InterfaceGeometry<Triangle3D6<Node>>();

    constexpr auto new_geometry_id = 1;
    const auto     p_new_geometry =
        geometry.Create(new_geometry_id, CreateNodesForSixPlusSixNoded3DSurfaceInterface());

    KRATOS_EXPECT_NE(p_new_geometry, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<const InterfaceGeometry<Triangle3D6<Node>>*>(p_new_geometry.get()), nullptr);
    KRATOS_EXPECT_EQ(p_new_geometry->PointsNumber(), 12);
    KRATOS_EXPECT_EQ(p_new_geometry->Id(), new_geometry_id);
    KRATOS_EXPECT_EQ(p_new_geometry->LocalSpaceDimension(), 2);
    KRATOS_EXPECT_EQ(p_new_geometry->WorkingSpaceDimension(), 3);
    KRATOS_EXPECT_EQ(p_new_geometry->GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Triangle);
    KRATOS_EXPECT_EQ(p_new_geometry->GetGeometryOrderType(),
                     GeometryData::KratosGeometryOrderType::Kratos_Quadratic_Order);
}

KRATOS_TEST_CASE_IN_SUITE(CreatingInterfaceWithThreeNodesThrows, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    PointerVector<Node> nodes;
    nodes.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(2, 5.0, 0.0, 0.0));
    nodes.push_back(Kratos::make_intrusive<Node>(3, -1.0, 0.0, 0.0));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(InterfaceGeometry<Line2D3<Node>>{nodes},
                                      "Number of nodes must be 2+2, 3+3, 6+6, 4+4 or 8+8")
    constexpr auto geometry_id = 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((InterfaceGeometry<Line2D3<Node>>{geometry_id, nodes}),
                                      "Number of nodes must be 2+2, 3+3, 6+6, 4+4 or 8+8")
}

KRATOS_TEST_CASE_IN_SUITE(MidGeometryOf2Plus2LineInterfaceIsDefinedByMidPointsOfNodePairs,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto     interface_geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    constexpr auto unused_part_index  = std::size_t{0};
    const auto&    r_mid_geometry     = interface_geometry.GetGeometryPart(unused_part_index);

    const auto expected_mid_points = std::vector<Point>{{-0.5, 0.1, 0.0}, {6.0, 0.1, 0.0}};
    ExpectPointsAreNear(r_mid_geometry.Points(), expected_mid_points);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi_start = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto xi_end   = array_1d<double, 3>{1.0, 0.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_end), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_end), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry  = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi_start  = array_1d<double, 3>{-1.0, 0.0, 0.0};
    const auto xi_end    = array_1d<double, 3>{1.0, 0.0, 0.0};
    const auto xi_middle = array_1d<double, 3>{0.0, 0.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_start), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, xi_middle), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_end), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, xi_middle), 0.0);

    // Values for third shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, xi_start), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, xi_end), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, xi_middle), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectShapeFunctionValuesInNodes_ForThreePlusThreeNodedSurfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry                     = CreateThreePlusThreeNoded3DSurfaceInterfaceGeometry();
    const auto local_coord_first_node_pair  = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto local_coord_second_node_pair = array_1d<double, 3>{1.0, 0.0, 0.0};
    const auto local_coord_third_node_pair  = array_1d<double, 3>{0.0, 1.0, 0.0};

    // Values for first shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, local_coord_first_node_pair), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, local_coord_second_node_pair), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(0, local_coord_third_node_pair), 0.0);

    // Values for second shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, local_coord_first_node_pair), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, local_coord_second_node_pair), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(1, local_coord_third_node_pair), 0.0);

    // Values for third shape function
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, local_coord_first_node_pair), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, local_coord_second_node_pair), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.ShapeFunctionValue(2, local_coord_third_node_pair), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, xi);

    // Note that the shape function values are evaluated per nodal pair!
    Vector expected_result{2};
    expected_result <<= 0.25, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllShapeFunctionValuesAtPosition_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Vector result;
    geometry.ShapeFunctionsValues(result, xi);

    // Note that the shape function values are evaluated per nodal pair!
    Vector expected_result{3};
    expected_result <<= -0.125, 0.375, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_result, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllLocalGradientsAtPosition_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.ShapeFunctionsLocalGradients(result, xi);

    Matrix expected_result(2, 1);
    expected_result <<= -0.5, 0.5;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectAllLocalGradientsAtPosition_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.ShapeFunctionsLocalGradients(result, xi);

    Matrix expected_result(3, 1);
    expected_result <<= 0.0, 1.0, -1.0;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectJacobian_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.Jacobian(result, xi);

    Matrix expected_result(2, 1);
    expected_result <<= 3.25, 0.0;
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(result, expected_result, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectJacobian_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    geometry.Jacobian(result, xi);

    Matrix expected_result(2, 1);
    expected_result <<= 3.0, -0.1;
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(result, expected_result, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectDeterminantOfJacobian_ForTwoPlusTwoNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DeterminantOfJacobian(xi), 3.25, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectDeterminantOfJacobian_ForThreePlusThreeNodedGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    KRATOS_EXPECT_RELATIVE_NEAR(geometry.DeterminantOfJacobian(xi), std::sqrt(9.01), Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_ReturnsCorrectWorkingSpaceDimension, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    KRATOS_EXPECT_EQ(geometry.WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingInverseJacobian, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();
    const auto xi       = array_1d<double, 3>{0.5, 0.0, 0.0};

    Matrix result;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.InverseOfJacobian(result, xi),
        "Inverse of Jacobian is not implemented for the interface geometry")
}

KRATOS_TEST_CASE_IN_SUITE(TwoPlusTwoLineInterfaceGeometry_LengthDomainSizeAndAreaReturnTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    constexpr auto expected_length = 6.5;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.Length(), expected_length, Defaults::relative_tolerance)
    KRATOS_EXPECT_DOUBLE_EQ(geometry.DomainSize(), geometry.Length());

    // The underlying line geometry defines the area to be the same as the length.
    // Although this might be counterintuitive, this behavior propagates to the
    // interface geometry.
    KRATOS_EXPECT_DOUBLE_EQ(geometry.Area(), geometry.Length());
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeNodedLineInterfaceGeometry_LengthDomainSizeAndAreaReturnTheLengthOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    // This number was not calculated by hand, meaning this unit test is a regression test.
    constexpr auto expected_length = 6.504159;
    KRATOS_EXPECT_RELATIVE_NEAR(geometry.Length(), expected_length, Defaults::relative_tolerance)
    KRATOS_EXPECT_DOUBLE_EQ(geometry.DomainSize(), geometry.Length());
    KRATOS_EXPECT_DOUBLE_EQ(geometry.Area(), geometry.Length());
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeNodedSurfaceInterfaceGeometry_AreaAndDomainSizeReturnsTheAreaOfUnderlyingLineGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded3DSurfaceInterfaceGeometry();

    KRATOS_EXPECT_DOUBLE_EQ(geometry.Area(), 0.5);
    KRATOS_EXPECT_DOUBLE_EQ(geometry.DomainSize(), geometry.Area());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalCoordinatesAreCorrectlyMappedToLocalCoordinate_ForTwoPlusTwoNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    const auto global_position = array_1d<double, 3>{4.375, 0.1, 0.0};
    auto       result          = array_1d<double, 3>{};
    result                     = geometry.PointLocalCoordinates(result, global_position);

    const auto expected_result = array_1d<double, 3>{0.5, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(result, expected_result, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GlobalCoordinatesAreCorrectlyMappedToLocalCoordinate_ForThreePlusThreeNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    const auto global_position = array_1d<double, 3>{4.5625, 0.175, 0.0};
    auto       result          = array_1d<double, 3>{};
    result                     = geometry.PointLocalCoordinates(result, global_position);

    const auto expected_result = array_1d<double, 3>{0.5, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(result, expected_result, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GetLocalCoordinatesOfAllNodesOfThreePlusThreeNodedLineInterfaceGeometry,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    Matrix result;
    geometry.PointsLocalCoordinates(result);

    Matrix expected_result{3, 1};
    expected_result <<= -1.0, 1.0, 0.0;
    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TwoPlusTwoLineInterfaceGeometryHasTwoEdgesWithOppositeOrientations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    const auto edges = geometry.GenerateEdges();

    KRATOS_EXPECT_EQ(edges.size(), 2);
    AssertNodeIdsOfGeometry(edges[0], {1, 2});
    AssertNodeIdsOfGeometry(edges[1], {4, 3});
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeLineInterfaceGeometryHasTwoEdgesWithOppositeOrientations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    const auto edges = geometry.GenerateEdges();

    KRATOS_EXPECT_EQ(edges.size(), 2);
    AssertNodeIdsOfGeometry(edges[0], {1, 2, 3});
    AssertNodeIdsOfGeometry(edges[1], {5, 4, 6});
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeLineInterfaceGeometryReturnsTwoLinesForBoundariesEntities,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();

    const auto edges = geometry.GenerateBoundariesEntities();

    KRATOS_EXPECT_EQ(edges.size(), 2);
    KRATOS_EXPECT_TRUE(std::all_of(edges.begin(), edges.end(), [](const auto& r_edge) {
        return r_edge.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear;
    }))
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeSurfaceInterfaceGeometryThrowsWhenCallingGenerateEdges,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded3DSurfaceInterfaceGeometry();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.GenerateEdges(),
                                      "Edges can only be generated for line geometries. "
                                      "This is a surface interface geometry, which does not "
                                      "support edges.");
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeSurfaceInterfaceGeometryHasTwoPlanesWithOppositeOrientations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded3DSurfaceInterfaceGeometry();

    const auto faces = geometry.GenerateFaces();

    KRATOS_EXPECT_EQ(faces.size(), 2);
    AssertNodeIdsOfGeometry(faces[0], {1, 2, 3});
    AssertNodeIdsOfGeometry(faces[1], {4, 6, 5});
}

KRATOS_TEST_CASE_IN_SUITE(SixPlusSixSurfaceInterfaceGeometryHasTwoPlanesWithOppositeOrientations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixPlusSixNoded3DSurfaceInterfaceGeometry();

    const auto faces = geometry.GenerateFaces();

    KRATOS_EXPECT_EQ(faces.size(), 2);
    AssertNodeIdsOfGeometry(faces[0], {1, 2, 3, 4, 5, 6});
    AssertNodeIdsOfGeometry(faces[1], {7, 9, 8, 12, 11, 10});
}

KRATOS_TEST_CASE_IN_SUITE(SixPlusSixSurfaceInterfaceGeometryReturnsTwoTrianglesForBoundariesEntities,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateSixPlusSixNoded3DSurfaceInterfaceGeometry();

    const auto faces = geometry.GenerateBoundariesEntities();

    KRATOS_EXPECT_EQ(faces.size(), 2);
    KRATOS_EXPECT_TRUE(std::all_of(faces.begin(), faces.end(), [](const auto& r_face) {
        return r_face.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle;
    }))
}

KRATOS_TEST_CASE_IN_SUITE(FourPlusFourSurfaceInterfaceGeometryHasTwoPlanesWithOppositeOrientations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateFourPlusFourNoded3DSurfaceInterfaceGeometry();

    const auto faces = geometry.GenerateFaces();

    KRATOS_EXPECT_EQ(faces.size(), 2);
    AssertNodeIdsOfGeometry(faces[0], {1, 2, 3, 4});
    AssertNodeIdsOfGeometry(faces[1], {5, 8, 7, 6});
}

KRATOS_TEST_CASE_IN_SUITE(EightPlusEightSurfaceInterfaceGeometryHasTwoPlanesWithOppositeOrientations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateEightPlusEightNoded3DSurfaceInterfaceGeometry();

    const auto faces = geometry.GenerateFaces();

    KRATOS_EXPECT_EQ(faces.size(), 2);
    AssertNodeIdsOfGeometry(faces[0], {1, 2, 3, 4, 5, 6, 7, 8});
    AssertNodeIdsOfGeometry(faces[1], {9, 12, 11, 10, 16, 15, 14, 13});
}

KRATOS_TEST_CASE_IN_SUITE(EightPlusEightSurfaceInterfaceGeometryReturnsTwoQuadrilateralsForBoundariesEntities,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateEightPlusEightNoded3DSurfaceInterfaceGeometry();

    const auto faces = geometry.GenerateBoundariesEntities();

    KRATOS_EXPECT_EQ(faces.size(), 2);
    KRATOS_EXPECT_TRUE(std::all_of(faces.begin(), faces.end(), [](const auto& r_face) {
        return r_face.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral;
    }))
}

KRATOS_TEST_CASE_IN_SUITE(TwoPlusTwoLineInterfaceGeometryThrowsWhenCallingGenerateFaces,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateTwoPlusTwoNoded2DLineInterfaceGeometry();

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.GenerateFaces(),
                                      "Faces can only be generated for surface geometries. "
                                      "This is a line interface geometry, which does not "
                                      "support faces.");
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingVolume, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Volume(),
                                      "Calling base class 'Volume' method instead of derived class "
                                      "one. Please check the definition of derived class. ")
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometry_Throws_WhenCallingFunctionsRelatedToIntegrationPoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto                          geometry = CreateThreePlusThreeNoded2DLineInterfaceGeometry();
    Geometry<Node>::JacobiansType dummy_jacobian;
    Geometry<Node>::ShapeFunctionsGradientsType dummy_shape_functions_gradients;
    Vector                                      dummy_vector;
    Matrix                                      dummy_matrix;
    IndexType                                   dummy_index = 0;
    GeometryData::IntegrationMethod dummy_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    Geometry<Node>::IntegrationPointsArrayType dummy_integration_points;
    IntegrationInfo dummy_integration_info(dummy_index, dummy_integration_method);
    Geometry<Node>::GeometriesArrayType     dummy_geometries;
    std::vector<Node::CoordinatesArrayType> dummy_coordinates;
    const std::string                       message =
        "This Geometry type does not support functionality related to integration schemes.\n";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Normal(dummy_index), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Normal(dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.UnitNormal(dummy_index), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.UnitNormal(dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.Jacobian(dummy_jacobian, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.Jacobian(dummy_jacobian, dummy_integration_method, dummy_matrix), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.Jacobian(dummy_matrix, dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.Jacobian(dummy_matrix, dummy_index, dummy_integration_method, dummy_matrix), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.DeterminantOfJacobian(dummy_vector, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.DeterminantOfJacobian(dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.InverseOfJacobian(dummy_jacobian, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.InverseOfJacobian(dummy_matrix, dummy_index, dummy_integration_method), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.ShapeFunctionsIntegrationPointsGradients(
                                          dummy_shape_functions_gradients, dummy_integration_method),
                                      message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.ShapeFunctionsIntegrationPointsGradients(
                                          dummy_shape_functions_gradients, dummy_vector, dummy_integration_method),
                                      message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.ShapeFunctionsIntegrationPointsGradients(
            dummy_shape_functions_gradients, dummy_vector, dummy_integration_method, dummy_matrix),
        message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.CreateIntegrationPoints(dummy_integration_points, dummy_integration_info), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry.GetDefaultIntegrationInfo(), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.CreateQuadraturePointGeometries(dummy_geometries, dummy_index,
                                                 dummy_integration_points, dummy_integration_info),
        message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.CreateQuadraturePointGeometries(dummy_geometries, dummy_index, dummy_integration_info), message)
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geometry.GlobalSpaceDerivatives(dummy_coordinates, dummy_index, dummy_index), message)
}
} // namespace Kratos::Testing
