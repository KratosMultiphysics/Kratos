//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

// System includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/multipatch_utilities.h"

namespace Kratos::Testing
{
namespace
{
void GenerateMultipatchTestMesh(ModelPart& rModelPart, const std::string& rElementName)
{
    Quadrilateral2D4<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(2, 0.0, 4.0, 0.0),
        Kratos::make_intrusive<Node>(3, 4.0, 4.0, 0.0),
        Kratos::make_intrusive<Node>(4, 4.0, 0.0, 0.0));
    Parameters settings(R"({
        "element_name": "",
        "number_of_divisions": [4, 4, 1],
        "create_skin_sub_model_part": false
    })");
    settings["element_name"].SetString(rElementName);
    StructuredMeshGeneratorProcess(geometry, rModelPart, settings).Execute();
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(FemMultipatchActivationAndCachedEdges, KratosCoreFastSuite)
{
    for (const std::string element_name : {"Element2D3N", "Element2D4N"}) {
        Model model;
        auto& r_background = model.CreateModelPart("Background");
        GenerateMultipatchTestMesh(r_background, element_name);
        const bool triangular = element_name == "Element2D3N";
        const std::size_t factor = triangular ? 2 : 1;
        KRATOS_EXPECT_EQ(r_background.NumberOfElements(), 16 * factor);
        KRATOS_EXPECT_EQ(r_background.NumberOfNodes(), 25);
        for (const auto& r_element : r_background.Elements()) {
            KRATOS_EXPECT_NEAR(r_element.GetGeometry().Area(), 1.0 / factor, 1e-12);
        }

        FemMultipatchUtilities utility(r_background);
        auto owners = utility.UpdateBackgroundActivation({{1.0, 1.0, 2.0, 2.0}}, 1e-10);
        KRATOS_EXPECT_EQ(owners.size(), factor);
        auto p_retained = r_background.pGetElement(owners.begin()->first);
        KRATOS_EXPECT_FALSE(p_retained->Is(ACTIVE));
        const auto edges = utility.ExtractSurrogateEdges(owners);
        KRATOS_EXPECT_EQ(edges.size(), 1);
        KRATOS_EXPECT_EQ(edges.at(1).size(), 4);
        double twice_signed_area = 0.0;
        for (const auto& r_nodes : edges.at(1)) {
            const auto& r_first = r_background.GetNode(r_nodes[0]);
            const auto& r_second = r_background.GetNode(r_nodes[1]);
            twice_signed_area += r_first.X() * r_second.Y() - r_second.X() * r_first.Y();
        }
        // Following active elements makes the hole boundary clockwise.
        KRATOS_EXPECT_NEAR(twice_signed_area, -2.0, 1e-12);

        owners = utility.UpdateBackgroundActivation({}, 1e-10);
        KRATOS_EXPECT_TRUE(utility.ExtractSurrogateEdges(owners).empty());
        KRATOS_EXPECT_TRUE(p_retained->Is(ACTIVE));
        KRATOS_EXPECT_EQ(r_background.pGetElement(p_retained->Id()), p_retained);

        // Reuse the cached topology for a different aligned inactive region.
        owners = utility.UpdateBackgroundActivation({{2.0, 2.0, 3.0, 3.0}}, 1e-10);
        KRATOS_EXPECT_EQ(owners.size(), factor);
        KRATOS_EXPECT_TRUE(p_retained->Is(ACTIVE));
        KRATOS_EXPECT_EQ(utility.ExtractSurrogateEdges(owners).at(1).size(), 4);
        KRATOS_EXPECT_EQ(r_background.NumberOfElements(), 16 * factor);
    }
}

KRATOS_TEST_CASE_IN_SUITE(FemMultipatchQuadrilateralGaussPoints, KratosCoreFastSuite)
{
    Model model;
    auto& r_background = model.CreateModelPart("Background");
    GenerateMultipatchTestMesh(r_background, "Element2D4N");
    FemMultipatchUtilities utility(r_background);
    utility.UpdateBackgroundActivation({{1.0, 1.0, 2.0, 2.0}}, 1e-10);
    const auto active_points = utility.GetGaussPoints(r_background);
    KRATOS_EXPECT_EQ(active_points.size(), 15 * 4);
    const auto all_points = utility.GetGaussPoints(r_background, false);
    KRATOS_EXPECT_EQ(all_points.size(), 16 * 4);
    const double gauss_offset = 0.5 / std::sqrt(3.0);
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_EXPECT_EQ(all_points[i].first, 1);
        KRATOS_EXPECT_NEAR(std::abs(all_points[i].second[0] - 0.5), gauss_offset, 1e-12);
        KRATOS_EXPECT_NEAR(std::abs(all_points[i].second[1] - 0.5), gauss_offset, 1e-12);
        KRATOS_EXPECT_NEAR(all_points[i].second[2], 0.0, 1e-12);
    }
}

KRATOS_TEST_CASE_IN_SUITE(FemMultipatchTopologyRebuild, KratosCoreFastSuite)
{
    Model model;
    auto& r_background = model.CreateModelPart("Background");
    GenerateMultipatchTestMesh(r_background, "Element2D4N");
    FemMultipatchUtilities utility(r_background);
    r_background.CreateNewNode(26, 5.0, 0.0, 0.0);
    r_background.CreateNewNode(27, 6.0, 0.0, 0.0);
    r_background.CreateNewNode(28, 5.0, 1.0, 0.0);
    auto p_triangle = r_background.CreateNewElement(
        "Element2D3N", 17, {26, 27, 28}, r_background.pGetProperties(0));
    utility.RebuildTopology();
    utility.UpdateBackgroundActivation({}, 1e-10);
    KRATOS_EXPECT_TRUE(p_triangle->Is(ACTIVE));
    KRATOS_EXPECT_EQ(utility.GetGaussPoints(r_background).size(), 65);
}

KRATOS_TEST_CASE_IN_SUITE(FemMultipatchNonManifoldEdge, KratosCoreFastSuite)
{
    Model model;
    auto& r_background = model.CreateModelPart("Background");
    r_background.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_background.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_background.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_background.CreateNewNode(4, 0.0, -1.0, 0.0);
    r_background.CreateNewNode(5, 0.5, 1.0, 0.0);
    auto p_properties = r_background.CreateNewProperties(0);
    r_background.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    r_background.CreateNewElement("Element2D3N", 2, {2, 1, 4}, p_properties);
    r_background.CreateNewElement("Element2D3N", 3, {1, 2, 5}, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(FemMultipatchUtilities utility(r_background),
                                    "Non-manifold background edge");
}

} // namespace Kratos::Testing
