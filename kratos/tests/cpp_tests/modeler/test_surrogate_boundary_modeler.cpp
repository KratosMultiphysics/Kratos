//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// System includes
#include <array>
#include <cmath>

// Project includes
#include "includes/kratos_application.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "modeler/surrogate_boundary_modeler.h"
#include "testing/testing.h"
#include "utilities/structured_surrogate_boundary_utility.h"

namespace Kratos::Testing {

namespace {

void AddRectangularSkinLoop(
    ModelPart& rSkin,
    Properties::Pointer pProperties,
    std::size_t& rNextNodeId,
    std::size_t& rNextConditionId,
    const double XMin,
    const double YMin,
    const double XMax,
    const double YMax)
{
    const std::array<std::array<double, 2>, 4> coordinates{{
        {{XMin, YMin}}, {{XMax, YMin}}, {{XMax, YMax}}, {{XMin, YMax}}
    }};
    std::array<std::size_t, 4> node_ids;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        node_ids[i] = rNextNodeId++;
        rSkin.CreateNewNode(node_ids[i], coordinates[i][0], coordinates[i][1], 0.0);
    }
    for (std::size_t i = 0; i < node_ids.size(); ++i) {
        rSkin.CreateNewCondition(
            "LineCondition2D2N",
            rNextConditionId++,
            std::vector<ModelPart::IndexType>{node_ids[i], node_ids[(i + 1) % node_ids.size()]},
            pProperties);
    }
}

array_1d<double, 3> ElementCenter(const Element& rElement)
{
    array_1d<double, 3> center{0.0, 0.0, 0.0};
    for (const auto& r_node : rElement.GetGeometry()) {
        center += r_node.Coordinates();
    }
    center /= rElement.GetGeometry().size();
    return center;
}

bool IsStrictlyInsideRectangle(
    const array_1d<double, 3>& rPoint,
    const double XMin,
    const double YMin,
    const double XMax,
    const double YMax)
{
    return rPoint[0] > XMin && rPoint[0] < XMax &&
           rPoint[1] > YMin && rPoint[1] < YMax;
}

bool ProjectionTouchesRectangle(
    const Matrix& rProjection,
    const double XMin,
    const double YMin,
    const double XMax,
    const double YMax)
{
    constexpr double tolerance = 1.0e-12;
    for (std::size_t i = 0; i < rProjection.size1(); ++i) {
        const double x = rProjection(i, 0);
        const double y = rProjection(i, 1);
        const bool on_vertical_side =
            (std::abs(x - XMin) < tolerance || std::abs(x - XMax) < tolerance) &&
            y >= YMin - tolerance && y <= YMax + tolerance;
        const bool on_horizontal_side =
            (std::abs(y - YMin) < tolerance || std::abs(y - YMax) < tolerance) &&
            x >= XMin - tolerance && x <= XMax + tolerance;
        if (on_vertical_side || on_horizontal_side) {
            return true;
        }
    }
    return false;
}

} // unnamed namespace

KRATOS_TEST_CASE_IN_SUITE(StructuredSurrogateBoundaryTopologySupports2DAnd3D, KratosCoreFastSuite)
{
    using Utility2D = StructuredSurrogateBoundaryUtility<2>;
    const Utility2D::IndexArrayType size_2d{3, 3};
    std::vector<bool> active_2d(9, false);
    active_2d[Utility2D::LinearIndex({1, 1}, size_2d)] = true;
    KRATOS_EXPECT_EQ(Utility2D::ExtractBoundaryFaces(active_2d, size_2d).size(), 4);

    using Utility3D = StructuredSurrogateBoundaryUtility<3>;
    const Utility3D::IndexArrayType size_3d{3, 3, 3};
    std::vector<bool> active_3d(27, false);
    active_3d[Utility3D::LinearIndex({1, 1, 1}, size_3d)] = true;
    KRATOS_EXPECT_EQ(Utility3D::ExtractBoundaryFaces(active_3d, size_3d).size(), 6);
}

KRATOS_TEST_CASE_IN_SUITE(SurrogateBoundaryModelerCreatesFaceBased2DMesh, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    r_skin.CreateNewNode(1, 0.75, 0.75, 0.0);
    r_skin.CreateNewNode(2, 1.25, 0.75, 0.0);
    r_skin.CreateNewNode(3, 1.25, 1.25, 0.0);
    r_skin.CreateNewNode(4, 0.75, 1.25, 0.0);
    r_skin.CreateNewCondition("LineCondition2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_skin_properties);
    r_skin.CreateNewCondition("LineCondition2D2N", 2, std::vector<ModelPart::IndexType>{2, 3}, p_skin_properties);
    r_skin.CreateNewCondition("LineCondition2D2N", 3, std::vector<ModelPart::IndexType>{3, 4}, p_skin_properties);
    r_skin.CreateNewCondition("LineCondition2D2N", 4, std::vector<ModelPart::IndexType>{4, 1}, p_skin_properties);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [2.0, 2.0, 0.0],
        "number_of_elements" : [8, 8, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0
    })");
    SurrogateBoundaryModeler modeler(model, settings);
    modeler.SetupModelPart();

    const auto& r_background = model.GetModelPart("Background");
    const auto& r_volume = r_background.GetSubModelPart("FluidDomain");
    const auto& r_surrogate = r_background.GetSubModelPart("SurrogateBoundary");
    KRATOS_EXPECT_GT(r_volume.NumberOfElements(), 0);
    KRATOS_EXPECT_LT(r_volume.NumberOfElements(), 64);
    KRATOS_EXPECT_GT(r_surrogate.NumberOfConditions(), 0);

    for (const auto& r_condition : r_surrogate.Conditions()) {
        KRATOS_EXPECT_EQ(r_condition.GetGeometry().size(), 4);
        KRATOS_EXPECT_EQ(r_condition.GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
        const Matrix& r_face_coordinates = r_condition.GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
        KRATOS_EXPECT_EQ(r_face_coordinates.size1(), 2);
        KRATOS_EXPECT_EQ(r_face_coordinates.size2(), 3);
        const Matrix& r_projection = r_condition.GetValue(SURROGATE_BOUNDARY_PROJECTION);
        KRATOS_EXPECT_EQ(r_projection.size1(), 2);
        KRATOS_EXPECT_EQ(r_projection.size2(), 3);

        const array_1d<double, 3> center{
            0.5 * (r_face_coordinates(0, 0) + r_face_coordinates(1, 0)),
            0.5 * (r_face_coordinates(0, 1) + r_face_coordinates(1, 1)),
            0.0};
        const double tangent_x = r_face_coordinates(1, 0) - r_face_coordinates(0, 0);
        const double tangent_y = r_face_coordinates(1, 1) - r_face_coordinates(0, 1);
        array_1d<double, 3> normal{tangent_y, -tangent_x, 0.0};
        normal /= norm_2(normal);
        const array_1d<double, 3> shift{
            0.5 * (r_projection(0, 0) + r_projection(1, 0)) - center[0],
            0.5 * (r_projection(0, 1) + r_projection(1, 1)) - center[1],
            0.0};
        KRATOS_EXPECT_GT(inner_prod(normal, shift), 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SurrogateBoundaryModelerSupportsMultipleHoles, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    std::size_t next_node_id = 1;
    std::size_t next_condition_id = 1;
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        0.75, 0.65, 1.35, 1.35);
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        2.65, 0.65, 3.25, 1.35);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [4.0, 2.0, 0.0],
        "number_of_elements" : [20, 10, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0
    })");
    SurrogateBoundaryModeler modeler(model, settings);
    modeler.SetupModelPart();

    const auto& r_background = model.GetModelPart("Background");
    const auto& r_volume = r_background.GetSubModelPart("FluidDomain");
    const auto& r_surrogate = r_background.GetSubModelPart("SurrogateBoundary");
    KRATOS_EXPECT_GT(r_volume.NumberOfElements(), 0);
    KRATOS_EXPECT_LT(r_volume.NumberOfElements(), 200);

    for (const auto& r_element : r_volume.Elements()) {
        const auto center = ElementCenter(r_element);
        KRATOS_EXPECT_FALSE(IsStrictlyInsideRectangle(center, 0.75, 0.65, 1.35, 1.35));
        KRATOS_EXPECT_FALSE(IsStrictlyInsideRectangle(center, 2.65, 0.65, 3.25, 1.35));
    }

    std::size_t first_hole_projection_count = 0;
    std::size_t second_hole_projection_count = 0;
    for (const auto& r_condition : r_surrogate.Conditions()) {
        const Matrix& r_projection = r_condition.GetValue(SURROGATE_BOUNDARY_PROJECTION);
        first_hole_projection_count += ProjectionTouchesRectangle(
            r_projection, 0.75, 0.65, 1.35, 1.35);
        second_hole_projection_count += ProjectionTouchesRectangle(
            r_projection, 2.65, 0.65, 3.25, 1.35);
    }
    KRATOS_EXPECT_GT(first_hole_projection_count, 0);
    KRATOS_EXPECT_GT(second_hole_projection_count, 0);
}

KRATOS_TEST_CASE_IN_SUITE(SurrogateBoundaryModelerSupportsOuterBoundaryAndHole, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    std::size_t next_node_id = 1;
    std::size_t next_condition_id = 1;
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        0.4, 0.4, 2.6, 2.6);
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        1.2, 1.2, 1.8, 1.8);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [3.0, 3.0, 0.0],
        "number_of_elements" : [24, 24, 1],
        "physical_domain" : "inside",
        "lambda" : 1.0
    })");
    SurrogateBoundaryModeler modeler(model, settings);
    modeler.SetupModelPart();

    const auto& r_background = model.GetModelPart("Background");
    const auto& r_volume = r_background.GetSubModelPart("FluidDomain");
    const auto& r_surrogate = r_background.GetSubModelPart("SurrogateBoundary");
    KRATOS_EXPECT_GT(r_volume.NumberOfElements(), 0);
    KRATOS_EXPECT_LT(r_volume.NumberOfElements(), 24 * 24);

    for (const auto& r_element : r_volume.Elements()) {
        const auto center = ElementCenter(r_element);
        KRATOS_EXPECT_TRUE(IsStrictlyInsideRectangle(center, 0.4, 0.4, 2.6, 2.6));
        KRATOS_EXPECT_FALSE(IsStrictlyInsideRectangle(center, 1.2, 1.2, 1.8, 1.8));
    }

    std::size_t outer_projection_count = 0;
    std::size_t hole_projection_count = 0;
    for (const auto& r_condition : r_surrogate.Conditions()) {
        const Matrix& r_projection = r_condition.GetValue(SURROGATE_BOUNDARY_PROJECTION);
        outer_projection_count += ProjectionTouchesRectangle(
            r_projection, 0.4, 0.4, 2.6, 2.6);
        hole_projection_count += ProjectionTouchesRectangle(
            r_projection, 1.2, 1.2, 1.8, 1.8);
    }
    KRATOS_EXPECT_GT(outer_projection_count, 0);
    KRATOS_EXPECT_GT(hole_projection_count, 0);
}

} // namespace Kratos::Testing
