//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli


// System includes
#include <array>
#include <cmath>
#include <vector>

// Project includes
#include "includes/kratos_application.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "modeler/gap_surrogate_boundary_modeler.h"
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

void AddCircularSkinLoop(
    ModelPart& rSkin,
    Properties::Pointer pProperties,
    std::size_t& rNextNodeId,
    std::size_t& rNextConditionId,
    const double CenterX,
    const double CenterY,
    const double Radius,
    const std::size_t NumberOfSegments)
{
    const double pi = std::acos(-1.0);
    std::vector<std::size_t> node_ids(NumberOfSegments);
    for (std::size_t i = 0; i < NumberOfSegments; ++i) {
        const double angle = 2.0 * pi * i / NumberOfSegments;
        node_ids[i] = rNextNodeId++;
        rSkin.CreateNewNode(
            node_ids[i],
            CenterX + Radius * std::cos(angle),
            CenterY + Radius * std::sin(angle),
            0.0);
    }
    for (std::size_t i = 0; i < NumberOfSegments; ++i) {
        rSkin.CreateNewCondition(
            "LineCondition2D2N",
            rNextConditionId++,
            std::vector<ModelPart::IndexType>{
                node_ids[i], node_ids[(i + 1) % NumberOfSegments]},
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
    modeler.PrepareGeometryModel(); // Standard pre-SetupModelPart lifecycle call.
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

KRATOS_TEST_CASE_IN_SUITE(GapSurrogateBoundaryModelerCreatesOneGapElementPerFace, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    std::size_t next_node_id = 1;
    std::size_t next_condition_id = 1;
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        0.73, 0.67, 1.31, 1.37);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [2.0, 2.0, 0.0],
        "number_of_elements" : [12, 12, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "gap_element_name" : "Element2D4N"
    })");
    GapSurrogateBoundaryModeler modeler(model, settings);
    modeler.SetupModelPart();

    const auto& r_background = model.GetModelPart("Background");
    const auto& r_surrogate =
        r_background.GetSubModelPart("SurrogateBoundary");
    const auto& r_gap_elements =
        r_background.GetSubModelPart("GapElements");
    KRATOS_EXPECT_GT(r_surrogate.NumberOfConditions(), 0);
    KRATOS_EXPECT_EQ(
        r_gap_elements.NumberOfElements(),
        r_surrogate.NumberOfConditions());

    for (const auto& r_condition : r_surrogate.Conditions()) {
        const Matrix& r_boundary_quadrature =
            r_condition.GetValue(SURROGATE_BOUNDARY_PROJECTION);
        KRATOS_EXPECT_EQ(r_boundary_quadrature.size1(), 2);
        KRATOS_EXPECT_EQ(r_boundary_quadrature.size2(), 7);
        KRATOS_EXPECT_GT(r_boundary_quadrature(0, 6), 0.0);
        KRATOS_EXPECT_GT(r_boundary_quadrature(1, 6), 0.0);
    }
    for (const auto& r_element : r_gap_elements.Elements()) {
        KRATOS_EXPECT_EQ(
            r_element.GetGeometry().GetGeometryType(),
            GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry);
        KRATOS_EXPECT_EQ(
            r_element.GetGeometry().IntegrationPointsNumber(
                GeometryData::IntegrationMethod::GI_GAUSS_2),
            4);
        KRATOS_EXPECT_EQ(r_element.GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
        KRATOS_EXPECT_EQ(r_element.GetValue(NEIGHBOUR_CONDITIONS).size(), 1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(
    GapSurrogateBoundaryModelerOmitsZeroMeasureTangencyPatches,
    KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    std::size_t next_node_id = 1;
    std::size_t next_condition_id = 1;

    // R/h=8 and the cardinal polygon vertices lie exactly on grid nodes. At
    // each grid tangency the surrogate topology contains a radial spike whose
    // two endpoints project to the same point. The corresponding ruled patch
    // has zero area and must not be instantiated as a volume element.
    AddCircularSkinLoop(
        r_skin,
        p_skin_properties,
        next_node_id,
        next_condition_id,
        0.0,
        0.0,
        0.5,
        128);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [-1.0, -1.0, 0.0],
        "upper_point" : [1.0, 1.0, 0.0],
        "number_of_elements" : [32, 32, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "gap_element_name" : "Element2D4N"
    })");
    GapSurrogateBoundaryModeler modeler(model, settings);
    modeler.SetupModelPart();

    const auto& r_background = model.GetModelPart("Background");
    const auto& r_surrogate =
        r_background.GetSubModelPart("SurrogateBoundary");
    const auto& r_gap_elements =
        r_background.GetSubModelPart("GapElements");
    KRATOS_EXPECT_GT(r_gap_elements.NumberOfElements(), 0);
    KRATOS_EXPECT_LT(
        r_gap_elements.NumberOfElements(),
        r_surrogate.NumberOfConditions());

    std::size_t number_of_zero_measure_boundary_pieces = 0;
    for (const auto& r_condition : r_surrogate.Conditions()) {
        const Matrix& r_boundary_quadrature =
            r_condition.GetValue(SURROGATE_BOUNDARY_PROJECTION);
        KRATOS_EXPECT_EQ(r_boundary_quadrature.size1(), 2);
        KRATOS_EXPECT_EQ(r_boundary_quadrature.size2(), 7);
        if (r_boundary_quadrature(0, 6) <= 1.0e-12 &&
            r_boundary_quadrature(1, 6) <= 1.0e-12) {
            ++number_of_zero_measure_boundary_pieces;
        }
    }
    KRATOS_EXPECT_GT(number_of_zero_measure_boundary_pieces, 0);
}

KRATOS_TEST_CASE_IN_SUITE(SurrogateBoundaryModelerUpdatesMovingSkinOnPersistentGrid, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    std::size_t next_node_id = 1;
    std::size_t next_condition_id = 1;
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        0.55, 0.65, 1.15, 1.35);

    ModelPart& r_background = model.CreateModelPart("Background");
    r_background.AddNodalSolutionStepVariable(VELOCITY);
    r_background.SetBufferSize(3);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [3.0, 2.0, 0.0],
        "number_of_elements" : [15, 10, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "is_moving" : true
    })");
    SurrogateBoundaryModeler modeler(model, settings);
    modeler.PrepareGeometryModel(); // Initial lifecycle call is a no-op.
    modeler.SetupModelPart();

    auto& r_volume = r_background.GetSubModelPart("FluidDomain");
    auto& r_grid = r_background.GetSubModelPart("BackgroundGrid");
    const std::size_t initial_active_elements = r_volume.NumberOfElements();
    std::vector<std::size_t> initial_active_ids;
    initial_active_ids.reserve(initial_active_elements);
    for (const auto& r_element : r_volume.Elements()) {
        initial_active_ids.push_back(r_element.Id());
    }

    KRATOS_EXPECT_EQ(r_grid.NumberOfElements(), 150);
    KRATOS_EXPECT_EQ(r_grid.NumberOfNodes(), 176);
    KRATOS_EXPECT_EQ(r_background.NumberOfElements(), 150);

    Node& r_history_node = *r_grid.NodesBegin();
    r_history_node.FastGetSolutionStepValue(VELOCITY, 1) =
        array_1d<double, 3>{1.25, -0.5, 0.0};
    const std::size_t history_node_id = r_history_node.Id();

    for (auto& r_node : r_skin.Nodes()) {
        r_node.X() += 0.8;
    }
    modeler.PrepareGeometryModel();

    std::vector<std::size_t> updated_active_ids;
    updated_active_ids.reserve(r_volume.NumberOfElements());
    for (const auto& r_element : r_volume.Elements()) {
        updated_active_ids.push_back(r_element.Id());
    }
    KRATOS_EXPECT_FALSE(initial_active_ids == updated_active_ids);
    KRATOS_EXPECT_EQ(r_volume.NumberOfElements(), initial_active_elements);
    KRATOS_EXPECT_EQ(r_background.NumberOfElements(), 150);
    KRATOS_EXPECT_EQ(r_grid.NumberOfElements(), 150);
    const array_1d<double, 3> expected_history{1.25, -0.5, 0.0};
    KRATOS_EXPECT_VECTOR_NEAR(
        r_background.GetNode(history_node_id).FastGetSolutionStepValue(VELOCITY, 1),
        expected_history,
        1.0e-15);

    const auto& r_surrogate =
        r_background.GetSubModelPart("SurrogateBoundary");
    KRATOS_EXPECT_GT(r_surrogate.NumberOfConditions(), 0);
    KRATOS_EXPECT_EQ(
        r_background.NumberOfConditions(),
        r_surrogate.NumberOfConditions());
    for (const auto& r_condition : r_surrogate.Conditions()) {
        const auto& r_neighbours = r_condition.GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_EXPECT_EQ(r_neighbours.size(), 1);
        KRATOS_EXPECT_TRUE(r_volume.HasElement(r_neighbours[0].Id()));
    }
}

KRATOS_TEST_CASE_IN_SUITE(GapSurrogateBoundaryModelerRebuildsMovingGapWithoutAccumulation, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_skin = model.CreateModelPart("Skin");
    auto p_skin_properties = r_skin.CreateNewProperties(1);
    std::size_t next_node_id = 1;
    std::size_t next_condition_id = 1;
    AddRectangularSkinLoop(
        r_skin, p_skin_properties, next_node_id, next_condition_id,
        0.63, 0.67, 1.17, 1.33);

    Parameters settings(R"({
        "input_model_part_name" : "Skin",
        "output_model_part_name" : "Background",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [2.5, 2.0, 0.0],
        "number_of_elements" : [15, 12, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "gap_element_name" : "Element2D4N",
        "is_moving" : true
    })");
    GapSurrogateBoundaryModeler modeler(model, settings);
    modeler.PrepareGeometryModel();
    modeler.SetupModelPart();

    auto& r_background = model.GetModelPart("Background");
    auto& r_grid = r_background.GetSubModelPart("BackgroundGrid");
    auto& r_surrogate = r_background.GetSubModelPart("SurrogateBoundary");
    auto& r_gap = r_background.GetSubModelPart("GapElements");
    KRATOS_EXPECT_EQ(r_grid.NumberOfElements(), 180);
    KRATOS_EXPECT_EQ(r_gap.NumberOfElements(), r_surrogate.NumberOfConditions());
    KRATOS_EXPECT_EQ(
        r_background.NumberOfElements(),
        r_grid.NumberOfElements() + r_gap.NumberOfElements());

    for (auto& r_node : r_skin.Nodes()) {
        r_node.X() += 0.5;
    }
    modeler.PrepareGeometryModel();

    KRATOS_EXPECT_EQ(r_grid.NumberOfElements(), 180);
    KRATOS_EXPECT_EQ(r_gap.NumberOfElements(), r_surrogate.NumberOfConditions());
    KRATOS_EXPECT_EQ(
        r_background.NumberOfElements(),
        r_grid.NumberOfElements() + r_gap.NumberOfElements());
    for (const auto& r_element : r_gap.Elements()) {
        KRATOS_EXPECT_EQ(r_element.GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
        KRATOS_EXPECT_EQ(r_element.GetValue(NEIGHBOUR_CONDITIONS).size(), 1);
    }
}

} // namespace Kratos::Testing
