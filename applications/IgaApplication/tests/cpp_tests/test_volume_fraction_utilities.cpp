//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes
#include <array>
#include <cmath>
#include <vector>

// Project includes
#include "containers/model.h"
#include "includes/global_variables.h"
#include "testing/testing.h"
#include "custom_utilities/volume_fraction_utilities.h"

namespace Kratos::Testing
{

namespace
{

using IndexType = std::size_t;

/// Appends a closed polygonal loop of LineCondition2D2N to the skin model part
void CreateSkinLoop(
    ModelPart& rSkinModelPart,
    const std::vector<std::array<double, 2>>& rVertices)
{
    if (!rSkinModelPart.HasProperties(0)) {
        rSkinModelPart.CreateNewProperties(0);
    }
    Properties::Pointer p_properties = rSkinModelPart.pGetProperties(0);

    const IndexType first_node_id = rSkinModelPart.NumberOfNodes() + 1;
    const IndexType first_condition_id = rSkinModelPart.NumberOfConditions() + 1;
    const IndexType number_of_vertices = rVertices.size();
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        rSkinModelPart.CreateNewNode(first_node_id + i, rVertices[i][0], rVertices[i][1], 0.0);
    }
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        const std::vector<IndexType> node_ids{first_node_id + i, first_node_id + (i + 1) % number_of_vertices};
        rSkinModelPart.CreateNewCondition("LineCondition2D2N", first_condition_id + i, node_ids, p_properties);
    }
}

/// Vertices of a circle approximated by a regular polygon, counter-clockwise
std::vector<std::array<double, 2>> CircleVertices(
    const double CenterX,
    const double CenterY,
    const double Radius,
    const IndexType NumberOfVertices)
{
    std::vector<std::array<double, 2>> vertices(NumberOfVertices);
    for (IndexType k = 0; k < NumberOfVertices; ++k) {
        const double angle = 2.0 * Globals::Pi * static_cast<double>(k) / static_cast<double>(NumberOfVertices);
        vertices[k] = {CenterX + Radius * std::cos(angle), CenterY + Radius * std::sin(angle)};
    }
    return vertices;
}

} // namespace

// Domain on a 4 x 4 grid over [0, 1]^2, with analytical volume fractions:
// - Outer skin (given clockwise): square [0.125, 0.875]^2 with its lower-right corner cut by the line y = x - 0.5,
//   which crosses span boundaries and passes through the span corner (0.75, 0.25). Fractions are multiples of 1/8.
// - Inner skin: two circular holes centred at span corners, so each takes a quarter disk from its four spans:
//   radius 0.125 at (0.5, 0.5) removes pi/16 of each span, radius 0.0625 at (0.25, 0.5) removes pi/64.
KRATOS_TEST_CASE_IN_SUITE(VolumeFractionUtilitiesDomainWithHoles, KratosIgaFastSuite)
{
    Model model;
    ModelPart& r_outer_skin = model.CreateModelPart("outer_skin");
    ModelPart& r_inner_skin = model.CreateModelPart("inner_skin");
    CreateSkinLoop(r_outer_skin, {{{0.125, 0.125}}, {{0.125, 0.875}}, {{0.875, 0.875}}, {{0.875, 0.375}}, {{0.625, 0.125}}});
    CreateSkinLoop(r_inner_skin, CircleVertices(0.5, 0.5, 0.125, 1024));
    CreateSkinLoop(r_inner_skin, CircleVertices(0.25, 0.5, 0.0625, 1024));

    const std::vector<double> spans{0.0, 0.25, 0.5, 0.75, 1.0};
    const Matrix volume_fractions = VolumeFractionUtilities::ComputeKnotSpanVolumeFractions(
        r_outer_skin, r_inner_skin, spans, spans);

    // Entry (i, j): i-th span in x, j-th span in y
    const double large_hole = Globals::Pi / 16.0;
    const double small_hole = Globals::Pi / 64.0;
    const std::array<std::array<double, 4>, 4> expected_volume_fractions{{
        {0.25,  0.5 - small_hole,              0.5 - small_hole,              0.25},
        {0.5,   1.0 - large_hole - small_hole, 1.0 - large_hole - small_hole, 0.5 },
        {0.375, 1.0 - large_hole,              1.0 - large_hole,              0.5 },
        {0.0,   0.375,                         0.5,                           0.25}}};

    KRATOS_EXPECT_EQ(volume_fractions.size1(), 4);
    KRATOS_EXPECT_EQ(volume_fractions.size2(), 4);
    for (IndexType i = 0; i < 4; ++i) {
        for (IndexType j = 0; j < 4; ++j) {
            KRATOS_EXPECT_NEAR(volume_fractions(i, j), expected_volume_fractions[i][j], 1.0e-5);
        }
    }
}

} // namespace Kratos::Testing
