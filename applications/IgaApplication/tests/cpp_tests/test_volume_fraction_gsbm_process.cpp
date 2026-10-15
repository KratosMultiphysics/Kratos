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
#include <vector>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "iga_application_variables.h"
#include "includes/kratos_components.h"

namespace Kratos::Testing
{

namespace
{

using IndexType = std::size_t;

/// Appends a closed polygonal loop of LineCondition2D2N to the skin
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
        rSkinModelPart.CreateNewNode(
            first_node_id + i, rVertices[i][0], rVertices[i][1], 0.0
        );
    }
    for (IndexType i = 0; i < number_of_vertices; ++i) {
        const std::vector<IndexType> node_ids{
            first_node_id + i, first_node_id + (i + 1) % number_of_vertices
        };
        rSkinModelPart.CreateNewCondition(
            "LineCondition2D2N", first_condition_id + i, node_ids, p_properties
        );
    }
}

} // namespace

// Background patch [0, 1] x [0, 0.75] with 4 x 3 knot spans of size 0.25 and
// degrees (2, 3), so any mix-up of u and v fails. The domain is the rectangle
// [0.125, 0.875] x [0.125, 0.625] with a square hole [0.3125, 0.4375]^2 inside
// span (1, 1), which gives volume fractions 0.25, 0.5, 0.75 and 1. Each element
// and facet side is checked against the span found from its center, shifted
// off the knot line towards its own side for facets.
KRATOS_TEST_CASE_IN_SUITE(VolumeFractionGsbmProcess, KratosIgaFastSuite)
{
    Model model;
    model.CreateModelPart("IgaModelPart");

    NurbsGeometryModeler nurbs_modeler(model, Parameters(R"({
        "model_part_name"     : "IgaModelPart",
        "geometry_name"       : "background",
        "lower_point_xyz"     : [0.0, 0.0, 0.0],
        "upper_point_xyz"     : [1.0, 0.75, 0.0],
        "lower_point_uvw"     : [0.0, 0.0, 0.0],
        "upper_point_uvw"     : [1.0, 0.75, 0.0],
        "polynomial_order"    : [2, 3],
        "number_of_knot_spans": [4, 3]
    })"));
    nurbs_modeler.SetupGeometryModel();

    const Modeler& r_registered_modeler = KratosComponents<Modeler>::Get("IgaModelerGsbm");
    const auto p_gsbm_modeler = r_registered_modeler.Create(model, Parameters(R"({
        "analysis_model_part_name": "IgaModelPart",
        "geometry_name"           : "background",
        "element_condition_list"  : [{
            "type"                            : "element",
            "name"                            : "GenericElement",
            "iga_model_part"                  : "Domain",
            "shape_function_derivatives_order": 4
        }, {
            "type"                            : "condition",
            "name"                            : "GenericCondition",
            "iga_model_part"                  : "Facets",
            "shape_function_derivatives_order": 4,
            "gsbm_parameters"                 : {}
        }]
    })"));
    p_gsbm_modeler->SetupModelPart();

    CreateSkinLoop(
        model.CreateModelPart("OuterSkin"),
        {{0.125, 0.125}, {0.875, 0.125}, {0.875, 0.625}, {0.125, 0.625}}
    );
    CreateSkinLoop(
        model.CreateModelPart("InnerSkin"),
        {{0.3125, 0.3125}, {0.4375, 0.3125}, {0.4375, 0.4375}, {0.3125, 0.4375}}
    );

    VolumeFractionGsbmProcess gsbm_process(model, Parameters(R"({
        "model_part_name"                   : "IgaModelPart",
        "geometry_name"                     : "background",
        "elements_model_part_name"          : "IgaModelPart.Domain",
        "facets_model_part_name"            : "IgaModelPart.Facets",
        "skin_model_part_outer_initial_name": "OuterSkin",
        "skin_model_part_inner_initial_name": "InnerSkin"
    })"));
    gsbm_process.ExecuteInitialize();

    // Volume fraction of span (i, j): the rectangle covers half of the border
    // spans, and the hole a quarter of span (1, 1)
    const std::array<std::array<double, 3>, 4> expected_volume_fractions{{
        {0.25, 0.5 , 0.25},
        {0.5 , 0.75, 0.5 },
        {0.5 , 1.0 , 0.5 },
        {0.25, 0.5 , 0.25}
    }};
    const auto expected_volume_fraction = [&](const array_1d<double, 3>& rPoint) {
        const IndexType i = static_cast<IndexType>(rPoint[0] / 0.25);
        const IndexType j = static_cast<IndexType>(rPoint[1] / 0.25);
        return expected_volume_fractions[i][j];
    };
    const double tolerance = 1.0e-10;

    // Elements: 4 x 3 spans with (2 + 1) x (3 + 1) Gauss points
    const ModelPart& r_domain = model.GetModelPart("IgaModelPart.Domain");
    KRATOS_EXPECT_EQ(r_domain.NumberOfElements(), 144);
    for (const auto& r_element : r_domain.Elements()) {
        const array_1d<double, 3> gauss_point = r_element.GetGeometry().Center().Coordinates();
        KRATOS_EXPECT_NEAR(
            r_element.GetValue(VOLUME_FRACTION), expected_volume_fraction(gauss_point), tolerance
        );
    }

    // Facets: 3 vertical lines x 3 spans x (3 + 1) points and 2 horizontal
    // lines x 4 spans x (2 + 1) points
    const ModelPart& r_facets = model.GetModelPart("IgaModelPart.Facets");
    KRATOS_EXPECT_EQ(r_facets.NumberOfConditions(), 60);
    const double offset = 1.0e-6;
    const IndexType point_index = 0;
    for (const auto& r_facet : r_facets.Conditions()) {
        // Part 0 is the - side (behind the normal), part 1 the + side
        for (IndexType side = 0; side < 2; ++side) {
            const auto& r_side = r_facet.GetGeometry().GetGeometryPart(side);
            const double sign = side == 0 ? -1.0 : 1.0;
            const array_1d<double, 3> point = r_side.Center().Coordinates()
                + sign * offset * r_side.UnitNormal(point_index);
            KRATOS_EXPECT_NEAR(
                r_side.GetValue(VOLUME_FRACTION), expected_volume_fraction(point), tolerance
            );
        }
    }
}

} // namespace Kratos::Testing
