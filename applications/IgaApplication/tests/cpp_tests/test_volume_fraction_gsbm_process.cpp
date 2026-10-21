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
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "geometries/nurbs_curve_geometry.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"
#include "includes/kratos_components.h"

namespace Kratos::Testing
{

namespace
{

using IndexType = std::size_t;
using NurbsCurveType = NurbsCurveGeometry<2, PointerVector<Node>>;

/// Appends a closed loop of linear NURBS curves, one per edge with its layer,
/// to the skin, as ImportNurbsSbmModeler does
void CreateNurbsSkinLoop(
    ModelPart& rSkinModelPart,
    const std::vector<std::array<double, 2>>& rVertices,
    const std::vector<std::string>& rLayers)
{
    const IndexType first_id = rSkinModelPart.NumberOfGeometries();
    const IndexType number_of_vertices = rVertices.size();
    Vector knots(2);
    knots[0] = 0.0;
    knots[1] = 1.0;

    for (IndexType i = 0; i < number_of_vertices; ++i) {
        const auto& r_start = rVertices[i];
        const auto& r_end = rVertices[(i + 1) % number_of_vertices];
        PointerVector<Node> control_points;
        control_points.push_back(Kratos::make_intrusive<Node>(0, r_start[0], r_start[1], 0.0));
        control_points.push_back(Kratos::make_intrusive<Node>(0, r_end[0], r_end[1], 0.0));

        auto p_curve = Kratos::make_shared<NurbsCurveType>(control_points, 1, knots);
        p_curve->SetId(first_id + i);
        p_curve->SetValue(IDENTIFIER, rLayers[i]);
        rSkinModelPart.AddGeometry(p_curve);
    }
}

/// Adds a counterclockwise ellipse to the skin, as one rational quadratic
/// curve with 9 control points: the affine image of the exact circle
void CreateNurbsEllipse(
    ModelPart& rSkinModelPart,
    const array_1d<double, 3>& rCenter,
    const std::array<double, 2>& rSemiAxes,
    const std::string& rLayer)
{
    const std::array<std::array<double, 2>, 9> directions{{
        {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}
    }};
    PointerVector<Node> control_points;
    for (const auto& r_direction : directions) {
        control_points.push_back(Kratos::make_intrusive<Node>(0,
            rCenter[0] + rSemiAxes[0] * r_direction[0],
            rCenter[1] + rSemiAxes[1] * r_direction[1], 0.0
        ));
    }

    const double w = std::sqrt(2.0) / 2.0;
    const std::array<double, 10> knot_values{0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1};
    const std::array<double, 9> weight_values{1, w, 1, w, 1, w, 1, w, 1};
    Vector knots(knot_values.size());
    std::copy(knot_values.begin(), knot_values.end(), knots.begin());
    Vector weights(weight_values.size());
    std::copy(weight_values.begin(), weight_values.end(), weights.begin());

    auto p_curve = Kratos::make_shared<NurbsCurveType>(control_points, 2, knots, weights);
    p_curve->SetId(rSkinModelPart.NumberOfGeometries());
    p_curve->SetValue(IDENTIFIER, rLayer);
    rSkinModelPart.AddGeometry(p_curve);
}

/// Closest point to rPoint on the boundary of the rectangle [rLower, rUpper]
array_1d<double, 3> ClosestPointOnRectangleBoundary(
    const array_1d<double, 3>& rPoint,
    const std::array<double, 2>& rLower,
    const std::array<double, 2>& rUpper)
{
    array_1d<double, 3> closest_point = rPoint;
    const bool is_inside = rLower[0] < rPoint[0] && rPoint[0] < rUpper[0]
        && rLower[1] < rPoint[1] && rPoint[1] < rUpper[1];

    // Outside: clamp onto the rectangle, a corner if outside in both directions
    if (!is_inside) {
        for (IndexType dir = 0; dir < 2; ++dir) {
            closest_point[dir] = std::clamp(rPoint[dir], rLower[dir], rUpper[dir]);
        }
        return closest_point;
    }

    // Inside: move to the nearest side
    double min_distance = std::numeric_limits<double>::max();
    for (IndexType dir = 0; dir < 2; ++dir) {
        for (const double bound : {rLower[dir], rUpper[dir]}) {
            if (std::abs(rPoint[dir] - bound) < min_distance) {
                min_distance = std::abs(rPoint[dir] - bound);
                closest_point = rPoint;
                closest_point[dir] = bound;
            }
        }
    }
    return closest_point;
}

} // namespace

// Background patch [0, 1] x [0, 0.75] with 4 x 3 knot spans of size 0.25 and
// degrees (2, 3), so any mix-up of u and v fails. The domain is the rectangle
// [0.125, 0.875] x [0.125, 0.625] with a square hole [0.3125, 0.4375]^2 inside
// span (1, 1), which gives volume fractions 0.25, 0.5, 0.75 and 1. Each element
// and facet side is checked against the span found from its center, shifted
// off the knot line towards its own side for facets. Cut facets must get a node
// at the closest skin point, with the outward normal and tangent of its side
// and zero curvature. Both loops are counterclockwise, as SBM expects. The
// skins are linear NURBS curves, tessellated exactly, with a layer per side of
// the rectangle and one for the hole, so the projection data is per layer.
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

    const std::vector<std::array<double, 2>> outer_vertices{
        {0.125, 0.125}, {0.875, 0.125}, {0.875, 0.625}, {0.125, 0.625}
    };
    const std::vector<std::array<double, 2>> inner_vertices{
        {0.3125, 0.3125}, {0.4375, 0.3125}, {0.4375, 0.4375}, {0.3125, 0.4375}
    };
    const std::vector<std::string> outer_layers{"bottom", "right", "top", "left"};
    CreateNurbsSkinLoop(model.CreateModelPart("OuterSkin"), outer_vertices, outer_layers);
    CreateNurbsSkinLoop(
        model.CreateModelPart("InnerSkin"), inner_vertices, {"hole", "hole", "hole", "hole"}
    );

    VolumeFractionGsbmProcess gsbm_process(model, Parameters(R"({
        "model_part_name"                   : "IgaModelPart",
        "geometry_name"                     : "background",
        "elements_model_part_name"          : "IgaModelPart.Domain",
        "skeleton_model_part_name"          : "IgaModelPart.Facets",
        "skin_model_part_outer_initial_name": "OuterSkin",
        "skin_model_part_inner_initial_name": "InnerSkin"
    })"));
    gsbm_process.ExecuteInitialize();

    // Straight curves: one segment each, sharing the corner nodes
    for (const std::string name : {"OuterSkin.tessellation", "InnerSkin.tessellation"}) {
        KRATOS_EXPECT_EQ(model.GetModelPart(name).NumberOfConditions(), 4);
        KRATOS_EXPECT_EQ(model.GetModelPart(name).NumberOfNodes(), 4);
    }

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

    // Projection data of a node, per layer, checking that the node lies on its
    // layer's side of the rectangle
    const auto find_projection_data = [&](const Node& rNode, const bool IsInner)
        -> const ModelPart* {
        if (IsInner) {
            return &model.GetModelPart("InnerSkin.hole.projection_data");
        }
        const std::array<IndexType, 4> side_directions{1, 0, 1, 0};
        const std::array<double, 4> side_coordinates{0.125, 0.875, 0.625, 0.125};
        const ModelPart* p_projection_data = nullptr;
        for (IndexType side = 0; side < 4; ++side) {
            const std::string name = "OuterSkin." + outer_layers[side] + ".projection_data";
            if (model.HasModelPart(name) && model.GetModelPart(name).HasNode(rNode.Id())) {
                KRATOS_EXPECT_EQ(p_projection_data, nullptr);  // in one layer only
                p_projection_data = &model.GetModelPart(name);
                KRATOS_EXPECT_NEAR(
                    rNode[side_directions[side]], side_coordinates[side], tolerance
                );
            }
        }
        return p_projection_data;
    };

    // Facets: 3 vertical lines x 3 spans x (3 + 1) points and 2 horizontal
    // lines x 4 spans x (2 + 1) points
    const ModelPart& r_facets = model.GetModelPart("IgaModelPart.Facets");
    KRATOS_EXPECT_EQ(r_facets.NumberOfConditions(), 60);
    const std::array<double, 2> outer_lower{0.125, 0.125};
    const std::array<double, 2> outer_upper{0.875, 0.625};
    const std::array<double, 2> inner_lower{0.3125, 0.3125};
    const std::array<double, 2> inner_upper{0.4375, 0.4375};
    const double offset = 1.0e-6;
    const IndexType point_index = 0;
    IndexType number_of_cut_facets = 0;
    IndexType number_of_projection_nodes = 0;
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

        // Uncut facets ([[alpha]] = 0) get no closest skin point
        const auto& r_minus = r_facet.GetGeometry().GetGeometryPart(0);
        const auto& r_plus = r_facet.GetGeometry().GetGeometryPart(1);
        if (r_minus.GetValue(VOLUME_FRACTION) == r_plus.GetValue(VOLUME_FRACTION)) {
            KRATOS_EXPECT_FALSE(r_facet.Has(NEIGHBOUR_NODES));
            continue;
        }
        ++number_of_cut_facets;

        // Cut facets: one node at the closest point of the two skins
        const array_1d<double, 3> center = r_minus.Center().Coordinates();
        const array_1d<double, 3> outer_point =
            ClosestPointOnRectangleBoundary(center, outer_lower, outer_upper);
        const array_1d<double, 3> inner_point =
            ClosestPointOnRectangleBoundary(center, inner_lower, inner_upper);
        const bool is_inner = norm_2(inner_point - center) < norm_2(outer_point - center);
        const array_1d<double, 3>& r_expected_point = is_inner ? inner_point : outer_point;

        const auto& r_projection_nodes = r_facet.GetValue(NEIGHBOUR_NODES);
        KRATOS_EXPECT_EQ(r_projection_nodes.size(), 1);
        const Node& r_projection_node = r_projection_nodes[0];
        KRATOS_EXPECT_VECTOR_NEAR(r_projection_node.Coordinates(), r_expected_point, tolerance);
        KRATOS_EXPECT_EQ(r_facet.GetValue(IDENTIFIER), is_inner ? "inner" : "outer");
        const ModelPart* p_projections = find_projection_data(r_projection_node, is_inner);
        KRATOS_EXPECT_NE(p_projections, nullptr);
        KRATOS_EXPECT_TRUE(p_projections->HasNode(r_projection_node.Id()));
        ++number_of_projection_nodes;

        // Unit axis normal out of the domain at the node's side: away from the
        // rectangle, into the hole
        const array_1d<double, 3>& r_normal = r_projection_node.GetValue(NORMAL);
        const IndexType dir = std::abs(r_normal[0]) > 0.5 ? 0 : 1;
        KRATOS_EXPECT_NEAR(std::abs(r_normal[dir]), 1.0, tolerance);
        KRATOS_EXPECT_NEAR(r_normal[1 - dir], 0.0, tolerance);
        const bool is_upper_side = (r_normal[dir] > 0.0) != is_inner;
        const auto& r_side_bounds = is_upper_side
            ? (is_inner ? inner_upper : outer_upper)
            : (is_inner ? inner_lower : outer_lower);
        KRATOS_EXPECT_NEAR(r_projection_node[dir], r_side_bounds[dir], tolerance);

        // Tangent with the domain on its left, and straight sides
        const array_1d<double, 3>& r_tangent = r_projection_node.GetValue(LOCAL_TANGENT);
        KRATOS_EXPECT_NEAR(r_tangent[0], -r_normal[1], tolerance);
        KRATOS_EXPECT_NEAR(r_tangent[1], r_normal[0], tolerance);
        KRATOS_EXPECT_NEAR(r_projection_node.GetValue(CURVATURE), 0.0, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(
            r_projection_node.GetValue(CURVATURE_DERIVATIVES), ZeroVector(2), tolerance
        );
    }

    // All facets are cut except the 8 on u = 0.5 in the bottom and top rows of
    // spans, where the volume fraction is 0.5 on both sides
    KRATOS_EXPECT_EQ(number_of_cut_facets, 52);
    KRATOS_EXPECT_EQ(number_of_projection_nodes, 52);
}

// Ellipse with semi-axes 0.35 and 0.25 and a circular hole of radius 0.1, both
// centred in [0, 1]^2 with 8 x 8 spans of degree 2, as exact rational curves
// given counterclockwise. The volume fractions must add up to the area between
// them.
// The closest points must be exact, and carry the outward normal, the
// curvature and its arc-length derivatives of the skin there.
KRATOS_TEST_CASE_IN_SUITE(VolumeFractionGsbmProcessNurbsEllipse, KratosIgaFastSuite)
{
    Model model;
    model.CreateModelPart("IgaModelPart");

    NurbsGeometryModeler nurbs_modeler(model, Parameters(R"({
        "model_part_name"     : "IgaModelPart",
        "geometry_name"       : "background",
        "lower_point_xyz"     : [0.0, 0.0, 0.0],
        "upper_point_xyz"     : [1.0, 1.0, 0.0],
        "lower_point_uvw"     : [0.0, 0.0, 0.0],
        "upper_point_uvw"     : [1.0, 1.0, 0.0],
        "polynomial_order"    : [2, 2],
        "number_of_knot_spans": [8, 8]
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
            "shape_function_derivatives_order": 2
        }, {
            "type"                            : "condition",
            "name"                            : "GenericCondition",
            "iga_model_part"                  : "Facets",
            "shape_function_derivatives_order": 3,
            "gsbm_parameters"                 : {}
        }]
    })"));
    p_gsbm_modeler->SetupModelPart();

    const array_1d<double, 3> center{0.5, 0.5, 0.0};
    const std::array<double, 2> ellipse_axes{0.35, 0.25};
    const std::array<double, 2> hole_axes{0.1, 0.1};
    CreateNurbsEllipse(model.CreateModelPart("OuterSkin"), center, ellipse_axes, "ellipse");
    CreateNurbsEllipse(model.CreateModelPart("InnerSkin"), center, hole_axes, "hole");

    VolumeFractionGsbmProcess gsbm_process(model, Parameters(R"({
        "model_part_name"                   : "IgaModelPart",
        "geometry_name"                     : "background",
        "elements_model_part_name"          : "IgaModelPart.Domain",
        "skeleton_model_part_name"          : "IgaModelPart.Facets",
        "skin_model_part_outer_initial_name": "OuterSkin",
        "skin_model_part_inner_initial_name": "InnerSkin"
    })"));
    gsbm_process.ExecuteInitialize();

    // Area between the skins, from the volume fractions at the Gauss points:
    // the inscribed polygons miss about their perimeter times the tolerance
    double area = 0.0;
    for (const auto& r_element : model.GetModelPart("IgaModelPart.Domain").Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        area += r_element.GetValue(VOLUME_FRACTION)
            * r_geometry.IntegrationPoints()[0].Weight() * r_geometry.DeterminantOfJacobian(0);
    }
    const double expected_area = Globals::Pi
        * (ellipse_axes[0] * ellipse_axes[1] - hole_axes[0] * hole_axes[1]);
    KRATOS_EXPECT_NEAR(area, expected_area, 1.0e-7);

    const double tolerance = 1.0e-10;
    IndexType number_of_cut_facets = 0;
    for (const auto& r_facet : model.GetModelPart("IgaModelPart.Facets").Conditions()) {
        if (!r_facet.Has(NEIGHBOUR_NODES)) {
            continue;
        }
        ++number_of_cut_facets;
        const bool is_inner = r_facet.GetValue(IDENTIFIER) == "inner";
        const Node& r_node = r_facet.GetValue(NEIGHBOUR_NODES)[0];
        const ModelPart& r_projections = model.GetModelPart(
            is_inner ? "InnerSkin.hole.projection_data" : "OuterSkin.ellipse.projection_data"
        );
        KRATOS_EXPECT_TRUE(r_projections.HasNode(r_node.Id()));

        // On the skin, at (a cos \theta, b sin \theta) from the center
        const double a = is_inner ? hole_axes[0] : ellipse_axes[0];
        const double b = is_inner ? hole_axes[1] : ellipse_axes[1];
        const double cos_theta = (r_node.X() - center[0]) / a;
        const double sin_theta = (r_node.Y() - center[1]) / b;
        KRATOS_EXPECT_NEAR(cos_theta * cos_theta + sin_theta * sin_theta, 1.0, tolerance);

        // Normal along the gradient of the ellipse's equation, into the hole,
        // and the facet point along it
        const double sign = is_inner ? -1.0 : 1.0;  // +1 if \tau is counterclockwise
        array_1d<double, 3> gradient = ZeroVector(3);
        gradient[0] = cos_theta / a;
        gradient[1] = sin_theta / b;
        const array_1d<double, 3> expected_normal = sign * gradient / norm_2(gradient);
        const array_1d<double, 3>& r_normal = r_node.GetValue(NORMAL);
        const array_1d<double, 3>& r_tangent = r_node.GetValue(LOCAL_TANGENT);
        KRATOS_EXPECT_VECTOR_NEAR(r_normal, expected_normal, tolerance);
        KRATOS_EXPECT_NEAR(r_tangent[0], -r_normal[1], tolerance);
        KRATOS_EXPECT_NEAR(r_tangent[1], r_normal[0], tolerance);
        const array_1d<double, 3> point =
            r_facet.GetGeometry().GetGeometryPart(0).Center().Coordinates();
        KRATOS_EXPECT_NEAR(inner_prod(point - r_node.Coordinates(), r_tangent), 0.0, tolerance);

        // Counterclockwise: \kappa = ab / D^(3/2), D = a^2 sin^2 + b^2 cos^2,
        // d/ds = D^(-1/2) d/d\theta; clockwise, \kappa and \kappa'' flip
        const double d = a * a * sin_theta * sin_theta + b * b * cos_theta * cos_theta;
        const double d_theta = 2.0 * (a * a - b * b) * sin_theta * cos_theta;
        const double d_theta_theta =
            2.0 * (a * a - b * b) * (cos_theta * cos_theta - sin_theta * sin_theta);
        const double kappa = a * b / std::pow(d, 1.5);
        const double kappa_s = -1.5 * a * b * d_theta / std::pow(d, 3);
        const double kappa_ss = -1.5 * a * b / std::sqrt(d)
            * (d_theta_theta / std::pow(d, 3) - 3.0 * d_theta * d_theta / std::pow(d, 4));
        const Vector& r_curvature_derivatives = r_node.GetValue(CURVATURE_DERIVATIVES);
        KRATOS_EXPECT_NEAR(r_node.GetValue(CURVATURE), sign * kappa, 1.0e-8);
        KRATOS_EXPECT_NEAR(r_curvature_derivatives[0], kappa_s, 1.0e-8);
        KRATOS_EXPECT_NEAR(r_curvature_derivatives[1], sign * kappa_ss, 1.0e-8);
    }

    // One node per cut facet, on either skin
    KRATOS_EXPECT_GT(number_of_cut_facets, 0);
    KRATOS_EXPECT_EQ(
        model.GetModelPart("OuterSkin.ellipse.projection_data").NumberOfNodes()
            + model.GetModelPart("InnerSkin.hole.projection_data").NumberOfNodes(),
        number_of_cut_facets
    );
}

} // namespace Kratos::Testing
