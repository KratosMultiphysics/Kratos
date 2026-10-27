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
#include <string>
#include <utility>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_modelers/nurbs_geometry_modeler.h"
#include "geometries/nurbs_surface_geometry.h"
#include "iga_application_variables.h"
#include "includes/kratos_components.h"

namespace Kratos::Testing
{

// Background patch [0, 1.5] x [0, 2] with 3 x 2 knot spans of sizes 0.5 x 1
// and degrees (2, 3), so any mix-up of u and v fails. Checks the elements on
// all spans with their knot span indices and, for every facet on the 3
// interior knot lines, its geometry, the knot span indices and sizes of each
// side, and that each side is the one-sided limit of the B-spline basis. A
// last entry creates more conditions on the same skeleton, sharing geometries.
KRATOS_TEST_CASE_IN_SUITE(IgaModelerGsbmElementsAndFacets, KratosIgaFastSuite)
{
    using IndexType = std::size_t;
    using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;

    Model model;
    model.CreateModelPart("IgaModelPart");

    NurbsGeometryModeler nurbs_modeler(model, Parameters(R"({
        "model_part_name"     : "IgaModelPart",
        "geometry_name"       : "background",
        "lower_point_xyz"     : [0.0, 0.0, 0.0],
        "upper_point_xyz"     : [1.5, 2.0, 0.0],
        "lower_point_uvw"     : [0.0, 0.0, 0.0],
        "upper_point_uvw"     : [1.5, 2.0, 0.0],
        "polynomial_order"    : [2, 3],
        "number_of_knot_spans": [3, 2]
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
        }, {
            "type"                            : "condition",
            "name"                            : "GenericCondition",
            "iga_model_part"                  : "SharedFacets",
            "gsbm_parameters"                 : {"skeleton_model_part": "Facets"}
        }]
    })"));
    p_gsbm_modeler->SetupModelPart();

    const double tolerance = 1.0e-10;

    // Elements: 3 x 2 spans with (2 + 1) x (3 + 1) Gauss points, covering the
    // area of the patch
    const ModelPart& r_domain = model.GetModelPart("IgaModelPart.Domain");
    KRATOS_EXPECT_EQ(r_domain.NumberOfElements(), 72);
    double area = 0.0;
    for (const auto& r_element : r_domain.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        area += r_geometry.IntegrationPoints()[0].Weight() * r_geometry.DeterminantOfJacobian(0);

        // Knot span indices: the span containing the Gauss point
        const auto& r_knot_span_indices = r_element.GetValue(KNOT_SPAN_INDICES);
        const array_1d<double, 3> gauss_point = r_geometry.Center().Coordinates();
        KRATOS_EXPECT_EQ(r_knot_span_indices.size(), 2);
        KRATOS_EXPECT_EQ(r_knot_span_indices[0], static_cast<IndexType>(gauss_point[0] / 0.5));
        KRATOS_EXPECT_EQ(r_knot_span_indices[1], static_cast<IndexType>(gauss_point[1] / 1.0));
    }
    KRATOS_EXPECT_NEAR(area, 3.0, tolerance);

    // Facets: 2 vertical lines x 2 spans x (3 + 1) points and 1 horizontal
    // line x 3 spans x (2 + 1) points
    const auto& r_surface = dynamic_cast<const NurbsSurfaceType&>(
        model.GetModelPart("IgaModelPart").GetGeometry("background")
    );
    const ModelPart& r_facets = model.GetModelPart("IgaModelPart.Facets");
    KRATOS_EXPECT_EQ(r_facets.NumberOfConditions(), 25);

    // Reference basis evaluated just off the knot line
    NurbsSurfaceShapeFunction reference(
        r_surface.PolynomialDegreeU(), r_surface.PolynomialDegreeV(), 3
    );
    const double offset = 1.0e-9;
    const double limit_tolerance = 1.0e-6;
    const IndexType point_index = 0;

    IndexType number_of_vertical_facets = 0;
    IndexType number_of_horizontal_facets = 0;
    double length = 0.0;
    for (const auto& r_facet : r_facets.Conditions()) {
        const auto& r_minus = r_facet.GetGeometry().GetGeometryPart(0);
        const auto& r_plus = r_facet.GetGeometry().GetGeometryPart(1);
        length += r_minus.IntegrationPoints()[0].Weight() * r_minus.DeterminantOfJacobian(0);

        // Both sides are the same point on a knot line, with the normal
        // pointing from the - to the + side
        const array_1d<double, 3> center = r_minus.Center().Coordinates();
        const array_1d<double, 3> center_plus = r_plus.Center().Coordinates();
        KRATOS_EXPECT_VECTOR_NEAR(center_plus, center, tolerance);
        const array_1d<double, 3> normal = r_minus.UnitNormal(point_index);
        const array_1d<double, 3> normal_plus = r_plus.UnitNormal(point_index);
        KRATOS_EXPECT_VECTOR_NEAR(normal_plus, normal, tolerance);
        const bool is_vertical = std::abs(normal[0] - 1.0) < tolerance;
        if (is_vertical) {
            ++number_of_vertical_facets;
            KRATOS_EXPECT_TRUE(
                std::abs(center[0] - 0.5) < tolerance || std::abs(center[0] - 1.0) < tolerance
            );
        } else {
            ++number_of_horizontal_facets;
            KRATOS_EXPECT_NEAR(normal[1], 1.0, tolerance);
            KRATOS_EXPECT_NEAR(center[1], 1.0, tolerance);
        }

        for (const double side : {-1.0, 1.0}) {
            const auto& r_side = side < 0.0 ? r_minus : r_plus;
            const auto& r_point = r_side.IntegrationPoints()[0];

            // Knot span sizes of the side's span, in u and v
            const Vector& r_knot_span_sizes = r_side.GetValue(KNOT_SPAN_SIZES);
            KRATOS_EXPECT_EQ(r_knot_span_sizes.size(), 2);
            KRATOS_EXPECT_NEAR(r_knot_span_sizes[0], 0.5, tolerance);
            KRATOS_EXPECT_NEAR(r_knot_span_sizes[1], 1.0, tolerance);

            // Knot span indices of the side's span: the one just off the knot
            // line, on that side
            const auto& r_knot_span_indices = r_side.GetValue(KNOT_SPAN_INDICES);
            const double x_off = center[0] + side * offset * normal[0];
            const double y_off = center[1] + side * offset * normal[1];
            KRATOS_EXPECT_EQ(r_knot_span_indices.size(), 2);
            KRATOS_EXPECT_EQ(r_knot_span_indices[0], static_cast<IndexType>(x_off / 0.5));
            KRATOS_EXPECT_EQ(r_knot_span_indices[1], static_cast<IndexType>(y_off / 1.0));
            reference.ComputeBSplineShapeFunctionValues(r_surface.KnotsU(), r_surface.KnotsV(),
                r_point[0] + side * offset * normal[0], r_point[1] + side * offset * normal[1]);

            // Same control points, in the same order
            const auto control_point_indices = reference.ControlPointIndices(
                r_surface.NumberOfControlPointsU(), r_surface.NumberOfControlPointsV());
            KRATOS_EXPECT_EQ(r_side.size(), control_point_indices.size());
            for (IndexType i = 0; i < r_side.size(); ++i) {
                KRATOS_EXPECT_EQ(
                    r_side[i].Id(), r_surface.pGetPoint(control_point_indices[i])->Id()
                );
            }

            // Same shape functions and derivatives up to order 3
            const auto method = r_side.GetDefaultIntegrationMethod();
            const Matrix& r_N = r_side.ShapeFunctionsValues();
            for (IndexType i = 0; i < r_side.size(); ++i) {
                KRATOS_EXPECT_NEAR(r_N(0, i), reference(i, 0), limit_tolerance);
            }
            for (IndexType order = 1; order <= 3; ++order) {
                const Matrix& r_derivatives =
                    r_side.ShapeFunctionDerivatives(order, point_index, method);
                for (IndexType i = 0; i < r_side.size(); ++i) {
                    for (IndexType k = 0; k <= order; ++k) {
                        const IndexType row =
                            NurbsSurfaceShapeFunction::IndexOfShapeFunctionRow(order - k, k);
                        KRATOS_EXPECT_NEAR(r_derivatives(i, k), reference(i, row), limit_tolerance);
                    }
                }
            }
        }
    }
    KRATOS_EXPECT_EQ(number_of_vertical_facets, 16);
    KRATOS_EXPECT_EQ(number_of_horizontal_facets, 9);
    KRATOS_EXPECT_NEAR(length, 2.0 * 2.0 + 1.0 * 1.5, tolerance);

    // Shared skeleton: a new condition per facet, on the same geometry
    const ModelPart& r_shared_facets = model.GetModelPart("IgaModelPart.SharedFacets");
    KRATOS_EXPECT_EQ(r_shared_facets.NumberOfConditions(), r_facets.NumberOfConditions());
    KRATOS_EXPECT_EQ(r_shared_facets.NumberOfNodes(), r_facets.NumberOfNodes());
    auto it_facet = r_facets.ConditionsBegin();
    for (const auto& r_condition : r_shared_facets.Conditions()) {
        KRATOS_EXPECT_NE(r_condition.Id(), it_facet->Id());
        KRATOS_EXPECT_EQ(&r_condition.GetGeometry(), &it_facet->GetGeometry());
        ++it_facet;
    }
}

// Three patches on [0, 1] x [0, 0.75] with 4 x 3 knot spans and degrees
// (3, 2), (2, 3) and (2, 2), one per field as for divergence-conforming
// Stokes. An entry with "geometry_names" creates one element per Gauss point
// of a common rule, 4 x 4 per span, enough for degree 3: each geometry has a
// part per patch, in the given order and at the same point, and the knot span
// indices of that point. The elements' model part holds every patch's nodes.
KRATOS_TEST_CASE_IN_SUITE(IgaModelerGsbmElementsMixedSpaces, KratosIgaFastSuite)
{
    using IndexType = std::size_t;

    Model model;
    const ModelPart& r_iga_model_part = model.CreateModelPart("IgaModelPart");

    const std::array<std::pair<std::string, std::array<int, 2>>, 3> patches{{
        {"VelocityX", {3, 2}}, {"VelocityY", {2, 3}}, {"Pressure", {2, 2}}
    }};
    for (const auto& [r_name, r_degrees] : patches) {
        Parameters parameters(R"({
            "model_part_name"     : "IgaModelPart",
            "geometry_name"       : "",
            "lower_point_xyz"     : [0.0, 0.0, 0.0],
            "upper_point_xyz"     : [1.0, 0.75, 0.0],
            "lower_point_uvw"     : [0.0, 0.0, 0.0],
            "upper_point_uvw"     : [1.0, 0.75, 0.0],
            "polynomial_order"    : [0, 0],
            "number_of_knot_spans": [4, 3]
        })");
        parameters["geometry_name"].SetString(r_name);
        parameters["polynomial_order"][0].SetInt(r_degrees[0]);
        parameters["polynomial_order"][1].SetInt(r_degrees[1]);
        NurbsGeometryModeler(model, parameters).SetupGeometryModel();
    }

    const Modeler& r_registered_modeler = KratosComponents<Modeler>::Get("IgaModelerGsbm");
    r_registered_modeler.Create(model, Parameters(R"({
        "analysis_model_part_name": "IgaModelPart",
        "geometry_name"           : "Pressure",
        "element_condition_list"  : [{
            "type"                            : "element",
            "name"                            : "GenericElement",
            "iga_model_part"                  : "Domain",
            "shape_function_derivatives_order": 2,
            "geometry_names"                  : ["VelocityX", "VelocityY", "Pressure"]
        }]
    })"))->SetupModelPart();

    // 4 x 3 knot spans with 4 x 4 Gauss points each
    const ModelPart& r_domain = model.GetModelPart("IgaModelPart.Domain");
    KRATOS_EXPECT_EQ(r_domain.NumberOfElements(), 192);

    // Control points per patch: (4 + p_u) x (3 + p_v)
    KRATOS_EXPECT_EQ(r_domain.NumberOfNodes(), 35 + 36 + 30);

    const double tolerance = 1.0e-12;
    for (const auto& r_element : r_domain.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        KRATOS_EXPECT_EQ(r_geometry.NumberOfGeometryParts(), 3);
        const array_1d<double, 3> center = r_geometry.GetGeometryPart(0).Center().Coordinates();

        for (IndexType part = 0; part < 3; ++part) {
            // A quadrature point of the part's patch, at the same point, with
            // its (p_u + 1) x (p_v + 1) nonzero basis functions
            const auto& r_part = r_geometry.GetGeometryPart(part);
            const auto& [r_name, r_degrees] = patches[part];
            KRATOS_EXPECT_EQ(&r_part.GetGeometryParent(0), &r_iga_model_part.GetGeometry(r_name));
            KRATOS_EXPECT_EQ(r_part.size(), (r_degrees[0] + 1) * (r_degrees[1] + 1));
            KRATOS_EXPECT_VECTOR_NEAR(r_part.Center().Coordinates(), center, tolerance);
        }

        // Knot span of the point, on the element's own geometry
        const auto& r_indices = r_geometry.GetValue(KNOT_SPAN_INDICES);
        KRATOS_EXPECT_EQ(r_indices.size(), 2);
        KRATOS_EXPECT_EQ(r_indices[0], static_cast<IndexType>(center[0] / 0.25));
        KRATOS_EXPECT_EQ(r_indices[1], static_cast<IndexType>(center[1] / 0.25));
    }
}

} // namespace Kratos::Testing
