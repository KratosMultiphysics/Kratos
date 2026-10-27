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
#include <iterator>
#include <string>
#include <vector>

// External includes

// Project includes
#include "custom_modelers/iga_modeler_gsbm.h"
#include "geometries/coupling_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "iga_application_variables.h"
#include "integration/integration_point_utilities.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{

namespace
{

using IndexType = IgaModelerGsbm::IndexType;
using SizeType = IgaModelerGsbm::SizeType;
using GeometryType = IgaModelerGsbm::GeometryType;
using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;

/**
 * @brief Creates a one-sided quadrature point on a line facet of a 2D patch,
 * storing in it the control points and the basis values and derivatives of
 * the given knot span, which is specific to the facet side (+ or -).
 *
 * @param rSurface The background B-spline surface.
 * @param rShapeFunctionContainer Container reused for all points, sized for
 *        NumberOfShapeFunctionDerivatives.
 * @param NumberOfShapeFunctionDerivatives Shape function derivatives are
 *        stored up to order this value - 1.
 * @param rSpans Knot span indices in u and v used to evaluate the basis.
 * @param rCoordinates Parameter coordinates (u, v) of the point.
 * @param Weight Integration weight of the point.
 * @param rTangent Facet tangent (u, v), which defines the Normal() of the
 *        point.
 * @return The quadrature point geometry, with the surface as parent.
 */
GeometryType::Pointer CreateOneSidedQuadraturePoint(
    NurbsSurfaceType& rSurface,
    NurbsSurfaceShapeFunction& rShapeFunctionContainer,
    const SizeType NumberOfShapeFunctionDerivatives,
    const std::array<IndexType, 2>& rSpans,
    const std::array<double, 2>& rCoordinates,
    const double Weight,
    const std::array<double, 2>& rTangent)
{
    // Fill container with basis values
    rShapeFunctionContainer.ComputeBSplineShapeFunctionValuesAtSpan(
        rSurface.KnotsU(), rSurface.KnotsV(),
        rSpans[0], rSpans[1],
        rCoordinates[0], rCoordinates[1]
    );

    // Control points supporting the chosen span
    const SizeType num_cps = rShapeFunctionContainer.NumberOfNonzeroControlPoints();
    GeometryType::PointsArrayType control_points(num_cps);

    const std::vector<int> control_point_indices = rShapeFunctionContainer.ControlPointIndices(
        rSurface.NumberOfControlPointsU(),
        rSurface.NumberOfControlPointsV()
    );

    for (IndexType i = 0; i < num_cps; ++i) {
        control_points(i) = rSurface.pGetPoint(control_point_indices[i]);
    }

    // Shape function values
    Matrix N(1, num_cps);

    for (IndexType i = 0; i < num_cps; ++i) {
        N(0, i) = rShapeFunctionContainer(i, 0);
    }

    // Shape function derivatives
    DenseVector<Matrix> shape_function_derivatives(
        NumberOfShapeFunctionDerivatives - 1
    );

    IndexType shape_derivative_index = 1;
    for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n) {
        shape_function_derivatives[n].resize(num_cps, n + 2);
        for (IndexType k = 0; k < n + 2; ++k) {
            for (IndexType i = 0; i < num_cps; ++i) {
                shape_function_derivatives[n](i, k) =
                    rShapeFunctionContainer(i, shape_derivative_index + k);
            }
        }
        shape_derivative_index += n + 2;
    }

    // Integration point (u, v) with its weight, and the basis evaluated there
    GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
        rSurface.GetDefaultIntegrationMethod(),
        IntegrationPoint<3>(rCoordinates[0], rCoordinates[1], 0.0, Weight),
        N, shape_function_derivatives
    );

    return CreateQuadraturePointsUtility<Node>::CreateQuadraturePointCurveOnSurface(
        data_container, control_points, rTangent[0], rTangent[1], &rSurface
    );
}

/**
 * @brief Finds the knot span of a quadrature point, per direction.
 *
 * @param rSpans Knot span boundaries per direction.
 * @param rQuadraturePoint The quadrature point geometry.
 * @return The index of the knot span, from 0, per direction.
 */
std::vector<IndexType> FindKnotSpanIndices(
    const std::vector<std::vector<double>>& rSpans,
    const GeometryType& rQuadraturePoint)
{
    const auto& r_point = rQuadraturePoint.IntegrationPoints()[0];
    std::vector<IndexType> knot_span_indices(rSpans.size());

    for (IndexType dir = 0; dir < rSpans.size(); ++dir) {
        knot_span_indices[dir] =
            std::upper_bound(rSpans[dir].begin(), rSpans[dir].end(), r_point[dir])
            - rSpans[dir].begin() - 1;
    }

    return knot_span_indices;
}

} // namespace

void IgaModelerGsbm::SetupModelPart()
{
    KRATOS_TRY

    for (const std::string key : {"analysis_model_part_name", "geometry_name"}) {
        KRATOS_ERROR_IF(mParameters[key].GetString().empty()) << "IgaModelerGsbm: \"" << key
            << "\" must be provided." << std::endl;
    }

    ModelPart& r_analysis_model_part = mpModel->GetModelPart(
        mParameters["analysis_model_part_name"].GetString()
    );

    const Parameters iga_physics_parameters = mParameters["element_condition_list"];
    KRATOS_ERROR_IF(iga_physics_parameters.size() == 0)
        << "IgaModelerGsbm: \"element_condition_list\" is empty." << std::endl;

    CreateIntegrationDomain(r_analysis_model_part, iga_physics_parameters);

    KRATOS_CATCH("")
}

void IgaModelerGsbm::CreateIntegrationDomain(
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    for (IndexType i = 0; i < rPhysicsParameters.size(); ++i) {
        CreateIntegrationDomainPerUnit(rModelPart, rPhysicsParameters[i]);
    }
}

void IgaModelerGsbm::CreateIntegrationDomainPerUnit(
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    const std::array<std::string, 3> required_keys{"type", "name", "iga_model_part"};
    for (const std::string& r_key : required_keys) {
        KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has(r_key)) << "IgaModelerGsbm: \"" << r_key
            << "\" needs to be specified in every entry of \"element_condition_list\"."
            << std::endl;
    }

    // Entries on the skeleton of an earlier entry create no quadrature points
    const bool reuses_skeleton = rPhysicsParameters.Has("gsbm_parameters")
        && rPhysicsParameters["gsbm_parameters"].Has("skeleton_model_part");
    KRATOS_ERROR_IF(!reuses_skeleton && !rPhysicsParameters.Has("shape_function_derivatives_order"))
        << "IgaModelerGsbm: \"shape_function_derivatives_order\" needs to be specified in "
        << "every entry of \"element_condition_list\" that creates quadrature points."
        << std::endl;

    const std::string sub_model_part_name = rPhysicsParameters["iga_model_part"].GetString();
    ModelPart& r_sub_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
        ? rModelPart.GetSubModelPart(sub_model_part_name)
        : rModelPart.CreateSubModelPart(sub_model_part_name);

    GeometriesArrayType geometry_list;
    GetGeometryList(geometry_list, rModelPart, rPhysicsParameters);

    if (rPhysicsParameters.Has("geometry_names"))
        // Elements with one space per field: one patch per field, same knots
        CreateQuadraturePointGeometriesMixed(geometry_list, r_sub_model_part, rPhysicsParameters);
    else if (!rPhysicsParameters.Has("gsbm_parameters"))
        // Elements on the patch, or body-fitted conditions on one or more BReps
        CreateQuadraturePointGeometries(geometry_list, r_sub_model_part, rPhysicsParameters);
    else if (!reuses_skeleton)
        // GSBM facets: always a single geometry, the background patch
        CreateQuadraturePointGeometriesGsbm(
            geometry_list[0], r_sub_model_part, rPhysicsParameters
        );
    else
        // More GSBM conditions on the skeleton of an earlier entry
        CreateConditionsOnGsbmSkeleton(geometry_list, r_sub_model_part, rPhysicsParameters);
}

void IgaModelerGsbm::GetGeometryList(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    const std::string type = rPhysicsParameters["type"].GetString();

    if (rPhysicsParameters.Has("gsbm_parameters")
        && rPhysicsParameters["gsbm_parameters"].Has("skeleton_model_part")) {
        // GSBM conditions on an existing skeleton: the earlier entry's facets
        const std::string skeleton_name =
            rPhysicsParameters["gsbm_parameters"]["skeleton_model_part"].GetString();
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart(skeleton_name)
            && rModelPart.GetSubModelPart(skeleton_name).NumberOfConditions() > 0)
            << "IgaModelerGsbm: no skeleton found in \"" << skeleton_name << "\". The entry "
            << "that creates it must come earlier in \"element_condition_list\"."
            << std::endl;
        for (auto& r_facet : rModelPart.GetSubModelPart(skeleton_name).Conditions()) {
            rGeometryList.push_back(r_facet.pGetGeometry());
        }
    } else if (type == "element" && rPhysicsParameters.Has("geometry_names")) {
        // Elements with one space per field: the patch of each field
        const Parameters geometry_names = rPhysicsParameters["geometry_names"];
        for (IndexType i = 0; i < geometry_names.size(); ++i) {
            const std::string geometry_name = geometry_names[i].GetString();
            KRATOS_ERROR_IF_NOT(rModelPart.HasGeometry(geometry_name))
                << "IgaModelerGsbm: geometry \"" << geometry_name << "\" not found in "
                << rModelPart.FullName() << "." << std::endl;
            rGeometryList.push_back(rModelPart.pGetGeometry(geometry_name));
        }
    } else if (type == "element" || rPhysicsParameters.Has("gsbm_parameters")) {
        // Elements and GSBM facets: the background patch
        const std::string& r_geometry_name = mParameters["geometry_name"].GetString();
        KRATOS_ERROR_IF_NOT(rModelPart.HasGeometry(r_geometry_name))
            << "IgaModelerGsbm: geometry \"" << r_geometry_name << "\" not found in "
            << rModelPart.FullName() << "." << std::endl;
        rGeometryList.push_back(rModelPart.pGetGeometry(r_geometry_name));
    } else if (type == "condition") {
        // Body-fitted conditions: the given boundary geometries of the patch
        KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("brep_ids"))
            << "IgaModelerGsbm: missing \"brep_ids\" in a body-fitted condition entry."
            << std::endl;
        for (IndexType i = 0; i < rPhysicsParameters["brep_ids"].size(); ++i) {
            const IndexType brep_id = rPhysicsParameters["brep_ids"][i].GetInt();
            rGeometryList.push_back(rModelPart.pGetGeometry(brep_id));
        }
    } else {
        KRATOS_ERROR << "IgaModelerGsbm: type \"" << type << "\" not defined. "
            << "Available types are \"element\" and \"condition\"." << std::endl;
    }
}

void IgaModelerGsbm::CreateQuadraturePointGeometries(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters) const
{
    const std::string type = rParameters["type"].GetString();
    const std::string name = rParameters["name"].GetString();

    const int shape_function_derivatives_order =
        rParameters["shape_function_derivatives_order"].GetInt();
    KRATOS_ERROR_IF(shape_function_derivatives_order < 2) << "IgaModelerGsbm: "
        << "shape_function_derivatives_order must be at least 2." << std::endl;

    for (auto& r_geometry : rGeometryList) {
        const SizeType dimension = r_geometry.LocalSpaceDimension();

        // Gauss points per span and direction
        IntegrationInfo integration_info = r_geometry.GetDefaultIntegrationInfo();
        if (rParameters.Has("number_of_integration_points_per_span")) {
            for (IndexType dir = 0; dir < dimension; ++dir) {
                const int num = rParameters["number_of_integration_points_per_span"].GetInt();
                integration_info.SetNumberOfIntegrationPointsPerSpan(dir, num);
            }
        }

        // Quadrature points of every knot span
        GeometriesArrayType geometries;
        r_geometry.CreateQuadraturePointGeometries(
            geometries, shape_function_derivatives_order, integration_info
        );

        if (type == "element") {
            // Knot span boundaries
            std::vector<std::vector<double>> spans(dimension);
            for (IndexType dir = 0; dir < dimension; ++dir) {
                r_geometry.SpansLocalSpace(spans[dir], dir);
            }

            for (auto& r_quadrature_point : geometries) {
                // Index of the Gauss point's knot span, per direction
                r_quadrature_point.SetValue(
                    KNOT_SPAN_INDICES, FindKnotSpanIndices(spans, r_quadrature_point)
                );
            }

            CreateElements(geometries.ptr_begin(), geometries.ptr_end(), rModelPart, name);
            rModelPart.AddNodes(r_geometry.Points().ptr_begin(), r_geometry.Points().ptr_end());
        } else if (type == "condition") {
            const ModelPart& r_parent_model_part = rModelPart.GetParentModelPart();
            KRATOS_ERROR_IF_NOT(r_parent_model_part.Has(KNOT_SPAN_SIZES))
                << "IgaModelerGsbm: KNOT_SPAN_SIZES not found in "
                << r_parent_model_part.FullName() << "." << std::endl;

            for (auto& r_quadrature_point : geometries) {
                // Uniform knot span sizes of the patch (h of Nitsche's penalty)
                r_quadrature_point.SetValue(
                    KNOT_SPAN_SIZES, r_parent_model_part.GetValue(KNOT_SPAN_SIZES)
                );
            }

            CreateConditions(geometries.ptr_begin(), geometries.ptr_end(), rModelPart, name);
        } else {
            KRATOS_ERROR << "IgaModelerGsbm: type \"" << type << "\" not defined. "
                << "Available types are \"element\" and \"condition\"." << std::endl;
        }

        KRATOS_INFO_IF("IgaModelerGsbm", mEchoLevel > 0)
            << "Created " << geometries.size() << " " << type << "s of type " << name
            << " in " << rModelPart.FullName() << "." << std::endl;
    }
}

void IgaModelerGsbm::CreateQuadraturePointGeometriesMixed(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters) const
{
    const std::string name = rParameters["name"].GetString();
    KRATOS_ERROR_IF(rParameters["type"].GetString() != "element") << "IgaModelerGsbm: "
        << "\"geometry_names\" is only available for elements." << std::endl;

    const int shape_function_derivatives_order =
        rParameters["shape_function_derivatives_order"].GetInt();
    KRATOS_ERROR_IF(shape_function_derivatives_order < 2) << "IgaModelerGsbm: "
        << "shape_function_derivatives_order must be at least 2." << std::endl;

    // One 2D B-spline patch per field
    std::vector<NurbsSurfaceType*> surfaces;
    for (auto& r_geometry : rGeometryList) {
        const auto p_surface = dynamic_cast<NurbsSurfaceType*>(&r_geometry);
        KRATOS_ERROR_IF(p_surface == nullptr) << "IgaModelerGsbm: geometry #"
            << r_geometry.Id() << " of \"geometry_names\" is not a 2D NURBS surface."
            << std::endl;
        surfaces.push_back(p_surface);
    }

    // Knot span boundaries, the same in every patch
    std::vector<std::vector<double>> spans(2);
    for (IndexType dir = 0; dir < 2; ++dir) {
        surfaces[0]->SpansLocalSpace(spans[dir], dir);
        for (const auto p_surface : surfaces) {
            std::vector<double> surface_spans;
            p_surface->SpansLocalSpace(surface_spans, dir);
            const bool is_same = surface_spans.size() == spans[dir].size()
                && std::equal(surface_spans.begin(), surface_spans.end(), spans[dir].begin(),
                    [](const double A, const double B) { return std::abs(A - B) < 1.0e-12; });
            KRATOS_ERROR_IF_NOT(is_same) << "IgaModelerGsbm: the patches of "
                << "\"geometry_names\" must have the same knot spans." << std::endl;
        }
    }

    // Gauss points per span and direction: enough for the highest degree
    std::array<SizeType, 2> number_of_points{0, 0};
    for (IndexType dir = 0; dir < 2; ++dir) {
        for (const auto p_surface : surfaces) {
            number_of_points[dir] =
                std::max(number_of_points[dir], p_surface->PolynomialDegree(dir) + 1);
        }
        if (rParameters.Has("number_of_integration_points_per_span")) {
            number_of_points[dir] = rParameters["number_of_integration_points_per_span"].GetInt();
        }
    }

    // One Gauss rule for all patches, so their quadrature points coincide
    GeometryType::IntegrationPointsArrayType integration_points;
    surfaces[0]->CreateIntegrationPoints(
        integration_points, number_of_points[0], number_of_points[1]
    );

    // Quadrature points of every patch
    std::vector<GeometriesArrayType> quadrature_points(surfaces.size());
    for (IndexType k = 0; k < surfaces.size(); ++k) {
        IntegrationInfo integration_info = surfaces[k]->GetDefaultIntegrationInfo();
        surfaces[k]->CreateQuadraturePointGeometries(
            quadrature_points[k], shape_function_derivatives_order, integration_points,
            integration_info
        );
    }

    // Per Gauss point: one geometry part per field, in the given order
    GeometriesArrayType geometries;
    for (IndexType i = 0; i < integration_points.size(); ++i) {
        std::vector<GeometryType::Pointer> parts;
        for (const auto& r_quadrature_points : quadrature_points) {
            parts.push_back(r_quadrature_points(i));
        }
        auto p_geometry = Kratos::make_shared<CouplingGeometry<Node>>(parts);

        // Index of the Gauss point's knot span, per direction
        p_geometry->SetValue(KNOT_SPAN_INDICES, FindKnotSpanIndices(spans, *parts[0]));
        geometries.push_back(p_geometry);
    }

    CreateElements(geometries.ptr_begin(), geometries.ptr_end(), rModelPart, name);
    for (const auto p_surface : surfaces) {
        rModelPart.AddNodes(p_surface->Points().ptr_begin(), p_surface->Points().ptr_end());
    }

    KRATOS_INFO_IF("IgaModelerGsbm", mEchoLevel > 0)
        << "Created " << geometries.size() << " elements of type " << name << " on "
        << surfaces.size() << " patches in " << rModelPart.FullName() << "." << std::endl;
}

void IgaModelerGsbm::CreateQuadraturePointGeometriesGsbm(
    GeometryType& rBackground,
    ModelPart& rModelPart,
    const Parameters rParameters) const
{
    KRATOS_ERROR_IF(rParameters["type"].GetString() != "condition")
        << "IgaModelerGsbm: entries with \"gsbm_parameters\" create the facet "
        << "conditions, so they require type \"condition\"." << std::endl;

    const std::string name = rParameters["name"].GetString();

    const int shape_function_derivatives_order =
        rParameters["shape_function_derivatives_order"].GetInt();
    KRATOS_ERROR_IF(shape_function_derivatives_order < 2)
        << "IgaModelerGsbm: \"shape_function_derivatives_order\" must be at least 2, "
        << "so that the first derivatives needed for the normals and Jacobians "
        << "are stored." << std::endl;
    const SizeType number_of_shape_function_derivatives =
        static_cast<SizeType>(shape_function_derivatives_order);

    const SizeType dimension = rBackground.LocalSpaceDimension();
    IntegrationInfo integration_info = rBackground.GetDefaultIntegrationInfo();
    GeometriesArrayType quadrature_geometries;

    // Fill Gauss points per span and direction
    if (rParameters.Has("number_of_integration_points_per_span")) {
        for (IndexType dir = 0; dir < dimension; ++dir) {
            const int num = rParameters["number_of_integration_points_per_span"].GetInt();
            integration_info.SetNumberOfIntegrationPointsPerSpan(dir, num);
        }
    }

    // 2D case: facets on the interior knot lines
    if (dimension == 2) {
        const auto p_surface = dynamic_cast<NurbsSurfaceType*>(&rBackground);

        KRATOS_ERROR_IF(p_surface == nullptr)
            << "IgaModelerGsbm: the background patch is not a NURBS surface." << std::endl;
        KRATOS_ERROR_IF(p_surface->IsRational())
            << "IgaModelerGsbm: the background patch is rational. "
            << "Only B-spline surfaces are supported." << std::endl;

        const std::array<SizeType, 2> polynomial_degrees{
            p_surface->PolynomialDegreeU(), p_surface->PolynomialDegreeV()
        };
        const std::array<const Vector*, 2> knots{&p_surface->KnotsU(), &p_surface->KnotsV()};
        const SizeType max_polynomial_degree =
            std::max(polynomial_degrees[0], polynomial_degrees[1]);

        // Derivatives up to order p are needed by the ghost penalty and the
        // Taylor expansions of the facet terms
        KRATOS_ERROR_IF(number_of_shape_function_derivatives < max_polynomial_degree + 1)
            << "IgaModelerGsbm: \"shape_function_derivatives_order\" of the GSBM "
            << "facets must be at least p + 1 = " << max_polynomial_degree + 1
            << ", so that the p-th derivatives are stored." << std::endl;

        // Knot span boundaries, exactly as stored in the knot vectors
        std::array<std::vector<double>, 2> spans;
        for (IndexType dir = 0; dir < 2; ++dir) {
            p_surface->SpansLocalSpace(spans[dir], dir);

            // Kratos knot vectors have multiplicity p at the ends, so any
            // extra knot is a repeated interior one
            const SizeType expected_number_of_knots =
                spans[dir].size() + 2 * polynomial_degrees[dir] - 2;
            KRATOS_ERROR_IF(knots[dir]->size() != expected_number_of_knots)
                << "IgaModelerGsbm: the background patch has repeated interior "
                << "knots, which are not supported (the ghost penalty only "
                << "controls the jump of the p-th derivative)." << std::endl;
        }

        // Span indices, first knot index at p - 1 (open knot vector)
        const auto knot_span_index =
            [&polynomial_degrees](const IndexType Direction, const IndexType KnotIndex) {
                return KnotIndex - (polynomial_degrees[Direction] - 1);
            };

        NurbsSurfaceShapeFunction shape_function_container(
            polynomial_degrees[0], polynomial_degrees[1],
            number_of_shape_function_derivatives
        );

        // Facet tangents per direction, so that Normal() of both sides is
        // n-, pointing from - to + (+x or +y)
        const std::array<std::array<double, 2>, 2> tangents{
            {{0.0, -1.0}, {1.0, 0.0}}
        };

        for (IndexType dir = 0; dir < 2; ++dir) {
            // Facets on the knot lines normal to direction dir,
            // integrated along the other direction
            const IndexType other = 1 - dir;

            // Gauss points along the other direction
            IntegrationPointUtilities::IntegrationPointsArrayType line_points;
            IntegrationPointUtilities::CreateIntegrationPoints1DGauss(
                line_points, spans[other],
                integration_info.GetNumberOfIntegrationPointsPerSpan(other)
            );

            for (IndexType k = 1; k + 1 < spans[dir].size(); ++k) {
                // Loop over interior knot lines (skip patch boundary)
                const double knot = spans[dir][k];

                // Spans on each side of the knot line: lower (left/below)
                // and upper (right/above)
                const IndexType lower_span = NurbsUtilities::GetLowerSpan(
                    polynomial_degrees[dir], *knots[dir], knot
                );
                const IndexType upper_span = NurbsUtilities::GetUpperSpan(
                    polynomial_degrees[dir], *knots[dir], knot
                );

                for (const auto& r_line_point : line_points) {
                    // Loop through the 1D Gauss rules along the knot lines
                    const double coordinate = r_line_point[0];
                    const IndexType span_other = NurbsUtilities::GetLowerSpan(
                        polynomial_degrees[other], *knots[other], coordinate
                    );

                    // Parametric coordinates of the facet point
                    std::array<double, 2> uv;
                    uv[dir] = knot;
                    uv[other] = coordinate;

                    // Knot indices of the - span
                    std::array<IndexType, 2> spans_minus;
                    spans_minus[dir] = lower_span;
                    spans_minus[other] = span_other;

                    // Knot indices of the + span
                    std::array<IndexType, 2> spans_plus = spans_minus;
                    spans_plus[dir] = upper_span;

                    // Same point evaluated from each side
                    auto p_minus = CreateOneSidedQuadraturePoint(
                        *p_surface, shape_function_container,
                        number_of_shape_function_derivatives,
                        spans_minus, uv, r_line_point.Weight(), tangents[dir]
                    );
                    auto p_plus = CreateOneSidedQuadraturePoint(
                        *p_surface, shape_function_container,
                        number_of_shape_function_derivatives,
                        spans_plus, uv, r_line_point.Weight(), tangents[dir]
                    );

                    // Knot span indices and sizes of each side's span
                    std::vector<IndexType> knot_span_indices_minus(2);
                    std::vector<IndexType> knot_span_indices_plus(2);
                    Vector knot_span_sizes_minus(2);
                    Vector knot_span_sizes_plus(2);

                    for (IndexType d = 0; d < 2; ++d) {
                        knot_span_indices_minus[d] = knot_span_index(d, spans_minus[d]);
                        knot_span_indices_plus[d] = knot_span_index(d, spans_plus[d]);
                        knot_span_sizes_minus[d] = (*knots[d])[spans_minus[d] + 1]
                                                 - (*knots[d])[spans_minus[d]];
                        knot_span_sizes_plus[d] = (*knots[d])[spans_plus[d] + 1]
                                                 - (*knots[d])[spans_plus[d]];
                    }

                    p_minus->SetValue(KNOT_SPAN_INDICES, knot_span_indices_minus);
                    p_plus->SetValue(KNOT_SPAN_INDICES, knot_span_indices_plus);
                    p_minus->SetValue(KNOT_SPAN_SIZES, knot_span_sizes_minus);
                    p_plus->SetValue(KNOT_SPAN_SIZES, knot_span_sizes_plus);

                    quadrature_geometries.push_back(
                        Kratos::make_shared<CouplingGeometry<Node>>(p_minus, p_plus)
                    );
                }
            }
        }
    } else {
        // TODO: 3D case
        KRATOS_ERROR << "IgaModelerGsbm: 3D case not implemented yet." << std::endl;
    }

    CreateConditions(
        quadrature_geometries.ptr_begin(), quadrature_geometries.ptr_end(),
        rModelPart, name
    );

    if (!quadrature_geometries.empty()) {
        rModelPart.AddNodes(rBackground.Points().ptr_begin(), rBackground.Points().ptr_end());
    }

    KRATOS_INFO_IF("IgaModelerGsbm", mEchoLevel > 0)
        << "Created " << quadrature_geometries.size() << " facet conditions of type "
        << name << " in " << rModelPart.FullName() << "." << std::endl;
}

void IgaModelerGsbm::CreateConditionsOnGsbmSkeleton(
    GeometriesArrayType& rFacets,
    ModelPart& rModelPart,
    const Parameters rParameters) const
{
    KRATOS_ERROR_IF(rParameters["type"].GetString() != "condition")
        << "IgaModelerGsbm: entries on an existing GSBM skeleton create conditions, so "
        << "they require type \"condition\"." << std::endl;

    const std::string name = rParameters["name"].GetString();

    // Same geometries, so alpha and d set on them are shared by both entries
    CreateConditions(rFacets.ptr_begin(), rFacets.ptr_end(), rModelPart, name);

    // Control points of the skeleton, as in the entry that created it
    ModelPart& r_skeleton_model_part = rModelPart.GetParentModelPart().GetSubModelPart(
        rParameters["gsbm_parameters"]["skeleton_model_part"].GetString()
    );
    rModelPart.AddNodes(r_skeleton_model_part.NodesBegin(), r_skeleton_model_part.NodesEnd());

    KRATOS_INFO_IF("IgaModelerGsbm", mEchoLevel > 0)
        << "Created " << rFacets.size() << " conditions of type " << name
        << " on the skeleton of " << r_skeleton_model_part.FullName() << " in "
        << rModelPart.FullName() << "." << std::endl;
}

void IgaModelerGsbm::CreateElements(
    GeometriesArrayType::ptr_iterator GeometriesBegin,
    GeometriesArrayType::ptr_iterator GeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rElementName) const
{
    KRATOS_ERROR_IF_NOT(KratosComponents<Element>::Has(rElementName))
        << "IgaModelerGsbm: element \"" << rElementName << "\" is not registered." << std::endl;

    // Reference element: Id=0, placeholder geometry and no properties
    const Element& r_reference_element = KratosComponents<Element>::Get(rElementName);

    const ModelPart& r_root_model_part = rModelPart.GetRootModelPart();
    IndexType id = r_root_model_part.NumberOfElements() > 0
        // If there are registered elements, continue from the last one
        ? r_root_model_part.Elements().back().Id() + 1
        // Else, start from 1
        : 1;

    ModelPart::ElementsContainerType new_elements;
    new_elements.reserve(std::distance(GeometriesBegin, GeometriesEnd));

    for (auto it = GeometriesBegin; it != GeometriesEnd; ++it) {
        // One element per geometry, without properties (assigned later)
        const auto p_new_element = r_reference_element.Create(id++, *it, Properties::Pointer());
        new_elements.push_back(p_new_element);
    }

    rModelPart.AddElements(new_elements.begin(), new_elements.end());
}

void IgaModelerGsbm::CreateConditions(
    GeometriesArrayType::ptr_iterator GeometriesBegin,
    GeometriesArrayType::ptr_iterator GeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rConditionName) const
{
    KRATOS_ERROR_IF_NOT(KratosComponents<Condition>::Has(rConditionName))
        << "IgaModelerGsbm: condition \"" << rConditionName << "\" is not registered." << std::endl;

    // Reference condition: Id=0, placeholder geometry and no properties
    const Condition& r_reference_condition = KratosComponents<Condition>::Get(rConditionName);

    const ModelPart& r_root_model_part = rModelPart.GetRootModelPart();
    IndexType id = r_root_model_part.NumberOfConditions() > 0
        // If there are registered conditions, continue from the last one
        ? r_root_model_part.Conditions().back().Id() + 1
        // Else, start from 1
        : 1;

    ModelPart::ConditionsContainerType new_conditions;
    new_conditions.reserve(std::distance(GeometriesBegin, GeometriesEnd));

    for (auto it = GeometriesBegin; it != GeometriesEnd; ++it) {
        // One condition per geometry, without properties (assigned later)
        const auto p_new_condition = r_reference_condition.Create(id++, *it, Properties::Pointer());
        new_conditions.push_back(p_new_condition);
    }

    rModelPart.AddConditions(new_conditions.begin(), new_conditions.end());
}

const Parameters IgaModelerGsbm::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level"               : 0,
        "analysis_model_part_name" : "",
        "geometry_name"            : "",
        "element_condition_list"   : []
    })");
}

} // namespace Kratos
