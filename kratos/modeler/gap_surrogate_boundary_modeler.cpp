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
//

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>

// Project includes
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "modeler/gap_surrogate_boundary_modeler.h"
#include "utilities/gap_quadrature_utility.h"

namespace Kratos {

namespace {

struct GapFaceData
{
    Condition* pBoundaryCondition = nullptr;
    Element* pParentElement = nullptr;
    std::array<array_1d<double, 3>, 2> SurrogatePoints;
    std::array<array_1d<double, 3>, 2> TruePoints;
    std::array<std::size_t, 2> SurrogateNodeIds;
    array_1d<double, 3> PatchCenter = ZeroVector(3);
    bool HasPositiveMeasure = true;
};

struct IncidentFace
{
    std::size_t FaceIndex = 0;
    std::size_t EndpointIndex = 0;
};

Parameters MakeBaseModelerParameters(const Parameters& rParameters)
{
    Parameters base_parameters = rParameters.Clone();
    const std::array<std::string, 5> gap_only_entries{{
        "gap_element_sub_model_part_name",
        "gap_interface_sub_model_part_name",
        "gap_element_name",
        "gap_interface_condition_name",
        "gap_integration_order"}};
    for (const auto& r_entry : gap_only_entries) {
        base_parameters.RemoveValue(r_entry);
    }

    // Endpoint projections are the four corners of the ruled gap patch. This
    // internal setting is forced by the derived modeler, irrespective of a
    // value supplied by the caller.
    if (base_parameters.Has("project_surrogate_face_end_points")) {
        base_parameters["project_surrogate_face_end_points"].SetBool(true);
    } else {
        base_parameters.AddBool("project_surrogate_face_end_points", true);
    }
    return base_parameters;
}

std::size_t FindMaximumElementId(const ModelPart& rModelPart)
{
    std::size_t maximum_id = 0;
    for (const auto& r_element : rModelPart.Elements()) {
        maximum_id = std::max(maximum_id, r_element.Id());
    }
    return maximum_id;
}

std::size_t FindMaximumConditionId(const ModelPart& rModelPart)
{
    std::size_t maximum_id = 0;
    for (const auto& r_condition : rModelPart.Conditions()) {
        maximum_id = std::max(maximum_id, r_condition.Id());
    }
    return maximum_id;
}

array_1d<double, 3> MatrixRowToPoint(
    const Matrix& rMatrix,
    const std::size_t Row)
{
    return array_1d<double, 3>{
        rMatrix(Row, 0), rMatrix(Row, 1), rMatrix(Row, 2)};
}

double Distance(const array_1d<double, 3>& rA, const array_1d<double, 3>& rB)
{
    return norm_2(rA - rB);
}

std::size_t FindParentNodeIdAtPoint(
    const Geometry<Node>& rParentGeometry,
    const array_1d<double, 3>& rPoint,
    const double Tolerance)
{
    for (const auto& r_node : rParentGeometry) {
        if (Distance(r_node.Coordinates(), rPoint) <= Tolerance) {
            return r_node.Id();
        }
    }
    KRATOS_ERROR << "GapSurrogateBoundaryModeler could not identify the grid "
                 << "node at surrogate endpoint " << rPoint << " in parent "
                 << "element #" << rParentGeometry.Id() << "." << std::endl;
}

Point InterpolatePatchPoint(
    const GapFaceData& rFace,
    const double Xi,
    const double Eta)
{
    const double first_xi = 0.5 * (1.0 - Xi);
    const double second_xi = 0.5 * (1.0 + Xi);
    const array_1d<double, 3> surrogate =
        first_xi * rFace.SurrogatePoints[0] +
        second_xi * rFace.SurrogatePoints[1];
    const array_1d<double, 3> true_point =
        first_xi * rFace.TruePoints[0] +
        second_xi * rFace.TruePoints[1];
    const array_1d<double, 3> global_coordinates =
        0.5 * (1.0 - Eta) * surrogate +
        0.5 * (1.0 + Eta) * true_point;
    return Point(global_coordinates[0], global_coordinates[1], global_coordinates[2]);
}

double CalculateSignedPatchJacobian(
    const GapFaceData& rFace,
    const double Xi,
    const double Eta)
{
    const array_1d<double, 3> surrogate_tangent =
        0.5 * (rFace.SurrogatePoints[1] - rFace.SurrogatePoints[0]);
    const array_1d<double, 3> true_tangent =
        0.5 * (rFace.TruePoints[1] - rFace.TruePoints[0]);
    const array_1d<double, 3> derivative_xi =
        0.5 * (1.0 - Eta) * surrogate_tangent +
        0.5 * (1.0 + Eta) * true_tangent;

    const double first_xi = 0.5 * (1.0 - Xi);
    const double second_xi = 0.5 * (1.0 + Xi);
    const array_1d<double, 3> surrogate =
        first_xi * rFace.SurrogatePoints[0] +
        second_xi * rFace.SurrogatePoints[1];
    const array_1d<double, 3> true_point =
        first_xi * rFace.TruePoints[0] +
        second_xi * rFace.TruePoints[1];
    const array_1d<double, 3> derivative_eta = 0.5 * (true_point - surrogate);

    return derivative_xi[0] * derivative_eta[1] -
           derivative_xi[1] * derivative_eta[0];
}

bool ValidateGapPatchAndCheckHasPositiveMeasure(
    const GapFaceData& rFace,
    const double CoordinateTolerance)
{
    const double maximum_side_length = std::max({
        Distance(rFace.SurrogatePoints[0], rFace.SurrogatePoints[1]),
        Distance(rFace.TruePoints[0], rFace.TruePoints[1]),
        CoordinateTolerance});
    const double jacobian_tolerance =
        CoordinateTolerance * maximum_side_length;
    const auto& r_line_integration_points =
        LineGaussLegendreIntegrationPoints2::IntegrationPoints();

    bool has_nonzero_measure_point = false;
    bool has_zero_measure_point = false;
    double first_nonzero_signed_jacobian = 0.0;
    for (const auto& r_eta_point : r_line_integration_points) {
        for (const auto& r_xi_point : r_line_integration_points) {
            const double signed_jacobian = CalculateSignedPatchJacobian(
                rFace, r_xi_point.X(), r_eta_point.X());
            if (std::abs(signed_jacobian) <= jacobian_tolerance) {
                has_zero_measure_point = true;
                continue;
            }

            if (!has_nonzero_measure_point) {
                first_nonzero_signed_jacobian = signed_jacobian;
                has_nonzero_measure_point = true;
            } else {
                KRATOS_ERROR_IF(
                    first_nonzero_signed_jacobian * signed_jacobian < 0.0)
                    << "GapSurrogateBoundaryModeler found a folded ruled gap "
                    << "patch for surrogate condition #"
                    << rFace.pBoundaryCondition->Id()
                    << ". Check the endpoint projections." << std::endl;
            }
        }
    }

    // A zero-measure patch is a legitimate topological limit. For example,
    // when a smooth skin is tangent to the Cartesian grid, a one-cell-wide
    // surrogate spike can have both endpoint projections at the same tangent
    // point. It owns no gap volume and therefore must not create an element,
    // but the face remains in the topology: its parent can still be the
    // volume-side neighbour of an adjacent, non-zero gap patch.
    if (!has_nonzero_measure_point) {
        return false;
    }

    KRATOS_ERROR_IF(has_zero_measure_point)
        << "GapSurrogateBoundaryModeler found a partially degenerate ruled "
        << "gap patch for surrogate condition #"
        << rFace.pBoundaryCondition->Id()
        << ". Check the endpoint projections." << std::endl;
    return true;
}

Geometry<Node>::Pointer CreateGapQuadratureGeometry(
    const GapFaceData& rFace)
{
    Geometry<Node>::IntegrationPointsArrayType points;
    for (const auto& r_eta : LineGaussLegendreIntegrationPoints2::IntegrationPoints()) {
        for (const auto& r_xi : LineGaussLegendreIntegrationPoints2::IntegrationPoints()) {
            const Point point = InterpolatePatchPoint(rFace, r_xi.X(), r_eta.X());
            const double weight = std::abs(CalculateSignedPatchJacobian(
                rFace, r_xi.X(), r_eta.X())) * r_xi.Weight() * r_eta.Weight();
            points.emplace_back(point.X(), point.Y(), 0.0, weight);
        }
    }
    return GapQuadratureUtility::Create(rFace.pParentElement->GetGeometry(), points);
}

void SetExactBoundaryQuadrature(GapFaceData& rFace, const double Tolerance)
{
    const auto true_tangent = rFace.TruePoints[1] - rFace.TruePoints[0];
    const double true_length = std::hypot(true_tangent[0], true_tangent[1]);

    const auto surrogate_tangent =
        rFace.SurrogatePoints[1] - rFace.SurrogatePoints[0];
    const double surrogate_length = std::hypot(
        surrogate_tangent[0], surrogate_tangent[1]);
    KRATOS_ERROR_IF(surrogate_length <= Tolerance)
        << "GapSurrogateBoundaryModeler found a zero-length surrogate face."
        << std::endl;

    array_1d<double, 3> surrogate_normal{
        surrogate_tangent[1] / surrogate_length,
        -surrogate_tangent[0] / surrogate_length,
        0.0};
    array_1d<double, 3> true_normal = surrogate_normal;
    if (true_length > Tolerance) {
        true_normal[0] = true_tangent[1] / true_length;
        true_normal[1] = -true_tangent[0] / true_length;
        if (inner_prod(true_normal, surrogate_normal) < 0.0) {
            true_normal *= -1.0;
        }
    }

    const auto& r_integration_points =
        LineGaussLegendreIntegrationPoints2::IntegrationPoints();
    Matrix boundary_quadrature(r_integration_points.size(), 7);
    for (std::size_t g = 0; g < r_integration_points.size(); ++g) {
        const double local_coordinate = r_integration_points[g].X();
        const double first_weight = 0.5 * (1.0 - local_coordinate);
        const double second_weight = 0.5 * (1.0 + local_coordinate);
        const auto global_point =
            first_weight * rFace.TruePoints[0] +
            second_weight * rFace.TruePoints[1];
        for (std::size_t d = 0; d < 3; ++d) {
            boundary_quadrature(g, d) = global_point[d];
            boundary_quadrature(g, 3 + d) = true_normal[d];
        }
        boundary_quadrature(g, 6) =
            0.5 * true_length * r_integration_points[g].Weight();
    }

    // Both surrogate endpoints can project to the same true point at a skin
    // vertex or a grid tangency. The ruled patch is either triangular with
    // positive area or entirely collapsed (handled by the topology code), but
    // its physical-boundary side has zero measure in both cases. Encode a
    // zero-weight boundary rule so the CFD condition contributes nothing on
    // this collapsed side.
    rFace.pBoundaryCondition->SetValue(
        SURROGATE_BOUNDARY_PROJECTION, boundary_quadrature);
}

Condition::NodesArrayType CollectUniqueParentNodes(
    const Element& rFirstParent,
    const Element& rSecondParent)
{
    Condition::NodesArrayType nodes;
    std::unordered_map<std::size_t, bool> inserted_node_ids;
    for (const Element* p_parent : {&rFirstParent, &rSecondParent}) {
        for (std::size_t n = 0; n < p_parent->GetGeometry().size(); ++n) {
            const auto p_node = p_parent->GetGeometry().pGetPoint(n);
            if (inserted_node_ids.emplace(p_node->Id(), true).second) {
                nodes.push_back(p_node);
            }
        }
    }
    return nodes;
}

bool ParentExtensionsCoincideOnConnector(
    const Element& rFirstParent,
    const Element& rSecondParent,
    const array_1d<double, 3>& rConnectorStart,
    const array_1d<double, 3>& rConnectorEnd,
    const double Tolerance)
{
    std::vector<array_1d<double, 3>> shared_points;
    for (const auto& r_first_node : rFirstParent.GetGeometry()) {
        for (const auto& r_second_node : rSecondParent.GetGeometry()) {
            if (r_first_node.Id() == r_second_node.Id()) {
                shared_points.push_back(r_first_node.Coordinates());
            }
        }
    }
    if (shared_points.size() < 2) {
        return false;
    }

    // Two adjacent Q1 elements have the same trace (and the same linear
    // extrapolation) on the straight line containing their common edge. A
    // connector on that line therefore needs no artificial coupling entity.
    const auto edge = shared_points[1] - shared_points[0];
    const auto connector = rConnectorEnd - rConnectorStart;
    const auto offset = rConnectorStart - shared_points[0];
    const double edge_norm = std::hypot(edge[0], edge[1]);
    const double connector_norm = std::hypot(connector[0], connector[1]);
    const double direction_cross = std::abs(
        edge[0] * connector[1] - edge[1] * connector[0]);
    const double offset_cross = std::abs(
        edge[0] * offset[1] - edge[1] * offset[0]);
    return direction_cross <= Tolerance * edge_norm * connector_norm &&
           offset_cross <= Tolerance * edge_norm;
}

} // unnamed namespace

GapSurrogateBoundaryModeler::GapSurrogateBoundaryModeler(
    Model& rModel,
    Parameters ModelerParameters)
    : SurrogateBoundaryModeler(
          rModel, MakeBaseModelerParameters(ModelerParameters)),
      mpGapModel(&rModel),
      mGapParameters(ModelerParameters)
{
    mGapParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

const Parameters GapSurrogateBoundaryModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "echo_level" : 0,
        "domain_size" : 2,
        "input_model_part_name" : "SkinModelPart",
        "output_model_part_name" : "SurrogateModelPart",
        "volume_sub_model_part_name" : "FluidDomain",
        "surrogate_sub_model_part_name" : "SurrogateBoundary",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [1.0, 1.0, 0.0],
        "number_of_elements" : [10, 10, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "classification_samples_per_direction" : 3,
        "element_name" : "Element2D4N",
        "condition_name" : "LineCondition2D2N",
        "properties_id" : 1,
        "number_of_boundary_gauss_points" : 2,
        "project_surrogate_face_end_points" : true,
        "is_moving" : false,
        "background_sub_model_part_name" : "BackgroundGrid",
        "gap_element_sub_model_part_name" : "GapElements",
        "gap_interface_sub_model_part_name" : "GapInterfaces",
        "gap_element_name" : "",
        "gap_interface_condition_name" : "",
        "gap_integration_order" : 2
    })");
}

void GapSurrogateBoundaryModeler::SetupModelPart()
{
    KRATOS_ERROR_IF(mpGapModel == nullptr)
        << "GapSurrogateBoundaryModeler has no Model container." << std::endl;
    KRATOS_ERROR_IF(mGapParameters["domain_size"].GetInt() != 2)
        << "GapSurrogateBoundaryModeler currently supports domain_size=2 only."
        << std::endl;
    KRATOS_ERROR_IF(mGapParameters["gap_integration_order"].GetInt() != 2)
        << "GapSurrogateBoundaryModeler currently supports a 2x2 gap "
        << "quadrature only." << std::endl;
    KRATOS_ERROR_IF(mGapParameters["gap_element_name"].GetString().empty())
        << "GapSurrogateBoundaryModeler requires 'gap_element_name'."
        << std::endl;

    SurrogateBoundaryModeler::SetupModelPart();
    SetupGapModelPart2D();
}

void GapSurrogateBoundaryModeler::PrepareGeometryModel()
{
    // The standard Kratos lifecycle calls this stage once before
    // SetupModelPart. Match the base-modeler no-op on that initial call.
    if (!IsMovingModelInitialized()) {
        return;
    }
    KRATOS_ERROR_IF_NOT(mGapParameters["is_moving"].GetBool())
        << "GapSurrogateBoundaryModeler::PrepareGeometryModel requires "
        << "'is_moving' to be true." << std::endl;

    // Gap entities hold pointers to the old surrogate conditions and active
    // parent elements. Destroy this transient layer before the base modeler
    // replaces the active view and its boundary.
    ClearGapModelPart2D();
    SurrogateBoundaryModeler::PrepareGeometryModel();
    SetupGapModelPart2D();
}

void GapSurrogateBoundaryModeler::ClearGapModelPart2D()
{
    const std::string output_name =
        mGapParameters["output_model_part_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpGapModel->HasModelPart(output_name))
        << "GapSurrogateBoundaryModeler::SetupModelPart must be called before "
        << "PrepareGeometryModel." << std::endl;
    ModelPart& r_output = mpGapModel->GetModelPart(output_name);

    const std::string gap_elements_name =
        mGapParameters["gap_element_sub_model_part_name"].GetString();
    if (r_output.HasSubModelPart(gap_elements_name)) {
        ModelPart& r_gap_elements = r_output.GetSubModelPart(gap_elements_name);
        std::vector<std::size_t> element_ids;
        element_ids.reserve(r_gap_elements.NumberOfElements());
        for (const auto& r_element : r_gap_elements.Elements()) {
            element_ids.push_back(r_element.Id());
        }
        for (const std::size_t element_id : element_ids) {
            r_output.RemoveElementFromAllLevels(element_id);
        }
        r_gap_elements.Nodes().clear();
    }

    const std::string gap_interfaces_name =
        mGapParameters["gap_interface_sub_model_part_name"].GetString();
    if (r_output.HasSubModelPart(gap_interfaces_name)) {
        ModelPart& r_gap_interfaces =
            r_output.GetSubModelPart(gap_interfaces_name);
        std::vector<std::size_t> condition_ids;
        condition_ids.reserve(r_gap_interfaces.NumberOfConditions());
        for (const auto& r_condition : r_gap_interfaces.Conditions()) {
            condition_ids.push_back(r_condition.Id());
        }
        for (const std::size_t condition_id : condition_ids) {
            r_output.RemoveConditionFromAllLevels(condition_id);
        }
        r_gap_interfaces.Nodes().clear();
    }
}

void GapSurrogateBoundaryModeler::SetupGapModelPart2D()
{
    const std::string output_name =
        mGapParameters["output_model_part_name"].GetString();
    ModelPart& r_output = mpGapModel->GetModelPart(output_name);
    ModelPart& r_surrogate = r_output.GetSubModelPart(
        mGapParameters["surrogate_sub_model_part_name"].GetString());
    ModelPart& r_gap_elements = r_output.HasSubModelPart(
        mGapParameters["gap_element_sub_model_part_name"].GetString())
        ? r_output.GetSubModelPart(
              mGapParameters["gap_element_sub_model_part_name"].GetString())
        : r_output.CreateSubModelPart(
              mGapParameters["gap_element_sub_model_part_name"].GetString());
    ModelPart& r_gap_interfaces = r_output.HasSubModelPart(
        mGapParameters["gap_interface_sub_model_part_name"].GetString())
        ? r_output.GetSubModelPart(
              mGapParameters["gap_interface_sub_model_part_name"].GetString())
        : r_output.CreateSubModelPart(
              mGapParameters["gap_interface_sub_model_part_name"].GetString());

    const std::string gap_element_name =
        mGapParameters["gap_element_name"].GetString();
    KRATOS_ERROR_IF_NOT(KratosComponents<Element>::Has(gap_element_name))
        << "Gap element '" << gap_element_name << "' is not registered."
        << std::endl;
    const Element& r_gap_element_prototype =
        KratosComponents<Element>::Get(gap_element_name);

    const std::string interface_condition_name =
        mGapParameters["gap_interface_condition_name"].GetString();
    const bool create_interfaces = !interface_condition_name.empty();
    KRATOS_ERROR_IF(create_interfaces &&
                    !KratosComponents<Condition>::Has(interface_condition_name))
        << "Gap-interface condition '" << interface_condition_name
        << "' is not registered." << std::endl;

    const auto lower = mGapParameters["lower_point"].GetVector();
    const auto upper = mGapParameters["upper_point"].GetVector();
    const double domain_scale = std::max({
        upper[0] - lower[0], upper[1] - lower[1], 1.0});
    const double tolerance = 1.0e-12 * domain_scale;

    std::vector<GapFaceData> gap_faces;
    gap_faces.reserve(r_surrogate.NumberOfConditions());
    std::unordered_map<std::size_t, std::vector<IncidentFace>> incident_faces;

    for (auto& r_condition : r_surrogate.Conditions()) {
        auto& r_neighbours = r_condition.GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(r_neighbours.size() != 1)
            << "Surrogate condition #" << r_condition.Id()
            << " must have exactly one parent element." << std::endl;
        const Matrix& r_face_coordinates =
            r_condition.GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
        const Matrix& r_projections =
            r_condition.GetValue(SURROGATE_BOUNDARY_PROJECTION);
        KRATOS_ERROR_IF(r_face_coordinates.size1() != 2 ||
                        r_face_coordinates.size2() != 3)
            << "Surrogate condition #" << r_condition.Id()
            << " has invalid face coordinates." << std::endl;
        KRATOS_ERROR_IF(r_projections.size1() != 4 ||
                        r_projections.size2() != 3)
            << "Surrogate condition #" << r_condition.Id()
            << " requires two Gauss-point and two endpoint projections."
            << std::endl;

        GapFaceData face;
        face.pBoundaryCondition = &r_condition;
        face.pParentElement = &r_neighbours[0];
        face.SurrogatePoints[0] = MatrixRowToPoint(r_face_coordinates, 0);
        face.SurrogatePoints[1] = MatrixRowToPoint(r_face_coordinates, 1);
        face.TruePoints[0] = MatrixRowToPoint(r_projections, 2);
        face.TruePoints[1] = MatrixRowToPoint(r_projections, 3);
        for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
            face.SurrogateNodeIds[endpoint] = FindParentNodeIdAtPoint(
                face.pParentElement->GetGeometry(),
                face.SurrogatePoints[endpoint],
                tolerance);
        }
        face.PatchCenter = 0.25 * (
            face.SurrogatePoints[0] + face.SurrogatePoints[1] +
            face.TruePoints[0] + face.TruePoints[1]);

        const auto surrogate_tangent =
            face.SurrogatePoints[1] - face.SurrogatePoints[0];
        const double surrogate_length = std::hypot(
            surrogate_tangent[0], surrogate_tangent[1]);
        const array_1d<double, 3> outward_normal{
            surrogate_tangent[1] / surrogate_length,
            -surrogate_tangent[0] / surrogate_length,
            0.0};
        for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
            const auto gap_vector =
                face.TruePoints[endpoint] - face.SurrogatePoints[endpoint];
            KRATOS_ERROR_IF(inner_prod(outward_normal, gap_vector) < -tolerance)
                << "GapSurrogateBoundaryModeler currently supports additive "
                << "gaps only. Endpoint " << endpoint << " of condition #"
                << r_condition.Id() << " projects behind the oriented "
                << "surrogate face." << std::endl;
        }

        face.HasPositiveMeasure =
            ValidateGapPatchAndCheckHasPositiveMeasure(face, tolerance);
        SetExactBoundaryQuadrature(face, tolerance);
        const std::size_t face_index = gap_faces.size();
        gap_faces.push_back(face);
        for (std::size_t endpoint = 0; endpoint < 2; ++endpoint) {
            incident_faces[face.SurrogateNodeIds[endpoint]].push_back(
                IncidentFace{face_index, endpoint});
        }
    }

    // IDs belong to the root container and must be selected globally, also
    // when the requested output ModelPart is nested.
    ModelPart& r_root_output = r_output.GetRootModelPart();
    std::size_t next_element_id = FindMaximumElementId(r_root_output) + 1;
    std::size_t number_of_zero_measure_patches = 0;
    for (auto& r_face : gap_faces) {
        if (!r_face.HasPositiveMeasure) {
            ++number_of_zero_measure_patches;
            continue;
        }
        auto p_gap_geometry = CreateGapQuadratureGeometry(r_face);
        auto p_gap_element = r_gap_element_prototype.Create(
            next_element_id++,
            p_gap_geometry,
            r_face.pBoundaryCondition->pGetProperties());
        p_gap_element->SetValue(
            NEIGHBOUR_ELEMENTS,
            GlobalPointersVector<Element>({r_face.pParentElement}));
        // Keep the positive-area patch-to-face association explicit. This is
        // useful to update only the affected patch in a future moving-body
        // workflow without searching the full gap container. Zero-measure
        // faces intentionally have no corresponding element.
        p_gap_element->SetValue(
            NEIGHBOUR_CONDITIONS,
            GlobalPointersVector<Condition>({r_face.pBoundaryCondition}));
        r_output.AddElement(p_gap_element);
        r_gap_elements.AddElement(p_gap_element);
        for (std::size_t n = 0; n < p_gap_geometry->size(); ++n) {
            r_gap_elements.AddNode(p_gap_geometry->pGetPoint(n));
        }
    }

    std::size_t number_of_interface_conditions = 0;
    if (create_interfaces) {
        const Condition& r_interface_prototype =
            KratosComponents<Condition>::Get(interface_condition_name);
        std::size_t next_condition_id =
            FindMaximumConditionId(r_root_output) + 1;

        for (const auto& r_entry : incident_faces) {
            const auto& r_incident = r_entry.second;
            KRATOS_ERROR_IF(r_incident.size() != 2)
                << "GapSurrogateBoundaryModeler expects a closed manifold "
                << "surrogate boundary. Grid node #" << r_entry.first
                << " belongs to " << r_incident.size()
                << " surrogate faces instead of two." << std::endl;

            const GapFaceData& r_first = gap_faces[r_incident[0].FaceIndex];
            const GapFaceData& r_second = gap_faces[r_incident[1].FaceIndex];
            if (r_first.pParentElement->Id() == r_second.pParentElement->Id()) {
                // Both patches use the same polynomial extension, hence their
                // trace is already continuous and no coupling term is needed.
                continue;
            }

            const auto& r_first_true =
                r_first.TruePoints[r_incident[0].EndpointIndex];
            const auto& r_second_true =
                r_second.TruePoints[r_incident[1].EndpointIndex];
            KRATOS_ERROR_IF(Distance(r_first_true, r_second_true) > tolerance)
                << "Adjacent gap patches at surrogate node #" << r_entry.first
                << " do not share the same projected true endpoint." << std::endl;

            const auto& r_surrogate_point =
                r_first.SurrogatePoints[r_incident[0].EndpointIndex];
            // Materialize the uBLAS expression. With auto, its reference to
            // the temporary scalar 0.5 dangles after this statement.
            const array_1d<double, 3> true_point =
                0.5 * (r_first_true + r_second_true);
            if (Distance(r_surrogate_point, true_point) <= tolerance) {
                continue;
            }
            if (ParentExtensionsCoincideOnConnector(
                    *r_first.pParentElement,
                    *r_second.pParentElement,
                    r_surrogate_point,
                    true_point,
                    tolerance)) {
                continue;
            }

            Matrix connector_coordinates(2, 3);
            for (std::size_t d = 0; d < 3; ++d) {
                connector_coordinates(0, d) = r_surrogate_point[d];
                connector_coordinates(1, d) = true_point[d];
            }
            const array_1d<double, 3> connector =
                true_point - r_surrogate_point;
            array_1d<double, 3> interface_normal{
                connector[1], -connector[0], 0.0};
            if (inner_prod(
                    interface_normal,
                    r_second.PatchCenter - r_first.PatchCenter) < 0.0) {
                // Reversing the connector reverses the normal computed by the
                // CFD interface condition. The first neighbour is therefore
                // always the '+' side and the normal points towards side '-'.
                for (std::size_t d = 0; d < 3; ++d) {
                    std::swap(
                        connector_coordinates(0, d),
                        connector_coordinates(1, d));
                }
            }

            auto condition_nodes = CollectUniqueParentNodes(
                *r_first.pParentElement, *r_second.pParentElement);
            auto p_interface_geometry =
                Kratos::make_shared<Geometry<Node>>(condition_nodes);
            auto p_interface_condition = r_interface_prototype.Create(
                next_condition_id++,
                p_interface_geometry,
                r_first.pBoundaryCondition->pGetProperties());
            p_interface_condition->SetValue(
                NEIGHBOUR_ELEMENTS,
                GlobalPointersVector<Element>({
                    r_first.pParentElement, r_second.pParentElement}));
            p_interface_condition->SetValue(
                SURROGATE_BOUNDARY_FACE_COORDINATES,
                connector_coordinates);
            r_output.AddCondition(p_interface_condition);
            r_gap_interfaces.AddCondition(p_interface_condition);
            for (std::size_t n = 0; n < condition_nodes.size(); ++n) {
                r_gap_interfaces.AddNode(condition_nodes(n));
            }
            ++number_of_interface_conditions;
        }
    }

    KRATOS_INFO_IF(
        "GapSurrogateBoundaryModeler",
        mGapParameters["echo_level"].GetInt() > 0)
        << "Created " << r_gap_elements.NumberOfElements()
        << " ruled gap elements, omitted "
        << number_of_zero_measure_patches
        << " zero-measure patches, and created "
        << number_of_interface_conditions
        << " internal gap-interface conditions." << std::endl;
}

} // namespace Kratos
