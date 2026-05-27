//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

// Project includes
#include "nurbs_geometry_modeler_gap_sbm.h"
#include "custom_utilities/create_breps_sbm_utilities.h"
#include "custom_processes/snake_gap_sbm_process.h"
#include "iga_application_variables.h"

// System includes
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Kratos
{

namespace
{

constexpr double skin_point_tol = 1.0e-10;

enum class SegmentSelection
{
    All,
    CouplingOnly,
    ExcludeCoupling
};

using SkinCoordinatesType = array_1d<double, 3>;
using BrepCurveType = BrepCurve<PointerVector<Node>, PointerVector<Point>>;
using BrepCurveOnSurfaceSbmType = BrepCurveOnSurface<PointerVector<Node>, true, PointerVector<Point>>;
using BrepCurveOnSurfaceType = BrepCurveOnSurface<PointerVector<Node>, false, PointerVector<Point>>;
using BrepSurfaceType = BrepSurface<PointerVector<Node>, true, PointerVector<Point>>;
using BrepSurfaceTypeUnshifted = BrepSurface<PointerVector<Node>, false, PointerVector<Point>>;
using NurbsCurveType = NurbsCurveGeometry<3, PointerVector<Node>>;
using NurbsSurfaceType = NurbsGeometryModeler::NurbsSurfaceGeometryType;
using TrimmingCurveType = NurbsCurveGeometry<2, PointerVector<Point>>;

struct SkinSegment
{
    SkinCoordinatesType Start = ZeroVector(3);
    SkinCoordinatesType End = ZeroVector(3);
    std::string LayerName;
    std::string ConditionName;
    bool HasBrepId = false;
    int BrepId = 0;
    std::string BrepModelPartFullName;
};

struct PendingNodeMatch
{
    Node::Pointer pPendingNode;
    Node::Pointer pClosestCreatedNode;
    double ClosestDistance = std::numeric_limits<double>::max();
    bool HasExactMatch = false;
    bool IsResolved = false;
};

struct QuantizedCoord
{
    long long X = 0;
    long long Y = 0;
    long long Z = 0;

    bool operator==(const QuantizedCoord& rOther) const
    {
        return X == rOther.X && Y == rOther.Y && Z == rOther.Z;
    }
};

struct QuantizedCoordHasher
{
    std::size_t operator()(const QuantizedCoord& rValue) const noexcept
    {
        std::size_t seed = 0;
        seed ^= std::hash<long long>{}(rValue.X) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<long long>{}(rValue.Y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<long long>{}(rValue.Z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

bool IsCloseCoordinates(
    const SkinCoordinatesType& rA,
    const SkinCoordinatesType& rB,
    const double Tolerance = skin_point_tol)
{
    return norm_2(rA - rB) <= Tolerance * (1.0 + std::max({
        std::abs(rA[0]), std::abs(rA[1]), std::abs(rA[2]),
        std::abs(rB[0]), std::abs(rB[1]), std::abs(rB[2])}));
}

template<class TDataContainerType>
void SetBrepMetadata(
    TDataContainerType& rContainer,
    const int BrepId,
    const std::string& rBrepModelPartFullName)
{
    rContainer.SetValue(BREP_ID, BrepId);
    rContainer.SetValue(BREP_MODEL_PART_FULL_NAME, rBrepModelPartFullName);
}

ModelPart::IndexType NextLocalRefinementConnectorGeometryId(const ModelPart& rModelPart)
{
    constexpr ModelPart::IndexType local_refinement_connector_brep_id_base = 1000000000;

    ModelPart::IndexType next_id = local_refinement_connector_brep_id_base;
    for (const auto& r_geometry : rModelPart.Geometries()) {
        const auto id = r_geometry.Id();
        constexpr auto id_bits = sizeof(ModelPart::IndexType) * 8;
        const bool id_generated_from_string = id & (ModelPart::IndexType(1) << (id_bits - 1));
        const bool id_self_assigned = id & (ModelPart::IndexType(1) << (id_bits - 2));
        if (id_generated_from_string || id_self_assigned) {
            continue;
        }
        if (id >= local_refinement_connector_brep_id_base) {
            next_id = std::max(next_id, static_cast<ModelPart::IndexType>(id + 1));
        }
    }

    return next_id;
}

SkinSegment ReverseSegment(const SkinSegment& rSegment)
{
    SkinSegment reversed = rSegment;
    reversed.Start = rSegment.End;
    reversed.End = rSegment.Start;
    return reversed;
}

double ComputeSignedLoopArea(const std::vector<SkinSegment>& rSegments)
{
    double area = 0.0;
    for (const auto& r_segment : rSegments) {
        area += r_segment.Start[0] * r_segment.End[1] - r_segment.End[0] * r_segment.Start[1];
    }
    return 0.5 * area;
}

std::vector<SkinSegment> ExtractSegments(
    const ModelPart& rSkinSubModelPart,
    const SegmentSelection Selection)
{
    std::vector<SkinSegment> segments;
    segments.reserve(rSkinSubModelPart.NumberOfConditions());

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        KRATOS_ERROR_IF(r_condition.GetGeometry().PointsNumber() != 2)
            << "NurbsGeometryModelerGapSbm: expected line skin condition with 2 points, got "
            << r_condition.GetGeometry().PointsNumber() << " for condition #" << r_condition.Id() << std::endl;

        const std::string layer_name = r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
        const bool is_coupling_side = layer_name == "COUPLING_SIDE";
        if (Selection == SegmentSelection::CouplingOnly && !is_coupling_side) {
            continue;
        }
        if (Selection == SegmentSelection::ExcludeCoupling && is_coupling_side) {
            continue;
        }

        SkinSegment segment;
        segment.Start = r_condition.GetGeometry()[0].Coordinates();
        segment.End = r_condition.GetGeometry()[1].Coordinates();
        segment.LayerName = layer_name;
        if (r_condition.Has(CONDITION_NAME)) {
            segment.ConditionName = r_condition.GetValue(CONDITION_NAME);
        }
        if (r_condition.Has(BREP_ID)) {
            segment.HasBrepId = true;
            segment.BrepId = r_condition.GetValue(BREP_ID);
            if (r_condition.Has(BREP_MODEL_PART_FULL_NAME)) {
                segment.BrepModelPartFullName = r_condition.GetValue(BREP_MODEL_PART_FULL_NAME);
            }
        }
        segments.push_back(std::move(segment));
    }

    return segments;
}

bool HasCouplingSideSegments(const ModelPart& rSkinSubModelPart)
{
    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        if (r_condition.Has(LAYER_NAME) && r_condition.GetValue(LAYER_NAME) == "COUPLING_SIDE") {
            return true;
        }
    }

    return false;
}

bool HasCouplingSideSegmentsRecursive(const ModelPart& rSkinSubModelPart)
{
    if (HasCouplingSideSegments(rSkinSubModelPart)) {
        return true;
    }

    for (const auto& r_sub_model_part : rSkinSubModelPart.SubModelParts()) {
        if (HasCouplingSideSegmentsRecursive(r_sub_model_part)) {
            return true;
        }
    }

    return false;
}

void MarkCouplingSideConditionsForErase(
    ModelPart& rSkinSubModelPart,
    std::unordered_set<IndexType>& rCandidateNodeIds)
{
    for (auto& r_condition : rSkinSubModelPart.Conditions()) {
        const std::string layer_name = r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
        if (layer_name != "COUPLING_SIDE") {
            continue;
        }

        r_condition.Set(TO_ERASE, true);
        const auto& r_geometry = r_condition.GetGeometry();
        for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
            rCandidateNodeIds.insert(r_geometry[i].Id());
        }
    }

    for (auto& r_sub_model_part : rSkinSubModelPart.SubModelParts()) {
        MarkCouplingSideConditionsForErase(r_sub_model_part, rCandidateNodeIds);
    }
}

std::vector<SkinSegment> OrderSegmentsHeadToTail(const std::vector<SkinSegment>& rSegments)
{
    KRATOS_ERROR_IF(rSegments.empty())
        << "NurbsGeometryModelerGapSbm: cannot order an empty set of skin segments." << std::endl;

    std::vector<SkinSegment> ordered_segments;
    ordered_segments.reserve(rSegments.size());

    std::vector<bool> used(rSegments.size(), false);
    ordered_segments.push_back(rSegments.front());
    used.front() = true;

    SkinCoordinatesType current_end = ordered_segments.back().End;
    while (ordered_segments.size() < rSegments.size()) {
        bool found = false;
        for (std::size_t i = 0; i < rSegments.size(); ++i) {
            if (used[i]) {
                continue;
            }

            if (IsCloseCoordinates(rSegments[i].Start, current_end)) {
                ordered_segments.push_back(rSegments[i]);
                current_end = ordered_segments.back().End;
                used[i] = true;
                found = true;
                break;
            }

            if (IsCloseCoordinates(rSegments[i].End, current_end)) {
                ordered_segments.push_back(ReverseSegment(rSegments[i]));
                current_end = ordered_segments.back().End;
                used[i] = true;
                found = true;
                break;
            }
        }

        KRATOS_ERROR_IF_NOT(found)
            << "NurbsGeometryModelerGapSbm: failed to order skin segments into a closed loop." << std::endl;
    }

    KRATOS_ERROR_IF_NOT(IsCloseCoordinates(ordered_segments.back().End, ordered_segments.front().Start))
        << "NurbsGeometryModelerGapSbm: ordered skin segments do not form a closed loop." << std::endl;

    return ordered_segments;
}

void ResetSubModelPartContents(ModelPart& rModelPart)
{
    rModelPart.Nodes().clear();
    rModelPart.Elements().clear();
    rModelPart.Conditions().clear();
    rModelPart.Geometries().clear();
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        ResetSubModelPartContents(r_sub_model_part);
    }
}

ModelPart::IndexType NextNodeId(ModelPart& rModelPart)
{
    ModelPart::IndexType max_id = 0;
    for (const auto& r_node : rModelPart.GetRootModelPart().Nodes()) {
        max_id = std::max(max_id, r_node.Id());
    }
    return max_id + 1;
}

ModelPart::IndexType NextConditionId(ModelPart& rModelPart)
{
    ModelPart::IndexType max_id = 0;
    for (const auto& r_condition : rModelPart.GetRootModelPart().Conditions()) {
        max_id = std::max(max_id, r_condition.Id());
    }
    return max_id + 1;
}

NurbsSurfaceType::Pointer FindPatchSurfaceGeometry(ModelPart& rPatchModelPart)
{
    for (auto& r_geometry : rPatchModelPart.GetRootModelPart().Geometries()) {
        auto p_geometry = rPatchModelPart.GetRootModelPart().pGetGeometry(r_geometry.Id());
        if (auto p_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_geometry)) {
            return p_surface;
        }
        if (auto p_brep_curve_on_surface = std::dynamic_pointer_cast<BrepCurveOnSurfaceSbmType>(p_geometry)) {
            if (auto p_curve_on_surface = p_brep_curve_on_surface->pGetCurveOnSurface()) {
                auto p_host_surface = p_curve_on_surface->pGetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX);
                if (auto p_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_host_surface)) {
                    return p_surface;
                }
            }
        }
    }

    return nullptr;
}

NurbsSurfaceType::Pointer GetPatchSurfaceGeometryById(
    ModelPart& rPatchModelPart,
    const IndexType GeometryId)
{
    KRATOS_ERROR_IF_NOT(rPatchModelPart.HasGeometry(GeometryId))
        << "NurbsGeometryModelerGapSbm: patch model part '" << rPatchModelPart.FullName()
        << "' has no geometry with id " << GeometryId << "." << std::endl;

    auto p_geometry = rPatchModelPart.pGetGeometry(GeometryId);

    auto p_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_geometry);
    if (p_surface) {
        return p_surface;
    }

    if (auto p_brep_surface = std::dynamic_pointer_cast<BrepSurfaceType>(p_geometry)) {
        auto p_background_geometry =
            p_brep_surface->pGetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX);
        p_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_background_geometry);
    } else if (auto p_brep_surface = std::dynamic_pointer_cast<BrepSurfaceTypeUnshifted>(p_geometry)) {
        auto p_background_geometry =
            p_brep_surface->pGetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX);
        p_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_background_geometry);
    }

    KRATOS_ERROR_IF(!p_surface)
        << "NurbsGeometryModelerGapSbm: geometry #" << GeometryId << " in patch model part '"
        << rPatchModelPart.FullName()
        << "' is neither a NurbsSurfaceType nor a BrepSurface backed by NurbsSurfaceType." << std::endl;

    return p_surface;
}

TrimmingCurveType::Pointer CreateStraightTrimmingCurve(
    const Point::Pointer& pFirstPoint,
    const Point::Pointer& pSecondPoint,
    const Vector& rActiveRangeKnotVector)
{
    PointerVector<Point> control_points;
    control_points.push_back(pFirstPoint);
    control_points.push_back(pSecondPoint);

    const int polynomial_degree = 1;
    Vector knot_vector = ZeroVector(4);
    knot_vector[0] = rActiveRangeKnotVector[0];
    knot_vector[1] = rActiveRangeKnotVector[0];
    knot_vector[2] = rActiveRangeKnotVector[1];
    knot_vector[3] = rActiveRangeKnotVector[1];

    return Kratos::make_shared<TrimmingCurveType>(
        control_points,
        polynomial_degree,
        knot_vector);
}

int CreateStraightConnectorBrepId(
    ModelPart& rPatchModelPart,
    const NurbsSurfaceType::Pointer& pSurface,
    const SkinCoordinatesType& rStartPoint,
    const SkinCoordinatesType& rEndPoint)
{
    KRATOS_ERROR_IF_NOT(pSurface)
        << "NurbsGeometryModelerGapSbm: null host surface while creating a local-refinement connector BREP in '"
        << rPatchModelPart.FullName() << "'." << std::endl;

    NurbsSurfaceType::CoordinatesArrayType start_local_coordinates = ZeroVector(3);
    NurbsSurfaceType::CoordinatesArrayType end_local_coordinates = ZeroVector(3);
    const int start_projection_status =
        pSurface->ProjectionPointGlobalToLocalSpace(rStartPoint, start_local_coordinates, 1.0e-12);
    const int end_projection_status =
        pSurface->ProjectionPointGlobalToLocalSpace(rEndPoint, end_local_coordinates, 1.0e-12);

    KRATOS_ERROR_IF(start_projection_status == 0 || end_projection_status == 0)
        << "NurbsGeometryModelerGapSbm: failed to project local-refinement connector endpoints onto patch surface '"
        << rPatchModelPart.FullName() << "'." << std::endl;

    Vector active_range_knot_vector = ZeroVector(2);
    if (std::abs(start_local_coordinates[0] - end_local_coordinates[0]) <= skin_point_tol) {
        active_range_knot_vector[0] = start_local_coordinates[1];
        active_range_knot_vector[1] = end_local_coordinates[1];
    } else {
        active_range_knot_vector[0] = start_local_coordinates[0];
        active_range_knot_vector[1] = end_local_coordinates[0];
    }
    std::sort(active_range_knot_vector.begin(), active_range_knot_vector.end());

    auto p_first_point = Kratos::make_shared<Point>(
        start_local_coordinates[0],
        start_local_coordinates[1],
        0.0);
    auto p_second_point = Kratos::make_shared<Point>(
        end_local_coordinates[0],
        end_local_coordinates[1],
        0.0);
    auto p_trimming_curve = CreateStraightTrimmingCurve(
        p_first_point,
        p_second_point,
        active_range_knot_vector);

    const IndexType new_brep_id = NextLocalRefinementConnectorGeometryId(rPatchModelPart.GetRootModelPart());
    NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);
    const bool curve_direction = true;
    auto p_brep_curve_on_surface = Kratos::make_shared<BrepCurveOnSurfaceType>(
        pSurface,
        p_trimming_curve,
        brep_active_range,
        curve_direction);
    p_brep_curve_on_surface->SetId(new_brep_id);
    SetBrepMetadata(
        *p_brep_curve_on_surface,
        static_cast<int>(new_brep_id),
        rPatchModelPart.FullName());
    rPatchModelPart.AddGeometry(p_brep_curve_on_surface);

    return static_cast<int>(new_brep_id);
}

int CreateBoundaryConnectorBrepId(
    ModelPart& rPatchModelPart,
    const SkinCoordinatesType& rStartPoint,
    const SkinCoordinatesType& rEndPoint)
{
    auto p_surface = FindPatchSurfaceGeometry(rPatchModelPart);
    KRATOS_ERROR_IF_NOT(p_surface)
        << "NurbsGeometryModelerGapSbm: could not find a NurbsSurface geometry in patch model part '"
        << rPatchModelPart.FullName() << "' while creating a local-refinement connector BREP." << std::endl;

    return CreateStraightConnectorBrepId(
        rPatchModelPart,
        p_surface,
        rStartPoint,
        rEndPoint);
}

int CreateLocalRefinementClosureConnectorBrepId(
    ModelPart& rPatchModelPart,
    const SkinCoordinatesType& rStartPoint,
    const SkinCoordinatesType& rEndPoint)
{
    Vector active_range_knot_vector = ZeroVector(2);
    active_range_knot_vector[0] = 0.0;
    active_range_knot_vector[1] = 1.0;

    Node::Pointer p_first_node(new Node(1, rStartPoint));
    Node::Pointer p_second_node(new Node(2, rEndPoint));

    PointerVector<Node> control_points;
    control_points.push_back(p_first_node);
    control_points.push_back(p_second_node);

    const int polynomial_degree = 1;
    Vector knot_vector = ZeroVector(4);
    knot_vector[0] = active_range_knot_vector[0];
    knot_vector[1] = active_range_knot_vector[0];
    knot_vector[2] = active_range_knot_vector[1];
    knot_vector[3] = active_range_knot_vector[1];

    auto p_nurbs_curve = Kratos::make_shared<NurbsCurveType>(
        control_points,
        polynomial_degree,
        knot_vector);

    const IndexType new_brep_id = NextLocalRefinementConnectorGeometryId(rPatchModelPart.GetRootModelPart());

    auto p_brep_curve = Kratos::make_shared<BrepCurveType>(p_nurbs_curve);
    p_brep_curve->SetId(new_brep_id);
    SetBrepMetadata(
        *p_brep_curve,
        static_cast<int>(new_brep_id),
        rPatchModelPart.FullName());
    rPatchModelPart.AddGeometry(p_brep_curve);

    return static_cast<int>(new_brep_id);
}

QuantizedCoord MakeQuantizedCoord(const SkinCoordinatesType& rCoordinates)
{
    return QuantizedCoord{
        static_cast<long long>(std::llround(rCoordinates[0] / skin_point_tol)),
        static_cast<long long>(std::llround(rCoordinates[1] / skin_point_tol)),
        static_cast<long long>(std::llround(rCoordinates[2] / skin_point_tol))};
}

void RebuildSkinLoop(
    ModelPart& rSkinSubModelPart,
    const std::vector<SkinSegment>& rOrderedSegments,
    const bool IsInnerLoop)
{
    KRATOS_ERROR_IF(rOrderedSegments.empty())
        << "NurbsGeometryModelerGapSbm: cannot rebuild an empty skin loop." << std::endl;

    ResetSubModelPartContents(rSkinSubModelPart);
    if (!rSkinSubModelPart.RecursivelyHasProperties(0)) {
        rSkinSubModelPart.CreateNewProperties(0);
    }
    auto p_properties = rSkinSubModelPart.pGetProperties(0);

    ModelPart& r_interface_vertices = rSkinSubModelPart.HasSubModelPart("interface_vertices")
        ? rSkinSubModelPart.GetSubModelPart("interface_vertices")
        : rSkinSubModelPart.CreateSubModelPart("interface_vertices");
    ModelPart& r_loop_sub_model_part = rSkinSubModelPart.HasSubModelPart("0")
        ? rSkinSubModelPart.GetSubModelPart("0")
        : rSkinSubModelPart.CreateSubModelPart("0");

    ModelPart::IndexType next_node_id = NextNodeId(rSkinSubModelPart);
    ModelPart::IndexType next_condition_id = NextConditionId(rSkinSubModelPart);

    std::unordered_map<QuantizedCoord, Node::Pointer, QuantizedCoordHasher> node_lookup;
    std::vector<Node::Pointer> ordered_nodes;
    ordered_nodes.reserve(rOrderedSegments.size());

    auto get_or_create_node = [&](const SkinCoordinatesType& rCoordinates) {
        const auto key = MakeQuantizedCoord(rCoordinates);
        const auto it = node_lookup.find(key);
        if (it != node_lookup.end()) {
            return it->second;
        }

        auto p_node = rSkinSubModelPart.CreateNewNode(
            next_node_id++, rCoordinates[0], rCoordinates[1], rCoordinates[2]);
        r_interface_vertices.AddNode(p_node);
        r_loop_sub_model_part.AddNode(p_node);
        node_lookup.emplace(key, p_node);
        ordered_nodes.push_back(p_node);
        return p_node;
    };

    for (std::size_t i = 0; i < rOrderedSegments.size(); ++i) {
        const auto& r_segment = rOrderedSegments[i];
        auto p_node_1 = get_or_create_node(r_segment.Start);
        auto p_node_2 = (i + 1 == rOrderedSegments.size())
            ? get_or_create_node(rOrderedSegments.front().Start)
            : get_or_create_node(r_segment.End);

        ModelPart& r_layer_model_part = rSkinSubModelPart.HasSubModelPart(r_segment.LayerName)
            ? rSkinSubModelPart.GetSubModelPart(r_segment.LayerName)
            : rSkinSubModelPart.CreateSubModelPart(r_segment.LayerName);
        r_layer_model_part.AddNode(p_node_1);
        r_layer_model_part.AddNode(p_node_2);

        Condition::Pointer p_condition = rSkinSubModelPart.CreateNewCondition(
            "LineCondition2D2N",
            next_condition_id++,
            {{p_node_1->Id(), p_node_2->Id()}},
            p_properties);

        p_condition->SetValue(LAYER_NAME, "COUPLING_SIDE");
        if (!r_segment.ConditionName.empty()) {
            p_condition->SetValue(CONDITION_NAME, r_segment.ConditionName);
        }
        if (r_segment.HasBrepId) {
            SetBrepMetadata(*p_condition, r_segment.BrepId, r_segment.BrepModelPartFullName);
        }

        r_layer_model_part.AddCondition(p_condition);
        r_loop_sub_model_part.AddCondition(p_condition);
    }

    const std::size_t number_of_nodes = ordered_nodes.size();
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        auto p_previous = ordered_nodes[(i + number_of_nodes - 1) % number_of_nodes];
        auto p_current = ordered_nodes[i];
        auto p_next = ordered_nodes[(i + 1) % number_of_nodes];

        array_1d<double, 3> tangent = p_next->Coordinates() - p_previous->Coordinates();
        const double tangent_norm = norm_2(tangent);
        if (tangent_norm > 1.0e-16) {
            tangent /= tangent_norm;
        }

        Vector normal = ZeroVector(3);
        normal[0] = tangent[1];
        normal[1] = -tangent[0];

        p_current->SetValue(LOCAL_TANGENT, tangent);
        p_current->SetValue(NORMAL, normal);
        p_current->SetValue(CURVATURE, 0.0);
        p_current->SetValue(IDENTIFIER, IsInnerLoop ? "inner" : "outer");
    }

    for (auto& r_condition : rSkinSubModelPart.Conditions()) {
        const std::string layer_name = r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
        const std::string condition_name = r_condition.Has(CONDITION_NAME) ? r_condition.GetValue(CONDITION_NAME) : "";

        for (IndexType i = 0; i < r_condition.GetGeometry().PointsNumber(); ++i) {
            auto p_node = rSkinSubModelPart.pGetNode(r_condition.GetGeometry()[i].Id());
            auto connected_layers = p_node->GetValue(CONNECTED_LAYERS);
            auto connected_conditions = p_node->GetValue(CONNECTED_CONDITIONS);

            if (std::find(connected_layers.begin(), connected_layers.end(), layer_name) == connected_layers.end()) {
                connected_layers.push_back(layer_name);
                connected_conditions.push_back(condition_name);
            }

            p_node->SetValue(CONNECTED_LAYERS, connected_layers);
            p_node->SetValue(CONNECTED_CONDITIONS, connected_conditions);
        }
    }
}

void EnsurePropertiesAvailable(ModelPart& rModelPart)
{
    if (!rModelPart.RecursivelyHasProperties(0)) {
        rModelPart.CreateNewProperties(0);
    } else {
        rModelPart.pGetProperties(0);
    }
}

bool ParsePatchSkinName(
    const std::string& rSkinModelPartName,
    std::string& rPrefix,
    int& rPatchIndex)
{
    const std::string suffix = "_skin";
    if (rSkinModelPartName.size() <= suffix.size() ||
        rSkinModelPartName.substr(rSkinModelPartName.size() - suffix.size()) != suffix) {
        return false;
    }

    const std::string stem = rSkinModelPartName.substr(0, rSkinModelPartName.size() - suffix.size());
    std::size_t digit_pos = stem.size();
    while (digit_pos > 0 && std::isdigit(static_cast<unsigned char>(stem[digit_pos - 1])) != 0) {
        --digit_pos;
    }
    if (digit_pos == stem.size()) {
        return false;
    }

    rPrefix = stem.substr(0, digit_pos);
    rPatchIndex = std::stoi(stem.substr(digit_pos));
    return true;
}

std::string ExtractPatchSuffixFromSkinName(const std::string& rSkinModelPartName)
{
    const std::string suffix = "_skin";
    KRATOS_ERROR_IF(rSkinModelPartName.size() <= suffix.size() ||
        rSkinModelPartName.substr(rSkinModelPartName.size() - suffix.size()) != suffix)
        << "NurbsGeometryModelerGapSbm: expected skin model part name ending with '_skin', got '"
        << rSkinModelPartName << "'." << std::endl;

    const std::size_t suffix_pos = rSkinModelPartName.size() - suffix.size();
    const std::size_t patch_separator_pos = rSkinModelPartName.rfind('_', suffix_pos - 1);
    KRATOS_ERROR_IF(patch_separator_pos == std::string::npos)
        << "NurbsGeometryModelerGapSbm: could not extract patch suffix from skin model part name '"
        << rSkinModelPartName << "'." << std::endl;

    return rSkinModelPartName.substr(
        patch_separator_pos + 1,
        suffix_pos - patch_separator_pos - 1);
}

std::string FindPatchModelPartNameFromSkinName(
    Model& rModel,
    const std::string& rSkinModelPartName)
{
    const std::string patch_suffix = ExtractPatchSuffixFromSkinName(rSkinModelPartName);
    const std::string dotted_suffix = "." + patch_suffix;
    const std::string skin_suffix = "_" + patch_suffix + "_skin";

    std::string owning_model_part_suffix;
    const std::size_t skin_suffix_pos = rSkinModelPartName.rfind(skin_suffix);
    if (skin_suffix_pos != std::string::npos) {
        owning_model_part_suffix = rSkinModelPartName.substr(0, skin_suffix_pos);
        const std::string skin_prefix = "skin_";
        if (owning_model_part_suffix.substr(0, skin_prefix.size()) == skin_prefix) {
            owning_model_part_suffix = owning_model_part_suffix.substr(skin_prefix.size());
        }
    }

    std::string candidate_name;
    for (const auto& r_model_part_name : rModel.GetModelPartNames()) {
        const bool is_patch_model_part = r_model_part_name.size() > dotted_suffix.size() &&
            r_model_part_name.substr(r_model_part_name.size() - dotted_suffix.size()) == dotted_suffix;
        if (!is_patch_model_part) {
            continue;
        }

        if (!owning_model_part_suffix.empty()) {
            const std::string full_dotted_suffix = "." + owning_model_part_suffix + dotted_suffix;
            if (r_model_part_name.size() < full_dotted_suffix.size() ||
                r_model_part_name.substr(r_model_part_name.size() - full_dotted_suffix.size()) != full_dotted_suffix) {
                continue;
            }
        }

        KRATOS_ERROR_IF(!candidate_name.empty() && candidate_name != r_model_part_name)
            << "NurbsGeometryModelerGapSbm: found multiple candidate patch model parts for skin '"
            << rSkinModelPartName << "': '" << candidate_name << "' and '" << r_model_part_name << "'."
            << std::endl;
        candidate_name = r_model_part_name;
    }

    if (candidate_name.empty() && rModel.HasModelPart(patch_suffix)) {
        candidate_name = patch_suffix;
    }

    KRATOS_ERROR_IF(candidate_name.empty())
        << "NurbsGeometryModelerGapSbm: could not find the patch model part corresponding to skin '"
        << rSkinModelPartName << "'." << std::endl;

    return candidate_name;
}

std::vector<std::string> CollectRefinementSkinModelPartNames(
    Model& rModel,
    const std::string& rBaseSkinModelPartName)
{
    std::string prefix;
    int patch_index = 0;
    KRATOS_ERROR_IF_NOT(ParsePatchSkinName(rBaseSkinModelPartName, prefix, patch_index))
        << "NurbsGeometryModelerGapSbm: expected base skin model part name like 'Patch1_skin', got '"
        << rBaseSkinModelPartName << "'." << std::endl;

    std::vector<std::string> refinement_skin_names;
    for (int i = patch_index + 1; ; ++i) {
        const std::string candidate_name = prefix + std::to_string(i) + "_skin";
        if (!rModel.HasModelPart(candidate_name)) {
            break;
        }
        refinement_skin_names.push_back(candidate_name);
    }

    return refinement_skin_names;
}

std::vector<SkinSegment> ExtractCouplingSegmentsFromBrepCurves(
    const ModelPart& rSkinSubModelPart,
    ModelPart& rPatchModelPart)
{
    std::vector<SkinSegment> segments;
    std::unordered_set<int> visited_brep_ids;

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const std::string layer_name = r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
        if (layer_name != "COUPLING_SIDE") {
            continue;
        }

        SkinCoordinatesType global_start = ZeroVector(3);
        SkinCoordinatesType global_end = ZeroVector(3);

        SkinSegment segment;
        segment.LayerName = "COUPLING_SIDE";
        if (r_condition.Has(CONDITION_NAME)) {
            segment.ConditionName = r_condition.GetValue(CONDITION_NAME);
        }

        if (r_condition.Has(BREP_ID)) {
            const int brep_id = r_condition.GetValue(BREP_ID);
            if (!visited_brep_ids.insert(brep_id).second) {
                continue;
            }

            KRATOS_ERROR_IF_NOT(rPatchModelPart.HasGeometry(static_cast<ModelPart::IndexType>(brep_id)))
                << "NurbsGeometryModelerGapSbm: patch model part '" << rPatchModelPart.FullName()
                << "' has no BrepCurveOnSurface geometry with id " << brep_id << "." << std::endl;

            auto p_brep_geometry = rPatchModelPart.pGetGeometry(static_cast<ModelPart::IndexType>(brep_id));
            auto p_brep_curve_on_surface = std::dynamic_pointer_cast<BrepCurveOnSurfaceSbmType>(p_brep_geometry);
            KRATOS_ERROR_IF(!p_brep_curve_on_surface)
                << "NurbsGeometryModelerGapSbm: geometry #" << brep_id
                << " is not a BrepCurveOnSurfaceSbmType." << std::endl;

            NurbsInterval brep_domain_interval = p_brep_curve_on_surface->DomainInterval();
            SkinCoordinatesType local_start = ZeroVector(3);
            SkinCoordinatesType local_end = ZeroVector(3);
            local_start[0] = brep_domain_interval.GetT0();
            local_end[0] = brep_domain_interval.GetT1();

            p_brep_curve_on_surface->GlobalCoordinates(global_start, local_start);
            p_brep_curve_on_surface->GlobalCoordinates(global_end, local_end);

            segment.HasBrepId = true;
            segment.BrepId = brep_id;
            segment.BrepModelPartFullName = r_condition.Has(BREP_MODEL_PART_FULL_NAME)
                ? r_condition.GetValue(BREP_MODEL_PART_FULL_NAME)
                : rPatchModelPart.FullName();
        } else {
            KRATOS_ERROR_IF(r_condition.GetGeometry().PointsNumber() < 2)
                << "NurbsGeometryModelerGapSbm: COUPLING_SIDE condition #" << r_condition.Id()
                << " in '" << rSkinSubModelPart.FullName() << "' has fewer than 2 points." << std::endl;
            global_start = r_condition.GetGeometry()[0].Coordinates();
            global_end = r_condition.GetGeometry()[1].Coordinates();
            segment.HasBrepId = true;
            segment.BrepId = CreateBoundaryConnectorBrepId(rPatchModelPart, global_start, global_end);
            segment.BrepModelPartFullName = rPatchModelPart.FullName();
        }

        segment.Start = global_start;
        segment.End = global_end;
        segments.push_back(std::move(segment));
    }

    return segments;
}

bool HasLayer(
    const Node& rNode,
    const std::string& rLayerName)
{
    const auto& r_connected_layers = rNode.GetValue(CONNECTED_LAYERS);
    return std::find(r_connected_layers.begin(), r_connected_layers.end(), rLayerName) != r_connected_layers.end();
}

std::string TrimWhitespace(const std::string& rValue)
{
    const auto first = std::find_if_not(
        rValue.begin(),
        rValue.end(),
        [](const unsigned char c) { return std::isspace(c) != 0; });
    const auto last = std::find_if_not(
        rValue.rbegin(),
        rValue.rend(),
        [](const unsigned char c) { return std::isspace(c) != 0; }).base();

    if (first >= last) {
        return "";
    }

    return std::string(first, last);
}

std::string GetConditionNameForConnectedLayer(
    const Node& rNode,
    const std::string& rLayerName)
{
    const auto& r_connected_layers = rNode.GetValue(CONNECTED_LAYERS);
    const auto& r_connected_conditions = rNode.GetValue(CONNECTED_CONDITIONS);
    for (IndexType i = 0; i < r_connected_layers.size(); ++i) {
        if (r_connected_layers[i] != rLayerName) {
            continue;
        }

        KRATOS_ERROR_IF(i >= r_connected_conditions.size())
            << "NurbsGeometryModelerGapSbm: node #" << rNode.Id()
            << " has CONNECTED_LAYERS entry '" << rLayerName
            << "' without a matching CONNECTED_CONDITIONS entry." << std::endl;
        return r_connected_conditions[i];
    }

    return "";
}

bool NodeHasNonCouplingCondition(
    const ModelPart& rModelPart,
    const IndexType NodeId)
{
    for (const auto& r_condition : rModelPart.Conditions()) {
        const std::string layer_name = r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
        if (layer_name == "COUPLING_SIDE") {
            continue;
        }

        const auto& r_geometry = r_condition.GetGeometry();
        for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
            if (r_geometry[i].Id() == NodeId) {
                return true;
            }
        }
    }

    return false;
}

Geometry<Node>::Pointer GetBrepGeometryById(
    ModelPart& rPatchModelPart,
    const IndexType BrepId)
{
    if (rPatchModelPart.HasGeometry(BrepId)) {
        return rPatchModelPart.pGetGeometry(BrepId);
    }

    ModelPart& r_root_model_part = rPatchModelPart.GetRootModelPart();
    KRATOS_ERROR_IF_NOT(r_root_model_part.HasGeometry(BrepId))
        << "NurbsGeometryModelerGapSbm: patch model part '" << rPatchModelPart.FullName()
        << "' has no geometry with BREP_ID " << BrepId << "." << std::endl;

    return r_root_model_part.pGetGeometry(BrepId);
}

void SetCouplingNodeData(
    Node& rNode,
    const array_1d<double, 3>& rTangent,
    const Vector& rNormal,
    const std::string& rConditionName,
    const int BrepId,
    const std::string& rBrepModelPartFullName)
{
    rNode.SetValue(NORMAL, rNormal);
    rNode.SetValue(LOCAL_TANGENT, rTangent);
    rNode.SetValue(CURVATURE, 0.0);

    std::vector<std::string> connected_layers{"COUPLING_SIDE"};
    std::vector<std::string> connected_conditions{rConditionName};
    rNode.SetValue(CONNECTED_LAYERS, connected_layers);
    rNode.SetValue(CONNECTED_CONDITIONS, connected_conditions);
    SetBrepMetadata(rNode, BrepId, rBrepModelPartFullName);
}

void EnsureCouplingNodeMembership(
    ModelPart& rTargetLoop,
    ModelPart& rInterfaceVertices,
    ModelPart& rLayerModelPart,
    ModelPart& rLoopSubModelPart,
    const Node::Pointer& pNode)
{
    if (!rTargetLoop.HasNode(pNode->Id())) {
        rTargetLoop.AddNode(pNode);
    }
    if (!rInterfaceVertices.HasNode(pNode->Id())) {
        rInterfaceVertices.AddNode(pNode);
    }
    if (!rLayerModelPart.HasNode(pNode->Id())) {
        rLayerModelPart.AddNode(pNode);
    }
    if (!rLoopSubModelPart.HasNode(pNode->Id())) {
        rLoopSubModelPart.AddNode(pNode);
    }
}

void RebuildSkinLoopsForLocalRefinement(
    Model& rModel,
    ModelPart& rSkinModelPart,
    const std::string& rSkinModelPartName,
    const bool RebuildInnerLoop,
    const bool RebuildOuterLoop)
{
    const auto refinement_skin_names = CollectRefinementSkinModelPartNames(rModel, rSkinModelPartName);
    Vector aggregated_refinement_knot_span_sizes;
    bool has_aggregated_refinement_knot_span_sizes = false;
    KRATOS_ERROR_IF(refinement_skin_names.empty())
        << "NurbsGeometryModelerGapSbm: local refinement requested for '" << rSkinModelPartName
        << "' but no refinement patch skin model parts were found." << std::endl;

    auto rebuild_one_loop = [&](const std::string& rTargetLoopName, const std::string& rSurrogateBrepLayerName) {
        KRATOS_ERROR_IF_NOT(rSkinModelPart.HasSubModelPart(rTargetLoopName))
            << "NurbsGeometryModelerGapSbm: missing skin submodelpart '" << rTargetLoopName
            << "' in '" << rSkinModelPart.FullName() << "'." << std::endl;

        ModelPart& r_target_loop = rSkinModelPart.GetSubModelPart(rTargetLoopName);
        if (!HasCouplingSideSegmentsRecursive(r_target_loop)) {
            return;
        }

        KRATOS_ERROR_IF_NOT(r_target_loop.HasSubModelPart("interface_vertices"))
            << "NurbsGeometryModelerGapSbm: missing 'interface_vertices' in skin loop '"
            << r_target_loop.FullName() << "'." << std::endl;

        ModelPart& r_interface_vertices = r_target_loop.GetSubModelPart("interface_vertices");
        std::vector<Node::Pointer> pending_nodes;
        std::string coupling_condition_name;
        for (auto& r_node : r_interface_vertices.Nodes()) {
            if (!HasLayer(r_node, "COUPLING_SIDE")) {
                continue;
            }

            const std::string condition_name = GetConditionNameForConnectedLayer(r_node, "COUPLING_SIDE");
            KRATOS_ERROR_IF(condition_name.empty())
                << "NurbsGeometryModelerGapSbm: node #" << r_node.Id()
                << " is connected to COUPLING_SIDE but has an empty condition name." << std::endl;

            if (coupling_condition_name.empty()) {
                coupling_condition_name = condition_name;
            } else {
                KRATOS_ERROR_IF(coupling_condition_name != condition_name)
                    << "NurbsGeometryModelerGapSbm: inconsistent COUPLING_SIDE condition names in '"
                    << r_interface_vertices.FullName() << "': '" << coupling_condition_name
                    << "' and '" << condition_name << "'." << std::endl;
            }
            pending_nodes.push_back(r_interface_vertices.pGetNode(r_node.Id()));
        }

        KRATOS_ERROR_IF(pending_nodes.empty())
            << "NurbsGeometryModelerGapSbm: no pending interface vertices connected to COUPLING_SIDE were found in '"
            << r_interface_vertices.FullName() << "'." << std::endl;

        std::vector<PendingNodeMatch> pending_matches;
        pending_matches.reserve(pending_nodes.size());
        for (auto& p_pending_node : pending_nodes) {
            PendingNodeMatch match;
            match.pPendingNode = p_pending_node;
            pending_matches.push_back(std::move(match));
        }

        auto find_pending_match = [&](const IndexType NodeId) -> PendingNodeMatch* {
            for (auto& r_match : pending_matches) {
                if (r_match.pPendingNode->Id() == NodeId) {
                    return &r_match;
                }
            }
            return nullptr;
        };

        IndexType next_node_id = NextNodeId(rSkinModelPart);

        ModelPart& r_root_model_part = rSkinModelPart.GetRootModelPart();
        for (auto& r_condition : r_root_model_part.Conditions()) {
            r_condition.Set(TO_ERASE, false);
        }
        for (auto& r_node : r_root_model_part.Nodes()) {
            r_node.Set(TO_ERASE, false);
        }

        std::unordered_set<IndexType> candidate_node_ids;
        MarkCouplingSideConditionsForErase(r_target_loop, candidate_node_ids);

        r_root_model_part.RemoveConditionsFromAllLevels(TO_ERASE);

        for (const auto node_id : candidate_node_ids) {
            if (!r_root_model_part.HasNode(node_id)) {
                continue;
            }
            if (!NodeHasNonCouplingCondition(r_root_model_part, node_id)) {
                r_root_model_part.GetNode(node_id).Set(TO_ERASE, true);
            }
        }

        r_root_model_part.RemoveNodesFromAllLevels(TO_ERASE);

        ModelPart& r_layer_model_part = r_target_loop.HasSubModelPart("COUPLING_SIDE")
            ? r_target_loop.GetSubModelPart("COUPLING_SIDE")
            : r_target_loop.CreateSubModelPart("COUPLING_SIDE");
        ModelPart& r_loop_sub_model_part = r_target_loop.HasSubModelPart("0")
            ? r_target_loop.GetSubModelPart("0")
            : r_target_loop.CreateSubModelPart("0");

        EnsurePropertiesAvailable(rSkinModelPart);
        auto p_properties = rSkinModelPart.pGetProperties(0);
        IndexType next_condition_id = NextConditionId(rSkinModelPart);

        auto mark_exact_pending_match = [&](const Node::Pointer& pMatchedNode) {
            if (auto* p_match = find_pending_match(pMatchedNode->Id())) {
                p_match->HasExactMatch = true;
                p_match->IsResolved = true;
            }
        };

        auto find_or_create_node = [&](const Node& rSourceNode,
                                       const SkinCoordinatesType& rCoordinates,
                                       const array_1d<double, 3>& rTangent,
                                       const Vector& rNormal,
                                       const int BrepId,
                                       const std::string& rBrepModelPartFullName) {
            for (auto& p_node : pending_nodes) {
                if (IsCloseCoordinates(p_node->Coordinates(), rCoordinates)) {
                    mark_exact_pending_match(p_node);
                    p_node->SetValue(NORMAL, rNormal);
                    p_node->SetValue(LOCAL_TANGENT, rTangent);
                    p_node->SetValue(CURVATURE, 0.0);
                    SetBrepMetadata(*p_node, BrepId, rBrepModelPartFullName);
                    EnsureCouplingNodeMembership(
                        r_target_loop,
                        r_interface_vertices,
                        r_layer_model_part,
                        r_loop_sub_model_part,
                        p_node);
                    return p_node;
                }
            }

            auto p_node = rSkinModelPart.CreateNewNode(
                next_node_id++,
                rCoordinates[0],
                rCoordinates[1],
                rCoordinates[2]);
            SetCouplingNodeData(
                *p_node,
                rTangent,
                rNormal,
                coupling_condition_name,
                BrepId,
                rBrepModelPartFullName);
            EnsureCouplingNodeMembership(
                r_target_loop,
                r_interface_vertices,
                r_layer_model_part,
                r_loop_sub_model_part,
                p_node);
            pending_nodes.push_back(p_node);
            return p_node;
        };

        std::size_t number_of_inserted_conditions = 0;
        for (const auto& r_refinement_skin_name : refinement_skin_names) {
            ModelPart& r_refinement_patch_model_part = rModel.GetModelPart(
                FindPatchModelPartNameFromSkinName(rModel, r_refinement_skin_name));
            KRATOS_ERROR_IF_NOT(r_refinement_patch_model_part.HasSubModelPart("surrogate_outer"))
                << "NurbsGeometryModelerGapSbm: refinement patch model part '"
                << r_refinement_patch_model_part.FullName()
                << "' has no 'surrogate_outer' submodelpart." << std::endl;

            ModelPart& r_surrogate_outer = r_refinement_patch_model_part.GetSubModelPart("surrogate_outer");
            const Vector current_patch_knot_span_sizes = r_surrogate_outer.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
            if (current_patch_knot_span_sizes.size() >= 2) {
                if (!has_aggregated_refinement_knot_span_sizes) {
                    aggregated_refinement_knot_span_sizes = current_patch_knot_span_sizes;
                    has_aggregated_refinement_knot_span_sizes = true;
                } else {
                    const std::size_t component_count = std::min(
                        aggregated_refinement_knot_span_sizes.size(),
                        current_patch_knot_span_sizes.size());
                    for (std::size_t i = 0; i < component_count; ++i) {
                        aggregated_refinement_knot_span_sizes[i] = std::min(
                            aggregated_refinement_knot_span_sizes[i],
                            current_patch_knot_span_sizes[i]);
                    }
                }
            }

            std::unordered_map<IndexType, Node::Pointer> patch_component_nodes;
            std::unordered_map<IndexType, Node::Pointer> patch_component_source_nodes;
            std::unordered_map<IndexType, std::size_t> patch_component_node_degrees;

            auto refinement_box_corneres = r_surrogate_outer.GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);

            const double box_boundary_tolerance = 1.0e-12;
            auto is_on_box_border = [&](const auto& rCoordinates) {
                const bool on_left_or_right =
                    (std::abs(rCoordinates[0] - refinement_box_corneres[0][0]) < box_boundary_tolerance ||
                     std::abs(rCoordinates[0] - refinement_box_corneres[0][1]) < box_boundary_tolerance) &&
                    (rCoordinates[1] >= refinement_box_corneres[1][0] - box_boundary_tolerance &&
                     rCoordinates[1] <= refinement_box_corneres[1][1] + box_boundary_tolerance);

                const bool on_bottom_or_top =
                    (std::abs(rCoordinates[1] - refinement_box_corneres[1][0]) < box_boundary_tolerance ||
                     std::abs(rCoordinates[1] - refinement_box_corneres[1][1]) < box_boundary_tolerance) &&
                    (rCoordinates[0] >= refinement_box_corneres[0][0] - box_boundary_tolerance &&
                     rCoordinates[0] <= refinement_box_corneres[0][1] + box_boundary_tolerance);

                return on_left_or_right || on_bottom_or_top;
            };

            auto ReorderNodesOnBoxBoundary = [&](const Node::Pointer& pOpenEndpoint,
                                                 const Node::Pointer& pPendingNode) -> std::array<Node::Pointer, 2> {
                KRATOS_ERROR_IF_NOT(is_on_box_border(pOpenEndpoint->Coordinates()))
                    << "NurbsGeometryModelerGapSbm: open endpoint #" << pOpenEndpoint->Id()
                    << " is not on the local-refinement box boundary. Coordinates: "
                    << pOpenEndpoint->Coordinates() << std::endl;
                KRATOS_ERROR_IF_NOT(is_on_box_border(pPendingNode->Coordinates()))
                    << "NurbsGeometryModelerGapSbm: pending node #" << pPendingNode->Id()
                    << " is not on the local-refinement box boundary. Coordinates: "
                    << pPendingNode->Coordinates() << std::endl;

                const auto& r_open_coordinates = pOpenEndpoint->Coordinates();
                const auto& r_pending_coordinates = pPendingNode->Coordinates();
                const bool same_x =
                    std::abs(r_open_coordinates[0] - r_pending_coordinates[0]) < box_boundary_tolerance;
                const bool same_y =
                    std::abs(r_open_coordinates[1] - r_pending_coordinates[1]) < box_boundary_tolerance;

                KRATOS_ERROR_IF(same_x == same_y)
                    << "NurbsGeometryModelerGapSbm: cannot identify a unique common box boundary for nodes #"
                    << pOpenEndpoint->Id() << " and #" << pPendingNode->Id()
                    << ". Coordinates: " << r_open_coordinates << " and "
                    << r_pending_coordinates << std::endl;

                const bool open_first_by_smaller_coordinate =
                    same_x ? r_open_coordinates[1] < r_pending_coordinates[1]
                           : r_open_coordinates[0] < r_pending_coordinates[0];

                bool use_smaller_coordinate_first = false;
                if (same_x) {
                    const bool on_left =
                        std::abs(r_open_coordinates[0] - refinement_box_corneres[0][0]) < box_boundary_tolerance &&
                        std::abs(r_pending_coordinates[0] - refinement_box_corneres[0][0]) < box_boundary_tolerance;
                    const bool on_right =
                        std::abs(r_open_coordinates[0] - refinement_box_corneres[0][1]) < box_boundary_tolerance &&
                        std::abs(r_pending_coordinates[0] - refinement_box_corneres[0][1]) < box_boundary_tolerance;
                    KRATOS_ERROR_IF(!on_left && !on_right)
                        << "NurbsGeometryModelerGapSbm: nodes #" << pOpenEndpoint->Id()
                        << " and #" << pPendingNode->Id()
                        << " share x but are not on the left or right box boundary."
                        << std::endl;
                    use_smaller_coordinate_first = on_right;
                } else {
                    const bool on_bottom =
                        std::abs(r_open_coordinates[1] - refinement_box_corneres[1][0]) < box_boundary_tolerance &&
                        std::abs(r_pending_coordinates[1] - refinement_box_corneres[1][0]) < box_boundary_tolerance;
                    const bool on_top =
                        std::abs(r_open_coordinates[1] - refinement_box_corneres[1][1]) < box_boundary_tolerance &&
                        std::abs(r_pending_coordinates[1] - refinement_box_corneres[1][1]) < box_boundary_tolerance;
                    KRATOS_ERROR_IF(!on_bottom && !on_top)
                        << "NurbsGeometryModelerGapSbm: nodes #" << pOpenEndpoint->Id()
                        << " and #" << pPendingNode->Id()
                        << " share y but are not on the bottom or top box boundary."
                        << std::endl;
                    use_smaller_coordinate_first = on_bottom;
                }

                const bool open_is_second =
                    use_smaller_coordinate_first == open_first_by_smaller_coordinate;
                if (open_is_second) {
                    return {pPendingNode, pOpenEndpoint};
                }
                return {pOpenEndpoint, pPendingNode};
            };

            for (const auto& r_condition : r_surrogate_outer.Conditions()) {
                const std::string condition_layer_name = r_condition.Has(LAYER_NAME)
                    ? TrimWhitespace(r_condition.GetValue(LAYER_NAME))
                    : "";
                if (condition_layer_name != rSurrogateBrepLayerName) {
                    continue;
                }

                auto cond_center = r_condition.GetGeometry().Center();
                if (!is_on_box_border(cond_center)) {
                    continue;
                }

                KRATOS_ERROR_IF_NOT(r_condition.Has(BREP_ID))
                    << "NurbsGeometryModelerGapSbm: surrogate_outer condition #" << r_condition.Id()
                    << " in '" << r_surrogate_outer.FullName() << "' has no BREP_ID." << std::endl;

                const int brep_id = r_condition.GetValue(BREP_ID);
                const std::string brep_model_part_full_name = r_condition.Has(BREP_MODEL_PART_FULL_NAME)
                    ? r_condition.GetValue(BREP_MODEL_PART_FULL_NAME)
                    : r_refinement_patch_model_part.FullName();
                GetBrepGeometryById(
                    r_refinement_patch_model_part,
                    static_cast<IndexType>(brep_id));

                const auto& r_geometry = r_condition.GetGeometry();
                KRATOS_ERROR_IF(r_geometry.PointsNumber() < 2)
                    << "NurbsGeometryModelerGapSbm: surrogate_outer condition #" << r_condition.Id()
                    << " in '" << r_surrogate_outer.FullName() << "' has fewer than 2 nodes." << std::endl;

                SkinCoordinatesType coord_1 = r_geometry[1].Coordinates();
                SkinCoordinatesType coord_2 = r_geometry[0].Coordinates();
                array_1d<double, 3> tangent = coord_2 - coord_1;
                const double tangent_norm = norm_2(tangent);
                KRATOS_ERROR_IF(tangent_norm <= 1.0e-16)
                    << "NurbsGeometryModelerGapSbm: zero-length local-refinement coupling segment in '"
                    << r_surrogate_outer.FullName() << "' condition #" << r_condition.Id() << "." << std::endl;
                tangent /= tangent_norm;

                Vector normal = ZeroVector(3);
                normal[0] = tangent[1];
                normal[1] = -tangent[0];

                auto p_node_1 = find_or_create_node(
                    r_geometry[1],
                    coord_1,
                    tangent,
                    normal,
                    brep_id,
                    brep_model_part_full_name);
                auto p_node_2 = find_or_create_node(
                    r_geometry[0],
                    coord_2,
                    tangent,
                    normal,
                    brep_id,
                    brep_model_part_full_name);

                patch_component_nodes[p_node_1->Id()] = p_node_1;
                patch_component_nodes[p_node_2->Id()] = p_node_2;
                patch_component_source_nodes[p_node_1->Id()] = r_geometry.pGetPoint(1);
                patch_component_source_nodes[p_node_2->Id()] = r_geometry.pGetPoint(0);
                ++patch_component_node_degrees[p_node_1->Id()];
                ++patch_component_node_degrees[p_node_2->Id()];

                auto p_new_condition = rSkinModelPart.CreateNewCondition(
                    "LineCondition2D2N",
                    next_condition_id++,
                    {{p_node_1->Id(), p_node_2->Id()}},
                    p_properties);
                p_new_condition->SetValue(LAYER_NAME, "COUPLING_SIDE");
                p_new_condition->SetValue(CONDITION_NAME, coupling_condition_name);
                SetBrepMetadata(*p_new_condition, brep_id, brep_model_part_full_name);

                r_target_loop.AddCondition(p_new_condition);
                r_layer_model_part.AddCondition(p_new_condition);
                r_loop_sub_model_part.AddCondition(p_new_condition);
                ++number_of_inserted_conditions;
            }

            std::vector<Node::Pointer> open_component_endpoints;
            open_component_endpoints.reserve(patch_component_node_degrees.size());
            for (const auto& r_entry : patch_component_node_degrees) {
                if (r_entry.second != 1) {
                    continue;
                }

                auto it_node = patch_component_nodes.find(r_entry.first);
                if (it_node == patch_component_nodes.end()) {
                    continue;
                }

                PendingNodeMatch* p_match = find_pending_match(r_entry.first);
                if (p_match != nullptr && p_match->HasExactMatch) {
                    continue;
                }

                open_component_endpoints.push_back(it_node->second);
            }

            struct EndpointPendingCandidate
            {
                double Distance = std::numeric_limits<double>::max();
                IndexType EndpointId = 0;
                IndexType PendingNodeId = 0;
            };

            std::vector<EndpointPendingCandidate> endpoint_pending_candidates;
            for (const auto& p_open_endpoint : open_component_endpoints) {
                for (const auto& r_match : pending_matches) {
                    if (r_match.IsResolved) {
                        continue;
                    }

                    EndpointPendingCandidate candidate;
                    candidate.Distance = norm_2(
                        p_open_endpoint->Coordinates() - r_match.pPendingNode->Coordinates());
                    candidate.EndpointId = p_open_endpoint->Id();
                    candidate.PendingNodeId = r_match.pPendingNode->Id();
                    endpoint_pending_candidates.push_back(candidate);
                }
            }

            std::sort(
                endpoint_pending_candidates.begin(),
                endpoint_pending_candidates.end(),
                [](const EndpointPendingCandidate& rLeft, const EndpointPendingCandidate& rRight) {
                    return rLeft.Distance < rRight.Distance;
                });

            std::unordered_set<IndexType> paired_endpoint_ids;
            std::unordered_set<IndexType> paired_pending_node_ids;
            for (const auto& r_candidate : endpoint_pending_candidates) {
                if (paired_endpoint_ids.find(r_candidate.EndpointId) != paired_endpoint_ids.end()) {
                    continue;
                }
                if (paired_pending_node_ids.find(r_candidate.PendingNodeId) != paired_pending_node_ids.end()) {
                    continue;
                }

                auto it_endpoint = patch_component_nodes.find(r_candidate.EndpointId);
                PendingNodeMatch* p_match = find_pending_match(r_candidate.PendingNodeId);
                if (it_endpoint == patch_component_nodes.end() || p_match == nullptr || p_match->IsResolved) {
                    continue;
                }

                auto p_open_endpoint = it_endpoint->second;
                auto p_pending_node = p_match->pPendingNode;

                const auto ordered_nodes = ReorderNodesOnBoxBoundary(p_open_endpoint, p_pending_node);
                auto p_node_1 = ordered_nodes[0];
                auto p_node_2 = ordered_nodes[1];

                const int connector_brep_id = CreateLocalRefinementClosureConnectorBrepId(
                    r_refinement_patch_model_part,
                    p_node_1->Coordinates(),
                    p_node_2->Coordinates());
                (void)connector_brep_id;

                EnsureCouplingNodeMembership(
                    r_target_loop,
                    r_interface_vertices,
                    r_layer_model_part,
                    r_loop_sub_model_part,
                    p_pending_node);
                EnsureCouplingNodeMembership(
                    r_target_loop,
                    r_interface_vertices,
                    r_layer_model_part,
                    r_loop_sub_model_part,
                    p_open_endpoint);

                auto p_new_condition = rSkinModelPart.CreateNewCondition(
                    "LineCondition2D2N",
                    next_condition_id++,
                    {{p_node_1->Id(), p_node_2->Id()}},
                    p_properties);
                p_new_condition->SetValue(LAYER_NAME, "COUPLING_SIDE");
                p_new_condition->SetValue(CONDITION_NAME, coupling_condition_name);
                SetBrepMetadata(*p_new_condition, connector_brep_id, r_refinement_patch_model_part.FullName());
                
                auto it_source_node = patch_component_source_nodes.find(p_open_endpoint->Id());
                KRATOS_ERROR_IF(it_source_node == patch_component_source_nodes.end())
                    << "NurbsGeometryModelerGapSbm: missing surrogate source node for copied endpoint #"
                    << p_open_endpoint->Id() << " while rebuilding skin loop '"
                    << r_target_loop.FullName() << "'." << std::endl;

                const auto& r_surrogate_source_node = *(it_source_node->second);
                const bool has_neighbour_geometries = r_surrogate_source_node.Has(NEIGHBOUR_GEOMETRIES);
                if (has_neighbour_geometries) {
                    p_new_condition->SetValue(
                        NEIGHBOUR_GEOMETRIES,
                        r_surrogate_source_node.GetValue(NEIGHBOUR_GEOMETRIES));
                }
                else
                {
                    KRATOS_ERROR                        << "NurbsGeometryModelerGapSbm: neither node #" << p_node_1->Id() << " nor node #" << p_node_2->Id()
                        << " has NEIGHBOUR_GEOMETRIES data to copy after pairing." << std::endl;
                }
                   
                r_target_loop.AddCondition(p_new_condition);
                r_layer_model_part.AddCondition(p_new_condition);
                r_loop_sub_model_part.AddCondition(p_new_condition);
                ++number_of_inserted_conditions;

                p_match->IsResolved = true;
                paired_endpoint_ids.insert(r_candidate.EndpointId);
                paired_pending_node_ids.insert(r_candidate.PendingNodeId);
            }

            KRATOS_ERROR_IF(paired_endpoint_ids.size() != open_component_endpoints.size())
                << "NurbsGeometryModelerGapSbm: failed to close all local-refinement coupling endpoints for patch '"
                << r_refinement_patch_model_part.FullName() << "' while rebuilding skin loop '"
                << r_target_loop.FullName() << "'." << std::endl;
        }

        KRATOS_ERROR_IF(number_of_inserted_conditions == 0)
            << "NurbsGeometryModelerGapSbm: no surrogate_outer conditions with LAYER_NAME='"
            << rSurrogateBrepLayerName << "' were found while rebuilding skin loop '"
            << r_target_loop.FullName() << "'." << std::endl;

        for (const auto& r_match : pending_matches) {
            KRATOS_ERROR_IF_NOT(r_match.IsResolved)
                << "NurbsGeometryModelerGapSbm: pending interface vertex #"
                << r_match.pPendingNode->Id()
                << " was not resolved while rebuilding the local-refinement coupling loop '"
                << r_target_loop.FullName() << "'." << std::endl;
        }
    };

    if (RebuildInnerLoop) {
        rebuild_one_loop("inner", "COUPLING_CONDITION_INNER");
    }
    if (RebuildOuterLoop) {
        rebuild_one_loop("outer", "COUPLING_CONDITION_OUTER");
    }

    if (has_aggregated_refinement_knot_span_sizes) {
        rSkinModelPart.SetValue(KNOT_SPAN_SIZES, aggregated_refinement_knot_span_sizes);
        if (rSkinModelPart.HasSubModelPart("inner")) {
            rSkinModelPart.GetSubModelPart("inner").SetValue(
                KNOT_SPAN_SIZES,
                aggregated_refinement_knot_span_sizes);
        }
        if (rSkinModelPart.HasSubModelPart("outer")) {
            rSkinModelPart.GetSubModelPart("outer").SetValue(
                KNOT_SPAN_SIZES,
                aggregated_refinement_knot_span_sizes);
        }
    }
}

} // namespace

///@name Stages
///@{

///@}
///@name Private Operations
///@{
void NurbsGeometryModelerGapSbm::CreateAndAddRegularGrid2D(
    ModelPart& rModelPart, 
    const Point& rPointAXyz, 
    const Point& rPointBXyz,
    const Point& rPointAUvw, 
    const Point& rPointBUvw, 
    const std::size_t OrderU, 
    const std::size_t OrderV, 
    const std::size_t NumKnotSpansU, 
    const std::size_t NumKnotSpansV, 
    const bool AddSurfaceToModelPart)
{   

    // Call the CreateAndAddRegularGrid2D method of the base class NurbsGeometryModeler
    NurbsGeometryModeler::CreateAndAddRegularGrid2D(
        rModelPart,
        rPointAXyz,
        rPointBXyz,
        rPointAUvw,
        rPointBUvw,
        OrderU,
        OrderV,
        NumKnotSpansU,
        NumKnotSpansV,
        false);
        
    // Create the Domain/Iga Model Part
    const std::string iga_model_part_name = mParameters["model_part_name"].GetString();
    ModelPart& r_iga_model_part = mpModel->HasModelPart(iga_model_part_name)
                                ? mpModel->GetModelPart(iga_model_part_name)
                                : mpModel->CreateModelPart(iga_model_part_name);

    // compute unique_knot_vector_u
    Vector unique_knot_vector_u(2+(NumKnotSpansU-1));
    unique_knot_vector_u[0] = mKnotVectorU[0]; 
    unique_knot_vector_u[NumKnotSpansU] = mKnotVectorU[mKnotVectorU.size()-1];
    for (std::size_t i_knot_insertion = 0; i_knot_insertion < NumKnotSpansU-1; i_knot_insertion++) {
        unique_knot_vector_u[i_knot_insertion+1] = mInsertKnotsU[i_knot_insertion];
    }
    // compute unique_knot_vector_v
    Vector unique_knot_vector_v(2+(NumKnotSpansV-1));
    unique_knot_vector_v[0] = mKnotVectorV[0]; 
    unique_knot_vector_v[NumKnotSpansV] = mKnotVectorV[mKnotVectorV.size()-1];
    for (std::size_t i_knot_insertion = 0; i_knot_insertion < NumKnotSpansV-1; i_knot_insertion++) {
        unique_knot_vector_v[i_knot_insertion+1] = mInsertKnotsV[i_knot_insertion];
    }
    // Set the value of the knot vectors
    r_iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
    r_iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);

    // If neither skin_inner nor skin_outer exists
    if (!(mParameters.Has("skin_model_part_inner_initial_name") || mParameters.Has("skin_model_part_outer_initial_name"))){
        
        Vector knot_step_uv= ZeroVector(2);
        knot_step_uv[0] = std::abs(unique_knot_vector_u[std::ceil(unique_knot_vector_u.size()/2) +1] - unique_knot_vector_u[std::ceil(unique_knot_vector_u.size()/2)] ) ;
        knot_step_uv[1] = std::abs(unique_knot_vector_v[std::ceil(unique_knot_vector_v.size()/2) +1] - unique_knot_vector_v[std::ceil(unique_knot_vector_v.size()/2)] ) ;

        // saving the knot span sizes
        r_iga_model_part.SetValue(KNOT_SPAN_SIZES, knot_step_uv);

        // Create the breps for the outer sbm boundary
        CreateBrepsSbmUtilities<Node, Point, false> breps_sbm_utilities(mEchoLevel);
        breps_sbm_utilities.CreateSurrogateBoundary(mpSurface, rPointAUvw, rPointBUvw, rModelPart);

        KRATOS_WARNING("None of the 'skin_model_part_name' have not been defined ") << 
                        "in the nurbs_geometry_modeler_sbm in the project paramer json" << std::endl;
        return;
    }

    // Create the True Model Part -> contains all the true boundary features
    std::string skin_model_part_name;

    // Retrieve skin_model_part_inner_initial_name if it exists
    std::string skin_model_part_inner_initial_name = "skin_model_part_outer_initial_name";
    if (mParameters.Has("skin_model_part_inner_initial_name")) {
        skin_model_part_inner_initial_name = mParameters["skin_model_part_inner_initial_name"].GetString();
    }

    // Retrieve skin_model_part_outer_initial_name if it exists;
    std::string skin_model_part_outer_initial_name = "skin_model_part_outer_initial_name";
    if (mParameters.Has("skin_model_part_outer_initial_name")) {
        skin_model_part_outer_initial_name = mParameters["skin_model_part_outer_initial_name"].GetString();
    }

    // Create the surrogate sub model parts inner and outer
    KRATOS_ERROR_IF(r_iga_model_part.HasSubModelPart("surrogate_inner"))
        << "NurbsGeometryModelerGapSbm: submodelpart 'surrogate_inner' already exists in '"
        << r_iga_model_part.FullName() << "'." << std::endl;
    KRATOS_ERROR_IF(r_iga_model_part.HasSubModelPart("surrogate_outer"))
        << "NurbsGeometryModelerGapSbm: submodelpart 'surrogate_outer' already exists in '"
        << r_iga_model_part.FullName() << "'." << std::endl;
    ModelPart& r_surrogate_sub_model_part_inner = r_iga_model_part.CreateSubModelPart("surrogate_inner");
    ModelPart& r_surrogate_sub_model_part_outer = r_iga_model_part.CreateSubModelPart("surrogate_outer");
    EnsurePropertiesAvailable(r_surrogate_sub_model_part_inner);
    EnsurePropertiesAvailable(r_surrogate_sub_model_part_outer);
    
    if (mParameters.Has("skin_model_part_name"))
        skin_model_part_name = mParameters["skin_model_part_name"].GetString();
    else
        KRATOS_ERROR << "The skin_model_part name '" << skin_model_part_name << "' was not defined in the project parameters.\n" << std::endl;

    // inner
    auto& skin_inner_initial = mpModel->HasModelPart(skin_model_part_inner_initial_name)
        ? mpModel->GetModelPart(skin_model_part_inner_initial_name)
        : mpModel->CreateModelPart(skin_model_part_inner_initial_name);
    // outer
    auto& skin_outer_initial = mpModel->HasModelPart(skin_model_part_outer_initial_name)
        ? mpModel->GetModelPart(skin_model_part_outer_initial_name)
        : mpModel->CreateModelPart(skin_model_part_outer_initial_name);
    
    // Skin model part refined after Snake Process
    KRATOS_ERROR_IF(mpModel->HasModelPart(skin_model_part_name))
        << "NurbsGeometryModelerGapSbm: model part '" << skin_model_part_name
        << "' already exists. Refusing to overwrite." << std::endl;
    ModelPart& r_skin_model_part = mpModel->CreateModelPart(skin_model_part_name);
    KRATOS_ERROR_IF(r_skin_model_part.HasSubModelPart("inner"))
        << "NurbsGeometryModelerGapSbm: submodelpart 'inner' already exists in '"
        << r_skin_model_part.FullName() << "'." << std::endl;
    KRATOS_ERROR_IF(r_skin_model_part.HasSubModelPart("outer"))
        << "NurbsGeometryModelerGapSbm: submodelpart 'outer' already exists in '"
        << r_skin_model_part.FullName() << "'." << std::endl;
    r_skin_model_part.CreateSubModelPart("inner");
    r_skin_model_part.CreateSubModelPart("outer");
    EnsurePropertiesAvailable(r_skin_model_part);
    EnsurePropertiesAvailable(r_skin_model_part.GetSubModelPart("inner"));
    EnsurePropertiesAvailable(r_skin_model_part.GetSubModelPart("outer"));

    const auto& r_inner_sizes = skin_inner_initial.GetValue(KNOT_SPAN_SIZES);
    const auto& r_outer_sizes = skin_outer_initial.GetValue(KNOT_SPAN_SIZES);
    const bool inner_has_sizes = r_inner_sizes.size() >= 2;
    const bool outer_has_sizes = r_outer_sizes.size() >= 2;

    const Vector& r_fallback_sizes = inner_has_sizes ? r_inner_sizes : r_outer_sizes;
    r_skin_model_part.GetSubModelPart("inner").SetValue(KNOT_SPAN_SIZES, r_fallback_sizes); // pass span sizes to skin for consistent curve generation
    r_skin_model_part.GetSubModelPart("outer").SetValue(KNOT_SPAN_SIZES, r_fallback_sizes);
    r_skin_model_part.SetValue(KNOT_SPAN_SIZES, r_fallback_sizes);

    // Create the parameters for the SnakeSbmProcess
    Kratos::Parameters snake_parameters;
    snake_parameters.AddString("model_part_name", iga_model_part_name);
    snake_parameters.AddString("skin_model_part_name", skin_model_part_name);
    snake_parameters.AddDouble("echo_level", mEchoLevel); //FIXME:
    snake_parameters.AddString("skin_model_part_inner_initial_name", skin_model_part_inner_initial_name);
    snake_parameters.AddString("skin_model_part_outer_initial_name", skin_model_part_outer_initial_name);
    snake_parameters.AddString("gap_element_name", mParameters["gap_element_name"].GetString());
    snake_parameters.AddString("gap_interface_condition_name", mParameters["gap_interface_condition_name"].GetString());
    snake_parameters.AddString("gap_sbm_type", mParameters["gap_sbm_type"].GetString());
    if (mParameters.Has("lambda_inner"))
        snake_parameters.AddDouble("lambda_inner", mParameters["lambda_inner"].GetDouble());
    if (mParameters.Has("lambda_outer"))
        snake_parameters.AddDouble("lambda_outer", mParameters["lambda_outer"].GetDouble());
    if (mParameters.Has("number_of_inner_loops"))
        snake_parameters.AddDouble("number_of_inner_loops", mParameters["number_of_inner_loops"].GetInt());
    if (mParameters.Has("number_internal_divisions"))
        snake_parameters.AddDouble("number_internal_divisions", mParameters["number_internal_divisions"].GetInt());
    if (mParameters.Has("number_initial_points_if_importing_nurbs"))
        snake_parameters.AddInt("number_initial_points_if_importing_nurbs", mParameters["number_initial_points_if_importing_nurbs"].GetInt());
    if (mParameters.Has("gap_approximation_order"))
        snake_parameters.AddInt("gap_approximation_order", mParameters["gap_approximation_order"].GetInt());
    if (mParameters.Has("polynomial_order"))
        snake_parameters.AddVector("polynomial_order", mParameters["polynomial_order"].GetVector());
    if (mParameters.Has("use_for_multipatch"))
        snake_parameters.AddBool("use_for_multipatch", mParameters["use_for_multipatch"].GetBool());
    if (mParameters.Has("use_for_local_refinement"))
        snake_parameters.AddBool("use_for_local_refinement", mParameters["use_for_local_refinement"].GetBool());
    if (mParameters.Has("gap_relative_tolerance_for_subdivisions"))
        snake_parameters.AddDouble("gap_relative_tolerance_for_subdivisions", mParameters["gap_relative_tolerance_for_subdivisions"].GetDouble());
    if (mParameters.Has("number_of_interpolation_levels"))
        snake_parameters.AddInt("number_of_interpolation_levels", mParameters["number_of_interpolation_levels"].GetInt());
    if (mParameters.Has("create_surr_outer_from_surr_inner"))
        snake_parameters.AddBool("create_surr_outer_from_surr_inner", mParameters["create_surr_outer_from_surr_inner"].GetBool());
    if (mParameters.Has("create_surr_inner_from_surr_outer"))
        snake_parameters.AddBool("create_surr_inner_from_surr_outer", mParameters["create_surr_inner_from_surr_outer"].GetBool());


    // Create the surrogate_sub_model_part for inner and outer

    SnakeGapSbmProcess snake_sbm_process(*mpModel, snake_parameters);
    snake_sbm_process.ExecuteInitialize();

    // Create the breps for the outer sbm boundary
    CreateBrepsSbmUtilities<Node, Point, true> breps_sbm_utilities(mEchoLevel);
    breps_sbm_utilities.CreateSurrogateBoundary(
        mpSurface,
        r_surrogate_sub_model_part_inner,
        r_surrogate_sub_model_part_outer,
        rPointAUvw,
        rPointBUvw,
        r_iga_model_part);

    if (mParameters.Has("use_for_local_refinement") &&
        mParameters["use_for_local_refinement"].GetBool()) {
        RebuildSkinLoopsForLocalRefinement(
            *mpModel,
            r_skin_model_part,
            skin_model_part_name,
            skin_inner_initial.NumberOfNodes() > 0 || skin_inner_initial.NumberOfGeometries() > 0,
            skin_outer_initial.NumberOfNodes() > 0 || skin_outer_initial.NumberOfGeometries() > 0);
    }

    // if (mParameters.Has("use_for_local_refinement") &&
    //     mParameters["use_for_local_refinement"].GetBool()) 
    //     return;
    snake_sbm_process.Execute();

}

// 3D 
void NurbsGeometryModelerGapSbm::CreateAndAddRegularGrid3D( 
    ModelPart& rModelPart,
    const Point& rPointAXyz,
    const Point& rPointBXyz,
    const Point& rPointAUvw,
    const Point& rPointBUvw,
    const std::size_t OrderU,
    const std::size_t OrderV,
    const std::size_t OrderW,
    const std::size_t NumKnotSpansU,
    const std::size_t NumKnotSpansV,
    const std::size_t NumKnotSpansW,
    const bool AddVolumeToModelPart)
{   
    KRATOS_ERROR<< "CreateAndAddRegularGrid3D is not yet implemented for the NurbsGeometryModelerGapSbm."<<std::endl;

    // TODO: implement the 3D version of the Gap-SBM modeler (next PR)
    // // Call the CreateAndAddRegularGrid3D method of the base class NurbsGeometryModeler
    // NurbsGeometryModeler::CreateAndAddRegularGrid3D(
    //     rModelPart,
    //     rPointAXyz,
    //     rPointBXyz,
    //     rPointAUvw,
    //     rPointBUvw,
    //     OrderU,
    //     OrderV,
    //     OrderW,
    //     NumKnotSpansU,
    //     NumKnotSpansV,
    //     NumKnotSpansW,
    //     false);
             
    // // Create the Domain/Iga Model Part
    // const std::string iga_model_part_name = mParameters["model_part_name"].GetString();
    // ModelPart& r_iga_model_part = mpModel->HasModelPart(iga_model_part_name)
    //                             ? mpModel->GetModelPart(iga_model_part_name)
    //                             : mpModel->CreateModelPart(iga_model_part_name);

    // // Create the True Model Part -> contains all the true boundary features
    // std::string skin_model_part_name;
    // std::string skin_model_part_inner_initial_name = mParameters["skin_model_part_inner_initial_name"].GetString();
    // std::string skin_model_part_outer_initial_name = mParameters["skin_model_part_outer_initial_name"].GetString();

    // // Create the surrogate sub model parts inner and outer
    // // ModelPart& r_surrogate_sub_model_part_inner = r_iga_model_part.CreateSubModelPart("surrogate_inner");  // uncomment this line (next PR) 
    // // ModelPart& r_surrogate_sub_model_part_outer = r_iga_model_part.CreateSubModelPart("surrogate_outer");  // uncomment this line (next PR)

    // // If there is not neither skin_inner nor skin_outer throw a warning since you are using the sbm modeler
    // if (!(mParameters.Has("skin_model_part_inner_initial_name") || mParameters.Has("skin_model_part_outer_initial_name"))){
        
    //     // Create the breps for the outer sbm boundary
    //     CreateBrepsSbmUtilities<Node, Point> breps_sbm_utilities_3d(mEchoLevel);
    //     // TODO: NEXT PR CreateSurrogateBoundary with Volume
    //     // breps_sbm_utilities_3d.CreateSurrogateBoundary(mpVolume, rPointAUvw, rPointBUvw, rModelPart);

    //     KRATOS_WARNING("None of the 'skin_model_part_name' have not been defined ") << 
    //                     "in the nurbs_geometry_modeler_sbm in the project paramer json" << std::endl;
    //     return;
    // }
    
    // if (mParameters.Has("skin_model_part_name"))
    //     skin_model_part_name = mParameters["skin_model_part_name"].GetString();
    // else
    //     KRATOS_ERROR << "The skin_model_part name '" << skin_model_part_name << "' was not defined in the project parameters.\n" << std::endl;

    // // inner
    // mpModel->HasModelPart(skin_model_part_inner_initial_name)
    //     ? mpModel->GetModelPart(skin_model_part_inner_initial_name)
    //     : mpModel->CreateModelPart(skin_model_part_inner_initial_name);
    // // outer
    // mpModel->HasModelPart(skin_model_part_outer_initial_name)
    //     ? mpModel->GetModelPart(skin_model_part_outer_initial_name)
    //     : mpModel->CreateModelPart(skin_model_part_outer_initial_name);
    
    // // Skin model part refined after Snake Process
    // ModelPart& r_skin_model_part = mpModel->CreateModelPart(skin_model_part_name);
    // r_skin_model_part.CreateSubModelPart("inner");
    // r_skin_model_part.CreateSubModelPart("outer");
    
    
    // // compute unique_knot_vector_u
    // Vector unique_knot_vector_u(2+(NumKnotSpansU-1));
    // unique_knot_vector_u[0] = mKnotVectorU[0]; unique_knot_vector_u[NumKnotSpansU] = mKnotVectorU[mKnotVectorU.size()-1];
    // for (std::size_t i_knot_insertion = 0; i_knot_insertion < NumKnotSpansU-1; i_knot_insertion++) {
    //     unique_knot_vector_u[i_knot_insertion+1] = mInsertKnotsU[i_knot_insertion];
    // }

    // // compute unique_knot_vector_v
    // Vector unique_knot_vector_v(2+(NumKnotSpansV-1));
    // unique_knot_vector_v[0] = mKnotVectorV[0]; unique_knot_vector_v[NumKnotSpansV] = mKnotVectorV[mKnotVectorV.size()-1];
    // for (std::size_t i_knot_insertion = 0; i_knot_insertion < NumKnotSpansV-1; i_knot_insertion++) {
    //     unique_knot_vector_v[i_knot_insertion+1] = mInsertKnotsV[i_knot_insertion];
    // }

    // // compute unique_knot_vector_w
    // Vector unique_knot_vector_w(2+(NumKnotSpansW-1));
    // unique_knot_vector_w[0] = mKnotVectorW[0]; unique_knot_vector_w[NumKnotSpansW] = mKnotVectorW[mKnotVectorW.size()-1];
    // for (std::size_t i_knot_insertion = 0; i_knot_insertion < NumKnotSpansW-1; i_knot_insertion++) {
    //     unique_knot_vector_w[i_knot_insertion+1] = mInsertKnotsW[i_knot_insertion];
    // }

    // // Set the value of the knot vectors
    // r_iga_model_part.SetValue(KNOT_VECTOR_U, unique_knot_vector_u);
    // r_iga_model_part.SetValue(KNOT_VECTOR_V, unique_knot_vector_v);
    // r_iga_model_part.SetValue(KNOT_VECTOR_W, unique_knot_vector_w);

    // // Create the parameters for the SnakeSbmProcess
    // Kratos::Parameters snake_parameters;
    // snake_parameters.AddString("model_part_name", iga_model_part_name);
    // snake_parameters.AddString("skin_model_part_name", skin_model_part_name);
    // snake_parameters.AddDouble("echo_level", mEchoLevel);
    // snake_parameters.AddString("skin_model_part_inner_initial_name", skin_model_part_inner_initial_name);
    // snake_parameters.AddString("skin_model_part_outer_initial_name", skin_model_part_outer_initial_name);
    // snake_parameters.AddString("gap_element_name", mParameters["gap_element_name"].GetString());
    // snake_parameters.AddString("gap_interface_condition_name", mParameters["gap_interface_condition_name"].GetString());
    // snake_parameters.AddString("gap_sbm_type", mParameters["gap_sbm_type"].GetString());
    // if (mParameters.Has("lambda_inner"))
    //     snake_parameters.AddDouble("lambda_inner", mParameters["lambda_inner"].GetDouble());
    // if (mParameters.Has("lambda_outer"))
    //     snake_parameters.AddDouble("lambda_outer", mParameters["lambda_outer"].GetDouble());
    // if (mParameters.Has("number_of_inner_loops"))
    //     snake_parameters.AddDouble("number_of_inner_loops", mParameters["number_of_inner_loops"].GetInt());
    // if (mParameters.Has("gap_approximation_order"))
    //     snake_parameters.AddInt("gap_approximation_order", mParameters["gap_approximation_order"].GetInt());
    // if (mParameters.Has("gap_relative_tolerance_for_subdivisions"))
    //     snake_parameters.AddDouble("gap_relative_tolerance_for_subdivisions", mParameters["gap_relative_tolerance_for_subdivisions"].GetDouble());
    // if (mParameters.Has("number_of_interpolation_levels"))
    //     snake_parameters.AddInt("number_of_interpolation_levels", mParameters["number_of_interpolation_levels"].GetInt());
    
    // if (mParameters.Has("polynomial_order"))
    //     snake_parameters.AddVector("polynomial_order", mParameters["polynomial_order"].GetVector());
    
    // KRATOS_ERROR << "The NurbsGeometryModelerGapSbm is not yet implemented for 3D. " 
    //     << "Please use the 2D version or implement the 3D version." << std::endl;

    // // TODO: NEXT PR SnakeSbmProcess in 3D
    // // // Create the surrogate_sub_model_part for inner and outer // TODO: extend this in 3D
    // // SnakeSbmProcess snake_sbm_process(*mpModel, snake_parameters);
    // // snake_sbm_process.Execute();

    // // Create the breps for the outer sbm boundary // TODO: extend this in 3D
    // CreateBrepsSbmUtilities<Node, Point> breps_sbm_utilities_3d(mEchoLevel);
    // // TODO: NEXT PR CreateSurrogateBoundary with Volume
    // // breps_sbm_utilities_3d.CreateSurrogateBoundary(mpVolume, r_surrogate_sub_model_part_inner, r_surrogate_sub_model_part_outer, rPointAUvw, rPointBUvw, r_iga_model_part);
}


const Parameters NurbsGeometryModelerGapSbm::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 1,
        "model_part_name" : "IgaModelPart",
        "lower_point_xyz": [0.0, 0.0, 0.0],
        "upper_point_xyz": [1.0, 1.0, 0.0],
        "lower_point_uvw": [0.0, 0.0, 0.0],
        "upper_point_uvw": [1.0, 1.0, 0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [10, 10],
        "number_of_inner_loops": 0,
        "number_initial_points_if_importing_nurbs": 100,
        "number_internal_divisions": 1,
        "gap_relative_tolerance_for_subdivisions": 0.1,
        "number_of_interpolation_levels": 3,
        "gap_approximation_order": 0,
        "gap_element_name": "",
        "gap_interface_condition_name": "",
        "gap_sbm_type": "interpolation",
        "create_surr_outer_from_surr_inner": false,
        "create_surr_inner_from_surr_outer": false,
        "use_for_multipatch": false,
        "use_for_local_refinement": false
    })");
}

const Parameters NurbsGeometryModelerGapSbm::GetValidParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 1,
        "model_part_name" : "IgaModelPart",
        "lower_point_xyz": [0.0, 0.0, 0.0],
        "upper_point_xyz": [1.0, 1.0, 0.0],
        "lower_point_uvw": [0.0, 0.0, 0.0],
        "upper_point_uvw": [1.0, 1.0, 0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [10, 10],
        "lambda_inner": 0.5,
        "lambda_outer": 0.5,
        "number_of_inner_loops": 0,
        "number_initial_points_if_importing_nurbs": 100,
        "number_internal_divisions": 1,
        "gap_relative_tolerance_for_subdivisions": 0.1,
        "number_of_interpolation_levels": 3,
        "gap_approximation_order": 0,
        "skin_model_part_inner_initial_name": "r_skin_model_part_inner_initial",
        "skin_model_part_outer_initial_name": "r_skin_model_part_outer_initial",
        "skin_model_part_name": "r_skin_model_part",
        "gap_element_name": "",
        "gap_interface_condition_name": "",
        "gap_sbm_type": "interpolation",
        "create_surr_outer_from_surr_inner": false,
        "create_surr_inner_from_surr_outer": false,
        "use_for_multipatch": false,
        "use_for_local_refinement": false
    })");
}

} // end namespace kratos
