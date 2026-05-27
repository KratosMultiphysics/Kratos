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
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "custom_modelers/local_refinement_modeler.h"
#include "custom_modelers/iga_modeler_sbm.h"
#include "custom_modelers/nurbs_geometry_modeler_gap_sbm.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_curve_geometry.h"
#include "iga_application_variables.h"
#include "includes/node.h"

namespace Kratos
{

namespace
{

using RectangleType = LocalRefinementModeler::RectangleType;
using IndexType = LocalRefinementModeler::IndexType;
using CurveGeometryType = NurbsCurveGeometry<2, PointerVector<Node>>;
using CurvePointerType = CurveGeometryType::Pointer;

constexpr std::size_t U_MIN = 0;
constexpr std::size_t U_MAX = 1;
constexpr std::size_t V_MIN = 2;
constexpr std::size_t V_MAX = 3;
constexpr double POINT_TOL = 1.0e-8;
constexpr double PARAM_TOL = 1.0e-10;

enum class Orientation
{
    Clockwise,
    CounterClockwise
};

struct Point2D
{
    double X = 0.0;
    double Y = 0.0;
};

struct CurveMetadata
{
    std::string Identifier;
    std::string ConditionName;
    bool HasBrepId = false;
    int BrepId = 0;
    std::string BrepModelPartFullName;
};

CurveMetadata MakeBoundaryConnectorMetadata(const std::string& rIdentifier)
{
    CurveMetadata metadata;
    metadata.ConditionName = "SolidCouplingCondition";
    metadata.Identifier = rIdentifier;
    return metadata;
}

struct SamplePoint
{
    double T = 0.0;
    Point2D Point;
};

struct IntersectionEvent
{
    double T = 0.0;
    Point2D Point;
};

struct RawCurvePiece
{
    bool IsInside = false;
    std::vector<Point2D> Points;
};

struct CurvePiece
{
    CurvePointerType pCurve;
    Point2D StartPoint;
    Point2D EndPoint;
};

struct CurveSplitResult
{
    std::vector<CurvePiece> Pieces;
    bool HasBoundaryIntersection = false;
};

ModelPart& GetOrCreateModelPartByName(Model& rModel, const std::string& rName)
{
    if (rModel.HasModelPart(rName)) {
        return rModel.GetModelPart(rName);
    }

    const auto dot_pos = rName.find('.');
    if (dot_pos == std::string::npos) {
        return rModel.CreateModelPart(rName);
    }

    const std::string root_name = rName.substr(0, dot_pos);
    ModelPart& r_root = rModel.HasModelPart(root_name)
        ? rModel.GetModelPart(root_name)
        : rModel.CreateModelPart(root_name);

    ModelPart* p_current = &r_root;
    std::size_t start = dot_pos + 1;
    while (start < rName.size()) {
        const std::size_t next_dot = rName.find('.', start);
        const std::string token = rName.substr(start, next_dot - start);
        if (!token.empty()) {
            p_current = p_current->HasSubModelPart(token)
                ? &p_current->GetSubModelPart(token)
                : &p_current->CreateSubModelPart(token);
        }
        if (next_dot == std::string::npos) {
            break;
        }
        start = next_dot + 1;
    }

    return *p_current;
}

RectangleType ClipRectangle(const RectangleType& rRectangle, const RectangleType& rBounds)
{
    RectangleType result;
    result[U_MIN] = std::max(rRectangle[U_MIN], rBounds[U_MIN]);
    result[U_MAX] = std::min(rRectangle[U_MAX], rBounds[U_MAX]);
    result[V_MIN] = std::max(rRectangle[V_MIN], rBounds[V_MIN]);
    result[V_MAX] = std::min(rRectangle[V_MAX], rBounds[V_MAX]);
    return result;
}

bool HasPositiveArea(const RectangleType& rRectangle)
{
    return (rRectangle[U_MAX] - rRectangle[U_MIN] > 0.0) &&
           (rRectangle[V_MAX] - rRectangle[V_MIN] > 0.0);
}

std::string RectangleToString(const RectangleType& rRectangle)
{
    std::ostringstream buffer;
    buffer << "[" << rRectangle[U_MIN] << ", " << rRectangle[U_MAX]
           << "] x [" << rRectangle[V_MIN] << ", " << rRectangle[V_MAX] << "]";
    return buffer.str();
}

double Distance(const Point2D& rA, const Point2D& rB)
{
    return std::hypot(rA.X - rB.X, rA.Y - rB.Y);
}

bool IsClose(const Point2D& rA, const Point2D& rB, const double Tolerance = POINT_TOL)
{
    return Distance(rA, rB) <= Tolerance * (1.0 + std::max({std::abs(rA.X), std::abs(rA.Y), std::abs(rB.X), std::abs(rB.Y)}));
}

bool IsInsideOrOn(const Point2D& rPoint, const RectangleType& rRect, const double Tolerance = POINT_TOL)
{
    return rPoint.X >= rRect[U_MIN] - Tolerance &&
           rPoint.X <= rRect[U_MAX] + Tolerance &&
           rPoint.Y >= rRect[V_MIN] - Tolerance &&
           rPoint.Y <= rRect[V_MAX] + Tolerance;
}

bool IsOnBoundary(const Point2D& rPoint, const RectangleType& rRect, const double Tolerance = POINT_TOL)
{
    if (!IsInsideOrOn(rPoint, rRect, Tolerance)) {
        return false;
    }

    return std::abs(rPoint.X - rRect[U_MIN]) <= Tolerance ||
           std::abs(rPoint.X - rRect[U_MAX]) <= Tolerance ||
           std::abs(rPoint.Y - rRect[V_MIN]) <= Tolerance ||
           std::abs(rPoint.Y - rRect[V_MAX]) <= Tolerance;
}

std::string MapSurrogateConditionLayerName(const std::string& rSkinLayerName)
{
    if (rSkinLayerName == "COUPLING_SIDE_OUTER") {
        return "COUPLING_CONDITION_OUTER";
    }
    if (rSkinLayerName == "COUPLING_SIDE_INNER") {
        return "COUPLING_CONDITION_INNER";
    }
    return rSkinLayerName;
}

std::string FindSkinConditionLayerName(
    const ModelPart& rSkinLoopModelPart,
    const IndexType ProjectionNodeId0,
    const IndexType ProjectionNodeId1)
{
    for (const auto& r_condition : rSkinLoopModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.PointsNumber() < 2) {
            continue;
        }

        const IndexType node_id_0 = r_geometry[0].Id();
        const IndexType node_id_1 = r_geometry[1].Id();
        const bool same_orientation =
            node_id_0 == ProjectionNodeId0 && node_id_1 == ProjectionNodeId1;
        const bool opposite_orientation =
            node_id_0 == ProjectionNodeId1 && node_id_1 == ProjectionNodeId0;
        if (!same_orientation && !opposite_orientation) {
            continue;
        }

        return r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
    }

    return "";
}

int ResolveBrepIdFromConditionGeometryOrNodes(const Condition& rCondition)
{
    if (rCondition.Has(BREP_ID)) {
        return rCondition.GetValue(BREP_ID);
    }

    const auto& r_geometry = rCondition.GetGeometry();
    if (r_geometry.Has(BREP_ID)) {
        return r_geometry.GetValue(BREP_ID);
    }

    int resolved_brep_id = 0;
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        const auto& r_node = r_geometry[i];
        if (!r_node.Has(BREP_ID)) {
            continue;
        }

        const int node_brep_id = r_node.GetValue(BREP_ID);
        if (resolved_brep_id == 0) {
            resolved_brep_id = node_brep_id;
        } else {
            KRATOS_ERROR_IF(resolved_brep_id != node_brep_id)
                << "LocalRefinementModeler: conflicting node BREP_ID values on condition #"
                << rCondition.Id() << ": " << resolved_brep_id << " vs " << node_brep_id << std::endl;
        }
    }

    return resolved_brep_id;
}

std::string FindSkinConditionLayerNameByBrepId(
    const ModelPart& rSkinLoopModelPart,
    const int BrepId)
{
    if (BrepId == 0) {
        return "";
    }

    for (const auto& r_condition : rSkinLoopModelPart.Conditions()) {
        const int skin_brep_id = ResolveBrepIdFromConditionGeometryOrNodes(r_condition);
        if (skin_brep_id != BrepId) {
            continue;
        }

        return r_condition.Has(LAYER_NAME) ? r_condition.GetValue(LAYER_NAME) : "";
    }

    return "";
}

std::string FindCommonConnectedLayerName(
    const Node& rNode0,
    const Node& rNode1)
{
    if (!rNode0.Has(CONNECTED_LAYERS) || !rNode1.Has(CONNECTED_LAYERS)) {
        return "";
    }

    const auto& r_layers_0 = rNode0.GetValue(CONNECTED_LAYERS);
    const auto& r_layers_1 = rNode1.GetValue(CONNECTED_LAYERS);
    for (const auto& r_layer_0 : r_layers_0) {
        if (std::find(r_layers_1.begin(), r_layers_1.end(), r_layer_0) != r_layers_1.end()) {
            return r_layer_0;
        }
    }

    return "";
}

void MarkSurrogateConditionLayerNames(
    ModelPart& rPatchModelPart,
    const ModelPart& rPatchSkinModelPart)
{
    auto mark_conditions = [&](ModelPart& rSurrogateModelPart, const ModelPart& rSkinLoopModelPart) {
        for (auto& r_condition : rSurrogateModelPart.Conditions()) {
            const auto& r_geometry = r_condition.GetGeometry();
            if (r_geometry.PointsNumber() < 2) {
                continue;
            }

            const auto& r_node_0 = r_geometry[0];
            const auto& r_node_1 = r_geometry[1];
            std::string skin_layer_name;

            if (r_node_0.Has(PROJECTION_NODE_ID) &&
                r_node_1.Has(PROJECTION_NODE_ID) &&
                rSkinLoopModelPart.HasNode(r_node_0.GetValue(PROJECTION_NODE_ID)) &&
                rSkinLoopModelPart.HasNode(r_node_1.GetValue(PROJECTION_NODE_ID))) {
                const auto& r_projection_node_0 =
                    rSkinLoopModelPart.GetNode(r_node_0.GetValue(PROJECTION_NODE_ID));
                const auto& r_projection_node_1 =
                    rSkinLoopModelPart.GetNode(r_node_1.GetValue(PROJECTION_NODE_ID));
                const double knot_span_size = std::max(
                    norm_2(r_node_0.Coordinates() - r_node_1.Coordinates()),
                    POINT_TOL);

                const bool projections_are_close =
                    norm_2(r_node_0.Coordinates() - r_projection_node_0.Coordinates()) <= 1.0e-1 * knot_span_size &&
                    norm_2(r_node_1.Coordinates() - r_projection_node_1.Coordinates()) <= 1.0e-1 * knot_span_size;

                if (projections_are_close) {
                    skin_layer_name = FindSkinConditionLayerName(
                        rSkinLoopModelPart,
                        r_node_0.GetValue(PROJECTION_NODE_ID),
                        r_node_1.GetValue(PROJECTION_NODE_ID));
                }
            }

            if (skin_layer_name.empty()) {
                skin_layer_name = FindSkinConditionLayerNameByBrepId(
                    rSkinLoopModelPart,
                    ResolveBrepIdFromConditionGeometryOrNodes(r_condition));
            }

            if (skin_layer_name.empty()) {
                skin_layer_name = FindCommonConnectedLayerName(r_node_0, r_node_1);
            }

            if (!skin_layer_name.empty()) {
                r_condition.SetValue(LAYER_NAME, MapSurrogateConditionLayerName(skin_layer_name));
            }
        }
    };

    if (rPatchModelPart.HasSubModelPart("surrogate_outer") &&
        rPatchSkinModelPart.HasSubModelPart("outer")) {
        mark_conditions(
            rPatchModelPart.GetSubModelPart("surrogate_outer"),
            rPatchSkinModelPart.GetSubModelPart("outer"));
    }
    if (rPatchModelPart.HasSubModelPart("surrogate_inner") &&
        rPatchSkinModelPart.HasSubModelPart("inner")) {
        mark_conditions(
            rPatchModelPart.GetSubModelPart("surrogate_inner"),
            rPatchSkinModelPart.GetSubModelPart("inner"));
    }
}

Point2D ClampPointToBoundary(const Point2D& rPoint, const RectangleType& rRect)
{
    Point2D result = rPoint;
    if (std::abs(result.X - rRect[U_MIN]) <= POINT_TOL) {
        result.X = rRect[U_MIN];
    } else if (std::abs(result.X - rRect[U_MAX]) <= POINT_TOL) {
        result.X = rRect[U_MAX];
    }

    if (std::abs(result.Y - rRect[V_MIN]) <= POINT_TOL) {
        result.Y = rRect[V_MIN];
    } else if (std::abs(result.Y - rRect[V_MAX]) <= POINT_TOL) {
        result.Y = rRect[V_MAX];
    }

    return result;
}

Point2D EvaluateCurvePoint(const CurveGeometryType& rCurve, const double T)
{
    CurveGeometryType::CoordinatesArrayType local_coordinates = ZeroVector(3);
    CurveGeometryType::CoordinatesArrayType global_coordinates = ZeroVector(3);
    local_coordinates[0] = T;
    rCurve.GlobalCoordinates(global_coordinates, local_coordinates);

    return Point2D{global_coordinates[0], global_coordinates[1]};
}

CurveMetadata ExtractMetadata(const CurveGeometryType& rCurve)
{
    CurveMetadata metadata;

    if (rCurve.Has(IDENTIFIER)) {
        metadata.Identifier = rCurve.GetValue(IDENTIFIER);
    }
    if (rCurve.Has(CONDITION_NAME)) {
        metadata.ConditionName = rCurve.GetValue(CONDITION_NAME);
    }
    if (rCurve.Has(BREP_ID)) {
        metadata.HasBrepId = true;
        metadata.BrepId = rCurve.GetValue(BREP_ID);
    }
    if (rCurve.Has(BREP_MODEL_PART_FULL_NAME)) {
        metadata.BrepModelPartFullName = rCurve.GetValue(BREP_MODEL_PART_FULL_NAME);
    }

    return metadata;
}

void AssignMetadata(CurveGeometryType& rCurve, const CurveMetadata& rMetadata)
{
    if (!rMetadata.Identifier.empty()) {
        rCurve.SetValue(IDENTIFIER, rMetadata.Identifier);
    }
    if (!rMetadata.ConditionName.empty()) {
        rCurve.SetValue(CONDITION_NAME, rMetadata.ConditionName);
    }
    if (rMetadata.HasBrepId) {
        rCurve.SetValue(BREP_ID, rMetadata.BrepId);
    }
    if (!rMetadata.BrepModelPartFullName.empty()) {
        rCurve.SetValue(BREP_MODEL_PART_FULL_NAME, rMetadata.BrepModelPartFullName);
    }
}

Vector CreateChordLengthKnots(const std::vector<Point2D>& rPoints)
{
    const std::size_t number_of_points = rPoints.size();
    Vector knots(number_of_points);

    if (number_of_points == 0) {
        return knots;
    }

    knots[0] = 0.0;
    if (number_of_points == 1) {
        return knots;
    }

    std::vector<double> cumulative_length(number_of_points, 0.0);
    for (std::size_t i = 1; i < number_of_points; ++i) {
        cumulative_length[i] = cumulative_length[i - 1] + Distance(rPoints[i - 1], rPoints[i]);
    }

    const double total_length = cumulative_length.back();
    if (total_length <= POINT_TOL) {
        for (std::size_t i = 1; i < number_of_points; ++i) {
            knots[i] = static_cast<double>(i) / static_cast<double>(number_of_points - 1);
        }
        return knots;
    }

    for (std::size_t i = 1; i < number_of_points; ++i) {
        knots[i] = cumulative_length[i] / total_length;
    }

    return knots;
}

void DeduplicateConsecutivePoints(std::vector<Point2D>& rPoints)
{
    std::vector<Point2D> unique_points;
    unique_points.reserve(rPoints.size());

    for (const auto& r_point : rPoints) {
        if (unique_points.empty() || !IsClose(unique_points.back(), r_point)) {
            unique_points.push_back(r_point);
        }
    }

    rPoints.swap(unique_points);
}

CurvePointerType CreatePolylineCurve(
    const std::vector<Point2D>& rPoints,
    const CurveMetadata& rMetadata)
{
    KRATOS_ERROR_IF(rPoints.size() < 2)
        << "LocalRefinementModeler: at least two points are needed to create a polyline NURBS curve." << std::endl;

    PointerVector<Node> control_points;
    for (const auto& r_point : rPoints) {
        control_points.push_back(Kratos::make_intrusive<Node>(0, r_point.X, r_point.Y, 0.0));
    }

    Vector knot_vector = CreateChordLengthKnots(rPoints);
    auto p_curve = Kratos::make_shared<CurveGeometryType>(control_points, 1, knot_vector);
    AssignMetadata(*p_curve, rMetadata);
    return p_curve;
}

CurvePointerType CloneCurve(const CurveGeometryType& rCurve)
{
    PointerVector<Node> control_points;
    for (IndexType i = 0; i < rCurve.PointsNumber(); ++i) {
        control_points.push_back(Kratos::make_intrusive<Node>(0, rCurve[i].X(), rCurve[i].Y(), rCurve[i].Z()));
    }

    CurvePointerType p_curve;
    if (rCurve.IsRational()) {
        p_curve = Kratos::make_shared<CurveGeometryType>(
            control_points,
            rCurve.PolynomialDegree(0),
            rCurve.Knots(),
            rCurve.Weights());
    } else {
        p_curve = Kratos::make_shared<CurveGeometryType>(
            control_points,
            rCurve.PolynomialDegree(0),
            rCurve.Knots());
    }

    AssignMetadata(*p_curve, ExtractMetadata(rCurve));
    return p_curve;
}

int GetMinimumDegree(const Vector& rDegrees)
{
    KRATOS_ERROR_IF(rDegrees.size() == 0)
        << "LocalRefinementModeler: polynomial_order must contain at least one entry." << std::endl;

    int minimum_degree = static_cast<int>(rDegrees[0]);
    for (IndexType i = 1; i < rDegrees.size(); ++i) {
        minimum_degree = std::min(minimum_degree, static_cast<int>(rDegrees[i]));
    }
    return minimum_degree;
}

void SetOrAddIntValue(
    Parameters& rParameters,
    const std::string& rName,
    const int Value)
{
    if (rParameters.Has(rName)) {
        rParameters[rName].SetInt(Value);
    } else {
        rParameters.AddEmptyValue(rName).SetInt(Value);
    }
}

void SetOrAddDoubleValue(
    Parameters& rParameters,
    const std::string& rName,
    const double Value)
{
    if (rParameters.Has(rName)) {
        rParameters[rName].SetDouble(Value);
    } else {
        rParameters.AddEmptyValue(rName).SetDouble(Value);
    }
}

void SetOrAddBoolValue(
    Parameters& rParameters,
    const std::string& rName,
    const bool Value)
{
    if (rParameters.Has(rName)) {
        rParameters[rName].SetBool(Value);
    } else {
        rParameters.AddEmptyValue(rName).SetBool(Value);
    }
}

void SetOrAddStringValue(
    Parameters& rParameters,
    const std::string& rName,
    const std::string& rValue)
{
    if (rParameters.Has(rName)) {
        rParameters[rName].SetString(rValue);
    } else {
        rParameters.AddEmptyValue(rName).SetString(rValue);
    }
}

void SetDefaultVectorIfMissing(
    Parameters& rParameters,
    const std::string& rName,
    const Vector& rDefaultValue)
{
    if (!rParameters.Has(rName)) {
        rParameters.AddEmptyValue(rName).SetVector(rDefaultValue);
        return;
    }

    const Vector current_value = rParameters[rName].GetVector();
    if (current_value.size() == 0) {
        rParameters[rName].SetVector(rDefaultValue);
    }
}

Vector GetVectorOrDefault(
    const Parameters& rParameters,
    const std::string& rName,
    const Vector& rDefaultValue)
{
    if (!rParameters.Has(rName)) {
        return rDefaultValue;
    }

    const Vector value = rParameters[rName].GetVector();
    return value.size() == 0 ? rDefaultValue : value;
}

int GetSamplingPointsPerCurve(const Parameters& rParameters)
{
    if (rParameters.Has("geometry_parameters") &&
        rParameters["geometry_parameters"].Has("number_initial_points_if_importing_nurbs")) {
        return rParameters["geometry_parameters"]["number_initial_points_if_importing_nurbs"].GetInt();
    }
    return 200;
}

std::vector<SamplePoint> SampleCurve(
    const CurveGeometryType& rCurve,
    const int RequestedPoints)
{
    std::vector<double> spans;
    rCurve.SpansLocalSpace(spans);

    const int number_of_spans = static_cast<int>(std::max<std::size_t>(1, spans.size() > 1 ? spans.size() - 1 : 1));
    const int points_per_span = std::max(3, static_cast<int>(std::ceil(static_cast<double>(std::max(2, RequestedPoints)) / number_of_spans)));

    std::vector<SamplePoint> samples;
    samples.reserve(static_cast<std::size_t>(number_of_spans * points_per_span));

    for (int i_span = 0; i_span < number_of_spans; ++i_span) {
        const double span_start = spans.empty() ? rCurve.DomainInterval().GetT0() : spans[i_span];
        const double span_end = spans.empty() ? rCurve.DomainInterval().GetT1() : spans[i_span + 1];
        const int number_of_local_points = (i_span == number_of_spans - 1) ? points_per_span : points_per_span - 1;

        for (int i_point = 0; i_point < number_of_local_points; ++i_point) {
            const double local_fraction = static_cast<double>(i_point) / static_cast<double>(points_per_span - 1);
            const double t = span_start + (span_end - span_start) * local_fraction;
            samples.push_back(SamplePoint{t, EvaluateCurvePoint(rCurve, t)});
        }
    }

    if (samples.empty()) {
        const auto interval = rCurve.DomainInterval();
        samples.push_back(SamplePoint{interval.GetT0(), EvaluateCurvePoint(rCurve, interval.GetT0())});
        samples.push_back(SamplePoint{interval.GetT1(), EvaluateCurvePoint(rCurve, interval.GetT1())});
    } else {
        const double t_end = spans.empty() ? rCurve.DomainInterval().GetT1() : spans.back();
        const Point2D end_point = EvaluateCurvePoint(rCurve, t_end);
        if (!IsClose(samples.back().Point, end_point) || std::abs(samples.back().T - t_end) > PARAM_TOL) {
            samples.push_back(SamplePoint{t_end, end_point});
        }
    }

    return samples;
}

std::vector<std::pair<double, Point2D>> ComputeSegmentBoundaryIntersections(
    const Point2D& rPointA,
    const Point2D& rPointB,
    const RectangleType& rRect)
{
    std::vector<std::pair<double, Point2D>> intersections;
    const double dx = rPointB.X - rPointA.X;
    const double dy = rPointB.Y - rPointA.Y;

    auto add_intersection = [&](const double S, const Point2D& rPoint) {
        if (S <= POINT_TOL || S >= 1.0 - POINT_TOL) {
            return;
        }
        if (!IsInsideOrOn(rPoint, rRect)) {
            return;
        }

        const Point2D clamped_point = ClampPointToBoundary(rPoint, rRect);
        for (const auto& r_existing : intersections) {
            if (std::abs(r_existing.first - S) <= PARAM_TOL || IsClose(r_existing.second, clamped_point)) {
                return;
            }
        }

        intersections.emplace_back(S, clamped_point);
    };

    if (std::abs(dx) > POINT_TOL) {
        const double s_u_min = (rRect[U_MIN] - rPointA.X) / dx;
        add_intersection(s_u_min, Point2D{rRect[U_MIN], rPointA.Y + s_u_min * dy});

        const double s_u_max = (rRect[U_MAX] - rPointA.X) / dx;
        add_intersection(s_u_max, Point2D{rRect[U_MAX], rPointA.Y + s_u_max * dy});
    }

    if (std::abs(dy) > POINT_TOL) {
        const double s_v_min = (rRect[V_MIN] - rPointA.Y) / dy;
        add_intersection(s_v_min, Point2D{rPointA.X + s_v_min * dx, rRect[V_MIN]});

        const double s_v_max = (rRect[V_MAX] - rPointA.Y) / dy;
        add_intersection(s_v_max, Point2D{rPointA.X + s_v_max * dx, rRect[V_MAX]});
    }

    std::sort(
        intersections.begin(),
        intersections.end(),
        [](const auto& rLeft, const auto& rRight) { return rLeft.first < rRight.first; });

    return intersections;
}

std::vector<IntersectionEvent> CollectBoundaryIntersections(
    const CurveGeometryType& rCurve,
    const std::vector<SamplePoint>& rSamples,
    const RectangleType& rRect)
{
    std::vector<IntersectionEvent> intersections;

    for (std::size_t i = 0; i + 1 < rSamples.size(); ++i) {
        const auto segment_intersections = ComputeSegmentBoundaryIntersections(
            rSamples[i].Point,
            rSamples[i + 1].Point,
            rRect);

        for (const auto& r_intersection : segment_intersections) {
            const double t = rSamples[i].T + (rSamples[i + 1].T - rSamples[i].T) * r_intersection.first;
            const IntersectionEvent event{t, r_intersection.second};

            if (!intersections.empty()) {
                const auto& r_last = intersections.back();
                if (std::abs(r_last.T - event.T) <= PARAM_TOL || IsClose(r_last.Point, event.Point)) {
                    continue;
                }
            }

            intersections.push_back(event);
        }
    }

    if (!std::is_sorted(intersections.begin(), intersections.end(), [](const auto& rLeft, const auto& rRight) { return rLeft.T < rRight.T; })) {
        std::sort(
            intersections.begin(),
            intersections.end(),
            [](const auto& rLeft, const auto& rRight) { return rLeft.T < rRight.T; });
    }

    std::vector<IntersectionEvent> unique_intersections;
    unique_intersections.reserve(intersections.size());
    for (const auto& r_event : intersections) {
        if (unique_intersections.empty() ||
            (std::abs(unique_intersections.back().T - r_event.T) > PARAM_TOL &&
             !IsClose(unique_intersections.back().Point, r_event.Point))) {
            unique_intersections.push_back(r_event);
        }
    }

    return unique_intersections;
}

CurveSplitResult SplitCurveByRectangle(
    const CurveGeometryType& rCurve,
    const RectangleType& rRect,
    const int NumberOfSamplePoints)
{
    CurveSplitResult result;
    const CurveMetadata metadata = ExtractMetadata(rCurve);
    const auto samples = SampleCurve(rCurve, NumberOfSamplePoints);
    const auto intersections = CollectBoundaryIntersections(rCurve, samples, rRect);
    result.HasBoundaryIntersection = !intersections.empty();

    const Point2D start_point = samples.front().Point;
    const Point2D end_point = samples.back().Point;

    if (intersections.empty()) {
        CurvePiece piece;
        piece.pCurve = CloneCurve(rCurve);
        piece.StartPoint = start_point;
        piece.EndPoint = end_point;
        result.Pieces.push_back(std::move(piece));
        return result;
    }

    std::vector<double> break_parameters;
    std::vector<Point2D> break_points;
    break_parameters.reserve(intersections.size() + 2);
    break_points.reserve(intersections.size() + 2);

    break_parameters.push_back(samples.front().T);
    break_points.push_back(start_point);
    for (const auto& r_intersection : intersections) {
        break_parameters.push_back(r_intersection.T);
        break_points.push_back(r_intersection.Point);
    }
    break_parameters.push_back(samples.back().T);
    break_points.push_back(end_point);

    std::vector<RawCurvePiece> raw_pieces;
    raw_pieces.reserve(break_parameters.size());

    for (std::size_t i_piece = 0; i_piece + 1 < break_parameters.size(); ++i_piece) {
        const double t0 = break_parameters[i_piece];
        const double t1 = break_parameters[i_piece + 1];
        if (t1 - t0 <= PARAM_TOL) {
            continue;
        }

        RawCurvePiece raw_piece;
        raw_piece.Points.push_back(break_points[i_piece]);

        for (const auto& r_sample : samples) {
            if (r_sample.T > t0 + PARAM_TOL && r_sample.T < t1 - PARAM_TOL) {
                raw_piece.Points.push_back(r_sample.Point);
            }
        }

        raw_piece.Points.push_back(break_points[i_piece + 1]);
        DeduplicateConsecutivePoints(raw_piece.Points);
        if (raw_piece.Points.size() < 2 || Distance(raw_piece.Points.front(), raw_piece.Points.back()) <= POINT_TOL) {
            continue;
        }

        const double t_mid = 0.5 * (t0 + t1);
        raw_piece.IsInside = IsInsideOrOn(EvaluateCurvePoint(rCurve, t_mid), rRect);
        raw_pieces.push_back(std::move(raw_piece));
    }

    if (raw_pieces.empty()) {
        return result;
    }

    std::vector<RawCurvePiece> merged_pieces;
    merged_pieces.reserve(raw_pieces.size());
    for (auto& r_piece : raw_pieces) {
        if (!merged_pieces.empty() && merged_pieces.back().IsInside == r_piece.IsInside) {
            auto& r_last = merged_pieces.back();
            r_last.Points.insert(r_last.Points.end(), r_piece.Points.begin() + 1, r_piece.Points.end());
            DeduplicateConsecutivePoints(r_last.Points);
        } else {
            merged_pieces.push_back(std::move(r_piece));
        }
    }

    result.Pieces.reserve(merged_pieces.size());
    for (const auto& r_raw_piece : merged_pieces) {
        CurvePiece piece;
        piece.StartPoint = r_raw_piece.Points.front();
        piece.EndPoint = r_raw_piece.Points.back();

        if (merged_pieces.size() == 1 &&
            IsClose(piece.StartPoint, start_point) &&
            IsClose(piece.EndPoint, end_point)) {
            piece.pCurve = CloneCurve(rCurve);
        } else {
            piece.pCurve = CreatePolylineCurve(r_raw_piece.Points, metadata);
        }

        result.Pieces.push_back(std::move(piece));
    }

    return result;
}

double RectanglePerimeter(const RectangleType& rRect)
{
    return 2.0 * ((rRect[U_MAX] - rRect[U_MIN]) + (rRect[V_MAX] - rRect[V_MIN]));
}

double BoundaryCoordinateCW(const Point2D& rPoint, const RectangleType& rRect)
{
    KRATOS_ERROR_IF_NOT(IsOnBoundary(rPoint, rRect))
        << "LocalRefinementModeler: point (" << rPoint.X << ", " << rPoint.Y
        << ") is not on refinement region boundary " << RectangleToString(rRect) << std::endl;

    const double width = rRect[U_MAX] - rRect[U_MIN];
    const double height = rRect[V_MAX] - rRect[V_MIN];

    if (std::abs(rPoint.Y - rRect[V_MIN]) <= POINT_TOL) {
        return rPoint.X - rRect[U_MIN];
    }
    if (std::abs(rPoint.X - rRect[U_MAX]) <= POINT_TOL) {
        return width + (rPoint.Y - rRect[V_MIN]);
    }
    if (std::abs(rPoint.Y - rRect[V_MAX]) <= POINT_TOL) {
        return width + height + (rRect[U_MAX] - rPoint.X);
    }
    return 2.0 * width + height + (rRect[V_MAX] - rPoint.Y);
}

double BoundaryCoordinateCCW(const Point2D& rPoint, const RectangleType& rRect)
{
    KRATOS_ERROR_IF_NOT(IsOnBoundary(rPoint, rRect))
        << "LocalRefinementModeler: point (" << rPoint.X << ", " << rPoint.Y
        << ") is not on refinement region boundary " << RectangleToString(rRect) << std::endl;

    const double width = rRect[U_MAX] - rRect[U_MIN];
    const double height = rRect[V_MAX] - rRect[V_MIN];

    if (std::abs(rPoint.X - rRect[U_MIN]) <= POINT_TOL) {
        return rPoint.Y - rRect[V_MIN];
    }
    if (std::abs(rPoint.Y - rRect[V_MAX]) <= POINT_TOL) {
        return height + (rPoint.X - rRect[U_MIN]);
    }
    if (std::abs(rPoint.X - rRect[U_MAX]) <= POINT_TOL) {
        return height + width + (rRect[V_MAX] - rPoint.Y);
    }
    return 2.0 * height + width + (rRect[U_MAX] - rPoint.X);
}

Point2D PointFromBoundaryCoordinateCW(double Coordinate, const RectangleType& rRect)
{
    const double width = rRect[U_MAX] - rRect[U_MIN];
    const double height = rRect[V_MAX] - rRect[V_MIN];
    const double perimeter = RectanglePerimeter(rRect);

    Coordinate = std::fmod(Coordinate, perimeter);
    if (Coordinate < 0.0) {
        Coordinate += perimeter;
    }

    if (Coordinate <= width + POINT_TOL) {
        return Point2D{rRect[U_MIN] + Coordinate, rRect[V_MIN]};
    }
    if (Coordinate <= width + height + POINT_TOL) {
        return Point2D{rRect[U_MAX], rRect[V_MIN] + (Coordinate - width)};
    }
    if (Coordinate <= 2.0 * width + height + POINT_TOL) {
        return Point2D{rRect[U_MAX] - (Coordinate - width - height), rRect[V_MAX]};
    }

    return Point2D{rRect[U_MIN], rRect[V_MAX] - (Coordinate - 2.0 * width - height)};
}

Point2D PointFromBoundaryCoordinateCCW(double Coordinate, const RectangleType& rRect)
{
    const double width = rRect[U_MAX] - rRect[U_MIN];
    const double height = rRect[V_MAX] - rRect[V_MIN];
    const double perimeter = RectanglePerimeter(rRect);

    Coordinate = std::fmod(Coordinate, perimeter);
    if (Coordinate < 0.0) {
        Coordinate += perimeter;
    }

    if (Coordinate <= height + POINT_TOL) {
        return Point2D{rRect[U_MIN], rRect[V_MIN] + Coordinate};
    }
    if (Coordinate <= height + width + POINT_TOL) {
        return Point2D{rRect[U_MIN] + (Coordinate - height), rRect[V_MAX]};
    }
    if (Coordinate <= 2.0 * height + width + POINT_TOL) {
        return Point2D{rRect[U_MAX], rRect[V_MAX] - (Coordinate - height - width)};
    }

    return Point2D{rRect[U_MAX] - (Coordinate - 2.0 * height - width), rRect[V_MIN]};
}

std::vector<double> BoundaryBreaks(const RectangleType& rRect, const Orientation BoundaryOrientation)
{
    const double width = rRect[U_MAX] - rRect[U_MIN];
    const double height = rRect[V_MAX] - rRect[V_MIN];

    if (BoundaryOrientation == Orientation::Clockwise) {
        return {width, width + height, 2.0 * width + height, 2.0 * (width + height)};
    }

    return {height, height + width, 2.0 * height + width, 2.0 * (width + height)};
}

double BoundaryCoordinate(
    const Point2D& rPoint,
    const RectangleType& rRect,
    const Orientation BoundaryOrientation)
{
    return BoundaryOrientation == Orientation::Clockwise
        ? BoundaryCoordinateCW(rPoint, rRect)
        : BoundaryCoordinateCCW(rPoint, rRect);
}

Point2D PointFromBoundaryCoordinate(
    const double Coordinate,
    const RectangleType& rRect,
    const Orientation BoundaryOrientation)
{
    return BoundaryOrientation == Orientation::Clockwise
        ? PointFromBoundaryCoordinateCW(Coordinate, rRect)
        : PointFromBoundaryCoordinateCCW(Coordinate, rRect);
}

Orientation OppositeOrientation(const Orientation BoundaryOrientation)
{
    return BoundaryOrientation == Orientation::Clockwise
        ? Orientation::CounterClockwise
        : Orientation::Clockwise;
}

std::vector<CurvePiece> BuildBoundaryConnectorPieces(
    const Point2D& rStartPoint,
    const Point2D& rEndPoint,
    const RectangleType& rRect,
    const Orientation BoundaryOrientation,
    const CurveMetadata& rBoundaryConnectorMetadata)
{
    std::vector<CurvePiece> pieces;
    if (IsClose(rStartPoint, rEndPoint)) {
        return pieces;
    }

    CurveMetadata boundary_metadata = rBoundaryConnectorMetadata;
    boundary_metadata.HasBrepId = false;
    boundary_metadata.BrepId = 0;
    if (boundary_metadata.ConditionName.empty()) {
        boundary_metadata.ConditionName = "SolidCouplingCondition";
    }
    if (boundary_metadata.Identifier.empty()) {
        boundary_metadata.Identifier = "COUPLING_SIDE";
    }

    const double perimeter = RectanglePerimeter(rRect);
    double start_coordinate = BoundaryCoordinate(rStartPoint, rRect, BoundaryOrientation);
    double end_coordinate = BoundaryCoordinate(rEndPoint, rRect, BoundaryOrientation);
    if (end_coordinate < start_coordinate + POINT_TOL) {
        end_coordinate += perimeter;
    }

    std::vector<double> coordinates;
    coordinates.push_back(start_coordinate);
    const auto boundary_breaks = BoundaryBreaks(rRect, BoundaryOrientation);
    const int first_period = static_cast<int>(std::floor(start_coordinate / perimeter));
    const int last_period = static_cast<int>(std::ceil(end_coordinate / perimeter));
    for (int period = first_period; period <= last_period; ++period) {
        const double offset = period * perimeter;
        for (const double break_value : boundary_breaks) {
            const double shifted_break_value = break_value + offset;
            if (shifted_break_value > start_coordinate + POINT_TOL &&
                shifted_break_value < end_coordinate - POINT_TOL) {
                coordinates.push_back(shifted_break_value);
            }
        }
    }
    coordinates.push_back(end_coordinate);
    std::sort(coordinates.begin(), coordinates.end());
    coordinates.erase(
        std::unique(
            coordinates.begin(),
            coordinates.end(),
            [](const double A, const double B) {
                return std::abs(A - B) <= POINT_TOL;
            }),
        coordinates.end());

    for (std::size_t i = 0; i + 1 < coordinates.size(); ++i) {
        const Point2D point_a = PointFromBoundaryCoordinate(coordinates[i], rRect, BoundaryOrientation);
        const Point2D point_b = PointFromBoundaryCoordinate(coordinates[i + 1], rRect, BoundaryOrientation);
        if (IsClose(point_a, point_b)) {
            continue;
        }

        CurvePiece piece;
        piece.StartPoint = point_a;
        piece.EndPoint = point_b;
        piece.pCurve = CreatePolylineCurve({point_a, point_b}, boundary_metadata);
        pieces.push_back(std::move(piece));
    }

    return pieces;
}

Point2D CurveStartPoint(const CurveGeometryType& rCurve)
{
    const auto interval = rCurve.DomainInterval();
    return EvaluateCurvePoint(rCurve, interval.GetT0());
}

Point2D CurveEndPoint(const CurveGeometryType& rCurve)
{
    const auto interval = rCurve.DomainInterval();
    return EvaluateCurvePoint(rCurve, interval.GetT1());
}

std::vector<CurvePointerType> OrderConnectedCurves(const ModelPart& rSourceModelPart)
{
    std::vector<CurvePointerType> curves;
    curves.reserve(rSourceModelPart.NumberOfGeometries());

    for (const auto& r_geometry : rSourceModelPart.Geometries()) {
        auto p_curve = std::dynamic_pointer_cast<CurveGeometryType>(
            rSourceModelPart.pGetGeometry(r_geometry.Id()));
        KRATOS_ERROR_IF_NOT(p_curve)
            << "LocalRefinementModeler: geometry #" << r_geometry.Id()
            << " in model part '" << rSourceModelPart.FullName()
            << "' is not a NurbsCurveGeometry<2, PointerVector<Node>>." << std::endl;
        curves.push_back(p_curve);
    }

    if (curves.empty()) {
        return curves;
    }

    std::vector<Point2D> start_points(curves.size());
    std::vector<Point2D> end_points(curves.size());
    for (std::size_t i = 0; i < curves.size(); ++i) {
        start_points[i] = CurveStartPoint(*curves[i]);
        end_points[i] = CurveEndPoint(*curves[i]);
    }

    std::vector<CurvePointerType> ordered_curves;
    ordered_curves.reserve(curves.size());
    std::vector<bool> used(curves.size(), false);

    ordered_curves.push_back(curves.front());
    used[0] = true;
    Point2D current_end = end_points.front();

    for (std::size_t i = 1; i < curves.size(); ++i) {
        bool found_next = false;
        for (std::size_t j = 1; j < curves.size(); ++j) {
            if (used[j]) {
                continue;
            }
            if (IsClose(current_end, start_points[j])) {
                ordered_curves.push_back(curves[j]);
                used[j] = true;
                current_end = end_points[j];
                found_next = true;
                break;
            }
        }

        KRATOS_ERROR_IF_NOT(found_next)
            << "LocalRefinementModeler: could not reorder curves in model part '"
            << rSourceModelPart.FullName()
            << "'. Curves must already be oriented head-to-tail." << std::endl;
    }

    KRATOS_ERROR_IF_NOT(IsClose(current_end, start_points.front()))
        << "LocalRefinementModeler: the ordered curves in model part '"
        << rSourceModelPart.FullName()
        << "' do not form a closed loop." << std::endl;

    return ordered_curves;
}

void ValidateLoopClosure(const std::vector<CurvePiece>& rPieces, const std::string& rLoopLabel)
{
    if (rPieces.empty()) {
        return;
    }

    for (std::size_t i = 0; i + 1 < rPieces.size(); ++i) {
        KRATOS_ERROR_IF_NOT(IsClose(rPieces[i].EndPoint, rPieces[i + 1].StartPoint))
            << "LocalRefinementModeler: loop '" << rLoopLabel
            << "' has disconnected consecutive curves." << std::endl;
    }

    KRATOS_ERROR_IF_NOT(IsClose(rPieces.back().EndPoint, rPieces.front().StartPoint))
        << "LocalRefinementModeler: loop '" << rLoopLabel
        << "' is not closed after reconstruction." << std::endl;
}

std::vector<CurvePiece> RebuildLoop(
    const std::vector<CurvePiece>& rPieces,
    const RectangleType& rRect,
    const Orientation BetweenPieceOrientation,
    const CurveMetadata& rBoundaryConnectorMetadata,
    const std::string& rLoopLabel)
{
    std::vector<CurvePiece> rebuilt_pieces;
    if (rPieces.empty()) {
        return rebuilt_pieces;
    }

    rebuilt_pieces.reserve(rPieces.size() + 8);
    rebuilt_pieces.push_back(rPieces.front());

    const Point2D first_start_point = rPieces.front().StartPoint;
    Point2D current_end_point = rPieces.front().EndPoint;

    for (std::size_t i = 1; i < rPieces.size(); ++i) {
        KRATOS_ERROR_IF(IsClose(current_end_point, first_start_point))
            << "LocalRefinementModeler: loop '" << rLoopLabel
            << "' closed before consuming all curve pieces." << std::endl;

        const auto& r_next_piece = rPieces[i];
        if (!IsClose(current_end_point, r_next_piece.StartPoint)) {
            KRATOS_ERROR_IF_NOT(IsOnBoundary(current_end_point, rRect) && IsOnBoundary(r_next_piece.StartPoint, rRect))
                << "LocalRefinementModeler: loop '" << rLoopLabel
                << "' needs a boundary connector, but at least one end point is not on the refinement-region boundary."
                << std::endl;

            auto connector_pieces = BuildBoundaryConnectorPieces(
                current_end_point,
                r_next_piece.StartPoint,
                rRect,
                BetweenPieceOrientation,
                rBoundaryConnectorMetadata);

            KRATOS_ERROR_IF(connector_pieces.empty())
                << "LocalRefinementModeler: loop '" << rLoopLabel
                << "' requires a connector, but no refinement-boundary segment could be generated." << std::endl;

            current_end_point = connector_pieces.back().EndPoint;
            rebuilt_pieces.insert(
                rebuilt_pieces.end(),
                connector_pieces.begin(),
                connector_pieces.end());
        }

        KRATOS_ERROR_IF_NOT(IsClose(current_end_point, r_next_piece.StartPoint))
            << "LocalRefinementModeler: loop '" << rLoopLabel
            << "' is still disconnected after adding refinement-boundary connectors." << std::endl;

        rebuilt_pieces.push_back(r_next_piece);
        current_end_point = r_next_piece.EndPoint;
    }

    if (!IsClose(current_end_point, first_start_point)) {
        KRATOS_ERROR_IF_NOT(IsOnBoundary(current_end_point, rRect) && IsOnBoundary(first_start_point, rRect))
            << "LocalRefinementModeler: loop '" << rLoopLabel
            << "' is open and cannot be closed because one end point is not on the refinement-region boundary."
            << std::endl;

        auto closing_pieces = BuildBoundaryConnectorPieces(
            current_end_point,
            first_start_point,
            rRect,
            OppositeOrientation(BetweenPieceOrientation),
            rBoundaryConnectorMetadata);

        KRATOS_ERROR_IF(closing_pieces.empty())
            << "LocalRefinementModeler: loop '" << rLoopLabel
            << "' could not be closed on the refinement-region boundary." << std::endl;

        rebuilt_pieces.insert(
            rebuilt_pieces.end(),
            closing_pieces.begin(),
            closing_pieces.end());
    }

    ValidateLoopClosure(rebuilt_pieces, rLoopLabel);
    return rebuilt_pieces;
}

IndexType NextGeometryId(ModelPart& rModelPart)
{
    IndexType max_id = 0;
    for (const auto& r_geometry : rModelPart.GetRootModelPart().Geometries()) {
        max_id = std::max(max_id, static_cast<IndexType>(r_geometry.Id()));
    }
    return max_id + 1;
}

void FillSubModelPartWithCurves(
    ModelPart& rTargetSubModelPart,
    const std::vector<CurvePiece>& rPieces)
{
    IndexType next_geometry_id = NextGeometryId(rTargetSubModelPart);
    for (const auto& r_piece : rPieces) {
        r_piece.pCurve->SetId(next_geometry_id++);
        rTargetSubModelPart.AddGeometry(r_piece.pCurve);
    }
}

} // namespace

LocalRefinementModeler::LocalRefinementModeler(
    Model& rModel,
    const Parameters ModelParameters)
    : Modeler(rModel, ModelParameters)
    , mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

Modeler::Pointer LocalRefinementModeler::Create(
    Model& rModel,
    const Parameters ModelParameters) const
{
    return Kratos::make_shared<LocalRefinementModeler>(rModel, ModelParameters);
}

void LocalRefinementModeler::SetupModelPart()
{
    KRATOS_ERROR_IF_NOT(mParameters.Has("geometry_parameters"))
        << "LocalRefinementModeler: missing 'geometry_parameters' block." << std::endl;
    KRATOS_ERROR_IF(mParameters["geometry_modeler_type"].GetString() != "gap-sbm")
        << "LocalRefinementModeler: only 'gap-sbm' is supported." << std::endl;

    mEchoLevel = static_cast<SizeType>(mParameters["echo_level"].GetInt());
    GenerateRefinementRegions();

    const std::string root_model_part_name = mParameters["model_part_name"].GetString();
    KRATOS_ERROR_IF(root_model_part_name.empty())
        << "LocalRefinementModeler: 'model_part_name' must be specified." << std::endl;
    GetOrCreateModelPartByName(*mpModel, root_model_part_name);

    const Parameters& r_geometry_parameters = mParameters["geometry_parameters"];
    const bool has_inner = r_geometry_parameters.Has("skin_model_part_inner_initial_name");
    const bool has_outer = r_geometry_parameters.Has("skin_model_part_outer_initial_name");

    KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 0)
        << "SetupModelPart begin | refinement_regions=" << mRefinementRegions.size() << std::endl;

    for (IndexType i_region = 0; i_region < mRefinementRegions.size(); ++i_region) {
        const auto& r_rect = mRefinementRegions[i_region];
        KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 1)
            << "Processing refinement region " << i_region + 1
            << " -> " << RectangleToString(r_rect) << std::endl;

        if (has_inner) {
            const std::string inner_skin_name = r_geometry_parameters["skin_model_part_inner_initial_name"].GetString();
            ProcessSkinModelPart(inner_skin_name, "updated_inner_initial", r_rect, i_region + 1, false);
        }

        if (has_outer) {
            const std::string outer_skin_name = r_geometry_parameters["skin_model_part_outer_initial_name"].GetString();
            ProcessSkinModelPart(outer_skin_name, "updated_outer_initial", r_rect, i_region + 1, true);
        }
    }

    RunGapSbmPatchModelers();

    KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 0)
        << "SetupModelPart end" << std::endl;
}

void LocalRefinementModeler::RunGapSbmPatchModelers()
{
    const Parameters& r_geometry_parameters = mParameters["geometry_parameters"];
    const bool has_inner = r_geometry_parameters.Has("skin_model_part_inner_initial_name");
    const bool has_outer = r_geometry_parameters.Has("skin_model_part_outer_initial_name");

    KRATOS_ERROR_IF(has_inner && has_outer)
        << "LocalRefinementModeler: simultaneous inner/outer initial skin handling is not implemented yet." << std::endl;
    KRATOS_ERROR_IF_NOT(has_inner || has_outer)
        << "LocalRefinementModeler: either 'skin_model_part_inner_initial_name' or 'skin_model_part_outer_initial_name' must be provided." << std::endl;

    const std::string source_skin_model_part_name = has_inner
        ? r_geometry_parameters["skin_model_part_inner_initial_name"].GetString()
        : r_geometry_parameters["skin_model_part_outer_initial_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(source_skin_model_part_name))
        << "LocalRefinementModeler: source skin model part '" << source_skin_model_part_name << "' does not exist." << std::endl;

    const std::string patch_prefix = mParameters["child_patch_prefix"].GetString();

    for (IndexType i_region = 0; i_region < mRefinementRegions.size(); ++i_region) {
        ModelPart& r_source_skin_model_part = mpModel->GetModelPart(source_skin_model_part_name);
        const std::string refinement_sub_model_part_name =
            "Refinement_Region_" + std::to_string(i_region + 1) + "_outer_initial";

        KRATOS_ERROR_IF_NOT(r_source_skin_model_part.HasSubModelPart(refinement_sub_model_part_name))
            << "LocalRefinementModeler: submodelpart '" << refinement_sub_model_part_name
            << "' was not created in '" << r_source_skin_model_part.FullName() << "'." << std::endl;

        const ModelPart& r_refinement_outer_initial =
            r_source_skin_model_part.GetSubModelPart(refinement_sub_model_part_name);
        const std::string patch_suffix = patch_prefix + std::to_string(i_region + 2);
        const std::string patch_skin_model_part_name = BuildPatchSkinModelPartName(patch_suffix);
        const Parameters& r_region_parameters =
            mParameters["refinement_regions"][mRefinementRegionSourceIndices[i_region]];

        RunGapSbmPatch(
            mRefinementRegions[i_region],
            patch_suffix,
            patch_skin_model_part_name,
            "",
            r_refinement_outer_initial.FullName(),
            &r_region_parameters);

        ModelPart& r_refinement_patch_model_part =
            mpModel->GetModelPart(mParameters["model_part_name"].GetString()).GetSubModelPart(patch_suffix);
        ModelPart& r_refinement_patch_skin_model_part =
            mpModel->GetModelPart(patch_skin_model_part_name);
        MarkSurrogateConditionLayerNames(
            r_refinement_patch_model_part,
            r_refinement_patch_skin_model_part);
    }

    std::string base_inner_initial_name;
    std::string base_outer_initial_name;

    ModelPart& r_source_skin_model_part = mpModel->GetModelPart(source_skin_model_part_name);
    if (has_inner) {
        if (r_source_skin_model_part.HasSubModelPart("updated_inner_initial")) {
            base_inner_initial_name = r_source_skin_model_part.GetSubModelPart("updated_inner_initial").FullName();
        } else {
            base_inner_initial_name = r_source_skin_model_part.FullName();
        }
    } else {
        if (r_source_skin_model_part.HasSubModelPart("updated_outer_initial")) {
            base_outer_initial_name = r_source_skin_model_part.GetSubModelPart("updated_outer_initial").FullName();
        } else {
            base_outer_initial_name = r_source_skin_model_part.FullName();
        }
    }

    RunGapSbmPatch(
        mBaseRect,
        patch_prefix + std::to_string(1),
        BuildPatchSkinModelPartName(patch_prefix + std::to_string(1)),
        base_inner_initial_name,
        base_outer_initial_name,
        nullptr);
}

std::string LocalRefinementModeler::BuildPatchSkinModelPartName(const std::string& rPatchSuffix) const
{
    const Parameters& r_geometry_parameters = mParameters["geometry_parameters"];
    std::string skin_root_name = mParameters["model_part_name"].GetString();
    if (r_geometry_parameters.Has("skin_model_part_name")) {
        const std::string candidate_name = r_geometry_parameters["skin_model_part_name"].GetString();
        if (!candidate_name.empty()) {
            skin_root_name = candidate_name;
        }
    }

    return skin_root_name + "_" + rPatchSuffix + "_skin";
}

const Parameters LocalRefinementModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name" : "",
        "base_domain" : {
            "lower_point_uvw" : [0.0, 0.0, 0.0],
            "upper_point_uvw" : [1.0, 1.0, 0.0]
        },
        "refinement_regions" : [],
        "child_patch_prefix" : "Patch",
        "geometry_parameters" : {},
        "analysis_parameters" : {},
        "skin_coupling_model_part_name": "skin_coupling_model_part",
        "coupling_conditions_name": "",
        "geometry_modeler_type": "gap-sbm",
        "echo_level" : 0
    })");
}

void LocalRefinementModeler::GenerateRefinementRegions()
{
    mRefinementRegions.clear();
    mRefinementRegionSourceIndices.clear();

    const auto& base_lower = mParameters["base_domain"]["lower_point_uvw"].GetVector();
    const auto& base_upper = mParameters["base_domain"]["upper_point_uvw"].GetVector();

    mBaseRect = RectangleType{base_lower[0], base_upper[0], base_lower[1], base_upper[1]};
    KRATOS_ERROR_IF_NOT(HasPositiveArea(mBaseRect))
        << "LocalRefinementModeler: base_domain has non-positive area." << std::endl;

    const auto& refinement_regions = mParameters["refinement_regions"];
    for (IndexType i_region = 0; i_region < refinement_regions.size(); ++i_region) {
        const auto& lower = refinement_regions[i_region]["lower_point_uvw"].GetVector();
        const auto& upper = refinement_regions[i_region]["upper_point_uvw"].GetVector();

        RectangleType region{lower[0], upper[0], lower[1], upper[1]};
        RectangleType clipped_region = ClipRectangle(region, mBaseRect);
        if (!HasPositiveArea(clipped_region)) {
            KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 1)
                << "Skipping refinement region " << i_region + 1
                << " because it does not overlap the base domain." << std::endl;
            continue;
        }

        mRefinementRegions.push_back(clipped_region);
        mRefinementRegionSourceIndices.push_back(i_region);
    }
}

ModelPart& LocalRefinementModeler::CreateOrResetSubModelPart(
    ModelPart& rParentModelPart,
    const std::string& rName) const
{
    ModelPart& r_sub_model_part = rParentModelPart.HasSubModelPart(rName)
        ? rParentModelPart.GetSubModelPart(rName)
        : rParentModelPart.CreateSubModelPart(rName);

    r_sub_model_part.Nodes().clear();
    r_sub_model_part.Elements().clear();
    r_sub_model_part.Conditions().clear();
    r_sub_model_part.Geometries().clear();
    return r_sub_model_part;
}

void LocalRefinementModeler::ProcessSkinModelPart(
    const std::string& rSkinModelPartName,
    const std::string& rUpdatedSubModelPartName,
    const RectType& rRect,
    const IndexType RegionIndex,
    const bool ReverseBorderOrientation) const
{
    KRATOS_ERROR_IF(rSkinModelPartName.empty())
        << "LocalRefinementModeler: empty skin model part name provided." << std::endl;
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(rSkinModelPartName))
        << "LocalRefinementModeler: model part '" << rSkinModelPartName << "' does not exist." << std::endl;

    ModelPart& r_skin_model_part = mpModel->GetModelPart(rSkinModelPartName);
    const ModelPart& r_source_model_part =
        r_skin_model_part.HasSubModelPart(rUpdatedSubModelPartName)
            ? r_skin_model_part.GetSubModelPart(rUpdatedSubModelPartName)
            : r_skin_model_part;

    const auto ordered_curves = OrderConnectedCurves(r_source_model_part);
    std::vector<CurvePiece> in_pieces;
    std::vector<CurvePiece> out_pieces;
    const int number_of_sample_points = GetSamplingPointsPerCurve(mParameters);
    bool has_any_boundary_intersection = false;

    for (const auto& p_curve : ordered_curves) {
        const auto split_result = SplitCurveByRectangle(*p_curve, rRect, number_of_sample_points);
        has_any_boundary_intersection = has_any_boundary_intersection || split_result.HasBoundaryIntersection;

        if (split_result.Pieces.empty()) {
            continue;
        }

        for (const auto& r_piece : split_result.Pieces) {
            const Point2D midpoint = EvaluateCurvePoint(*r_piece.pCurve, r_piece.pCurve->DomainInterval().GetT0() +
                0.5 * (r_piece.pCurve->DomainInterval().GetT1() - r_piece.pCurve->DomainInterval().GetT0()));

            if (IsInsideOrOn(midpoint, rRect)) {
                in_pieces.push_back(r_piece);
            } else {
                out_pieces.push_back(r_piece);
            }
        }
    }

    KRATOS_ERROR_IF_NOT(has_any_boundary_intersection)
        << "LocalRefinementModeler: refinement regions whose boundary does not intersect the immersed curve are not implemented yet."
        << std::endl;

    const Orientation in_orientation = ReverseBorderOrientation
        ? Orientation::CounterClockwise
        : Orientation::CounterClockwise;
    const Orientation out_orientation = ReverseBorderOrientation
        ? Orientation::CounterClockwise
        : Orientation::Clockwise;

    const CurveMetadata refinement_boundary_metadata = MakeBoundaryConnectorMetadata(
        ReverseBorderOrientation ? "COUPLING_SIDE_OUTER" : "COUPLING_SIDE_INNER");
    const CurveMetadata base_boundary_metadata = MakeBoundaryConnectorMetadata("COUPLING_SIDE");

    auto rebuilt_in_pieces = RebuildLoop(
        in_pieces,
        rRect,
        in_orientation,
        refinement_boundary_metadata,
        rSkinModelPartName + "::region_" + std::to_string(RegionIndex) + "::in");

    auto rebuilt_out_pieces = RebuildLoop(
        out_pieces,
        rRect,
        out_orientation,
        base_boundary_metadata,
        rSkinModelPartName + "::region_" + std::to_string(RegionIndex) + "::out");

    // KRATOS_WATCH(rebuilt_in_pieces.size())
    // for (IndexType i_piece = 0; i_piece < rebuilt_in_pieces.size(); ++i_piece) {
    //     KRATOS_WATCH(i_piece)
    //     KRATOS_WATCH(rebuilt_in_pieces[i_piece].StartPoint.X)
    //     KRATOS_WATCH(rebuilt_in_pieces[i_piece].StartPoint.Y)
    //     KRATOS_WATCH(rebuilt_in_pieces[i_piece].EndPoint.X)
    //     KRATOS_WATCH(rebuilt_in_pieces[i_piece].EndPoint.Y)
    //     KRATOS_WATCH(rebuilt_in_pieces[i_piece].pCurve->GetValue(IDENTIFIER))
    // }

    // KRATOS_WATCH(rebuilt_out_pieces.size())
    // for (IndexType i_piece = 0; i_piece < rebuilt_out_pieces.size(); ++i_piece) {
    //     KRATOS_WATCH(i_piece)
    //     KRATOS_WATCH(rebuilt_out_pieces[i_piece].StartPoint.X)
    //     KRATOS_WATCH(rebuilt_out_pieces[i_piece].StartPoint.Y)
    //     KRATOS_WATCH(rebuilt_out_pieces[i_piece].EndPoint.X)
    //     KRATOS_WATCH(rebuilt_out_pieces[i_piece].EndPoint.Y)
    //     KRATOS_WATCH(rebuilt_out_pieces[i_piece].pCurve->GetValue(IDENTIFIER))
    // }

    const std::string refinement_sub_model_part_name =
        "Refinement_Region_" + std::to_string(RegionIndex) + "_outer_initial";

    ModelPart& r_refinement_sub_model_part =
        CreateOrResetSubModelPart(r_skin_model_part, refinement_sub_model_part_name);
    ModelPart& r_updated_sub_model_part =
        CreateOrResetSubModelPart(r_skin_model_part, rUpdatedSubModelPartName);

    FillSubModelPartWithCurves(r_refinement_sub_model_part, rebuilt_in_pieces);
    FillSubModelPartWithCurves(r_updated_sub_model_part, rebuilt_out_pieces);

    // if (r_source_model_part.Has(KNOT_SPAN_SIZES)) { //FIXME: probaly useless
    //     r_refinement_sub_model_part.SetValue(KNOT_SPAN_SIZES, r_source_model_part.GetValue(KNOT_SPAN_SIZES));
    //     r_updated_sub_model_part.SetValue(KNOT_SPAN_SIZES, r_source_model_part.GetValue(KNOT_SPAN_SIZES));
    // } else if (r_skin_model_part.Has(KNOT_SPAN_SIZES)) {
    //     r_refinement_sub_model_part.SetValue(KNOT_SPAN_SIZES, r_skin_model_part.GetValue(KNOT_SPAN_SIZES));
    //     r_updated_sub_model_part.SetValue(KNOT_SPAN_SIZES, r_skin_model_part.GetValue(KNOT_SPAN_SIZES));
    // }

    KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 1)
        << "Processed skin '" << rSkinModelPartName
        << "' for refinement region " << RegionIndex
        << " | source_geometries=" << r_source_model_part.NumberOfGeometries()
        << " | in_curves=" << rebuilt_in_pieces.size()
        << " | out_curves=" << rebuilt_out_pieces.size() << std::endl;
}

void LocalRefinementModeler::RunGapSbmPatch(
    const RectType& rRect,
    const std::string& rPatchSuffix,
    const std::string& rPatchSkinModelPartName,
    const std::string& rInnerInitialSkinModelPartName,
    const std::string& rOuterInitialSkinModelPartName,
    const Parameters* pRegionParameters) const
{
    const Parameters& r_geometry_parameters = mParameters["geometry_parameters"];
    const std::string root_model_part_name = mParameters["model_part_name"].GetString();
    ModelPart& r_root_model_part = mpModel->GetModelPart(root_model_part_name);
    ModelPart& r_patch_model_part = CreateOrResetSubModelPart(r_root_model_part, rPatchSuffix);
    const std::string patch_full_name = r_patch_model_part.FullName();

    Parameters patch_geometry = r_geometry_parameters.Clone();
    SetOrAddStringValue(patch_geometry, "model_part_name", patch_full_name);
    SetOrAddStringValue(patch_geometry, "skin_model_part_name", rPatchSkinModelPartName);
    SetOrAddBoolValue(patch_geometry, "use_for_multipatch", false);
    SetOrAddBoolValue(patch_geometry, "use_for_local_refinement", pRegionParameters == nullptr);

    if (patch_geometry.Has("skin_model_part_inner_initial_name")) {
        patch_geometry.RemoveValue("skin_model_part_inner_initial_name");
    }
    if (patch_geometry.Has("skin_model_part_outer_initial_name")) {
        patch_geometry.RemoveValue("skin_model_part_outer_initial_name");
    }

    if (!rInnerInitialSkinModelPartName.empty()) {
        patch_geometry.AddEmptyValue("skin_model_part_inner_initial_name").SetString(rInnerInitialSkinModelPartName);
    } else if (!rOuterInitialSkinModelPartName.empty()) {
        patch_geometry.AddEmptyValue("skin_model_part_outer_initial_name").SetString(rOuterInitialSkinModelPartName);
    } else {
        KRATOS_ERROR << "LocalRefinementModeler: patch '" << patch_full_name
                     << "' has neither an inner nor an outer initial skin model part." << std::endl;
    }

    const Vector base_lower_uvw = mParameters["base_domain"]["lower_point_uvw"].GetVector();
    const Vector base_upper_uvw = mParameters["base_domain"]["upper_point_uvw"].GetVector();
    const Vector base_lower_xyz = GetVectorOrDefault(r_geometry_parameters, "lower_point_xyz", base_lower_uvw);
    const Vector base_upper_xyz = GetVectorOrDefault(r_geometry_parameters, "upper_point_xyz", base_upper_uvw);

    SetDefaultVectorIfMissing(patch_geometry, "lower_point_uvw", base_lower_uvw);
    SetDefaultVectorIfMissing(patch_geometry, "upper_point_uvw", base_upper_uvw);
    SetDefaultVectorIfMissing(patch_geometry, "lower_point_xyz", base_lower_xyz);
    SetDefaultVectorIfMissing(patch_geometry, "upper_point_xyz", base_upper_xyz);

    Vector patch_lower_uvw = patch_geometry["lower_point_uvw"].GetVector();
    Vector patch_upper_uvw = patch_geometry["upper_point_uvw"].GetVector();
    patch_lower_uvw[0] = rRect[U_MIN];
    patch_lower_uvw[1] = rRect[V_MIN];
    patch_upper_uvw[0] = rRect[U_MAX];
    patch_upper_uvw[1] = rRect[V_MAX];
    patch_geometry["lower_point_uvw"].SetVector(patch_lower_uvw);
    patch_geometry["upper_point_uvw"].SetVector(patch_upper_uvw);
    const double base_u_length = base_upper_uvw[0] - base_lower_uvw[0];
    const double base_v_length = base_upper_uvw[1] - base_lower_uvw[1];

    KRATOS_ERROR_IF(base_u_length <= 0.0 || base_v_length <= 0.0)
        << "LocalRefinementModeler: invalid base geometry parameter-space extents." << std::endl;

    const double u_fraction_min = (rRect[U_MIN] - base_lower_uvw[0]) / base_u_length;
    const double u_fraction_max = (rRect[U_MAX] - base_lower_uvw[0]) / base_u_length;
    const double v_fraction_min = (rRect[V_MIN] - base_lower_uvw[1]) / base_v_length;
    const double v_fraction_max = (rRect[V_MAX] - base_lower_uvw[1]) / base_v_length;

    Vector patch_lower_xyz = patch_geometry["lower_point_xyz"].GetVector();
    Vector patch_upper_xyz = patch_geometry["upper_point_xyz"].GetVector();
    patch_lower_xyz[0] = base_lower_xyz[0] + u_fraction_min * (base_upper_xyz[0] - base_lower_xyz[0]);
    patch_upper_xyz[0] = base_lower_xyz[0] + u_fraction_max * (base_upper_xyz[0] - base_lower_xyz[0]);
    patch_lower_xyz[1] = base_lower_xyz[1] + v_fraction_min * (base_upper_xyz[1] - base_lower_xyz[1]);
    patch_upper_xyz[1] = base_lower_xyz[1] + v_fraction_max * (base_upper_xyz[1] - base_lower_xyz[1]);
    patch_lower_xyz[2] = base_lower_xyz[2];
    patch_upper_xyz[2] = base_upper_xyz[2];
    patch_geometry["lower_point_xyz"].SetVector(patch_lower_xyz);
    patch_geometry["upper_point_xyz"].SetVector(patch_upper_xyz);

    Vector patch_spans = patch_geometry["number_of_knot_spans"].GetVector();
    if (pRegionParameters != nullptr && pRegionParameters->Has("number_of_knot_spans")) {
        patch_spans = (*pRegionParameters)["number_of_knot_spans"].GetVector();
    } else if (pRegionParameters != nullptr) {
        const Vector base_spans = r_geometry_parameters["number_of_knot_spans"].GetVector();
        patch_spans[0] = std::max(
            1,
            static_cast<int>(std::round((rRect[U_MAX] - rRect[U_MIN]) / (mBaseRect[U_MAX] - mBaseRect[U_MIN]) * base_spans[0])));
        patch_spans[1] = std::max(
            1,
            static_cast<int>(std::round((rRect[V_MAX] - rRect[V_MIN]) / (mBaseRect[V_MAX] - mBaseRect[V_MIN]) * base_spans[1])));
    }
    patch_geometry["number_of_knot_spans"].SetVector(patch_spans);

    if (pRegionParameters != nullptr && pRegionParameters->Has("polynomial_order")) {
        patch_geometry["polynomial_order"].SetVector((*pRegionParameters)["polynomial_order"].GetVector());
    }

    const int gap_approximation_order =
        (pRegionParameters != nullptr && pRegionParameters->Has("gap_approximation_order"))
            ? (*pRegionParameters)["gap_approximation_order"].GetInt()
            : GetMinimumDegree(patch_geometry["polynomial_order"].GetVector());
    SetOrAddIntValue(patch_geometry, "gap_approximation_order", gap_approximation_order);

    int number_internal_divisions = gap_approximation_order;
    if (pRegionParameters != nullptr) {
        if (pRegionParameters->Has("number_internal_divisions")) {
            number_internal_divisions = (*pRegionParameters)["number_internal_divisions"].GetInt();
        }
    } else if (patch_geometry.Has("number_internal_divisions")) {
        number_internal_divisions = patch_geometry["number_internal_divisions"].GetInt();
    }
    SetOrAddIntValue(patch_geometry, "number_internal_divisions", number_internal_divisions);

    // KRATOS_WATCH(number_internal_divisions)
    // KRATOS_WATCH("******************************************************************************")

    if (pRegionParameters != nullptr && pRegionParameters->Has("gap_relative_tolerance_for_subdivisions")) {
        SetOrAddDoubleValue(
            patch_geometry,
            "gap_relative_tolerance_for_subdivisions",
            (*pRegionParameters)["gap_relative_tolerance_for_subdivisions"].GetDouble());
    }

    Vector patch_span_sizes(2);
    patch_span_sizes[0] = (rRect[U_MAX] - rRect[U_MIN]) / patch_spans[0];
    patch_span_sizes[1] = (rRect[V_MAX] - rRect[V_MIN]) / patch_spans[1];
    r_patch_model_part.SetValue(KNOT_SPAN_SIZES, patch_span_sizes);

    std::vector<Vector> patch_parameter_corners(2);
    patch_parameter_corners[0].resize(2);
    patch_parameter_corners[1].resize(2);
    patch_parameter_corners[0][0] = rRect[U_MIN];
    patch_parameter_corners[0][1] = rRect[U_MAX];
    patch_parameter_corners[1][0] = rRect[V_MIN];
    patch_parameter_corners[1][1] = rRect[V_MAX];
    r_patch_model_part.SetValue(PARAMETER_SPACE_CORNERS, patch_parameter_corners);

    Matrix patch_parameter_corners_matrix(2, 2, 0.0);
    patch_parameter_corners_matrix(0, 0) = rRect[U_MIN];
    patch_parameter_corners_matrix(0, 1) = rRect[U_MAX];
    patch_parameter_corners_matrix(1, 0) = rRect[V_MIN];
    patch_parameter_corners_matrix(1, 1) = rRect[V_MAX];
    r_patch_model_part.SetValue(PATCH_PARAMETER_SPACE_CORNERS, patch_parameter_corners_matrix);

    KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 1)
        << "Running NurbsGeometryModelerGapSbm for patch '" << patch_full_name
        << "' | skin='" << rPatchSkinModelPartName << "'" << std::endl;

    NurbsGeometryModelerGapSbm geometry_modeler(*mpModel, patch_geometry);
    geometry_modeler.SetupGeometryModel();
    geometry_modeler.PrepareGeometryModel();
    geometry_modeler.SetupModelPart();

    KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_parameters"))
        << "LocalRefinementModeler: missing 'analysis_parameters' block." << std::endl;

    Parameters patch_analysis = mParameters["analysis_parameters"].Clone();
    KRATOS_ERROR_IF_NOT(patch_analysis.Has("analysis_model_part_name"))
        << "LocalRefinementModeler: analysis_parameters must define 'analysis_model_part_name'." << std::endl;
    patch_analysis["analysis_model_part_name"].SetString(patch_full_name);

    if (patch_analysis.Has("element_condition_list") && patch_analysis["element_condition_list"].IsArray()) {
        Parameters condition_list = patch_analysis["element_condition_list"];
        Parameters reduced_elements("[]");

        for (IndexType i = 0; i < condition_list.size(); ++i) {
            const Parameters e = condition_list[i];
            if (e.Has("type") && e["type"].IsString() && e["type"].GetString() == "element") {
                Parameters elem = e.Clone();
                if (!elem.Has("geometry_type")) {
                    elem.AddEmptyValue("geometry_type").SetString("GeometrySurface");
                }
                reduced_elements.Append(elem);
            }
        }

        for (IndexType i = 0; i < condition_list.size(); ++i) {
            Parameters c = condition_list[i].Clone();
            if (!(c.Has("type") && c["type"].IsString() && c["type"].GetString() == "condition")) {
                continue;
            }
            if (!c.Has("sbm_parameters")) {
                continue;
            }
            if (!c.Has("geometry_type")) {
                c.AddEmptyValue("geometry_type").SetString("SurfaceEdge");
            }
            reduced_elements.Append(c);
        }

        patch_analysis["element_condition_list"] = reduced_elements;
    }

    if (!patch_analysis.Has("skin_model_part_name")) {
        patch_analysis.AddEmptyValue("skin_model_part_name");
    }
    patch_analysis["skin_model_part_name"].SetString(rPatchSkinModelPartName);

    KRATOS_INFO_IF("LocalRefinementModeler", mEchoLevel > 1)
        << "Running IgaModelerSbm for patch '" << patch_full_name
        << "' | analysis_skin='" << rPatchSkinModelPartName << "'" << std::endl;

    IgaModelerSbm iga_modeler(*mpModel, patch_analysis);
    iga_modeler.SetupModelPart();
}

} // namespace Kratos
