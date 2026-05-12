//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

// Project includes
#include "custom_processes/iga_contact_process_gap_sbm_mortar.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "containers/pointer_vector.h"
#include "custom_utilities/iga_sbm_utilities.h"
#include "geometries/brep_curve.h"
#include "geometries/nurbs_curve_geometry.h"
#include "iga_application_variables.h"
#include "includes/define.h"
#include "integration/integration_point_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/intersection_utilities.h"
#include "utilities/entities_utilities.h"

namespace Kratos
{

namespace
{
using NodeType = Node;
using GeometryType = Geometry<NodeType>;
using GeometryPointerType = GeometryType::Pointer;
using CoordinatesArrayType = GeometryType::CoordinatesArrayType;
using IntegrationPointsArrayType = GeometryType::IntegrationPointsArrayType;
using NurbsCurveGeometryType = NurbsCurveGeometry<3, PointerVector<NodeType>>;

using PointType = NodeType;
using PointTypePointer = NodeType::Pointer;
using PointVector = std::vector<PointType::Pointer>;
using PointIterator = std::vector<PointType::Pointer>::iterator;
using DistanceVector = std::vector<double>;
using DistanceIterator = std::vector<double>::iterator;
using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;

constexpr double parameter_tolerance = 1.0e-8;
constexpr double normal_tolerance = 1.0e-12;

struct CurveData
{
    GeometryPointerType pReferenceCurve = nullptr;
    GeometryPointerType pDeformedCurve = nullptr;
    std::vector<Geometry<Node>::Pointer> ReferenceGeometries;
};

struct OverlapData
{
    bool IsValid = false;

    double MasterDeformedBegin = 0.0;
    double MasterDeformedEnd = 0.0;
    double MasterReferenceBegin = 0.0;
    double MasterReferenceEnd = 0.0;

    double SlaveDeformedBegin = 0.0;
    double SlaveDeformedEnd = 0.0;
    double SlaveReferenceBegin = 0.0;
    double SlaveReferenceEnd = 0.0;

    CoordinatesArrayType MasterReferencePointBegin = ZeroVector(3);
    CoordinatesArrayType MasterReferencePointEnd = ZeroVector(3);
    CoordinatesArrayType SlaveReferencePointBegin = ZeroVector(3);
    CoordinatesArrayType SlaveReferencePointEnd = ZeroVector(3);
};

struct CandidateOverlapData
{
    const CurveData* pSlaveCurve = nullptr;
    OverlapData Overlap;
};

using CandidatePointData = std::pair<double, PointType::Pointer>;

struct EndpointProjectionData
{
    bool IsValid = false;
    const CurveData* pMasterCurve = nullptr;
    double MasterDeformedParameter = 0.0;
    double MasterReferenceParameter = 0.0;
    CoordinatesArrayType MasterDeformedPoint = ZeroVector(3);
    CoordinatesArrayType MasterReferencePoint = ZeroVector(3);
    double DistanceSquared = std::numeric_limits<double>::max();
};

struct SlaveProjectionData
{
    bool IsValid = false;
    const CurveData* pSlaveCurve = nullptr;
    double SlaveDeformedParameter = 0.0;
    double SlaveReferenceParameter = 0.0;
    CoordinatesArrayType SlaveDeformedPoint = ZeroVector(3);
    CoordinatesArrayType SlaveReferencePoint = ZeroVector(3);
    array_1d<double, 3> SlaveNormal = ZeroVector(3);
    double DistanceSquared = std::numeric_limits<double>::max();
};

struct CouplingEntity
{
    const CurveData* pSlaveCurve = nullptr;
    const CurveData* pMasterCurve = nullptr;
    double MasterDeformedBegin = 0.0;
    double MasterDeformedEnd = 0.0;
    double MasterReferenceBegin = 0.0;
    double MasterReferenceEnd = 0.0;
};

bool HaveCompatibleNormals(
    const GeometryType& rMasterDeformedCurve,
    const GeometryType& rSlaveDeformedCurve);

bool GetCurveDomain(
    const GeometryType& rCurveGeometry,
    double& rDomainBegin,
    double& rDomainEnd)
{
    Vector domain_interval;
    rCurveGeometry.DomainInterval(domain_interval);
    if (domain_interval.size() >= 2) {
        rDomainBegin = domain_interval[0];
        rDomainEnd = domain_interval[1];
        return true;
    }

    std::vector<double> spans;
    rCurveGeometry.SpansLocalSpace(spans);
    if (spans.size() >= 2) {
        rDomainBegin = spans.front();
        rDomainEnd = spans.back();
        return true;
    }

    return false;
}

double ClampToDomain(
    const GeometryType& rCurveGeometry,
    const double Parameter)
{
    double domain_begin = 0.0;
    double domain_end = 0.0;
    KRATOS_ERROR_IF_NOT(GetCurveDomain(rCurveGeometry, domain_begin, domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve curve domain." << std::endl;

    const double min_parameter = std::min(domain_begin, domain_end);
    const double max_parameter = std::max(domain_begin, domain_end);
    return std::max(min_parameter, std::min(max_parameter, Parameter));
}

CoordinatesArrayType MakeLocalCoordinate(const double Parameter)
{
    CoordinatesArrayType local_coordinates = ZeroVector(3);
    local_coordinates[0] = Parameter;
    return local_coordinates;
}

CoordinatesArrayType EvaluateCurvePoint(
    const GeometryType& rCurveGeometry,
    const double Parameter)
{
    CoordinatesArrayType global_coordinates = ZeroVector(3);
    const auto local_coordinates = MakeLocalCoordinate(Parameter);
    rCurveGeometry.GlobalCoordinates(global_coordinates, local_coordinates);
    return global_coordinates;
}

bool ProjectPointOntoCurve(
    const GeometryType& rCurveGeometry,
    const CoordinatesArrayType& rPointGlobalCoordinates,
    const double InitialGuess,
    double& rProjectedParameter,
    CoordinatesArrayType& rProjectedPoint,
    double& rDistanceSquared,
    int& rIsInside)
{
    CoordinatesArrayType local_coordinates = MakeLocalCoordinate(InitialGuess);
    rCurveGeometry.ProjectionPointGlobalToLocalSpace(rPointGlobalCoordinates, local_coordinates);
    rProjectedParameter = ClampToDomain(rCurveGeometry, local_coordinates[0]);
    local_coordinates[0] = rProjectedParameter;

    rCurveGeometry.GlobalCoordinates(rProjectedPoint, local_coordinates);

    const auto diff = rProjectedPoint - rPointGlobalCoordinates;
    rDistanceSquared = inner_prod(diff, diff);

    rIsInside = rCurveGeometry.IsInsideLocalSpace(local_coordinates, parameter_tolerance);
    return rIsInside != 0;
}

void MapDeformedPointToReferenceCurve(
    const GeometryType& rReferenceCurve,
    const std::vector<Geometry<Node>::Pointer>& rReferenceGeometries,
    const double DeformedParameter,
    const CoordinatesArrayType& rDeformedPoint,
    double& rReferenceParameter,
    CoordinatesArrayType& rReferencePoint,
    double& rConsistencyDistanceSquared)
{
    KRATOS_ERROR_IF(rReferenceGeometries.empty())
        << "::[IgaContactProcessGapSbmMortar]:: Missing surrogate reference geometries for parameter mapping."
        << std::endl;

    rReferenceParameter = ClampToDomain(rReferenceCurve, DeformedParameter);
    rReferencePoint = EvaluateCurvePoint(rReferenceCurve, rReferenceParameter);

    auto p_reference_node = Kratos::make_intrusive<NodeType>(
        0,
        rReferencePoint[0],
        rReferencePoint[1],
        rReferencePoint[2]);
    p_reference_node->SetValue(NEIGHBOUR_GEOMETRIES, rReferenceGeometries);

    CoordinatesArrayType reference_point_deformed = rReferencePoint;
    IgaSbmUtilities::GetDeformedPosition(*p_reference_node, reference_point_deformed);

    const auto difference = reference_point_deformed - rDeformedPoint;
    rConsistencyDistanceSquared = inner_prod(difference, difference);
}

CoordinatesArrayType EvaluateDeformedReferenceCurvePoint(
    const GeometryType& rReferenceCurve,
    const std::vector<Geometry<Node>::Pointer>& rReferenceGeometries,
    const double ReferenceParameter)
{
    KRATOS_ERROR_IF(rReferenceGeometries.empty())
        << "::[IgaContactProcessGapSbmMortar]:: Missing surrogate reference geometries for deformed evaluation."
        << std::endl;

    CoordinatesArrayType reference_point = EvaluateCurvePoint(rReferenceCurve, ReferenceParameter);
    auto p_reference_node = Kratos::make_intrusive<NodeType>(
        0,
        reference_point[0],
        reference_point[1],
        reference_point[2]);
    p_reference_node->SetValue(NEIGHBOUR_GEOMETRIES, rReferenceGeometries);

    CoordinatesArrayType reference_point_deformed = reference_point;
    IgaSbmUtilities::GetDeformedPosition(*p_reference_node, reference_point_deformed);
    return reference_point_deformed;
}

std::vector<double> BuildClippedCurveSpans(
    const GeometryType& rCurveGeometry,
    const double IntervalBegin,
    const double IntervalEnd);

bool ProjectPointOntoCurveAlongNormal(
    const CoordinatesArrayType& rSlavePointDeformed,
    const array_1d<double, 3>& rSlaveNormal,
    const CurveData& rMasterCurve,
    const double ProjectionDistanceLimitSquared,
    EndpointProjectionData& rBestProjection)
{
    rBestProjection = EndpointProjectionData();

    if (norm_2(rSlaveNormal) <= normal_tolerance) {
        return false;
    }

    double master_domain_begin = 0.0;
    double master_domain_end = 0.0;
    KRATOS_ERROR_IF_NOT(GetCurveDomain(*rMasterCurve.pDeformedCurve, master_domain_begin, master_domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve master deformed domain."
        << std::endl;

    const auto master_spans = BuildClippedCurveSpans(
        *rMasterCurve.pDeformedCurve,
        master_domain_begin,
        master_domain_end);

    if (master_spans.size() < 2) {
        return false;
    }

    CoordinatesArrayType ray_direction = rSlaveNormal;
    const auto master_center = rMasterCurve.pDeformedCurve->Center().Coordinates();
    if (inner_prod(master_center - rSlavePointDeformed, ray_direction) < 0.0) {
        ray_direction *= -1.0;
    }

    const double ray_length = norm_2(master_center - rSlavePointDeformed) + 2.0 * std::sqrt(ProjectionDistanceLimitSquared) + 1.0e-8;
    const CoordinatesArrayType ray_end = rSlavePointDeformed + ray_direction * ray_length;

    const auto try_accept_intersection = [
        &rSlavePointDeformed,
        &ray_direction,
        ProjectionDistanceLimitSquared,
        &rMasterCurve,
        &rBestProjection](
            const CoordinatesArrayType& rIntersectionPoint,
            const double MasterDeformedParameter,
            const double DistanceSquared) {
        if (DistanceSquared >= rBestProjection.DistanceSquared) {
            return;
        }

        const auto ray_vector = rIntersectionPoint - rSlavePointDeformed;
        const double distance_along_ray = inner_prod(ray_vector, ray_direction);
        if (distance_along_ray < -parameter_tolerance) {
            return;
        }

        double master_reference_parameter = 0.0;
        CoordinatesArrayType master_reference_point = ZeroVector(3);
        double consistency_distance_sq = 0.0;
        MapDeformedPointToReferenceCurve(
            *rMasterCurve.pReferenceCurve,
            rMasterCurve.ReferenceGeometries,
            MasterDeformedParameter,
            rIntersectionPoint,
            master_reference_parameter,
            master_reference_point,
            consistency_distance_sq);

        if (consistency_distance_sq > ProjectionDistanceLimitSquared) {
            return;
        }

        rBestProjection.IsValid = true;
        rBestProjection.pMasterCurve = &rMasterCurve;
        rBestProjection.MasterDeformedParameter = MasterDeformedParameter;
        rBestProjection.MasterReferenceParameter = master_reference_parameter;
        rBestProjection.MasterDeformedPoint = rIntersectionPoint;
        rBestProjection.MasterReferencePoint = master_reference_point;
        rBestProjection.DistanceSquared = DistanceSquared;
    };

    for (std::size_t i_span = 0; i_span + 1 < master_spans.size(); ++i_span) {
        const double span_begin = master_spans[i_span];
        const double span_end = master_spans[i_span + 1];

        const CoordinatesArrayType segment_begin = EvaluateCurvePoint(*rMasterCurve.pDeformedCurve, span_begin);
        const CoordinatesArrayType segment_end = EvaluateCurvePoint(*rMasterCurve.pDeformedCurve, span_end);

        CoordinatesArrayType intersection_point = ZeroVector(3);
        const int intersection_type = IntersectionUtilities::ComputeLineLineIntersection(
            rSlavePointDeformed,
            ray_end,
            segment_begin,
            segment_end,
            intersection_point);

        if (intersection_type == 0) {
            continue;
        }

        if (intersection_type == 2) {
            try_accept_intersection(
                segment_begin,
                span_begin,
                inner_prod(segment_begin - rSlavePointDeformed, segment_begin - rSlavePointDeformed));
            try_accept_intersection(
                segment_end,
                span_end,
                inner_prod(segment_end - rSlavePointDeformed, segment_end - rSlavePointDeformed));
            continue;
        }

        const auto segment_vector = segment_end - segment_begin;
        const double segment_length_sq = inner_prod(segment_vector, segment_vector);
        if (segment_length_sq <= normal_tolerance) {
            continue;
        }

        double segment_alpha = inner_prod(intersection_point - segment_begin, segment_vector) / segment_length_sq;
        segment_alpha = std::max(0.0, std::min(1.0, segment_alpha));

        const double master_deformed_parameter = span_begin + segment_alpha * (span_end - span_begin);
        const double distance_squared = inner_prod(
            intersection_point - rSlavePointDeformed,
            intersection_point - rSlavePointDeformed);

        try_accept_intersection(intersection_point, master_deformed_parameter, distance_squared);
    }

    return rBestProjection.IsValid;
}

bool ComputeCurveNormal(
    const GeometryType& rCurveGeometry,
    const double Parameter,
    array_1d<double, 3>& rNormal)
{
    Matrix jacobian = ZeroMatrix(3, 1);
    const auto local_coordinates = MakeLocalCoordinate(Parameter);

    if (rCurveGeometry.HasGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)) {
        const auto p_background_curve =
            rCurveGeometry.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        p_background_curve->Jacobian(jacobian, local_coordinates);
    } else {
        rCurveGeometry.Jacobian(jacobian, local_coordinates);
    }

    const double tangent_x = jacobian(0, 0);
    const double tangent_y = jacobian(1, 0);
    const double tangent_norm = std::sqrt(tangent_x * tangent_x + tangent_y * tangent_y);
    if (tangent_norm <= normal_tolerance) {
        return false;
    }

    rNormal = ZeroVector(3);
    rNormal[0] = tangent_y / tangent_norm;
    rNormal[1] = -tangent_x / tangent_norm;
    return true;
}

double ComputeMaxCenterSearchRadius(
    const PointVector& rCenterPoints,
    const double MinimumRadius)
{
    if (rCenterPoints.empty()) {
        return MinimumRadius;
    }

    array_1d<double, 3> min_coords;
    array_1d<double, 3> max_coords;
    for (IndexType i = 0; i < 3; ++i) {
        min_coords[i] = std::numeric_limits<double>::max();
        max_coords[i] = -std::numeric_limits<double>::max();
    }

    for (const auto& p_center_point : rCenterPoints) {
        for (IndexType i = 0; i < 3; ++i) {
            min_coords[i] = std::min(min_coords[i], p_center_point->Coordinates()[i]);
            max_coords[i] = std::max(max_coords[i], p_center_point->Coordinates()[i]);
        }
    }

    const auto diagonal = max_coords - min_coords;
    return std::max(
        std::sqrt(inner_prod(diagonal, diagonal)) + std::numeric_limits<double>::epsilon(),
        MinimumRadius);
}

void SearchNearestCurveCenterCandidates(
    const CoordinatesArrayType& rQueryPoint,
    DynamicBins& rBins,
    const std::size_t MaxConsideredNeighbours,
    const double InitialSearchRadius,
    const double MaxSearchRadius,
    PointVector& rBinsResults,
    DistanceVector& rBinsDistances,
    std::vector<CandidatePointData>& rCandidatePoints)
{
    PointType query_point(
        0,
        rQueryPoint[0],
        rQueryPoint[1],
        rQueryPoint[2]);

    rCandidatePoints.clear();
    if (rBinsResults.empty()) {
        return;
    }

    if (MaxConsideredNeighbours == 1) {
        DynamicBins::PointerType p_nearest = nullptr;
        double nearest_distance = std::numeric_limits<double>::max();
        rBins.SearchNearestPoint(query_point, p_nearest, nearest_distance);
        if (p_nearest) {
            rCandidatePoints.emplace_back(nearest_distance, p_nearest);
        }
        return;
    }

    std::size_t obtained_results = 0;
    double search_radius = InitialSearchRadius;
    while (obtained_results < MaxConsideredNeighbours && search_radius <= MaxSearchRadius) {
        obtained_results = rBins.SearchInRadius(
            query_point,
            search_radius,
            rBinsResults.begin(),
            rBinsDistances.begin(),
            rBinsResults.size());

        if (obtained_results >= MaxConsideredNeighbours || search_radius >= MaxSearchRadius) {
            break;
        }

        search_radius = std::min(search_radius * 2.0, MaxSearchRadius);
    }

    for (std::size_t i = 0; i < obtained_results; ++i) {
        rCandidatePoints.emplace_back(rBinsDistances[i], rBinsResults[i]);
    }

    if (rCandidatePoints.size() > MaxConsideredNeighbours) {
        auto nth = rCandidatePoints.begin() + MaxConsideredNeighbours;
        std::nth_element(
            rCandidatePoints.begin(),
            nth,
            rCandidatePoints.end(),
            [](const auto& rA, const auto& rB) {
                return rA.first < rB.first;
            });
        rCandidatePoints.resize(MaxConsideredNeighbours);
    }

    std::sort(
        rCandidatePoints.begin(),
        rCandidatePoints.end(),
        [](const auto& rA, const auto& rB) {
            return rA.first < rB.first;
        });
}

bool FindBestMasterProjectionForSlavePoint(
    const CoordinatesArrayType& rSlavePointDeformed,
    const array_1d<double, 3>& rSlaveNormal,
    const CurveData& rSlaveCurve,
    const std::vector<CurveData>& rMasterCurveData,
    DynamicBins& rMasterBins,
    const std::size_t MaxConsideredNeighbours,
    const double InitialSearchRadius,
    const double MaxSearchRadius,
    const double ProjectionDistanceLimitSquared,
    PointVector& rBinsResults,
    DistanceVector& rBinsDistances,
    std::vector<CandidatePointData>& rCandidatePoints,
    EndpointProjectionData& rBestProjection)
{
    rBestProjection = EndpointProjectionData();

    SearchNearestCurveCenterCandidates(
        rSlavePointDeformed,
        rMasterBins,
        MaxConsideredNeighbours,
        InitialSearchRadius,
        MaxSearchRadius,
        rBinsResults,
        rBinsDistances,
        rCandidatePoints);

    for (const auto& r_candidate : rCandidatePoints) {
        const IndexType candidate_id = r_candidate.second->Id();
        KRATOS_ERROR_IF(candidate_id == 0 || candidate_id > rMasterCurveData.size())
            << "::[IgaContactProcessGapSbmMortar]:: invalid master candidate id "
            << candidate_id << std::endl;

        const auto& r_master_curve = rMasterCurveData[candidate_id - 1];
        if (!HaveCompatibleNormals(*r_master_curve.pDeformedCurve, *rSlaveCurve.pDeformedCurve)) {
            continue;
        }

        EndpointProjectionData candidate_projection;
        if (!ProjectPointOntoCurveAlongNormal(
                rSlavePointDeformed,
                rSlaveNormal,
                r_master_curve,
                ProjectionDistanceLimitSquared,
                candidate_projection)) {
            continue;
        }

        if (candidate_projection.DistanceSquared >= rBestProjection.DistanceSquared) {
            continue;
        }

        rBestProjection = candidate_projection;
    }

    return rBestProjection.IsValid;
}

bool FindBestSlaveProjectionForMasterPoint(
    const CoordinatesArrayType& rMasterPointDeformed,
    const CurveData& rMasterCurve,
    const std::vector<CurveData>& rSlaveCurveData,
    DynamicBins& rSlaveBins,
    const std::size_t MaxConsideredNeighbours,
    const double InitialSearchRadius,
    const double MaxSearchRadius,
    const double ProjectionDistanceLimitSquared,
    PointVector& rBinsResults,
    DistanceVector& rBinsDistances,
    std::vector<CandidatePointData>& rCandidatePoints,
    SlaveProjectionData& rBestProjection)
{
    rBestProjection = SlaveProjectionData();

    SearchNearestCurveCenterCandidates(
        rMasterPointDeformed,
        rSlaveBins,
        MaxConsideredNeighbours,
        InitialSearchRadius,
        MaxSearchRadius,
        rBinsResults,
        rBinsDistances,
        rCandidatePoints);

    for (const auto& r_candidate : rCandidatePoints) {
        const IndexType candidate_id = r_candidate.second->Id();
        KRATOS_ERROR_IF(candidate_id == 0 || candidate_id > rSlaveCurveData.size())
            << "::[IgaContactProcessGapSbmMortar]:: invalid slave candidate id "
            << candidate_id << std::endl;

        const auto& r_slave_curve = rSlaveCurveData[candidate_id - 1];
        if (!HaveCompatibleNormals(*rMasterCurve.pDeformedCurve, *r_slave_curve.pDeformedCurve)) {
            continue;
        }

        double slave_domain_begin = 0.0;
        double slave_domain_end = 0.0;
        KRATOS_ERROR_IF_NOT(GetCurveDomain(*r_slave_curve.pDeformedCurve, slave_domain_begin, slave_domain_end))
            << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve slave deformed domain."
            << std::endl;

        double slave_parameter = 0.5 * (slave_domain_begin + slave_domain_end);
        CoordinatesArrayType slave_projection = ZeroVector(3);
        double distance_squared = 0.0;
        int is_inside = 0;
        ProjectPointOntoCurve(
            *r_slave_curve.pDeformedCurve,
            rMasterPointDeformed,
            slave_parameter,
            slave_parameter,
            slave_projection,
            distance_squared,
            is_inside);

        if (is_inside == 0 || distance_squared > ProjectionDistanceLimitSquared ||
            distance_squared >= rBestProjection.DistanceSquared) {
            continue;
        }

        double slave_reference_parameter = 0.0;
        CoordinatesArrayType slave_reference_point = ZeroVector(3);
        double consistency_distance_sq = 0.0;
        MapDeformedPointToReferenceCurve(
            *r_slave_curve.pReferenceCurve,
            r_slave_curve.ReferenceGeometries,
            slave_parameter,
            slave_projection,
            slave_reference_parameter,
            slave_reference_point,
            consistency_distance_sq);

        if (consistency_distance_sq > ProjectionDistanceLimitSquared) {
            continue;
        }

        array_1d<double, 3> slave_normal = ZeroVector(3);
        ComputeCurveNormal(*r_slave_curve.pReferenceCurve, slave_reference_parameter, slave_normal);

        rBestProjection.IsValid = true;
        rBestProjection.pSlaveCurve = &r_slave_curve;
        rBestProjection.SlaveDeformedParameter = slave_parameter;
        rBestProjection.SlaveReferenceParameter = slave_reference_parameter;
        rBestProjection.SlaveDeformedPoint = slave_projection;
        rBestProjection.SlaveReferencePoint = slave_reference_point;
        rBestProjection.SlaveNormal = slave_normal;
        rBestProjection.DistanceSquared = distance_squared;
    }

    return rBestProjection.IsValid;
}

void PushUniqueParameter(
    std::vector<double>& rParameters,
    const double Parameter)
{
    const auto it = std::find_if(
        rParameters.begin(),
        rParameters.end(),
        [Parameter](const double ExistingParameter) {
            return std::abs(ExistingParameter - Parameter) <= parameter_tolerance;
        });

    if (it == rParameters.end()) {
        rParameters.push_back(Parameter);
    }
}

Vector MakeVector2(
    const double Value0,
    const double Value1)
{
    Vector result = ZeroVector(2);
    result[0] = Value0;
    result[1] = Value1;
    return result;
}

Vector MakeVector3(const CoordinatesArrayType& rCoordinates)
{
    Vector result = ZeroVector(3);
    result[0] = rCoordinates[0];
    result[1] = rCoordinates[1];
    result[2] = rCoordinates[2];
    return result;
}

double ComputeCharacteristicLength(
    const CoordinatesArrayType& rMasterPoint1,
    const CoordinatesArrayType& rMasterPoint2,
    const CoordinatesArrayType& rSlavePoint1,
    const CoordinatesArrayType& rSlavePoint2)
{
    array_1d<double, 3> center = ZeroVector(3);
    center += rMasterPoint1;
    center += rMasterPoint2;
    center += rSlavePoint1;
    center += rSlavePoint2;
    center *= 0.25;

    double radius_squared = 0.0;
    const auto update_radius = [&center, &radius_squared](const CoordinatesArrayType& rPoint) {
        const auto difference = center - rPoint;
        const double distance_squared = inner_prod(difference, difference);
        radius_squared = std::max(radius_squared, distance_squared);
    };

    update_radius(rMasterPoint1);
    update_radius(rMasterPoint2);
    update_radius(rSlavePoint1);
    update_radius(rSlavePoint2);

    return std::sqrt(radius_squared + 1.0e-15);
}

std::vector<double> BuildClippedCurveSpans(
    const GeometryType& rCurveGeometry,
    const double IntervalBegin,
    const double IntervalEnd)
{
    const GeometryType* p_sampling_geometry = &rCurveGeometry;
    if (rCurveGeometry.HasGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)) {
        const auto p_background_curve =
            rCurveGeometry.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        p_sampling_geometry = p_background_curve.get();
    }

    std::vector<double> curve_spans;
    p_sampling_geometry->SpansLocalSpace(curve_spans);
    KRATOS_ERROR_IF(curve_spans.size() < 2)
        << "::[IgaContactProcessGapSbmMortar]:: Invalid span sampling for overlap quadrature."
        << std::endl;

    const double interval_min = std::min(IntervalBegin, IntervalEnd);
    const double interval_max = std::max(IntervalBegin, IntervalEnd);

    std::vector<double> clipped_spans;
    clipped_spans.reserve(curve_spans.size() + 2);
    PushUniqueParameter(clipped_spans, interval_min);
    for (const double span_parameter : curve_spans) {
        if (span_parameter > interval_min + parameter_tolerance &&
            span_parameter < interval_max - parameter_tolerance) {
            PushUniqueParameter(clipped_spans, span_parameter);
        }
    }
    PushUniqueParameter(clipped_spans, interval_max);
    std::sort(clipped_spans.begin(), clipped_spans.end());

    return clipped_spans;
}

IntegrationPointsArrayType BuildIntegrationPointsForInterval(
    const GeometryType& rCurveGeometry,
    const double IntervalBegin,
    const double IntervalEnd,
    const std::size_t NumberOfPoints)
{
    KRATOS_ERROR_IF(NumberOfPoints == 0)
        << "::[IgaContactProcessGapSbmMortar]:: NumberOfPoints must be positive." << std::endl;
    KRATOS_ERROR_IF(NumberOfPoints > IntegrationPointUtilities::s_gauss_legendre.size())
        << "::[IgaContactProcessGapSbmMortar]:: Unsupported number of Gauss points: "
        << NumberOfPoints << std::endl;

    const double delta = IntervalEnd - IntervalBegin;
    const double abs_delta = std::abs(delta);
    KRATOS_ERROR_IF(abs_delta <= parameter_tolerance)
        << "::[IgaContactProcessGapSbmMortar]:: Degenerate overlap interval." << std::endl;

    const auto clipped_spans = BuildClippedCurveSpans(rCurveGeometry, IntervalBegin, IntervalEnd);
    KRATOS_ERROR_IF(clipped_spans.size() < 2)
        << "::[IgaContactProcessGapSbmMortar]:: Invalid clipped spans for overlap quadrature." << std::endl;

    const auto& gauss_table = IntegrationPointUtilities::s_gauss_legendre[NumberOfPoints - 1];
    const std::size_t number_of_spans = clipped_spans.size() - 1;
    IntegrationPointsArrayType integration_points(number_of_spans * NumberOfPoints);

    std::size_t point_counter = 0;
    for (std::size_t i_span = 0; i_span < number_of_spans; ++i_span) {
        const double span_begin = clipped_spans[i_span];
        const double span_end = clipped_spans[i_span + 1];
        const double span_delta = span_end - span_begin;
        const double span_length = std::abs(span_delta);

        KRATOS_ERROR_IF(span_length <= parameter_tolerance)
            << "::[IgaContactProcessGapSbmMortar]:: Degenerate clipped span in overlap quadrature."
            << std::endl;

        for (std::size_t i_gp = 0; i_gp < NumberOfPoints; ++i_gp) {
            const double xi = gauss_table[i_gp][0];
            const double wi = gauss_table[i_gp][1];
            integration_points[point_counter++] =
                IntegrationPoint<1>(span_begin + span_delta * xi, wi * span_length);
        }
    }

    return integration_points;
}

void RescaleIntegrationPointWeightsByCurveMetric(
    const GeometryType& rCurveGeometry,
    IntegrationPointsArrayType& rIntegrationPoints)
{
    const GeometryType* p_metric_geometry = &rCurveGeometry;
    if (rCurveGeometry.HasGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)) {
        const auto p_background_curve =
            rCurveGeometry.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        p_metric_geometry = p_background_curve.get();
    }

    for (auto& r_integration_point : rIntegrationPoints) {
        Matrix jacobian = ZeroMatrix(3, 1);
        const auto local_coordinates = MakeLocalCoordinate(r_integration_point.X());
        p_metric_geometry->Jacobian(jacobian, local_coordinates);

        const double tangent_norm = std::sqrt(
            jacobian(0, 0) * jacobian(0, 0) +
            jacobian(1, 0) * jacobian(1, 0) +
            jacobian(2, 0) * jacobian(2, 0));

        r_integration_point.SetWeight(r_integration_point.Weight() * tangent_norm);
    }
}

GeometryPointerType BuildDeformedCurve(
    const GeometryType& rCurveGeometry,
    const std::vector<Geometry<Node>::Pointer>& rReferenceGeometries,
    std::size_t& rNextAuxiliaryNodeId)
{
    KRATOS_ERROR_IF(rReferenceGeometries.empty())
        << "::[IgaContactProcessGapSbmMortar]:: Missing surrogate reference geometries." << std::endl;

    const GeometryType* p_sampling_geometry = &rCurveGeometry;
    if (rCurveGeometry.HasGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)) {
        const auto p_background_curve =
            rCurveGeometry.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        p_sampling_geometry = p_background_curve.get();
    }

    if (const auto* p_nurbs_curve =
            dynamic_cast<const NurbsCurveGeometryType*>(p_sampling_geometry)) {
        PointerVector<NodeType> deformed_control_points;
        deformed_control_points.reserve(p_nurbs_curve->PointsNumber());

        for (IndexType i = 0; i < p_nurbs_curve->PointsNumber(); ++i) {
            const auto& r_control_point = p_nurbs_curve->GetPoint(i);

            auto p_reference_node = Kratos::make_intrusive<NodeType>(
                rNextAuxiliaryNodeId++,
                r_control_point.X(),
                r_control_point.Y(),
                r_control_point.Z());
            p_reference_node->SetValue(NEIGHBOUR_GEOMETRIES, rReferenceGeometries);

            CoordinatesArrayType deformed_control_point = r_control_point.Coordinates();
            IgaSbmUtilities::GetDeformedPosition(*p_reference_node, deformed_control_point);

            deformed_control_points.push_back(Kratos::make_intrusive<NodeType>(
                rNextAuxiliaryNodeId++,
                deformed_control_point[0],
                deformed_control_point[1],
                deformed_control_point[2]));
        }

        NurbsCurveGeometryType::Pointer p_deformed_curve_geometry;
        const auto& r_weights = p_nurbs_curve->Weights();
        if (r_weights.size() == p_nurbs_curve->PointsNumber() && r_weights.size() > 0) {
            p_deformed_curve_geometry = Kratos::make_shared<NurbsCurveGeometryType>(
                deformed_control_points,
                p_nurbs_curve->PolynomialDegree(0),
                p_nurbs_curve->Knots(),
                r_weights);
        } else {
            p_deformed_curve_geometry = Kratos::make_shared<NurbsCurveGeometryType>(
                deformed_control_points,
                p_nurbs_curve->PolynomialDegree(0),
                p_nurbs_curve->Knots());
        }

        return Kratos::make_shared<BrepCurve<PointerVector<NodeType>, PointerVector<Point>>>(
            p_deformed_curve_geometry);
    }

    std::vector<double> curve_spans;
    p_sampling_geometry->SpansLocalSpace(curve_spans);
    KRATOS_ERROR_IF(curve_spans.size() < 2)
        << "::[IgaContactProcessGapSbmMortar]:: Invalid span sampling for curve geometry #"
        << rCurveGeometry.Id() << std::endl;

    std::vector<double> sample_parameters;
    sample_parameters.reserve(curve_spans.size() * 2 - 1);
    sample_parameters.push_back(curve_spans.front());
    for (std::size_t i = 0; i + 1 < curve_spans.size(); ++i) {
        const double left = curve_spans[i];
        const double right = curve_spans[i + 1];
        sample_parameters.push_back(0.5 * (left + right));
        sample_parameters.push_back(right);
    }

    PointerVector<NodeType> deformed_curve_points;
    deformed_curve_points.reserve(sample_parameters.size());
    Vector deformed_curve_knots(sample_parameters.size());

    for (std::size_t i = 0; i < sample_parameters.size(); ++i) {
        const CoordinatesArrayType sample_global =
            EvaluateCurvePoint(rCurveGeometry, sample_parameters[i]);

        auto p_sample_node = Kratos::make_intrusive<NodeType>(
            rNextAuxiliaryNodeId++,
            sample_global[0],
            sample_global[1],
            sample_global[2]);
        p_sample_node->SetValue(NEIGHBOUR_GEOMETRIES, rReferenceGeometries);

        CoordinatesArrayType sample_global_deformed = sample_global;
        IgaSbmUtilities::GetDeformedPosition(*p_sample_node, sample_global_deformed);

        deformed_curve_points.push_back(Kratos::make_intrusive<NodeType>(
            rNextAuxiliaryNodeId++,
            sample_global_deformed[0],
            sample_global_deformed[1],
            sample_global_deformed[2]));
        deformed_curve_knots[i] = sample_parameters[i];
    }

    auto p_deformed_curve_geometry = Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<NodeType>>>(
        deformed_curve_points,
        1,
        deformed_curve_knots);

    return Kratos::make_shared<BrepCurve<PointerVector<NodeType>, PointerVector<Point>>>(
        p_deformed_curve_geometry);
}

CurveData BuildCurveData(
    const GeometryPointerType& pCurveGeometry,
    std::size_t& rNextAuxiliaryNodeId)
{
    KRATOS_ERROR_IF_NOT(pCurveGeometry)
        << "::[IgaContactProcessGapSbmMortar]:: null curve geometry pointer." << std::endl;
    KRATOS_ERROR_IF_NOT(pCurveGeometry->Has(NEIGHBOUR_GEOMETRIES))
        << "::[IgaContactProcessGapSbmMortar]:: curve geometry #" << pCurveGeometry->Id()
        << " missing NEIGHBOUR_GEOMETRIES." << std::endl;

    CurveData curve_data;
    curve_data.pReferenceCurve = pCurveGeometry;
    curve_data.ReferenceGeometries = pCurveGeometry->GetValue(NEIGHBOUR_GEOMETRIES);
    KRATOS_ERROR_IF(curve_data.ReferenceGeometries.empty())
        << "::[IgaContactProcessGapSbmMortar]:: curve geometry #" << pCurveGeometry->Id()
        << " has empty NEIGHBOUR_GEOMETRIES." << std::endl;

    curve_data.pDeformedCurve = BuildDeformedCurve(
        *pCurveGeometry,
        curve_data.ReferenceGeometries,
        rNextAuxiliaryNodeId);

    return curve_data;
}

bool HaveCompatibleNormals(
    const GeometryType& rMasterDeformedCurve,
    const GeometryType& rSlaveDeformedCurve)
{
    double master_domain_begin = 0.0;
    double master_domain_end = 0.0;
    double slave_domain_begin = 0.0;
    double slave_domain_end = 0.0;

    KRATOS_ERROR_IF_NOT(GetCurveDomain(rMasterDeformedCurve, master_domain_begin, master_domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve master deformed domain." << std::endl;
    KRATOS_ERROR_IF_NOT(GetCurveDomain(rSlaveDeformedCurve, slave_domain_begin, slave_domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve slave deformed domain." << std::endl;

    array_1d<double, 3> master_normal = ZeroVector(3);
    array_1d<double, 3> slave_normal = ZeroVector(3);
    const bool has_master_normal =
        ComputeCurveNormal(rMasterDeformedCurve, 0.5 * (master_domain_begin + master_domain_end), master_normal);
    const bool has_slave_normal =
        ComputeCurveNormal(rSlaveDeformedCurve, 0.5 * (slave_domain_begin + slave_domain_end), slave_normal);

    if (!has_master_normal || !has_slave_normal) {
        return true;
    }

    return inner_prod(master_normal, slave_normal) < 0.0;
}

bool BuildOverlapDataFromMasterSegment(
    const GeometryType& rMasterDeformedCurve,
    const GeometryType& rMasterReferenceCurve,
    const std::vector<Geometry<Node>::Pointer>& rMasterReferenceGeometries,
    const double MasterSegmentBegin,
    const double MasterSegmentEnd,
    const GeometryType& rSlaveDeformedCurve,
    const GeometryType& rSlaveReferenceCurve,
    const std::vector<Geometry<Node>::Pointer>& rSlaveReferenceGeometries,
    const double ProjectionDistanceLimitSquared,
    OverlapData& rOverlapData)
{
    if (std::abs(MasterSegmentEnd - MasterSegmentBegin) <= parameter_tolerance) {
        return false;
    }

    double slave_domain_begin = 0.0;
    double slave_domain_end = 0.0;
    KRATOS_ERROR_IF_NOT(GetCurveDomain(rSlaveDeformedCurve, slave_domain_begin, slave_domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve slave deformed domain." << std::endl;

    const CoordinatesArrayType master_overlap_deformed_point_begin =
        EvaluateCurvePoint(rMasterDeformedCurve, MasterSegmentBegin);
    const CoordinatesArrayType master_overlap_deformed_point_end =
        EvaluateCurvePoint(rMasterDeformedCurve, MasterSegmentEnd);

    double slave_overlap_begin = 0.5 * (slave_domain_begin + slave_domain_end);
    double slave_overlap_end = slave_overlap_begin;
    CoordinatesArrayType slave_projection_begin = ZeroVector(3);
    CoordinatesArrayType slave_projection_end = ZeroVector(3);
    double slave_distance_begin_sq = 0.0;
    double slave_distance_end_sq = 0.0;
    int slave_inside_begin = 0;
    int slave_inside_end = 0;

    ProjectPointOntoCurve(
        rSlaveDeformedCurve,
        master_overlap_deformed_point_begin,
        slave_overlap_begin,
        slave_overlap_begin,
        slave_projection_begin,
        slave_distance_begin_sq,
        slave_inside_begin);
    ProjectPointOntoCurve(
        rSlaveDeformedCurve,
        master_overlap_deformed_point_end,
        slave_overlap_end,
        slave_overlap_end,
        slave_projection_end,
        slave_distance_end_sq,
        slave_inside_end);

    if (slave_inside_begin == 0 || slave_inside_end == 0) {
        return false;
    }
    if (slave_distance_begin_sq > ProjectionDistanceLimitSquared ||
        slave_distance_end_sq > ProjectionDistanceLimitSquared) {
        return false;
    }
    if (std::abs(slave_overlap_end - slave_overlap_begin) <= parameter_tolerance) {
        return false;
    }

    double master_reference_begin = 0.0;
    double master_reference_end = 0.0;
    double master_reference_begin_distance_sq = 0.0;
    double master_reference_end_distance_sq = 0.0;
    CoordinatesArrayType master_reference_point_begin = ZeroVector(3);
    CoordinatesArrayType master_reference_point_end = ZeroVector(3);
    MapDeformedPointToReferenceCurve(
        rMasterReferenceCurve,
        rMasterReferenceGeometries,
        MasterSegmentBegin,
        master_overlap_deformed_point_begin,
        master_reference_begin,
        master_reference_point_begin,
        master_reference_begin_distance_sq);
    MapDeformedPointToReferenceCurve(
        rMasterReferenceCurve,
        rMasterReferenceGeometries,
        MasterSegmentEnd,
        master_overlap_deformed_point_end,
        master_reference_end,
        master_reference_point_end,
        master_reference_end_distance_sq);

    double slave_reference_begin = 0.0;
    double slave_reference_end = 0.0;
    double slave_reference_begin_distance_sq = 0.0;
    double slave_reference_end_distance_sq = 0.0;
    CoordinatesArrayType slave_reference_point_begin = ZeroVector(3);
    CoordinatesArrayType slave_reference_point_end = ZeroVector(3);
    MapDeformedPointToReferenceCurve(
        rSlaveReferenceCurve,
        rSlaveReferenceGeometries,
        slave_overlap_begin,
        slave_projection_begin,
        slave_reference_begin,
        slave_reference_point_begin,
        slave_reference_begin_distance_sq);
    MapDeformedPointToReferenceCurve(
        rSlaveReferenceCurve,
        rSlaveReferenceGeometries,
        slave_overlap_end,
        slave_projection_end,
        slave_reference_end,
        slave_reference_point_end,
        slave_reference_end_distance_sq);

    if (master_reference_begin_distance_sq > ProjectionDistanceLimitSquared ||
        master_reference_end_distance_sq > ProjectionDistanceLimitSquared ||
        slave_reference_begin_distance_sq > ProjectionDistanceLimitSquared ||
        slave_reference_end_distance_sq > ProjectionDistanceLimitSquared) {
        return false;
    }
    if (std::abs(master_reference_end - master_reference_begin) <= parameter_tolerance ||
        std::abs(slave_reference_end - slave_reference_begin) <= parameter_tolerance) {
        return false;
    }

    rOverlapData.IsValid = true;
    rOverlapData.MasterDeformedBegin = MasterSegmentBegin;
    rOverlapData.MasterDeformedEnd = MasterSegmentEnd;
    rOverlapData.MasterReferenceBegin = master_reference_begin;
    rOverlapData.MasterReferenceEnd = master_reference_end;
    rOverlapData.SlaveDeformedBegin = slave_overlap_begin;
    rOverlapData.SlaveDeformedEnd = slave_overlap_end;
    rOverlapData.SlaveReferenceBegin = slave_reference_begin;
    rOverlapData.SlaveReferenceEnd = slave_reference_end;
    rOverlapData.MasterReferencePointBegin = master_reference_point_begin;
    rOverlapData.MasterReferencePointEnd = master_reference_point_end;
    rOverlapData.SlaveReferencePointBegin = slave_reference_point_begin;
    rOverlapData.SlaveReferencePointEnd = slave_reference_point_end;

    return true;
}

bool ComputeOverlapOnMaster(
    const GeometryType& rMasterDeformedCurve,
    const GeometryType& rMasterReferenceCurve,
    const std::vector<Geometry<Node>::Pointer>& rMasterReferenceGeometries,
    const GeometryType& rSlaveDeformedCurve,
    const GeometryType& rSlaveReferenceCurve,
    const std::vector<Geometry<Node>::Pointer>& rSlaveReferenceGeometries,
    const double ProjectionDistanceLimitSquared,
    OverlapData& rOverlapData)
{
    double master_domain_begin = 0.0;
    double master_domain_end = 0.0;
    double slave_domain_begin = 0.0;
    double slave_domain_end = 0.0;

    KRATOS_ERROR_IF_NOT(GetCurveDomain(rMasterDeformedCurve, master_domain_begin, master_domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve master deformed domain." << std::endl;
    KRATOS_ERROR_IF_NOT(GetCurveDomain(rSlaveDeformedCurve, slave_domain_begin, slave_domain_end))
        << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve slave deformed domain." << std::endl;

    const CoordinatesArrayType master_point_begin =
        EvaluateCurvePoint(rMasterDeformedCurve, master_domain_begin);
    const CoordinatesArrayType master_point_end =
        EvaluateCurvePoint(rMasterDeformedCurve, master_domain_end);
    const CoordinatesArrayType slave_point_begin =
        EvaluateCurvePoint(rSlaveDeformedCurve, slave_domain_begin);
    const CoordinatesArrayType slave_point_end =
        EvaluateCurvePoint(rSlaveDeformedCurve, slave_domain_end);

    std::vector<double> master_overlap_parameters;
    master_overlap_parameters.reserve(4);

    auto project_master_endpoint = [&](const double MasterParameter, const CoordinatesArrayType& rMasterPoint) {
        double slave_parameter = 0.5 * (slave_domain_begin + slave_domain_end);
        CoordinatesArrayType slave_projection = ZeroVector(3);
        double distance_squared = 0.0;
        int is_inside = 0;
        ProjectPointOntoCurve(
            rSlaveDeformedCurve,
            rMasterPoint,
            slave_parameter,
            slave_parameter,
            slave_projection,
            distance_squared,
            is_inside);

        if (is_inside != 0 && distance_squared <= ProjectionDistanceLimitSquared) {
            PushUniqueParameter(master_overlap_parameters, MasterParameter);
        }
    };

    auto project_slave_endpoint = [&](const CoordinatesArrayType& rSlavePoint) {
        double master_parameter = 0.5 * (master_domain_begin + master_domain_end);
        CoordinatesArrayType master_projection = ZeroVector(3);
        double distance_squared = 0.0;
        int is_inside = 0;
        ProjectPointOntoCurve(
            rMasterDeformedCurve,
            rSlavePoint,
            master_parameter,
            master_parameter,
            master_projection,
            distance_squared,
            is_inside);

        if (is_inside != 0 && distance_squared <= ProjectionDistanceLimitSquared) {
            PushUniqueParameter(master_overlap_parameters, master_parameter);
        }
    };

    project_master_endpoint(master_domain_begin, master_point_begin);
    project_master_endpoint(master_domain_end, master_point_end);
    project_slave_endpoint(slave_point_begin);
    project_slave_endpoint(slave_point_end);

    if (master_overlap_parameters.size() < 2) {
        return false;
    }

    std::sort(master_overlap_parameters.begin(), master_overlap_parameters.end());
    const double master_overlap_begin = master_overlap_parameters.front();
    const double master_overlap_end = master_overlap_parameters.back();

    if (std::abs(master_overlap_end - master_overlap_begin) <= parameter_tolerance) {
        return false;
    }

    return BuildOverlapDataFromMasterSegment(
        rMasterDeformedCurve,
        rMasterReferenceCurve,
        rMasterReferenceGeometries,
        master_overlap_begin,
        master_overlap_end,
        rSlaveDeformedCurve,
        rSlaveReferenceCurve,
        rSlaveReferenceGeometries,
        ProjectionDistanceLimitSquared,
        rOverlapData);
}

} // namespace

IgaContactProcessGapSbmMortar::IgaContactProcessGapSbmMortar(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mpModel(&rModel)
    , mParameters(ThisParameters)
{
    mEchoLevel = mParameters.Has("echo_level") ? mParameters["echo_level"].GetInt() : 0;

    KRATOS_ERROR_IF_NOT(ThisParameters.Has("analysis_model_part_name"))
        << "::[IgaContactProcessGapSbmMortar]:: Missing \"analysis_model_part_name\" parameter." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("contact_sub_model_part_name"))
        << "::[IgaContactProcessGapSbmMortar]:: Missing \"contact_sub_model_part_name\" parameter." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("contact_parameters"))
        << "::[IgaContactProcessGapSbmMortar]:: Missing \"contact_parameters\" parameter." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters["contact_parameters"].Has("slave_model_part"))
        << "::[IgaContactProcessGapSbmMortar]:: Missing \"contact_parameters/slave_model_part\" parameter." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters["contact_parameters"].Has("master_model_part"))
        << "::[IgaContactProcessGapSbmMortar]:: Missing \"contact_parameters/master_model_part\" parameter." << std::endl;

    std::string slave_model_part_name =
        mParameters["contact_parameters"]["slave_model_part"]["sub_model_part_name"].GetString();
    slave_model_part_name += ".";
    slave_model_part_name += mParameters["contact_parameters"]["slave_model_part"]["layer_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_model_part_name))
        << "ERROR: SLAVE MODEL PART " << slave_model_part_name << " NOT CREATED" << std::endl;
    mrSlaveModelPart = &(mpModel->GetModelPart(slave_model_part_name));
    mpPropSlave = mrSlaveModelPart->pGetProperties(
        mParameters["contact_parameters"]["slave_model_part"]["property_id"].GetInt());

    std::string master_model_part_name =
        mParameters["contact_parameters"]["master_model_part"]["sub_model_part_name"].GetString();
    master_model_part_name += ".";
    master_model_part_name += mParameters["contact_parameters"]["master_model_part"]["layer_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_model_part_name))
        << "ERROR: MASTER MODEL PART " << master_model_part_name << " NOT CREATED" << std::endl;
    mrMasterModelPart = &(mpModel->GetModelPart(master_model_part_name));
    mpPropMaster = mrMasterModelPart->pGetProperties(
        mParameters["contact_parameters"]["master_model_part"]["property_id"].GetInt());

    const std::string contact_model_part_name =
        mParameters["analysis_model_part_name"].GetString() +
        ".ContactInterface." +
        mParameters["contact_sub_model_part_name"].GetString();

    if (mpModel->HasModelPart(contact_model_part_name)) {
        mpContactModelPart = &(mpModel->GetModelPart(contact_model_part_name));
    } else {
        mpContactModelPart = &(mpModel->CreateModelPart(contact_model_part_name));
    }
}

void IgaContactProcessGapSbmMortar::ClearCreatedPairingConditions()
{
    if (mCreatedPairingConditionIds.empty()) {
        return;
    }

    auto& r_root_model_part = mpContactModelPart->GetRootModelPart();
    for (const auto condition_id : mCreatedPairingConditionIds) {
        if (r_root_model_part.HasCondition(condition_id)) {
            r_root_model_part.RemoveCondition(condition_id);
        }
    }

    mCreatedPairingConditionIds.clear();
}

void IgaContactProcessGapSbmMortar::Execute()
{
    ClearCreatedPairingConditions();

    ModelPart& r_contact_sub_model_part = mpContactModelPart->HasSubModelPart("contact")
        ? mpContactModelPart->GetSubModelPart("contact")
        : mpContactModelPart->CreateSubModelPart("contact");

    ModelPart& r_slave_curve_model_part = mrSlaveModelPart->HasSubModelPart("contact")
        ? mrSlaveModelPart->GetSubModelPart("contact")
        : *mrSlaveModelPart;
    ModelPart& r_master_curve_model_part = mrMasterModelPart->HasSubModelPart("contact")
        ? mrMasterModelPart->GetSubModelPart("contact")
        : *mrMasterModelPart;

    if (r_slave_curve_model_part.NumberOfGeometries() == 0 || r_master_curve_model_part.NumberOfGeometries() == 0) {
        KRATOS_WARNING("IgaContactProcessGapSbmMortar")
            << "Master or slave contact model part has no geometries. Nothing to do." << std::endl;
        return;
    }

    if (!mpContactConditionProperties) {
        mpContactConditionProperties = Kratos::make_shared<Properties>(*mpPropMaster);
        if (mpContactConditionProperties->GetSubProperties().size() == 0) {
            mpContactConditionProperties->AddSubProperties(mpPropMaster);
            if (mpPropSlave != mpPropMaster) {
                mpContactConditionProperties->AddSubProperties(mpPropSlave);
            }
        }
    }

    PointVector slave_center_points;
    std::vector<CurveData> slave_curve_data;
    slave_curve_data.reserve(r_slave_curve_model_part.NumberOfGeometries());

    PointVector master_center_points;
    std::vector<CurveData> master_curve_data;
    master_curve_data.reserve(r_master_curve_model_part.NumberOfGeometries());

    std::size_t next_auxiliary_node_id = mpContactModelPart->GetRootModelPart().NumberOfNodes() > 0
        ? mpContactModelPart->GetRootModelPart().Nodes().back().Id() + 1
        : 1;

    for (auto& r_slave_geometry : r_slave_curve_model_part.Geometries()) {
        slave_curve_data.push_back(
            BuildCurveData(
                r_slave_curve_model_part.pGetGeometry(r_slave_geometry.Id()),
                next_auxiliary_node_id));
        const auto& r_center = slave_curve_data.back().pDeformedCurve->Center();
        slave_center_points.push_back(Kratos::make_intrusive<PointType>(
            static_cast<IndexType>(slave_curve_data.size()),
            r_center.X(),
            r_center.Y(),
            r_center.Z()));
    }

    for (auto& r_master_geometry : r_master_curve_model_part.Geometries()) {
        master_curve_data.push_back(
            BuildCurveData(
                r_master_curve_model_part.pGetGeometry(r_master_geometry.Id()),
                next_auxiliary_node_id));
        const auto& r_center = master_curve_data.back().pDeformedCurve->Center();
        master_center_points.push_back(Kratos::make_intrusive<PointType>(
            static_cast<IndexType>(master_curve_data.size()),
            r_center.X(),
            r_center.Y(),
            r_center.Z()));
    }

    DynamicBins slave_bins(slave_center_points.begin(), slave_center_points.end());
    DynamicBins master_bins(master_center_points.begin(), master_center_points.end());

    const SizeType number_of_considered_neighbours = static_cast<SizeType>(
        mParameters["numbered_considered_neighbours"].GetInt());
    KRATOS_ERROR_IF(number_of_considered_neighbours == 0)
        << "::[IgaContactProcessGapSbmMortar]:: \"numbered_considered_neighbours\" must be > 0."
        << std::endl;

    const SizeType max_considered_slave_neighbours =
        std::min(number_of_considered_neighbours, slave_center_points.size());
    const SizeType max_considered_master_neighbours =
        std::min(number_of_considered_neighbours, master_center_points.size());

    if (max_considered_slave_neighbours == 0 || max_considered_master_neighbours == 0) {
        KRATOS_WARNING("IgaContactProcessGapSbmMortar")
            << "No slave or master center points available for contact search. Nothing to do." << std::endl;
        return;
    }

    const Vector master_knot_span_sizes =
        mrMasterModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    const double projection_distance_limit =
        master_knot_span_sizes[0] * mParameters["projection_distance_scale"].GetDouble();

    KRATOS_WATCH(projection_distance_limit)
    exit(0);
    const double projection_distance_fallback =
        master_knot_span_sizes[0] * mParameters["projection_distance_fallback_scale"].GetDouble();
    const double projection_distance_limit_sq = projection_distance_limit * projection_distance_limit;
    const double projection_distance_fallback_sq = projection_distance_fallback * projection_distance_fallback;

    const double max_slave_bins_search_radius =
        ComputeMaxCenterSearchRadius(slave_center_points, projection_distance_limit);
    const double max_master_bins_search_radius =
        ComputeMaxCenterSearchRadius(master_center_points, projection_distance_limit);

    PointVector slave_bins_results(slave_center_points.size());
    DistanceVector slave_bins_distances(slave_center_points.size());
    PointVector master_bins_results(master_center_points.size());
    DistanceVector master_bins_distances(master_center_points.size());
    std::vector<CandidatePointData> candidate_points;
    candidate_points.reserve(std::max(slave_center_points.size(), master_center_points.size()));

    const std::string& r_condition_name = mParameters["contact_condition_name"].GetString();
    KRATOS_ERROR_IF_NOT(KratosComponents<Condition>::Has(r_condition_name))
        << "::[IgaContactProcessGapSbmMortar]:: condition \"" << r_condition_name
        << "\" is not registered." << std::endl;
    const Condition& r_reference_condition = KratosComponents<Condition>::Get(r_condition_name);

    const SizeType number_of_shape_function_derivatives =
        static_cast<SizeType>(mParameters["shape_function_derivatives_order"].GetInt());
    const SizeType number_of_integration_points_per_span =
        static_cast<SizeType>(mParameters["number_of_integration_points_per_span"].GetInt());
    KRATOS_ERROR_IF(number_of_shape_function_derivatives == 0)
        << "::[IgaContactProcessGapSbmMortar]:: \"shape_function_derivatives_order\" must be > 0."
        << std::endl;
    KRATOS_ERROR_IF(number_of_integration_points_per_span == 0)
        << "::[IgaContactProcessGapSbmMortar]:: \"number_of_integration_points_per_span\" must be > 0."
        << std::endl;

    ModelPart::ConditionsContainerType new_contact_conditions;
    IndexType next_condition_id = mpContactModelPart->GetRootModelPart().Conditions().empty()
        ? 1
        : mpContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;

    std::vector<CouplingEntity> coupling_entities;
    coupling_entities.reserve(slave_curve_data.size() * 2);

    const auto append_coupling_entity = [&coupling_entities](
                                           const CurveData* pSlaveCurve,
                                           const CurveData* pMasterCurve,
                                           const double MasterDeformedBegin,
                                           const double MasterDeformedEnd,
                                           const double MasterReferenceBegin,
                                           const double MasterReferenceEnd) {
        if (!pSlaveCurve || !pMasterCurve) {
            return;
        }
        if (std::abs(MasterDeformedEnd - MasterDeformedBegin) <= parameter_tolerance ||
            std::abs(MasterReferenceEnd - MasterReferenceBegin) <= parameter_tolerance) {
            return;
        }

        CouplingEntity coupling_entity;
        coupling_entity.pSlaveCurve = pSlaveCurve;
        coupling_entity.pMasterCurve = pMasterCurve;
        if (MasterDeformedBegin <= MasterDeformedEnd) {
            coupling_entity.MasterDeformedBegin = MasterDeformedBegin;
            coupling_entity.MasterDeformedEnd = MasterDeformedEnd;
            coupling_entity.MasterReferenceBegin = MasterReferenceBegin;
            coupling_entity.MasterReferenceEnd = MasterReferenceEnd;
        } else {
            coupling_entity.MasterDeformedBegin = MasterDeformedEnd;
            coupling_entity.MasterDeformedEnd = MasterDeformedBegin;
            coupling_entity.MasterReferenceBegin = MasterReferenceEnd;
            coupling_entity.MasterReferenceEnd = MasterReferenceBegin;
        }
        coupling_entities.push_back(coupling_entity);
    };

    for (const auto& r_slave_curve : slave_curve_data) {
        double slave_domain_begin = 0.0;
        double slave_domain_end = 0.0;
        KRATOS_ERROR_IF_NOT(GetCurveDomain(*r_slave_curve.pDeformedCurve, slave_domain_begin, slave_domain_end))
            << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve slave deformed domain." << std::endl;

        const CoordinatesArrayType slave_deformed_point_begin =
            EvaluateCurvePoint(*r_slave_curve.pDeformedCurve, slave_domain_begin);
        const CoordinatesArrayType slave_deformed_point_end =
            EvaluateCurvePoint(*r_slave_curve.pDeformedCurve, slave_domain_end);

        EndpointProjectionData endpoint_projection_begin;
        EndpointProjectionData endpoint_projection_end;

        array_1d<double, 3> slave_normal_begin = ZeroVector(3);
        array_1d<double, 3> slave_normal_end = ZeroVector(3);
        const bool has_slave_normal_begin = ComputeCurveNormal(*r_slave_curve.pDeformedCurve, slave_domain_begin, slave_normal_begin);
        const bool has_slave_normal_end = ComputeCurveNormal(*r_slave_curve.pDeformedCurve, slave_domain_end, slave_normal_end);

        if (!has_slave_normal_begin || !has_slave_normal_end) {
            continue;
        }

        const bool has_begin_projection = FindBestMasterProjectionForSlavePoint(
            slave_deformed_point_begin,
            slave_normal_begin,
            r_slave_curve,
            master_curve_data,
            master_bins,
            max_considered_master_neighbours,
            projection_distance_limit,
            max_master_bins_search_radius,
            projection_distance_limit_sq,
            master_bins_results,
            master_bins_distances,
            candidate_points,
            endpoint_projection_begin);
        const bool has_end_projection = FindBestMasterProjectionForSlavePoint(
            slave_deformed_point_end,
            slave_normal_end,
            r_slave_curve,
            master_curve_data,
            master_bins,
            max_considered_master_neighbours,
            projection_distance_limit,
            max_master_bins_search_radius,
            projection_distance_limit_sq,
            master_bins_results,
            master_bins_distances,
            candidate_points,
            endpoint_projection_end);

        if (!has_begin_projection || !has_end_projection) {
            continue;
        }

        if (endpoint_projection_begin.pMasterCurve == endpoint_projection_end.pMasterCurve) {
            if (endpoint_projection_begin.MasterDeformedParameter <= endpoint_projection_end.MasterDeformedParameter) {
                append_coupling_entity(
                    &r_slave_curve,
                    endpoint_projection_begin.pMasterCurve,
                    endpoint_projection_begin.MasterDeformedParameter,
                    endpoint_projection_end.MasterDeformedParameter,
                    endpoint_projection_begin.MasterReferenceParameter,
                    endpoint_projection_end.MasterReferenceParameter);
            } else {
                append_coupling_entity(
                    &r_slave_curve,
                    endpoint_projection_begin.pMasterCurve,
                    endpoint_projection_end.MasterDeformedParameter,
                    endpoint_projection_begin.MasterDeformedParameter,
                    endpoint_projection_end.MasterReferenceParameter,
                    endpoint_projection_begin.MasterReferenceParameter);
            }
            continue;
        }

        const auto create_entity_from_endpoint = [&](const EndpointProjectionData& rEndpointProjection,
                                                     const CoordinatesArrayType& rOppositeSlavePoint) {
            double master_domain_begin = 0.0;
            double master_domain_end = 0.0;
            KRATOS_ERROR_IF_NOT(GetCurveDomain(
                *rEndpointProjection.pMasterCurve->pDeformedCurve,
                master_domain_begin,
                master_domain_end))
                << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve master deformed domain." << std::endl;

            const CoordinatesArrayType master_point_begin =
                EvaluateCurvePoint(*rEndpointProjection.pMasterCurve->pDeformedCurve, master_domain_begin);
            const CoordinatesArrayType master_point_end =
                EvaluateCurvePoint(*rEndpointProjection.pMasterCurve->pDeformedCurve, master_domain_end);

            const auto delta_begin = master_point_begin - rOppositeSlavePoint;
            const auto delta_end = master_point_end - rOppositeSlavePoint;
            const double distance_begin_sq = inner_prod(delta_begin, delta_begin);
            const double distance_end_sq = inner_prod(delta_end, delta_end);

            if (distance_begin_sq <= distance_end_sq) {
                append_coupling_entity(
                    &r_slave_curve,
                    rEndpointProjection.pMasterCurve,
                    master_domain_begin,
                    rEndpointProjection.MasterDeformedParameter,
                    master_domain_begin,
                    rEndpointProjection.MasterReferenceParameter);
            } else {
                append_coupling_entity(
                    &r_slave_curve,
                    rEndpointProjection.pMasterCurve,
                    rEndpointProjection.MasterDeformedParameter,
                    master_domain_end,
                    rEndpointProjection.MasterReferenceParameter,
                    master_domain_end);
            }
        };

        create_entity_from_endpoint(endpoint_projection_begin, slave_deformed_point_end);
        create_entity_from_endpoint(endpoint_projection_end, slave_deformed_point_begin);
    }

    for (const auto& r_master_curve : master_curve_data) {
        double master_domain_begin = 0.0;
        double master_domain_end = 0.0;
        KRATOS_ERROR_IF_NOT(GetCurveDomain(*r_master_curve.pReferenceCurve, master_domain_begin, master_domain_end))
            << "::[IgaContactProcessGapSbmMortar]:: Could not retrieve master reference domain." << std::endl;

        std::vector<std::pair<double, double>> covered_intervals;
        covered_intervals.reserve(coupling_entities.size());
        for (const auto& r_coupling_entity : coupling_entities) {
            if (r_coupling_entity.pMasterCurve != &r_master_curve) {
                continue;
            }

            const double clamped_begin = std::max(
                master_domain_begin,
                std::min(master_domain_end, r_coupling_entity.MasterReferenceBegin));
            const double clamped_end = std::max(
                master_domain_begin,
                std::min(master_domain_end, r_coupling_entity.MasterReferenceEnd));

            if (std::abs(clamped_end - clamped_begin) <= parameter_tolerance) {
                continue;
            }

            covered_intervals.emplace_back(
                std::min(clamped_begin, clamped_end),
                std::max(clamped_begin, clamped_end));
        }

        std::sort(
            covered_intervals.begin(),
            covered_intervals.end(),
            [](const auto& rA, const auto& rB) {
                return rA.first < rB.first;
            });

        std::vector<std::pair<double, double>> missing_intervals;
        double current_position = master_domain_begin;

        if (covered_intervals.empty()) {
            missing_intervals.emplace_back(master_domain_begin, master_domain_end);
        } else {
            for (const auto& r_interval : covered_intervals) {
                if (r_interval.first > current_position + parameter_tolerance) {
                    missing_intervals.emplace_back(current_position, r_interval.first);
                }
                current_position = std::max(current_position, r_interval.second);
            }

            if (master_domain_end > current_position + parameter_tolerance) {
                missing_intervals.emplace_back(current_position, master_domain_end);
            }
        }

        for (const auto& r_missing_interval : missing_intervals) {
            if (std::abs(r_missing_interval.second - r_missing_interval.first) <= parameter_tolerance) {
                continue;
            }

            const double master_reference_midpoint = 0.5 * (r_missing_interval.first + r_missing_interval.second);
            const CoordinatesArrayType master_midpoint_deformed = EvaluateDeformedReferenceCurvePoint(
                *r_master_curve.pReferenceCurve,
                r_master_curve.ReferenceGeometries,
                master_reference_midpoint);

            SlaveProjectionData slave_projection;
            const bool has_slave_projection = FindBestSlaveProjectionForMasterPoint(
                master_midpoint_deformed,
                r_master_curve,
                slave_curve_data,
                slave_bins,
                max_considered_slave_neighbours,
                projection_distance_limit,
                max_slave_bins_search_radius,
                projection_distance_limit_sq,
                slave_bins_results,
                slave_bins_distances,
                candidate_points,
                slave_projection);

            if (!has_slave_projection) {
                continue;
            }

            append_coupling_entity(
                slave_projection.pSlaveCurve,
                &r_master_curve,
                r_missing_interval.first,
                r_missing_interval.second,
                r_missing_interval.first,
                r_missing_interval.second);
        }
    }

    for (const auto& r_coupling_entity : coupling_entities) {
        const auto& r_master_curve = *r_coupling_entity.pMasterCurve;
        // const double master_reference_midpoint = 0.5 * (
        //     r_coupling_entity.MasterReferenceBegin + r_coupling_entity.MasterReferenceEnd);
        // const auto master_reference_midpoint_point =
        //     EvaluateCurvePoint(*r_master_curve.pReferenceCurve, master_reference_midpoint);

        // KRATOS_WATCH("Processing coupling entity:");
        // KRATOS_WATCH("  Master curve ID: " << r_master_curve.pReferenceCurve->Id());
        // KRATOS_WATCH("  Master deformed begin: " << r_coupling_entity.MasterDeformedBegin);
        // KRATOS_WATCH("  Master deformed end: " << r_coupling_entity.MasterDeformedEnd);
        // KRATOS_WATCH("  Master reference begin: " << r_coupling_entity.MasterReferenceBegin);
        // KRATOS_WATCH("  Master reference end: " << r_coupling_entity.MasterReferenceEnd);
        // KRATOS_WATCH("  Master reference midpoint: " << master_reference_midpoint_point);

        IntegrationInfo master_integration_info = r_master_curve.pReferenceCurve->GetDefaultIntegrationInfo();
        master_integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);

        const auto master_integration_points = BuildIntegrationPointsForInterval(
            *r_master_curve.pReferenceCurve,
            r_coupling_entity.MasterReferenceBegin,
            r_coupling_entity.MasterReferenceEnd,
            number_of_integration_points_per_span);

        auto master_integration_points_rescaled = master_integration_points;
        RescaleIntegrationPointWeightsByCurveMetric(
            *r_master_curve.pReferenceCurve,
            master_integration_points_rescaled);

        GeometryType::GeometriesArrayType master_qp_geometries;
        r_master_curve.pReferenceCurve->CreateQuadraturePointGeometries(
            master_qp_geometries,
            number_of_shape_function_derivatives,
            master_integration_points_rescaled,
            master_integration_info);

        for (IndexType i_gp = 0; i_gp < master_qp_geometries.size(); ++i_gp) {
            auto p_master_qp_geometry = master_qp_geometries(i_gp);
            p_master_qp_geometry->SetValue(NEIGHBOUR_GEOMETRIES, r_master_curve.ReferenceGeometries);

            auto p_master_gp_node = Kratos::make_intrusive<NodeType>(
                next_auxiliary_node_id++,
                p_master_qp_geometry->Center().X(),
                p_master_qp_geometry->Center().Y(),
                p_master_qp_geometry->Center().Z());
            p_master_gp_node->SetValue(NEIGHBOUR_GEOMETRIES, r_master_curve.ReferenceGeometries);

            CoordinatesArrayType master_gp_deformed = p_master_qp_geometry->Center().Coordinates();
            IgaSbmUtilities::GetDeformedPosition(*p_master_gp_node, master_gp_deformed);

            SlaveProjectionData slave_projection;
            const bool has_slave_projection = FindBestSlaveProjectionForMasterPoint(
                master_gp_deformed,
                r_master_curve,
                slave_curve_data,
                slave_bins,
                max_considered_slave_neighbours,
                projection_distance_limit,
                max_slave_bins_search_radius,
                projection_distance_limit_sq,
                slave_bins_results,
                slave_bins_distances,
                candidate_points,
                slave_projection);

            if (!has_slave_projection) {
                continue;
            }

            OverlapData overlap_data;
            const bool has_overlap_data = BuildOverlapDataFromMasterSegment(
                *r_master_curve.pDeformedCurve,
                *r_master_curve.pReferenceCurve,
                r_master_curve.ReferenceGeometries,
                r_coupling_entity.MasterDeformedBegin,
                r_coupling_entity.MasterDeformedEnd,
                *slave_projection.pSlaveCurve->pDeformedCurve,
                *slave_projection.pSlaveCurve->pReferenceCurve,
                slave_projection.pSlaveCurve->ReferenceGeometries,
                projection_distance_limit_sq,
                overlap_data);

            if (!has_overlap_data) {
                continue;
            }

            const double characteristic_length = ComputeCharacteristicLength(
                overlap_data.MasterReferencePointBegin,
                overlap_data.MasterReferencePointEnd,
                overlap_data.SlaveReferencePointBegin,
                overlap_data.SlaveReferencePointEnd);
            Vector characteristic_length_vector = ZeroVector(3);
            characteristic_length_vector[0] = characteristic_length;

            std::vector<Vector> parameter_space_corners;
            parameter_space_corners.reserve(4);
            parameter_space_corners.push_back(
                MakeVector2(overlap_data.MasterDeformedBegin, overlap_data.MasterDeformedEnd));
            parameter_space_corners.push_back(
                MakeVector2(overlap_data.MasterReferenceBegin, overlap_data.MasterReferenceEnd));
            parameter_space_corners.push_back(
                MakeVector2(overlap_data.SlaveDeformedBegin, overlap_data.SlaveDeformedEnd));
            parameter_space_corners.push_back(
                MakeVector2(overlap_data.SlaveReferenceBegin, overlap_data.SlaveReferenceEnd));

            auto p_projection_node = Kratos::make_intrusive<NodeType>(
                next_auxiliary_node_id++,
                slave_projection.SlaveReferencePoint[0],
                slave_projection.SlaveReferencePoint[1],
                slave_projection.SlaveReferencePoint[2]);
            p_projection_node->SetValue(NEIGHBOUR_GEOMETRIES, slave_projection.pSlaveCurve->ReferenceGeometries);
            p_projection_node->SetValue(NORMAL, slave_projection.SlaveNormal);
            p_master_qp_geometry->SetValue(PROJECTION_NODE, p_projection_node);

            Condition::Pointer p_new_contact_condition = r_reference_condition.Create(
                next_condition_id++,
                p_master_qp_geometry,
                mpContactConditionProperties);

            p_new_contact_condition->SetValue(IDENTIFIER, "MASTER");
            p_new_contact_condition->SetValue(NEIGHBOUR_GEOMETRIES, r_master_curve.ReferenceGeometries);
            p_new_contact_condition->SetValue(KNOT_SPAN_SIZES, master_knot_span_sizes);
            p_new_contact_condition->SetValue(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_vector);
            p_new_contact_condition->SetValue(PARAMETER_SPACE_CORNERS, parameter_space_corners);
            p_new_contact_condition->SetValue(PROJECTION_NODE_COORDINATES, MakeVector3(slave_projection.SlaveReferencePoint));
            p_new_contact_condition->SetValue(PROJECTION_NODE_ID, static_cast<int>(slave_projection.pSlaveCurve->pReferenceCurve->Id()));
            p_new_contact_condition->SetValue(
                INTEGRATION_POINTS_MASTER,
                std::vector<Vector>{MakeVector3(p_master_qp_geometry->Center().Coordinates())});
            p_new_contact_condition->SetValue(
                INTEGRATION_POINTS_SLAVE,
                std::vector<Vector>{MakeVector3(slave_projection.SlaveReferencePoint)});
            p_new_contact_condition->SetValue(BREP_ID, static_cast<int>(r_master_curve.pReferenceCurve->Id()));

            new_contact_conditions.push_back(p_new_contact_condition);
            mCreatedPairingConditionIds.push_back(p_new_contact_condition->Id());
        }
    }

    if (!new_contact_conditions.empty()) {
        r_contact_sub_model_part.AddConditions(
            new_contact_conditions.begin(),
            new_contact_conditions.end());
        EntitiesUtilities::InitializeEntities<Condition>(r_contact_sub_model_part);
    }
}

} // namespace Kratos
