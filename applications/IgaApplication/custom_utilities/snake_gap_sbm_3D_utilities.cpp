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
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

// Project includes
#include "custom_utilities/snake_gap_sbm_3D_utilities.h"
#include "custom_processes/snake_gap_sbm_process.h"
#include "custom_utilities/snake_gap_sbm_3D_utilities.h"

namespace Kratos
{

namespace
{
//FIXME: Verify if this tolerance is appropriate for all gap volume scales.
constexpr double MinimumGapVolumeForQuadrature = 1.0e-14;

struct ParameterSpaceBounds
{
    double MinU = 0.0;
    double MaxU = 0.0;
    double MinV = 0.0;
    double MaxV = 0.0;
    double MinW = 0.0;
    double MaxW = 0.0;
};

ParameterSpaceBounds ReadParameterSpaceBounds3D(
    const ModelPart& rModelPart,
    const char* pCaller)
{
    ParameterSpaceBounds bounds;

    if (rModelPart.Has(PATCH_PARAMETER_SPACE_CORNERS)) {
        const Matrix& r_patch_corners = rModelPart.GetValue(PATCH_PARAMETER_SPACE_CORNERS);

        KRATOS_ERROR_IF(r_patch_corners.size1() < 3 || r_patch_corners.size2() < 2)
            << "[" << pCaller << "] PATCH_PARAMETER_SPACE_CORNERS must be at least 3x2.\n";

        bounds.MinU = r_patch_corners(0, 0);
        bounds.MaxU = r_patch_corners(0, 1);
        bounds.MinV = r_patch_corners(1, 0);
        bounds.MaxV = r_patch_corners(1, 1);
        bounds.MinW = r_patch_corners(2, 0);
        bounds.MaxW = r_patch_corners(2, 1);

        return bounds;
    }

    KRATOS_ERROR_IF_NOT(rModelPart.Has(PARAMETER_SPACE_CORNERS))
        << "[" << pCaller << "] ModelPart '" << rModelPart.Name()
        << "' does not have PATCH_PARAMETER_SPACE_CORNERS nor PARAMETER_SPACE_CORNERS.\n";

    const auto& r_parameter_space_corners = rModelPart.GetValue(PARAMETER_SPACE_CORNERS);

    KRATOS_ERROR_IF(r_parameter_space_corners.size() < 3)
        << "[" << pCaller << "] PARAMETER_SPACE_CORNERS must contain at least three vectors.\n";

    KRATOS_ERROR_IF(r_parameter_space_corners[0].size() < 2 ||
                    r_parameter_space_corners[1].size() < 2 ||
                    r_parameter_space_corners[2].size() < 2)
        << "[" << pCaller << "] PARAMETER_SPACE_CORNERS entries must contain [min, max] values."
        << " If the model part stores a multipatch box list, set PATCH_PARAMETER_SPACE_CORNERS on the patch model part.\n";

    bounds.MinU = r_parameter_space_corners[0][0];
    bounds.MaxU = r_parameter_space_corners[0][1];
    bounds.MinV = r_parameter_space_corners[1][0];
    bounds.MaxV = r_parameter_space_corners[1][1];
    bounds.MinW = r_parameter_space_corners[2][0];
    bounds.MaxW = r_parameter_space_corners[2][1];

    return bounds;
}

std::string GapSpanTypeToString(const SnakeGapSbm3DUtilities::GapSpanType Type)
{
    switch (Type) {
        case SnakeGapSbm3DUtilities::GapSpanType::Type1:
            return "TYPE_1";
        case SnakeGapSbm3DUtilities::GapSpanType::Type2:
            return "TYPE_2";
        case SnakeGapSbm3DUtilities::GapSpanType::Type3:
            return "TYPE_3";
        default:
            return "UNDEFINED";
    }
}

int GapSpanTypePriority(const SnakeGapSbm3DUtilities::GapSpanType Type)
{
    switch (Type) {
        case SnakeGapSbm3DUtilities::GapSpanType::Type1:
            return 1;
        case SnakeGapSbm3DUtilities::GapSpanType::Type2:
            return 2;
        case SnakeGapSbm3DUtilities::GapSpanType::Type3:
            return 3;
        default:
            return 100;
    }
}

double CalculateTetrahedronVolume(
    const Node& rNode0,
    const Node& rNode1,
    const Node& rNode2,
    const Node& rNode3)
{
    const array_1d<double, 3> edge_0 =
        rNode1.Coordinates() - rNode0.Coordinates();
    const array_1d<double, 3> edge_1 =
        rNode2.Coordinates() - rNode0.Coordinates();
    const array_1d<double, 3> edge_2 =
        rNode3.Coordinates() - rNode0.Coordinates();

    return std::abs(
        inner_prod(edge_0, MathUtils<double>::CrossProduct(edge_1, edge_2))) / 6.0;
}

array_1d<double, 3> ClosestPointOnTriangle(
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rA,
    const array_1d<double, 3>& rB,
    const array_1d<double, 3>& rC)
{
    const array_1d<double, 3> ab = rB - rA;
    const array_1d<double, 3> ac = rC - rA;
    const array_1d<double, 3> ap = rPoint - rA;
    const double d1 = inner_prod(ab, ap);
    const double d2 = inner_prod(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        return rA;
    }

    const array_1d<double, 3> bp = rPoint - rB;
    const double d3 = inner_prod(ab, bp);
    const double d4 = inner_prod(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        return rB;
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        array_1d<double, 3> closest_point = rA;
        noalias(closest_point) += v * ab;
        return closest_point;
    }

    const array_1d<double, 3> cp = rPoint - rC;
    const double d5 = inner_prod(ab, cp);
    const double d6 = inner_prod(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        return rC;
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        array_1d<double, 3> closest_point = rA;
        noalias(closest_point) += w * ac;
        return closest_point;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        array_1d<double, 3> closest_point = rB;
        noalias(closest_point) += w * (rC - rB);
        return closest_point;
    }

    const double denominator = 1.0 / (va + vb + vc);
    const double v = vb * denominator;
    const double w = vc * denominator;
    array_1d<double, 3> closest_point = rA;
    noalias(closest_point) += v * ab;
    noalias(closest_point) += w * ac;
    return closest_point;
}

} // unnamed namespace

SnakeGapSbm3DUtilities::SnakeGapSbm3DUtilities(const int EchoLevel)
    : mEchoLevel(EchoLevel)
{
}

void SnakeGapSbm3DUtilities::SetGapApproximationOrder(
    const std::size_t GapApproximationOrder)
{
    mGapApproximationOrder = GapApproximationOrder == 0 ? 1 : GapApproximationOrder;

    KRATOS_ERROR_IF(mGapApproximationOrder > 2)
        << "[SnakeGapSbm3DUtilities::SetGapApproximationOrder] "
        << "Only gap_approximation_order <= 2 is currently implemented "
        << "for 3D curved gap geometries. Requested order: "
        << mGapApproximationOrder << ".\n";
}

void SnakeGapSbm3DUtilities::SetStoreGapDebugGeometries(
    const bool StoreGapDebugGeometries)
{
    mStoreGapDebugGeometries = StoreGapDebugGeometries;
}

SnakeGapSbm3DUtilities::KnotSpanGridInfo
SnakeGapSbm3DUtilities::CreateKnotSpanGridInfo(
    const ModelPart& rBackgroundModelPart) const
{
    const char* p_caller = "SnakeGapSbm3DUtilities::CreateKnotSpanGridInfo";

    KRATOS_ERROR_IF_NOT(rBackgroundModelPart.Has(KNOT_SPAN_SIZES))
        << "[" << p_caller << "] ModelPart '" << rBackgroundModelPart.Name()
        << "' does not have KNOT_SPAN_SIZES.\n";

    const Vector& r_knot_span_sizes = rBackgroundModelPart.GetValue(KNOT_SPAN_SIZES);

    KRATOS_ERROR_IF(r_knot_span_sizes.size() < 3)
        << "[" << p_caller << "] KNOT_SPAN_SIZES must contain at least three entries.\n";

    KnotSpanGridInfo info;

    info.SpanSizeX = r_knot_span_sizes[0];
    info.SpanSizeY = r_knot_span_sizes[1];
    info.SpanSizeZ = r_knot_span_sizes[2];

    KRATOS_ERROR_IF(info.SpanSizeX <= 0.0 || info.SpanSizeY <= 0.0 || info.SpanSizeZ <= 0.0)
        << "[" << p_caller << "] Knot span sizes must be positive. Got "
        << info.SpanSizeX << ", " << info.SpanSizeY << ", " << info.SpanSizeZ << ".\n";

    const auto bounds = ReadParameterSpaceBounds3D(rBackgroundModelPart, p_caller);

    info.MinU = bounds.MinU;
    info.MaxU = bounds.MaxU;
    info.MinV = bounds.MinV;
    info.MaxV = bounds.MaxV;
    info.MinW = bounds.MinW;
    info.MaxW = bounds.MaxW;

    KRATOS_ERROR_IF(info.MaxU <= info.MinU || info.MaxV <= info.MinV || info.MaxW <= info.MinW)
        << "[" << p_caller << "] Invalid parameter space bounds: "
        << "U=[" << info.MinU << ", " << info.MaxU << "], "
        << "V=[" << info.MinV << ", " << info.MaxV << "], "
        << "W=[" << info.MinW << ", " << info.MaxW << "].\n";

    info.NumberOfSpansX = ComputeSpanCount(info.MaxU - info.MinU, info.SpanSizeX, "u", p_caller);
    info.NumberOfSpansY = ComputeSpanCount(info.MaxV - info.MinV, info.SpanSizeY, "v", p_caller);
    info.NumberOfSpansZ = ComputeSpanCount(info.MaxW - info.MinW, info.SpanSizeZ, "w", p_caller);

    return info;
}

std::size_t SnakeGapSbm3DUtilities::ComputeSpanCount(
    const double DomainLength,
    const double SpanSize,
    const char* pDirection,
    const char* pCaller) const
{
    const double spans_real = DomainLength / SpanSize;
    const double rounded_spans = std::round(spans_real);

    constexpr double absolute_tolerance = 1.0e-10;
    constexpr double relative_tolerance = 1.0e-9;

    const double tolerance = std::max(
        absolute_tolerance,
        relative_tolerance * std::max(1.0, std::abs(rounded_spans)));

    KRATOS_ERROR_IF(rounded_spans <= 0.0)
        << "[" << pCaller << "] Non-positive number of knot spans in " << pDirection << ".\n";

    KRATOS_ERROR_IF(std::abs(spans_real - rounded_spans) > tolerance)
        << "[" << pCaller << "] Non-integer number of knot spans in " << pDirection
        << " (" << std::setprecision(17) << spans_real << ")"
        << " from domain length " << DomainLength
        << " and span size " << SpanSize << ".\n";

    return static_cast<std::size_t>(rounded_spans);
}

std::size_t SnakeGapSbm3DUtilities::ComputeSpanIndex(
    const double Coordinate,
    const double MinValue,
    const double MaxValue,
    const double SpanSize,
    const std::size_t NumberOfSpans,
    const double Tolerance,
    const char* pCaller) const
{
    double clamped_coordinate = Coordinate;

    if (Coordinate < MinValue) {
        KRATOS_ERROR_IF(Coordinate < MinValue - Tolerance)
            << "[" << pCaller << "] Coordinate " << Coordinate
            << " below minimum parameter range " << MinValue << ".\n";
        clamped_coordinate = MinValue;
    } else if (Coordinate > MaxValue) {
        KRATOS_ERROR_IF(Coordinate > MaxValue + Tolerance)
            << "[" << pCaller << "] Coordinate " << Coordinate
            << " above maximum parameter range " << MaxValue << ".\n";
        clamped_coordinate = MaxValue;
    }

    std::size_t span_index = static_cast<std::size_t>(
        std::floor((clamped_coordinate - MinValue) / SpanSize + Tolerance));

    if (span_index >= NumberOfSpans) {
        span_index = NumberOfSpans - 1;
    }

    return span_index;
}

int SnakeGapSbm3DUtilities::ComputeGridLineIndex(
    const double Coordinate,
    const double MinValue,
    const double SpanSize,
    const std::size_t NumberOfSpans) const
{
    const double relative_coordinate = (Coordinate - MinValue) / SpanSize;
    int grid_index = static_cast<int>(std::round(relative_coordinate));

    if (grid_index < 0) {
        grid_index = 0;
    } else if (grid_index > static_cast<int>(NumberOfSpans)) {
        grid_index = static_cast<int>(NumberOfSpans);
    }

    return grid_index;
}

std::vector<std::size_t> SnakeGapSbm3DUtilities::ComputeNodeSpanIndices(
    const double Coordinate,
    const double MinValue,
    const double MaxValue,
    const double SpanSize,
    const std::size_t NumberOfSpans,
    const double Tolerance,
    const double BoundaryTolerance,
    const char* pCaller) const
{
    std::vector<std::size_t> span_indices;

    const std::size_t base_span_index = ComputeSpanIndex(
        Coordinate,
        MinValue,
        MaxValue,
        SpanSize,
        NumberOfSpans,
        Tolerance,
        pCaller);

    span_indices.push_back(base_span_index);

    const double clamped_coordinate = std::max(MinValue, std::min(MaxValue, Coordinate));
    const double relative_grid_coordinate = (clamped_coordinate - MinValue) / SpanSize;
    const double nearest_grid_line = std::round(relative_grid_coordinate);
    const double distance_to_grid_line = std::abs(relative_grid_coordinate - nearest_grid_line) * SpanSize;

    if (distance_to_grid_line <= BoundaryTolerance) {
        const int grid_line_index = static_cast<int>(nearest_grid_line);

        if (grid_line_index > 0) {
            span_indices.push_back(static_cast<std::size_t>(grid_line_index - 1));
        }

        if (grid_line_index < static_cast<int>(NumberOfSpans)) {
            span_indices.push_back(static_cast<std::size_t>(grid_line_index));
        }
    }

    std::sort(span_indices.begin(), span_indices.end());
    span_indices.erase(std::unique(span_indices.begin(), span_indices.end()), span_indices.end());

    return span_indices;
}

std::pair<std::size_t, std::size_t>
SnakeGapSbm3DUtilities::ComputeSpanIndexRange(
    const double MinCoordinate,
    const double MaxCoordinate,
    const double MinValue,
    const double MaxValue,
    const double SpanSize,
    const std::size_t NumberOfSpans,
    const double Tolerance,
    const double BoundaryTolerance,
    const char* rCaller) const
{
    KRATOS_ERROR_IF(NumberOfSpans == 0)
        << "[" << rCaller << "] NumberOfSpans is zero.\n";

    KRATOS_ERROR_IF(SpanSize <= 0.0)
        << "[" << rCaller << "] SpanSize must be positive. Got "
        << SpanSize << ".\n";

    KRATOS_ERROR_IF(MaxCoordinate < MinCoordinate)
        << "[" << rCaller << "] Invalid coordinate range: ["
        << MinCoordinate << ", " << MaxCoordinate << "].\n";

    if (MaxCoordinate < MinValue - BoundaryTolerance ||
        MinCoordinate > MaxValue + BoundaryTolerance) {
        KRATOS_WARNING("SnakeGapSbm3DUtilities")
            << "[" << rCaller << "] Bounding box range ["
            << MinCoordinate << ", " << MaxCoordinate
            << "] is outside domain ["
            << MinValue << ", " << MaxValue << "]. Clamping.\n";
    }

    const double expanded_min =
        std::max(MinValue, MinCoordinate - BoundaryTolerance);

    const double expanded_max =
        std::min(MaxValue, MaxCoordinate + BoundaryTolerance);

    const auto first_indices = ComputeNodeSpanIndices(
        expanded_min,
        MinValue,
        MaxValue,
        SpanSize,
        NumberOfSpans,
        Tolerance,
        BoundaryTolerance,
        rCaller);

    const auto last_indices = ComputeNodeSpanIndices(
        expanded_max,
        MinValue,
        MaxValue,
        SpanSize,
        NumberOfSpans,
        Tolerance,
        BoundaryTolerance,
        rCaller);

    KRATOS_ERROR_IF(first_indices.empty() || last_indices.empty())
        << "[" << rCaller << "] Failed to compute span range for ["
        << MinCoordinate << ", " << MaxCoordinate << "].\n";

    std::size_t first_index = NumberOfSpans - 1;
    std::size_t last_index = 0;

    for (const auto index : first_indices) {
        first_index = std::min(first_index, index);
        last_index = std::max(last_index, index);
    }

    for (const auto index : last_indices) {
        first_index = std::min(first_index, index);
        last_index = std::max(last_index, index);
    }

    KRATOS_ERROR_IF(first_index > last_index)
        << "[" << rCaller << "] Invalid computed span range: "
        << first_index << " > " << last_index << ".\n";

    return std::make_pair(first_index, last_index);
}

std::size_t SnakeGapSbm3DUtilities::FlattenSpanColumnIndex(
    const KnotSpanGridInfo& rGridInfo,
    const std::size_t J,
    const std::size_t K) const
{
    return J * rGridInfo.NumberOfSpansZ + K;
}

std::size_t SnakeGapSbm3DUtilities::FindSpanNnzIndex(
    const KnotSpanSkinBinsCSR& rSkinBins,
    const SpanKey3D& rSpan) const
{
    const auto& r_grid_info = rSkinBins.GridInfo;

    if (!IsSpanInsideDomain(rSpan, r_grid_info)) {
        return static_cast<std::size_t>(-1);
    }

    const std::size_t span_i = static_cast<std::size_t>(rSpan.I);
    const std::size_t span_j = static_cast<std::size_t>(rSpan.J);
    const std::size_t span_k = static_cast<std::size_t>(rSpan.K);

    const std::size_t column_index = FlattenSpanColumnIndex(r_grid_info, span_j, span_k);

    const auto row_begin = rSkinBins.Occupancy.index1_data()[span_i];
    const auto row_end = rSkinBins.Occupancy.index1_data()[span_i + 1];

    for (auto nnz_index = row_begin; nnz_index < row_end; ++nnz_index) {
        if (rSkinBins.Occupancy.index2_data()[nnz_index] == column_index) {
            return nnz_index;
        }
    }

    return static_cast<std::size_t>(-1);
}

bool SnakeGapSbm3DUtilities::IsSpanInsideDomain(
    const SpanKey3D& rSpan,
    const KnotSpanGridInfo& rGridInfo) const
{
    return rSpan.I >= 0 &&
           rSpan.J >= 0 &&
           rSpan.K >= 0 &&
           rSpan.I < static_cast<int>(rGridInfo.NumberOfSpansX) &&
           rSpan.J < static_cast<int>(rGridInfo.NumberOfSpansY) &&
           rSpan.K < static_cast<int>(rGridInfo.NumberOfSpansZ);
}

array_1d<double, 3> SnakeGapSbm3DUtilities::SpanCenter(
    const SpanKey3D& rSpan,
    const KnotSpanGridInfo& rGridInfo) const
{
    array_1d<double, 3> center = ZeroVector(3);

    center[0] = rGridInfo.MinU + (static_cast<double>(rSpan.I) + 0.5) * rGridInfo.SpanSizeX;
    center[1] = rGridInfo.MinV + (static_cast<double>(rSpan.J) + 0.5) * rGridInfo.SpanSizeY;
    center[2] = rGridInfo.MinW + (static_cast<double>(rSpan.K) + 0.5) * rGridInfo.SpanSizeZ;

    return center;
}

void SnakeGapSbm3DUtilities::ConditionBoundingBox(
    const Geometry<Node>& rGeometry,
    array_1d<double, 3>& rMinPoint,
    array_1d<double, 3>& rMaxPoint) const
{
    rMinPoint[0] = std::numeric_limits<double>::max();
    rMinPoint[1] = std::numeric_limits<double>::max();
    rMinPoint[2] = std::numeric_limits<double>::max();

    rMaxPoint[0] = -std::numeric_limits<double>::max();
    rMaxPoint[1] = -std::numeric_limits<double>::max();
    rMaxPoint[2] = -std::numeric_limits<double>::max();

    for (IndexType point_index = 0; point_index < rGeometry.size(); ++point_index) {
        const auto& r_coordinates = rGeometry[point_index].Coordinates();

        for (IndexType axis = 0; axis < 3; ++axis) {
            rMinPoint[axis] = std::min(rMinPoint[axis], r_coordinates[axis]);
            rMaxPoint[axis] = std::max(rMaxPoint[axis], r_coordinates[axis]);
        }
    }
}

std::size_t ComputeSpanCount3D(
    const double domain_length,
    const double span_size,
    const char* pDirection,
    const char* pCaller)
{
    const double spans_real = domain_length / span_size;
    const double rounded_spans = std::round(spans_real);
    constexpr double absolute_tolerance = 1.0e-10;
    constexpr double relative_tolerance = 1.0e-9;
    const double tolerance = std::max(
        absolute_tolerance,
        relative_tolerance * std::max(1.0, std::abs(rounded_spans)));

    KRATOS_ERROR_IF(rounded_spans <= 0.0)
        << "[" << pCaller << "] Non-positive number of knot spans in " << pDirection << ".\n";
    KRATOS_ERROR_IF(std::abs(spans_real - rounded_spans) > tolerance)
        << "[" << pCaller << "] Non-integer number of knot spans in " << pDirection
        << " (" << std::setprecision(17) << spans_real << ")"
        << " computed from domain length " << domain_length
        << " and span size " << span_size << ".\n";

    return static_cast<std::size_t>(rounded_spans);
}

SnakeGapSbm3DUtilities::KnotSpanSkinBinsCSR
SnakeGapSbm3DUtilities::BuildSkinBinsPerKnotSpan3D(
    const ModelPart& rSkinSubModelPart,
    const KnotSpanGridInfo& rGridInfo) const
{
    const char* p_caller = "SnakeGapSbm3DUtilities::BuildSkinBinsPerKnotSpan3D";

    KnotSpanSkinBinsCSR skin_bins;
    skin_bins.GridInfo = rGridInfo;

    const std::size_t number_of_flattened_yz_spans =
        rGridInfo.NumberOfSpansY * rGridInfo.NumberOfSpansZ;

    skin_bins.Occupancy.resize(
        rGridInfo.NumberOfSpansX,
        number_of_flattened_yz_spans,
        false);

    if (rSkinSubModelPart.NumberOfNodes() == 0 &&
        rSkinSubModelPart.NumberOfConditions() == 0) {
        return skin_bins;
    }

    const double min_span_size = std::min({
        rGridInfo.SpanSizeX,
        rGridInfo.SpanSizeY,
        rGridInfo.SpanSizeZ});

    const double tolerance = 1.0e-10;
    const double boundary_tolerance = 1.0e-6 * min_span_size;

    std::vector<std::vector<std::size_t>> column_indices_per_row(rGridInfo.NumberOfSpansX);

    auto add_span_to_pattern = [&](const std::size_t SpanIndexX,
                                   const std::size_t SpanIndexY,
                                   const std::size_t SpanIndexZ) {
        column_indices_per_row[SpanIndexX].push_back(
            FlattenSpanColumnIndex(rGridInfo, SpanIndexY, SpanIndexZ));
    };

    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const auto span_indices_x = ComputeNodeSpanIndices(
            r_node.X(),
            rGridInfo.MinU,
            rGridInfo.MaxU,
            rGridInfo.SpanSizeX,
            rGridInfo.NumberOfSpansX,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_indices_y = ComputeNodeSpanIndices(
            r_node.Y(),
            rGridInfo.MinV,
            rGridInfo.MaxV,
            rGridInfo.SpanSizeY,
            rGridInfo.NumberOfSpansY,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_indices_z = ComputeNodeSpanIndices(
            r_node.Z(),
            rGridInfo.MinW,
            rGridInfo.MaxW,
            rGridInfo.SpanSizeZ,
            rGridInfo.NumberOfSpansZ,
            tolerance,
            boundary_tolerance,
            p_caller);

        for (const std::size_t span_index_x : span_indices_x) {
            for (const std::size_t span_index_y : span_indices_y) {
                for (const std::size_t span_index_z : span_indices_z) {
                    add_span_to_pattern(span_index_x, span_index_y, span_index_z);
                }
            }
        }
    }

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.size() == 0) {
            continue;
        }

        array_1d<double, 3> min_point;
        array_1d<double, 3> max_point;
        ConditionBoundingBox(r_geometry, min_point, max_point);

        const auto span_range_x = ComputeSpanIndexRange(
            min_point[0],
            max_point[0],
            rGridInfo.MinU,
            rGridInfo.MaxU,
            rGridInfo.SpanSizeX,
            rGridInfo.NumberOfSpansX,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_range_y = ComputeSpanIndexRange(
            min_point[1],
            max_point[1],
            rGridInfo.MinV,
            rGridInfo.MaxV,
            rGridInfo.SpanSizeY,
            rGridInfo.NumberOfSpansY,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_range_z = ComputeSpanIndexRange(
            min_point[2],
            max_point[2],
            rGridInfo.MinW,
            rGridInfo.MaxW,
            rGridInfo.SpanSizeZ,
            rGridInfo.NumberOfSpansZ,
            tolerance,
            boundary_tolerance,
            p_caller);

        for (std::size_t span_index_x = span_range_x.first; span_index_x <= span_range_x.second; ++span_index_x) {
            for (std::size_t span_index_y = span_range_y.first; span_index_y <= span_range_y.second; ++span_index_y) {
                for (std::size_t span_index_z = span_range_z.first; span_index_z <= span_range_z.second; ++span_index_z) {
                    add_span_to_pattern(span_index_x, span_index_y, span_index_z);
                }
            }
        }
    }

    std::size_t number_of_non_zero_entries = 0;
    for (auto& r_column_indices : column_indices_per_row) {
        std::sort(r_column_indices.begin(), r_column_indices.end());
        r_column_indices.erase(
            std::unique(r_column_indices.begin(), r_column_indices.end()),
            r_column_indices.end());

        number_of_non_zero_entries += r_column_indices.size();
    }

    skin_bins.Occupancy.reserve(number_of_non_zero_entries);

    std::vector<std::unordered_map<std::size_t, std::size_t>> occupancy_index_lookup(rGridInfo.NumberOfSpansX);
    for (auto& r_lookup : occupancy_index_lookup) {
        r_lookup.reserve(8);
    }

    std::size_t nnz_counter = 0;
    for (std::size_t span_index_x = 0; span_index_x < rGridInfo.NumberOfSpansX; ++span_index_x) {
        for (const std::size_t span_column_index : column_indices_per_row[span_index_x]) {
            skin_bins.Occupancy.push_back(span_index_x, span_column_index, 0.0);
            occupancy_index_lookup[span_index_x].emplace(span_column_index, nnz_counter++);
        }
    }

    std::vector<std::vector<NodePointerType>> nodes_per_non_zero(number_of_non_zero_entries);
    std::vector<std::vector<ConditionPointerType>> conditions_per_non_zero(number_of_non_zero_entries);

    auto find_non_zero_index = [&](const std::size_t SpanIndexX,
                                   const std::size_t SpanIndexY,
                                   const std::size_t SpanIndexZ) {
        const std::size_t span_column_index =
            FlattenSpanColumnIndex(rGridInfo, SpanIndexY, SpanIndexZ);

        const auto it = occupancy_index_lookup[SpanIndexX].find(span_column_index);
        return it != occupancy_index_lookup[SpanIndexX].end() ?
            it->second :
            static_cast<std::size_t>(-1);
    };

    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const auto span_indices_x = ComputeNodeSpanIndices(
            r_node.X(),
            rGridInfo.MinU,
            rGridInfo.MaxU,
            rGridInfo.SpanSizeX,
            rGridInfo.NumberOfSpansX,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_indices_y = ComputeNodeSpanIndices(
            r_node.Y(),
            rGridInfo.MinV,
            rGridInfo.MaxV,
            rGridInfo.SpanSizeY,
            rGridInfo.NumberOfSpansY,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_indices_z = ComputeNodeSpanIndices(
            r_node.Z(),
            rGridInfo.MinW,
            rGridInfo.MaxW,
            rGridInfo.SpanSizeZ,
            rGridInfo.NumberOfSpansZ,
            tolerance,
            boundary_tolerance,
            p_caller);

        for (const std::size_t span_index_x : span_indices_x) {
            for (const std::size_t span_index_y : span_indices_y) {
                for (const std::size_t span_index_z : span_indices_z) {
                    const std::size_t non_zero_index =
                        find_non_zero_index(span_index_x, span_index_y, span_index_z);

                    KRATOS_DEBUG_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
                        << "[" << p_caller << "] Node nonzero not found in CSR pattern.\n";

                    nodes_per_non_zero[non_zero_index].push_back(
                        rSkinSubModelPart.pGetNode(r_node.Id()));

                    skin_bins.Occupancy.value_data()[non_zero_index] += 1.0;
                }
            }
        }
    }

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.size() == 0) {
            continue;
        }

        array_1d<double, 3> min_point;
        array_1d<double, 3> max_point;
        ConditionBoundingBox(r_geometry, min_point, max_point);

        const auto span_range_x = ComputeSpanIndexRange(
            min_point[0],
            max_point[0],
            rGridInfo.MinU,
            rGridInfo.MaxU,
            rGridInfo.SpanSizeX,
            rGridInfo.NumberOfSpansX,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_range_y = ComputeSpanIndexRange(
            min_point[1],
            max_point[1],
            rGridInfo.MinV,
            rGridInfo.MaxV,
            rGridInfo.SpanSizeY,
            rGridInfo.NumberOfSpansY,
            tolerance,
            boundary_tolerance,
            p_caller);

        const auto span_range_z = ComputeSpanIndexRange(
            min_point[2],
            max_point[2],
            rGridInfo.MinW,
            rGridInfo.MaxW,
            rGridInfo.SpanSizeZ,
            rGridInfo.NumberOfSpansZ,
            tolerance,
            boundary_tolerance,
            p_caller);

        for (std::size_t span_index_x = span_range_x.first; span_index_x <= span_range_x.second; ++span_index_x) {
            for (std::size_t span_index_y = span_range_y.first; span_index_y <= span_range_y.second; ++span_index_y) {
                for (std::size_t span_index_z = span_range_z.first; span_index_z <= span_range_z.second; ++span_index_z) {
                    const std::size_t non_zero_index =
                        find_non_zero_index(span_index_x, span_index_y, span_index_z);

                    KRATOS_DEBUG_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
                        << "[" << p_caller << "] Condition nonzero not found in CSR pattern.\n";

                    conditions_per_non_zero[non_zero_index].push_back(
                        rSkinSubModelPart.pGetCondition(r_condition.Id()));

                    skin_bins.Occupancy.value_data()[non_zero_index] += 1.0;
                }
            }
        }
    }

    skin_bins.CellDataByNnz.resize(number_of_non_zero_entries);

    for (std::size_t nnz_index = 0; nnz_index < number_of_non_zero_entries; ++nnz_index) {
        auto& r_cell_data = skin_bins.CellDataByNnz[nnz_index];

        auto& r_nodes = nodes_per_non_zero[nnz_index];
        std::sort(
            r_nodes.begin(),
            r_nodes.end(),
            [](const NodePointerType& pA, const NodePointerType& pB) {
                return pA->Id() < pB->Id();
            });

        r_nodes.erase(
            std::unique(
                r_nodes.begin(),
                r_nodes.end(),
                [](const NodePointerType& pA, const NodePointerType& pB) {
                    return pA->Id() == pB->Id();
                }),
            r_nodes.end());

        r_cell_data.Nodes = std::move(r_nodes);

        auto& r_conditions = conditions_per_non_zero[nnz_index];
        std::sort(
            r_conditions.begin(),
            r_conditions.end(),
            [](const ConditionPointerType& pA, const ConditionPointerType& pB) {
                return pA->Id() < pB->Id();
            });

        r_conditions.erase(
            std::unique(
                r_conditions.begin(),
                r_conditions.end(),
                [](const ConditionPointerType& pA, const ConditionPointerType& pB) {
                    return pA->Id() == pB->Id();
                }),
            r_conditions.end());

        r_cell_data.Conditions = std::move(r_conditions);
    }

    return skin_bins;
}

SnakeGapSbm3DUtilities::GridPointKey3D
SnakeGapSbm3DUtilities::ComputeGridPointKey(
    const array_1d<double, 3>& rPoint,
    const KnotSpanGridInfo& rGridInfo) const
{
    GridPointKey3D key;

    key.I = ComputeGridLineIndex(
        rPoint[0],
        rGridInfo.MinU,
        rGridInfo.SpanSizeX,
        rGridInfo.NumberOfSpansX);

    key.J = ComputeGridLineIndex(
        rPoint[1],
        rGridInfo.MinV,
        rGridInfo.SpanSizeY,
        rGridInfo.NumberOfSpansY);

    key.K = ComputeGridLineIndex(
        rPoint[2],
        rGridInfo.MinW,
        rGridInfo.SpanSizeZ,
        rGridInfo.NumberOfSpansZ);

    return key;
}

int SnakeGapSbm3DUtilities::ComputeFaceNormalAxis(
    const array_1d<double, 3>& rNormal) const
{
    int normal_axis = 0;

    for (int axis = 1; axis < 3; ++axis) {
        if (std::abs(rNormal[axis]) > std::abs(rNormal[normal_axis])) {
            normal_axis = axis;
        }
    }

    return normal_axis;
}

int SnakeGapSbm3DUtilities::ComputeFaceNormalSign(
    const array_1d<double, 3>& rNormal,
    const int NormalAxis) const
{
    return rNormal[NormalAxis] >= 0.0 ? 1 : -1;
}

SnakeGapSbm3DUtilities::SpanKey3D
SnakeGapSbm3DUtilities::ComputeExternalSpanFromSurrogateFace(
    const SurrogateFaceData& rFaceData,
    const KnotSpanGridInfo& rGridInfo) const
{
    SpanKey3D external_span;

    external_span.I = std::min({
        rFaceData.GridNodes[0].I,
        rFaceData.GridNodes[1].I,
        rFaceData.GridNodes[2].I,
        rFaceData.GridNodes[3].I});

    external_span.J = std::min({
        rFaceData.GridNodes[0].J,
        rFaceData.GridNodes[1].J,
        rFaceData.GridNodes[2].J,
        rFaceData.GridNodes[3].J});

    external_span.K = std::min({
        rFaceData.GridNodes[0].K,
        rFaceData.GridNodes[1].K,
        rFaceData.GridNodes[2].K,
        rFaceData.GridNodes[3].K});

    const int fixed_grid_index =
        rFaceData.NormalAxis == 0 ? rFaceData.GridNodes[0].I :
        rFaceData.NormalAxis == 1 ? rFaceData.GridNodes[0].J :
                                    rFaceData.GridNodes[0].K;

    if (rFaceData.NormalAxis == 0) {
        external_span.I = rFaceData.NormalSign > 0 ? fixed_grid_index : fixed_grid_index - 1;
    } else if (rFaceData.NormalAxis == 1) {
        external_span.J = rFaceData.NormalSign > 0 ? fixed_grid_index : fixed_grid_index - 1;
    } else {
        external_span.K = rFaceData.NormalSign > 0 ? fixed_grid_index : fixed_grid_index - 1;
    }

    KRATOS_ERROR_IF_NOT(IsSpanInsideDomain(external_span, rGridInfo))
        << "[SnakeGapSbm3DUtilities::ComputeExternalSpanFromSurrogateFace] External span "
        << SpanToString(external_span) << " from surrogate condition #"
        << rFaceData.pCondition->Id() << " is outside the knot-span grid.\n";

    return external_span;
}

SnakeGapSbm3DUtilities::SpanKey3D
SnakeGapSbm3DUtilities::ComputeActiveSpanFromSurrogateFace(
    const SurrogateFaceData& rFaceData) const
{
    SpanKey3D active_span = rFaceData.ExternalSpan;

    if (rFaceData.NormalAxis == 0) {
        active_span.I -= rFaceData.NormalSign;
    } else if (rFaceData.NormalAxis == 1) {
        active_span.J -= rFaceData.NormalSign;
    } else {
        active_span.K -= rFaceData.NormalSign;
    }

    return active_span;
}

SnakeGapSbm3DUtilities::SurrogateFaceData
SnakeGapSbm3DUtilities::CreateSurrogateFaceData(
    const Condition& rCondition,
    const KnotSpanGridInfo& rGridInfo) const
{
    const char* p_caller = "SnakeGapSbm3DUtilities::CreateSurrogateFaceData";

    const auto& r_geometry = rCondition.GetGeometry();

    KRATOS_ERROR_IF(r_geometry.PointsNumber() != 4)
        << "[" << p_caller << "] Surrogate condition #" << rCondition.Id()
        << " must have exactly four nodes. Got " << r_geometry.PointsNumber() << ".\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(NORMAL))
        << "[" << p_caller << "] Surrogate condition #" << rCondition.Id()
        << " does not have NORMAL.\n";

    const Vector& r_normal_vector = rCondition.GetValue(NORMAL);

    KRATOS_ERROR_IF(r_normal_vector.size() < 3)
        << "[" << p_caller << "] NORMAL on surrogate condition #" << rCondition.Id()
        << " must have three components.\n";

    SurrogateFaceData face_data;
    face_data.pCondition = &rCondition;

    for (IndexType axis = 0; axis < 3; ++axis) {
        face_data.Normal[axis] = r_normal_vector[axis];
    }

    for (IndexType local_index = 0; local_index < 4; ++local_index) {
        face_data.Nodes[local_index] =
            rCondition.pGetGeometry()->pGetPoint(local_index);

        face_data.GridNodes[local_index] =
            ComputeGridPointKey(r_geometry[local_index].Coordinates(), rGridInfo);
    }

    const double normal_norm = norm_2(face_data.Normal);

    KRATOS_ERROR_IF(normal_norm <= 1.0e-16)
        << "[" << p_caller << "] Zero NORMAL on surrogate condition #"
        << rCondition.Id() << ".\n";

    face_data.Normal /= normal_norm;

    face_data.NormalAxis = ComputeFaceNormalAxis(face_data.Normal);
    face_data.NormalSign = ComputeFaceNormalSign(face_data.Normal, face_data.NormalAxis);

    face_data.ExternalSpan = ComputeExternalSpanFromSurrogateFace(face_data, rGridInfo);
    face_data.ActiveSpan = ComputeActiveSpanFromSurrogateFace(face_data);

    KRATOS_ERROR_IF_NOT(IsSpanInsideDomain(face_data.ActiveSpan, rGridInfo))
        << "[" << p_caller << "] Active span " << SpanToString(face_data.ActiveSpan)
        << " from surrogate condition #" << rCondition.Id()
        << " is outside the knot-span grid.\n";

    return face_data;
}

std::vector<SnakeGapSbm3DUtilities::SurrogateFaceData>
SnakeGapSbm3DUtilities::BuildSurrogateFaceDataVector(
    const ModelPart& rSurrogateSubModelPart,
    const KnotSpanGridInfo& rGridInfo) const
{
    std::vector<SurrogateFaceData> surrogate_faces;
    surrogate_faces.reserve(rSurrogateSubModelPart.NumberOfConditions());

    for (const auto& r_condition : rSurrogateSubModelPart.Conditions()) {
        surrogate_faces.push_back(
            CreateSurrogateFaceData(r_condition, rGridInfo));
    }

    return surrogate_faces;
}

std::set<SnakeGapSbm3DUtilities::SpanKey3D>
SnakeGapSbm3DUtilities::ExtractActiveSpans(
    const std::vector<SurrogateFaceData>& rSurrogateFaces) const
{
    std::set<SpanKey3D> active_spans;

    for (const auto& r_face_data : rSurrogateFaces) {
        active_spans.insert(r_face_data.ActiveSpan);
    }

    return active_spans;
}

void SnakeGapSbm3DUtilities::AddUniqueActiveSpan(
    std::vector<SpanKey3D>& rSpans,
    const SpanKey3D& rSpan) const
{
    if (std::find(rSpans.begin(), rSpans.end(), rSpan) == rSpans.end()) {
        rSpans.push_back(rSpan);
    }
}

void SnakeGapSbm3DUtilities::AddUniqueConditionId(
    std::vector<IndexType>& rConditionIds,
    const IndexType ConditionId) const
{
    if (std::find(rConditionIds.begin(), rConditionIds.end(), ConditionId) == rConditionIds.end()) {
        rConditionIds.push_back(ConditionId);
    }
}

void SnakeGapSbm3DUtilities::RegisterExternalSpanCandidate(
    ExternalSpanDataMap& rExternalSpans,
    const SpanKey3D& rExternalSpan,
    const GapSpanType Type,
    const SpanKey3D& rAdjacentActiveSpan,
    const IndexType AdjacentSurrogateConditionId) const
{
    auto& r_data = rExternalSpans[rExternalSpan];

    r_data.Key = rExternalSpan;

    if (GapSpanTypePriority(Type) < GapSpanTypePriority(r_data.Type)) {
        r_data.Type = Type;
    }

    AddUniqueActiveSpan(r_data.AdjacentActiveSpans, rAdjacentActiveSpan);

    if (AdjacentSurrogateConditionId != 0) {
        AddUniqueConditionId(r_data.AdjacentSurrogateConditionIds, AdjacentSurrogateConditionId);
    }
}

SnakeGapSbm3DUtilities::ExternalSpanDataMap
SnakeGapSbm3DUtilities::ClassifyExternalSpans(
    const std::set<SpanKey3D>& rActiveSpans,
    const std::vector<SurrogateFaceData>& rSurrogateFaces,
    const KnotSpanGridInfo& rGridInfo,
    const ActiveSpanPredicateType& rIsSpanActive) const
{
    ExternalSpanDataMap external_spans;

    for (const auto& r_face_data : rSurrogateFaces) {
        RegisterExternalSpanCandidate(
            external_spans,
            r_face_data.ExternalSpan,
            GapSpanType::Type1,
            r_face_data.ActiveSpan,
            r_face_data.pCondition != nullptr ? r_face_data.pCondition->Id() : 0);
    }

    for (const auto& r_active_span : rActiveSpans) {
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int dk = -1; dk <= 1; ++dk) {
                    if (di == 0 && dj == 0 && dk == 0) {
                        continue;
                    }

                    SpanKey3D candidate_span{
                        r_active_span.I + di,
                        r_active_span.J + dj,
                        r_active_span.K + dk};

                    if (!IsSpanInsideDomain(candidate_span, rGridInfo)) {
                        continue;
                    }

                    if (rActiveSpans.find(candidate_span) != rActiveSpans.end()) {
                        continue;
                    }

                    if (rIsSpanActive && rIsSpanActive(candidate_span)) {
                        continue;
                    }

                    const int manhattan_distance =
                        std::abs(di) + std::abs(dj) + std::abs(dk);

                    const GapSpanType type =
                        manhattan_distance == 1 ? GapSpanType::Type1 :
                        manhattan_distance == 2 ? GapSpanType::Type2 :
                                                   GapSpanType::Type3;

                    RegisterExternalSpanCandidate(
                        external_spans,
                        candidate_span,
                        type,
                        r_active_span,
                        0);
                }
            }
        }
    }

    return external_spans;
}

SnakeGapSbm3DUtilities::NodePointerType
SnakeGapSbm3DUtilities::FindOrCreateProjectionNodeInSpan(
    ModelPart& rSkinSubModelPart,
    const SpanKey3D& rSpan,
    const array_1d<double, 3>& rReferencePoint,
    const KnotSpanSkinBinsCSR& rSkinBins,
    const KnotSpanGridInfo& rGridInfo) const
{
    const std::size_t nnz_index = FindSpanNnzIndex(rSkinBins, rSpan);

    if (nnz_index != static_cast<std::size_t>(-1)) {
        const auto& r_cell_data = rSkinBins.CellDataByNnz[nnz_index];
        NodePointerType p_best_node;
        double best_distance_squared = std::numeric_limits<double>::max();

        for (const auto& p_node : r_cell_data.Nodes) {
            const array_1d<double, 3> delta = p_node->Coordinates() - rReferencePoint;
            const double distance_squared = inner_prod(delta, delta);

            if (distance_squared < best_distance_squared) {
                best_distance_squared = distance_squared;
                p_best_node = p_node;
            }
        }

        if (p_best_node) {
            return p_best_node;
        }
    }

    struct ConditionDistanceData
    {
        ConditionPointerType pCondition;
        double DistanceSquared = std::numeric_limits<double>::max();
    };

    auto create_auxiliary_projection_node = [&](
        const array_1d<double, 3>& rPoint)
    {
        IndexType new_node_id = GetNextAuxiliarySkinNodeId(rSkinSubModelPart);
        while (rSkinSubModelPart.GetRootModelPart().HasNode(new_node_id)) {
            ++new_node_id;
        }

        return rSkinSubModelPart.CreateNewNode(
            new_node_id,
            rPoint[0],
            rPoint[1],
            rPoint[2]);
    };

    auto try_create_projection_node_from_conditions = [&](
        const std::vector<ConditionPointerType>& rConditions,
        const char* pCandidateSource) -> NodePointerType
    {
        std::vector<ConditionDistanceData> condition_candidates;
        condition_candidates.reserve(rConditions.size());

        for (const auto& p_condition : rConditions) {
            const auto& r_geometry = p_condition->GetGeometry();
            if (r_geometry.size() == 0) {
                continue;
            }

            array_1d<double, 3> condition_center = ZeroVector(3);
            for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
                noalias(condition_center) += r_geometry[point_index].Coordinates();
            }
            condition_center /= static_cast<double>(r_geometry.size());

            const array_1d<double, 3> delta =
                condition_center - rReferencePoint;

            ConditionDistanceData candidate;
            candidate.pCondition = p_condition;
            candidate.DistanceSquared = inner_prod(delta, delta);
            condition_candidates.push_back(candidate);
        }

        std::sort(
            condition_candidates.begin(),
            condition_candidates.end(),
            [](const ConditionDistanceData& rA, const ConditionDistanceData& rB) {
                return rA.DistanceSquared < rB.DistanceSquared;
            });

        for (std::size_t condition_index = 0;
             condition_index < condition_candidates.size();
             ++condition_index) {
            const auto& r_candidate = condition_candidates[condition_index];
            array_1d<double, 3> intersection_point = ZeroVector(3);

            if (!FindTriangleBoxIntersectionPoint(
                    r_candidate.pCondition->GetGeometry(),
                    rSpan,
                    rGridInfo,
                    intersection_point)) {
                continue;
            }

            ClampPointInsideSpan(intersection_point, rSpan, rGridInfo);

            auto p_auxiliary_node =
                create_auxiliary_projection_node(intersection_point);

            KRATOS_INFO_IF("SnakeGapSbm3DUtilities", mEchoLevel > 2)
                << "Created auxiliary skin projection node "
                << p_auxiliary_node->Id()
                << " from skin condition "
                << r_candidate.pCondition->Id()
                << " intersecting external span "
                << SpanToString(rSpan)
                << " after testing "
                << condition_index + 1
                << " "
                << pCandidateSource
                << " condition candidates.\n";

            return p_auxiliary_node;
        }

        return nullptr;
    };

    if (nnz_index != static_cast<std::size_t>(-1)) {
        const auto& r_cell_data = rSkinBins.CellDataByNnz[nnz_index];
        if (!r_cell_data.Conditions.empty()) {
            const auto p_projection_node = try_create_projection_node_from_conditions(
                r_cell_data.Conditions,
                "local");

            if (p_projection_node) {
                return p_projection_node;
            }
        }
    }

    std::vector<ConditionPointerType> all_skin_conditions;
    all_skin_conditions.reserve(rSkinSubModelPart.NumberOfConditions());
    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        all_skin_conditions.push_back(rSkinSubModelPart.pGetCondition(r_condition.Id()));
    }

    const auto p_projection_node = try_create_projection_node_from_conditions(
        all_skin_conditions,
        "global");

    if (p_projection_node) {
        return p_projection_node;
    }

    std::vector<ConditionDistanceData> condition_candidates;
    condition_candidates.reserve(all_skin_conditions.size());
    for (const auto& p_condition : all_skin_conditions) {
        const auto& r_geometry = p_condition->GetGeometry();
        if (r_geometry.size() == 0) {
            continue;
        }

        array_1d<double, 3> condition_center = ZeroVector(3);
        for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
            noalias(condition_center) += r_geometry[point_index].Coordinates();
        }
        condition_center /= static_cast<double>(r_geometry.size());

        const array_1d<double, 3> delta = condition_center - rReferencePoint;

        ConditionDistanceData candidate;
        candidate.pCondition = p_condition;
        candidate.DistanceSquared = inner_prod(delta, delta);
        condition_candidates.push_back(candidate);
    }

    std::sort(
        condition_candidates.begin(),
        condition_candidates.end(),
        [](const ConditionDistanceData& rA, const ConditionDistanceData& rB) {
            return rA.DistanceSquared < rB.DistanceSquared;
        });

    NodePointerType p_closest_condition_projection_node;
    IndexType closest_condition_id = 0;
    array_1d<double, 3> closest_skin_point = ZeroVector(3);
    double closest_skin_distance_squared = std::numeric_limits<double>::max();

    for (const auto& r_candidate : condition_candidates) {
        const auto& r_geometry = r_candidate.pCondition->GetGeometry();
        if (r_geometry.PointsNumber() < 3) {
            continue;
        }

        auto update_closest_point = [&](
            const array_1d<double, 3>& rPoint0,
            const array_1d<double, 3>& rPoint1,
            const array_1d<double, 3>& rPoint2)
        {
            const array_1d<double, 3> candidate_point =
                ClosestPointOnTriangle(
                    rReferencePoint,
                    rPoint0,
                    rPoint1,
                    rPoint2);
            const array_1d<double, 3> delta =
                candidate_point - rReferencePoint;
            const double distance_squared = inner_prod(delta, delta);
            if (distance_squared < closest_skin_distance_squared) {
                closest_skin_distance_squared = distance_squared;
                closest_skin_point = candidate_point;
                closest_condition_id = r_candidate.pCondition->Id();
            }
        };

        update_closest_point(
            r_geometry[0].Coordinates(),
            r_geometry[1].Coordinates(),
            r_geometry[2].Coordinates());

        if (r_geometry.PointsNumber() == 4) {
            update_closest_point(
                r_geometry[0].Coordinates(),
                r_geometry[2].Coordinates(),
                r_geometry[3].Coordinates());
        }
    }

    if (std::isfinite(closest_skin_distance_squared)) {
        array_1d<double, 3> closest_span_point = closest_skin_point;
        ClampPointInsideSpan(closest_span_point, rSpan, rGridInfo);

        p_closest_condition_projection_node =
            create_auxiliary_projection_node(closest_span_point);

        KRATOS_WARNING("SnakeGapSbm3DUtilities")
            << "[FindOrCreateProjectionNodeInSpan] No skin condition "
            << "intersects external span " << SpanToString(rSpan)
            << ". Created auxiliary projection node "
            << p_closest_condition_projection_node->Id()
            << " by clamping the closest point on skin condition "
            << closest_condition_id
            << " into the external span box. Distance from span center to "
            << "closest skin point: " << std::sqrt(closest_skin_distance_squared)
            << ", clamp distance: "
            << norm_2(closest_span_point - closest_skin_point) << "\n";

        return p_closest_condition_projection_node;
    }

    std::ostringstream tested_conditions_buffer;
    const std::size_t number_of_preferred_conditions =
        std::min<std::size_t>(2, condition_candidates.size());
    for (std::size_t condition_index = 0;
         condition_index < number_of_preferred_conditions;
         ++condition_index) {
        tested_conditions_buffer
            << " " << condition_candidates[condition_index].pCondition->Id();
    }

    KRATOS_ERROR
        << "[SnakeGapSbm3DUtilities::FindOrCreateProjectionNodeInSpan] "
        << "No skin node exists in external span "
        << SpanToString(rSpan)
        << " and none of the skin conditions intersects the span box.\n"
        << "  span center: " << rReferencePoint << "\n"
        << "  has local span bin: "
        << (nnz_index != static_cast<std::size_t>(-1) ? "yes" : "no") << "\n"
        << "  number of skin conditions: "
        << all_skin_conditions.size() << "\n"
        << "  two closest condition ids:"
        << tested_conditions_buffer.str() << "\n";

    return nullptr;
}

void SnakeGapSbm3DUtilities::AddInteriorSkinNodesForLargeTriangles(
    ModelPart& rSkinSubModelPart,
    const KnotSpanGridInfo& rGridInfo) const
{
    const double minimum_span_size = std::min({
        rGridInfo.SpanSizeX,
        rGridInfo.SpanSizeY,
        rGridInfo.SpanSizeZ});
    const double target_spacing = minimum_span_size / 2.0;

    if (target_spacing <= 0.0) {
        return;
    }

    const double coincident_node_tolerance =
        1.0e-12 * std::max(1.0, minimum_span_size);
    const double coincident_node_tolerance_squared =
        coincident_node_tolerance * coincident_node_tolerance;

    IndexType next_node_id = GetNextAuxiliarySkinNodeId(rSkinSubModelPart);
    std::size_t number_of_refined_triangles = 0;
    std::size_t number_of_added_nodes = 0;

    auto point_exists = [&](
        const array_1d<double, 3>& rPoint)
    {
        for (const auto& r_node : rSkinSubModelPart.Nodes()) {
            const array_1d<double, 3> delta =
                r_node.Coordinates() - rPoint;
            if (inner_prod(delta, delta) <= coincident_node_tolerance_squared) {
                return true;
            }
        }

        return false;
    };

    auto create_node_if_new = [&](
        const array_1d<double, 3>& rPoint)
    {
        if (point_exists(rPoint)) {
            return;
        }

        while (rSkinSubModelPart.GetRootModelPart().HasNode(next_node_id)) {
            ++next_node_id;
        }

        rSkinSubModelPart.CreateNewNode(
            next_node_id++,
            rPoint[0],
            rPoint[1],
            rPoint[2]);
        ++number_of_added_nodes;
    };

    auto add_triangle_centroid_nodes = [&](
        const array_1d<double, 3>& rPoint0,
        const array_1d<double, 3>& rPoint1,
        const array_1d<double, 3>& rPoint2)
    {
        const double edge_length_01 = norm_2(rPoint1 - rPoint0);
        const double edge_length_12 = norm_2(rPoint2 - rPoint1);
        const double edge_length_20 = norm_2(rPoint0 - rPoint2);
        const double maximum_edge_length =
            std::max({edge_length_01, edge_length_12, edge_length_20});

        if (maximum_edge_length <= target_spacing) {
            return;
        }

        const int number_of_subdivisions = std::max(
            2,
            static_cast<int>(std::ceil(maximum_edge_length / target_spacing)));

        auto barycentric_point = [&](
            const double Weight0,
            const double Weight1,
            const double Weight2)
        {
            array_1d<double, 3> point = ZeroVector(3);
            noalias(point) += Weight0 * rPoint0;
            noalias(point) += Weight1 * rPoint1;
            noalias(point) += Weight2 * rPoint2;
            return point;
        };

        for (int i = 0; i < number_of_subdivisions; ++i) {
            for (int j = 0; j < number_of_subdivisions - i; ++j) {
                const int k = number_of_subdivisions - i - j;

                const double lower_weight_0 =
                    (static_cast<double>(i) + 1.0 / 3.0) /
                    static_cast<double>(number_of_subdivisions);
                const double lower_weight_1 =
                    (static_cast<double>(j) + 1.0 / 3.0) /
                    static_cast<double>(number_of_subdivisions);
                const double lower_weight_2 =
                    (static_cast<double>(k) - 2.0 / 3.0) /
                    static_cast<double>(number_of_subdivisions);

                create_node_if_new(barycentric_point(
                    lower_weight_0,
                    lower_weight_1,
                    lower_weight_2));

                if (j < number_of_subdivisions - i - 1) {
                    const double upper_weight_0 =
                        (static_cast<double>(i) + 2.0 / 3.0) /
                        static_cast<double>(number_of_subdivisions);
                    const double upper_weight_1 =
                        (static_cast<double>(j) + 2.0 / 3.0) /
                        static_cast<double>(number_of_subdivisions);
                    const double upper_weight_2 =
                        (static_cast<double>(k) - 4.0 / 3.0) /
                        static_cast<double>(number_of_subdivisions);

                    create_node_if_new(barycentric_point(
                        upper_weight_0,
                        upper_weight_1,
                        upper_weight_2));
                }
            }
        }

        ++number_of_refined_triangles;
    };

    const std::vector<IndexType> condition_ids = [&]() {
        std::vector<IndexType> ids;
        ids.reserve(rSkinSubModelPart.NumberOfConditions());
        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            ids.push_back(r_condition.Id());
        }
        return ids;
    }();

    for (const auto condition_id : condition_ids) {
        const auto& r_geometry =
            rSkinSubModelPart.GetCondition(condition_id).GetGeometry();

        if (r_geometry.PointsNumber() < 3) {
            continue;
        }

        const array_1d<double, 3> point_0 = r_geometry[0].Coordinates();
        const array_1d<double, 3> point_1 = r_geometry[1].Coordinates();
        const array_1d<double, 3> point_2 = r_geometry[2].Coordinates();

        add_triangle_centroid_nodes(point_0, point_1, point_2);

        if (r_geometry.PointsNumber() == 4) {
            const array_1d<double, 3> point_3 =
                r_geometry[3].Coordinates();
            add_triangle_centroid_nodes(point_0, point_2, point_3);
        }
    }

    KRATOS_INFO_IF("SnakeGapSbm3DUtilities", mEchoLevel > 1)
        << "Added " << number_of_added_nodes
        << " interior skin nodes from "
        << number_of_refined_triangles
        << " large skin triangles.\n";
}

void SnakeGapSbm3DUtilities::AssignProjectionNodesToExternalSpans(
    ExternalSpanDataMap& rExternalSpans,
    ModelPart& rSkinSubModelPart,
    const KnotSpanSkinBinsCSR& rSkinBins,
    const KnotSpanGridInfo& rGridInfo) const
{
    for (auto& r_entry : rExternalSpans) {
        auto& r_span_data = r_entry.second;

        if (r_span_data.Type != GapSpanType::Type1 &&
            r_span_data.Type != GapSpanType::Type2 &&
            r_span_data.Type != GapSpanType::Type3) {
            continue;
        }

        const auto center = SpanCenter(r_span_data.Key, rGridInfo);

        const auto p_projection_node = FindOrCreateProjectionNodeInSpan(
            rSkinSubModelPart,
            r_span_data.Key,
            center,
            rSkinBins,
            rGridInfo);

        if (p_projection_node) {
            r_span_data.ProjectionNodeId = p_projection_node->Id();
            r_span_data.pProjectionNode = p_projection_node;
        } else {
            KRATOS_ERROR
                << "[SnakeGapSbm3DUtilities::AssignProjectionNodesToExternalSpans] "
                << "No skin projection node found in external span "
                << SpanToString(r_span_data.Key)
                << " of " << GapSpanTypeToString(r_span_data.Type)
                << " with center at " << center << ".\n";
        }
    }
}

void SnakeGapSbm3DUtilities::EnsureSkinNodesInExternalSpans(
    ModelPart& rSkinSubModelPart,
    const ExternalSpanDataMap& rExternalSpans,
    const KnotSpanGridInfo& rGridInfo) const
{
    for (const auto& r_entry : rExternalSpans) {
        const auto& r_span_data = r_entry.second;

        if (r_span_data.Type != GapSpanType::Type1 &&
            r_span_data.Type != GapSpanType::Type2 &&
            r_span_data.Type != GapSpanType::Type3) {
            continue;
        }

        const auto& r_span = r_span_data.Key;

        if (!IsSpanInsideDomain(r_span, rGridInfo)) {
            continue;
        }

        if (SpanHasSkinNode(rSkinSubModelPart, r_span, rGridInfo)) {
            continue;
        }

        bool found_intersecting_condition = false;
        array_1d<double, 3> auxiliary_point = ZeroVector(3);

        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            if (FindConditionSpanIntersectionPoint(
                    r_condition.GetGeometry(),
                    r_span,
                    rGridInfo,
                    auxiliary_point)) {
                found_intersecting_condition = true;
                break;
            }
        }

        if (!found_intersecting_condition) {
            array_1d<double, 3> box_min = ZeroVector(3);
            array_1d<double, 3> box_max = ZeroVector(3);
            array_1d<double, 3> sphere_center = ZeroVector(3);
            sphere_center[0] = 1.0; sphere_center[1] = 1.0; sphere_center[2] = 1.0; // Dummy non-zero center to avoid warnings in SpanBox.
            SpanBox(r_span, rGridInfo, box_min, box_max);

            KRATOS_WATCH(box_min)
            KRATOS_WATCH(box_max)
            KRATOS_WATCH(norm_2(box_max - sphere_center))
            KRATOS_WATCH(norm_2(box_min - sphere_center))
            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "External span " << SpanToString(r_span)
                << " has no skin node and no intersecting skin condition. "
                << "No auxiliary projection node was created.\n";
            continue;
        }

        ClampPointInsideSpan(auxiliary_point, r_span, rGridInfo);

        const IndexType new_node_id = GetNextAuxiliarySkinNodeId(rSkinSubModelPart);

        rSkinSubModelPart.CreateNewNode(
            new_node_id,
            auxiliary_point[0],
            auxiliary_point[1],
            auxiliary_point[2]);

        KRATOS_INFO_IF("SnakeGapSbm3DUtilities", mEchoLevel > 3)
            << "Added auxiliary skin node " << new_node_id
            << " in external span " << SpanToString(r_span)
            << " at " << auxiliary_point << "\n";
    }
}

SnakeGapSbm3DUtilities::ExternalSpanDataMap
SnakeGapSbm3DUtilities::InitializeExternalSpanData(
    ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart,
    const ActiveSpanPredicateType& rIsSpanActive) const
{
    const auto& r_background_model_part = rSurrogateSubModelPart.GetParentModelPart();

    const auto grid_info = CreateKnotSpanGridInfo(r_background_model_part);

    const auto surrogate_faces = BuildSurrogateFaceDataVector(
        rSurrogateSubModelPart,
        grid_info);

    const auto active_spans = ExtractActiveSpans(surrogate_faces);

    auto external_spans = ClassifyExternalSpans(
        active_spans,
        surrogate_faces,
        grid_info,
        rIsSpanActive);

    // AddInteriorSkinNodesForLargeTriangles( //FIXME:
    //     rSkinSubModelPart,
    //     grid_info);

    // Build bins after additional skin nodes have been created.
    const auto skin_bins = BuildSkinBinsPerKnotSpan3D(
        rSkinSubModelPart,
        grid_info);

    AssignProjectionNodesToExternalSpans(
        external_spans,
        rSkinSubModelPart,
        skin_bins,
        grid_info);

    PrintExternalSpanSummary(external_spans);

    return external_spans;
}

std::string SnakeGapSbm3DUtilities::SpanToString(
    const SpanKey3D& rSpan) const
{
    std::ostringstream buffer;
    buffer << "(" << rSpan.I << "," << rSpan.J << "," << rSpan.K << ")";
    return buffer.str();
}

void SnakeGapSbm3DUtilities::PrintExternalSpanSummary(
    const ExternalSpanDataMap& rExternalSpans) const
{
    if (mEchoLevel == 0) {
        return;
    }

    std::size_t number_type_1 = 0;
    std::size_t number_type_2 = 0;
    std::size_t number_type_3 = 0;
    std::size_t number_type_1_with_projection = 0;
    std::size_t number_type_2_with_projection = 0;
    std::size_t number_type_3_with_projection = 0;
    std::size_t number_missing_projection = 0;

    for (const auto& r_entry : rExternalSpans) {
        const auto& r_data = r_entry.second;

        if (r_data.Type == GapSpanType::Type1) {
            ++number_type_1;
            if (r_data.HasProjectionNode()) {
                ++number_type_1_with_projection;
            } else {
                ++number_missing_projection;
            }
        } else if (r_data.Type == GapSpanType::Type2) {
            ++number_type_2;
            if (r_data.HasProjectionNode()) {
                ++number_type_2_with_projection;
            } else {
                ++number_missing_projection;
            }
        } else if (r_data.Type == GapSpanType::Type3) {
            ++number_type_3;
            if (r_data.HasProjectionNode()) {
                ++number_type_3_with_projection;
            } else {
                ++number_missing_projection;
            }
        }
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "External span initialization summary:\n"
        << "  type 1 spans:                 " << number_type_1 << "\n"
        << "  type 2 spans:                 " << number_type_2 << "\n"
        << "  type 3 spans:                 " << number_type_3 << "\n"
        << "  type 1 with projection node:  " << number_type_1_with_projection << "\n"
        << "  type 2 with projection node:  " << number_type_2_with_projection << "\n"
        << "  type 3 with projection node:  " << number_type_3_with_projection << "\n"
        << "  missing projection nodes:     " << number_missing_projection << "\n";
}


void SnakeGapSbmProcess::CreateSbmExtendedGeometries3D()
{
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 1)
        << "Echo level: " << mEchoLevel
        << " | Gap type: " << mGapSbmType
        << " | Approx order: " << mGapApproximationOrder
        << " | Internal divisions: " << mInternalDivisions
        << " | Rel tol subdivisions: " << mGapRelativeToleranceForSubdivisions
        << " | Interp levels: " << mNumberOfInterpolationLevels << std::endl;
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0 || mpSkinModelPartInnerInitial->NumberOfGeometries()>0) 
    {
        const auto& r_surrogate_sub_model_part_inner = mpIgaModelPart->GetSubModelPart("surrogate_inner");
        auto& r_skin_sub_model_part_inner = mpSkinModelPart->GetSubModelPart("inner");

        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the inner skin." << std::endl;
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 2)
            << "Inner skin nodes: " << r_skin_sub_model_part_inner.NumberOfNodes()
            << " | conditions: " << r_skin_sub_model_part_inner.NumberOfConditions()
            << " | surrogate conditions: " << r_surrogate_sub_model_part_inner.NumberOfConditions() << std::endl;
        CreateSbmExtendedGeometries3D<true>(r_skin_sub_model_part_inner, r_surrogate_sub_model_part_inner);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the inner skin." << std::endl;
    }
    if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) 
    {
        const auto& r_surrogate_sub_model_part_outer = mpIgaModelPart->GetSubModelPart("surrogate_outer");
        auto& r_skin_sub_model_part_outer = mpSkinModelPart->GetSubModelPart("outer");
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the outer skin." << std::endl;
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 2)
            << "Outer skin nodes: " << r_skin_sub_model_part_outer.NumberOfNodes()
            << " | conditions: " << r_skin_sub_model_part_outer.NumberOfConditions()
            << " | surrogate conditions: " << r_surrogate_sub_model_part_outer.NumberOfConditions() << std::endl;
        CreateSbmExtendedGeometries3D<false>(r_skin_sub_model_part_outer, r_surrogate_sub_model_part_outer);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the outer skin." << std::endl;
    }
}

// type 1 element methods
Vector SnakeGapSbm3DUtilities::CreateOpenUnitKnotVectorDegree1() const
{
    return CreateOpenUnitKnotVector(1);
}

Vector SnakeGapSbm3DUtilities::CreateOpenUnitKnotVector(
    const std::size_t PolynomialDegree) const
{
    KRATOS_ERROR_IF(PolynomialDegree == 0)
        << "[SnakeGapSbm3DUtilities::CreateOpenUnitKnotVector] "
        << "Polynomial degree must be positive.\n";

    Vector knot_vector = ZeroVector(2 * (PolynomialDegree + 1));
    for (std::size_t i = PolynomialDegree + 1;
         i < 2 * (PolynomialDegree + 1);
         ++i) {
        knot_vector[i] = 1.0;
    }

    return knot_vector;
}

SnakeGapSbm3DUtilities::CanonicalFaceKey3D
SnakeGapSbm3DUtilities::MakeCanonicalFaceKey3D(
    const IndexType NodeId0,
    const IndexType NodeId1,
    const IndexType NodeId2) const
{
    CanonicalFaceKey3D key;
    key.NodeIds = {{NodeId0, NodeId1, NodeId2}};
    std::sort(key.NodeIds.begin(), key.NodeIds.end());
    return key;
}

SnakeGapSbm3DUtilities::CanonicalFaceKey3D
SnakeGapSbm3DUtilities::MakeCanonicalFaceKey3D(
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pNode2) const
{
    KRATOS_ERROR_IF_NOT(pNode0)
        << "[SnakeGapSbm3DUtilities::MakeCanonicalFaceKey3D] pNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pNode1)
        << "[SnakeGapSbm3DUtilities::MakeCanonicalFaceKey3D] pNode1 is null.\n";
    KRATOS_ERROR_IF_NOT(pNode2)
        << "[SnakeGapSbm3DUtilities::MakeCanonicalFaceKey3D] pNode2 is null.\n";

    std::array<const NodeType*, 3> node_pointers = {{
        pNode0.get(),
        pNode1.get(),
        pNode2.get()}};

    std::sort(
        node_pointers.begin(),
        node_pointers.end(),
        [](const NodeType* pNodeA, const NodeType* pNodeB) {
            return std::less<const NodeType*>{}(pNodeA, pNodeB);
        });

    CanonicalFaceKey3D key;
    key.NodePointers = node_pointers;
    key.NodeIds = {{
        node_pointers[0]->Id(),
        node_pointers[1]->Id(),
        node_pointers[2]->Id()}};
    return key;
}

IndexType SnakeGapSbm3DUtilities::GetNextGeometryId(
    const ModelPart& rRootModelPart) const
{
    IndexType max_id = 0;

    for (const auto& r_geometry : rRootModelPart.Geometries()) {
        max_id = std::max(max_id, r_geometry.Id());
    }

    return max_id + 1;
}

ModelPart& SnakeGapSbm3DUtilities::GetOrCreateSubModelPart(
    ModelPart& rRootModelPart,
    const std::string& rSubModelPartName) const
{
    if (rRootModelPart.HasSubModelPart(rSubModelPartName)) {
        return rRootModelPart.GetSubModelPart(rSubModelPartName);
    }

    return rRootModelPart.CreateSubModelPart(rSubModelPartName);
}

void SnakeGapSbm3DUtilities::ClearLateralFaceRegistry()
{
    mLateralFaceRegistry.clear();
}

const SnakeGapSbm3DUtilities::LateralFaceRegistry&
SnakeGapSbm3DUtilities::GetLateralFaceRegistry() const
{
    return mLateralFaceRegistry;
}

std::size_t SnakeGapSbm3DUtilities::ComputeNumberOfShapeFunctionsDerivatives(
    const ModelPart& rIgaModelPart) const
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::ComputeNumberOfShapeFunctionsDerivatives";

    KRATOS_ERROR_IF_NOT(rIgaModelPart.HasGeometry(1))
        << "[" << p_caller << "] IgaModelPart does not contain geometry #1.\n";

    auto p_volume = rIgaModelPart.pGetGeometry(1);

    auto p_nurbs_volume = std::dynamic_pointer_cast<NurbsVolumeType>(
        p_volume->pGetGeometryPart(
            Geometry<typename PointerVector<Node>::value_type>::BACKGROUND_GEOMETRY_INDEX));

    KRATOS_ERROR_IF_NOT(p_nurbs_volume)
        << "[" << p_caller << "] geometry #1 does not expose a NurbsVolumeType "
        << "as BACKGROUND_GEOMETRY.\n";

    const std::size_t brep_degree = p_nurbs_volume->PolynomialDegree(0);

    return 3 * brep_degree + 1;
}

bool SnakeGapSbm3DUtilities::IsInnerSurrogateLoop(
    const ModelPart& rSurrogateSubModelPart) const
{
    return rSurrogateSubModelPart.Name().find("inner") != std::string::npos;
}

IndexType SnakeGapSbm3DUtilities::ComputeStartingBrepId(
    const ModelPart& rIgaModelPart,
    const ModelPart& rSurrogateSubModelPart) const
{
    if (!IsInnerSurrogateLoop(rSurrogateSubModelPart)) {
        return 2;
    }

    if (rIgaModelPart.HasSubModelPart("surrogate_outer") &&
        rIgaModelPart.GetSubModelPart("surrogate_outer").NumberOfConditions() > 0) {
        return 2 + rIgaModelPart.GetSubModelPart("surrogate_outer").NumberOfConditions();
    }

    return 8;
}

std::size_t SnakeGapSbm3DUtilities::ComputeBrepPatchLoopSize(
    const ModelPart& rSurrogateSubModelPart) const
{
    if (!IsInnerSurrogateLoop(rSurrogateSubModelPart)) {
        return rSurrogateSubModelPart.NumberOfConditions();
    }

    if (rSurrogateSubModelPart.NumberOfElements() == 0) {
        return rSurrogateSubModelPart.NumberOfConditions();
    }

    const IndexType element_id = rSurrogateSubModelPart.ElementsBegin()->Id();

    const IndexType first_condition_id =
        rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[0].Id();

    const IndexType last_condition_id =
        rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[1].Id();

    return last_condition_id - first_condition_id + 1;
}

SnakeGapSbm3DUtilities::BrepPatchData
SnakeGapSbm3DUtilities::ComputeBrepPatchData(
    const ModelPart& rIgaModelPart,
    const IndexType BrepId) const
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::ComputeBrepPatchData";

    BrepPatchData brep_patch_data;

    KRATOS_ERROR_IF_NOT(rIgaModelPart.HasGeometry(BrepId))
        << "[" << p_caller << "] geometry #" << BrepId
        << " was not found in model part '" << rIgaModelPart.Name() << "'.\n";

    brep_patch_data.pBrepGeometry = rIgaModelPart.pGetGeometry(BrepId);

    brep_patch_data.pBrepSurface =
        std::dynamic_pointer_cast<BrepSurfaceOnVolumeType>(
            brep_patch_data.pBrepGeometry);

    KRATOS_ERROR_IF_NOT(brep_patch_data.pBrepSurface)
        << "[" << p_caller << "] geometry #"
        << brep_patch_data.pBrepGeometry->Id()
        << " is not a BrepSurfaceOnVolumeType.\n";

    const NurbsInterval brep_domain_interval_u =
        brep_patch_data.pBrepSurface->DomainIntervalU();

    const NurbsInterval brep_domain_interval_v =
        brep_patch_data.pBrepSurface->DomainIntervalV();

    array_1d<double, 3> vertex_00_local_coords = ZeroVector(3);
    array_1d<double, 3> vertex_01_local_coords = ZeroVector(3);
    array_1d<double, 3> vertex_10_local_coords = ZeroVector(3);
    array_1d<double, 3> vertex_11_local_coords = ZeroVector(3);

    vertex_00_local_coords[0] = brep_domain_interval_u.GetT0();
    vertex_00_local_coords[1] = brep_domain_interval_v.GetT0();

    vertex_01_local_coords[0] = brep_domain_interval_u.GetT0();
    vertex_01_local_coords[1] = brep_domain_interval_v.GetT1();

    vertex_10_local_coords[0] = brep_domain_interval_u.GetT1();
    vertex_10_local_coords[1] = brep_domain_interval_v.GetT0();

    vertex_11_local_coords[0] = brep_domain_interval_u.GetT1();
    vertex_11_local_coords[1] = brep_domain_interval_v.GetT1();

    brep_patch_data.pBrepSurface->GlobalCoordinates(
        brep_patch_data.Vertex00,
        vertex_00_local_coords);

    brep_patch_data.pBrepSurface->GlobalCoordinates(
        brep_patch_data.Vertex01,
        vertex_01_local_coords);

    brep_patch_data.pBrepSurface->GlobalCoordinates(
        brep_patch_data.Vertex10,
        vertex_10_local_coords);

    brep_patch_data.pBrepSurface->GlobalCoordinates(
        brep_patch_data.Vertex11,
        vertex_11_local_coords);

    brep_patch_data.MiddlePointLocalCoordinates[0] =
        0.5 * (brep_domain_interval_u.GetT0() + brep_domain_interval_u.GetT1());

    brep_patch_data.MiddlePointLocalCoordinates[1] =
        0.5 * (brep_domain_interval_v.GetT0() + brep_domain_interval_v.GetT1());

    brep_patch_data.pBrepSurface->GlobalCoordinates(
        brep_patch_data.MiddlePoint,
        brep_patch_data.MiddlePointLocalCoordinates);

    return brep_patch_data;
}

std::vector<SnakeGapSbm3DUtilities::BrepPatchData>
SnakeGapSbm3DUtilities::BuildBrepPatchDataVector(
    const ModelPart& rIgaModelPart,
    const ModelPart& rSurrogateSubModelPart) const
{
    const IndexType starting_brep_id =
        ComputeStartingBrepId(rIgaModelPart, rSurrogateSubModelPart);

    const std::size_t number_of_brep_patches =
        ComputeBrepPatchLoopSize(rSurrogateSubModelPart);

    std::vector<BrepPatchData> brep_patch_data_list;
    brep_patch_data_list.reserve(number_of_brep_patches);

    for (std::size_t i = 0; i < number_of_brep_patches; ++i) {
        brep_patch_data_list.push_back(
            ComputeBrepPatchData(
                rIgaModelPart,
                starting_brep_id + static_cast<IndexType>(i)));
    }

    return brep_patch_data_list;
}

const SnakeGapSbm3DUtilities::BrepPatchData&
SnakeGapSbm3DUtilities::FindBrepPatchMatchingCondition(
    const Condition& rCondition,
    const std::vector<BrepPatchData>& rBrepPatchDataList) const
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::FindBrepPatchMatchingCondition";

    KRATOS_ERROR_IF(rBrepPatchDataList.empty())
        << "[" << p_caller << "] Empty BREP patch data list.\n";

    constexpr double distance_tolerance = 1.0e-8;

    const auto& r_condition_geometry = rCondition.GetGeometry();
    const array_1d<double, 3> condition_center = r_condition_geometry.Center();

    for (const auto& r_brep_patch_data : rBrepPatchDataList) {
        const array_1d<double, 3> center_difference =
            condition_center - r_brep_patch_data.MiddlePoint;

        const double distance_squared =
            inner_prod(center_difference, center_difference);

        if (distance_squared <= distance_tolerance * distance_tolerance) {
            return r_brep_patch_data;
        }
    }

    KRATOS_ERROR
        << "[" << p_caller << "] Failed to find BREP patch matching surrogate condition #"
        << rCondition.Id()
        << " center=" << condition_center << ".\n";

    return rBrepPatchDataList.front();
}

Geometry<Node>::Pointer
SnakeGapSbm3DUtilities::CreateSurrogateFaceNeighbourGeometry(
    const BrepPatchData& rBrepPatchData,
    const std::size_t NumberOfShapeFunctionsDerivatives) const
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::CreateSurrogateFaceNeighbourGeometry";

    KRATOS_ERROR_IF_NOT(rBrepPatchData.pBrepSurface)
        << "[" << p_caller << "] Null BREP surface.\n";

    IntegrationPoint<2> integration_point(
        rBrepPatchData.MiddlePointLocalCoordinates[0],
        rBrepPatchData.MiddlePointLocalCoordinates[1],
        1.0);

    IntegrationPointsArrayType surrogate_integration_points;
    surrogate_integration_points.push_back(integration_point);

    IntegrationInfo integration_info =
        rBrepPatchData.pBrepSurface->GetDefaultIntegrationInfo();

    GeometriesArrayType quadrature_point_list;

    rBrepPatchData.pBrepSurface->CreateQuadraturePointGeometries(
        quadrature_point_list,
        NumberOfShapeFunctionsDerivatives,
        surrogate_integration_points,
        integration_info);

    KRATOS_ERROR_IF(quadrature_point_list.size() == 0)
        << "[" << p_caller << "] Failed to create face centre NEIGHBOUR_GEOMETRY "
        << "from BREP surface #" << rBrepPatchData.pBrepSurface->Id() << ".\n";

    return quadrature_point_list(0);
}

SnakeGapSbm3DUtilities::NurbsSurfaceType::Pointer
SnakeGapSbm3DUtilities::CreateType1CollapsedLateralCoonsSurface(
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pApexNode) const
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::CreateType1CollapsedLateralCoonsSurface";

    KRATOS_ERROR_IF_NOT(pNode0)
        << "[" << p_caller << "] First base node is null.\n";
    KRATOS_ERROR_IF_NOT(pNode1)
        << "[" << p_caller << "] Second base node is null.\n";
    KRATOS_ERROR_IF_NOT(pApexNode)
        << "[" << p_caller << "] Apex node is null.\n";

    PointerVector<NodeType> control_points;

    // Tensor-product order:
    // P00 = S0, P10 = S1, P01 = P, P11 = P.
    control_points.push_back(pNode0);
    control_points.push_back(pNode1);
    control_points.push_back(pApexNode);
    control_points.push_back(pApexNode);

    return Kratos::make_shared<NurbsSurfaceType>(
        control_points,
        std::size_t(1),
        std::size_t(1),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1());
}

SnakeGapSbm3DUtilities::NurbsVolumeType::Pointer
SnakeGapSbm3DUtilities::CreateType1CollapsedPyramidVolume(
    const SurrogateFaceData& rFaceData,
    const NodePointerType& pApexNode) const
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::CreateType1CollapsedPyramidVolume";

    KRATOS_ERROR_IF_NOT(pApexNode)
        << "[" << p_caller << "] Apex node is null.\n";

    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_ERROR_IF_NOT(rFaceData.Nodes[i])
            << "[" << p_caller << "] Null base node at local index " << i
            << " for surrogate condition #" << rFaceData.pCondition->Id() << ".\n";
    }

    PointerVector<NodeType> control_points;

    // Bottom surrogate face in tensor-product order.
    control_points.push_back(rFaceData.Nodes[0]); // P000
    control_points.push_back(rFaceData.Nodes[1]); // P100
    control_points.push_back(rFaceData.Nodes[3]); // P010
    control_points.push_back(rFaceData.Nodes[2]); // P110

    // Collapsed skin-side top face.
    control_points.push_back(pApexNode); // P001
    control_points.push_back(pApexNode); // P101
    control_points.push_back(pApexNode); // P011
    control_points.push_back(pApexNode); // P111

    return Kratos::make_shared<NurbsVolumeType>(
        control_points,
        std::size_t(1),
        std::size_t(1),
        std::size_t(1),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1());
}

SnakeGapSbm3DUtilities::NurbsVolumeType::Pointer
SnakeGapSbm3DUtilities::CreateType2CollapsedEdgeVolume(
    const NodePointerType& pEdgeNode0,
    const NodePointerType& pEdgeNode1,
    const std::vector<NodePointerType>& rSkinEdgeControlNodes) const
{
    KRATOS_ERROR_IF_NOT(pEdgeNode0)
        << "[CreateType2CollapsedEdgeVolume] pEdgeNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pEdgeNode1)
        << "[CreateType2CollapsedEdgeVolume] pEdgeNode1 is null.\n";
    KRATOS_ERROR_IF(
        rSkinEdgeControlNodes.size() != mGapApproximationOrder + 1)
        << "[CreateType2CollapsedEdgeVolume] Unexpected number of skin edge "
        << "control nodes: " << rSkinEdgeControlNodes.size()
        << ". Expected " << mGapApproximationOrder + 1 << ".\n";

    PointerVector<NodeType> control_points;

    // Tensor-product order: u along the surrogate edge, v along the
    // skin-skin edge, w from surrogate to skin. The surrogate edge is
    // duplicated for every v-control and the skin curve is collapsed in u.
    for (std::size_t v = 0; v <= mGapApproximationOrder; ++v) {
        control_points.push_back(pEdgeNode0);
        control_points.push_back(pEdgeNode1);
    }

    for (const auto& p_skin_control_node : rSkinEdgeControlNodes) {
        KRATOS_ERROR_IF_NOT(p_skin_control_node)
            << "[CreateType2CollapsedEdgeVolume] Null skin edge control node.\n";
        control_points.push_back(p_skin_control_node);
        control_points.push_back(p_skin_control_node);
    }

    return Kratos::make_shared<NurbsVolumeType>(
        control_points,
        1,
        mGapApproximationOrder,
        1,
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVector(mGapApproximationOrder),
        CreateOpenUnitKnotVectorDegree1());
}

SnakeGapSbm3DUtilities::NurbsVolumeType::Pointer
SnakeGapSbm3DUtilities::CreateType2LinearCollapsedEdgeVolume(
    const NodePointerType& pEdgeNode0,
    const NodePointerType& pEdgeNode1,
    const NodePointerType& pSkinNode0,
    const NodePointerType& pSkinNode1) const
{
    KRATOS_ERROR_IF_NOT(pEdgeNode0)
        << "[CreateType2LinearCollapsedEdgeVolume] pEdgeNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pEdgeNode1)
        << "[CreateType2LinearCollapsedEdgeVolume] pEdgeNode1 is null.\n";
    KRATOS_ERROR_IF_NOT(pSkinNode0)
        << "[CreateType2LinearCollapsedEdgeVolume] pSkinNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pSkinNode1)
        << "[CreateType2LinearCollapsedEdgeVolume] pSkinNode1 is null.\n";

    PointerVector<NodeType> control_points;

    // Tensor-product order consistent with CreateType2CollapsedEdgeVolume:
    //
    // u: along the surrogate edge / collapsed skin edge
    // v: along the skin-skin edge
    // w: from surrogate to skin
    //
    // w = 0: surrogate side, duplicated for the two v-control positions
    control_points.push_back(pEdgeNode0); // u0, v0, w0
    control_points.push_back(pEdgeNode1); // u1, v0, w0
    control_points.push_back(pEdgeNode0); // u0, v1, w0
    control_points.push_back(pEdgeNode1); // u1, v1, w0

    // w = 1: skin side, collapsed in u
    control_points.push_back(pSkinNode0); // u0, v0, w1
    control_points.push_back(pSkinNode0); // u1, v0, w1
    control_points.push_back(pSkinNode1); // u0, v1, w1
    control_points.push_back(pSkinNode1); // u1, v1, w1

    return Kratos::make_shared<NurbsVolumeType>(
        control_points,
        std::size_t(1),
        std::size_t(1),
        std::size_t(1),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1());
}

SnakeGapSbm3DUtilities::Type3CollapsedCornerData
SnakeGapSbm3DUtilities::CreateType3CollapsedCornerVolume(
    const ModelPart& rSkinSubModelPart,
    const double MaximumProjectionDistance,
    const NodePointerType& pSurrogateNode,
    const NodePointerType& pProjectionNode0,
    const NodePointerType& pProjectionNode1,
    const NodePointerType& pProjectionNode2) const
{
    KRATOS_ERROR_IF_NOT(pSurrogateNode)
        << "[CreateType3CollapsedCornerVolume] pSurrogateNode is null.\n";
    KRATOS_ERROR_IF_NOT(pProjectionNode0)
        << "[CreateType3CollapsedCornerVolume] pProjectionNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pProjectionNode1)
        << "[CreateType3CollapsedCornerVolume] pProjectionNode1 is null.\n";
    KRATOS_ERROR_IF_NOT(pProjectionNode2)
        << "[CreateType3CollapsedCornerVolume] pProjectionNode2 is null.\n";
    KRATOS_ERROR_IF(MaximumProjectionDistance <= 0.0)
        << "[CreateType3CollapsedCornerVolume] Maximum projection distance "
        << "must be positive. Got " << MaximumProjectionDistance << ".\n";

    if (mGapApproximationOrder == 1) {
        Type3CollapsedCornerData result;
        result.ProjectionNodes = {{
            pProjectionNode0,
            pProjectionNode1,
            pProjectionNode2}};

        PointerVector<NodeType> control_points;
        control_points.push_back(pSurrogateNode);
        control_points.push_back(pProjectionNode0);
        control_points.push_back(pProjectionNode1);
        control_points.push_back(pProjectionNode1);

        control_points.push_back(pProjectionNode2);
        control_points.push_back(pProjectionNode2);
        control_points.push_back(pProjectionNode2);
        control_points.push_back(pProjectionNode2);

        result.pVolume = Kratos::make_shared<NurbsVolumeType>(
            control_points,
            1,
            1,
            1,
            CreateOpenUnitKnotVectorDegree1(),
            CreateOpenUnitKnotVectorDegree1(),
            CreateOpenUnitKnotVectorDegree1());

        result.pTopSurface = CreateCollapsedTriangleSurface(
            pProjectionNode0,
            pProjectionNode1,
            pProjectionNode2);

        return result;
    }

    KRATOS_ERROR_IF(mGapApproximationOrder != 2)
        << "[CreateType3CollapsedCornerVolume] Type3 curved corner volumes "
        << "require gap_approximation_order == 2. Requested order: "
        << mGapApproximationOrder << ".\n";

    struct Type3JacobianQuality
    {
        bool AllFinite = true;
        bool StrictlyPositive = true;
        double MinDeterminant = std::numeric_limits<double>::max();
        double MaxDeterminant = -std::numeric_limits<double>::max();
    };

    auto evaluate_jacobian_quality = [](
        const NurbsVolumeType& rVolume)
    {
        IntegrationPointsArrayType integration_points(27);
        auto integration_point_it = integration_points.begin();

        IntegrationPointUtilities::IntegrationPoints3D(
            integration_point_it,
            3,
            3,
            3,
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            1.0);

        constexpr double determinant_tolerance = 1.0e-14;
        Type3JacobianQuality quality;

        for (const auto& r_integration_point : integration_points) {
            CoordinatesArrayType local_coordinates = ZeroVector(3);
            local_coordinates[0] = r_integration_point[0];
            local_coordinates[1] = r_integration_point[1];
            local_coordinates[2] = r_integration_point[2];

            Matrix jacobian;
            rVolume.Jacobian(jacobian, local_coordinates);

            const double determinant = MathUtils<double>::Det(jacobian);
            if (!std::isfinite(determinant)) {
                quality.AllFinite = false;
                quality.StrictlyPositive = false;
                continue;
            }

            quality.MinDeterminant =
                std::min(quality.MinDeterminant, determinant);
            quality.MaxDeterminant =
                std::max(quality.MaxDeterminant, determinant);

            if (determinant <= determinant_tolerance) {
                quality.StrictlyPositive = false;
            }
        }

        if (quality.MinDeterminant == std::numeric_limits<double>::max()) {
            quality.AllFinite = false;
            quality.StrictlyPositive = false;
            quality.MinDeterminant = -std::numeric_limits<double>::max();
            quality.MaxDeterminant = -std::numeric_limits<double>::max();
        }

        return quality;
    };

    auto create_top_surface = [this](
        const std::vector<NodePointerType>& rTopControlPoints)
    {
        KRATOS_ERROR_IF(rTopControlPoints.size() != 9)
            << "[CreateType3CollapsedCornerVolume] Type3 top face requires "
            << "nine quadratic control points. Got "
            << rTopControlPoints.size() << ".\n";

        PointerVector<NodeType> surface_control_points;
        for (const auto& p_control_point : rTopControlPoints) {
            KRATOS_ERROR_IF_NOT(p_control_point)
                << "[CreateType3CollapsedCornerVolume] Null top face "
                << "control point.\n";
            surface_control_points.push_back(p_control_point);
        }

        return Kratos::make_shared<NurbsSurfaceType>(
            surface_control_points,
            std::size_t(2),
            std::size_t(2),
            CreateOpenUnitKnotVector(2),
            CreateOpenUnitKnotVector(2));
    };

    auto create_curved_type3_volume = [this,
                                      &rSkinSubModelPart,
                                      &pSurrogateNode,
                                      &MaximumProjectionDistance,
                                      &create_top_surface](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const bool ReverseRadialDirection)
        -> Type3CollapsedCornerData
    {
        Type3CollapsedCornerData result;
        result.ProjectionNodes = {{pNode0, pNode1, pNode2}};

        auto get_oriented_edge = [this, &pSurrogateNode](
            const NodePointerType& pFirst,
            const NodePointerType& pSecond)
        {
            const auto p_cached_control_nodes =
                FindCachedSkinEdgeControlNodes(pFirst, pSecond);
            KRATOS_ERROR_IF_NOT(p_cached_control_nodes)
                << "[CreateType3CollapsedCornerVolume] Missing curved "
                << "skin-skin edge required by type3 top face.\n"
                << "  first projection node: " << pFirst->Id()
                << " coords: " << pFirst->Coordinates() << "\n"
                << "  second projection node: " << pSecond->Id()
                << " coords: " << pSecond->Coordinates() << "\n"
                << "  surrogate node: " << pSurrogateNode->Id() << "\n";

            auto control_nodes = *p_cached_control_nodes;
            if (control_nodes.front()->Id() != pFirst->Id()) {
                std::reverse(control_nodes.begin(), control_nodes.end());
            }
            return control_nodes;
        };

        const auto edge_01 = get_oriented_edge(pNode0, pNode1);
        const auto edge_12 = get_oriented_edge(pNode1, pNode2);
        const auto edge_02 = get_oriented_edge(pNode0, pNode2);

        array_1d<double, 3> coons_middle = ZeroVector(3);
        noalias(coons_middle) += 0.25 * edge_01[1]->Coordinates();
        noalias(coons_middle) += 0.25 * edge_02[1]->Coordinates();
        noalias(coons_middle) += 0.25 * edge_12[1]->Coordinates();
        noalias(coons_middle) += 0.25 * pNode2->Coordinates();

        array_1d<double, 3> middle_point = ZeroVector(3);
        double middle_projection_distance = std::numeric_limits<double>::max();
        const bool projected_middle_point =
            ProjectPointToClosestSkinBoundary(
                rSkinSubModelPart,
                coons_middle,
                middle_point,
                middle_projection_distance);

        KRATOS_ERROR_IF_NOT(projected_middle_point)
            << "[CreateType3CollapsedCornerVolume] Could not project type3 "
            << "top-face middle control point to the skin.\n"
            << "  surrogate node: " << pSurrogateNode->Id() << "\n"
            << "  projection nodes: "
            << pNode0->Id() << ", "
            << pNode1->Id() << ", "
            << pNode2->Id() << "\n"
            << "  initial middle point: " << coons_middle << "\n";

        KRATOS_ERROR_IF(middle_projection_distance > MaximumProjectionDistance)
            << "[CreateType3CollapsedCornerVolume] Type3 top-face middle "
            << "projection is farther than the allowed distance.\n"
            << "  surrogate node: " << pSurrogateNode->Id() << "\n"
            << "  projection nodes: "
            << pNode0->Id() << ", "
            << pNode1->Id() << ", "
            << pNode2->Id() << "\n"
            << "  initial middle point: " << coons_middle << "\n"
            << "  projected middle point: " << middle_point << "\n"
            << "  projection distance: " << middle_projection_distance << "\n"
            << "  maximum allowed distance: "
            << MaximumProjectionDistance << "\n";

        array_1d<double, 3> middle_control_point =
            4.0 * middle_point;
        noalias(middle_control_point) -= 0.25 * edge_01[0]->Coordinates();
        noalias(middle_control_point) -= 0.50 * edge_01[1]->Coordinates();
        noalias(middle_control_point) -= 0.25 * edge_01[2]->Coordinates();
        noalias(middle_control_point) -= 0.50 * edge_02[1]->Coordinates();
        noalias(middle_control_point) -= 0.50 * edge_12[1]->Coordinates();
        noalias(middle_control_point) -= pNode2->Coordinates();

        NodePointerType p_middle_node(new Node(0, middle_control_point));

        std::vector<NodePointerType> top_control_points;
        top_control_points.push_back(edge_01[0]);
        top_control_points.push_back(edge_01[1]);
        top_control_points.push_back(edge_01[2]);

        top_control_points.push_back(edge_02[1]);
        top_control_points.push_back(p_middle_node);
        top_control_points.push_back(edge_12[1]);

        top_control_points.push_back(pNode2);
        top_control_points.push_back(pNode2);
        top_control_points.push_back(pNode2);

        result.pTopSurface = create_top_surface(top_control_points);

        PointerVector<NodeType> control_points;
        for (const auto& r_top_control_point : top_control_points) {
            if (ReverseRadialDirection) {
                control_points.push_back(r_top_control_point);
                control_points.push_back(pSurrogateNode);
            } else {
                control_points.push_back(pSurrogateNode);
                control_points.push_back(r_top_control_point);
            }
        }

        result.pVolume = Kratos::make_shared<NurbsVolumeType>(
            control_points,
            std::size_t(1),
            std::size_t(2),
            std::size_t(2),
            CreateOpenUnitKnotVectorDegree1(),
            CreateOpenUnitKnotVector(2),
            CreateOpenUnitKnotVector(2));

        return result;
    };

    const std::array<std::array<NodePointerType, 3>, 6> node_permutations = {{
        {{pProjectionNode0, pProjectionNode1, pProjectionNode2}},
        {{pProjectionNode0, pProjectionNode2, pProjectionNode1}},
        {{pProjectionNode1, pProjectionNode0, pProjectionNode2}},
        {{pProjectionNode1, pProjectionNode2, pProjectionNode0}},
        {{pProjectionNode2, pProjectionNode0, pProjectionNode1}},
        {{pProjectionNode2, pProjectionNode1, pProjectionNode0}}
    }};

    auto select_valid_curved_volume = [&]()
        -> Type3CollapsedCornerData
    {
        Type3CollapsedCornerData best_finite_candidate;
        Type3JacobianQuality best_finite_quality;
        bool has_best_finite_candidate = false;

        auto consider_candidate = [&](
            const Type3CollapsedCornerData& rCandidateData)
        {
            if (!rCandidateData.pVolume) {
                return false;
            }

            const Type3JacobianQuality quality =
                evaluate_jacobian_quality(*rCandidateData.pVolume);

            if (!quality.AllFinite) {
                return false;
            }

            if (quality.StrictlyPositive) {
                best_finite_candidate = rCandidateData;
                best_finite_quality = quality;
                return true;
            }

            if (!has_best_finite_candidate ||
                quality.MinDeterminant >
                best_finite_quality.MinDeterminant) {
                best_finite_candidate = rCandidateData;
                best_finite_quality = quality;
                has_best_finite_candidate = true;
            }

            return false;
        };

        for (const auto& r_nodes : node_permutations) {
            auto candidate_data = create_curved_type3_volume(
                r_nodes[0],
                r_nodes[1],
                r_nodes[2],
                false);

            if (consider_candidate(candidate_data)) {
                return candidate_data;
            }

            candidate_data = create_curved_type3_volume(
                r_nodes[0],
                r_nodes[1],
                r_nodes[2],
                true);

            if (consider_candidate(candidate_data)) {
                return candidate_data;
            }
        }

        if (has_best_finite_candidate) {
            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "[CreateType3CollapsedCornerVolume] Could not find a "
                << "strictly positive quadratic type3 volume. Using the least "
                << "inverted finite curved candidate instead. No linear "
                << "fallback is used.\n"
                << "  surrogate node: " << pSurrogateNode->Id() << "\n"
                << "  projection nodes: "
                << pProjectionNode0->Id() << ", "
                << pProjectionNode1->Id() << ", "
                << pProjectionNode2->Id() << "\n"
                << "  minimum determinant on selection grid: "
                << best_finite_quality.MinDeterminant << "\n"
                << "  maximum determinant on selection grid: "
                << best_finite_quality.MaxDeterminant << "\n";

            return best_finite_candidate;
        }

        return Type3CollapsedCornerData();
    };

    auto curved_volume_data = select_valid_curved_volume();
    KRATOS_ERROR_IF_NOT(curved_volume_data.pVolume)
        << "[CreateType3CollapsedCornerVolume] Could not create a valid "
        << "quadratic type3 volume from curved type2 skin edges without "
        << "inverting the volume. No linear fallback is allowed.\n"
        << "  surrogate node: " << pSurrogateNode->Id() << "\n"
        << "  projection nodes: "
        << pProjectionNode0->Id() << ", "
        << pProjectionNode1->Id() << ", "
        << pProjectionNode2->Id() << "\n";

    KRATOS_ERROR_IF_NOT(curved_volume_data.pTopSurface)
        << "[CreateType3CollapsedCornerVolume] Valid type3 volume has null "
        << "top surface.\n";

    return curved_volume_data;
}

SnakeGapSbm3DUtilities::Type3CollapsedCornerData
SnakeGapSbm3DUtilities::CreateType3LinearCollapsedCornerVolume(
    const NodePointerType& pSurrogateNode,
    const NodePointerType& pProjectionNode0,
    const NodePointerType& pProjectionNode1,
    const NodePointerType& pProjectionNode2,
    const bool CreateTopSurface) const
{
    KRATOS_ERROR_IF_NOT(pSurrogateNode)
        << "[CreateType3LinearCollapsedCornerVolume] pSurrogateNode is null.\n";
    KRATOS_ERROR_IF_NOT(pProjectionNode0)
        << "[CreateType3LinearCollapsedCornerVolume] pProjectionNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pProjectionNode1)
        << "[CreateType3LinearCollapsedCornerVolume] pProjectionNode1 is null.\n";
    KRATOS_ERROR_IF_NOT(pProjectionNode2)
        << "[CreateType3LinearCollapsedCornerVolume] pProjectionNode2 is null.\n";

    Type3CollapsedCornerData result;
    result.ProjectionNodes = {{
        pProjectionNode0,
        pProjectionNode1,
        pProjectionNode2}};

    PointerVector<NodeType> control_points;
    control_points.push_back(pSurrogateNode);
    control_points.push_back(pProjectionNode0);
    control_points.push_back(pProjectionNode1);
    control_points.push_back(pProjectionNode1);

    control_points.push_back(pProjectionNode2);
    control_points.push_back(pProjectionNode2);
    control_points.push_back(pProjectionNode2);
    control_points.push_back(pProjectionNode2);

    result.pVolume = Kratos::make_shared<NurbsVolumeType>(
        control_points,
        std::size_t(1),
        std::size_t(1),
        std::size_t(1),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVectorDegree1());

    if (CreateTopSurface) {
        result.pTopSurface = CreateCollapsedTriangleSurface(
            pProjectionNode0,
            pProjectionNode1,
            pProjectionNode2);
    }

    return result;
}

SnakeGapSbm3DUtilities::Type3CollapsedCornerData
SnakeGapSbm3DUtilities::CreateType3CollapsedCornerVolumeFromTopControlNodes(
    const NodePointerType& pSurrogateNode,
    const std::array<NodePointerType, 3>& rProjectionNodes,
    const std::vector<NodePointerType>& rTopControlNodes,
    const bool CreateTopSurface) const
{
    KRATOS_ERROR_IF_NOT(pSurrogateNode)
        << "[CreateType3CollapsedCornerVolumeFromTopControlNodes] "
        << "pSurrogateNode is null.\n";
    KRATOS_ERROR_IF(rTopControlNodes.size() != 9)
        << "[CreateType3CollapsedCornerVolumeFromTopControlNodes] "
        << "Quadratic type3 top face requires nine control nodes. Got "
        << rTopControlNodes.size() << ".\n";

    Type3CollapsedCornerData result;
    result.ProjectionNodes = rProjectionNodes;

    PointerVector<NodeType> volume_control_points;
    PointerVector<NodeType> surface_control_points;

    for (const auto& p_top_control_node : rTopControlNodes) {
        KRATOS_ERROR_IF_NOT(p_top_control_node)
            << "[CreateType3CollapsedCornerVolumeFromTopControlNodes] "
            << "Null top control node.\n";

        if (CreateTopSurface) {
            surface_control_points.push_back(p_top_control_node);
        }
        volume_control_points.push_back(pSurrogateNode);
        volume_control_points.push_back(p_top_control_node);
    }

    if (CreateTopSurface) {
        result.pTopSurface = Kratos::make_shared<NurbsSurfaceType>(
            surface_control_points,
            std::size_t(2),
            std::size_t(2),
            CreateOpenUnitKnotVector(2),
            CreateOpenUnitKnotVector(2));
    }

    result.pVolume = Kratos::make_shared<NurbsVolumeType>(
        volume_control_points,
        std::size_t(1),
        std::size_t(2),
        std::size_t(2),
        CreateOpenUnitKnotVectorDegree1(),
        CreateOpenUnitKnotVector(2),
        CreateOpenUnitKnotVector(2));

    return result;
}

void SnakeGapSbm3DUtilities::OrientTriangleFaceAwayFromOppositeNode(
    NodePointerType& rpNode0,
    NodePointerType& rpNode1,
    NodePointerType& rpNode2,
    const NodePointerType& pOppositeNode) const
{
    KRATOS_ERROR_IF_NOT(rpNode0)
        << "[SnakeGapSbm3DUtilities::OrientTriangleFaceAwayFromOppositeNode] "
        << "First face node is null.\n";
    KRATOS_ERROR_IF_NOT(rpNode1)
        << "[SnakeGapSbm3DUtilities::OrientTriangleFaceAwayFromOppositeNode] "
        << "Second face node is null.\n";
    KRATOS_ERROR_IF_NOT(rpNode2)
        << "[SnakeGapSbm3DUtilities::OrientTriangleFaceAwayFromOppositeNode] "
        << "Third face node is null.\n";
    KRATOS_ERROR_IF_NOT(pOppositeNode)
        << "[SnakeGapSbm3DUtilities::OrientTriangleFaceAwayFromOppositeNode] "
        << "Opposite tetrahedron node is null.\n";

    const array_1d<double, 3> edge_01 =
        rpNode1->Coordinates() - rpNode0->Coordinates();
    const array_1d<double, 3> edge_02 =
        rpNode2->Coordinates() - rpNode0->Coordinates();
    const array_1d<double, 3> face_normal =
        MathUtils<double>::CrossProduct(edge_01, edge_02);
    const array_1d<double, 3> vector_to_opposite_node =
        pOppositeNode->Coordinates() - rpNode0->Coordinates();

    if (inner_prod(face_normal, vector_to_opposite_node) > 0.0) {
        std::swap(rpNode1, rpNode2);
    }
}

void SnakeGapSbm3DUtilities::RegisterType3LateralFace(
    ModelPart& rDebugModelPart,
    IndexType& rNextGeometryId,
    const IndexType SurrogateConditionId,
    const SpanKey3D& rExternalSpan,
    NodePointerType pNode0,
    NodePointerType pNode1,
    NodePointerType pNode2,
    const NodePointerType& pOppositeNode,
    const Geometry<Node>::Pointer& pNeighbourGeometry,
    const bool HasType1NeighbourPath,
    const bool HasNeighbourActiveSpan,
    const SpanKey3D& rNeighbourActiveSpan)
{
    OrientTriangleFaceAwayFromOppositeNode(
        pNode0,
        pNode1,
        pNode2,
        pOppositeNode);

    auto p_surface = GetOrCreateLateralFaceSurface(
        rDebugModelPart,
        rNextGeometryId,
        pNode0,
        pNode1,
        pNode2,
        pNeighbourGeometry);

    RegisterLateralFaceOccurrence(
        3,
        SurrogateConditionId,
        rExternalSpan,
        pNode0,
        pNode1,
        pNode2,
        pOppositeNode,
        p_surface,
        pNeighbourGeometry,
        HasType1NeighbourPath,
        HasNeighbourActiveSpan,
        rNeighbourActiveSpan);
}

void SnakeGapSbm3DUtilities::CheckType3QuadraturePointGeometries(
    const GeometriesArrayType& rQuadraturePointGeometries,
    const NodePointerType& pSurrogateNode,
    const NodePointerType& pProjectionNode0,
    const NodePointerType& pProjectionNode1,
    const NodePointerType& pProjectionNode2) const
{
    KRATOS_ERROR_IF(rQuadraturePointGeometries.size() == 0)
        << "[SnakeGapSbm3DUtilities::CheckType3QuadraturePointGeometries] "
        << "Type 3 volume has no quadrature point geometries.\n";

    constexpr double negative_determinant_tolerance = 1.0e-2;

    for (IndexType point_index = 0;
         point_index < rQuadraturePointGeometries.size();
         ++point_index) {
        const auto& p_quadrature_geometry =
            rQuadraturePointGeometries(point_index);

        KRATOS_ERROR_IF_NOT(p_quadrature_geometry)
            << "[SnakeGapSbm3DUtilities::CheckType3QuadraturePointGeometries] "
            << "Null quadrature point geometry at index "
            << point_index << ".\n";

        const double determinant_jacobian =
            p_quadrature_geometry->DeterminantOfJacobian(
                0,
                p_quadrature_geometry->GetDefaultIntegrationMethod());

        KRATOS_ERROR_IF(
            !std::isfinite(determinant_jacobian) ||
            determinant_jacobian < -negative_determinant_tolerance)
            << "[SnakeGapSbm3DUtilities::CheckType3QuadraturePointGeometries] "
            << "Invalid type 3 tetrahedron. Determinant of Jacobian must be finite "
            << "and not smaller than the accepted negative tolerance.\n"
            << "  quadrature point index: " << point_index << "\n"
            << "  determinant: " << determinant_jacobian << "\n"
            << "  negative determinant tolerance: "
            << negative_determinant_tolerance << "\n"
            << "  surrogate node: " << pSurrogateNode->Id() << "\n"
            << "  projection node 0: " << pProjectionNode0->Id() << "\n"
            << "  projection node 1: " << pProjectionNode1->Id() << "\n"
            << "  projection node 2: "
            << pProjectionNode2->Id() << "\n";
    }
}

std::array<SnakeGapSbm3DUtilities::SpanKey3D, 4>
SnakeGapSbm3DUtilities::GetSurroundingSpansOfGridEdge(
    const GridPointKey3D& rNode0,
    const GridPointKey3D& rNode1) const
{
    std::array<SpanKey3D, 4> spans;

    if (rNode0.I != rNode1.I && rNode0.J == rNode1.J && rNode0.K == rNode1.K) {
        const int i = std::min(rNode0.I, rNode1.I);
        const int j = rNode0.J;
        const int k = rNode0.K;
        spans = {{{i, j - 1, k - 1}, {i, j, k - 1}, {i, j - 1, k}, {i, j, k}}};
    } else if (rNode0.J != rNode1.J && rNode0.I == rNode1.I && rNode0.K == rNode1.K) {
        const int i = rNode0.I;
        const int j = std::min(rNode0.J, rNode1.J);
        const int k = rNode0.K;
        spans = {{{i - 1, j, k - 1}, {i, j, k - 1}, {i - 1, j, k}, {i, j, k}}};
    } else if (rNode0.K != rNode1.K && rNode0.I == rNode1.I && rNode0.J == rNode1.J) {
        const int i = rNode0.I;
        const int j = rNode0.J;
        const int k = std::min(rNode0.K, rNode1.K);
        spans = {{{i - 1, j - 1, k}, {i, j - 1, k}, {i - 1, j, k}, {i, j, k}}};
    } else {
        KRATOS_ERROR << "[SnakeGapSbm3DUtilities::GetSurroundingSpansOfGridEdge] "
                     << "The grid nodes do not define an axis-aligned grid edge.\n";
    }

    return spans;
}

bool SnakeGapSbm3DUtilities::IsValidSpan(
    const SpanKey3D& rSpan,
    const KnotSpanGridInfo& rGridInfo) const
{
    return rSpan.I >= 0 && rSpan.J >= 0 && rSpan.K >= 0 &&
           rSpan.I < static_cast<int>(rGridInfo.NumberOfSpansX) &&
           rSpan.J < static_cast<int>(rGridInfo.NumberOfSpansY) &&
           rSpan.K < static_cast<int>(rGridInfo.NumberOfSpansZ);
}

bool SnakeGapSbm3DUtilities::AreFaceAdjacentSpans(
    const SpanKey3D& rA,
    const SpanKey3D& rB) const
{
    return std::abs(rA.I - rB.I) +
           std::abs(rA.J - rB.J) +
           std::abs(rA.K - rB.K) == 1;
}

double SnakeGapSbm3DUtilities::MaximumSkinProjectionDistance(
    const KnotSpanGridInfo& rGridInfo) const
{
    return 3.5 * std::max({
        rGridInfo.SpanSizeX,
        rGridInfo.SpanSizeY,
        rGridInfo.SpanSizeZ});
}

bool SnakeGapSbm3DUtilities::ProjectPointToSkinBoundaryAlongDirection(
    const ModelPart& rSkinSubModelPart,
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rDirection,
    const double MaximumDistance,
    array_1d<double, 3>& rProjectionPoint) const
{
    const double direction_norm = norm_2(rDirection);
    KRATOS_ERROR_IF(direction_norm <= std::numeric_limits<double>::epsilon())
        << "[SnakeGapSbm3DUtilities::ProjectPointToSkinBoundaryAlongDirection] "
        << "Projection direction has zero norm.\n";

    array_1d<double, 3> direction = rDirection / direction_norm;

    const double point_on_skin_tolerance =
        1.0e-10 * std::max(1.0, MaximumDistance);
    const double ray_tolerance =
        1.0e-12 * std::max(1.0, MaximumDistance);

    bool found_projection = false;
    double best_distance = std::numeric_limits<double>::max();
    double best_closest_point_distance = std::numeric_limits<double>::max();
    array_1d<double, 3> best_closest_point = ZeroVector(3);

    auto consider_point_on_triangle = [&](
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC)
    {
        const array_1d<double, 3> closest_point =
            ClosestPointOnTriangle(rPoint, rA, rB, rC);
        const double distance = norm_2(closest_point - rPoint);
        if (distance < best_closest_point_distance) {
            best_closest_point_distance = distance;
            best_closest_point = closest_point;
        }
        if (distance <= point_on_skin_tolerance &&
            distance < best_distance) {
            best_distance = distance;
            rProjectionPoint = closest_point;
            found_projection = true;
        }
    };

    auto ray_triangle_intersection = [&](
        const array_1d<double, 3>& rRayDirection,
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC)
    {
        const array_1d<double, 3> edge_1 = rB - rA;
        const array_1d<double, 3> edge_2 = rC - rA;
        const array_1d<double, 3> p_vector =
            MathUtils<double>::CrossProduct(rRayDirection, edge_2);
        const double determinant = inner_prod(edge_1, p_vector);
        if (std::abs(determinant) <= ray_tolerance) {
            return;
        }

        const double inverse_determinant = 1.0 / determinant;
        const array_1d<double, 3> t_vector = rPoint - rA;
        const double u = inner_prod(t_vector, p_vector) * inverse_determinant;
        if (u < -ray_tolerance || u > 1.0 + ray_tolerance) {
            return;
        }

        const array_1d<double, 3> q_vector =
            MathUtils<double>::CrossProduct(t_vector, edge_1);
        const double v =
            inner_prod(rRayDirection, q_vector) * inverse_determinant;
        if (v < -ray_tolerance || u + v > 1.0 + ray_tolerance) {
            return;
        }

        const double distance =
            inner_prod(edge_2, q_vector) * inverse_determinant;
        if (distance <= point_on_skin_tolerance ||
            distance > MaximumDistance ||
            distance >= best_distance) {
            return;
        }

        best_distance = distance;
        rProjectionPoint = rPoint + distance * rRayDirection;
        found_projection = true;
    };

    auto consider_triangle = [&](
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC)
    {
        array_1d<double, 3> opposite_direction = direction;
        opposite_direction *= -1.0;

        consider_point_on_triangle(rA, rB, rC);
        ray_triangle_intersection(direction, rA, rB, rC);
        ray_triangle_intersection(opposite_direction, rA, rB, rC);
    };

    if (!mSkinProjectionTriangleData.empty()) {
        const double candidate_radius =
            MaximumDistance + ray_tolerance;
        for (const auto& r_triangle_data : mSkinProjectionTriangleData) {
            const double center_distance =
                norm_2(r_triangle_data.Center - rPoint);
            if (center_distance >
                candidate_radius + r_triangle_data.Radius) {
                continue;
            }

            consider_triangle(
                r_triangle_data.A,
                r_triangle_data.B,
                r_triangle_data.C);
        }
    } else {
        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            const auto& r_geometry = r_condition.GetGeometry();
            if (r_geometry.PointsNumber() < 3) {
                continue;
            }

            consider_triangle(
                r_geometry[0].Coordinates(),
                r_geometry[1].Coordinates(),
                r_geometry[2].Coordinates());

            if (r_geometry.PointsNumber() == 4) {
                consider_triangle(
                    r_geometry[0].Coordinates(),
                    r_geometry[2].Coordinates(),
                    r_geometry[3].Coordinates());
            }
        }
    }

    if (found_projection && best_distance <= MaximumDistance) {
        return true;
    }

    if (best_closest_point_distance <= MaximumDistance) {
        rProjectionPoint = best_closest_point;
        return true;
    }

    return false;
}

bool SnakeGapSbm3DUtilities::ProjectPointToClosestSkinBoundary(
    const ModelPart& rSkinSubModelPart,
    const array_1d<double, 3>& rPoint,
    array_1d<double, 3>& rProjectionPoint,
    double& rProjectionDistance) const
{
    bool found_projection = false;
    rProjectionDistance = std::numeric_limits<double>::max();
    rProjectionPoint = ZeroVector(3);

    auto consider_triangle = [&](
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC)
    {
        const array_1d<double, 3> closest_point =
            ClosestPointOnTriangle(rPoint, rA, rB, rC);
        const double distance = norm_2(closest_point - rPoint);
        if (distance < rProjectionDistance) {
            rProjectionDistance = distance;
            rProjectionPoint = closest_point;
            found_projection = true;
        }
    };

    if (!mSkinProjectionTriangleData.empty()) {
        for (const auto& r_triangle_data : mSkinProjectionTriangleData) {
            if (std::isfinite(rProjectionDistance) &&
                norm_2(r_triangle_data.Center - rPoint) >
                rProjectionDistance + r_triangle_data.Radius) {
                continue;
            }

            consider_triangle(
                r_triangle_data.A,
                r_triangle_data.B,
                r_triangle_data.C);
        }
    } else {
        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            const auto& r_geometry = r_condition.GetGeometry();
            if (r_geometry.PointsNumber() < 3) {
                continue;
            }

            consider_triangle(
                r_geometry[0].Coordinates(),
                r_geometry[1].Coordinates(),
                r_geometry[2].Coordinates());

            if (r_geometry.PointsNumber() == 4) {
                consider_triangle(
                    r_geometry[0].Coordinates(),
                    r_geometry[2].Coordinates(),
                    r_geometry[3].Coordinates());
            }
        }
    }

    return found_projection;
}

SnakeGapSbm3DUtilities::NodePointerType
SnakeGapSbm3DUtilities::CreateOrReuseAuxiliaryControlNode(
    ModelPart& rSkinSubModelPart,
    const array_1d<double, 3>& rPoint) const
{
    ModelPart& r_root_model_part = rSkinSubModelPart.GetRootModelPart();
    const double point_tolerance = 1.0e-12 * std::max({
        1.0,
        std::abs(rPoint[0]),
        std::abs(rPoint[1]),
        std::abs(rPoint[2])});

    for (auto& r_node : r_root_model_part.Nodes()) {
        if (norm_2(r_node.Coordinates() - rPoint) <= point_tolerance) {
            return r_root_model_part.pGetNode(r_node.Id());
        }
    }

    IndexType new_node_id = GetNextAuxiliarySkinNodeId(rSkinSubModelPart);
    while (r_root_model_part.HasNode(new_node_id)) {
        ++new_node_id;
    }

    return r_root_model_part.CreateNewNode(
        new_node_id,
        rPoint[0],
        rPoint[1],
        rPoint[2]);
}

std::vector<SnakeGapSbm3DUtilities::NodePointerType>
SnakeGapSbm3DUtilities::CreateStraightSkinEdgeControlNodes(
    ModelPart& rSkinSubModelPart,
    const NodePointerType& pSkinNode0,
    const NodePointerType& pSkinNode1) const
{
    static_cast<void>(rSkinSubModelPart);

    KRATOS_ERROR_IF_NOT(pSkinNode0)
        << "[SnakeGapSbm3DUtilities::CreateStraightSkinEdgeControlNodes] "
        << "First skin node is null.\n";
    KRATOS_ERROR_IF_NOT(pSkinNode1)
        << "[SnakeGapSbm3DUtilities::CreateStraightSkinEdgeControlNodes] "
        << "Second skin node is null.\n";

    std::vector<NodePointerType> control_nodes;
    control_nodes.reserve(mGapApproximationOrder + 1);
    control_nodes.push_back(pSkinNode0);

    for (std::size_t i = 1; i < mGapApproximationOrder; ++i) {
        const double t =
            static_cast<double>(i) /
            static_cast<double>(mGapApproximationOrder);
        array_1d<double, 3> point =
            (1.0 - t) * pSkinNode0->Coordinates() +
            t * pSkinNode1->Coordinates();
        control_nodes.push_back(NodePointerType(new Node(0, point)));
    }

    control_nodes.push_back(pSkinNode1);

    return control_nodes;
}

std::vector<SnakeGapSbm3DUtilities::NodePointerType>
SnakeGapSbm3DUtilities::GetOrCreateFinalSkinEdgeControlNodes(
    ModelPart& rSkinSubModelPart,
    const NodePointerType& pSkinNode0,
    const NodePointerType& pSkinNode1,
    const double MaximumProjectionDistance,
    const std::string& rContext)
{
    KRATOS_ERROR_IF_NOT(pSkinNode0)
        << "[SnakeGapSbm3DUtilities::GetOrCreateFinalSkinEdgeControlNodes] "
        << "First skin node is null.\n";
    KRATOS_ERROR_IF_NOT(pSkinNode1)
        << "[SnakeGapSbm3DUtilities::GetOrCreateFinalSkinEdgeControlNodes] "
        << "Second skin node is null.\n";

    const bool reverse_orientation = pSkinNode1->Id() < pSkinNode0->Id();
    const SkinEdgeKey key(
        std::min(pSkinNode0->Id(), pSkinNode1->Id()),
        std::max(pSkinNode0->Id(), pSkinNode1->Id()));

    const auto registry_it = mCurvedEdgeRegistry.find(key);
    if (registry_it != mCurvedEdgeRegistry.end()) {
        auto control_nodes = registry_it->second.CurrentControlNodes;
        if (reverse_orientation) {
            std::reverse(control_nodes.begin(), control_nodes.end());
        }
        return control_nodes;
    }

    const NodePointerType p_canonical_node_0 =
        reverse_orientation ? pSkinNode1 : pSkinNode0;
    const NodePointerType p_canonical_node_1 =
        reverse_orientation ? pSkinNode0 : pSkinNode1;

    CurvedEdgeData edge_data;
    edge_data.Key = key;
    edge_data.LinearControlNodes =
        CreateStraightSkinEdgeControlNodes(
            rSkinSubModelPart,
            p_canonical_node_0,
            p_canonical_node_1);
    edge_data.CurvedControlNodes = edge_data.LinearControlNodes;
    edge_data.CurrentControlNodes = edge_data.LinearControlNodes;

    if (mGapApproximationOrder > 1) {
        bool created_curved_edge = true;
        std::vector<NodePointerType> curved_control_nodes;
        curved_control_nodes.reserve(mGapApproximationOrder + 1);
        curved_control_nodes.push_back(p_canonical_node_0);

        const array_1d<double, 3> edge_vector =
            p_canonical_node_1->Coordinates() -
            p_canonical_node_0->Coordinates();

        array_1d<double, 3> projection_direction = ZeroVector(3);
        bool has_projection_direction = false;
        array_1d<double, 3> unused_closest_point = ZeroVector(3);
        double unused_projection_distance = 0.0;
        if (ProjectPointToClosestSkinBoundary(
                rSkinSubModelPart,
                0.5 * (p_canonical_node_0->Coordinates() +
                       p_canonical_node_1->Coordinates()),
                unused_closest_point,
                unused_projection_distance)) {
            projection_direction =
                unused_closest_point -
                0.5 * (p_canonical_node_0->Coordinates() +
                       p_canonical_node_1->Coordinates());
            has_projection_direction =
                norm_2(projection_direction) >
                std::numeric_limits<double>::epsilon();
        }

        if (!has_projection_direction) {
            array_1d<double, 3> fallback_axis = ZeroVector(3);
            fallback_axis[0] = 1.0;
            if (std::abs(edge_vector[0]) >
                std::abs(edge_vector[1])) {
                fallback_axis[0] = 0.0;
                fallback_axis[1] = 1.0;
            }
            projection_direction =
                MathUtils<double>::CrossProduct(edge_vector, fallback_axis);
            has_projection_direction =
                norm_2(projection_direction) >
                std::numeric_limits<double>::epsilon();
        }

        if (has_projection_direction) {
            for (std::size_t i = 1; i < mGapApproximationOrder; ++i) {
                const double t =
                    static_cast<double>(i) /
                    static_cast<double>(mGapApproximationOrder);
                const array_1d<double, 3> linear_point =
                    (1.0 - t) * p_canonical_node_0->Coordinates() +
                    t * p_canonical_node_1->Coordinates();

                array_1d<double, 3> projected_point = ZeroVector(3);
                if (!ProjectPointToSkinBoundaryAlongDirection(
                        rSkinSubModelPart,
                        linear_point,
                        projection_direction,
                        MaximumProjectionDistance,
                        projected_point)) {
                    created_curved_edge = false;
                    break;
                }

                if (mGapApproximationOrder == 2) {
                    const array_1d<double, 3> bezier_control_point =
                        2.0 * projected_point -
                        0.5 * (p_canonical_node_0->Coordinates() +
                               p_canonical_node_1->Coordinates());
                    curved_control_nodes.push_back(
                        NodePointerType(new Node(0, bezier_control_point)));
                } else {
                    curved_control_nodes.push_back(
                        NodePointerType(new Node(0, projected_point)));
                }
            }
        } else {
            created_curved_edge = false;
        }

        if (created_curved_edge) {
            curved_control_nodes.push_back(p_canonical_node_1);
            edge_data.CurvedControlNodes = curved_control_nodes;
            edge_data.CurrentControlNodes = curved_control_nodes;
            edge_data.Alpha = 1.0;
            edge_data.IsCurved = true;
        } else {
            mLinearSkinEdges.insert(key);
            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "[GetOrCreateFinalSkinEdgeControlNodes] Falling back to "
                << "a linear skin edge after high-order projection failed.\n"
                << "  context: " << rContext << "\n"
                << "  edge ids: " << key.first << " - " << key.second << "\n";
        }
    }

    mCurvedEdgeRegistry.emplace(key, edge_data);
    mSkinEdgeControlNodes[key] = edge_data.CurrentControlNodes;

    auto control_nodes = edge_data.CurrentControlNodes;
    if (reverse_orientation) {
        std::reverse(control_nodes.begin(), control_nodes.end());
    }
    return control_nodes;
}

std::vector<SnakeGapSbm3DUtilities::NodePointerType>
SnakeGapSbm3DUtilities::GetOrCreateQuadraticSkinEdgeControlNodes(
    ModelPart& rSkinSubModelPart,
    const NodePointerType& pSkinNode0,
    const NodePointerType& pSkinNode1,
    const double MaximumProjectionDistance,
    const std::string& rContext)
{
    KRATOS_ERROR_IF_NOT(pSkinNode0)
        << "[SnakeGapSbm3DUtilities::GetOrCreateQuadraticSkinEdgeControlNodes] "
        << "First skin node is null.\n";
    KRATOS_ERROR_IF_NOT(pSkinNode1)
        << "[SnakeGapSbm3DUtilities::GetOrCreateQuadraticSkinEdgeControlNodes] "
        << "Second skin node is null.\n";

    if (mGapApproximationOrder == 1) {
        return {pSkinNode0, pSkinNode1};
    }

    KRATOS_ERROR_IF(mGapApproximationOrder != 2)
        << "[SnakeGapSbm3DUtilities::GetOrCreateQuadraticSkinEdgeControlNodes] "
        << "Only quadratic skin edges are implemented. Requested order: "
        << mGapApproximationOrder << ".\n";

    const bool reverse_orientation = pSkinNode1->Id() < pSkinNode0->Id();
    const auto key = std::make_pair(
        std::min(pSkinNode0->Id(), pSkinNode1->Id()),
        std::max(pSkinNode0->Id(), pSkinNode1->Id()));

    const auto cached_it = mSkinEdgeControlNodes.find(key);
    if (cached_it != mSkinEdgeControlNodes.end()) {
        std::vector<NodePointerType> control_nodes = cached_it->second;
        if (reverse_orientation) {
            std::reverse(control_nodes.begin(), control_nodes.end());
        }
        return control_nodes;
    }

    const array_1d<double, 3> midpoint =
        0.5 * (pSkinNode0->Coordinates() + pSkinNode1->Coordinates());

    array_1d<double, 3> projected_midpoint = ZeroVector(3);
    double projection_distance = std::numeric_limits<double>::max();
    const bool closest_projection_succeeded =
        ProjectPointToClosestSkinBoundary(
            rSkinSubModelPart,
            midpoint,
            projected_midpoint,
            projection_distance);

    // KRATOS_INFO("SnakeGapSbm3DUtilities")
    //     << "[GetOrCreateQuadraticSkinEdgeControlNodes] "
    //     << "Closest midpoint skin projection for quadratic skin edge.\n"
    //     << "  context: " << rContext << "\n"
    //     << "  edge ids: " << pSkinNode0->Id()
    //     << " - " << pSkinNode1->Id() << "\n"
    //     << "  P0: " << pSkinNode0->Coordinates() << "\n"
    //     << "  P1: " << pSkinNode1->Coordinates() << "\n"
    //     << "  M: " << midpoint << "\n"
    //     << "  Q: " << projected_midpoint << "\n"
    //     << "  projection distance: " << projection_distance << "\n"
    //     << "  closest projection succeeded: "
    //     << closest_projection_succeeded << "\n";

    KRATOS_ERROR_IF_NOT(closest_projection_succeeded)
        << "[SnakeGapSbm3DUtilities::GetOrCreateQuadraticSkinEdgeControlNodes] "
        << "Could not find closest skin projection for quadratic skin edge "
        << "interpolation point.\n"
        << "  context: " << rContext << "\n"
        << "  edge ids: " << pSkinNode0->Id()
        << " - " << pSkinNode1->Id() << "\n"
        << "  P0: " << pSkinNode0->Coordinates() << "\n"
        << "  P1: " << pSkinNode1->Coordinates() << "\n"
        << "  M: " << midpoint << "\n"
        << "  Q: " << projected_midpoint << "\n"
        << "  projection distance: " << projection_distance << "\n"
        << "  closest projection succeeded: "
        << closest_projection_succeeded << "\n";

    KRATOS_ERROR_IF(projection_distance > MaximumProjectionDistance)
        << "[SnakeGapSbm3DUtilities::GetOrCreateQuadraticSkinEdgeControlNodes] "
        << "Closest skin projection for quadratic skin edge interpolation "
        << "point is farther than the allowed distance.\n"
        << "  context: " << rContext << "\n"
        << "  edge ids: " << pSkinNode0->Id()
        << " - " << pSkinNode1->Id() << "\n"
        << "  P0: " << pSkinNode0->Coordinates() << "\n"
        << "  P1: " << pSkinNode1->Coordinates() << "\n"
        << "  M: " << midpoint << "\n"
        << "  Q: " << projected_midpoint << "\n"
        << "  projection distance: " << projection_distance << "\n"
        << "  maximum allowed distance: " << MaximumProjectionDistance << "\n"
        << "  closest projection succeeded: "
        << closest_projection_succeeded << "\n";

    array_1d<double, 3> bezier_control_point =
        2.0 * projected_midpoint -
        0.5 * (pSkinNode0->Coordinates() + pSkinNode1->Coordinates());

    // const array_1d<double, 3> bezier_control_point =
    //     0.5 * (pSkinNode0->Coordinates() + pSkinNode1->Coordinates());

    auto p_control_node =
        CreateOrReuseAuxiliaryControlNode(rSkinSubModelPart, bezier_control_point);

    std::vector<NodePointerType> canonical_control_nodes;
    if (reverse_orientation) {
        canonical_control_nodes = {pSkinNode1, p_control_node, pSkinNode0};
    } else {
        canonical_control_nodes = {pSkinNode0, p_control_node, pSkinNode1};
    }

    mSkinEdgeControlNodes.emplace(key, canonical_control_nodes);

    std::vector<NodePointerType> control_nodes = canonical_control_nodes;
    if (reverse_orientation) {
        std::reverse(control_nodes.begin(), control_nodes.end());
    }
    return control_nodes;
}

const std::vector<SnakeGapSbm3DUtilities::NodePointerType>*
SnakeGapSbm3DUtilities::FindCachedSkinEdgeControlNodes(
    const NodePointerType& pNode0,
    const NodePointerType& pNode1) const
{
    if (!pNode0 || !pNode1 || mGapApproximationOrder == 1) {
        return nullptr;
    }

    const auto key = std::make_pair(
        std::min(pNode0->Id(), pNode1->Id()),
        std::max(pNode0->Id(), pNode1->Id()));

    const auto cached_it = mSkinEdgeControlNodes.find(key);
    if (cached_it == mSkinEdgeControlNodes.end()) {
        return nullptr;
    }

    return &cached_it->second;
}

void SnakeGapSbm3DUtilities::RegisterLateralFaceOccurrence(
    const int GapType,
    const IndexType SurrogateConditionId,
    const SpanKey3D& rExternalSpan,
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pNode2,
    const NodePointerType& pOppositeNode,
    const NurbsSurfaceType::Pointer& pSurfaceGeometry,
    const Geometry<Node>::Pointer& pNeighbourGeometry,
    const bool HasType1NeighbourPath,
    const bool HasNeighbourActiveSpan,
    const SpanKey3D& rNeighbourActiveSpan)
{
    KRATOS_ERROR_IF_NOT(pNode0)
        << "[RegisterLateralFaceOccurrence] pNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pNode1)
        << "[RegisterLateralFaceOccurrence] pNode1 is null.\n";
    KRATOS_ERROR_IF_NOT(pNode2)
        << "[RegisterLateralFaceOccurrence] pNode2 is null.\n";
    KRATOS_ERROR_IF_NOT(pOppositeNode)
        << "[RegisterLateralFaceOccurrence] opposite node is null.\n";
    KRATOS_ERROR_IF_NOT(pSurfaceGeometry)
        << "[RegisterLateralFaceOccurrence] surface geometry is null.\n";
    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[RegisterLateralFaceOccurrence] neighbour geometry is null.\n";

    const auto key = MakeCanonicalFaceKey3D(
        pNode0,
        pNode1,
        pNode2);

    LateralFaceOccurrence occurrence;
    occurrence.GapType = GapType;
    occurrence.SurrogateConditionId = SurrogateConditionId;
    occurrence.ExternalSpan = rExternalSpan;
    occurrence.pGeometry = pSurfaceGeometry;
    occurrence.pNeighbourGeometry = pNeighbourGeometry;
    occurrence.pOppositeNode = pOppositeNode;
    occurrence.FaceNodes = {{pNode0, pNode1, pNode2}};
    occurrence.HasType1NeighbourPath = HasType1NeighbourPath;
    occurrence.HasNeighbourActiveSpan = HasNeighbourActiveSpan;
    occurrence.NeighbourActiveSpan = rNeighbourActiveSpan;

    auto& r_occurrences = mLateralFaceRegistry[key];
    KRATOS_ERROR_IF(r_occurrences.size() >= 2)
        << "[SnakeGapSbm3DUtilities::RegisterLateralFaceOccurrence] "
        << "Non-manifold lateral face for node ids "
        << key.NodeIds[0] << ", "
        << key.NodeIds[1] << ", "
        << key.NodeIds[2] << ".\n";

    r_occurrences.push_back(occurrence);

    if (r_occurrences.size() == 2) {
        for (const auto& r_closed_occurrence : r_occurrences) {
            KRATOS_ERROR_IF_NOT(r_closed_occurrence.pGeometry)
                << "[SnakeGapSbm3DUtilities::RegisterLateralFaceOccurrence] "
                << "Closed lateral face occurrence has null surface geometry.\n";
        }
    }
}

SnakeGapSbm3DUtilities::NurbsSurfaceType::Pointer
SnakeGapSbm3DUtilities::GetOrCreateLateralFaceSurface(
    ModelPart& rDebugModelPart,
    IndexType& rNextGeometryId,
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pNode2,
    const Geometry<Node>::Pointer& pNeighbourGeometry)
{
    KRATOS_ERROR_IF_NOT(pNode0)
        << "[GetOrCreateLateralFaceSurface] pNode0 is null.\n";
    KRATOS_ERROR_IF_NOT(pNode1)
        << "[GetOrCreateLateralFaceSurface] pNode1 is null.\n";
    KRATOS_ERROR_IF_NOT(pNode2)
        << "[GetOrCreateLateralFaceSurface] pNode2 is null.\n";
    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[GetOrCreateLateralFaceSurface] neighbour geometry is null.\n";

    const auto key = MakeCanonicalFaceKey3D(
        pNode0,
        pNode1,
        pNode2);

    const auto registry_it = mLateralFaceRegistry.find(key);
    if (registry_it != mLateralFaceRegistry.end()) {
        for (const auto& r_occurrence : registry_it->second) {
            if (r_occurrence.pGeometry) {
                AddUniqueNeighbourGeometry(
                    *r_occurrence.pGeometry,
                    pNeighbourGeometry);
                return r_occurrence.pGeometry;
            }
        }
    }

    auto p_surface = CreateCollapsedTriangleSurface(
        pNode0,
        pNode1,
        pNode2);

    p_surface->SetId(rNextGeometryId++);
    AddNeighbourGeometry(*p_surface, pNeighbourGeometry);
    if (mStoreGapDebugGeometries) {
        rDebugModelPart.AddGeometry(p_surface);
    }

    return p_surface;
}

SnakeGapSbm3DUtilities::Type2CreationResult
SnakeGapSbm3DUtilities::CreateType2GapGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    Type1CreationResult& rType1CreationResult,
    const ExternalSpanDataMap& rExternalSpans,
    const KnotSpanGridInfo& rGridInfo,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives)
{
    static_cast<void>(NumberOfShapeFunctionsDerivatives);

    Type2CreationResult result;

    struct Type2VolumeKey {
        SurrogateEdgeKey3D EdgeKey;
        SpanKey3D SpanA;
        SpanKey3D SpanB;

        bool operator<(const Type2VolumeKey& rOther) const
        {
            if (!(EdgeKey == rOther.EdgeKey)) {
                return EdgeKey.NodeIds < rOther.EdgeKey.NodeIds;
            }

            return std::tie(SpanA, SpanB) < std::tie(rOther.SpanA, rOther.SpanB);
        }
    };

    struct ActiveEdgeData {
        SurrogateEdgeKey3D EdgeKey;
        NodePointerType pEdgeNode0;
        NodePointerType pEdgeNode1;
        GridPointKey3D EdgeGridNode0;
        GridPointKey3D EdgeGridNode1;
        Geometry<Node>::Pointer pNeighbourGeometry;
        IndexType SurrogateConditionId = 0;
        SpanKey3D ReferenceExternalSpan;
        std::map<SpanKey3D, std::size_t> Type1FaceIndexBySpan;
    };

    struct ExternalSpanCandidate {
        SpanKey3D Span;
        const ExternalSpanData* pData = nullptr;
        NodePointerType pProjectionNode;
        std::size_t Type1FaceIndex = std::numeric_limits<std::size_t>::max();

        bool HasType1Face() const
        {
            return Type1FaceIndex != std::numeric_limits<std::size_t>::max();
        }
    };

    ModelPart& r_gap_type2_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType2Debug");

    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    result.Summary.NumberOfCandidateEdges =
        rType1CreationResult.Type1LateralFacesByEdge.size();

    std::set<Type2VolumeKey> created_volume_keys;
    std::size_t number_type1_type1_pairs = 0;
    std::size_t number_type1_type2_pairs = 0;
    std::size_t number_type2_type2_pairs = 0;
    std::size_t number_pairs_with_unbacked_type1_candidate = 0;
    std::size_t number_skipped_unbacked_type1_candidate_pairs = 0;
    double total_type2_tetra_volume = 0.0;
    double total_type2_quadrature_weight = 0.0;
    double max_type2_relative_weight_error = 0.0;
    std::size_t number_pairs_with_common_active_span = 0;
    std::size_t number_pairs_without_common_active_span = 0;
    std::size_t number_type1_type1_pairs_without_common_active_span = 0;
    std::size_t number_type1_type2_pairs_without_common_active_span = 0;
    std::size_t number_type2_neighbour_selected_from_reused_face = 0;
    std::size_t number_type2_reused_lateral_face_candidates = 0;
    std::size_t number_type2_reused_lateral_face_candidates_with_common_active_span = 0;
    std::size_t number_type2_neighbour_selected_from_s0_projection_face = 0;
    std::size_t number_type2_neighbour_selected_from_s1_projection_face = 0;
    std::size_t number_type2_neighbour_selected_from_candidate_a_edge_face = 0;
    std::size_t number_type2_neighbour_selected_from_candidate_b_edge_face = 0;

    std::vector<ActiveEdgeData> active_edges;
    active_edges.reserve(rType1CreationResult.Type1LateralFacesByEdge.size());
    for (const auto& r_edge_entry : rType1CreationResult.Type1LateralFacesByEdge) {
        KRATOS_ERROR_IF(r_edge_entry.second.empty())
            << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
            << "Active edge has no type1 lateral faces.\n";

        const auto reference_face_index = r_edge_entry.second.front();
        KRATOS_ERROR_IF(reference_face_index >= rType1CreationResult.Type1LateralFaces.size())
            << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
            << "Invalid type1 lateral face index for active edge.\n";

        const auto& r_reference_face =
            rType1CreationResult.Type1LateralFaces[reference_face_index];

        ActiveEdgeData edge_data;
        edge_data.EdgeKey = r_edge_entry.first;
        edge_data.pEdgeNode0 = r_reference_face.pEdgeNode0;
        edge_data.pEdgeNode1 = r_reference_face.pEdgeNode1;
        edge_data.EdgeGridNode0 = r_reference_face.EdgeGridNode0;
        edge_data.EdgeGridNode1 = r_reference_face.EdgeGridNode1;
        edge_data.pNeighbourGeometry = r_reference_face.pNeighbourGeometry;
        edge_data.SurrogateConditionId = r_reference_face.SurrogateConditionId;
        edge_data.ReferenceExternalSpan = r_reference_face.ExternalSpan;

        KRATOS_ERROR_IF_NOT(edge_data.pEdgeNode0)
            << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
            << "Active edge has null first node.\n";
        KRATOS_ERROR_IF_NOT(edge_data.pEdgeNode1)
            << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
            << "Active edge has null second node.\n";
        KRATOS_ERROR_IF_NOT(edge_data.pNeighbourGeometry)
            << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
            << "Active edge has null neighbour geometry.\n";

        for (const auto type1_face_index : r_edge_entry.second) {
            KRATOS_ERROR_IF(type1_face_index >= rType1CreationResult.Type1LateralFaces.size())
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Invalid type1 lateral face index for active edge.\n";

            const auto& r_type1_face =
                rType1CreationResult.Type1LateralFaces[type1_face_index];
            KRATOS_ERROR_IF_NOT(r_type1_face.pProjectionNode)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Type 1 lateral face has null projection node.\n";

            edge_data.Type1FaceIndexBySpan.emplace(
                r_type1_face.ExternalSpan,
                type1_face_index);
        }

        active_edges.push_back(std::move(edge_data));
    }

    for (const auto& r_edge_data : active_edges) {
        auto p_neighbour_geometry = r_edge_data.pNeighbourGeometry;

        const auto surrounding_spans = GetSurroundingSpansOfGridEdge(
            r_edge_data.EdgeGridNode0,
            r_edge_data.EdgeGridNode1);

        std::vector<ExternalSpanCandidate> candidates;
        candidates.reserve(surrounding_spans.size());
        for (const auto& r_span : surrounding_spans) {
            if (!IsValidSpan(r_span, rGridInfo)) {
                continue;
            }

            const auto external_span_it = rExternalSpans.find(r_span);
            if (external_span_it == rExternalSpans.end()) {
                continue;
            }

            const auto& r_span_data = external_span_it->second;
            if (r_span_data.Type == GapSpanType::Type3 ||
                r_span_data.Type == GapSpanType::Undefined) {
                continue;
            }

            if (!r_span_data.HasProjectionNode()) {
                ++result.Summary.NumberOfSkippedEdges;
                continue;
            }

            ExternalSpanCandidate candidate;
            candidate.Span = r_span;
            candidate.pData = &r_span_data;
            candidate.pProjectionNode = r_span_data.pProjectionNode;

            const auto type1_face_it =
                r_edge_data.Type1FaceIndexBySpan.find(r_span);
            if (type1_face_it != r_edge_data.Type1FaceIndexBySpan.end()) {
                candidate.Type1FaceIndex = type1_face_it->second;
            }

            candidates.push_back(candidate);
        }

        for (std::size_t i = 0; i < candidates.size(); ++i) {
            for (std::size_t j = i + 1; j < candidates.size(); ++j) {
            const auto& r_candidate_a = candidates[i];
            const auto& r_candidate_b = candidates[j];

            if (!AreFaceAdjacentSpans(r_candidate_a.Span, r_candidate_b.Span)) {
                continue;
            }

            if (r_candidate_a.pData->ProjectionNodeId == r_candidate_b.pData->ProjectionNodeId) {
                ++result.Summary.NumberOfSkippedEdges;
                continue;
            }

            const bool candidate_a_is_type1 =
                r_candidate_a.pData->Type == GapSpanType::Type1;
            const bool candidate_b_is_type1 =
                r_candidate_b.pData->Type == GapSpanType::Type1;

            const bool has_unbacked_type1_candidate =
                (candidate_a_is_type1 && !r_candidate_a.HasType1Face()) ||
                (candidate_b_is_type1 && !r_candidate_b.HasType1Face());

            auto p_projection_node_a = r_candidate_a.pProjectionNode;
            auto p_projection_node_b = r_candidate_b.pProjectionNode;
            KRATOS_ERROR_IF_NOT(p_projection_node_a)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Projection node #" << r_candidate_a.pData->ProjectionNodeId
                << " from external span " << SpanToString(r_candidate_a.Span)
                << " has no stored node pointer.\n";
            KRATOS_ERROR_IF_NOT(p_projection_node_b)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Projection node #" << r_candidate_b.pData->ProjectionNodeId
                << " from external span " << SpanToString(r_candidate_b.Span)
                << " has no stored node pointer.\n";

            SpanKey3D span_a = r_candidate_a.Span;
            SpanKey3D span_b = r_candidate_b.Span;
            if (span_b < span_a) {
                std::swap(span_a, span_b);
            }

            const Type2VolumeKey volume_key{
                r_edge_data.EdgeKey,
                span_a,
                span_b};

            if (created_volume_keys.find(volume_key) != created_volume_keys.end()) {
                continue;
            }

            auto get_candidate_surrogate_condition_id = [&](
                const ExternalSpanCandidate& rCandidate)
            {
                if (rCandidate.HasType1Face()) {
                    return rType1CreationResult.Type1LateralFaces[
                        rCandidate.Type1FaceIndex].SurrogateConditionId;
                }

                return r_edge_data.SurrogateConditionId;
            };

            auto span_list_contains = [](
                const std::vector<SpanKey3D>& rSpans,
                const SpanKey3D& rSpan)
            {
                return std::find(rSpans.begin(), rSpans.end(), rSpan) !=
                       rSpans.end();
            };

            bool candidates_have_common_active_span = false;
            for (const auto& r_active_span :
                 r_candidate_a.pData->AdjacentActiveSpans) {
                if (span_list_contains(
                        r_candidate_b.pData->AdjacentActiveSpans,
                        r_active_span)) {
                    candidates_have_common_active_span = true;
                    break;
                }
            }

            if (candidates_have_common_active_span) {
                ++number_pairs_with_common_active_span;
            } else {
                ++number_pairs_without_common_active_span;
                if (candidate_a_is_type1 && candidate_b_is_type1) {
                    ++number_type1_type1_pairs_without_common_active_span;
                } else if (candidate_a_is_type1 || candidate_b_is_type1) {
                    ++number_type1_type2_pairs_without_common_active_span;
                }
            }

            auto p_type2_volume = CreateType2LinearCollapsedEdgeVolume(
                r_edge_data.pEdgeNode0,
                r_edge_data.pEdgeNode1,
                p_projection_node_a,
                p_projection_node_b);

            const array_1d<double, 3> volume_center =
                p_type2_volume->Center().Coordinates();

            const auto& p_s0 = r_edge_data.pEdgeNode0;
            const auto& p_s1 = r_edge_data.pEdgeNode1;
            const auto& p_pa = p_projection_node_a;
            const auto& p_pb = p_projection_node_b;

            p_neighbour_geometry = nullptr;
            bool selected_neighbour_has_type1_path = false;
            bool selected_neighbour_has_active_span = false;
            SpanKey3D selected_neighbour_active_span;
            std::string selected_neighbour_source;
            CanonicalFaceKey3D selected_neighbour_face_key;
            double selected_neighbour_distance =
                std::numeric_limits<double>::max();
            std::size_t number_reused_faces_for_volume = 0;
            std::size_t number_reused_faces_with_common_active_span_for_volume = 0;

            auto select_neighbour_from_reused_face = [&](
                const char* pFaceName,
                const NodePointerType& pNode0,
                const NodePointerType& pNode1,
                const NodePointerType& pNode2)
            {
                const auto face_key = MakeCanonicalFaceKey3D(
                    pNode0,
                    pNode1,
                    pNode2);

                const auto face_registry_it =
                    mLateralFaceRegistry.find(face_key);
                if (face_registry_it == mLateralFaceRegistry.end()) {
                    return;
                }

                const auto& r_existing_occurrences =
                    face_registry_it->second;
                KRATOS_ERROR_IF(r_existing_occurrences.empty())
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Reused lateral face has no occurrences.\n";
                KRATOS_ERROR_IF(r_existing_occurrences.size() >= 2)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Type2 volume would reuse an already closed lateral face.\n"
                    << "  face: " << pFaceName << "\n"
                    << "  face node ids: "
                    << face_key.NodeIds[0] << ", "
                    << face_key.NodeIds[1] << ", "
                    << face_key.NodeIds[2] << "\n"
                    << "  occurrences: " << r_existing_occurrences.size() << "\n"
                    << "  edge nodes: "
                    << r_edge_data.EdgeKey.NodeIds[0] << " - "
                    << r_edge_data.EdgeKey.NodeIds[1] << "\n"
                    << "  spans: " << SpanToString(r_candidate_a.Span)
                    << " - " << SpanToString(r_candidate_b.Span) << "\n";

                const auto& r_source_occurrence =
                    r_existing_occurrences.front();
                KRATOS_ERROR_IF_NOT(r_source_occurrence.pNeighbourGeometry)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Reused lateral face has null neighbour geometry.\n"
                    << "  face: " << pFaceName << "\n";

                ++number_reused_faces_for_volume;
                ++number_type2_reused_lateral_face_candidates;

                const bool source_has_common_active_span =
                    r_source_occurrence.HasNeighbourActiveSpan &&
                    span_list_contains(
                        r_candidate_a.pData->AdjacentActiveSpans,
                        r_source_occurrence.NeighbourActiveSpan) &&
                    span_list_contains(
                        r_candidate_b.pData->AdjacentActiveSpans,
                        r_source_occurrence.NeighbourActiveSpan);

                if (!source_has_common_active_span) {
                    return;
                }

                ++number_reused_faces_with_common_active_span_for_volume;
                ++number_type2_reused_lateral_face_candidates_with_common_active_span;

                const double distance =
                    norm_2(
                        r_source_occurrence.pNeighbourGeometry->Center().Coordinates() -
                        volume_center);

                if (!p_neighbour_geometry ||
                    distance < selected_neighbour_distance) {
                    p_neighbour_geometry =
                        r_source_occurrence.pNeighbourGeometry;
                    selected_neighbour_distance = distance;
                    selected_neighbour_has_type1_path =
                        r_source_occurrence.HasType1NeighbourPath;
                    selected_neighbour_has_active_span =
                        r_source_occurrence.HasNeighbourActiveSpan;
                    selected_neighbour_active_span =
                        r_source_occurrence.NeighbourActiveSpan;
                    selected_neighbour_source = pFaceName;
                    selected_neighbour_face_key = face_key;
                }
            };

            select_neighbour_from_reused_face(
                "S0-Pa-Pb",
                p_s0,
                p_pa,
                p_pb);
            select_neighbour_from_reused_face(
                "S1-Pb-Pa",
                p_s1,
                p_pb,
                p_pa);
            select_neighbour_from_reused_face(
                "S0-S1-Pa",
                p_s0,
                p_s1,
                p_pa);
            select_neighbour_from_reused_face(
                "S1-S0-Pb",
                p_s1,
                p_s0,
                p_pb);

            KRATOS_ERROR_IF(number_reused_faces_for_volume == 0)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Type2 volume has no reused lateral face. "
                << "All lateral faces would be new.\n"
                << "  edge nodes: "
                << r_edge_data.EdgeKey.NodeIds[0] << " - "
                << r_edge_data.EdgeKey.NodeIds[1] << "\n"
                << "  span 0: " << SpanToString(r_candidate_a.Span)
                << " type: " << GapSpanTypeToString(r_candidate_a.pData->Type)
                << " projection node: " << p_projection_node_a->Id() << "\n"
                << "  span 1: " << SpanToString(r_candidate_b.Span)
                << " type: " << GapSpanTypeToString(r_candidate_b.pData->Type)
                << " projection node: " << p_projection_node_b->Id() << "\n";

            KRATOS_ERROR_IF(number_reused_faces_with_common_active_span_for_volume == 0)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Type2 volume has reused lateral faces, but none carries "
                << "a neighbour active span common to both external spans.\n"
                << "  reused faces: " << number_reused_faces_for_volume << "\n"
                << "  edge nodes: "
                << r_edge_data.EdgeKey.NodeIds[0] << " - "
                << r_edge_data.EdgeKey.NodeIds[1] << "\n"
                << "  span 0: " << SpanToString(r_candidate_a.Span)
                << " type: " << GapSpanTypeToString(r_candidate_a.pData->Type)
                << " projection node: " << p_projection_node_a->Id() << "\n"
                << "  span 1: " << SpanToString(r_candidate_b.Span)
                << " type: " << GapSpanTypeToString(r_candidate_b.pData->Type)
                << " projection node: " << p_projection_node_b->Id() << "\n";

            KRATOS_ERROR_IF_NOT(selected_neighbour_has_type1_path)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Selected neighbour geometry for type2 volume does not "
                << "have a continuous path to a type1 volume.\n"
                << "  selected source face: " << selected_neighbour_source << "\n"
                << "  selected source face node ids: "
                << selected_neighbour_face_key.NodeIds[0] << ", "
                << selected_neighbour_face_key.NodeIds[1] << ", "
                << selected_neighbour_face_key.NodeIds[2] << "\n"
                << "  edge nodes: "
                << r_edge_data.EdgeKey.NodeIds[0] << " - "
                << r_edge_data.EdgeKey.NodeIds[1] << "\n"
                << "  span 0: " << SpanToString(r_candidate_a.Span) << "\n"
                << "  span 1: " << SpanToString(r_candidate_b.Span) << "\n";

            KRATOS_ERROR_IF_NOT(selected_neighbour_has_active_span)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Selected neighbour geometry for type2 volume has no "
                << "active span provenance.\n"
                << "  selected source face: " << selected_neighbour_source << "\n";

            ++number_type2_neighbour_selected_from_reused_face;
            if (selected_neighbour_source == "S0-Pa-Pb") {
                ++number_type2_neighbour_selected_from_s0_projection_face;
            } else if (selected_neighbour_source == "S1-Pb-Pa") {
                ++number_type2_neighbour_selected_from_s1_projection_face;
            } else if (selected_neighbour_source == "S0-S1-Pa") {
                ++number_type2_neighbour_selected_from_candidate_a_edge_face;
            } else if (selected_neighbour_source == "S1-S0-Pb") {
                ++number_type2_neighbour_selected_from_candidate_b_edge_face;
            }

            std::vector<Geometry<Node>::Pointer> type2_neighbour_geometries;
            AddUniqueNeighbourGeometry(
                type2_neighbour_geometries,
                p_neighbour_geometry);

            const double type2_volume = CalculateTetrahedronVolume(
                *r_edge_data.pEdgeNode0,
                *r_edge_data.pEdgeNode1,
                *p_projection_node_a,
                *p_projection_node_b);
            const bool create_type2_quadrature =
                type2_volume >= MinimumGapVolumeForQuadrature;

            if (create_type2_quadrature) {
                const IntegrationPointsArrayType type2_integration_points =
                    CreateCoonsVolumeGaussPoints(
                        IntegrationOrder,
                        *p_type2_volume);
                double type2_weight_sum = 0.0;
                for (const auto& r_integration_point : type2_integration_points) {
                    type2_weight_sum += r_integration_point.Weight();
                }

                total_type2_tetra_volume += type2_volume;
                total_type2_quadrature_weight += type2_weight_sum;
                max_type2_relative_weight_error = std::max(
                    max_type2_relative_weight_error,
                    std::abs(type2_weight_sum - type2_volume) /
                    std::max(type2_volume, MinimumGapVolumeForQuadrature));
            } else {
                ++result.Summary.NumberOfSkippedEdges;
            }

            auto orient_face_away_from_volume = [&volume_center](
                NodePointerType& rpNode0,
                NodePointerType& rpNode1,
                NodePointerType& rpNode2)
            {
                const array_1d<double, 3> edge_01 =
                    rpNode1->Coordinates() - rpNode0->Coordinates();
                const array_1d<double, 3> edge_02 =
                    rpNode2->Coordinates() - rpNode0->Coordinates();
                const array_1d<double, 3> face_normal =
                    MathUtils<double>::CrossProduct(edge_01, edge_02);

                array_1d<double, 3> face_center = ZeroVector(3);
                noalias(face_center) += rpNode0->Coordinates();
                noalias(face_center) += rpNode1->Coordinates();
                noalias(face_center) += rpNode2->Coordinates();
                face_center /= 3.0;

                const array_1d<double, 3> outward_vector =
                    face_center - volume_center;

                if (inner_prod(face_normal, outward_vector) < 0.0) {
                    std::swap(rpNode1, rpNode2);
                }
            };

            auto register_type2_face = [&](
                const IndexType SurrogateConditionId,
                const SpanKey3D& rOccurrenceSpan,
                const NodePointerType& pNode0,
                const NodePointerType& pNode1,
                const NodePointerType& pNode2,
                const NodePointerType& pOppositeNode,
                const Geometry<Node>::Pointer& pFaceNeighbourGeometry)
            {
                KRATOS_ERROR_IF_NOT(pFaceNeighbourGeometry)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Trying to register a type2 face with null neighbour geometry.\n";

                const auto face_key = MakeCanonicalFaceKey3D(
                    pNode0,
                    pNode1,
                    pNode2);

                const auto face_registry_it =
                    mLateralFaceRegistry.find(face_key);
                const std::size_t previous_occurrences =
                    face_registry_it == mLateralFaceRegistry.end() ?
                    0 : face_registry_it->second.size();

                KRATOS_ERROR_IF(previous_occurrences >= 2)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Trying to add a third occurrence to a type2 face.\n"
                    << "  face node ids: "
                    << face_key.NodeIds[0] << ", "
                    << face_key.NodeIds[1] << ", "
                    << face_key.NodeIds[2] << "\n"
                    << "  occurrences: " << previous_occurrences << "\n"
                    << "  external span: " << SpanToString(rOccurrenceSpan)
                    << "\n";

                NurbsSurfaceType::Pointer p_surface;
                if (previous_occurrences == 0) {
                    p_surface = CreateCollapsedTriangleSurface(
                        pNode0,
                        pNode1,
                        pNode2);
                    p_surface->SetId(next_geometry_id++);
                    AddNeighbourGeometry(*p_surface, pFaceNeighbourGeometry);
                    if (mStoreGapDebugGeometries) {
                        r_gap_type2_debug.AddGeometry(p_surface);
                    }
                } else {
                    p_surface = face_registry_it->second.front().pGeometry;
                    KRATOS_ERROR_IF_NOT(p_surface)
                        << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                        << "Existing type2 face occurrence has null geometry.\n";
                    AddUniqueNeighbourGeometry(
                        *p_surface,
                        pFaceNeighbourGeometry);
                }

                RegisterLateralFaceOccurrence(
                    2,
                    SurrogateConditionId,
                    rOccurrenceSpan,
                    pNode0,
                    pNode1,
                    pNode2,
                    pOppositeNode,
                    p_surface,
                    pFaceNeighbourGeometry,
                    selected_neighbour_has_type1_path,
                    selected_neighbour_has_active_span,
                    selected_neighbour_active_span);

                return std::make_pair(p_surface, previous_occurrences);
            };

            created_volume_keys.insert(volume_key);

            if (candidate_a_is_type1 && candidate_b_is_type1) {
                ++number_type1_type1_pairs;
            } else if (candidate_a_is_type1 || candidate_b_is_type1) {
                ++number_type1_type2_pairs;
            } else {
                ++number_type2_type2_pairs;
            }

            if (has_unbacked_type1_candidate) {
                ++number_pairs_with_unbacked_type1_candidate;
            }

            if (r_candidate_a.HasType1Face()) {
                rType1CreationResult.Type1LateralFaces[
                    r_candidate_a.Type1FaceIndex].IsClosedByType2 = true;
            }
            if (r_candidate_b.HasType1Face()) {
                rType1CreationResult.Type1LateralFaces[
                    r_candidate_b.Type1FaceIndex].IsClosedByType2 = true;
            }

            p_type2_volume->SetId(next_geometry_id++);
            for (const auto& p_current_neighbour_geometry : type2_neighbour_geometries) {
                AddNeighbourGeometry(
                    *p_type2_volume,
                    p_current_neighbour_geometry);
            }
            if (mStoreGapDebugGeometries) {
                r_gap_type2_debug.AddGeometry(p_type2_volume);
            }

            if (create_type2_quadrature) {
                Type2VolumeQuadratureData type2_volume_data;

                type2_volume_data.NeighbourGeometries =
                    type2_neighbour_geometries;

                type2_volume_data.CharacteristicLength =
                    std::max(
                        norm_2(r_edge_data.pEdgeNode0->Coordinates() - r_edge_data.pEdgeNode1->Coordinates()),
                        norm_2(p_projection_node_a->Coordinates() - p_projection_node_b->Coordinates()));

                type2_volume_data.EdgeKey = r_edge_data.EdgeKey;
                type2_volume_data.pEdgeNode0 = r_edge_data.pEdgeNode0;
                type2_volume_data.pEdgeNode1 = r_edge_data.pEdgeNode1;
                type2_volume_data.pProjectionNode0 = p_projection_node_a;
                type2_volume_data.pProjectionNode1 = p_projection_node_b;
                type2_volume_data.ProjectionNodeId0 = p_projection_node_a->Id();
                type2_volume_data.ProjectionNodeId1 = p_projection_node_b->Id();
                type2_volume_data.ExternalSpan0 = r_candidate_a.Span;
                type2_volume_data.ExternalSpan1 = r_candidate_b.Span;
                type2_volume_data.FirstType1LateralFaceIndex =
                    r_candidate_a.HasType1Face() ? r_candidate_a.Type1FaceIndex : 0;
                type2_volume_data.SecondType1LateralFaceIndex =
                    r_candidate_b.HasType1Face() ? r_candidate_b.Type1FaceIndex : 0;
                type2_volume_data.RequiresQuadrature = true;
                type2_volume_data.pVolumeGeometry = p_type2_volume;

                result.VolumeQuadratureDataList.emplace_back(
                    std::move(type2_volume_data));

                ++result.Summary.NumberOfCreatedVolumes;
            }

            // Open face at S0: S0 - Pa - Pb
            {
                auto p_face_node_0 = r_edge_data.pEdgeNode0;
                auto p_face_node_1 = p_projection_node_a;
                auto p_face_node_2 = p_projection_node_b;
                orient_face_away_from_volume(
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2);

                const auto registration_data = register_type2_face(
                    r_edge_data.SurrogateConditionId,
                    r_candidate_a.Span,
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2,
                    r_edge_data.pEdgeNode1,
                    p_neighbour_geometry);

                if (registration_data.second == 0) {
                    Type2OpenFaceData open_face_data;
                    open_face_data.pSurrogateNode = p_face_node_0;
                    open_face_data.pProjectionNode0 = p_projection_node_a;
                    open_face_data.pProjectionNode1 = p_projection_node_b;
                    open_face_data.pOppositeNode = r_edge_data.pEdgeNode1;
                    open_face_data.pSurfaceGeometry = registration_data.first;
                    open_face_data.NeighbourGeometries.clear();
                    open_face_data.NeighbourGeometries.push_back(p_neighbour_geometry);
                    open_face_data.ParentEdgeKey = r_edge_data.EdgeKey;
                    open_face_data.ExternalSpan0 = r_candidate_a.Span;
                    open_face_data.ExternalSpan1 = r_candidate_b.Span;

                    result.OpenFaceDataList.emplace_back(std::move(open_face_data));
                    ++result.Summary.NumberOfOpenFaces;
                }
            }

            // Open face at S1: S1 - Pb - Pa
            {
                auto p_face_node_0 = r_edge_data.pEdgeNode1;
                auto p_face_node_1 = p_projection_node_b;
                auto p_face_node_2 = p_projection_node_a;
                orient_face_away_from_volume(
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2);

                const auto registration_data = register_type2_face(
                    r_edge_data.SurrogateConditionId,
                    r_candidate_b.Span,
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2,
                    r_edge_data.pEdgeNode0,
                    p_neighbour_geometry);

                if (registration_data.second == 0) {
                    Type2OpenFaceData open_face_data;
                    open_face_data.pSurrogateNode = p_face_node_0;
                    open_face_data.pProjectionNode0 = p_projection_node_a;
                    open_face_data.pProjectionNode1 = p_projection_node_b;
                    open_face_data.pOppositeNode = r_edge_data.pEdgeNode0;
                    open_face_data.pSurfaceGeometry = registration_data.first;
                    open_face_data.NeighbourGeometries.clear();
                    open_face_data.NeighbourGeometries.push_back(p_neighbour_geometry);
                    open_face_data.ParentEdgeKey = r_edge_data.EdgeKey;
                    open_face_data.ExternalSpan0 = r_candidate_a.Span;
                    open_face_data.ExternalSpan1 = r_candidate_b.Span;

                    result.OpenFaceDataList.emplace_back(std::move(open_face_data));
                    ++result.Summary.NumberOfOpenFaces;
                }
            }

            {
                auto p_face_node_0 = p_s0;
                auto p_face_node_1 = p_s1;
                auto p_face_node_2 = p_pa;
                orient_face_away_from_volume(
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2);

                register_type2_face(
                    get_candidate_surrogate_condition_id(r_candidate_a),
                    r_candidate_a.Span,
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2,
                    p_pb,
                    p_neighbour_geometry);
            }
            {
                auto p_face_node_0 = p_s1;
                auto p_face_node_1 = p_s0;
                auto p_face_node_2 = p_pb;
                orient_face_away_from_volume(
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2);

                register_type2_face(
                    get_candidate_surrogate_condition_id(r_candidate_b),
                    r_candidate_b.Span,
                    p_face_node_0,
                    p_face_node_1,
                    p_face_node_2,
                    p_pa,
                    p_neighbour_geometry);
            }
            }
        }
    }

    std::size_t number_registry_open_type1 = 0;
    std::size_t number_registry_open_type2 = 0;
    std::size_t number_registry_closed_type1_type1 = 0;
    std::size_t number_registry_closed_type1_type2 = 0;
    std::size_t number_registry_closed_type2_type2 = 0;
    std::size_t number_registry_closed_same_neighbour = 0;
    std::size_t number_registry_closed_different_neighbours = 0;
    std::size_t number_registry_shared_faces_with_complete_neighbours = 0;

    for (auto& r_registry_entry : mLateralFaceRegistry) {
        auto& r_occurrences = r_registry_entry.second;
        if (r_occurrences.size() != 2) {
            continue;
        }

        for (const auto& r_target_occurrence : r_occurrences) {
            KRATOS_ERROR_IF_NOT(r_target_occurrence.pNeighbourGeometry)
                << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                << "Shared lateral face occurrence has null neighbour geometry.\n";

            for (const auto& r_surface_occurrence : r_occurrences) {
                KRATOS_ERROR_IF_NOT(r_surface_occurrence.pGeometry)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Shared lateral face occurrence has null surface geometry.\n";

                AddUniqueNeighbourGeometry(
                    *r_surface_occurrence.pGeometry,
                    r_target_occurrence.pNeighbourGeometry);
            }
        }
    }

    for (const auto& r_registry_entry : mLateralFaceRegistry) {
        const auto& r_occurrences = r_registry_entry.second;

        if (r_occurrences.size() == 1) {
            if (r_occurrences.front().GapType == 1) {
                ++number_registry_open_type1;
            } else if (r_occurrences.front().GapType == 2) {
                ++number_registry_open_type2;
            }
        } else if (r_occurrences.size() == 2) {
            const int gap_type_0 = r_occurrences[0].GapType;
            const int gap_type_1 = r_occurrences[1].GapType;

            for (const auto& r_occurrence : r_occurrences) {
                KRATOS_ERROR_IF(
                    r_occurrence.GapType == 2 &&
                    !r_occurrence.HasType1NeighbourPath)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Shared type2 lateral face occurrence has no "
                    << "continuous neighbour path to a type1 volume.\n"
                    << "  face node ids: "
                    << r_registry_entry.first.NodeIds[0] << ", "
                    << r_registry_entry.first.NodeIds[1] << ", "
                    << r_registry_entry.first.NodeIds[2] << "\n"
                    << "  external span: "
                    << SpanToString(r_occurrence.ExternalSpan) << "\n";

                KRATOS_ERROR_IF_NOT(r_occurrence.pGeometry)
                    << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                    << "Shared lateral face occurrence has null surface geometry.\n";

                const auto stored_neighbour_geometries =
                    r_occurrence.pGeometry->Has(NEIGHBOUR_GEOMETRIES) ?
                    r_occurrence.pGeometry->GetValue(NEIGHBOUR_GEOMETRIES) :
                    NeighbourGeometriesVectorType();

                for (const auto& r_other_occurrence : r_occurrences) {
                    KRATOS_ERROR_IF_NOT(ContainsNeighbourGeometry(
                        stored_neighbour_geometries,
                        r_other_occurrence.pNeighbourGeometry))
                        << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
                        << "Shared lateral face geometry does not contain "
                        << "both occurrence neighbour geometries.\n"
                        << "  face node ids: "
                        << r_registry_entry.first.NodeIds[0] << ", "
                        << r_registry_entry.first.NodeIds[1] << ", "
                        << r_registry_entry.first.NodeIds[2] << "\n"
                        << "  checked occurrence external span: "
                        << SpanToString(r_occurrence.ExternalSpan) << "\n"
                        << "  missing neighbour from external span: "
                        << SpanToString(r_other_occurrence.ExternalSpan) << "\n";
                }
            }

            ++number_registry_shared_faces_with_complete_neighbours;

            if (gap_type_0 == 1 && gap_type_1 == 1) {
                ++number_registry_closed_type1_type1;
            } else if ((gap_type_0 == 1 && gap_type_1 == 2) ||
                       (gap_type_0 == 2 && gap_type_1 == 1)) {
                ++number_registry_closed_type1_type2;
            } else if (gap_type_0 == 2 && gap_type_1 == 2) {
                ++number_registry_closed_type2_type2;
            }

            const auto neighbour_geometries =
                CollectUniqueNeighbourGeometries(r_occurrences);
            if (neighbour_geometries.size() <= 1) {
                ++number_registry_closed_same_neighbour;
            } else {
                ++number_registry_closed_different_neighbours;
            }
        }
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type 2 creation summary:\n"
        << "  candidate active edges: " << result.Summary.NumberOfCandidateEdges << "\n"
        << "  created volumes: " << result.Summary.NumberOfCreatedVolumes << "\n"
        << "  skipped candidates: " << result.Summary.NumberOfSkippedEdges << "\n"
        << "  open faces: " << result.Summary.NumberOfOpenFaces << "\n"
        << "  created type1/type1 pairs: " << number_type1_type1_pairs << "\n"
        << "  created type1/type2 pairs: " << number_type1_type2_pairs << "\n"
        << "  created type2/type2 pairs: " << number_type2_type2_pairs << "\n"
        << "  pairs with type1 candidate without lateral face on edge: "
        << number_pairs_with_unbacked_type1_candidate << "\n"
        << "  skipped pairs with type1 candidate without lateral face on edge: "
        << number_skipped_unbacked_type1_candidate_pairs << "\n"
        << "  total type2 tetra volume: " << total_type2_tetra_volume << "\n"
        << "  total type2 quadrature weight: " << total_type2_quadrature_weight << "\n"
        << "  max type2 relative weight error: " << max_type2_relative_weight_error << "\n"
        << "  pairs with common active span: "
        << number_pairs_with_common_active_span << "\n"
        << "  pairs without common active span: "
        << number_pairs_without_common_active_span << "\n"
        << "  type1/type1 pairs without common active span: "
        << number_type1_type1_pairs_without_common_active_span << "\n"
        << "  type1/type2 pairs without common active span: "
        << number_type1_type2_pairs_without_common_active_span << "\n"
        << "  type2 neighbours selected from reused faces: "
        << number_type2_neighbour_selected_from_reused_face << "\n"
        << "  type2 reused lateral face candidates: "
        << number_type2_reused_lateral_face_candidates << "\n"
        << "  type2 reused lateral face candidates with common active span: "
        << number_type2_reused_lateral_face_candidates_with_common_active_span << "\n"
        << "  type2 neighbours selected from S0-Pa-Pb faces: "
        << number_type2_neighbour_selected_from_s0_projection_face << "\n"
        << "  type2 neighbours selected from S1-Pb-Pa faces: "
        << number_type2_neighbour_selected_from_s1_projection_face << "\n"
        << "  type2 neighbours selected from S0-S1-Pa faces: "
        << number_type2_neighbour_selected_from_candidate_a_edge_face << "\n"
        << "  type2 neighbours selected from S1-S0-Pb faces: "
        << number_type2_neighbour_selected_from_candidate_b_edge_face << "\n"
        << "  registry open type1 faces: " << number_registry_open_type1 << "\n"
        << "  registry open type2 faces: " << number_registry_open_type2 << "\n"
        << "  registry closed type1/type1 faces: "
        << number_registry_closed_type1_type1 << "\n"
        << "  registry closed type1/type2 faces: "
        << number_registry_closed_type1_type2 << "\n"
        << "  registry closed type2/type2 faces: "
        << number_registry_closed_type2_type2 << "\n"
        << "  registry closed same-neighbour faces: "
        << number_registry_closed_same_neighbour << "\n"
        << "  registry closed different-neighbour faces: "
        << number_registry_closed_different_neighbours << "\n"
        << "  registry shared faces with complete neighbours: "
        << number_registry_shared_faces_with_complete_neighbours << "\n"
        << "  GapType2Debug geometries: " << r_gap_type2_debug.NumberOfGeometries() << "\n";

    KRATOS_ERROR_IF(number_registry_open_type1 != 0)
        << "[SnakeGapSbm3DUtilities::CreateType2GapGeometries] "
        << "Type2 creation left open type1 lateral faces. "
        << "Every type1 lateral face must be closed by a type2 volume.\n"
        << "  registry open type1 faces: " << number_registry_open_type1 << "\n";

    return result;
}

bool SnakeGapSbm3DUtilities::HaveCommonActiveSpan(
    const std::vector<SpanKey3D>& rFirst,
    const std::vector<SpanKey3D>& rSecond,
    const std::vector<SpanKey3D>& rThird) const
{
    for (const auto& r_active_span : rFirst) {
        if (std::find(rSecond.begin(), rSecond.end(), r_active_span) !=
                rSecond.end() &&
            std::find(rThird.begin(), rThird.end(), r_active_span) !=
                rThird.end()) {
            return true;
        }
    }

    return false;
}

#if 1
SnakeGapSbm3DUtilities::Type3CreationResult
SnakeGapSbm3DUtilities::CreateType3GapGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    const ExternalSpanDataMap& rExternalSpans,
    const KnotSpanGridInfo& rGridInfo,
    const Type2CreationResult& rType2CreationResult,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives)
{
    static_cast<void>(IntegrationOrder);
    static_cast<void>(NumberOfShapeFunctionsDerivatives);

    Type3CreationResult result;
    result.Summary.NumberOfCandidateOpenFaces =
        rType2CreationResult.OpenFaceDataList.size();

    ModelPart& r_gap_type3_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType3Debug");
    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    std::map<IndexType, NodePointerType> surrogate_nodes_by_id;
    std::map<IndexType, LateralFaceOccurrence>
        reference_occurrence_by_surrogate_node_id;

    for (const auto& r_open_face_data : rType2CreationResult.OpenFaceDataList) {
        KRATOS_ERROR_IF_NOT(r_open_face_data.pSurrogateNode)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face has null surrogate node.\n";

        surrogate_nodes_by_id.emplace(
            r_open_face_data.pSurrogateNode->Id(),
            r_open_face_data.pSurrogateNode);

        const auto face_key = MakeCanonicalFaceKey3D(
            r_open_face_data.pSurrogateNode,
            r_open_face_data.pProjectionNode0,
            r_open_face_data.pProjectionNode1);
        const auto registry_it = mLateralFaceRegistry.find(face_key);
        if (registry_it != mLateralFaceRegistry.end() &&
            registry_it->second.size() == 1 &&
            registry_it->second.front().GapType == 2) {
            reference_occurrence_by_surrogate_node_id.emplace(
                r_open_face_data.pSurrogateNode->Id(),
                registry_it->second.front());
        }
    }

    struct DirectionalExternalSpans {
        std::array<std::vector<SpanKey3D>, 6> Spans;

        std::vector<SpanKey3D>& XNegative() { return Spans[0]; }
        std::vector<SpanKey3D>& XPositive() { return Spans[1]; }
        std::vector<SpanKey3D>& YNegative() { return Spans[2]; }
        std::vector<SpanKey3D>& YPositive() { return Spans[3]; }
        std::vector<SpanKey3D>& ZNegative() { return Spans[4]; }
        std::vector<SpanKey3D>& ZPositive() { return Spans[5]; }
    };

    struct Type3VolumeKey {
        IndexType SurrogateNodeId = 0;
        std::array<IndexType, 3> ProjectionNodeIds = {{0, 0, 0}};

        bool operator<(const Type3VolumeKey& rOther) const
        {
            return std::tie(SurrogateNodeId, ProjectionNodeIds) <
                   std::tie(rOther.SurrogateNodeId, rOther.ProjectionNodeIds);
        }

        static Type3VolumeKey Create(
            const NodePointerType& pSurrogateNode,
            const NodePointerType& pProjectionNode0,
            const NodePointerType& pProjectionNode1,
            const NodePointerType& pProjectionNode2)
        {
            Type3VolumeKey key;
            key.SurrogateNodeId = pSurrogateNode->Id();
            key.ProjectionNodeIds = {{
                pProjectionNode0->Id(),
                pProjectionNode1->Id(),
                pProjectionNode2->Id()}};
            std::sort(
                key.ProjectionNodeIds.begin(),
                key.ProjectionNodeIds.end());
            return key;
        }
    };

    std::set<Type3VolumeKey> created_volume_keys;
    std::size_t number_skipped_degenerate_candidate_volumes = 0;
    std::size_t number_skipped_empty_quadrature_volumes = 0;

    const std::array<const char*, 6> direction_names = {{
        "x negative",
        "x positive",
        "y negative",
        "y positive",
        "z negative",
        "z positive"}};

    auto span_list_to_string = [this](
        const std::vector<SpanKey3D>& rSpans)
    {
        std::ostringstream buffer;
        for (const auto& r_span : rSpans) {
            buffer << " " << SpanToString(r_span);
        }
        return buffer.str();
    };

    auto span_centers_share_side = [&](
        const SpanKey3D& rSpan0,
        const SpanKey3D& rSpan1)
    {
        const array_1d<double, 3> center_0 = SpanCenter(rSpan0, rGridInfo);
        const array_1d<double, 3> center_1 = SpanCenter(rSpan1, rGridInfo);
        std::size_t number_of_equal_coordinates = 0;
        for (IndexType axis = 0; axis < 3; ++axis) {
            if (center_0[axis] == center_1[axis]) {
                ++number_of_equal_coordinates;
            }
        }

        return number_of_equal_coordinates == 2;
    };

    auto calculate_characteristic_length = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const NodePointerType& pNode3)
    {
        const std::array<NodePointerType, 4> nodes = {{
            pNode0,
            pNode1,
            pNode2,
            pNode3}};

        double characteristic_length = 0.0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            for (std::size_t j = i + 1; j < nodes.size(); ++j) {
                characteristic_length = std::max(
                    characteristic_length,
                    norm_2(nodes[i]->Coordinates() - nodes[j]->Coordinates()));
            }
        }

        return characteristic_length;
    };

    auto signed_face_side = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const NodePointerType& pPoint)
    {
        const array_1d<double, 3> edge_01 =
            pNode1->Coordinates() - pNode0->Coordinates();
        const array_1d<double, 3> edge_02 =
            pNode2->Coordinates() - pNode0->Coordinates();
        const array_1d<double, 3> face_normal =
            MathUtils<double>::CrossProduct(edge_01, edge_02);

        return inner_prod(face_normal, pPoint->Coordinates() - pNode0->Coordinates());
    };

    auto face_side_tolerance = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const NodePointerType& pPoint0,
        const NodePointerType& pPoint1)
    {
        double characteristic_length = std::max({
            norm_2(pNode0->Coordinates() - pNode1->Coordinates()),
            norm_2(pNode1->Coordinates() - pNode2->Coordinates()),
            norm_2(pNode2->Coordinates() - pNode0->Coordinates()),
            norm_2(pPoint0->Coordinates() - pNode0->Coordinates()),
            norm_2(pPoint1->Coordinates() - pNode0->Coordinates())});

        return 1.0e-12 * std::max(
            1.0,
            characteristic_length * characteristic_length * characteristic_length);
    };

    for (const auto& r_surrogate_node_entry : surrogate_nodes_by_id) {
        const IndexType surrogate_node_id = r_surrogate_node_entry.first;
        const NodePointerType p_surrogate_node = r_surrogate_node_entry.second;
        const auto& r_surrogate_coordinates = p_surrogate_node->Coordinates();
        const GridPointKey3D surrogate_grid_node =
            ComputeGridPointKey(r_surrogate_coordinates, rGridInfo);

        std::vector<SpanKey3D> incident_external_spans;
        for (int di = -1; di <= 0; ++di) {
            for (int dj = -1; dj <= 0; ++dj) {
                for (int dk = -1; dk <= 0; ++dk) {
                    const SpanKey3D incident_span{
                        surrogate_grid_node.I + di,
                        surrogate_grid_node.J + dj,
                        surrogate_grid_node.K + dk};

                    if (!IsValidSpan(incident_span, rGridInfo)) {
                        continue;
                    }

                    if (rExternalSpans.find(incident_span) !=
                        rExternalSpans.end()) {
                        incident_external_spans.push_back(incident_span);
                    }
                }
            }
        }

        if (incident_external_spans.size() < 4) {
            continue;
        }

        DirectionalExternalSpans directional_spans;
        for (const auto& r_span : incident_external_spans) {
            const array_1d<double, 3> center =
                SpanCenter(r_span, rGridInfo);

            if (center[0] < r_surrogate_coordinates[0]) {
                directional_spans.XNegative().push_back(r_span);
            } else if (center[0] > r_surrogate_coordinates[0]) {
                directional_spans.XPositive().push_back(r_span);
            }

            if (center[1] < r_surrogate_coordinates[1]) {
                directional_spans.YNegative().push_back(r_span);
            } else if (center[1] > r_surrogate_coordinates[1]) {
                directional_spans.YPositive().push_back(r_span);
            }

            if (center[2] < r_surrogate_coordinates[2]) {
                directional_spans.ZNegative().push_back(r_span);
            } else if (center[2] > r_surrogate_coordinates[2]) {
                directional_spans.ZPositive().push_back(r_span);
            }
        }

        std::size_t number_of_full_direction_vectors = 0;
        for (std::size_t direction_index = 0;
             direction_index < directional_spans.Spans.size();
             ++direction_index) {
            const auto& r_direction_spans =
                directional_spans.Spans[direction_index];

            KRATOS_ERROR_IF(r_direction_spans.size() > 4)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "More than four external knot spans found in one "
                << "surrogate-node direction.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  direction: " << direction_names[direction_index] << "\n"
                << "  number of spans: " << r_direction_spans.size() << "\n"
                << "  spans:" << span_list_to_string(r_direction_spans)
                << "\n";

            if (r_direction_spans.size() == 4) {
                ++number_of_full_direction_vectors;
            }
        }

        if (number_of_full_direction_vectors == 0) {
            continue;
        }

        // if (number_of_full_direction_vectors > 2) continue; //FIXME: useful for debug

        KRATOS_ERROR_IF(number_of_full_direction_vectors > 3)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "More than three surrogate-node direction vectors contain "
            << "four external knot spans.\n"
            << "  surrogate node: " << surrogate_node_id << "\n"
            << "  number of full direction vectors: "
            << number_of_full_direction_vectors << "\n";

        const auto reference_occurrence_it =
            reference_occurrence_by_surrogate_node_id.find(surrogate_node_id);
        KRATOS_ERROR_IF(
            reference_occurrence_it ==
            reference_occurrence_by_surrogate_node_id.end())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Cannot create type3 closure because no reference type2 "
            << "occurrence was found for the surrogate node.\n"
            << "  surrogate node: " << surrogate_node_id << "\n";

        const auto& r_reference_occurrence =
            reference_occurrence_it->second;
        KRATOS_ERROR_IF_NOT(r_reference_occurrence.pNeighbourGeometry)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Reference type2 occurrence for type3 creation has null "
            << "neighbour geometry.\n"
            << "  surrogate node: " << surrogate_node_id << "\n";
        KRATOS_ERROR_IF_NOT(r_reference_occurrence.HasType1NeighbourPath)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Reference type2 occurrence for type3 creation has no "
            << "continuous neighbour path to a type1 volume.\n"
            << "  surrogate node: " << surrogate_node_id << "\n";
        KRATOS_ERROR_IF_NOT(r_reference_occurrence.HasNeighbourActiveSpan)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Reference type2 occurrence for type3 creation has no "
            << "active span provenance.\n"
            << "  surrogate node: " << surrogate_node_id << "\n";

        auto direction_has_registered_face = [&](
            const std::size_t DirectionIndex)
        {
            const auto& r_direction_spans =
                directional_spans.Spans[DirectionIndex];
            if (r_direction_spans.size() != 4) {
                return false;
            }

            std::array<NodePointerType, 4> direction_projection_nodes;
            for (std::size_t i = 0; i < r_direction_spans.size(); ++i) {
                const auto span_it = rExternalSpans.find(r_direction_spans[i]);
                KRATOS_ERROR_IF(span_it == rExternalSpans.end())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Selected external span missing from external span map "
                    << "while ordering type3 directions.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[DirectionIndex] << "\n"
                    << "  span: " << SpanToString(r_direction_spans[i]) << "\n";
                KRATOS_ERROR_IF_NOT(span_it->second.HasProjectionNode())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Selected external span has no projection node while "
                    << "ordering type3 directions.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[DirectionIndex] << "\n"
                    << "  span: " << SpanToString(r_direction_spans[i]) << "\n";

                direction_projection_nodes[i] = span_it->second.pProjectionNode;
            }

            for (std::size_t i = 0; i < r_direction_spans.size(); ++i) {
                for (std::size_t j = i + 1; j < r_direction_spans.size(); ++j) {
                    if (!span_centers_share_side(
                            r_direction_spans[i],
                            r_direction_spans[j])) {
                        continue;
                    }

                    const auto face_key = MakeCanonicalFaceKey3D(
                        p_surrogate_node,
                        direction_projection_nodes[i],
                        direction_projection_nodes[j]);
                    if (mLateralFaceRegistry.find(face_key) !=
                        mLateralFaceRegistry.end()) {
                        return true;
                    }
                }
            }

            return false;
        };

        std::vector<std::size_t> full_direction_indices;
        full_direction_indices.reserve(number_of_full_direction_vectors);
        for (std::size_t direction_index = 0;
             direction_index < directional_spans.Spans.size();
             ++direction_index) {
            if (directional_spans.Spans[direction_index].size() == 4) {
                full_direction_indices.push_back(direction_index);
            }
        }

        std::vector<std::size_t> remaining_direction_indices =
            full_direction_indices;

        while (!remaining_direction_indices.empty()) {
            const auto registered_direction_it = std::find_if(
                remaining_direction_indices.begin(),
                remaining_direction_indices.end(),
                direction_has_registered_face);

            KRATOS_ERROR_IF(
                registered_direction_it == remaining_direction_indices.end())
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "No remaining full type3 direction has a registered "
                << "lateral face. Cannot assign a neighbour geometry to the "
                << "next type3 volume without inventing one.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  remaining full direction vectors: "
                << remaining_direction_indices.size() << "\n";

            const auto direction_index = *registered_direction_it;
            remaining_direction_indices.erase(registered_direction_it);

            const auto& r_selected_direction_spans =
                directional_spans.Spans[direction_index];

            std::array<const ExternalSpanData*, 4> p_span_data = {{
                nullptr, nullptr, nullptr, nullptr}};
            std::array<NodePointerType, 4> projection_nodes;
            for (std::size_t i = 0; i < r_selected_direction_spans.size(); ++i) {
                const auto span_it =
                    rExternalSpans.find(r_selected_direction_spans[i]);
                KRATOS_ERROR_IF(span_it == rExternalSpans.end())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Selected external span missing from external span map.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n"
                    << "  span: " << SpanToString(r_selected_direction_spans[i]) << "\n";
                KRATOS_ERROR_IF_NOT(span_it->second.HasProjectionNode())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Selected external span has no projection node.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n"
                    << "  span: " << SpanToString(r_selected_direction_spans[i]) << "\n";

                p_span_data[i] = &span_it->second;
                projection_nodes[i] = span_it->second.pProjectionNode;
            }

            std::array<std::vector<std::size_t>, 4> adjacency;
            std::size_t number_of_cycle_edges = 0;
            for (std::size_t i = 0; i < r_selected_direction_spans.size(); ++i) {
                for (std::size_t j = i + 1; j < r_selected_direction_spans.size(); ++j) {
                    if (!span_centers_share_side(
                        r_selected_direction_spans[i],
                        r_selected_direction_spans[j])) {
                        continue;
                    }

                    adjacency[i].push_back(j);
                    adjacency[j].push_back(i);
                    ++number_of_cycle_edges;

                    const auto face_key = MakeCanonicalFaceKey3D(
                        p_surrogate_node,
                        projection_nodes[i],
                        projection_nodes[j]);
                    const auto face_registry_it =
                        mLateralFaceRegistry.find(face_key);

                    if (face_registry_it != mLateralFaceRegistry.end()) {
                        KRATOS_ERROR_IF(face_registry_it->second.size() > 2)
                            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                            << "Projection-projection-surrogate face already exists "
                            << "with more than two occurrences.\n"
                            << "  surrogate node: " << surrogate_node_id << "\n"
                            << "  direction: " << direction_names[direction_index] << "\n"
                            << "  span 0: " << SpanToString(r_selected_direction_spans[i]) << "\n"
                            << "  span 1: " << SpanToString(r_selected_direction_spans[j]) << "\n"
                            << "  face node ids: "
                            << face_key.NodeIds[0] << ", "
                            << face_key.NodeIds[1] << ", "
                            << face_key.NodeIds[2] << "\n"
                            << "  occurrences: "
                            << face_registry_it->second.size() << "\n";
                    }
                }
            }

            KRATOS_ERROR_IF(number_of_cycle_edges != 4)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "The four external knot spans in a direction do not define "
                << "a simple four-edge cycle.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  direction: " << direction_names[direction_index] << "\n"
                << "  number of side-sharing pairs: "
                << number_of_cycle_edges << "\n"
                << "  spans:" << span_list_to_string(r_selected_direction_spans)
                << "\n";
            for (std::size_t i = 0; i < adjacency.size(); ++i) {
                KRATOS_ERROR_IF(adjacency[i].size() != 2)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "The four external knot spans in a direction are not a "
                    << "simple cycle: every span must have two side neighbours.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n"
                    << "  span: " << SpanToString(r_selected_direction_spans[i]) << "\n"
                    << "  neighbours: " << adjacency[i].size() << "\n";
            }

            std::array<std::size_t, 4> ordered_indices = {{0, 0, 0, 0}};
            ordered_indices[0] = 0;
            ordered_indices[1] = adjacency[0][0];
            for (std::size_t order_index = 2; order_index < 4; ++order_index) {
                const std::size_t previous_index =
                    ordered_indices[order_index - 2];
                const std::size_t current_index =
                    ordered_indices[order_index - 1];
                ordered_indices[order_index] =
                    adjacency[current_index][0] == previous_index ?
                    adjacency[current_index][1] :
                    adjacency[current_index][0];
            }
            KRATOS_ERROR_IF(
                adjacency[ordered_indices[3]][0] != ordered_indices[0] &&
                adjacency[ordered_indices[3]][1] != ordered_indices[0])
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Ordered side-sharing spans do not close into a cycle.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  direction: " << direction_names[direction_index] << "\n";

            using Type3TetraProjectionIndices =
                std::array<std::size_t, 3>;
            using Type3DirectionTriangulation =
                std::array<Type3TetraProjectionIndices, 2>;

            const std::array<Type3DirectionTriangulation, 2>
                triangulation_candidates = {{
                    Type3DirectionTriangulation{{
                        Type3TetraProjectionIndices{{ordered_indices[0], ordered_indices[1], ordered_indices[2]}},
                        Type3TetraProjectionIndices{{ordered_indices[0], ordered_indices[2], ordered_indices[3]}}}},
                    Type3DirectionTriangulation{{
                        Type3TetraProjectionIndices{{ordered_indices[0], ordered_indices[1], ordered_indices[3]}},
                        Type3TetraProjectionIndices{{ordered_indices[1], ordered_indices[2], ordered_indices[3]}}}}}};

            auto existing_face_is_compatible = [&](
                const NodePointerType& pNode0,
                const NodePointerType& pNode1,
                const NodePointerType& pNode2,
                const NodePointerType& pCandidateOppositeNode,
                std::size_t& rNumberOfExistingFaces)
            {
                const auto face_key = MakeCanonicalFaceKey3D(
                    pNode0,
                    pNode1,
                    pNode2);
                const auto face_registry_it =
                    mLateralFaceRegistry.find(face_key);
                if (face_registry_it == mLateralFaceRegistry.end()) {
                    return true;
                }

                ++rNumberOfExistingFaces;

                if (face_registry_it->second.size() >= 2) {
                    return false;
                }

                for (const auto& r_occurrence : face_registry_it->second) {
                    KRATOS_ERROR_IF_NOT(r_occurrence.pOppositeNode)
                        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                        << "Registered lateral face occurrence has no opposite "
                        << "node, so type3 cannot determine the already occupied "
                        << "side of the face.\n"
                        << "  surrogate node: " << surrogate_node_id << "\n"
                        << "  direction: " << direction_names[direction_index] << "\n"
                        << "  face node ids: "
                        << face_key.NodeIds[0] << ", "
                        << face_key.NodeIds[1] << ", "
                        << face_key.NodeIds[2] << "\n";

                    const double existing_side = signed_face_side(
                        pNode0,
                        pNode1,
                        pNode2,
                        r_occurrence.pOppositeNode);
                    const double candidate_side = signed_face_side(
                        pNode0,
                        pNode1,
                        pNode2,
                        pCandidateOppositeNode);
                    const double side_tolerance = face_side_tolerance(
                        pNode0,
                        pNode1,
                        pNode2,
                        r_occurrence.pOppositeNode,
                        pCandidateOppositeNode);

                    if (std::abs(existing_side) > side_tolerance &&
                        std::abs(candidate_side) > side_tolerance &&
                        existing_side * candidate_side > 0.0) {
                        return false;
                    }
                }

                return true;
            };

            auto tetra_existing_face_count = [&](
                const Type3TetraProjectionIndices& rTetraIndices,
                bool& rIsCompatible)
            {
                rIsCompatible = true;
                std::size_t number_of_existing_faces = 0;

                auto check_face = [&](
                    const NodePointerType& pNode0,
                    const NodePointerType& pNode1,
                    const NodePointerType& pNode2,
                    const NodePointerType& pCandidateOppositeNode)
                {
                    if (!existing_face_is_compatible(
                            pNode0,
                            pNode1,
                            pNode2,
                            pCandidateOppositeNode,
                            number_of_existing_faces)) {
                        rIsCompatible = false;
                    }
                };

                const auto& p_projection_node_0 =
                    projection_nodes[rTetraIndices[0]];
                const auto& p_projection_node_1 =
                    projection_nodes[rTetraIndices[1]];
                const auto& p_projection_node_2 =
                    projection_nodes[rTetraIndices[2]];

                check_face(
                    p_surrogate_node,
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2);
                check_face(
                    p_surrogate_node,
                    p_projection_node_1,
                    p_projection_node_2,
                    p_projection_node_0);
                check_face(
                    p_surrogate_node,
                    p_projection_node_2,
                    p_projection_node_0,
                    p_projection_node_1);
                check_face(
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2,
                    p_surrogate_node);

                return number_of_existing_faces;
            };

            auto triangulation_has_valid_internal_face = [&](
                const Type3DirectionTriangulation& rTriangulation)
            {
                std::array<std::size_t, 2> common_indices = {{0, 0}};
                std::size_t number_of_common_indices = 0;
                std::size_t opposite_index_0 = 0;
                std::size_t opposite_index_1 = 0;

                for (const auto index_0 : rTriangulation[0]) {
                    const bool is_common =
                        std::find(
                            rTriangulation[1].begin(),
                            rTriangulation[1].end(),
                            index_0) != rTriangulation[1].end();
                    if (is_common) {
                        common_indices[number_of_common_indices++] = index_0;
                    } else {
                        opposite_index_0 = index_0;
                    }
                }
                for (const auto index_1 : rTriangulation[1]) {
                    if (std::find(
                            rTriangulation[0].begin(),
                            rTriangulation[0].end(),
                            index_1) == rTriangulation[0].end()) {
                        opposite_index_1 = index_1;
                        break;
                    }
                }

                KRATOS_ERROR_IF(number_of_common_indices != 2)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Invalid type3 triangulation candidate: the two tetrahedra "
                    << "must share exactly two projection nodes.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n";

                const auto& p_common_node_0 =
                    projection_nodes[common_indices[0]];
                const auto& p_common_node_1 =
                    projection_nodes[common_indices[1]];
                const auto& p_opposite_node_0 =
                    projection_nodes[opposite_index_0];
                const auto& p_opposite_node_1 =
                    projection_nodes[opposite_index_1];

                const double side_0 = signed_face_side(
                    p_surrogate_node,
                    p_common_node_0,
                    p_common_node_1,
                    p_opposite_node_0);
                const double side_1 = signed_face_side(
                    p_surrogate_node,
                    p_common_node_0,
                    p_common_node_1,
                    p_opposite_node_1);
                const double side_tolerance = face_side_tolerance(
                    p_surrogate_node,
                    p_common_node_0,
                    p_common_node_1,
                    p_opposite_node_0,
                    p_opposite_node_1);

                return std::abs(side_0) <= side_tolerance ||
                       std::abs(side_1) <= side_tolerance ||
                       side_0 * side_1 < 0.0;
            };

            Type3DirectionTriangulation tetra_projection_indices;
            bool found_valid_triangulation = false;
            std::size_t best_triangulation_existing_face_count = 0;

            for (auto triangulation_candidate : triangulation_candidates) {
                if (!triangulation_has_valid_internal_face(triangulation_candidate)) {
                    continue;
                }

                bool first_tetra_is_compatible = false;
                bool second_tetra_is_compatible = false;
                const std::size_t first_existing_face_count =
                    tetra_existing_face_count(
                        triangulation_candidate[0],
                        first_tetra_is_compatible);
                const std::size_t second_existing_face_count =
                    tetra_existing_face_count(
                        triangulation_candidate[1],
                        second_tetra_is_compatible);

                if (!first_tetra_is_compatible ||
                    !second_tetra_is_compatible) {
                    continue;
                }

                if (first_existing_face_count == 0 &&
                    second_existing_face_count > 0) {
                    std::swap(
                        triangulation_candidate[0],
                        triangulation_candidate[1]);
                }

                if (first_existing_face_count == 0 &&
                    second_existing_face_count == 0) {
                    continue;
                }

                const std::size_t total_existing_face_count =
                    first_existing_face_count + second_existing_face_count;
                if (!found_valid_triangulation ||
                    total_existing_face_count >
                    best_triangulation_existing_face_count) {
                    tetra_projection_indices = triangulation_candidate;
                    best_triangulation_existing_face_count =
                        total_existing_face_count;
                    found_valid_triangulation = true;
                }
            }

            KRATOS_ERROR_IF_NOT(found_valid_triangulation)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Could not find a type3 triangulation that stays on the "
                << "free side of already registered faces.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  direction: " << direction_names[direction_index] << "\n"
                << "  spans:" << span_list_to_string(r_selected_direction_spans)
                << "\n";

            auto register_type3_face_if_needed = [&](
                const NodePointerType& pNode0,
                const NodePointerType& pNode1,
                const NodePointerType& pNode2,
                const NodePointerType& pOppositeNode,
                const Geometry<Node>::Pointer& pFaceNeighbourGeometry,
                const bool HasType1NeighbourPath,
                const bool HasNeighbourActiveSpan,
                const SpanKey3D& rNeighbourActiveSpan)
            {
                const auto face_key = MakeCanonicalFaceKey3D(
                    pNode0,
                    pNode1,
                    pNode2);
                const auto face_registry_it =
                    mLateralFaceRegistry.find(face_key);

                if (face_registry_it != mLateralFaceRegistry.end()) {
                    KRATOS_ERROR_IF(face_registry_it->second.size() > 2)
                        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                        << "Non-manifold type3 closure face.\n"
                        << "  surrogate node: " << surrogate_node_id << "\n"
                        << "  direction: " << direction_names[direction_index] << "\n"
                        << "  face node ids: "
                        << face_key.NodeIds[0] << ", "
                        << face_key.NodeIds[1] << ", "
                        << face_key.NodeIds[2] << "\n"
                        << "  occurrences: "
                        << face_registry_it->second.size() << "\n";

                    if (face_registry_it->second.size() == 2) {
                        return;
                    }
                }

                RegisterType3LateralFace(
                    r_gap_type3_debug,
                    next_geometry_id,
                    r_reference_occurrence.SurrogateConditionId,
                    r_reference_occurrence.ExternalSpan,
                    pNode0,
                    pNode1,
                    pNode2,
                    pOppositeNode,
                    pFaceNeighbourGeometry,
                    HasType1NeighbourPath,
                    HasNeighbourActiveSpan,
                    rNeighbourActiveSpan);
            };

            auto register_type3_top_face_if_needed = [&](
                const NodePointerType& pNode0,
                const NodePointerType& pNode1,
                const NodePointerType& pNode2,
                const NodePointerType& pOppositeNode,
                const NurbsSurfaceType::Pointer& pTopSurface,
                const Geometry<Node>::Pointer& pFaceNeighbourGeometry,
                const bool HasType1NeighbourPath,
                const bool HasNeighbourActiveSpan,
                const SpanKey3D& rNeighbourActiveSpan)
            {
                KRATOS_ERROR_IF_NOT(pTopSurface)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Type3 top face has null curved surface.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n";

                const auto face_key = MakeCanonicalFaceKey3D(
                    pNode0,
                    pNode1,
                    pNode2);
                const auto face_registry_it =
                    mLateralFaceRegistry.find(face_key);

                if (face_registry_it != mLateralFaceRegistry.end()) {
                    KRATOS_ERROR_IF(face_registry_it->second.size() > 2)
                        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                        << "Non-manifold type3 top closure face.\n"
                        << "  surrogate node: " << surrogate_node_id << "\n"
                        << "  direction: " << direction_names[direction_index] << "\n"
                        << "  face node ids: "
                        << face_key.NodeIds[0] << ", "
                        << face_key.NodeIds[1] << ", "
                        << face_key.NodeIds[2] << "\n"
                        << "  occurrences: "
                        << face_registry_it->second.size() << "\n";

                    if (face_registry_it->second.size() == 2) {
                        return;
                    }

                    for (const auto& r_occurrence : face_registry_it->second) {
                        KRATOS_ERROR_IF_NOT(r_occurrence.pGeometry)
                            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                            << "Existing type3 top face occurrence has null "
                            << "surface geometry.\n"
                            << "  surrogate node: " << surrogate_node_id << "\n"
                            << "  face node ids: "
                            << face_key.NodeIds[0] << ", "
                            << face_key.NodeIds[1] << ", "
                            << face_key.NodeIds[2] << "\n";
                    }
                }

                pTopSurface->SetId(next_geometry_id++);
                AddNeighbourGeometry(*pTopSurface, pFaceNeighbourGeometry);
                if (mStoreGapDebugGeometries) {
                    r_gap_type3_debug.AddGeometry(pTopSurface);
                }

                RegisterLateralFaceOccurrence(
                    3,
                    r_reference_occurrence.SurrogateConditionId,
                    r_reference_occurrence.ExternalSpan,
                    pNode0,
                    pNode1,
                    pNode2,
                    pOppositeNode,
                    pTopSurface,
                    pFaceNeighbourGeometry,
                    HasType1NeighbourPath,
                    HasNeighbourActiveSpan,
                    rNeighbourActiveSpan);
            };

            for (const auto& r_tetra_indices : tetra_projection_indices) {
                NodePointerType p_projection_node_0 =
                    projection_nodes[r_tetra_indices[0]];
                NodePointerType p_projection_node_1 =
                    projection_nodes[r_tetra_indices[1]];
                NodePointerType p_projection_node_2 =
                    projection_nodes[r_tetra_indices[2]];

                const array_1d<double, 3> edge_0 =
                    p_projection_node_0->Coordinates() -
                    p_surrogate_node->Coordinates();
                const array_1d<double, 3> edge_1 =
                    p_projection_node_1->Coordinates() -
                    p_surrogate_node->Coordinates();
                const array_1d<double, 3> edge_2 =
                    p_projection_node_2->Coordinates() -
                    p_surrogate_node->Coordinates();
                double determinant = inner_prod(
                    edge_0,
                    MathUtils<double>::CrossProduct(edge_1, edge_2));
                const double type3_volume = std::abs(determinant) / 6.0;
                const bool create_type3_quadrature =
                    type3_volume >= MinimumGapVolumeForQuadrature;

                if (determinant < 0.0) {
                    std::swap(p_projection_node_1, p_projection_node_2);
                    determinant = -determinant;
                }

                auto type3_geometry_data = CreateType3LinearCollapsedCornerVolume(
                    p_surrogate_node,
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2);

                p_projection_node_0 = type3_geometry_data.ProjectionNodes[0];
                p_projection_node_1 = type3_geometry_data.ProjectionNodes[1];
                p_projection_node_2 = type3_geometry_data.ProjectionNodes[2];

                const Type3VolumeKey volume_key = Type3VolumeKey::Create(
                    p_surrogate_node,
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2);
                KRATOS_ERROR_IF(
                    created_volume_keys.find(volume_key) !=
                    created_volume_keys.end())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Duplicated type3 closure tetra.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  projection nodes: "
                    << volume_key.ProjectionNodeIds[0] << ", "
                    << volume_key.ProjectionNodeIds[1] << ", "
                    << volume_key.ProjectionNodeIds[2] << "\n";

                auto p_type3_volume = type3_geometry_data.pVolume;

                const array_1d<double, 3> volume_center =
                    p_type3_volume->Center().Coordinates();

                Geometry<Node>::Pointer p_selected_neighbour_geometry = nullptr;
                bool selected_neighbour_has_type1_path = false;
                bool selected_neighbour_has_active_span = false;
                SpanKey3D selected_neighbour_active_span;
                std::string selected_neighbour_source;
                CanonicalFaceKey3D selected_neighbour_face_key;
                double selected_neighbour_distance =
                    std::numeric_limits<double>::max();
                std::size_t number_reused_faces_for_volume = 0;
                std::size_t number_reused_faces_with_type1_path = 0;

                auto select_neighbour_from_reused_face = [&](
                    const char* pFaceName,
                    const NodePointerType& pNode0,
                    const NodePointerType& pNode1,
                    const NodePointerType& pNode2)
                {
                    const auto face_key = MakeCanonicalFaceKey3D(
                        pNode0,
                        pNode1,
                        pNode2);

                    const auto face_registry_it =
                        mLateralFaceRegistry.find(face_key);
                    if (face_registry_it == mLateralFaceRegistry.end()) {
                        return;
                    }

                    const auto& r_existing_occurrences =
                        face_registry_it->second;
                    KRATOS_ERROR_IF(r_existing_occurrences.empty())
                        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                        << "Reused type3 lateral face has no occurrences.\n"
                        << "  face: " << pFaceName << "\n";
                    KRATOS_ERROR_IF(r_existing_occurrences.size() > 2)
                        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                        << "Non-manifold lateral face before type3 neighbour "
                        << "selection.\n"
                        << "  face: " << pFaceName << "\n"
                        << "  face node ids: "
                        << face_key.NodeIds[0] << ", "
                        << face_key.NodeIds[1] << ", "
                        << face_key.NodeIds[2] << "\n"
                        << "  occurrences: " << r_existing_occurrences.size()
                        << "\n"
                        << "  surrogate node: " << surrogate_node_id << "\n"
                        << "  projection nodes: "
                        << p_projection_node_0->Id() << ", "
                        << p_projection_node_1->Id() << ", "
                        << p_projection_node_2->Id() << "\n";

                    ++number_reused_faces_for_volume;

                    for (const auto& r_source_occurrence :
                         r_existing_occurrences) {
                        KRATOS_ERROR_IF_NOT(r_source_occurrence.pNeighbourGeometry)
                            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                            << "Reused type3 lateral face has null neighbour "
                            << "geometry.\n"
                            << "  face: " << pFaceName << "\n";

                        if (!r_source_occurrence.HasType1NeighbourPath) {
                            continue;
                        }

                        ++number_reused_faces_with_type1_path;

                        const double distance =
                            norm_2(
                                r_source_occurrence.pNeighbourGeometry->Center().Coordinates() -
                                volume_center);

                        if (!p_selected_neighbour_geometry ||
                            distance < selected_neighbour_distance) {
                            p_selected_neighbour_geometry =
                                r_source_occurrence.pNeighbourGeometry;
                            selected_neighbour_distance = distance;
                            selected_neighbour_has_type1_path =
                                r_source_occurrence.HasType1NeighbourPath;
                            selected_neighbour_has_active_span =
                                r_source_occurrence.HasNeighbourActiveSpan;
                            selected_neighbour_active_span =
                                r_source_occurrence.NeighbourActiveSpan;
                            selected_neighbour_source = pFaceName;
                            selected_neighbour_face_key = face_key;
                        }
                    }
                };

                select_neighbour_from_reused_face(
                    "S-P0-P1",
                    p_surrogate_node,
                    p_projection_node_0,
                    p_projection_node_1);
                select_neighbour_from_reused_face(
                    "S-P1-P2",
                    p_surrogate_node,
                    p_projection_node_1,
                    p_projection_node_2);
                select_neighbour_from_reused_face(
                    "S-P2-P0",
                    p_surrogate_node,
                    p_projection_node_2,
                    p_projection_node_0);
                select_neighbour_from_reused_face(
                    "P0-P1-P2",
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2);

                KRATOS_ERROR_IF(number_reused_faces_for_volume == 0)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Type3 volume has no reused lateral face. "
                    << "All lateral faces would be new.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n"
                    << "  projection nodes: "
                    << p_projection_node_0->Id() << ", "
                    << p_projection_node_1->Id() << ", "
                    << p_projection_node_2->Id() << "\n";

                KRATOS_ERROR_IF(number_reused_faces_with_type1_path == 0)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Type3 volume has reused lateral faces, but none has "
                    << "a continuous neighbour path to a type1 volume.\n"
                    << "  reused faces: " << number_reused_faces_for_volume
                    << "\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n"
                    << "  projection nodes: "
                    << p_projection_node_0->Id() << ", "
                    << p_projection_node_1->Id() << ", "
                    << p_projection_node_2->Id() << "\n";

                KRATOS_ERROR_IF_NOT(p_selected_neighbour_geometry)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Type3 volume did not select a neighbour geometry.\n"
                    << "  surrogate node: " << surrogate_node_id << "\n"
                    << "  direction: " << direction_names[direction_index] << "\n"
                    << "  projection nodes: "
                    << p_projection_node_0->Id() << ", "
                    << p_projection_node_1->Id() << ", "
                    << p_projection_node_2->Id() << "\n";

                KRATOS_ERROR_IF_NOT(selected_neighbour_has_active_span)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Selected neighbour geometry for type3 volume has no "
                    << "active span provenance.\n"
                    << "  selected source face: "
                    << selected_neighbour_source << "\n"
                    << "  selected source face node ids: "
                    << selected_neighbour_face_key.NodeIds[0] << ", "
                    << selected_neighbour_face_key.NodeIds[1] << ", "
                    << selected_neighbour_face_key.NodeIds[2] << "\n";

                created_volume_keys.insert(volume_key);
                p_type3_volume->SetId(next_geometry_id++);
                AddNeighbourGeometry(
                    *p_type3_volume,
                    p_selected_neighbour_geometry);
                if (mStoreGapDebugGeometries) {
                    r_gap_type3_debug.AddGeometry(p_type3_volume);
                }

                register_type3_face_if_needed(
                    p_surrogate_node,
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2,
                    p_selected_neighbour_geometry,
                    selected_neighbour_has_type1_path,
                    selected_neighbour_has_active_span,
                    selected_neighbour_active_span);
                register_type3_face_if_needed(
                    p_surrogate_node,
                    p_projection_node_1,
                    p_projection_node_2,
                    p_projection_node_0,
                    p_selected_neighbour_geometry,
                    selected_neighbour_has_type1_path,
                    selected_neighbour_has_active_span,
                    selected_neighbour_active_span);
                register_type3_face_if_needed(
                    p_surrogate_node,
                    p_projection_node_2,
                    p_projection_node_0,
                    p_projection_node_1,
                    p_selected_neighbour_geometry,
                    selected_neighbour_has_type1_path,
                    selected_neighbour_has_active_span,
                    selected_neighbour_active_span);
                register_type3_top_face_if_needed(
                    p_projection_node_0,
                    p_projection_node_1,
                    p_projection_node_2,
                    p_surrogate_node,
                    type3_geometry_data.pTopSurface,
                    p_selected_neighbour_geometry,
                    selected_neighbour_has_type1_path,
                    selected_neighbour_has_active_span,
                    selected_neighbour_active_span);

                if (!create_type3_quadrature) {
                    ++result.Summary.NumberOfSkippedFaces;
                    ++number_skipped_degenerate_candidate_volumes;
                    continue;
                }

                Type3VolumeQuadratureData type3_volume_data;
                type3_volume_data.NeighbourGeometries.push_back(
                    p_selected_neighbour_geometry);
                type3_volume_data.CharacteristicLength =
                    calculate_characteristic_length(
                        p_surrogate_node,
                        p_projection_node_0,
                        p_projection_node_1,
                        p_projection_node_2);
                type3_volume_data.SurrogateNodeId = surrogate_node_id;
                type3_volume_data.ProjectionNodeId0 =
                    p_projection_node_0->Id();
                type3_volume_data.ProjectionNodeId1 =
                    p_projection_node_1->Id();
                type3_volume_data.ProjectionNodeId2 =
                    p_projection_node_2->Id();
                type3_volume_data.pSurrogateNode = p_surrogate_node;
                type3_volume_data.pProjectionNode0 = p_projection_node_0;
                type3_volume_data.pProjectionNode1 = p_projection_node_1;
                type3_volume_data.pProjectionNode2 = p_projection_node_2;
                type3_volume_data.RequiresQuadrature = true;
                type3_volume_data.pVolumeGeometry = p_type3_volume;
                result.VolumeQuadratureDataList.push_back(
                    std::move(type3_volume_data));
                ++result.Summary.NumberOfCreatedVolumes;
            }

        }

    }
    std::size_t number_type3_shared_faces = 0;
    std::size_t number_type3_shared_faces_with_complete_neighbours = 0;

    for (auto& r_registry_entry : mLateralFaceRegistry) {
        auto& r_occurrences = r_registry_entry.second;
        if (r_occurrences.size() != 2) {
            continue;
        }

        bool has_type3_occurrence = false;
        for (const auto& r_occurrence : r_occurrences) {
            if (r_occurrence.GapType == 3) {
                has_type3_occurrence = true;
                break;
            }
        }

        if (!has_type3_occurrence) {
            continue;
        }

        ++number_type3_shared_faces;

        for (const auto& r_target_occurrence : r_occurrences) {
            KRATOS_ERROR_IF_NOT(r_target_occurrence.pNeighbourGeometry)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Shared type3 lateral face occurrence has null neighbour "
                << "geometry.\n"
                << "  face node ids: "
                << r_registry_entry.first.NodeIds[0] << ", "
                << r_registry_entry.first.NodeIds[1] << ", "
                << r_registry_entry.first.NodeIds[2] << "\n";

            for (const auto& r_surface_occurrence : r_occurrences) {
                KRATOS_ERROR_IF_NOT(r_surface_occurrence.pGeometry)
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Shared type3 lateral face occurrence has null surface "
                    << "geometry.\n"
                    << "  face node ids: "
                    << r_registry_entry.first.NodeIds[0] << ", "
                    << r_registry_entry.first.NodeIds[1] << ", "
                    << r_registry_entry.first.NodeIds[2] << "\n";

                AddUniqueNeighbourGeometry(
                    *r_surface_occurrence.pGeometry,
                    r_target_occurrence.pNeighbourGeometry);
            }
        }
    }

    for (const auto& r_registry_entry : mLateralFaceRegistry) {
        const auto& r_occurrences = r_registry_entry.second;
        if (r_occurrences.size() != 2) {
            continue;
        }

        bool has_type3_occurrence = false;
        for (const auto& r_occurrence : r_occurrences) {
            if (r_occurrence.GapType == 3) {
                has_type3_occurrence = true;
                break;
            }
        }

        if (!has_type3_occurrence) {
            continue;
        }

        for (const auto& r_occurrence : r_occurrences) {
            KRATOS_ERROR_IF(
                (r_occurrence.GapType == 2 ||
                 r_occurrence.GapType == 3) &&
                !r_occurrence.HasType1NeighbourPath)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Shared type2/type3 lateral face occurrence has no "
                << "continuous neighbour path to a type1 volume.\n"
                << "  face node ids: "
                << r_registry_entry.first.NodeIds[0] << ", "
                << r_registry_entry.first.NodeIds[1] << ", "
                << r_registry_entry.first.NodeIds[2] << "\n"
                << "  gap type: " << r_occurrence.GapType << "\n"
                << "  external span: "
                << SpanToString(r_occurrence.ExternalSpan) << "\n";

            KRATOS_ERROR_IF_NOT(r_occurrence.pGeometry)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Shared type3 lateral face occurrence has null surface "
                << "geometry.\n"
                << "  face node ids: "
                << r_registry_entry.first.NodeIds[0] << ", "
                << r_registry_entry.first.NodeIds[1] << ", "
                << r_registry_entry.first.NodeIds[2] << "\n";

            const auto stored_neighbour_geometries =
                r_occurrence.pGeometry->Has(NEIGHBOUR_GEOMETRIES) ?
                r_occurrence.pGeometry->GetValue(NEIGHBOUR_GEOMETRIES) :
                NeighbourGeometriesVectorType();

            for (const auto& r_other_occurrence : r_occurrences) {
                KRATOS_ERROR_IF_NOT(ContainsNeighbourGeometry(
                    stored_neighbour_geometries,
                    r_other_occurrence.pNeighbourGeometry))
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Shared type3 lateral face geometry does not contain "
                    << "both occurrence neighbour geometries.\n"
                    << "  face node ids: "
                    << r_registry_entry.first.NodeIds[0] << ", "
                    << r_registry_entry.first.NodeIds[1] << ", "
                    << r_registry_entry.first.NodeIds[2] << "\n"
                    << "  checked occurrence external span: "
                    << SpanToString(r_occurrence.ExternalSpan) << "\n"
                    << "  missing neighbour from external span: "
                    << SpanToString(r_other_occurrence.ExternalSpan) << "\n";
            }
        }

        ++number_type3_shared_faces_with_complete_neighbours;
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type 3 creation summary:\n"
        << "  candidate open type2 faces: "
        << result.Summary.NumberOfCandidateOpenFaces << "\n"
        << "  created volumes: " << result.Summary.NumberOfCreatedVolumes << "\n"
        << "  skipped volumes total: "
        << result.Summary.NumberOfSkippedFaces << "\n"
        << "  skipped degenerate volumes: "
        << number_skipped_degenerate_candidate_volumes << "\n"
        << "  skipped volumes with empty quadrature: "
        << number_skipped_empty_quadrature_volumes << "\n"
        << "  shared faces containing type3: "
        << number_type3_shared_faces << "\n"
        << "  shared type3 faces with complete neighbours: "
        << number_type3_shared_faces_with_complete_neighbours << "\n";

    return result;
}

#endif

#if 0
SnakeGapSbm3DUtilities::Type3CreationResult
SnakeGapSbm3DUtilities::CreateType3GapGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    const ExternalSpanDataMap& rExternalSpans,
    const KnotSpanGridInfo& rGridInfo,
    const Type2CreationResult& rType2CreationResult,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives)
{
    Type3CreationResult result;
    result.Summary.NumberOfCandidateOpenFaces =
        rType2CreationResult.OpenFaceDataList.size();

    ModelPart& r_gap_type3_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType3Debug");

    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    struct Type3OpenFaceCandidate {
        NodePointerType pSurrogateNode;
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;
        Geometry<Node>::Pointer pNeighbourGeometry;
        IndexType SurrogateConditionId = 0;
        SpanKey3D ExternalSpan;
        SpanKey3D ExternalSpan0;
        SpanKey3D ExternalSpan1;
    };

    std::vector<Type3OpenFaceCandidate> open_face_candidates;
    open_face_candidates.reserve(rType2CreationResult.OpenFaceDataList.size());

    std::map<IndexType, std::vector<NodePointerType>>
        incident_projection_nodes_by_surrogate_node_id;

    for (const auto& r_open_face_data : rType2CreationResult.OpenFaceDataList) {
        KRATOS_ERROR_IF_NOT(r_open_face_data.pSurrogateNode)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face has null surrogate node.\n";
        KRATOS_ERROR_IF_NOT(r_open_face_data.pProjectionNode0)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face has null first projection node.\n";
        KRATOS_ERROR_IF_NOT(r_open_face_data.pProjectionNode1)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face has null second projection node.\n";
        KRATOS_ERROR_IF_NOT(AreFaceAdjacentSpans(
            r_open_face_data.ExternalSpan0,
            r_open_face_data.ExternalSpan1))
            << "Trying to create type3 from non-adjacent type2 open face spans";

        const auto key = MakeCanonicalFaceKey3D(
            r_open_face_data.pSurrogateNode->Id(),
            r_open_face_data.pProjectionNode0->Id(),
            r_open_face_data.pProjectionNode1->Id());

        const auto registry_it = mLateralFaceRegistry.find(key);
        if (registry_it == mLateralFaceRegistry.end()) {
            ++result.Summary.NumberOfSkippedFaces;
            continue;
        }

        const auto& r_occurrences = registry_it->second;
        if (r_occurrences.size() > 2) {
            KRATOS_ERROR << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                         << "Non-manifold lateral face for node ids "
                         << key.NodeIds[0] << ", "
                         << key.NodeIds[1] << ", "
                         << key.NodeIds[2] << ". Occurrences: "
                         << r_occurrences.size() << ".\n";
        }

        if (r_occurrences.size() != 1 || r_occurrences.front().GapType != 2) {
            ++result.Summary.NumberOfSkippedFaces;
            continue;
        }

        const auto& r_reference_occurrence = r_occurrences.front();
        KRATOS_ERROR_IF_NOT(r_reference_occurrence.pGeometry)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face occurrence has null geometry.\n";
        KRATOS_ERROR_IF_NOT(r_reference_occurrence.pNeighbourGeometry)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face occurrence has null neighbour geometry.\n";

        Type3OpenFaceCandidate candidate;
        candidate.pSurrogateNode = r_open_face_data.pSurrogateNode;
        candidate.pProjectionNode0 = r_open_face_data.pProjectionNode0;
        candidate.pProjectionNode1 = r_open_face_data.pProjectionNode1;
        candidate.pNeighbourGeometry = r_reference_occurrence.pNeighbourGeometry;
        candidate.SurrogateConditionId = r_reference_occurrence.SurrogateConditionId;
        candidate.ExternalSpan = r_reference_occurrence.ExternalSpan;
        candidate.ExternalSpan0 = r_open_face_data.ExternalSpan0;
        candidate.ExternalSpan1 = r_open_face_data.ExternalSpan1;

        open_face_candidates.push_back(std::move(candidate));

        const IndexType surrogate_node_id = r_open_face_data.pSurrogateNode->Id();
        auto& r_incident_nodes =
            incident_projection_nodes_by_surrogate_node_id[surrogate_node_id];

        auto add_incident_node = [&r_incident_nodes](
            const NodePointerType& pProjectionNode)
        {
            const auto node_it = std::find_if(
                r_incident_nodes.begin(),
                r_incident_nodes.end(),
                [&pProjectionNode](const NodePointerType& pExistingNode) {
                    return pExistingNode->Id() == pProjectionNode->Id();
                });

            if (node_it == r_incident_nodes.end()) {
                r_incident_nodes.push_back(pProjectionNode);
            }
        };

        add_incident_node(r_open_face_data.pProjectionNode0);
        add_incident_node(r_open_face_data.pProjectionNode1);
    }

#if 0
    const double point_tolerance =
        1.0e-10 * std::max({
            rGridInfo.SpanSizeX,
            rGridInfo.SpanSizeY,
            rGridInfo.SpanSizeZ});

    auto find_existing_skin_node_at_point = [&](
        const array_1d<double, 3>& rPoint) -> NodePointerType
    {
        for (const auto& r_node : rSkinSubModelPart.Nodes()) {
            if (norm_2(r_node.Coordinates() - rPoint) <= point_tolerance) {
                return rSkinSubModelPart.pGetNode(r_node.Id());
            }
        }

        return NodePointerType();
    };

    auto create_or_reuse_skin_node = [&](
        const array_1d<double, 3>& rPoint,
        bool& rCreatedNode) -> NodePointerType
    {
        if (auto p_existing_node = find_existing_skin_node_at_point(rPoint)) {
            rCreatedNode = false;
            return p_existing_node;
        }

        const IndexType new_node_id = GetNextAuxiliarySkinNodeId(rSkinSubModelPart);
        auto p_node = rSkinSubModelPart.CreateNewNode(
            new_node_id,
            rPoint[0],
            rPoint[1],
            rPoint[2]);

        rCreatedNode = true;
        return p_node;
    };

    auto condition_overlaps_span = [&](
        const Geometry<Node>& rGeometry,
        const array_1d<double, 3>& rSpanMin,
        const array_1d<double, 3>& rSpanMax) -> bool
    {
        array_1d<double, 3> condition_min;
        array_1d<double, 3> condition_max;
        ConditionBoundingBox(rGeometry, condition_min, condition_max);

        for (IndexType axis = 0; axis < 3; ++axis) {
            if (condition_max[axis] < rSpanMin[axis] - point_tolerance ||
                condition_min[axis] > rSpanMax[axis] + point_tolerance) {
                return false;
            }
        }

        return true;
    };

    auto squared_distance_3d = [](
        const array_1d<double, 3>& rPoint0,
        const array_1d<double, 3>& rPoint1)
    {
        return inner_prod(rPoint0 - rPoint1, rPoint0 - rPoint1);
    };

    auto closest_point_on_triangle_3d = [](
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC) -> array_1d<double, 3>
    {
        const array_1d<double, 3> ab = rB - rA;
        const array_1d<double, 3> ac = rC - rA;
        const array_1d<double, 3> ap = rPoint - rA;
        const double d1 = inner_prod(ab, ap);
        const double d2 = inner_prod(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0) {
            return rA;
        }

        const array_1d<double, 3> bp = rPoint - rB;
        const double d3 = inner_prod(ab, bp);
        const double d4 = inner_prod(ac, bp);
        if (d3 >= 0.0 && d4 <= d3) {
            return rB;
        }

        const double vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
            const double v = d1 / (d1 - d3);
            array_1d<double, 3> closest_point = rA;
            noalias(closest_point) += v * ab;
            return closest_point;
        }

        const array_1d<double, 3> cp = rPoint - rC;
        const double d5 = inner_prod(ab, cp);
        const double d6 = inner_prod(ac, cp);
        if (d6 >= 0.0 && d5 <= d6) {
            return rC;
        }

        const double vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
            const double w = d2 / (d2 - d6);
            array_1d<double, 3> closest_point = rA;
            noalias(closest_point) += w * ac;
            return closest_point;
        }

        const double va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
            const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            array_1d<double, 3> closest_point = rB;
            noalias(closest_point) += w * (rC - rB);
            return closest_point;
        }

        const double denominator = 1.0 / (va + vb + vc);
        const double v = vb * denominator;
        const double w = vc * denominator;
        array_1d<double, 3> closest_point = rA;
        noalias(closest_point) += v * ab;
        noalias(closest_point) += w * ac;
        return closest_point;
    };

    auto find_closest_skin_point_in_span = [&](
        const array_1d<double, 3>& rReferencePoint,
        const SpanKey3D& rSpan,
        array_1d<double, 3>& rClosestPoint) -> bool
    {
        array_1d<double, 3> span_min = ZeroVector(3);
        array_1d<double, 3> span_max = ZeroVector(3);
        SpanBox(rSpan, rGridInfo, span_min, span_max);

        bool found_candidate = false;
        double best_distance_squared = std::numeric_limits<double>::max();

        auto consider_point = [&](
            const array_1d<double, 3>& rCandidatePoint)
        {
            if (!IsPointInsideBox(
                    rCandidatePoint,
                    span_min,
                    span_max,
                    point_tolerance)) {
                return;
            }

            const double distance_squared =
                squared_distance_3d(rCandidatePoint, rReferencePoint);
            if (distance_squared < best_distance_squared) {
                best_distance_squared = distance_squared;
                rClosestPoint = rCandidatePoint;
                found_candidate = true;
            }
        };

        auto consider_triangle = [&](
            const array_1d<double, 3>& rA,
            const array_1d<double, 3>& rB,
            const array_1d<double, 3>& rC)
        {
            consider_point(closest_point_on_triangle_3d(
                rReferencePoint,
                rA,
                rB,
                rC));
        };

        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            const auto& r_geometry = r_condition.GetGeometry();
            if (r_geometry.PointsNumber() < 3 ||
                !condition_overlaps_span(r_geometry, span_min, span_max)) {
                continue;
            }

            consider_triangle(
                r_geometry[0].Coordinates(),
                r_geometry[1].Coordinates(),
                r_geometry[2].Coordinates());

            if (r_geometry.PointsNumber() == 4) {
                consider_triangle(
                    r_geometry[0].Coordinates(),
                    r_geometry[2].Coordinates(),
                    r_geometry[3].Coordinates());
            }
        }

        if (found_candidate) {
            return true;
        }

        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            array_1d<double, 3> intersection_point = ZeroVector(3);
            if (FindConditionSpanIntersectionPoint(
                    r_condition.GetGeometry(),
                    rSpan,
                    rGridInfo,
                    intersection_point)) {
                consider_point(intersection_point);
            }
        }

        for (const auto& r_node : rSkinSubModelPart.Nodes()) {
            consider_point(r_node.Coordinates());
        }

        return found_candidate;
    };

    auto ray_triangle_intersection_3d = [](
        const array_1d<double, 3>& rOrigin,
        const array_1d<double, 3>& rDirection,
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC,
        double& rDistance,
        array_1d<double, 3>& rIntersectionPoint)
    {
        const double tolerance = 1.0e-12;
        const array_1d<double, 3> edge_1 = rB - rA;
        const array_1d<double, 3> edge_2 = rC - rA;
        const array_1d<double, 3> p_vector =
            MathUtils<double>::CrossProduct(rDirection, edge_2);
        const double determinant = inner_prod(edge_1, p_vector);
        if (std::abs(determinant) <= tolerance) {
            return false;
        }

        const double inverse_determinant = 1.0 / determinant;
        const array_1d<double, 3> t_vector = rOrigin - rA;
        const double u = inner_prod(t_vector, p_vector) * inverse_determinant;
        if (u < -tolerance || u > 1.0 + tolerance) {
            return false;
        }

        const array_1d<double, 3> q_vector =
            MathUtils<double>::CrossProduct(t_vector, edge_1);
        const double v = inner_prod(rDirection, q_vector) * inverse_determinant;
        if (v < -tolerance || u + v > 1.0 + tolerance) {
            return false;
        }

        const double distance = inner_prod(edge_2, q_vector) * inverse_determinant;
        if (distance <= tolerance) {
            return false;
        }

        rDistance = distance;
        rIntersectionPoint = rOrigin + distance * rDirection;
        return true;
    };

    auto find_ray_skin_intersection = [&](
        const array_1d<double, 3>& rOrigin,
        const array_1d<double, 3>& rDirection,
        array_1d<double, 3>& rIntersectionPoint) -> bool
    {
        bool found_intersection = false;
        double best_distance = std::numeric_limits<double>::max();

        auto consider_triangle = [&](
            const array_1d<double, 3>& rA,
            const array_1d<double, 3>& rB,
            const array_1d<double, 3>& rC)
        {
            double distance = 0.0;
            array_1d<double, 3> intersection_point = ZeroVector(3);

            if (ray_triangle_intersection_3d(
                    rOrigin,
                    rDirection,
                    rA,
                    rB,
                    rC,
                    distance,
                    intersection_point) &&
                distance < best_distance) {
                best_distance = distance;
                rIntersectionPoint = intersection_point;
                found_intersection = true;
            }
        };

        for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
            const auto& r_geometry = r_condition.GetGeometry();
            if (r_geometry.PointsNumber() < 3) {
                continue;
            }

            consider_triangle(
                r_geometry[0].Coordinates(),
                r_geometry[1].Coordinates(),
                r_geometry[2].Coordinates());

            if (r_geometry.PointsNumber() == 4) {
                consider_triangle(
                    r_geometry[0].Coordinates(),
                    r_geometry[2].Coordinates(),
                    r_geometry[3].Coordinates());
            }
        }

        return found_intersection;
    };

    struct CornerProjectionData {
        NodePointerType pNode;
        SpanKey3D Type3Span;
        bool HasType3Span = false;
    };

    std::unordered_map<IndexType, CornerProjectionData> corner_projection_by_surrogate_node_id;

    auto get_corner_projection_for_surrogate_node = [&](
        const NodePointerType& pSurrogateNode) -> CornerProjectionData
    {
        const IndexType surrogate_node_id = pSurrogateNode->Id();

        const auto cache_it =
            corner_projection_by_surrogate_node_id.find(surrogate_node_id);
        if (cache_it != corner_projection_by_surrogate_node_id.end()) {
            ++result.Summary.NumberOfReusedCornerProjectionNodes;
            return cache_it->second;
        }

        const GridPointKey3D grid_node =
            ComputeGridPointKey(pSurrogateNode->Coordinates(), rGridInfo);

        std::vector<SpanKey3D> type3_spans;
        for (int di = -1; di <= 0; ++di) {
            for (int dj = -1; dj <= 0; ++dj) {
                for (int dk = -1; dk <= 0; ++dk) {
                    const SpanKey3D span{
                        grid_node.I + di,
                        grid_node.J + dj,
                        grid_node.K + dk};

                    if (!IsValidSpan(span, rGridInfo)) {
                        continue;
                    }

                    const auto external_span_it = rExternalSpans.find(span);
                    if (external_span_it != rExternalSpans.end() &&
                        external_span_it->second.Type == GapSpanType::Type3) {
                        type3_spans.push_back(span);
                    }
                }
            }
        }

        KRATOS_ERROR_IF(type3_spans.size() > 1)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Found more than one type 3 span incident to surrogate node #"
            << surrogate_node_id << ".\n";

        bool created_node = false;
        CornerProjectionData data;

        if (type3_spans.size() == 1) {
            array_1d<double, 3> closest_point = ZeroVector(3);
            KRATOS_ERROR_IF_NOT(find_closest_skin_point_in_span(
                pSurrogateNode->Coordinates(),
                type3_spans.front(),
                closest_point))
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Failed to find a skin projection point for surrogate node #"
                << surrogate_node_id
                << " in type 3 span " << SpanToString(type3_spans.front()) << ".\n";

            data.pNode = create_or_reuse_skin_node(
                closest_point,
                created_node);
            data.Type3Span = type3_spans.front();
            data.HasType3Span = true;
        } else {
            auto normal = pSurrogateNode->GetValue(NORMAL);
            const double normal_norm = norm_2(normal);
            KRATOS_ERROR_IF(normal_norm <= std::numeric_limits<double>::epsilon())
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Normal is zero for surrogate node #"
                << surrogate_node_id << ".\n";

            array_1d<double, 3> direction = normal / normal_norm;

            array_1d<double, 3> intersection_point = ZeroVector(3);
            bool intersection_found = find_ray_skin_intersection(
                pSurrogateNode->Coordinates(),
                direction,
                intersection_point);

            if (!intersection_found) {
                direction *= -1.0;
                KRATOS_ERROR_IF_NOT(find_ray_skin_intersection(
                    pSurrogateNode->Coordinates(),
                    direction,
                    intersection_point))
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Failed to intersect the true boundary from surrogate node #"
                    << surrogate_node_id
                    << " along both mean projection directions.\n";
            }

            data.pNode = create_or_reuse_skin_node(
                intersection_point,
                created_node);
            data.HasType3Span = false;
        }


        if (created_node) {
            ++result.Summary.NumberOfCreatedCornerProjectionNodes;
        } else {
            ++result.Summary.NumberOfReusedCornerProjectionNodes;
        }

        corner_projection_by_surrogate_node_id.emplace(
            surrogate_node_id,
            data);

        return data;
    };

#endif

#if 0
    std::map<IndexType, NodePointerType>
        corner_projection_nodes_by_surrogate_node_id;

    auto get_corner_projection_for_surrogate_node = [&](
        const NodePointerType& pSurrogateNode) -> NodePointerType
    {
        const IndexType surrogate_node_id = pSurrogateNode->Id();
        const auto cache_it =
            corner_projection_nodes_by_surrogate_node_id.find(surrogate_node_id);

        if (cache_it != corner_projection_nodes_by_surrogate_node_id.end()) {
            ++result.Summary.NumberOfReusedCornerProjectionNodes;
            return cache_it->second;
        }

        const auto incident_nodes_it =
            incident_projection_nodes_by_surrogate_node_id.find(surrogate_node_id);
        KRATOS_ERROR_IF(
            incident_nodes_it ==
            incident_projection_nodes_by_surrogate_node_id.end())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "No incident type 2 projection nodes found for surrogate node #"
            << surrogate_node_id << ".\n";

        const GridPointKey3D surrogate_grid_node =
            ComputeGridPointKey(pSurrogateNode->Coordinates(), rGridInfo);

        auto p_corner_projection_node =
            GetOrCreateType3CornerProjectionNode(
                rRootModelPart,
                rSkinSubModelPart,
                pSurrogateNode,
                surrogate_grid_node,
                incident_nodes_it->second,
                rExternalSpans,
                rGridInfo,
                result.Summary);

        corner_projection_nodes_by_surrogate_node_id.emplace(
            surrogate_node_id,
            p_corner_projection_node);

        return p_corner_projection_node;
    };
#endif

    auto calculate_type3_characteristic_length = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const NodePointerType& pNode3)
    {
        const std::array<NodePointerType, 4> nodes = {{
            pNode0,
            pNode1,
            pNode2,
            pNode3
        }};

        double characteristic_length = 0.0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            for (std::size_t j = i + 1; j < nodes.size(); ++j) {
                characteristic_length = std::max(
                    characteristic_length,
                    norm_2(nodes[i]->Coordinates() - nodes[j]->Coordinates()));
            }
        }

        return characteristic_length;
    };

    for (const auto& r_candidate : open_face_candidates) {
        const auto open_face_key = MakeCanonicalFaceKey3D(
            r_candidate.pSurrogateNode->Id(),
            r_candidate.pProjectionNode0->Id(),
            r_candidate.pProjectionNode1->Id());

        const auto registry_it = mLateralFaceRegistry.find(open_face_key);
        if (registry_it == mLateralFaceRegistry.end() ||
            registry_it->second.size() != 1 ||
            registry_it->second.front().GapType != 2) {
            ++result.Summary.NumberOfSkippedFaces;
            continue;
        }

        const auto corner_projection_data =
            get_corner_projection_for_surrogate_node(r_candidate.pSurrogateNode);
        const auto p_corner_projection_node = corner_projection_data.pNode;

        const array_1d<double, 3> edge_0 =
            r_candidate.pProjectionNode0->Coordinates() -
            r_candidate.pSurrogateNode->Coordinates();
        const array_1d<double, 3> edge_1 =
            r_candidate.pProjectionNode1->Coordinates() -
            r_candidate.pSurrogateNode->Coordinates();
        const array_1d<double, 3> corner_edge =
            p_corner_projection_node->Coordinates() -
            r_candidate.pSurrogateNode->Coordinates();
        const array_1d<double, 3> cross_product =
            MathUtils<double>::CrossProduct(edge_1, corner_edge);
        const double tetrahedron_volume =
            std::abs(inner_prod(edge_0, cross_product)) / 6.0;
        const double minimum_span_size = std::min({
            rGridInfo.SpanSizeX,
            rGridInfo.SpanSizeY,
            rGridInfo.SpanSizeZ});
        const double volume_tolerance =
            1.0e-12 * minimum_span_size * minimum_span_size * minimum_span_size;

        if (tetrahedron_volume <= volume_tolerance) {
            ++result.Summary.NumberOfSkippedFaces;
            continue;
        }

        auto p_type3_volume = CreateType3CollapsedCornerVolume(
            r_candidate.pSurrogateNode,
            r_candidate.pProjectionNode0,
            r_candidate.pProjectionNode1,
            p_corner_projection_node);

        GeometriesArrayType volume_quadrature_point_geometries =
            CreateAndTagVolumeQuadraturePointGeometries(
                p_type3_volume,
                r_candidate.pNeighbourGeometry,
                IntegrationOrder,
                NumberOfShapeFunctionsDerivatives);

        if (volume_quadrature_point_geometries.size() == 0) {
            ++result.Summary.NumberOfSkippedFaces;
            continue;
        }

        CheckType3QuadraturePointGeometries(
            volume_quadrature_point_geometries,
            r_candidate.pSurrogateNode,
            r_candidate.pProjectionNode0,
            r_candidate.pProjectionNode1,
            p_corner_projection_node);

        p_type3_volume->SetId(next_geometry_id++);
        AddNeighbourGeometry(*p_type3_volume, r_candidate.pNeighbourGeometry);
        if (mStoreGapDebugGeometries) {
            r_gap_type3_debug.AddGeometry(p_type3_volume);
        }

        RegisterType3LateralFace(
            r_gap_type3_debug,
            next_geometry_id,
            r_candidate.SurrogateConditionId,
            r_candidate.ExternalSpan,
            r_candidate.pSurrogateNode,
            r_candidate.pProjectionNode0,
            r_candidate.pProjectionNode1,
            p_corner_projection_node,
            r_candidate.pNeighbourGeometry);

        RegisterType3LateralFace(
            r_gap_type3_debug,
            next_geometry_id,
            r_candidate.SurrogateConditionId,
            r_candidate.ExternalSpan,
            r_candidate.pSurrogateNode,
            r_candidate.pProjectionNode0,
            p_corner_projection_node,
            r_candidate.pProjectionNode1,
            r_candidate.pNeighbourGeometry);

        RegisterType3LateralFace(
            r_gap_type3_debug,
            next_geometry_id,
            r_candidate.SurrogateConditionId,
            r_candidate.ExternalSpan,
            r_candidate.pSurrogateNode,
            p_corner_projection_node,
            r_candidate.pProjectionNode1,
            r_candidate.pProjectionNode0,
            r_candidate.pNeighbourGeometry);

        RegisterType3LateralFace(
            r_gap_type3_debug,
            next_geometry_id,
            r_candidate.SurrogateConditionId,
            r_candidate.ExternalSpan,
            r_candidate.pProjectionNode0,
            r_candidate.pProjectionNode1,
            p_corner_projection_node,
            r_candidate.pSurrogateNode,
            r_candidate.pNeighbourGeometry);

        Type3VolumeQuadratureData type3_volume_data;
        type3_volume_data.VolumeQuadraturePointGeometries =
            std::move(volume_quadrature_point_geometries);
        type3_volume_data.NeighbourGeometries.push_back(
            r_candidate.pNeighbourGeometry);
        type3_volume_data.CharacteristicLength =
            calculate_type3_characteristic_length(
                r_candidate.pSurrogateNode,
                r_candidate.pProjectionNode0,
                r_candidate.pProjectionNode1,
                p_corner_projection_node);
        type3_volume_data.SurrogateNodeId =
            r_candidate.pSurrogateNode->Id();
        type3_volume_data.ProjectionNodeId0 =
            r_candidate.pProjectionNode0->Id();
        type3_volume_data.ProjectionNodeId1 =
            r_candidate.pProjectionNode1->Id();
        type3_volume_data.CornerProjectionNodeId =
            p_corner_projection_node->Id();

        result.VolumeQuadratureDataList.emplace_back(
            std::move(type3_volume_data));

        ++result.Summary.NumberOfCreatedVolumes;
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type 3 creation summary:\n"
        << "  candidate open faces: " << result.Summary.NumberOfCandidateOpenFaces << "\n"
        << "  created volumes: " << result.Summary.NumberOfCreatedVolumes << "\n"
        << "  skipped faces: " << result.Summary.NumberOfSkippedFaces << "\n"
        << "  created corner projection nodes: "
        << result.Summary.NumberOfCreatedCornerProjectionNodes << "\n"
        << "  reused corner projection nodes: "
        << result.Summary.NumberOfReusedCornerProjectionNodes << "\n"
        << "  GapType3Debug geometries: " << r_gap_type3_debug.NumberOfGeometries() << "\n";

    return result;
}

#endif

void SnakeGapSbm3DUtilities::UpdateLateralFaceRegistryGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart)
{
    ModelPart& r_gap_type2_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType2Debug");
    ModelPart& r_gap_type3_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType3Debug");
    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    auto is_skin_node = [&](
        const NodePointerType& pNode)
    {
        return pNode && rSkinSubModelPart.HasNode(pNode->Id());
    };

    auto make_edge_key = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1)
    {
        return SkinEdgeKey(
            std::min(pNode0->Id(), pNode1->Id()),
            std::max(pNode0->Id(), pNode1->Id()));
    };

    auto make_top_key = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2)
    {
        return SkinTopFaceKey(
            pNode0->Id(),
            pNode1->Id(),
            pNode2->Id());
    };

    for (auto& r_top_entry : mCurvedTopFaceRegistry) {
        auto& r_top_data = r_top_entry.second;
        if (!r_top_data.IsCurved) {
            continue;
        }

        if (r_top_data.pSurfaceGeometry) {
            continue;
        }

        PointerVector<NodeType> control_points;
        for (const auto& p_control_node : r_top_data.CurrentControlNodes) {
            KRATOS_ERROR_IF_NOT(p_control_node)
                << "[UpdateLateralFaceRegistryGeometries] Null top-face "
                << "control node.\n";
            control_points.push_back(p_control_node);
        }

        r_top_data.pSurfaceGeometry = Kratos::make_shared<NurbsSurfaceType>(
            control_points,
            std::size_t(2),
            std::size_t(2),
            CreateOpenUnitKnotVector(2),
            CreateOpenUnitKnotVector(2));
        r_top_data.pSurfaceGeometry->SetId(next_geometry_id++);
        if (mStoreGapDebugGeometries) {
            r_gap_type3_debug.AddGeometry(r_top_data.pSurfaceGeometry);
        }
    }

    for (auto& r_registry_entry : mLateralFaceRegistry) {
        auto& r_occurrences = r_registry_entry.second;
        if (r_occurrences.empty()) {
            continue;
        }

        const auto& r_reference_occurrence = r_occurrences.front();
        std::array<NodePointerType, 3> face_nodes =
            r_reference_occurrence.FaceNodes;

        if (!face_nodes[0] || !face_nodes[1] || !face_nodes[2]) {
            continue;
        }

        NurbsSurfaceType::Pointer p_final_surface;
        const bool node_0_is_skin = is_skin_node(face_nodes[0]);
        const bool node_1_is_skin = is_skin_node(face_nodes[1]);
        const bool node_2_is_skin = is_skin_node(face_nodes[2]);
        const std::size_t number_of_skin_nodes =
            static_cast<std::size_t>(node_0_is_skin) +
            static_cast<std::size_t>(node_1_is_skin) +
            static_cast<std::size_t>(node_2_is_skin);

        if (number_of_skin_nodes == 3) {
            const auto top_key = make_top_key(
                face_nodes[0],
                face_nodes[1],
                face_nodes[2]);
            const auto top_it = mCurvedTopFaceRegistry.find(top_key);
            if (top_it != mCurvedTopFaceRegistry.end() &&
                top_it->second.IsCurved) {
                p_final_surface = top_it->second.pSurfaceGeometry;
            }
        }

        if (!p_final_surface && number_of_skin_nodes == 2) {
            NodePointerType p_skin_node_0;
            NodePointerType p_skin_node_1;
            if (node_0_is_skin && node_1_is_skin) {
                p_skin_node_0 = face_nodes[0];
                p_skin_node_1 = face_nodes[1];
            } else if (node_1_is_skin && node_2_is_skin) {
                p_skin_node_0 = face_nodes[1];
                p_skin_node_1 = face_nodes[2];
            } else {
                p_skin_node_0 = face_nodes[2];
                p_skin_node_1 = face_nodes[0];
            }

            const auto edge_it = mCurvedEdgeRegistry.find(
                make_edge_key(p_skin_node_0, p_skin_node_1));
            if (edge_it != mCurvedEdgeRegistry.end() &&
                edge_it->second.IsCurved) {
                p_final_surface = CreateCollapsedTriangleSurface(
                    face_nodes[0],
                    face_nodes[1],
                    face_nodes[2]);
                p_final_surface->SetId(next_geometry_id++);
                if (mStoreGapDebugGeometries) {
                    bool has_type3_occurrence = false;
                    for (const auto& r_occurrence : r_occurrences) {
                        if (r_occurrence.GapType == 3) {
                            has_type3_occurrence = true;
                            break;
                        }
                    }
                    if (has_type3_occurrence) {
                        r_gap_type3_debug.AddGeometry(p_final_surface);
                    } else {
                        r_gap_type2_debug.AddGeometry(p_final_surface);
                    }
                }
            }
        }

        if (!p_final_surface) {
            continue;
        }

        for (const auto& r_occurrence : r_occurrences) {
            AddUniqueNeighbourGeometry(
                *p_final_surface,
                r_occurrence.pNeighbourGeometry);
        }

        for (auto& r_occurrence : r_occurrences) {
            r_occurrence.pGeometry = p_final_surface;
        }
    }
}

void SnakeGapSbm3DUtilities::FinalizeType2AndType3GapGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    Type2CreationResult& rType2CreationResult,
    Type3CreationResult& rType3CreationResult,
    const KnotSpanGridInfo& rGridInfo,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives)
{
    mCurvedEdgeRegistry.clear();
    mCurvedTopFaceRegistry.clear();
    mSkinEdgeControlNodes.clear();
    mLinearSkinEdges.clear();

    const std::size_t estimated_type2_volumes =
        rType2CreationResult.VolumeQuadratureDataList.size();
    const std::size_t estimated_type3_volumes =
        rType3CreationResult.VolumeQuadratureDataList.size();
    const std::size_t estimated_skin_edges =
        estimated_type2_volumes + 3 * estimated_type3_volumes;
    mCurvedEdgeRegistry.reserve(estimated_skin_edges);
    mSkinEdgeControlNodes.reserve(estimated_skin_edges);
    mCurvedTopFaceRegistry.reserve(estimated_type3_volumes);

    mSkinProjectionTriangleData.clear();
    mSkinProjectionTriangleData.reserve(2 * rSkinSubModelPart.NumberOfConditions());

    auto add_skin_projection_triangle = [&](
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB,
        const array_1d<double, 3>& rC)
    {
        SkinProjectionTriangleData triangle_data;
        triangle_data.A = rA;
        triangle_data.B = rB;
        triangle_data.C = rC;
        triangle_data.Center = (rA + rB + rC) / 3.0;
        triangle_data.Radius = std::max({
            norm_2(rA - triangle_data.Center),
            norm_2(rB - triangle_data.Center),
            norm_2(rC - triangle_data.Center)});
        mSkinProjectionTriangleData.push_back(std::move(triangle_data));
    };

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.PointsNumber() < 3) {
            continue;
        }

        add_skin_projection_triangle(
            r_geometry[0].Coordinates(),
            r_geometry[1].Coordinates(),
            r_geometry[2].Coordinates());

        if (r_geometry.PointsNumber() == 4) {
            add_skin_projection_triangle(
                r_geometry[0].Coordinates(),
                r_geometry[2].Coordinates(),
                r_geometry[3].Coordinates());
        }
    }

    const double maximum_projection_distance =
        0.75 * std::min({
            rGridInfo.SpanSizeX,
            rGridInfo.SpanSizeY,
            rGridInfo.SpanSizeZ});

    auto make_edge_key = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1)
    {
        KRATOS_ERROR_IF_NOT(pNode0)
            << "[FinalizeType2AndType3GapGeometries] First edge node is null.\n";
        KRATOS_ERROR_IF_NOT(pNode1)
            << "[FinalizeType2AndType3GapGeometries] Second edge node is null.\n";

        return SkinEdgeKey(
            std::min(pNode0->Id(), pNode1->Id()),
            std::max(pNode0->Id(), pNode1->Id()));
    };

    auto make_top_key = [](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2)
    {
        KRATOS_ERROR_IF_NOT(pNode0)
            << "[FinalizeType2AndType3GapGeometries] First top node is null.\n";
        KRATOS_ERROR_IF_NOT(pNode1)
            << "[FinalizeType2AndType3GapGeometries] Second top node is null.\n";
        KRATOS_ERROR_IF_NOT(pNode2)
            << "[FinalizeType2AndType3GapGeometries] Third top node is null.\n";

        return SkinTopFaceKey(
            pNode0->Id(),
            pNode1->Id(),
            pNode2->Id());
    };

    std::unordered_map<
        SkinEdgeKey,
        std::vector<SkinTopFaceKey>,
        SkinEdgeKeyHash> top_faces_by_edge;
    top_faces_by_edge.reserve(3 * estimated_type3_volumes);

    auto get_oriented_edge_controls = [&](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1)
    {
        auto controls = GetOrCreateFinalSkinEdgeControlNodes(
            rSkinSubModelPart,
            pNode0,
            pNode1,
            maximum_projection_distance,
            "finalization");
        KRATOS_ERROR_IF(controls.size() != mGapApproximationOrder + 1)
            << "[FinalizeType2AndType3GapGeometries] Unexpected edge control "
            << "node count: " << controls.size() << ". Expected "
            << mGapApproximationOrder + 1 << ".\n";
        return controls;
    };

    for (auto& r_type2_data : rType2CreationResult.VolumeQuadratureDataList) {
        if (!r_type2_data.RequiresQuadrature) {
            continue;
        }

        get_oriented_edge_controls(
            r_type2_data.pProjectionNode0,
            r_type2_data.pProjectionNode1);
    }

    auto is_edge_curved = [&](
        const SkinEdgeKey& rKey)
    {
        const auto edge_it = mCurvedEdgeRegistry.find(rKey);
        return edge_it != mCurvedEdgeRegistry.end() &&
               edge_it->second.IsCurved;
    };

    auto create_quadratic_top_controls = [&](
        const NodePointerType& pNode0,
        const NodePointerType& pNode1,
        const NodePointerType& pNode2,
        const double TopAlpha)
    {
        KRATOS_ERROR_IF(mGapApproximationOrder != 2)
            << "[FinalizeType2AndType3GapGeometries] Top-face curving "
            << "currently requires quadratic gap geometry. Requested order: "
            << mGapApproximationOrder << ".\n";

        const auto edge_01 = get_oriented_edge_controls(pNode0, pNode1);
        const auto edge_12 = get_oriented_edge_controls(pNode1, pNode2);
        const auto edge_02 = get_oriented_edge_controls(pNode0, pNode2);

        array_1d<double, 3> coons_middle = ZeroVector(3);
        noalias(coons_middle) += 0.25 * edge_01[1]->Coordinates();
        noalias(coons_middle) += 0.25 * edge_02[1]->Coordinates();
        noalias(coons_middle) += 0.25 * edge_12[1]->Coordinates();
        noalias(coons_middle) += 0.25 * pNode2->Coordinates();

        array_1d<double, 3> interpolation_point = coons_middle;

        if (TopAlpha > 0.0) {
            const array_1d<double, 3> edge_a =
                pNode1->Coordinates() - pNode0->Coordinates();
            const array_1d<double, 3> edge_b =
                pNode2->Coordinates() - pNode0->Coordinates();
            const array_1d<double, 3> face_normal =
                MathUtils<double>::CrossProduct(edge_a, edge_b);

            array_1d<double, 3> projected_point = ZeroVector(3);
            if (norm_2(face_normal) >
                std::numeric_limits<double>::epsilon() &&
                ProjectPointToSkinBoundaryAlongDirection(
                    rSkinSubModelPart,
                    coons_middle,
                    face_normal,
                    maximum_projection_distance,
                    projected_point)) {
                interpolation_point += TopAlpha * (
                    projected_point - coons_middle);
            }
        }

        array_1d<double, 3> middle_control_point =
            4.0 * interpolation_point;
        noalias(middle_control_point) -= 0.25 * edge_01[0]->Coordinates();
        noalias(middle_control_point) -= 0.50 * edge_01[1]->Coordinates();
        noalias(middle_control_point) -= 0.25 * edge_01[2]->Coordinates();
        noalias(middle_control_point) -= 0.50 * edge_02[1]->Coordinates();
        noalias(middle_control_point) -= 0.50 * edge_12[1]->Coordinates();
        noalias(middle_control_point) -= pNode2->Coordinates();

        NodePointerType p_middle_node(new Node(0, middle_control_point));

        std::vector<NodePointerType> top_control_nodes;
        top_control_nodes.reserve(9);
        top_control_nodes.push_back(edge_01[0]);
        top_control_nodes.push_back(edge_01[1]);
        top_control_nodes.push_back(edge_01[2]);
        top_control_nodes.push_back(edge_02[1]);
        top_control_nodes.push_back(p_middle_node);
        top_control_nodes.push_back(edge_12[1]);
        top_control_nodes.push_back(pNode2);
        top_control_nodes.push_back(pNode2);
        top_control_nodes.push_back(pNode2);

        return top_control_nodes;
    };

    auto top_has_curvature = [&](
        const CurvedTopFaceData& rTopData)
    {
        return rTopData.Alpha > 0.0 ||
               is_edge_curved(make_edge_key(
                   rTopData.CornerNodes[0],
                   rTopData.CornerNodes[1])) ||
               is_edge_curved(make_edge_key(
                   rTopData.CornerNodes[1],
                   rTopData.CornerNodes[2])) ||
               is_edge_curved(make_edge_key(
                   rTopData.CornerNodes[0],
                   rTopData.CornerNodes[2]));
    };

    auto refresh_top_face = [&](
        CurvedTopFaceData& rTopData)
    {
        rTopData.CurrentControlNodes = create_quadratic_top_controls(
            rTopData.CornerNodes[0],
            rTopData.CornerNodes[1],
            rTopData.CornerNodes[2],
            rTopData.Alpha);
        rTopData.IsCurved = top_has_curvature(rTopData);
        rTopData.pSurfaceGeometry = NurbsSurfaceType::Pointer();
    };

    auto register_top_face_incidence = [&](
        const CurvedTopFaceData& rTopData)
    {
        top_faces_by_edge[make_edge_key(
            rTopData.CornerNodes[0],
            rTopData.CornerNodes[1])].push_back(rTopData.Key);
        top_faces_by_edge[make_edge_key(
            rTopData.CornerNodes[1],
            rTopData.CornerNodes[2])].push_back(rTopData.Key);
        top_faces_by_edge[make_edge_key(
            rTopData.CornerNodes[0],
            rTopData.CornerNodes[2])].push_back(rTopData.Key);
    };

    auto refresh_top_faces_for_edge = [&](
        const SkinEdgeKey& rEdgeKey)
    {
        const auto incidence_it = top_faces_by_edge.find(rEdgeKey);
        if (incidence_it == top_faces_by_edge.end()) {
            return;
        }

        for (const auto& r_top_key : incidence_it->second) {
            auto top_it = mCurvedTopFaceRegistry.find(r_top_key);
            if (top_it != mCurvedTopFaceRegistry.end()) {
                refresh_top_face(top_it->second);
            }
        }
    };

    for (auto& r_type3_data : rType3CreationResult.VolumeQuadratureDataList) {
        if (!r_type3_data.RequiresQuadrature) {
            continue;
        }

        get_oriented_edge_controls(
            r_type3_data.pProjectionNode0,
            r_type3_data.pProjectionNode1);
        get_oriented_edge_controls(
            r_type3_data.pProjectionNode1,
            r_type3_data.pProjectionNode2);
        get_oriented_edge_controls(
            r_type3_data.pProjectionNode0,
            r_type3_data.pProjectionNode2);

        if (mGapApproximationOrder == 2) {
            const SkinTopFaceKey top_key = make_top_key(
                r_type3_data.pProjectionNode0,
                r_type3_data.pProjectionNode1,
                r_type3_data.pProjectionNode2);

            if (mCurvedTopFaceRegistry.find(top_key) ==
                mCurvedTopFaceRegistry.end()) {
                CurvedTopFaceData top_data;
                top_data.Key = top_key;
                top_data.CornerNodes = {{
                    r_type3_data.pProjectionNode0,
                    r_type3_data.pProjectionNode1,
                    r_type3_data.pProjectionNode2}};
                top_data.Alpha = 1.0;
                top_data.IsCurved = true;
                refresh_top_face(top_data);
                register_top_face_incidence(top_data);
                mCurvedTopFaceRegistry.emplace(top_key, std::move(top_data));
            }
        }
    }

    auto blend_control_nodes = [&](
        const std::vector<NodePointerType>& rLinearControlNodes,
        const std::vector<NodePointerType>& rCurvedControlNodes,
        const double Alpha)
    {
        KRATOS_ERROR_IF(rLinearControlNodes.size() != rCurvedControlNodes.size())
            << "[FinalizeType2AndType3GapGeometries] Cannot blend control "
            << "node vectors with different sizes.\n";

        std::vector<NodePointerType> blended_control_nodes;
        blended_control_nodes.reserve(rLinearControlNodes.size());

        for (std::size_t i = 0; i < rLinearControlNodes.size(); ++i) {
            KRATOS_ERROR_IF_NOT(rLinearControlNodes[i])
                << "[FinalizeType2AndType3GapGeometries] Null linear control node.\n";
            KRATOS_ERROR_IF_NOT(rCurvedControlNodes[i])
                << "[FinalizeType2AndType3GapGeometries] Null curved control node.\n";

            if (i == 0 || i + 1 == rLinearControlNodes.size()) {
                blended_control_nodes.push_back(rLinearControlNodes[i]);
                continue;
            }

            array_1d<double, 3> point = rLinearControlNodes[i]->Coordinates();
            noalias(point) += Alpha * (
                rCurvedControlNodes[i]->Coordinates() -
                rLinearControlNodes[i]->Coordinates());
            blended_control_nodes.push_back(NodePointerType(new Node(0, point)));
        }

        return blended_control_nodes;
    };

    auto reduce_edge_curvature = [&](const SkinEdgeKey& rKey)
    {
        auto edge_it = mCurvedEdgeRegistry.find(rKey);
        if (edge_it == mCurvedEdgeRegistry.end()) {
            return false;
        }

        auto& r_edge_data = edge_it->second;
        if (r_edge_data.Alpha <= 0.0) {
            return false;
        }

        r_edge_data.Alpha *= 0.5;
        if (r_edge_data.Alpha < 1.0e-6) {
            r_edge_data.Alpha = 0.0;
        }
        r_edge_data.IsCurved = r_edge_data.Alpha > 0.0;
        r_edge_data.CurrentControlNodes = blend_control_nodes(
            r_edge_data.LinearControlNodes,
            r_edge_data.CurvedControlNodes,
            r_edge_data.Alpha);
        mSkinEdgeControlNodes[rKey] = r_edge_data.CurrentControlNodes;
        if (!r_edge_data.IsCurved) {
            mLinearSkinEdges.insert(rKey);
        }
        refresh_top_faces_for_edge(rKey);

        return true;
    };

    auto reduce_top_curvature = [&](const SkinTopFaceKey& rKey)
    {
        auto top_it = mCurvedTopFaceRegistry.find(rKey);
        if (top_it == mCurvedTopFaceRegistry.end()) {
            return false;
        }

        auto& r_top_data = top_it->second;
        if (r_top_data.Alpha <= 0.0) {
            return false;
        }

        r_top_data.Alpha *= 0.5;
        if (r_top_data.Alpha < 1.0e-6) {
            r_top_data.Alpha = 0.0;
        }
        r_top_data.IsCurved = r_top_data.Alpha > 0.0;
        refresh_top_face(r_top_data);

        return true;
    };

    const std::size_t mapping_sampling_order =
        std::max<std::size_t>(3, mGapApproximationOrder + 2);
    IntegrationPointsArrayType mapping_sampling_points(
        mapping_sampling_order *
        mapping_sampling_order *
        mapping_sampling_order);
    auto mapping_sampling_point_it = mapping_sampling_points.begin();

    IntegrationPointUtilities::IntegrationPoints3D(
        mapping_sampling_point_it,
        mapping_sampling_order,
        mapping_sampling_order,
        mapping_sampling_order,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0);

    std::vector<CoordinatesArrayType> mapping_sampling_coordinates;
    mapping_sampling_coordinates.reserve(mapping_sampling_points.size());
    for (const auto& r_sampling_point : mapping_sampling_points) {
        CoordinatesArrayType local_coordinates = ZeroVector(3);
        local_coordinates[0] = r_sampling_point[0];
        local_coordinates[1] = r_sampling_point[1];
        local_coordinates[2] = r_sampling_point[2];
        mapping_sampling_coordinates.push_back(local_coordinates);
    }

    CoordinatesArrayType center_coordinates = ZeroVector(3);
    center_coordinates[0] = 0.5;
    center_coordinates[1] = 0.5;
    center_coordinates[2] = 0.5;

    auto is_volume_mapping_valid = [&](
        const NurbsVolumeType& rVolume)
    {
        Matrix center_jacobian;
        rVolume.Jacobian(center_jacobian, center_coordinates);
        const double center_determinant = MathUtils<double>::Det(center_jacobian);
        if (!std::isfinite(center_determinant) ||
            std::abs(center_determinant) <= 1.0e-14) {
            return false;
        }

        const double orientation_sign =
            center_determinant < 0.0 ? -1.0 : 1.0;

        for (const auto& r_local_coordinates : mapping_sampling_coordinates) {
            Matrix jacobian;
            rVolume.Jacobian(jacobian, r_local_coordinates);

            const double determinant =
                orientation_sign * MathUtils<double>::Det(jacobian);
            if (!std::isfinite(determinant) || determinant <= 1.0e-14) {
                return false;
            }
        }

        return true;
    };

    auto create_final_type2_volume = [&](
        const Type2VolumeQuadratureData& rType2Data)
    {
        if (mGapApproximationOrder == 1 ||
            !is_edge_curved(make_edge_key(
                rType2Data.pProjectionNode0,
                rType2Data.pProjectionNode1))) {
            return CreateType2LinearCollapsedEdgeVolume(
                rType2Data.pEdgeNode0,
                rType2Data.pEdgeNode1,
                rType2Data.pProjectionNode0,
                rType2Data.pProjectionNode1);
        }

        const auto edge_control_nodes = get_oriented_edge_controls(
            rType2Data.pProjectionNode0,
            rType2Data.pProjectionNode1);

        return CreateType2CollapsedEdgeVolume(
            rType2Data.pEdgeNode0,
            rType2Data.pEdgeNode1,
            edge_control_nodes);
    };

    auto create_final_type3_volume_data = [&](
        const Type3VolumeQuadratureData& rType3Data)
    {
        if (mGapApproximationOrder == 1) {
            return CreateType3LinearCollapsedCornerVolume(
                rType3Data.pSurrogateNode,
                rType3Data.pProjectionNode0,
                rType3Data.pProjectionNode1,
                rType3Data.pProjectionNode2,
                false);
        }

        const SkinTopFaceKey top_key = make_top_key(
            rType3Data.pProjectionNode0,
            rType3Data.pProjectionNode1,
            rType3Data.pProjectionNode2);
        const auto top_it = mCurvedTopFaceRegistry.find(top_key);
        KRATOS_ERROR_IF(top_it == mCurvedTopFaceRegistry.end())
            << "[FinalizeType2AndType3GapGeometries] Missing final top-face "
            << "registry entry for type3 volume.\n";

        if (!top_it->second.IsCurved) {
            return CreateType3LinearCollapsedCornerVolume(
                rType3Data.pSurrogateNode,
                rType3Data.pProjectionNode0,
                rType3Data.pProjectionNode1,
                rType3Data.pProjectionNode2,
                false);
        }

        return CreateType3CollapsedCornerVolumeFromTopControlNodes(
            rType3Data.pSurrogateNode,
            std::array<NodePointerType, 3>{{
                rType3Data.pProjectionNode0,
                rType3Data.pProjectionNode1,
                rType3Data.pProjectionNode2}},
            top_it->second.CurrentControlNodes,
            false);
    };

    bool all_mappings_valid_after_damping = false;
    for (std::size_t pass = 0; pass < 24; ++pass) {
        bool all_valid = true;
        std::unordered_set<SkinEdgeKey, SkinEdgeKeyHash> edge_keys_to_reduce;
        std::unordered_set<SkinTopFaceKey, SkinTopFaceKeyHash> top_keys_to_reduce;

        for (const auto& r_type2_data :
             rType2CreationResult.VolumeQuadratureDataList) {
            if (!r_type2_data.RequiresQuadrature) {
                continue;
            }

            const auto p_volume = create_final_type2_volume(r_type2_data);
            if (is_volume_mapping_valid(*p_volume)) {
                continue;
            }

            all_valid = false;
            edge_keys_to_reduce.insert(make_edge_key(
                r_type2_data.pProjectionNode0,
                r_type2_data.pProjectionNode1));
        }

        for (const auto& r_type3_data :
             rType3CreationResult.VolumeQuadratureDataList) {
            if (!r_type3_data.RequiresQuadrature) {
                continue;
            }

            const auto type3_volume_data =
                create_final_type3_volume_data(r_type3_data);
            if (is_volume_mapping_valid(*type3_volume_data.pVolume)) {
                continue;
            }

            all_valid = false;
            edge_keys_to_reduce.insert(make_edge_key(
                r_type3_data.pProjectionNode0,
                r_type3_data.pProjectionNode1));
            edge_keys_to_reduce.insert(make_edge_key(
                r_type3_data.pProjectionNode1,
                r_type3_data.pProjectionNode2));
            edge_keys_to_reduce.insert(make_edge_key(
                r_type3_data.pProjectionNode0,
                r_type3_data.pProjectionNode2));
            top_keys_to_reduce.insert(make_top_key(
                r_type3_data.pProjectionNode0,
                r_type3_data.pProjectionNode1,
                r_type3_data.pProjectionNode2));
        }

        if (all_valid) {
            all_mappings_valid_after_damping = true;
            break;
        }

        bool reduced_curvature = false;
        for (const auto& r_edge_key : edge_keys_to_reduce) {
            reduced_curvature |= reduce_edge_curvature(r_edge_key);
        }
        for (const auto& r_top_key : top_keys_to_reduce) {
            reduced_curvature |= reduce_top_curvature(r_top_key);
        }

        if (!reduced_curvature) {
            break;
        }
    }

    if (!all_mappings_valid_after_damping) {
        for (const auto& r_type2_data :
             rType2CreationResult.VolumeQuadratureDataList) {
            if (!r_type2_data.RequiresQuadrature) {
                continue;
            }

            const auto p_volume = create_final_type2_volume(r_type2_data);
            KRATOS_ERROR_IF_NOT(is_volume_mapping_valid(*p_volume))
                << "[FinalizeType2AndType3GapGeometries] Type2 volume mapping "
                << "remains invalid after high-order damping/fallback.\n"
                << "  surrogate edge nodes: "
                << r_type2_data.pEdgeNode0->Id() << ", "
                << r_type2_data.pEdgeNode1->Id() << "\n"
                << "  skin projection nodes: "
                << r_type2_data.pProjectionNode0->Id() << ", "
                << r_type2_data.pProjectionNode1->Id() << "\n";
        }

        for (const auto& r_type3_data :
             rType3CreationResult.VolumeQuadratureDataList) {
            if (!r_type3_data.RequiresQuadrature) {
                continue;
            }

            const auto type3_volume_data =
                create_final_type3_volume_data(r_type3_data);
            KRATOS_ERROR_IF_NOT(is_volume_mapping_valid(*type3_volume_data.pVolume))
                << "[FinalizeType2AndType3GapGeometries] Type3 volume mapping "
                << "remains invalid after high-order damping/fallback.\n"
                << "  surrogate node: "
                << r_type3_data.pSurrogateNode->Id() << "\n"
                << "  skin projection nodes: "
                << r_type3_data.pProjectionNode0->Id() << ", "
                << r_type3_data.pProjectionNode1->Id() << ", "
                << r_type3_data.pProjectionNode2->Id() << "\n";
        }
    }

    UpdateLateralFaceRegistryGeometries(
        rRootModelPart,
        rSkinSubModelPart);

    ModelPart& r_gap_type2_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType2Debug");
    ModelPart& r_gap_type3_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType3Debug");
    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    for (auto& r_type2_data : rType2CreationResult.VolumeQuadratureDataList) {
        if (!r_type2_data.RequiresQuadrature) {
            continue;
        }

        auto p_volume = create_final_type2_volume(r_type2_data);
        p_volume->SetId(next_geometry_id++);
        for (const auto& p_neighbour_geometry : r_type2_data.NeighbourGeometries) {
            AddUniqueNeighbourGeometry(*p_volume, p_neighbour_geometry);
        }
        if (mStoreGapDebugGeometries) {
            r_gap_type2_debug.AddGeometry(p_volume);
        }

        r_type2_data.pVolumeGeometry = p_volume;
        r_type2_data.VolumeQuadraturePointGeometries =
            CreateAndTagVolumeQuadraturePointGeometries(
                p_volume,
                r_type2_data.NeighbourGeometries.front(),
                IntegrationOrder,
                NumberOfShapeFunctionsDerivatives);
    }

    for (auto& r_type3_data : rType3CreationResult.VolumeQuadratureDataList) {
        if (!r_type3_data.RequiresQuadrature) {
            continue;
        }

        auto type3_volume_data = create_final_type3_volume_data(r_type3_data);
        auto p_volume = type3_volume_data.pVolume;
        p_volume->SetId(next_geometry_id++);
        for (const auto& p_neighbour_geometry : r_type3_data.NeighbourGeometries) {
            AddUniqueNeighbourGeometry(*p_volume, p_neighbour_geometry);
        }
        if (mStoreGapDebugGeometries) {
            r_gap_type3_debug.AddGeometry(p_volume);
        }

        r_type3_data.pVolumeGeometry = p_volume;
        r_type3_data.VolumeQuadraturePointGeometries =
            CreateAndTagVolumeQuadraturePointGeometries(
                p_volume,
                r_type3_data.NeighbourGeometries.front(),
                IntegrationOrder,
                NumberOfShapeFunctionsDerivatives);

        CheckType3QuadraturePointGeometries(
            r_type3_data.VolumeQuadraturePointGeometries,
            r_type3_data.pSurrogateNode,
            r_type3_data.pProjectionNode0,
            r_type3_data.pProjectionNode1,
            r_type3_data.pProjectionNode2);
    }

    mSkinProjectionTriangleData.clear();
}


#if 0
SnakeGapSbm3DUtilities::Type3CreationResult
SnakeGapSbm3DUtilities::CreateType3GapGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    const ExternalSpanDataMap& rExternalSpans,
    const KnotSpanGridInfo& rGridInfo,
    const Type2CreationResult& rType2CreationResult,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives)
{
    static_cast<void>(rSkinSubModelPart);

    Type3CreationResult result;
    result.Summary.NumberOfCandidateOpenFaces =
        rType2CreationResult.OpenFaceDataList.size();

    ModelPart& r_gap_type3_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType3Debug");
    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    struct LocalSpanEdgeKey {
        SpanKey3D Span0;
        SpanKey3D Span1;

        bool operator<(const LocalSpanEdgeKey& rOther) const
        {
            return std::tie(Span0, Span1) <
                   std::tie(rOther.Span0, rOther.Span1);
        }
    };

    struct LocalOpenFace {
        NodePointerType pSurrogateNode;
        SpanKey3D SpanA;
        SpanKey3D SpanB;
        NodePointerType pProjectionNodeA;
        NodePointerType pProjectionNodeB;
        LateralFaceOccurrence Type2Occurrence;
    };

    struct Type3TriangleCandidate {
        SpanKey3D SpanA;
        SpanKey3D SpanB;
        SpanKey3D SpanC;
        NodePointerType pProjectionNodeA;
        NodePointerType pProjectionNodeB;
        NodePointerType pProjectionNodeC;
        std::array<LocalSpanEdgeKey, 3> Edges;
        std::vector<std::size_t> CoveredOpenEdgeIndices;
        bool HasType3Span = false;
    };

    struct Type3VolumeKey {
        IndexType SurrogateNodeId = 0;
        std::array<IndexType, 3> ProjectionNodeIds = {0, 0, 0};

        bool operator<(const Type3VolumeKey& rOther) const
        {
            return std::tie(SurrogateNodeId, ProjectionNodeIds) <
                   std::tie(rOther.SurrogateNodeId, rOther.ProjectionNodeIds);
        }
    };

    std::map<IndexType, std::map<LocalSpanEdgeKey, LocalOpenFace>>
        open_faces_by_surrogate_node_id;
    std::set<Type3VolumeKey> created_volume_keys;

    // Freeze the type 2 open-face graph before type 3 registrations modify it.
    for (const auto& r_open_face_data : rType2CreationResult.OpenFaceDataList) {
        KRATOS_ERROR_IF_NOT(r_open_face_data.pSurrogateNode)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Type 2 open face has null surrogate node.\n";

        const auto span_a_it =
            rExternalSpans.find(r_open_face_data.ExternalSpan0);
        const auto span_b_it =
            rExternalSpans.find(r_open_face_data.ExternalSpan1);
        KRATOS_ERROR_IF(span_a_it == rExternalSpans.end())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Missing external span A "
            << SpanToString(r_open_face_data.ExternalSpan0) << ".\n";
        KRATOS_ERROR_IF(span_b_it == rExternalSpans.end())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Missing external span B "
            << SpanToString(r_open_face_data.ExternalSpan1) << ".\n";
        KRATOS_ERROR_IF_NOT(span_a_it->second.HasProjectionNode())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "External span A has no projection node.\n";
        KRATOS_ERROR_IF_NOT(span_b_it->second.HasProjectionNode())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "External span B has no projection node.\n";

        SpanKey3D span_a = r_open_face_data.ExternalSpan0;
        SpanKey3D span_b = r_open_face_data.ExternalSpan1;
        NodePointerType p_projection_node_a =
            span_a_it->second.pProjectionNode;
        NodePointerType p_projection_node_b =
            span_b_it->second.pProjectionNode;

        const auto open_face_key = MakeCanonicalFaceKey3D(
            r_open_face_data.pSurrogateNode->Id(),
            p_projection_node_a->Id(),
            p_projection_node_b->Id());
        const auto open_face_registry_it =
            mLateralFaceRegistry.find(open_face_key);

        if (open_face_registry_it == mLateralFaceRegistry.end() ||
            open_face_registry_it->second.size() != 1 ||
            open_face_registry_it->second.front().GapType != 2) {
            continue;
        }

        const LateralFaceOccurrence type2_occurrence =
            open_face_registry_it->second.front();
        KRATOS_ERROR_IF_NOT(type2_occurrence.pNeighbourGeometry)
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Open type 2 face has null neighbour geometry.\n";

        LocalSpanEdgeKey local_edge_key;
        local_edge_key.Span0 = span_a;
        local_edge_key.Span1 = span_b;
        if (local_edge_key.Span1 < local_edge_key.Span0) {
            std::swap(local_edge_key.Span0, local_edge_key.Span1);
        }

        auto& r_local_edges = open_faces_by_surrogate_node_id[
            r_open_face_data.pSurrogateNode->Id()];
        if (r_local_edges.find(local_edge_key) != r_local_edges.end()) {
            continue;
        }

        LocalOpenFace local_open_face;
        local_open_face.pSurrogateNode =
            r_open_face_data.pSurrogateNode;
        local_open_face.SpanA = span_a;
        local_open_face.SpanB = span_b;
        local_open_face.pProjectionNodeA = p_projection_node_a;
        local_open_face.pProjectionNodeB = p_projection_node_b;
        local_open_face.Type2Occurrence = type2_occurrence;
        r_local_edges.emplace(
            local_edge_key,
            std::move(local_open_face));
    }

    const double minimum_span_size = std::min({
        rGridInfo.SpanSizeX,
        rGridInfo.SpanSizeY,
        rGridInfo.SpanSizeZ});
    const double determinant_tolerance =
        6.0e-12 * minimum_span_size * minimum_span_size *
        minimum_span_size;

    for (auto& r_node_open_faces : open_faces_by_surrogate_node_id) {
        auto& r_open_edges = r_node_open_faces.second;
        KRATOS_ERROR_IF(r_open_edges.empty())
            << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
            << "Empty local type 2 face graph for surrogate node "
            << r_node_open_faces.first << ".\n";

        const NodePointerType p_surrogate_node =
            r_open_edges.begin()->second.pSurrogateNode;
        const GridPointKey3D surrogate_grid_node =
            ComputeGridPointKey(
                p_surrogate_node->Coordinates(),
                rGridInfo);

        std::map<SpanKey3D, const ExternalSpanData*>
            incident_external_spans;
        std::vector<SpanKey3D> incident_span_keys;

        for (int di = -1; di <= 0; ++di) {
            for (int dj = -1; dj <= 0; ++dj) {
                for (int dk = -1; dk <= 0; ++dk) {
                    const SpanKey3D incident_span{
                        surrogate_grid_node.I + di,
                        surrogate_grid_node.J + dj,
                        surrogate_grid_node.K + dk};

                    if (!IsValidSpan(incident_span, rGridInfo)) {
                        continue;
                    }

                    const auto external_span_it =
                        rExternalSpans.find(incident_span);
                    if (external_span_it == rExternalSpans.end() ||
                        !external_span_it->second.HasProjectionNode()) {
                        continue;
                    }

                    incident_external_spans.emplace(
                        incident_span,
                        &external_span_it->second);
                    incident_span_keys.push_back(incident_span);
                }
            }
        }

        for (const auto& r_open_edge_entry : r_open_edges) {
            const auto& r_open_face = r_open_edge_entry.second;
            KRATOS_ERROR_IF(
                incident_external_spans.find(r_open_face.SpanA) ==
                    incident_external_spans.end() ||
                incident_external_spans.find(r_open_face.SpanB) ==
                    incident_external_spans.end())
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Type 2 open face references external spans that are not "
                << "incident to surrogate node " << p_surrogate_node->Id()
                << ": " << SpanToString(r_open_face.SpanA) << " and "
                << SpanToString(r_open_face.SpanB) << ".\n";
        }

        std::vector<LocalSpanEdgeKey> open_edge_keys;
        open_edge_keys.reserve(r_open_edges.size());
        std::map<LocalSpanEdgeKey, std::size_t> open_edge_index_by_key;
        for (const auto& r_open_edge_entry : r_open_edges) {
            const std::size_t open_edge_index = open_edge_keys.size();
            open_edge_keys.push_back(r_open_edge_entry.first);
            open_edge_index_by_key.emplace(
                r_open_edge_entry.first,
                open_edge_index);
        }

        auto make_local_edge_key = [](
            const SpanKey3D& rSpan0,
            const SpanKey3D& rSpan1)
        {
            LocalSpanEdgeKey edge_key;
            edge_key.Span0 = rSpan0;
            edge_key.Span1 = rSpan1;
            if (edge_key.Span1 < edge_key.Span0) {
                std::swap(edge_key.Span0, edge_key.Span1);
            }
            return edge_key;
        };

        auto calculate_oriented_tetra_determinant = [&p_surrogate_node](
            const NodePointerType& pProjectionNodeA,
            const NodePointerType& pProjectionNodeB,
            const NodePointerType& pProjectionNodeC)
        {
            const array_1d<double, 3> edge_a =
                pProjectionNodeA->Coordinates() -
                p_surrogate_node->Coordinates();
            const array_1d<double, 3> edge_b =
                pProjectionNodeB->Coordinates() -
                p_surrogate_node->Coordinates();
            const array_1d<double, 3> edge_c =
                pProjectionNodeC->Coordinates() -
                p_surrogate_node->Coordinates();
            const array_1d<double, 3> cross_bc =
                MathUtils<double>::CrossProduct(edge_b, edge_c);
            return inner_prod(edge_a, cross_bc);
        };

        std::vector<Type3TriangleCandidate> triangle_candidates;
        for (std::size_t first = 0; first < incident_span_keys.size(); ++first) {
            for (std::size_t second = first + 1;
                 second < incident_span_keys.size();
                 ++second) {
                for (std::size_t third = second + 1;
                     third < incident_span_keys.size();
                     ++third) {
                    Type3TriangleCandidate candidate;
                    candidate.SpanA = incident_span_keys[first];
                    candidate.SpanB = incident_span_keys[second];
                    candidate.SpanC = incident_span_keys[third];

                    const auto span_a_it =
                        incident_external_spans.find(candidate.SpanA);
                    const auto span_b_it =
                        incident_external_spans.find(candidate.SpanB);
                    const auto span_c_it =
                        incident_external_spans.find(candidate.SpanC);
                    KRATOS_ERROR_IF(
                        span_a_it == incident_external_spans.end() ||
                        span_b_it == incident_external_spans.end() ||
                        span_c_it == incident_external_spans.end())
                        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                        << "Internal error while building local type 3 "
                        << "triangulation candidates.\n";

                    candidate.pProjectionNodeA =
                        span_a_it->second->pProjectionNode;
                    candidate.pProjectionNodeB =
                        span_b_it->second->pProjectionNode;
                    candidate.pProjectionNodeC =
                        span_c_it->second->pProjectionNode;
                    candidate.Edges = {{
                        make_local_edge_key(candidate.SpanA, candidate.SpanB),
                        make_local_edge_key(candidate.SpanB, candidate.SpanC),
                        make_local_edge_key(candidate.SpanC, candidate.SpanA)}};
                    candidate.HasType3Span =
                        span_a_it->second->Type == GapSpanType::Type3 ||
                        span_b_it->second->Type == GapSpanType::Type3 ||
                        span_c_it->second->Type == GapSpanType::Type3;

                    for (const auto& r_edge_key : candidate.Edges) {
                        const auto open_edge_index_it =
                            open_edge_index_by_key.find(r_edge_key);
                        if (open_edge_index_it !=
                            open_edge_index_by_key.end()) {
                            candidate.CoveredOpenEdgeIndices.push_back(
                                open_edge_index_it->second);
                        }
                    }

                    if (candidate.CoveredOpenEdgeIndices.empty()) {
                        continue;
                    }

                    const double determinant =
                        calculate_oriented_tetra_determinant(
                            candidate.pProjectionNodeA,
                            candidate.pProjectionNodeB,
                            candidate.pProjectionNodeC);
                    if (std::abs(determinant) <= determinant_tolerance) {
                        continue;
                    }

                    triangle_candidates.push_back(std::move(candidate));
                }
            }
        }

        std::sort(
            triangle_candidates.begin(),
            triangle_candidates.end(),
            [](const Type3TriangleCandidate& rLeft,
               const Type3TriangleCandidate& rRight)
            {
                if (rLeft.HasType3Span != rRight.HasType3Span) {
                    return rLeft.HasType3Span;
                }
                if (rLeft.CoveredOpenEdgeIndices.size() !=
                    rRight.CoveredOpenEdgeIndices.size()) {
                    return rLeft.CoveredOpenEdgeIndices.size() >
                           rRight.CoveredOpenEdgeIndices.size();
                }

                const std::array<SpanKey3D, 3> left_spans = {{
                    rLeft.SpanA,
                    rLeft.SpanB,
                    rLeft.SpanC}};
                const std::array<SpanKey3D, 3> right_spans = {{
                    rRight.SpanA,
                    rRight.SpanB,
                    rRight.SpanC}};
                return left_spans < right_spans;
            });

        std::vector<std::vector<std::size_t>>
            candidate_indices_by_open_edge(open_edge_keys.size());
        for (std::size_t candidate_index = 0;
             candidate_index < triangle_candidates.size();
             ++candidate_index) {
            for (const std::size_t open_edge_index :
                 triangle_candidates[candidate_index].CoveredOpenEdgeIndices) {
                candidate_indices_by_open_edge[open_edge_index].push_back(
                    candidate_index);
            }
        }

        const std::size_t invalid_candidate_index =
            std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> open_edge_covered_by(
            open_edge_keys.size(),
            invalid_candidate_index);
        std::map<LocalSpanEdgeKey, std::size_t> selected_edge_usage;
        std::vector<std::size_t> selected_triangle_indices;

        auto can_select_candidate = [&](
            const Type3TriangleCandidate& rCandidate)
        {
            for (const std::size_t open_edge_index :
                 rCandidate.CoveredOpenEdgeIndices) {
                if (open_edge_covered_by[open_edge_index] !=
                    invalid_candidate_index) {
                    return false;
                }
            }

            for (const auto& r_edge_key : rCandidate.Edges) {
                const bool is_type2_open_edge =
                    open_edge_index_by_key.find(r_edge_key) !=
                    open_edge_index_by_key.end();
                const std::size_t maximum_usage =
                    is_type2_open_edge ? 1u : 2u;
                const auto usage_it = selected_edge_usage.find(r_edge_key);
                const std::size_t current_usage =
                    usage_it == selected_edge_usage.end() ?
                    0u : usage_it->second;
                if (current_usage >= maximum_usage) {
                    return false;
                }
            }

            return true;
        };

        auto add_candidate_usage = [&](
            const Type3TriangleCandidate& rCandidate)
        {
            for (const std::size_t open_edge_index :
                 rCandidate.CoveredOpenEdgeIndices) {
                open_edge_covered_by[open_edge_index] =
                    selected_triangle_indices.back();
            }
            for (const auto& r_edge_key : rCandidate.Edges) {
                ++selected_edge_usage[r_edge_key];
            }
        };

        auto remove_candidate_usage = [&](
            const Type3TriangleCandidate& rCandidate)
        {
            for (const std::size_t open_edge_index :
                 rCandidate.CoveredOpenEdgeIndices) {
                open_edge_covered_by[open_edge_index] =
                    invalid_candidate_index;
            }
            for (const auto& r_edge_key : rCandidate.Edges) {
                auto usage_it = selected_edge_usage.find(r_edge_key);
                KRATOS_ERROR_IF(usage_it == selected_edge_usage.end())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Internal type 3 triangulation usage bookkeeping "
                    << "error.\n";
                if (--usage_it->second == 0) {
                    selected_edge_usage.erase(usage_it);
                }
            }
        };

        std::function<bool()> select_triangles = [&]()
        {
            std::size_t selected_open_edge_index =
                invalid_candidate_index;
            std::vector<std::size_t> usable_candidate_indices;

            for (std::size_t open_edge_index = 0;
                 open_edge_index < open_edge_keys.size();
                 ++open_edge_index) {
                if (open_edge_covered_by[open_edge_index] !=
                    invalid_candidate_index) {
                    continue;
                }

                std::vector<std::size_t> current_usable_candidate_indices;
                for (const std::size_t candidate_index :
                     candidate_indices_by_open_edge[open_edge_index]) {
                    if (can_select_candidate(
                        triangle_candidates[candidate_index])) {
                        current_usable_candidate_indices.push_back(
                            candidate_index);
                    }
                }

                if (current_usable_candidate_indices.empty()) {
                    return false;
                }

                if (selected_open_edge_index == invalid_candidate_index ||
                    current_usable_candidate_indices.size() <
                    usable_candidate_indices.size()) {
                    selected_open_edge_index = open_edge_index;
                    usable_candidate_indices =
                        std::move(current_usable_candidate_indices);
                }
            }

            if (selected_open_edge_index == invalid_candidate_index) {
                return true;
            }

            for (const std::size_t candidate_index :
                 usable_candidate_indices) {
                const auto& r_candidate =
                    triangle_candidates[candidate_index];
                selected_triangle_indices.push_back(candidate_index);
                add_candidate_usage(r_candidate);
                if (select_triangles()) {
                    return true;
                }
                remove_candidate_usage(r_candidate);
                selected_triangle_indices.pop_back();
            }

            return false;
        };

        if (!select_triangles()) {
            std::ostringstream error_message;
            error_message
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Could not build a local type 3 triangulation that closes "
                << "all initial type 2 open faces exactly once.\n"
                << "  surrogate node: " << p_surrogate_node->Id() << "\n"
                << "  open edges:\n";
            for (const auto& r_open_edge_entry : r_open_edges) {
                const auto& r_open_face = r_open_edge_entry.second;
                error_message
                    << "    " << SpanToString(r_open_face.SpanA)
                    << " -- " << SpanToString(r_open_face.SpanB)
                    << " projection nodes "
                    << r_open_face.pProjectionNodeA->Id()
                    << ", "
                    << r_open_face.pProjectionNodeB->Id()
                    << "\n";
            }
            KRATOS_ERROR << error_message.str();
        }

        for (const auto& r_usage_entry : selected_edge_usage) {
            const auto& r_edge_key = r_usage_entry.first;
            const std::size_t usage = r_usage_entry.second;

            const bool is_initial_type2_open_edge =
                open_edge_index_by_key.find(r_edge_key) != open_edge_index_by_key.end();

            if (is_initial_type2_open_edge) {
                KRATOS_ERROR_IF(usage != 1)
                    << "Initial type2 open edge not used exactly once.\n"
                    << "  edge: " << SpanToString(r_edge_key.Span0)
                    << " -- " << SpanToString(r_edge_key.Span1)
                    << "  usage: " << usage << "\n";
            } else {
                KRATOS_ERROR_IF(usage != 2)
                    << "Internal type3 edge used only once. This creates an artificial "
                    << "open type3 face and will generate a wrong Dirichlet condition.\n"
                    << "  edge: " << SpanToString(r_edge_key.Span0)
                    << " -- " << SpanToString(r_edge_key.Span1)
                    << "  usage: " << usage << "\n";
            }
        }

        auto get_open_face_for_edge = [&](
            const LocalSpanEdgeKey& rEdgeKey) -> const LocalOpenFace*
        {
            const auto open_edge_it = r_open_edges.find(rEdgeKey);
            if (open_edge_it == r_open_edges.end()) {
                return nullptr;
            }
            return &open_edge_it->second;
        };

        auto get_face_occurrence = [&](
            const LocalSpanEdgeKey& rEdgeKey,
            const LateralFaceOccurrence& rFallbackOccurrence,
            const SpanKey3D& rFallbackExternalSpan)
        {
            const auto p_open_face = get_open_face_for_edge(rEdgeKey);
            if (p_open_face) {
                return p_open_face->Type2Occurrence;
            }

            LateralFaceOccurrence occurrence = rFallbackOccurrence;
            occurrence.ExternalSpan = rFallbackExternalSpan;
            return occurrence;
        };

        for (const std::size_t selected_triangle_index :
             selected_triangle_indices) {
            const auto& r_candidate =
                triangle_candidates[selected_triangle_index];

            KRATOS_ERROR_IF(r_candidate.CoveredOpenEdgeIndices.empty())
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Selected type 3 triangle does not cover an initial "
                << "type 2 open edge.\n";

            const auto primary_open_edge_key =
                open_edge_keys[r_candidate.CoveredOpenEdgeIndices.front()];
            const auto primary_open_edge_it =
                r_open_edges.find(primary_open_edge_key);
            KRATOS_ERROR_IF(primary_open_edge_it == r_open_edges.end())
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Internal error: selected type 3 triangle references a "
                << "missing open type 2 edge.\n";

            const LateralFaceOccurrence primary_occurrence =
                primary_open_edge_it->second.Type2Occurrence;
            KRATOS_ERROR_IF_NOT(primary_occurrence.pNeighbourGeometry)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Type 2 graph edge has null neighbour geometry.\n";

            SpanKey3D span_a = r_candidate.SpanA;
            SpanKey3D span_b = r_candidate.SpanB;
            SpanKey3D span_c = r_candidate.SpanC;
            NodePointerType p_projection_node_a =
                r_candidate.pProjectionNodeA;
            NodePointerType p_projection_node_b =
                r_candidate.pProjectionNodeB;
            NodePointerType p_projection_node_c =
                r_candidate.pProjectionNodeC;

            double determinant =
                calculate_oriented_tetra_determinant(
                    p_projection_node_a,
                    p_projection_node_b,
                    p_projection_node_c);

            if (std::abs(determinant) <= determinant_tolerance) {
                ++result.Summary.NumberOfSkippedFaces;
                KRATOS_WARNING("SnakeGapSbm3DUtilities")
                    << "Skipping degenerate type3 tetra at surrogate node "
                    << p_surrogate_node->Id()
                    << " spans "
                    << SpanToString(span_a) << ", "
                    << SpanToString(span_b) << ", "
                    << SpanToString(span_c)
                    << " projection nodes "
                    << p_projection_node_a->Id() << ", "
                    << p_projection_node_b->Id() << ", "
                    << p_projection_node_c->Id()
                    << " determinant = " << determinant << "\n";
                continue;
            }

            if (determinant < 0.0) {
                std::swap(p_projection_node_a, p_projection_node_b);
                std::swap(span_a, span_b);
                determinant = -determinant;
            }

            Type3VolumeKey volume_key;
            volume_key.SurrogateNodeId = p_surrogate_node->Id();
            volume_key.ProjectionNodeIds = {{
                p_projection_node_a->Id(),
                p_projection_node_b->Id(),
                p_projection_node_c->Id()}};
            std::sort(
                volume_key.ProjectionNodeIds.begin(),
                volume_key.ProjectionNodeIds.end());

            KRATOS_ERROR_IF(created_volume_keys.find(volume_key) !=
                            created_volume_keys.end())
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "The local type 3 triangulation selected a duplicated "
                << "tetrahedron at surrogate node "
                << p_surrogate_node->Id() << " with projection node ids "
                << volume_key.ProjectionNodeIds[0] << ", "
                << volume_key.ProjectionNodeIds[1] << ", "
                << volume_key.ProjectionNodeIds[2] << ".\n";

            auto p_type3_volume =
                CreateType3CollapsedCornerVolume(
                    p_surrogate_node,
                    p_projection_node_a,
                    p_projection_node_b,
                    p_projection_node_c);

            GeometriesArrayType volume_quadrature_point_geometries =
                CreateAndTagVolumeQuadraturePointGeometries(
                    p_type3_volume,
                    primary_occurrence.pNeighbourGeometry,
                    IntegrationOrder,
                    NumberOfShapeFunctionsDerivatives);

            KRATOS_ERROR_IF(
                volume_quadrature_point_geometries.size() == 0)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "A non-degenerate type 3 tetrahedron generated no "
                << "quadrature point geometries.\n";

            CheckType3QuadraturePointGeometries(
                volume_quadrature_point_geometries,
                p_surrogate_node,
                p_projection_node_a,
                p_projection_node_b,
                p_projection_node_c);

            std::vector<Geometry<Node>::Pointer> neighbour_geometries;
            for (const std::size_t open_edge_index :
                 r_candidate.CoveredOpenEdgeIndices) {
                const auto open_edge_it =
                    r_open_edges.find(open_edge_keys[open_edge_index]);
                KRATOS_ERROR_IF(open_edge_it == r_open_edges.end())
                    << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                    << "Internal error: selected type 3 triangle references "
                    << "a missing open type 2 edge.\n";
                AddUniqueNeighbourGeometry(
                    neighbour_geometries,
                    open_edge_it->second.Type2Occurrence.pNeighbourGeometry);
            }

            created_volume_keys.insert(volume_key);

            p_type3_volume->SetId(next_geometry_id++);
            for (const auto& p_neighbour_geometry : neighbour_geometries) {
                AddNeighbourGeometry(*p_type3_volume, p_neighbour_geometry);
            }
            if (mStoreGapDebugGeometries) {
                r_gap_type3_debug.AddGeometry(p_type3_volume);
            }

            const auto edge_ab = make_local_edge_key(span_a, span_b);
            const auto edge_bc = make_local_edge_key(span_b, span_c);
            const auto edge_ca = make_local_edge_key(span_c, span_a);

            const auto occurrence_ab =
                get_face_occurrence(edge_ab, primary_occurrence, span_a);
            const auto occurrence_bc =
                get_face_occurrence(edge_bc, primary_occurrence, span_b);
            const auto occurrence_ca =
                get_face_occurrence(edge_ca, primary_occurrence, span_c);

            RegisterType3LateralFace(
                r_gap_type3_debug,
                next_geometry_id,
                occurrence_ab.SurrogateConditionId,
                occurrence_ab.ExternalSpan,
                p_surrogate_node,
                p_projection_node_a,
                p_projection_node_b,
                p_projection_node_c,
                occurrence_ab.pNeighbourGeometry);

            RegisterType3LateralFace(
                r_gap_type3_debug,
                next_geometry_id,
                occurrence_bc.SurrogateConditionId,
                occurrence_bc.ExternalSpan,
                p_surrogate_node,
                p_projection_node_b,
                p_projection_node_c,
                p_projection_node_a,
                occurrence_bc.pNeighbourGeometry);

            RegisterType3LateralFace(
                r_gap_type3_debug,
                next_geometry_id,
                occurrence_ca.SurrogateConditionId,
                occurrence_ca.ExternalSpan,
                p_surrogate_node,
                p_projection_node_c,
                p_projection_node_a,
                p_projection_node_b,
                occurrence_ca.pNeighbourGeometry);

            RegisterType3LateralFace(
                r_gap_type3_debug,
                next_geometry_id,
                primary_occurrence.SurrogateConditionId,
                span_c,
                p_projection_node_a,
                p_projection_node_b,
                p_projection_node_c,
                p_surrogate_node,
                primary_occurrence.pNeighbourGeometry);

            Type3VolumeQuadratureData type3_volume_data;
            type3_volume_data.VolumeQuadraturePointGeometries =
                std::move(volume_quadrature_point_geometries);
            type3_volume_data.NeighbourGeometries =
                std::move(neighbour_geometries);
            type3_volume_data.CharacteristicLength = std::max({
                norm_2(
                    p_projection_node_a->Coordinates() -
                    p_surrogate_node->Coordinates()),
                norm_2(
                    p_projection_node_b->Coordinates() -
                    p_surrogate_node->Coordinates()),
                norm_2(
                    p_projection_node_c->Coordinates() -
                    p_surrogate_node->Coordinates())});
            type3_volume_data.SurrogateNodeId =
                p_surrogate_node->Id();
            type3_volume_data.ProjectionNodeId0 =
                p_projection_node_a->Id();
            type3_volume_data.ProjectionNodeId1 =
                p_projection_node_b->Id();
            type3_volume_data.ProjectionNodeId2 =
                p_projection_node_c->Id();

            result.VolumeQuadratureDataList.push_back(
                std::move(type3_volume_data));
            ++result.Summary.NumberOfCreatedVolumes;
        }
    }

    for (const auto& r_node_open_faces : open_faces_by_surrogate_node_id) {
        for (const auto& r_open_edge_entry : r_node_open_faces.second) {
            const auto& r_open_face = r_open_edge_entry.second;
            const auto face_key = MakeCanonicalFaceKey3D(
                r_open_face.pSurrogateNode->Id(),
                r_open_face.pProjectionNodeA->Id(),
                r_open_face.pProjectionNodeB->Id());
            const auto registry_it = mLateralFaceRegistry.find(face_key);

            std::size_t number_of_type2_occurrences = 0;
            std::size_t number_of_type3_occurrences = 0;
            std::size_t number_of_occurrences = 0;
            if (registry_it != mLateralFaceRegistry.end()) {
                number_of_occurrences = registry_it->second.size();
                for (const auto& r_occurrence : registry_it->second) {
                    if (r_occurrence.GapType == 2) {
                        ++number_of_type2_occurrences;
                    } else if (r_occurrence.GapType == 3) {
                        ++number_of_type3_occurrences;
                    }
                }
            }

            KRATOS_ERROR_IF(
                number_of_occurrences != 2 ||
                number_of_type2_occurrences != 1 ||
                number_of_type3_occurrences != 1)
                << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
                << "Initial type 2 open face was not closed by exactly one "
                << "type 3 occurrence.\n"
                << "  surrogate node: "
                << r_open_face.pSurrogateNode->Id() << "\n"
                << "  span A: " << SpanToString(r_open_face.SpanA) << "\n"
                << "  span B: " << SpanToString(r_open_face.SpanB) << "\n"
                << "  projection node A: "
                << r_open_face.pProjectionNodeA->Id() << "\n"
                << "  projection node B: "
                << r_open_face.pProjectionNodeB->Id() << "\n"
                << "  registry occurrences: "
                << number_of_occurrences << "\n"
                << "  GapType 2 occurrences: "
                << number_of_type2_occurrences << "\n"
                << "  GapType 3 occurrences: "
                << number_of_type3_occurrences << "\n";
        }
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type 3 creation summary:\n"
        << "  candidate open faces: "
        << result.Summary.NumberOfCandidateOpenFaces << "\n"
        << "  created volumes: "
        << result.Summary.NumberOfCreatedVolumes << "\n"
        << "  skipped faces: "
        << result.Summary.NumberOfSkippedFaces << "\n"
        << "  GapType3Debug geometries: "
        << r_gap_type3_debug.NumberOfGeometries() << "\n";


    //debug
    std::unordered_set<IndexType> surrogate_node_ids;

    for (const auto& r_open_face_data : rType2CreationResult.OpenFaceDataList) {
        surrogate_node_ids.insert(r_open_face_data.pSurrogateNode->Id());
    }
    std::size_t open_type3_faces = 0;
    std::size_t open_type3_lateral_faces = 0;

    for (const auto& r_entry : mLateralFaceRegistry) {
        const auto& r_face_key = r_entry.first;
        const auto& r_occurrences = r_entry.second;

        if (r_occurrences.size() != 1) {
            continue;
        }

        const auto& r_occurrence = r_occurrences.front();

        if (r_occurrence.GapType != 3) {
            continue;
        }

        ++open_type3_faces;

        const bool has_surrogate_node =
            surrogate_node_ids.find(r_face_key.NodeIds[0]) != surrogate_node_ids.end() ||
            surrogate_node_ids.find(r_face_key.NodeIds[1]) != surrogate_node_ids.end() ||
            surrogate_node_ids.find(r_face_key.NodeIds[2]) != surrogate_node_ids.end();

        if (has_surrogate_node) {
            ++open_type3_lateral_faces;

            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "[CreateType3GapGeometries] Open type3 lateral face detected. "
                << "This face contains a surrogate node and will later receive a "
                << "Dirichlet condition, but it is probably an internal face left "
                << "open by an incomplete type3 triangulation.\n"
                << "  face node ids: "
                << r_face_key.NodeIds[0] << ", "
                << r_face_key.NodeIds[1] << ", "
                << r_face_key.NodeIds[2] << "\n"
                << "  occurrence gap type: " << r_occurrence.GapType << "\n";
        }
    }

    KRATOS_ERROR_IF(open_type3_lateral_faces > 0)
        << "[SnakeGapSbm3DUtilities::CreateType3GapGeometries] "
        << "There are open type3 lateral faces containing surrogate nodes. "
        << "This means the type3 triangulation is not closed and Dirichlet BCs "
        << "would be applied on internal faces. Number of bad faces: "
        << open_type3_lateral_faces << "\n";

    return result;
}

#endif

#if 0
SnakeGapSbm3DUtilities::Type3CreationResult
SnakeGapSbm3DUtilities::CreateType3GapGeometries(
    ModelPart& rRootModelPart,
    ModelPart& rSkinSubModelPart,
    const ExternalSpanDataMap& rExternalSpans,
    const KnotSpanGridInfo& rGridInfo,
    const Type2CreationResult& rType2CreationResult,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives)
{
    static_cast<void>(rSkinSubModelPart);

    Type3CreationResult result;
    result.Summary.NumberOfCandidateOpenFaces =
        rType2CreationResult.OpenFaceDataList.size();

    ModelPart& r_gap_type3_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType3Debug");

    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    struct LocalOpenFace
    {
        NodePointerType pSurrogateNode;
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;
        CanonicalFaceKey3D RegistryKey;
    };

    struct Type3VolumeKey
    {
        IndexType SurrogateNodeId = 0;
        std::array<IndexType, 3> ProjectionNodeIds = {{0, 0, 0}};

        bool operator<(const Type3VolumeKey& rOther) const
        {
            return std::tie(SurrogateNodeId, ProjectionNodeIds) <
                   std::tie(rOther.SurrogateNodeId, rOther.ProjectionNodeIds);
        }

        static Type3VolumeKey Create(
            const NodePointerType& pSurrogateNode,
            const NodePointerType& pProjectionNode0,
            const NodePointerType& pProjectionNode1,
            const NodePointerType& pProjectionNode2)
        {
            Type3VolumeKey key;
            key.SurrogateNodeId = pSurrogateNode->Id();
            key.ProjectionNodeIds = {{
                pProjectionNode0->Id(),
                pProjectionNode1->Id(),
                pProjectionNode2->Id()}};
            std::sort(
                key.ProjectionNodeIds.begin(),
                key.ProjectionNodeIds.end());
            return key;
        }
    };

    struct LocalSideFace
    {
        NodePointerType pProjectionNode0;
        NodePointerType pProjectionNode1;
        NodePointerType pOppositeProjectionNode;
    };

    struct Point2D
    {
        double X = 0.0;
        double Y = 0.0;
    };

    struct LocalTriangle
    {
        std::size_t Index0 = 0;
        std::size_t Index1 = 0;
        std::size_t Index2 = 0;
    };

    struct PlannedTriangleCandidate
    {
        LocalTriangle Triangle;
        std::array<CanonicalFaceKey3D, 4> FaceKeys;
        Type3VolumeKey VolumeKey;
    };

    struct PolygonUtilities
    {
        double Cross(
            const Point2D& rA,
            const Point2D& rB,
            const Point2D& rC) const
        {
            return (rB.X - rA.X) * (rC.Y - rA.Y) -
                   (rB.Y - rA.Y) * (rC.X - rA.X);
        }

        bool IsPointOnSegment(
            const Point2D& rA,
            const Point2D& rB,
            const Point2D& rPoint,
            const double Tolerance) const
        {
            if (std::abs(Cross(rA, rB, rPoint)) > Tolerance) {
                return false;
            }

            return rPoint.X >= std::min(rA.X, rB.X) - Tolerance &&
                   rPoint.X <= std::max(rA.X, rB.X) + Tolerance &&
                   rPoint.Y >= std::min(rA.Y, rB.Y) - Tolerance &&
                   rPoint.Y <= std::max(rA.Y, rB.Y) + Tolerance;
        }

        bool DoSegmentsIntersect(
            const Point2D& rA0,
            const Point2D& rA1,
            const Point2D& rB0,
            const Point2D& rB1,
            const double Tolerance) const
        {
            const double c0 = Cross(rA0, rA1, rB0);
            const double c1 = Cross(rA0, rA1, rB1);
            const double c2 = Cross(rB0, rB1, rA0);
            const double c3 = Cross(rB0, rB1, rA1);

            if (((c0 > Tolerance && c1 < -Tolerance) ||
                 (c0 < -Tolerance && c1 > Tolerance)) &&
                ((c2 > Tolerance && c3 < -Tolerance) ||
                 (c2 < -Tolerance && c3 > Tolerance))) {
                return true;
            }

            return IsPointOnSegment(rA0, rA1, rB0, Tolerance) ||
                   IsPointOnSegment(rA0, rA1, rB1, Tolerance) ||
                   IsPointOnSegment(rB0, rB1, rA0, Tolerance) ||
                   IsPointOnSegment(rB0, rB1, rA1, Tolerance);
        }

        bool IsPointInTriangle(
            const Point2D& rA,
            const Point2D& rB,
            const Point2D& rC,
            const Point2D& rPoint,
            const double Tolerance) const
        {
            const double c0 = Cross(rA, rB, rPoint);
            const double c1 = Cross(rB, rC, rPoint);
            const double c2 = Cross(rC, rA, rPoint);

            return c0 >= -Tolerance &&
                   c1 >= -Tolerance &&
                   c2 >= -Tolerance;
        }
    };

    const double minimum_span_size = std::min({
        rGridInfo.SpanSizeX,
        rGridInfo.SpanSizeY,
        rGridInfo.SpanSizeZ});

    const double determinant_tolerance =
        6.0e-12 *
        minimum_span_size *
        minimum_span_size *
        minimum_span_size;

    const double projection_tolerance =
        1.0e-12 * std::max(1.0, minimum_span_size);
    const double projection_area_tolerance =
        projection_tolerance * projection_tolerance;
    const PolygonUtilities polygon_utilities;

    std::set<Type3VolumeKey> created_volume_keys;
    std::map<IndexType, NodePointerType> surrogate_nodes_by_id;
    std::unordered_set<
        CanonicalFaceKey3D,
        CanonicalFaceKey3DHasher> initial_type2_open_face_keys;

    for (const auto& r_open_face_data : rType2CreationResult.OpenFaceDataList) {
        KRATOS_ERROR_IF_NOT(r_open_face_data.pSurrogateNode)
            << "[CreateType3GapGeometries] Null surrogate node in type2 open face data.\n";
        KRATOS_ERROR_IF_NOT(r_open_face_data.pProjectionNode0)
            << "[CreateType3GapGeometries] Null first projection node in type2 open face data.\n";
        KRATOS_ERROR_IF_NOT(r_open_face_data.pProjectionNode1)
            << "[CreateType3GapGeometries] Null second projection node in type2 open face data.\n";

        surrogate_nodes_by_id.emplace(
            r_open_face_data.pSurrogateNode->Id(),
            r_open_face_data.pSurrogateNode);
        initial_type2_open_face_keys.insert(MakeCanonicalFaceKey3D(
            r_open_face_data.pSurrogateNode->Id(),
            r_open_face_data.pProjectionNode0->Id(),
            r_open_face_data.pProjectionNode1->Id()));
    }

    for (const auto& r_surrogate_node_entry : surrogate_nodes_by_id) {
        const IndexType surrogate_node_id = r_surrogate_node_entry.first;
        const NodePointerType p_surrogate_node = r_surrogate_node_entry.second;

        std::vector<LocalOpenFace> open_faces;
        for (const auto& r_registry_entry : mLateralFaceRegistry) {
            const auto& r_face_key = r_registry_entry.first;
            const auto& r_occurrences = r_registry_entry.second;
            if (initial_type2_open_face_keys.find(r_face_key) ==
                initial_type2_open_face_keys.end()) {
                continue;
            }

            if (r_occurrences.size() != 1) {
                continue;
            }

            if (r_occurrences.front().GapType != 2) {
                continue;
            }

            const std::size_t number_of_surrogate_nodes =
                static_cast<std::size_t>(r_face_key.NodeIds[0] == surrogate_node_id) +
                static_cast<std::size_t>(r_face_key.NodeIds[1] == surrogate_node_id) +
                static_cast<std::size_t>(r_face_key.NodeIds[2] == surrogate_node_id);
            if (number_of_surrogate_nodes != 1) {
                continue;
            }

            std::array<IndexType, 2> projection_node_ids = {{0, 0}};
            std::size_t next_projection_node_index = 0;
            for (const auto node_id : r_face_key.NodeIds) {
                if (node_id != surrogate_node_id) {
                    projection_node_ids[next_projection_node_index++] = node_id;
                }
            }

            KRATOS_ERROR_IF_NOT(rSkinSubModelPart.HasNode(projection_node_ids[0]))
                << "[CreateType3GapGeometries] Open registry face connected "
                << "to surrogate node references a projection node that does "
                << "not exist in the skin model part.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  missing projection node: "
                << projection_node_ids[0] << "\n"
                << "  face node ids: "
                << r_face_key.NodeIds[0] << ", "
                << r_face_key.NodeIds[1] << ", "
                << r_face_key.NodeIds[2] << "\n";
            KRATOS_ERROR_IF_NOT(rSkinSubModelPart.HasNode(projection_node_ids[1]))
                << "[CreateType3GapGeometries] Open registry face connected "
                << "to surrogate node references a projection node that does "
                << "not exist in the skin model part.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  missing projection node: "
                << projection_node_ids[1] << "\n"
                << "  face node ids: "
                << r_face_key.NodeIds[0] << ", "
                << r_face_key.NodeIds[1] << ", "
                << r_face_key.NodeIds[2] << "\n";

            LocalOpenFace local_open_face;
            local_open_face.pSurrogateNode = p_surrogate_node;
            local_open_face.pProjectionNode0 =
                rSkinSubModelPart.pGetNode(projection_node_ids[0]);
            local_open_face.pProjectionNode1 =
                rSkinSubModelPart.pGetNode(projection_node_ids[1]);
            local_open_face.RegistryKey = r_face_key;
            open_faces.push_back(std::move(local_open_face));
        }

        bool has_incident_type3_span = false;
        const GridPointKey3D surrogate_grid_node =
            ComputeGridPointKey(
                p_surrogate_node->Coordinates(),
                rGridInfo);

        for (int di = -1; di <= 0 && !has_incident_type3_span; ++di) {
            for (int dj = -1; dj <= 0 && !has_incident_type3_span; ++dj) {
                for (int dk = -1; dk <= 0; ++dk) {
                    const SpanKey3D incident_span{
                        surrogate_grid_node.I + di,
                        surrogate_grid_node.J + dj,
                        surrogate_grid_node.K + dk};

                    if (!IsValidSpan(incident_span, rGridInfo)) {
                        continue;
                    }

                    const auto external_span_it =
                        rExternalSpans.find(incident_span);
                    if (external_span_it == rExternalSpans.end()) {
                        continue;
                    }

                    has_incident_type3_span =
                        external_span_it->second.Type == GapSpanType::Type3 &&
                        external_span_it->second.HasProjectionNode();
                    if (has_incident_type3_span) {
                        break;
                    }
                }
            }
        }

        if (open_faces.empty() && !has_incident_type3_span) {
            continue;
        }

        if (has_incident_type3_span) {
            // TODO: seed the local type3 closure from the incident external
            // type3 span. For now this branch is skipped for debugging.
            result.Summary.NumberOfSkippedFaces += open_faces.size();

            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "[CreateType3GapGeometries] Skipping surrogate node with "
                << "incident type3 span. TODO seed logic.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  open faces: " << open_faces.size() << "\n";

            continue;
        }

        std::map<IndexType, NodePointerType> cycle_nodes_by_id;
        std::map<IndexType, std::vector<IndexType>> adjacency;
        for (const auto& r_open_face : open_faces) {
            const IndexType node_id_0 = r_open_face.pProjectionNode0->Id();
            const IndexType node_id_1 = r_open_face.pProjectionNode1->Id();
            cycle_nodes_by_id.emplace(node_id_0, r_open_face.pProjectionNode0);
            cycle_nodes_by_id.emplace(node_id_1, r_open_face.pProjectionNode1);
            adjacency[node_id_0].push_back(node_id_1);
            adjacency[node_id_1].push_back(node_id_0);
        }

        std::ostringstream graph_dump;
        graph_dump << "  surrogate node: " << surrogate_node_id << "\n";
        graph_dump << "  open faces:\n";
        for (const auto& r_open_face : open_faces) {
            graph_dump
                << "    V-" << r_open_face.pProjectionNode0->Id()
                << "-" << r_open_face.pProjectionNode1->Id()
                << " registry face ids "
                << r_open_face.RegistryKey.NodeIds[0] << ", "
                << r_open_face.RegistryKey.NodeIds[1] << ", "
                << r_open_face.RegistryKey.NodeIds[2] << "\n";
        }
        graph_dump << "  graph nodes:\n";
        for (const auto& r_adjacency_entry : adjacency) {
            graph_dump << "    " << r_adjacency_entry.first << " degree "
                       << r_adjacency_entry.second.size() << " neighbours";
            for (const auto neighbour_id : r_adjacency_entry.second) {
                graph_dump << " " << neighbour_id;
            }
            graph_dump << "\n";
        }

        for (const auto& r_adjacency_entry : adjacency) {
            KRATOS_ERROR_IF(r_adjacency_entry.second.size() != 2)
                << "[CreateType3GapGeometries] Open projection-node graph is "
                << "not a simple cycle: every node must have degree two.\n"
                << graph_dump.str();
        }

        std::set<IndexType> visited_projection_node_ids;
        std::vector<std::vector<IndexType>> connected_components;
        for (const auto& r_node_entry : cycle_nodes_by_id) {
            const IndexType seed_node_id = r_node_entry.first;
            if (visited_projection_node_ids.find(seed_node_id) !=
                visited_projection_node_ids.end()) {
                continue;
            }

            std::vector<IndexType> component_node_ids;
            std::vector<IndexType> pending_node_ids;
            pending_node_ids.push_back(seed_node_id);
            visited_projection_node_ids.insert(seed_node_id);

            while (!pending_node_ids.empty()) {
                const IndexType node_id = pending_node_ids.back();
                pending_node_ids.pop_back();
                component_node_ids.push_back(node_id);

                for (const auto neighbour_id : adjacency[node_id]) {
                    if (visited_projection_node_ids.insert(neighbour_id).second) {
                        pending_node_ids.push_back(neighbour_id);
                    }
                }
            }

            connected_components.push_back(std::move(component_node_ids));
        }

        for (const auto& r_component_node_ids : connected_components) {
            std::set<IndexType> component_node_id_set(
                r_component_node_ids.begin(),
                r_component_node_ids.end());
            std::size_t number_of_component_edges = 0;
            for (const auto node_id : r_component_node_ids) {
                for (const auto neighbour_id : adjacency[node_id]) {
                    if (node_id < neighbour_id &&
                        component_node_id_set.find(neighbour_id) !=
                            component_node_id_set.end()) {
                        ++number_of_component_edges;
                    }
                }
            }

            KRATOS_ERROR_IF(r_component_node_ids.size() < 3)
                << "[CreateType3GapGeometries] Open projection-node graph "
                << "component has fewer than three nodes.\n"
                << graph_dump.str();

            if (number_of_component_edges != r_component_node_ids.size()) {
                std::ostringstream component_dump;
                component_dump << "  component nodes:";
                for (const auto node_id : r_component_node_ids) {
                    component_dump << " " << node_id;
                }
                component_dump << "\n";

                KRATOS_ERROR
                    << "[CreateType3GapGeometries] Open projection-node graph "
                    << "component is not a simple cycle because number of edges "
                    << "and nodes differ.\n"
                    << component_dump.str()
                    << graph_dump.str();
            }

            CanonicalFaceKey3D component_reference_face_key;
            bool found_component_reference_face = false;
            for (const auto& r_open_face : open_faces) {
                if (component_node_id_set.find(r_open_face.pProjectionNode0->Id()) !=
                        component_node_id_set.end() &&
                    component_node_id_set.find(r_open_face.pProjectionNode1->Id()) !=
                        component_node_id_set.end()) {
                    component_reference_face_key = r_open_face.RegistryKey;
                    found_component_reference_face = true;
                    break;
                }
            }
            KRATOS_ERROR_IF_NOT(found_component_reference_face)
                << "[CreateType3GapGeometries] Could not find a registry "
                << "reference face for type3 cycle component.\n"
                << graph_dump.str();

        std::vector<IndexType> ordered_projection_node_ids;
        ordered_projection_node_ids.reserve(r_component_node_ids.size());
        IndexType previous_node_id = 0;
        IndexType current_node_id = r_component_node_ids.front();

        while (ordered_projection_node_ids.size() < r_component_node_ids.size()) {
            ordered_projection_node_ids.push_back(current_node_id);
            const auto& r_neighbours = adjacency[current_node_id];
            const IndexType next_node_id =
                r_neighbours[0] == previous_node_id ? r_neighbours[1] : r_neighbours[0];
            previous_node_id = current_node_id;
            current_node_id = next_node_id;

            KRATOS_ERROR_IF(
                std::find(
                    ordered_projection_node_ids.begin(),
                    ordered_projection_node_ids.end(),
                    current_node_id) != ordered_projection_node_ids.end() &&
                ordered_projection_node_ids.size() < r_component_node_ids.size())
                << "[CreateType3GapGeometries] Open projection-node graph is "
                << "not separable into simple cycles.\n"
                << graph_dump.str();
        }

        KRATOS_ERROR_IF(current_node_id != ordered_projection_node_ids.front())
            << "[CreateType3GapGeometries] Open projection-node graph "
            << "component is not closed by a simple cycle.\n"
            << graph_dump.str();

        array_1d<double, 3> centroid = ZeroVector(3);
        for (const auto node_id : ordered_projection_node_ids) {
            noalias(centroid) += cycle_nodes_by_id[node_id]->Coordinates();
        }
        centroid /= static_cast<double>(ordered_projection_node_ids.size());

        array_1d<double, 3> normal = ZeroVector(3);
        for (std::size_t i = 0; i < ordered_projection_node_ids.size(); ++i) {
            const auto& r_point_i =
                cycle_nodes_by_id[ordered_projection_node_ids[i]]->Coordinates();
            const auto& r_point_j =
                cycle_nodes_by_id[
                    ordered_projection_node_ids[
                        (i + 1) % ordered_projection_node_ids.size()]]->Coordinates();
            normal[0] += (r_point_i[1] - r_point_j[1]) * (r_point_i[2] + r_point_j[2]);
            normal[1] += (r_point_i[2] - r_point_j[2]) * (r_point_i[0] + r_point_j[0]);
            normal[2] += (r_point_i[0] - r_point_j[0]) * (r_point_i[1] + r_point_j[1]);
        }

        const double normal_norm = norm_2(normal);
        KRATOS_ERROR_IF(normal_norm <= projection_tolerance)
            << "[CreateType3GapGeometries] Could not define a local projection "
            << "plane for type3 cycle.\n"
            << graph_dump.str();
        normal /= normal_norm;

        array_1d<double, 3> axis_x =
            cycle_nodes_by_id[ordered_projection_node_ids.front()]->Coordinates() -
            centroid;
        const double axis_x_norm = norm_2(axis_x);
        KRATOS_ERROR_IF(axis_x_norm <= projection_tolerance)
            << "[CreateType3GapGeometries] Could not define local projection "
            << "axis for type3 cycle.\n"
            << graph_dump.str();
        axis_x /= axis_x_norm;
        array_1d<double, 3> axis_y =
            MathUtils<double>::CrossProduct(normal, axis_x);

        std::vector<Point2D> projected_points;
        projected_points.reserve(ordered_projection_node_ids.size());
        for (const auto node_id : ordered_projection_node_ids) {
            const array_1d<double, 3> point_vector =
                cycle_nodes_by_id[node_id]->Coordinates() - centroid;
            Point2D point;
            point.X = inner_prod(point_vector, axis_x);
            point.Y = inner_prod(point_vector, axis_y);
            projected_points.push_back(point);
        }

        double signed_area = 0.0;
        for (std::size_t i = 0; i < projected_points.size(); ++i) {
            const auto& r_point_i = projected_points[i];
            const auto& r_point_j =
                projected_points[(i + 1) % projected_points.size()];
            signed_area += r_point_i.X * r_point_j.Y -
                           r_point_j.X * r_point_i.Y;
        }
        signed_area *= 0.5;

        KRATOS_ERROR_IF(std::abs(signed_area) <= projection_area_tolerance)
            << "[CreateType3GapGeometries] Projected type3 cycle has near-zero "
            << "area.\n"
            << graph_dump.str();

        if (signed_area < 0.0) {
            std::reverse(
                ordered_projection_node_ids.begin(),
                ordered_projection_node_ids.end());
            std::reverse(projected_points.begin(), projected_points.end());
        }

        for (std::size_t i = 0; i < projected_points.size(); ++i) {
            const std::size_t i_next = (i + 1) % projected_points.size();
            for (std::size_t j = i + 1; j < projected_points.size(); ++j) {
                const std::size_t j_next = (j + 1) % projected_points.size();
                if (i == j_next || i_next == j) {
                    continue;
                }
                if (i == 0 && j_next == 0) {
                    continue;
                }

                KRATOS_ERROR_IF(polygon_utilities.DoSegmentsIntersect(
                    projected_points[i],
                    projected_points[i_next],
                    projected_points[j],
                    projected_points[j_next],
                    projection_tolerance))
                    << "[CreateType3GapGeometries] Projected type3 cycle "
                    << "self-intersects.\n"
                    << graph_dump.str();
            }
        }

        const auto make_polygon_edge_key =
            [](const std::size_t Index0, const std::size_t Index1) {
                std::array<std::size_t, 2> edge_key{{Index0, Index1}};
                std::sort(edge_key.begin(), edge_key.end());
                return edge_key;
            };

        std::set<std::array<std::size_t, 2>> original_boundary_edges;
        for (std::size_t i = 0; i < projected_points.size(); ++i) {
            original_boundary_edges.insert(make_polygon_edge_key(
                i,
                (i + 1) % projected_points.size()));
        }

        std::unordered_map<
            CanonicalFaceKey3D,
            std::size_t,
            CanonicalFaceKey3DHasher> planned_face_occurrences;
        std::set<Type3VolumeKey> planned_volume_keys;
        std::vector<PlannedTriangleCandidate> planned_candidates;
        std::vector<LocalTriangle> planned_triangles;
        std::ostringstream ear_reject_dump;

        const auto registry_occurrences =
            [&](const CanonicalFaceKey3D& rFaceKey) -> std::size_t {
                const auto registry_it = mLateralFaceRegistry.find(rFaceKey);
                return registry_it == mLateralFaceRegistry.end() ?
                    0 : registry_it->second.size();
            };

        const auto planned_occurrences =
            [&](const CanonicalFaceKey3D& rFaceKey) -> std::size_t {
                const auto planned_it = planned_face_occurrences.find(rFaceKey);
                return planned_it == planned_face_occurrences.end() ?
                    0 : planned_it->second;
            };

        const auto append_face_occurrence_dump =
            [&](const std::array<CanonicalFaceKey3D, 4>& rFaceKeys) {
                for (const auto& r_face_key : rFaceKeys) {
                    ear_reject_dump
                        << "    face ids "
                        << r_face_key.NodeIds[0] << " "
                        << r_face_key.NodeIds[1] << " "
                        << r_face_key.NodeIds[2]
                        << " registry_occ="
                        << registry_occurrences(r_face_key)
                        << " planned_occ="
                        << planned_occurrences(r_face_key)
                        << "\n";
                }
            };

        const auto append_ear_reject =
            [&](const LocalTriangle& rTriangle,
                const std::string& rReason) {
                ear_reject_dump
                    << "  EAR_REJECT "
                    << ordered_projection_node_ids[rTriangle.Index0] << " "
                    << ordered_projection_node_ids[rTriangle.Index1] << " "
                    << ordered_projection_node_ids[rTriangle.Index2]
                    << " reason=" << rReason << "\n";
            };

        const auto make_candidate_face_keys =
            [&](const LocalTriangle& rTriangle) {
                const auto node_id_0 =
                    ordered_projection_node_ids[rTriangle.Index0];
                const auto node_id_1 =
                    ordered_projection_node_ids[rTriangle.Index1];
                const auto node_id_2 =
                    ordered_projection_node_ids[rTriangle.Index2];

                return std::array<CanonicalFaceKey3D, 4>{{
                    MakeCanonicalFaceKey3D(
                        surrogate_node_id,
                        node_id_0,
                        node_id_1),
                    MakeCanonicalFaceKey3D(
                        surrogate_node_id,
                        node_id_1,
                        node_id_2),
                    MakeCanonicalFaceKey3D(
                        surrogate_node_id,
                        node_id_2,
                        node_id_0),
                    MakeCanonicalFaceKey3D(
                        node_id_0,
                        node_id_1,
                        node_id_2)}};
            };

        const auto append_overlap_dump =
            [&](const LocalTriangle& rCandidateTriangle,
                const LocalTriangle& rPlannedTriangle) {
                const std::array<IndexType, 3> candidate_node_ids{{
                    ordered_projection_node_ids[rCandidateTriangle.Index0],
                    ordered_projection_node_ids[rCandidateTriangle.Index1],
                    ordered_projection_node_ids[rCandidateTriangle.Index2]}};
                const std::array<IndexType, 3> planned_node_ids{{
                    ordered_projection_node_ids[rPlannedTriangle.Index0],
                    ordered_projection_node_ids[rPlannedTriangle.Index1],
                    ordered_projection_node_ids[rPlannedTriangle.Index2]}};

                std::size_t shared_node_count = 1; // The surrogate node V.
                for (const auto candidate_node_id : candidate_node_ids) {
                    if (std::find(
                            planned_node_ids.begin(),
                            planned_node_ids.end(),
                            candidate_node_id) != planned_node_ids.end()) {
                        ++shared_node_count;
                    }
                }

                ear_reject_dump
                    << "    overlap_with_tetra "
                    << surrogate_node_id << " "
                    << planned_node_ids[0] << " "
                    << planned_node_ids[1] << " "
                    << planned_node_ids[2] << "\n"
                    << "    shared_node_count=" << shared_node_count << "\n";
            };

        const auto is_point_strictly_in_triangle =
            [&](const Point2D& rPoint, const LocalTriangle& rTriangle) {
                return polygon_utilities.Cross(
                           projected_points[rTriangle.Index0],
                           projected_points[rTriangle.Index1],
                           rPoint) > projection_tolerance &&
                       polygon_utilities.Cross(
                           projected_points[rTriangle.Index1],
                           projected_points[rTriangle.Index2],
                           rPoint) > projection_tolerance &&
                       polygon_utilities.Cross(
                           projected_points[rTriangle.Index2],
                           projected_points[rTriangle.Index0],
                           rPoint) > projection_tolerance;
            };

        const auto do_segments_properly_intersect =
            [&](const Point2D& rA0,
                const Point2D& rA1,
                const Point2D& rB0,
                const Point2D& rB1) {
                const double c0 = polygon_utilities.Cross(rA0, rA1, rB0);
                const double c1 = polygon_utilities.Cross(rA0, rA1, rB1);
                const double c2 = polygon_utilities.Cross(rB0, rB1, rA0);
                const double c3 = polygon_utilities.Cross(rB0, rB1, rA1);

                return ((c0 > projection_tolerance && c1 < -projection_tolerance) ||
                        (c0 < -projection_tolerance && c1 > projection_tolerance)) &&
                       ((c2 > projection_tolerance && c3 < -projection_tolerance) ||
                        (c2 < -projection_tolerance && c3 > projection_tolerance));
            };

        const auto triangles_have_area_overlap =
            [&](const LocalTriangle& rCandidateTriangle,
                const LocalTriangle& rPlannedTriangle) {
                const std::array<std::size_t, 3> candidate_indices{{
                    rCandidateTriangle.Index0,
                    rCandidateTriangle.Index1,
                    rCandidateTriangle.Index2}};
                const std::array<std::size_t, 3> planned_indices{{
                    rPlannedTriangle.Index0,
                    rPlannedTriangle.Index1,
                    rPlannedTriangle.Index2}};

                for (std::size_t i = 0; i < 3; ++i) {
                    const auto candidate_i0 = candidate_indices[i];
                    const auto candidate_i1 = candidate_indices[(i + 1) % 3];
                    for (std::size_t j = 0; j < 3; ++j) {
                        const auto planned_i0 = planned_indices[j];
                        const auto planned_i1 = planned_indices[(j + 1) % 3];
                        if (do_segments_properly_intersect(
                            projected_points[candidate_i0],
                            projected_points[candidate_i1],
                            projected_points[planned_i0],
                            projected_points[planned_i1])) {
                            return true;
                        }
                    }
                }

                for (const auto candidate_index : candidate_indices) {
                    if (is_point_strictly_in_triangle(
                        projected_points[candidate_index],
                        rPlannedTriangle)) {
                        return true;
                    }
                }

                for (const auto planned_index : planned_indices) {
                    if (is_point_strictly_in_triangle(
                        projected_points[planned_index],
                        rCandidateTriangle)) {
                        return true;
                    }
                }

                return false;
            };

        const auto build_candidate =
            [&](const LocalTriangle& rTriangle,
                const bool IsNewDiagonal,
                PlannedTriangleCandidate& rCandidate) {
                const auto node_id_0 =
                    ordered_projection_node_ids[rTriangle.Index0];
                const auto node_id_1 =
                    ordered_projection_node_ids[rTriangle.Index1];
                const auto node_id_2 =
                    ordered_projection_node_ids[rTriangle.Index2];

                const double projected_area_cross = polygon_utilities.Cross(
                    projected_points[rTriangle.Index0],
                    projected_points[rTriangle.Index1],
                    projected_points[rTriangle.Index2]);

                if (projected_area_cross <= projection_area_tolerance) {
                    append_ear_reject(
                        rTriangle,
                        "non_positive_projected_area cross=" +
                        std::to_string(projected_area_cross));
                    append_face_occurrence_dump(
                        make_candidate_face_keys(rTriangle));
                    return false;
                }

                const NodePointerType p_node_0 = cycle_nodes_by_id[node_id_0];
                const NodePointerType p_node_1 = cycle_nodes_by_id[node_id_1];
                const NodePointerType p_node_2 = cycle_nodes_by_id[node_id_2];

                const array_1d<double, 3> edge_0 =
                    p_node_0->Coordinates() -
                    p_surrogate_node->Coordinates();
                const array_1d<double, 3> edge_1 =
                    p_node_1->Coordinates() -
                    p_surrogate_node->Coordinates();
                const array_1d<double, 3> edge_2 =
                    p_node_2->Coordinates() -
                    p_surrogate_node->Coordinates();
                const double determinant = inner_prod(
                    edge_0,
                    MathUtils<double>::CrossProduct(edge_1, edge_2));

                if (determinant <= determinant_tolerance) {
                    append_ear_reject(
                        rTriangle,
                        "non_positive_determinant det=" +
                        std::to_string(determinant));
                    append_face_occurrence_dump(
                        make_candidate_face_keys(rTriangle));
                    return false;
                }

                rCandidate.Triangle = rTriangle;
                rCandidate.FaceKeys = make_candidate_face_keys(rTriangle);
                rCandidate.VolumeKey = Type3VolumeKey::Create(
                    p_surrogate_node,
                    p_node_0,
                    p_node_1,
                    p_node_2);

                if (created_volume_keys.find(rCandidate.VolumeKey) !=
                        created_volume_keys.end() ||
                    planned_volume_keys.find(rCandidate.VolumeKey) !=
                        planned_volume_keys.end()) {
                    append_ear_reject(rTriangle, "duplicated_tetra");
                    append_face_occurrence_dump(rCandidate.FaceKeys);
                    return false;
                }

                for (const auto& r_face_key : rCandidate.FaceKeys) {
                    if (registry_occurrences(r_face_key) +
                        planned_occurrences(r_face_key) + 1 > 2) {
                        append_ear_reject(
                            rTriangle,
                            "face_occurrence_limit");
                        append_face_occurrence_dump(rCandidate.FaceKeys);
                        return false;
                    }
                }

                if (IsNewDiagonal) {
                    const auto& r_new_diagonal_face_key =
                        rCandidate.FaceKeys[2];
                    if (registry_occurrences(r_new_diagonal_face_key) != 0 ||
                        planned_occurrences(r_new_diagonal_face_key) != 0) {
                        append_ear_reject(
                            rTriangle,
                            "new_diagonal_face_already_exists");
                        append_face_occurrence_dump(rCandidate.FaceKeys);
                        return false;
                    }
                }

                for (const auto& r_planned_triangle : planned_triangles) {
                    if (triangles_have_area_overlap(rTriangle, r_planned_triangle)) {
                        append_ear_reject(
                            rTriangle,
                            "overlap_with_planned_tetra");
                        append_face_occurrence_dump(rCandidate.FaceKeys);
                        append_overlap_dump(rTriangle, r_planned_triangle);
                        return false;
                    }
                }

                return true;
            };

        const auto apply_candidate =
            [&](const PlannedTriangleCandidate& rCandidate) {
                planned_candidates.push_back(rCandidate);
                planned_triangles.push_back(rCandidate.Triangle);
                planned_volume_keys.insert(rCandidate.VolumeKey);
                for (const auto& r_face_key : rCandidate.FaceKeys) {
                    ++planned_face_occurrences[r_face_key];
                }
            };

        const auto undo_candidate =
            [&](const PlannedTriangleCandidate& rCandidate) {
                planned_candidates.pop_back();
                planned_triangles.pop_back();
                planned_volume_keys.erase(rCandidate.VolumeKey);
                for (const auto& r_face_key : rCandidate.FaceKeys) {
                    auto planned_it =
                        planned_face_occurrences.find(r_face_key);
                    KRATOS_ERROR_IF(
                        planned_it == planned_face_occurrences.end() ||
                        planned_it->second == 0)
                        << "[CreateType3GapGeometries] Invalid planned face "
                        << "occurrence rollback.\n";

                    --planned_it->second;
                    if (planned_it->second == 0) {
                        planned_face_occurrences.erase(planned_it);
                    }
                }
            };

        std::function<bool(const std::vector<std::size_t>&)>
            plan_triangulation =
                [&](const std::vector<std::size_t>& rPolygonIndices) {
                    if (rPolygonIndices.size() < 3) {
                        return false;
                    }

                    for (std::size_t i = 0; i < rPolygonIndices.size(); ++i) {
                        const std::size_t previous_index =
                            rPolygonIndices[
                                (i + rPolygonIndices.size() - 1) %
                                rPolygonIndices.size()];
                        const std::size_t current_index = rPolygonIndices[i];
                        const std::size_t next_index =
                            rPolygonIndices[
                                (i + 1) % rPolygonIndices.size()];

                        const LocalTriangle candidate_triangle{
                            previous_index,
                            current_index,
                            next_index};

                        bool contains_other_point = false;
                        IndexType contained_node_id = 0;
                        for (const auto polygon_index : rPolygonIndices) {
                            if (polygon_index == previous_index ||
                                polygon_index == current_index ||
                                polygon_index == next_index) {
                                continue;
                            }

                            if (polygon_utilities.IsPointInTriangle(
                                projected_points[previous_index],
                                projected_points[current_index],
                                projected_points[next_index],
                                projected_points[polygon_index],
                                projection_tolerance)) {
                                contains_other_point = true;
                                contained_node_id =
                                    ordered_projection_node_ids[polygon_index];
                                break;
                            }
                        }

                        if (contains_other_point) {
                            append_ear_reject(
                                candidate_triangle,
                                "contains_projection_point point_id=" +
                                std::to_string(contained_node_id));
                            append_face_occurrence_dump(
                                make_candidate_face_keys(candidate_triangle));
                            continue;
                        }

                        const bool is_new_diagonal =
                            rPolygonIndices.size() > 3 &&
                            original_boundary_edges.find(make_polygon_edge_key(
                                previous_index,
                                next_index)) == original_boundary_edges.end();

                        PlannedTriangleCandidate candidate;
                        if (!build_candidate(
                            candidate_triangle,
                            is_new_diagonal,
                            candidate)) {
                            continue;
                        }

                        apply_candidate(candidate);

                        if (rPolygonIndices.size() == 3) {
                            return true;
                        }

                        std::vector<std::size_t> remaining_polygon_indices =
                            rPolygonIndices;
                        remaining_polygon_indices.erase(
                            remaining_polygon_indices.begin() + i);

                        if (plan_triangulation(remaining_polygon_indices)) {
                            return true;
                        }

                        undo_candidate(candidate);
                    }

                    return false;
                };

        std::vector<std::size_t> polygon_indices(projected_points.size());
        for (std::size_t i = 0; i < polygon_indices.size(); ++i) {
            polygon_indices[i] = i;
        }

        KRATOS_ERROR_IF_NOT(plan_triangulation(polygon_indices))
            << "[CreateType3GapGeometries] Could not find a valid type3 "
            << "triangulation for the projected cycle without exceeding two "
            << "face occurrences or overlapping planned tetrahedra.\n"
            << graph_dump.str()
            << "  ear rejections:\n"
            << ear_reject_dump.str();

        std::vector<LocalTriangle> local_triangles;
        local_triangles.reserve(planned_candidates.size());
        for (const auto& r_candidate : planned_candidates) {
            local_triangles.push_back(r_candidate.Triangle);
        }

        for (const auto& r_triangle : local_triangles) {
            NodePointerType p_projection_node_0 =
                cycle_nodes_by_id[
                    ordered_projection_node_ids[r_triangle.Index0]];
            NodePointerType p_projection_node_1 =
                cycle_nodes_by_id[
                    ordered_projection_node_ids[r_triangle.Index1]];
            NodePointerType p_projection_node_2 =
                cycle_nodes_by_id[
                    ordered_projection_node_ids[r_triangle.Index2]];

            const array_1d<double, 3> edge_0 =
                p_projection_node_0->Coordinates() -
                p_surrogate_node->Coordinates();
            const array_1d<double, 3> edge_1 =
                p_projection_node_1->Coordinates() -
                p_surrogate_node->Coordinates();
            const array_1d<double, 3> edge_2 =
                p_projection_node_2->Coordinates() -
                p_surrogate_node->Coordinates();
            double determinant = inner_prod(
                edge_0,
                MathUtils<double>::CrossProduct(edge_1, edge_2));

            KRATOS_ERROR_IF(std::abs(determinant) <= determinant_tolerance)
                << "[CreateType3GapGeometries] Degenerate type3 tetra.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  projection nodes: "
                << p_projection_node_0->Id() << ", "
                << p_projection_node_1->Id() << ", "
                << p_projection_node_2->Id() << "\n"
                << "  determinant: " << determinant << "\n";

            if (determinant < 0.0) {
                std::swap(p_projection_node_1, p_projection_node_2);
                determinant = -determinant;
            }

            const Type3VolumeKey volume_key = Type3VolumeKey::Create(
                p_surrogate_node,
                p_projection_node_0,
                p_projection_node_1,
                p_projection_node_2);

            KRATOS_ERROR_IF(
                created_volume_keys.find(volume_key) !=
                created_volume_keys.end())
                << "[CreateType3GapGeometries] Duplicated type3 tetra.\n"
                << "  surrogate node: " << surrogate_node_id << "\n"
                << "  projection nodes: "
                << p_projection_node_0->Id() << ", "
                << p_projection_node_1->Id() << ", "
                << p_projection_node_2->Id() << "\n";

            created_volume_keys.insert(volume_key);

            const auto first_open_face_it =
                mLateralFaceRegistry.find(component_reference_face_key);
            KRATOS_ERROR_IF(
                first_open_face_it == mLateralFaceRegistry.end() ||
                first_open_face_it->second.empty())
                << "[CreateType3GapGeometries] Could not retrieve a reference "
                << "open face occurrence from the registry.\n"
                << graph_dump.str();
            const LateralFaceOccurrence reference_occurrence =
                first_open_face_it->second.front();
            Geometry<Node>::Pointer p_type3_neighbour_geometry =
                reference_occurrence.pNeighbourGeometry;

            KRATOS_ERROR_IF_NOT(p_type3_neighbour_geometry)
                << "[CreateType3GapGeometries] Null neighbour geometry in "
                << "reference open face.\n";

            auto p_type3_volume = CreateType3CollapsedCornerVolume(
                p_surrogate_node,
                p_projection_node_0,
                p_projection_node_1,
                p_projection_node_2);

            GeometriesArrayType volume_quadrature_point_geometries =
                CreateAndTagVolumeQuadraturePointGeometries(
                    p_type3_volume,
                    p_type3_neighbour_geometry,
                    IntegrationOrder,
                    NumberOfShapeFunctionsDerivatives);

            KRATOS_ERROR_IF(volume_quadrature_point_geometries.size() == 0)
                << "[CreateType3GapGeometries] Type3 tetra generated no "
                << "quadrature point geometries.\n";

            CheckType3QuadraturePointGeometries(
                volume_quadrature_point_geometries,
                p_surrogate_node,
                p_projection_node_0,
                p_projection_node_1,
                p_projection_node_2);

            p_type3_volume->SetId(next_geometry_id++);
            AddNeighbourGeometry(
                *p_type3_volume,
                p_type3_neighbour_geometry);
            if (mStoreGapDebugGeometries) {
                r_gap_type3_debug.AddGeometry(p_type3_volume);
            }

            const std::array<LocalSideFace, 3> lateral_faces = {{
                {p_projection_node_0, p_projection_node_1, p_projection_node_2},
                {p_projection_node_1, p_projection_node_2, p_projection_node_0},
                {p_projection_node_2, p_projection_node_0, p_projection_node_1}}};

            for (const auto& r_lateral_face : lateral_faces) {
                const auto canonical_face_key = MakeCanonicalFaceKey3D(
                    p_surrogate_node->Id(),
                    r_lateral_face.pProjectionNode0->Id(),
                    r_lateral_face.pProjectionNode1->Id());
                const auto registry_it =
                    mLateralFaceRegistry.find(canonical_face_key);

                LateralFaceOccurrence registration_occurrence =
                    reference_occurrence;
                const bool face_already_exists =
                    registry_it != mLateralFaceRegistry.end();

                if (face_already_exists) {
                    KRATOS_ERROR_IF(registry_it->second.size() >= 2)
                        << "[CreateType3GapGeometries] Trying to add a third "
                        << "occurrence to a lateral face.\n"
                        << "  face node ids: "
                        << canonical_face_key.NodeIds[0] << ", "
                        << canonical_face_key.NodeIds[1] << ", "
                        << canonical_face_key.NodeIds[2] << "\n"
                        << "  occurrences: "
                        << registry_it->second.size() << "\n";

                    registration_occurrence =
                        registry_it->second.front();
                }

                RegisterType3LateralFace(
                    r_gap_type3_debug,
                    next_geometry_id,
                    registration_occurrence.SurrogateConditionId,
                    registration_occurrence.ExternalSpan,
                    p_surrogate_node,
                    r_lateral_face.pProjectionNode0,
                    r_lateral_face.pProjectionNode1,
                    r_lateral_face.pOppositeProjectionNode,
                    p_type3_neighbour_geometry);
            }

            const auto top_face_key = MakeCanonicalFaceKey3D(
                p_projection_node_0->Id(),
                p_projection_node_1->Id(),
                p_projection_node_2->Id());
            const auto top_registry_it =
                mLateralFaceRegistry.find(top_face_key);

            LateralFaceOccurrence top_registration_occurrence =
                reference_occurrence;
            if (top_registry_it != mLateralFaceRegistry.end()) {
                KRATOS_ERROR_IF(top_registry_it->second.size() >= 2)
                    << "[CreateType3GapGeometries] Trying to add a third "
                    << "occurrence to a top face.\n"
                    << "  top face node ids: "
                    << top_face_key.NodeIds[0] << ", "
                    << top_face_key.NodeIds[1] << ", "
                    << top_face_key.NodeIds[2] << "\n"
                    << "  occurrences: "
                    << top_registry_it->second.size() << "\n";

                top_registration_occurrence =
                    top_registry_it->second.front();
            }

            RegisterType3LateralFace(
                r_gap_type3_debug,
                next_geometry_id,
                top_registration_occurrence.SurrogateConditionId,
                top_registration_occurrence.ExternalSpan,
                p_projection_node_0,
                p_projection_node_1,
                p_projection_node_2,
                p_surrogate_node,
                p_type3_neighbour_geometry);

            Type3VolumeQuadratureData type3_volume_data;
            type3_volume_data.VolumeQuadraturePointGeometries =
                std::move(volume_quadrature_point_geometries);
            type3_volume_data.NeighbourGeometries.push_back(
                p_type3_neighbour_geometry);

            type3_volume_data.CharacteristicLength = std::max({
                norm_2(
                    p_projection_node_0->Coordinates() -
                    p_surrogate_node->Coordinates()),
                norm_2(
                    p_projection_node_1->Coordinates() -
                    p_surrogate_node->Coordinates()),
                norm_2(
                    p_projection_node_2->Coordinates() -
                    p_surrogate_node->Coordinates())});

            type3_volume_data.SurrogateNodeId =
                p_surrogate_node->Id();
            type3_volume_data.ProjectionNodeId0 =
                p_projection_node_0->Id();
            type3_volume_data.ProjectionNodeId1 =
                p_projection_node_1->Id();
            type3_volume_data.ProjectionNodeId2 =
                p_projection_node_2->Id();

            result.VolumeQuadratureDataList.push_back(
                std::move(type3_volume_data));

            ++result.Summary.NumberOfCreatedVolumes;
        }
        }
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type 3 creation summary:\n"
        << "  candidate OpenFaceDataList faces: "
        << result.Summary.NumberOfCandidateOpenFaces << "\n"
        << "  created volumes: "
        << result.Summary.NumberOfCreatedVolumes << "\n"
        << "  skipped faces: "
        << result.Summary.NumberOfSkippedFaces << "\n"
        << "  GapType3Debug geometries: "
        << r_gap_type3_debug.NumberOfGeometries() << "\n";

    return result;
}

#endif

void SnakeGapSbm3DUtilities::AddNeighbourGeometry(
    Geometry<Node>& rGeometry,
    const Geometry<Node>::Pointer& pNeighbourGeometry) const
{
    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::AddNeighbourGeometry] "
        << "Trying to add null neighbour geometry.\n";

    NeighbourGeometriesVectorType neighbour_geometries;

    if (rGeometry.Has(NEIGHBOUR_GEOMETRIES)) {
        neighbour_geometries = rGeometry.GetValue(NEIGHBOUR_GEOMETRIES);
    }

    neighbour_geometries.push_back(pNeighbourGeometry);

    rGeometry.SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
}

void SnakeGapSbm3DUtilities::AddUniqueNeighbourGeometry(
    Geometry<Node>& rGeometry,
    const Geometry<Node>::Pointer& pCandidateNeighbourGeometry) const
{
    KRATOS_ERROR_IF_NOT(pCandidateNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::AddUniqueNeighbourGeometry] "
        << "Trying to add null neighbour geometry.\n";

    NeighbourGeometriesVectorType neighbour_geometries;

    if (rGeometry.Has(NEIGHBOUR_GEOMETRIES)) {
        neighbour_geometries = rGeometry.GetValue(NEIGHBOUR_GEOMETRIES);
    }

    if (!ContainsNeighbourGeometry(
            neighbour_geometries,
            pCandidateNeighbourGeometry)) {
        neighbour_geometries.push_back(pCandidateNeighbourGeometry);
        rGeometry.SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
    } else if (!rGeometry.Has(NEIGHBOUR_GEOMETRIES)) {
        rGeometry.SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
    }
}

void SnakeGapSbm3DUtilities::SetType1GeometryData(
    Geometry<Node>& rGeometry,
    const IndexType ProjectionNodeId,
    const Geometry<Node>::Pointer& pNeighbourGeometry,
    const std::string& rIdentifier) const
{
    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::SetType1GeometryData] "
        << "NEIGHBOUR_GEOMETRY is null.\n";

    rGeometry.SetValue(IDENTIFIER, rIdentifier);
    rGeometry.SetValue(PROJECTION_NODE_ID, ProjectionNodeId);
    AddNeighbourGeometry(rGeometry, pNeighbourGeometry);
}

void SnakeGapSbm3DUtilities::CreateAndTagSurfaceQuadraturePointGeometries(
    const NurbsSurfaceType::Pointer& pSurfaceGeometry,
    const Geometry<Node>::Pointer& pNeighbourGeometry) const
{
    KRATOS_ERROR_IF_NOT(pSurfaceGeometry)
        << "[SnakeGapSbm3DUtilities::CreateAndTagSurfaceQuadraturePointGeometries] "
        << "Surface geometry is null.\n";

    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::CreateAndTagSurfaceQuadraturePointGeometries] "
        << "Neighbour geometry is null.\n";

    IntegrationPointsArrayType integration_points;

    pSurfaceGeometry->CreateIntegrationPoints(
        integration_points,
        2,
        2);

    IntegrationInfo integration_info = pSurfaceGeometry->GetDefaultIntegrationInfo();

    GeometriesArrayType quadrature_point_geometries;

    pSurfaceGeometry->CreateQuadraturePointGeometries(
        quadrature_point_geometries,
        1,
        integration_points,
        integration_info);

    for (auto& p_quadrature_geometry : quadrature_point_geometries) {
        AddNeighbourGeometry(p_quadrature_geometry, pNeighbourGeometry);
    }
}

SnakeGapSbm3DUtilities::IntegrationPointsArrayType
SnakeGapSbm3DUtilities::CreateCoonsVolumeGaussPoints(
    const std::size_t IntegrationOrder,
    const NurbsVolumeType& rGapVolume) const
{
    KRATOS_ERROR_IF(IntegrationOrder == 0)
        << "[SnakeGapSbm3DUtilities::CreateCoonsVolumeGaussPoints] "
        << "Integration order must be positive.\n";

    IntegrationPointsArrayType integration_points(
        IntegrationOrder * IntegrationOrder * IntegrationOrder);

    auto integration_point_it = integration_points.begin();

    IntegrationPointUtilities::IntegrationPoints3D(
        integration_point_it,
        IntegrationOrder,
        IntegrationOrder,
        IntegrationOrder,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0);

    constexpr double determinant_tolerance = 1.0e-14;

    std::vector<double> determinant_values;
    determinant_values.reserve(integration_points.size());

    std::size_t number_of_positive_det = 0;
    std::size_t number_of_negative_det = 0;

    for (auto& r_integration_point : integration_points) {
        CoordinatesArrayType local_coordinates = ZeroVector(3);
        local_coordinates[0] = r_integration_point[0];
        local_coordinates[1] = r_integration_point[1];
        local_coordinates[2] = r_integration_point[2];

        Matrix jacobian;
        rGapVolume.Jacobian(jacobian, local_coordinates);

        const double det_jacobian = MathUtils<double>::Det(jacobian);
        determinant_values.push_back(det_jacobian);

        if (det_jacobian > determinant_tolerance) {
            ++number_of_positive_det;
        } else if (det_jacobian < -determinant_tolerance) {
            ++number_of_negative_det;
        }
    }

    CoordinatesArrayType center_local_coordinates = ZeroVector(3);
    center_local_coordinates[0] = 0.5;
    center_local_coordinates[1] = 0.5;
    center_local_coordinates[2] = 0.5;

    Matrix center_jacobian;
    rGapVolume.Jacobian(center_jacobian, center_local_coordinates);

    const double center_det_jacobian =
        MathUtils<double>::Det(center_jacobian);

    double orientation_sign = center_det_jacobian < 0.0 ? -1.0 : 1.0;

    if (number_of_negative_det > 0 && number_of_positive_det == 0) {
        orientation_sign = -1.0;
    } else if (number_of_positive_det > 0 && number_of_negative_det == 0) {
        orientation_sign = 1.0;
    }

    for (std::size_t point_index = 0; point_index < integration_points.size(); ++point_index) {
        auto& r_integration_point = integration_points[point_index];

        const double oriented_det_jacobian =
            orientation_sign * determinant_values[point_index];

        if (oriented_det_jacobian < -1.0e-1) {
            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "[CreateCoonsVolumeGaussPoints] Inconsistent Coons volume orientation "
                << "at local point " << r_integration_point.Coordinates()
                << ". Oriented determinant = " << oriented_det_jacobian << ".\n";
        }

        // r_integration_point.SetWeight(
        //     r_integration_point.Weight() * std::max(oriented_det_jacobian, 0.0));

        r_integration_point.SetWeight(
            r_integration_point.Weight() * oriented_det_jacobian); //FIXME:
    }

    return integration_points;
}

double SnakeGapSbm3DUtilities::CalculateType1CharacteristicLength(
    const SurrogateFaceData& rFaceData,
    const NodePointerType& pApexNode) const
{
    KRATOS_ERROR_IF_NOT(pApexNode)
        << "[SnakeGapSbm3DUtilities::CalculateType1CharacteristicLength] "
        << "Apex node is null.\n";

    std::array<array_1d<double, 3>, 5> points;

    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_ERROR_IF_NOT(rFaceData.Nodes[i])
            << "[SnakeGapSbm3DUtilities::CalculateType1CharacteristicLength] "
            << "Null base node at local index " << i << ".\n";

        points[i] = rFaceData.Nodes[i]->Coordinates();
    }

    points[4] = pApexNode->Coordinates();

    double max_distance = 0.0;

    for (std::size_t i = 0; i < points.size(); ++i) {
        for (std::size_t j = i + 1; j < points.size(); ++j) {
            max_distance = std::max(
                max_distance,
                norm_2(points[i] - points[j]));
        }
    }

    return max_distance;
}

SnakeGapSbm3DUtilities::GeometriesArrayType
SnakeGapSbm3DUtilities::CreateAndTagVolumeQuadraturePointGeometries(
    const NurbsVolumeType::Pointer& pVolumeGeometry,
    const Geometry<Node>::Pointer& pNeighbourGeometry,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives) const
{
    KRATOS_ERROR_IF_NOT(pVolumeGeometry)
        << "[SnakeGapSbm3DUtilities::CreateAndTagVolumeQuadraturePointGeometries] "
        << "Volume geometry is null.\n";

    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::CreateAndTagVolumeQuadraturePointGeometries] "
        << "Neighbour geometry is null.\n";

    IntegrationPointsArrayType volume_integration_points =
        CreateCoonsVolumeGaussPoints(
            IntegrationOrder,
            *pVolumeGeometry);

    double weight_sum = 0.0;
    for (const auto& r_integration_point : volume_integration_points) {
        weight_sum += r_integration_point.Weight();
    }

    GeometriesArrayType volume_quadrature_point_geometries;

    if (weight_sum <= 1.0e-16) {
        KRATOS_WARNING("SnakeGapSbm3DUtilities")
            << "[CreateAndTagVolumeQuadraturePointGeometries] "
            << "Skipping volume quadrature point geometries because weight sum is nearly zero.\n";

        return volume_quadrature_point_geometries;
    }

    IntegrationInfo volume_integration_info(
        {IntegrationOrder, IntegrationOrder, IntegrationOrder},
        {IntegrationInfo::QuadratureMethod::CUSTOM,
         IntegrationInfo::QuadratureMethod::CUSTOM,
         IntegrationInfo::QuadratureMethod::CUSTOM});

    pVolumeGeometry->CreateQuadraturePointGeometries(
        volume_quadrature_point_geometries,
        NumberOfShapeFunctionsDerivatives,
        volume_integration_points,
        volume_integration_info);

    KRATOS_ERROR_IF(volume_quadrature_point_geometries.size() == 0)
        << "[SnakeGapSbm3DUtilities::CreateAndTagVolumeQuadraturePointGeometries] "
        << "Failed to create volume quadrature point geometries.\n";

    return volume_quadrature_point_geometries;
}

void SnakeGapSbm3DUtilities::RegisterType1LateralFace(
    const SurrogateFaceData& rFaceData,
    const IndexType LocalEdgeIndex,
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pApexNode,
    const NurbsSurfaceType::Pointer& pLateralSurface,
    const Geometry<Node>::Pointer& pNeighbourGeometry)
{
    KRATOS_ERROR_IF_NOT(pNode0)
        << "[SnakeGapSbm3DUtilities::RegisterType1LateralFace] First node is null.\n";
    KRATOS_ERROR_IF_NOT(pNode1)
        << "[SnakeGapSbm3DUtilities::RegisterType1LateralFace] Second node is null.\n";
    KRATOS_ERROR_IF_NOT(pApexNode)
        << "[SnakeGapSbm3DUtilities::RegisterType1LateralFace] Apex node is null.\n";
    KRATOS_ERROR_IF_NOT(pLateralSurface)
        << "[SnakeGapSbm3DUtilities::RegisterType1LateralFace] Surface is null.\n";
    KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::RegisterType1LateralFace] Neighbour geometry is null.\n";

    const auto key = MakeCanonicalFaceKey3D(
        pNode0,
        pNode1,
        pApexNode);

    LateralFaceOccurrence occurrence;
    occurrence.GapType = 1;
    occurrence.SurrogateConditionId = rFaceData.pCondition->Id();
    occurrence.ExternalSpan = rFaceData.ExternalSpan;
    occurrence.pGeometry = pLateralSurface;
    occurrence.pNeighbourGeometry = pNeighbourGeometry;
    occurrence.pOppositeNode = pApexNode;
    occurrence.FaceNodes = {{pNode0, pNode1, pApexNode}};
    occurrence.HasType1NeighbourPath = true;
    occurrence.HasNeighbourActiveSpan = true;
    occurrence.NeighbourActiveSpan = rFaceData.ActiveSpan;

    mLateralFaceRegistry[key].push_back(occurrence);
}

SnakeGapSbm3DUtilities::Type1CreationSummary
SnakeGapSbm3DUtilities::ComputeType1CreationSummary(
    const std::size_t NumberOfPyramids) const
{
    Type1CreationSummary summary;
    summary.NumberOfPyramids = NumberOfPyramids;

    for (const auto& r_entry : mLateralFaceRegistry) {
        const auto number_of_occurrences = r_entry.second.size();

        summary.NumberOfLateralFaces += number_of_occurrences;

        if (number_of_occurrences == 1) {
            ++summary.NumberOfOpenFaces;
        } else if (number_of_occurrences == 2) {
            ++summary.NumberOfClosedFaces;
        } else {
            ++summary.NumberOfNonManifoldFaces;
        }
    }

    return summary;
}

void SnakeGapSbm3DUtilities::PrintType1CreationSummary(
    const Type1CreationSummary& rSummary) const
{
    if (mEchoLevel == 0) {
        return;
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type 1 geometry creation summary:\n"
        << "  type 1 pyramids created:          " << rSummary.NumberOfPyramids << "\n"
        << "  type 1 lateral face occurrences:  " << rSummary.NumberOfLateralFaces << "\n"
        << "  open lateral faces:               " << rSummary.NumberOfOpenFaces << "\n"
        << "  closed lateral faces:             " << rSummary.NumberOfClosedFaces << "\n"
        << "  non-manifold lateral faces:       " << rSummary.NumberOfNonManifoldFaces << "\n";
}

SnakeGapSbm3DUtilities::NurbsSurfaceType::Pointer
SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurface(
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pNode2) const
{
    KRATOS_ERROR_IF_NOT(pNode0)
        << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurface] pNode0 is null.\n";

    KRATOS_ERROR_IF_NOT(pNode1)
        << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurface] pNode1 is null.\n";

    KRATOS_ERROR_IF_NOT(pNode2)
        << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurface] pNode2 is null.\n";

    if (mGapApproximationOrder > 1) {
        const bool has_edge_01 = FindCachedSkinEdgeControlNodes(pNode0, pNode1);
        const bool has_edge_12 = FindCachedSkinEdgeControlNodes(pNode1, pNode2);
        const bool has_edge_20 = FindCachedSkinEdgeControlNodes(pNode2, pNode0);
        const std::size_t number_of_curved_edges =
            static_cast<std::size_t>(has_edge_01) +
            static_cast<std::size_t>(has_edge_12) +
            static_cast<std::size_t>(has_edge_20);

        if (number_of_curved_edges == 3) {
            return CreateCollapsedTriangleSurfaceWithCurvedTop(
                pNode0,
                pNode1,
                pNode2);
        }

        if (number_of_curved_edges == 1 && has_edge_01) {
            auto skin_edge_control_nodes =
                *FindCachedSkinEdgeControlNodes(pNode0, pNode1);
            if (skin_edge_control_nodes.front()->Id() != pNode0->Id()) {
                std::reverse(
                    skin_edge_control_nodes.begin(),
                    skin_edge_control_nodes.end());
            }
            return CreateCollapsedTriangleSurfaceWithCurvedEdge(
                pNode2,
                skin_edge_control_nodes);
        }

        if (number_of_curved_edges == 1 && has_edge_12) {
            auto skin_edge_control_nodes =
                *FindCachedSkinEdgeControlNodes(pNode1, pNode2);
            if (skin_edge_control_nodes.front()->Id() != pNode1->Id()) {
                std::reverse(
                    skin_edge_control_nodes.begin(),
                    skin_edge_control_nodes.end());
            }
            return CreateCollapsedTriangleSurfaceWithCurvedEdge(
                pNode0,
                skin_edge_control_nodes);
        }

        if (number_of_curved_edges == 1 && has_edge_20) {
            auto skin_edge_control_nodes =
                *FindCachedSkinEdgeControlNodes(pNode2, pNode0);
            if (skin_edge_control_nodes.front()->Id() != pNode2->Id()) {
                std::reverse(
                    skin_edge_control_nodes.begin(),
                    skin_edge_control_nodes.end());
            }
            return CreateCollapsedTriangleSurfaceWithCurvedEdge(
                pNode1,
                skin_edge_control_nodes);
        }
    }

    PointerVector<NodeType> control_points;

    control_points.push_back(pNode0);
    control_points.push_back(pNode1);
    control_points.push_back(pNode2);
    control_points.push_back(pNode2);

    Vector knot_vector = ZeroVector(2);
    knot_vector[0] = 0.0;
    knot_vector[1] = 1.0;

    return Kratos::make_shared<NurbsSurfaceType>(
        control_points,
        1,
        1,
        knot_vector,
        knot_vector);
}

SnakeGapSbm3DUtilities::NurbsSurfaceType::Pointer
SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurfaceWithCurvedEdge(
    const NodePointerType& pApexNode,
    const std::vector<NodePointerType>& rSkinEdgeControlNodes) const
{
    KRATOS_ERROR_IF_NOT(pApexNode)
        << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurfaceWithCurvedEdge] "
        << "Apex node is null.\n";
    KRATOS_ERROR_IF(rSkinEdgeControlNodes.size() != mGapApproximationOrder + 1)
        << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurfaceWithCurvedEdge] "
        << "Unexpected number of skin edge control nodes: "
        << rSkinEdgeControlNodes.size() << ". Expected "
        << mGapApproximationOrder + 1 << ".\n";

    PointerVector<NodeType> control_points;

    for (const auto& p_control_node : rSkinEdgeControlNodes) {
        KRATOS_ERROR_IF_NOT(p_control_node)
            << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurfaceWithCurvedEdge] "
            << "Null skin edge control node.\n";
        control_points.push_back(p_control_node);
    }

    for (std::size_t i = 0; i <= mGapApproximationOrder; ++i) {
        control_points.push_back(pApexNode);
    }

    return Kratos::make_shared<NurbsSurfaceType>(
        control_points,
        mGapApproximationOrder,
        1,
        CreateOpenUnitKnotVector(mGapApproximationOrder),
        CreateOpenUnitKnotVectorDegree1());
}

SnakeGapSbm3DUtilities::NurbsSurfaceType::Pointer
SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurfaceWithCurvedTop(
    const NodePointerType& pNode0,
    const NodePointerType& pNode1,
    const NodePointerType& pNode2) const
{
    KRATOS_ERROR_IF(mGapApproximationOrder != 2)
        << "[SnakeGapSbm3DUtilities::CreateCollapsedTriangleSurfaceWithCurvedTop] "
        << "Only quadratic collapsed top faces are implemented. Requested order: "
        << mGapApproximationOrder << ".\n";

    auto get_oriented_edge = [this](
        const NodePointerType& pFirst,
        const NodePointerType& pSecond)
    {
        auto control_nodes = *FindCachedSkinEdgeControlNodes(pFirst, pSecond);
        if (control_nodes.front()->Id() != pFirst->Id()) {
            std::reverse(control_nodes.begin(), control_nodes.end());
        }
        return control_nodes;
    };

    const auto edge_01 = get_oriented_edge(pNode0, pNode1);
    const auto edge_12 = get_oriented_edge(pNode1, pNode2);
    const auto edge_02 = get_oriented_edge(pNode0, pNode2);

    array_1d<double, 3> middle_point = ZeroVector(3);
    noalias(middle_point) += 0.5 * edge_01[1]->Coordinates();
    noalias(middle_point) += 0.5 * edge_02[1]->Coordinates();
    noalias(middle_point) += 0.5 * edge_12[1]->Coordinates();
    noalias(middle_point) -= 0.25 * pNode0->Coordinates();
    noalias(middle_point) -= 0.25 * pNode1->Coordinates();

    NodePointerType p_middle_node(new Node(0, middle_point));

    PointerVector<NodeType> control_points;

    control_points.push_back(edge_01[0]);
    control_points.push_back(edge_01[1]);
    control_points.push_back(edge_01[2]);

    control_points.push_back(edge_02[1]);
    control_points.push_back(p_middle_node);
    control_points.push_back(edge_12[1]);

    control_points.push_back(pNode2);
    control_points.push_back(pNode2);
    control_points.push_back(pNode2);

    return Kratos::make_shared<NurbsSurfaceType>(
        control_points,
        std::size_t(2),
        std::size_t(2),
        CreateOpenUnitKnotVector(2),
        CreateOpenUnitKnotVector(2));
}

SnakeGapSbm3DUtilities::Type1CreationResult
SnakeGapSbm3DUtilities::CreateType1GapGeometries(
    ModelPart& rRootModelPart,
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart,
    const ExternalSpanDataMap& rExternalSpans,
    const std::size_t GapVolumeIntegrationOrder)
{
    const char* p_caller =
        "SnakeGapSbm3DUtilities::CreateType1GapGeometries";

    ClearLateralFaceRegistry();

    ModelPart& r_gap_type1_debug = GetOrCreateSubModelPart(
        rRootModelPart,
        "GapType1Debug");

    const auto grid_info = CreateKnotSpanGridInfo(
        rSurrogateSubModelPart.GetParentModelPart());

    const auto surrogate_faces = BuildSurrogateFaceDataVector(
        rSurrogateSubModelPart,
        grid_info);

    IndexType next_geometry_id = GetNextGeometryId(rRootModelPart);

    std::size_t number_of_type1_pyramids = 0;

    Type1CreationResult result;

    const ModelPart& r_iga_model_part = rSurrogateSubModelPart.GetRootModelPart();

    // const std::size_t number_of_shape_functions_derivatives =
    //     ComputeNumberOfShapeFunctionsDerivatives(r_iga_model_part);
    
    const std::size_t  number_of_shape_functions_derivatives = 9; //FIXME:

    const std::vector<BrepPatchData> brep_patch_data_list =
        BuildBrepPatchDataVector(
            r_iga_model_part,
            rSurrogateSubModelPart);
    
    result.VolumeQuadratureDataList.reserve(surrogate_faces.size());
    mLateralFaceRegistry.reserve(4 * surrogate_faces.size());

    for (const auto& r_face_data : surrogate_faces) {
        const auto external_span_it = rExternalSpans.find(r_face_data.ExternalSpan);

        KRATOS_ERROR_IF(external_span_it == rExternalSpans.end())
            << "[" << p_caller << "] External span "
            << SpanToString(r_face_data.ExternalSpan)
            << " from surrogate condition #" << r_face_data.pCondition->Id()
            << " was not found in ExternalSpanDataMap.\n";

        const auto& r_external_span_data = external_span_it->second;

        if (r_external_span_data.Type != GapSpanType::Type1) {
            continue;
        }

        KRATOS_ERROR_IF_NOT(r_external_span_data.HasProjectionNode())
            << "[" << p_caller << "] Type 1 span "
            << SpanToString(r_face_data.ExternalSpan)
            << " from surrogate condition #" << r_face_data.pCondition->Id()
            << " has no PROJECTION_NODE_ID.\n";

        const IndexType projection_node_id = r_external_span_data.ProjectionNodeId;

        KRATOS_ERROR_IF_NOT(rSkinSubModelPart.HasNode(projection_node_id))
            << "[" << p_caller << "] Projection node #" << projection_node_id
            << " for type 1 span " << SpanToString(r_face_data.ExternalSpan)
            << " is not in skin sub model part '" << rSkinSubModelPart.Name() << "'.\n";

        const auto p_apex_node = rSkinSubModelPart.pGetNode(projection_node_id);

        const auto& r_brep_patch_data = FindBrepPatchMatchingCondition(
            *r_face_data.pCondition,
            brep_patch_data_list);
        
        auto p_neighbour_geometry = CreateSurrogateFaceNeighbourGeometry(
            r_brep_patch_data,
            number_of_shape_functions_derivatives);

        auto p_volume_geometry = CreateType1CollapsedPyramidVolume(
            r_face_data,
            p_apex_node);

        p_volume_geometry->SetId(next_geometry_id++);
        SetType1GeometryData(
            *p_volume_geometry,
            projection_node_id,
            p_neighbour_geometry,
            "TYPE_1_VOLUME");

        if (mStoreGapDebugGeometries) {
            r_gap_type1_debug.AddGeometry(p_volume_geometry);
        }

        GeometriesArrayType volume_quadrature_point_geometries =
        CreateAndTagVolumeQuadraturePointGeometries(
            p_volume_geometry,
            p_neighbour_geometry,
            GapVolumeIntegrationOrder,
            number_of_shape_functions_derivatives);

        if (volume_quadrature_point_geometries.size() > 0) {
            Type1VolumeQuadratureData type1_volume_data;

            type1_volume_data.VolumeQuadraturePointGeometries = std::move(volume_quadrature_point_geometries);

            type1_volume_data.NeighbourGeometries.push_back(p_neighbour_geometry);

            type1_volume_data.CharacteristicLength =
                CalculateType1CharacteristicLength(
                    r_face_data,
                    p_apex_node);

            type1_volume_data.BaseNodes = r_face_data.Nodes;
            type1_volume_data.pProjectionNode = p_apex_node;

            type1_volume_data.ProjectionNodeId = projection_node_id;
            type1_volume_data.SurrogateConditionId = r_face_data.pCondition->Id();
            type1_volume_data.ExternalSpan = r_face_data.ExternalSpan;

            result.VolumeQuadratureDataList.emplace_back(std::move(type1_volume_data));
        }

        const std::array<std::array<IndexType, 2>, 4> edge_node_indices = {{
            {{0, 1}},
            {{1, 2}},
            {{2, 3}},
            {{3, 0}}
        }};

        const array_1d<double, 3> volume_center = p_volume_geometry->Center().Coordinates();

        for (IndexType local_edge_index = 0; local_edge_index < 4; ++local_edge_index) {
            const auto local_node_0 = edge_node_indices[local_edge_index][0];
            const auto local_node_1 = edge_node_indices[local_edge_index][1];

            auto p_node_0 = r_face_data.Nodes[local_node_0];
            auto p_node_1 = r_face_data.Nodes[local_node_1];
            auto grid_node_0 = r_face_data.GridNodes[local_node_0];
            auto grid_node_1 = r_face_data.GridNodes[local_node_1];

            const array_1d<double, 3> edge_vector =
                p_node_1->Coordinates() - p_node_0->Coordinates();
            const array_1d<double, 3> apex_vector =
                p_apex_node->Coordinates() - p_node_0->Coordinates();
            const array_1d<double, 3> lateral_normal =
                MathUtils<double>::CrossProduct(edge_vector, apex_vector);

            array_1d<double, 3> face_center = ZeroVector(3);
            noalias(face_center) += p_node_0->Coordinates();
            noalias(face_center) += p_node_1->Coordinates();
            noalias(face_center) += p_apex_node->Coordinates();
            face_center /= 3.0;
            const array_1d<double, 3> outward_vector =
                face_center - volume_center;

            if (inner_prod(lateral_normal, outward_vector) < 0.0) {
                std::swap(p_node_0, p_node_1);
                std::swap(grid_node_0, grid_node_1);
            }

            auto p_lateral_surface = CreateType1CollapsedLateralCoonsSurface(
                p_node_0,
                p_node_1,
                p_apex_node);

            p_lateral_surface->SetId(next_geometry_id++);

            // type 2 debug
            Type1LateralFaceData lateral_face_data;

            lateral_face_data.EdgeKey = SurrogateEdgeKey3D(
                p_node_0->Id(),
                p_node_1->Id());

            lateral_face_data.pEdgeNode0 = p_node_0;
            lateral_face_data.pEdgeNode1 = p_node_1;
            lateral_face_data.pProjectionNode = p_apex_node;

            lateral_face_data.EdgeGridNode0 = grid_node_0;
            lateral_face_data.EdgeGridNode1 = grid_node_1;

            lateral_face_data.pSurfaceGeometry = p_lateral_surface;
            lateral_face_data.pNeighbourGeometry = p_neighbour_geometry;

            lateral_face_data.SurrogateConditionId = r_face_data.pCondition->Id();
            lateral_face_data.ActiveSpan = r_face_data.ActiveSpan;
            lateral_face_data.ExternalSpan = r_face_data.ExternalSpan;

            const std::size_t lateral_face_index =
                result.Type1LateralFaces.size();

            result.Type1LateralFaces.push_back(lateral_face_data);

            result.Type1LateralFacesByEdge[lateral_face_data.EdgeKey].push_back(
                lateral_face_index);
            //---------------------------------------------

            SetType1GeometryData(
                *p_lateral_surface,
                projection_node_id,
                p_neighbour_geometry,
                "TYPE_1_LATERAL_SURFACE");

            if (mStoreGapDebugGeometries) {
                r_gap_type1_debug.AddGeometry(p_lateral_surface);
            }

            CreateAndTagSurfaceQuadraturePointGeometries(
                p_lateral_surface,
                p_neighbour_geometry);

            RegisterType1LateralFace(
                r_face_data,
                local_edge_index,
                p_node_0,
                p_node_1,
                p_apex_node,
                p_lateral_surface,
                p_neighbour_geometry);
        }

        ++number_of_type1_pyramids;
    }

    result.Summary = ComputeType1CreationSummary(number_of_type1_pyramids);

    PrintType1CreationSummary(result.Summary);

    return result;
}

// create quadrature point for surfaces
SnakeGapSbm3DUtilities::IntegrationPointsArrayType
SnakeGapSbm3DUtilities::CreateCoonsSurfaceGaussPoints(
    const std::size_t IntegrationOrder,
    const NurbsSurfaceType& rGapSurface) const
{
    KRATOS_ERROR_IF(IntegrationOrder == 0)
        << "[SnakeGapSbm3DUtilities::CreateCoonsSurfaceGaussPoints] "
        << "Integration order must be positive.\n";

    IntegrationPointsArrayType integration_points(
        IntegrationOrder * IntegrationOrder);

    auto integration_point_it = integration_points.begin();

    IntegrationPointUtilities::IntegrationPoints2D(
        integration_point_it,
        IntegrationOrder,
        IntegrationOrder,
        0.0,
        1.0,
        0.0,
        1.0);

    double raw_weight_sum = 0.0;

    for (auto& r_integration_point : integration_points) {
        CoordinatesArrayType local_coordinates = ZeroVector(3);
        local_coordinates[0] = r_integration_point[0];
        local_coordinates[1] = r_integration_point[1];

        Matrix jacobian;
        rGapSurface.Jacobian(jacobian, local_coordinates);

        KRATOS_ERROR_IF(jacobian.size1() < 3 || jacobian.size2() < 2)
            << "[SnakeGapSbm3DUtilities::CreateCoonsSurfaceGaussPoints] "
            << "Invalid surface Jacobian size: "
            << jacobian.size1() << " x " << jacobian.size2() << ".\n";

        array_1d<double, 3> tangent_u = ZeroVector(3);
        array_1d<double, 3> tangent_v = ZeroVector(3);

        for (IndexType i = 0; i < 3; ++i) {
            tangent_u[i] = jacobian(i, 0);
            tangent_v[i] = jacobian(i, 1);
        }

        const array_1d<double, 3> normal =
            MathUtils<double>::CrossProduct(tangent_u, tangent_v);

        const double surface_jacobian = norm_2(normal);

        r_integration_point.SetWeight(
            r_integration_point.Weight() * surface_jacobian);

        raw_weight_sum += r_integration_point.Weight();
    }

    const double exact_area = rGapSurface.Area();

    return integration_points;
}

SnakeGapSbm3DUtilities::GeometriesArrayType
SnakeGapSbm3DUtilities::CreateSurfaceQuadraturePointGeometries(
    const NurbsSurfaceType::Pointer& pSurfaceGeometry,
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives) const
{
    KRATOS_ERROR_IF_NOT(pSurfaceGeometry)
        << "[SnakeGapSbm3DUtilities::CreateSurfaceQuadraturePointGeometries] "
        << "Surface geometry is null.\n";

    IntegrationPointsArrayType surface_integration_points =
        CreateCoonsSurfaceGaussPoints(
            IntegrationOrder,
            *pSurfaceGeometry);

    GeometriesArrayType surface_quadrature_point_geometries;

    double weight_sum = 0.0;
    for (const auto& r_integration_point : surface_integration_points) {
        weight_sum += r_integration_point.Weight();
    }

    if (weight_sum <= 1.0e-16) {
        KRATOS_WARNING("SnakeGapSbm3DUtilities")
            << "[CreateSurfaceQuadraturePointGeometries] "
            << "Skipping surface quadrature point geometries because weight sum is nearly zero.\n";

        return surface_quadrature_point_geometries;
    }

    IntegrationInfo surface_integration_info(
        {IntegrationOrder, IntegrationOrder},
        {IntegrationInfo::QuadratureMethod::CUSTOM,
         IntegrationInfo::QuadratureMethod::CUSTOM});

    pSurfaceGeometry->CreateQuadraturePointGeometries(
        surface_quadrature_point_geometries,
        NumberOfShapeFunctionsDerivatives,
        surface_integration_points,
        surface_integration_info);

    KRATOS_ERROR_IF(surface_quadrature_point_geometries.size() == 0)
        << "[SnakeGapSbm3DUtilities::CreateSurfaceQuadraturePointGeometries] "
        << "Failed to create surface quadrature point geometries.\n";

    return surface_quadrature_point_geometries;
}

bool SnakeGapSbm3DUtilities::ContainsNeighbourGeometry(
    const std::vector<Geometry<Node>::Pointer>& rNeighbourGeometries,
    const Geometry<Node>::Pointer& pCandidateNeighbourGeometry) const
{
    if (!pCandidateNeighbourGeometry) {
        return false;
    }

    for (const auto& p_existing_neighbour_geometry : rNeighbourGeometries) {
        if (p_existing_neighbour_geometry.get() == pCandidateNeighbourGeometry.get()) {
            return true;
        }
    }

    return false;
}

void SnakeGapSbm3DUtilities::AddUniqueNeighbourGeometry(
    std::vector<Geometry<Node>::Pointer>& rNeighbourGeometries,
    const Geometry<Node>::Pointer& pCandidateNeighbourGeometry) const
{
    KRATOS_ERROR_IF_NOT(pCandidateNeighbourGeometry)
        << "[SnakeGapSbm3DUtilities::AddUniqueNeighbourGeometry] "
        << "Trying to add null neighbour geometry.\n";

    if (!ContainsNeighbourGeometry(rNeighbourGeometries, pCandidateNeighbourGeometry)) {
        rNeighbourGeometries.push_back(pCandidateNeighbourGeometry);
    }
}

std::vector<Geometry<Node>::Pointer>
SnakeGapSbm3DUtilities::CollectUniqueNeighbourGeometries(
    const std::vector<LateralFaceOccurrence>& rOccurrences) const
{
    std::vector<Geometry<Node>::Pointer> neighbour_geometries;
    neighbour_geometries.reserve(rOccurrences.size());

    for (const auto& r_occurrence : rOccurrences) {
        AddUniqueNeighbourGeometry(
            neighbour_geometries,
            r_occurrence.pNeighbourGeometry);
    }

    return neighbour_geometries;
}

double SnakeGapSbm3DUtilities::CalculateSurfaceCharacteristicLength(
    const Geometry<Node>& rGeometry) const
{
    double characteristic_length = 0.0;

    for (IndexType i = 0; i < rGeometry.size(); ++i) {
        for (IndexType j = i + 1; j < rGeometry.size(); ++j) {
            characteristic_length = std::max(
                characteristic_length,
                norm_2(rGeometry[i].Coordinates() - rGeometry[j].Coordinates()));
        }
    }

    return characteristic_length;
}

std::vector<SnakeGapSbm3DUtilities::LateralSurfaceQuadratureData>
SnakeGapSbm3DUtilities::CreateOpenLateralSurfaceQuadratureData(
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives) const
{
    std::vector<LateralSurfaceQuadratureData> data_list;

    for (const auto& r_entry : mLateralFaceRegistry) {
        const auto& r_occurrences = r_entry.second;

        KRATOS_ERROR_IF(r_occurrences.size() > 2)
            << "[SnakeGapSbm3DUtilities::CreateOpenLateralSurfaceQuadratureData] "
            << "Non-manifold lateral face with "
            << r_occurrences.size() << " occurrences.\n";

        if (r_occurrences.size() != 1) {
            continue;
        }

        const auto& r_occurrence = r_occurrences.front();

        KRATOS_ERROR_IF_NOT(r_occurrence.pGeometry)
            << "[SnakeGapSbm3DUtilities::CreateOpenLateralSurfaceQuadratureData] "
            << "Null lateral surface geometry.\n";

        LateralSurfaceQuadratureData data;

        data.SurfaceQuadraturePointGeometries =
            CreateSurfaceQuadraturePointGeometries(
                r_occurrence.pGeometry,
                IntegrationOrder,
                NumberOfShapeFunctionsDerivatives);

        if (data.SurfaceQuadraturePointGeometries.size() == 0) {
            continue;
        }

        AddUniqueNeighbourGeometry(
            data.NeighbourGeometries,
            r_occurrence.pNeighbourGeometry);

        data.CharacteristicLength =
            CalculateSurfaceCharacteristicLength(*r_occurrence.pGeometry);

        data.GapType = r_occurrence.GapType;
        data.SurrogateConditionId = r_occurrence.SurrogateConditionId;
        data.ExternalSpan = r_occurrence.ExternalSpan;

        data_list.push_back(std::move(data));
    }

    return data_list;
}

std::vector<SnakeGapSbm3DUtilities::LateralSurfaceQuadratureData>
SnakeGapSbm3DUtilities::CreateInterfaceLateralSurfaceQuadratureData(
    const std::size_t IntegrationOrder,
    const std::size_t NumberOfShapeFunctionsDerivatives) const
{
    std::vector<LateralSurfaceQuadratureData> data_list;

    for (const auto& r_entry : mLateralFaceRegistry) {
        const auto& r_occurrences = r_entry.second;

        KRATOS_ERROR_IF(r_occurrences.size() > 2)
            << "[SnakeGapSbm3DUtilities::CreateInterfaceLateralSurfaceQuadratureData] "
            << "Non-manifold lateral face with "
            << r_occurrences.size() << " occurrences.\n";

        if (r_occurrences.size() != 2) {
            continue;
        }

        std::vector<Geometry<Node>::Pointer> neighbour_geometries =
            CollectUniqueNeighbourGeometries(r_occurrences);

        KRATOS_ERROR_IF(neighbour_geometries.size() == 0)
            << "[SnakeGapSbm3DUtilities::CreateInterfaceLateralSurfaceQuadratureData] "
            << "No neighbour geometries found for lateral face with "
            << r_occurrences.size() << " occurrences.\n";

        if (neighbour_geometries.size() <= 1) {
            continue;
        }

        const auto& r_reference_occurrence = r_occurrences.front();

        KRATOS_ERROR_IF_NOT(r_reference_occurrence.pGeometry)
            << "[SnakeGapSbm3DUtilities::CreateInterfaceLateralSurfaceQuadratureData] "
            << "Null lateral surface geometry.\n";

        LateralSurfaceQuadratureData data;

        data.SurfaceQuadraturePointGeometries =
            CreateSurfaceQuadraturePointGeometries(
                r_reference_occurrence.pGeometry,
                IntegrationOrder,
                NumberOfShapeFunctionsDerivatives);

        if (data.SurfaceQuadraturePointGeometries.size() == 0) {
            continue;
        }

        data.NeighbourGeometries = std::move(neighbour_geometries);

        data.CharacteristicLength =
            CalculateSurfaceCharacteristicLength(*r_reference_occurrence.pGeometry);

        data.GapType = r_reference_occurrence.GapType;
        data.SurrogateConditionId = r_reference_occurrence.SurrogateConditionId;
        data.ExternalSpan = r_reference_occurrence.ExternalSpan;

        data_list.push_back(std::move(data));
    }

    return data_list;
}


template <bool TIsInnerLoop>
void SnakeGapSbmProcess::CreateSbmExtendedGeometries3D(
    ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart)
{
    SnakeGapSbm3DUtilities utilities(mEchoLevel);
    utilities.SetGapApproximationOrder(mGapApproximationOrder);
    const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    KRATOS_ERROR_IF(knot_span_sizes.size() < 3)
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: KNOT_SPAN_SIZES must contain at least three entries." << std::endl;

    const double min_knot_span_size = std::min(knot_span_sizes[0], std::min(knot_span_sizes[1], knot_span_sizes[2]));
    const double max_knot_span_size = std::max(knot_span_sizes[0], std::max(knot_span_sizes[1], knot_span_sizes[2]));
    KRATOS_ERROR_IF(min_knot_span_size <= 0.0)
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: knot span sizes must be positive." << std::endl;

    const double area_tol = 1.0e-12 * min_knot_span_size * min_knot_span_size;
    const double volume_tol = 1.0e-12 * min_knot_span_size * min_knot_span_size * min_knot_span_size;
    const double distance_tol = 1.0e-12 * min_knot_span_size;
    const double determinant_tol = 1.0e-14 * std::max(1.0, max_knot_span_size * max_knot_span_size * max_knot_span_size);

    const std::size_t gap_volume_integration_order = mGapInterpolationOrder+1; //TODO: make this an input parameter


    ModelPart& r_root_model_part = mpIgaModelPart->GetRootModelPart();
    const bool store_gap_debug_geometries =
        mThisParameters.Has("store_gap_debug_geometries")
            ? mThisParameters["store_gap_debug_geometries"].GetBool()
            : false;
    utilities.SetStoreGapDebugGeometries(store_gap_debug_geometries);

    ModelPart* p_gap_type1_debug = nullptr;
    ModelPart* p_gap_type2_debug = nullptr;
    ModelPart* p_gap_type3_debug = nullptr;
    IndexType next_debug_geometry_id = 1;
    std::unordered_map<const NodeType*, Node::Pointer> debug_node_map;
    IndexType next_debug_node_id = 1;   

    if (store_gap_debug_geometries) {
        p_gap_type1_debug = r_root_model_part.HasSubModelPart("GapType1Debug")
            ? &r_root_model_part.GetSubModelPart("GapType1Debug")
            : &r_root_model_part.CreateSubModelPart("GapType1Debug");
        p_gap_type2_debug = r_root_model_part.HasSubModelPart("GapType2Debug")
            ? &r_root_model_part.GetSubModelPart("GapType2Debug")
            : &r_root_model_part.CreateSubModelPart("GapType2Debug");
        p_gap_type3_debug = r_root_model_part.HasSubModelPart("GapType3Debug")
            ? &r_root_model_part.GetSubModelPart("GapType3Debug")
            : &r_root_model_part.CreateSubModelPart("GapType3Debug");

        for (const auto& r_geometry : r_root_model_part.Geometries()) {
            next_debug_geometry_id = std::max(next_debug_geometry_id, r_geometry.Id() + 1);
        }
        for (const auto& r_node : r_root_model_part.Nodes()) {
            next_debug_node_id = std::max(next_debug_node_id, r_node.Id() + 1);
        }
    }

    auto get_next_debug_geometry_id = [&]() {
        while (r_root_model_part.HasGeometry(next_debug_geometry_id)) {
            ++next_debug_geometry_id;
        }
        return next_debug_geometry_id++;
    };

    auto get_next_debug_node_id = [&]() {
        while (r_root_model_part.HasNode(next_debug_node_id)) {
            ++next_debug_node_id;
        }
        return next_debug_node_id++;
    };
    
    auto get_or_create_debug_node = [&](ModelPart& rDebugPart, const Node::Pointer& pOriginalNode) -> Node::Pointer {
        KRATOS_ERROR_IF_NOT(pOriginalNode)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null original node in debug geometry."
            << std::endl;
    
        const NodeType* p_key = pOriginalNode.get();
        auto it = debug_node_map.find(p_key);
        if (it != debug_node_map.end()) {
            if (!rDebugPart.HasNode(it->second->Id())) {
                rDebugPart.AddNode(it->second);
            }
            return it->second;
        }
    
        if (r_root_model_part.HasNode(pOriginalNode->Id())) {
            const auto& r_root_node_same_id = r_root_model_part.GetNode(pOriginalNode->Id());
            const double coordinate_mismatch =
                norm_2(r_root_node_same_id.Coordinates() - pOriginalNode->Coordinates());
        }
    
        const IndexType new_debug_node_id = get_next_debug_node_id();
        auto p_debug_node = r_root_model_part.CreateNewNode(
            new_debug_node_id,
            pOriginalNode->X(),
            pOriginalNode->Y(),
            pOriginalNode->Z());
    
        rDebugPart.AddNode(p_debug_node);
        debug_node_map.emplace(p_key, p_debug_node);
        return p_debug_node;
    };

    auto p_volume = mpIgaModelPart->pGetGeometry(1);
    KRATOS_ERROR_IF_NOT(p_volume)
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: geometry with id 1 was not found." << std::endl;

    auto p_nurbs_volume = std::dynamic_pointer_cast<NurbsVolumeType>(
        p_volume->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    KRATOS_ERROR_IF_NOT(p_nurbs_volume)
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: geometry with id 1 does not expose a NurbsVolumeType background geometry." << std::endl;

    const auto bounds = ReadParameterSpaceBounds3D(
        rSurrogateSubModelPart.GetParentModelPart(),
        "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D");

    const int number_of_spans_x = static_cast<int>(ComputeSpanCount3D(
        bounds.MaxU - bounds.MinU,
        knot_span_sizes[0],
        "u",
        "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D"));
    
    const int number_of_spans_y = static_cast<int>(ComputeSpanCount3D(
        bounds.MaxV - bounds.MinV,
        knot_span_sizes[1],
        "v",
        "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D"));
    
    const int number_of_spans_z = static_cast<int>(ComputeSpanCount3D(
        bounds.MaxW - bounds.MinW,
        knot_span_sizes[2],
        "w",
        "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D"));

    const std::size_t brep_degree = p_nurbs_volume->PolynomialDegree(0);
    const std::size_t number_of_shape_functions_derivatives = 6 * brep_degree + 1;
    // if (mGapApproximationOrder == 0) {
    //     mGapApproximationOrder = brep_degree;
    // }

    const auto grid_info = utilities.CreateKnotSpanGridInfo(
        rSurrogateSubModelPart.GetParentModelPart());

    //-------------------------------------------------
    const double span_tolerance =
        1.0e-10 * std::max({knot_span_sizes[0], knot_span_sizes[1], knot_span_sizes[2]});
    
    auto is_point_inside_box = [&](const array_1d<double, 3>& rPoint,
                                   const array_1d<double, 3>& rBoxMin,
                                   const array_1d<double, 3>& rBoxMax) {
        return rPoint[0] >= rBoxMin[0] - span_tolerance && rPoint[0] <= rBoxMax[0] + span_tolerance &&
               rPoint[1] >= rBoxMin[1] - span_tolerance && rPoint[1] <= rBoxMax[1] + span_tolerance &&
               rPoint[2] >= rBoxMin[2] - span_tolerance && rPoint[2] <= rBoxMax[2] + span_tolerance;
    };
    //----------------------------------------------------

    PointVector skin_condition_points;
    skin_condition_points.reserve(3 * rSkinSubModelPart.NumberOfConditions());

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();

        KRATOS_ERROR_IF(r_geometry.PointsNumber() < 3)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "IsPointInsideSkinBoundary3D requires triangular skin conditions. "
            << "Condition #" << r_condition.Id()
            << " has " << r_geometry.PointsNumber() << " points."
            << std::endl;

        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            skin_condition_points.push_back(Kratos::make_intrusive<PointType>(
                r_condition.Id(),
                r_geometry[i_node].X(),
                r_geometry[i_node].Y(),
                r_geometry[i_node].Z()));
        }
    }

    DynamicBins skin_points_bin(
        skin_condition_points.begin(),
        skin_condition_points.end());
    
    // Use the same fake-Gauss-point criterion as SnakeSbmProcess::MarkKnotSpansAvailable3D.
    const double lambda = TIsInnerLoop ? mLambdaInner : mLambdaOuter;
    auto is_span_active_from_skin_center =
    [&](const SnakeGapSbm3DUtilities::SpanKey3D& rSpan) -> bool
    {
        const int num_fake_gauss_points = 2;
        const int total_fake_gauss_points =
            num_fake_gauss_points * num_fake_gauss_points * num_fake_gauss_points;

        int number_of_inside_gaussian_points = 0;

        const double tolerance = min_knot_span_size / 1.0e8;

        for (IndexType i_GPx = 0; i_GPx < num_fake_gauss_points; ++i_GPx) {
            const double x_coord =
                (static_cast<double>(rSpan.I) * knot_span_sizes[0] + tolerance) +
                (knot_span_sizes[0] - 2.0 * tolerance) /
                    (num_fake_gauss_points - 1) * i_GPx +
                bounds.MinU;

            for (IndexType i_GPy = 0; i_GPy < num_fake_gauss_points; ++i_GPy) {
                const double y_coord =
                    (static_cast<double>(rSpan.J) * knot_span_sizes[1] + tolerance) +
                    (knot_span_sizes[1] - 2.0 * tolerance) /
                        (num_fake_gauss_points - 1) * i_GPy +
                    bounds.MinV;

                for (IndexType i_GPz = 0; i_GPz < num_fake_gauss_points; ++i_GPz) {
                    const double z_coord =
                        (static_cast<double>(rSpan.K) * knot_span_sizes[2] + tolerance) +
                        (knot_span_sizes[2] - 2.0 * tolerance) /
                            (num_fake_gauss_points - 1) * i_GPz +
                        bounds.MinW;

                    PointType gauss_point(0, x_coord, y_coord, z_coord);

                    if (this->IsPointInsideSkinBoundary3D(
                            gauss_point,
                            skin_points_bin,
                            rSkinSubModelPart)) {
                        ++number_of_inside_gaussian_points;
                    }
                }
            }
        }

        return number_of_inside_gaussian_points >=
               lambda * total_fake_gauss_points;
    };

    //FIXME: DEBUG, REMOVE
    for (const auto& r_condition : rSurrogateSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();

        for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
            const auto& r_node = r_geometry[i];

            PointType point(0, r_node.X(), r_node.Y(), r_node.Z());

            if (!this->IsPointInsideSkinBoundary3D(
                    point,
                    skin_points_bin,
                    rSkinSubModelPart)) {

                KRATOS_ERROR
                    << "[SnakeGapSbmProcess] Surrogate boundary contains a node outside the skin.\n"
                    << "  surrogate condition id: " << r_condition.Id() << "\n"
                    << "  local node index: " << i << "\n"
                    << "  node id: " << r_node.Id() << "\n"
                    << "  coordinates: " << r_node.Coordinates() << "\n";
            }
        }
    }

    const auto external_spans = utilities.InitializeExternalSpanData(
        rSkinSubModelPart,
        rSurrogateSubModelPart,
        is_span_active_from_skin_center);

    
    auto type1_creation_result = utilities.CreateType1GapGeometries(
        r_root_model_part,
        rSkinSubModelPart,
        rSurrogateSubModelPart,
        external_spans,
        gap_volume_integration_order);

    
    IndexType id_element = 1;
    for (const auto& r_element : r_root_model_part.Elements()) {
        id_element = std::max(id_element, r_element.Id() + 1);
    }

    auto p_properties = mpIgaModelPart->pGetProperties(1);

    for (auto& r_type1_data : type1_creation_result.VolumeQuadratureDataList) {
        this->CreateElements(
            r_type1_data.VolumeQuadraturePointGeometries.ptr_begin(),
            r_type1_data.VolumeQuadraturePointGeometries.ptr_end(),
            *mpGapElementsSubModelPart,
            std::string("GapSbmSolidElement"),
            id_element,
            p_properties,
            r_type1_data.NeighbourGeometries,
            r_type1_data.CharacteristicLength);
    }

    std::map<SnakeGapSbm3DUtilities::SpanKey3D, std::size_t> type1_lateral_faces_per_span;
    for (const auto& r_type1_lateral_face : type1_creation_result.Type1LateralFaces) {
        ++type1_lateral_faces_per_span[r_type1_lateral_face.ExternalSpan];
    }

    std::vector<SnakeGapSbm3DUtilities::SpanKey3D> missing_type1_pyramid_spans;
    for (const auto& r_external_span_entry : external_spans) {
        const auto& r_external_span_data = r_external_span_entry.second;
        if (r_external_span_data.Type != SnakeGapSbm3DUtilities::GapSpanType::Type1) {
            continue;
        }

        if (type1_lateral_faces_per_span.find(r_external_span_data.Key) ==
            type1_lateral_faces_per_span.end()) {
            missing_type1_pyramid_spans.push_back(r_external_span_data.Key);
        }
    }

    if (!missing_type1_pyramid_spans.empty()) {
        std::ostringstream buffer;
        buffer
            << "[SnakeGapSbmProcess::CreateSbmExtendedGeometries3D] "
            << "Missing type1 pyramid creation for "
            << missing_type1_pyramid_spans.size()
            << " external TYPE_1 knot spans after type1 element creation.\n";

        for (const auto& r_span : missing_type1_pyramid_spans) {
            buffer << "  missing TYPE_1 span: "
                   << utilities.SpanToString(r_span) << "\n";
        }

        KRATOS_ERROR << buffer.str();
    }

    auto type2_creation_result = utilities.CreateType2GapGeometries(
        r_root_model_part,
        rSkinSubModelPart,
        type1_creation_result,
        external_spans,
        grid_info,
        gap_volume_integration_order,
        number_of_shape_functions_derivatives);

    auto type3_creation_result = utilities.CreateType3GapGeometries(
        r_root_model_part,
        rSkinSubModelPart,
        external_spans,
        grid_info,
        type2_creation_result,
        gap_volume_integration_order,
        number_of_shape_functions_derivatives);

    KRATOS_WATCH("PRE ELEMENTS TYPE 2 & 3")
    utilities.FinalizeType2AndType3GapGeometries(
        r_root_model_part,
        rSkinSubModelPart,
        type2_creation_result,
        type3_creation_result,
        grid_info,
        gap_volume_integration_order,
        number_of_shape_functions_derivatives);

        KRATOS_WATCH("POST ELEMENTS TYPE 2 & 3")

    for (auto& r_type2_data : type2_creation_result.VolumeQuadratureDataList) {
        if (r_type2_data.VolumeQuadraturePointGeometries.size() == 0) {
            continue;
        }

        this->CreateElements(
            r_type2_data.VolumeQuadraturePointGeometries.ptr_begin(),
            r_type2_data.VolumeQuadraturePointGeometries.ptr_end(),
            *mpGapElementsSubModelPart,
            std::string("GapSbmSolidElement"),
            id_element,
            p_properties,
            r_type2_data.NeighbourGeometries,
            r_type2_data.CharacteristicLength);
    }

    for (auto& r_type3_data : type3_creation_result.VolumeQuadratureDataList) {
        if (r_type3_data.VolumeQuadraturePointGeometries.size() == 0) {
            continue;
        }

        this->CreateElements(
            r_type3_data.VolumeQuadraturePointGeometries.ptr_begin(),
            r_type3_data.VolumeQuadraturePointGeometries.ptr_end(),
            *mpGapElementsSubModelPart,
            std::string("GapSbmSolidElement"),
            id_element,
            p_properties,
            r_type3_data.NeighbourGeometries,
            r_type3_data.CharacteristicLength);
    }

    std::size_t number_of_unclosed_type2_faces = 0;
    std::size_t number_of_type2_type3_closed_faces = 0;
    std::size_t number_of_already_closed_non_type3_faces = 0;
    std::size_t number_of_non_manifold_faces = 0;

    auto make_process_canonical_face_key = [](
        const SnakeGapSbm3DUtilities::NodePointerType& pNode0,
        const SnakeGapSbm3DUtilities::NodePointerType& pNode1,
        const SnakeGapSbm3DUtilities::NodePointerType& pNode2)
    {
        KRATOS_ERROR_IF_NOT(pNode0)
            << "[SnakeGapSbmProcess::CreateSbmExtendedGeometries3D] "
            << "Null first node while creating canonical face key.\n";
        KRATOS_ERROR_IF_NOT(pNode1)
            << "[SnakeGapSbmProcess::CreateSbmExtendedGeometries3D] "
            << "Null second node while creating canonical face key.\n";
        KRATOS_ERROR_IF_NOT(pNode2)
            << "[SnakeGapSbmProcess::CreateSbmExtendedGeometries3D] "
            << "Null third node while creating canonical face key.\n";

        std::array<const NodeType*, 3> node_pointers = {{pNode0.get(), pNode1.get(), pNode2.get()}};
        std::sort(
            node_pointers.begin(),
            node_pointers.end(),
            [](const NodeType* pNodeA, const NodeType* pNodeB) {
                return std::less<const NodeType*>{}(pNodeA, pNodeB);
            });

        SnakeGapSbm3DUtilities::CanonicalFaceKey3D key;
        key.NodePointers = node_pointers;
        key.NodeIds = {{node_pointers[0]->Id(), node_pointers[1]->Id(), node_pointers[2]->Id()}};

        return key;
    };

    for (const auto& r_open_face_data : type2_creation_result.OpenFaceDataList) {
        const auto key = make_process_canonical_face_key(
            r_open_face_data.pSurrogateNode,
            r_open_face_data.pProjectionNode0,
            r_open_face_data.pProjectionNode1);

        const auto it = utilities.GetLateralFaceRegistry().find(key);

        KRATOS_ERROR_IF(it == utilities.GetLateralFaceRegistry().end())
            << "Missing type2 face in registry.\n";

        const auto& r_occurrences = it->second;

        std::size_t number_of_type1_occurrences = 0;
        std::size_t number_of_type2_occurrences = 0;
        std::size_t number_of_type3_occurrences = 0;

        for (const auto& r_occurrence : r_occurrences) {
            if (r_occurrence.GapType == 1) {
                ++number_of_type1_occurrences;
            } else if (r_occurrence.GapType == 2) {
                ++number_of_type2_occurrences;
            } else if (r_occurrence.GapType == 3) {
                ++number_of_type3_occurrences;
            }
        }

        KRATOS_ERROR_IF(number_of_type2_occurrences == 0)
            << "OpenFaceDataList face has no type2 occurrence.\n"
            << "  surrogate node = " << r_open_face_data.pSurrogateNode->Id() << "\n"
            << "  spans = "
            << utilities.SpanToString(r_open_face_data.ExternalSpan0)
            << " - "
            << utilities.SpanToString(r_open_face_data.ExternalSpan1) << "\n";

        if (r_occurrences.size() == 1) {
            ++number_of_unclosed_type2_faces;
        } else if (r_occurrences.size() == 2) {
            if (number_of_type3_occurrences == 1) {
                ++number_of_type2_type3_closed_faces;
            } else {
                ++number_of_already_closed_non_type3_faces;
            }

        } else {
            ++number_of_non_manifold_faces;

            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "Non-manifold type2 face at surrogate node "
                << r_open_face_data.pSurrogateNode->Id()
                << " spans "
                << utilities.SpanToString(r_open_face_data.ExternalSpan0)
                << " - "
                << utilities.SpanToString(r_open_face_data.ExternalSpan1)
                << " occurrences = " << r_occurrences.size()
                << ", type1 = " << number_of_type1_occurrences
                << ", type2 = " << number_of_type2_occurrences
                << ", type3 = " << number_of_type3_occurrences << "\n";
        }
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "Type2 open-face closure audit:\n"
        << "  type2-type3 closed faces: "
        << number_of_type2_type3_closed_faces << "\n"
        << "  already closed without type3: "
        << number_of_already_closed_non_type3_faces << "\n"
        << "  remaining open type2 faces: "
        << number_of_unclosed_type2_faces << "\n"
        << "  non-manifold faces: "
        << number_of_non_manifold_faces << "\n";

    KRATOS_ERROR_IF(number_of_non_manifold_faces > 0)
        << "Non-manifold lateral faces detected after type3 creation: "
        << number_of_non_manifold_faces << "\n";

    //debug //FIXME:
    const std::size_t gap_surface_integration_order = mGapInterpolationOrder+1;

    auto open_lateral_surface_data_list = utilities.CreateOpenLateralSurfaceQuadratureData(
                                                    gap_surface_integration_order,
                                                    number_of_shape_functions_derivatives);
                                                    

    KRATOS_ERROR_IF_NOT(mpGapConditionsSubModelPart)
    << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
    << "mpGapConditionsSubModelPart is not initialized.\n";

    ModelPart& r_gap_conditions_model_part = *mpGapConditionsSubModelPart;
    
    std::size_t id_condition = 1;
    
    if (r_gap_conditions_model_part.GetRootModelPart().NumberOfConditions() > 0) {
        id_condition =
            r_gap_conditions_model_part.GetRootModelPart().Conditions().back().Id() + 1;
    }

    struct SkinTriangleProjectionData
    {
        array_1d<double, 3> Point0 = ZeroVector(3);
        array_1d<double, 3> Point1 = ZeroVector(3);
        array_1d<double, 3> Point2 = ZeroVector(3);
        array_1d<double, 3> Center = ZeroVector(3);
        double Radius = 0.0;
    };

    std::vector<SkinTriangleProjectionData> skin_triangle_projection_data;
    PointVector skin_triangle_center_points;
    skin_triangle_projection_data.reserve(rSkinSubModelPart.NumberOfConditions());
    skin_triangle_center_points.reserve(rSkinSubModelPart.NumberOfConditions());

    double max_skin_triangle_radius = 0.0;

    auto append_skin_triangle_projection_data = [&](
        const array_1d<double, 3>& rPoint0,
        const array_1d<double, 3>& rPoint1,
        const array_1d<double, 3>& rPoint2)
    {
        SkinTriangleProjectionData data;
        data.Point0 = rPoint0;
        data.Point1 = rPoint1;
        data.Point2 = rPoint2;
        noalias(data.Center) = (rPoint0 + rPoint1 + rPoint2) / 3.0;

        data.Radius = std::max({
            norm_2(rPoint0 - data.Center),
            norm_2(rPoint1 - data.Center),
            norm_2(rPoint2 - data.Center)});

        const IndexType triangle_id =
            static_cast<IndexType>(skin_triangle_projection_data.size() + 1);

        skin_triangle_projection_data.push_back(data);
        skin_triangle_center_points.push_back(Kratos::make_intrusive<PointType>(
            triangle_id,
            data.Center[0],
            data.Center[1],
            data.Center[2]));

        max_skin_triangle_radius =
            std::max(max_skin_triangle_radius, data.Radius);
    };

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();

        KRATOS_ERROR_IF(r_geometry.PointsNumber() < 3)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "Open lateral projection requires skin conditions with at least "
            << "three points. Condition #" << r_condition.Id()
            << " has " << r_geometry.PointsNumber() << " points.\n";

        for (IndexType triangle_index = 1;
             triangle_index + 1 < r_geometry.PointsNumber();
             ++triangle_index) {
            append_skin_triangle_projection_data(
                r_geometry[0].Coordinates(),
                r_geometry[triangle_index].Coordinates(),
                r_geometry[triangle_index + 1].Coordinates());
        }
    }

    KRATOS_ERROR_IF(skin_triangle_center_points.empty())
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
        << "No skin triangles available for open lateral projection.\n";

    DynamicBins skin_triangle_center_bins(
        skin_triangle_center_points.begin(),
        skin_triangle_center_points.end());

    PointVector candidate_skin_triangle_centers(
        skin_triangle_center_points.size());
    DistanceVector candidate_skin_triangle_center_distances(
        skin_triangle_center_points.size());

    ModelPart& r_skin_root_model_part = rSkinSubModelPart.GetRootModelPart();

    NodePointerContainerType skin_root_node_pointers;
    skin_root_node_pointers.reserve(r_skin_root_model_part.NumberOfNodes());
    for (const auto& r_node : r_skin_root_model_part.Nodes()) {
        skin_root_node_pointers.push_back(
            r_skin_root_model_part.pGetNode(r_node.Id()));
    }

    NodeBinsType skin_root_node_bins(
        skin_root_node_pointers.begin(),
        skin_root_node_pointers.end());

    NodePointerContainerType coincident_skin_root_nodes(8);
    const double open_lateral_projection_node_tolerance =
        1.0e-12 * std::max(1.0, max_knot_span_size);
    const double open_lateral_projection_node_tolerance_squared =
        open_lateral_projection_node_tolerance *
        open_lateral_projection_node_tolerance;

    NodePointerContainerType created_open_lateral_projection_nodes;

    IndexType next_open_lateral_projection_node_id = 1;
    for (const auto& r_node : r_skin_root_model_part.Nodes()) {
        next_open_lateral_projection_node_id =
            std::max(next_open_lateral_projection_node_id, r_node.Id() + 1);
    }

    auto find_closest_skin_projection = [&](
        const Geometry<Node>& rQuadratureGeometry) -> array_1d<double, 3>
    {
        const auto quadrature_center = rQuadratureGeometry.Center();

        array_1d<double, 3> quadrature_center_coordinates = ZeroVector(3);
        quadrature_center_coordinates[0] = quadrature_center.X();
        quadrature_center_coordinates[1] = quadrature_center.Y();
        quadrature_center_coordinates[2] = quadrature_center.Z();

        PointType center_search_point(
            0,
            quadrature_center_coordinates[0],
            quadrature_center_coordinates[1],
            quadrature_center_coordinates[2]);

        const auto p_nearest_center =
            skin_triangle_center_bins.SearchNearestPoint(center_search_point);

        KRATOS_ERROR_IF_NOT(p_nearest_center)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "Could not find a nearest skin triangle center.\n";

        KRATOS_ERROR_IF(p_nearest_center->Id() == 0 ||
                        p_nearest_center->Id() > skin_triangle_projection_data.size())
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "Invalid skin triangle center id " << p_nearest_center->Id()
            << " while projecting an open lateral quadrature geometry.\n";

        array_1d<double, 3> closest_point = ZeroVector(3);
        double closest_distance_squared = std::numeric_limits<double>::max();

        auto update_closest_projection = [&](
            const SkinTriangleProjectionData& rSkinTriangle)
        {
            const array_1d<double, 3> candidate_point =
                ClosestPointOnTriangle(
                    quadrature_center_coordinates,
                    rSkinTriangle.Point0,
                    rSkinTriangle.Point1,
                    rSkinTriangle.Point2);

            const array_1d<double, 3> delta =
                candidate_point - quadrature_center_coordinates;
            const double distance_squared = inner_prod(delta, delta);

            if (distance_squared < closest_distance_squared) {
                closest_distance_squared = distance_squared;
                closest_point = candidate_point;
            }
        };

        update_closest_projection(
            skin_triangle_projection_data[p_nearest_center->Id() - 1]);

        const double search_tolerance =
            1.0e-12 * std::max(1.0, max_knot_span_size);
        const double search_radius =
            std::sqrt(closest_distance_squared) +
            max_skin_triangle_radius +
            search_tolerance;

        const std::size_t number_of_candidate_centers =
            skin_triangle_center_bins.SearchInRadius(
                center_search_point,
                search_radius,
                candidate_skin_triangle_centers.begin(),
                candidate_skin_triangle_center_distances.begin(),
                candidate_skin_triangle_centers.size());

        for (std::size_t i_center = 0;
             i_center < number_of_candidate_centers;
             ++i_center) {
            const auto p_candidate_center =
                candidate_skin_triangle_centers[i_center];

            KRATOS_ERROR_IF_NOT(p_candidate_center)
                << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
                << "Null skin triangle center returned by DynamicBins.\n";

            KRATOS_ERROR_IF(p_candidate_center->Id() == 0 ||
                            p_candidate_center->Id() > skin_triangle_projection_data.size())
                << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
                << "Invalid skin triangle center id "
                << p_candidate_center->Id()
                << " returned by DynamicBins.\n";

            update_closest_projection(
                skin_triangle_projection_data[p_candidate_center->Id() - 1]);
        }

        KRATOS_ERROR_IF_NOT(std::isfinite(closest_distance_squared))
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "Failed to project an open lateral quadrature center onto the skin.\n";

        return closest_point;
    };

    auto find_or_create_open_lateral_projection_node = [&](
        const array_1d<double, 3>& rProjectionPoint) -> NodeType::Pointer
    {
        auto add_to_skin_sub_model_part_if_needed = [&](
            const NodeType::Pointer& pNode) -> NodeType::Pointer
        {
            KRATOS_ERROR_IF_NOT(pNode)
                << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
                << "Trying to reuse a null projected skin node.\n";

            if (!rSkinSubModelPart.HasNode(pNode->Id())) {
                rSkinSubModelPart.AddNode(pNode);
            }

            return pNode;
        };

        for (const auto& p_created_node : created_open_lateral_projection_nodes) {
            const array_1d<double, 3> delta =
                p_created_node->Coordinates() - rProjectionPoint;
            if (inner_prod(delta, delta) <=
                open_lateral_projection_node_tolerance_squared) {
                return add_to_skin_sub_model_part_if_needed(p_created_node);
            }
        }

        PointType projected_point_to_search(
            0,
            rProjectionPoint[0],
            rProjectionPoint[1],
            rProjectionPoint[2]);

        const std::size_t number_of_coincident_root_nodes =
            skin_root_node_bins.SearchInRadius(
                projected_point_to_search,
                open_lateral_projection_node_tolerance,
                coincident_skin_root_nodes.begin(),
                coincident_skin_root_nodes.size());

        for (std::size_t i_node = 0;
             i_node < number_of_coincident_root_nodes;
             ++i_node) {
            const auto p_existing_node = coincident_skin_root_nodes[i_node];
            KRATOS_ERROR_IF_NOT(p_existing_node)
                << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
                << "Null node returned while searching coincident skin nodes.\n";

            const array_1d<double, 3> delta =
                p_existing_node->Coordinates() - rProjectionPoint;
            if (inner_prod(delta, delta) <=
                open_lateral_projection_node_tolerance_squared) {
                return add_to_skin_sub_model_part_if_needed(p_existing_node);
            }
        }

        while (r_skin_root_model_part.HasNode(next_open_lateral_projection_node_id)) {
            ++next_open_lateral_projection_node_id;
        }

        auto p_projection_node = rSkinSubModelPart.CreateNewNode(
            next_open_lateral_projection_node_id++,
            rProjectionPoint[0],
            rProjectionPoint[1],
            rProjectionPoint[2]);

        created_open_lateral_projection_nodes.push_back(p_projection_node);

        return p_projection_node;
    };

    auto set_projection_node = [](
        Geometry<Node>& rQuadratureGeometry,
        const NodeType::Pointer& pProjectionNode)
    {
        KRATOS_ERROR_IF_NOT(pProjectionNode)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "Trying to set a null projection node.\n";

        rQuadratureGeometry.SetValue(PROJECTION_NODE, pProjectionNode);
    };

    std::size_t number_of_open_lateral_projection_nodes = 0;
    double max_open_lateral_projection_distance = 0.0;
    
    for (auto& r_lateral_data : open_lateral_surface_data_list) {
        if (r_lateral_data.SurfaceQuadraturePointGeometries.size() == 0) {
            continue;
        }

        for (auto& r_quadrature_geometry :
             r_lateral_data.SurfaceQuadraturePointGeometries) {
            const array_1d<double, 3> projected_point =
                find_closest_skin_projection(r_quadrature_geometry);

            const auto quadrature_center = r_quadrature_geometry.Center();
            array_1d<double, 3> quadrature_center_coordinates = ZeroVector(3);
            quadrature_center_coordinates[0] = quadrature_center.X();
            quadrature_center_coordinates[1] = quadrature_center.Y();
            quadrature_center_coordinates[2] = quadrature_center.Z();

            max_open_lateral_projection_distance = std::max(
                max_open_lateral_projection_distance,
                norm_2(projected_point - quadrature_center_coordinates));

            const auto p_projection_node =
                find_or_create_open_lateral_projection_node(projected_point);

            set_projection_node(
                r_quadrature_geometry,
                p_projection_node);
            ++number_of_open_lateral_projection_nodes;
        }
    
        this->CreateConditions(
            r_lateral_data.SurfaceQuadraturePointGeometries.ptr_begin(),
            r_lateral_data.SurfaceQuadraturePointGeometries.ptr_end(),
            r_gap_conditions_model_part,
            std::string("GapSbmSolidCondition"),
            id_condition,
            PropertiesPointerType(),
            knot_span_sizes,
            r_lateral_data.NeighbourGeometries,
            r_lateral_data.CharacteristicLength);
    }

    KRATOS_INFO_IF("SnakeGapSbm3DUtilities", mEchoLevel > 2)
        << "Set " << number_of_open_lateral_projection_nodes
        << " PROJECTION_NODE values on open lateral quadrature geometries.\n";

    KRATOS_INFO_IF("SnakeGapSbm3DUtilities", mEchoLevel >= 0)
        << "Maximum open lateral quadrature-to-skin projection distance: "
        << max_open_lateral_projection_distance << "\n";


    // //TODO: interface conditions on lateral surfaces
    auto interface_lateral_surface_data_list =
                                                utilities.CreateInterfaceLateralSurfaceQuadratureData(
                                                    gap_surface_integration_order,
                                                    number_of_shape_functions_derivatives);

    for (auto& r_interface_data : interface_lateral_surface_data_list) {
        this->CreateConditions(
            r_interface_data.SurfaceQuadraturePointGeometries.ptr_begin(),
            r_interface_data.SurfaceQuadraturePointGeometries.ptr_end(),
            *mpGapInterfaceSubModelPart,
            mGapInterfaceConditionName,
            id_condition,
            PropertiesPointerType(),
            knot_span_sizes,
            r_interface_data.NeighbourGeometries,
            r_interface_data.CharacteristicLength);
    }

    // KRATOS_WATCH("a ver")
    // exit(0);

}

} // namespace Kratos
