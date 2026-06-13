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
#include <iomanip>
#include <limits>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

// Project includes
#include "custom_utilities/snake_gap_sbm_3D_utilities.h"
#include "custom_processes/snake_gap_sbm_process.h"
#include "custom_utilities/snake_gap_sbm_3D_utilities.h"

namespace Kratos
{

namespace
{
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

} // unnamed namespace

SnakeGapSbm3DUtilities::SnakeGapSbm3DUtilities(const int EchoLevel)
    : mEchoLevel(EchoLevel)
{
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

    if (nnz_index == static_cast<std::size_t>(-1)) {
        return nullptr;
    }

    const auto& r_cell_data = rSkinBins.CellDataByNnz[nnz_index];

    if (!r_cell_data.Nodes.empty()) {
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

        return p_best_node;
    }

    if (!r_cell_data.Conditions.empty()) {
        bool found_intersection = false;
        array_1d<double, 3> auxiliary_point = ZeroVector(3);

        for (const auto& p_condition : r_cell_data.Conditions) {
            if (FindTriangleBoxIntersectionPoint(
                    p_condition->GetGeometry(),
                    rSpan,
                    rGridInfo,
                    auxiliary_point)) {
                found_intersection = true;
                break;
            }
        }

        if (!found_intersection) {
            return nullptr;
        }

        ClampPointInsideSpan(auxiliary_point, rSpan, rGridInfo);

        const IndexType new_node_id = GetNextAuxiliarySkinNodeId(rSkinSubModelPart);

        auto p_auxiliary_node = rSkinSubModelPart.CreateNewNode(
            new_node_id,
            auxiliary_point[0],
            auxiliary_point[1],
            auxiliary_point[2]);

        return p_auxiliary_node;
    }

    return nullptr;
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
            r_span_data.Type != GapSpanType::Type2) {
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
        } else {
            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "No skin node found in external span " << SpanToString(r_span_data.Key)
                << " of " << GapSpanTypeToString(r_span_data.Type)
                << "with center at " << center
                << ". PROJECTION_NODE_ID remains zero.\n";
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

        // I would ensure Type1 and Type2 only.
        if (r_span_data.Type != GapSpanType::Type1 &&
            r_span_data.Type != GapSpanType::Type2) {
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

    // FIXME: is it necessary??
    // New step: create auxiliary skin nodes where needed.
    EnsureSkinNodesInExternalSpans(
        rSkinSubModelPart,
        external_spans,
        grid_info);

    // Build bins after auxiliary nodes have been created.
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
        }
    }

    KRATOS_INFO("SnakeGapSbm3DUtilities")
        << "External span initialization summary:\n"
        << "  type 1 spans:                 " << number_type_1 << "\n"
        << "  type 2 spans:                 " << number_type_2 << "\n"
        << "  type 3 spans:                 " << number_type_3 << "\n"
        << "  type 1 with projection node:  " << number_type_1_with_projection << "\n"
        << "  type 2 with projection node:  " << number_type_2_with_projection << "\n"
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
    Vector knot_vector = ZeroVector(4);
    knot_vector[0] = 0.0;
    knot_vector[1] = 0.0;
    knot_vector[2] = 1.0;
    knot_vector[3] = 1.0;
    return knot_vector;
}

SnakeGapSbm3DUtilities::CanonicalFaceKey3D
SnakeGapSbm3DUtilities::MakeCanonicalFaceKey3D(
    const IndexType NodeId0,
    const IndexType NodeId1,
    const IndexType NodeId2) const
{
    CanonicalFaceKey3D key{{NodeId0, NodeId1, NodeId2}};
    std::sort(key.NodeIds.begin(), key.NodeIds.end());
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

        if (oriented_det_jacobian < -1.0e-12) {
            KRATOS_WARNING("SnakeGapSbm3DUtilities")
                << "[CreateCoonsVolumeGaussPoints] Inconsistent Coons volume orientation "
                << "at local point " << r_integration_point.Coordinates()
                << ". Oriented determinant = " << oriented_det_jacobian << ".\n";
        }

        r_integration_point.SetWeight(
            r_integration_point.Weight() * std::max(oriented_det_jacobian, 0.0));
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
        pNode0->Id(),
        pNode1->Id(),
        pApexNode->Id());

    LateralFaceOccurrence occurrence;
    occurrence.GapType = 1;
    occurrence.SurrogateConditionId = rFaceData.pCondition->Id();
    occurrence.ExternalSpan = rFaceData.ExternalSpan;
    occurrence.pGeometry = pLateralSurface;
    occurrence.pNeighbourGeometry = pNeighbourGeometry;

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

SnakeGapSbm3DUtilities::Type1CreationResult
SnakeGapSbm3DUtilities::CreateType1GapGeometries(
    ModelPart& rRootModelPart,
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart,
    const ExternalSpanDataMap& rExternalSpans)
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
    
    const std::size_t  number_of_shape_functions_derivatives = 10; //FIXME:

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

        r_gap_type1_debug.AddGeometry(p_volume_geometry);

        // TODO: fix with mGapApproximationOrder
        const std::size_t gap_volume_integration_order = 4;

        GeometriesArrayType volume_quadrature_point_geometries =
        CreateAndTagVolumeQuadraturePointGeometries(
            p_volume_geometry,
            p_neighbour_geometry,
            gap_volume_integration_order,
            number_of_shape_functions_derivatives);

        if (volume_quadrature_point_geometries.size() > 0) {
            Type1VolumeQuadratureData type1_volume_data;

            type1_volume_data.VolumeQuadraturePointGeometries = std::move(volume_quadrature_point_geometries);

            type1_volume_data.NeighbourGeometries.push_back(p_neighbour_geometry);

            type1_volume_data.CharacteristicLength =
                CalculateType1CharacteristicLength(
                    r_face_data,
                    p_apex_node);

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
            }

            auto p_lateral_surface = CreateType1CollapsedLateralCoonsSurface(
                p_node_0,
                p_node_1,
                p_apex_node);

            p_lateral_surface->SetId(next_geometry_id++);
            SetType1GeometryData(
                *p_lateral_surface,
                projection_node_id,
                p_neighbour_geometry,
                "TYPE_1_LATERAL_SURFACE");

            r_gap_type1_debug.AddGeometry(p_lateral_surface);

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

    if (raw_weight_sum > 1.0e-16 && exact_area > 1.0e-16) {
        const double scale_factor = exact_area / raw_weight_sum;

        for (auto& r_integration_point : integration_points) {
            r_integration_point.SetWeight(
                r_integration_point.Weight() * scale_factor);
        }
    } else {
        KRATOS_WARNING("SnakeGapSbm3DUtilities")
            << "[CreateCoonsSurfaceGaussPoints] Surface has nearly zero area. "
            << "raw_weight_sum=" << raw_weight_sum
            << ", exact_area=" << exact_area << ".\n";
    }

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

        if (r_occurrences.size() != 2) {
            continue;
        }

        std::vector<Geometry<Node>::Pointer> neighbour_geometries =
            CollectUniqueNeighbourGeometries(r_occurrences);

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

    ModelPart& r_root_model_part = mpIgaModelPart->GetRootModelPart();
    const bool store_gap_debug_geometries =
        mThisParameters.Has("store_gap_debug_geometries")
            ? mThisParameters["store_gap_debug_geometries"].GetBool()
            : false;

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
    skin_condition_points.reserve(rSkinSubModelPart.NumberOfConditions());

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();

        KRATOS_ERROR_IF(r_geometry.PointsNumber() < 3)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: "
            << "IsPointInsideSkinBoundary3D requires triangular skin conditions. "
            << "Condition #" << r_condition.Id()
            << " has " << r_geometry.PointsNumber() << " points."
            << std::endl;

        skin_condition_points.push_back(Kratos::make_intrusive<PointType>(
            r_condition.Id(),
            r_geometry[0].X(),
            r_geometry[0].Y(),
            r_geometry[0].Z()));
    }

    DynamicBins skin_points_bin(
        skin_condition_points.begin(),
        skin_condition_points.end());
    
    // lambda to identify interior and exterior knot spans based on the position of their centers with respect to the skin boundary 
    auto is_span_active_from_skin_center =
    [&](const SnakeGapSbm3DUtilities::SpanKey3D& rSpan) -> bool
    {
        const double x0 = bounds.MinU + static_cast<double>(rSpan.I) * knot_span_sizes[0];
        const double x1 = x0 + knot_span_sizes[0];

        const double y0 = bounds.MinV + static_cast<double>(rSpan.J) * knot_span_sizes[1];
        const double y1 = y0 + knot_span_sizes[1];

        const double z0 = bounds.MinW + static_cast<double>(rSpan.K) * knot_span_sizes[2];
        const double z1 = z0 + knot_span_sizes[2];

        const double xc = 0.5 * (x0 + x1);
        const double yc = 0.5 * (y0 + y1);
        const double zc = 0.5 * (z0 + z1);

        std::array<PointType, 21> check_points = {
            // center
            PointType(0, xc, yc, zc),

            // corners
            PointType(0, x0, y0, z0),
            PointType(0, x1, y0, z0),
            PointType(0, x0, y1, z0),
            PointType(0, x1, y1, z0),

            PointType(0, x0, y0, z1),
            PointType(0, x1, y0, z1),
            PointType(0, x0, y1, z1),
            PointType(0, x1, y1, z1),

            // edge midpoints, z = z0
            PointType(0, xc, y0, z0),
            PointType(0, xc, y1, z0),
            PointType(0, x0, yc, z0),
            PointType(0, x1, yc, z0),

            // edge midpoints, z = z1
            PointType(0, xc, y0, z1),
            PointType(0, xc, y1, z1),
            PointType(0, x0, yc, z1),
            PointType(0, x1, yc, z1),

            // vertical edge midpoints
            PointType(0, x0, y0, zc),
            PointType(0, x1, y0, zc),
            PointType(0, x0, y1, zc),
            PointType(0, x1, y1, zc)
        };

        for (const auto& r_point : check_points) {
            if (!this->IsPointInsideSkinBoundary3D(
                    r_point,
                    skin_points_bin,
                    rSkinSubModelPart)) {
                return false;
            }
        }

        return true;
    };

    const auto external_spans = utilities.InitializeExternalSpanData(
        rSkinSubModelPart,
        rSurrogateSubModelPart,
        is_span_active_from_skin_center);

    
    auto type1_creation_result = utilities.CreateType1GapGeometries(
        r_root_model_part,
        rSkinSubModelPart,
        rSurrogateSubModelPart,
        external_spans);

    
    IndexType id_element = 1;
    for (const auto& r_element : r_root_model_part.Elements()) {
        id_element = std::max(id_element, r_element.Id() + 1);
    }

    for (auto& r_type1_data : type1_creation_result.VolumeQuadratureDataList) {
        this->CreateElements(
            r_type1_data.VolumeQuadraturePointGeometries.ptr_begin(),
            r_type1_data.VolumeQuadraturePointGeometries.ptr_end(),
            *mpGapElementsSubModelPart,
            std::string("GapSbmSolidElement"),
            id_element,
            PropertiesPointerType(),
            r_type1_data.NeighbourGeometries,
            r_type1_data.CharacteristicLength);
    }


    //debug //FIXME:
    const std::size_t gap_surface_integration_order = 4;

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
    
    for (auto& r_lateral_data : open_lateral_surface_data_list) {
        if (r_lateral_data.SurfaceQuadraturePointGeometries.size() == 0) {
            continue;
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
