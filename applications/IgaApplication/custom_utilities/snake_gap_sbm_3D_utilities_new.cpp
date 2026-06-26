// //    |  /           |
// //    ' /   __| _` | __|  _ \   __|
// //    . \  |   (   | |   (   |\__ `
// //   _|\_\_|  \__,_|\__|\___/ ____/
// //                   Multi-Physics
// //
// //  License:         BSD License
// //                   Kratos default license: kratos/license.txt
// //

// // System includes
// #include <algorithm>
// #include <array>
// #include <cmath>
// #include <iomanip>
// #include <limits>
// #include <map>
// #include <set>
// #include <sstream>
// #include <tuple>
// #include <unordered_map>
// #include <unordered_set>
// #include <utility>

// // Project includes
// #include "custom_utilities/snake_gap_sbm_3D_utilities.h"

// namespace Kratos
// {

// namespace
// {
//     struct ParameterSpaceBounds
//     {
//         double MinU = 0.0;
//         double MaxU = 0.0;
//         double MinV = 0.0;
//         double MaxV = 0.0;
//         double MinW = 0.0;
//         double MaxW = 0.0;
//     };

//     ParameterSpaceBounds ReadParameterSpaceBounds3D(
//         const ModelPart& rModelPart,
//         const char* pCaller)
//     {
//         ParameterSpaceBounds bounds;

//         if (rModelPart.Has(PATCH_PARAMETER_SPACE_CORNERS)) {
//             const Matrix& r_patch_corners = rModelPart.GetValue(PATCH_PARAMETER_SPACE_CORNERS);
//             KRATOS_ERROR_IF(r_patch_corners.size1() < 3 || r_patch_corners.size2() < 2)
//                 << "[" << pCaller << "] PATCH_PARAMETER_SPACE_CORNERS must be at least 3x2.\n";

//             bounds.MinU = r_patch_corners(0, 0);
//             bounds.MaxU = r_patch_corners(0, 1);
//             bounds.MinV = r_patch_corners(1, 0);
//             bounds.MaxV = r_patch_corners(1, 1);
//             bounds.MinW = r_patch_corners(2, 0);
//             bounds.MaxW = r_patch_corners(2, 1);
//             return bounds;
//         }

//         const auto& r_parameter_space_corners = rModelPart.GetValue(PARAMETER_SPACE_CORNERS);
//         KRATOS_ERROR_IF(r_parameter_space_corners.size() < 3)
//             << "[" << pCaller << "] PARAMETER_SPACE_CORNERS must contain at least three vectors.\n";
//         KRATOS_ERROR_IF(r_parameter_space_corners[0].size() < 2 ||
//                         r_parameter_space_corners[1].size() < 2 ||
//                         r_parameter_space_corners[2].size() < 2)
//             << "[" << pCaller << "] PARAMETER_SPACE_CORNERS entries must contain [min, max] values."
//             << " If the model part stores a multipatch box list, use PATCH_PARAMETER_SPACE_CORNERS on the patch model part instead.\n";

//         bounds.MinU = r_parameter_space_corners[0][0];
//         bounds.MaxU = r_parameter_space_corners[0][1];
//         bounds.MinV = r_parameter_space_corners[1][0];
//         bounds.MaxV = r_parameter_space_corners[1][1];
//         bounds.MinW = r_parameter_space_corners[2][0];
//         bounds.MaxW = r_parameter_space_corners[2][1];
//         return bounds;
//     }

//     std::size_t ComputeSpanCount3D(
//         const double domain_length,
//         const double span_size,
//         const char* pDirection,
//         const char* pCaller)
//     {
//         const double spans_real = domain_length / span_size;
//         const double rounded_spans = std::round(spans_real);
//         constexpr double absolute_tolerance = 1.0e-10;
//         constexpr double relative_tolerance = 1.0e-9;
//         const double tolerance = std::max(
//             absolute_tolerance,
//             relative_tolerance * std::max(1.0, std::abs(rounded_spans)));

//         KRATOS_ERROR_IF(rounded_spans <= 0.0)
//             << "[" << pCaller << "] Non-positive number of knot spans in " << pDirection << ".\n";
//         KRATOS_ERROR_IF(std::abs(spans_real - rounded_spans) > tolerance)
//             << "[" << pCaller << "] Non-integer number of knot spans in " << pDirection
//             << " (" << std::setprecision(17) << spans_real << ")"
//             << " computed from domain length " << domain_length
//             << " and span size " << span_size << ".\n";

//         return static_cast<std::size_t>(rounded_spans);
//     }

//     struct SurrogateSegmentKey
//     {
//         IndexType FirstNodeId = 0;
//         IndexType SecondNodeId = 0;

//         SurrogateSegmentKey() = default;

//         SurrogateSegmentKey(const IndexType NodeIdA, const IndexType NodeIdB)
//             : FirstNodeId(std::min(NodeIdA, NodeIdB))
//             , SecondNodeId(std::max(NodeIdA, NodeIdB))
//         {
//         }

//         bool operator==(const SurrogateSegmentKey& rOther) const
//         {
//             return FirstNodeId == rOther.FirstNodeId && SecondNodeId == rOther.SecondNodeId;
//         }
//     };

//     struct SurrogateSegmentKeyHasher
//     {
//         std::size_t operator()(const SurrogateSegmentKey& rKey) const
//         {
//             const std::size_t first_hash = std::hash<IndexType>{}(rKey.FirstNodeId);
//             const std::size_t second_hash = std::hash<IndexType>{}(rKey.SecondNodeId);
//             return first_hash ^ (second_hash + 0x9e3779b9 + (first_hash << 6) + (first_hash >> 2));
//         }
//     };

//     struct SegmentNeighbourRecord
//     {
//         IndexType ConditionId = 0;
//         int Orientation = 1;
//         Geometry<Node>::Pointer pNeighbourGeometry;
//     };

//     struct SurrogateSegmentData
//     {
//         IndexType SegmentId = 0;
//         IndexType ProjectionSurfaceId = 0;
//         Node::Pointer pCanonicalFirstNode;
//         Node::Pointer pCanonicalSecondNode;
//         Geometry<Node>::Pointer pSegmentGeometry;
//         SnakeGapSbmProcess::NurbsSurfaceType::Pointer pProjectionSurface;
//         std::vector<SegmentNeighbourRecord> NeighbourRecords;
//     };

//     struct SurrogateFaceSegmentData
//     {
//         std::array<IndexType, 4> SegmentIds = {0, 0, 0, 0};
//         std::array<int, 4> Orientations = {1, 1, 1, 1};
//     };
// } // unnamed namespace

// SnakeGapSbmProcess::KnotSpanSkinBinsCSR
// SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D(
//     const ModelPart& rSkinSubModelPart,
//     const ModelPart& rSurrogateSubModelPart) const
// {
//     KnotSpanSkinBinsCSR knot_span_data;

//     const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();
//     const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);
//     KRATOS_ERROR_IF(knot_span_sizes.size() < 3)
//         << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] KNOT_SPAN_SIZES must have at least three entries.\n";

//     const auto bounds = ReadParameterSpaceBounds3D(
//         r_parent_model_part,
//         "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D");

//     const double span_size_x = knot_span_sizes[0];
//     const double span_size_y = knot_span_sizes[1];
//     const double span_size_z = knot_span_sizes[2];
//     KRATOS_ERROR_IF(span_size_x <= 0.0 || span_size_y <= 0.0 || span_size_z <= 0.0)
//         << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] Knot span sizes must be positive.\n";

//     const double domain_length_u = bounds.MaxU - bounds.MinU;
//     const double domain_length_v = bounds.MaxV - bounds.MinV;
//     const double domain_length_w = bounds.MaxW - bounds.MinW;
//     KRATOS_ERROR_IF(domain_length_u <= 0.0 || domain_length_v <= 0.0 || domain_length_w <= 0.0)
//         << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] Invalid parameter space extents.\n";

//     const std::size_t number_of_spans_x = ComputeSpanCount3D(
//         domain_length_u,
//         span_size_x,
//         "u",
//         "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D");
//     const std::size_t number_of_spans_y = ComputeSpanCount3D(
//         domain_length_v,
//         span_size_y,
//         "v",
//         "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D");
//     const std::size_t number_of_spans_z = ComputeSpanCount3D(
//         domain_length_w,
//         span_size_z,
//         "w",
//         "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D");

//     knot_span_data.NumberOfSpansX = number_of_spans_x;
//     knot_span_data.NumberOfSpansY = number_of_spans_y;
//     knot_span_data.NumberOfSpansZ = number_of_spans_z;
//     knot_span_data.MinU = bounds.MinU;
//     knot_span_data.MaxU = bounds.MaxU;
//     knot_span_data.MinV = bounds.MinV;
//     knot_span_data.MaxV = bounds.MaxV;
//     knot_span_data.MinW = bounds.MinW;
//     knot_span_data.MaxW = bounds.MaxW;
//     knot_span_data.SpanSizeX = span_size_x;
//     knot_span_data.SpanSizeY = span_size_y;
//     knot_span_data.SpanSizeZ = span_size_z;

//     const std::size_t number_of_flattened_yz_spans = number_of_spans_y * number_of_spans_z;
//     if (rSkinSubModelPart.NumberOfNodes() == 0 && rSkinSubModelPart.NumberOfConditions() == 0) {
//         knot_span_data.Occupancy.resize(number_of_spans_x, number_of_flattened_yz_spans, false);
//         return knot_span_data;
//     }

//     const double min_knot_span_size = std::min({span_size_x, span_size_y, span_size_z});
//     const double tolerance = 1.0e-10;
//     const double boundary_tolerance = 1.0e-10 * min_knot_span_size;
//     const auto compute_span_index = [tolerance](
//         double coordinate,
//         double min_value,
//         double max_value,
//         double span_size,
//         std::size_t span_count) {
//         double clamped = coordinate;
//         if (coordinate < min_value) {
//             KRATOS_ERROR_IF(coordinate < min_value - tolerance)
//                 << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] coordinate below minimum parameter range.\n";
//             clamped = min_value;
//         } else if (coordinate > max_value) {
//             KRATOS_ERROR_IF(coordinate > max_value + tolerance)
//                 << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] coordinate above maximum parameter range.\n";
//             clamped = max_value;
//         }
//         std::size_t span_index = static_cast<std::size_t>(
//             std::floor((clamped - min_value) / span_size + tolerance));
//         if (span_index >= span_count) {
//             span_index = span_count - 1;
//         }
//         return span_index;
//     };

//     const auto compute_node_span_indices = [&](double coordinate,
//                                                double min_value,
//                                                double max_value,
//                                                double span_size,
//                                                std::size_t span_count) {
//         std::vector<std::size_t> span_indices;
//         const std::size_t base_span_index = compute_span_index(
//             coordinate,
//             min_value,
//             max_value,
//             span_size,
//             span_count);
//         span_indices.push_back(base_span_index);

//         const double clamped_coordinate = std::max(min_value, std::min(max_value, coordinate));
//         const double relative_grid_coordinate = (clamped_coordinate - min_value) / span_size;
//         const double nearest_grid_line = std::round(relative_grid_coordinate);
//         const double distance_to_grid_line = std::abs(relative_grid_coordinate - nearest_grid_line) * span_size;

//         if (distance_to_grid_line <= boundary_tolerance) {
//             const int grid_line_index = static_cast<int>(nearest_grid_line);
//             if (grid_line_index > 0) {
//                 span_indices.push_back(static_cast<std::size_t>(grid_line_index - 1));
//             }
//             if (grid_line_index < static_cast<int>(span_count)) {
//                 span_indices.push_back(static_cast<std::size_t>(grid_line_index));
//             }
//         }

//         std::sort(span_indices.begin(), span_indices.end());
//         span_indices.erase(std::unique(span_indices.begin(), span_indices.end()), span_indices.end());
//         return span_indices;
//     };

//     std::vector<std::vector<std::size_t>> column_indices_per_row(number_of_spans_x);
//     auto add_span_to_pattern = [&](const std::size_t SpanIndexX,
//                                    const std::size_t SpanIndexY,
//                                    const std::size_t SpanIndexZ) {
//         column_indices_per_row[SpanIndexX].push_back(
//             FlattenSpanColumnIndex(knot_span_data, SpanIndexY, SpanIndexZ));
//     };

//     for (const auto& r_node : rSkinSubModelPart.Nodes()) {
//         const auto span_indices_x = compute_node_span_indices(r_node.X(), bounds.MinU, bounds.MaxU, span_size_x, number_of_spans_x);
//         const auto span_indices_y = compute_node_span_indices(r_node.Y(), bounds.MinV, bounds.MaxV, span_size_y, number_of_spans_y);
//         const auto span_indices_z = compute_node_span_indices(r_node.Z(), bounds.MinW, bounds.MaxW, span_size_z, number_of_spans_z);
//         for (const std::size_t span_index_x : span_indices_x) {
//             for (const std::size_t span_index_y : span_indices_y) {
//                 for (const std::size_t span_index_z : span_indices_z) {
//                     add_span_to_pattern(span_index_x, span_index_y, span_index_z);
//                 }
//             }
//         }
//     }

//     for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
//         const auto& r_geometry = r_condition.GetGeometry();
//         if (r_geometry.size() == 0) {
//             continue;
//         }
//         array_1d<double, 3> centroid = ZeroVector(3);
//         for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
//             centroid += r_geometry[point_index].Coordinates();
//         }
//         centroid /= static_cast<double>(r_geometry.size());
//         const std::size_t span_index_x = compute_span_index(centroid[0], bounds.MinU, bounds.MaxU, span_size_x, number_of_spans_x);
//         const std::size_t span_index_y = compute_span_index(centroid[1], bounds.MinV, bounds.MaxV, span_size_y, number_of_spans_y);
//         const std::size_t span_index_z = compute_span_index(centroid[2], bounds.MinW, bounds.MaxW, span_size_z, number_of_spans_z);
//         add_span_to_pattern(span_index_x, span_index_y, span_index_z);
//     }

//     std::size_t number_of_non_zero_entries = 0;
//     for (auto& r_column_indices : column_indices_per_row) {
//         std::sort(r_column_indices.begin(), r_column_indices.end());
//         r_column_indices.erase(std::unique(r_column_indices.begin(), r_column_indices.end()), r_column_indices.end());
//         number_of_non_zero_entries += static_cast<std::size_t>(r_column_indices.size());
//     }

//     auto& r_occupancy_matrix = knot_span_data.Occupancy;
//     r_occupancy_matrix.resize(number_of_spans_x, number_of_flattened_yz_spans, false);
//     r_occupancy_matrix.reserve(number_of_non_zero_entries);

//     std::vector<std::unordered_map<std::size_t, std::size_t>> occupancy_index_lookup(number_of_spans_x);
//     for (auto& r_lookup : occupancy_index_lookup) {
//         r_lookup.reserve(8);
//     }

//     std::size_t nnz_counter = 0;
//     for (std::size_t span_index_x = 0; span_index_x < number_of_spans_x; ++span_index_x) {
//         for (const std::size_t span_column_index : column_indices_per_row[span_index_x]) {
//             r_occupancy_matrix.push_back(span_index_x, span_column_index, 0.0);
//             occupancy_index_lookup[span_index_x].emplace(span_column_index, nnz_counter++);
//         }
//     }

//     std::vector<NodePointerContainerType> nodes_per_non_zero(number_of_non_zero_entries);
//     std::vector<ConditionPointerContainerType> conditions_per_non_zero(number_of_non_zero_entries);

//     auto find_non_zero_index = [&](const std::size_t SpanIndexX,
//                                    const std::size_t SpanIndexY,
//                                    const std::size_t SpanIndexZ) {
//         const std::size_t span_column_index = FlattenSpanColumnIndex(knot_span_data, SpanIndexY, SpanIndexZ);
//         const auto nnz_it = occupancy_index_lookup[SpanIndexX].find(span_column_index);
//         return nnz_it != occupancy_index_lookup[SpanIndexX].end() ? nnz_it->second : static_cast<std::size_t>(-1);
//     };

//     for (const auto& r_node : rSkinSubModelPart.Nodes()) {
//         const auto span_indices_x = compute_node_span_indices(r_node.X(), bounds.MinU, bounds.MaxU, span_size_x, number_of_spans_x);
//         const auto span_indices_y = compute_node_span_indices(r_node.Y(), bounds.MinV, bounds.MaxV, span_size_y, number_of_spans_y);
//         const auto span_indices_z = compute_node_span_indices(r_node.Z(), bounds.MinW, bounds.MaxW, span_size_z, number_of_spans_z);
//         for (const std::size_t span_index_x : span_indices_x) {
//             for (const std::size_t span_index_y : span_indices_y) {
//                 for (const std::size_t span_index_z : span_indices_z) {
//                     const std::size_t non_zero_index = find_non_zero_index(span_index_x, span_index_y, span_index_z);
//                     KRATOS_DEBUG_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
//                         << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] node nonzero not found in CSR pattern.\n";

//                     nodes_per_non_zero[non_zero_index].push_back(rSkinSubModelPart.pGetNode(r_node.Id()));
//                     r_occupancy_matrix.value_data()[non_zero_index] += 1.0;
//                 }
//             }
//         }
//     }

//     for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
//         const auto& r_geometry = r_condition.GetGeometry();
//         if (r_geometry.size() == 0) {
//             continue;
//         }
//         array_1d<double, 3> centroid = ZeroVector(3);
//         for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
//             centroid += r_geometry[point_index].Coordinates();
//         }
//         centroid /= static_cast<double>(r_geometry.size());
//         const std::size_t span_index_x = compute_span_index(centroid[0], bounds.MinU, bounds.MaxU, span_size_x, number_of_spans_x);
//         const std::size_t span_index_y = compute_span_index(centroid[1], bounds.MinV, bounds.MaxV, span_size_y, number_of_spans_y);
//         const std::size_t span_index_z = compute_span_index(centroid[2], bounds.MinW, bounds.MaxW, span_size_z, number_of_spans_z);
//         const std::size_t non_zero_index = find_non_zero_index(span_index_x, span_index_y, span_index_z);
//         KRATOS_DEBUG_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
//             << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix3D] condition nonzero not found in CSR pattern.\n";

//         conditions_per_non_zero[non_zero_index].push_back(rSkinSubModelPart.pGetCondition(r_condition.Id()));
//         r_occupancy_matrix.value_data()[non_zero_index] += 1.0;
//     }

//     knot_span_data.CellBinsByNnz.resize(number_of_non_zero_entries);
//     for (std::size_t k = 0; k < number_of_non_zero_entries; ++k) {
//         auto& r_cell_bins = knot_span_data.CellBinsByNnz[k];
//         auto& r_nodes = nodes_per_non_zero[k];
//         if (!r_nodes.empty()) {
//             std::sort(r_nodes.begin(), r_nodes.end(),
//                       [](const Node::Pointer& a, const Node::Pointer& b) {
//                           return a->Id() < b->Id();
//                       });
//             r_nodes.erase(std::unique(r_nodes.begin(), r_nodes.end(),
//                                       [](const Node::Pointer& a, const Node::Pointer& b) {
//                                           return a->Id() == b->Id();
//                                       }),
//                           r_nodes.end());

//             r_cell_bins.Nodes = std::move(r_nodes);
//             r_cell_bins.pNodeBins = std::make_unique<NodeBinsType>(
//                 r_cell_bins.Nodes.begin(), r_cell_bins.Nodes.end());
//             r_cell_bins.HasNodeBins = true;
//         }

//         auto& r_conditions = conditions_per_non_zero[k];
//         if (!r_conditions.empty()) {
//             std::sort(r_conditions.begin(), r_conditions.end(),
//                       [](const ConditionPointerType& a, const ConditionPointerType& b) {
//                           return a->Id() < b->Id();
//                       });
//             r_conditions.erase(std::unique(r_conditions.begin(), r_conditions.end(),
//                                            [](const ConditionPointerType& a, const ConditionPointerType& b) {
//                                                return a->Id() == b->Id();
//                                            }),
//                                r_conditions.end());

//             r_cell_bins.Conditions = std::move(r_conditions);
//             r_cell_bins.ConditionBins = BinsObjectDynamic<ConditionConfigure>(
//                 r_cell_bins.Conditions.begin(), r_cell_bins.Conditions.end());
//             r_cell_bins.HasConditionBins = true;
//         }
//     }

//     return knot_span_data;
// }

// array_1d<double, 3> SnakeGapSbmProcess::CalculateAxisAlignedBrepFaceNormal3D(
//     const CoordinatesArrayType& rVertex00,
//     const CoordinatesArrayType& rVertex01,
//     const CoordinatesArrayType& rVertex10,
//     const CoordinatesArrayType& rVertex11,
//     const NodeType& rFirstSurrogateNode,
//     const double MaxKnotSpanSize) const
// {
//     KRATOS_ERROR_IF(MaxKnotSpanSize <= 0.0)
//         << "::[SnakeGapSbmProcess]::CalculateAxisAlignedBrepFaceNormal3D: MaxKnotSpanSize must be positive."
//         << " Got " << MaxKnotSpanSize << std::endl;

//     const double coordinate_tolerance = 1.0e-12 * MaxKnotSpanSize;
//     IndexType normal_axis = 3;
//     for (IndexType axis = 0; axis < 3; ++axis) {
//         const double min_coordinate = std::min(
//             std::min(rVertex00[axis], rVertex01[axis]),
//             std::min(rVertex10[axis], rVertex11[axis]));
//         const double max_coordinate = std::max(
//             std::max(rVertex00[axis], rVertex01[axis]),
//             std::max(rVertex10[axis], rVertex11[axis]));

//         if (max_coordinate - min_coordinate <= coordinate_tolerance) {
//             KRATOS_ERROR_IF(normal_axis != 3)
//                 << "::[SnakeGapSbmProcess]::CalculateAxisAlignedBrepFaceNormal3D: more than one common coordinate found for BREP face."
//                 << " Vertices: " << rVertex00 << ", " << rVertex01 << ", "
//                 << rVertex10 << ", " << rVertex11
//                 << " | tolerance: " << coordinate_tolerance << std::endl;
//             normal_axis = axis;
//         }
//     }

//     KRATOS_ERROR_IF(normal_axis == 3)
//         << "::[SnakeGapSbmProcess]::CalculateAxisAlignedBrepFaceNormal3D: no common coordinate found for BREP face."
//         << " Vertices: " << rVertex00 << ", " << rVertex01 << ", "
//         << rVertex10 << ", " << rVertex11
//         << " | tolerance: " << coordinate_tolerance << std::endl;

//     array_1d<double, 3> normal = ZeroVector(3);
//     normal[normal_axis] = 1.0;

//     KRATOS_ERROR_IF_NOT(rFirstSurrogateNode.Has(NORMAL))
//         << "::[SnakeGapSbmProcess]::CalculateAxisAlignedBrepFaceNormal3D: NORMAL not set for surrogate node "
//         << rFirstSurrogateNode.Id() << std::endl;

//     return normal;
// }

// SnakeGapSbmProcess::NurbsSurfaceType::Pointer
// SnakeGapSbmProcess::CreateProjectionSurfaceFromSurrogateSegment(
//     const Node::Pointer& pSurrogateNode1,
//     const Node::Pointer& pSurrogateNode2,
//     const ModelPart& rSkinSubModelPart) const
// {
//     KRATOS_ERROR_IF_NOT(pSurrogateNode1)
//         << "::[SnakeGapSbmProcess]::CreateProjectionSurfaceFromSurrogateSegment: first surrogate node is null."
//         << std::endl;
//     KRATOS_ERROR_IF_NOT(pSurrogateNode2)
//         << "::[SnakeGapSbmProcess]::CreateProjectionSurfaceFromSurrogateSegment: second surrogate node is null."
//         << std::endl;
//     KRATOS_ERROR_IF_NOT(pSurrogateNode1->Has(PROJECTION_NODE_ID))
//         << "::[SnakeGapSbmProcess]::CreateProjectionSurfaceFromSurrogateSegment: PROJECTION_NODE_ID not set for surrogate node "
//         << pSurrogateNode1->Id() << "." << std::endl;
//     KRATOS_ERROR_IF_NOT(pSurrogateNode2->Has(PROJECTION_NODE_ID))
//         << "::[SnakeGapSbmProcess]::CreateProjectionSurfaceFromSurrogateSegment: PROJECTION_NODE_ID not set for surrogate node "
//         << pSurrogateNode2->Id() << "." << std::endl;

//     const IndexType projection_id_1 = pSurrogateNode1->GetValue(PROJECTION_NODE_ID);
//     const IndexType projection_id_2 = pSurrogateNode2->GetValue(PROJECTION_NODE_ID);
//     KRATOS_ERROR_IF(projection_id_1 == 0 || !rSkinSubModelPart.HasNode(projection_id_1))
//         << "::[SnakeGapSbmProcess]::CreateProjectionSurfaceFromSurrogateSegment: invalid projection node id "
//         << projection_id_1 << " for surrogate node " << pSurrogateNode1->Id() << "." << std::endl;
//     KRATOS_ERROR_IF(projection_id_2 == 0 || !rSkinSubModelPart.HasNode(projection_id_2))
//         << "::[SnakeGapSbmProcess]::CreateProjectionSurfaceFromSurrogateSegment: invalid projection node id "
//         << projection_id_2 << " for surrogate node " << pSurrogateNode2->Id() << "." << std::endl;

//     const auto p_skin_node_1 = rSkinSubModelPart.pGetNode(projection_id_1);
//     const auto p_skin_node_2 = rSkinSubModelPart.pGetNode(projection_id_2);

//     const array_1d<double, 3> S1 = pSurrogateNode1->Coordinates();
//     const array_1d<double, 3> S2 = pSurrogateNode2->Coordinates();
//     const array_1d<double, 3> K1 = p_skin_node_1->Coordinates();
//     const array_1d<double, 3> K2 = p_skin_node_2->Coordinates();

//     // Corner convention for the Coons patch:
//     // P00 = S1, P10 = S2, P01 = K1, P11 = K2.
//     const auto B0 = [&](const double U) -> array_1d<double, 3> {
//         return (1.0 - U) * S1 + U * S2;
//     };
//     const auto B1 = [&](const double U) -> array_1d<double, 3> {
//         // This boundary is straight for now; it can later be replaced by a reconstructed 3D skin curve.
//         return (1.0 - U) * K1 + U * K2;
//     };
//     const auto L0 = [&](const double V) -> array_1d<double, 3> {
//         return (1.0 - V) * S1 + V * K1;
//     };
//     const auto L1 = [&](const double V) -> array_1d<double, 3> {
//         return (1.0 - V) * S2 + V * K2;
//     };
//     const auto bilinear_blend = [&](const double U, const double V) -> array_1d<double, 3> {
//         return (1.0 - U) * (1.0 - V) * S1
//             + U * (1.0 - V) * S2
//             + (1.0 - U) * V * K1
//             + U * V * K2;
//     };
//     const auto coons_point = [&](const double U, const double V) -> array_1d<double, 3> {
//         return (1.0 - V) * B0(U)
//             + V * B1(U)
//             + (1.0 - U) * L0(V)
//             + U * L1(V)
//             - bilinear_blend(U, V);
//     };

//     const array_1d<double, 3> P00 = coons_point(0.0, 0.0);
//     const array_1d<double, 3> P10 = coons_point(1.0, 0.0);
//     const array_1d<double, 3> P01 = coons_point(0.0, 1.0);
//     const array_1d<double, 3> P11 = coons_point(1.0, 1.0);

//     PointerVector<NodeType> control_points;
//     control_points.push_back(NodeType::Pointer(new NodeType(1, P00)));
//     control_points.push_back(NodeType::Pointer(new NodeType(2, P10)));
//     control_points.push_back(NodeType::Pointer(new NodeType(3, P01)));
//     control_points.push_back(NodeType::Pointer(new NodeType(4, P11)));

//     const std::size_t polynomial_degree_u = 1;
//     const std::size_t polynomial_degree_v = 1;

//     Vector knot_vector_u = ZeroVector(4);
//     knot_vector_u[2] = 1.0;
//     knot_vector_u[3] = 1.0;

//     Vector knot_vector_v = ZeroVector(4);
//     knot_vector_v[2] = 1.0;
//     knot_vector_v[3] = 1.0;

//     return Kratos::make_shared<NurbsSurfaceType>(
//         control_points,
//         polynomial_degree_u,
//         polynomial_degree_v,
//         knot_vector_u,
//         knot_vector_v);
// }

// SnakeGapSbmProcess::NurbsSurfaceType::Pointer
// SnakeGapSbmProcess::CreateGapCoonsSurface(
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
//     const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
//     const bool ReverseSurfaceNormal) const
// {
//     KRATOS_ERROR_IF_NOT(pProjectionSurface00_01)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsSurface: projection surface 00-01 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface00_10)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsSurface: projection surface 00-10 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface01_11)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsSurface: projection surface 01-11 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface10_11)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsSurface: projection surface 10-11 is null." << std::endl;

//     auto create_unit_knot_vector = []() {
//         Vector knot_vector = ZeroVector(4);
//         knot_vector[2] = 1.0;
//         knot_vector[3] = 1.0;
//         return knot_vector;
//     };

//     auto create_surface_from_points = [&](const array_1d<double, 3>& rP00,
//                                           const array_1d<double, 3>& rP10,
//                                           const array_1d<double, 3>& rP01,
//                                           const array_1d<double, 3>& rP11) {
//         PointerVector<NodeType> control_points;
//         control_points.push_back(NodeType::Pointer(new NodeType(1, rP00)));
//         control_points.push_back(NodeType::Pointer(new NodeType(2, rP10)));
//         control_points.push_back(NodeType::Pointer(new NodeType(3, rP01)));
//         control_points.push_back(NodeType::Pointer(new NodeType(4, rP11)));

//         return Kratos::make_shared<NurbsSurfaceType>(
//             control_points,
//             std::size_t(1),
//             std::size_t(1),
//             create_unit_knot_vector(),
//             create_unit_knot_vector());
//     };

//     auto surface_point = [](const NurbsSurfaceType& rSurface,
//                             const double U,
//                             const double V) -> array_1d<double, 3> {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = U;
//         local_coordinates[1] = V;

//         array_1d<double, 3> global_coordinates = ZeroVector(3);
//         rSurface.GlobalCoordinates(global_coordinates, local_coordinates);
//         return global_coordinates;
//     };

//     // The top skin-side boundary curves are extracted from the lateral faces.
//     // Later these lambdas can forward genuinely curved skin-side traces.
//     const auto top_b0 = [&](const double U) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface00_10, U, 1.0);
//     };
//     const auto top_b1 = [&](const double U) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface01_11, U, 1.0);
//     };
//     const auto top_l0 = [&](const double V) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface00_01, V, 1.0);
//     };
//     const auto top_l1 = [&](const double V) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface10_11, V, 1.0);
//     };

//     const array_1d<double, 3> P001 = top_b0(0.0);
//     const array_1d<double, 3> P101 = top_b0(1.0);
//     const array_1d<double, 3> P011 = top_b1(0.0);
//     const array_1d<double, 3> P111 = top_b1(1.0);

//     auto check_matching_top_edge_endpoint = [](const array_1d<double, 3>& rFirstPoint,
//                                                const array_1d<double, 3>& rSecondPoint,
//                                                const char* pEndpointName) {
//         const array_1d<double, 3> difference = rFirstPoint - rSecondPoint;
//         KRATOS_ERROR_IF(norm_2(difference) > 1.0e-10)
//             << "::[SnakeGapSbmProcess]::CreateGapCoonsSurface: inconsistent top skin boundary endpoint "
//             << pEndpointName << ". First point: " << rFirstPoint
//             << " second point: " << rSecondPoint << "." << std::endl;
//     };

//     check_matching_top_edge_endpoint(P001, top_l0(0.0), "P001");
//     check_matching_top_edge_endpoint(P011, top_l0(1.0), "P011");
//     check_matching_top_edge_endpoint(P101, top_l1(0.0), "P101");
//     check_matching_top_edge_endpoint(P111, top_l1(1.0), "P111");

//     const auto top_coons_point = [&](const double U, const double V) -> array_1d<double, 3> {
//         const array_1d<double, 3> bilinear_blend =
//             (1.0 - U) * (1.0 - V) * P001
//             + U * (1.0 - V) * P101
//             + (1.0 - U) * V * P011
//             + U * V * P111;

//         return (1.0 - V) * top_b0(U)
//             + V * top_b1(U)
//             + (1.0 - U) * top_l0(V)
//             + U * top_l1(V)
//             - bilinear_blend;
//     };

//     const array_1d<double, 3> top_point_00 = top_coons_point(0.0, 0.0);
//     const array_1d<double, 3> top_point_10 = top_coons_point(1.0, 0.0);
//     const array_1d<double, 3> top_point_01 = top_coons_point(0.0, 1.0);
//     const array_1d<double, 3> top_point_11 = top_coons_point(1.0, 1.0);

//     if (ReverseSurfaceNormal) {
//         return create_surface_from_points(
//             top_point_00,
//             top_point_01,
//             top_point_10,
//             top_point_11);
//     }

//     return create_surface_from_points(
//         top_point_00,
//         top_point_10,
//         top_point_01,
//         top_point_11);
// }

// SnakeGapSbmProcess::NurbsVolumeType::Pointer
// SnakeGapSbmProcess::CreateGapCoonsVolume(
//     const Node::Pointer& pNode00,
//     const Node::Pointer& pNode01,
//     const Node::Pointer& pNode10,
//     const Node::Pointer& pNode11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
//     const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
//     const NurbsSurfaceType::Pointer& pTopSkinFace) const
// {
//     KRATOS_ERROR_IF_NOT(pNode00)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: node 00 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pNode01)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: node 01 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pNode10)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: node 10 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pNode11)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: node 11 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface00_01)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: projection surface 00-01 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface00_10)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: projection surface 00-10 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface01_11)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: projection surface 01-11 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pProjectionSurface10_11)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: projection surface 10-11 is null." << std::endl;
//     KRATOS_ERROR_IF_NOT(pTopSkinFace)
//         << "::[SnakeGapSbmProcess]::CreateGapCoonsVolume: top skin face is null." << std::endl;

//     auto create_unit_knot_vector = []() {
//         Vector knot_vector = ZeroVector(4);
//         knot_vector[2] = 1.0;
//         knot_vector[3] = 1.0;
//         return knot_vector;
//     };

//     PointerVector<NodeType> bottom_control_points;
//     bottom_control_points.push_back(pNode00);
//     bottom_control_points.push_back(pNode10);
//     bottom_control_points.push_back(pNode01);
//     bottom_control_points.push_back(pNode11);

//     auto p_bottom_surrogate_face = Kratos::make_shared<NurbsSurfaceType>(
//         bottom_control_points,
//         std::size_t(1),
//         std::size_t(1),
//         create_unit_knot_vector(),
//         create_unit_knot_vector());

//     auto surface_point = [](const NurbsSurfaceType& rSurface,
//                             const double U,
//                             const double V) -> array_1d<double, 3> {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = U;
//         local_coordinates[1] = V;

//         array_1d<double, 3> global_coordinates = ZeroVector(3);
//         rSurface.GlobalCoordinates(global_coordinates, local_coordinates);
//         return global_coordinates;
//     };

//     const auto F_u0 = [&](const double V, const double W) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface00_01, V, W);
//     };
//     const auto F_u1 = [&](const double V, const double W) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface10_11, V, W);
//     };
//     const auto F_v0 = [&](const double U, const double W) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface00_10, U, W);
//     };
//     const auto F_v1 = [&](const double U, const double W) -> array_1d<double, 3> {
//         return surface_point(*pProjectionSurface01_11, U, W);
//     };
//     const auto F_w0 = [&](const double U, const double V) -> array_1d<double, 3> {
//         return surface_point(*p_bottom_surrogate_face, U, V);
//     };
//     const auto F_w1 = [&](const double U, const double V) -> array_1d<double, 3> {
//         return surface_point(*pTopSkinFace, U, V);
//     };

//     const array_1d<double, 3> P000 = F_w0(0.0, 0.0);
//     const array_1d<double, 3> P100 = F_w0(1.0, 0.0);
//     const array_1d<double, 3> P010 = F_w0(0.0, 1.0);
//     const array_1d<double, 3> P110 = F_w0(1.0, 1.0);
//     const array_1d<double, 3> P001 = F_w1(0.0, 0.0);
//     const array_1d<double, 3> P101 = F_w1(1.0, 0.0);
//     const array_1d<double, 3> P011 = F_w1(0.0, 1.0);
//     const array_1d<double, 3> P111 = F_w1(1.0, 1.0);

//     const auto coons_volume_point = [&](const double U,
//                                         const double V,
//                                         const double W) -> array_1d<double, 3> {
//         const array_1d<double, 3> face_blend =
//             (1.0 - U) * F_u0(V, W) + U * F_u1(V, W)
//             + (1.0 - V) * F_v0(U, W) + V * F_v1(U, W)
//             + (1.0 - W) * F_w0(U, V) + W * F_w1(U, V);

//         const array_1d<double, 3> uv_edge_blend =
//             (1.0 - U) * (1.0 - V) * F_u0(0.0, W)
//             + U * (1.0 - V) * F_u1(0.0, W)
//             + (1.0 - U) * V * F_u0(1.0, W)
//             + U * V * F_u1(1.0, W);

//         const array_1d<double, 3> uw_edge_blend =
//             (1.0 - U) * (1.0 - W) * F_w0(0.0, V)
//             + U * (1.0 - W) * F_w0(1.0, V)
//             + (1.0 - U) * W * F_w1(0.0, V)
//             + U * W * F_w1(1.0, V);

//         const array_1d<double, 3> vw_edge_blend =
//             (1.0 - V) * (1.0 - W) * F_w0(U, 0.0)
//             + V * (1.0 - W) * F_w0(U, 1.0)
//             + (1.0 - V) * W * F_w1(U, 0.0)
//             + V * W * F_w1(U, 1.0);

//         const array_1d<double, 3> vertex_blend =
//             (1.0 - U) * (1.0 - V) * (1.0 - W) * P000
//             + U * (1.0 - V) * (1.0 - W) * P100
//             + (1.0 - U) * V * (1.0 - W) * P010
//             + U * V * (1.0 - W) * P110
//             + (1.0 - U) * (1.0 - V) * W * P001
//             + U * (1.0 - V) * W * P101
//             + (1.0 - U) * V * W * P011
//             + U * V * W * P111;

//         return face_blend - uv_edge_blend - uw_edge_blend - vw_edge_blend + vertex_blend;
//     };

//     PointerVector<NodeType> volume_control_points;
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(1, coons_volume_point(0.0, 0.0, 0.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(2, coons_volume_point(1.0, 0.0, 0.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(3, coons_volume_point(0.0, 1.0, 0.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(4, coons_volume_point(1.0, 1.0, 0.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(5, coons_volume_point(0.0, 0.0, 1.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(6, coons_volume_point(1.0, 0.0, 1.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(7, coons_volume_point(0.0, 1.0, 1.0))));
//     volume_control_points.push_back(NodeType::Pointer(new NodeType(8, coons_volume_point(1.0, 1.0, 1.0))));

//     return Kratos::make_shared<NurbsVolumeType>(
//         volume_control_points,
//         std::size_t(1),
//         std::size_t(1),
//         std::size_t(1),
//         create_unit_knot_vector(),
//         create_unit_knot_vector(),
//         create_unit_knot_vector());
// }

// SnakeGapSbmProcess::IntegrationPointsArrayType
// SnakeGapSbmProcess::CreateCoonsVolumeGaussPoints(
//     const std::size_t Order,
//     const NurbsVolumeType& rGapVolume) const
// {
//     KRATOS_ERROR_IF(Order == 0)
//         << "::[SnakeGapSbmProcess]::CreateCoonsVolumeGaussPoints: integration order must be positive."
//         << std::endl;

//     IntegrationPointsArrayType integration_points(Order * Order * Order);
//     auto integration_point_it = integration_points.begin();
//     IntegrationPointUtilities::IntegrationPoints3D(
//         integration_point_it,
//         Order,
//         Order,
//         Order,
//         0.0,
//         1.0,
//         0.0,
//         1.0,
//         0.0,
//         1.0);

//     constexpr double determinant_tolerance = 1.0e-14;
//     std::vector<double> determinant_values;
//     determinant_values.reserve(integration_points.size());

//     std::size_t number_of_positive_det = 0;
//     std::size_t number_of_negative_det = 0;
//     std::size_t number_of_almost_zero_det = 0;

//     for (auto& r_integration_point : integration_points) {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = r_integration_point[0];
//         local_coordinates[1] = r_integration_point[1];
//         local_coordinates[2] = r_integration_point[2];

//         Matrix jacobian;
//         rGapVolume.Jacobian(jacobian, local_coordinates);
//         const double det_jacobian = MathUtils<double>::Det(jacobian);
//         determinant_values.push_back(det_jacobian);

//         if (det_jacobian > determinant_tolerance) {
//             ++number_of_positive_det;
//         } else if (det_jacobian < -determinant_tolerance) {
//             ++number_of_negative_det;
//         } else {
//             ++number_of_almost_zero_det;
//         }
//     }

//     CoordinatesArrayType center_local_coordinates = ZeroVector(3);
//     center_local_coordinates[0] = 0.5;
//     center_local_coordinates[1] = 0.5;
//     center_local_coordinates[2] = 0.5;

//     Matrix center_jacobian;
//     rGapVolume.Jacobian(center_jacobian, center_local_coordinates);
//     const double center_det_jacobian = MathUtils<double>::Det(center_jacobian);

//     // A globally negative determinant is a consistent parametric orientation flip.
//     // Mixed signs remain an error because they indicate a folded volume.
//     double orientation_sign = center_det_jacobian < 0.0 ? -1.0 : 1.0;
//     if (number_of_negative_det > 0 && number_of_positive_det == 0) {
//         orientation_sign = -1.0;
//     } else if (number_of_positive_det > 0 && number_of_negative_det == 0) {
//         orientation_sign = 1.0;
//     }

//     for (std::size_t point_index = 0; point_index < integration_points.size(); ++point_index) {
//         auto& r_integration_point = integration_points[point_index];
//         const double oriented_det_jacobian = orientation_sign * determinant_values[point_index];
//         if (oriented_det_jacobian < -1.0e-12) {
//             KRATOS_WARNING("CreateCoonsVolumeGaussPoints")
//                 << "::[SnakeGapSbmProcess]::CreateCoonsVolumeGaussPoints: inconsistent Coons volume orientation at local point "
//                 << r_integration_point.Coordinates() << ". Oriented determinant: " << oriented_det_jacobian
//                 << ". This is a local folding or inconsistent boundary orientation, not a pure global orientation flip."
//                 << std::endl;
//         }

//         r_integration_point.SetWeight(
//             r_integration_point.Weight() * std::max(oriented_det_jacobian, 0.0));
//     }

//     return integration_points;
// }

// void SnakeGapSbmProcess::CreateGapElementsAndConditions(
//     const Node::Pointer& pNode00,
//     const Node::Pointer& pNode01,
//     const Node::Pointer& pNode10,
//     const Node::Pointer& pNode11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
//     const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
//     const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
//     const IntegrationParameters& rIntegrationParameters,
//     ModelPart& rSkinSubModelPart) const
// {
//     auto p_top_surface_canonical = CreateGapCoonsSurface(
//         pProjectionSurface00_01,
//         pProjectionSurface00_10,
//         pProjectionSurface01_11,
//         pProjectionSurface10_11,
//         false);

//     KRATOS_ERROR_IF_NOT(p_top_surface_canonical)
//         << "::[SnakeGapSbmProcess]:: Failed to create the canonical 3D gap top surface."
//         << std::endl;

//     auto p_top_surface_reversed = CreateGapCoonsSurface(
//         pProjectionSurface00_01,
//         pProjectionSurface00_10,
//         pProjectionSurface01_11,
//         pProjectionSurface10_11,
//         true);

//     KRATOS_ERROR_IF_NOT(p_top_surface_reversed)
//         << "::[SnakeGapSbmProcess]:: Failed to create the reversed 3D gap top surface."
//         << std::endl;

//     auto p_gap_coons_volume = CreateGapCoonsVolume(
//         pNode00,
//         pNode01,
//         pNode10,
//         pNode11,
//         pProjectionSurface00_01,
//         pProjectionSurface00_10,
//         pProjectionSurface01_11,
//         pProjectionSurface10_11,
//         p_top_surface_canonical);

//     KRATOS_ERROR_IF_NOT(p_gap_coons_volume)
//         << "::[SnakeGapSbmProcess]:: Failed to create a 3D gap Coons volume."
//         << std::endl;

//     auto surface_point = [](const NurbsSurfaceType& rSurface,
//                             const double U,
//                             const double V) {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = U;
//         local_coordinates[1] = V;

//         array_1d<double, 3> global_coordinates = ZeroVector(3);
//         rSurface.GlobalCoordinates(global_coordinates, local_coordinates);
//         return global_coordinates;
//     };

//     auto volume_point = [](const NurbsVolumeType& rVolume,
//                            const double U,
//                            const double V,
//                            const double W) {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = U;
//         local_coordinates[1] = V;
//         local_coordinates[2] = W;

//         array_1d<double, 3> global_coordinates = ZeroVector(3);
//         rVolume.GlobalCoordinates(global_coordinates, local_coordinates);
//         return global_coordinates;
//     };

//     auto surface_normal = [](const NurbsSurfaceType& rSurface,
//                              const double U,
//                              const double V) {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = U;
//         local_coordinates[1] = V;

//         array_1d<double, 3> normal = rSurface.Normal(local_coordinates);
//         const double normal_norm = norm_2(normal);
//         if (normal_norm > 1.0e-16) {
//             normal /= normal_norm;
//         }
//         return normal;
//     };

//     struct TopSurfaceOrientationInfo
//     {
//         double Score = 0.0;
//         double MinDot = std::numeric_limits<double>::max();
//         double MaxDot = -std::numeric_limits<double>::max();
//         std::size_t PositiveCount = 0;
//         std::size_t NegativeCount = 0;
//         std::size_t AlmostZeroCount = 0;
//         std::size_t DegenerateCount = 0;
//         std::size_t SampleCount = 0;
//         bool MixedSigns = false;
//         bool GloballyInverted = false;
//         bool HasDegenerateSamples = false;
//     };

//     auto check_top_surface_orientation = [&](const NurbsSurfaceType& rSurface,
//                                              const NurbsVolumeType& rVolume) {
//         TopSurfaceOrientationInfo info;
//         const std::array<double, 3> sample_coordinates = {
//             0.1127016653792583,
//             0.5,
//             0.8872983346207417};
//         const array_1d<double, 3> volume_center = volume_point(rVolume, 0.5, 0.5, 0.5);

//         for (const double u : sample_coordinates) {
//             for (const double v : sample_coordinates) {
//                 ++info.SampleCount;

//                 const array_1d<double, 3> top_point = surface_point(rSurface, u, v);
//                 array_1d<double, 3> outward_reference = top_point - volume_center;
//                 const double reference_norm = norm_2(outward_reference);
//                 array_1d<double, 3> normal = surface_normal(rSurface, u, v);
//                 const double normal_norm = norm_2(normal);

//                 if (reference_norm <= 1.0e-16 || normal_norm <= 1.0e-16) {
//                     ++info.DegenerateCount;
//                     continue;
//                 }

//                 outward_reference /= reference_norm;
//                 normal /= normal_norm;

//                 const double normal_dot_outward = inner_prod(normal, outward_reference);
//                 info.Score += normal_dot_outward;
//                 info.MinDot = std::min(info.MinDot, normal_dot_outward);
//                 info.MaxDot = std::max(info.MaxDot, normal_dot_outward);

//                 constexpr double orientation_tolerance = 1.0e-12;
//                 if (normal_dot_outward > orientation_tolerance) {
//                     ++info.PositiveCount;
//                 } else if (normal_dot_outward < -orientation_tolerance) {
//                     ++info.NegativeCount;
//                 } else {
//                     ++info.AlmostZeroCount;
//                 }
//             }
//         }

//         if (info.PositiveCount + info.NegativeCount + info.AlmostZeroCount == 0) {
//             info.MinDot = 0.0;
//             info.MaxDot = 0.0;
//         }

//         info.MixedSigns = info.PositiveCount > 0 && info.NegativeCount > 0;
//         info.GloballyInverted =
//             info.NegativeCount > 0 && info.PositiveCount == 0 && info.AlmostZeroCount == 0;
//         info.HasDegenerateSamples = info.DegenerateCount > 0 || info.AlmostZeroCount > 0;
//         return info;
//     };

//     const TopSurfaceOrientationInfo canonical_orientation =
//         check_top_surface_orientation(*p_top_surface_canonical, *p_gap_coons_volume);
//     const TopSurfaceOrientationInfo reversed_orientation =
//         check_top_surface_orientation(*p_top_surface_reversed, *p_gap_coons_volume);
//     const double canonical_score = canonical_orientation.Score;
//     const double reversed_score = reversed_orientation.Score;
//     const double selected_score = std::max(canonical_score, reversed_score);

//     const bool use_reversed_top_surface = reversed_score > canonical_score;
//     NurbsSurfaceType::Pointer p_gap_coons_surface =
//         use_reversed_top_surface ? p_top_surface_reversed : p_top_surface_canonical;
//     const TopSurfaceOrientationInfo& r_selected_orientation =
//         use_reversed_top_surface ? reversed_orientation : canonical_orientation;

//     KRATOS_WARNING_IF("CreateGapElementsAndConditions", selected_score < -1.0e-12)
//         << "Selected top skin surface still appears to point inward."
//         << " canonical_score=" << canonical_score
//         << " reversed_score=" << reversed_score
//         << std::endl;
//     KRATOS_WARNING_IF("CreateGapElementsAndConditions", r_selected_orientation.MixedSigns)
//         << "Selected top skin surface has mixed normal orientation samples and may be folded."
//         << " selected=" << (use_reversed_top_surface ? "reversed" : "canonical")
//         << " score=" << r_selected_orientation.Score
//         << " min_dot=" << r_selected_orientation.MinDot
//         << " max_dot=" << r_selected_orientation.MaxDot
//         << " positive=" << r_selected_orientation.PositiveCount
//         << " negative=" << r_selected_orientation.NegativeCount
//         << " almost_zero=" << r_selected_orientation.AlmostZeroCount
//         << " degenerate=" << r_selected_orientation.DegenerateCount
//         << std::endl;
//     KRATOS_WARNING_IF("CreateGapElementsAndConditions", r_selected_orientation.HasDegenerateSamples)
//         << "Selected top skin surface has degenerate or nearly tangent orientation samples."
//         << " selected=" << (use_reversed_top_surface ? "reversed" : "canonical")
//         << " score=" << r_selected_orientation.Score
//         << " min_dot=" << r_selected_orientation.MinDot
//         << " max_dot=" << r_selected_orientation.MaxDot
//         << " positive=" << r_selected_orientation.PositiveCount
//         << " negative=" << r_selected_orientation.NegativeCount
//         << " almost_zero=" << r_selected_orientation.AlmostZeroCount
//         << " degenerate=" << r_selected_orientation.DegenerateCount
//         << std::endl;

//     KRATOS_INFO_IF("CreateGapElementsAndConditions", mEchoLevel > 2)
//         << "Top skin surface orientation diagnostics:"
//         << " canonical score=" << canonical_orientation.Score
//         << " min_dot=" << canonical_orientation.MinDot
//         << " max_dot=" << canonical_orientation.MaxDot
//         << " positive=" << canonical_orientation.PositiveCount
//         << " negative=" << canonical_orientation.NegativeCount
//         << " almost_zero=" << canonical_orientation.AlmostZeroCount
//         << " degenerate=" << canonical_orientation.DegenerateCount
//         << " | reversed score=" << reversed_orientation.Score
//         << " min_dot=" << reversed_orientation.MinDot
//         << " max_dot=" << reversed_orientation.MaxDot
//         << " positive=" << reversed_orientation.PositiveCount
//         << " negative=" << reversed_orientation.NegativeCount
//         << " almost_zero=" << reversed_orientation.AlmostZeroCount
//         << " degenerate=" << reversed_orientation.DegenerateCount
//         << " | selected=" << (use_reversed_top_surface ? "reversed" : "canonical")
//         << std::endl;

//     const std::size_t surface_integration_order = static_cast<std::size_t>(2 * mGapInterpolationOrder + 1);
//     IntegrationPointsArrayType surface_integration_points;
//     GeometriesArrayType surface_quadrature_point_list;

//     p_gap_coons_surface->CreateIntegrationPoints(
//         surface_integration_points,
//         surface_integration_order,
//         surface_integration_order);

//     for (auto& r_integration_point : surface_integration_points) {
//         const double determinant_jacobian = p_gap_coons_surface->DeterminantOfJacobian(r_integration_point);
//         r_integration_point.SetWeight(r_integration_point.Weight() * determinant_jacobian);
//     }

//     IntegrationInfo surface_integration_info = p_gap_coons_surface->GetDefaultIntegrationInfo();
//     surface_integration_info.SetNumberOfIntegrationPointsPerSpan(0, 2*surface_integration_order+1);
//     surface_integration_info.SetNumberOfIntegrationPointsPerSpan(1, 2*surface_integration_order+1);

//     p_gap_coons_surface->CreateQuadraturePointGeometries(
//         surface_quadrature_point_list,
//         rIntegrationParameters.NumberOfShapeFunctionsDerivatives,
//         surface_integration_points,
//         surface_integration_info);

//     if (surface_quadrature_point_list.size() > 0) {
//         KRATOS_ERROR_IF_NOT(mpGapConditionsSubModelPart)
//             << "::[SnakeGapSbmProcess]::CreateGapElementsAndConditions: mpGapConditionsSubModelPart is not initialized."
//             << std::endl;
//         ModelPart& r_gap_conditions_model_part = *mpGapConditionsSubModelPart;

//         std::vector<Geometry<Node>::Pointer> neighbour_geometries_cond{
//             rSurrogateBrepMiddleGeometry};

//         PropertiesPointerType p_props = r_gap_conditions_model_part.pGetProperties(0);
//         KRATOS_ERROR_IF_NOT(p_props)
//             << "::[SnakeGapSbmProcess]::CreateGapElementsAndConditions: failed to get properties 0 for gap conditions."
//             << std::endl;

//         std::size_t id_cond = 1;
//         if (r_gap_conditions_model_part.GetRootModelPart().Conditions().size() > 0) {
//             id_cond = r_gap_conditions_model_part.GetRootModelPart().Conditions().back().Id() + 1;
//         }

//         const double characteristic_condition_length = norm_2(
//             pNode00->Coordinates() - pNode11->Coordinates()) / 2.0;

//         this->CreateConditions(
//             surface_quadrature_point_list.ptr_begin(),
//             surface_quadrature_point_list.ptr_end(),
//             r_gap_conditions_model_part,
//             mGapConditionName,
//             id_cond,
//             p_props,
//             rIntegrationParameters.rKnotSpanSizes,
//             neighbour_geometries_cond,
//             characteristic_condition_length);
//     }

//     auto get_projected_coordinates = [&](const Node::Pointer& p_surrogate_node) {
//         const IndexType projection_id = p_surrogate_node->GetValue(PROJECTION_NODE_ID);
//         return rSkinSubModelPart.pGetNode(projection_id)->Coordinates();
//     };

//     const array_1d<double,3> projected_point_00 = get_projected_coordinates(pNode00);
//     const array_1d<double,3> projected_point_01 = get_projected_coordinates(pNode01);
//     const array_1d<double,3> projected_point_10 = get_projected_coordinates(pNode10);
//     const array_1d<double,3> projected_point_11 = get_projected_coordinates(pNode11);

//     const std::size_t integration_order = mGapInterpolationOrder + 1;
//     IntegrationPointsArrayType diagnostic_volume_integration_points(
//         integration_order * integration_order * integration_order);
//     auto diagnostic_integration_point_it = diagnostic_volume_integration_points.begin();
//     IntegrationPointUtilities::IntegrationPoints3D(
//         diagnostic_integration_point_it,
//         integration_order,
//         integration_order,
//         integration_order,
//         0.0,
//         1.0,
//         0.0,
//         1.0,
//         0.0,
//         1.0);

//     const bool is_folded_gap_volume = CheckGapVolumeOrientation(
//         *p_gap_coons_volume,
//         *pNode00,
//         *pNode01,
//         *pNode10,
//         *pNode11,
//         projected_point_00,
//         projected_point_01,
//         projected_point_10,
//         projected_point_11,
//         pProjectionSurface00_01,
//         pProjectionSurface00_10,
//         pProjectionSurface01_11,
//         pProjectionSurface10_11,
//         diagnostic_volume_integration_points);

//     IntegrationPointsArrayType volume_integration_points =
//         CreateCoonsVolumeGaussPoints(integration_order, *p_gap_coons_volume);

//     double volume_integration_weight_sum = 0.0;
//     for (const auto& r_integration_point : volume_integration_points) {
//         volume_integration_weight_sum += r_integration_point.Weight();
//     }

//     if (volume_integration_weight_sum < 1.0e-14) {
//         return;
//     }

//     if (is_folded_gap_volume) {
//         KRATOS_ERROR
//             << "Folded 3D gap element for surrogate nodes "
//             << pNode00->Coordinates() << ", "
//             << pNode01->Coordinates() << ", "
//             << pNode10->Coordinates() << ", "
//             << pNode11->Coordinates()
//             << " with projections "
//             << projected_point_00 << ", "
//             << projected_point_01 << ", "
//             << projected_point_10 << ", "
//             << projected_point_11
//             << " and volume weight sum " << volume_integration_weight_sum
//             << "." << std::endl;
//     }

//     GeometriesArrayType volume_quadrature_point_list;
//     IntegrationInfo volume_integration_info(
//         {integration_order, integration_order, integration_order},
//         {IntegrationInfo::QuadratureMethod::CUSTOM,
//          IntegrationInfo::QuadratureMethod::CUSTOM,
//          IntegrationInfo::QuadratureMethod::CUSTOM});

//     p_gap_coons_volume->CreateQuadraturePointGeometries(
//         volume_quadrature_point_list,
//         rIntegrationParameters.NumberOfShapeFunctionsDerivatives,
//         volume_integration_points,
//         volume_integration_info);

//     if (volume_quadrature_point_list.size() == 0) {
//         return;
//     }

//     const double characteristic_length = CalculateGapVolumeElementCharacteristicLength(
//         *pNode00,
//         *pNode01,
//         *pNode10,
//         *pNode11,
//         projected_point_00,
//         projected_point_01,
//         projected_point_10,
//         projected_point_11);
//     KRATOS_ERROR_IF(characteristic_length <= 0.0)
//         << "Characteristic length for the 3D gap element must be positive." << std::endl;

//     IndexType id_element = 1;
//     if (mpGapElementsSubModelPart->GetRootModelPart().Elements().size() > 0) {
//         id_element = mpGapElementsSubModelPart->GetRootModelPart().Elements().back().Id() + 1;
//     }

//     std::vector<Geometry<Node>::Pointer> neighbour_geometries_gap_volume{
//         rSurrogateBrepMiddleGeometry};

//     this->CreateElements(
//         volume_quadrature_point_list.ptr_begin(),
//         volume_quadrature_point_list.ptr_end(),
//         *mpGapElementsSubModelPart,
//         std::string("GapSbmSolidElement"),
//         id_element,
//         PropertiesPointerType(),
//         neighbour_geometries_gap_volume,
//         characteristic_length);
// }

// bool SnakeGapSbmProcess::CheckGapVolumeOrientation(
//     const NurbsVolumeType& rGapVolume,
//     const NodeType& rNode00,
//     const NodeType& rNode01,
//     const NodeType& rNode10,
//     const NodeType& rNode11,
//     const array_1d<double,3>& rProjectedPoint00,
//     const array_1d<double,3>& rProjectedPoint01,
//     const array_1d<double,3>& rProjectedPoint10,
//     const array_1d<double,3>& rProjectedPoint11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_01,
//     const NurbsSurfaceType::Pointer& pProjectionSurface00_10,
//     const NurbsSurfaceType::Pointer& pProjectionSurface01_11,
//     const NurbsSurfaceType::Pointer& pProjectionSurface10_11,
//     const IntegrationPointsArrayType& rIntegrationPoints) const
// {
//     constexpr double determinant_tolerance = 1.0e-14;

//     double min_det_jacobian = std::numeric_limits<double>::max();
//     double max_det_jacobian = -std::numeric_limits<double>::max();
//     std::size_t number_of_positive_det = 0;
//     std::size_t number_of_negative_det = 0;
//     std::size_t number_of_almost_zero_det = 0;

//     CoordinatesArrayType center_local_coordinates = ZeroVector(3);
//     center_local_coordinates[0] = 0.5;
//     center_local_coordinates[1] = 0.5;
//     center_local_coordinates[2] = 0.5;

//     Matrix center_jacobian;
//     rGapVolume.Jacobian(center_jacobian, center_local_coordinates);
//     const double center_det_jacobian = MathUtils<double>::Det(center_jacobian);

//     for (const auto& r_integration_point : rIntegrationPoints) {
//         CoordinatesArrayType local_coordinates = ZeroVector(3);
//         local_coordinates[0] = r_integration_point[0];
//         local_coordinates[1] = r_integration_point[1];
//         local_coordinates[2] = r_integration_point[2];

//         Matrix jacobian;
//         rGapVolume.Jacobian(jacobian, local_coordinates);
//         const double det_jacobian = MathUtils<double>::Det(jacobian);

//         min_det_jacobian = std::min(min_det_jacobian, det_jacobian);
//         max_det_jacobian = std::max(max_det_jacobian, det_jacobian);

//         if (det_jacobian > determinant_tolerance) {
//             ++number_of_positive_det;
//         } else if (det_jacobian < -determinant_tolerance) {
//             ++number_of_negative_det;
//         } else {
//             ++number_of_almost_zero_det;
//         }
//     }

//     const bool all_positive =
//         number_of_positive_det > 0 && number_of_negative_det == 0 && number_of_almost_zero_det == 0;
//     const bool all_negative =
//         number_of_negative_det > 0 && number_of_positive_det == 0 && number_of_almost_zero_det == 0;
//     const bool mixed_signs = number_of_positive_det > 0 && number_of_negative_det > 0;
//     const bool has_degenerate_points = number_of_almost_zero_det > 0;
//     const bool should_print = mEchoLevel > 2 || all_negative || mixed_signs || has_degenerate_points;

//     if (!should_print) {
//         return false;
//     }

//     std::string classification;
//     if (all_positive) {
//         classification = "orientation is valid";
//     } else if (all_negative) {
//         classification = "volume is globally inverted; integration weights are corrected by one consistent orientation sign";
//     } else if (mixed_signs) {
//         classification = "volume is locally folded; likely invalid Coons volume or inconsistent boundary curves";
//     } else if (has_degenerate_points) {
//         classification = "volume is degenerate or nearly degenerate";
//     } else {
//         classification = "no determinant samples were available";
//     }

//     KRATOS_WARNING_IF("CheckGapVolumeOrientation", all_negative || mixed_signs || has_degenerate_points)
//         << classification
//         << " center_detJ=" << center_det_jacobian
//         << " min_detJ=" << min_det_jacobian
//         << " max_detJ=" << max_det_jacobian
//         << " positive_detJ=" << number_of_positive_det
//         << " negative_detJ=" << number_of_negative_det
//         << " near_zero_detJ=" << number_of_almost_zero_det
//         << std::endl;
//     KRATOS_INFO_IF("CheckGapVolumeOrientation", mEchoLevel > 2 && all_positive)
//         << classification
//         << " center_detJ=" << center_det_jacobian
//         << " min_detJ=" << min_det_jacobian
//         << " max_detJ=" << max_det_jacobian
//         << std::endl;
//     return mixed_signs;
// }

// double SnakeGapSbmProcess::CalculateGapVolumeElementCharacteristicLength(
//     const NodeType& rNode00,
//     const NodeType& rNode01,
//     const NodeType& rNode10,
//     const NodeType& rNode11,
//     const array_1d<double,3>& rProjectedPoint00,
//     const array_1d<double,3>& rProjectedPoint01,
//     const array_1d<double,3>& rProjectedPoint10,
//     const array_1d<double,3>& rProjectedPoint11) const
// {
//     array_1d<double,3> center = ZeroVector(3);
//     center += rNode00.Coordinates();
//     center += rNode01.Coordinates();
//     center += rNode10.Coordinates();
//     center += rNode11.Coordinates();
//     center += rProjectedPoint00;
//     center += rProjectedPoint01;
//     center += rProjectedPoint10;
//     center += rProjectedPoint11;
//     center *= 0.125;

//     double radius_squared = 0.0;
//     const auto update_radius = [&center, &radius_squared](const array_1d<double,3>& rPoint) {
//         const array_1d<double,3> difference = center - rPoint;
//         const double distance_squared = inner_prod(difference, difference);
//         if (distance_squared > radius_squared) {
//             radius_squared = distance_squared;
//         }
//     };

//     update_radius(rNode00.Coordinates());
//     update_radius(rNode01.Coordinates());
//     update_radius(rNode10.Coordinates());
//     update_radius(rNode11.Coordinates());
//     update_radius(rProjectedPoint00);
//     update_radius(rProjectedPoint01);
//     update_radius(rProjectedPoint10);
//     update_radius(rProjectedPoint11);

//     radius_squared += 1.0e-15;
//     return std::sqrt(radius_squared);
// }

// void SnakeGapSbmProcess::CreateSbmExtendedGeometries3D()
// {
//     mEchoLevel = mThisParameters["echo_level"].GetInt();
//     KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 1)
//         << "Echo level: " << mEchoLevel
//         << " | Gap type: " << mGapSbmType
//         << " | Approx order: " << mGapApproximationOrder
//         << " | Internal divisions: " << mInternalDivisions
//         << " | Rel tol subdivisions: " << mGapRelativeToleranceForSubdivisions
//         << " | Interp levels: " << mNumberOfInterpolationLevels << std::endl;
//     if (mpSkinModelPartInnerInitial->NumberOfNodes()>0 || mpSkinModelPartInnerInitial->NumberOfGeometries()>0) 
//     {
//         const auto& r_surrogate_sub_model_part_inner = mpIgaModelPart->GetSubModelPart("surrogate_inner");
//         auto& r_skin_sub_model_part_inner = mpSkinModelPart->GetSubModelPart("inner");

//         KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
//             << "Creating the extended SBM geometries for the inner skin." << std::endl;
//         KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 2)
//             << "Inner skin nodes: " << r_skin_sub_model_part_inner.NumberOfNodes()
//             << " | conditions: " << r_skin_sub_model_part_inner.NumberOfConditions()
//             << " | surrogate conditions: " << r_surrogate_sub_model_part_inner.NumberOfConditions() << std::endl;
//         CreateSbmExtendedGeometries3D<true>(r_skin_sub_model_part_inner, r_surrogate_sub_model_part_inner);
//         KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
//             << "Finished creating the extended SBM geometries for the inner skin." << std::endl;
//     }
//     if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) 
//     {
//         const auto& r_surrogate_sub_model_part_outer = mpIgaModelPart->GetSubModelPart("surrogate_outer");
//         auto& r_skin_sub_model_part_outer = mpSkinModelPart->GetSubModelPart("outer");
//         KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
//             << "Creating the extended SBM geometries for the outer skin." << std::endl;
//         KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 2)
//             << "Outer skin nodes: " << r_skin_sub_model_part_outer.NumberOfNodes()
//             << " | conditions: " << r_skin_sub_model_part_outer.NumberOfConditions()
//             << " | surrogate conditions: " << r_surrogate_sub_model_part_outer.NumberOfConditions() << std::endl;
//         CreateSbmExtendedGeometries3D<false>(r_skin_sub_model_part_outer, r_surrogate_sub_model_part_outer);
//         KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
//             << "Finished creating the extended SBM geometries for the outer skin." << std::endl;
//     }
// }


// template <bool TIsInnerLoop>
// void SnakeGapSbmProcess::CreateSbmExtendedGeometries3D(
//     ModelPart& rSkinSubModelPart,
//     const ModelPart& rSurrogateSubModelPart)
// {
//     const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
//     KRATOS_ERROR_IF(knot_span_sizes.size() < 3)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: KNOT_SPAN_SIZES must contain at least three entries." << std::endl;

//     const double min_knot_span_size = std::min(knot_span_sizes[0], std::min(knot_span_sizes[1], knot_span_sizes[2]));
//     const double max_knot_span_size = std::max(knot_span_sizes[0], std::max(knot_span_sizes[1], knot_span_sizes[2]));
//     KRATOS_ERROR_IF(min_knot_span_size <= 0.0)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: knot span sizes must be positive." << std::endl;

//     const double area_tol = 1.0e-12 * min_knot_span_size * min_knot_span_size;
//     const double volume_tol = 1.0e-12 * min_knot_span_size * min_knot_span_size * min_knot_span_size;
//     const double distance_tol = 1.0e-12 * min_knot_span_size;
//     const double determinant_tol = 1.0e-14 * std::max(1.0, max_knot_span_size * max_knot_span_size * max_knot_span_size);

//     ModelPart& r_root_model_part = mpIgaModelPart->GetRootModelPart();
//     const bool store_gap_debug_geometries =
//         mThisParameters.Has("store_gap_debug_geometries")
//             ? mThisParameters["store_gap_debug_geometries"].GetBool()
//             : false;

//     ModelPart* p_gap_type1_debug = nullptr;
//     ModelPart* p_gap_type2_debug = nullptr;
//     ModelPart* p_gap_type3_debug = nullptr;
//     IndexType next_debug_geometry_id = 1;
//     std::unordered_map<const NodeType*, Node::Pointer> debug_node_map;
//     IndexType next_debug_node_id = 1;   

//     if (store_gap_debug_geometries) {
//         p_gap_type1_debug = r_root_model_part.HasSubModelPart("GapType1Debug")
//             ? &r_root_model_part.GetSubModelPart("GapType1Debug")
//             : &r_root_model_part.CreateSubModelPart("GapType1Debug");
//         p_gap_type2_debug = r_root_model_part.HasSubModelPart("GapType2Debug")
//             ? &r_root_model_part.GetSubModelPart("GapType2Debug")
//             : &r_root_model_part.CreateSubModelPart("GapType2Debug");
//         p_gap_type3_debug = r_root_model_part.HasSubModelPart("GapType3Debug")
//             ? &r_root_model_part.GetSubModelPart("GapType3Debug")
//             : &r_root_model_part.CreateSubModelPart("GapType3Debug");

//         for (const auto& r_geometry : r_root_model_part.Geometries()) {
//             next_debug_geometry_id = std::max(next_debug_geometry_id, r_geometry.Id() + 1);
//         }
//         for (const auto& r_node : r_root_model_part.Nodes()) {
//             next_debug_node_id = std::max(next_debug_node_id, r_node.Id() + 1);
//         }
//     }

//     auto get_next_debug_geometry_id = [&]() {
//         while (r_root_model_part.HasGeometry(next_debug_geometry_id)) {
//             ++next_debug_geometry_id;
//         }
//         return next_debug_geometry_id++;
//     };

//     auto get_next_debug_node_id = [&]() {
//         while (r_root_model_part.HasNode(next_debug_node_id)) {
//             ++next_debug_node_id;
//         }
//         return next_debug_node_id++;
//     };
    
//     auto get_or_create_debug_node = [&](ModelPart& rDebugPart, const Node::Pointer& pOriginalNode) -> Node::Pointer {
//         KRATOS_ERROR_IF_NOT(pOriginalNode)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null original node in debug geometry."
//             << std::endl;
    
//         const NodeType* p_key = pOriginalNode.get();
//         auto it = debug_node_map.find(p_key);
//         if (it != debug_node_map.end()) {
//             if (!rDebugPart.HasNode(it->second->Id())) {
//                 rDebugPart.AddNode(it->second);
//             }
//             return it->second;
//         }
    
//         if (r_root_model_part.HasNode(pOriginalNode->Id())) {
//             const auto& r_root_node_same_id = r_root_model_part.GetNode(pOriginalNode->Id());
//             const double coordinate_mismatch =
//                 norm_2(r_root_node_same_id.Coordinates() - pOriginalNode->Coordinates());
//         }
    
//         const IndexType new_debug_node_id = get_next_debug_node_id();
//         auto p_debug_node = r_root_model_part.CreateNewNode(
//             new_debug_node_id,
//             pOriginalNode->X(),
//             pOriginalNode->Y(),
//             pOriginalNode->Z());
    
//         rDebugPart.AddNode(p_debug_node);
//         debug_node_map.emplace(p_key, p_debug_node);
//         return p_debug_node;
//     };

//     auto p_volume = mpIgaModelPart->pGetGeometry(1);
//     KRATOS_ERROR_IF_NOT(p_volume)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: geometry with id 1 was not found." << std::endl;

//     auto p_nurbs_volume = std::dynamic_pointer_cast<NurbsVolumeType>(
//         p_volume->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
//     KRATOS_ERROR_IF_NOT(p_nurbs_volume)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: geometry with id 1 does not expose a NurbsVolumeType background geometry." << std::endl;

//     const auto bounds = ReadParameterSpaceBounds3D(
//         rSurrogateSubModelPart.GetParentModelPart(),
//         "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D");

//     const int number_of_spans_x = static_cast<int>(ComputeSpanCount3D(
//         bounds.MaxU - bounds.MinU,
//         knot_span_sizes[0],
//         "u",
//         "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D"));
    
//     const int number_of_spans_y = static_cast<int>(ComputeSpanCount3D(
//         bounds.MaxV - bounds.MinV,
//         knot_span_sizes[1],
//         "v",
//         "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D"));
    
//     const int number_of_spans_z = static_cast<int>(ComputeSpanCount3D(
//         bounds.MaxW - bounds.MinW,
//         knot_span_sizes[2],
//         "w",
//         "SnakeGapSbmProcess::CreateSbmExtendedGeometries3D"));

//     const std::size_t brep_degree = p_nurbs_volume->PolynomialDegree(0);
//     const std::size_t number_of_shape_functions_derivatives = 3 * brep_degree + 1;
//     if (mGapApproximationOrder == 0) {
//         mGapApproximationOrder = brep_degree;
//     }

//     struct SpanKey3D
//     {
//         int I = 0;
//         int J = 0;
//         int K = 0;

//         bool operator<(const SpanKey3D& rOther) const
//         {
//             return std::tie(I, J, K) < std::tie(rOther.I, rOther.J, rOther.K);
//         }

//         bool operator==(const SpanKey3D& rOther) const
//         {
//             return I == rOther.I && J == rOther.J && K == rOther.K;
//         }
//     };

//     struct GridPointKey3D
//     {
//         int I = 0;
//         int J = 0;
//         int K = 0;
//     };

//     struct VertexKey
//     {
//         int Kind = 0; // 0: surrogate, 1: skin
//         IndexType Id = 0;

//         bool operator<(const VertexKey& rOther) const
//         {
//             return std::tie(Kind, Id) < std::tie(rOther.Kind, rOther.Id);
//         }

//         bool operator==(const VertexKey& rOther) const
//         {
//             return Kind == rOther.Kind && Id == rOther.Id;
//         }
//     };

//     struct LateralFaceKey
//     {
//         std::vector<VertexKey> Vertices;

//         bool operator<(const LateralFaceKey& rOther) const
//         {
//             return Vertices < rOther.Vertices;
//         }
//     };

//     struct LateralFaceData
//     {
//         LateralFaceKey Key;
//         NurbsSurfaceType::Pointer pGeometry;
//         std::vector<Geometry<Node>::Pointer> NeighbourGeometries;
//         std::vector<IndexType> AdjacentGapVolumeIds;
//         bool IsTriangle = false;
//         bool IsDegenerate = false;
//     };

//     struct ProjectionFaceData
//     {
//         Node::Pointer pS0;
//         Node::Pointer pS1;
//         Node::Pointer pK0;
//         Node::Pointer pK1;
//         SpanKey3D ExternalSpan;
//         SpanKey3D ActiveSpan;
//         std::vector<Geometry<Node>::Pointer> CandidateNeighbourGeometries;
//         int NormalAxis = -1;
//         int NormalSign = 0;
//         bool IsType1 = false;
//     };

//     struct SurrogateEdgeData
//     {
//         Node::Pointer pS0;
//         Node::Pointer pS1;
//         std::vector<ProjectionFaceData> ProjectionFaces;
//         std::set<SpanKey3D> ExternalSpanIds;
//     };

//     struct LocalProjectionData
//     {
//         Node::Pointer pSurrogateNode;
//         Node::Pointer pSkinNode;
//         SpanKey3D ExternalSpan;
//         std::array<int, 3> SignTriple{{0, 0, 0}};
//         array_1d<double, 3> Vector = ZeroVector(3);
//     };

//     struct SurrogateNodeData
//     {
//         Node::Pointer pNode;
//         GridPointKey3D Grid;
//         std::vector<LocalProjectionData> ProjectionSegments;
//         std::vector<Geometry<Node>::Pointer> IncidentFaceNeighbourGeometries;
//     };

//     struct SurrogateFaceData
//     {
//         const Condition* pCondition = nullptr;
//         std::array<Node::Pointer, 4> Nodes; // S00, S10, S11, S01
//         std::array<GridPointKey3D, 4> GridNodes;
//         array_1d<double, 3> OutwardNormal = ZeroVector(3);
//         SpanKey3D ExternalSpan;
//         SpanKey3D ActiveSpan;
//         Geometry<Node>::Pointer pNeighbourGeometry;
//         int NormalAxis = -1;
//         int NormalSign = 0;
//     };

//     struct GapVolumeData
//     {
//         int Type = 0;
//         NurbsVolumeType::Pointer pVolumeGeometry;
//         NurbsSurfaceType::Pointer pTopSurfaceGeometry;
//         Geometry<Node>::Pointer pChosenNeighbourGeometry;
//         std::vector<LateralFaceKey> LateralFaceKeys;
//     };

//     struct VolumeDiagnostics
//     {
//         double MinDeterminant = std::numeric_limits<double>::max();
//         double MaxDeterminant = -std::numeric_limits<double>::max();
//         double WeightSum = 0.0;
//         std::size_t PositiveDeterminants = 0;
//         std::size_t NegativeDeterminants = 0;
//         std::size_t NearZeroDeterminants = 0;

//         bool HasMixedSigns() const
//         {
//             return PositiveDeterminants > 0 && NegativeDeterminants > 0;
//         }

//         bool HasConsistentNegativeSign() const
//         {
//             return NegativeDeterminants > 0 && PositiveDeterminants == 0;
//         }
//     };

//     struct GapTypeDiagnostics
//     {
//         std::size_t Volumes = 0;
//         std::size_t TopConditions = 0;
//         std::size_t CollapsedTops = 0;
//         double MinWeight = std::numeric_limits<double>::max();
//         double MaxWeight = 0.0;
//         double MinDeterminant = std::numeric_limits<double>::max();

//         void UpdateVolume(const double WeightSum, const double MinDet)
//         {
//             ++Volumes;
//             MinWeight = std::min(MinWeight, WeightSum);
//             MaxWeight = std::max(MaxWeight, WeightSum);
//             MinDeterminant = std::min(MinDeterminant, MinDet);
//         }

//         double PrintableMinWeight() const
//         {
//             return MinWeight == std::numeric_limits<double>::max() ? 0.0 : MinWeight;
//         }

//         double PrintableMinDeterminant() const
//         {
//             return MinDeterminant == std::numeric_limits<double>::max() ? 0.0 : MinDeterminant;
//         }
//     };

//     auto span_to_string = [](const SpanKey3D& rSpan) {
//         std::ostringstream buffer;
//         buffer << "(" << rSpan.I << "," << rSpan.J << "," << rSpan.K << ")";
//         return buffer.str();
//     };

//     auto get_span_component = [](const SpanKey3D& rSpan, const int Axis) {
//         return Axis == 0 ? rSpan.I : (Axis == 1 ? rSpan.J : rSpan.K);
//     };

//     auto add_to_span = [](SpanKey3D& rSpan, const int Axis, const int Value) {
//         if (Axis == 0) {
//             rSpan.I += Value;
//         } else if (Axis == 1) {
//             rSpan.J += Value;
//         } else {
//             rSpan.K += Value;
//         }
//     };

//     auto is_span_inside_domain = [&](const SpanKey3D& rSpan) {
//         return rSpan.I >= 0 && rSpan.J >= 0 && rSpan.K >= 0 &&
//                rSpan.I < number_of_spans_x &&
//                rSpan.J < number_of_spans_y &&
//                rSpan.K < number_of_spans_z;
//     };

//     auto span_to_box = [&](const SpanKey3D& rSpan) {
//         std::array<array_1d<double, 3>, 2> box;
    
//         box[0] = ZeroVector(3);
//         box[1] = ZeroVector(3);
    
//         box[0][0] = bounds.MinU + rSpan.I * knot_span_sizes[0];
//         box[0][1] = bounds.MinV + rSpan.J * knot_span_sizes[1];
//         box[0][2] = bounds.MinW + rSpan.K * knot_span_sizes[2];
    
//         box[1][0] = box[0][0] + knot_span_sizes[0];
//         box[1][1] = box[0][1] + knot_span_sizes[1];
//         box[1][2] = box[0][2] + knot_span_sizes[2];
    
//         return box;
//     };
    
//     const double span_tolerance =
//         1.0e-10 * std::max({knot_span_sizes[0], knot_span_sizes[1], knot_span_sizes[2]});
    
//     auto is_point_inside_box = [&](const array_1d<double, 3>& rPoint,
//                                    const array_1d<double, 3>& rBoxMin,
//                                    const array_1d<double, 3>& rBoxMax) {
//         return rPoint[0] >= rBoxMin[0] - span_tolerance && rPoint[0] <= rBoxMax[0] + span_tolerance &&
//                rPoint[1] >= rBoxMin[1] - span_tolerance && rPoint[1] <= rBoxMax[1] + span_tolerance &&
//                rPoint[2] >= rBoxMin[2] - span_tolerance && rPoint[2] <= rBoxMax[2] + span_tolerance;
//     };
    
//     auto point_span_index = [&](const double Coordinate,
//                                 const double MinValue,
//                                 const double SpanSize,
//                                 const int NumberOfSpans) {
//         int index = static_cast<int>(std::floor((Coordinate - MinValue) / SpanSize + 1.0e-12));
//         index = std::max(0, std::min(index, NumberOfSpans - 1));
//         return index;
//     };
    
//     auto node_is_in_span = [&](const Node& rNode, const SpanKey3D& rSpan) {
//         const int i = point_span_index(rNode.X(), bounds.MinU, knot_span_sizes[0], number_of_spans_x);
//         const int j = point_span_index(rNode.Y(), bounds.MinV, knot_span_sizes[1], number_of_spans_y);
//         const int k = point_span_index(rNode.Z(), bounds.MinW, knot_span_sizes[2], number_of_spans_z);
    
//         return i == rSpan.I && j == rSpan.J && k == rSpan.K;
//     };
    
//     auto span_has_skin_node = [&](const SpanKey3D& rSpan) {
//         for (const auto& r_node : rSkinSubModelPart.Nodes()) {
//             if (node_is_in_span(r_node, rSpan)) {
//                 return true;
//             }
//         }
//         return false;
//     };
    
//     auto clip_segment_with_box = [&](const array_1d<double, 3>& rA,
//                                      const array_1d<double, 3>& rB,
//                                      const array_1d<double, 3>& rBoxMin,
//                                      const array_1d<double, 3>& rBoxMax,
//                                      array_1d<double, 3>& rPointInsideBox) {
//         double t_min = 0.0;
//         double t_max = 1.0;
    
//         const array_1d<double, 3> direction = rB - rA;
    
//         for (IndexType axis = 0; axis < 3; ++axis) {
//             if (std::abs(direction[axis]) <= span_tolerance) {
//                 if (rA[axis] < rBoxMin[axis] - span_tolerance ||
//                     rA[axis] > rBoxMax[axis] + span_tolerance) {
//                     return false;
//                 }
//             } else {
//                 double t1 = (rBoxMin[axis] - rA[axis]) / direction[axis];
//                 double t2 = (rBoxMax[axis] - rA[axis]) / direction[axis];
    
//                 if (t1 > t2) {
//                     std::swap(t1, t2);
//                 }
    
//                 t_min = std::max(t_min, t1);
//                 t_max = std::min(t_max, t2);
    
//                 if (t_min > t_max + span_tolerance) {
//                     return false;
//                 }
//             }
//         }
    
//         const double t = 0.5 * (t_min + t_max);
//         rPointInsideBox = rA + t * direction;
    
//         return is_point_inside_box(rPointInsideBox, rBoxMin, rBoxMax);
//     };
    
//     auto is_point_inside_triangle = [&](const array_1d<double, 3>& rPoint,
//                                         const array_1d<double, 3>& rA,
//                                         const array_1d<double, 3>& rB,
//                                         const array_1d<double, 3>& rC) {
//         const array_1d<double, 3> v0 = rC - rA;
//         const array_1d<double, 3> v1 = rB - rA;
//         const array_1d<double, 3> v2 = rPoint - rA;
    
//         const double dot00 = inner_prod(v0, v0);
//         const double dot01 = inner_prod(v0, v1);
//         const double dot02 = inner_prod(v0, v2);
//         const double dot11 = inner_prod(v1, v1);
//         const double dot12 = inner_prod(v1, v2);
    
//         const double denominator = dot00 * dot11 - dot01 * dot01;
//         if (std::abs(denominator) <= std::numeric_limits<double>::epsilon()) {
//             return false;
//         }
    
//         const double inv_denominator = 1.0 / denominator;
//         const double u = (dot11 * dot02 - dot01 * dot12) * inv_denominator;
//         const double v = (dot00 * dot12 - dot01 * dot02) * inv_denominator;
    
//         return u >= -1.0e-10 && v >= -1.0e-10 && u + v <= 1.0 + 1.0e-10;
//     };
    
//     auto find_triangle_box_intersection_point = [&](const Geometry<Node>& rTriangle,
//                                                     const SpanKey3D& rSpan,
//                                                     array_1d<double, 3>& rPointInsideSpan) {
//         if (rTriangle.PointsNumber() != 3) {
//             return false;
//         }
    
//         const auto box = span_to_box(rSpan);
//         const array_1d<double, 3>& r_box_min = box[0];
//         const array_1d<double, 3>& r_box_max = box[1];
    
//         const array_1d<double, 3> a = rTriangle[0].Coordinates();
//         const array_1d<double, 3> b = rTriangle[1].Coordinates();
//         const array_1d<double, 3> c = rTriangle[2].Coordinates();
    
//         if (is_point_inside_box(a, r_box_min, r_box_max)) {
//             rPointInsideSpan = a;
//             return true;
//         }
    
//         if (is_point_inside_box(b, r_box_min, r_box_max)) {
//             rPointInsideSpan = b;
//             return true;
//         }
    
//         if (is_point_inside_box(c, r_box_min, r_box_max)) {
//             rPointInsideSpan = c;
//             return true;
//         }
    
//         if (clip_segment_with_box(a, b, r_box_min, r_box_max, rPointInsideSpan)) {
//             return true;
//         }
    
//         if (clip_segment_with_box(b, c, r_box_min, r_box_max, rPointInsideSpan)) {
//             return true;
//         }
    
//         if (clip_segment_with_box(c, a, r_box_min, r_box_max, rPointInsideSpan)) {
//             return true;
//         }
    
//         const array_1d<double, 3> ab = b - a;
//         const array_1d<double, 3> ac = c - a;
    
//         array_1d<double, 3> normal;
//         MathUtils<double>::CrossProduct(normal, ab, ac);
    
//         const double normal_norm = norm_2(normal);
//         if (normal_norm <= std::numeric_limits<double>::epsilon()) {
//             return false;
//         }
    
//         normal /= normal_norm;
    
//         const std::array<array_1d<double, 3>, 8> box_corners = {{
//             array_1d<double, 3>{r_box_min[0], r_box_min[1], r_box_min[2]},
//             array_1d<double, 3>{r_box_max[0], r_box_min[1], r_box_min[2]},
//             array_1d<double, 3>{r_box_min[0], r_box_max[1], r_box_min[2]},
//             array_1d<double, 3>{r_box_max[0], r_box_max[1], r_box_min[2]},
//             array_1d<double, 3>{r_box_min[0], r_box_min[1], r_box_max[2]},
//             array_1d<double, 3>{r_box_max[0], r_box_min[1], r_box_max[2]},
//             array_1d<double, 3>{r_box_min[0], r_box_max[1], r_box_max[2]},
//             array_1d<double, 3>{r_box_max[0], r_box_max[1], r_box_max[2]}
//         }};
    
//         for (const auto& r_corner : box_corners) {
//             const double signed_distance = inner_prod(normal, r_corner - a);
//             const array_1d<double, 3> projected_point = r_corner - signed_distance * normal;
    
//             if (is_point_inside_box(projected_point, r_box_min, r_box_max) &&
//                 is_point_inside_triangle(projected_point, a, b, c)) {
//                 rPointInsideSpan = projected_point;
//                 return true;
//             }
//         }
    
//         return false;
//     };
    
//     auto get_next_auxiliary_skin_node_id = [&]() {
//         auto& r_root_model_part = rSkinSubModelPart.GetRootModelPart();
    
//         IndexType next_id = 1;
//         if (r_root_model_part.NumberOfNodes() > 0) {
//             next_id = (r_root_model_part.NodesEnd() - 1)->Id() + 1;
//         }
    
//         while (r_root_model_part.HasNode(next_id)) {
//             ++next_id;
//         }
    
//         return next_id;
//     };
    
//     auto clamp_point_inside_span = [&](array_1d<double, 3>& rPoint,
//                                        const SpanKey3D& rSpan) {
//         const auto box = span_to_box(rSpan);
    
//         rPoint[0] = std::max(box[0][0] + span_tolerance, std::min(box[1][0] - span_tolerance, rPoint[0]));
//         rPoint[1] = std::max(box[0][1] + span_tolerance, std::min(box[1][1] - span_tolerance, rPoint[1]));
//         rPoint[2] = std::max(box[0][2] + span_tolerance, std::min(box[1][2] - span_tolerance, rPoint[2]));
//     };
    
//     auto ensure_skin_nodes_in_external_spans = [&](const std::set<SpanKey3D>& rExternalSpans) {
//         for (const auto& r_external_span : rExternalSpans) {
//             if (!is_span_inside_domain(r_external_span)) {
//                 continue;
//             }
    
//             if (span_has_skin_node(r_external_span)) {
//                 continue;
//             }
    
//             bool found_intersecting_triangle = false;
//             array_1d<double, 3> auxiliary_point = ZeroVector(3);
    
//             for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
//                 const auto& r_geometry = r_condition.GetGeometry();
    
//                 if (find_triangle_box_intersection_point(
//                         r_geometry,
//                         r_external_span,
//                         auxiliary_point)) {
//                     found_intersecting_triangle = true;
//                     break;
//                 }
//             }
    
//             KRATOS_ERROR_IF_NOT(found_intersecting_triangle)
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: external span "
//                 << span_to_string(r_external_span)
//                 << " has no skin node and no intersecting skin triangle was found." << std::endl;
    
//             clamp_point_inside_span(auxiliary_point, r_external_span);
    
//             const IndexType new_node_id = get_next_auxiliary_skin_node_id();
    
//             rSkinSubModelPart.CreateNewNode(
//                 new_node_id,
//                 auxiliary_point[0],
//                 auxiliary_point[1],
//                 auxiliary_point[2]);
    
//             KRATOS_INFO_IF("CreateSbmExtendedGeometries3D", mEchoLevel > 1)
//                 << "Added auxiliary skin node " << new_node_id
//                 << " in external span " << span_to_string(r_external_span)
//                 << " at " << auxiliary_point << std::endl;
//         }
//     };

//     auto compute_grid_index = [&](const array_1d<double, 3>& rPoint) {
//         GridPointKey3D key;
//         key.I = static_cast<int>(std::llround((rPoint[0] - bounds.MinU) / knot_span_sizes[0]));
//         key.J = static_cast<int>(std::llround((rPoint[1] - bounds.MinV) / knot_span_sizes[1]));
//         key.K = static_cast<int>(std::llround((rPoint[2] - bounds.MinW) / knot_span_sizes[2]));
//         return key;
//     };

//     auto span_sign_around_node = [](const SpanKey3D& rSpan, const GridPointKey3D& rGrid) {
//         std::array<int, 3> signs{{-1, -1, -1}};
//         signs[0] = rSpan.I >= rGrid.I ? 1 : -1;
//         signs[1] = rSpan.J >= rGrid.J ? 1 : -1;
//         signs[2] = rSpan.K >= rGrid.K ? 1 : -1;
//         return signs;
//     };

//     auto make_vertex_key = [](const Node::Pointer& pNode, const int Kind) {
//         KRATOS_ERROR_IF_NOT(pNode)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null node in lateral face key." << std::endl;
//         return VertexKey{Kind, pNode->Id()};
//     };

//     auto make_face_key = [&](const std::vector<std::pair<Node::Pointer, int>>& rNodes) {
//         LateralFaceKey key;
//         key.Vertices.reserve(rNodes.size());
//         for (const auto& r_node_data : rNodes) {
//             key.Vertices.push_back(make_vertex_key(r_node_data.first, r_node_data.second));
//         }
//         std::sort(key.Vertices.begin(), key.Vertices.end());
//         key.Vertices.erase(std::unique(key.Vertices.begin(), key.Vertices.end()), key.Vertices.end());
//         return key;
//     };

//     auto create_unit_knot_vector = []() {
//         Vector knot_vector = ZeroVector(4);
//         knot_vector[2] = 1.0;
//         knot_vector[3] = 1.0;
//         return knot_vector;
//     };

//     auto create_linear_surface = [&](const std::vector<Node::Pointer>& rOrderedNodes) {
//         KRATOS_ERROR_IF(rOrderedNodes.size() != 3 && rOrderedNodes.size() != 4)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: linear surface requires three or four nodes." << std::endl;

//         PointerVector<NodeType> control_points;
//         if (rOrderedNodes.size() == 3) {
//             control_points.push_back(rOrderedNodes[0]);
//             control_points.push_back(rOrderedNodes[1]);
//             control_points.push_back(rOrderedNodes[2]);
//             control_points.push_back(rOrderedNodes[2]);
//         } else {
//             control_points.push_back(rOrderedNodes[0]); // P00
//             control_points.push_back(rOrderedNodes[1]); // P10
//             control_points.push_back(rOrderedNodes[3]); // P01
//             control_points.push_back(rOrderedNodes[2]); // P11
//         }

//         return Kratos::make_shared<NurbsSurfaceType>(
//             control_points,
//             std::size_t(1),
//             std::size_t(1),
//             create_unit_knot_vector(),
//             create_unit_knot_vector());
//     };

//     auto create_linear_volume = [&](const std::array<Node::Pointer, 8>& rControlNodes) {
//         PointerVector<NodeType> control_points;
//         for (const auto& p_node : rControlNodes) {
//             KRATOS_ERROR_IF_NOT(p_node)
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null control point in gap volume." << std::endl;
//             control_points.push_back(p_node);
//         }

//         return Kratos::make_shared<NurbsVolumeType>(
//             control_points,
//             std::size_t(1),
//             std::size_t(1),
//             std::size_t(1),
//             create_unit_knot_vector(),
//             create_unit_knot_vector(),
//             create_unit_knot_vector());
//     };

//     auto compute_face_area = [&](const std::vector<Node::Pointer>& rOrderedNodes) {
//         if (rOrderedNodes.size() < 3) {
//             return 0.0;
//         }

//         const array_1d<double, 3> p0 = rOrderedNodes[0]->Coordinates();
//         double area = 0.0;
//         for (std::size_t i = 1; i + 1 < rOrderedNodes.size(); ++i) {
//             const array_1d<double, 3> a = rOrderedNodes[i]->Coordinates() - p0;
//             const array_1d<double, 3> b = rOrderedNodes[i + 1]->Coordinates() - p0;
//             area += 0.5 * norm_2(MathUtils<double>::CrossProduct(a, b));
//         }
//         return area;
//     };

//     auto compute_polygon_normal = [](const std::vector<Node::Pointer>& rOrderedNodes) {
//         array_1d<double, 3> normal = ZeroVector(3);
//         if (rOrderedNodes.size() < 3) {
//             return normal;
//         }

//         for (std::size_t i = 0; i < rOrderedNodes.size(); ++i) {
//             const auto& r_current_point = rOrderedNodes[i]->Coordinates();
//             const auto& r_next_point = rOrderedNodes[(i + 1) % rOrderedNodes.size()]->Coordinates();
//             normal += MathUtils<double>::CrossProduct(r_current_point, r_next_point);
//         }
//         return normal;
//     };

//     auto compress_ordered_face_nodes = [](const std::vector<Node::Pointer>& rOrderedNodes) {
//         std::vector<Node::Pointer> compressed_nodes;
//         compressed_nodes.reserve(rOrderedNodes.size());
//         for (const auto& p_node : rOrderedNodes) {
//             if (!p_node) {
//                 continue;
//             }
//             if (compressed_nodes.empty() || compressed_nodes.back()->Id() != p_node->Id()) {
//                 compressed_nodes.push_back(p_node);
//             }
//         }
//         if (compressed_nodes.size() > 1 &&
//             compressed_nodes.front()->Id() == compressed_nodes.back()->Id()) {
//             compressed_nodes.pop_back();
//         }
//         return compressed_nodes;
//     };

//     auto is_ordered_face_non_degenerate = [&](const std::vector<Node::Pointer>& rOrderedNodes,
//                                               const array_1d<double, 3>* pExpectedNormal,
//                                               std::string& rFailureReason) {
//         const std::vector<Node::Pointer> face_nodes = compress_ordered_face_nodes(rOrderedNodes);
//         if (face_nodes.size() < 3) {
//             rFailureReason = "less than three distinct ordered top nodes";
//             return false;
//         }

//         for (std::size_t i = 0; i < face_nodes.size(); ++i) {
//             const array_1d<double, 3> edge =
//                 face_nodes[(i + 1) % face_nodes.size()]->Coordinates() - face_nodes[i]->Coordinates();
//             if (norm_2(edge) <= distance_tol) {
//                 rFailureReason = "zero-length top edge";
//                 return false;
//             }
//         }

//         const double area = compute_face_area(face_nodes);
//         if (area <= area_tol) {
//             rFailureReason = "top area below tolerance";
//             return false;
//         }

//         array_1d<double, 3> reference_normal = compute_polygon_normal(face_nodes);
//         const double reference_normal_norm = norm_2(reference_normal);
//         if (reference_normal_norm <= area_tol) {
//             rFailureReason = "top polygon normal below tolerance";
//             return false;
//         }

//         if (pExpectedNormal) {
//             const double signed_orientation = inner_prod(reference_normal, *pExpectedNormal);
//             if (signed_orientation <= area_tol) {
//                 rFailureReason = "top orientation is inconsistent with the expected normal";
//                 return false;
//             }
//             reference_normal = *pExpectedNormal;
//         } else {
//             reference_normal /= reference_normal_norm;
//         }

//         for (std::size_t i = 0; i < face_nodes.size(); ++i) {
//             const array_1d<double, 3> edge_to_next =
//                 face_nodes[(i + 1) % face_nodes.size()]->Coordinates() - face_nodes[i]->Coordinates();
//             const array_1d<double, 3> edge_to_previous =
//                 face_nodes[(i + face_nodes.size() - 1) % face_nodes.size()]->Coordinates() - face_nodes[i]->Coordinates();
//             const double corner_measure =
//                 inner_prod(MathUtils<double>::CrossProduct(edge_to_next, edge_to_previous), reference_normal);
//             if (corner_measure <= area_tol) {
//                 rFailureReason = "top face is concave, inverted, or has a nearly collinear corner";
//                 return false;
//             }
//         }

//         return true;
//     };

//     auto is_top_quad_sane = [&](const std::array<Node::Pointer, 4>& rBottomNodes,
//                                 const std::array<Node::Pointer, 4>& rTopNodes,
//                                 const array_1d<double, 3>& rOutwardNormal,
//                                 std::string& rFailureReason) {
//         for (std::size_t i = 0; i < rTopNodes.size(); ++i) {
//             KRATOS_ERROR_IF_NOT(rBottomNodes[i])
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null bottom node in type-1 top validation." << std::endl;
//             KRATOS_ERROR_IF_NOT(rTopNodes[i])
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null top node in type-1 top validation." << std::endl;

//             const array_1d<double, 3> projection_vector =
//                 rTopNodes[i]->Coordinates() - rBottomNodes[i]->Coordinates();
//             if (inner_prod(projection_vector, rOutwardNormal) <= distance_tol) {
//                 rFailureReason = "projection vector is not outward";
//                 return false;
//             }
//         }

//         std::set<IndexType> unique_top_node_ids;
//         for (const auto& p_node : rTopNodes) {
//             unique_top_node_ids.insert(p_node->Id());
//         }
//         if (unique_top_node_ids.size() != rTopNodes.size()) {
//             rFailureReason = "top quad has repeated nodes";
//             return false;
//         }

//         return is_ordered_face_non_degenerate(
//             {rTopNodes[0], rTopNodes[1], rTopNodes[2], rTopNodes[3]},
//             &rOutwardNormal,
//             rFailureReason);
//     };

//     auto find_problematic_top_quad_corners = [&](const std::array<Node::Pointer, 4>& rBottomNodes,
//                                                  const std::array<Node::Pointer, 4>& rTopNodes,
//                                                  const array_1d<double, 3>& rOutwardNormal) {
//         std::vector<std::size_t> problematic_corners;
//         for (std::size_t i = 0; i < rTopNodes.size(); ++i) {
//             const array_1d<double, 3> projection_vector =
//                 rTopNodes[i]->Coordinates() - rBottomNodes[i]->Coordinates();
//             const array_1d<double, 3> edge_to_next =
//                 rTopNodes[(i + 1) % rTopNodes.size()]->Coordinates() - rTopNodes[i]->Coordinates();
//             const array_1d<double, 3> edge_to_previous =
//                 rTopNodes[(i + rTopNodes.size() - 1) % rTopNodes.size()]->Coordinates() - rTopNodes[i]->Coordinates();
//             const double corner_measure =
//                 inner_prod(MathUtils<double>::CrossProduct(edge_to_next, edge_to_previous), rOutwardNormal);
//             if (inner_prod(projection_vector, rOutwardNormal) <= distance_tol ||
//                 norm_2(edge_to_next) <= distance_tol ||
//                 norm_2(edge_to_previous) <= distance_tol ||
//                 corner_measure <= area_tol) {
//                 problematic_corners.push_back(i);
//             }
//         }
//         if (problematic_corners.empty()) {
//             problematic_corners = {0, 1, 2, 3};
//         }
//         return problematic_corners;
//     };

//     auto collect_unique_face_nodes = [](const std::array<Node::Pointer, 4>& rFaceNodes) {
//         std::vector<Node::Pointer> unique_nodes;
//         unique_nodes.reserve(4);
//         std::unordered_set<IndexType> used_node_ids;
//         for (const auto& p_node : rFaceNodes) {
//             if (p_node && used_node_ids.insert(p_node->Id()).second) {
//                 unique_nodes.push_back(p_node);
//             }
//         }
//         return unique_nodes;
//     };

//     auto add_debug_face_geometry = [&](const int GapType, const std::array<Node::Pointer, 4>& rFaceNodes) {
//         if (!store_gap_debug_geometries) {
//             return;
//         }

//         ModelPart* p_debug_part = nullptr;
//         if (GapType == 1) {
//             p_debug_part = p_gap_type1_debug;
//         } else if (GapType == 2) {
//             p_debug_part = p_gap_type2_debug;
//         } else if (GapType == 3) {
//             p_debug_part = p_gap_type3_debug;
//         }
//         if (!p_debug_part) {
//             return;
//         }

//         std::vector<Node::Pointer> face_nodes = collect_unique_face_nodes(rFaceNodes);
//         if (face_nodes.size() < 3 || face_nodes.size() > 4) {
//             return;
//         }
//         if (compute_face_area(face_nodes) <= area_tol) {
//             return;
//         }

//         std::vector<Node::Pointer> debug_face_nodes;
//         debug_face_nodes.reserve(face_nodes.size());
//         for (const auto& p_node : face_nodes) {
//             debug_face_nodes.push_back(get_or_create_debug_node(*p_debug_part, p_node));
//         }

//         std::vector<IndexType> debug_node_ids;
//         debug_node_ids.reserve(debug_face_nodes.size());
//         for (const auto& p_debug_node : debug_face_nodes) {
//             debug_node_ids.push_back(p_debug_node->Id());
//         }

//         const IndexType geometry_id = get_next_debug_geometry_id();
//         if (debug_node_ids.size() == 3) {
//             p_debug_part->CreateNewGeometry("Triangle3D3", geometry_id, debug_node_ids);
//         } else {
//             p_debug_part->CreateNewGeometry("Quadrilateral3D4", geometry_id, debug_node_ids);
//         }
//     };

//     auto store_gap_debug_faces = [&](const int GapType, const std::array<Node::Pointer, 8>& rVolumeNodes) {
//         if (!store_gap_debug_geometries) {
//             return;
//         }

//         const std::array<std::array<std::size_t, 4>, 6> face_indices = {{
//             {{0, 1, 3, 2}},
//             {{4, 5, 7, 6}},
//             {{0, 1, 5, 4}},
//             {{1, 3, 7, 5}},
//             {{2, 3, 7, 6}},
//             {{0, 2, 6, 4}}
//         }};

//         for (const auto& r_face : face_indices) {
//             add_debug_face_geometry(
//                 GapType,
//                 {{rVolumeNodes[r_face[0]], rVolumeNodes[r_face[1]], rVolumeNodes[r_face[2]], rVolumeNodes[r_face[3]]}});
//         }
//     };

//     KnotSpanSkinBinsCSR skin_bins_per_knot_span;

//     auto characteristic_length = [&](const std::vector<Node::Pointer>& rNodes) {
//         KRATOS_ERROR_IF(rNodes.empty())
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: cannot compute characteristic length from an empty node list." << std::endl;

//         array_1d<double, 3> center = ZeroVector(3);
//         for (const auto& p_node : rNodes) {
//             center += p_node->Coordinates();
//         }
//         center /= static_cast<double>(rNodes.size());

//         double radius_squared = 0.0;
//         for (const auto& p_node : rNodes) {
//             const array_1d<double, 3> diff = p_node->Coordinates() - center;
//             radius_squared = std::max(radius_squared, inner_prod(diff, diff));
//         }

//         const double length = std::sqrt(radius_squared);
//         KRATOS_ERROR_IF(length <= distance_tol)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: characteristic length is not positive." << std::endl;
//         return length;
//     };
//     //------------------------------------------------------------------------------------
//     auto find_closest_skin_node_in_external_span = [&](const array_1d<double, 3>& rPoint,
//         const SpanKey3D& rExternalSpan) -> Node::Pointer {
//         KRATOS_ERROR_IF_NOT(is_span_inside_domain(rExternalSpan))
//         << "::[SnakeGapSbmProcess]::FindClosestSkinNodeInExternalSpan: external span "
//         << span_to_string(rExternalSpan) << " is outside the knot-span grid." << std::endl;

//         const std::size_t nnz_index = FindSpanNnzIndex(
//         skin_bins_per_knot_span,
//         static_cast<std::size_t>(rExternalSpan.I),
//         static_cast<std::size_t>(rExternalSpan.J),
//         static_cast<std::size_t>(rExternalSpan.K));
//         KRATOS_ERROR_IF(nnz_index == static_cast<std::size_t>(-1) ||
//         nnz_index >= skin_bins_per_knot_span.CellBinsByNnz.size() ||
//         skin_bins_per_knot_span.CellBinsByNnz[nnz_index].Nodes.empty())
//         << "::[SnakeGapSbmProcess]::FindClosestSkinNodeInExternalSpan: no skin node exists in external span "
//         << span_to_string(rExternalSpan) << "." << std::endl;

//         Node::Pointer p_best_node;
//         double best_distance_squared = std::numeric_limits<double>::max();
//         for (const auto& p_skin_node : skin_bins_per_knot_span.CellBinsByNnz[nnz_index].Nodes) {
//         const array_1d<double, 3> diff = p_skin_node->Coordinates() - rPoint;
//         const double distance_squared = inner_prod(diff, diff);
//         if (distance_squared < best_distance_squared) {
//         best_distance_squared = distance_squared;
//         p_best_node = p_skin_node;
//         }
//         }

//         KRATOS_ERROR_IF_NOT(p_best_node)
//         << "::[SnakeGapSbmProcess]::FindClosestSkinNodeInExternalSpan: failed to select a skin node in span "
//         << span_to_string(rExternalSpan) << "." << std::endl;
//         return p_best_node;
//         };

//         auto find_closest_skin_node_to_entity_in_external_span = [&](const std::vector<Node::Pointer>& rSurrogateEntityNodes,
//                         const SpanKey3D& rExternalSpan) {
//         KRATOS_ERROR_IF(rSurrogateEntityNodes.empty())
//         << "::[SnakeGapSbmProcess]::FindClosestSkinNodeToEntityInExternalSpan: empty surrogate entity." << std::endl;
//         array_1d<double, 3> centroid = ZeroVector(3);
//         for (const auto& p_node : rSurrogateEntityNodes) {
//         centroid += p_node->Coordinates();
//         }
//         centroid /= static_cast<double>(rSurrogateEntityNodes.size());
//         return find_closest_skin_node_in_external_span(centroid, rExternalSpan);
//         };

//         auto get_or_create_closest_skin_node_in_external_span =
//         [&](const array_1d<double, 3>& rPoint,
//             const SpanKey3D& rExternalSpan,
//             const std::string& rCallerInfo) -> Node::Pointer
//     {
//         if (!span_has_skin_node(rExternalSpan)) {
//             bool found_intersecting_triangle = false;
//             array_1d<double, 3> auxiliary_point = ZeroVector(3);
    
//             for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
//                 const auto& r_geometry = r_condition.GetGeometry();
    
//                 if (find_triangle_box_intersection_point(
//                         r_geometry,
//                         rExternalSpan,
//                         auxiliary_point)) {
//                     found_intersecting_triangle = true;
//                     break;
//                 }
//             }
    
//             KRATOS_ERROR_IF_NOT(found_intersecting_triangle)
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: external span "
//                 << span_to_string(rExternalSpan)
//                 << " has no skin node and no intersecting skin triangle was found."
//                 << " caller_info=" << rCallerInfo
//                 << " query_point=" << rPoint
//                 << std::endl;
    
//             clamp_point_inside_span(auxiliary_point, rExternalSpan);
    
//             const IndexType new_node_id = get_next_auxiliary_skin_node_id();
    
//             auto p_auxiliary_node = rSkinSubModelPart.CreateNewNode(
//                 new_node_id,
//                 auxiliary_point[0],
//                 auxiliary_point[1],
//                 auxiliary_point[2]);
    
//             return p_auxiliary_node;
//         }
    
//         return find_closest_skin_node_in_external_span(rPoint, rExternalSpan);
//     };

//     auto distinct_neighbours = [](const std::vector<Geometry<Node>::Pointer>& rNeighbours) {
//         std::vector<Geometry<Node>::Pointer> neighbours;
//         std::set<IndexType> ids;
//         for (const auto& p_neighbour : rNeighbours) {
//             if (p_neighbour && ids.insert(p_neighbour->Id()).second) {
//                 neighbours.push_back(p_neighbour);
//             }
//         }
//         return neighbours;
//     };

//     auto choose_neighbour_geometry = [&](const std::vector<Geometry<Node>::Pointer>& rCandidates,
//                                          const array_1d<double, 3>& rCentroid) {
//         const auto candidates = distinct_neighbours(rCandidates);
//         KRATOS_ERROR_IF(candidates.empty())
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: no candidate NEIGHBOUR_GEOMETRY was available." << std::endl;

//         Geometry<Node>::Pointer p_best_neighbour;
//         double best_distance_squared = std::numeric_limits<double>::max();
//         for (const auto& p_neighbour : candidates) {
//             const array_1d<double, 3> diff = p_neighbour->Center() - rCentroid;
//             const double distance_squared = inner_prod(diff, diff);
//             if (distance_squared < best_distance_squared) {
//                 best_distance_squared = distance_squared;
//                 p_best_neighbour = p_neighbour;
//             }
//         }
//         return p_best_neighbour;
//     };

//     auto scan_volume = [&](const NurbsVolumeType& rVolume, const std::size_t Order) {
//         VolumeDiagnostics diagnostics;
//         IntegrationPointsArrayType integration_points(Order * Order * Order);
//         auto integration_point_it = integration_points.begin();
//         IntegrationPointUtilities::IntegrationPoints3D(
//             integration_point_it,
//             Order,
//             Order,
//             Order,
//             0.0,
//             1.0,
//             0.0,
//             1.0,
//             0.0,
//             1.0);

//         for (const auto& r_integration_point : integration_points) {
//             CoordinatesArrayType local_coordinates = ZeroVector(3);
//             local_coordinates[0] = r_integration_point[0];
//             local_coordinates[1] = r_integration_point[1];
//             local_coordinates[2] = r_integration_point[2];

//             Matrix jacobian;
//             rVolume.Jacobian(jacobian, local_coordinates);
//             const double det_jacobian = MathUtils<double>::Det(jacobian);

//             diagnostics.MinDeterminant = std::min(diagnostics.MinDeterminant, det_jacobian);
//             diagnostics.MaxDeterminant = std::max(diagnostics.MaxDeterminant, det_jacobian);
//             if (det_jacobian > determinant_tol) {
//                 ++diagnostics.PositiveDeterminants;
//             } else if (det_jacobian < -determinant_tol) {
//                 ++diagnostics.NegativeDeterminants;
//             } else {
//                 ++diagnostics.NearZeroDeterminants;
//             }
//             diagnostics.WeightSum += det_jacobian * r_integration_point.Weight();
//         }
//         return diagnostics;
//     };

//     auto swap_volume_parametric_directions = [](std::array<Node::Pointer, 8>& rNodes) {
//         std::swap(rNodes[1], rNodes[2]);
//         std::swap(rNodes[5], rNodes[6]);
//     };

//     std::vector<GapVolumeData> gap_volumes;
//     std::map<LateralFaceKey, LateralFaceData> lateral_face_registry;
//     std::map<std::pair<IndexType, IndexType>, SurrogateEdgeData> surrogate_edges;
//     std::map<IndexType, SurrogateNodeData> surrogate_nodes;
//     std::map<IndexType, std::vector<std::size_t>> face_indices_by_node;
//     std::set<SpanKey3D> all_external_spans;
//     std::set<SpanKey3D> type1_external_spans;
//     std::set<std::pair<IndexType, SpanKey3D>> type3_projection_keys;

//     std::size_t number_of_external_spans = 0;
//     std::size_t number_of_type2_projection_faces = 0;
//     std::size_t number_of_type3_projection_segments = 0;
//     std::size_t number_of_top_conditions = 0;
//     std::size_t number_of_collapsed_top_faces = 0;
//     std::size_t number_of_interface_conditions = 0;
//     std::size_t number_of_non_manifold_lateral_faces = 0;
//     double min_volume_weight = std::numeric_limits<double>::max();
//     double max_volume_weight = 0.0;
//     double min_determinant = std::numeric_limits<double>::max();
//     std::array<GapTypeDiagnostics, 4> gap_type_diagnostics;

//     KRATOS_ERROR_IF_NOT(mpGapConditionsSubModelPart)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: mpGapConditionsSubModelPart is not initialized." << std::endl;
//     KRATOS_ERROR_IF(mGapConditionName.empty())
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: mGapConditionName is empty." << std::endl;

//     ModelPart& r_top_conditions_model_part = *mpGapConditionsSubModelPart;
//     PropertiesPointerType p_top_condition_properties = r_top_conditions_model_part.pGetProperties(0);
//     KRATOS_ERROR_IF_NOT(p_top_condition_properties)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to get properties 0 for gap conditions."
//         << std::endl;

//     auto register_lateral_face = [&](const std::vector<Node::Pointer>& rOrderedNodes,
//                                      const std::vector<int>& rNodeKinds,
//                                      const Geometry<Node>::Pointer& pNeighbourGeometry,
//                                      const IndexType GapVolumeId) {
//         KRATOS_ERROR_IF(rOrderedNodes.size() != rNodeKinds.size())
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: inconsistent lateral face node-kind data." << std::endl;
//         KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: lateral face registered without NEIGHBOUR_GEOMETRY." << std::endl;

//         std::vector<std::pair<Node::Pointer, int>> key_nodes;
//         key_nodes.reserve(rOrderedNodes.size());
//         for (std::size_t i = 0; i < rOrderedNodes.size(); ++i) {
//             key_nodes.push_back({rOrderedNodes[i], rNodeKinds[i]});
//         }
//         const LateralFaceKey key = make_face_key(key_nodes);

//         auto registry_it = lateral_face_registry.find(key);
//         if (registry_it == lateral_face_registry.end()) {
//             LateralFaceData data;
//             data.Key = key;
//             data.IsTriangle = key.Vertices.size() == 3;
//             data.IsDegenerate = key.Vertices.size() < 3 || compute_face_area(rOrderedNodes) <= area_tol;
//             if (!data.IsDegenerate) {
//                 data.pGeometry = create_linear_surface(rOrderedNodes);
//             }
//             registry_it = lateral_face_registry.emplace(key, std::move(data)).first;
//         }

//         auto& r_data = registry_it->second;
//         r_data.AdjacentGapVolumeIds.push_back(GapVolumeId);
//         auto neighbours = r_data.NeighbourGeometries;
//         neighbours.push_back(pNeighbourGeometry);
//         r_data.NeighbourGeometries = distinct_neighbours(neighbours);
//         return key;
//     };

//     auto create_top_condition = [&](const int GapType,
//                                     const NurbsSurfaceType::Pointer& pSurface,
//                                     const std::vector<Node::Pointer>& rTopNodes,
//                                     const Geometry<Node>::Pointer& pNeighbourGeometry) {
//         KRATOS_ERROR_IF_NOT(pSurface)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: top surface is null." << std::endl;
//         KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: top condition requested without NEIGHBOUR_GEOMETRY." << std::endl;

//         IntegrationPointsArrayType surface_integration_points;
//         GeometriesArrayType surface_quadrature_point_list;
//         const std::size_t surface_integration_order = static_cast<std::size_t>(2 * mGapInterpolationOrder + 1);
//         pSurface->CreateIntegrationPoints(
//             surface_integration_points,
//             surface_integration_order,
//             surface_integration_order);

//         double surface_weight_sum = 0.0;
//         for (auto& r_integration_point : surface_integration_points) {
//             const double determinant_jacobian = pSurface->DeterminantOfJacobian(r_integration_point);
//             r_integration_point.SetWeight(r_integration_point.Weight() * determinant_jacobian);
//             surface_weight_sum += r_integration_point.Weight();
//         }
//         if (surface_weight_sum <= area_tol) {
//             return false;
//         }

//         IntegrationInfo surface_integration_info = pSurface->GetDefaultIntegrationInfo();
//         surface_integration_info.SetNumberOfIntegrationPointsPerSpan(0, 2 * surface_integration_order + 1);
//         surface_integration_info.SetNumberOfIntegrationPointsPerSpan(1, 2 * surface_integration_order + 1);

//         pSurface->CreateQuadraturePointGeometries(
//             surface_quadrature_point_list,
//             number_of_shape_functions_derivatives,
//             surface_integration_points,
//             surface_integration_info);
//         if (surface_quadrature_point_list.size() == 0) {
//             return false;
//         }
//         for (auto geometry_it = surface_quadrature_point_list.ptr_begin();
//              geometry_it != surface_quadrature_point_list.ptr_end();
//              ++geometry_it) {
//             KRATOS_ERROR_IF_NOT(*geometry_it)
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null top quadrature geometry."
//                 << " gap_type=" << GapType
//                 << " condition_name=" << mGapConditionName
//                 << " surface_weight_sum=" << surface_weight_sum
//                 << "." << std::endl;
//         }

//         std::size_t id_condition = 1;
//         if (r_top_conditions_model_part.GetRootModelPart().Conditions().size() > 0) {
//             id_condition = r_top_conditions_model_part.GetRootModelPart().Conditions().back().Id() + 1;
//         }

//         this->CreateConditions(
//             surface_quadrature_point_list.ptr_begin(),
//             surface_quadrature_point_list.ptr_end(),
//             r_top_conditions_model_part,
//             mGapConditionName,
//             id_condition,
//             p_top_condition_properties,
//             knot_span_sizes,
//             std::vector<Geometry<Node>::Pointer>{pNeighbourGeometry},
//             characteristic_length(rTopNodes));
//         number_of_top_conditions += surface_quadrature_point_list.size();
//         gap_type_diagnostics[GapType].TopConditions += surface_quadrature_point_list.size();
//         return true;
//     };

//     auto is_volume_candidate_sane = [&](std::array<Node::Pointer, 8> volume_nodes,
//                                         const std::size_t GapType,
//                                         const SpanKey3D& rDiagnosticSpan,
//                                         std::string& rFailureReason) {
//         const std::size_t integration_order = mGapInterpolationOrder + 1;
//         auto p_gap_volume = create_linear_volume(volume_nodes);
//         VolumeDiagnostics diagnostics = scan_volume(*p_gap_volume, integration_order);
//         if (diagnostics.HasConsistentNegativeSign()) {
//             swap_volume_parametric_directions(volume_nodes);
//             p_gap_volume = create_linear_volume(volume_nodes);
//             diagnostics = scan_volume(*p_gap_volume, integration_order);
//         }

//         if (diagnostics.HasMixedSigns() ||
//             diagnostics.HasConsistentNegativeSign() ||
//             diagnostics.WeightSum <= volume_tol) {
//             std::ostringstream message;
//             message << "type=" << GapType
//                     << " span=" << span_to_string(rDiagnosticSpan)
//                     << " min_det=" << diagnostics.MinDeterminant
//                     << " max_det=" << diagnostics.MaxDeterminant
//                     << " weight_sum=" << diagnostics.WeightSum
//                     << " positive_det=" << diagnostics.PositiveDeterminants
//                     << " negative_det=" << diagnostics.NegativeDeterminants
//                     << " near_zero_det=" << diagnostics.NearZeroDeterminants;
//             rFailureReason = message.str();
//             return false;
//         }
//         return true;
//     };

//     auto create_gap_volume = [&](const int GapType,
//                                  const SpanKey3D& rDiagnosticSpan,
//                                  std::array<Node::Pointer, 8> volume_nodes,
//                                  const std::vector<Node::Pointer>& rCharacteristicNodes,
//                                  const Geometry<Node>::Pointer& pNeighbourGeometry,
//                                  const std::string& rDiagnosticLabel) -> std::optional<IndexType>{
//         KRATOS_ERROR_IF_NOT(pNeighbourGeometry)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: gap volume has no NEIGHBOUR_GEOMETRY." << std::endl;

//         const std::size_t integration_order = mGapInterpolationOrder + 1;
//         auto p_gap_volume = create_linear_volume(volume_nodes);
//         VolumeDiagnostics diagnostics = scan_volume(*p_gap_volume, integration_order);
//         if (diagnostics.HasConsistentNegativeSign()) {
//             swap_volume_parametric_directions(volume_nodes);
//             p_gap_volume = create_linear_volume(volume_nodes);
//             diagnostics = scan_volume(*p_gap_volume, integration_order);
//         }

//         if (diagnostics.WeightSum < volume_tol) 
//             return std::nullopt;
            
//         KRATOS_ERROR_IF(diagnostics.HasMixedSigns() ||
//                         diagnostics.HasConsistentNegativeSign() ||
//                         diagnostics.WeightSum <= volume_tol)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: invalid gap volume."
//             << " type=" << GapType
//             << " span=" << span_to_string(rDiagnosticSpan)
//             << " label=" << rDiagnosticLabel
//             << " neighbour_geometry_id=" << pNeighbourGeometry->Id()
//             << " min_det=" << diagnostics.MinDeterminant
//             << " max_det=" << diagnostics.MaxDeterminant
//             << " weight_sum=" << diagnostics.WeightSum
//             << " positive_det=" << diagnostics.PositiveDeterminants
//             << " negative_det=" << diagnostics.NegativeDeterminants
//             << " near_zero_det=" << diagnostics.NearZeroDeterminants
//             << "volume_nodes=[" << volume_nodes[0]->Coordinates() << "," << volume_nodes[1]->Coordinates() << "," << volume_nodes[2]->Coordinates() << "," << volume_nodes[3]->Coordinates()
//              << "," << volume_nodes[4]->Coordinates() << "," << volume_nodes[5]->Coordinates() << "," << volume_nodes[6]->Coordinates() << "," << volume_nodes[7]->Coordinates()
//             << "." << std::endl;

//         IntegrationPointsArrayType volume_integration_points(integration_order * integration_order * integration_order);
//         auto integration_point_it = volume_integration_points.begin();
//         IntegrationPointUtilities::IntegrationPoints3D(
//             integration_point_it,
//             integration_order,
//             integration_order,
//             integration_order,
//             0.0,
//             1.0,
//             0.0,
//             1.0,
//             0.0,
//             1.0);

//         double volume_weight_sum = 0.0;
//         for (auto& r_integration_point : volume_integration_points) {
//             CoordinatesArrayType local_coordinates = ZeroVector(3);
//             local_coordinates[0] = r_integration_point[0];
//             local_coordinates[1] = r_integration_point[1];
//             local_coordinates[2] = r_integration_point[2];

//             Matrix jacobian;
//             p_gap_volume->Jacobian(jacobian, local_coordinates);
//             const double det_jacobian = MathUtils<double>::Det(jacobian);
//             KRATOS_ERROR_IF(det_jacobian < -determinant_tol)
//                 << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: negative determinant after orientation repair."
//                 << " type=" << GapType
//                 << " span=" << span_to_string(rDiagnosticSpan)
//                 << " det=" << det_jacobian << "." << std::endl;
//             const double positive_det_jacobian =
//                 std::abs(det_jacobian) <= determinant_tol ? 0.0 : det_jacobian;
//             const double weight = r_integration_point.Weight() * positive_det_jacobian;
//             r_integration_point.SetWeight(weight);
//             volume_weight_sum += weight;
//         }

//         KRATOS_ERROR_IF(volume_weight_sum <= volume_tol)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: volume integration weight is not positive."
//             << " type=" << GapType
//             << " span=" << span_to_string(rDiagnosticSpan)
//             << " weight_sum=" << volume_weight_sum << "." << std::endl;

//         store_gap_debug_faces(GapType, volume_nodes);

//         GeometriesArrayType volume_quadrature_point_list;
//         IntegrationInfo volume_integration_info(
//             {integration_order, integration_order, integration_order},
//             {IntegrationInfo::QuadratureMethod::CUSTOM,
//              IntegrationInfo::QuadratureMethod::CUSTOM,
//              IntegrationInfo::QuadratureMethod::CUSTOM});

//         p_gap_volume->CreateQuadraturePointGeometries(
//             volume_quadrature_point_list,
//             number_of_shape_functions_derivatives,
//             volume_integration_points,
//             volume_integration_info);
//         KRATOS_ERROR_IF(volume_quadrature_point_list.size() == 0)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to create volume quadrature geometries."
//             << " type=" << GapType
//             << " span=" << span_to_string(rDiagnosticSpan) << "." << std::endl;

//         IndexType id_element = 1;
//         if (mpGapElementsSubModelPart->GetRootModelPart().Elements().size() > 0) {
//             id_element = mpGapElementsSubModelPart->GetRootModelPart().Elements().back().Id() + 1;
//         }

//         this->CreateElements(
//             volume_quadrature_point_list.ptr_begin(),
//             volume_quadrature_point_list.ptr_end(),
//             *mpGapElementsSubModelPart,
//             mGapElementName,
//             id_element,
//             PropertiesPointerType(),
//             std::vector<Geometry<Node>::Pointer>{pNeighbourGeometry},
//             characteristic_length(rCharacteristicNodes));

//         min_volume_weight = std::min(min_volume_weight, volume_weight_sum);
//         max_volume_weight = std::max(max_volume_weight, volume_weight_sum);
//         min_determinant = std::min(min_determinant, diagnostics.MinDeterminant);
//         gap_type_diagnostics[GapType].UpdateVolume(volume_weight_sum, diagnostics.MinDeterminant);

//         GapVolumeData gap_volume_data;
//         gap_volume_data.Type = GapType;
//         gap_volume_data.pVolumeGeometry = p_gap_volume;
//         gap_volume_data.pChosenNeighbourGeometry = pNeighbourGeometry;
//         gap_volumes.push_back(gap_volume_data);
//         return static_cast<IndexType>(gap_volumes.size());
//     };

//     struct BrepPatchData
//     {
//         GeometryType::Pointer pBrepGeometry;
//         BrepSurfaceOnVolumeType::Pointer pBrepSurface;
//         CoordinatesArrayType Vertex00 = ZeroVector(3);
//         CoordinatesArrayType Vertex01 = ZeroVector(3);
//         CoordinatesArrayType Vertex10 = ZeroVector(3);
//         CoordinatesArrayType Vertex11 = ZeroVector(3);
//         CoordinatesArrayType MiddlePoint = ZeroVector(3);
//         CoordinatesArrayType MiddlePointLocalCoordinates = ZeroVector(3);
//     };

//     IndexType starting_brep_id = 2;
//     std::size_t size_surrogate_loop = rSurrogateSubModelPart.NumberOfConditions();
//     if constexpr (TIsInnerLoop) {
//         if (rSurrogateSubModelPart.NumberOfElements() > 0) {
//             const IndexType element_id = rSurrogateSubModelPart.ElementsBegin()->Id();
//             const IndexType first_condition_id = rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[0].Id();
//             const IndexType last_condition_id = rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[1].Id();
//             size_surrogate_loop = last_condition_id - first_condition_id + 1;
//         }
//         starting_brep_id = mpIgaModelPart->HasSubModelPart("surrogate_outer") &&
//                            mpIgaModelPart->GetSubModelPart("surrogate_outer").NumberOfConditions() > 0
//             ? 2 + mpIgaModelPart->GetSubModelPart("surrogate_outer").NumberOfConditions()
//             : 8;
//     }

//     auto compute_brep_patch_data = [&](const IndexType BrepId) {
//         BrepPatchData brep_patch_data;
//         brep_patch_data.pBrepGeometry = mpIgaModelPart->pGetGeometry(BrepId);
//         KRATOS_ERROR_IF_NOT(brep_patch_data.pBrepGeometry)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: geometry with id "
//             << BrepId << " was not found." << std::endl;

//         brep_patch_data.pBrepSurface =
//             std::dynamic_pointer_cast<BrepSurfaceOnVolumeType>(brep_patch_data.pBrepGeometry);
//         KRATOS_ERROR_IF_NOT(brep_patch_data.pBrepSurface)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: geometry with id "
//             << brep_patch_data.pBrepGeometry->Id() << " is not a BrepSurfaceOnVolumeType." << std::endl;

//         NurbsInterval brep_domain_interval_u = brep_patch_data.pBrepSurface->DomainIntervalU();
//         NurbsInterval brep_domain_interval_v = brep_patch_data.pBrepSurface->DomainIntervalV();

//         CoordinatesArrayType vertex_00_local_coords = ZeroVector(3);
//         CoordinatesArrayType vertex_01_local_coords = ZeroVector(3);
//         CoordinatesArrayType vertex_10_local_coords = ZeroVector(3);
//         CoordinatesArrayType vertex_11_local_coords = ZeroVector(3);

//         vertex_00_local_coords[0] = brep_domain_interval_u.GetT0();
//         vertex_00_local_coords[1] = brep_domain_interval_v.GetT0();
//         vertex_01_local_coords[0] = brep_domain_interval_u.GetT0();
//         vertex_01_local_coords[1] = brep_domain_interval_v.GetT1();
//         vertex_10_local_coords[0] = brep_domain_interval_u.GetT1();
//         vertex_10_local_coords[1] = brep_domain_interval_v.GetT0();
//         vertex_11_local_coords[0] = brep_domain_interval_u.GetT1();
//         vertex_11_local_coords[1] = brep_domain_interval_v.GetT1();

//         brep_patch_data.pBrepSurface->GlobalCoordinates(brep_patch_data.Vertex00, vertex_00_local_coords);
//         brep_patch_data.pBrepSurface->GlobalCoordinates(brep_patch_data.Vertex01, vertex_01_local_coords);
//         brep_patch_data.pBrepSurface->GlobalCoordinates(brep_patch_data.Vertex10, vertex_10_local_coords);
//         brep_patch_data.pBrepSurface->GlobalCoordinates(brep_patch_data.Vertex11, vertex_11_local_coords);

//         brep_patch_data.MiddlePointLocalCoordinates[0] =
//             0.5 * (brep_domain_interval_u.GetT0() + brep_domain_interval_u.GetT1());
//         brep_patch_data.MiddlePointLocalCoordinates[1] =
//             0.5 * (brep_domain_interval_v.GetT0() + brep_domain_interval_v.GetT1());
//         brep_patch_data.pBrepSurface->GlobalCoordinates(
//             brep_patch_data.MiddlePoint,
//             brep_patch_data.MiddlePointLocalCoordinates);
//         return brep_patch_data;
//     };

//     std::vector<BrepPatchData> brep_patch_data_list;
//     brep_patch_data_list.reserve(size_surrogate_loop);
//     for (std::size_t i = 0; i < size_surrogate_loop; ++i) {
//         brep_patch_data_list.push_back(compute_brep_patch_data(starting_brep_id + i));
//     }

//     auto find_brep_patch_matching_condition = [&](const Condition& rCondition) -> const BrepPatchData& {
//         const auto& r_condition_geometry = rCondition.GetGeometry();
//         for (const auto& r_brep_patch_data : brep_patch_data_list) {
//             const array_1d<double, 3> center_difference =
//                 r_condition_geometry.Center() - r_brep_patch_data.MiddlePoint;
//             if (inner_prod(center_difference, center_difference) <= distance_tol * distance_tol) {
//                 return r_brep_patch_data;
//             }
//         }

//         KRATOS_ERROR << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to find BREP patch matching surrogate condition #"
//             << rCondition.Id() << " center=" << r_condition_geometry.Center() << "." << std::endl;
//         return brep_patch_data_list.front();
//     };

//     auto match_condition_node = [&](const Condition& rCondition,
//                                     const CoordinatesArrayType& rVertex,
//                                     const char* pVertexName) {
//         const auto& r_geometry = rCondition.GetGeometry();
//         double best_distance_squared = std::numeric_limits<double>::max();
//         std::size_t best_index = r_geometry.PointsNumber();

//         for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
//             const array_1d<double, 3> diff = r_geometry[i].Coordinates() - rVertex;
//             const double distance_squared = inner_prod(diff, diff);
//             if (distance_squared < best_distance_squared) {
//                 best_distance_squared = distance_squared;
//                 best_index = i;
//             }
//         }

//         KRATOS_ERROR_IF(best_index == r_geometry.PointsNumber() || best_distance_squared > distance_tol * distance_tol)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to match surrogate vertex "
//             << pVertexName << " to condition #" << rCondition.Id()
//             << " vertex=" << rVertex
//             << " best_distance=" << std::sqrt(best_distance_squared) << "." << std::endl;
//         return r_geometry.pGetPoint(best_index);
//     };

//     auto create_face_neighbour_geometry = [&](const BrepPatchData& rBrepPatchData) {
//         IntegrationPoint<2> integration_point(
//             rBrepPatchData.MiddlePointLocalCoordinates[0],
//             rBrepPatchData.MiddlePointLocalCoordinates[1],
//             1.0);
//         IntegrationPointsArrayType surrogate_integration_points;
//         surrogate_integration_points.push_back(integration_point);

//         IntegrationInfo integration_info = rBrepPatchData.pBrepSurface->GetDefaultIntegrationInfo();
//         GeometriesArrayType quadrature_point_list;
//         rBrepPatchData.pBrepSurface->CreateQuadraturePointGeometries(
//             quadrature_point_list,
//             number_of_shape_functions_derivatives,
//             surrogate_integration_points,
//             integration_info);

//         KRATOS_ERROR_IF(quadrature_point_list.size() == 0)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to create face centre NEIGHBOUR_GEOMETRY." << std::endl;
//         return quadrature_point_list(0);
//     };

//     std::vector<SurrogateFaceData> surrogate_faces;
//     surrogate_faces.reserve(rSurrogateSubModelPart.NumberOfConditions());

//     for (const auto& r_surrogate_condition : rSurrogateSubModelPart.Conditions()) {
//         KRATOS_ERROR_IF(r_surrogate_condition.GetGeometry().PointsNumber() != 4)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: surrogate condition #"
//             << r_surrogate_condition.Id() << " must have exactly four nodes." << std::endl;

//         const auto& r_brep_patch_data = find_brep_patch_matching_condition(r_surrogate_condition);
//         SurrogateFaceData face_data;
//         face_data.pCondition = &r_surrogate_condition;
//         face_data.Nodes = {{
//             match_condition_node(r_surrogate_condition, r_brep_patch_data.Vertex00, "00"),
//             match_condition_node(r_surrogate_condition, r_brep_patch_data.Vertex10, "10"),
//             match_condition_node(r_surrogate_condition, r_brep_patch_data.Vertex11, "11"),
//             match_condition_node(r_surrogate_condition, r_brep_patch_data.Vertex01, "01")}};

//         Vector normal_vector = r_surrogate_condition.GetValue(NORMAL);
//         KRATOS_ERROR_IF(normal_vector.size() < 3)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: NORMAL on surrogate condition #"
//             << r_surrogate_condition.Id() << " must have three components." << std::endl;
//         for (int d = 0; d < 3; ++d) {
//             face_data.OutwardNormal[d] = normal_vector[d];
//         }
//         const double normal_norm = norm_2(face_data.OutwardNormal);
//         KRATOS_ERROR_IF(normal_norm <= 1.0e-16)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: zero NORMAL on surrogate condition #"
//             << r_surrogate_condition.Id() << "." << std::endl;
//         face_data.OutwardNormal /= normal_norm;

//         face_data.NormalAxis = 0;
//         for (int d = 1; d < 3; ++d) {
//             if (std::abs(face_data.OutwardNormal[d]) > std::abs(face_data.OutwardNormal[face_data.NormalAxis])) {
//                 face_data.NormalAxis = d;
//             }
//         }
//         face_data.NormalSign = face_data.OutwardNormal[face_data.NormalAxis] >= 0.0 ? 1 : -1;
//         face_data.pNeighbourGeometry = create_face_neighbour_geometry(r_brep_patch_data);

//         for (std::size_t i = 0; i < face_data.Nodes.size(); ++i) {
//             face_data.GridNodes[i] = compute_grid_index(face_data.Nodes[i]->Coordinates());
//         }

//         const int fixed_grid_index =
//             face_data.NormalAxis == 0 ? face_data.GridNodes[0].I :
//             face_data.NormalAxis == 1 ? face_data.GridNodes[0].J :
//                                         face_data.GridNodes[0].K;
//         face_data.ExternalSpan.I = std::min(std::min(face_data.GridNodes[0].I, face_data.GridNodes[1].I),
//                                             std::min(face_data.GridNodes[2].I, face_data.GridNodes[3].I));
//         face_data.ExternalSpan.J = std::min(std::min(face_data.GridNodes[0].J, face_data.GridNodes[1].J),
//                                             std::min(face_data.GridNodes[2].J, face_data.GridNodes[3].J));
//         face_data.ExternalSpan.K = std::min(std::min(face_data.GridNodes[0].K, face_data.GridNodes[1].K),
//                                             std::min(face_data.GridNodes[2].K, face_data.GridNodes[3].K));
//         add_to_span(face_data.ExternalSpan, face_data.NormalAxis, face_data.NormalSign > 0 ? fixed_grid_index - get_span_component(face_data.ExternalSpan, face_data.NormalAxis) : fixed_grid_index - 1 - get_span_component(face_data.ExternalSpan, face_data.NormalAxis));
//         face_data.ActiveSpan = face_data.ExternalSpan;
//         add_to_span(face_data.ActiveSpan, face_data.NormalAxis, -face_data.NormalSign);

//         if (!is_span_inside_domain(face_data.ExternalSpan)) {
//             KRATOS_INFO_IF("CreateSbmExtendedGeometries3D", mEchoLevel > 1)
//                 << "Skipping surrogate face condition #" << r_surrogate_condition.Id()
//                 << " because its external span " << span_to_string(face_data.ExternalSpan)
//                 << " is outside the background grid." << std::endl;
//             continue;
//         }

//         type1_external_spans.insert(face_data.ExternalSpan);
//         all_external_spans.insert(face_data.ExternalSpan);
//         const std::size_t face_index = surrogate_faces.size();
//         for (const auto& p_node : face_data.Nodes) {
//             auto& r_node_data = surrogate_nodes[p_node->Id()];
//             r_node_data.pNode = p_node;
//             r_node_data.Grid = compute_grid_index(p_node->Coordinates());
//             r_node_data.IncidentFaceNeighbourGeometries.push_back(face_data.pNeighbourGeometry);
//             face_indices_by_node[p_node->Id()].push_back(face_index);
//         }
//         surrogate_faces.push_back(face_data);
//     }

//     ensure_skin_nodes_in_external_spans(type1_external_spans);

//     skin_bins_per_knot_span =
//         CreateSkinBinsPerKnotSpanMatrix3D(rSkinSubModelPart, rSurrogateSubModelPart);


//     auto edge_key = [](const Node::Pointer& pA, const Node::Pointer& pB) {
//         return std::make_pair(std::min(pA->Id(), pB->Id()), std::max(pA->Id(), pB->Id()));
//     };

//     auto add_projection_face_to_edge = [&](const Node::Pointer& pS0,
//                                            const Node::Pointer& pS1,
//                                            const Node::Pointer& pK0,
//                                            const Node::Pointer& pK1,
//                                            const SpanKey3D& rExternalSpan,
//                                            const SpanKey3D& rActiveSpan,
//                                            const std::vector<Geometry<Node>::Pointer>& rCandidateNeighbours,
//                                            const int NormalAxis,
//                                            const int NormalSign,
//                                            const bool IsType1) {
//         auto key = edge_key(pS0, pS1);
//         auto& r_edge_data = surrogate_edges[key];
//         const bool local_matches_canonical = pS0->Id() <= pS1->Id();
//         r_edge_data.pS0 = local_matches_canonical ? pS0 : pS1;
//         r_edge_data.pS1 = local_matches_canonical ? pS1 : pS0;

//         ProjectionFaceData projection_face;
//         projection_face.pS0 = r_edge_data.pS0;
//         projection_face.pS1 = r_edge_data.pS1;
//         projection_face.pK0 = local_matches_canonical ? pK0 : pK1;
//         projection_face.pK1 = local_matches_canonical ? pK1 : pK0;
//         projection_face.ExternalSpan = rExternalSpan;
//         projection_face.ActiveSpan = rActiveSpan;
//         projection_face.CandidateNeighbourGeometries = distinct_neighbours(rCandidateNeighbours);
//         projection_face.NormalAxis = NormalAxis;
//         projection_face.NormalSign = NormalSign;
//         projection_face.IsType1 = IsType1;

//         if (r_edge_data.ExternalSpanIds.insert(rExternalSpan).second) {
//             r_edge_data.ProjectionFaces.push_back(projection_face);
//         }
//     };

//     auto add_projection_segment_to_node = [&](const Node::Pointer& pSurrogateNode,
//         const Node::Pointer& pSkinNode,
//         const SpanKey3D& rExternalSpan) {
//         KRATOS_ERROR_IF_NOT(pSurrogateNode)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null surrogate node in projection segment."
//         << std::endl;
//         KRATOS_ERROR_IF_NOT(pSkinNode)
//         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: null skin node in projection segment."
//         << std::endl;

//         auto& r_node_data = surrogate_nodes[pSurrogateNode->Id()];
//         r_node_data.pNode = pSurrogateNode;
//         r_node_data.Grid = compute_grid_index(pSurrogateNode->Coordinates());

//         for (const auto& r_existing_projection : r_node_data.ProjectionSegments) {
//         if (r_existing_projection.ExternalSpan == rExternalSpan) {
//         return;
//         }
//         }

//         LocalProjectionData projection_data;
//         projection_data.pSurrogateNode = pSurrogateNode;
//         projection_data.pSkinNode = pSkinNode;
//         projection_data.ExternalSpan = rExternalSpan;
//         projection_data.SignTriple = span_sign_around_node(rExternalSpan, r_node_data.Grid);
//         projection_data.Vector = pSkinNode->Coordinates() - pSurrogateNode->Coordinates();

//         r_node_data.ProjectionSegments.push_back(projection_data);
//     };

//     for (const auto& r_face : surrogate_faces) {
//         const std::array<Node::Pointer, 4> original_skin_nodes{{get_or_create_closest_skin_node_in_external_span(
//                                                                         r_face.Nodes[0]->Coordinates(),
//                                                                         r_face.ExternalSpan,
//                                                                         "type1 face_condition_id=" + std::to_string(r_face.pCondition->Id()) +
//                                                                         " normal_axis=" + std::to_string(r_face.NormalAxis) +
//                                                                         " normal_sign=" + std::to_string(r_face.NormalSign) +
//                                                                         " active_span=" + span_to_string(r_face.ActiveSpan)),
//                                                                 get_or_create_closest_skin_node_in_external_span(
//                                                                     r_face.Nodes[1]->Coordinates(),
//                                                                     r_face.ExternalSpan,
//                                                                     "type1 face_condition_id=" + std::to_string(r_face.pCondition->Id()) +
//                                                                     " normal_axis=" + std::to_string(r_face.NormalAxis) +
//                                                                     " normal_sign=" + std::to_string(r_face.NormalSign) +
//                                                                     " active_span=" + span_to_string(r_face.ActiveSpan)),
//                                                                 get_or_create_closest_skin_node_in_external_span(
//                                                                     r_face.Nodes[2]->Coordinates(),
//                                                                     r_face.ExternalSpan,
//                                                                     "type1 face_condition_id=" + std::to_string(r_face.pCondition->Id()) +
//                                                                     " normal_axis=" + std::to_string(r_face.NormalAxis) +
//                                                                     " normal_sign=" + std::to_string(r_face.NormalSign) +
//                                                                     " active_span=" + span_to_string(r_face.ActiveSpan)),
//                                                                 get_or_create_closest_skin_node_in_external_span(
//                                                                     r_face.Nodes[3]->Coordinates(),
//                                                                     r_face.ExternalSpan,
//                                                                     "type1 face_condition_id=" + std::to_string(r_face.pCondition->Id()) +
//                                                                     " normal_axis=" + std::to_string(r_face.NormalAxis) +
//                                                                     " normal_sign=" + std::to_string(r_face.NormalSign) +
//                                                                     " active_span=" + span_to_string(r_face.ActiveSpan))}};

//         const std::array<Node::Pointer, 4> bottom_nodes{{
//             r_face.Nodes[0],
//             r_face.Nodes[1],
//             r_face.Nodes[2],
//             r_face.Nodes[3]}};

//         auto make_type1_volume_nodes = [&](const std::array<Node::Pointer, 4>& rTopNodes) {
//             return std::array<Node::Pointer, 8>{{
//                 r_face.Nodes[0],
//                 r_face.Nodes[1],
//                 r_face.Nodes[3],
//                 r_face.Nodes[2],
//                 rTopNodes[0],
//                 rTopNodes[1],
//                 rTopNodes[3],
//                 rTopNodes[2]}};
//         };

//         std::vector<std::array<Node::Pointer, 4>> top_candidates;
//         auto append_top_candidate = [&](const std::array<Node::Pointer, 4>& rCandidate) {
//             const bool already_present = std::any_of(
//                 top_candidates.begin(),
//                 top_candidates.end(),
//                 [&](const std::array<Node::Pointer, 4>& rExistingCandidate) {
//                     for (std::size_t i = 0; i < rCandidate.size(); ++i) {
//                         if (rCandidate[i]->Id() != rExistingCandidate[i]->Id()) {
//                             return false;
//                         }
//                     }
//                     return true;
//                 });
//             if (!already_present) {
//                 top_candidates.push_back(rCandidate);
//             }
//         };

//         append_top_candidate(original_skin_nodes);
//         std::string original_top_failure;
//         if (!is_top_quad_sane(bottom_nodes, original_skin_nodes, r_face.OutwardNormal, original_top_failure)) {
//             const auto problematic_corners =
//                 find_problematic_top_quad_corners(bottom_nodes, original_skin_nodes, r_face.OutwardNormal);

//             for (const std::size_t corner_index : problematic_corners) {
//                 auto candidate = original_skin_nodes;
//                 candidate[corner_index] = original_skin_nodes[(corner_index + 1) % original_skin_nodes.size()];
//                 append_top_candidate(candidate);

//                 candidate = original_skin_nodes;
//                 candidate[corner_index] =
//                     original_skin_nodes[(corner_index + original_skin_nodes.size() - 1) % original_skin_nodes.size()];
//                 append_top_candidate(candidate);
//             }

//             for (std::size_t i = 0; i < problematic_corners.size(); ++i) {
//                 for (std::size_t j = i + 1; j < problematic_corners.size(); ++j) {
//                     auto candidate = original_skin_nodes;
//                     const std::size_t first_corner = problematic_corners[i];
//                     const std::size_t second_corner = problematic_corners[j];
//                     candidate[first_corner] = original_skin_nodes[(first_corner + 1) % original_skin_nodes.size()];
//                     candidate[second_corner] =
//                         original_skin_nodes[(second_corner + original_skin_nodes.size() - 1) % original_skin_nodes.size()];
//                     append_top_candidate(candidate);
//                 }
//             }
//         }

//         const auto p_full_collapse_node = find_closest_skin_node_to_entity_in_external_span(
//             {r_face.Nodes[0], r_face.Nodes[1], r_face.Nodes[2], r_face.Nodes[3]},
//             r_face.ExternalSpan);
//         append_top_candidate({{
//             p_full_collapse_node,
//             p_full_collapse_node,
//             p_full_collapse_node,
//             p_full_collapse_node}});

//         std::array<Node::Pointer, 4> skin_nodes = original_skin_nodes;
//         bool selected_type1_candidate = false;
//         bool selected_top_is_sane_quad = false;
//         std::string last_candidate_failure;
//         for (const auto& r_candidate_top_nodes : top_candidates) {
//             std::string top_failure;
//             const bool candidate_is_sane_quad =
//                 is_top_quad_sane(bottom_nodes, r_candidate_top_nodes, r_face.OutwardNormal, top_failure);
//             std::set<IndexType> candidate_top_node_ids;
//             for (const auto& p_node : r_candidate_top_nodes) {
//                 candidate_top_node_ids.insert(p_node->Id());
//             }
//             const bool candidate_has_collapsed_corners = candidate_top_node_ids.size() < r_candidate_top_nodes.size();
//             const bool candidate_is_collapsed_fallback =
//                 !candidate_is_sane_quad &&
//                 candidate_has_collapsed_corners &&
//                 (compute_face_area({r_candidate_top_nodes[0], r_candidate_top_nodes[1], r_candidate_top_nodes[2], r_candidate_top_nodes[3]}) <= area_tol ||
//                  is_ordered_face_non_degenerate(
//                      {r_candidate_top_nodes[0], r_candidate_top_nodes[1], r_candidate_top_nodes[2], r_candidate_top_nodes[3]},
//                      &r_face.OutwardNormal,
//                      top_failure));
//             if (!candidate_is_sane_quad && !candidate_is_collapsed_fallback) {
//                 last_candidate_failure = top_failure;
//                 continue;
//             }

//             std::string volume_failure;
//             if (!is_volume_candidate_sane(
//                     make_type1_volume_nodes(r_candidate_top_nodes),
//                     1,
//                     r_face.ExternalSpan,
//                     volume_failure)) {
//                 last_candidate_failure = volume_failure;
//                 continue;
//             }

//             skin_nodes = r_candidate_top_nodes;
//             selected_top_is_sane_quad = candidate_is_sane_quad;
//             selected_type1_candidate = true;
//             break;
//         }

//         KRATOS_ERROR_IF_NOT(selected_type1_candidate)
//             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to create a sane type-1 top fallback."
//             << " face_condition_id=" << r_face.pCondition->Id()
//             << " external_span=" << span_to_string(r_face.ExternalSpan)
//             << " original_top_failure=" << original_top_failure
//             << " last_candidate_failure=" << last_candidate_failure
//             << "." << std::endl;

//         const bool collapsed_top_face = !selected_top_is_sane_quad;
//         if (collapsed_top_face) {
//             ++number_of_collapsed_top_faces;
//             ++gap_type_diagnostics[1].CollapsedTops;
//         }

//         const std::array<Node::Pointer, 8> volume_nodes = make_type1_volume_nodes(skin_nodes);

//         const auto gap_volume_id = create_gap_volume(
//             1,
//             r_face.ExternalSpan,
//             volume_nodes,
//             {r_face.Nodes[0], r_face.Nodes[1], r_face.Nodes[2], r_face.Nodes[3],
//             skin_nodes[0], skin_nodes[1], skin_nodes[2], skin_nodes[3]},
//             r_face.pNeighbourGeometry,
//             "type-1 face condition #" + std::to_string(r_face.pCondition->Id()));
        
//         if (!gap_volume_id) {
//             continue;
//         }

//         if (!collapsed_top_face) {
//             auto p_top_surface = create_linear_surface({skin_nodes[0], skin_nodes[1], skin_nodes[2], skin_nodes[3]});
//             create_top_condition(
//                 1,
//                 p_top_surface,
//                 {skin_nodes[0], skin_nodes[1], skin_nodes[2], skin_nodes[3]},
//                 r_face.pNeighbourGeometry);
//         }

//         register_lateral_face({r_face.Nodes[0], r_face.Nodes[1], skin_nodes[1], skin_nodes[0]}, {0, 0, 1, 1}, r_face.pNeighbourGeometry, *gap_volume_id);
//         register_lateral_face({r_face.Nodes[1], r_face.Nodes[2], skin_nodes[2], skin_nodes[1]}, {0, 0, 1, 1}, r_face.pNeighbourGeometry, *gap_volume_id);
//         register_lateral_face({r_face.Nodes[2], r_face.Nodes[3], skin_nodes[3], skin_nodes[2]}, {0, 0, 1, 1}, r_face.pNeighbourGeometry, *gap_volume_id);
//         register_lateral_face({r_face.Nodes[3], r_face.Nodes[0], skin_nodes[0], skin_nodes[3]}, {0, 0, 1, 1}, r_face.pNeighbourGeometry, *gap_volume_id);

//         add_projection_face_to_edge(r_face.Nodes[0], r_face.Nodes[1], skin_nodes[0], skin_nodes[1],
//                                     r_face.ExternalSpan, r_face.ActiveSpan, {r_face.pNeighbourGeometry},
//                                     r_face.NormalAxis, r_face.NormalSign, true);
//         add_projection_face_to_edge(r_face.Nodes[1], r_face.Nodes[2], skin_nodes[1], skin_nodes[2],
//                                     r_face.ExternalSpan, r_face.ActiveSpan, {r_face.pNeighbourGeometry},
//                                     r_face.NormalAxis, r_face.NormalSign, true);
//         add_projection_face_to_edge(r_face.Nodes[2], r_face.Nodes[3], skin_nodes[2], skin_nodes[3],
//                                     r_face.ExternalSpan, r_face.ActiveSpan, {r_face.pNeighbourGeometry},
//                                     r_face.NormalAxis, r_face.NormalSign, true);
//         add_projection_face_to_edge(r_face.Nodes[3], r_face.Nodes[0], skin_nodes[3], skin_nodes[0],
//                                     r_face.ExternalSpan, r_face.ActiveSpan, {r_face.pNeighbourGeometry},
//                                     r_face.NormalAxis, r_face.NormalSign, true);

//         add_projection_segment_to_node(r_face.Nodes[0], skin_nodes[0], r_face.ExternalSpan);
//         add_projection_segment_to_node(r_face.Nodes[1], skin_nodes[1], r_face.ExternalSpan);
//         add_projection_segment_to_node(r_face.Nodes[2], skin_nodes[2], r_face.ExternalSpan);
//         add_projection_segment_to_node(r_face.Nodes[3], skin_nodes[3], r_face.ExternalSpan);
//     }

//     return;


//     // TYPE 2 ELEMENTS -- CONSTRUCTION MISSING FACES FROM EDGE PAIRS

//     for (auto& r_edge_entry : surrogate_edges) {
//         auto& r_edge_data = r_edge_entry.second;
//         const std::size_t initial_projection_face_count = r_edge_data.ProjectionFaces.size();
//         for (std::size_t i = 0; i < initial_projection_face_count; ++i) {
//             const auto& r_first_face = r_edge_data.ProjectionFaces[i];
//             if (!r_first_face.IsType1) {
//                 continue;
//             }
//             for (std::size_t j = i + 1; j < initial_projection_face_count; ++j) {
//                 const auto& r_second_face = r_edge_data.ProjectionFaces[j];
//                 if (!r_second_face.IsType1 ||
//                     r_first_face.NormalAxis == r_second_face.NormalAxis ||
//                     !(r_first_face.ActiveSpan == r_second_face.ActiveSpan)) {
//                     continue;
//                 }

//                 SpanKey3D type2_span = r_first_face.ActiveSpan;
//                 add_to_span(type2_span, r_first_face.NormalAxis, r_first_face.NormalSign);
//                 add_to_span(type2_span, r_second_face.NormalAxis, r_second_face.NormalSign);
//                 if (!is_span_inside_domain(type2_span) || type1_external_spans.find(type2_span) != type1_external_spans.end()) {
//                     continue;
//                 }
//                 all_external_spans.insert(type2_span);

//                 ensure_skin_nodes_in_external_spans(std::set<SpanKey3D>{type2_span});

//                 skin_bins_per_knot_span =
//                     CreateSkinBinsPerKnotSpanMatrix3D(rSkinSubModelPart, rSurrogateSubModelPart);

//                 std::vector<Geometry<Node>::Pointer> candidates = r_first_face.CandidateNeighbourGeometries;
//                 candidates.insert(candidates.end(), r_second_face.CandidateNeighbourGeometries.begin(), r_second_face.CandidateNeighbourGeometries.end());

//                 auto p_k0_type2 = get_or_create_closest_skin_node_in_external_span(
//                     r_edge_data.pS0->Coordinates(),
//                     type2_span,
//                     "type2 edge=" + std::to_string(r_edge_data.pS0->Id()) + "-" +
//                     std::to_string(r_edge_data.pS1->Id()) +
//                     " from_spans=" + span_to_string(r_first_face.ExternalSpan) + "," +
//                     span_to_string(r_second_face.ExternalSpan) +
//                     " active_span=" + span_to_string(r_first_face.ActiveSpan));
                
//                 auto p_k1_type2 = get_or_create_closest_skin_node_in_external_span(
//                     r_edge_data.pS1->Coordinates(),
//                     type2_span,
//                     "type2 edge=" + std::to_string(r_edge_data.pS0->Id()) + "-" +
//                     std::to_string(r_edge_data.pS1->Id()) +
//                     " from_spans=" + span_to_string(r_first_face.ExternalSpan) + "," +
//                     span_to_string(r_second_face.ExternalSpan) +
//                     " active_span=" + span_to_string(r_first_face.ActiveSpan));

//                 add_projection_face_to_edge(
//                     r_edge_data.pS0,
//                     r_edge_data.pS1,
//                     p_k0_type2,
//                     p_k1_type2,
//                     type2_span,
//                     r_first_face.ActiveSpan,
//                     candidates,
//                     -1,
//                     0,
//                     false);
                
//                 add_projection_segment_to_node(r_edge_data.pS0, p_k0_type2, type2_span);
//                 add_projection_segment_to_node(r_edge_data.pS1, p_k1_type2, type2_span);
//                 ++number_of_type2_projection_faces;
//             }
//         }
//     }

//     std::set<std::tuple<IndexType, IndexType, SpanKey3D, SpanKey3D>> created_type2_pairs;
//     for (const auto& r_edge_entry : surrogate_edges) {
//         const auto& r_edge_data = r_edge_entry.second;
//         for (std::size_t i = 0; i < r_edge_data.ProjectionFaces.size(); ++i) {
//             for (std::size_t j = i + 1; j < r_edge_data.ProjectionFaces.size(); ++j) {
//                 const auto& r_face_a = r_edge_data.ProjectionFaces[i];
//                 const auto& r_face_b = r_edge_data.ProjectionFaces[j];
//                 if (r_face_a.IsType1 && r_face_b.IsType1) {
//                     continue;
//                 }

//                 const int manhattan_distance =
//                     std::abs(r_face_a.ExternalSpan.I - r_face_b.ExternalSpan.I) +
//                     std::abs(r_face_a.ExternalSpan.J - r_face_b.ExternalSpan.J) +
//                     std::abs(r_face_a.ExternalSpan.K - r_face_b.ExternalSpan.K);
//                 if (manhattan_distance != 1) {
//                     continue;
//                 }

//                 SpanKey3D first_span = r_face_a.ExternalSpan;
//                 SpanKey3D second_span = r_face_b.ExternalSpan;
//                 if (second_span < first_span) {
//                     std::swap(first_span, second_span);
//                 }
//                 const auto pair_key = std::make_tuple(r_edge_data.pS0->Id(), r_edge_data.pS1->Id(), first_span, second_span);
//                 if (!created_type2_pairs.insert(pair_key).second) {
//                     continue;
//                 }

//                 std::vector<Geometry<Node>::Pointer> candidates = r_face_a.CandidateNeighbourGeometries;
//                 candidates.insert(candidates.end(), r_face_b.CandidateNeighbourGeometries.begin(), r_face_b.CandidateNeighbourGeometries.end());
//                 const array_1d<double, 3> top_centroid =
//                     0.25 * (r_face_a.pK0->Coordinates() + r_face_a.pK1->Coordinates() +
//                             r_face_b.pK0->Coordinates() + r_face_b.pK1->Coordinates());
//                 auto p_chosen_neighbour = choose_neighbour_geometry(candidates, top_centroid);

//                 Node::Pointer p_k0_a = r_face_a.pK0;
//                 Node::Pointer p_k1_a = r_face_a.pK1;
//                 Node::Pointer p_k0_b = r_face_b.pK0;
//                 Node::Pointer p_k1_b = r_face_b.pK1;
//                 auto make_type2_volume_nodes = [&]() {
//                     return std::array<Node::Pointer, 8>{{
//                         r_edge_data.pS0,
//                         r_edge_data.pS1,
//                         r_edge_data.pS0,
//                         r_edge_data.pS1,
//                         p_k0_a,
//                         p_k1_a,
//                         p_k0_b,
//                         p_k1_b}};
//                 };

//                 std::vector<Node::Pointer> top_nodes{p_k0_a, p_k1_a, p_k1_b, p_k0_b};

//                 const std::array<Node::Pointer, 8> volume_nodes = make_type2_volume_nodes();

//                 auto are_spans_face_adjacent = [](const SpanKey3D& a, const SpanKey3D& b) {
//                     const int di = std::abs(a.I - b.I);
//                     const int dj = std::abs(a.J - b.J);
//                     const int dk = std::abs(a.K - b.K);
                
//                     return (di + dj + dk) == 1;
//                 };
//                 KRATOS_INFO("Type2PairDebug")
//                         << "Creating type-2 candidate"
//                         << " edge=" << r_edge_data.pS0->Id() << "-" << r_edge_data.pS1->Id() << "\n"
//                         << " span_a=" << span_to_string(r_face_a.ExternalSpan)
//                         << " is_type1_a=" << r_face_a.IsType1 << "\n"
//                         << " span_b=" << span_to_string(r_face_b.ExternalSpan)
//                         << " is_type1_b=" << r_face_b.IsType1 << "\n"
//                         << " adjacent=" << are_spans_face_adjacent(r_face_a.ExternalSpan, r_face_b.ExternalSpan) << "\n"
//                         << "S0=" << r_edge_data.pS0->Coordinates() << "\n"
//                         << "S1=" << r_edge_data.pS1->Coordinates() << "\n"
//                         << "A K0=" << r_face_a.pK0->Coordinates()
//                         << " K1=" << r_face_a.pK1->Coordinates() << "\n"
//                         << "B K0=" << r_face_b.pK0->Coordinates()
//                         << " K1=" << r_face_b.pK1->Coordinates() << "\n"
//                         << std::endl;
//                 const auto gap_volume_id = create_gap_volume(
//                     2,
//                     r_face_a.ExternalSpan,
//                     volume_nodes,
//                     {r_edge_data.pS0, r_edge_data.pS1, p_k0_a, p_k1_a, p_k0_b, p_k1_b},
//                     p_chosen_neighbour,
//                     "type-2 edge " + std::to_string(r_edge_data.pS0->Id()) + "-" + std::to_string(r_edge_data.pS1->Id()));
                
//                 if (!gap_volume_id) {
//                     continue;
//                 }

//                 if (compute_face_area(top_nodes) > area_tol) {
//                     auto p_top_surface = create_linear_surface(top_nodes);
//                     create_top_condition(2, p_top_surface, top_nodes, p_chosen_neighbour);
//                 }

//                 register_lateral_face({r_edge_data.pS0, r_edge_data.pS1, p_k1_a, p_k0_a}, {0, 0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//                 register_lateral_face({r_edge_data.pS0, r_edge_data.pS1, p_k1_b, p_k0_b}, {0, 0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//                 register_lateral_face({r_edge_data.pS0, p_k0_a, p_k0_b}, {0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//                 register_lateral_face({r_edge_data.pS1, p_k1_a, p_k1_b}, {0, 1, 1}, p_chosen_neighbour, *gap_volume_id);

//             }
//         }
//     }

//     return;

//     // for (const auto& r_node_faces_entry : face_indices_by_node) {
//     //     auto& r_node_data = surrogate_nodes[r_node_faces_entry.first];
//     //     const auto& r_face_indices = r_node_faces_entry.second;
//     //     for (std::size_t i = 0; i < r_face_indices.size(); ++i) {
//     //         for (std::size_t j = i + 1; j < r_face_indices.size(); ++j) {
//     //             for (std::size_t k = j + 1; k < r_face_indices.size(); ++k) {
//     //                 const auto& r_face_i = surrogate_faces[r_face_indices[i]];
//     //                 const auto& r_face_j = surrogate_faces[r_face_indices[j]];
//     //                 const auto& r_face_k = surrogate_faces[r_face_indices[k]];
//     //                 if (r_face_i.NormalAxis == r_face_j.NormalAxis ||
//     //                     r_face_i.NormalAxis == r_face_k.NormalAxis ||
//     //                     r_face_j.NormalAxis == r_face_k.NormalAxis ||
//     //                     !(r_face_i.ActiveSpan == r_face_j.ActiveSpan) ||
//     //                     !(r_face_i.ActiveSpan == r_face_k.ActiveSpan)) {
//     //                     continue;
//     //                 }

//     //                 SpanKey3D type3_span = r_face_i.ActiveSpan;
//     //                 add_to_span(type3_span, r_face_i.NormalAxis, r_face_i.NormalSign);
//     //                 add_to_span(type3_span, r_face_j.NormalAxis, r_face_j.NormalSign);
//     //                 add_to_span(type3_span, r_face_k.NormalAxis, r_face_k.NormalSign);
//     //                 if (!is_span_inside_domain(type3_span)) {
//     //                     continue;
//     //                 }

//     //                 const auto projection_key = std::make_pair(r_node_data.pNode->Id(), type3_span);
//     //                 if (!type3_projection_keys.insert(projection_key).second) {
//     //                     continue;
//     //                 }

//     //                 LocalProjectionData projection_data;
//     //                 projection_data.pSurrogateNode = r_node_data.pNode;
//     //                 // projection_data.pSkinNode = get_or_create_closest_skin_node_in_external_span(
//     //                 //     r_node_data.pNode->Coordinates(),
//     //                 //     type3_span);
                    
//     //                 projection_data.pSkinNode = get_or_create_closest_skin_node_in_external_span(
//     //                         r_node_data.pNode->Coordinates(),
//     //                         type3_span,
//     //                         "type3 node=" + std::to_string(r_node_data.pNode->Id()) +
//     //                         " active_span=" + span_to_string(r_face_i.ActiveSpan));
//     //                 projection_data.ExternalSpan = type3_span;
//     //                 projection_data.SignTriple = span_sign_around_node(type3_span, r_node_data.Grid);
//     //                 projection_data.Vector = projection_data.pSkinNode->Coordinates() - r_node_data.pNode->Coordinates();
//     //                 r_node_data.ProjectionSegments.push_back(projection_data);
//     //                 all_external_spans.insert(type3_span);
//     //                 ++number_of_type3_projection_segments;
//     //             }
//     //         }
//     //     }
//     // }

//     // std::set<std::tuple<IndexType, int, int>> created_type3_groups;
//     // auto find_projection_by_signs = [](const std::vector<LocalProjectionData>& rProjections,
//     //                                    const std::array<int, 3>& rSigns) -> const LocalProjectionData* {
//     //     for (const auto& r_projection : rProjections) {
//     //         if (r_projection.SignTriple == rSigns) {
//     //             return &r_projection;
//     //         }
//     //     }
//     //     return nullptr;
//     // };

//     // for (auto& r_node_entry : surrogate_nodes) {
//     //     auto& r_node_data = r_node_entry.second;
//     //     if (r_node_data.ProjectionSegments.size() < 4) {
//     //         continue;
//     //     }

//     //     for (int fixed_axis = 0; fixed_axis < 3; ++fixed_axis) {
//     //         for (const int fixed_sign : {-1, 1}) {
//     //             const auto group_key = std::make_tuple(r_node_data.pNode->Id(), fixed_axis, fixed_sign);
//     //             if (!created_type3_groups.insert(group_key).second) {
//     //                 continue;
//     //             }

//     //             std::array<std::array<int, 3>, 4> ordered_signs;
//     //             if (fixed_axis == 0) {
//     //                 ordered_signs = {{{fixed_sign, -1, -1}, {fixed_sign, 1, -1}, {fixed_sign, 1, 1}, {fixed_sign, -1, 1}}};
//     //             } else if (fixed_axis == 1) {
//     //                 ordered_signs = {{{-1, fixed_sign, -1}, {1, fixed_sign, -1}, {1, fixed_sign, 1}, {-1, fixed_sign, 1}}};
//     //             } else {
//     //                 ordered_signs = {{{-1, -1, fixed_sign}, {1, -1, fixed_sign}, {1, 1, fixed_sign}, {-1, 1, fixed_sign}}};
//     //             }

//     //             std::array<const LocalProjectionData*, 4> projections{{nullptr, nullptr, nullptr, nullptr}};
//     //             bool complete_group = true;
//     //             for (std::size_t i = 0; i < ordered_signs.size(); ++i) {
//     //                 projections[i] = find_projection_by_signs(r_node_data.ProjectionSegments, ordered_signs[i]);
//     //                 if (!projections[i]) {
//     //                     complete_group = false;
//     //                     break;
//     //                 }
//     //             }
//     //             if (!complete_group) {
//     //                 continue;
//     //             }

//     //             std::vector<Geometry<Node>::Pointer> candidates;
//     //             const std::array<std::vector<Node::Pointer>, 4> lateral_triangles{{
//     //                 {r_node_data.pNode, projections[0]->pSkinNode, projections[1]->pSkinNode},
//     //                 {r_node_data.pNode, projections[1]->pSkinNode, projections[2]->pSkinNode},
//     //                 {r_node_data.pNode, projections[2]->pSkinNode, projections[3]->pSkinNode},
//     //                 {r_node_data.pNode, projections[3]->pSkinNode, projections[0]->pSkinNode}}};
//     //             for (const auto& r_triangle_nodes : lateral_triangles) {
//     //                 const auto key = make_face_key({{r_triangle_nodes[0], 0}, {r_triangle_nodes[1], 1}, {r_triangle_nodes[2], 1}});
//     //                 const auto registry_it = lateral_face_registry.find(key);
//     //                 if (registry_it != lateral_face_registry.end()) {
//     //                     candidates.insert(candidates.end(),
//     //                                       registry_it->second.NeighbourGeometries.begin(),
//     //                                       registry_it->second.NeighbourGeometries.end());
//     //                 }
//     //             }
//     //             if (candidates.empty()) {
//     //                 candidates = r_node_data.IncidentFaceNeighbourGeometries;
//     //             }

//     //             const array_1d<double, 3> top_centroid =
//     //                 0.25 * (projections[0]->pSkinNode->Coordinates() +
//     //                         projections[1]->pSkinNode->Coordinates() +
//     //                         projections[2]->pSkinNode->Coordinates() +
//     //                         projections[3]->pSkinNode->Coordinates());
//     //             auto p_chosen_neighbour = choose_neighbour_geometry(candidates, top_centroid);

//     //             std::array<Node::Pointer, 4> top_skin_nodes{{
//     //                 projections[0]->pSkinNode,
//     //                 projections[1]->pSkinNode,
//     //                 projections[2]->pSkinNode,
//     //                 projections[3]->pSkinNode}};
//     //             auto make_type3_volume_nodes = [&]() {
//     //                 return std::array<Node::Pointer, 8>{{
//     //                     r_node_data.pNode,
//     //                     r_node_data.pNode,
//     //                     r_node_data.pNode,
//     //                     r_node_data.pNode,
//     //                     top_skin_nodes[0],
//     //                     top_skin_nodes[1],
//     //                     top_skin_nodes[3],
//     //                     top_skin_nodes[2]}};
//     //             };

//     //             std::vector<Node::Pointer> top_nodes{
//     //                 top_skin_nodes[0],
//     //                 top_skin_nodes[1],
//     //                 top_skin_nodes[2],
//     //                 top_skin_nodes[3]};
//     //             std::string type3_top_failure;
//     //             std::string type3_volume_failure;
//     //             bool type3_top_is_valid =
//     //                 is_ordered_face_non_degenerate(top_nodes, nullptr, type3_top_failure);
//     //             bool type3_volume_is_valid =
//     //                 is_volume_candidate_sane(make_type3_volume_nodes(), 3, projections[0]->ExternalSpan, type3_volume_failure);

//     //             bool collapsed_type3_top = false;
//     //             if (!type3_top_is_valid || !type3_volume_is_valid) {
//     //                 collapsed_type3_top = true;
//     //                 // const auto p_collapse_node = get_or_create_closest_skin_node_in_external_span(
//     //                 //     r_node_data.pNode->Coordinates(),
//     //                 //     projections[0]->ExternalSpan);

//     //                 const auto p_collapse_node = get_or_create_closest_skin_node_in_external_span(
//     //                     r_node_data.pNode->Coordinates(),
//     //                     projections[0]->ExternalSpan,
//     //                     "type3 collapsed node=" + std::to_string(r_node_data.pNode->Id()) +
//     //                     " fixed_axis=" + std::to_string(fixed_axis) +
//     //                     " fixed_sign=" + std::to_string(fixed_sign) +
//     //                     " span=" + span_to_string(projections[0]->ExternalSpan));
//     //                 top_skin_nodes = {{
//     //                     p_collapse_node,
//     //                     p_collapse_node,
//     //                     p_collapse_node,
//     //                     p_collapse_node}};
//     //                 top_nodes = {
//     //                     top_skin_nodes[0],
//     //                     top_skin_nodes[1],
//     //                     top_skin_nodes[2],
//     //                     top_skin_nodes[3]};
//     //                 type3_volume_failure.clear();
//     //                 type3_volume_is_valid =
//     //                     is_volume_candidate_sane(make_type3_volume_nodes(), 3, projections[0]->ExternalSpan, type3_volume_failure);

//     //                 KRATOS_ERROR_IF_NOT(type3_volume_is_valid)
//     //                     << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: failed to create a sane collapsed type-3 volume."
//     //                     << " node=" << r_node_data.pNode->Id()
//     //                     << " span=" << span_to_string(projections[0]->ExternalSpan)
//     //                     << " top_failure=" << type3_top_failure
//     //                     << " volume_failure=" << type3_volume_failure
//     //                     << "." << std::endl;

//     //                 ++number_of_collapsed_top_faces;
//     //                 ++gap_type_diagnostics[3].CollapsedTops;
//     //             }

//     //             const std::array<Node::Pointer, 8> volume_nodes = make_type3_volume_nodes();

//             //     const auto gap_volume_id = create_gap_volume(
//             //         3,
//             //         projections[0]->ExternalSpan,
//             //         volume_nodes,
//             //         {r_node_data.pNode,
//             //          top_skin_nodes[0],
//             //          top_skin_nodes[1],
//             //          top_skin_nodes[2],
//             //          top_skin_nodes[3]},
//             //         p_chosen_neighbour,
//             //         "type-3 node " + std::to_string(r_node_data.pNode->Id()));
                
//             //     if (!gap_volume_id) {
//             //         continue;
//             //     }

//             //     if (!collapsed_type3_top && compute_face_area(top_nodes) > area_tol) {
//             //         auto p_top_surface = create_linear_surface(top_nodes);
//             //         create_top_condition(3, p_top_surface, top_nodes, p_chosen_neighbour);
//             //     }

//             //     const std::array<std::vector<Node::Pointer>, 4> selected_lateral_triangles{{
//             //         {r_node_data.pNode, top_skin_nodes[0], top_skin_nodes[1]},
//             //         {r_node_data.pNode, top_skin_nodes[1], top_skin_nodes[2]},
//             //         {r_node_data.pNode, top_skin_nodes[2], top_skin_nodes[3]},
//             //         {r_node_data.pNode, top_skin_nodes[3], top_skin_nodes[0]}}};
//             //     register_lateral_face(selected_lateral_triangles[0], {0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//             //     register_lateral_face(selected_lateral_triangles[1], {0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//             //     register_lateral_face(selected_lateral_triangles[2], {0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//             //     register_lateral_face(selected_lateral_triangles[3], {0, 1, 1}, p_chosen_neighbour, *gap_volume_id);
//             // }
//     //     }
//     // }

//     // number_of_external_spans = all_external_spans.size();

//     // for (const auto& r_registry_entry : lateral_face_registry) {
//     //     const auto& r_lateral_face_data = r_registry_entry.second;
//     //     const auto neighbours = distinct_neighbours(r_lateral_face_data.NeighbourGeometries);
//     //     KRATOS_ERROR_IF(neighbours.empty())
//     //         << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: lateral face has no NEIGHBOUR_GEOMETRY." << std::endl;
//     //     if (neighbours.size() == 1) {
//     //         continue;
//     //     }
//     //     if (neighbours.size() > 2) {
//     //         ++number_of_non_manifold_lateral_faces;
//     //         KRATOS_ERROR
//     //             << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries3D: non-manifold lateral face with "
//     //             << neighbours.size() << " distinct NEIGHBOUR_GEOMETRY entries." << std::endl;
//     //     }
//     //     if (neighbours[0]->Id() == neighbours[1]->Id() || r_lateral_face_data.IsDegenerate) {
//     //         continue;
//     //     }

//     //     IntegrationPointsArrayType interface_integration_points;
//     //     GeometriesArrayType interface_quadrature_point_list;
//     //     const std::size_t interface_integration_order = static_cast<std::size_t>(2 * mGapInterpolationOrder + 1);
//     //     r_lateral_face_data.pGeometry->CreateIntegrationPoints(
//     //         interface_integration_points,
//     //         interface_integration_order,
//     //         interface_integration_order);

//     //     double interface_weight_sum = 0.0;
//     //     for (auto& r_integration_point : interface_integration_points) {
//     //         const double determinant_jacobian =
//     //             r_lateral_face_data.pGeometry->DeterminantOfJacobian(r_integration_point);
//     //         r_integration_point.SetWeight(r_integration_point.Weight() * determinant_jacobian);
//     //         interface_weight_sum += r_integration_point.Weight();
//     //     }
//     //     if (interface_weight_sum <= area_tol) {
//     //         continue;
//     //     }

//     //     IntegrationInfo interface_integration_info = r_lateral_face_data.pGeometry->GetDefaultIntegrationInfo();
//     //     interface_integration_info.SetNumberOfIntegrationPointsPerSpan(0, 2 * interface_integration_order + 1);
//     //     interface_integration_info.SetNumberOfIntegrationPointsPerSpan(1, 2 * interface_integration_order + 1);
//     //     r_lateral_face_data.pGeometry->CreateQuadraturePointGeometries(
//     //         interface_quadrature_point_list,
//     //         number_of_shape_functions_derivatives,
//     //         interface_integration_points,
//     //         interface_integration_info);
//     //     if (interface_quadrature_point_list.size() == 0) {
//     //         continue;
//     //     }

//     //     std::size_t id_interface_condition = 1;
//     //     if (mpGapInterfaceSubModelPart->GetRootModelPart().Conditions().size() > 0) {
//     //         id_interface_condition = mpGapInterfaceSubModelPart->GetRootModelPart().Conditions().back().Id() + 1;
//     //     }

//     //     this->CreateConditions(
//     //         interface_quadrature_point_list.ptr_begin(),
//     //         interface_quadrature_point_list.ptr_end(),
//     //         *mpGapInterfaceSubModelPart,
//     //         mGapInterfaceConditionName,
//     //         id_interface_condition,
//     //         PropertiesPointerType(),
//     //         knot_span_sizes,
//     //         neighbours,
//     //         std::max(std::sqrt(interface_weight_sum), distance_tol));
//     //     number_of_interface_conditions += interface_quadrature_point_list.size();
//     // }

//     // if (min_volume_weight == std::numeric_limits<double>::max()) {
//     //     min_volume_weight = 0.0;
//     // }
//     // if (min_determinant == std::numeric_limits<double>::max()) {
//     //     min_determinant = 0.0;
//     // }

//     KRATOS_INFO_IF("CreateSbmExtendedGeometries3D", mEchoLevel > 0)
//         << "3D Gap-SBM local projection summary:"
//         << " external_spans=" << number_of_external_spans
//         << " type1_spans=" << type1_external_spans.size()
//         << " type1_volumes=" << gap_type_diagnostics[1].Volumes
//         << " type1_top_conditions=" << gap_type_diagnostics[1].TopConditions
//         << " type1_collapsed_tops=" << gap_type_diagnostics[1].CollapsedTops
//         << " type1_min_weight=" << gap_type_diagnostics[1].PrintableMinWeight()
//         << " type1_max_weight=" << gap_type_diagnostics[1].MaxWeight
//         << " type1_min_det=" << gap_type_diagnostics[1].PrintableMinDeterminant()
//         << " type2_projection_faces=" << number_of_type2_projection_faces
//         << " type2_volumes=" << gap_type_diagnostics[2].Volumes
//         << " type2_top_conditions=" << gap_type_diagnostics[2].TopConditions
//         << " type2_collapsed_tops=" << gap_type_diagnostics[2].CollapsedTops
//         << " type2_min_weight=" << gap_type_diagnostics[2].PrintableMinWeight()
//         << " type2_max_weight=" << gap_type_diagnostics[2].MaxWeight
//         << " type2_min_det=" << gap_type_diagnostics[2].PrintableMinDeterminant()
//         << " type3_projection_segments=" << number_of_type3_projection_segments
//         << " type3_pyramids=" << gap_type_diagnostics[3].Volumes
//         << " type3_top_conditions=" << gap_type_diagnostics[3].TopConditions
//         << " type3_collapsed_tops=" << gap_type_diagnostics[3].CollapsedTops
//         << " type3_min_weight=" << gap_type_diagnostics[3].PrintableMinWeight()
//         << " type3_max_weight=" << gap_type_diagnostics[3].MaxWeight
//         << " type3_min_det=" << gap_type_diagnostics[3].PrintableMinDeterminant()
//         << " top_conditions=" << number_of_top_conditions
//         << " collapsed_top_faces=" << number_of_collapsed_top_faces
//         << " lateral_faces=" << lateral_face_registry.size()
//         << " interface_conditions=" << number_of_interface_conditions
//         << " non_manifold_lateral_faces=" << number_of_non_manifold_lateral_faces
//         << " min_volume_weight=" << min_volume_weight
//         << " max_volume_weight=" << max_volume_weight
//         << " min_det=" << min_determinant
//         << std::endl;

//     return;
// }
    

// } // namespace Kratos
