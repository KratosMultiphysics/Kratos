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
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "custom_processes/snake_gap_sbm_process.h"

namespace Kratos
{

namespace
{
    // Store the minimum and maximum coordinates of the two parametric directions.
    struct ParameterSpaceBounds
    {
        double MinU = 0.0;
        double MaxU = 0.0;
        double MinV = 0.0;
        double MaxV = 0.0;
    };

    // Read the two-dimensional parametric bounds from the model part.
    ParameterSpaceBounds ReadParameterSpaceBounds(
        const ModelPart& rModelPart,
        const char* pCaller)
    {
        ParameterSpaceBounds bounds;

        const auto& r_parameter_space_corners = rModelPart.GetValue(PARAMETER_SPACE_CORNERS);
        KRATOS_ERROR_IF(r_parameter_space_corners.size() < 2)
            << "[" << pCaller << "] PARAMETER_SPACE_CORNERS must contain at least two vectors.\n";
        KRATOS_ERROR_IF(r_parameter_space_corners[0].size() < 2 || r_parameter_space_corners[1].size() < 2)
            << "[" << pCaller << "] PARAMETER_SPACE_CORNERS entries must contain [min, max] values."
            << "\n";

        bounds.MinU = r_parameter_space_corners[0][0];
        bounds.MaxU = r_parameter_space_corners[0][1];
        bounds.MinV = r_parameter_space_corners[1][0];
        bounds.MaxV = r_parameter_space_corners[1][1];
        return bounds;
    }

    // Validate uniform knot spacing and recover the integer number of spans.
    std::size_t ComputeSpanCount(
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

    // Maintain the parallel layer and condition-name lists attached to each skin node.
    void AddLayerConditionMetadata(
        Node& rNode,
        const std::string& rLayerName,
        const std::string& rConditionName)
    {
        auto connected_layers = rNode.GetValue(CONNECTED_LAYERS);
        auto connected_conditions = rNode.GetValue(CONNECTED_CONDITIONS);
        if (connected_conditions.size() < connected_layers.size()) {
            connected_conditions.resize(connected_layers.size());
        }

        const auto layer_it = std::find(connected_layers.begin(), connected_layers.end(), rLayerName);
        if (layer_it == connected_layers.end()) {
            connected_layers.push_back(rLayerName);
            connected_conditions.push_back(rConditionName);
        } else {
            const std::size_t layer_index = static_cast<std::size_t>(std::distance(connected_layers.begin(), layer_it));
            if (connected_conditions[layer_index].empty() && !rConditionName.empty()) {
                connected_conditions[layer_index] = rConditionName;
            }
        }

        rNode.SetValue(CONNECTED_LAYERS, connected_layers);
        rNode.SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }

    // Return the condition associated with a layer, if available.
    std::string GetConditionNameForLayer(
        const Node& rNode,
        const std::string& rLayerName)
    {
        const auto& connected_layers = rNode.GetValue(CONNECTED_LAYERS);
        const auto& connected_conditions = rNode.GetValue(CONNECTED_CONDITIONS);
        const auto layer_it = std::find(connected_layers.begin(), connected_layers.end(), rLayerName);
        if (layer_it == connected_layers.end()) {
            return "";
        }

        const std::size_t layer_index = static_cast<std::size_t>(std::distance(connected_layers.begin(), layer_it));
        if (layer_index >= connected_conditions.size()) {
            return "";
        }

        return connected_conditions[layer_index];
    }

    // Store the common layer and matching condition selected for a node pair.
    struct CommonLayerInfo
    {
        std::string LayerName;
        std::string ConditionName;
        bool Found = false;
    };

    // Find shared layer metadata, retaining the original last-match selection.
    CommonLayerInfo FindCommonLayerInfo(
        const Node& rFirstNode,
        const Node& rSecondNode)
    {
        const auto& connected_layers_1 = rFirstNode.GetValue(CONNECTED_LAYERS);
        const auto& connected_layers_2 = rSecondNode.GetValue(CONNECTED_LAYERS);
        std::unordered_set<std::string> layer_lookup(connected_layers_2.begin(), connected_layers_2.end());

        CommonLayerInfo layer_info;

        for (const auto& r_layer_name : connected_layers_1) {
            if (layer_lookup.find(r_layer_name) == layer_lookup.end()) {
                continue;
            }

            layer_info.LayerName = r_layer_name;
            layer_info.ConditionName = GetConditionNameForLayer(rFirstNode, r_layer_name);
            if (layer_info.ConditionName.empty()) {
                layer_info.ConditionName = GetConditionNameForLayer(rSecondNode, r_layer_name);
            }
            layer_info.Found = true;
        }
        return layer_info;
    }

    // Find an unused regular geometry ID, ignoring Kratos special IDs.
    ModelPart::IndexType NextRegularGeometryId(const ModelPart& rModelPart)
    {
        ModelPart::IndexType next_id = 0;
        for (const auto& r_geometry : rModelPart.Geometries()) {
            const auto id = r_geometry.Id();
            constexpr auto id_bits = sizeof(ModelPart::IndexType) * 8;
            const bool id_generated_from_string = id & (ModelPart::IndexType(1) << (id_bits - 1));
            const bool id_self_assigned = id & (ModelPart::IndexType(1) << (id_bits - 2));
            if (id_generated_from_string || id_self_assigned) {
                continue;
            }
            next_id = std::max(next_id, static_cast<ModelPart::IndexType>(id));
        }
        return next_id + 1;
    }
} // unnamed namespace

//-------------------------------

// Construct the gap and skin quadrature geometries for the inner and outer loops.
void SnakeGapSbmProcess::CreateSbmExtendedGeometries2D()
{
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 1)
        << "Echo level: " << mEchoLevel
        << " | Gap type: " << mGapSbmType
        << " | Approx order: " << mGapApproximationOrder
        << " | Internal divisions: " << mInternalDivisions
        << " | Rel tol subdivisions: " << mGapRelativeToleranceForSubdivisions
        << " | Interp levels: " << mNumberOfInterpolationLevels << std::endl;
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0 || mpSkinModelPartInnerInitial->NumberOfGeometries()>0) {
        const auto& r_surrogate_sub_model_part_inner = mpIgaModelPart->GetSubModelPart("surrogate_inner");
        auto& r_skin_sub_model_part_inner = mpSkinModelPart->GetSubModelPart("inner");

        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the inner skin." << std::endl;
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 2)
            << "Inner skin nodes: " << r_skin_sub_model_part_inner.NumberOfNodes()
            << " | conditions: " << r_skin_sub_model_part_inner.NumberOfConditions()
            << " | surrogate conditions: " << r_surrogate_sub_model_part_inner.NumberOfConditions() << std::endl;
        CreateSbmExtendedGeometries2D<true>(r_skin_sub_model_part_inner, r_surrogate_sub_model_part_inner);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the inner skin." << std::endl;
    }
    if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) {
        const auto& r_surrogate_sub_model_part_outer = mpIgaModelPart->GetSubModelPart("surrogate_outer");
        auto& r_skin_sub_model_part_outer = mpSkinModelPart->GetSubModelPart("outer");
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the outer skin." << std::endl;
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 2)
            << "Outer skin nodes: " << r_skin_sub_model_part_outer.NumberOfNodes()
            << " | conditions: " << r_skin_sub_model_part_outer.NumberOfConditions()
            << " | surrogate conditions: " << r_surrogate_sub_model_part_outer.NumberOfConditions() << std::endl;
        CreateSbmExtendedGeometries2D<false>(r_skin_sub_model_part_outer, r_surrogate_sub_model_part_outer);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the outer skin." << std::endl;
    }
}

template <bool TIsInnerLoop>
void SnakeGapSbmProcess::CreateSbmExtendedGeometries2D(
    ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart)
{
    // Get the mesh sizes from the surrogate model part
    const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 1)
        << "Knot span sizes: " << knot_span_sizes << std::endl;

    auto p_surface = mpIgaModelPart->pGetGeometry(1);
    KRATOS_ERROR_IF_NOT(p_surface)
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries: background geometry '"
        << mBackgroundSurfaceGeometryName << "' was not found." << std::endl;

    auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_surface);
    if (!p_nurbs_surface) {
        p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
            p_surface->pGetGeometryPart(
                Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    }
    KRATOS_ERROR_IF_NOT(p_nurbs_surface)
        << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries: background geometry '"
        << mBackgroundSurfaceGeometryName << "'"
        << " does not expose a NurbsSurfaceType background geometry." << std::endl;

    // Build knot-span lookups and set projections from surrogate to skin
    const auto& skin_bins_per_knot_span = CreateSkinBinsPerKnotSpanMatrix2D(rSkinSubModelPart, rSurrogateSubModelPart);

    SetSurrogateToSkinProjections2D<TIsInnerLoop>(rSurrogateSubModelPart, rSkinSubModelPart,
                                                skin_bins_per_knot_span);

    std::size_t brep_degree = p_nurbs_surface->PolynomialDegree(0);
    std::size_t number_of_shape_functions_derivatives = 2 * brep_degree + 1;

    if (mGapApproximationOrder == 0)
        mGapApproximationOrder = brep_degree;

    IndexType first_condition_id;
    IndexType last_condition_id;
    IndexType starting_brep_id;
    std::size_t size_surrogate_loop;

    if constexpr (TIsInnerLoop) {
        IndexType element_id = rSurrogateSubModelPart.ElementsBegin()->Id();
        first_condition_id = rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[0].Id();
        last_condition_id = rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[1].Id();
        size_surrogate_loop = last_condition_id - first_condition_id + 1;
        if (mpIgaModelPart->GetSubModelPart("surrogate_outer").NumberOfConditions() > 0)
        {
            starting_brep_id = 2 + mpIgaModelPart->GetSubModelPart("surrogate_outer").NumberOfConditions(); //1 surface + outer surrogate loop
        }
        else
            starting_brep_id = 6; //1 surface + 4 external boundaries
    }
    else {
        size_surrogate_loop = rSurrogateSubModelPart.NumberOfConditions();
        first_condition_id = rSurrogateSubModelPart.ConditionsBegin()->Id();
        last_condition_id = first_condition_id + size_surrogate_loop - 1;
        starting_brep_id = 2; //1 surface
    }

    struct QuantizedCoord
    {
        long long X;
        long long Y;
        long long Z;

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
            seed ^= std::hash<long long>{}(rValue.X) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            seed ^= std::hash<long long>{}(rValue.Y) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            seed ^= std::hash<long long>{}(rValue.Z) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            return seed;
        }
    };

    // Quantization is only used to find nearby hash buckets. A geometric
    // distance check below decides whether a node is an actual match.
    const double coordinate_tolerance = 1.0e-12;
    auto quantize_coord = [coordinate_tolerance](const CoordinatesArrayType& rCoords) {
        return QuantizedCoord{
            static_cast<long long>(std::llround(rCoords[0] / coordinate_tolerance)),
            static_cast<long long>(std::llround(rCoords[1] / coordinate_tolerance)),
            static_cast<long long>(std::llround(rCoords[2] / coordinate_tolerance))};
    };

    std::unordered_map<QuantizedCoord, std::vector<Node::Pointer>, QuantizedCoordHasher> surrogate_node_lookup;
    surrogate_node_lookup.reserve(rSurrogateSubModelPart.NumberOfNodes() * 2);
    for (auto& r_surrogate_node : rSurrogateSubModelPart.Nodes()) {
        surrogate_node_lookup[quantize_coord(r_surrogate_node.Coordinates())].push_back(&r_surrogate_node);
    }

    const double coordinate_tolerance_sq = coordinate_tolerance * coordinate_tolerance;
    auto find_surrogate_node = [&](const CoordinatesArrayType& rCoordinates) {
        const auto center_key = quantize_coord(rCoordinates);
        Node::Pointer p_best_node;
        double best_distance_sq = coordinate_tolerance_sq;
        for (long long dx = -1; dx <= 1; ++dx) {
            for (long long dy = -1; dy <= 1; ++dy) {
                for (long long dz = -1; dz <= 1; ++dz) {
                    const QuantizedCoord key{center_key.X + dx, center_key.Y + dy, center_key.Z + dz};
                    const auto bucket_it = surrogate_node_lookup.find(key);
                    if (bucket_it == surrogate_node_lookup.end()) {
                        continue;
                    }
                    for (const auto& p_candidate : bucket_it->second) {
                        const auto difference = p_candidate->Coordinates() - rCoordinates;
                        const double distance_sq = inner_prod(difference, difference);
                        if (distance_sq <= best_distance_sq) {
                            best_distance_sq = distance_sq;
                            p_best_node = p_candidate;
                        }
                    }
                }
            }
        }
        return p_best_node;
    };

    for (std::size_t j = 0; j < size_surrogate_loop; ++j) {
        auto p_brep_geometry = mpIgaModelPart->pGetGeometry(starting_brep_id + j);
        KRATOS_ERROR_IF_NOT(p_brep_geometry)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries: geometry with id "
            << (starting_brep_id + j) << " was not found." << std::endl;

        auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

        KRATOS_ERROR_IF_NOT(p_brep_curve_on_surface_surrogate1_surrogate2)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries: geometry with id "
            << p_brep_geometry->Id() << " is not a BrepCurveOnSurfaceType." << std::endl;

        NurbsInterval brep_domain_interval = p_brep_curve_on_surface_surrogate1_surrogate2->DomainInterval();
        CoordinatesArrayType surrogate_vertex_1 = ZeroVector(3);
        CoordinatesArrayType surrogate_vertex_1_local_coords = ZeroVector(3);
        CoordinatesArrayType surrogate_vertex_2 = ZeroVector(3);
        CoordinatesArrayType surrogate_vertex_2_local_coords = ZeroVector(3);
        surrogate_vertex_1_local_coords[0] = brep_domain_interval.GetT0();
        surrogate_vertex_2_local_coords[0] = brep_domain_interval.GetT1();

        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

        const array_1d<double, 3> vector_cond = surrogate_vertex_2 - surrogate_vertex_1;
        array_1d<double, 3> normal_cond;
        normal_cond[0] = vector_cond[1];
        normal_cond[1] = -vector_cond[0];
        normal_cond[2] =  0.0;

        // Use the BREP midpoint as the surrogate reference for the gap DOFs.
        CoordinatesArrayType surrogate_middle_point_local_coords = ZeroVector(3);
        surrogate_middle_point_local_coords[0] = 0.5 * (surrogate_vertex_1_local_coords[0] + surrogate_vertex_2_local_coords[0]);

        IntegrationPoint<1> integration_point(surrogate_middle_point_local_coords[0]);
        IntegrationPointsArrayType surrogate_integration_points_list;
        surrogate_integration_points_list.push_back(integration_point);

        IntegrationInfo integration_info = p_brep_curve_on_surface_surrogate1_surrogate2->GetDefaultIntegrationInfo();

        IntegrationParameters integration_parameters(
            number_of_shape_functions_derivatives,
            integration_info,
            knot_span_sizes);
        integration_parameters.pSkinBinsPerSpan = &skin_bins_per_knot_span;

        GeometriesArrayType quadrature_point_list;
        p_brep_curve_on_surface_surrogate1_surrogate2->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives,
                                                            surrogate_integration_points_list, integration_info);

        KRATOS_ERROR_IF(quadrature_point_list.size() == 0)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries: no quadrature geometries were created for BREP geometry #"
            << p_brep_geometry->Id() << std::endl;

        GeometryType::Pointer surrogate_brep_middle_geometry = quadrature_point_list(0);

        // Store the surrogate middle geometry for the lateral Breps
        Node::Pointer p_surrogate_node_1;
        Node::Pointer p_surrogate_node_2;

        p_surrogate_node_1 = find_surrogate_node(surrogate_vertex_1);
        p_surrogate_node_2 = find_surrogate_node(surrogate_vertex_2);

        if (p_surrogate_node_1) {
            auto& r_neighbour_geometries = p_surrogate_node_1->GetValue(NEIGHBOUR_GEOMETRIES);
            r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);
        }
        if (p_surrogate_node_2) {
            auto& r_neighbour_geometries = p_surrogate_node_2->GetValue(NEIGHBOUR_GEOMETRIES);
            r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);
        }

        // Fall back to a full scan if the spatial hash did not find a match.
        if (!p_surrogate_node_1 || !p_surrogate_node_2) {
            for (auto& r_surrogate_node : rSurrogateSubModelPart.Nodes()) {
                const double distance_sq_1 = inner_prod(r_surrogate_node.Coordinates() - surrogate_vertex_1,
                                                        r_surrogate_node.Coordinates() - surrogate_vertex_1);
                const double distance_sq_2 = inner_prod(r_surrogate_node.Coordinates() - surrogate_vertex_2,
                                                        r_surrogate_node.Coordinates() - surrogate_vertex_2);

                if (!p_surrogate_node_1 && distance_sq_1 <= coordinate_tolerance_sq) {
                    p_surrogate_node_1 = &r_surrogate_node;
                    auto& r_neighbour_geometries =
                        p_surrogate_node_1->GetValue(NEIGHBOUR_GEOMETRIES);
                    r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);
                }
                if (!p_surrogate_node_2 && distance_sq_2 <= coordinate_tolerance_sq)
                {
                    p_surrogate_node_2 = &r_surrogate_node;
                    auto& r_neighbour_geometries =
                        p_surrogate_node_2->GetValue(NEIGHBOUR_GEOMETRIES);
                    r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);
                }
                if (p_surrogate_node_1 && p_surrogate_node_2) {
                    break;
                }
            }
        }

        KRATOS_ERROR_IF(!p_surrogate_node_1 || !p_surrogate_node_2)
            << "::[SnakeGapSbmProcess]::CreateSbmExtendedGeometries: failed to match surrogate end nodes "
            << "for BREP geometry #" << p_brep_curve_on_surface_surrogate1_surrogate2->Id()
            << ". Vertex 1: " << surrogate_vertex_1
            << " | Vertex 2: " << surrogate_vertex_2
            << " | coordinate matches: "
            << (find_surrogate_node(surrogate_vertex_1) ? "yes" : "no") << ", "
            << (find_surrogate_node(surrogate_vertex_2) ? "yes" : "no")
            << std::endl;

        //------------------------------------------------------------------
        // check layer projection compatibility and build projection chain if necessary

        const IndexType projection_id_1_base = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
        const IndexType projection_id_2_base = p_surrogate_node_2->GetValue(PROJECTION_NODE_ID);

        KRATOS_ERROR_IF(projection_id_1_base == 0 || projection_id_2_base == 0)
            << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: PROJECTION_NODE_ID not set for surrogate nodes "
            << p_surrogate_node_1->Id() << " and/or " << p_surrogate_node_2->Id()
            << p_surrogate_node_1->Coordinates() << " and/or " << p_surrogate_node_2->Coordinates()
            << " | ids=" << projection_id_1_base << "," << projection_id_2_base << std::endl;

        auto have_common_layer = [](const NodeType& rFirstProjection, const NodeType& rSecondProjection) {
            const auto& layers_1 = rFirstProjection.GetValue(CONNECTED_LAYERS);
            const auto& layers_2 = rSecondProjection.GetValue(CONNECTED_LAYERS);
            std::unordered_set<std::string> layer_lookup(layers_1.begin(), layers_1.end());
            for (const auto& layer : layers_2) {
                if (layer_lookup.find(layer) != layer_lookup.end()) {
                    return true;
                }
            }
            return false;
        };

        // Create temporary surrogate nodes with the metadata of their skin projections.
        auto create_projected_clone = [&](const CoordinatesArrayType& rCoordinates, const IndexType projection_id) {
            Node::Pointer p_clone = Kratos::make_intrusive<Node>(0, rCoordinates);
            KRATOS_ERROR_IF(projection_id == 0 || !rSkinSubModelPart.HasNode(projection_id))
                << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: invalid PROJECTION_NODE_ID="
                << projection_id << " | coords=" << rCoordinates << std::endl;
            const auto& r_projection_node = rSkinSubModelPart.GetNode(projection_id);
            p_clone->SetValue(PROJECTION_NODE_ID, projection_id);
            p_clone->SetValue(CONNECTED_LAYERS, r_projection_node.GetValue(CONNECTED_LAYERS));
            p_clone->SetValue(CONNECTED_CONDITIONS, r_projection_node.GetValue(CONNECTED_CONDITIONS));
            return p_clone;
        };

        auto process_pair = [&](const Node::Pointer& rp_surrogate_node_1, const Node::Pointer& rp_surrogate_node_2) {
            const auto common_layer_info = FindCommonLayerInfo(*rp_surrogate_node_1, *rp_surrogate_node_2);
            const std::string common_layer_name = common_layer_info.LayerName;
            const std::string common_condition_name = common_layer_info.ConditionName;
            KRATOS_ERROR_IF(common_layer_info.Found && common_condition_name.empty())
                << ":::[SnakeGapSbmProcess]::: Missing condition metadata for common layer '"
                << common_layer_name << "' between surrogate nodes "
                << rp_surrogate_node_1->Id() << " and " << rp_surrogate_node_2->Id() << "." << std::endl;

            // Connect the original projections across layers before the normal checks.
            std::vector<Node::Pointer> chained_segment_nodes{rp_surrogate_node_1};
            const auto& p_first_chain_node = rp_surrogate_node_1;
            const auto& p_second_chain_node = rp_surrogate_node_2;
            const IndexType first_projection_id = p_first_chain_node->GetValue(PROJECTION_NODE_ID);
            const IndexType second_projection_id = p_second_chain_node->GetValue(PROJECTION_NODE_ID);

            KRATOS_ERROR_IF(first_projection_id == 0 || second_projection_id == 0)
                << "::[SnakeGapSbmProcess]:: Missing PROJECTION_NODE_ID while building corrected projection chain."
                << std::endl;

            const auto& r_first_projection = rSkinSubModelPart.GetNode(first_projection_id);
            const auto& r_second_projection = rSkinSubModelPart.GetNode(second_projection_id);
            if (!have_common_layer(r_first_projection, r_second_projection)) {

                const auto projection_chain = AssessProjectionsFeasibility(
                    rSkinSubModelPart,
                    p_first_chain_node,
                    p_second_chain_node);

                KRATOS_ERROR_IF(projection_chain.size() < 2)
                    << "::[SnakeGapSbmProcess]:: Empty projection chain detected for corrected surrogate nodes "
                    << p_first_chain_node->Id() << " and " << p_second_chain_node->Id() << std::endl;

                const double division_count = static_cast<double>(projection_chain.size() - 1);
                for (std::size_t chain_index = 1; chain_index + 1 < projection_chain.size(); ++chain_index) {
                    const double t = static_cast<double>(chain_index) / division_count;
                    array_1d<double, 3> interpolated_coordinates =
                        (1.0 - t) * p_first_chain_node->Coordinates()
                        + t * p_second_chain_node->Coordinates();
                    chained_segment_nodes.push_back(
                        create_projected_clone(interpolated_coordinates, projection_chain[chain_index]));
                }
            }
            chained_segment_nodes.push_back(p_second_chain_node);

            std::vector<Node::Pointer> surrogate_segment_nodes;
            for (std::size_t i = 0; i + 1 < chained_segment_nodes.size(); ++i) {
                Node::Pointer p_first_for_fallback = chained_segment_nodes[i];
                Node::Pointer p_second_for_fallback = chained_segment_nodes[i + 1];
                const auto pair_layer_info = FindCommonLayerInfo(*p_first_for_fallback, *p_second_for_fallback);
                const bool is_coupling_layer = pair_layer_info.LayerName == "COUPLING_SIDE" ||
                    pair_layer_info.LayerName == "COUPLING_SIDE_OUTER" ||
                    pair_layer_info.LayerName == "COUPLING_SIDE_INNER";

                const IndexType proj_id_first_base  = p_first_for_fallback->GetValue(PROJECTION_NODE_ID);
                const IndexType proj_id_second_base = p_second_for_fallback->GetValue(PROJECTION_NODE_ID);

                array_1d<double, 3> v_proj_1_base = rSkinSubModelPart.GetNode(proj_id_first_base).Coordinates() - p_first_for_fallback->Coordinates();
                array_1d<double, 3> v_proj_2_base = rSkinSubModelPart.GetNode(proj_id_second_base).Coordinates() - p_second_for_fallback->Coordinates();

                // Check whether the current skin projections lie on the expected side.
                double dot1_base = v_proj_1_base[0]*normal_cond[0] + v_proj_1_base[1]*normal_cond[1];
                double dot2_base = v_proj_2_base[0]*normal_cond[0] + v_proj_2_base[1]*normal_cond[1];

                if (i == 0) {
                    surrogate_segment_nodes.push_back(p_first_for_fallback);
                }

                // Correct the first projection if it lies on the wrong side.
                if (!is_coupling_layer && dot1_base < -1.0e-12) {
                    Vector direction(3);
                    direction[0] = 0.5 * (vector_cond[0] - normal_cond[0]);
                    direction[1] = 0.5 * (vector_cond[1] - normal_cond[1]);
                    direction[2] = 0.5 * (vector_cond[2] - normal_cond[2]);

                    if constexpr (!TIsInnerLoop) {
                        direction *= -1.0;
                    }

                    const double direction_norm = norm_2(direction);
                    KRATOS_ERROR_IF(direction_norm < 1e-13) << "Error: zero direction vector.\n";
                    direction = direction / direction_norm;
                    const auto split_pair = FindClosestPairInLayerWithNormalDirection<TIsInnerLoop>(
                        p_first_for_fallback->Coordinates(),
                        "",
                        rSkinSubModelPart,
                        integration_parameters.rKnotSpanSizes,
                        skin_bins_per_knot_span,
                        direction, nullptr, nullptr, TIsInnerLoop);

                    // Only one should satisfy the inner product condition
                    array_1d<double, 3> v_pair_1 = rSkinSubModelPart.GetNode(split_pair.first).Coordinates() - p_first_for_fallback->Coordinates();
                    array_1d<double, 3> v_pair_2 = rSkinSubModelPart.GetNode(split_pair.second).Coordinates() - p_second_for_fallback->Coordinates();

                    const double dot_pair1 = v_pair_1[0]*normal_cond[0] + v_pair_1[1]*normal_cond[1];
                    const double dot_pair2 = v_pair_2[0]*normal_cond[0] + v_pair_2[1]*normal_cond[1];

                    Node::Pointer p_split_node;
                    if (dot_pair1 > 0) {
                        p_split_node = rSkinSubModelPart.pGetNode(split_pair.first);
                    } else if (dot_pair2 > 0) {
                        p_split_node = rSkinSubModelPart.pGetNode(split_pair.second);
                    } else {
                        p_split_node = rSkinSubModelPart.pGetNode(split_pair.first);
                        KRATOS_WARNING("Neither candidate split nodes satisfy the inner product condition.\n");
                    }

                    if (pair_layer_info.Found) {
                        AddLayerConditionMetadata(*p_split_node, pair_layer_info.LayerName, pair_layer_info.ConditionName);
                    }

                    surrogate_segment_nodes.push_back(
                        create_projected_clone(p_first_for_fallback->Coordinates(), p_split_node->Id()));

                }
                if (!is_coupling_layer && dot2_base < -1.0e-12) {
                    Vector direction(3);
                    direction[0] = 0.5 * (vector_cond[0] + normal_cond[0]);
                    direction[1] = 0.5 * (vector_cond[1] + normal_cond[1]);
                    direction[2] = 0.5 * (vector_cond[2] + normal_cond[2]);

                    if constexpr (!TIsInnerLoop) {
                        direction *= -1;
                    }

                    const double direction_norm = norm_2(direction);
                    KRATOS_ERROR_IF(direction_norm < 1e-13) << "Error: zero direction vector.\n";
                    direction = -direction/direction_norm;

                    const auto split_pair = FindClosestPairInLayerWithNormalDirection<TIsInnerLoop>(
                        p_second_for_fallback->Coordinates(),
                        "",
                        rSkinSubModelPart,
                        integration_parameters.rKnotSpanSizes,
                        skin_bins_per_knot_span,
                        direction, nullptr, nullptr, TIsInnerLoop);

                    // Only one satisfy inner product condition
                    array_1d<double, 3> v_pair_1 = rSkinSubModelPart.GetNode(split_pair.first).Coordinates() - p_first_for_fallback->Coordinates();
                    array_1d<double, 3> v_pair_2 = rSkinSubModelPart.GetNode(split_pair.second).Coordinates() - p_second_for_fallback->Coordinates();

                    const double dot_pair1 = v_pair_1[0]*normal_cond[0] + v_pair_1[1]*normal_cond[1];
                    const double dot_pair2 = v_pair_2[0]*normal_cond[0] + v_pair_2[1]*normal_cond[1];

                    Node::Pointer p_split_node;
                    if (dot_pair1 > 0) {
                        // Use first of the pair
                        p_split_node = rSkinSubModelPart.pGetNode(split_pair.first);
                    } else if (dot_pair2 > 0) {
                        // Use second of the pair
                        p_split_node = rSkinSubModelPart.pGetNode(split_pair.second);
                    } else {
                        p_split_node = rSkinSubModelPart.pGetNode(split_pair.first);
                        KRATOS_WARNING("Neither candidate split nodes satisfy the inner product condition.\n");
                    }

                    if (pair_layer_info.Found) {
                        AddLayerConditionMetadata(*p_split_node, pair_layer_info.LayerName, pair_layer_info.ConditionName);
                    }

                    surrogate_segment_nodes.push_back(
                        create_projected_clone(p_second_for_fallback->Coordinates(), p_split_node->Id()));

                }

                surrogate_segment_nodes.push_back(p_second_for_fallback);
            }

            const std::size_t segment_count =
                surrogate_segment_nodes.size() > 0 ? surrogate_segment_nodes.size() - 1 : 0;
            KRATOS_ERROR_IF(segment_count == 0)
                << "::[SnakeGapSbmProcess]:: Missing surrogate segment definition for condition." << std::endl;

            for (std::size_t segment = 0; segment < segment_count; ++segment) {
                const auto& p_first = surrogate_segment_nodes[segment];
                const auto& p_second = surrogate_segment_nodes[segment + 1];

                CreateGapAndSkinQuadraturePoints<TIsInnerLoop>(
                    integration_parameters,
                    p_nurbs_surface,
                    p_first,
                    p_second,
                    surrogate_brep_middle_geometry,
                    *mpIgaModelPart,
                    rSkinSubModelPart);
            }

        };

        process_pair(p_surrogate_node_1, p_surrogate_node_2);
    }
    //---------------------------------------------------------------------------
    bool is_entering = false;
    for (IndexType i_cond_id = first_condition_id; i_cond_id <= last_condition_id; ++i_cond_id) {
        is_entering = !is_entering;

        const auto& surrogate_condition = rSurrogateSubModelPart.pGetCondition(i_cond_id);

        const auto& p_surrogate_1 = surrogate_condition->GetGeometry()(0);
        const auto& p_surrogate_2 = surrogate_condition->GetGeometry()(1);

        const auto surrogate_middle_point = surrogate_condition->GetGeometry().Center();

        const bool is_first_surrogate_already_computed = (p_surrogate_1->GetValue(ACTIVATION_LEVEL) == 1);
        const bool is_second_surrogate_already_computed = (p_surrogate_2->GetValue(ACTIVATION_LEVEL) == 1);

        // Project the surrogate vertices to the skin boundary
        // search the projection of the first vertex
        if (!is_first_surrogate_already_computed) {
            p_surrogate_1->SetValue(ACTIVATION_LEVEL, 1.0);
            const IndexType id_closest_true_node = p_surrogate_1->GetValue(PROJECTION_NODE_ID);

            const auto skin_vertex_1 = rSkinSubModelPart.pGetNode(id_closest_true_node);

            Vector active_range_knot_vector = ZeroVector(2);
            active_range_knot_vector[0] = 0;
            active_range_knot_vector[1] = 1;

            // surrogate_1 - skin_1
            // create the brep connecting vertex and closest true point
            Point surrogate_1(*p_surrogate_1);
            Point skin_1(*skin_vertex_1);

            double characteristic_condition_length = norm_2(surrogate_1-skin_1)/2;

            Node::Pointer p_first_point = Kratos::make_intrusive<Node>(1, surrogate_1);
            Node::Pointer p_second_point = Kratos::make_intrusive<Node>(2, skin_1);

            if (is_entering) {
                // change the order to preserve the anticlockwise orientation
                Node::Pointer p_temp_pointer = p_first_point;
                p_first_point = p_second_point;
                p_second_point = p_temp_pointer;
            }

            auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_first_point, p_second_point, active_range_knot_vector);
            auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);

            IntegrationInfo brep_integration_info_surrogate1_skin1 = p_brep_curve_surrogate1_skin1->GetDefaultIntegrationInfo();

            brep_integration_info_surrogate1_skin1.SetNumberOfIntegrationPointsPerSpan(0,2*mGapInterpolationOrder+1);

            IntegrationPointsArrayType brep_integration_points_list_surrogate1_skin1;
            GeometriesArrayType brep_quadrature_point_list_surrogate1_skin1;

            p_brep_curve_surrogate1_skin1->CreateIntegrationPoints(brep_integration_points_list_surrogate1_skin1, brep_integration_info_surrogate1_skin1);

            const double brep_curve_surrogate1_skin1_length = norm_2(skin_1 - surrogate_1);
            for (auto& integration_point : brep_integration_points_list_surrogate1_skin1) {
                integration_point.SetWeight(integration_point.Weight() * brep_curve_surrogate1_skin1_length);
            }

            p_brep_curve_surrogate1_skin1->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate1_skin1, number_of_shape_functions_derivatives,
                                                                brep_integration_points_list_surrogate1_skin1, brep_integration_info_surrogate1_skin1);

            std::size_t id = 1;
            if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
                id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

            auto& neighbour_geometries = p_surrogate_1->GetValue(NEIGHBOUR_GEOMETRIES);

            if (norm_2(surrogate_middle_point-neighbour_geometries[0]->Center()) < 1e-12) {
                // neighbour_geometries already disposed in the correct way
            } else if (norm_2(surrogate_middle_point-neighbour_geometries[1]->Center()) < 1e-12) {
                // neighbour_geometries not disposed in the correct way: reverse the entries
                std::swap(neighbour_geometries[0], neighbour_geometries[1]);

                p_surrogate_1->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
            } else {
                KRATOS_ERROR << "::[SnakeSbmProcess]:: The surrogate middle point is not close to any of the neighbour geometries." << std::endl;
            }
            this->CreateConditions(
                brep_quadrature_point_list_surrogate1_skin1.ptr_begin(), brep_quadrature_point_list_surrogate1_skin1.ptr_end(),
                *mpGapInterfaceSubModelPart, mGapInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries, characteristic_condition_length);

        }

        if (!is_second_surrogate_already_computed) {
            p_surrogate_2->SetValue(ACTIVATION_LEVEL, 1.0);

            const IndexType id_closest_true_node = p_surrogate_2->GetValue(PROJECTION_NODE_ID);
            const auto skin_vertex_2 = rSkinSubModelPart.pGetNode(id_closest_true_node);

            Vector active_range_knot_vector = ZeroVector(2);
            active_range_knot_vector[0] = 0;
            active_range_knot_vector[1] = 1;

            // surrogate_2 - skin_2
            // create the brep connecting vertex and closest true point
            Point surrogate_2(*p_surrogate_2);
            Point skin_2(*skin_vertex_2);

            double characteristic_condition_length = norm_2(surrogate_2-skin_2)/2;

            Node::Pointer p_first_point = Kratos::make_intrusive<Node>(1, surrogate_2);
            Node::Pointer p_second_point = Kratos::make_intrusive<Node>(2, skin_2);

            // check if the interface length is zero -> not necessary to create quadrature points there
            double min_span_size = std::min(knot_span_sizes[0], knot_span_sizes[1]);
            if (norm_2(surrogate_2-skin_2) < 1e-10*min_span_size) continue;

            if (!is_entering) {
                // change the order to preserve the anticlockwise orientation
                Node::Pointer p_temp_pointer = p_first_point;
                p_first_point = p_second_point;
                p_second_point = p_temp_pointer;
            }

            // change the order to preserve the anticlockwise orientation
            auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(p_first_point, p_second_point, active_range_knot_vector);
            auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);

            IntegrationInfo brep_integration_info_surrogate2_skin2 = p_brep_curve_surrogate2_skin2->GetDefaultIntegrationInfo();

            brep_integration_info_surrogate2_skin2.SetNumberOfIntegrationPointsPerSpan(0,2*mGapInterpolationOrder+1);

            IntegrationPointsArrayType brep_integration_points_list_surrogate2_skin2;
            GeometriesArrayType brep_quadrature_point_list_surrogate2_skin2;

            p_brep_curve_surrogate2_skin2->CreateIntegrationPoints(brep_integration_points_list_surrogate2_skin2, brep_integration_info_surrogate2_skin2);

            const double brep_curve_surrogate2_skin2_length = norm_2(skin_2 - surrogate_2);
            for (auto& integration_point : brep_integration_points_list_surrogate2_skin2) {
                integration_point.SetWeight(integration_point.Weight() * brep_curve_surrogate2_skin2_length);
            }

            p_brep_curve_surrogate2_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate2_skin2, number_of_shape_functions_derivatives,
                                                                brep_integration_points_list_surrogate2_skin2, brep_integration_info_surrogate2_skin2);

            std::size_t id = 1;
            if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
                id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

            auto& neighbour_geometries = p_surrogate_2->GetValue(NEIGHBOUR_GEOMETRIES);

            if (norm_2(surrogate_middle_point-neighbour_geometries[0]->Center()) < 1e-12) {
                // neighbour_geometries already disposed in the correct way

            } else if (norm_2(surrogate_middle_point-neighbour_geometries[1]->Center()) < 1e-12) {
                // neighbour_geometries not disposed in the correct way: reverse the entries
                std::swap(neighbour_geometries[0], neighbour_geometries[1]);

                p_surrogate_2->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
            } else {
                KRATOS_ERROR << "::[SnakeSbmProcess]:: The surrogate middle point is not close to any of the neighbour geometries." << std::endl;
            }

            this->CreateConditions(
                brep_quadrature_point_list_surrogate2_skin2.ptr_begin(), brep_quadrature_point_list_surrogate2_skin2.ptr_end(),
                *mpGapInterfaceSubModelPart, mGapInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries, characteristic_condition_length);

        }
    }

}

//-------------------------------



// Project a skin-segment midpoint and check the admissible side of the surrogate.
template <bool TIsInnerLoop>
SnakeGapSbmProcess::ProjectionResult SnakeGapSbmProcess::ComputeSingleProjection(
    const ModelPart& rSkinSubModelPart,
    const CoordinatesArrayType& rSurrogateMidPoint,
    const IndexType LeftProjectionId,
    const IndexType RightProjectionId,
    const std::string& rCommonLayerName,
    const std::string& rCommonConditionName,
    const array_1d<double, 3>& rNormalCond,
    const IntegrationParameters& rIntegrationParameters,
    const KnotSpanSkinBinsCSR& rSkinBinsPerKnotSpan,
    const std::vector<IndexType>& rConsideredNodeIds) const
{
    ProjectionResult projection_result;

    const auto& r_left_skin_node = rSkinSubModelPart.GetNode(LeftProjectionId);
    const auto& r_right_skin_node = rSkinSubModelPart.GetNode(RightProjectionId);

    CoordinatesArrayType skin_mid_point_coords = 0.5 * (r_left_skin_node.Coordinates() + r_right_skin_node.Coordinates());

    Vector normal_direction = r_right_skin_node.Coordinates() - r_left_skin_node.Coordinates();
    const double direction_norm = norm_2(normal_direction);
    if (direction_norm > 1.0e-16) {
        normal_direction /= direction_norm;
    }

    if constexpr (!TIsInnerLoop) {
        normal_direction *= -1;
    }

    const double tmp = normal_direction[0];
    normal_direction[0] = -normal_direction[1];
    normal_direction[1] = tmp;
    const double normal_norm = norm_2(normal_direction);
    if (normal_norm > 1.0e-16) {
        normal_direction /= normal_norm;
    }

    projection_result.SuccessfulProjection = true;

    IndexType id_skin_node = LeftProjectionId;
    if (norm_2(r_left_skin_node.Coordinates() - r_right_skin_node.Coordinates()) >= rIntegrationParameters.rKnotSpanSizes[0] / 1e5) {
        id_skin_node = FindClosestNodeInLayerWithDirection<TIsInnerLoop>(
            skin_mid_point_coords,
            rCommonLayerName,
            rSkinSubModelPart,
            rIntegrationParameters.rKnotSpanSizes,
            rSkinBinsPerKnotSpan,
            normal_direction,
            &rConsideredNodeIds);

        auto candidate_surrogate_to_skin_node = rSkinSubModelPart.GetNode(id_skin_node).Coordinates() - rSurrogateMidPoint;
        const double dot_product = inner_prod(candidate_surrogate_to_skin_node, rNormalCond);

        if (dot_product < -1e-12) {
            const auto split_pair = FindClosestPairInLayerWithNormalDirection<TIsInnerLoop>(
                skin_mid_point_coords,
                rCommonLayerName,
                rSkinSubModelPart,
                rIntegrationParameters.rKnotSpanSizes,
                rSkinBinsPerKnotSpan,
                normal_direction, nullptr, nullptr, false, true);

            auto candidate_surrogate_to_skin_node_1 = rSkinSubModelPart.GetNode(split_pair.first).Coordinates() - rSurrogateMidPoint;
            auto candidate_surrogate_to_skin_node_2 = rSkinSubModelPart.GetNode(split_pair.second).Coordinates() - rSurrogateMidPoint;
            const double dot_prod_1 = inner_prod(candidate_surrogate_to_skin_node_1, rNormalCond);
            const double dot_prod_2 = inner_prod(candidate_surrogate_to_skin_node_2, rNormalCond);
            if (dot_prod_1 >= 0 && dot_prod_2 < 0) {
                id_skin_node = split_pair.first;
            } else if (dot_prod_2 >= 0 && dot_prod_1 < 0) {
                id_skin_node = split_pair.second;
            } else {
                projection_result.SuccessfulProjection = false;

                if (dot_prod_1 > dot_prod_2 )
                    id_skin_node = split_pair.first;
                else
                    id_skin_node = split_pair.second;
            }
        }
    }

    projection_result.ProjectionId = id_skin_node;
    auto p_skin_node = rSkinSubModelPart.pGetNode(id_skin_node);
    AddLayerConditionMetadata(*p_skin_node, rCommonLayerName, rCommonConditionName);
    projection_result.ProjectionPoint = rSkinSubModelPart.GetNode(id_skin_node).Coordinates();

    return projection_result;
}

// Build sparse knot-span occupancy and local search bins for skin nodes and conditions.
SnakeGapSbmProcess::KnotSpanSkinBinsCSR
SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D(
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart) const
{
    KnotSpanSkinBinsCSR knot_span_data;

    const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();

    const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);
    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] KNOT_SPAN_SIZES must have at least two entries.\n";
    const auto bounds = ReadParameterSpaceBounds(
        r_parent_model_part,
        "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D");

    const double span_size_x = knot_span_sizes[0];
    const double span_size_y = knot_span_sizes[1];
    KRATOS_ERROR_IF(span_size_x <= 0.0 || span_size_y <= 0.0)
        << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] Knot span sizes must be positive.\n";

    const double min_u = bounds.MinU;
    const double max_u = bounds.MaxU;
    const double min_v = bounds.MinV;
    const double max_v = bounds.MaxV;

    const double domain_length_u = max_u - min_u;
    const double domain_length_v = max_v - min_v;
    KRATOS_ERROR_IF(domain_length_u <= 0.0 || domain_length_v <= 0.0)
        << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] Invalid parameter space extents.\n";

    const std::size_t number_of_spans_x = ComputeSpanCount(
        domain_length_u,
        span_size_x,
        "u",
        "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D");
    const std::size_t number_of_spans_y = ComputeSpanCount(
        domain_length_v,
        span_size_y,
        "v",
        "SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D");

    knot_span_data.NumberOfSpansX = number_of_spans_x;
    knot_span_data.NumberOfSpansY = number_of_spans_y;
    knot_span_data.NumberOfSpansZ = 1;
    knot_span_data.MinU = min_u;
    knot_span_data.MaxU = max_u;
    knot_span_data.MinV = min_v;
    knot_span_data.MaxV = max_v;
    knot_span_data.MinW = 0.0;
    knot_span_data.MaxW = 0.0;
    knot_span_data.SpanSizeX = span_size_x;
    knot_span_data.SpanSizeY = span_size_y;
    knot_span_data.SpanSizeZ = 0.0;

    if (rSkinSubModelPart.NumberOfNodes() == 0 && rSkinSubModelPart.NumberOfConditions() == 0) {
        knot_span_data.Occupancy.resize(number_of_spans_x, number_of_spans_y, false);
        return knot_span_data;
    }

    const double tolerance = 1.0e-10;
    const double max_u_with_tolerance = max_u + tolerance;
    const double max_v_with_tolerance = max_v + tolerance;
    const double min_u_with_tolerance = min_u - tolerance;
    const double min_v_with_tolerance = min_v - tolerance;

    const auto compute_span_index = [tolerance](
        double coordinate,
        double min_value,
        double max_value,
        double span_size,
        std::size_t span_count) {
        double clamped = coordinate;
        if (coordinate < min_value) {
            KRATOS_ERROR_IF(coordinate < min_value - tolerance)
                << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] coordinate below minimum parameter range.\n";
            clamped = min_value;
        } else if (coordinate > max_value) {
            KRATOS_ERROR_IF(coordinate > max_value + tolerance)
                << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] coordinate above maximum parameter range.\n";
            clamped = max_value;
        }
        const double relative = (clamped - min_value) / span_size;
        std::size_t idx = static_cast<std::size_t>(std::floor(relative + tolerance));
        if (idx >= span_count) idx = span_count - 1;
        return idx;
    };

    const auto get_touched_span_range = [tolerance](
        const double MinCoordinate,
        const double MaxCoordinate,
        const double MinDomainCoordinate,
        const double MaxDomainCoordinate,
        const double SpanSize,
        const std::size_t SpanCount) {
        const double min_clamped = std::max(MinDomainCoordinate, std::min(MaxDomainCoordinate, MinCoordinate));
        const double max_clamped = std::max(MinDomainCoordinate, std::min(MaxDomainCoordinate, MaxCoordinate));
        const double min_relative = (min_clamped - MinDomainCoordinate) / SpanSize;
        const double max_relative = (max_clamped - MinDomainCoordinate) / SpanSize;
        std::size_t first_span = static_cast<std::size_t>(std::floor(min_relative + tolerance));
        std::size_t last_span = static_cast<std::size_t>(std::floor(max_relative + tolerance));
        if (first_span > 0 && std::abs(min_relative - std::round(min_relative)) <= tolerance) {
            --first_span;
        }
        first_span = std::min(first_span, SpanCount - 1);
        last_span = std::min(last_span, SpanCount - 1);
        return std::make_pair(first_span, last_span);
    };
    const auto for_each_condition_span = [&](const Condition& rCondition, const auto& rCallback) {
        const auto& r_geometry = rCondition.GetGeometry();
        if (r_geometry.size() == 0) {
            return;
        }
        double min_condition_x = std::numeric_limits<double>::max();
        double max_condition_x = std::numeric_limits<double>::lowest();
        double min_condition_y = std::numeric_limits<double>::max();
        double max_condition_y = std::numeric_limits<double>::lowest();
        for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
            const auto& r_coordinates = r_geometry[point_index].Coordinates();
            min_condition_x = std::min(min_condition_x, r_coordinates[0]);
            max_condition_x = std::max(max_condition_x, r_coordinates[0]);
            min_condition_y = std::min(min_condition_y, r_coordinates[1]);
            max_condition_y = std::max(max_condition_y, r_coordinates[1]);
        }
        const auto x_range = get_touched_span_range(
            min_condition_x, max_condition_x, min_u, max_u, span_size_x, number_of_spans_x);
        const auto y_range = get_touched_span_range(
            min_condition_y, max_condition_y, min_v, max_v, span_size_y, number_of_spans_y);
        for (std::size_t span_index_x = x_range.first; span_index_x <= x_range.second; ++span_index_x) {
            for (std::size_t span_index_y = y_range.first; span_index_y <= y_range.second; ++span_index_y) {
                rCallback(span_index_x, span_index_y);
            }
        }
    };

    // --- Pass 1: gather unique columns per row for sparsity pattern ---
    std::vector<std::vector<std::size_t>> column_indices_per_row(number_of_spans_x);

    auto add_span_to_pattern = [&](const std::size_t SpanIndexX, const std::size_t SpanIndexY) {
        column_indices_per_row[SpanIndexX].push_back(
            FlattenSpanColumnIndex(knot_span_data, SpanIndexY));
    };

    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const double u = r_node.X();
        const double v = r_node.Y();

        KRATOS_ERROR_IF(u < min_u_with_tolerance || u > max_u_with_tolerance)
            << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] node " << r_node.Id()
            << " has u-parameter " << u << " outside [" << min_u << ", " << max_u << "].\n";
        KRATOS_ERROR_IF(v < min_v_with_tolerance || v > max_v_with_tolerance)
            << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] node " << r_node.Id()
            << " has v-parameter " << v << " outside [" << min_v << ", " << max_v << "].\n";

        const std::size_t span_index_x = compute_span_index(u, min_u, max_u, span_size_x, number_of_spans_x);
        const std::size_t span_index_y = compute_span_index(v, min_v, max_v, span_size_y, number_of_spans_y);
        add_span_to_pattern(span_index_x, span_index_y);
    }

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        for_each_condition_span(r_condition, add_span_to_pattern);
    }

    std::size_t number_of_non_zero_entries = 0;
    for (auto& r_column_indices : column_indices_per_row) {
        std::sort(r_column_indices.begin(), r_column_indices.end());
        r_column_indices.erase(std::unique(r_column_indices.begin(), r_column_indices.end()), r_column_indices.end());
        number_of_non_zero_entries += static_cast<std::size_t>(r_column_indices.size());
    }

    // --- Materialize CSR matrix with zero values and build nnz slots ---
    auto& r_occupancy_matrix = knot_span_data.Occupancy;
    r_occupancy_matrix.resize(number_of_spans_x, number_of_spans_y, false);
    r_occupancy_matrix.reserve(number_of_non_zero_entries);
    std::vector<std::unordered_map<std::size_t, std::size_t>> occupancy_index_lookup(number_of_spans_x);
    for (auto& r_lookup : occupancy_index_lookup) {
        r_lookup.reserve(8);
    }

    std::size_t nnz_counter = 0;
    for (std::size_t span_index_x = 0; span_index_x < number_of_spans_x; ++span_index_x) {
        for (const std::size_t span_column_index : column_indices_per_row[span_index_x]) {
            r_occupancy_matrix.push_back(span_index_x, span_column_index, 0.0);
            occupancy_index_lookup[span_index_x].emplace(span_column_index, nnz_counter++);
        }
    }

    r_occupancy_matrix.complete_index1_data();

    std::vector<NodePointerContainerType> nodes_per_non_zero(number_of_non_zero_entries);
    std::vector<ConditionPointerContainerType> conditions_per_non_zero(number_of_non_zero_entries);

    auto find_non_zero_index = [&](const std::size_t SpanIndexX, const std::size_t SpanIndexY) {
        const std::size_t span_column_index = FlattenSpanColumnIndex(knot_span_data, SpanIndexY);
        const auto nnz_it = occupancy_index_lookup[SpanIndexX].find(span_column_index);
        return nnz_it != occupancy_index_lookup[SpanIndexX].end() ? nnz_it->second : static_cast<std::size_t>(-1);
    };

    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const std::size_t span_index_x = compute_span_index(r_node.X(), min_u, max_u, span_size_x, number_of_spans_x);
        const std::size_t span_index_y = compute_span_index(r_node.Y(), min_v, max_v, span_size_y, number_of_spans_y);

        const std::size_t non_zero_index = find_non_zero_index(span_index_x, span_index_y);
        KRATOS_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
            << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] node nonzero (ix,iy) not found in CSR pattern.\n";

        nodes_per_non_zero[non_zero_index].push_back(rSkinSubModelPart.pGetNode(r_node.Id()));
        r_occupancy_matrix.value_data()[non_zero_index] += 1.0;
    }

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        for_each_condition_span(r_condition, [&](const std::size_t SpanIndexX, const std::size_t SpanIndexY) {
            const std::size_t non_zero_index = find_non_zero_index(SpanIndexX, SpanIndexY);
            KRATOS_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
                << "[SnakeGapSbmProcess::CreateSkinBinsPerKnotSpanMatrix2D] condition nonzero (ix,iy) not found in CSR pattern.\n";

            conditions_per_non_zero[non_zero_index].push_back(rSkinSubModelPart.pGetCondition(r_condition.Id()));
            r_occupancy_matrix.value_data()[non_zero_index] += 1.0;
        });
    }

    // Build local search structures for each occupied span.
    knot_span_data.CellBinsByNnz.resize(number_of_non_zero_entries);
    for (std::size_t k = 0; k < number_of_non_zero_entries; ++k) {
        auto& r_cell_bins = knot_span_data.CellBinsByNnz[k];
        auto& r_nodes = nodes_per_non_zero[k];

        auto& r_conditions = conditions_per_non_zero[k];
        if (!r_nodes.empty()) {
            std::sort(r_nodes.begin(), r_nodes.end(),
                      [](const Node::Pointer& a, const Node::Pointer& b) {
                          return a->Id() < b->Id();
                      });
            r_nodes.erase(std::unique(r_nodes.begin(), r_nodes.end(),
                                      [](const Node::Pointer& a, const Node::Pointer& b) {
                                          return a->Id() == b->Id();
                                      }),
                          r_nodes.end());

            r_cell_bins.Nodes = std::move(r_nodes);
            r_cell_bins.pNodeBins = std::make_unique<NodeBinsType>(
                r_cell_bins.Nodes.begin(), r_cell_bins.Nodes.end());
            r_cell_bins.HasNodeBins = true;
        }

        if (!r_conditions.empty()) {
            std::sort(r_conditions.begin(), r_conditions.end(),
                      [](const ConditionPointerType& a, const ConditionPointerType& b) {
                          return a->Id() < b->Id();
                      });
            r_conditions.erase(std::unique(r_conditions.begin(), r_conditions.end(),
                                           [](const ConditionPointerType& a, const ConditionPointerType& b) {
                                               return a->Id() == b->Id();
                                           }),
                               r_conditions.end());

            r_cell_bins.Conditions = std::move(r_conditions);
            r_cell_bins.ConditionBins = BinsObjectDynamic<ConditionConfigure>(
                r_cell_bins.Conditions.begin(), r_cell_bins.Conditions.end());
            r_cell_bins.HasConditionBins = true;
        }
    }

    return knot_span_data;
}

// Connect projections on different layers through shared interface vertices.
std::vector<IndexType> SnakeGapSbmProcess::AssessProjectionsFeasibility(
    const ModelPart& rSkinSubModelPart,
    Node::Pointer pSurrogateNode1,
    Node::Pointer pSurrogateNode2)
{
    const IndexType projection_id_1 = pSurrogateNode1->GetValue(PROJECTION_NODE_ID);
    const IndexType projection_id_2 = pSurrogateNode2->GetValue(PROJECTION_NODE_ID);

    const auto& r_projection_node_1 = rSkinSubModelPart.GetNode(projection_id_1);
    const auto& r_projection_node_2 = rSkinSubModelPart.GetNode(projection_id_2);

    const auto& connected_layers_surrogate_node_1 = r_projection_node_1.GetValue(CONNECTED_LAYERS);
    const auto& connected_layers_surrogate_node_2 = r_projection_node_2.GetValue(CONNECTED_LAYERS);

    auto share_layer = [](const std::vector<std::string>& rLayers1, const std::vector<std::string>& rLayers2) {
        for (const auto& layer : rLayers1) {
            if (std::find(rLayers2.begin(), rLayers2.end(), layer) != rLayers2.end()) {
                return true;
            }
        }
        return false;
    };

    // Start from the first projection and insert interface vertices only when required.
    std::vector<IndexType> projection_chain{projection_id_1};
    if (share_layer(connected_layers_surrogate_node_1, connected_layers_surrogate_node_2)) {
        if (projection_id_2 != projection_id_1) {
            projection_chain.push_back(projection_id_2);
        }
        return projection_chain;
    }

    const ModelPart::NodesContainerType* p_interface_nodes = nullptr;
    if (rSkinSubModelPart.HasSubModelPart("interface_vertices")) {
        p_interface_nodes = &rSkinSubModelPart.GetSubModelPart("interface_vertices").Nodes();
    } else {
        p_interface_nodes = &rSkinSubModelPart.Nodes();
    }

    std::unordered_set<IndexType> visited_ids;
    visited_ids.insert(projection_id_1);

    IndexType current_projection_id = projection_id_1;
    const auto surrogate_target_coordinates = pSurrogateNode2->Coordinates();
    const std::size_t max_iterations = p_interface_nodes->size() + 2;
    std::size_t iter = 0;

    while (!share_layer(rSkinSubModelPart.GetNode(current_projection_id).GetValue(CONNECTED_LAYERS),
                        connected_layers_surrogate_node_2))
    {
        iter++;
        KRATOS_ERROR_IF(iter > max_iterations)
            << "::[SnakeGapSbmProcess]::AssessProjectionsFeasibility: Unable to connect layers from projection nodes "
            << projection_id_1 << " to " << projection_id_2 << " without cycles." << std::endl;

        const auto& current_layers = rSkinSubModelPart.GetNode(current_projection_id).GetValue(CONNECTED_LAYERS);

        IndexType best_candidate_id = std::numeric_limits<IndexType>::max();
        double best_distance = std::numeric_limits<double>::max();

        for (const auto& r_interface_node : *p_interface_nodes) {
            const IndexType candidate_id = r_interface_node.Id();
            if (visited_ids.find(candidate_id) != visited_ids.end()) {
                continue;
            }

            const auto& candidate_layers = r_interface_node.GetValue(CONNECTED_LAYERS);
            if (!share_layer(current_layers, candidate_layers)) {
                continue;
            }

            array_1d<double, 3> diff = r_interface_node.Coordinates();
            diff -= surrogate_target_coordinates;
            const double distance = norm_2(diff);

            if (distance < best_distance) {
                best_distance = distance;
                best_candidate_id = candidate_id;
            }
        }

        KRATOS_ERROR_IF(best_candidate_id == std::numeric_limits<IndexType>::max())
            << "::[SnakeGapSbmProcess]::AssessProjectionsFeasibility: no interface node shares layers with projection node "
            << current_projection_id << " while targeting surrogate node " << pSurrogateNode2->Id() << std::endl;

        projection_chain.push_back(best_candidate_id);
        visited_ids.insert(best_candidate_id);
        current_projection_id = best_candidate_id;

        if (share_layer(rSkinSubModelPart.GetNode(current_projection_id).GetValue(CONNECTED_LAYERS),
                        connected_layers_surrogate_node_2)) {
            if (current_projection_id != projection_id_2) {
                projection_chain.push_back(projection_id_2);
            }
            break;
        }
    }

    if (projection_chain.back() != projection_id_2) {
        projection_chain.push_back(projection_id_2);
    }

    auto common_layer_between = [&](IndexType id1, IndexType id2) {
        const auto& layers_1 = rSkinSubModelPart.GetNode(id1).GetValue(CONNECTED_LAYERS);
        const auto& layers_2 = rSkinSubModelPart.GetNode(id2).GetValue(CONNECTED_LAYERS);
        for (const auto& layer : layers_1) {
            if (std::find(layers_2.begin(), layers_2.end(), layer) != layers_2.end()) {
                return layer;
            }
        }
        return std::string();
    };

    // Remove intermediate vertices when both adjacent segments belong to the same layer.
    if (projection_chain.size() >= 3) {
        std::size_t i = 0;
        while (i + 2 < projection_chain.size()) {
            const auto layer_1 = common_layer_between(projection_chain[i], projection_chain[i + 1]);
            const auto layer_2 = common_layer_between(projection_chain[i + 1], projection_chain[i + 2]);
            if (!layer_1.empty() && !layer_2.empty() && layer_1 == layer_2) {
                projection_chain.erase(projection_chain.begin() + static_cast<std::ptrdiff_t>(i + 1));
                if (i > 0) {
                    --i;
                }
                continue;
            }
            ++i;
        }
    }

    return projection_chain;
}

// Reconstruct the gap and its skin boundary, with optional recursive subdivision.
template <bool TIsInnerLoop>
void SnakeGapSbmProcess::CreateGapAndSkinQuadraturePoints(
    IntegrationParameters& rIntegrationParameters,
    NurbsSurfaceType::Pointer& pNurbsSurface,
    const Node::Pointer& pSurrogateNode1,
    const Node::Pointer& pSurrogateNode2,
    const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
    ModelPart& rIgaModelPart,
    ModelPart& rSkinSubModelPart)
{
    const std::size_t gap_quadrature_order = std::max(
        mGapInterpolationOrder, mGapApproximationOrder);

    const IndexType id_closest_true_node = pSurrogateNode1->GetValue(PROJECTION_NODE_ID);

    const auto& p_skin_node_1 = rSkinSubModelPart.pGetNode(id_closest_true_node);

    // search the projection of the second vertex
    const IndexType id_closest_true_node_2 = pSurrogateNode2->GetValue(PROJECTION_NODE_ID);

    const auto& p_skin_node_2 = rSkinSubModelPart.pGetNode(id_closest_true_node_2);

    const auto common_layer_info = FindCommonLayerInfo(*pSurrogateNode1, *pSurrogateNode2);
    KRATOS_ERROR_IF(!common_layer_info.Found) << ":::[SnakeGapSbmProcess]::: No common layer found between the two surrogate nodes "
                                    << pSurrogateNode1->Id() << " and " << pSurrogateNode2->Id() <<  "\n"
                                    << pSurrogateNode1->Coordinates() << pSurrogateNode2->Coordinates() << "\n"
                                    << pSurrogateNode1->GetValue(CONNECTED_LAYERS)
                                    << pSurrogateNode2->GetValue(CONNECTED_LAYERS) << std::endl;

    const std::string& layer_name = common_layer_info.LayerName;
    const std::string& condition_name = common_layer_info.ConditionName;
    KRATOS_ERROR_IF(condition_name.empty())
        << ":::[SnakeGapSbmProcess]::: Missing condition metadata for common layer '"
        << layer_name << "' between surrogate nodes "
        << pSurrogateNode1->Id() << " and " << pSurrogateNode2->Id() << "." << std::endl;

    KRATOS_ERROR_IF(!rIntegrationParameters.pSkinBinsPerSpan)
        << "::[SnakeGapSbmProcess]::CreateGapAndSkinQuadraturePoints: skin span bins are not initialized." << std::endl;

    const auto& r_skin_bins_per_span = *rIntegrationParameters.pSkinBinsPerSpan;

    ModelPart& r_layer_model_part = rIgaModelPart.HasSubModelPart(layer_name) ?
                                    rIgaModelPart.GetSubModelPart(layer_name) :
                                    rIgaModelPart.CreateSubModelPart(layer_name);

    Vector active_range_knot_vector = ZeroVector(2);
    active_range_knot_vector[0] = 0;
    active_range_knot_vector[1] = 1;

    // surrogate_1 - skin_1
    // create the brep connecting vertex and closest true point
    Point skin_1(*p_skin_node_1);
    Point skin_2(*p_skin_node_2);

    Node::Pointer p_skin1_brep_point = Kratos::make_intrusive<Node>(2, skin_1);
    Node::Pointer p_skin2_brep_point = Kratos::make_intrusive<Node>(2, skin_2);

    auto p_nurbs_curve_skin1_skin2_linear = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);
    auto p_brep_curve_skin1_skin2_linear = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2_linear);

    std::vector<Geometry<Node>::Pointer> neighbour_geometries_skin1_skin2;
    neighbour_geometries_skin1_skin2.push_back(rSurrogateBrepMiddleGeometry);

    const double subdivision_threshold = rIntegrationParameters.rKnotSpanSizes[0] * mGapRelativeToleranceForSubdivisions;
    const double subdivision_threshold_sq = subdivision_threshold * subdivision_threshold;
    std::size_t max_internal_divisions = mInternalDivisions > 0 ? static_cast<std::size_t>(mInternalDivisions) : 0;

    // Reconstruct the skin curve from either recursive projections or ordered skin samples.
    auto build_skin_curve = [&](IndexType left_id,
                                IndexType right_id,
                                std::vector<IndexType>& rInterpolationNodes) {
        rInterpolationNodes.clear();
        const auto& r_skin_left = rSkinSubModelPart.GetNode(left_id);
        const auto& r_skin_right = rSkinSubModelPart.GetNode(right_id);

        auto p_left_point = Kratos::make_intrusive<Node>(2, r_skin_left.Coordinates());
        auto p_right_point = Kratos::make_intrusive<Node>(2, r_skin_right.Coordinates());

        auto p_nurbs_curve_skin = this->CreateBrepCurve(p_left_point, p_right_point, active_range_knot_vector);

        int p = mGapApproximationOrder;

        if (p == 1 && mGapSbmType != "default") {
            rInterpolationNodes = {left_id, right_id};
            return p_nurbs_curve_skin;
        }

        bool linear_interpolation_is_enough = false;
        std::vector<array_1d<double, 3>> interpolation_points;

        if (mGapSbmType != "default") {
            std::size_t number_interpolation_cuts = mNumberOfInterpolationLevels;
            std::size_t segment_count = number_interpolation_cuts == 0 ? 1 : static_cast<std::size_t>(1) << number_interpolation_cuts;
            const IndexType invalid_projection_id = std::numeric_limits<IndexType>::max();

            std::vector<IndexType> interpolation_projection_ids(segment_count + 1, invalid_projection_id);
            interpolation_projection_ids[0] = left_id;
            interpolation_projection_ids[segment_count] = right_id;
            const std::size_t recursion_segment_count = segment_count;

            auto compute_recursive_interpolation = [&](auto&& self,
                                                       std::size_t depth,
                                                       std::size_t left_index,
                                                       std::size_t right_index,
                                                       IndexType left_node_id,
                                                       IndexType right_node_id) -> void
            {
                if (depth == 0) {
                    return;
                }

                const std::size_t mid_index = (left_index + right_index) / 2;
                if (interpolation_projection_ids[mid_index] != invalid_projection_id) {
                    if (depth > 1) {
                        self(self, depth - 1, left_index, mid_index, left_node_id, interpolation_projection_ids[mid_index]);
                        self(self, depth - 1, mid_index, right_index, interpolation_projection_ids[mid_index], right_node_id);
                    }
                    return;
                }

                const auto& r_left_skin_node = rSkinSubModelPart.GetNode(left_node_id);
                const auto& r_right_skin_node = rSkinSubModelPart.GetNode(right_node_id);

                CoordinatesArrayType skin_mid_point_coords = 0.5 * (r_left_skin_node.Coordinates() + r_right_skin_node.Coordinates());

                const array_1d<double, 3> tangent = r_right_skin_node.Coordinates() - r_left_skin_node.Coordinates();
                const double tangent_norm = norm_2(tangent);

                if (tangent_norm >= rIntegrationParameters.rKnotSpanSizes[0] / 1.0e5) {
                    array_1d<double, 3> normal_skin_segment;
                    normal_skin_segment[0] = tangent[1];
                    normal_skin_segment[1] = -tangent[0];
                    normal_skin_segment[2] = 0.0;

                    const auto projection = ComputeSingleProjection<TIsInnerLoop>(
                        rSkinSubModelPart,
                        skin_mid_point_coords,
                        left_node_id,
                        right_node_id,
                        layer_name,
                        condition_name,
                        normal_skin_segment,
                        rIntegrationParameters,
                        r_skin_bins_per_span,
                        interpolation_projection_ids);
                    const IndexType id_skin_node = projection.ProjectionId;
                    const bool added_new_node = std::find(
                        interpolation_projection_ids.begin(), interpolation_projection_ids.end(), id_skin_node)
                        == interpolation_projection_ids.end();
                    interpolation_projection_ids[mid_index] = id_skin_node;

                    if (added_new_node && depth > 1) {
                        self(self, depth - 1, left_index, mid_index, left_node_id, id_skin_node);
                        self(self, depth - 1, mid_index, right_index, id_skin_node, right_node_id);
                    }
                    return;
                }

                interpolation_projection_ids[mid_index] = left_node_id;
            };

            if (segment_count > 1) {
                compute_recursive_interpolation(compute_recursive_interpolation, number_interpolation_cuts, 0, recursion_segment_count,
                                                left_id, right_id);
            }

            std::unordered_set<IndexType> seen_projection_ids;

            auto new_end = std::remove_if(
                interpolation_projection_ids.begin(),
                interpolation_projection_ids.end(),
                [&](IndexType projection_id) {
                    if (projection_id == invalid_projection_id) {
                        return true;
                    }
                    return !seen_projection_ids.insert(projection_id).second;
                });

            interpolation_projection_ids.erase(new_end, interpolation_projection_ids.end());

            segment_count = interpolation_projection_ids.size();
            interpolation_points.reserve(segment_count);

            rInterpolationNodes = interpolation_projection_ids;

            for (std::size_t i = 0; i < segment_count; ++i) {
                const IndexType projection_id = interpolation_projection_ids[i];
                KRATOS_ERROR_IF(projection_id == invalid_projection_id)
                    << "::[SnakeGapSbmProcess]:: Missing interpolation projection id at subdivision index " << i << std::endl;

                const auto& r_skin_node = rSkinSubModelPart.GetNode(projection_id);
                interpolation_points.push_back(r_skin_node.Coordinates());
            }

            const double ridge = 1e-11;
            if (!(norm_2(r_skin_right.Coordinates() - r_skin_left.Coordinates()) > rIntegrationParameters.rKnotSpanSizes[0]/50 && p>1)) {
                linear_interpolation_is_enough = true;
            }

            if (!linear_interpolation_is_enough && int(segment_count) >= p+1) {
                p_nurbs_curve_skin = FitBezierUV_LS_Generic(interpolation_points, p, ridge);
            }
        }
        else //default case
        {
            const IndexType first_id = rSkinSubModelPart.NodesBegin()->Id();
            const IndexType last_id  = first_id + rSkinSubModelPart.NumberOfNodes() - 1;
            bool check_enough_samples = false;

            auto next_id = [&](IndexType id){ return (id < last_id) ? (id + 1) : first_id; };
            auto previous_id = [&](IndexType id){ return (id > first_id) ? (id - 1) : last_id; };

            auto advance_to_next_skin_id = [&](IndexType current_id) {
                if constexpr (TIsInnerLoop) {
                    return previous_id(current_id);
                } else {
                    return next_id(current_id);
                }
            };

            std::size_t id_distance = 0;
            IndexType current_id = left_id;
            ModelPart::SizeType iter = 0;
            const ModelPart::SizeType max_iterations = rSkinSubModelPart.NumberOfNodes();
            while (current_id != right_id) {
                current_id = advance_to_next_skin_id(current_id);

                ++id_distance;
                ++iter;
                KRATOS_ERROR_IF(iter > max_iterations)
                    << "::[SnakeGapSbmProcess]:: id_distance loop between node IDs "
                    << left_id << " and " << right_id
                    << " exceeded the number of skin nodes (" << max_iterations << ").\n";
            }

            if (id_distance > static_cast<std::size_t>(2*p+1)) {
                check_enough_samples = true;
            }

            if (norm_2(r_skin_right.Coordinates() - r_skin_left.Coordinates()) > rIntegrationParameters.rKnotSpanSizes[0]/50 && check_enough_samples && p>1) {
                p_nurbs_curve_skin = FitUV_BetweenSkinNodes_Generic<TIsInnerLoop>(
                    rSkinSubModelPart, *pNurbsSurface, left_id, right_id, p, /*ridge=*/1e-14);
            }
        }

        return p_nurbs_curve_skin;
    };

    // Generate skin quadrature and attach the corresponding surrogate reference geometry.
    auto create_skin_conditions = [&](const BrepCurveType::Pointer& p_brep_curve_skin,
                                      const std::vector<IndexType>& rInterpolationNodes,
                                      const Node::Pointer& p_left_skin,
                                      const Node::Pointer& p_right_skin) {
        IntegrationPointsArrayType brep_integration_points_list_skin;
        GeometriesArrayType brep_quadrature_point_list_skin;

        const int number_of_integration_points = static_cast<int>(2 * gap_quadrature_order + 1);
        rIntegrationParameters.CurveIntegrationInfo.SetNumberOfIntegrationPointsPerSpan(0, number_of_integration_points);

        p_brep_curve_skin->CreateIntegrationPoints(brep_integration_points_list_skin, rIntegrationParameters.CurveIntegrationInfo);

        for (auto& integration_point : brep_integration_points_list_skin) {
            std::vector<CoordinatesArrayType> ders;
            p_brep_curve_skin->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->GlobalSpaceDerivatives(ders, integration_point.Coordinates(), 1);
            const double weight_correcting_factor = std::sqrt(ders[1][0]*ders[1][0] + ders[1][1]*ders[1][1]);
            integration_point.SetWeight(integration_point.Weight() * weight_correcting_factor);
        }

        if (norm_2(p_right_skin->Coordinates() - p_left_skin->Coordinates()) > 1e-12) {
            p_brep_curve_skin->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin, rIntegrationParameters.NumberOfShapeFunctionsDerivatives,
                                                               brep_integration_points_list_skin, rIntegrationParameters.CurveIntegrationInfo);
        }

        for (auto& quad_geom : brep_quadrature_point_list_skin) {
            quad_geom.SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries_skin1_skin2);
        }

        if (mGapSbmType == "sbm") {
            for (auto& quad_geom: brep_quadrature_point_list_skin) {
                CoordinatesArrayType quadrature_point_coords = quad_geom.Center().Coordinates();
                auto normal = quad_geom.Normal(0);
                normal = -normal / MathUtils<double>::Norm(normal);

                const auto split_pair = FindClosestPairInLayerWithNormalDirection<TIsInnerLoop>(
                    quadrature_point_coords,
                    layer_name,
                    rSkinSubModelPart,
                    rIntegrationParameters.rKnotSpanSizes,
                    r_skin_bins_per_span,
                    normal);

                KRATOS_ERROR_IF(split_pair.first == std::numeric_limits<IndexType>::max() ||
                                split_pair.second == std::numeric_limits<IndexType>::max())
                    << "::[SnakeGapSbmProcess]::CreateGapAndSkinQuadraturePoints: invalid skin pair for layer "
                    << layer_name << " at quadrature point " << quadrature_point_coords << std::endl;

                auto p_skin_node_1_pair = rSkinSubModelPart.pGetNode(split_pair.first);
                auto p_skin_node_2_pair = rSkinSubModelPart.pGetNode(split_pair.second);

                const array_1d<double, 3> segment = p_skin_node_2_pair->Coordinates() - p_skin_node_1_pair->Coordinates();
                const double segment_length_sq = inner_prod(segment, segment);
                array_1d<double, 3> closest_point = p_skin_node_1_pair->Coordinates();
                if (segment_length_sq > 1e-16) {
                    const array_1d<double, 3> to_point = quadrature_point_coords - p_skin_node_1_pair->Coordinates();
                    double t = inner_prod(to_point, segment) / segment_length_sq;
                    t = std::max(0.0, std::min(1.0, t));
                    closest_point = p_skin_node_1_pair->Coordinates() + t * segment;
                }

                ModelPart& r_skin_sub_model_part = mpSkinModelPart->GetSubModelPart(rSkinSubModelPart.Name());
                ModelPart& r_skin_layer_sub_model_part = r_skin_sub_model_part.HasSubModelPart(layer_name) ?
                                                        r_skin_sub_model_part.GetSubModelPart(layer_name) :
                                                        r_skin_sub_model_part.CreateSubModelPart(layer_name);

                const auto& r_skin_root_model_part = r_skin_layer_sub_model_part.GetRootModelPart();
                const IndexType new_skin_node_id = r_skin_root_model_part.NumberOfNodes() > 0 ?
                                                (r_skin_root_model_part.NodesEnd() - 1)->Id() + 1 : 1;
                NodeType::Pointer p_skin_node = r_skin_layer_sub_model_part.CreateNewNode(
                    new_skin_node_id, closest_point[0], closest_point[1], closest_point[2]);
                const auto& r_normal_1_pair = p_skin_node_1_pair->GetValue(NORMAL);
                const auto& r_normal_2_pair = p_skin_node_2_pair->GetValue(NORMAL);
                array_1d<double, 3> avg_normal_pair = 0.5 * (r_normal_1_pair + r_normal_2_pair);
                avg_normal_pair = avg_normal_pair/norm_2(avg_normal_pair);
                p_skin_node->SetValue(NORMAL, avg_normal_pair);

                AddLayerConditionMetadata(*p_skin_node, layer_name, condition_name);
                NodePointerVector empty_vector;

                if constexpr (TIsInnerLoop)
                    p_skin_node->SetValue(IDENTIFIER, "inner");
                else
                    p_skin_node->SetValue(IDENTIFIER, "outer");

                empty_vector.push_back(p_skin_node); // Just it_node-plane neighbours
                NodePointerVector neighbour_nodes;

                neighbour_nodes.push_back(p_skin_node);
                quad_geom.SetValue(NEIGHBOUR_NODES, neighbour_nodes);
            }
        }

        const double characteristic_condition_length = CalculateGapElementCharacteristicLength(
            pSurrogateNode1->Coordinates(), pSurrogateNode2->Coordinates(),
            p_left_skin->Coordinates(), p_right_skin->Coordinates());

        std::size_t id = 1;
        if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
            id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

        // needed for GapSBM ContactMechanics -----------------------
        if (condition_name != "GapSbmContactCondition")
        {
            this->CreateConditions(
                brep_quadrature_point_list_skin.ptr_begin(), brep_quadrature_point_list_skin.ptr_end(),
                r_layer_model_part, condition_name, id, PropertiesPointerType(), rIntegrationParameters.rKnotSpanSizes,
                neighbour_geometries_skin1_skin2, characteristic_condition_length);

            ModelPart& r_contact_sub_model_part = r_layer_model_part.HasSubModelPart("contact") ?
                                                    r_layer_model_part.GetSubModelPart("contact") :
                                                    r_layer_model_part.CreateSubModelPart("contact");
            auto& r_root_model_part = r_layer_model_part.GetRootModelPart();
            p_brep_curve_skin->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries_skin1_skin2);
            if (p_brep_curve_skin->Id() == 0 || r_root_model_part.HasGeometry(p_brep_curve_skin->Id())) {
                const IndexType new_geom_id = NextRegularGeometryId(r_root_model_part);
                p_brep_curve_skin->SetId(new_geom_id);
            }
            r_contact_sub_model_part.AddGeometry(p_brep_curve_skin);
        }
        //------------------------------------------------------------
    };

    // Integrate the signed gap area and create its quadrature-point elements.
    auto create_area_elements = [&](const BrepCurveType& rB0,
                                    const BrepCurveType& rL0,
                                    const BrepCurveType& rL1,
                                    const BrepCurveType& rB1,
                                    const array_1d<double, 3>& rP00,
                                    const array_1d<double, 3>& rP01,
                                    const array_1d<double, 3>& rP10,
                                    const array_1d<double, 3>& rP11,
                                    const double contribution_sign) {
        IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
            (gap_quadrature_order + 1), /*Order*/
            rB0,
            rL0,
            rL1,
            rB1,
            rP00,
            rP01,
            rP10,
            rP11);

        for (auto& integration_point : surface_integration_points) {
            integration_point.SetWeight(integration_point.Weight() * contribution_sign);
        }

        double surface_integration_weight_sum = 0.0;
        for (const auto& r_integration_point : surface_integration_points) {
            surface_integration_weight_sum += r_integration_point.Weight();
        }

        if (std::abs(surface_integration_weight_sum) < 1e-14) return;

        GeometriesArrayType surface_quadrature_point_list;
        IntegrationInfo surface_integration_info = pNurbsSurface->GetDefaultIntegrationInfo();
        pNurbsSurface->CreateQuadraturePointGeometries(surface_quadrature_point_list, rIntegrationParameters.NumberOfShapeFunctionsDerivatives,
                                                       surface_integration_points, surface_integration_info);

        double characteristic_length = CalculateGapElementCharacteristicLength(rP00, rP10, rP01, rP11);
        KRATOS_ERROR_IF(characteristic_length <= 0.0)
            << "Characteristic length for the gap element must be positive." << std::endl;

        IndexType id_element = 1;
        if (mpGapElementsSubModelPart->GetRootModelPart().Elements().size() > 0)
            id_element = mpGapElementsSubModelPart->GetRootModelPart().Elements().back().Id() + 1;

        this->CreateElements(
            surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
            *mpGapElementsSubModelPart, mGapElementName, id_element, PropertiesPointerType(), neighbour_geometries_skin1_skin2, characteristic_length);
    };

    auto p_nurbs_curve_surrogate1_surrogate2 = this->CreateBrepCurve(pSurrogateNode1, pSurrogateNode2, active_range_knot_vector);
    auto p_brep_curve_surrogate1_surrogate2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_surrogate2);
    auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(pSurrogateNode1, p_skin1_brep_point, active_range_knot_vector);
    auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(pSurrogateNode2, p_skin2_brep_point, active_range_knot_vector);
    auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);
    auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);

    const array_1d<double, 3> initial_skin_segment = p_skin_node_2->Coordinates() - p_skin_node_1->Coordinates();
    array_1d<double, 3> initial_normal_skin_segment;
    initial_normal_skin_segment[0] = initial_skin_segment[1];
    initial_normal_skin_segment[1] = -initial_skin_segment[0];
    initial_normal_skin_segment[2] = 0.0;

    const CoordinatesArrayType initial_skin_mid_point = EvaluateGlobalCurvePoint(*p_brep_curve_skin1_skin2_linear, 0.5);
    std::vector<IndexType> subdivision_projection_ids{p_skin_node_1->Id(), p_skin_node_2->Id()};
    // No subdivision projection is needed for a linear gap without subdivisions.
    const auto initial_projection_info = mGapApproximationOrder == 1 && max_internal_divisions == 0
        ? ProjectionResult{p_skin_node_1->Id(), initial_skin_mid_point, true}
        : ComputeSingleProjection<TIsInnerLoop>(
        rSkinSubModelPart,
        initial_skin_mid_point,
        p_skin_node_1->Id(),
        p_skin_node_2->Id(),
        layer_name,
        condition_name,
        initial_normal_skin_segment,
        rIntegrationParameters,
        r_skin_bins_per_span,
        subdivision_projection_ids);

    const array_1d<double, 3> initial_projection_delta = initial_projection_info.ProjectionPoint - initial_skin_mid_point;
    const double initial_distance_to_projection_sq = inner_prod(initial_projection_delta, initial_projection_delta);
    const bool initial_reached_max_divisions = max_internal_divisions == 0;
    const bool initial_degenerate_split = initial_projection_info.ProjectionId == p_skin_node_1->Id() || initial_projection_info.ProjectionId == p_skin_node_2->Id();

    if ((initial_distance_to_projection_sq <= subdivision_threshold_sq || initial_reached_max_divisions || initial_degenerate_split)) {
        std::vector<IndexType> interpolation_nodes_id;
        auto p_nurbs_curve_skin = build_skin_curve(p_skin_node_1->Id(), p_skin_node_2->Id(), interpolation_nodes_id);
        auto p_brep_curve_skin = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin);

        create_skin_conditions(p_brep_curve_skin, interpolation_nodes_id, p_skin_node_1, p_skin_node_2);

        if constexpr (TIsInnerLoop) {
            create_area_elements(
                *p_brep_curve_surrogate1_surrogate2,
                *p_brep_curve_surrogate1_skin1,
                *p_brep_curve_surrogate2_skin2,
                *p_brep_curve_skin,
                pSurrogateNode1->Coordinates(),
                p_skin_node_1->Coordinates(),
                pSurrogateNode2->Coordinates(),
                p_skin_node_2->Coordinates(),
                1.0);
        } else {
            auto p_b0_reversed = ReverseBezierUV_Generic(p_nurbs_curve_surrogate1_surrogate2);
            auto p_b1_reversed = ReverseBezierUV_Generic(p_nurbs_curve_skin);
            auto p_brep_curve_b0 = Kratos::make_shared<BrepCurveType>(p_b0_reversed);
            auto p_brep_curve_b1 = Kratos::make_shared<BrepCurveType>(p_b1_reversed);

            create_area_elements(
                *p_brep_curve_b0,
                *p_brep_curve_surrogate2_skin2,
                *p_brep_curve_surrogate1_skin1,
                *p_brep_curve_b1,
                pSurrogateNode2->Coordinates(),
                p_skin_node_2->Coordinates(),
                pSurrogateNode1->Coordinates(),
                p_skin_node_1->Coordinates(),
                1.0);
        }
    } else {
        if constexpr (TIsInnerLoop) {
            create_area_elements(
                *p_brep_curve_surrogate1_surrogate2,
                *p_brep_curve_surrogate1_skin1,
                *p_brep_curve_surrogate2_skin2,
                *p_brep_curve_skin1_skin2_linear,
                pSurrogateNode1->Coordinates(),
                p_skin_node_1->Coordinates(),
                pSurrogateNode2->Coordinates(),
                p_skin_node_2->Coordinates(),
                1.0);
        } else {
            auto p_b0_reversed = ReverseBezierUV_Generic(p_nurbs_curve_surrogate1_surrogate2);
            auto p_b1_reversed = ReverseBezierUV_Generic(p_nurbs_curve_skin1_skin2_linear);
            auto p_brep_curve_b0 = Kratos::make_shared<BrepCurveType>(p_b0_reversed);
            auto p_brep_curve_b1 = Kratos::make_shared<BrepCurveType>(p_b1_reversed);

            create_area_elements(
                *p_brep_curve_b0,
                *p_brep_curve_surrogate2_skin2,
                *p_brep_curve_surrogate1_skin1,
                *p_brep_curve_b1,
                pSurrogateNode2->Coordinates(),
                p_skin_node_2->Coordinates(),
                pSurrogateNode1->Coordinates(),
                p_skin_node_1->Coordinates(),
                1.0);
        }
    }

    // Recursively add signed correction patches until the subdivision criterion is met.
    auto process_piece = [&](auto&& self,
                             const Node::Pointer& p_left_skin,
                             const Node::Pointer& p_right_skin,
                             const BrepCurveType::Pointer& p_brep_curve_linear,
                             std::size_t division_level,
                             std::size_t positive_count,
                             std::size_t negative_count) -> void
    {
        KRATOS_INFO_IF("SnakeGapSbmProcess", mEchoLevel > 3)
            << "process_piece: div=" << division_level
            << " left_id=" << p_left_skin->Id()
            << " right_id=" << p_right_skin->Id()
            << " pos=" << positive_count
            << " neg=" << negative_count << std::endl;

        const array_1d<double, 3> skin_segment = p_right_skin->Coordinates() - p_left_skin->Coordinates();
        array_1d<double, 3> normal_skin_segment;
        normal_skin_segment[0] = skin_segment[1];
        normal_skin_segment[1] = -skin_segment[0];
        normal_skin_segment[2] = 0.0;

        const CoordinatesArrayType skin_mid_point = EvaluateGlobalCurvePoint(*p_brep_curve_linear, 0.5);

        const auto projection_info = mGapApproximationOrder == 1 && max_internal_divisions == 0
            ? ProjectionResult{p_left_skin->Id(), skin_mid_point, true}
            : ComputeSingleProjection<TIsInnerLoop>(
            rSkinSubModelPart,
            skin_mid_point,
            p_left_skin->Id(),
            p_right_skin->Id(),
            layer_name,
            condition_name,
            normal_skin_segment,
            rIntegrationParameters,
            r_skin_bins_per_span,
            subdivision_projection_ids);

        const array_1d<double, 3> projection_delta = projection_info.ProjectionPoint - skin_mid_point;
        const double distance_to_projection_sq = inner_prod(projection_delta, projection_delta);

        const bool reached_max_divisions = division_level >= max_internal_divisions;
        const bool degenerate_split = std::find(
            subdivision_projection_ids.begin(), subdivision_projection_ids.end(), projection_info.ProjectionId) !=
            subdivision_projection_ids.end();

        const double dot_product = inner_prod(projection_delta, normal_skin_segment);
        const double new_piece_contribution_sign = (dot_product >= 0.0) ? 1.0 : -1.0;

        KRATOS_INFO_IF("SnakeGapSbmProcess", mEchoLevel > 3)
            << "process_piece: proj_id=" << projection_info.ProjectionId
            << " dist2=" << distance_to_projection_sq
            << " dot=" << dot_product
            << " sign=" << new_piece_contribution_sign << std::endl;

        if ((distance_to_projection_sq <= subdivision_threshold_sq || reached_max_divisions || degenerate_split)) {
            KRATOS_INFO_IF("SnakeGapSbmProcess", mEchoLevel > 3)
                << "process_piece: finalize_piece"
                << " reached_max=" << reached_max_divisions
                << " degenerate=" << degenerate_split
                << std::endl;

            std::vector<IndexType> interpolation_nodes_id;
            auto p_nurbs_curve_skin = build_skin_curve(p_left_skin->Id(), p_right_skin->Id(), interpolation_nodes_id);
            auto p_brep_curve_skin = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin);

            create_skin_conditions(p_brep_curve_skin, interpolation_nodes_id, p_left_skin, p_right_skin);

            const auto p_new_skin_node = rSkinSubModelPart.pGetNode(projection_info.ProjectionId);

            auto p_nurbs_curve_l0 = this->CreateBrepCurve(p_left_skin, p_new_skin_node, active_range_knot_vector);
            auto p_nurbs_curve_l1 = this->CreateBrepCurve(p_right_skin, p_new_skin_node, active_range_knot_vector);

            auto p_brep_curve_l0 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_l0);
            auto p_brep_curve_l1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_l1);

            auto p_brep_curve_b0 = p_brep_curve_linear;
            auto p_brep_curve_b1 = p_brep_curve_skin;

            if constexpr (!TIsInnerLoop) {
                auto p_b0_background = std::dynamic_pointer_cast<NurbsCurveGeometryType>(
                    p_brep_curve_b0->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
                auto p_b1_background = std::dynamic_pointer_cast<NurbsCurveGeometryType>(
                    p_brep_curve_b1->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
                auto p_b0_reversed = ReverseBezierUV_Generic(p_b0_background);
                auto p_b1_reversed = ReverseBezierUV_Generic(p_b1_background);
                p_brep_curve_b0 = Kratos::make_shared<BrepCurveType>(p_b0_reversed);
                p_brep_curve_b1 = Kratos::make_shared<BrepCurveType>(p_b1_reversed);
            }

            const auto& skin_left_coords = p_left_skin->Coordinates();
            const auto& skin_right_coords = p_right_skin->Coordinates();
            const auto& skin_new_coords = p_new_skin_node->Coordinates();
            create_area_elements(
                *p_brep_curve_b0,
                *p_brep_curve_l0,
                *p_brep_curve_l1,
                *p_brep_curve_b1,
                skin_left_coords,
                skin_new_coords,
                skin_right_coords,
                skin_new_coords,
                new_piece_contribution_sign);
            return;
        } else {
            KRATOS_INFO_IF("SnakeGapSbmProcess", mEchoLevel > 3)
                << "process_piece: subdivide"
                << " reached_max=" << reached_max_divisions
                << std::endl;

            const auto p_new_skin_node = rSkinSubModelPart.pGetNode(projection_info.ProjectionId);

            subdivision_projection_ids.push_back(projection_info.ProjectionId);
            std::size_t next_positive_count = positive_count + (new_piece_contribution_sign > 0.0 ? 1 : 0);
            std::size_t next_negative_count = negative_count + (new_piece_contribution_sign < 0.0 ? 1 : 0);

            auto p_nurbs_curve_l0 = this->CreateBrepCurve(p_left_skin, p_new_skin_node, active_range_knot_vector);
            auto p_nurbs_curve_l1 = this->CreateBrepCurve(p_right_skin, p_new_skin_node, active_range_knot_vector);
            auto p_nurbs_curve_b1 = this->CreateBrepCurve(p_new_skin_node, p_new_skin_node, active_range_knot_vector);

            auto p_brep_curve_l0 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_l0);
            auto p_brep_curve_l1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_l1);
            auto p_brep_curve_b1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_b1);

            auto p_brep_curve_b0 = p_brep_curve_linear;
            if constexpr (!TIsInnerLoop) {
                auto p_b0_background = std::dynamic_pointer_cast<NurbsCurveGeometryType>(
                    p_brep_curve_b0->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
                auto p_b0_reversed = ReverseBezierUV_Generic(p_b0_background);
                p_brep_curve_b0 = Kratos::make_shared<BrepCurveType>(p_b0_reversed);
            }

            const auto& skin_left_coords = p_left_skin->Coordinates();
            const auto& skin_right_coords = p_right_skin->Coordinates();
            const auto& skin_new_coords = p_new_skin_node->Coordinates();

            create_area_elements(
                *p_brep_curve_b0,
                *p_brep_curve_l0,
                *p_brep_curve_l1,
                *p_brep_curve_b1,
                skin_left_coords,
                skin_new_coords,
                skin_right_coords,
                skin_new_coords,
                new_piece_contribution_sign);

            KRATOS_INFO_IF("SnakeGapSbmProcess", mEchoLevel > 3)
                << "process_piece: recurse"
                << " next_div=" << (division_level + 1)
                << " next_pos=" << next_positive_count
                << " next_neg=" << next_negative_count << std::endl;

            auto p_nurbs_left_linear = this->CreateBrepCurve(p_left_skin, p_new_skin_node, active_range_knot_vector);
            auto p_nurbs_right_linear = this->CreateBrepCurve(p_new_skin_node, p_right_skin, active_range_knot_vector);
            auto p_brep_left_linear = Kratos::make_shared<BrepCurveType>(p_nurbs_left_linear);
            auto p_brep_right_linear = Kratos::make_shared<BrepCurveType>(p_nurbs_right_linear);

            self(self, p_left_skin, p_new_skin_node, p_brep_left_linear, division_level + 1, next_positive_count, next_negative_count);
            self(self, p_new_skin_node, p_right_skin, p_brep_right_linear, division_level + 1, next_positive_count, next_negative_count);
        }
    };

    if (!(initial_distance_to_projection_sq <= subdivision_threshold_sq || initial_reached_max_divisions || initial_degenerate_split)
       ) {
            process_piece(process_piece, p_skin_node_1, p_skin_node_2, p_brep_curve_skin1_skin2_linear, 0, 0, 0);
    }

}

// Instantiate registered conditions and attach surrogate reference data.
void SnakeGapSbmProcess::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rConditionName,
    std::size_t& rIdCounter,
    PropertiesPointerType pProperties,
    const Vector KnotSpanSizes,
    const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries,
    const double CharacteristicLength) const
{
    KRATOS_ERROR_IF(CharacteristicLength <= 0.0)
        << "Characteristic length for gap conditions must be positive." << std::endl;

    const std::string and_token = "_AND_";
    std::vector<std::string> condition_names;
    const bool has_and_token = rConditionName.find(and_token) != std::string::npos;
    if (has_and_token) {
        size_t start_pos = 0;
        while (true) {
            const size_t end_pos = rConditionName.find(and_token, start_pos);
            const std::string name = rConditionName.substr(start_pos, end_pos - start_pos);
            KRATOS_ERROR_IF(name.empty()) << "Empty condition name token in " << rConditionName << std::endl;
            condition_names.push_back(name);
            if (end_pos == std::string::npos) {
                break;
            }
            start_pos = end_pos + and_token.size();
        }
    }

    if (has_and_token) {
        for (const auto& condition_name : condition_names) {
            CreateConditions(
                rGeometriesBegin,
                rGeometriesEnd,
                rModelPart,
                condition_name,
                rIdCounter,
                pProperties,
                KnotSpanSizes,
                pSurrogateReferenceGeometries,
                CharacteristicLength);
        }
        return;
    }

    KRATOS_ERROR_IF_NOT(KratosComponents<Condition>::Has(rConditionName))
        << "Condition \"" << rConditionName << "\" is not registered."
        << std::endl;
    const Condition& reference_condition =
        KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 3)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    array_1d<double, 3> characteristic_length_vector = ZeroVector(3);
    characteristic_length_vector[0] = CharacteristicLength;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
        auto p_condition = reference_condition.Create(rIdCounter, *it, pProperties);
        new_condition_list.push_back(p_condition);

        if ((*it)->Has(BREP_ID)) {
            p_condition->SetValue(BREP_ID, (*it)->GetValue(BREP_ID));
        }

        // Set knot span sizes to the condition
        p_condition->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);

        p_condition->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);
        p_condition->SetValue(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_vector);
        if (rConditionName == "GapSbmContactCondition" ||
            rConditionName == "GapSbmALMContactCondition") {
            p_condition->SetValue(IDENTIFIER, "INACTIVE");
        }

        rIdCounter++;
    }

    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}

// Instantiate gap elements and attach surrogate reference data.
void SnakeGapSbmProcess::CreateElements(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rElementName,
    std::size_t& rIdCounter,
    PropertiesPointerType pProperties,
    const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries,
    const double CharacteristicLength) const
{
    KRATOS_ERROR_IF(!KratosComponents<Element>::Has(rElementName))
        << rElementName << " not registered." << std::endl;

    const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

    ElementsContainerType new_element_list;

    KRATOS_INFO_IF("CreateElements", mEchoLevel > 3)
        << "Creating elements of type " << rElementName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    std::size_t num_elements = std::distance(rGeometriesBegin, rGeometriesEnd);
    new_element_list.reserve(num_elements);

    array_1d<double, 3> characteristic_length_vector = ZeroVector(3);
    characteristic_length_vector[0] = CharacteristicLength;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
        auto p_element = rReferenceElement.Create(rIdCounter, *it, pProperties);
        new_element_list.push_back(p_element);

        p_element->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);
        p_element->SetValue(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_vector);
        rIdCounter++;
    }

    rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
}

// Search the skin bins for the closest node satisfying the layer and direction criteria.
template <bool TIsInnerLoop>
IndexType SnakeGapSbmProcess::FindClosestNodeInLayerWithDirection(
    const array_1d<double, 3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
    const Vector& rDirection,
    const std::vector<IndexType>* pConsideredNodeIds) const
{
    KRATOS_ERROR_IF(rSkinBinsPerSpan.NumberOfSpansX == 0 ||
                    rSkinBinsPerSpan.NumberOfSpansY == 0)
        << "::[SnakeGapSbmProcess]::FindClosestNodeInLayerWithDirection: skin span bins are empty." << std::endl;

    Vector direction = rDirection;
    constexpr double dir_tol = 1.0e-12;
    double direction_norm = norm_2(direction);
    if (direction_norm < dir_tol) {
        direction = ZeroVector(3);
        direction[0] = 1.0;
        direction_norm = 1.0;
    }
    direction /= direction_norm;

    if constexpr (TIsInnerLoop) {
        direction *= -1.0;
    }

    const double span_size_x = (rSkinBinsPerSpan.SpanSizeX > 0.0) ? rSkinBinsPerSpan.SpanSizeX
                           : (rKnotSpanSizes.size() > 0 ? rKnotSpanSizes[0] : 1.0);
    const double span_size_y = (rSkinBinsPerSpan.SpanSizeY > 0.0) ? rSkinBinsPerSpan.SpanSizeY
                           : (rKnotSpanSizes.size() > 1 ? rKnotSpanSizes[1] : span_size_x);
    double reference_span_size = std::max(span_size_x, span_size_y);
    if (rKnotSpanSizes.size() > 0) {
        reference_span_size = 0.0;
        for (std::size_t i = 0; i < std::min<std::size_t>(2, rKnotSpanSizes.size()); ++i) {
            reference_span_size = std::max(reference_span_size, rKnotSpanSizes[i]);
        }
    }
    if (reference_span_size <= dir_tol) {
        reference_span_size = 1.0;
    }

    const double min_u = rSkinBinsPerSpan.MinU;
    const double max_u = rSkinBinsPerSpan.MaxU;
    const double min_v = rSkinBinsPerSpan.MinV;
    const double max_v = rSkinBinsPerSpan.MaxV;

    auto compute_index = [](double coordinate, double min_value, double max_value, double span_size, std::size_t span_count){
        double clamped = coordinate;
        if (clamped < min_value) {
            clamped = min_value;
        } else if (clamped > max_value) {
            clamped = max_value;
        }

        if (span_count == 0) {
            return std::size_t(0);
        }

        double relative = (clamped - min_value) / span_size;
        std::size_t index = static_cast<std::size_t>(std::floor(relative + 1.0e-12));
        if (index >= span_count) {
            index = span_count - 1;
        }
        return index;
    };

    const std::size_t number_spans_x = rSkinBinsPerSpan.NumberOfSpansX;
    const std::size_t number_spans_y = rSkinBinsPerSpan.NumberOfSpansY;
    const std::size_t base_ix = compute_index(rStartPoint[0], min_u, max_u, span_size_x, number_spans_x);
    const std::size_t base_iy = compute_index(rStartPoint[1], min_v, max_v, span_size_y, number_spans_y);
    const int max_search_level = static_cast<int>(std::max(number_spans_x, number_spans_y));

    ConditionPointerContainerType candidate_conditions;
    candidate_conditions.reserve(2000);
    const bool filter_by_layer = !rLayer.empty();
    auto node_has_requested_layer = [&](Node::Pointer p_candidate_node) {
        if (!filter_by_layer) {
            return true;
        }
        const auto& candidate_layers = p_candidate_node->GetValue(CONNECTED_LAYERS);
        return std::find(candidate_layers.begin(), candidate_layers.end(), rLayer) != candidate_layers.end();
    };
    const IndexType invalid_node_id = std::numeric_limits<IndexType>::max();
    auto already_considered = [&](const IndexType NodeId) {
        return pConsideredNodeIds != nullptr &&
            std::find(pConsideredNodeIds->begin(), pConsideredNodeIds->end(), NodeId) != pConsideredNodeIds->end();
    };
    auto search_with_direction = [&](const array_1d<double, 3>& rSearchDirection) {
        IndexType best_node_id = std::numeric_limits<IndexType>::max();
        double best_intersection_distance = std::numeric_limits<double>::max();
        double best_intersection_node_distance = std::numeric_limits<double>::max();

        auto update_candidate = [&](Node::Pointer p_candidate_node,
                                    const CoordinatesArrayType& rIntersectionPoint,
                                    const bool ConditionMatchesLayer){
            if (!ConditionMatchesLayer && !node_has_requested_layer(p_candidate_node)) {
                return;
            }
            array_1d<double, 3> diff_intersection = rIntersectionPoint;
            diff_intersection -= rStartPoint;
            const double distance_intersection = norm_2(diff_intersection);

            const double distance_tolerance = 1.0e-12 * std::max(1.0, distance_intersection);
            if (distance_intersection < best_intersection_distance - distance_tolerance) {
                best_intersection_distance = distance_intersection;
                best_intersection_node_distance = std::numeric_limits<double>::max();
                best_node_id = invalid_node_id;
            } else if (distance_intersection > best_intersection_distance + distance_tolerance) {
                return;
            }

            const double distance_intersection_node = norm_2(p_candidate_node->Coordinates() - rIntersectionPoint);
            const bool candidate_considered = already_considered(p_candidate_node->Id());
            const bool best_considered = already_considered(best_node_id);
            if (best_node_id == invalid_node_id ||
                (best_considered && !candidate_considered) ||
                (best_considered == candidate_considered &&
                 distance_intersection_node < best_intersection_node_distance - distance_tolerance)) {
                best_intersection_distance = std::min(best_intersection_distance, distance_intersection);
                best_intersection_node_distance = distance_intersection_node;
                best_node_id = p_candidate_node->Id();
            }
        };

        bool found_intersection = false;
        for (int level = 0; level < max_search_level && !found_intersection; ++level) {
            const std::size_t extension = static_cast<std::size_t>(level+1);
            const double current_length = reference_span_size * static_cast<double>(extension);

            array_1d<double, 3> segment_start = rStartPoint - rSearchDirection * reference_span_size*0.5;
            array_1d<double, 3> segment_end = rStartPoint;
            segment_end += rSearchDirection * current_length;

            const std::size_t min_ix = (base_ix > extension) ? base_ix - extension : 0;
            const std::size_t max_ix = std::min<std::size_t>(base_ix + extension, number_spans_x > 0 ? number_spans_x - 1 : 0);
            const std::size_t min_iy = (base_iy > extension) ? base_iy - extension : 0;
            const std::size_t max_iy = std::min<std::size_t>(base_iy + extension, number_spans_y > 0 ? number_spans_y - 1 : 0);

            candidate_conditions.clear();

            Node::Pointer start_node = Kratos::make_intrusive<Node>(0, segment_start[0], segment_start[1], segment_start[2]);
            Node::Pointer end_node   = Kratos::make_intrusive<Node>(0, segment_end[0],   segment_end[1],   segment_end[2]);

            Condition::NodesArrayType ray_nodes;
            ray_nodes.reserve(2);
            ray_nodes.push_back(start_node);
            ray_nodes.push_back(end_node);
            ConditionPointerType p_ray_condition = Kratos::make_intrusive<Condition>(0, ray_nodes);

            for (std::size_t ix = min_ix; ix <= max_ix; ++ix) {
                if (ix >= rSkinBinsPerSpan.NumberOfSpansX) continue;
                for (std::size_t iy = min_iy; iy <= max_iy; ++iy) {
                    if (iy >= rSkinBinsPerSpan.NumberOfSpansY) continue;
                    const std::size_t k = FindSpanNnzIndex(rSkinBinsPerSpan, ix, iy, 0);
                    if (k == static_cast<std::size_t>(-1)) continue;
                    if (k >= rSkinBinsPerSpan.CellBinsByNnz.size()) continue;
                    const auto& r_cell_bins = rSkinBinsPerSpan.CellBinsByNnz[k];
                    if (!r_cell_bins.HasConditionBins) continue;
                    auto& r_bins = const_cast<BinsObjectDynamic<ConditionConfigure>&>(r_cell_bins.ConditionBins);
                    r_bins.SearchObjects(p_ray_condition, candidate_conditions);
                }
            }

            if (candidate_conditions.empty()) {
                continue;
            }

            for (const auto& p_condition : candidate_conditions) {
                const auto p_geometry = p_condition->pGetGeometry();
                if (p_geometry->size() < 2) {
                    continue;
                }

                const auto& node_a = p_geometry->pGetPoint(0);
                const auto& node_b = p_geometry->pGetPoint(1);
                const bool condition_matches_layer =
                    !filter_by_layer || (p_condition->Has(LAYER_NAME) && p_condition->GetValue(LAYER_NAME) == rLayer);
                if (!condition_matches_layer && !node_has_requested_layer(node_a) && !node_has_requested_layer(node_b)) {
                    continue;
                }

                CoordinatesArrayType intersection_point;
                const bool intersects = SegmentsIntersect(start_node, end_node, node_a, node_b, intersection_point);

                if (intersects) {
                    update_candidate(node_a, intersection_point, condition_matches_layer);
                    update_candidate(node_b, intersection_point, condition_matches_layer);
                    found_intersection = true;
                }
            }

            if (best_node_id != invalid_node_id) {
                return std::make_pair(best_node_id, best_intersection_distance);
            }
        }

        return std::make_pair(best_node_id, best_intersection_distance);
    };

    const auto forward_candidate = search_with_direction(direction);
    const array_1d<double, 3> opposite_direction = -direction;
    const auto opposite_candidate = search_with_direction(opposite_direction);
    auto best_candidate = forward_candidate;
    const double distance_tolerance = 1.0e-12 * std::max(1.0, opposite_candidate.second);
    if (opposite_candidate.first != invalid_node_id &&
        (forward_candidate.first == invalid_node_id ||
         (already_considered(forward_candidate.first) && !already_considered(opposite_candidate.first)) ||
         (already_considered(forward_candidate.first) == already_considered(opposite_candidate.first) &&
          opposite_candidate.second < forward_candidate.second - distance_tolerance))) {
        best_candidate = opposite_candidate;
    }

    if (best_candidate.first != invalid_node_id) {
        if (rLayer == "COUPLING_SIDE" && best_candidate.second > reference_span_size / 2.0) {
            KRATOS_WARNING("SnakeGapSbmProcess")
                << "::[SnakeGapSbmProcess]::FindClosestNodeInLayerWithDirection: projection on COUPLING_SIDE farther than reference_span_size/2."
                << " distance=" << best_candidate.second
                << ", threshold=" << (reference_span_size / 2.0)
                << ", start_point=" << rStartPoint
                << ", selected_node=" << best_candidate.first
                << std::endl;
        }
        return best_candidate.first;
    }

    KRATOS_ERROR << "::[SnakeGapSbmProcess]::FindClosestNodeInLayerWithDirection: no node found for layer "
                 << (rLayer.empty() ? std::string("<any>") : rLayer) << " starting from point " << rStartPoint
                 << ". Skin sub model part '" << rSkinSubModelPart.Name() << "' has "
                 << rSkinSubModelPart.NumberOfNodes() << " nodes and "
                 << rSkinSubModelPart.NumberOfConditions() << " conditions."
                 << " Search parameters: base_ix=" << base_ix << ", base_iy=" << base_iy
                 << ", span_size_x=" << span_size_x << ", span_size_y=" << span_size_y
                 << ", reference_span_size=" << reference_span_size
                 << ", direction=" << direction << std::endl;

    return 0;
}

// Assign skin projections to surrogate nodes and check for crossing segments.
template <bool TIsInnerLoop>
void SnakeGapSbmProcess::SetSurrogateToSkinProjections2D(
    const ModelPart& rSurrogateSubModelPart,
    const ModelPart& rSkinSubModelPart,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan)
{
    const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();
    const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);

    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: KNOT_SPAN_SIZES must contain at least two entries." << std::endl;

    KRATOS_ERROR_IF(rSkinBinsPerSpan.NumberOfSpansX == 0 ||
                    rSkinBinsPerSpan.NumberOfSpansY == 0)
        << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: skin span bins are empty." << std::endl;

    const double span_size_x = (rSkinBinsPerSpan.SpanSizeX > 0.0) ? rSkinBinsPerSpan.SpanSizeX : knot_span_sizes[0];
    const double span_size_y = (rSkinBinsPerSpan.SpanSizeY > 0.0) ? rSkinBinsPerSpan.SpanSizeY : knot_span_sizes[1];
    KRATOS_ERROR_IF(span_size_x <= 0.0 || span_size_y <= 0.0)
        << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: invalid knot span sizes."
        << " span_size_x=" << span_size_x << ", span_size_y=" << span_size_y << std::endl;

    const double min_u = rSkinBinsPerSpan.MinU;
    const double max_u = rSkinBinsPerSpan.MaxU;
    const double min_v = rSkinBinsPerSpan.MinV;
    const double max_v = rSkinBinsPerSpan.MaxV;
    const std::size_t span_count_x = rSkinBinsPerSpan.NumberOfSpansX;
    const std::size_t span_count_y = rSkinBinsPerSpan.NumberOfSpansY;
    double min_knot_span_size = std::min(span_size_x, span_size_y);

    const double span_probe_length = 1.0e-2 * min_knot_span_size;

    auto clamp_coordinate = [](double coordinate, double min_value, double max_value) {
        if (coordinate < min_value) {
            return min_value;
        }
        if (coordinate > max_value) {
            return max_value;
        }
        return coordinate;
    };

    auto compute_span_index = [&](double coordinate,
                                  double min_value,
                                  double max_value,
                                  double span_size,
                                  std::size_t span_count) {
        if (span_count == 0) {
            return std::size_t(0);
        }

        const double clamped = clamp_coordinate(coordinate, min_value, max_value);
        const double relative = (clamped - min_value) / span_size;
        std::size_t index = static_cast<std::size_t>(std::floor(relative + 1.0e-12));
        if (index >= span_count) {
            index = span_count - 1;
        }
        return index;
    };

    struct SpanId
    {
        std::size_t X = 0;
        std::size_t Y = 0;

        bool operator==(const SpanId& rOther) const
        {
            return X == rOther.X && Y == rOther.Y;
        }
    };

    auto add_span_id = [](std::vector<SpanId>& rSpanIds,
                          const std::size_t SpanIndexX,
                          const std::size_t SpanIndexY) {
        const SpanId span_id{SpanIndexX, SpanIndexY};
        if (std::find(rSpanIds.begin(), rSpanIds.end(), span_id) == rSpanIds.end()) {
            rSpanIds.push_back(span_id);
        }
    };

    auto add_span_id_at_point = [&](std::vector<SpanId>& rSpanIds,
                                    const array_1d<double, 3>& rPoint) {
        add_span_id(
            rSpanIds,
            compute_span_index(rPoint[0], min_u, max_u, span_size_x, span_count_x),
            compute_span_index(rPoint[1], min_v, max_v, span_size_y, span_count_y));
    };

    // Select the closest admissible skin node using directional span probes.
    auto select_candidate = [&](const Node::Pointer& pSurrogateNode) {
        KRATOS_ERROR_IF_NOT(pSurrogateNode->Has(NORMAL))
            << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: surrogate node "
            << pSurrogateNode->Id() << " has no NORMAL value." << std::endl;

        array_1d<double, 3> normal_direction = pSurrogateNode->GetValue(NORMAL);
        const double normal_norm = norm_2(normal_direction);
        KRATOS_ERROR_IF(normal_norm <= 1.0e-16)
            << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: surrogate node "
            << pSurrogateNode->Id() << " has a zero NORMAL value." << std::endl;
        normal_direction /= normal_norm;

        array_1d<double, 3> reference_axis = ZeroVector(3);
        reference_axis[2] = 1.0;
        if (std::abs(inner_prod(normal_direction, reference_axis)) > 0.9) {
            reference_axis = ZeroVector(3);
            reference_axis[0] = 1.0;
        }
        array_1d<double, 3> tangent_direction = MathUtils<double>::CrossProduct(normal_direction, reference_axis);
        const double tangent_norm = norm_2(tangent_direction);
        KRATOS_ERROR_IF(tangent_norm <= 1.0e-16)
            << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: surrogate node "
            << pSurrogateNode->Id() << " has a degenerate tangent from NORMAL " << normal_direction << std::endl;
        tangent_direction /= tangent_norm;

        std::vector<SpanId> span_ids;
        span_ids.reserve(3);
        const array_1d<double, 3> normal_probe = pSurrogateNode->Coordinates() + span_probe_length * normal_direction;
        add_span_id_at_point(span_ids, normal_probe);
        add_span_id_at_point(span_ids, normal_probe + span_probe_length * tangent_direction);
        add_span_id_at_point(span_ids, normal_probe - span_probe_length * tangent_direction);

        IndexType best_any_id = std::numeric_limits<IndexType>::max();
        double best_any_distance = std::numeric_limits<double>::max();
        auto consider_candidate = [&](const Node& r_candidate_node) {
            const array_1d<double, 3> diff = r_candidate_node.Coordinates() - pSurrogateNode->Coordinates();
            const double distance_sq = inner_prod(diff, diff);

            if (distance_sq < best_any_distance) {
                best_any_distance = distance_sq;
                best_any_id = r_candidate_node.Id();
            }
        };

        std::unordered_set<IndexType> visited_candidate_ids;
        visited_candidate_ids.reserve(16);

        auto search_span = [&](const SpanId& r_span_id) {
            const std::size_t ix = r_span_id.X;
            const std::size_t iy = r_span_id.Y;
            if (ix >= rSkinBinsPerSpan.NumberOfSpansX || iy >= rSkinBinsPerSpan.NumberOfSpansY) {
                return;
            }

            const std::size_t k = FindSpanNnzIndex(rSkinBinsPerSpan, ix, iy, 0);
            if (k == static_cast<std::size_t>(-1)) {
                return;
            }
            if (k >= rSkinBinsPerSpan.CellBinsByNnz.size()) {
                return;
            }

            const auto& r_cell_bins = rSkinBinsPerSpan.CellBinsByNnz[k];
            if (!r_cell_bins.HasNodeBins) {
                return;
            }

            for (const auto& p_candidate_node : r_cell_bins.Nodes) {
                if (!p_candidate_node) {
                    continue;
                }
                if (!visited_candidate_ids.insert(p_candidate_node->Id()).second) {
                    continue;
                }
                consider_candidate(*p_candidate_node);
            }
        };

        // Search the normal span first to preserve its priority.
        if (!span_ids.empty()) {
            search_span(span_ids.front());
        }

        std::vector<SpanId> opposite_normal_span_ids;
        opposite_normal_span_ids.reserve(1);
        add_span_id_at_point(
            opposite_normal_span_ids,
            pSurrogateNode->Coordinates() - span_probe_length * normal_direction);

        const IndexType normal_any_id = best_any_id;
        const double normal_any_distance = best_any_distance;
        best_any_id = std::numeric_limits<IndexType>::max();
        best_any_distance = std::numeric_limits<double>::max();

        // Compare both tangent spans and the span opposite to the normal.
        for (std::size_t span_index = 1; span_index < span_ids.size(); ++span_index) {
            opposite_normal_span_ids.push_back(span_ids[span_index]);
        }
        for (const auto& r_span_id : opposite_normal_span_ids) {
            search_span(r_span_id);
        }

        const double alternative_distance_limit = (1.0 - 1.0e-5) * min_knot_span_size;
        const double alternative_distance_limit_sq = alternative_distance_limit * alternative_distance_limit;
        if (normal_any_id != std::numeric_limits<IndexType>::max() &&
            (best_any_distance >= normal_any_distance || best_any_distance >= alternative_distance_limit_sq)) {
            best_any_id = normal_any_id;
            best_any_distance = normal_any_distance;
        }
        // At a surrogate corner, the closest skin nodes can lie in an adjacent
        // span that the short normal/tangent probes do not reach.
        if (best_any_id == std::numeric_limits<IndexType>::max()) {
            const std::size_t base_ix = compute_span_index(
                pSurrogateNode->X(), min_u, max_u, span_size_x, span_count_x);
            const std::size_t base_iy = compute_span_index(
                pSurrogateNode->Y(), min_v, max_v, span_size_y, span_count_y);
            for (std::size_t ix = base_ix > 0 ? base_ix - 1 : 0;
                 ix <= std::min(base_ix + 1, span_count_x - 1); ++ix) {
                for (std::size_t iy = base_iy > 0 ? base_iy - 1 : 0;
                     iy <= std::min(base_iy + 1, span_count_y - 1); ++iy) {
                    const std::size_t k = FindSpanNnzIndex(rSkinBinsPerSpan, ix, iy, 0);
                    if (k == static_cast<std::size_t>(-1) || k >= rSkinBinsPerSpan.CellBinsByNnz.size()) {
                        continue;
                    }
                    for (const auto& p_candidate_node : rSkinBinsPerSpan.CellBinsByNnz[k].Nodes) {
                        if (p_candidate_node) {
                            consider_candidate(*p_candidate_node);
                        }
                    }
                }
            }
        }

        if (best_any_id == std::numeric_limits<IndexType>::max() &&
            rSkinSubModelPart.HasSubModelPart("interface_vertices")) {
            const auto& r_interface_vertices = rSkinSubModelPart.GetSubModelPart("interface_vertices");
            for (const auto& r_interface_node : r_interface_vertices.Nodes()) {
                consider_candidate(r_interface_node);
            }
        }

        return best_any_id;
    };

    // Assign each surrogate node once, even when shared by several conditions.
    for (auto& r_surrogate_condition : rSurrogateSubModelPart.Conditions()) {
        for (std::size_t i = 0; i < 2; ++i) {
            auto p_surrogate_node = r_surrogate_condition.pGetGeometry()->pGetPoint(i);
            if (p_surrogate_node->Has(PROJECTION_NODE_ID)) {
                continue;
            }
            const IndexType skin_node_id = select_candidate(p_surrogate_node);
            KRATOS_ERROR_IF(skin_node_id == std::numeric_limits<IndexType>::max())
                << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections2D: no skin node candidate found for surrogate node "
                << p_surrogate_node->Id() << " at " << p_surrogate_node->Coordinates() << std::endl;

            const auto& r_skin_node = rSkinSubModelPart.GetNode(skin_node_id);
            p_surrogate_node->SetValue(PROJECTION_NODE_ID, skin_node_id);
            p_surrogate_node->SetValue(CONNECTED_LAYERS, r_skin_node.GetValue(CONNECTED_LAYERS));
            p_surrogate_node->SetValue(CONNECTED_CONDITIONS, r_skin_node.GetValue(CONNECTED_CONDITIONS));
        }
    }

    // Detect crossing projections and report an invalid gap reconstruction.
    std::size_t count_intersections = 0;
    for (auto& r_surrogate_condition : rSurrogateSubModelPart.Conditions()) {
        Node::Pointer p_surrogate_node_1 = r_surrogate_condition.pGetGeometry()->pGetPoint(0);
        Node::Pointer p_surrogate_node_2 = r_surrogate_condition.pGetGeometry()->pGetPoint(1);

        IndexType projection_node_id_1 = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
        IndexType projection_node_id_2 = p_surrogate_node_2->GetValue(PROJECTION_NODE_ID);

        Node::Pointer p_projection_node_1 = rSkinSubModelPart.pGetNode(projection_node_id_1);
        Node::Pointer p_projection_node_2 = rSkinSubModelPart.pGetNode(projection_node_id_2);

        if (projection_node_id_1 == projection_node_id_2) continue;

        if (SegmentsIntersect(p_surrogate_node_1, p_projection_node_1,
            p_surrogate_node_2, p_projection_node_2))
        {
            p_surrogate_node_1->SetValue(PROJECTION_NODE_ID, projection_node_id_2);
            p_surrogate_node_2->SetValue(PROJECTION_NODE_ID, projection_node_id_1);

            count_intersections += 1;
        }
    }
    KRATOS_ERROR_IF(count_intersections > 0) <<
        "::[SnakeGapSbmProcess]:: Intersections found between projections in the gap reconstruction. Please check the data." << std::endl;

}

// Find the skin segment intersected by the directional projection search.
template <bool TIsInnerLoop>
std::pair<SnakeGapSbmProcess::IndexType, SnakeGapSbmProcess::IndexType> SnakeGapSbmProcess::FindClosestPairInLayerWithNormalDirection(
    const array_1d<double, 3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
    const Vector& rDirection,
    std::vector<IndexType>* pClosestIntersectionCandidateIds,
    const std::vector<IndexType>* pConsideredNodeIds,
    const bool ReverseNormalFilter,
    const bool IgnoreNormalFilter) const
{

    KRATOS_ERROR_IF(rSkinBinsPerSpan.NumberOfSpansX == 0 ||
                    rSkinBinsPerSpan.NumberOfSpansY == 0)
        << "::[SnakeGapSbmProcess]::FindClosestPairInLayerWithNormalDirection: skin span bins are empty." << std::endl;

    Vector direction = rDirection;
    constexpr double dir_tol = 1.0e-12;
    double direction_norm = norm_2(direction);
    if (direction_norm < dir_tol) {
        direction = ZeroVector(3);
        direction[0] = 1.0;
        direction_norm = 1.0;
    }
    direction /= direction_norm;

    if constexpr (TIsInnerLoop) {
        direction *= -1.0;
    }

    const double span_size_x = (rSkinBinsPerSpan.SpanSizeX > 0.0) ? rSkinBinsPerSpan.SpanSizeX
                           : (rKnotSpanSizes.size() > 0 ? rKnotSpanSizes[0] : 1.0);
    const double span_size_y = (rSkinBinsPerSpan.SpanSizeY > 0.0) ? rSkinBinsPerSpan.SpanSizeY
                           : (rKnotSpanSizes.size() > 1 ? rKnotSpanSizes[1] : span_size_x);
    double reference_span_size = std::max(span_size_x, span_size_y);
    if (rKnotSpanSizes.size() > 0) {
        reference_span_size = 0.0;
        for (std::size_t i = 0; i < rKnotSpanSizes.size(); ++i) {
            reference_span_size = std::max(reference_span_size, rKnotSpanSizes[i]);
        }
    }
    if (reference_span_size <= dir_tol) {
        reference_span_size = 1.0;
    }

    const double min_u = rSkinBinsPerSpan.MinU;
    const double max_u = rSkinBinsPerSpan.MaxU;
    const double min_v = rSkinBinsPerSpan.MinV;
    const double max_v = rSkinBinsPerSpan.MaxV;

    auto compute_index = [](double coordinate, double min_value, double max_value, double span_size, std::size_t span_count){
        double clamped = coordinate;
        if (clamped < min_value) {
            clamped = min_value;
        } else if (clamped > max_value) {
            clamped = max_value;
        }

        if (span_count == 0) {
            return std::size_t(0);
        }

        double relative = (clamped - min_value) / span_size;
        std::size_t index = static_cast<std::size_t>(std::floor(relative + 1.0e-12));
        if (index >= span_count) {
            index = span_count - 1;
        }
        return index;
    };

    std::pair<IndexType, IndexType> best_pair = {std::numeric_limits<IndexType>::max(), std::numeric_limits<IndexType>::max()};
    double best_intersection_distance = std::numeric_limits<double>::max();
    std::vector<IndexType> closest_intersection_candidate_ids;
    auto update_candidate = [&](Node::Pointer p_node_a, Node::Pointer p_node_b, const CoordinatesArrayType& rIntersectionPoint){

        array_1d<double, 3> diff_intersection = rIntersectionPoint;
        diff_intersection -= rStartPoint;
        const double distance_intersection = norm_2(diff_intersection);
        const double distance_tolerance = 1.0e-12 * std::max(1.0, distance_intersection);
        if (distance_intersection < best_intersection_distance - distance_tolerance) {
            best_intersection_distance = distance_intersection;
            best_pair = {p_node_a->Id(), p_node_b->Id()};
            closest_intersection_candidate_ids.clear();
            closest_intersection_candidate_ids.push_back(p_node_a->Id());
            closest_intersection_candidate_ids.push_back(p_node_b->Id());
        } else if (std::abs(distance_intersection - best_intersection_distance) <= distance_tolerance) {
            best_intersection_distance = std::min(best_intersection_distance, distance_intersection);
            if (closest_intersection_candidate_ids.empty()) {
                best_pair = {p_node_a->Id(), p_node_b->Id()};
            }
            for (const IndexType candidate_id : {p_node_a->Id(), p_node_b->Id()}) {
                if (std::find(closest_intersection_candidate_ids.begin(), closest_intersection_candidate_ids.end(), candidate_id) ==
                    closest_intersection_candidate_ids.end()) {
                    closest_intersection_candidate_ids.push_back(candidate_id);
                }
            }
        }
    };

    const std::size_t number_spans_x = rSkinBinsPerSpan.NumberOfSpansX;
    const std::size_t number_spans_y = rSkinBinsPerSpan.NumberOfSpansY;
    const std::size_t base_ix = compute_index(rStartPoint[0], min_u, max_u, span_size_x, number_spans_x);
    const std::size_t base_iy = compute_index(rStartPoint[1], min_v, max_v, span_size_y, number_spans_y);

    const int max_search_level = static_cast<int>(std::max(number_spans_x, number_spans_y));

    ConditionPointerContainerType candidate_conditions;
    candidate_conditions.reserve(2000);

    const bool filter_by_layer = !rLayer.empty();
    auto node_has_requested_layer = [&](Node::Pointer p_candidate_node) {
        if (!filter_by_layer) {
            return true;
        }
        const auto& candidate_layers = p_candidate_node->GetValue(CONNECTED_LAYERS);
        return std::find(candidate_layers.begin(), candidate_layers.end(), rLayer) != candidate_layers.end();
    };

    auto condition_matches_requested_layer = [&](const ConditionPointerType& p_condition) {
        if (!filter_by_layer) {
            return true;
        }
        const bool condition_matches =
            p_condition->Has(LAYER_NAME) && p_condition->GetValue(LAYER_NAME) == rLayer;
        if (condition_matches) {
            return true;
        }

        const auto p_geometry = p_condition->pGetGeometry();
        if (p_geometry->size() < 2) {
            return false;
        }

        return node_has_requested_layer(p_geometry->pGetPoint(0)) || node_has_requested_layer(p_geometry->pGetPoint(1));
    };

    auto search_with_direction = [&](const array_1d<double, 3>& rSearchDirection) {
        bool found_intersection = false;
        for (int level = 0; level < max_search_level && !found_intersection; ++level) {
            const std::size_t extension = static_cast<std::size_t>(level+1);
            const double current_length = reference_span_size * static_cast<double>(extension);
            array_1d<double, 3> segment_start = rStartPoint - rSearchDirection * reference_span_size*0.5;
            array_1d<double, 3> segment_end = rStartPoint;
            segment_end += rSearchDirection * current_length;

            const std::size_t min_ix = (base_ix > extension) ? base_ix - extension : 0;
            const std::size_t max_ix = std::min<std::size_t>(base_ix + extension, number_spans_x > 0 ? number_spans_x - 1 : 0);
            const std::size_t min_iy = (base_iy > extension) ? base_iy - extension : 0;
            const std::size_t max_iy = std::min<std::size_t>(base_iy + extension, number_spans_y > 0 ? number_spans_y - 1 : 0);

            candidate_conditions.clear();

            Node::Pointer start_node = Kratos::make_intrusive<Node>(0, segment_start[0], segment_start[1], segment_start[2]);
            Node::Pointer end_node   = Kratos::make_intrusive<Node>(0, segment_end[0],   segment_end[1],   segment_end[2]);
            Condition::NodesArrayType ray_nodes;
            ray_nodes.reserve(2);
            ray_nodes.push_back(start_node);
            ray_nodes.push_back(end_node);
            ConditionPointerType p_ray_condition = Kratos::make_intrusive<Condition>(0, ray_nodes);

            for (std::size_t ix = min_ix; ix <= max_ix; ++ix) {
                if (ix >= rSkinBinsPerSpan.NumberOfSpansX) continue;
                for (std::size_t iy = min_iy; iy <= max_iy; ++iy) {
                    if (iy >= rSkinBinsPerSpan.NumberOfSpansY) continue;
                    const std::size_t k = FindSpanNnzIndex(rSkinBinsPerSpan, ix, iy, 0);
                    if (k == static_cast<std::size_t>(-1)) continue;
                    if (k >= rSkinBinsPerSpan.CellBinsByNnz.size()) continue;
                    const auto& r_cell_bins = rSkinBinsPerSpan.CellBinsByNnz[k];
                    if (!r_cell_bins.HasConditionBins) continue;
                    auto& r_bins = const_cast<BinsObjectDynamic<ConditionConfigure>&>(r_cell_bins.ConditionBins);
                    r_bins.SearchObjects(p_ray_condition, candidate_conditions);
                }
            }

            if (candidate_conditions.empty()) {
                continue;
            }

            for (const auto& p_condition : candidate_conditions) {
                if (!condition_matches_requested_layer(p_condition)) {
                    continue;
                }

                const auto p_geometry = p_condition->pGetGeometry();
                if (p_geometry->size() < 2) {
                    continue;
                }

                const auto& node_a = p_geometry->pGetPoint(0);
                const auto& node_b = p_geometry->pGetPoint(1);

                CoordinatesArrayType intersection_point;
                bool intersects = SegmentsIntersect(start_node, end_node, node_a, node_b, intersection_point);

                if (intersects) {
                    const array_1d<double, 3> intersection_delta = intersection_point - rStartPoint;
                    if (inner_prod(intersection_delta, rSearchDirection) < -dir_tol) {
                        continue;
                    }
                }

                const auto& normal_a = node_a->GetValue(NORMAL);
                const auto& normal_b = node_b->GetValue(NORMAL);
                const double normal_sign = ReverseNormalFilter ? -1.0 : 1.0;
                if (!IgnoreNormalFilter &&
                    normal_sign * inner_prod(normal_a, rSearchDirection) < 0.0 &&
                    normal_sign * inner_prod(normal_b, rSearchDirection) < 0.0) {
                    continue;
                }

                if (intersects) {
                    update_candidate(node_a, node_b, intersection_point);
                    found_intersection = true;
                }
            }

            if (best_pair.first != std::numeric_limits<IndexType>::max()) {
                return true;
            }
        }

        return found_intersection;
    };

    search_with_direction(direction);
    const auto forward_pair = best_pair;
    const double forward_distance = best_intersection_distance;
    auto forward_candidate_ids = std::move(closest_intersection_candidate_ids);

    best_pair = {std::numeric_limits<IndexType>::max(), std::numeric_limits<IndexType>::max()};
    best_intersection_distance = std::numeric_limits<double>::max();
    closest_intersection_candidate_ids.clear();
    const array_1d<double, 3> opposite_direction = -direction;
    search_with_direction(opposite_direction);

    auto has_unconsidered_candidate = [&](const std::vector<IndexType>& rCandidateIds) {
        return std::any_of(rCandidateIds.begin(), rCandidateIds.end(), [&](const IndexType CandidateId) {
            return pConsideredNodeIds == nullptr ||
                std::find(pConsideredNodeIds->begin(), pConsideredNodeIds->end(), CandidateId) == pConsideredNodeIds->end();
        });
    };
    const bool forward_has_new_candidate = has_unconsidered_candidate(forward_candidate_ids);
    const bool opposite_has_new_candidate = has_unconsidered_candidate(closest_intersection_candidate_ids);
    const double distance_tolerance = 1.0e-12 * std::max(1.0, forward_distance);
    // Keep the closest geometric intersection even when its nodes were already
    // used. In that case the caller stops subdividing this skin segment.
    if (forward_pair.first != std::numeric_limits<IndexType>::max() &&
        (best_pair.first == std::numeric_limits<IndexType>::max() ||
         forward_distance < best_intersection_distance - distance_tolerance ||
         (std::abs(forward_distance - best_intersection_distance) <= distance_tolerance &&
          (forward_has_new_candidate || !opposite_has_new_candidate)))) {
        best_pair = forward_pair;
        best_intersection_distance = forward_distance;
        closest_intersection_candidate_ids = std::move(forward_candidate_ids);
    }

    if (best_pair.first != std::numeric_limits<IndexType>::max()) {
        if (pClosestIntersectionCandidateIds != nullptr) {
            *pClosestIntersectionCandidateIds = closest_intersection_candidate_ids;
        }
        if (rLayer == "COUPLING_SIDE" && best_intersection_distance > reference_span_size / 2.0) {
            KRATOS_WARNING("SnakeGapSbmProcess")
                << "::[SnakeGapSbmProcess]::FindClosestPairInLayerWithNormalDirection: projection on COUPLING_SIDE farther than reference_span_size/2."
                << " distance=" << best_intersection_distance
                << ", threshold=" << (reference_span_size / 2.0)
                << ", start_point=" << rStartPoint
                << ", selected_pair=(" << best_pair.first << ", " << best_pair.second << ")"
                << std::endl;
        }
        return best_pair;
    }

    KRATOS_ERROR << "::[SnakeGapSbmProcess]::FindClosestPairInLayerWithNormalDirection: no intersection found for layer "
                 << (rLayer.empty() ? std::string("<any>") : rLayer) << " starting from point " << rStartPoint
                 << ". Skin sub model part '" << rSkinSubModelPart.Name() << "' has "
                 << rSkinSubModelPart.NumberOfNodes() << " nodes and "
                 << rSkinSubModelPart.NumberOfConditions() << " conditions."
                 << " Search parameters: base_ix=" << base_ix << ", base_iy=" << base_iy
                 << ", span_size_x=" << span_size_x << ", span_size_y=" << span_size_y
                 << ", reference_span_size=" << reference_span_size
                 << ", direction=" << direction << std::endl;

    return best_pair;
}

// Associate the surrogate reference geometry with intermediate skin nodes.
void SnakeGapSbmProcess::AttachSurrogateMiddleGeometryToSkinNodes(
    const ModelPart& rSkinSubModelPart,
    const std::vector<IndexType>& rSegmentProjectionIds,
    GeometryType::Pointer pSurrogateMiddleGeometry,
    const NodeType::Pointer& pProjectionNode1,
    const NodeType::Pointer& pProjectionNode2)
{
    // Also attach the surrogate middle geometry to all intermediate projection skin nodes
    std::unordered_set<IndexType> unique_projection_ids(
        rSegmentProjectionIds.begin(),
        rSegmentProjectionIds.end());

    unique_projection_ids.erase(std::numeric_limits<IndexType>::max());

    for (const auto id_proj : unique_projection_ids) {
        const auto p_skin_node = rSkinSubModelPart.pGetNode(id_proj);
        auto& r_neigh = p_skin_node->GetValue(NEIGHBOUR_GEOMETRIES);
        bool already_present = false;
        for (const auto& p_g : r_neigh) {
            if (p_g.get() == pSurrogateMiddleGeometry.get()) {
                already_present = true;
                break;
            }
        }
        if (!already_present) {
            r_neigh.push_back(pSurrogateMiddleGeometry);
        }
    }
}

// Merge neighbour geometries at the coincident endpoints of a closed skin.
void SnakeGapSbmProcess::SynchronizeEndSkinNodeNeighbourGeometries(
    const ModelPart& rSkinSubModelPart)
{
    // Assign same neighbour geometries to first and last node of the skin
    if (rSkinSubModelPart.NumberOfNodes() < 2) {
        return;
    }

    auto it_first = rSkinSubModelPart.NodesBegin();
    auto it_last = rSkinSubModelPart.NodesEnd();
    --it_last; // last valid iterator

    auto& r_first_node = *it_first;
    auto& r_last_node  = *it_last;

    auto& neigh_first = r_first_node.GetValue(NEIGHBOUR_GEOMETRIES);
    auto& neigh_last  = r_last_node.GetValue(NEIGHBOUR_GEOMETRIES);

    // Push back last's neighbours into first (avoid duplicates)
    for (const auto& p_g : neigh_last) {
        bool present = false;
        for (const auto& p_h : neigh_first) {
            if (p_g.get() == p_h.get()) {
                present = true;
                break;
            }
        }
        if (!present) {
            neigh_first.push_back(p_g);
        }
    }

    // Push back first's neighbours into last (avoid duplicates)
    for (const auto& p_g : neigh_first) {
        bool present = false;
        for (const auto& p_h : neigh_last) {
            if (p_g.get() == p_h.get()) {
                present = true;
                break;
            }
        }
        if (!present) {
            neigh_last.push_back(p_g);
        }
    }
}



//-------------------------------

// ------------------------------------------------------------------
// Evaluate curve point at parameter t in global coordinates
// ------------------------------------------------------------------
array_1d<double, 3> SnakeGapSbmProcess::EvaluateGlobalCurvePoint(
    const GeometryType& rCurve,
    const double T)
{
    GeometryType::CoordinatesArrayType local_coordinates;
    for (std::size_t i = 0; i < local_coordinates.size(); ++i) {
        local_coordinates[i] = 0.0;
    }
    local_coordinates[0] = T;

    array_1d<double, 3> global_point;
    for (std::size_t i = 0; i < 3; ++i) {
        global_point[i] = 0.0;
    }
    rCurve.GlobalCoordinates(global_point, local_coordinates);
    return global_point;
}

std::vector<SnakeGapSbmProcess::IndexType> SnakeGapSbmProcess::CollectOrderedSkinNodePath(
    const ModelPart& rSkinSubModelPart,
    IndexType StartNodeId,
    IndexType EndNodeId,
    bool ReverseConditionOrientation,
    const std::string& rLayerName,
    const std::string& rConditionName) const
{
    if (StartNodeId == EndNodeId) {
        return {StartNodeId};
    }

    struct PathEdge
    {
        IndexType NeighborId = 0;
        int OrientationScore = 0;
    };

    auto add_or_update_edge = [](
        std::unordered_map<IndexType, std::vector<PathEdge>>& rAdjacency,
        const IndexType FromId,
        const IndexType ToId,
        const int ScoreDelta) {
        auto& r_edges = rAdjacency[FromId];
        for (auto& r_edge : r_edges) {
            if (r_edge.NeighborId == ToId) {
                r_edge.OrientationScore += ScoreDelta;
                return;
            }
        }
        r_edges.push_back(PathEdge{ToId, ScoreDelta});
    };

    std::unordered_map<IndexType, std::vector<PathEdge>> adjacency;
    std::unordered_set<IndexType> matched_node_ids;
    std::size_t matched_condition_count = 0;

    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        if (!rLayerName.empty()) {
            if (!r_condition.Has(LAYER_NAME) || r_condition.GetValue(LAYER_NAME) != rLayerName) {
                continue;
            }
        }
        if (!rConditionName.empty()) {
            if (!r_condition.Has(CONDITION_NAME) || r_condition.GetValue(CONDITION_NAME) != rConditionName) {
                continue;
            }
        }

        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.PointsNumber() < 2) {
            continue;
        }

        const IndexType node_id_0 = r_geometry[0].Id();
        const IndexType node_id_1 = r_geometry[r_geometry.PointsNumber() - 1].Id();

        matched_node_ids.insert(node_id_0);
        matched_node_ids.insert(node_id_1);

        const int forward_score = ReverseConditionOrientation ? -1 : 1;
        const int backward_score = -forward_score;
        add_or_update_edge(adjacency, node_id_0, node_id_1, forward_score);
        add_or_update_edge(adjacency, node_id_1, node_id_0, backward_score);
        ++matched_condition_count;
    }

    KRATOS_ERROR_IF(matched_condition_count == 0)
        << "::[SnakeGapSbmProcess]::CollectOrderedSkinNodePath: no skin conditions matched layer='"
        << rLayerName << "' condition='" << rConditionName << "' in '"
        << rSkinSubModelPart.FullName() << "'." << std::endl;

    KRATOS_ERROR_IF(matched_node_ids.find(StartNodeId) == matched_node_ids.end())
        << "::[SnakeGapSbmProcess]::CollectOrderedSkinNodePath: start node #"
        << StartNodeId << " does not belong to the filtered skin conditions." << std::endl;
    KRATOS_ERROR_IF(matched_node_ids.find(EndNodeId) == matched_node_ids.end())
        << "::[SnakeGapSbmProcess]::CollectOrderedSkinNodePath: end node #"
        << EndNodeId << " does not belong to the filtered skin conditions." << std::endl;

    const auto& r_end_coords = rSkinSubModelPart.GetNode(EndNodeId).Coordinates();
    auto compare_neighbors = [&](IndexType IdA, IndexType IdB) {
        const auto& r_coords_a = rSkinSubModelPart.GetNode(IdA).Coordinates();
        const auto& r_coords_b = rSkinSubModelPart.GetNode(IdB).Coordinates();

        const double dist_a = norm_2(r_coords_a - r_end_coords);
        const double dist_b = norm_2(r_coords_b - r_end_coords);
        if (std::abs(dist_a - dist_b) > 1.0e-15) {
            return dist_a < dist_b;
        }
        if (std::abs(r_coords_a[0] - r_coords_b[0]) > 1.0e-15) {
            return r_coords_a[0] < r_coords_b[0];
        }
        if (std::abs(r_coords_a[1] - r_coords_b[1]) > 1.0e-15) {
            return r_coords_a[1] < r_coords_b[1];
        }
        return r_coords_a[2] < r_coords_b[2];
    };

    for (auto& r_entry : adjacency) {
        std::stable_sort(
            r_entry.second.begin(),
            r_entry.second.end(),
            [&](const PathEdge& rEdgeA, const PathEdge& rEdgeB) {
                if (rEdgeA.OrientationScore != rEdgeB.OrientationScore) {
                    return rEdgeA.OrientationScore > rEdgeB.OrientationScore;
                }
                return compare_neighbors(rEdgeA.NeighborId, rEdgeB.NeighborId);
            });
    }

    std::unordered_set<IndexType> visited_ids;
    visited_ids.reserve(matched_node_ids.size());
    visited_ids.insert(StartNodeId);

    std::vector<IndexType> current_path{StartNodeId};
    std::vector<IndexType> best_path;
    int best_orientation_score = std::numeric_limits<int>::min();

    auto is_better_candidate = [&](const int CandidateScore, const std::vector<IndexType>& rCandidatePath) {
        if (best_path.empty()) {
            return true;
        }
        if (CandidateScore != best_orientation_score) {
            return CandidateScore > best_orientation_score;
        }
        if (rCandidatePath.size() != best_path.size()) {
            return rCandidatePath.size() < best_path.size();
        }
        return false;
    };

    auto explore_paths = [&](auto&& rSelf, const IndexType CurrentId, const int CurrentScore) -> void {
        if (CurrentId == EndNodeId) {
            if (is_better_candidate(CurrentScore, current_path)) {
                best_orientation_score = CurrentScore;
                best_path = current_path;
            }
            return;
        }

        const auto adjacency_it = adjacency.find(CurrentId);
        if (adjacency_it == adjacency.end()) {
            return;
        }

        for (const auto& r_edge : adjacency_it->second) {
            const IndexType neighbor_id = r_edge.NeighborId;
            if (!visited_ids.insert(neighbor_id).second) {
                continue;
            }

            current_path.push_back(neighbor_id);
            rSelf(rSelf, neighbor_id, CurrentScore + r_edge.OrientationScore);
            current_path.pop_back();
            visited_ids.erase(neighbor_id);
        }
    };

    explore_paths(explore_paths, StartNodeId, 0);

    KRATOS_ERROR_IF(best_path.empty())
        << "::[SnakeGapSbmProcess]::CollectOrderedSkinNodePath: no topological path was found from node #"
        << StartNodeId << " to node #" << EndNodeId
        << " for layer='" << rLayerName << "' condition='" << rConditionName << "'."
        << std::endl;

    return best_path;
}

// Compute the radius of the gap corners about their centroid.
double SnakeGapSbmProcess::CalculateGapElementCharacteristicLength(
    const array_1d<double, 3>& rSurrogatePoint1,
    const array_1d<double, 3>& rSurrogatePoint2,
    const array_1d<double, 3>& rSkinPoint1,
    const array_1d<double, 3>& rSkinPoint2) const
{
    array_1d<double, 3> center = ZeroVector(3);
    center += rSurrogatePoint1;
    center += rSurrogatePoint2;
    center += rSkinPoint1;
    center += rSkinPoint2;
    center *= 0.25;

    double radius_squared = 0.0;
    const auto update_radius = [&center, &radius_squared](const array_1d<double, 3>& rPoint){
        array_1d<double, 3> difference = center - rPoint;
        const double distance_squared = inner_prod(difference, difference);
        if (distance_squared > radius_squared) {
            radius_squared = distance_squared;
        }
    };

    update_radius(rSurrogatePoint1);
    update_radius(rSurrogatePoint2);
    update_radius(rSkinPoint1);
    update_radius(rSkinPoint2);

    radius_squared += 1E-15; //numerical safety to avoid 0

    return std::sqrt(radius_squared);
}

// ------------------------------------------------------------------
// Coons patch mapping X(ξ,η)
// ------------------------------------------------------------------
array_1d<double, 3> SnakeGapSbmProcess::CoonsPoint(
    const double Xi,
    const double Eta,
    const BrepCurveType& rB0,
    const BrepCurveType& rL0,
    const BrepCurveType& rL1,
    const BrepCurveType& rB1,
    const array_1d<double, 3>&     rP00,
    const array_1d<double, 3>&     rP01,
    const array_1d<double, 3>&     rP10,
    const array_1d<double, 3>&     rP11)
{
    const array_1d<double, 3> brep_B0_at_xi = EvaluateGlobalCurvePoint(rB0, Xi);
    const array_1d<double, 3> brep_B1_at_xi = EvaluateGlobalCurvePoint(rB1, Xi);
    const array_1d<double, 3> brep_L0_at_eta = EvaluateGlobalCurvePoint(rL0, Eta);
    const array_1d<double, 3> brep_L1_at_eta = EvaluateGlobalCurvePoint(rL1, Eta);

    const double om_xi  = 1.0 - Xi;
    const double om_eta = 1.0 - Eta;

    return  om_xi * brep_L0_at_eta + Xi * brep_L1_at_eta
          + om_eta * brep_B0_at_xi + Eta * brep_B1_at_xi
          - om_xi * om_eta * rP00
          - Xi    * om_eta * rP10
          - Xi    * Eta    * rP11
          - om_xi * Eta    * rP01;
}

// Evaluate the analytical derivative of the Coons mapping.
array_1d<double, 3> SnakeGapSbmProcess::CoonsDerivative(
    const double Xi,
    const double Eta,
    const bool with_respect_to_xi,
    const BrepCurveType& rB0,
    const BrepCurveType& rL0,
    const BrepCurveType& rL1,
    const BrepCurveType& rB1,
    const array_1d<double, 3>&     rP00,
    const array_1d<double, 3>&     rP01,
    const array_1d<double, 3>&     rP10,
    const array_1d<double, 3>&     rP11)
{

    // Evaluate dX/ds at normalized s (chain rule: dX/ds = dX/dt * dt/ds = dX/dt * (±len))
    auto EvalDerAt = [&](const BrepCurveType& rCurve, double t) -> array_1d<double, 3> {
        GeometryType::CoordinatesArrayType lc;
        lc.resize(rCurve.LocalSpaceDimension(), false);
        lc[0] = t;

        std::vector<CoordinatesArrayType> D; // D[0]=X, D[1]=dX/dt
        rCurve.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->GlobalSpaceDerivatives(D, lc, 1);

        array_1d<double, 3> dXds = ZeroVector(3);
        dXds[0] = D[1][0];
        dXds[1] = D[1][1];
        dXds[2] = D[1][2];
        return dXds;
    };

    // Boundary values and derivatives (w.r.t. ξ or η in [0,1])

    const array_1d<double, 3> brep_B0_at_xi = EvaluateGlobalCurvePoint(rB0, Xi);
    const array_1d<double, 3> brep_B1_at_xi = EvaluateGlobalCurvePoint(rB1, Xi);
    const array_1d<double, 3> brep_L0_at_eta = EvaluateGlobalCurvePoint(rL0, Eta);
    const array_1d<double, 3> brep_L1_at_eta = EvaluateGlobalCurvePoint(rL1, Eta);

    const array_1d<double, 3> dB0_at_xi = EvalDerAt(rB0, Xi);   // dB0_at_xi/dξ (xi)
    const array_1d<double, 3> dB1_at_xi = EvalDerAt(rB1, Xi);   // dB1_at_xi/dξ (xi)
    const array_1d<double, 3> dL0_at_eta = EvalDerAt(rL0, Eta);  // dL0_at_eta/dη (eta)
    const array_1d<double, 3> dL1_at_eta = EvalDerAt(rL1, Eta);  // dL1_at_eta/dη (eta)

    array_1d<double, 3> dS = ZeroVector(3);
    if (with_respect_to_xi) {
        dS  = (1.0 - Eta) * dB0_at_xi + Eta * dB1_at_xi;     // boundary derivative blend
        dS += -brep_L0_at_eta + brep_L1_at_eta;                          // side linear term
        dS += (1.0 - Eta) * rP00 - (1.0 - Eta) * rP10 + Eta * rP01 - Eta * rP11; // corner correction
    } else {
        dS  = -brep_B0_at_xi + brep_B1_at_xi;                           // bottom/top linear term
        dS += (1.0 - Xi) * dL0_at_eta + Xi * dL1_at_eta;        // side derivative blend
        dS += (1.0 - Xi) * rP00 + Xi * rP10 - (1.0 - Xi) * rP01 - Xi * rP11;     // corner correction
    }
    return dS;
}

SnakeGapSbmProcess::IntegrationPointsArrayType
SnakeGapSbmProcess::CreateCoonsPatchGaussPoints(
    const std::size_t             Order,
    const BrepCurveType& rB0,
    const BrepCurveType& rL0,
    const BrepCurveType& rL1,
    const BrepCurveType& rB1,
    const array_1d<double, 3>&     rP00,
    const array_1d<double, 3>&     rP01,
    const array_1d<double, 3>&     rP10,
    const array_1d<double, 3>&     rP11) const
{
    IntegrationPointsArrayType gp_list;
    gp_list.reserve(Order * Order);

    // 1-D Gauss nodes / weights on [0,1]
    IntegrationPointUtilities::IntegrationPointsArrayType integration_points(Order);
    auto integration_points_it = integration_points.begin();
    IntegrationPointUtilities::IntegrationPoints1D(integration_points_it, Order, 0.0, 1.0);

    std::vector<double> xi(Order);
    std::vector<double> w(Order);
    for (std::size_t i = 0; i < Order; ++i) {
        xi[i] = integration_points[i].X();
        w[i] = integration_points[i].Weight();
    }

    // --- build a reference normal to define the positive orientation ----
    array_1d<double, 3> a = rP10 - rP00;
    array_1d<double, 3> b = rP01 - rP00;
    // --- reference normal from Coons derivatives at center --------------
    const double xi_c  = 0.5, eta_c = 0.5;
    const array_1d<double, 3> dXi_c  = CoonsDerivative(xi_c, eta_c, true ,
                                                        rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
    const array_1d<double, 3> dEta_c = CoonsDerivative(xi_c, eta_c, false,
                                                        rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
    array_1d<double, 3> n_ref = MathUtils<double>::CrossProduct(dXi_c, dEta_c);

    // normalize (if non-degenerate)
    const double nrm = norm_2(n_ref);
    if (nrm > 0.0) n_ref /= nrm;

    for (std::size_t i=0; i<Order; ++i)
    for (std::size_t j=0; j<Order; ++j) {
        const double xi_i  = xi[i];
        const double eta_j = xi[j];
        double w_ij  = w[i] * w[j];

        // --- Jacobian (signed) ------------------------------------------
        const array_1d<double, 3> dXi  = CoonsDerivative(xi_i, eta_j, true ,
                                                          rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
        const array_1d<double, 3> dEta = CoonsDerivative(xi_i, eta_j, false,
                                                          rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
        const array_1d<double, 3> cr   = MathUtils<double>::CrossProduct(dXi, dEta);

        // Signed Jacobian: projection onto reference normal (can be < 0)
        const double jac_signed = inner_prod(cr, n_ref);

        // --- Global coordinates -----------------------------------------
        const array_1d<double, 3> X =
            CoonsPoint(xi_i, eta_j, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);

        array_1d<double, 3> e1 = rP10 - rP00;
        array_1d<double, 3> e2 = X - rP00;

        // Store global coords + signed weight
        gp_list.emplace_back( IntegrationPoint<3>( X[0], X[1], X[2],
                                                   w_ij * jac_signed ) );
    }
    return gp_list;
}

template void SnakeGapSbmProcess::SetSurrogateToSkinProjections2D<true>(
    const ModelPart& rSurrogateSubModelPart,
    const ModelPart& rSkinSubModelPart,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan);

template void SnakeGapSbmProcess::SetSurrogateToSkinProjections2D<false>(
    const ModelPart& rSurrogateSubModelPart,
    const ModelPart& rSkinSubModelPart,
    const KnotSpanSkinBinsCSR& rSkinBinsPerSpan);

} // namespace Kratos
