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
#include <filesystem>
#include <fstream>
#include <iomanip>
// External includes

// Project includes
#include "snake_gap_sbm_process.h"
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"

namespace Kratos
{

namespace
{
using CurvePointerType = SnakeGapSbmProcess::NurbsCurveGeometryType::Pointer;

void AppendPolynomialCurveInfo(const CurvePointerType& pCurve, const IndexType skin1Id, const IndexType skin2Id)
{
    if (!pCurve) {
        KRATOS_WARNING("SnakeGapSbmProcess") << "AppendPolynomialCurveInfo called with null curve pointer." << std::endl;
        return;
    }

    std::error_code mkdir_error;
    std::filesystem::create_directories("txt_files", mkdir_error);
    if (mkdir_error) {
        KRATOS_WARNING("SnakeGapSbmProcess") << "Unable to create directory txt_files: " << mkdir_error.message() << std::endl;
    }

    std::ofstream out("txt_files/polynomial_curves.txt", std::ios::app);
    if (!out.is_open()) {
        KRATOS_WARNING("SnakeGapSbmProcess") << "Unable to open txt_files/polynomial_curves.txt for writing." << std::endl;
        return;
    }

    out << std::setprecision(15);
    out << "curve skin1_id=" << skin1Id << " skin2_id=" << skin2Id;
    out << " degree=" << pCurve->PolynomialDegree(0);

    const auto& r_knots = pCurve->Knots();
    out << " knots=[";
    for (std::size_t i = 0; i < r_knots.size(); ++i) {
        if (i > 0) out << ' ';
        out << r_knots[i];
    }
    out << "]";

    const auto& r_weights = pCurve->Weights();
    const std::size_t number_of_control_points = pCurve->size();
    out << " weights=[";
    if (r_weights.size() > 0) {
        for (std::size_t i = 0; i < r_weights.size(); ++i) {
            if (i > 0) out << ' ';
            out << r_weights[i];
        }
    } else {
        for (std::size_t i = 0; i < number_of_control_points; ++i) {
            if (i > 0) out << ' ';
            out << 1.0;
        }
    }
    out << "]";

    out << " control_points=[";
    for (std::size_t i = 0; i < number_of_control_points; ++i) {
        const auto& r_cp = (*pCurve)[i];
        if (i > 0) out << ' ';
        out << "(" << r_cp.Id() << "," << r_cp.X() << "," << r_cp.Y() << "," << r_cp.Z() << ")";
    }
    out << "]\n";
}
} // unnamed namespace

SnakeGapSbmProcess::SnakeGapSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    SnakeSbmProcess(rModel, ThisParameters)
{
    
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("gap_element_name")) << "::[SnakeGapSbmProcess]::" 
                    << "Missing \"gap_element_name\" section." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("gap_interface_condition_name")) << "::[SnakeGapSbmProcess]::" 
                    << "Missing \"gap_interface_condition_name\" section." << std::endl;

    ThisParameters.AddMissingParameters(this->GetDefaultParameters());
    
    mpGapElementsSubModelPart = &(mpIgaModelPart->CreateSubModelPart("GapElements"));
    mpGapInterfaceSubModelPart = &(mpIgaModelPart->CreateSubModelPart("GapInterfaces"));
    mGapElementName = ThisParameters["gap_element_name"].GetString();
    mGapInterfaceConditionName = ThisParameters["gap_interface_condition_name"].GetString();
    mGapSbmType = ThisParameters["gap_sbm_type"].GetString(); 

    // TODO: in future PR, implement also the "sbm" type
    if (mGapSbmType != "default" && mGapSbmType != "interpolation" /*&& mGapSbmType != "sbm"*/) {
        KRATOS_ERROR << "::[SnakeGapSbmProcess]::"
                     << "The gap_sbm_type \"" << mGapSbmType << "\" is not supported. Available options are: "
                     << "default, interpolation." << std::endl;
    }

    if (ThisParameters.Has("gap_approximation_order"))
        mGapApproximationOrder = ThisParameters["gap_approximation_order"].GetInt();
    if (ThisParameters.Has("use_for_multipatch")) 
        mUseForMultipatch = ThisParameters["use_for_multipatch"].GetBool();
    if (ThisParameters.Has("polynomial_order"))
        mGapInterpolationOrder = ThisParameters["polynomial_order"].GetVector()[0];

    mLambdaInner = 0.0;
    mLambdaOuter = 1.0;
    mInternalDivisions = ThisParameters["number_internal_divisions"].GetInt();

    if (mUseForMultipatch) mLambdaInner = 0.001;
}

SnakeGapSbmProcess::KnotSpanIdsCSR
SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix(
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart) const
{
    KnotSpanIdsCSR knot_span_data;

    // --- Read parameters from parent model part (as in your code) ---
    const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();

    const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);
    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] KNOT_SPAN_SIZES must have at least two entries.\n";

    const auto& parameter_space_corners = r_parent_model_part.GetValue(PARAMETER_SPACE_CORNERS);
    KRATOS_ERROR_IF(parameter_space_corners.size() < 2)
        << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] PARAMETER_SPACE_CORNERS must have at least two vectors.\n";
    KRATOS_ERROR_IF(parameter_space_corners[0].size() < 2 || parameter_space_corners[1].size() < 2)
        << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] PARAMETER_SPACE_CORNERS vectors must contain min and max values.\n";

    const double span_size_x = knot_span_sizes[0];
    const double span_size_y = knot_span_sizes[1];
    KRATOS_ERROR_IF(span_size_x <= 0.0 || span_size_y <= 0.0)
        << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Knot span sizes must be positive.\n";

    const double min_u = parameter_space_corners[0][0];
    const double max_u = parameter_space_corners[0][1];
    const double min_v = parameter_space_corners[1][0];
    const double max_v = parameter_space_corners[1][1];

    const double domain_length_u = max_u - min_u;
    const double domain_length_v = max_v - min_v;
    KRATOS_ERROR_IF(domain_length_u <= 0.0 || domain_length_v <= 0.0)
        << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Invalid parameter space extents.\n";

    const double number_of_spans_x_real = domain_length_u / span_size_x;
    const double number_of_spans_y_real = domain_length_v / span_size_y;

    const auto compute_span_count = [](double spans_real) -> std::size_t {
        constexpr double tolerance = 1.0e-10;
        std::size_t span_count = static_cast<std::size_t>(std::round(spans_real));
        KRATOS_ERROR_IF(span_count <= 0)
            << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Non-positive number of knot spans.\n";
        KRATOS_ERROR_IF(std::abs(spans_real - static_cast<double>(span_count)) > tolerance)
            << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Non-integer number of knot spans (" << spans_real << ").\n";
        return span_count;
    };

    const std::size_t number_of_spans_x = compute_span_count(number_of_spans_x_real);
    const std::size_t number_of_spans_y = compute_span_count(number_of_spans_y_real);

    // --- Fill metadata ---
    knot_span_data.NumberOfSpansX = number_of_spans_x;
    knot_span_data.NumberOfSpansY = number_of_spans_y;
    knot_span_data.MinU = min_u;
    knot_span_data.MaxU = max_u;
    knot_span_data.MinV = min_v;
    knot_span_data.MaxV = max_v;
    knot_span_data.SpanSizeX = span_size_x;
    knot_span_data.SpanSizeY = span_size_y;

    // Early exit: no nodes
    if (rSkinSubModelPart.NumberOfNodes() == 0) {
        knot_span_data.Occupancy.resize(number_of_spans_x, number_of_spans_y, false);
        return knot_span_data;
    }

    const double tolerance = 1.0e-12;
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
                << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] coordinate below minimum parameter range.\n";
            clamped = min_value;
        } else if (coordinate > max_value) {
            KRATOS_ERROR_IF(coordinate > max_value + tolerance)
                << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] coordinate above maximum parameter range.\n";
            clamped = max_value;
        }
        const double relative = (clamped - min_value) / span_size;
        std::size_t idx = static_cast<std::size_t>(std::floor(relative + tolerance));
        if (idx >= span_count) idx = span_count - 1;
        return idx;
    };

    // --- Pass 1: gather unique columns per row for sparsity pattern ---
    std::vector<std::vector<std::size_t>> column_indices_per_row(number_of_spans_x);

    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const double u = r_node.X();
        const double v = r_node.Y();

        KRATOS_ERROR_IF(u < min_u_with_tolerance || u > max_u_with_tolerance)
            << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] node " << r_node.Id()
            << " has u-parameter " << u << " outside [" << min_u << ", " << max_u << "].\n";
        KRATOS_ERROR_IF(v < min_v_with_tolerance || v > max_v_with_tolerance)
            << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] node " << r_node.Id()
            << " has v-parameter " << v << " outside [" << min_v << ", " << max_v << "].\n";

        const std::size_t span_index_x = compute_span_index(u, min_u, max_u, span_size_x, number_of_spans_x);
        const std::size_t span_index_y = compute_span_index(v, min_v, max_v, span_size_y, number_of_spans_y);
        column_indices_per_row[span_index_x].push_back(span_index_y);
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

    for (std::size_t span_index_x = 0; span_index_x < number_of_spans_x; ++span_index_x) {
        for (const std::size_t span_index_y : column_indices_per_row[span_index_x]) {
            r_occupancy_matrix.push_back(span_index_x, span_index_y, 0.0);  // keep strictly increasing (i,j)
        }
    }

    // Temporary per-nnz buckets for node ids
    std::vector<std::vector<IndexType>> node_ids_per_non_zero(number_of_non_zero_entries);

    // --- Pass 2: fill tmp lists and counts in A.value_data() ---
    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const std::size_t span_index_x = compute_span_index(r_node.X(), min_u, max_u, span_size_x, number_of_spans_x);
        const std::size_t span_index_y = compute_span_index(r_node.Y(), min_v, max_v, span_size_y, number_of_spans_y);

        const std::size_t non_zero_index = FindNnzIndex(r_occupancy_matrix, span_index_x, span_index_y);
        KRATOS_DEBUG_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1))
            << "[SnakeGapSbmProcess::CreateSkinNodesPerKnotSpanMatrix] nonzero (ix,iy) not found in CSR pattern.\n";

        node_ids_per_non_zero[non_zero_index].push_back(r_node.Id());
        r_occupancy_matrix.value_data()[non_zero_index] += 1.0;  // optional: store count per cell
    }

    // --- Pack payload to compact pool ---
    CommitPayload(node_ids_per_non_zero, knot_span_data);

    return knot_span_data;
}

// === Builder: conditions -> KnotSpanIdsCSR ===
SnakeGapSbmProcess::KnotSpanIdsCSR
SnakeGapSbmProcess::CreateSkinConditionsPerKnotSpanMatrix(
    const ModelPart& rSkinSubModelPart,
    const SnakeGapSbmProcess::KnotSpanIdsCSR& rReferenceMatrix) const
{
    KnotSpanIdsCSR knot_span_data;

    // Copy metadata from reference
    knot_span_data.NumberOfSpansX = rReferenceMatrix.NumberOfSpansX;
    knot_span_data.NumberOfSpansY = rReferenceMatrix.NumberOfSpansY;
    knot_span_data.MinU = rReferenceMatrix.MinU;
    knot_span_data.MaxU = rReferenceMatrix.MaxU;
    knot_span_data.MinV = rReferenceMatrix.MinV;
    knot_span_data.MaxV = rReferenceMatrix.MaxV;
    knot_span_data.SpanSizeX = rReferenceMatrix.SpanSizeX;
    knot_span_data.SpanSizeY = rReferenceMatrix.SpanSizeY;

    const std::size_t number_of_spans_x = knot_span_data.NumberOfSpansX;
    const std::size_t number_of_spans_y = knot_span_data.NumberOfSpansY;
    if (number_of_spans_x == 0 || number_of_spans_y == 0) {
        knot_span_data.Occupancy.resize(0, 0, false);
        return knot_span_data;
    }
    if (rSkinSubModelPart.NumberOfConditions() == 0) {
        knot_span_data.Occupancy.resize(number_of_spans_x, number_of_spans_y, false);
        return knot_span_data;
    }

    const double min_u = knot_span_data.MinU;
    const double max_u = knot_span_data.MaxU;
    const double min_v = knot_span_data.MinV;
    const double max_v = knot_span_data.MaxV;
    const double span_size_u = knot_span_data.SpanSizeX;
    const double span_size_v = knot_span_data.SpanSizeY;
    const double tolerance = 1e-12;

    const auto compute_span_index = [tolerance](
        double coordinate,
        double min_value,
        double max_value,
        double span_size,
        std::size_t span_count) {
        double clamped_coordinate = coordinate;
        if (coordinate < min_value) {
            KRATOS_ERROR_IF(coordinate < min_value - tolerance) << "coord < min\n";
            clamped_coordinate = min_value;
        } else if (coordinate > max_value) {
            KRATOS_ERROR_IF(coordinate > max_value + tolerance) << "coord > max\n";
            clamped_coordinate = max_value;
        }
        std::size_t span_index = static_cast<std::size_t>(
            std::floor((clamped_coordinate - min_value) / span_size + tolerance));
        if (span_index >= span_count) {
            span_index = span_count - 1;
        }
        return span_index;
    };

    // Sparsity pattern
    std::vector<std::vector<std::size_t>> column_indices_per_row(number_of_spans_x);
    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.size() == 0) {
            continue;
        }
        array_1d<double, 3> centroid = ZeroVector(3);
        for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
            centroid += r_geometry[point_index].Coordinates();
        }
        centroid /= static_cast<double>(r_geometry.size());
        const std::size_t span_index_x = compute_span_index(centroid[0], min_u, max_u, span_size_u, number_of_spans_x);
        const std::size_t span_index_y = compute_span_index(centroid[1], min_v, max_v, span_size_v, number_of_spans_y);
        column_indices_per_row[span_index_x].push_back(span_index_y);
    }

    std::size_t number_of_non_zero_entries = 0;
    for (auto& r_column_indices : column_indices_per_row) {
        std::sort(r_column_indices.begin(), r_column_indices.end());
        r_column_indices.erase(std::unique(r_column_indices.begin(), r_column_indices.end()), r_column_indices.end());
        number_of_non_zero_entries += r_column_indices.size();
    }

    auto& r_occupancy_matrix = knot_span_data.Occupancy;
    r_occupancy_matrix.resize(number_of_spans_x, number_of_spans_y, false);
    r_occupancy_matrix.reserve(number_of_non_zero_entries);
    for (std::size_t span_index_x = 0; span_index_x < number_of_spans_x; ++span_index_x) {
        for (std::size_t span_index_y : column_indices_per_row[span_index_x]) {
            r_occupancy_matrix.push_back(span_index_x, span_index_y, 0.0);
        }
    }

    // Fill
    std::vector<std::vector<IndexType>> condition_ids_per_non_zero(number_of_non_zero_entries);
    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geometry = r_condition.GetGeometry();
        if (r_geometry.size() == 0) {
            continue;
        }
        array_1d<double, 3> centroid = ZeroVector(3);
        for (IndexType point_index = 0; point_index < r_geometry.size(); ++point_index) {
            centroid += r_geometry[point_index].Coordinates();
        }
        centroid /= static_cast<double>(r_geometry.size());
        const std::size_t span_index_x = compute_span_index(centroid[0], min_u, max_u, span_size_u, number_of_spans_x);
        const std::size_t span_index_y = compute_span_index(centroid[1], min_v, max_v, span_size_v, number_of_spans_y);
        const std::size_t non_zero_index = FindNnzIndex(r_occupancy_matrix, span_index_x, span_index_y);
        KRATOS_DEBUG_ERROR_IF(non_zero_index == static_cast<std::size_t>(-1)) << "pattern miss\n";
        condition_ids_per_non_zero[non_zero_index].push_back(r_condition.Id());
        r_occupancy_matrix.value_data()[non_zero_index] += 1.0;
    }

    CommitPayload(condition_ids_per_non_zero, knot_span_data);
    return knot_span_data;
}

void SnakeGapSbmProcess::CreateSbmExtendedGeometries()
{
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0 || mpSkinModelPartInnerInitial->NumberOfGeometries()>0) 
    {
        const auto& r_surrogate_sub_model_part_inner = mpIgaModelPart->GetSubModelPart("surrogate_inner");
        const auto& r_skin_sub_model_part_inner = mpSkinModelPart->GetSubModelPart("inner");

        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the inner skin." << std::endl;
        CreateSbmExtendedGeometries<true>(r_skin_sub_model_part_inner, r_surrogate_sub_model_part_inner);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the inner skin." << std::endl;
    }
    if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) 
    {
        const auto& r_surrogate_sub_model_part_outer = mpIgaModelPart->GetSubModelPart("surrogate_outer");
        const auto& r_skin_sub_model_part_outer = mpSkinModelPart->GetSubModelPart("outer");
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the outer skin." << std::endl;
        CreateSbmExtendedGeometries<false>(r_skin_sub_model_part_outer, r_surrogate_sub_model_part_outer);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the outer skin." << std::endl;
    }
}

template <bool TIsInnerLoop>
void SnakeGapSbmProcess::CreateSbmExtendedGeometries(
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart)
{
    // Get the mesh sizes from the surrogate model part
    const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    double knot_span_reference_size = knot_span_sizes[0];
    if (knot_span_sizes[1] > knot_span_reference_size) {
        knot_span_reference_size = knot_span_sizes[1];
    }
    if (knot_span_sizes.size() > 2) {
        if (knot_span_sizes[2] > knot_span_reference_size) {
            knot_span_reference_size = knot_span_sizes[2];
        }
    }

    auto p_surface = mpIgaModelPart->pGetGeometry(1);    
    auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
                            p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // Build knot-span lookups and set projections from surrogate to skin
    const auto& skin_nodes_per_knot_span = CreateSkinNodesPerKnotSpanMatrix(rSkinSubModelPart, rSurrogateSubModelPart);
    const auto& skin_conditions_per_knot_span = CreateSkinConditionsPerKnotSpanMatrix(rSkinSubModelPart, skin_nodes_per_knot_span);

    SetSurrogateToSkinProjections<TIsInnerLoop>(rSurrogateSubModelPart, rSkinSubModelPart, skin_nodes_per_knot_span);
    // Loop over the nodes of the surrogate sub model part
    IndexType element_id = rSurrogateSubModelPart.ElementsBegin()->Id();
    std::size_t brep_degree = p_nurbs_surface->PolynomialDegree(0);
    std::size_t number_of_shape_functions_derivatives = 2 * brep_degree + 1;

    if (mGapApproximationOrder == 0)
        mGapApproximationOrder = brep_degree;

    IndexType first_condition_id;
    IndexType last_condition_id;
    IndexType starting_brep_id;
    std::size_t size_surrogate_loop;

    if constexpr (TIsInnerLoop)  {
        first_condition_id = rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[0].Id();
        last_condition_id = rSurrogateSubModelPart.pGetElement(element_id)->GetGeometry()[1].Id();
        size_surrogate_loop = last_condition_id - first_condition_id + 1;
        if (mpSkinModelPartOuterInitial->NumberOfNodes() > 0 || mpSkinModelPartOuterInitial->NumberOfGeometries() > 0) 
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

    for (std::size_t j = 0; j < size_surrogate_loop; ++j) {
        auto p_brep_geometry = mpIgaModelPart->pGetGeometry(starting_brep_id + j);
        auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

        KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeGapSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
                                            << " is not a BrepCurveOnSurfaceType." << std::endl;

        NurbsInterval brep_domain_interval = p_brep_curve_on_surface_surrogate1_surrogate2->DomainInterval();
        CoordinatesArrayType surrogate_vertex_1 = ZeroVector(3); 
        CoordinatesArrayType surrogate_vertex_1_local_coords = ZeroVector(3);
        CoordinatesArrayType surrogate_vertex_2 = ZeroVector(3); 
        CoordinatesArrayType surrogate_vertex_2_local_coords = ZeroVector(3);
        surrogate_vertex_1_local_coords[0] = brep_domain_interval.GetT0();
        surrogate_vertex_2_local_coords[0] = brep_domain_interval.GetT1();

        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

        // retrieve middle point of the brep (used as surrogate reference for the dof)
        CoordinatesArrayType surrogate_middle_point = ZeroVector(3); 
        CoordinatesArrayType surrogate_middle_point_local_coords = ZeroVector(3);
        surrogate_middle_point_local_coords[0] = 0.5 * (surrogate_vertex_1_local_coords[0] + surrogate_vertex_2_local_coords[0]);
        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_middle_point, surrogate_middle_point_local_coords);

        IntegrationPoint<1> integration_point(surrogate_middle_point_local_coords[0]);
        IntegrationPointsArrayType surrogate_integration_points_list;
        surrogate_integration_points_list.push_back(integration_point);

        IntegrationInfo integration_info = p_brep_curve_on_surface_surrogate1_surrogate2->GetDefaultIntegrationInfo();

        IntegrationParameters integration_parameters(
            number_of_shape_functions_derivatives, 
            integration_info, 
            knot_span_sizes);
        integration_parameters.pSkinNodesPerSpan = &skin_nodes_per_knot_span;
        integration_parameters.pSkinConditionsPerSpan = &skin_conditions_per_knot_span;

        GeometriesArrayType quadrature_point_list;
        p_brep_curve_on_surface_surrogate1_surrogate2->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives, 
                                                            surrogate_integration_points_list, integration_info);

        GeometryType::Pointer surrogate_brep_middle_geometry = quadrature_point_list(0);

        // Store the surrogate middle geometry for the lateral Breps
        const double tol = 1.0e-12;
        bool check_cond_1 = false;
        bool check_cond_2 = false;
        Node::Pointer p_surrogate_node_1;
        Node::Pointer p_surrogate_node_2;
        for (auto& r_surrogate_node : rSurrogateSubModelPart.Nodes())
        {
            // Check if this node coincides with surrogate_vertex_1
            if (norm_2(r_surrogate_node.Coordinates() - surrogate_vertex_1) < tol)
            {
                // Direct reference to the stored neighbour geometries
                auto& r_neighbour_geometries =
                    r_surrogate_node.GetValue(NEIGHBOUR_GEOMETRIES);
                
                p_surrogate_node_1 = &r_surrogate_node;

                // Append the new middle geometry; no SetValue call required
                r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);

                check_cond_1 = true; // Only one node should satisfy the condition
            }
            else if (norm_2(r_surrogate_node.Coordinates() - surrogate_vertex_2) < tol)
            {
                // Direct reference to the stored neighbour geometries
                auto& r_neighbour_geometries =
                    r_surrogate_node.GetValue(NEIGHBOUR_GEOMETRIES);

                p_surrogate_node_2 = &r_surrogate_node;

                // Append the new middle geometry; no SetValue call required
                r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);

                check_cond_2 = true; // Only one node should satisfy the condition
            }
            if (check_cond_1 && check_cond_2) {
                // If both conditions are satisfied, break the loop
                break;
            }
        }

        array_1d<double,3> vector_cond = -p_surrogate_node_1->Coordinates() + p_surrogate_node_2->Coordinates();
        array_1d<double,3> normal_cond;
        normal_cond[0] = vector_cond[1];
        normal_cond[1] = -vector_cond[0];
        normal_cond[2] =  0.0;

        const double t0 = brep_domain_interval.GetT0();
        const double t1 = brep_domain_interval.GetT1();

        //------------------------------------------------------------------
        // 3. Loop over sub-intervals of the upper curve
        //------------------------------------------------------------------
        Node::Pointer p_first_node = p_surrogate_node_1;
        Node::Pointer p_second_node = nullptr;

        auto connected_layers_1 = p_surrogate_node_1->GetValue(CONNECTED_LAYERS);
        auto connected_layers_2 = p_surrogate_node_2->GetValue(CONNECTED_LAYERS);
        const auto projection_id_1 = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
        const auto projection_id_2 = p_surrogate_node_2->GetValue(PROJECTION_NODE_ID);
        const auto projection_node_1 = rSkinSubModelPart.pGetNode(projection_id_1);
        const auto projection_node_2 = rSkinSubModelPart.pGetNode(projection_id_2);
        std::string common_layer_name = "";
        std::string common_condition_name = "";
        IndexType condition_count = 0;
        bool layer_found = false;
        // Find the common layer between the two surrogate nodes
        for (auto& layer : connected_layers_1) {
            for (auto& layer_2 : connected_layers_2) {
                if (layer == layer_2) {
                    common_layer_name = layer;
                    common_condition_name = p_surrogate_node_1->GetValue(CONNECTED_CONDITIONS)[condition_count];
                    layer_found = true;
                    break;                  
                }
            }
            condition_count++;
        }

        std::size_t current_condition_internal_divisions = mInternalDivisions;

        if (norm_2(projection_node_2->Coordinates() - projection_node_1->Coordinates()) < knot_span_reference_size/1e5) 
            current_condition_internal_divisions = 0;

        const std::size_t subdivision_depth = current_condition_internal_divisions;
        const std::size_t segment_count = subdivision_depth == 0 ? 1 : static_cast<std::size_t>(1) << subdivision_depth;
    
        const double dt_subdivision = (t1 - t0) / static_cast<double>(segment_count);

        std::vector<Node::Pointer> surrogate_segment_nodes(segment_count + 1);
        surrogate_segment_nodes[0] = p_surrogate_node_1;
        surrogate_segment_nodes[segment_count] = p_surrogate_node_2;

        for (std::size_t k = 1; k < segment_count; ++k) {
            CoordinatesArrayType local_coords = ZeroVector(3);
            CoordinatesArrayType global_coords = ZeroVector(3);
            local_coords[0] = t0 + dt_subdivision * static_cast<double>(k);
            p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(global_coords, local_coords);

            Node::Pointer p_new_surrogate_node = Node::Pointer(new Node(0, global_coords));

            auto& r_new_connected_layers = p_new_surrogate_node->GetValue(CONNECTED_LAYERS);
            if (layer_found) {
                r_new_connected_layers.push_back(common_layer_name);
            }
            auto& r_new_connected_conditions = p_new_surrogate_node->GetValue(CONNECTED_CONDITIONS);
            if (!common_condition_name.empty()) {
                r_new_connected_conditions.push_back(common_condition_name);
            }

            surrogate_segment_nodes[k] = p_new_surrogate_node;
        }

        std::vector<IndexType> segment_projection_ids(segment_count + 1, std::numeric_limits<IndexType>::max());
        segment_projection_ids[0] = projection_id_1;
        segment_projection_ids[segment_count] = projection_id_2;

        auto compute_recursive_projection = [&](auto&& self,
                                                std::size_t depth,
                                                std::size_t left_index,
                                                std::size_t right_index) -> void
        {
            if (depth == 0) {
                return;
            }

            const std::size_t mid_index = (left_index + right_index) / 2;
            if (segment_projection_ids[mid_index] != std::numeric_limits<IndexType>::max()) {
                if (depth > 1) {
                    self(self, depth - 1, left_index, mid_index);
                    self(self, depth - 1, mid_index, right_index);
                }
                return;
            }

            const auto& r_left_skin_node = rSkinSubModelPart.GetNode(segment_projection_ids[left_index]);
            const auto& r_right_skin_node = rSkinSubModelPart.GetNode(segment_projection_ids[right_index]);

            CoordinatesArrayType skin_mid_point_coords = 0.5 * (r_left_skin_node.Coordinates() + r_right_skin_node.Coordinates());

            Vector normal_direction = r_right_skin_node.Coordinates() - r_left_skin_node.Coordinates();
            normal_direction /= norm_2(normal_direction);

            double temp = normal_direction[0];
            normal_direction[0] = -normal_direction[1];
            normal_direction[1] = temp;
            const double normal_norm = norm_2(normal_direction);
            if (normal_norm > 1.0e-16) {
                normal_direction /= normal_norm;
            }
            IndexType id_skin_node = -1;
            if (norm_2(r_left_skin_node.Coordinates()-r_right_skin_node.Coordinates()) < integration_parameters.rKnotSpanSizes[0]/1e5)
                id_skin_node = segment_projection_ids[left_index];
            
            else
            {   
                id_skin_node = FindClosestNodeInLayerWithDirection<TIsInnerLoop>(
                    skin_mid_point_coords,
                    common_layer_name,
                    rSkinSubModelPart,
                    integration_parameters.rKnotSpanSizes,
                    skin_conditions_per_knot_span,
                    normal_direction);

            
                // choose if dot product is < 0
                auto candidate_surrogate_to_skin_node = rSkinSubModelPart.GetNode(id_skin_node).Coordinates() - surrogate_segment_nodes[mid_index]->Coordinates();
                const double dot_product = inner_prod(candidate_surrogate_to_skin_node, normal_cond);
                if (dot_product < -1e-12) {
                    // Avoid that the projection cross the surrogate condition
                    const auto split_pair = FindClosestPairInLayerWithNormalDirection<TIsInnerLoop>(
                        skin_mid_point_coords,
                        common_layer_name,
                        rSkinSubModelPart,
                        integration_parameters.rKnotSpanSizes,
                        skin_conditions_per_knot_span,
                        normal_direction);
                    
                    // choose positive inner prod
                    auto candidate_surrogate_to_skin_node_1 = rSkinSubModelPart.GetNode(split_pair.first).Coordinates() - surrogate_segment_nodes[mid_index]->Coordinates();
                    auto candidate_surrogate_to_skin_node_2 = rSkinSubModelPart.GetNode(split_pair.second).Coordinates() - surrogate_segment_nodes[mid_index]->Coordinates();
                    const double dot_prod_1 = inner_prod(candidate_surrogate_to_skin_node_1, normal_cond);
                    const double dot_prod_2 = inner_prod(candidate_surrogate_to_skin_node_2, normal_cond);
                    if (dot_prod_1 >= 0 && dot_prod_2 < 0) {
                        id_skin_node = split_pair.first;
                    } else if (dot_prod_2 >= 0 && dot_prod_1 < 0) {
                        id_skin_node = split_pair.second;
                    } else {
                        KRATOS_ERROR << "Both candidate split nodes have negative inner product with normal direction.\n";
                    }
                }
            }

            segment_projection_ids[mid_index] = id_skin_node;

            if (depth > 1) {
                self(self, depth - 1, left_index, mid_index);
                self(self, depth - 1, mid_index, right_index);
            }
        };

        if (segment_count > 1) {
            compute_recursive_projection(compute_recursive_projection, subdivision_depth, 0, segment_count);
        }

        for (std::size_t k = 1; k < segment_count; ++k) {
            const IndexType projection_id = segment_projection_ids[k];
            KRATOS_ERROR_IF(projection_id == std::numeric_limits<IndexType>::max())
                << "::[SnakeGapSbmProcess]:: Missing projection id for surrogate subdivision node at index "
                << k << std::endl;
            surrogate_segment_nodes[k]->SetValue(PROJECTION_NODE_ID, projection_id);
        }

        for (std::size_t segment = 0; segment < segment_count; ++segment) {
            auto p_first = surrogate_segment_nodes[segment];
            auto p_second = surrogate_segment_nodes[segment + 1];

            // Use inner products
            const CoordinatesArrayType center = surrogate_brep_middle_geometry->Center();
            const IndexType proj_id_first  = p_first->GetValue(PROJECTION_NODE_ID);
            const IndexType proj_id_second = p_second->GetValue(PROJECTION_NODE_ID);

            array_1d<double,3> v_proj_1 = rSkinSubModelPart.GetNode(proj_id_first).Coordinates() - p_first->Coordinates();
            array_1d<double,3> v_proj_2 = rSkinSubModelPart.GetNode(proj_id_second).Coordinates() - p_second->Coordinates();

            // v_proj_1 /= norm_2(v_proj_1); //FIXME:
            // v_proj_2 /= norm_2(v_proj_2);

            const double dot1 = v_proj_1[0]*normal_cond[0] + v_proj_1[1]*normal_cond[1];
            const double dot2 = v_proj_2[0]*normal_cond[0] + v_proj_2[1]*normal_cond[1]; 

            // dot1
            if (dot1 < -1.0e-12) {
                Vector direction(3);
                direction[0] = vector_cond[0];
                direction[1] = vector_cond[1];
                direction[2] = vector_cond[2];
                KRATOS_ERROR_IF(norm_2(direction) < 1e-13) << "Error: zero direction vector.\n";
                direction = direction / norm_2(direction);
                
                // if constexpr(TIsInnerLoop)
                //     direction *=-1;
                const auto split_pair = FindClosestPairInLayerWithNormalDirection<true>(
                    p_first->Coordinates(),
                    common_layer_name,
                    rSkinSubModelPart,
                    integration_parameters.rKnotSpanSizes,
                    skin_conditions_per_knot_span,
                    direction);
                                    
                // Only one should satisfy the inner product condition
                array_1d<double,3> v_pair_1 = rSkinSubModelPart.GetNode(split_pair.first).Coordinates() - p_first->Coordinates();
                array_1d<double,3> v_pair_2 = rSkinSubModelPart.GetNode(split_pair.second).Coordinates() - p_second->Coordinates();

                const double dot_pair1 = v_pair_1[0]*normal_cond[0] + v_pair_1[1]*normal_cond[1];
                const double dot_pair2 = v_pair_2[0]*normal_cond[0] + v_pair_2[1]*normal_cond[1];

                Node::Pointer p_split_node;
                if (dot_pair1 > 0) {
                    p_split_node = rSkinSubModelPart.pGetNode(split_pair.first);
                } else if (dot_pair2 > 0) {
                    p_split_node = rSkinSubModelPart.pGetNode(split_pair.second);
                } else {
                    KRATOS_ERROR << "Neither candidate split nodes satisfy the inner product condition.\n";
                }

                auto p_fake_surrogate_node = Node::Pointer(new Node(0, p_first->Coordinates()));

                p_fake_surrogate_node->SetValue(PROJECTION_NODE_ID, p_split_node->Id());
                const auto& connected_layers = p_split_node->GetValue(CONNECTED_LAYERS);
                const auto& connected_conditions = p_split_node->GetValue(CONNECTED_CONDITIONS);
                p_fake_surrogate_node->SetValue(CONNECTED_LAYERS, connected_layers);
                p_fake_surrogate_node->SetValue(CONNECTED_CONDITIONS, connected_conditions);

                // Create split segments: [p_first -> p_split_node] and [p_split_node -> p_second]
                CreateGapAndSkinQuadraturePoints<TIsInnerLoop>(
                    integration_parameters,
                    p_nurbs_surface,
                    p_first,
                    p_fake_surrogate_node,
                    surrogate_brep_middle_geometry,
                    *mpIgaModelPart,
                    rSkinSubModelPart);

                p_first = p_fake_surrogate_node;
            }
            if (dot2 < -1.0e-12) {
                Vector direction(3);
                direction[0] = vector_cond[0];
                direction[1] = vector_cond[1];
                direction[2] = vector_cond[2];
                KRATOS_ERROR_IF(norm_2(direction) < 1e-13) << "Error: zero direction vector.\n";
                direction = direction/norm_2(direction);

                // if constexpr(TIsInnerLoop) //FIXME: check this
                //     direction *=-1;

                const auto split_pair = FindClosestPairInLayerWithNormalDirection<false>(
                    p_second->Coordinates(),
                    common_layer_name,
                    rSkinSubModelPart,
                    integration_parameters.rKnotSpanSizes,
                    skin_conditions_per_knot_span,
                    direction);
                

                // Only one satisfy inner product condition
                array_1d<double,3> v_pair_1 = rSkinSubModelPart.GetNode(split_pair.first).Coordinates() - p_first->Coordinates();
                array_1d<double,3> v_pair_2 = rSkinSubModelPart.GetNode(split_pair.second).Coordinates() - p_second->Coordinates();

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
                    KRATOS_ERROR << "Neither candidate split nodes satisfy the inner product condition.\n";
                }

                auto p_fake_surrogate_node = Node::Pointer(new Node(0, p_second->Coordinates()));

                p_fake_surrogate_node->SetValue(PROJECTION_NODE_ID, p_split_node->Id());
                const auto& connected_layers = p_split_node->GetValue(CONNECTED_LAYERS);
                const auto& connected_conditions = p_split_node->GetValue(CONNECTED_CONDITIONS);
                p_fake_surrogate_node->SetValue(CONNECTED_LAYERS, connected_layers);
                p_fake_surrogate_node->SetValue(CONNECTED_CONDITIONS, connected_conditions);

                CreateGapAndSkinQuadraturePoints<TIsInnerLoop>(
                    integration_parameters,
                    p_nurbs_surface,
                    p_fake_surrogate_node,
                    p_second,
                    surrogate_brep_middle_geometry,
                    *mpIgaModelPart,
                    rSkinSubModelPart);

                p_second = p_fake_surrogate_node;
            } 
            // Fallback: standard behavior
            CreateGapAndSkinQuadraturePoints<TIsInnerLoop>(
                integration_parameters,
                p_nurbs_surface,
                p_first,
                p_second,
                surrogate_brep_middle_geometry,
                *mpIgaModelPart,
                rSkinSubModelPart);
        }

        if (mUseForMultipatch) {
            // Attach the surrogate middle geometry to all associated skin nodes
            AttachSurrogateMiddleGeometryToSkinNodes(
                rSkinSubModelPart,
                segment_projection_ids,
                surrogate_brep_middle_geometry,
                projection_node_1,
                projection_node_2);
        }
    }

    if (mUseForMultipatch) {
        // Ensure first and last skin nodes share the same neighbour geometries
        SynchronizeEndSkinNodeNeighbourGeometries(rSkinSubModelPart);
    }

    //---------------------------------------------------------------------------
    bool is_entering = false;
    for (IndexType i_cond_id = first_condition_id; i_cond_id <= last_condition_id; ++i_cond_id)
    {
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
            
            Node::Pointer p_first_point = Node::Pointer(new Node(1, surrogate_1));
            Node::Pointer p_second_point = Node::Pointer(new Node(2, skin_1));

            // check if the interface length is zero -> not necessary to create quadrature points there
            if (norm_2(surrogate_1-skin_1) < 1e-13) continue;

            if (is_entering) {
                // change the order to preserve the anticlockwise orientation
                Node::Pointer p_temp_pointer = p_first_point;
                p_first_point = p_second_point;
                p_second_point = p_temp_pointer;
            } 

            auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_first_point, p_second_point, active_range_knot_vector);
            auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
        
            IntegrationInfo brep_integration_info_surrogate1_skin1 = p_brep_curve_surrogate1_skin1->GetDefaultIntegrationInfo();

            brep_integration_info_surrogate1_skin1.SetNumberOfIntegrationPointsPerSpan(0,mGapInterpolationOrder+1);

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
                const auto& temp_geom = neighbour_geometries[0];
                neighbour_geometries[0] = neighbour_geometries[1];
                neighbour_geometries[1] = temp_geom;

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

            Node::Pointer p_first_point = Node::Pointer(new Node(1, surrogate_2));
            Node::Pointer p_second_point = Node::Pointer(new Node(2, skin_2));

            // check if the interface length is zero -> not necessary to create quadrature points there
            if (norm_2(surrogate_2-skin_2) < 1e-13) continue;

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

            brep_integration_info_surrogate2_skin2.SetNumberOfIntegrationPointsPerSpan(0,mGapInterpolationOrder+1);

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
                const auto& temp_geom = neighbour_geometries[0];
                neighbour_geometries[0] = neighbour_geometries[1];
                neighbour_geometries[1] = temp_geom;

                p_surrogate_2->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
            } else {
                KRATOS_ERROR << "::[SnakeSbmProcess]:: The surrogate middle point is not close to any of the neighbour geometries." << std::endl;
            }

            this->CreateConditions(
                brep_quadrature_point_list_surrogate2_skin2.ptr_begin(), brep_quadrature_point_list_surrogate2_skin2.ptr_end(),
                *mpGapInterfaceSubModelPart, mGapInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries, characteristic_condition_length);

        }
    }
    
    // For the inner loop, create multipatch coupling conditions on the skin
    if constexpr (TIsInnerLoop) {
        if (mUseForMultipatch) {
            CreateInnerSkinMultipatchCouplingConditions(
                rSkinSubModelPart,
                knot_span_sizes,
                p_nurbs_surface);
        }
    }
}

template <bool TIsInnerLoop>
void SnakeGapSbmProcess::CreateGapAndSkinQuadraturePoints(
    IntegrationParameters& rIntegrationParameters,
    NurbsSurfaceType::Pointer& pNurbsSurface,
    const Node::Pointer& pSurrogateNode1, 
    const Node::Pointer& pSurrogateNode2, 
    const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
    ModelPart& rIgaModelPart,
    const ModelPart& rSkinSubModelPart)
{
    const IndexType id_closest_true_node = pSurrogateNode1->GetValue(PROJECTION_NODE_ID);

    const auto& p_skin_node_1 = rSkinSubModelPart.pGetNode(id_closest_true_node);
    
    // search the projection of the second vertex
    const IndexType id_closest_true_node_2 = pSurrogateNode2->GetValue(PROJECTION_NODE_ID);

    const auto& p_skin_node_2 = rSkinSubModelPart.pGetNode(id_closest_true_node_2);

    // retrieve condition name for the skin condition
    auto connected_layers_1 = pSurrogateNode1->GetValue(CONNECTED_LAYERS);
    auto connected_layers_2 = pSurrogateNode2->GetValue(CONNECTED_LAYERS);
    std::string layer_name = "";   
    std::string condition_name = "";
    IndexType condition_count = 0;
    bool layer_found = false;
    // Find the common layer between the two surrogate nodes
    for (auto& layer : connected_layers_1) {
        for (auto& layer_2 : connected_layers_2) {
            if (layer == layer_2) {
                layer_name = layer;
                condition_name = pSurrogateNode1->GetValue(CONNECTED_CONDITIONS)[condition_count];
                layer_found = true;
                break;                  
            }
        }
        condition_count++;
    }

    KRATOS_ERROR_IF(!layer_found) << ":::[SnakeGapSbmProcess]::: No common layer found between the two surrogate nodes "
                                    << pSurrogateNode1->Id() << " and " << pSurrogateNode2->Id() <<  "\n"
                                    << pSurrogateNode1->Coordinates() << pSurrogateNode2->Coordinates() << "\n"
                                    << connected_layers_1 << connected_layers_2 << std::endl;

    KRATOS_ERROR_IF(!rIntegrationParameters.pSkinConditionsPerSpan)
        << "::[SnakeGapSbmProcess]::CreateGapAndSkinQuadraturePoints: skin condition span matrix is not initialized." << std::endl;

    const auto& r_skin_conditions_per_span = *rIntegrationParameters.pSkinConditionsPerSpan;

    ModelPart& r_layer_model_part = rIgaModelPart.HasSubModelPart(layer_name) ? 
                                    rIgaModelPart.GetSubModelPart(layer_name) : 
                                    rIgaModelPart.CreateSubModelPart(layer_name);       
                                    
    Vector active_range_knot_vector = ZeroVector(2);
    active_range_knot_vector[0] = 0;
    active_range_knot_vector[1] = 1;

    // surrogate_1 - skin_1
    // create the brep connecting vertex and closest true point
    Point skin_1(*p_skin_node_1);
    
    Node::Pointer p_skin1_brep_point = Node::Pointer(new Node(2, skin_1));
    
    auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(pSurrogateNode1, p_skin1_brep_point, active_range_knot_vector);
    auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
                                            
    // surrogate_2 - skin_2
    // create the brep connecting vertex and closest true point
    Point skin_2(*p_skin_node_2);
    
    Node::Pointer p_skin2_brep_point = Node::Pointer(new Node(2, skin_2));

    auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(pSurrogateNode2, p_skin2_brep_point, active_range_knot_vector);
    auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);      

    // skin_1 - skin_2
    int p = mGapApproximationOrder;

    auto p_nurbs_curve_skin1_skin2 = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);

    Vector temp_interpolation_nodes_id;
    std::vector<IndexType> interpolation_nodes_id;
    std::vector<array_1d<double,3>> interpolation_points;

    if (mGapSbmType == "interpolation")
    {
        std::size_t number_interpolation_cuts = (p-1);
        if (!(norm_2(skin_2 - skin_1) > rIntegrationParameters.rKnotSpanSizes[0]/50 && p>1))
            number_interpolation_cuts = 0;


        std::size_t segment_count = number_interpolation_cuts == 0 ? 1 : static_cast<std::size_t>(1) << number_interpolation_cuts;
        const IndexType invalid_projection_id = std::numeric_limits<IndexType>::max();

        std::vector<IndexType> interpolation_projection_ids(segment_count + 1, invalid_projection_id);
        interpolation_projection_ids[0] = id_closest_true_node;
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

            Vector normal_direction = ZeroVector(3);
            if (tangent_norm > 1.0e-16) {
                normal_direction = tangent / tangent_norm;

                const double temp = normal_direction[0];
                normal_direction[0] = -normal_direction[1];
                normal_direction[1] = temp;

                const double normal_norm = norm_2(normal_direction);
                if (normal_norm > 1.0e-16) {
                    normal_direction /= normal_norm;

                    if constexpr(!TIsInnerLoop)
                        normal_direction*=-1;
                }
            }

            IndexType id_skin_node = invalid_projection_id;
            if (tangent_norm >= rIntegrationParameters.rKnotSpanSizes[0] / 1.0e5) {
                const auto split_pair = FindClosestPairInLayerWithNormalDirection<TIsInnerLoop>(
                    skin_mid_point_coords,
                    layer_name,
                    rSkinSubModelPart,
                    rIntegrationParameters.rKnotSpanSizes,
                    r_skin_conditions_per_span,
                    normal_direction);

                std::vector<std::pair<double, IndexType>> candidates;
                auto consider_candidate = [&](IndexType candidate_id) {
                    if (candidate_id == std::numeric_limits<IndexType>::max()) {
                        return;
                    }
                    const auto& r_candidate_node = rSkinSubModelPart.GetNode(candidate_id);
                    const double distance = norm_2(r_candidate_node.Coordinates() - skin_mid_point_coords);
                    candidates.emplace_back(distance, candidate_id);
                };
                consider_candidate(split_pair.first);
                consider_candidate(split_pair.second);
                std::sort(candidates.begin(), candidates.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });

                const auto already_inserted = [&](IndexType candidate_id) {
                    return std::find(interpolation_projection_ids.begin(), interpolation_projection_ids.end(), candidate_id)
                        != interpolation_projection_ids.end();
                };

                bool added_new_node = false;
                for (const auto& [distance, candidate_id] : candidates) {
                    if (!already_inserted(candidate_id)) {
                        id_skin_node = candidate_id;
                        added_new_node = true;
                    }
                }

                if (!added_new_node && !candidates.empty()) {
                    id_skin_node = candidates.front().second;
                } else if (candidates.empty()) {
                    id_skin_node = left_node_id;
                }

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
                                            id_closest_true_node, id_closest_true_node_2);
        }

        interpolation_projection_ids[recursion_segment_count] = id_closest_true_node_2;

        std::unordered_set<IndexType> seen_projection_ids;

        auto new_end = std::remove_if(
            interpolation_projection_ids.begin(),
            interpolation_projection_ids.end(),
            [&](IndexType projection_id) {
                if (projection_id == invalid_projection_id) {
                    return true; // tieni gli invalidi (se li vuoi tenere)
                }
                // se l'inserimento fallisce (second == false), è un duplicato -> rimuovi
                return !seen_projection_ids.insert(projection_id).second;
            });
        
        interpolation_projection_ids.erase(new_end, interpolation_projection_ids.end());

        segment_count = interpolation_projection_ids.size();

        interpolation_nodes_id.reserve(segment_count);
        temp_interpolation_nodes_id.resize(segment_count);
        interpolation_points.reserve(segment_count);

        for (std::size_t i = 0; i < segment_count; ++i) {
            const IndexType projection_id = interpolation_projection_ids[i];
            KRATOS_ERROR_IF(projection_id == invalid_projection_id)
                << "::[SnakeGapSbmProcess]:: Missing interpolation projection id at subdivision index " << i << std::endl;

            const auto& r_skin_node = rSkinSubModelPart.GetNode(projection_id);
            interpolation_nodes_id.push_back(projection_id);
            temp_interpolation_nodes_id[i] = projection_id;
            interpolation_points.push_back(r_skin_node.Coordinates());
        }

        const double ridge = 1e-11;

        if (norm_2(skin_2 - skin_1) > rIntegrationParameters.rKnotSpanSizes[0]/50 && p>1 && segment_count >= p+1)
            p_nurbs_curve_skin1_skin2 = FitBezierUV_LS_Generic(interpolation_points, p, ridge);
    }
    else if (mGapSbmType == "default")
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
        IndexType current_id = id_closest_true_node;
        ModelPart::SizeType iter = 0;
        const ModelPart::SizeType max_iterations = rSkinSubModelPart.NumberOfNodes();
        while (current_id != id_closest_true_node_2) {
            current_id = advance_to_next_skin_id(current_id);
            ++id_distance;
            ++iter;
            KRATOS_ERROR_IF(iter > max_iterations)
                << "::[SnakeGapSbmProcess]:: id_distance loop between node IDs "
                << id_closest_true_node << " and " << id_closest_true_node_2
                << " exceeded the number of skin nodes (" << max_iterations << ").\n";
        }

        if (id_distance > static_cast<std::size_t>(2*p+1)) {
            check_enough_samples = true;
        }

        if (norm_2(skin_2 - skin_1) > rIntegrationParameters.rKnotSpanSizes[0]/50 && check_enough_samples && p>1)
        {
            p_nurbs_curve_skin1_skin2 = FitUV_BetweenSkinNodes_Generic<TIsInnerLoop>(
                rSkinSubModelPart, *pNurbsSurface, id_closest_true_node, id_closest_true_node_2, p, /*ridge=*/1e-14);
        }
    }

    AppendPolynomialCurveInfo(p_nurbs_curve_skin1_skin2, p_skin_node_1->Id(), p_skin_node_2->Id());

    auto p_brep_curve_skin1_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2);    
    
    IntegrationPointsArrayType brep_integration_points_list_skin1_skin2;
    GeometriesArrayType brep_quadrature_point_list_skin1_skin2;
    rIntegrationParameters.CurveIntegrationInfo.SetNumberOfIntegrationPointsPerSpan(0, ((mGapInterpolationOrder+1)));

    p_brep_curve_skin1_skin2->CreateIntegrationPoints(brep_integration_points_list_skin1_skin2, rIntegrationParameters.CurveIntegrationInfo);

    // fix the value of the weights
    for (auto& integration_point : brep_integration_points_list_skin1_skin2) {
        std::vector<CoordinatesArrayType> ders;
        p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->GlobalSpaceDerivatives(ders, integration_point.Coordinates(), 1);
        const double weight_correcting_factor = std::sqrt(ders[1][0]*ders[1][0] +ders[1][1]*ders[1][1]);

        integration_point.SetWeight(integration_point.Weight() * weight_correcting_factor);
    }

    if (norm_2(skin_2 - skin_1) > 1e-12)
        p_brep_curve_skin1_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin1_skin2, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                                  brep_integration_points_list_skin1_skin2, rIntegrationParameters.CurveIntegrationInfo);
    
    if (mGapSbmType == "interpolation")
        for (auto& quad_geom : brep_quadrature_point_list_skin1_skin2)
            quad_geom.SetValue(INTERPOLATION_NODES_ID, interpolation_nodes_id);

    std::size_t id = 1;
    if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
        id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;
    
    std::vector<Geometry<Node>::Pointer> neighbour_geometries_skin1_skin2;
    neighbour_geometries_skin1_skin2.push_back(rSurrogateBrepMiddleGeometry);
    
    // TODO: next PR --- GapSbmType == SBM
    // // Find the closest point on the true boundary for each of the quadrature geometries
    // for (auto& quad_geom: brep_quadrature_point_list_skin1_skin2)
    // {
    //     CoordinatesArrayType quadrature_point_coords = quad_geom.Center().Coordinates();
    //     DynamicBinsPointerType p_quadrature_point = DynamicBinsPointerType(new PointType(1, quadrature_point_coords[0], quadrature_point_coords[1], quadrature_point_coords[2]));
        
    //     auto normal = quad_geom.Normal(0);
    //     normal = normal / MathUtils<double>::Norm(normal);

    //     // IndexType id_skin_node = FindClosestNodeInLayer(p_quadrature_point, rBinSearchParameters, layer_name, rSkinSubModelPart);
        
    //     IndexType id_skin_node = FindClosestNodeInLayerWithDirection<TIsInnerLoop>(quadrature_point_coords, layer_name, rSkinSubModelPart, rIntegrationParameters.rKnotSpanSizes,
    //         r_skin_conditions_per_span, normal);
    //     NodeType::Pointer p_skin_node = rSkinSubModelPart.pGetNode(id_skin_node);
    //     NodePointerVector empty_vector;

    //     if constexpr(TIsInnerLoop)
    //         p_skin_node->SetValue(IDENTIFIER, "inner");
    //     else 
    //         p_skin_node->SetValue(IDENTIFIER, "outer");
    //     empty_vector.push_back(p_skin_node); // Just it_node-plane neighbours
    //     quad_geom.SetValue(NEIGHBOUR_NODES, empty_vector);
    // }
    const double characteristic_condition_length = CalculateGapElementCharacteristicLength(
        pSurrogateNode1->Coordinates(), pSurrogateNode2->Coordinates(), p_skin_node_1->Coordinates(), p_skin_node_2->Coordinates());
    
    if (!mUseForMultipatch)
        this->CreateConditions(
            brep_quadrature_point_list_skin1_skin2.ptr_begin(), brep_quadrature_point_list_skin1_skin2.ptr_end(),
            r_layer_model_part, condition_name, id, PropertiesPointerType(), rIntegrationParameters.rKnotSpanSizes, neighbour_geometries_skin1_skin2, characteristic_condition_length);

    // check for void elements/true coincident with surrogate boundary
    if (
        (std::abs(pSurrogateNode1->X() - p_skin_node_1->X()) < 1e-14 &&
         std::abs(pSurrogateNode1->X() - p_skin_node_2->X()) < 1e-14 &&
         std::abs(pSurrogateNode1->X() - pSurrogateNode2->X()) < 1e-14)  ||
        (std::abs(pSurrogateNode1->Y() - p_skin_node_1->Y()) < 1e-14 &&
         std::abs(pSurrogateNode1->Y() - p_skin_node_2->Y()) < 1e-14 &&
         std::abs(pSurrogateNode1->Y() - pSurrogateNode2->Y()) < 1e-14))
        return;
    
    // surrogate_1 - surrogate_2
    // create ONLY the brep connecting the two vertices
    auto p_nurbs_curve_surrogate1_surrogate2 = this->CreateBrepCurve(pSurrogateNode1, pSurrogateNode2, active_range_knot_vector);
    auto p_brep_curve_surrogate1_surrogate2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_surrogate2);      

    // Creation integration points on the cut elements
    GeometriesArrayType surface_quadrature_point_list;
    IntegrationInfo surface_integration_info = pNurbsSurface->GetDefaultIntegrationInfo();
    double characteristic_length = 0.0;
    if constexpr (TIsInnerLoop)  
    {

        array_1d<double,3> surrogate_1 = pSurrogateNode1->Coordinates();
        array_1d<double,3> surrogate_2 = pSurrogateNode2->Coordinates();
        array_1d<double,3> skin_1 = p_skin_node_1->Coordinates();
        array_1d<double,3> skin_2 = p_skin_node_2->Coordinates();

        IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
            (mGapInterpolationOrder+1), /*Order*/
            *p_brep_curve_surrogate1_surrogate2,   // B0
            *p_brep_curve_surrogate1_skin1,       // L0
            *p_brep_curve_surrogate2_skin2,       // L1
            *p_brep_curve_skin1_skin2,            // B1
            surrogate_1,  // P00
            skin_1,       // P01
            surrogate_2,  // P10
            skin_2);      // P11
        
        characteristic_length = CalculateGapElementCharacteristicLength(surrogate_1, surrogate_2, skin_1, skin_2);
        pNurbsSurface->CreateQuadraturePointGeometries(surface_quadrature_point_list, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
            surface_integration_points, surface_integration_info);

    }
    else
    {
        auto p_nurbs_curve_skin2_skin1 = ReverseBezierUV_Generic(p_nurbs_curve_skin1_skin2);
        auto p_brep_curve_skin2_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin2_skin1);  

        auto p_nurbs_curve_surrogate2_surrogate1 = ReverseBezierUV_Generic(p_nurbs_curve_surrogate1_surrogate2);
        auto p_brep_curve_surrogate2_surrogate1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_surrogate1);  

        array_1d<double,3> surrogate_1 = pSurrogateNode1->Coordinates();
        array_1d<double,3> surrogate_2 = pSurrogateNode2->Coordinates();
        array_1d<double,3> skin_1 = p_skin_node_1->Coordinates();
        array_1d<double,3> skin_2 = p_skin_node_2->Coordinates();

        IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
            (mGapInterpolationOrder+1), /*Order*/
            *p_brep_curve_surrogate2_surrogate1,   // B0
            *p_brep_curve_surrogate2_skin2,       // L0
            *p_brep_curve_surrogate1_skin1,       // L1
            *p_brep_curve_skin2_skin1,            // B1
            surrogate_2,  // P00
            skin_2,       // P01
            surrogate_1,  // P10
            skin_1);      // P11

        characteristic_length = CalculateGapElementCharacteristicLength(surrogate_1, surrogate_2, skin_1, skin_2);
        pNurbsSurface->CreateQuadraturePointGeometries(surface_quadrature_point_list, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                        surface_integration_points, surface_integration_info);
    }
    KRATOS_ERROR_IF(characteristic_length <= 0.0)
        << "Characteristic length for the gap element must be positive." << std::endl;

    IndexType id_element = 1;
    if (mpGapElementsSubModelPart->GetRootModelPart().Elements().size() > 0)
        id_element = mpGapElementsSubModelPart->GetRootModelPart().Elements().back().Id() + 1;

    this->CreateElements(
        surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
        *mpGapElementsSubModelPart, mGapElementName, id_element, PropertiesPointerType(), neighbour_geometries_skin1_skin2, characteristic_length);
}

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
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    KRATOS_ERROR_IF(CharacteristicLength <= 0.0)
        << "Characteristic length for gap conditions must be positive." << std::endl;

    IndexType geometry_count = 0;
    array_1d<double,3> characteristic_length_vector = ZeroVector(3);
    characteristic_length_vector[0] = CharacteristicLength;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_condition_list.push_back(
            reference_condition.Create(rIdCounter, (*it), pProperties));
        
        // Set knot span sizes to the condition
        new_condition_list.GetContainer()[geometry_count]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);

        new_condition_list.GetContainer()[geometry_count]->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);
        new_condition_list.GetContainer()[geometry_count]->SetValue(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_vector);

        rIdCounter++;
        geometry_count++;
    }

    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}


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

    KRATOS_INFO_IF("CreateElements", mEchoLevel > 2)
        << "Creating elements of type " << rElementName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    std::size_t num_elements = std::distance(rGeometriesBegin, rGeometriesEnd);
    new_element_list.reserve(num_elements);
    IndexType geometry_count = 0;

    array_1d<double,3> characteristic_length_vector = ZeroVector(3);
    characteristic_length_vector[0] = CharacteristicLength;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_element_list.push_back(
            rReferenceElement.Create(rIdCounter, (*it), pProperties));

        new_element_list.GetContainer()[geometry_count]->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);
        new_element_list.GetContainer()[geometry_count]->SetValue(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_length_vector);
        rIdCounter++;
        geometry_count++;
    }

    rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
}

double SnakeGapSbmProcess::CalculateGapElementCharacteristicLength(
    const array_1d<double,3>& rSurrogatePoint1,
    const array_1d<double,3>& rSurrogatePoint2,
    const array_1d<double,3>& rSkinPoint1,
    const array_1d<double,3>& rSkinPoint2) const
{
    array_1d<double,3> center = ZeroVector(3);
    center += rSurrogatePoint1;
    center += rSurrogatePoint2;
    center += rSkinPoint1;
    center += rSkinPoint2;
    center *= 0.25;

    double radius_squared = 0.0;
    const auto update_radius = [&center, &radius_squared](const array_1d<double,3>& rPoint){
        array_1d<double,3> difference = center - rPoint;
        const double distance_squared = inner_prod(difference, difference);
        if (distance_squared > radius_squared) {
            radius_squared = distance_squared;
        }
    };

    update_radius(rSurrogatePoint1);
    update_radius(rSurrogatePoint2);
    update_radius(rSkinPoint1);
    update_radius(rSkinPoint2);

    return std::sqrt(radius_squared);
}


// ------------------------------------------------------------------
// Evaluate curve point at parameter t in global coordinates
// ------------------------------------------------------------------
array_1d<double,3> SnakeGapSbmProcess::GlobalPoint(
    const GeometryType& rCurve,
    const double        T)
{
    GeometryType::CoordinatesArrayType local_coordinates;
    for (std::size_t i = 0; i < local_coordinates.size(); ++i) {
        local_coordinates[i] = 0.0;
    }
    local_coordinates[0] = T;

    array_1d<double,3> global_point;
    for (std::size_t i = 0; i < 3; ++i) {
        global_point[i] = 0.0;
    }
    rCurve.GlobalCoordinates(global_point, local_coordinates);
    return global_point;
}

// ------------------------------------------------------------------
// Coons patch mapping X(ξ,η)
// ------------------------------------------------------------------
array_1d<double,3> SnakeGapSbmProcess::CoonsPoint(
    const double                  Xi,
    const double                  Eta,
    const BrepCurveType&          rB0,
    const BrepCurveType&          rL0,
    const BrepCurveType&          rL1,
    const BrepCurveType&          rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11)
{
    const array_1d<double,3> brep_B0_at_xi = GlobalPoint(rB0, Xi);
    const array_1d<double,3> brep_B1_at_xi = GlobalPoint(rB1, Xi);
    const array_1d<double,3> brep_L0_at_eta = GlobalPoint(rL0, Eta);
    const array_1d<double,3> brep_L1_at_eta = GlobalPoint(rL1, Eta);

    const double om_xi  = 1.0 - Xi;
    const double om_eta = 1.0 - Eta;

    return  om_xi * brep_L0_at_eta + Xi * brep_L1_at_eta
          + om_eta * brep_B0_at_xi + Eta * brep_B1_at_xi
          - om_xi * om_eta * rP00
          - Xi    * om_eta * rP10
          - Xi    * Eta    * rP11
          - om_xi * Eta    * rP01;
}

// ------------------------------------------------------------------
// Finite-difference derivative ∂X/∂ξ or ∂X/∂η
// ------------------------------------------------------------------

array_1d<double,3> SnakeGapSbmProcess::CoonsDerivativeFD(
    const double                  Xi,
    const double                  Eta,
    const bool                    WithRespectToXi,
    const BrepCurveType&          rB0,
    const BrepCurveType&          rL0,
    const BrepCurveType&          rL1,
    const BrepCurveType&          rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11,
    const double                  Step)
{
    if (WithRespectToXi) {
        return 0.5/Step *
               ( CoonsPoint(Xi+Step, Eta, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) -
                 CoonsPoint(Xi-Step, Eta, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) );
    } else {
        return 0.5/Step *
               ( CoonsPoint(Xi, Eta+Step, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) -
                 CoonsPoint(Xi, Eta-Step, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) );
    }
}

array_1d<double,3> SnakeGapSbmProcess::CoonsDerivative(
    const double                  Xi,
    const double                  Eta,
    const bool                    WithRespectToXi,
    const BrepCurveType&          rB0,
    const BrepCurveType&          rL0,
    const BrepCurveType&          rL1,
    const BrepCurveType&          rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11)
{

    // Evaluate dX/ds at normalized s (chain rule: dX/ds = dX/dt * dt/ds = dX/dt * (±len))
    auto EvalDerAt = [&](const BrepCurveType& rCurve, double t) -> array_1d<double,3> {
        GeometryType::CoordinatesArrayType lc;
        lc.resize(rCurve.LocalSpaceDimension(), false);
        lc[0] = t;

        std::vector<CoordinatesArrayType> D; // D[0]=X, D[1]=dX/dt
        rCurve.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->GlobalSpaceDerivatives(D, lc, 1);

        array_1d<double,3> dXds = ZeroVector(3);
        dXds[0] = D[1][0];
        dXds[1] = D[1][1];
        dXds[2] = D[1][2];
        return dXds;
    };

    // Boundary values and derivatives (w.r.t. ξ or η in [0,1])

    const array_1d<double,3> brep_B0_at_xi = GlobalPoint(rB0, Xi);
    const array_1d<double,3> brep_B1_at_xi = GlobalPoint(rB1, Xi);
    const array_1d<double,3> brep_L0_at_eta = GlobalPoint(rL0, Eta);
    const array_1d<double,3> brep_L1_at_eta = GlobalPoint(rL1, Eta);

    const array_1d<double,3> dB0_at_xi = EvalDerAt(rB0, Xi);   // dB0_at_xi/dξ (xi)
    const array_1d<double,3> dB1_at_xi = EvalDerAt(rB1, Xi);   // dB1_at_xi/dξ (xi)
    const array_1d<double,3> dL0_at_eta = EvalDerAt(rL0, Eta);  // dL0_at_eta/dη (eta)
    const array_1d<double,3> dL1_at_eta = EvalDerAt(rL1, Eta);  // dL1_at_eta/dη (eta)

    array_1d<double,3> dS = ZeroVector(3);
    if (WithRespectToXi) {
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

//TODO: Extend for Optimal boundary.
SnakeGapSbmProcess::IntegrationPointsArrayType
SnakeGapSbmProcess::CreateCoonsPatchGaussPoints(
    const std::size_t             Order,
    const BrepCurveType&           rB0,
    const BrepCurveType&           rL0,
    const BrepCurveType&           rL1,
    const BrepCurveType&           rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11) const
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
    array_1d<double,3> a = rP10 - rP00;
    array_1d<double,3> b = rP01 - rP00;
    // --- reference normal from Coons derivatives at center --------------
    const double xi_c  = 0.5, eta_c = 0.5;
    const array_1d<double,3> dXi_c  = CoonsDerivative(xi_c, eta_c, true ,
                                                        rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
    const array_1d<double,3> dEta_c = CoonsDerivative(xi_c, eta_c, false,
                                                        rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
    array_1d<double,3> n_ref = MathUtils<double>::CrossProduct(dXi_c, dEta_c);


    // normalize (if non-degenerate)
    const double nrm = norm_2(n_ref);
    if (nrm > 0.0) n_ref /= nrm;


    for (std::size_t i=0; i<Order; ++i)
    for (std::size_t j=0; j<Order; ++j)
    {
        const double xi_i  = xi[i];
        const double eta_j = xi[j];
        double w_ij  = w[i] * w[j];

        // --- Jacobian (signed) ------------------------------------------
        const array_1d<double,3> dXi  = CoonsDerivative(xi_i, eta_j, true ,
                                                          rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
        const array_1d<double,3> dEta = CoonsDerivative(xi_i, eta_j, false,
                                                          rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
        const array_1d<double,3> cr   = MathUtils<double>::CrossProduct(dXi, dEta);

        // Signed Jacobian: projection onto reference normal (can be < 0)
        const double jac_signed = inner_prod(cr, n_ref);

        // --- Global coordinates -----------------------------------------
        const array_1d<double,3> X =
            CoonsPoint(xi_i, eta_j, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);

        
        array_1d<double,3> e1 = rP10 - rP00;
        array_1d<double,3> e2 = X - rP00;


        // Store global coords + signed weight
        gp_list.emplace_back( IntegrationPoint<3>( X[0], X[1], X[2],
                                                   w_ij * jac_signed ) );
    }
    return gp_list;
}

template <bool TIsInnerLoop>
void SnakeGapSbmProcess::SetSurrogateToSkinProjections(
    const ModelPart& rSurrogateSubModelPart,
    const ModelPart& rSkinSubModelPart,
    const KnotSpanIdsCSR& rSkinNodesPerSpan)
{
    const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();
    const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);
    const auto& parameter_space_corners = r_parent_model_part.GetValue(PARAMETER_SPACE_CORNERS);

    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections: KNOT_SPAN_SIZES must contain at least two entries." << std::endl;
    KRATOS_ERROR_IF(parameter_space_corners.size() < 2)
        << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections: PARAMETER_SPACE_CORNERS must contain two vectors." << std::endl;
    KRATOS_ERROR_IF(parameter_space_corners[0].size() < 2 || parameter_space_corners[1].size() < 2)
        << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections: PARAMETER_SPACE_CORNERS entries must contain min and max values." << std::endl;

    const double span_size_x = knot_span_sizes[0];
    const double span_size_y = knot_span_sizes[1];
    const double min_u = parameter_space_corners[0][0];
    const double max_u = parameter_space_corners[0][1];
    const double min_v = parameter_space_corners[1][0];
    const double max_v = parameter_space_corners[1][1];

    const std::size_t span_count_x = rSkinNodesPerSpan.NumberOfSpansX;
    const std::size_t span_count_y = rSkinNodesPerSpan.NumberOfSpansY;

    auto clamp_coordinate = [](double coordinate, double min_value, double max_value) {
        if (coordinate < min_value) {
            return min_value;
        }
        if (coordinate > max_value) {
            return max_value;
        }
        return coordinate;
    };

    constexpr double tol = 1.0e-8;
    auto compute_candidate_indices = [&](double coordinate,
                                         double min_value,
                                         double max_value,
                                         double span_size,
                                         std::size_t span_count) -> std::vector<std::size_t>
    {
        std::vector<std::size_t> indices;
        if (span_count == 0) {
            return indices;
        }

        const double clamped = clamp_coordinate(coordinate, min_value, max_value);
        const double relative = (clamped - min_value) / span_size;

        const double relative_minus = std::max(0.0, relative - tol);
        const double relative_plus = std::min(static_cast<double>(span_count), relative + tol);

        std::size_t index_min = static_cast<std::size_t>(std::floor(relative_minus));
        if (index_min >= span_count) {
            index_min = span_count - 1;
        }

        std::size_t index_max = static_cast<std::size_t>(std::floor(relative_plus));
        if (index_max >= span_count) {
            index_max = span_count - 1;
        }

        indices.push_back(index_min);
        if (index_max != index_min) {
            indices.push_back(index_max);
        }

        return indices;
    };

    auto gather_candidate_node_ids = [&](const Node& rSurrogateNode) {
        std::vector<IndexType> candidate_ids;
        const auto x_indices = compute_candidate_indices(rSurrogateNode.X(), min_u, max_u, span_size_x, span_count_x);
        const auto y_indices = compute_candidate_indices(rSurrogateNode.Y(), min_v, max_v, span_size_y, span_count_y);

        // Helper to append ids of nonzero k
        auto append_ids_by_k = [&](const KnotSpanIdsCSR& S, std::size_t k) {
            const auto view = CellIdsByK(S, k);
            if (view.size == 0) {
                return;
            }

            const IndexType* const pool_begin = S.pool.data();
            const IndexType* const pool_end = pool_begin + S.pool.size();
            const IndexType* const first = view.data;
            const IndexType* const last = view.data + view.size;

            KRATOS_DEBUG_ERROR_IF(first < pool_begin || last > pool_end)
                << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections detected an invalid "
                << "pointer range for CSR entry k = " << k
                << " (pool bounds [" << static_cast<const void*>(pool_begin) << ", "
                << static_cast<const void*>(pool_end) << "), "
                << "requested [" << static_cast<const void*>(first) << ", "
                << static_cast<const void*>(last) << "))." << std::endl;

            candidate_ids.insert(candidate_ids.end(), first, last);
        };

        const auto& S = rSkinNodesPerSpan;                 // KnotSpanIdsCSR
        const auto& A = S.Occupancy;

        // Query only cells (ix,iy)
        for (const std::size_t ix : x_indices) {
            if (ix >= S.NumberOfSpansX) continue;
            for (const std::size_t iy : y_indices) {
                if (iy >= S.NumberOfSpansY) continue;

                const std::size_t k = FindNnzIndex(A, ix, iy);
                if (k != static_cast<std::size_t>(-1)) {
                    append_ids_by_k(S, k);
                }
            }
        }

        // If empty, collect all ids from all nonzeros
        if (candidate_ids.empty()) {
            const std::size_t nnz = static_cast<std::size_t>(A.value_data().size());
            for (std::size_t k = 0; k < nnz; ++k) {
                append_ids_by_k(S, k);
            }
        }

        std::sort(candidate_ids.begin(), candidate_ids.end());
        candidate_ids.erase(std::unique(candidate_ids.begin(), candidate_ids.end()), candidate_ids.end());

        return candidate_ids;
    };

    struct CandidateSelectionResult
    {
        IndexType Id = std::numeric_limits<IndexType>::max();
        double Distance = std::numeric_limits<double>::max();
        bool MatchesForcedLayers = false;
    };

    auto select_candidate = [&](Node::Pointer pSurrogateNode, const std::vector<std::string>& forced_layers) {
        CandidateSelectionResult result;

        const auto candidate_ids = gather_candidate_node_ids(*pSurrogateNode);

        IndexType fallback_id = std::numeric_limits<IndexType>::max();
        double fallback_distance = std::numeric_limits<double>::max();

        auto evaluate_candidate = [&](const Node& r_candidate_node) {
            const IndexType candidate_id = r_candidate_node.Id();
            const auto& candidate_layers = r_candidate_node.GetValue(CONNECTED_LAYERS);

            bool matches_forced = false;
            for (const auto& forced_layer : forced_layers) {
                if (std::find(candidate_layers.begin(), candidate_layers.end(), forced_layer) != candidate_layers.end()) {
                    matches_forced = true;
                    break;
                }
            }

            array_1d<double, 3> difference = r_candidate_node.Coordinates();
            difference -= pSurrogateNode->Coordinates();
            const double distance = norm_2(difference);

            if (distance < fallback_distance) {
                fallback_distance = distance;
                fallback_id = candidate_id;
            }

            if (distance < result.Distance) {
                result.Id = r_candidate_node.Id();
                result.Distance = distance;
                result.MatchesForcedLayers = matches_forced;
            }
        };

        for (const IndexType candidate_id : candidate_ids) {
            evaluate_candidate(rSkinSubModelPart.GetNode(candidate_id));
        }

        auto process_fallback_nodes = [&]() {
            const ModelPart::NodesContainerType* p_nodes = nullptr;
            if (rSkinSubModelPart.HasSubModelPart("interface_vertices")) {
                p_nodes = &rSkinSubModelPart.GetSubModelPart("interface_vertices").Nodes();
            } else {
                p_nodes = &rSkinSubModelPart.Nodes();
            }

            for (const auto& r_node : *p_nodes) {
                evaluate_candidate(r_node);
            }
        };

        if (candidate_ids.empty() || result.Id == std::numeric_limits<IndexType>::max()) {
            process_fallback_nodes();
        }

        if (result.Id == std::numeric_limits<IndexType>::max() && fallback_id != std::numeric_limits<IndexType>::max()) {
            result.Id = fallback_id;
            result.Distance = fallback_distance;
            result.MatchesForcedLayers = forced_layers.empty();
        }
        return result;
    };

    bool is_entering = false;
    auto has_projection = [](Node& rNode){
                return rNode.Has(PROJECTION_NODE_ID);
            };

    for (auto& r_surrogate_condition : rSurrogateSubModelPart.Conditions())
    {
        is_entering = !is_entering;

        Node::Pointer p_surrogate_node_1 = r_surrogate_condition.pGetGeometry()->pGetPoint(0);
        Node::Pointer p_surrogate_node_2 = r_surrogate_condition.pGetGeometry()->pGetPoint(1);

        const bool has_proj_1 = has_projection(*p_surrogate_node_1);
        const bool has_proj_2 = has_projection(*p_surrogate_node_2);

        IndexType skin_node_id_1 = std::numeric_limits<IndexType>::max();
        IndexType skin_node_id_2 = std::numeric_limits<IndexType>::max();

        if (has_proj_1) {
            skin_node_id_1 = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
        }
        if (has_proj_2) {
            skin_node_id_2 = p_surrogate_node_2->GetValue(PROJECTION_NODE_ID);
        }

        if (has_proj_1 && has_proj_2) {
            const auto node_1_connected_layers = p_surrogate_node_1->GetValue(CONNECTED_LAYERS);
            const auto node_2_connected_layers = p_surrogate_node_2->GetValue(CONNECTED_LAYERS);

            bool have_common_layer = false;
            for (const auto& layer1 : node_1_connected_layers) {
                if (std::find(node_2_connected_layers.begin(), node_2_connected_layers.end(), layer1) != node_2_connected_layers.end()) {
                    have_common_layer = true;
                    break;
                }
            }

            if (have_common_layer) {
                continue;
            }

            AssestProjectionsFeasibility(rSkinSubModelPart, p_surrogate_node_1, p_surrogate_node_2);
            continue;
        }

        std::vector<std::string> forced_layers;

        if (has_proj_1 && !has_proj_2) {
            forced_layers = p_surrogate_node_1->GetValue(CONNECTED_LAYERS);
        } else if (!has_proj_1 && has_proj_2) {
            forced_layers = p_surrogate_node_2->GetValue(CONNECTED_LAYERS);
        }

        if (!has_proj_1) {
            const auto selection = select_candidate(p_surrogate_node_1, forced_layers);
            KRATOS_ERROR_IF(selection.Id == std::numeric_limits<IndexType>::max())
                << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections: no skin node candidate found for surrogate node "
                << p_surrogate_node_1->Id() << " at " << p_surrogate_node_1->Coordinates() << std::endl;

            skin_node_id_1 = selection.Id;

            const auto& connected_layers = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_LAYERS);
            const auto& connected_conditions = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_CONDITIONS);

            p_surrogate_node_1->SetValue(PROJECTION_NODE_ID, skin_node_id_1);
            p_surrogate_node_1->SetValue(CONNECTED_LAYERS, connected_layers);
            p_surrogate_node_1->SetValue(CONNECTED_CONDITIONS, connected_conditions);

            if (!selection.MatchesForcedLayers && !forced_layers.empty()) {

                AssestProjectionsFeasibility(rSkinSubModelPart, p_surrogate_node_1, p_surrogate_node_2);
            }

            forced_layers = connected_layers;
        }

        if (!has_proj_2) {
            const auto selection = select_candidate(p_surrogate_node_2, forced_layers);
            KRATOS_ERROR_IF(selection.Id == std::numeric_limits<IndexType>::max())
                << "::[SnakeGapSbmProcess]::SetSurrogateToSkinProjections: no skin node candidate found for surrogate node "
                << p_surrogate_node_2->Id() << " at " << p_surrogate_node_2->Coordinates() << std::endl;

            skin_node_id_2 = selection.Id;

            const auto& connected_layers = rSkinSubModelPart.GetNode(skin_node_id_2).GetValue(CONNECTED_LAYERS);
            const auto& connected_conditions = rSkinSubModelPart.GetNode(skin_node_id_2).GetValue(CONNECTED_CONDITIONS);

            p_surrogate_node_2->SetValue(PROJECTION_NODE_ID, skin_node_id_2);
            p_surrogate_node_2->SetValue(CONNECTED_LAYERS, connected_layers);
            p_surrogate_node_2->SetValue(CONNECTED_CONDITIONS, connected_conditions);

            if (!selection.MatchesForcedLayers && !forced_layers.empty()) {
                AssestProjectionsFeasibility(rSkinSubModelPart, p_surrogate_node_1, p_surrogate_node_2);
            }
        }
    }

    // correct the projections
    std::size_t count_intersections = 1;
    IndexType iter_check = 0;
    while (count_intersections > 0 && iter_check < 5)
    {
        iter_check +=1;
        count_intersections = 0;
        for (auto& r_surrogate_condition : rSurrogateSubModelPart.Conditions())
        {
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
    }

    KRATOS_ERROR_IF(iter_check == 5) << "::[SnakeSbmProcess]:: Maximum iteration reached when checking intersections between projections. Please check the input data." << std::endl;

}

void SnakeGapSbmProcess::AssestProjectionsFeasibility(
    const ModelPart& rSkinSubModelPart,
    Node::Pointer pSurrogateNode1, 
    Node::Pointer pSurrogateNode2)
{
    // TODO: manage the case where there is just one element between two separeted layers and no projections onto the middle layer
    // i.e. when the first surrogate node is projected onto layer A and the second surrogate node onto layer C.
    // -> need a vertex correction on both sides, not just one
    std::vector<std::string> forced_layers = pSurrogateNode1->GetValue(CONNECTED_LAYERS);
    std::vector<std::string> connected_layers_surrogate_node_2 = pSurrogateNode2->GetValue(CONNECTED_LAYERS);

    const ModelPart::NodesContainerType* p_interface_nodes = nullptr;
    if (rSkinSubModelPart.HasSubModelPart("interface_vertices")) {
        p_interface_nodes = &rSkinSubModelPart.GetSubModelPart("interface_vertices").Nodes();
    } else {
        p_interface_nodes = &rSkinSubModelPart.Nodes();
    }

    IndexType nearest_node_id = std::numeric_limits<IndexType>::max();
    double minimum_distance = std::numeric_limits<double>::max();

    for (const auto& r_interface_node : *p_interface_nodes) {
        const auto& candidate_layers = r_interface_node.GetValue(CONNECTED_LAYERS);
        bool matches_layers = false;
        for (const auto& forced_layer : forced_layers) {
            if (std::find(candidate_layers.begin(), candidate_layers.end(), forced_layer) != candidate_layers.end()) {
                matches_layers = true;
                break;
            }
        }
        if (!matches_layers && !forced_layers.empty()) {
            continue;
        }

        array_1d<double, 3> difference = r_interface_node.Coordinates();
        difference -= pSurrogateNode2->Coordinates();
        const double distance = norm_2(difference);
        if (distance < minimum_distance) {
            minimum_distance = distance;
            nearest_node_id = r_interface_node.Id();
        }
    }

    KRATOS_ERROR_IF(nearest_node_id == std::numeric_limits<IndexType>::max())
        << "::[SnakeGapSbmProcess]::AssestProjectionsFeasibility: no interface node satisfies the forced layers for surrogate node "
        << pSurrogateNode2->Id() << std::endl;

    const auto& r_projection_node = rSkinSubModelPart.GetNode(nearest_node_id);
    const auto& connected_layers = rSkinSubModelPart.GetNode(nearest_node_id).GetValue(CONNECTED_LAYERS);
    const auto& connected_conditions = rSkinSubModelPart.GetNode(nearest_node_id).GetValue(CONNECTED_CONDITIONS);

    // check if the found node matches also the layer of the other surrogate node projection
    bool matches_layers = false;
    for (const auto& layer : connected_layers_surrogate_node_2) {
        if (std::find(connected_layers.begin(), connected_layers.end(), layer) != connected_layers.end()) {
            matches_layers = true;
            break;
        }
    }

    KRATOS_ERROR_IF(!matches_layers)
        << "::[SnakeGapSbmProcess]::AssestProjectionsFeasibility: It wasn't possible to find a common layer for the nodes: " 
        << pSurrogateNode1->Coordinates() << ", " << pSurrogateNode1->Coordinates() 
        << connected_layers << connected_layers_surrogate_node_2
        << ".\n Just one element in the defined layer. Refine more the mesh" << std::endl;

    if (norm_2(r_projection_node.Coordinates()-pSurrogateNode1->Coordinates()) <
        norm_2(r_projection_node.Coordinates()-pSurrogateNode2->Coordinates()))
    {
        pSurrogateNode1->SetValue(PROJECTION_NODE_ID, nearest_node_id);
        pSurrogateNode1->SetValue(CONNECTED_LAYERS, connected_layers);
        pSurrogateNode1->SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }
    else
    {
        pSurrogateNode2->SetValue(PROJECTION_NODE_ID, nearest_node_id);
        pSurrogateNode2->SetValue(CONNECTED_LAYERS, connected_layers);
        pSurrogateNode2->SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }
    return;
}

template <bool TIsInnerLoop>
IndexType SnakeGapSbmProcess::FindClosestNodeInLayerWithDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection)
{
    KRATOS_ERROR_IF(rSkinConditionsPerSpan.NumberOfSpansX == 0 ||
                    rSkinConditionsPerSpan.NumberOfSpansY == 0)
        << "::[SnakeGapSbmProcess]::FindClosestNodeInLayerWithDirection: condition span matrix is empty." << std::endl;

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

    const double span_size_x = (rSkinConditionsPerSpan.SpanSizeX > 0.0) ? rSkinConditionsPerSpan.SpanSizeX
                           : (rKnotSpanSizes.size() > 0 ? rKnotSpanSizes[0] : 1.0);
    const double span_size_y = (rSkinConditionsPerSpan.SpanSizeY > 0.0) ? rSkinConditionsPerSpan.SpanSizeY
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

    const double min_u = rSkinConditionsPerSpan.MinU;
    const double max_u = rSkinConditionsPerSpan.MaxU;
    const double min_v = rSkinConditionsPerSpan.MinV;
    const double max_v = rSkinConditionsPerSpan.MaxV;

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

    IndexType best_node_id = std::numeric_limits<IndexType>::max();
    double best_intersection_distance = std::numeric_limits<double>::max();
    double best_intersection_node_distance = std::numeric_limits<double>::max();

    auto update_candidate = [&](Node::Pointer p_candidate_node, const CoordinatesArrayType& rIntersectionPoint){
        const auto& candidate_layers = p_candidate_node->GetValue(CONNECTED_LAYERS);
        if (std::find(candidate_layers.begin(), candidate_layers.end(), rLayer) == candidate_layers.end()) {
            return;
        }
        array_1d<double,3> diff_intersection = rIntersectionPoint;
        diff_intersection -= rStartPoint;
        const double distance_intersection = norm_2(diff_intersection);

        if (distance_intersection <= best_intersection_distance) {
            best_intersection_distance = distance_intersection;

            array_1d<double,3> diff_intersection_node = p_candidate_node->Coordinates();
            diff_intersection_node -= rIntersectionPoint;
            const double distance_intersection_node = norm_2(diff_intersection_node);
            if (distance_intersection_node < best_intersection_node_distance)
            {
                best_intersection_node_distance = distance_intersection_node;
                best_node_id = p_candidate_node->Id();
            }
            
        }
    };

    const std::size_t number_spans_x = rSkinConditionsPerSpan.NumberOfSpansX;
    const std::size_t number_spans_y = rSkinConditionsPerSpan.NumberOfSpansY;
    const std::size_t base_ix = compute_index(rStartPoint[0], min_u, max_u, span_size_x, number_spans_x);
    const std::size_t base_iy = compute_index(rStartPoint[1], min_v, max_v, span_size_y, number_spans_y);

    const int max_search_level = static_cast<int>(std::max(number_spans_x, number_spans_y));

    std::vector<IndexType> candidate_conditions;
    candidate_conditions.reserve(32);

    bool found_intersection = false;
    for (int level = 0; level < max_search_level && !found_intersection; ++level) {

        const std::size_t extension = static_cast<std::size_t>(level+1);
        const double current_length = reference_span_size * static_cast<double>(extension);

        array_1d<double,3> segment_start = rStartPoint - direction * reference_span_size/2;
        array_1d<double,3> segment_end = rStartPoint;
        segment_end += direction * current_length;

        const std::size_t min_ix = (base_ix > extension) ? base_ix - extension : 0;
        const std::size_t max_ix = std::min<std::size_t>(base_ix + extension, number_spans_x > 0 ? number_spans_x - 1 : 0);
        const std::size_t min_iy = (base_iy > extension) ? base_iy - extension : 0;
        const std::size_t max_iy = std::min<std::size_t>(base_iy + extension, number_spans_y > 0 ? number_spans_y - 1 : 0);

        candidate_conditions.clear();
        // Collect unique condition ids from CSR cells in [min_ix..max_ix] × [min_iy..max_iy]
        const auto& S = rSkinConditionsPerSpan;
        std::unordered_set<IndexType> seen;
        seen.reserve(256); // optional

        for (std::size_t ix = min_ix; ix <= max_ix; ++ix) {
            if (ix >= S.NumberOfSpansX) continue;
            for (std::size_t iy = min_iy; iy <= max_iy; ++iy) {
                if (iy >= S.NumberOfSpansY) continue;

                const std::size_t k = FindNnzIndex(S.Occupancy, ix, iy);
                if (k == static_cast<std::size_t>(-1)) continue;

                const auto ids = CellIdsByK(S, k);
                for (std::size_t t = 0; t < ids.size; ++t) {
                    const IndexType cid = ids.data[t];
                    if (seen.insert(cid).second) {
                        candidate_conditions.push_back(cid);
                    }
                }
            }
        }


        if (candidate_conditions.empty()) {
            continue;
        }

        // Example: framework factory
        Node::Pointer start_node = Node::Pointer(new Node(0, segment_start[0], segment_start[1], segment_start[2]));
        Node::Pointer end_node   = Node::Pointer(new Node(0, segment_end[0],   segment_end[1],   segment_end[2]));

        for (const IndexType condition_id : candidate_conditions) {
            const auto r_condition = rSkinSubModelPart.pGetCondition(condition_id);
            const auto p_geometry = r_condition->pGetGeometry();
            if (p_geometry->size() < 2) {
                continue;
            }

            const auto& node_a = p_geometry->pGetPoint(0);
            const auto& node_b = p_geometry->pGetPoint(1);

            CoordinatesArrayType intersection_point;
        
            const bool intersects = SegmentsIntersect(start_node, end_node, node_a, node_b, intersection_point);

            if (intersects) {
                best_intersection_node_distance = std::numeric_limits<double>::max();
                update_candidate(node_a, intersection_point);
                update_candidate(node_b, intersection_point);
                found_intersection = true;
            }
        }

        // // Fallback: try the opposite direction if nothing found at this level
        // if (!found_intersection) {
        //     array_1d<double,3> segment_end_neg = rStartPoint - direction * current_length;
        //     Node::Pointer end_node_neg = Node::Pointer(new Node(0, segment_end_neg[0], segment_end_neg[1], segment_end_neg[2]));
        //     for (const IndexType condition_id : candidate_conditions) {
        //         const auto r_condition = rSkinSubModelPart.pGetCondition(condition_id);
        //         const auto p_geometry = r_condition->pGetGeometry();
        //         if (p_geometry->size() < 2) {
        //             continue;
        //         }

        //         const auto& node_a = p_geometry->pGetPoint(0);
        //         const auto& node_b = p_geometry->pGetPoint(1);

        //         CoordinatesArrayType intersection_point;

        //         // bool intersects = SegmentsIntersectRay(start_node, end_node_neg, node_a, node_b, intersection_point);
        //         bool intersects = SegmentsIntersect(start_node, end_node_neg, node_a, node_b, intersection_point);

        //         if (intersects) {
        //             update_candidate(node_a, intersection_point);
        //             update_candidate(node_b, intersection_point);
        //             found_intersection = true;
        //         }
        //     }
        // }

        if (best_node_id != std::numeric_limits<IndexType>::max()) {
            return best_node_id;
        }
    }

    KRATOS_ERROR << "::[SnakeGapSbmProcess]::FindClosestNodeInLayerWithDirection: no node found for layer "
                 << rLayer << " starting from point " << rStartPoint
                 << ". Skin sub model part '" << rSkinSubModelPart.Name() << "' has "
                 << rSkinSubModelPart.NumberOfNodes() << " nodes and "
                 << rSkinSubModelPart.NumberOfConditions() << " conditions."
                 << " Search parameters: base_ix=" << base_ix << ", base_iy=" << base_iy
                 << ", span_size_x=" << span_size_x << ", span_size_y=" << span_size_y
                 << ", reference_span_size=" << reference_span_size
                 << ", direction=" << direction << std::endl;

    return 0;
}

template <bool TIsInnerLoop>
std::pair<SnakeGapSbmProcess::IndexType, SnakeGapSbmProcess::IndexType> SnakeGapSbmProcess::FindClosestPairInLayerWithNormalDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection)
{
    KRATOS_ERROR_IF(rSkinConditionsPerSpan.NumberOfSpansX == 0 ||
                    rSkinConditionsPerSpan.NumberOfSpansY == 0)
        << "::[SnakeGapSbmProcess]::FindClosestNodeInLayerWithDirection: condition span matrix is empty." << std::endl;

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

    const double span_size_x = (rSkinConditionsPerSpan.SpanSizeX > 0.0) ? rSkinConditionsPerSpan.SpanSizeX
                           : (rKnotSpanSizes.size() > 0 ? rKnotSpanSizes[0] : 1.0);
    const double span_size_y = (rSkinConditionsPerSpan.SpanSizeY > 0.0) ? rSkinConditionsPerSpan.SpanSizeY
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

    const double min_u = rSkinConditionsPerSpan.MinU;
    const double max_u = rSkinConditionsPerSpan.MaxU;
    const double min_v = rSkinConditionsPerSpan.MinV;
    const double max_v = rSkinConditionsPerSpan.MaxV;

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
    auto update_candidate = [&](Node::Pointer p_node_a, Node::Pointer p_node_b, const CoordinatesArrayType& rIntersectionPoint){
        const auto& candidate_layers_a = p_node_a->GetValue(CONNECTED_LAYERS);
        const auto& candidate_layers_b = p_node_b->GetValue(CONNECTED_LAYERS);
        if ((std::find(candidate_layers_a.begin(), candidate_layers_a.end(), rLayer) == candidate_layers_a.end()) &&
            (std::find(candidate_layers_b.begin(), candidate_layers_b.end(), rLayer) == candidate_layers_b.end())) {
            return;
        }
        array_1d<double,3> diff_intersection = rIntersectionPoint;
        diff_intersection -= rStartPoint;
        const double distance_intersection = norm_2(diff_intersection);
        if (distance_intersection < best_intersection_distance) {
            best_intersection_distance = distance_intersection;
            best_pair = {p_node_a->Id(), p_node_b->Id()};
        }
    };

    const std::size_t number_spans_x = rSkinConditionsPerSpan.NumberOfSpansX;
    const std::size_t number_spans_y = rSkinConditionsPerSpan.NumberOfSpansY;
    const std::size_t base_ix = compute_index(rStartPoint[0], min_u, max_u, span_size_x, number_spans_x);
    const std::size_t base_iy = compute_index(rStartPoint[1], min_v, max_v, span_size_y, number_spans_y);

    const int max_search_level = static_cast<int>(std::max(number_spans_x, number_spans_y));

    std::vector<IndexType> candidate_conditions;
    candidate_conditions.reserve(32);

    bool found_intersection = false;
    for (int level = 0; level < max_search_level && !found_intersection; ++level) {
        const std::size_t extension = static_cast<std::size_t>(level+1);
        const double current_length = reference_span_size * static_cast<double>(extension);

        // Cast a forward segment (ray-like) starting at rStartPoint along +direction
        array_1d<double,3> segment_start = rStartPoint - direction * reference_span_size/2;
        array_1d<double,3> segment_end = rStartPoint;
        segment_end += direction * current_length;

        const std::size_t min_ix = (base_ix > extension) ? base_ix - extension : 0;
        const std::size_t max_ix = std::min<std::size_t>(base_ix + extension, number_spans_x > 0 ? number_spans_x - 1 : 0);
        const std::size_t min_iy = (base_iy > extension) ? base_iy - extension : 0;
        const std::size_t max_iy = std::min<std::size_t>(base_iy + extension, number_spans_y > 0 ? number_spans_y - 1 : 0);

        candidate_conditions.clear();
        // Collect unique condition ids from CSR cells in [min_ix..max_ix] × [min_iy..max_iy]
        const auto& S = rSkinConditionsPerSpan;
        std::unordered_set<IndexType> seen;
        seen.reserve(256); // optional

        for (std::size_t ix = min_ix; ix <= max_ix; ++ix) {
            if (ix >= S.NumberOfSpansX) continue;
            for (std::size_t iy = min_iy; iy <= max_iy; ++iy) {
                if (iy >= S.NumberOfSpansY) continue;

                const std::size_t k = FindNnzIndex(S.Occupancy, ix, iy);
                if (k == static_cast<std::size_t>(-1)) continue;

                const auto ids = CellIdsByK(S, k);
                for (std::size_t t = 0; t < ids.size; ++t) {
                    const IndexType cid = ids.data[t];
                    if (seen.insert(cid).second) {
                        candidate_conditions.push_back(cid);
                    }
                }
            }
        }


        if (candidate_conditions.empty()) {
            continue;
        }

        // Example: framework factory
        Node::Pointer start_node = Node::Pointer(new Node(0, segment_start[0], segment_start[1], segment_start[2]));
        Node::Pointer end_node   = Node::Pointer(new Node(0, segment_end[0],   segment_end[1],   segment_end[2]));

        for (const IndexType condition_id : candidate_conditions) {
            const auto r_condition = rSkinSubModelPart.pGetCondition(condition_id);
            const auto p_geometry = r_condition->pGetGeometry();
            if (p_geometry->size() < 2) {
                continue;
            }

            const auto& node_a = p_geometry->pGetPoint(0);
            const auto& node_b = p_geometry->pGetPoint(1);

            CoordinatesArrayType intersection_point; 
            // bool intersects = SegmentsIntersectRay(start_node, end_node, node_a, node_b, intersection_point);
            bool intersects = SegmentsIntersect(start_node, end_node, node_a, node_b, intersection_point);

            if (intersects) {
                update_candidate(node_a, node_b, intersection_point);
                found_intersection = true;
            }
        }

        if (best_pair.first != std::numeric_limits<IndexType>::max()) {
            return best_pair;
        }
    }

    KRATOS_ERROR << "::[SnakeGapSbmProcess]::FindClosestPairInLayerWithNormalDirection: no intersection found for layer "
                 << rLayer << " starting from point " << rStartPoint
                 << ". Skin sub model part '" << rSkinSubModelPart.Name() << "' has "
                 << rSkinSubModelPart.NumberOfNodes() << " nodes and "
                 << rSkinSubModelPart.NumberOfConditions() << " conditions."
                 << " Search parameters: base_ix=" << base_ix << ", base_iy=" << base_iy
                 << ", span_size_x=" << span_size_x << ", span_size_y=" << span_size_y
                 << ", reference_span_size=" << reference_span_size
                 << ", direction=" << direction << std::endl;

    return best_pair;
}

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

    // Attach the surrogate middle geometry to the two projection (skin) nodes
    auto& r_neighbour_geometries_1 = pProjectionNode1->GetValue(NEIGHBOUR_GEOMETRIES);
    r_neighbour_geometries_1.push_back(pSurrogateMiddleGeometry);

    auto& r_neighbour_geometries_2 = pProjectionNode2->GetValue(NEIGHBOUR_GEOMETRIES);
    r_neighbour_geometries_2.push_back(pSurrogateMiddleGeometry);
}

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

void SnakeGapSbmProcess::CreateInnerSkinMultipatchCouplingConditions(
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    NurbsSurfaceType::Pointer& pNurbsSurface)
{
    // Ensure/obtain MultipatchCouplingConditions submodelpart and add conditions there
    ModelPart& r_multipatch_coupling = mpIgaModelPart->HasSubModelPart("MultipatchCouplingConditions")
        ? mpIgaModelPart->GetSubModelPart("MultipatchCouplingConditions")
        : mpIgaModelPart->CreateSubModelPart("MultipatchCouplingConditions");

    // Validate inner skin conditions share a common neighbour geometry
    for (const auto& r_condition : rSkinSubModelPart.Conditions()) {
        const auto& r_geom = r_condition.GetGeometry();
        KRATOS_ERROR_IF(r_geom.PointsNumber() < 2)
            << "::[SnakeGapSbmProcess]:: Skin condition #" << r_condition.Id()
            << " has less than 2 nodes." << std::endl;

        const auto& r_node_0 = r_geom[0];
        const auto& r_node_1 = r_geom[1];

        const double characteristich_length = norm_2(r_node_0-r_node_1)/2;

        const auto& neigh_0 = r_node_0.GetValue(NEIGHBOUR_GEOMETRIES);
        const auto& neigh_1 = r_node_1.GetValue(NEIGHBOUR_GEOMETRIES);

        // If either node has no neighbour geometries, this typically indicates
        // the skin coincides with the surrogate (no unique projection target).
        KRATOS_ERROR_IF(neigh_0.empty() || neigh_1.empty())
            << "::[SnakeGapSbmProcess]:: Empty NEIGHBOUR_GEOMETRIES detected on skin condition #" << r_condition.Id()
            << ". Node #" << r_node_0.Id() << " neigh_size=" << neigh_0.size()
            << ", Node #" << r_node_1.Id() << " neigh_size=" << neigh_1.size()
            << ". This usually happens when skin and surrogate coincide; adjust layers or skip coupling for this condition." << std::endl;

        bool has_common = false;
        Geometry<Node>::Pointer p_common_geometry;
        // Among common neighbour geometries, pick the one closest to the skin condition center
        const auto skin_center = r_geom.Center();
        double best_dist = std::numeric_limits<double>::max();
        std::unordered_set<const void*> visited;
        for (const auto& p_g0 : neigh_0) {
            const void* key = p_g0.get();
            if (!visited.insert(key).second) continue; // skip duplicates from neigh_0
            bool is_common = false;
            for (const auto& p_g1 : neigh_1) {
                if (p_g1.get() == key) { is_common = true; break; }
            }
            if (!is_common) continue;

            const double dist = norm_2(skin_center - p_g0->Center());
            if (dist < best_dist) {
                best_dist = dist;
                has_common = true;
                p_common_geometry = p_g0;
            }
        }

        KRATOS_ERROR_IF_NOT(has_common)
            << "::[SnakeGapSbmProcess]:: No common NEIGHBOUR_GEOMETRIES between nodes of skin condition #"
            << r_condition.Id() << std::endl;

        int brep_id = 0;
        // Retrieve BREP_ID of the skin condition identifying the curve on surface
        if (r_condition.Has(BREP_ID)) {
            brep_id = r_condition.GetValue(BREP_ID);
        } else {
            KRATOS_ERROR << "SnakeGapSbmProcess :: "
                << "Skin condition #" << r_condition.Id()
                << " has no BREP_ID set" << std::endl;
        }

        // Determine condition type from the NURBS curve CONDITION_NAME set in MultipatchModeler
        std::string condition_type_name = "";
        if (r_condition.Has(CONDITION_NAME)) {
            condition_type_name = r_condition.GetValue(CONDITION_NAME);
        } else {
            KRATOS_ERROR << "SnakeGapSbmProcess :: "
                << "Skin condition #" << r_condition.Id()
                << " has no CONDITION_NAME set" << std::endl;
        }

        // Create conditions on the skin from the BREP curve with id = brep_id
        auto& r_brep_model_part = mpIgaModelPart->GetParentModelPart();
        auto p_brep_geometry = r_brep_model_part.pGetGeometry(static_cast<IndexType>(brep_id));
        auto p_brep_curve_on_surface = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);
        KRATOS_ERROR_IF(!p_brep_curve_on_surface)
            << "SnakeGapSbmProcess :: Geometry #" << brep_id << " is not a BrepCurveOnSurfaceType." << std::endl;

        IntegrationInfo brep_integration_info = p_brep_curve_on_surface->GetDefaultIntegrationInfo();
        // Choose GP count as 2*max_degree+1, where max_degree is the max between
        // the current patch surface degree and the degree of the surface that hosts this BREP
        std::size_t surface_deg_u = 1, surface_deg_v = 1;
        if (pNurbsSurface) {
            surface_deg_u = pNurbsSurface->PolynomialDegree(0);
            surface_deg_v = pNurbsSurface->PolynomialDegree(1);
        }
        std::size_t host_deg_u = 0, host_deg_v = 0;
        // Retrieve underlying surface from the curve-on-surface
        if (auto p_curve_on_surface = p_brep_curve_on_surface->pGetCurveOnSurface()) {
            auto p_host_surface_geom = p_curve_on_surface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
            if (p_host_surface_geom) {
                if (auto p_host_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(p_host_surface_geom)) {
                    host_deg_u = p_host_surface->PolynomialDegree(0);
                    host_deg_v = p_host_surface->PolynomialDegree(1);
                }
            }
        }
        const std::size_t max_deg = std::max({surface_deg_u, surface_deg_v, host_deg_u, host_deg_v});
        const int points_per_span = static_cast<int>(2 * max_deg + 1);
        brep_integration_info.SetNumberOfIntegrationPointsPerSpan(0, points_per_span);

        IntegrationPointsArrayType brep_integration_points_list;
        GeometriesArrayType brep_quadrature_point_list;

        p_brep_curve_on_surface->CreateIntegrationPoints(
            brep_integration_points_list, brep_integration_info);

        p_brep_curve_on_surface->CreateQuadraturePointGeometries(
            brep_quadrature_point_list,
            /*num shape derivs*/ points_per_span,
            brep_integration_points_list,
            brep_integration_info);

        // Common neighbour geometry to attach
        std::vector<Geometry<Node>::Pointer> surrogate_refs;
        if (p_common_geometry) surrogate_refs.push_back(p_common_geometry);

        // Create the conditions with neighbour geometries = {p_common_geometry}
        std::size_t id = 1;
        if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
            id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

        this->CreateConditions(
            brep_quadrature_point_list.ptr_begin(), brep_quadrature_point_list.ptr_end(),
            r_multipatch_coupling, condition_type_name, id, PropertiesPointerType(),
            rKnotSpanSizes, surrogate_refs, characteristich_length);
    }
}


const Parameters SnakeGapSbmProcess::GetDefaultParameters() const
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
        "number_of_knot_spans" : [10, 10]
    })");
}

const Parameters SnakeGapSbmProcess::GetValidParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 0,
        "lower_point_xyz": [-0.5, -0.5,0.0],
        "upper_point_xyz": [0.5,0.5,0.0],
        "lower_point_uvw": [-0.5,-0.5,0.0],
        "upper_point_uvw": [0.5, 0.5,0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [7, 7],
        "number_of_inner_loops": 0,
        "number_initial_points_if_importing_nurbs": 1000,
        "number_internal_divisions": 0,
        "gap_sbm_type": "default",
        "lambda_inner" : 0.0,
        "lambda_outer" : 1.0,
        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",    
        "skin_model_part_inner_initial_name": "initial_skin_model_part_in",           
        "skin_model_part_name": "skin_model_part",
        "gap_element_name": "CutSbmSolidElement",
        "gap_interface_condition_name": "CutSbmSolidInterfaceCondition"
    })");
}

}  // namespace Kratos.
