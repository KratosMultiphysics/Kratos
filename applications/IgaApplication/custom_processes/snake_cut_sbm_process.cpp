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

// External includes

// Project includes
#include <algorithm>
#include <cmath>
#include <limits>
#include "snake_cut_sbm_process.h"
#include "iga_application_variables.h"

namespace Kratos
{

SnakeCutSbmProcess::SnakeCutSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    SnakeSbmProcess(rModel, ThisParameters)
{

    KRATOS_ERROR_IF_NOT(ThisParameters.Has("cut_element_name")) << "::[SnakeCutSbmProcess]::" 
                    << "Missing \"cut_element_name\" section." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("cut_interface_condition_name")) << "::[SnakeCutSbmProcess]::" 
                    << "Missing \"cut_interface_condition_name\" section." << std::endl;
    
    mpCutElementsSubModelPart = &(mpIgaModelPart->CreateSubModelPart("CutElements"));
    mpCutInterfaceSubModelPart = &(mpIgaModelPart->CreateSubModelPart("CutInterfaces"));
    mCutElementName = ThisParameters["cut_element_name"].GetString();
    mCutInterfaceConditionName = ThisParameters["cut_interface_condition_name"].GetString();
    mLambdaInner = 0.0;
    mLambdaOuter = 1.0;
    mInternalDivision = ThisParameters["number_internal_divisions"].GetInt();
}

SnakeCutSbmProcess::KnotSpanIdsCSR
SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix(
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart) const
{
    KnotSpanIdsCSR out;

    // --- Read parameters from parent model part (as in your code) ---
    const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();

    const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);
    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] KNOT_SPAN_SIZES must have at least two entries.\n";

    const auto& parameter_space_corners = r_parent_model_part.GetValue(PARAMETER_SPACE_CORNERS);
    KRATOS_ERROR_IF(parameter_space_corners.size() < 2)
        << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] PARAMETER_SPACE_CORNERS must have at least two vectors.\n";
    KRATOS_ERROR_IF(parameter_space_corners[0].size() < 2 || parameter_space_corners[1].size() < 2)
        << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] PARAMETER_SPACE_CORNERS vectors must contain min and max values.\n";

    const double span_size_x = knot_span_sizes[0];
    const double span_size_y = knot_span_sizes[1];
    KRATOS_ERROR_IF(span_size_x <= 0.0 || span_size_y <= 0.0)
        << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Knot span sizes must be positive.\n";

    const double min_u = parameter_space_corners[0][0];
    const double max_u = parameter_space_corners[0][1];
    const double min_v = parameter_space_corners[1][0];
    const double max_v = parameter_space_corners[1][1];

    const double domain_length_u = max_u - min_u;
    const double domain_length_v = max_v - min_v;
    KRATOS_ERROR_IF(domain_length_u <= 0.0 || domain_length_v <= 0.0)
        << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Invalid parameter space extents.\n";

    const double spans_x_real = domain_length_u / span_size_x;
    const double spans_y_real = domain_length_v / span_size_y;

    auto compute_span_count = [](double spans_real) -> SizeType {
        constexpr double tolerance = 1.0e-10;
        SizeType span_count = static_cast<SizeType>(std::round(spans_real));
        KRATOS_ERROR_IF(span_count <= 0)
            << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Non-positive number of knot spans.\n";
        KRATOS_ERROR_IF(std::abs(spans_real - static_cast<double>(span_count)) > tolerance)
            << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] Non-integer number of knot spans (" << spans_real << ").\n";
        return span_count;
    };

    const SizeType Nx = compute_span_count(spans_x_real);
    const SizeType Ny = compute_span_count(spans_y_real);

    // --- Fill metadata ---
    out.NumberOfSpansX = Nx;
    out.NumberOfSpansY = Ny;
    out.MinU = min_u; out.MaxU = max_u;
    out.MinV = min_v; out.MaxV = max_v;
    out.SpanSizeX = span_size_x; out.SpanSizeY = span_size_y;

    // Early exit: no nodes
    if (rSkinSubModelPart.NumberOfNodes() == 0) {
        out.Occupancy.resize(Nx, Ny, false);
        return out;
    }

    const double tol = 1.0e-12;
    const double max_u_with_tol = max_u + tol;
    const double max_v_with_tol = max_v + tol;
    const double min_u_with_tol = min_u - tol;
    const double min_v_with_tol = min_v - tol;

    auto compute_index = [tol](double coordinate, double min_value, double max_value, double span_size, SizeType span_count) {
        double clamped = coordinate;
        if (coordinate < min_value) {
            KRATOS_ERROR_IF(coordinate < min_value - tol)
                << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] coordinate below minimum parameter range.\n";
            clamped = min_value;
        } else if (coordinate > max_value) {
            KRATOS_ERROR_IF(coordinate > max_value + tol)
                << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] coordinate above maximum parameter range.\n";
            clamped = max_value;
        }
        double relative = (clamped - min_value) / span_size;
        SizeType idx = static_cast<SizeType>(std::floor(relative + tol));
        if (idx >= span_count) idx = span_count - 1;
        return idx;
    };

    // --- Pass 1: gather unique columns per row for sparsity pattern ---
    std::vector<std::vector<SizeType>> cols_per_row(Nx);

    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const double u = r_node.X();
        const double v = r_node.Y();

        KRATOS_ERROR_IF(u < min_u_with_tol || u > max_u_with_tol)
            << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] node " << r_node.Id()
            << " has u-parameter " << u << " outside [" << min_u << ", " << max_u << "].\n";
        KRATOS_ERROR_IF(v < min_v_with_tol || v > max_v_with_tol)
            << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] node " << r_node.Id()
            << " has v-parameter " << v << " outside [" << min_v << ", " << max_v << "].\n";

        const SizeType ix = compute_index(u, min_u, max_u, span_size_x, Nx);
        const SizeType iy = compute_index(v, min_v, max_v, span_size_y, Ny);
        cols_per_row[ix].push_back(iy);
    }

    SizeType nnz = 0;
    for (auto& cols : cols_per_row) {
        std::sort(cols.begin(), cols.end());
        cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        nnz += static_cast<SizeType>(cols.size());
    }

    // --- Materialize CSR matrix with zero values and build nnz slots ---
    auto& A = out.Occupancy;
    A.resize(Nx, Ny, false);
    A.reserve(nnz);

    for (SizeType i = 0; i < Nx; ++i) {
        for (const SizeType j : cols_per_row[i]) {
            A.push_back(i, j, 0.0);  // keep strictly increasing (i,j)
        }
    }

    // Temporary per-nnz buckets for node ids
    std::vector<std::vector<IndexType>> tmp(nnz);

    // --- Pass 2: fill tmp lists and counts in A.value_data() ---
    for (const auto& r_node : rSkinSubModelPart.Nodes()) {
        const SizeType ix = compute_index(r_node.X(), min_u, max_u, span_size_x, Nx);
        const SizeType iy = compute_index(r_node.Y(), min_v, max_v, span_size_y, Ny);

        const SizeType k = FindNnzIndex(A, ix, iy);
        KRATOS_DEBUG_ERROR_IF(k == static_cast<SizeType>(-1))
            << "[SnakeCutSbmProcess::CreateSkinNodesPerKnotSpanMatrix] nonzero (ix,iy) not found in CSR pattern.\n";

        tmp[k].push_back(r_node.Id());
        A.value_data()[k] += 1.0;  // optional: store count per cell
    }

    // --- Pack payload to compact pool ---
    CommitPayload(tmp, out);

    return out;
}

// === Builder: conditions -> KnotSpanIdsCSR ===
SnakeCutSbmProcess::KnotSpanIdsCSR
SnakeCutSbmProcess::CreateSkinConditionsPerKnotSpanMatrix(
    const ModelPart& rSkinSubModelPart,
    const SnakeCutSbmProcess::KnotSpanIdsCSR& rReferenceMatrix) const
{
    KnotSpanIdsCSR out;

    // Copy metadata from reference
    out.NumberOfSpansX = rReferenceMatrix.NumberOfSpansX;
    out.NumberOfSpansY = rReferenceMatrix.NumberOfSpansY;
    out.MinU = rReferenceMatrix.MinU; out.MaxU = rReferenceMatrix.MaxU;
    out.MinV = rReferenceMatrix.MinV; out.MaxV = rReferenceMatrix.MaxV;
    out.SpanSizeX = rReferenceMatrix.SpanSizeX; out.SpanSizeY = rReferenceMatrix.SpanSizeY;

    const SizeType Nx = out.NumberOfSpansX, Ny = out.NumberOfSpansY;
    if (Nx == 0 || Ny == 0) { out.Occupancy.resize(0, 0, false); return out; }
    if (rSkinSubModelPart.NumberOfConditions() == 0) { out.Occupancy.resize(Nx, Ny, false); return out; }

    const double u0 = out.MinU, u1 = out.MaxU, v0 = out.MinV, v1 = out.MaxV;
    const double hx = out.SpanSizeX, hy = out.SpanSizeY;
    const double tol = 1e-12;

    auto idx = [tol](double c, double mn, double mx, double h, SizeType n) {
        double cl = c;
        if (c < mn) { KRATOS_ERROR_IF(c < mn - tol) << "coord < min\n"; cl = mn; }
        else if (c > mx) { KRATOS_ERROR_IF(c > mx + tol) << "coord > max\n"; cl = mx; }
        SizeType k = static_cast<SizeType>(std::floor((cl - mn) / h + tol));
        if (k >= n) k = n - 1;
        return k;
    };

    // Sparsity pattern
    std::vector<std::vector<SizeType>> cols(Nx);
    for (const auto& r_cond : rSkinSubModelPart.Conditions()) {
        const auto& g = r_cond.GetGeometry();
        if (g.size() == 0) continue;
        array_1d<double,3> c = ZeroVector(3);
        for (IndexType i = 0; i < g.size(); ++i) c += g[i].Coordinates();
        c /= static_cast<double>(g.size());
        const SizeType ix = idx(c[0], u0, u1, hx, Nx);
        const SizeType iy = idx(c[1], v0, v1, hy, Ny);
        cols[ix].push_back(iy);
    }
    SizeType nnz = 0;
    for (auto& c : cols) { std::sort(c.begin(), c.end()); c.erase(std::unique(c.begin(), c.end()), c.end()); nnz += c.size(); }

    auto& A = out.Occupancy;
    A.resize(Nx, Ny, false);
    A.reserve(nnz);
    for (SizeType i = 0; i < Nx; ++i) for (SizeType j : cols[i]) A.push_back(i, j, 0.0);

    // Fill
    std::vector<std::vector<IndexType>> tmp(nnz);
    for (const auto& r_cond : rSkinSubModelPart.Conditions()) {
        const auto& g = r_cond.GetGeometry();
        if (g.size() == 0) continue;
        array_1d<double,3> c = ZeroVector(3);
        for (IndexType i = 0; i < g.size(); ++i) c += g[i].Coordinates();
        c /= static_cast<double>(g.size());
        const SizeType ix = idx(c[0], u0, u1, hx, Nx);
        const SizeType iy = idx(c[1], v0, v1, hy, Ny);
        const SizeType k = FindNnzIndex(A, ix, iy);
        KRATOS_DEBUG_ERROR_IF(k == static_cast<SizeType>(-1)) << "pattern miss\n";
        tmp[k].push_back(r_cond.Id());
        A.value_data()[k] += 1.0;
    }

    CommitPayload(tmp, out);
    return out;
}

void SnakeCutSbmProcess::CreateSbmExtendedGeometries()
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
void SnakeCutSbmProcess::CreateSbmExtendedGeometries(
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart)
{
    // get the data

    // Create the testBins for the search in radius
    PointVector points;
    for (auto &i_node : rSkinSubModelPart.Nodes()) {
        points.push_back(PointTypePointer(new PointType(i_node.Id(), i_node.X(), i_node.Y(), i_node.Z())));
    }
    // Get the mesh sizes from the surrogate model part
    const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    double knot_span_reference_size = knot_span_sizes[0];
    if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    const int domain_size = mpIgaModelPart->GetProcessInfo()[DOMAIN_SIZE];
    double search_radius;
    if (domain_size == 2) {
        search_radius = 2*std::sqrt(2.0) * knot_span_reference_size;
    } else {
        KRATOS_ERROR << "This method is only implemented for 2D (DOMAIN_SIZE == 2). "
                    << "Current DOMAIN_SIZE: " << domain_size << std::endl;
    }

    DynamicBins testBins(points.begin(), points.end());

    // Maximum number of results to be found in the search in radius
    const int number_of_results = 1e6; 

    ModelPart::NodesContainerType::ContainerType results(number_of_results);
    std::vector<double> list_of_distances(number_of_results);

    BinSearchParameters bin_search_parameters(
        testBins, 
        number_of_results, 
        results, 
        list_of_distances, 
        search_radius);

    auto p_surface = mpIgaModelPart->pGetGeometry(1);
    IndexType id_brep_curve_on_surface = (mpIgaModelPart->GeometriesEnd()-1)->Id() + 1;
    
    auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
                            p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // Build knot-span lookups and set projections from surrogate to skin
    
    const auto skin_nodes_per_knot_span = CreateSkinNodesPerKnotSpanMatrix(rSkinSubModelPart, rSurrogateSubModelPart);
    const auto skin_conditions_per_knot_span = CreateSkinConditionsPerKnotSpanMatrix(rSkinSubModelPart, skin_nodes_per_knot_span);

    
    SetSurrogateToSkinProjections<TIsInnerLoop>(rSurrogateSubModelPart, rSkinSubModelPart, skin_nodes_per_knot_span);
    // Loop over the nodes of the surrogate sub model part
    IndexType iel = 1;
    SizeType brep_degree = p_nurbs_surface->PolynomialDegree(0);
    SizeType number_of_shape_functions_derivatives = 2*brep_degree+1;

    mCutApproximationOrder = brep_degree;

    IndexType first_condition_id;
    IndexType last_condition_id;
    IndexType starting_brep_id;
    SizeType size_surrogate_loop;

    const bool is_inner = TIsInnerLoop;
    if constexpr (TIsInnerLoop)  {
        first_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[0].Id();
        last_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[1].Id();
        size_surrogate_loop = last_condition_id - first_condition_id + 1;
        if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) 
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

    for (SizeType j = 0; j < size_surrogate_loop; ++j) {
        auto p_brep_geometry = mpIgaModelPart->pGetGeometry(starting_brep_id + j);
        auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

        KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeCutSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
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

        // retrieve middle point of the brep
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

        // STORE THE SURROGATE BREP MIDDLE GEOMETRY FOR THE LATERAL BREPS
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

        const double t0 = brep_domain_interval.GetT0();
        const double t1 = brep_domain_interval.GetT1();
        const double dt = (t1 - t0) / static_cast<double>(mInternalDivision);

        //------------------------------------------------------------------
        // 3. Loop over sub-intervals of the upper curve
        //------------------------------------------------------------------
        Node::Pointer p_first_node = p_surrogate_node_1;
        Node::Pointer p_second_node = nullptr;

        // FIXME:
        auto connected_layers_1 = p_surrogate_node_1->GetValue(CONNECTED_LAYERS);
        auto connected_layers_2 = p_surrogate_node_2->GetValue(CONNECTED_LAYERS);
        const auto projection_id_1 = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
        const auto projection_id_2 = p_surrogate_node_2->GetValue(PROJECTION_NODE_ID);
        const auto projection_node_1 = rSkinSubModelPart.pGetNode(projection_id_1);
        const auto projection_node_2 = rSkinSubModelPart.pGetNode(projection_id_2);
        SizeType n_skin_nodes = rSkinSubModelPart.NumberOfNodes();
        IndexType first_node_id = rSkinSubModelPart.NodesBegin()->Id();
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

        SizeType current_internal_divisions = mInternalDivision;

        if (abs(projection_id_2 - projection_id_1) < mInternalDivision*10 || 
            norm_2(projection_node_2->Coordinates() - projection_node_1->Coordinates()) < knot_span_reference_size/5)
            current_internal_divisions = 1;

            const SizeType subdivision_depth = current_internal_divisions;
            const SizeType segment_count = subdivision_depth == 0 ? 1 : static_cast<SizeType>(1) << subdivision_depth;
        
            const double dt_subdivision = (t1 - t0) / static_cast<double>(segment_count);

            std::vector<Node::Pointer> surrogate_segment_nodes(segment_count + 1);
            surrogate_segment_nodes[0] = p_surrogate_node_1;
            surrogate_segment_nodes[segment_count] = p_surrogate_node_2;
    
            for (SizeType k = 1; k < segment_count; ++k) {
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
                                                    SizeType depth,
                                                    SizeType left_index,
                                                    SizeType right_index) -> void
            {
                if (depth == 0) {
                    return;
                }
    
                const SizeType mid_index = (left_index + right_index) / 2;
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
                normal_direction[0] = normal_direction[1];
                normal_direction[1] = -temp;
                const double normal_norm = norm_2(normal_direction);
                if (normal_norm > 1.0e-16) {
                    normal_direction /= normal_norm;
                }

                IndexType id_skin_node = -1;
                if (norm_2(r_left_skin_node.Coordinates()-r_right_skin_node.Coordinates()) < integration_parameters.KnotSpanSizes[0]/1e5)
                    id_skin_node = segment_projection_ids[left_index];
                
                else
                {
                    id_skin_node = FindClosestNodeInLayerWithDirection(
                    skin_mid_point_coords,
                    common_layer_name,
                    rSkinSubModelPart,
                    integration_parameters.KnotSpanSizes,
                    skin_conditions_per_knot_span,
                    normal_direction);
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
    
            for (SizeType k = 1; k < segment_count; ++k) {
                const IndexType projection_id = segment_projection_ids[k];
                KRATOS_ERROR_IF(projection_id == std::numeric_limits<IndexType>::max())
                    << "::[SnakeCutSbmProcess]:: Missing projection id for surrogate subdivision node at index "
                    << k << std::endl;
                surrogate_segment_nodes[k]->SetValue(PROJECTION_NODE_ID, projection_id);
            }
    
            for (SizeType segment = 0; segment < segment_count; ++segment) {
                auto p_first = surrogate_segment_nodes[segment];
                auto p_second = surrogate_segment_nodes[segment + 1];
    
                CreateCutAndSkinQuadraturePoints<TIsInnerLoop>(
                    integration_parameters,
                    bin_search_parameters,
                    p_nurbs_surface,
                    p_first,
                    p_second,
                    surrogate_brep_middle_geometry,
                    *mpIgaModelPart,
                    rSkinSubModelPart);
            }
        }

    //     for (SizeType d = 0; d < current_internal_divisions; ++d)
    //     {
    //         //---------------- 3.1 Current sub-interval --------------------
    //         const double sub_t0 = t0 +  d      * dt;
    //         const double sub_t1 = t0 + (d + 1) * dt;
    //         NurbsInterval sub_interval(sub_t0, sub_t1);

    //         surrogate_vertex_1_local_coords[0] = sub_t0;
    //         surrogate_vertex_2_local_coords[0] = sub_t1;

    //         p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
    //         p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

    //         if (d != current_internal_divisions - 1) {
    //             // create the new node 
    //             p_second_node = Node::Pointer(new Node(0, surrogate_vertex_2));

    //             auto second_node_connected_layers = p_second_node->GetValue(CONNECTED_LAYERS);
    //             second_node_connected_layers.push_back(common_layer_name);
    //             p_second_node->SetValue(CONNECTED_LAYERS, second_node_connected_layers);

    //             auto second_node_connected_conditions = p_second_node->GetValue(CONNECTED_CONDITIONS);
    //             second_node_connected_conditions.push_back(common_condition_name);
    //             p_second_node->SetValue(CONNECTED_CONDITIONS, second_node_connected_conditions);

    //             IndexType projection_node_id = (projection_id_1 * (current_internal_divisions -1 - d) + projection_id_2 * (d+1)) / (current_internal_divisions);

    //             projection_node_id -= first_node_id;
    //             if (abs(projection_id_2-projection_id_1) > n_skin_nodes/2) 
    //                 if (projection_id_2 > projection_id_1) 
    //                     projection_node_id = ((projection_id_1 + n_skin_nodes - first_node_id) * (current_internal_divisions -1 - d) + (projection_id_2-first_node_id) * (d+1)) / (current_internal_divisions);
    //                 else
    //                     projection_node_id = ((projection_id_1-first_node_id) * (current_internal_divisions -1 - d) + (projection_id_2 + n_skin_nodes - first_node_id) * (d+1)) / (current_internal_divisions);

    //             if (projection_node_id > n_skin_nodes) projection_node_id -= n_skin_nodes;
                
    //             projection_node_id += first_node_id;

    //             p_second_node->SetValue(PROJECTION_NODE_ID, projection_node_id);
    //         }
    //         else {
    //             // last division -> use the second surrogate node
    //             p_second_node = p_surrogate_node_2;
    //         }

    //         CreateCutAndSkinQuadraturePoints<TIsInnerLoop>(integration_parameters, bin_search_parameters,
    //                                     p_nurbs_surface, p_first_node, p_second_node,
    //                                     surrogate_brep_middle_geometry,
    //                                     *mpIgaModelPart, rSkinSubModelPart);    

    //         p_first_node = p_second_node;
    //     }
    // }
    
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
            DynamicBinsPointerType p_surrogate_vertex_1 = DynamicBinsPointerType(new PointType(1, p_surrogate_1->X(), p_surrogate_1->Y(), p_surrogate_1->Z()));

            const IndexType id_closest_true_node = p_surrogate_1->GetValue(PROJECTION_NODE_ID);
            
            const auto skin_vertex_1 = rSkinSubModelPart.pGetNode(id_closest_true_node);

            Vector active_range_knot_vector = ZeroVector(2);
            active_range_knot_vector[0] = 0;
            active_range_knot_vector[1] = 1;
            NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

            // surrogate_1 - skin_1
            // create the brep connecting vertex and closest true point
            Point surrogate_1(*p_surrogate_1);
            Point skin_1(*skin_vertex_1);
            
            Node::Pointer p_first_point = Node::Pointer(new Node(1, surrogate_1));
            Node::Pointer p_second_point = Node::Pointer(new Node(2, skin_1));

            if (is_entering) {
                // change the order to preserve the anticlockwise orientation
                Node::Pointer p_temp_pointer = p_first_point;
                p_first_point = p_second_point;
                p_second_point = p_temp_pointer;
            } 

            auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_first_point, p_second_point, active_range_knot_vector);
            auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
        
            IntegrationInfo brep_integration_info_surrogate1_skin1 = p_brep_curve_surrogate1_skin1->GetDefaultIntegrationInfo();

            brep_integration_info_surrogate1_skin1.SetNumberOfIntegrationPointsPerSpan(0,2*mCutApproximationOrder+1);

            IntegrationPointsArrayType brep_integration_points_list_surrogate1_skin1;
            GeometriesArrayType brep_quadrature_point_list_surrogate1_skin1;

            p_brep_curve_surrogate1_skin1->CreateIntegrationPoints(brep_integration_points_list_surrogate1_skin1, brep_integration_info_surrogate1_skin1);

            const double brep_curve_surrogate1_skin1_length = norm_2(skin_1 - surrogate_1);
            for (auto& integration_point : brep_integration_points_list_surrogate1_skin1) {
                integration_point.SetWeight(integration_point.Weight() * brep_curve_surrogate1_skin1_length);
            }

            p_brep_curve_surrogate1_skin1->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate1_skin1, number_of_shape_functions_derivatives, 
                                                                brep_integration_points_list_surrogate1_skin1, brep_integration_info_surrogate1_skin1);
            
            SizeType id = 1;
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
                *mpCutInterfaceSubModelPart, mCutInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries);

        }

        if (!is_second_surrogate_already_computed) {
            p_surrogate_2->SetValue(ACTIVATION_LEVEL, 1.0); 
           
            const IndexType id_closest_true_node = p_surrogate_2->GetValue(PROJECTION_NODE_ID);
            const auto skin_vertex_2 = rSkinSubModelPart.pGetNode(id_closest_true_node);


            Vector active_range_knot_vector = ZeroVector(2);
            active_range_knot_vector[0] = 0;
            active_range_knot_vector[1] = 1;
            NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

            // surrogate_2 - skin_2
            // create the brep connecting vertex and closest true point
            Point surrogate_2(*p_surrogate_2);
            Point skin_2(*skin_vertex_2);

            Node::Pointer p_first_point = Node::Pointer(new Node(1, surrogate_2));
            Node::Pointer p_second_point = Node::Pointer(new Node(2, skin_2));

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

            brep_integration_info_surrogate2_skin2.SetNumberOfIntegrationPointsPerSpan(0,2*mCutApproximationOrder+1);

            IntegrationPointsArrayType brep_integration_points_list_surrogate2_skin2;
            GeometriesArrayType brep_quadrature_point_list_surrogate2_skin2;

            p_brep_curve_surrogate2_skin2->CreateIntegrationPoints(brep_integration_points_list_surrogate2_skin2, brep_integration_info_surrogate2_skin2);

            const double brep_curve_surrogate2_skin2_length = norm_2(skin_2 - surrogate_2);
            for (auto& integration_point : brep_integration_points_list_surrogate2_skin2) {
                integration_point.SetWeight(integration_point.Weight() * brep_curve_surrogate2_skin2_length);
            }

            p_brep_curve_surrogate2_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate2_skin2, number_of_shape_functions_derivatives, 
                                                                brep_integration_points_list_surrogate2_skin2, brep_integration_info_surrogate2_skin2);

                      
            SizeType id = 1;
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
                *mpCutInterfaceSubModelPart, mCutInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries);

        }
    }
}

template <bool TIsInnerLoop>
void SnakeCutSbmProcess::CreateCutAndSkinQuadraturePoints(
    IntegrationParameters& rIntegrationParameters,
    BinSearchParameters& rBinSearchParameters,
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

    KRATOS_ERROR_IF(!layer_found) << ":::[SnakeCutSbmProcess]::: No common layer found between the two surrogate nodes "
                                    << pSurrogateNode1->Id() << " and " << pSurrogateNode2->Id() <<  "\n"
                                    << pSurrogateNode1->Coordinates() << pSurrogateNode2->Coordinates() << "\n"
                                    << connected_layers_1 << connected_layers_2 << std::endl;

    KRATOS_ERROR_IF(!rIntegrationParameters.pSkinConditionsPerSpan)
        << "::[SnakeCutSbmProcess]::CreateCutAndSkinQuadraturePoints: skin condition span matrix is not initialized." << std::endl;

    const auto& r_skin_conditions_per_span = *rIntegrationParameters.pSkinConditionsPerSpan;

    ModelPart& r_layer_model_part = rIgaModelPart.HasSubModelPart(layer_name) ? 
                                    rIgaModelPart.GetSubModelPart(layer_name) : 
                                    rIgaModelPart.CreateSubModelPart(layer_name);       
                                    

    Vector active_range_knot_vector = ZeroVector(2);
    active_range_knot_vector[0] = 0;
    active_range_knot_vector[1] = 1;
    NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

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
    // auto p_nurbs_curve_skin1_skin2 = FitQuadraticBezierBetween<TIsInnerLoop>(
    //     rSkinSubModelPart, id_closest_true_node, id_closest_true_node_2,
    //     p_skin1_brep_point, p_skin2_brep_point);
    

    int p = 1;
    p = mCutApproximationOrder;

    auto p_nurbs_curve_skin1_skin2 = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);

    Vector temp_interpolation_nodes_id;
    std::vector<IndexType> interpolation_nodes_id;
    std::vector<array_1d<double,3>> interpolation_points;

    // SizeType number_interpolation_cuts = p-1;
    // if (!(norm_2(skin_2 - skin_1) > rIntegrationParameters.KnotSpanSizes[0]/10 && p>1))
    //     number_interpolation_cuts = 0;


    // const SizeType segment_count = number_interpolation_cuts == 0 ? 1 : static_cast<SizeType>(1) << number_interpolation_cuts;
    // const IndexType invalid_projection_id = std::numeric_limits<IndexType>::max();

    // std::vector<IndexType> interpolation_projection_ids(segment_count + 1, invalid_projection_id);
    // interpolation_projection_ids[0] = id_closest_true_node;
    // interpolation_projection_ids[segment_count] = id_closest_true_node_2;

    // auto compute_recursive_interpolation = [&](auto&& self,
    //                                            SizeType depth,
    //                                            SizeType left_index,
    //                                            SizeType right_index) -> void
    // {
    //     if (depth == 0) {
    //         return;
    //     }

    //     const SizeType mid_index = (left_index + right_index) / 2;
    //     if (interpolation_projection_ids[mid_index] != invalid_projection_id) {
    //         if (depth > 1) {
    //             self(self, depth - 1, left_index, mid_index);
    //             self(self, depth - 1, mid_index, right_index);
    //         }
    //         return;
    //     }

    //     const auto& r_left_skin_node = rSkinSubModelPart.GetNode(interpolation_projection_ids[left_index]);
    //     const auto& r_right_skin_node = rSkinSubModelPart.GetNode(interpolation_projection_ids[right_index]);

    //     CoordinatesArrayType skin_mid_point_coords = 0.5 * (r_left_skin_node.Coordinates() + r_right_skin_node.Coordinates());

    //     const array_1d<double, 3> tangent = r_right_skin_node.Coordinates() - r_left_skin_node.Coordinates();
    //     const double tangent_norm = norm_2(tangent);

    //     Vector normal_direction = ZeroVector(3);
    //     if (tangent_norm > 1.0e-16) {
    //         normal_direction = tangent / tangent_norm;

    //         const double temp = normal_direction[0];
    //         normal_direction[0] = normal_direction[1];
    //         normal_direction[1] = -temp;

    //         const double normal_norm = norm_2(normal_direction);
    //         if (normal_norm > 1.0e-16) {
    //             normal_direction /= normal_norm;
    //         }
    //     }

    //     IndexType id_skin_node = interpolation_projection_ids[left_index];
    //     if (tangent_norm >= rIntegrationParameters.KnotSpanSizes[0] / 1.0e5) {
    //         id_skin_node = FindClosestNodeInLayerWithDirection(
    //             skin_mid_point_coords,
    //             layer_name,
    //             rSkinSubModelPart,
    //             rIntegrationParameters.KnotSpanSizes,
    //             r_skin_conditions_per_span,
    //             normal_direction);
    //     }

    //     interpolation_projection_ids[mid_index] = id_skin_node;

    //     if (depth > 1) {
    //         self(self, depth - 1, left_index, mid_index);
    //         self(self, depth - 1, mid_index, right_index);
    //     }
    // };

    // if (segment_count > 1) {
    //     compute_recursive_interpolation(compute_recursive_interpolation, number_interpolation_cuts, 0, segment_count);
    // }

    // interpolation_nodes_id.reserve(segment_count + 1);
    // temp_interpolation_nodes_id.resize(segment_count+1);
    // interpolation_points.reserve(segment_count + 1);

    // for (SizeType i = 0; i <= segment_count; ++i) {
    //     const IndexType projection_id = interpolation_projection_ids[i];
    //     KRATOS_ERROR_IF(projection_id == invalid_projection_id)
    //         << "::[SnakeCutSbmProcess]:: Missing interpolation projection id at subdivision index " << i << std::endl;

    //     const auto& r_skin_node = rSkinSubModelPart.GetNode(projection_id);
    //     interpolation_nodes_id.push_back(projection_id);
    //     temp_interpolation_nodes_id[i] = projection_id;
    //     interpolation_points.push_back(r_skin_node.Coordinates());
    // }

    // const double ridge = 1e-12;

    // if (norm_2(skin_2 - skin_1) > rIntegrationParameters.KnotSpanSizes[0]/10 && p>1)
    //     p_nurbs_curve_skin1_skin2 = FitBezierUV_LS_Generic(interpolation_points, p, ridge);

    if (norm_2(skin_2 - skin_1) > rIntegrationParameters.KnotSpanSizes[0]/10 && abs(id_closest_true_node-id_closest_true_node_2)> (2*p+1) && p>1)
        p_nurbs_curve_skin1_skin2 = FitUV_BetweenSkinNodes_Generic<TIsInnerLoop>(
            rSkinSubModelPart, *pNurbsSurface, id_closest_true_node, id_closest_true_node_2, p, /*ridge=*/1e-14);


    auto p_brep_curve_skin1_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2);    
    
    IntegrationPointsArrayType brep_integration_points_list_skin1_skin2;
    GeometriesArrayType brep_quadrature_point_list_skin1_skin2;
    rIntegrationParameters.CurveIntegrationInfo.SetNumberOfIntegrationPointsPerSpan(0, ((mCutApproximationOrder+1)*(mCutApproximationOrder+1)));

    p_brep_curve_skin1_skin2->CreateIntegrationPoints(brep_integration_points_list_skin1_skin2, rIntegrationParameters.CurveIntegrationInfo);

    const double p_brep_curve_skin1_skin2_length = (p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX))->Length();

    // FIXME: fix the value of the weights
    for (auto& integration_point : brep_integration_points_list_skin1_skin2) {
        std::vector<CoordinatesArrayType> ders;
        p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->GlobalSpaceDerivatives(ders, integration_point.Coordinates(), 1);
        const double weight_correcting_factor = std::sqrt(ders[1][0]*ders[1][0] +ders[1][1]*ders[1][1]);

        integration_point.SetWeight(integration_point.Weight() * weight_correcting_factor);
    }

    if (norm_2(skin_2 - skin_1) > 1e-12)
        p_brep_curve_skin1_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin1_skin2, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                                  brep_integration_points_list_skin1_skin2, rIntegrationParameters.CurveIntegrationInfo);
    
    
    for (auto& quad_geom : brep_quadrature_point_list_skin1_skin2)
    {
        quad_geom.SetValue(INTERPOLATION_NODES_ID, interpolation_nodes_id);
        // quad_geom.SetValue(TEMP_INTERPOLATION_NODES_ID, temp_interpolation_nodes_id);
    }

    SizeType id = 1;
    if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
        id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;
    
    std::vector<Geometry<Node>::Pointer> neighbour_geometries_skin1_skin2;
    neighbour_geometries_skin1_skin2.push_back(rSurrogateBrepMiddleGeometry);

    // Find the closest point on the true boundary for each of the quadrature geometries
    for (auto& quad_geom: brep_quadrature_point_list_skin1_skin2)
    {
        CoordinatesArrayType quadrature_point_coords = quad_geom.Center().Coordinates();
        DynamicBinsPointerType p_quadrature_point = DynamicBinsPointerType(new PointType(1, quadrature_point_coords[0], quadrature_point_coords[1], quadrature_point_coords[2]));
        
        auto normal = quad_geom.Normal(0);
        normal = normal / MathUtils<double>::Norm(normal);

        // IndexType id_skin_node = FindClosestNodeInLayer(p_quadrature_point, rBinSearchParameters, layer_name, rSkinSubModelPart);
        
        IndexType id_skin_node = FindClosestNodeInLayerWithDirection(quadrature_point_coords, layer_name, rSkinSubModelPart, rIntegrationParameters.KnotSpanSizes,
            r_skin_conditions_per_span, normal);
        NodeType::Pointer p_skin_node = rSkinSubModelPart.pGetNode(id_skin_node);
        NodePointerVector empty_vector;

        if constexpr(TIsInnerLoop)
            p_skin_node->SetValue(IDENTIFIER, "inner");
        else 
            p_skin_node->SetValue(IDENTIFIER, "outer");
        empty_vector.push_back(p_skin_node); // Just it_node-plane neighbours
        quad_geom.SetValue(NEIGHBOUR_NODES, empty_vector);
    }
    
    this->CreateConditions(
        brep_quadrature_point_list_skin1_skin2.ptr_begin(), brep_quadrature_point_list_skin1_skin2.ptr_end(),
        r_layer_model_part, condition_name, id, PropertiesPointerType(), rIntegrationParameters.KnotSpanSizes, neighbour_geometries_skin1_skin2);

    // FIXME: check for void elements/true coincident with surrogate boundary
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
    if constexpr (TIsInnerLoop)  
    {

        array_1d<double,3> surrogate_1 = pSurrogateNode1->Coordinates();
        array_1d<double,3> surrogate_2 = pSurrogateNode2->Coordinates();
        array_1d<double,3> skin_1 = p_skin_node_1->Coordinates();
        array_1d<double,3> skin_2 = p_skin_node_2->Coordinates();

        IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
            (mCutApproximationOrder+1)*2, /*Order*/
            *p_brep_curve_surrogate1_surrogate2,   // B0
            *p_brep_curve_surrogate1_skin1,       // L0
            *p_brep_curve_surrogate2_skin2,       // L1
            *p_brep_curve_skin1_skin2,            // B1
            surrogate_1,  // P00
            skin_1,       // P01
            surrogate_2,  // P10
            skin_2);      // P11

        
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
            (mCutApproximationOrder+1)*2, /*Order*/
            *p_brep_curve_surrogate2_surrogate1,   // B0
            *p_brep_curve_surrogate2_skin2,       // L0
            *p_brep_curve_surrogate1_skin1,       // L1
            *p_brep_curve_skin2_skin1,            // B1
            surrogate_2,  // P00
            skin_2,       // P01
            surrogate_1,  // P10
            skin_1);      // P11

        pNurbsSurface->CreateQuadraturePointGeometries(surface_quadrature_point_list, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                        surface_integration_points, surface_integration_info);
    }

    IndexType id_element = 1;
    if (mpCutElementsSubModelPart->GetRootModelPart().Elements().size() > 0)
        id_element = mpCutElementsSubModelPart->GetRootModelPart().Elements().back().Id() + 1;

    this->CreateElements(
        surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
        *mpCutElementsSubModelPart, mCutElementName, id_element, PropertiesPointerType(), neighbour_geometries_skin1_skin2);
}






bool SnakeCutSbmProcess::ProjectToSkinBoundary(
        const ModelPart* pSkinModelPart,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rProjectedPoint,
        CoordinatesArrayType& rProjectedPointLocal,
        std::vector<array_1d<double, 3>>& rCurveDerivatives,
        int nInitialGuesses)
{
    rProjectedPoint = ZeroVector(3);
    rCurveDerivatives.resize(3);
    bool is_projected_at_least_once = false;
    double best_distance = 1e12;
    std::vector<array_1d<double, 3>> best_curve_derivatives(2, ZeroVector(3));
    std::string best_layer_name = "";

    for (auto &i_curve : pSkinModelPart->Geometries())
    {   
        int nurbs_curve_id = i_curve.Id();
        auto p_nurbs_curve_geometry = pSkinModelPart->pGetGeometry(nurbs_curve_id);
        auto nurbs_curve_geometry = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_nurbs_curve_geometry);
        KRATOS_ERROR_IF(!nurbs_curve_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << nurbs_curve_id 
                            << " is not a NurbsCurveGeometryType." << std::endl;

        const double t0 = nurbs_curve_geometry->DomainInterval().GetT0();
        const double t1 = nurbs_curve_geometry->DomainInterval().GetT1();

        for (int i_guess = 0; i_guess < nInitialGuesses; ++i_guess) {
            CoordinatesArrayType projected_point_local = ZeroVector(3);
            CoordinatesArrayType projected_point = ZeroVector(3);
            std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));

            projected_point_local[0] = t0 + (t1 - t0) * double(i_guess) / (nInitialGuesses - 1);

            bool is_projected = nurbs_curve_geometry->ProjectionPointGlobalToLocalSpace(rPoint, projected_point_local, 1e-13);

            if (!is_projected) continue;

            nurbs_curve_geometry->GlobalCoordinates(projected_point, projected_point_local);

            double curr_distance = norm_2(rPoint - projected_point);

            if (curr_distance < best_distance) {
                best_distance = curr_distance;
                rProjectedPoint = projected_point;
                nurbs_curve_geometry->GlobalSpaceDerivatives(best_curve_derivatives, projected_point_local, 2);
                rProjectedPointLocal = projected_point_local;
                is_projected_at_least_once = true;
                best_layer_name = i_curve.GetValue(IDENTIFIER);


            }
        }
    }

    rCurveDerivatives = best_curve_derivatives;

    if (!is_projected_at_least_once)
    {
        KRATOS_WARNING("::[IgaContactProcessSbm]:: no projection found on the skin boundary")  
                    << " for the point: " << rPoint << std::endl;
    }

    return is_projected_at_least_once;
}


void SnakeCutSbmProcess::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rConditionName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const Vector KnotSpanSizes,
    const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const
{
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    IndexType geometry_count = 0;
    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_condition_list.push_back(
            reference_condition.Create(rIdCounter, (*it), pProperties));
        
        // Set knot span sizes to the condition
        new_condition_list.GetContainer()[geometry_count]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);

        new_condition_list.GetContainer()[geometry_count]->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);

        rIdCounter++;
        geometry_count++;
    }

    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}


void SnakeCutSbmProcess::CreateElements(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rElementName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const
{
    KRATOS_ERROR_IF(!KratosComponents<Element>::Has(rElementName))
        << rElementName << " not registered." << std::endl;

    const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

    ElementsContainerType new_element_list;

    KRATOS_INFO_IF("CreateElements", mEchoLevel > 2)
        << "Creating elements of type " << rElementName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    SizeType num_elements = std::distance(rGeometriesBegin, rGeometriesEnd);
    new_element_list.reserve(num_elements);
    IndexType geometry_count = 0;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_element_list.push_back(
            rReferenceElement.Create(rIdCounter, (*it), pProperties));

        new_element_list.GetContainer()[geometry_count]->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);
        rIdCounter++;
        geometry_count++;
    }

    rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
}


// ------------------------------------------------------------------
// 1-D Gauss–Legendre on [0,1]
// ------------------------------------------------------------------
void SnakeCutSbmProcess::GaussLegendreOnUnitInterval(
    const std::size_t      Order,
    std::vector<double>&   rXi,
    std::vector<double>&   rWeight)
{
    KRATOS_ERROR_IF(Order < 1 || Order > 10)
        << "Gauss order " << Order << " not implemented (1…10 supported)." << std::endl;

    // Gauss–Legendre nodes/weights on [-1,1] (rows: n=1..10).
    // Only the first 'n' entries of each row are used; the rest are zeros.
    static const double sXi[10][10] = {
        // n = 1
        {  0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 2
        { -0.5773502691896257,  +0.5773502691896257, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 3
        {  0.0, -0.7745966692414834, +0.7745966692414834, 0, 0, 0, 0, 0, 0, 0 },
        // n = 4
        { -0.3399810435848563, +0.3399810435848563, -0.8611363115940526, +0.8611363115940526, 0, 0, 0, 0, 0, 0 },
        // n = 5
        {  0.0, -0.5384693101056831, +0.5384693101056831, -0.9061798459386640, +0.9061798459386640, 0, 0, 0, 0, 0 },
        // n = 6
        { -0.2386191860831969, +0.2386191860831969, -0.6612093864662645, +0.6612093864662645,
          -0.9324695142031521, +0.9324695142031521, 0, 0, 0, 0 },
        // n = 7
        {  0.0, -0.4058451513773972, +0.4058451513773972, -0.7415311855993945, +0.7415311855993945,
          -0.9491079123427585, +0.9491079123427585, 0, 0, 0 },
        // n = 8
        { -0.1834346424956498, +0.1834346424956498, -0.5255324099163290, +0.5255324099163290,
          -0.7966664774136267, +0.7966664774136267, -0.9602898564975363, +0.9602898564975363, 0, 0 },
        // n = 9
        {  0.0, -0.3242534234038089, +0.3242534234038089, -0.6133714327005904, +0.6133714327005904,
          -0.8360311073266358, +0.8360311073266358, -0.9681602395076261, +0.9681602395076261, 0 },
        // n = 10
        { -0.1488743389816312, +0.1488743389816312, -0.4333953941292472, +0.4333953941292472,
          -0.6794095682990244, +0.6794095682990244, -0.8650633666889845, +0.8650633666889845,
          -0.9739065285171717, +0.9739065285171717 }
    };

    static const double sW[10][10] = {
        // n = 1
        {  2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 2
        {  1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 3
        {  0.8888888888888888, 0.5555555555555556, 0.5555555555555556, 0, 0, 0, 0, 0, 0, 0 },
        // n = 4
        {  0.6521451548625461, 0.6521451548625461, 0.3478548451374539, 0.3478548451374539, 0, 0, 0, 0, 0, 0 },
        // n = 5
        {  0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891, 0, 0, 0, 0, 0 },
        // n = 6
        {  0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.3607615730481386,
           0.1713244923791704, 0.1713244923791704, 0, 0, 0, 0 },
        // n = 7
        {  0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766,
           0.1294849661688697, 0.1294849661688697, 0, 0, 0 },
        // n = 8
        {  0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873,
           0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763, 0, 0 },
        // n = 9
        {  0.3302393550012598, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354,
           0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0 },
        // n = 10
        {  0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963,
           0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806,
           0.0666713443086881, 0.0666713443086881 }
    };

    const std::size_t n = Order;

    rXi.resize(n);
    rWeight.resize(n);

    for (std::size_t k = 0; k < n; ++k) {
        // Map node from [-1,1] to [0,1]
        rXi[k]     = 0.5 * (sXi[n-1][k] + 1.0);
        // Rescale weight: w' = w / 2
        rWeight[k] = 0.5 *  sW [n-1][k];
    }
}
// ------------------------------------------------------------------
// Global point on Brep curve
// ------------------------------------------------------------------
array_1d<double,3> SnakeCutSbmProcess::GlobalPoint(
    const GeometryType& rCurve,
    const double         T)
{
    CoordinatesArrayType local(3,0.0);
    local[0] = T;

    array_1d<double,3> P;
    rCurve.GlobalCoordinates(P, local);
    return P;
}

// ------------------------------------------------------------------
// Coons patch mapping X(ξ,η)
// ------------------------------------------------------------------
array_1d<double,3> SnakeCutSbmProcess::CoonsPoint(
    const double                  Xi,
    const double                  Eta,
    const GeometryType&          rB0,
    const GeometryType&          rL0,
    const GeometryType&          rL1,
    const GeometryType&          rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11)
{
    const array_1d<double,3> B0 = GlobalPoint(rB0, Xi);
    const array_1d<double,3> B1 = GlobalPoint(rB1, Xi);
    const array_1d<double,3> L0 = GlobalPoint(rL0, Eta);
    const array_1d<double,3> L1 = GlobalPoint(rL1, Eta);

    const double om_xi  = 1.0 - Xi;
    const double om_eta = 1.0 - Eta;

    return  om_xi * L0 + Xi * L1
          + om_eta * B0 + Eta * B1
          - om_xi * om_eta * rP00
          - Xi    * om_eta * rP10
          - Xi    * Eta    * rP11
          - om_xi * Eta    * rP01;
}

// ------------------------------------------------------------------
// Finite-difference derivative ∂X/∂ξ or ∂X/∂η
// ------------------------------------------------------------------
array_1d<double,3> SnakeCutSbmProcess::CoonsDerivativeFD(
    const double                  Xi,
    const double                  Eta,
    const bool                    WithRespectToXi,
    const GeometryType&          rB0,
    const GeometryType&          rL0,
    const GeometryType&          rL1,
    const GeometryType&          rB1,
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

// ------------------------------------------------------------------
// Build Gauss points on the curved quadrilateral (tensor-product)
// ------------------------------------------------------------------
// WORKING
// SnakeCutSbmProcess::IntegrationPointsArrayType
// SnakeCutSbmProcess::CreateCoonsPatchGaussPoints(
//     const std::size_t             Order,
//     const GeometryType&           rB0,
//     const GeometryType&           rL0,
//     const GeometryType&           rL1,
//     const GeometryType&           rB1,
//     const array_1d<double,3>&     rP00,
//     const array_1d<double,3>&     rP01,
//     const array_1d<double,3>&     rP10,
//     const array_1d<double,3>&     rP11) const
// {
//     IntegrationPointsArrayType gp_list;
//     gp_list.reserve(Order * Order);

//     // 1-D Gauss nodes / weights on [0,1]
//     std::vector<double> xi, w;
//     GaussLegendreOnUnitInterval(Order, xi, w);

//     for (std::size_t i=0; i<Order; ++i)
//     for (std::size_t j=0; j<Order; ++j)
//     {
//         const double xi_i  = xi[i];
//         const double eta_j = xi[j];          // same grid for η
//         const double w_ij  = w[i] * w[j];

//         // --- Jacobian ------------------------------------------------
//         const array_1d<double,3> dXi =
//             CoonsDerivativeFD(xi_i, eta_j, true ,
//                               rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);

//         const array_1d<double,3> dEta =
//             CoonsDerivativeFD(xi_i, eta_j, false,
//                               rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);

//         const array_1d<double,3> cross =
//             MathUtils<double>::CrossProduct(dXi, dEta);

//         const double jac = norm_2(cross);    // |X_ξ × X_η|

//         // --- Global coordinates --------------------------------------
//         const array_1d<double,3> X =
//             CoonsPoint(xi_i, eta_j,
//                        rB0,rL0,rL1,rB1,
//                        rP00,rP01,rP10,rP11);

//         // store *global* coords + weight (w·|J|)
//         gp_list.emplace_back( IntegrationPoint<3>( X[0], X[1], X[2],
//                                                    w_ij * jac ) );
//     }

//     return gp_list;
// }

// FIXME: trial for optimal boundary
SnakeCutSbmProcess::IntegrationPointsArrayType
SnakeCutSbmProcess::CreateCoonsPatchGaussPoints(
    const std::size_t             Order,
    const GeometryType&           rB0,
    const GeometryType&           rL0,
    const GeometryType&           rL1,
    const GeometryType&           rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11) const
{
    IntegrationPointsArrayType gp_list;
    gp_list.reserve(Order * Order);

    // 1-D Gauss nodes / weights on [0,1]
    std::vector<double> xi, w;
    GaussLegendreOnUnitInterval(Order, xi, w);

    // --- build a reference normal to define the positive orientation ----
    array_1d<double,3> a = rP10 - rP00;
    array_1d<double,3> b = rP01 - rP00;
    // --- reference normal from Coons derivatives at center --------------
    const double xi_c  = 0.5, eta_c = 0.5;
    const array_1d<double,3> dXi_c  = CoonsDerivativeFD(xi_c, eta_c, true ,
                                                        rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
    const array_1d<double,3> dEta_c = CoonsDerivativeFD(xi_c, eta_c, false,
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
        const array_1d<double,3> dXi  = CoonsDerivativeFD(xi_i, eta_j, true ,
                                                          rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);
        const array_1d<double,3> dEta = CoonsDerivativeFD(xi_i, eta_j, false,
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



void SnakeCutSbmProcess::BuildParabolicNurbsData(
    Node::Pointer         pNode0,
    Node::Pointer         pNodeM,
    Node::Pointer         pNode2,
    PointerVector<Node>&  rCtrlPtsPointerVector,
    Vector&               rKnots,
    Vector&               rWeights,
    const double          tm)
{
    KRATOS_ERROR_IF(tm <= 0.0 || tm >= 1.0)
        << "tm must be in (0,1). Given: " << tm << std::endl;

    // ----------------------------------------------------------------
    // Control node 0 and 2 are the given ones
    // Control node 1 is created on-the-fly (ID = 0 → virtual, change if needed)
    // ----------------------------------------------------------------
    const array_1d<double,3> P0 = pNode0->Coordinates();
    const array_1d<double,3> Pm = pNodeM->Coordinates();
    const array_1d<double,3> P2 = pNode2->Coordinates();

    const double t  = tm;
    const double denom = 2.0 * t * (1.0 - t);

    array_1d<double,3> P1 =
          ( Pm
          - std::pow(1.0 - t, 2) * P0
          - std::pow(    t, 2)   * P2 ) / denom;

    // ----------------------------------------------------------------
    // Create a new node for P1 (ID = 0  → not stored in any ModelPart)
    // If you need it in a ModelPart, use pModelPart->CreateNewNode(...)
    // ----------------------------------------------------------------
    Node::Pointer pNode1 = Node::Pointer(new Node(1, P1));

    // ----------------------------------------------------------------
    // Fill PointerVector<Node>  (size 3)
    // ----------------------------------------------------------------
    rCtrlPtsPointerVector.clear();
    rCtrlPtsPointerVector.reserve(3);
    rCtrlPtsPointerVector.push_back(pNode0);
    rCtrlPtsPointerVector.push_back(pNode1);
    rCtrlPtsPointerVector.push_back(pNode2);

    // ----------------------------------------------------------------
    // Knot vector [0 0 0 1 1 1]  (Vector of size 6)
    // ----------------------------------------------------------------
    if (rKnots.size() != 6) rKnots.resize(6);
    rKnots[0] = rKnots[1] = rKnots[2] = 0.0;
    rKnots[3] = rKnots[4] = rKnots[5] = 1.0;

    // ----------------------------------------------------------------
    // Weights  (all 1 for polynomial parabola)
    // ----------------------------------------------------------------
    if (rWeights.size() != 3) rWeights.resize(3);
    rWeights[0] = rWeights[1] = rWeights[2] = 1.0;
}

template <bool TIsInnerLoop>
void SnakeCutSbmProcess::SetSurrogateToSkinProjections(
    const ModelPart& rSurrogateSubModelPart,
    const ModelPart& rSkinSubModelPart,
    const KnotSpanIdsCSR& rSkinNodesPerSpan)
{
    const auto& r_parent_model_part = rSurrogateSubModelPart.GetParentModelPart();
    const Vector& knot_span_sizes = r_parent_model_part.GetValue(KNOT_SPAN_SIZES);
    const auto& parameter_space_corners = r_parent_model_part.GetValue(PARAMETER_SPACE_CORNERS);

    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "::[SnakeCutSbmProcess]::SetSurrogateToSkinProjections: KNOT_SPAN_SIZES must contain at least two entries." << std::endl;
    KRATOS_ERROR_IF(parameter_space_corners.size() < 2)
        << "::[SnakeCutSbmProcess]::SetSurrogateToSkinProjections: PARAMETER_SPACE_CORNERS must contain two vectors." << std::endl;
    KRATOS_ERROR_IF(parameter_space_corners[0].size() < 2 || parameter_space_corners[1].size() < 2)
        << "::[SnakeCutSbmProcess]::SetSurrogateToSkinProjections: PARAMETER_SPACE_CORNERS entries must contain min and max values." << std::endl;

    const double span_size_x = knot_span_sizes[0];
    const double span_size_y = knot_span_sizes[1];
    const double min_u = parameter_space_corners[0][0];
    const double max_u = parameter_space_corners[0][1];
    const double min_v = parameter_space_corners[1][0];
    const double max_v = parameter_space_corners[1][1];

    const SizeType span_count_x = rSkinNodesPerSpan.NumberOfSpansX;
    const SizeType span_count_y = rSkinNodesPerSpan.NumberOfSpansY;

    auto clamp_coordinate = [](double coordinate, double min_value, double max_value) {
        if (coordinate < min_value) {
            return min_value;
        }
        if (coordinate > max_value) {
            return max_value;
        }
        return coordinate;
    };

    constexpr double tol = 1.0e-12;

    auto compute_candidate_indices = [&](double coordinate,
                                         double min_value,
                                         double max_value,
                                         double span_size,
                                         SizeType span_count) -> std::vector<SizeType>
    {
        std::vector<SizeType> indices;
        if (span_count == 0) {
            return indices;
        }

        const double clamped = clamp_coordinate(coordinate, min_value, max_value);
        const double relative = (clamped - min_value) / span_size;

        const double relative_minus = std::max(0.0, relative - tol);
        const double relative_plus = std::min(static_cast<double>(span_count), relative + tol);

        SizeType index_min = static_cast<SizeType>(std::floor(relative_minus));
        if (index_min >= span_count) {
            index_min = span_count - 1;
        }

        SizeType index_max = static_cast<SizeType>(std::floor(relative_plus));
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
        auto append_ids_by_k = [&](const KnotSpanIdsCSR& S, SizeType k) {
            const auto view = CellIdsByK(S, k);
            if (view.size) {
                candidate_ids.insert(candidate_ids.end(), view.data, view.data + view.size);
            }
        };

        const auto& S = rSkinNodesPerSpan;                 // KnotSpanIdsCSR
        const auto& A = S.Occupancy;

        // Query only cells (ix,iy)
        for (const SizeType ix : x_indices) {
            if (ix >= S.NumberOfSpansX) continue;
            for (const SizeType iy : y_indices) {
                if (iy >= S.NumberOfSpansY) continue;

                const SizeType k = FindNnzIndex(A, ix, iy);
                if (k != static_cast<SizeType>(-1)) append_ids_by_k(S, k);
            }
        }

        // If empty, collect all ids from all nonzeros
        if (candidate_ids.empty()) {
            const SizeType nnz = static_cast<SizeType>(A.value_data().size());
            for (SizeType k = 0; k < nnz; ++k) {
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

        for (const IndexType candidate_id : candidate_ids) {
            const auto& r_candidate_node = rSkinSubModelPart.GetNode(candidate_id);
            const auto& candidate_layers = r_candidate_node.GetValue(CONNECTED_LAYERS);

            bool matches_forced = true;
            for (const auto& forced_layer : forced_layers) {
                if (std::find(candidate_layers.begin(), candidate_layers.end(), forced_layer) == candidate_layers.end()) {
                    matches_forced = false;
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

            if ((forced_layers.empty() || matches_forced) && distance < result.Distance) {
                result.Id = candidate_id;
                result.Distance = distance;
                result.MatchesForcedLayers = true;
            }
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
                << "::[SnakeCutSbmProcess]::SetSurrogateToSkinProjections: no skin node candidate found for surrogate node "
                << p_surrogate_node_1->Id() << " at " << p_surrogate_node_1->Coordinates() << std::endl;

            skin_node_id_1 = selection.Id;

            const auto& connected_layers = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_LAYERS);
            const auto& connected_conditions = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_CONDITIONS);

            p_surrogate_node_1->SetValue(PROJECTION_NODE_ID, skin_node_id_1);
            p_surrogate_node_1->SetValue(CONNECTED_LAYERS, connected_layers);
            p_surrogate_node_1->SetValue(CONNECTED_CONDITIONS, connected_conditions);

            forced_layers = connected_layers;
        }

        if (!has_proj_2) {
            const auto selection = select_candidate(p_surrogate_node_2, forced_layers);
            KRATOS_ERROR_IF(selection.Id == std::numeric_limits<IndexType>::max())
                << "::[SnakeCutSbmProcess]::SetSurrogateToSkinProjections: no skin node candidate found for surrogate node "
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
    SizeType count_intersections = 1;
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

            bool intersects_another_projection = false;
            
            if (SegmentsIntersect(p_surrogate_node_1, p_projection_node_1,
                p_surrogate_node_2, p_projection_node_2)) 
            {
                p_surrogate_node_1->SetValue(PROJECTION_NODE_ID, projection_node_id_2); 
                p_surrogate_node_2->SetValue(PROJECTION_NODE_ID, projection_node_id_1);

                KRATOS_WATCH(p_surrogate_node_1->Coordinates())
                KRATOS_WATCH(p_surrogate_node_2->Coordinates())
                KRATOS_WATCH("-------------------------")

                count_intersections += 1;
            }
        }
    }

    KRATOS_ERROR_IF(iter_check == 5) << "::[SnakeSbmProcess]:: Maximum iteration reached when checking intersections between projections. Please check the input data." << std::endl;

}

void SnakeCutSbmProcess::AssestProjectionsFeasibility(
    const ModelPart& rSkinSubModelPart,
    Node::Pointer pSurrogateNode1, 
    Node::Pointer pSurrogateNode2)
{
    std::vector<std::string> forced_layers = pSurrogateNode1->GetValue(CONNECTED_LAYERS);

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
        bool matches_layers = true;
        for (const auto& forced_layer : forced_layers) {
            if (std::find(candidate_layers.begin(), candidate_layers.end(), forced_layer) == candidate_layers.end()) {
                matches_layers = false;
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
        << "::[SnakeCutSbmProcess]::AssestProjectionsFeasibility: no interface node satisfies the forced layers for surrogate node "
        << pSurrogateNode2->Id() << std::endl;

    const IndexType id_closest_true_node_2 = nearest_node_id;
    const auto& r_projection_node = rSkinSubModelPart.GetNode(id_closest_true_node_2);

    if (norm_2(r_projection_node.Coordinates()-pSurrogateNode1->Coordinates()) <
        norm_2(r_projection_node.Coordinates()-pSurrogateNode2->Coordinates()))
    {
        pSurrogateNode1->SetValue(PROJECTION_NODE_ID, id_closest_true_node_2);

        auto connected_layers = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_LAYERS);
        auto connected_conditions = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_CONDITIONS);

        pSurrogateNode1->SetValue(CONNECTED_LAYERS, connected_layers);
        pSurrogateNode1->SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }
    else
    {
        pSurrogateNode2->SetValue(PROJECTION_NODE_ID, id_closest_true_node_2);

        auto connected_layers = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_LAYERS);
        auto connected_conditions = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_CONDITIONS);

        pSurrogateNode2->SetValue(CONNECTED_LAYERS, connected_layers);
        pSurrogateNode2->SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }
    return;
}

IndexType SnakeCutSbmProcess::FindClosestNodeInLayer(
    const DynamicBinsPointerType& rStartPoint,
    BinSearchParameters& rSearchParameters,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart)
{
    // Reset search parameters
    rSearchParameters.reset();

    // Perform the search
    SizeType obtained_results = rSearchParameters.TestBins.SearchInRadius(
        *rStartPoint,
        rSearchParameters.SearchRadius,
        rSearchParameters.Results.begin(),
        rSearchParameters.ListOfDistances.begin(),
        rSearchParameters.NumberOfResults);

    KRATOS_ERROR_IF(obtained_results == 0) 
        << "::[FindClosestNodeInLayer]:: No points found for projection of point: "
        << rStartPoint << std::endl;

    // Find nearest node id that matches the given layer
    IndexType nearest_node_id = -1;
    double minimum_distance = 1e14;

    for (IndexType k = 0; k < obtained_results; ++k) {
        double current_distance = rSearchParameters.ListOfDistances[k];
        IndexType current_id = rSearchParameters.Results[k]->Id();

        const auto& connected_layers = rSkinSubModelPart.GetNode(current_id).GetValue(CONNECTED_LAYERS);
        bool layer_match = std::find(connected_layers.begin(), connected_layers.end(), rLayer) != connected_layers.end();

        if (layer_match && current_distance < minimum_distance) {
            minimum_distance = current_distance;
            nearest_node_id = current_id;
        }
    }

    KRATOS_ERROR_IF(nearest_node_id == -1) 
        << "::[FindClosestNodeInLayer]:: No node found matching layer: " << rLayer << std::endl;

    return nearest_node_id;
}

IndexType SnakeCutSbmProcess::FindClosestNodeInLayerWithDirection(
    const array_1d<double,3>& rStartPoint,
    const std::string& rLayer,
    const ModelPart& rSkinSubModelPart,
    const Vector& rKnotSpanSizes,
    const KnotSpanIdsCSR& rSkinConditionsPerSpan,
    const Vector& rDirection)
{
    KRATOS_ERROR_IF(rSkinConditionsPerSpan.NumberOfSpansX == 0 ||
                    rSkinConditionsPerSpan.NumberOfSpansY == 0)
        << "::[SnakeCutSbmProcess]::FindClosestNodeInLayerWithDirection: condition span matrix is empty." << std::endl;

    Vector direction = rDirection;
    constexpr double dir_tol = 1.0e-12;
    double direction_norm = norm_2(direction);
    if (direction_norm < dir_tol) {
        direction = ZeroVector(3);
        direction[0] = 1.0;
        direction_norm = 1.0;
    }
    direction /= direction_norm;

    const double span_size_x = (rSkinConditionsPerSpan.SpanSizeX > 0.0) ? rSkinConditionsPerSpan.SpanSizeX
                           : (rKnotSpanSizes.size() > 0 ? rKnotSpanSizes[0] : 1.0);
    const double span_size_y = (rSkinConditionsPerSpan.SpanSizeY > 0.0) ? rSkinConditionsPerSpan.SpanSizeY
                           : (rKnotSpanSizes.size() > 1 ? rKnotSpanSizes[1] : span_size_x);

    double reference_span_size = std::max(span_size_x, span_size_y);
    if (rKnotSpanSizes.size() > 0) {
        reference_span_size = 0.0;
        for (SizeType i = 0; i < rKnotSpanSizes.size(); ++i) {
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

    auto compute_index = [](double coordinate, double min_value, double max_value, double span_size, SizeType span_count){
        double clamped = coordinate;
        if (clamped < min_value) {
            clamped = min_value;
        } else if (clamped > max_value) {
            clamped = max_value;
        }

        if (span_count == 0) {
            return SizeType(0);
        }

        double relative = (clamped - min_value) / span_size;
        SizeType index = static_cast<SizeType>(std::floor(relative + 1.0e-12));
        if (index >= span_count) {
            index = span_count - 1;
        }
        return index;
    };

    IndexType best_node_id = std::numeric_limits<IndexType>::max();
    double best_distance = std::numeric_limits<double>::max();

    auto update_candidate = [&](Node::Pointer p_candidate_node){
        const auto& candidate_layers = p_candidate_node->GetValue(CONNECTED_LAYERS);
        if (std::find(candidate_layers.begin(), candidate_layers.end(), rLayer) == candidate_layers.end()) {
            return;
        }

        array_1d<double,3> diff = p_candidate_node->Coordinates();
        diff -= rStartPoint;
        const double distance = norm_2(diff);

        if (distance < best_distance) {
            best_distance = distance;
            best_node_id = p_candidate_node->Id();
        }
    };

    const SizeType number_spans_x = rSkinConditionsPerSpan.NumberOfSpansX;
    const SizeType number_spans_y = rSkinConditionsPerSpan.NumberOfSpansY;
    const SizeType base_ix = compute_index(rStartPoint[0], min_u, max_u, span_size_x, number_spans_x);
    const SizeType base_iy = compute_index(rStartPoint[1], min_v, max_v, span_size_y, number_spans_y);

    const int max_search_level = static_cast<int>(std::max(number_spans_x, number_spans_y));

    std::vector<IndexType> candidate_conditions;
    candidate_conditions.reserve(32);

    bool found_intersection = false;
    for (int level = 0; level < max_search_level && !found_intersection; ++level) {
        const SizeType extension = static_cast<SizeType>(level + 1);
        const double current_length = reference_span_size * static_cast<double>(extension);

        array_1d<double,3> segment_start = rStartPoint;
        segment_start -= direction * reference_span_size/2;
        array_1d<double,3> segment_end = rStartPoint;
        segment_end += direction * current_length;

        const SizeType min_ix = (base_ix > extension) ? base_ix - extension : 0;
        const SizeType max_ix = std::min<SizeType>(base_ix + extension, number_spans_x > 0 ? number_spans_x - 1 : 0);
        const SizeType min_iy = (base_iy > extension) ? base_iy - extension : 0;
        const SizeType max_iy = std::min<SizeType>(base_iy + extension, number_spans_y > 0 ? number_spans_y - 1 : 0);

        candidate_conditions.clear();
        // Collect unique condition ids from CSR cells in [min_ix..max_ix] × [min_iy..max_iy]
        const auto& S = rSkinConditionsPerSpan;
        std::unordered_set<IndexType> seen;
        seen.reserve(256); // optional

        for (SizeType ix = min_ix; ix <= max_ix; ++ix) {
            if (ix >= S.NumberOfSpansX) continue;
            for (SizeType iy = min_iy; iy <= max_iy; ++iy) {
                if (iy >= S.NumberOfSpansY) continue;

                const SizeType k = FindNnzIndex(S.Occupancy, ix, iy);
                if (k == static_cast<SizeType>(-1)) continue;

                const auto ids = CellIdsByK(S, k);
                for (SizeType t = 0; t < ids.size; ++t) {
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

            const bool intersects = SegmentsIntersect(start_node, end_node, node_a, node_b);

            if (intersects) {
                update_candidate(node_a);
                update_candidate(node_b);
                found_intersection = true;
            }
        }

        if (best_node_id != std::numeric_limits<IndexType>::max()) {
            return best_node_id;
        }
    }

    KRATOS_ERROR << "::[SnakeCutSbmProcess]::FindClosestNodeInLayerWithDirection: no node found for layer "
                 << rLayer << " starting from point " << rStartPoint << std::endl;

    return 0;
}

// IndexType SnakeCutSbmProcess::FindClosestNodeInLayerWithDirection(
//     const array_1d<double,3>& rStartPoint,
//     const std::string& rLayer,
//     const ModelPart& rSkinSubModelPart,
//     const Vector& rKnotSpanSizes,
//     const KnotSpanIdsCSR& rSkinConditionsPerSpan,
//     const Vector& rDirection)
// {
//     // Preconditions
//     KRATOS_ERROR_IF(rSkinConditionsPerSpan.NumberOfSpansX == 0 ||
//                     rSkinConditionsPerSpan.NumberOfSpansY == 0)
//         << "::[SnakeCutSbmProcess]::FindClosestNodeInLayerWithDirection: empty span matrix." << std::endl;

//     // --- Normalize direction; fallback +x ---
//     Vector dir = rDirection;
//     constexpr double tiny = 1e-12;
//     double n = norm_2(dir);
//     if (n < tiny) { dir = ZeroVector(3); dir[0] = 1.0; n = 1.0; }
//     dir /= n;

//     // --- Span sizes and reference scale ---
//     const double span_x = (rSkinConditionsPerSpan.SpanSizeX > 0.0) ? rSkinConditionsPerSpan.SpanSizeX
//                        : (rKnotSpanSizes.size() > 0 ? rKnotSpanSizes[0] : 1.0);
//     const double span_y = (rSkinConditionsPerSpan.SpanSizeY > 0.0) ? rSkinConditionsPerSpan.SpanSizeY
//                        : (rKnotSpanSizes.size() > 1 ? rKnotSpanSizes[1] : span_x);

//     double ref_span = std::max(span_x, span_y);
//     if (rKnotSpanSizes.size() > 0) {
//         ref_span = 0.0;
//         for (SizeType i = 0; i < rKnotSpanSizes.size(); ++i) ref_span = std::max(ref_span, rKnotSpanSizes[i]);
//         if (ref_span <= tiny) ref_span = 1.0;
//     }
//     const double eps = 1e-9 * ref_span;

//     // --- Bounds ---
//     const double min_u = rSkinConditionsPerSpan.MinU;
//     const double max_u = rSkinConditionsPerSpan.MaxU;
//     const double min_v = rSkinConditionsPerSpan.MinV;
//     const double max_v = rSkinConditionsPerSpan.MaxV;

//     const SizeType NX = rSkinConditionsPerSpan.NumberOfSpansX;
//     const SizeType NY = rSkinConditionsPerSpan.NumberOfSpansY;

//     // --- Helpers ---
//     auto clamp_val = [](double x, double a, double b){ return x < a ? a : (x > b ? b : x); };
//     auto compute_index = [&](double c, double cmin, double cmax, double span, SizeType count){
//         if (count == 0) return SizeType(0);
//         double clamped = clamp_val(c, cmin, cmax);
//         SizeType idx = static_cast<SizeType>(std::floor((clamped - cmin) / span + 1e-12));
//         if (idx >= count) idx = count - 1;
//         return idx;
//     };
//     auto in_bounds = [&](SizeType i, SizeType j){ return i < NX && j < NY; };

//     // --- Start cell and point ---
//     SizeType ix = compute_index(rStartPoint[0], min_u, max_u, span_x, NX);
//     SizeType iy = compute_index(rStartPoint[1], min_v, max_v, span_y, NY);

//     array_1d<double,3> p = rStartPoint;                 // marching point
//     array_1d<double,3> seg_start = rStartPoint - dir*eps; // start slightly inside
//     const double max_ray_len = 2.0 * ( (max_u-min_u) + (max_v-min_v) ) + ref_span;

//     // --- Choose base step length along the dominant axis ---
//     const double ax = std::abs(dir[0]);
//     const double ay = std::abs(dir[1]);
//     double base_step = (ax >= ay) ? span_x/ std::max(ax, tiny) : span_y/ std::max(ay, tiny);
//     if (base_step <= 0.0 || !std::isfinite(base_step)) base_step = ref_span;

//     // --- Best candidate bookkeeping ---
//     IndexType best_node_id = std::numeric_limits<IndexType>::max();
//     double best_distance = std::numeric_limits<double>::max();
//     auto update_candidate = [&](const Node& rnode){
//         const auto& layers = rnode.GetValue(CONNECTED_LAYERS);
//         if (std::find(layers.begin(), layers.end(), rLayer) == layers.end()) return;
//         array_1d<double,3> d = rnode.Coordinates();
//         d -= rStartPoint;
//         const double dist = norm_2(d);
//         if (dist < best_distance) { best_distance = dist; best_node_id = rnode.Id(); }
//     };

//     // --- March: unlock only cells entered along dir ---
//     double traveled = 0.0;
//     SizeType safety = (NX + NY) * 4 + 16; // conservative cap

//     while (in_bounds(ix, iy) && traveled <= max_ray_len && safety--) {
//         // Process current cell (ix,iy) with segment [seg_start, seg_end)
//         // seg_end will be computed after deciding step

//         // Decide tentative next point
//         double step = base_step;
//         array_1d<double,3> p_next;

//         for (int backoff = 0; backoff < 8; ++backoff) {
//             p_next = p + dir * (step + eps); // small extension to include boundary
//             SizeType jx = compute_index(p_next[0], min_u, max_u, span_x, NX);
//             SizeType jy = compute_index(p_next[1], min_v, max_v, span_y, NY);

//             // Accept if change is 0 or exactly +/−1 on a single axis
//             const int dx = static_cast<int>(jx) - static_cast<int>(ix);
//             const int dy = static_cast<int>(jy) - static_cast<int>(iy);

//             const bool single_axis_step =
//                 ((dx == 0 && (dy == 1 || dy == -1)) ||
//                  (dy == 0 && (dx == 1 || dx == -1)) ||
//                  (dx == 0 && dy == 0)); // still inside

//             if (single_axis_step) {
//                 // ok, we will move to (jx,jy); seg_end = p_next limited by bbox
//                 // fetch and test conditions in current cell before moving
//                 if (ix < rSkinConditionsPerSpan.NumberOfSpansX && iy < rSkinConditionsPerSpan.NumberOfSpansY) {
//                     const auto& S = rSkinConditionsPerSpan;
//                     const SizeType k = FindNnzIndex(S.Occupancy, ix, iy);
//                     if (k != static_cast<SizeType>(-1)) {
//                         const auto ids = CellIdsByK(S, k);
                
//                         Node start_node(0, seg_start[0], seg_start[1], seg_start[2]);
//                         Node end_node  (0, p_next[0],    p_next[1],    p_next[2]);
                
//                         bool hit = false;
//                         for (SizeType t = 0; t < ids.size; ++t) {
//                             const IndexType cond_id = ids.data[t];
//                             const auto& cond = rSkinSubModelPart.GetCondition(cond_id);
//                             const auto& geo  = cond.GetGeometry();
//                             if (geo.size() < 2) continue;
//                             const Node& a = geo[0];
//                             const Node& b = geo[1];
//                             if (SegmentsIntersect(start_node, end_node, a, b)) {
//                                 update_candidate(a);
//                                 update_candidate(b);
//                                 hit = true;
//                             }
//                         }
//                         if (hit && best_node_id != std::numeric_limits<IndexType>::max()) {
//                             return best_node_id;
//                         }
//                     }
//                 }

//                 // Advance to next cell if changed
//                 seg_start = p_next; // next segment will start here
//                 traveled += step;
//                 p = p_next;

//                 ix = jx; iy = jy;
//                 break; // exit backoff loop and continue marching
//             }

//             // If diagonal or jump > 1, halve the step and retry
//             step *= 0.5;
//         }

//         // If backoff exhausted, take a minimal step to avoid stalling
//         if (!in_bounds(ix, iy)) break;
//         if (step < 1e-6 * base_step) {
//             p_next = p + dir * (1e-6 * base_step + eps);
//             seg_start = p_next;
//             traveled += 1e-6 * base_step;
//             p = p_next;

//             SizeType jx = compute_index(p[0], min_u, max_u, span_x, NX);
//             SizeType jy = compute_index(p[1], min_v, max_v, span_y, NY);
//             ix = jx; iy = jy;
//         }
//     }

//     KRATOS_ERROR << "::[SnakeCutSbmProcess]::FindClosestNodeInLayerWithDirection: no intersection for layer "
//                  << rLayer << " from point " << rStartPoint << std::endl;
//     return 0;
// }


void SnakeCutSbmProcess::FindClosestTruePointToSurrogateVertexByNurbs()
{
    // // get the data
    // mEchoLevel = mThisParameters["echo_level"].GetInt();
    // std::string iga_model_part_name = mThisParameters["model_part_name"].GetString();
    // std::string skin_model_part_name = mThisParameters["skin_model_part_name"].GetString();
    // std::string skin_model_part_inner_initial_name = mThisParameters["skin_model_part_inner_initial_name"].GetString();
    
    // // Loop over the surrogate vertex and call the search in radius or something like that
    // // TODO: Outer boyndary
    // const bool is_inner = true;
    
    // std::string rSurrogateSubModelPartName; 
    // std::string skin_sub_model_part_name; 
    // rSurrogateSubModelPartName = "surrogate_inner";
    // skin_sub_model_part_name = "inner";

    // ModelPart& rCutSbmSubModelPart = mpIgaModelPart->CreateSubModelPart("extended_sbm");

    // ModelPart& rSkinSubModelPart = mpSkinModelPart->GetSubModelPart(skin_sub_model_part_name);
    // ModelPart& rSurrogateSubModelPart = mpIgaModelPart->GetSubModelPart(rSurrogateSubModelPartName);

    // auto p_surface = mpIgaModelPart->pGetGeometry(1);
    // IndexType id_brep_curve_on_surface = (mpIgaModelPart->GeometriesEnd()-1)->Id() + 1;
    
    // typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceType;
    // auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
    //                         p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    // IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // // Get the mesh sizes from the surrogate model part
    // const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    // double knot_span_reference_size = knot_span_sizes[0];
    // if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    // if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    // std::string condition_name = "ExtendedSbmSolidCondition";
    // std::string element_name = "ExtendedSbmSolidElement";
    // // std::string element_name = "SolidElement";

    // // Loop over the nodes of the surrogate sub model part
    // IndexType iel = 1;
    // SizeType number_of_shape_functions_derivatives = 5;

    // IndexType first_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[0].Id();
    // IndexType last_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[1].Id();

    // SizeType size_surrogate_loop = last_condition_id - first_condition_id + 1;

    // for (SizeType j = 0; j < size_surrogate_loop; ++j) {
    //     auto p_brep_geometry = mpIgaModelPart->pGetGeometry(6 + j);
    //     auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

    //     KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeCutSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
    //                                         << " is not a BrepCurveOnSurfaceType." << std::endl;

    //     NurbsInterval brep_domain_interval = p_brep_curve_on_surface_surrogate1_surrogate2->DomainInterval();
    //     CoordinatesArrayType surrogate_vertex_1 = ZeroVector(3); 
    //     CoordinatesArrayType surrogate_vertex_1_local_coords = ZeroVector(3);
    //     CoordinatesArrayType surrogate_vertex_2 = ZeroVector(3); 
    //     CoordinatesArrayType surrogate_vertex_2_local_coords = ZeroVector(3);
    //     surrogate_vertex_1_local_coords[0] = brep_domain_interval.GetT0();
    //     surrogate_vertex_2_local_coords[0] = brep_domain_interval.GetT1();

    //     p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
    //     p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

    //     // retrieve middle point of the brep
    //     CoordinatesArrayType surrogate_middle_point = ZeroVector(3); 
    //     CoordinatesArrayType surrogate_middle_point_local_coords = ZeroVector(3);
    //     surrogate_middle_point_local_coords[0] = 0.5 * (surrogate_vertex_1_local_coords[0] + surrogate_vertex_2_local_coords[0]);
    //     p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_middle_point, surrogate_middle_point_local_coords);

    //     IntegrationPoint<1> integration_point(surrogate_middle_point_local_coords[0]);
    //     IntegrationPointsArrayType surrogate_integration_points_list;
    //     surrogate_integration_points_list.push_back(integration_point);

    //     IntegrationInfo integration_info = p_brep_curve_on_surface_surrogate1_surrogate2->GetDefaultIntegrationInfo();
    //     GeometriesArrayType quadrature_point_list;
    //     p_brep_curve_on_surface_surrogate1_surrogate2->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives, 
    //                                                         surrogate_integration_points_list, integration_info);

    //     GeometryType::Pointer surrogate_brep_middle_geometry = quadrature_point_list(0);

    //     CoordinatesArrayType skin_vertex_1 = ZeroVector(3);
    //     CoordinatesArrayType skin_vertex_local_1 = ZeroVector(3);
    //     std::vector<array_1d<double, 3>> curve_derivatives;
    //     bool is_projected = this->ProjectToSkinBoundary(&*mpSkinModelPartInnerInitial, surrogate_vertex_1, skin_vertex_1, skin_vertex_local_1, curve_derivatives, 50);
        
    //     // second vertex
    //     CoordinatesArrayType skin_vertex_2 = ZeroVector(3);
    //     CoordinatesArrayType skin_vertex_local_2 = ZeroVector(3);
    //     std::vector<array_1d<double, 3>> curve_derivatives_2;
    //     bool is_projected_2 = this->ProjectToSkinBoundary(&*mpSkinModelPartInnerInitial, surrogate_vertex_2, skin_vertex_2, skin_vertex_local_2, curve_derivatives_2, 50);

    //     Vector active_range_knot_vector = ZeroVector(2);
    //     active_range_knot_vector[0] = 0;
    //     active_range_knot_vector[1] = 1;
    //     NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

    //     // surrogate_1 - skin_1
    //     // create the brep connecting vertex and closest true point
    //     Point surrogate_1(surrogate_vertex_1);
    //     Point skin_1(skin_vertex_1);
        
    //     Node::Pointer p_surrogate1_brep_point = Node::Pointer(new Node(1, surrogate_1));
    //     Node::Pointer p_skin1_brep_point = Node::Pointer(new Node(2, skin_1));

    //     auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_surrogate1_brep_point, p_skin1_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
    
    //     // IntegrationInfo brep_integration_info_surrogate1_skin1 = p_brep_curve_surrogate1_skin1->GetDefaultIntegrationInfo();
    //     IntegrationPointsArrayType brep_integration_points_list_surrogate1_skin1;
    //     GeometriesArrayType brep_quadrature_point_list_surrogate1_skin1;

    //     p_brep_curve_surrogate1_skin1->CreateIntegrationPoints(brep_integration_points_list_surrogate1_skin1, integration_info);
    //     p_brep_curve_surrogate1_skin1->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate1_skin1, number_of_shape_functions_derivatives, 
    //                                                         brep_integration_points_list_surrogate1_skin1, integration_info);

                                                            
        
    //     SizeType id = 1;
    //     if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
    //         id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

    //     this->CreateConditions(
    //         brep_quadrature_point_list_surrogate1_skin1.ptr_begin(), brep_quadrature_point_list_surrogate1_skin1.ptr_end(),
    //         *mpIgaModelPart, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, "");
            

    //     // surrogate_2 - skin_2
    //     // create the brep connecting vertex and closest true point
    //     Point surrogate_2(surrogate_vertex_2);
    //     Point skin_2(skin_vertex_2);
        
    //     Node::Pointer p_surrogate2_brep_point = Node::Pointer(new Node(1, surrogate_2));
    //     Node::Pointer p_skin2_brep_point = Node::Pointer(new Node(2, skin_2));

    //     auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(p_surrogate2_brep_point, p_skin2_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);      
    
    //     IntegrationPointsArrayType brep_integration_points_list_surrogate2_skin2;
    //     GeometriesArrayType brep_quadrature_point_list_surrogate2_skin2;

    //     p_brep_curve_surrogate2_skin2->CreateIntegrationPoints(brep_integration_points_list_surrogate2_skin2, integration_info);
    //     p_brep_curve_surrogate2_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate2_skin2, number_of_shape_functions_derivatives, 
    //                                                 brep_integration_points_list_surrogate2_skin2, integration_info);
        
    //     if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
    //         id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

    //     this->CreateConditions(
    //         brep_quadrature_point_list_surrogate2_skin2.ptr_begin(), brep_quadrature_point_list_surrogate2_skin2.ptr_end(),
    //         *mpIgaModelPart, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, "");


    //     // skin_1 - skin_2
    //     // create the brep connecting vertex and closest true point
    //     auto p_nurbs_curve_geometry = mpSkinModelPartInnerInitial->pGetGeometry(0);
    //     auto nurbs_curve_geometry = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_nurbs_curve_geometry);
    //     auto t0 = nurbs_curve_geometry->DomainInterval().GetT0();
    //     auto t1 = nurbs_curve_geometry->DomainInterval().GetT1();
    //     CoordinatesArrayType skin_middle_point_local_coords = ZeroVector(3);
    //     skin_middle_point_local_coords = (skin_vertex_local_1+skin_vertex_local_2)*0.5;

    //     if (skin_vertex_local_1[0] < skin_vertex_local_2[0]) 
    //         skin_middle_point_local_coords[0] = 0.5*((t1-t0) + skin_vertex_local_1[0]+skin_vertex_local_2[0]);

    //         if (skin_middle_point_local_coords[0] > t1) 
    //             skin_middle_point_local_coords[0] -= t1-t0;
    //         else if (skin_middle_point_local_coords[0] > rSkinSubModelPart.Nodes().front().Id())
    //             skin_middle_point_local_coords[0] += t1-t0;

    //     CoordinatesArrayType skin_middle_point = ZeroVector(3);
    //     nurbs_curve_geometry->GlobalCoordinates(skin_middle_point, skin_middle_point_local_coords);
    //     Point p_skin_middle_point(skin_middle_point);
    //     Node::Pointer p_skin_middle_brep_point = Node::Pointer(new Node(2, p_skin_middle_point));

    //     PointerVector<Node> ctrl_pts;
    //     Vector               knots;
    //     Vector               weights;

    //     BuildParabolicNurbsData(p_skin1_brep_point, p_skin_middle_brep_point, p_skin2_brep_point,
    //                             ctrl_pts, knots, weights, 0.5);

    //     const double polynomial_degree = 2;
    //     typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_nurbs_curve_skin1_skin2(
    //             new NurbsCurveGeometry<3, PointerVector<Node>>(
    //                 ctrl_pts,
    //                 polynomial_degree,
    //                 knots,
    //                 weights)); 

    //     // auto p_nurbs_curve_skin1_skin2 = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_skin1_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2);      
    
    //     IntegrationPointsArrayType brep_integration_points_list_skin1_skin2;
    //     GeometriesArrayType brep_quadrature_point_list_skin1_skin2;

    //     p_brep_curve_skin1_skin2->CreateIntegrationPoints(brep_integration_points_list_skin1_skin2, integration_info);

    //     const double p_brep_curve_skin1_skin2_length = (p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX))->Length();
    //     for (auto& integration_point : brep_integration_points_list_skin1_skin2) {
    //         integration_point.SetWeight(integration_point.Weight() * p_brep_curve_skin1_skin2_length);
    //     }

    //     p_brep_curve_skin1_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin1_skin2, number_of_shape_functions_derivatives, 
    //                                                               brep_integration_points_list_skin1_skin2, integration_info);
        
        
    //     if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
    //         id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;


    //     std::string identifier_name = "TRUE_BOUNDARY";
    //     this->CreateConditions(
    //         brep_quadrature_point_list_skin1_skin2.ptr_begin(), brep_quadrature_point_list_skin1_skin2.ptr_end(),
    //         r_cut_sbm_sub_model_part, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, identifier_name);

        
    //     // surrogate_1 - surrogate_2
    //     // create ONLY the brep connecting the two vertices
    //     auto p_nurbs_curve_surrogate1_surrogate2 = this->CreateBrepCurve(p_surrogate1_brep_point, p_surrogate2_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_surrogate1_surrogate2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_surrogate2);      


        
    //     // Creation integration points on the cut elements

    //     IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
    //                         3, /*Order*/
    //                         *p_brep_curve_surrogate1_surrogate2,   // B0
    //                         *p_brep_curve_surrogate1_skin1,       // L0
    //                         *p_brep_curve_surrogate2_skin2,       // L1
    //                         *p_brep_curve_skin1_skin2,            // B1
    //                         surrogate_1,  // P00
    //                         skin_1,       // P01
    //                         surrogate_2,  // P10
    //                         skin_2);      // P11
        
    //     GeometriesArrayType surface_quadrature_point_list;

    //     p_nurbs_surface->CreateQuadraturePointGeometries(surface_quadrature_point_list, number_of_shape_functions_derivatives, 
    //                                                     surface_integration_points, surface_integration_info);

    //     IndexType id_element = 1;
    //     if (mpIgaModelPart->GetRootModelPart().Elements().size() > 0)
    //         id_element = mpIgaModelPart->GetRootModelPart().Elements().back().Id() + 1;

    //     this->CreateElements(
    //         surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
    //         *mpIgaModelPart, element_name, id_element, PropertiesPointerType(), surrogate_brep_middle_geometry);
    
    // }
    
}

}  // namespace Kratos.
