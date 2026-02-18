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
#include "snake_sbm_process.h"
#include "iga_application_variables.h"

namespace Kratos
{
SnakeSbmProcess::SnakeSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    Process(),
    mpModel(&rModel),
    mThisParameters(ThisParameters)
{
    mThisParameters.AddMissingParameters(this->GetDefaultParameters());

    mEchoLevel = mThisParameters["echo_level"].GetInt();
    mLambdaInner = mThisParameters["lambda_inner"].GetDouble();
    mLambdaOuter = mThisParameters["lambda_outer"].GetDouble();
    mNumberOfInnerLoops = mThisParameters["number_of_inner_loops"].GetInt();
    mNumberInitialPointsIfImportingNurbs = mThisParameters["number_initial_points_if_importing_nurbs"].GetInt();
    mCreateSurrOuterFromSurrInner = mThisParameters["create_surr_outer_from_surr_inner"].GetBool();

    std::string iga_model_part_name = mThisParameters["model_part_name"].GetString();
    std::string skin_model_part_inner_initial_name = mThisParameters["skin_model_part_inner_initial_name"].GetString();
    std::string skin_model_part_outer_initial_name = mThisParameters["skin_model_part_outer_initial_name"].GetString();
    std::string skin_model_part_name = mThisParameters["skin_model_part_name"].GetString();

    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(iga_model_part_name)) << "::[SnakeSbmProcess]::" 
                    << "Model Part \"iga_model_part\" does not exist. "<< std::endl;  
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(skin_model_part_inner_initial_name)) << "::[SnakeSbmProcess]::" 
                    << "Model Part \"skin_model_part_inner_initial\" does not exist. "<< std::endl;  
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(skin_model_part_outer_initial_name)) << "::[SnakeSbmProcess]::" 
                    << "Model Part \"skin_model_part_outer_initial\" does not exist. "<< std::endl;  
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(skin_model_part_name)) << "::[SnakeSbmProcess]::" 
                    << "Model Part \"skin_model_part\" does not exist. "<< std::endl;  

    mpIgaModelPart = &(mpModel->GetModelPart(iga_model_part_name));
    mpSkinModelPartInnerInitial = &(mpModel->GetModelPart(skin_model_part_inner_initial_name));
    mpSkinModelPartOuterInitial = &(mpModel->GetModelPart(skin_model_part_outer_initial_name));
    mpSkinModelPart = &(mpModel->GetModelPart(skin_model_part_name));

}


void SnakeSbmProcess::CreateTheSnakeCoordinates(bool RemoveIslands)
{   
    // Initilize the property of skin_model_part_in and out
    // skin model part may have nodes if imported from an stl file or geometries if imported from a nurbs file
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0 || mpSkinModelPartInnerInitial->NumberOfGeometries()>0) 
    {
        if (!mpSkinModelPartInnerInitial->HasProperties(0)) mpSkinModelPartInnerInitial->CreateNewProperties(0);
        if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
        // template argument IsInnerLoop set true
        CreateTheSnakeCoordinates<true>(*mpSkinModelPartInnerInitial, mNumberOfInnerLoops, mLambdaInner, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart, mNumberInitialPointsIfImportingNurbs, RemoveIslands);
            
    }

    if (mCreateSurrOuterFromSurrInner) {
        GenerateOuterInitialFromSurrogateInner();
    } else {
        // Normal case sbm from imported skin
        if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) {
            if (!mpSkinModelPartOuterInitial->HasProperties(0)) mpSkinModelPartOuterInitial->CreateNewProperties(0);
            if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
            // template argument IsInnerLoop set false
            CreateTheSnakeCoordinates<false>(*mpSkinModelPartOuterInitial, 1, mLambdaOuter, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart, mNumberInitialPointsIfImportingNurbs, false);
        }
    }
}    



void SnakeSbmProcess::GenerateOuterInitialFromSurrogateInner()
{
    ModelPart* pSkinModelPartOuterInitialFromOuter = nullptr;

    const std::string base_outer_name = mThisParameters["skin_model_part_outer_initial_name"].GetString();
    const std::string generated_outer_name = base_outer_name + std::string("_from_surrogate");

    if (mpModel->HasModelPart(generated_outer_name)) {
        pSkinModelPartOuterInitialFromOuter = &mpModel->GetModelPart(generated_outer_name);
        pSkinModelPartOuterInitialFromOuter->Clear();
    } else {
        pSkinModelPartOuterInitialFromOuter = &mpModel->CreateModelPart(generated_outer_name);
    }

    auto& r_root = mpIgaModelPart->GetRootModelPart();
    auto& r_inner = mpSkinModelPart->GetSubModelPart("inner");
    auto& r_out = *pSkinModelPartOuterInitialFromOuter;

    r_out.Clear();
    if (!r_out.HasProperties(0)) r_out.CreateNewProperties(0);
    auto p_props = r_out.pGetProperties(0);

    double cx = 0.0, cy = 0.0; std::size_t n_nodes = 0;
    for (const auto& r_node : r_inner.Nodes()) { cx += r_node.X(); cy += r_node.Y(); ++n_nodes; }
    if (n_nodes > 0) { cx /= static_cast<double>(n_nodes); cy /= static_cast<double>(n_nodes); }

    IndexType next_node_id = r_root.NumberOfNodes() + 1;
    IndexType next_cond_id = r_root.NumberOfConditions() + 1;

    const auto& knot_u = mpIgaModelPart->GetValue(KNOT_VECTOR_U);
    const auto& knot_v = mpIgaModelPart->GetValue(KNOT_VECTOR_V);
    double step_u = 1.0, step_v = 1.0;
    if (knot_u.size() >= 2) {
        const std::size_t iu = static_cast<std::size_t>(std::ceil(knot_u.size() / 2.0));
        if (iu + 1 < knot_u.size()) step_u = std::abs(knot_u[iu + 1] - knot_u[iu]);
        else step_u = std::abs(knot_u[iu] - knot_u[iu - 1]);
    }
    if (knot_v.size() >= 2) {
        const std::size_t iv = static_cast<std::size_t>(std::ceil(knot_v.size() / 2.0));
        if (iv + 1 < knot_v.size()) step_v = std::abs(knot_v[iv + 1] - knot_v[iv]);
        else step_v = std::abs(knot_v[iv] - knot_v[iv - 1]);
    }

    const double s = 2.0 * std::max(step_u, step_v);

    std::vector<std::array<double,4>> segs;
    segs.reserve(r_inner.NumberOfConditions());
    for (const auto& r_cond : r_inner.Conditions()) {
        const auto& g = r_cond.GetGeometry();
        if (g.size() != 2) continue;
        segs.push_back({g[0].X(), g[0].Y(), g[1].X(), g[1].Y()});
    }

    if (!segs.empty()) {
        struct Point2D {
            double x;
            double y;
        };

        std::vector<Point2D> candidate_points;
        candidate_points.reserve(segs.size());

        // First, compute the offset points without creating entities in the model part
        for (std::size_t i = 0; i < segs.size(); ++i) {
            const double x0 = segs[i][0];
            const double y0 = segs[i][1];
            const double x1 = segs[i][2];
            const double y1 = segs[i][3];

            double tx = x1 - x0; double ty = y1 - y0;
            const double tlen = std::hypot(tx, ty);
            if (tlen <= 1e-15) continue;
            tx /= tlen; ty /= tlen;

            double nx = -ty, ny = tx;
            const double mx = 0.5 * (x0 + x1);
            const double my = 0.5 * (y0 + y1);
            const double wx = mx - cx;
            const double wy = my - cy;
            if (nx * wx + ny * wy < 0.0) { nx = -nx; ny = -ny; }

            const double ax = x0 + s * nx;
            const double ay = y0 + s * ny;
            const double bx = x1 + s * nx;
            const double by = y1 + s * ny;

            if (i == 0) {
                candidate_points.push_back({ax, ay});
                candidate_points.push_back({bx, by});
            } else if (i + 1 < segs.size()) {
                candidate_points.push_back({bx, by});
            } else {
                // last segment: closing will be handled below
            }
        }

        if (!candidate_points.empty()) {

            auto polygon_area = [](const std::vector<Point2D>& r_points) {
                const std::size_t n = r_points.size();
                if (n < 3) {
                    return 0.0;
                }
                double area = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    const auto& p = r_points[i];
                    const auto& q = r_points[(i + 1) % n];
                    area += p.x * q.y - q.x * p.y;
                }
                return 0.5 * area;
            };

            std::vector<Point2D> polygon_points = candidate_points;

            // Build a convex hull of the offset polyline to remove possible
            // self-intersections that can appear for concave shapes.
            if (candidate_points.size() >= 3) {
                std::vector<Point2D> pts = candidate_points;
                std::sort(pts.begin(), pts.end(), [](const Point2D& a, const Point2D& b) {
                    if (a.x < b.x) return true;
                    if (a.x > b.x) return false;
                    return a.y < b.y;
                });
                pts.erase(std::unique(pts.begin(), pts.end(),
                    [](const Point2D& a, const Point2D& b) {
                        return std::abs(a.x - b.x) < 1e-12 && std::abs(a.y - b.y) < 1e-12;
                    }), pts.end());

                if (pts.size() >= 3) {
                    auto cross = [](const Point2D& O, const Point2D& A, const Point2D& B) {
                        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
                    };

                    std::vector<Point2D> hull;
                    hull.reserve(pts.size() * 2);

                    // Lower hull
                    for (const auto& p : pts) {
                        while (hull.size() >= 2 &&
                               cross(hull[hull.size() - 2], hull[hull.size() - 1], p) <= 0.0) {
                            hull.pop_back();
                        }
                        hull.push_back(p);
                    }

                    // Upper hull
                    const std::size_t lower_size = hull.size();
                    for (std::size_t i = pts.size(); i-- > 0;) {
                        const auto& p = pts[i];
                        while (hull.size() > lower_size &&
                               cross(hull[hull.size() - 2], hull[hull.size() - 1], p) <= 0.0) {
                            hull.pop_back();
                        }
                        hull.push_back(p);
                    }

                    if (!hull.empty()) {
                        // Last point is equal to the first one
                        hull.pop_back();
                    }

                    if (hull.size() >= 3) {
                        polygon_points = hull;

                        // Preserve the original orientation (clockwise/anticlockwise)
                        const double area_candidate = polygon_area(candidate_points);
                        const double area_hull = polygon_area(polygon_points);
                        if (area_candidate * area_hull < 0.0) {
                            std::reverse(polygon_points.begin(), polygon_points.end());
                        }
                    }
                }
            }

            if (polygon_points.size() >= 2) {
                const IndexType first_start_id = next_node_id++;
                r_out.CreateNewNode(first_start_id, polygon_points[0].x, polygon_points[0].y, 0.0);
                IndexType last_id = first_start_id;

                for (std::size_t i = 1; i < polygon_points.size(); ++i) {
                    const IndexType node_id = next_node_id++;
                    r_out.CreateNewNode(node_id, polygon_points[i].x, polygon_points[i].y, 0.0);
                    r_out.CreateNewCondition("LineCondition2D2N", next_cond_id++, {{last_id, node_id}}, p_props);
                    last_id = node_id;
                }

                if (polygon_points.size() > 2) {
                    r_out.CreateNewCondition("LineCondition2D2N", next_cond_id++, {{last_id, first_start_id}}, p_props);
                }
            }
        }
    }

    if (pSkinModelPartOuterInitialFromOuter->NumberOfNodes()>0 || pSkinModelPartOuterInitialFromOuter->NumberOfGeometries()>0) {
        if (!pSkinModelPartOuterInitialFromOuter->HasProperties(0)) pSkinModelPartOuterInitialFromOuter->CreateNewProperties(0);
        if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
        CreateTheSnakeCoordinates<false>(*pSkinModelPartOuterInitialFromOuter, 1, mLambdaOuter, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart, mNumberInitialPointsIfImportingNurbs, false);
    }
}

template <bool TIsInnerLoop>
void SnakeSbmProcess::CreateTheSnakeCoordinates(
    const ModelPart& rSkinModelPartInitial,
    const std::size_t NumberOfLoops,
    const double Lambda,
    IndexType EchoLevel,
    ModelPart& rIgaModelPart,
    ModelPart& rSkinModelPart,
    const int NumberInitialPointsIfImportingNurbs,
    bool RemoveIslands) 
{ 
    KRATOS_ERROR_IF(rIgaModelPart.GetValue(KNOT_VECTOR_U).size() == 0) << "::[SnakeSbmProcess]::" 
                << "The iga model part has KNOT_VECTOR_U of size 0" << std::endl;
    KRATOS_ERROR_IF(rIgaModelPart.GetValue(KNOT_VECTOR_V).size() == 0) << "::[SnakeSbmProcess]::" 
                << "The iga model part has KNOT_VECTOR_V of size 0" << std::endl;
    
    Vector knot_vector_u = rIgaModelPart.GetValue(KNOT_VECTOR_U);
    Vector knot_vector_v = rIgaModelPart.GetValue(KNOT_VECTOR_V);
    
    const bool is_inner = TIsInnerLoop;

    std::string surrogate_sub_model_part_name; 
    std::string skin_sub_model_part_name; 
    // ModelPart skin_sub_model_part; 
    if (is_inner)  {
        surrogate_sub_model_part_name = "surrogate_inner";
        skin_sub_model_part_name = "inner";
    }
    else {
        surrogate_sub_model_part_name = "surrogate_outer";
        skin_sub_model_part_name = "outer";
    }
    
    ModelPart& r_skin_sub_model_part = rSkinModelPart.GetSubModelPart(skin_sub_model_part_name);
    ModelPart& r_surrogate_sub_model_part = rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name);

    array_1d<double, 2> knot_step_uv(2);
    knot_step_uv[0] = std::abs(knot_vector_u[std::ceil(knot_vector_u.size()/2) +1]  - knot_vector_u[std::ceil(knot_vector_u.size()/2)] ) ;
    knot_step_uv[1] = std::abs(knot_vector_v[std::ceil(knot_vector_v.size()/2) +1]  - knot_vector_v[std::ceil(knot_vector_v.size()/2)] ) ;

    Vector mesh_sizes_uv(2);
    mesh_sizes_uv[0] = knot_step_uv[0]; 
    mesh_sizes_uv[1] = knot_step_uv[1];
    auto& r_surrogate_model_part = rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name);
    // Note that in here we are saving the knot span info in the parent model part database
    r_surrogate_model_part.GetParentModelPart().SetValue(KNOT_SPAN_SIZES, mesh_sizes_uv);

    array_1d<double, 2> starting_pos_uv;
    starting_pos_uv[0] = knot_vector_u[0];
    starting_pos_uv[1] = knot_vector_v[0];

    std::vector<Vector> parameter_external_coordinates(2);
    parameter_external_coordinates[0].resize(2);
    parameter_external_coordinates[1].resize(2);

    parameter_external_coordinates[0][0] = knot_vector_u[0];
    parameter_external_coordinates[1][0] = knot_vector_v[0];
    parameter_external_coordinates[0][1] = knot_vector_u[knot_vector_u.size()-1];
    parameter_external_coordinates[1][1] = knot_vector_v[knot_vector_v.size()-1];
    // Note that in here we are saving the knot span info in the parent model part database
    r_surrogate_model_part.GetParentModelPart().SetValue(PARAMETER_SPACE_CORNERS, parameter_external_coordinates);
    r_skin_sub_model_part.SetValue(PARAMETER_SPACE_CORNERS, parameter_external_coordinates);

    // Create the matrix of active/inactive knot spans, one for inner and one for outer loop
    std::vector<int> n_knot_spans_uv(2);
    n_knot_spans_uv[0] = knot_vector_u.size()-1; 
    n_knot_spans_uv[1] = knot_vector_v.size()-1;

    std::vector<std::vector<std::vector<int>>> knot_spans_available;
    knot_spans_available.reserve(NumberOfLoops);

    for (IndexType i = 0; i < NumberOfLoops; ++i) {
        std::vector<std::vector<int>> matrix; 
        matrix.reserve(n_knot_spans_uv[1]);
        for (int j = 0; j <= n_knot_spans_uv[1]-1; ++j) {
            std::vector<int> row(n_knot_spans_uv[0]); 
            matrix.push_back(row); 
        }
        knot_spans_available.push_back(matrix);

        // create one Skin subModelPart for each loop
        r_skin_sub_model_part.CreateSubModelPart(std::to_string(i)); 
    }
    
    // Optimized Snake -> for inner loops
    int id_matrix_knot_spans_available = 0;
    IndexType id_first_node;
    bool new_inner_loop = true;
    bool is_skin_nurbs = false;
    
    if (EchoLevel >  0)
    {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::",  is_inner) << "Inner :: Starting SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Starting SnakeStep" << std::endl;
    }
            
    if (rSkinModelPartInitial.NumberOfConditions() > 0) {

        auto p_skin_sub_model_part_loop = &(r_skin_sub_model_part.GetSubModelPart("0"));
        
        for (auto& r_condition : rSkinModelPartInitial.Conditions()) {
            if (new_inner_loop) {
                const std::string loop_sub_model_part_name = std::to_string(id_matrix_knot_spans_available);
                p_skin_sub_model_part_loop = &(r_skin_sub_model_part.GetSubModelPart(loop_sub_model_part_name));
                id_first_node = r_condition.GetGeometry()[0].Id();

                const double x_true_boundary0_loop = r_condition.GetGeometry()[0].X();
                const double y_true_boundary0_loop = r_condition.GetGeometry()[0].Y();
                const int id_new_node = rSkinModelPart.GetRootModelPart().NumberOfNodes() + 1;

                p_skin_sub_model_part_loop->CreateNewNode(id_new_node, x_true_boundary0_loop, y_true_boundary0_loop, 0.0);
                new_inner_loop = false;
            }
            // Collect the coordinates of the points of the current condition
            const auto& r_coords_true_boundary_1 = r_condition.GetGeometry()[0].Coordinates();
            const auto& r_coords_true_boundary_2 = r_condition.GetGeometry()[1].Coordinates();

            std::vector<std::vector<double>> condition_coordinates(2);
            condition_coordinates[0].resize(2);
            condition_coordinates[1].resize(2); 
            
            condition_coordinates[0][0] = r_coords_true_boundary_1[0];
            condition_coordinates[1][0] = r_coords_true_boundary_1[1];
            condition_coordinates[0][1] = r_coords_true_boundary_2[0];
            condition_coordinates[1][1] = r_coords_true_boundary_2[1];
            
            // Collect the intersections of the skin boundary with the knot values
            std::vector<std::vector<int>> knot_span_uv(2);
            knot_span_uv[0].resize(2);
            knot_span_uv[1].resize(2);

            knot_span_uv[0][0] = (r_coords_true_boundary_1[0] - starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_1st_point
            knot_span_uv[1][0] = (r_coords_true_boundary_1[1] - starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_1st_point
            knot_span_uv[0][1] = (r_coords_true_boundary_2[0] - starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_2nd_point
            knot_span_uv[1][1] = (r_coords_true_boundary_2[1] - starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_2nd_point

            // In the inner case : check is the immersed object is inside the rectangular domain
            if (is_inner && IsInside(knot_span_uv, n_knot_spans_uv))
                KRATOS_ERROR << "[SnakeSbmProcess]:: The skin boundary provided is bigger than the background geometry in the parameter space." << std::endl;
            
            // In the outer case : additional check knot_span_uv computation on the domain border 
            if (!is_inner)
            {
                if (knot_span_uv[0][0] == n_knot_spans_uv[0]) knot_span_uv[0][0]--;
                if (knot_span_uv[1][0] == n_knot_spans_uv[1]) knot_span_uv[1][0]--;
                if (knot_span_uv[0][1] == n_knot_spans_uv[0]) knot_span_uv[0][1]--;
                if (knot_span_uv[1][1] == n_knot_spans_uv[1]) knot_span_uv[1][1]--;
            }
            
            SnakeStep(id_matrix_knot_spans_available, knot_span_uv, condition_coordinates, knot_step_uv, starting_pos_uv, 
                        *p_skin_sub_model_part_loop, knot_spans_available);
            
            if (r_condition.GetGeometry()[1].Id() == id_first_node) {
                id_matrix_knot_spans_available++;
                new_inner_loop = true;
            }
        }
    }
    else if (rSkinModelPartInitial.Geometries().size()>0) // if the skin model part is defined by nurbs geometries
    {
        is_skin_nurbs = true;
        // number of sampling points per curve side
        const int number_initial_points_if_importing_nurbs = NumberInitialPointsIfImportingNurbs; 
        int first_node_id = r_skin_sub_model_part.GetRootModelPart().NumberOfNodes()+1;
        const std::size_t n_boundary_curves = rSkinModelPartInitial.NumberOfGeometries();

        // Reorder curves to form a single closed loop: each curve's start must match previous curve's end (within tol)
        const double tol = 1e-7;
        std::vector<IndexType> ordered_indices;
        ordered_indices.reserve(n_boundary_curves);
        std::vector<bool> used(n_boundary_curves, false);

        // Precompute start/end coordinates for each curve (local t=0 and t=1)
        std::vector<CoordinatesArrayType> starts(n_boundary_curves), ends(n_boundary_curves);
        for (IndexType i = 0; i < n_boundary_curves; ++i) {
            auto p_geom = rSkinModelPartInitial.pGetGeometry(i);
            NurbsCurveGeometryPointerType p_curve_i = std::dynamic_pointer_cast<Kratos::NurbsCurveGeometry<2, Kratos::PointerVector<Kratos::Node>>>(p_geom);
            KRATOS_ERROR_IF_NOT(p_curve_i) << "NURBS curve " << i << " not defined in the initial Model Part. Check the importNurbsSbmModeler." << std::endl;
            CoordinatesArrayType local0 = ZeroVector(3);
            CoordinatesArrayType local1 = ZeroVector(3); local1[0] = 1.0;
            starts[i].resize(3,false); ends[i].resize(3,false);
            p_curve_i->GlobalCoordinates(starts[i], local0);
            p_curve_i->GlobalCoordinates(ends[i], local1);
        }

        // Greedy ordering: start from 0 and find next whose start matches current end
        ordered_indices.push_back(0);
        used[0] = true;
        CoordinatesArrayType current_end = ends[0];
        for (IndexType k = 1; k < n_boundary_curves; ++k) {
            bool found = false;
            for (IndexType j = 1; j < n_boundary_curves; ++j) { // j=0 already used
                if (used[j]) continue;
                if (norm_2(current_end - starts[j]) <= tol) {
                    ordered_indices.push_back(j);
                    used[j] = true;
                    current_end = ends[j];
                    found = true;
                    break;
                }
            }
            KRATOS_ERROR_IF_NOT(found)
                << "[SnakeSbmProcess] Could not find the next NURBS curve to continue the single closed loop."
                << " Ensure curves are connected head-to-tail and no reversal is needed." << std::endl;
        }
        // Check closure: last end must match first start
        KRATOS_ERROR_IF(norm_2(current_end - starts[ordered_indices.front()]) > tol)
            << "[SnakeSbmProcess] The ordered NURBS curves do not form a closed loop (last end != first start)." << std::endl;

        bool new_inner_loop = true;

        //Create vertices sub model part
        std::string interface_sub_model_part_name = "interface_vertices";
        ModelPart& r_skin_interface_sub_model_part = r_skin_sub_model_part.HasSubModelPart(interface_sub_model_part_name) ? 
                                                    r_skin_sub_model_part.GetSubModelPart(interface_sub_model_part_name) : r_skin_sub_model_part.CreateSubModelPart(interface_sub_model_part_name);
                                                         
        for (IndexType i_ordered = 0; i_ordered < n_boundary_curves; i_ordered++) 
        {
            const IndexType i_boundary_curve = ordered_indices[i_ordered];
            NurbsCurveGeometryPointerType p_curve = std::dynamic_pointer_cast<Kratos::NurbsCurveGeometry<2, Kratos::PointerVector<Kratos::Node>>>(rSkinModelPartInitial.pGetGeometry(i_boundary_curve));
            if (!p_curve) 
                KRATOS_ERROR << "NURBS curve " << i_boundary_curve << " not defined in the initial Model Part. Check the importNurbsSbmModeler." << std::endl;

            // first point
            CoordinatesArrayType first_point_coords(3);
            Vector first_point_local_coord = ZeroVector(3);
            p_curve->GlobalCoordinates(first_point_coords, first_point_local_coord);

            std::string layer_name = p_curve->GetValue(IDENTIFIER);
            ModelPart& r_skin_layer_sub_model_part = r_skin_sub_model_part.HasSubModelPart(layer_name) ? 
                                                     r_skin_sub_model_part.GetSubModelPart(layer_name) : 
                                                     r_skin_sub_model_part.CreateSubModelPart(layer_name);

            // check the first point of the curve
            if (new_inner_loop) 
            {

                Node::Pointer p_node = new Node(first_node_id, first_point_coords[0], first_point_coords[1], first_point_coords[2]);

                // Create two nodes and two conditions for each skin condition

                //needed for the call to the assign_vector_variable_to_nodes_process
                // compute normal at the node coords
                std::vector<CoordinatesArrayType> global_space_derivatives;
                std::size_t derivative_order = 2;
                CoordinatesArrayType new_point_local_coord = ZeroVector(3); //first point at local coord zero
                p_curve->GlobalSpaceDerivatives(global_space_derivatives, new_point_local_coord, derivative_order);
                CoordinatesArrayType tangent_vector = global_space_derivatives[1];
                double tangent_magnitude = norm_2(tangent_vector);
                tangent_vector /= tangent_magnitude;
                Vector normal_vector = ZeroVector(3);
                normal_vector[0] = tangent_vector[1];
                normal_vector[1] = -tangent_vector[0];

                p_node->SetValue(NORMAL, normal_vector);
                p_node->SetValue(LOCAL_TANGENT, tangent_vector);

                // compute the curvature
                CoordinatesArrayType curve_first_derivative_vector = global_space_derivatives[1];
                CoordinatesArrayType curve_second_derivative_vector = global_space_derivatives[2];

                double curvature = norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector, curve_second_derivative_vector)) / pow(norm_2(curve_first_derivative_vector), 3);
                p_node->SetValue(CURVATURE, curvature);
        
                r_skin_layer_sub_model_part.AddNode(p_node);
                
                // add to the interface sub model part
                r_skin_interface_sub_model_part.AddNode(p_node);
                // cut sbm modifications 
                std::string condition_name = p_curve->GetValue(CONDITION_NAME);
                auto connected_layers = p_node->GetValue(CONNECTED_LAYERS);
                auto connected_condition_names = p_node->GetValue(CONNECTED_CONDITIONS);
                connected_layers.push_back(layer_name);
                connected_condition_names.push_back(condition_name);

                p_node->SetValue(CONNECTED_LAYERS, connected_layers);
                p_node->SetValue(CONNECTED_CONDITIONS, connected_condition_names);

                new_inner_loop = false;

            } else {
                const int last_node_id = r_skin_sub_model_part.GetRootModelPart().NumberOfNodes();
                Node& r_last_node = r_skin_sub_model_part.GetNode(last_node_id);
                KRATOS_ERROR_IF(norm_2(first_point_coords - r_last_node) > tol)
                                << "[SnakeSbmProcess] NURBS curves reordering failed: expected continuity between curves but points differ by > tol." << std::endl;


                // Create two nodes and two conditions for each skin condition
                // compute normal at the node coords
                std::vector<CoordinatesArrayType> global_space_derivatives;
                std::size_t derivative_order = 2;
                CoordinatesArrayType new_point_local_coord = ZeroVector(3); //first point at local coord zero
                p_curve->GlobalSpaceDerivatives(global_space_derivatives, new_point_local_coord, derivative_order);
                CoordinatesArrayType tangent_vector = global_space_derivatives[1];
                double tangent_magnitude = norm_2(tangent_vector);
                tangent_vector /= tangent_magnitude;
                Vector normal_vector = ZeroVector(3);
                normal_vector[0] = tangent_vector[1];
                normal_vector[1] = -tangent_vector[0];

                r_last_node.SetValue(NORMAL, normal_vector);
                r_last_node.SetValue(LOCAL_TANGENT, tangent_vector);
                
                // cut sbm modifications 
                std::string condition_name = p_curve->GetValue(CONDITION_NAME);
                auto connected_layers = r_last_node.GetValue(CONNECTED_LAYERS);
                auto connected_condition_names = r_last_node.GetValue(CONNECTED_CONDITIONS);
                connected_layers.push_back(layer_name);
                connected_condition_names.push_back(condition_name);

                r_last_node.SetValue(CONNECTED_LAYERS, connected_layers);
                r_last_node.SetValue(CONNECTED_CONDITIONS, connected_condition_names);

            }
            // add the specified number of points
            Vector second_point_local_coord = ZeroVector(3);
            CoordinatesArrayType second_point_coords(3);
            for (int i = 1; i < number_initial_points_if_importing_nurbs; i++)
            {
                second_point_local_coord[0] = (double) i/(number_initial_points_if_importing_nurbs-1);
                p_curve->GlobalCoordinates(second_point_coords, second_point_local_coord);
                //***********************************************************
                    // Collect the coordinates of the points
                std::vector<std::vector<double>> xy_coord_i_cond(2);
                xy_coord_i_cond[0].resize(2); xy_coord_i_cond[1].resize(2); 
                
                xy_coord_i_cond[0][0] = first_point_coords[0]; // x_true_boundary1
                xy_coord_i_cond[1][0] = first_point_coords[1]; // y_true_boundary1
                xy_coord_i_cond[0][1] = second_point_coords[0]; // x_true_boundary2
                xy_coord_i_cond[1][1] = second_point_coords[1]; // y_true_boundary2
                
                // Collect the intersections of the skin boundary with the knot values
                std::vector<std::vector<int>> knot_span_uv(2);
                knot_span_uv[0].resize(2); knot_span_uv[1].resize(2);

                knot_span_uv[0][0] = (first_point_coords[0]-starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_1st_point
                knot_span_uv[1][0] = (first_point_coords[1]-starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_1st_point
                knot_span_uv[0][1] = (second_point_coords[0]-starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_2nd_point
                knot_span_uv[1][1] = (second_point_coords[1]-starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_2nd_point

                if (is_inner &&
                            (knot_span_uv[0][0] < 0 || knot_span_uv[0][0] >= n_knot_spans_uv[0] ||
                            knot_span_uv[1][0] < 0 || knot_span_uv[1][0] >= n_knot_spans_uv[1] ||
                            knot_span_uv[0][1] < 0 || knot_span_uv[0][1] >= n_knot_spans_uv[0] ||
                            knot_span_uv[1][1] < 0 || knot_span_uv[1][1] >= n_knot_spans_uv[1]) )
                    KRATOS_ERROR << "[SnakeSbmUtilities]:: The skin boundary provided is bigger than the background geometry in the parameter space." << std::endl;

                // additional check knot_span_uv computation on the domain border [especially for outer boundary]
                if (knot_span_uv[0][0] == n_knot_spans_uv[0]) knot_span_uv[0][0]--; 
                if (knot_span_uv[1][0] == n_knot_spans_uv[1]) knot_span_uv[1][0]--;
                if (knot_span_uv[0][1] == n_knot_spans_uv[0]) knot_span_uv[0][1]--; 
                if (knot_span_uv[1][1] == n_knot_spans_uv[1]) knot_span_uv[1][1]--;

                std::vector<double> local_coords{first_point_local_coord[0], second_point_local_coord[0]};
                
                SnakeStepNurbs(id_matrix_knot_spans_available, knot_span_uv, xy_coord_i_cond, knot_step_uv, starting_pos_uv, local_coords,
                                p_curve, r_skin_layer_sub_model_part, knot_spans_available);
                
                first_point_local_coord = second_point_local_coord;
                first_point_coords = second_point_coords;
            }
            // check the last point of the curve
            if (norm_2(second_point_coords - r_skin_sub_model_part.GetNode(first_node_id)) < 1e-15)
            {
                const int last_node_id = r_skin_sub_model_part.GetRootModelPart().NumberOfNodes();
                Condition* p_closing_condition = nullptr;
                for (auto cond_it = r_skin_sub_model_part.ConditionsEnd(); cond_it != r_skin_sub_model_part.ConditionsBegin();) {
                    --cond_it;
                    if (int(cond_it->GetGeometry()[1].Id()) == last_node_id) {
                        p_closing_condition = &(*cond_it);
                        break;
                    }
                }
                KRATOS_ERROR_IF(!p_closing_condition) << "[SnakeSbmProcess] Could not locate closing condition." << std::endl;

                auto p_first_node = r_skin_sub_model_part.pGetNode(first_node_id);
                p_closing_condition->GetGeometry()(1) = p_first_node;

                // cut sbm modifications 
                Node& r_last_node = r_skin_sub_model_part.GetNode(first_node_id);
                std::string layer_name = p_curve->GetValue(IDENTIFIER);
                std::string condition_name = p_curve->GetValue(CONDITION_NAME);
                auto connected_layers = r_last_node.GetValue(CONNECTED_LAYERS);
                auto connected_condition_names = r_last_node.GetValue(CONNECTED_CONDITIONS);
                connected_layers.push_back(layer_name);
                connected_condition_names.push_back(condition_name);

                r_last_node.SetValue(CONNECTED_LAYERS, connected_layers);
                r_last_node.SetValue(CONNECTED_CONDITIONS, connected_condition_names);

                r_skin_sub_model_part.RemoveNodeFromAllLevels(last_node_id);

                first_node_id = r_skin_sub_model_part.GetRootModelPart().NumberOfNodes()+1;
                new_inner_loop = true;
                id_matrix_knot_spans_available++;
            } else // cut sbm modifications
            {
                const int last_node_id = r_skin_sub_model_part.GetRootModelPart().NumberOfNodes();
                Node& r_last_node = r_skin_sub_model_part.GetNode(last_node_id);
                r_skin_interface_sub_model_part.AddNode(&r_last_node);
            }
        }
    }
    else {
        KRATOS_ERROR << "::[SnakeSbmProcess]:: Reference Skin model part for SBM is empty." << std::endl;;
    }

    PointVector points;
    for (auto &i_cond : r_skin_sub_model_part.Conditions()) {
        points.push_back(Kratos::make_intrusive<PointType>(
            i_cond.Id(),
            i_cond.GetGeometry()[0].X(),
            i_cond.GetGeometry()[0].Y(),
            i_cond.GetGeometry()[0].Z()
        ));
    }
    DynamicBins points_bin(points.begin(), points.end());

    if (EchoLevel >  0) {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Ending SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Ending SnakeStep" << std::endl;
        
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Starting MarkKnotSpansAvailable" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Starting MarkKnotSpansAvailable" << std::endl;
    }

    for (IndexType id_inner_loop = 0; id_inner_loop < NumberOfLoops; id_inner_loop++) {
        // Mark the knot_spans_available's for inner and outer loops
        const std::string loop_sub_model_part_name = std::to_string(id_inner_loop);
        auto p_skin_sub_model_part_loop = &(r_skin_sub_model_part.GetSubModelPart(loop_sub_model_part_name));

        if (is_skin_nurbs)
            MarkKnotSpansAvailable(id_inner_loop, points_bin, r_skin_sub_model_part, Lambda, 
                n_knot_spans_uv, knot_step_uv, starting_pos_uv, knot_spans_available);   
        else
            MarkKnotSpansAvailable(id_inner_loop, points_bin, *p_skin_sub_model_part_loop, Lambda, 
                                    n_knot_spans_uv, knot_step_uv, starting_pos_uv, knot_spans_available);       
    
        if (EchoLevel >  0) {
            KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Ending MarkKnotSpansAvailable" << std::endl;
            KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Ending MarkKnotSpansAvailable" << std::endl;
        }

        if (RemoveIslands)
        {
            auto& plane = knot_spans_available[id_inner_loop]; 
            KeepLargestZeroIsland<TIsInnerLoop>(plane);
        }
        
        if (is_inner) {
            if (is_skin_nurbs)
                CreateSurrogateBuondaryFromSnakeInner(id_inner_loop, r_skin_sub_model_part, points_bin, n_knot_spans_uv, 
                    knot_vector_u, knot_vector_v, starting_pos_uv, knot_spans_available, r_surrogate_sub_model_part);
            else
                CreateSurrogateBuondaryFromSnakeInner(id_inner_loop, *p_skin_sub_model_part_loop, points_bin, n_knot_spans_uv, 
                                                        knot_vector_u, knot_vector_v, starting_pos_uv, knot_spans_available, r_surrogate_sub_model_part);
            
            if (EchoLevel >  0)
                KRATOS_INFO("::[SnakeSbmProcess]::") << "Inner :: Snake process has finished" << std::endl;
        }
        else {
            if (is_skin_nurbs)
                CreateSurrogateBuondaryFromSnakeOuter(id_inner_loop, r_skin_sub_model_part, points_bin, n_knot_spans_uv, knot_vector_u,
                                                    knot_vector_v, starting_pos_uv, knot_spans_available, r_surrogate_sub_model_part);
            else
                CreateSurrogateBuondaryFromSnakeOuter(id_inner_loop, *p_skin_sub_model_part_loop, points_bin, n_knot_spans_uv, knot_vector_u,
                                                    knot_vector_v, starting_pos_uv, knot_spans_available, r_surrogate_sub_model_part);
            
            if (EchoLevel >  0)
                KRATOS_INFO("::[SnakeSbmProcess]::") << "Outer :: Snake process has finished" << std::endl;
        }
    }

    if (EchoLevel >  0) {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Loop finished" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Loop finished" << std::endl;
    }
}


void SnakeSbmProcess::SnakeStep(
    const int IdMatrix, 
    const std::vector<std::vector<int>>& rKnotSpansUV, 
    const std::vector<std::vector<double>>& rConditionCoordinates, 
    const array_1d<double, 2>& rKnotStepUV, 
    const array_1d<double, 2>& rStartingPosition,
    ModelPart& rSkinModelPart, 
    std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable)
{
    bool is_splitted = false;

    if (rKnotSpansUV[0][0] != rKnotSpansUV[0][1] || rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
        // Check if we are jumping some cut knot spans. If yes we split the true segment
        if (std::abs(rKnotSpansUV[1][0]-rKnotSpansUV[1][1]) > 1 || std::abs(rKnotSpansUV[0][0]-rKnotSpansUV[0][1]) > 1 || 
                (rKnotSpansUV[0][0] != rKnotSpansUV[0][1] && rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) ) {
            is_splitted = true;

            // Split the segment and do it recursively
            const double x_true_boundary_split = (rConditionCoordinates[0][0] + rConditionCoordinates[0][1]) / 2.0;
            const double y_true_boundary_split = (rConditionCoordinates[1][0] + rConditionCoordinates[1][1]) / 2.0;
            int knot_span_u_point_split = (x_true_boundary_split - rStartingPosition[0]) / rKnotStepUV[0];
            int knot_span_v_point_split = (y_true_boundary_split - rStartingPosition[1]) / rKnotStepUV[1];

            if (knot_span_u_point_split == static_cast<int>(rKnotSpansAvailable[IdMatrix][0].size())) knot_span_u_point_split--;
            if (knot_span_v_point_split == static_cast<int>(rKnotSpansAvailable[IdMatrix].size())) knot_span_v_point_split--;

            // update xy_coord for the first split segment
            std::vector<std::vector<double>> condition_coordinates_split(2);
            condition_coordinates_split[0].resize(2);
            condition_coordinates_split[1].resize(2); 
            condition_coordinates_split[0][0] = rConditionCoordinates[0][0]; // x_true_boundary1
            condition_coordinates_split[1][0] = rConditionCoordinates[1][0]; // y_true_boundary1
            condition_coordinates_split[0][1] = x_true_boundary_split; // x_true_boundary_split
            condition_coordinates_split[1][1] = y_true_boundary_split; // y_true_boundary_split
            // update knot_span_uv for the first split segment
            std::vector<std::vector<int>> knot_span_uv_split(2);
            knot_span_uv_split[0].resize(2);
            knot_span_uv_split[1].resize(2); 
            knot_span_uv_split[0][0] = rKnotSpansUV[0][0]; // knot_span_u_1st_point
            knot_span_uv_split[1][0] = rKnotSpansUV[1][0]; // knot_span_v_1st_point
            knot_span_uv_split[0][1] = knot_span_u_point_split; // knot_span_u_point_split
            knot_span_uv_split[1][1] = knot_span_v_point_split; // knot_span_v_point_split
            
            // __We do it recursively first split__
            SnakeStep(IdMatrix, knot_span_uv_split, condition_coordinates_split, rKnotStepUV, rStartingPosition, 
                        rSkinModelPart, rKnotSpansAvailable);

            // update xy_coord for the second split segment
            condition_coordinates_split[0][0] = x_true_boundary_split; // x_true_boundary_split
            condition_coordinates_split[1][0] = y_true_boundary_split; // y_true_boundary_split
            condition_coordinates_split[0][1] = rConditionCoordinates[0][1]; // x_true_boundary2
            condition_coordinates_split[1][1] = rConditionCoordinates[1][1]; // y_true_boundary2
            // update knot_span_uv for the first split segment
            knot_span_uv_split[0][0] = knot_span_u_point_split; // knot_span_u_point_split
            knot_span_uv_split[1][0] = knot_span_v_point_split; // knot_span_v_point_split
            knot_span_uv_split[0][1] = rKnotSpansUV[0][1]; // knot_span_u_2nd_point
            knot_span_uv_split[1][1] = rKnotSpansUV[1][1]; // knot_span_v_2nd_point

            // __We do it recursively second split__
            SnakeStep(IdMatrix, knot_span_uv_split, condition_coordinates_split, rKnotStepUV, rStartingPosition, 
                        rSkinModelPart, rKnotSpansAvailable);
        }
        // Check if the true boundary crosses an u or a v knot value
        else if (rKnotSpansUV[0][0] != rKnotSpansUV[0][1]) { // u knot value is crossed
            // Find the "rKnotSpansAvailable" using the intersection
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][1]] = 2;
        }
        else if (rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) { // v knot value is crossed
            // Find the "rKnotSpansAvailable" using the intersection (Snake_coordinate classic -> External Boundary)
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][1]][rKnotSpansUV[0][0]] = 2;
        }
    }
    if (!is_splitted) {
        // Call the root model part for the Ids of the node
        const IndexType node_id_1 = (rSkinModelPart.GetRootModelPart().NodesEnd() - 1)->Id();
        const IndexType node_id_2 = node_id_1 + 1;
        // Create two nodes and two conditions for each skin condition
        rSkinModelPart.CreateNewNode(
            node_id_2,
            (rConditionCoordinates[0][0] + rConditionCoordinates[0][1]) / 2.0,
            (rConditionCoordinates[1][0] + rConditionCoordinates[1][1]) / 2.0,
            0.0);
        rSkinModelPart.CreateNewNode(
            node_id_2 + 1,
            rConditionCoordinates[0][1],
            rConditionCoordinates[1][1],
            0.0);
        auto p_cond_prop = rSkinModelPart.pGetProperties(0);
        auto p_cond1 = rSkinModelPart.CreateNewCondition("LineCondition2D2N", node_id_1, {{node_id_1, node_id_2}}, p_cond_prop );
        auto p_cond2 = rSkinModelPart.CreateNewCondition("LineCondition2D2N", node_id_2, {{node_id_2, node_id_2 + 1}}, p_cond_prop );
        rSkinModelPart.AddCondition(p_cond1);
        rSkinModelPart.AddCondition(p_cond2);
    }
}


void SnakeSbmProcess::SnakeStepNurbs(
            const int IdMatrix, 
            const std::vector<std::vector<int>>& rKnotSpansUV, 
            const std::vector<std::vector<double>>& rConditionCoordinates, 
            const array_1d<double, 2>& rKnotStepUV, 
            const array_1d<double, 2>& rStartingPosition,
            const std::vector<double>& rLocalCoords,
            const NurbsCurveGeometryPointerType& rpCurve,
            ModelPart& rSkinModelPart, 
            std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable)
{
    bool is_splitted = false;

    if (rKnotSpansUV[0][0] != rKnotSpansUV[0][1] || rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) { 
        // intersection between true and surrogate boundary
        // Check if we are jumping some cut knot spans. If yes we split the true segment
        if (std::abs(rKnotSpansUV[1][0]-rKnotSpansUV[1][1]) > 1 || std::abs(rKnotSpansUV[0][0]-rKnotSpansUV[0][1]) > 1 || 
                (rKnotSpansUV[0][0] != rKnotSpansUV[0][1] && rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) ) {
            is_splitted = true;
            
            // Split the segment and do it recursively
            Vector local_coords_split = ZeroVector(3);
            local_coords_split[0] = (rLocalCoords[0] + rLocalCoords[1]) / 2;
            CoordinatesArrayType xy_true_boundary_split;
            rpCurve->GlobalCoordinates(xy_true_boundary_split, local_coords_split);

            int knot_span_u_point_split = (xy_true_boundary_split[0]-rStartingPosition[0]) / rKnotStepUV[0] ;
            int knot_span_v_point_split = (xy_true_boundary_split[1]-rStartingPosition[1]) / rKnotStepUV[1] ;

            // check if it's exactly the same of the first or second point
            bool is_passing_through_diagonal = (std::abs(rKnotSpansUV[0][0] - rKnotSpansUV[0][1]) == 1)  && 
                                                (std::abs(rKnotSpansUV[1][0] - rKnotSpansUV[1][1]) == 1);
            
            if (is_passing_through_diagonal)
            {
                // additional check to avoid infinite loop: check if the splitted segment is too small
                const double split_segment_length = sqrt(std::pow((xy_true_boundary_split[0] - rConditionCoordinates[0][0]),2) + 
                                                         std::pow((xy_true_boundary_split[1] - rConditionCoordinates[1][0]),2)); 
                
                // exactly passing trough a diagonal vertex 
                const double minumum_length = std::min(rKnotStepUV[0]/5000000, rKnotStepUV[1]/5000000); //TODO: to assest in future PR
                if (split_segment_length <= minumum_length)
                {
                    KRATOS_WARNING("[SnakeSbmProcess] :: one skin segment is exactly passing trough a diagonal vertex");
                    // -> we mark an arbitrary knot span (the one with the x of the first point and the y of the second point) as cut.
                    knot_span_u_point_split = rKnotSpansUV[0][0];
                    knot_span_v_point_split = rKnotSpansUV[1][1];


                    // TODO:: check thorugh diagonal. Assest in future PR 
                    rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][0]] = 2;
                    rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][1]] = 2;
                    rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][1]][rKnotSpansUV[0][0]] = 2;
                }
            }

            if (knot_span_u_point_split == int (rKnotSpansAvailable[IdMatrix][0].size())) knot_span_u_point_split--;
            if (knot_span_v_point_split == int (rKnotSpansAvailable[IdMatrix].size())) knot_span_v_point_split--;

            // update xy_coord for the first split segment
            std::vector<std::vector<double>> condition_coordinates_split(2);
            condition_coordinates_split[0].resize(2); condition_coordinates_split[1].resize(2); 
            condition_coordinates_split[0][0] = rConditionCoordinates[0][0];
            condition_coordinates_split[1][0] = rConditionCoordinates[1][0];
            condition_coordinates_split[0][1] = xy_true_boundary_split[0];
            condition_coordinates_split[1][1] = xy_true_boundary_split[1];
            // update knot_span_uv for the first split segment
            std::vector<std::vector<int>> knot_span_uv_split(2);
            knot_span_uv_split[0].resize(2); knot_span_uv_split[1].resize(2); 
            knot_span_uv_split[0][0] = rKnotSpansUV[0][0];
            knot_span_uv_split[1][0] = rKnotSpansUV[1][0];
            knot_span_uv_split[0][1] = knot_span_u_point_split;
            knot_span_uv_split[1][1] = knot_span_v_point_split;

            std::vector<double> local_coords_split_segment1{rLocalCoords[0], local_coords_split[0]};
            
            // __We do it recursively first split__
            SnakeStepNurbs(IdMatrix, knot_span_uv_split, condition_coordinates_split, rKnotStepUV, rStartingPosition, local_coords_split_segment1,
                        rpCurve, rSkinModelPart, rKnotSpansAvailable);

            // update xy_coord for the second split segment
            condition_coordinates_split[0][0] = xy_true_boundary_split[0];
            condition_coordinates_split[1][0] = xy_true_boundary_split[1];
            condition_coordinates_split[0][1] = rConditionCoordinates[0][1];
            condition_coordinates_split[1][1] = rConditionCoordinates[1][1];
            // update knot_span_uv for the second split segment
            knot_span_uv_split[0][0] = knot_span_u_point_split;
            knot_span_uv_split[1][0] = knot_span_v_point_split;
            knot_span_uv_split[0][1] = rKnotSpansUV[0][1];
            knot_span_uv_split[1][1] = rKnotSpansUV[1][1];

            std::vector<double> local_coords_split_segment2{local_coords_split[0], rLocalCoords[1]};

            // __We do it recursively second split__
            SnakeStepNurbs(IdMatrix, knot_span_uv_split, condition_coordinates_split, rKnotStepUV, rStartingPosition, local_coords_split_segment2,
                        rpCurve, rSkinModelPart, rKnotSpansAvailable);
        }
        // Check if the true boundary crosses an u or a v knot value
        else if (rKnotSpansUV[0][0] != rKnotSpansUV[0][1]) { // u knot value is crossed
            // Find the "knot_spans_available" using the intersection
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][1]] = 2;

        }
        else if (rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) { // v knot value is crossed
            // Find the "knot_spans_available" using the intersection (Snake_coordinate classic -> External Boundary)
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][0]][rKnotSpansUV[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUV[1][1]][rKnotSpansUV[0][0]] = 2;
        }
    }
    if (!is_splitted) {
        // Call the root model part for the Ids of the node
        const IndexType node_id_1 = (rSkinModelPart.GetRootModelPart().NodesEnd() - 1)->Id();

        std::string layer_name = rpCurve->GetValue(IDENTIFIER);
        std::string condition_name = rpCurve->GetValue(CONDITION_NAME);

        // check if we are jumping to the next layer. If yes we have to dublicate the skin node
        if (!rSkinModelPart.HasNode(node_id_1))
        {
            auto p_node_from_root = rSkinModelPart.GetRootModelPart().pGetNode(node_id_1);

            rSkinModelPart.AddNode(p_node_from_root);
        }

        
        const IndexType node_id_2 = node_id_1 + 1;
        // Create two nodes and two conditions for each skin condition
        auto p_new_node = new Node(node_id_2, rConditionCoordinates[0][1], rConditionCoordinates[1][1], 0.0);
        // compute normal and tangent informations at the local coord of the point 
        std::vector<CoordinatesArrayType> global_space_derivatives;
        std::size_t derivative_order = 2;
        CoordinatesArrayType new_point_local_coord = ZeroVector(3);
        new_point_local_coord[0] = rLocalCoords[1];
        rpCurve->GlobalSpaceDerivatives(global_space_derivatives, new_point_local_coord, derivative_order);
        CoordinatesArrayType tangent_vector = global_space_derivatives[1];
        double tangent_magnitude = norm_2(tangent_vector);
        tangent_vector /= tangent_magnitude;
        Vector normal_vector = ZeroVector(3);
        normal_vector[0] = tangent_vector[1];
        normal_vector[1] = -tangent_vector[0];
        p_new_node->SetValue(NORMAL, normal_vector);
        p_new_node->SetValue(LOCAL_TANGENT, tangent_vector);


        // compute the curvature
        CoordinatesArrayType curve_first_derivative_vector = global_space_derivatives[1];
        CoordinatesArrayType curve_second_derivative_vector = global_space_derivatives[2];

        double curvature = norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector, curve_second_derivative_vector)) / pow(norm_2(curve_first_derivative_vector), 3);
        p_new_node->SetValue(CURVATURE, curvature);

        //cut sbm modifications
        auto connected_layers = p_new_node->GetValue(CONNECTED_LAYERS);
        auto connected_condition_names = p_new_node->GetValue(CONNECTED_CONDITIONS);
        connected_layers.push_back(layer_name);
        connected_condition_names.push_back(condition_name);

        p_new_node->SetValue(CONNECTED_LAYERS, connected_layers);
        p_new_node->SetValue(CONNECTED_CONDITIONS, connected_condition_names);

        rSkinModelPart.AddNode(p_new_node);

        Properties::Pointer p_cond_prop = rSkinModelPart.pGetProperties(0);
        Condition::Pointer p_cond = rSkinModelPart.CreateNewCondition("LineCondition2D2N", node_id_1, {{node_id_1, node_id_2}}, p_cond_prop );

        p_cond->SetValue(CONDITION_NAME, condition_name);
        p_cond->SetValue(LAYER_NAME, layer_name);
        // In NURBS case, tag the skin condition with its originating BREP_ID
        if (rpCurve) {
            if (rpCurve->Has(BREP_ID)) {
                p_cond->SetValue(BREP_ID, rpCurve->GetValue(BREP_ID));
            }
        }
        rSkinModelPart.AddCondition(p_cond);
    }
}


bool SnakeSbmProcess::IsPointInsideSkinBoundary(
    const Point& rPoint1, 
    DynamicBins& rPointsBin, 
    const ModelPart& rSkinModelPart)
{
    // Get the nearest point of the true boundary
    DynamicBinsPointerType p_point_to_search = DynamicBinsPointerType(new PointType(1, rPoint1.X(), rPoint1.Y(), 0.0));
    DynamicBinsPointerType p_nearest_point = rPointsBin.SearchNearestPoint(*p_point_to_search);
    
    // Get the closest Condition the initial_skin_model_part_in.Conditions
    IndexType id_1 = p_nearest_point->Id();
    const IndexType first_condition_id = rSkinModelPart.ConditionsBegin()->Id();
    const IndexType number_conditions = rSkinModelPart.NumberOfConditions();

    auto compute_cross_product_z = [&](IndexType condition_id) {
        const auto& r_condition = rSkinModelPart.GetCondition(condition_id);
        const auto& r_third_point = r_condition.GetGeometry()[1].Coordinates();

        array_1d<double,3> v_1;
        array_1d<double,3> v_2;

        v_1 = r_condition.GetGeometry()[0] - rPoint1;
        v_2 = r_third_point - rPoint1;

        array_1d<double,3> cross_product;
        MathUtils<double>::CrossProduct(cross_product, v_1, v_2);
        return cross_product[2];
    };

    const double cross_product_main_z = compute_cross_product_z(id_1);

    IndexType id_2;
    if (id_1 == first_condition_id) {
        id_2 = static_cast<IndexType>(first_condition_id + number_conditions - 1);
    } else {
        id_2 = id_1 - 1;
    }

    const double cross_product_previous_z = compute_cross_product_z(id_2);

    if (cross_product_main_z * cross_product_previous_z < 0.0) {
        const auto& r_condition_main = rSkinModelPart.GetCondition(id_1);
        const auto& r_condition_previous = rSkinModelPart.GetCondition(id_2);

        Node temp_point(0, rPoint1.X(), rPoint1.Y(), rPoint1.Z());
        const Node& r_main_second_node = r_condition_main.GetGeometry()[1];
        const Node& r_prev_first_node = r_condition_previous.GetGeometry()[0];
        const Node& r_prev_second_node = r_condition_previous.GetGeometry()[1];

        const double extension_length = 1e10;

        const auto& r_main_first_node = r_condition_main.GetGeometry()[0];
        array_1d<double, 3> main_extended_coordinates = r_main_second_node.Coordinates();
        array_1d<double, 3> main_direction = r_main_second_node - r_main_first_node;
        const double main_direction_norm = norm_2(main_direction);
        if (main_direction_norm != 0.0) {
            main_extended_coordinates += (main_direction / main_direction_norm) * extension_length;
        }
        const Node main_extended_node(0, main_extended_coordinates[0], main_extended_coordinates[1], main_extended_coordinates[2]);

        array_1d<double, 3> prev_extended_coordinates = r_prev_first_node.Coordinates();
        array_1d<double, 3> prev_direction = r_prev_first_node - r_prev_second_node;
        const double prev_direction_norm = norm_2(prev_direction);
        if (prev_direction_norm != 0.0) {
            prev_extended_coordinates += (prev_direction / prev_direction_norm) * extension_length;
        }
        const Node prev_extended_node(0, prev_extended_coordinates[0], prev_extended_coordinates[1], prev_extended_coordinates[2]);
        
        if (SegmentsIntersect(temp_point, main_extended_node, prev_extended_node, r_prev_second_node)) {
            return cross_product_previous_z > 0.0;
        }
    }

    return cross_product_main_z > 0.0;
}


/**
    * Marking process:
    *   1) We set to 2 all the cut knot spans
    *   2) We check the 8 neighbor knot spans and set them to 1|0  if inside|outside
    *   3) We check all the cut knot spans and set them    to 1|-1 if inside|outside 
    */
void SnakeSbmProcess::MarkKnotSpansAvailable(
    const int IdMatrix,
    DynamicBins& rPointsBin, 
    const ModelPart& rSkinModelPart,
    const double Lambda, 
    const std::vector<int>& rNumberKnotSpans, 
    const array_1d<double, 2>& rKnotStepUV,
    const Vector& rStartingPosition,
    std::vector<std::vector<std::vector<int>>> & rKnotSpansAvailable) 
{
    for (int i = 0; i < rNumberKnotSpans[1]; i++) {
        for (int j = 0; j < rNumberKnotSpans[0]; j++) {
            if (rKnotSpansAvailable[IdMatrix][i][j] == 2) {
                // Check the 8 neighbor knot spans -> Is there any completely inside? Note that we can just check 1 point.
                
                // right node
                if (i != rNumberKnotSpans[1]-1)
                    if (rKnotSpansAvailable[IdMatrix][i+1][j] == 0) { 
                        Point gauss_point = Point((j+0.5) * rKnotStepUV[0] + rStartingPosition[0], (i+1+0.5) * rKnotStepUV[1] +rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i+1][j] = 1;}
                    }
                // left node    
                if (i != 0)
                    if (rKnotSpansAvailable[IdMatrix][i-1][j] == 0) { 
                        Point gauss_point = Point((j+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i-1+0.5) * rKnotStepUV[1] + rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i-1][j] = 1;}
                    }
                // up node
                if (j != rNumberKnotSpans[0]-1)
                    if (rKnotSpansAvailable[IdMatrix][i][j+1] == 0) { 
                        Point gauss_point = Point((j+1+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i+0.5) * rKnotStepUV[1]+rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i][j+1] = 1;}
                    }
                //down node
                if (j != 0)
                    if (rKnotSpansAvailable[IdMatrix][i][j-1] == 0) { 
                        Point gauss_point = Point((j-1+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i+0.5) * rKnotStepUV[1]+rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i][j-1] = 1;}
                    } 

                // corner right-down node
                if (j != 0 && i != rNumberKnotSpans[1]-1)
                    if (rKnotSpansAvailable[IdMatrix][i+1][j-1] == 0) {
                        Point gauss_point = Point((j-1+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i+1+0.5) * rKnotStepUV[1]+rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i+1][j-1] = 1;}
                    }
                // corner left-down node
                if (j != 0 && i != 0)
                    if (rKnotSpansAvailable[IdMatrix][i-1][j-1] == 0) {
                        Point gauss_point = Point((j-1+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i-1+0.5) * rKnotStepUV[1]+rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i-1][j-1] = 1;}
                    }
                // corner right-up node
                if (j != rNumberKnotSpans[0]-1 && i != rNumberKnotSpans[1]-1)
                    if (rKnotSpansAvailable[IdMatrix][i+1][j+1] == 0) {
                        Point gauss_point = Point((j+1+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i+1+0.5) * rKnotStepUV[1]+rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i+1][j+1] = 1;}
                    }
                // corner left-up node
                if (j != rNumberKnotSpans[0]-1 && i != 0)
                    if (rKnotSpansAvailable[IdMatrix][i-1][j+1] == 0) {
                        Point gauss_point = Point((j+1+0.5) * rKnotStepUV[0]+rStartingPosition[0], (i-1+0.5) * rKnotStepUV[1]+rStartingPosition[1], 0);
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {rKnotSpansAvailable[IdMatrix][i-1][j+1] = 1;}
                    }

                // Create 49 "fake" gauss_points to check if the majority are inside or outside
                const int num_fake_gauss_points = 7;
                int number_of_inside_gaussian_points = 0;
                const double tollerance = rKnotStepUV[0]/1e8; // Tolerance to avoid numerical issues
                for (IndexType i_GPx = 0; i_GPx < num_fake_gauss_points; i_GPx++){
                    double x_coord = (j*rKnotStepUV[0]+tollerance) +
                                     (rKnotStepUV[0]-2*tollerance)/(num_fake_gauss_points-1)*(i_GPx) 
                                     + rStartingPosition[0];

                    // NOTE:: The v-knot spans are upside down in the matrix!!
                    for (IndexType i_GPy = 0; i_GPy < num_fake_gauss_points; i_GPy++) 
                    {
                        double y_coord = (i*rKnotStepUV[1]+tollerance) + 
                                         (rKnotStepUV[1]-2*tollerance)/(num_fake_gauss_points-1)*(i_GPy) 
                                        + rStartingPosition[1];
                        Point gauss_point = Point(x_coord, y_coord, 0);  // GAUSSIAN POINT
                        if (IsPointInsideSkinBoundary(gauss_point, rPointsBin, rSkinModelPart)) {
                            // Sum over the number of num_fake_gauss_points per knot span
                            number_of_inside_gaussian_points++;
                        }
                    }
                    
                }
            
                // Mark the knot span as available or not depending on the number of Gauss Points Inside/Outside
                if (number_of_inside_gaussian_points < Lambda*num_fake_gauss_points*num_fake_gauss_points) {
                    rKnotSpansAvailable[IdMatrix][i][j] = -1; // Cut knot spans that have been checked
                }
                else{
                    rKnotSpansAvailable[IdMatrix][i][j] = 1; // The knot span is considered DEACTIVE
                }
            }
        }
    }
}

/**
    * summary of knot_spans_available:
        " 1"  -> interior knot spans                                  
        "-1"  -> exterior knot spans well checked
        " 0"  -> exterior knot spans OR very interior knot spans (more 
                    than one ks away from surrogate boundary)
    */
void SnakeSbmProcess::CreateSurrogateBuondaryFromSnakeInner(
    const int IdMatrix, 
    const ModelPart& rSkinModelPartInner, 
    DynamicBins& rPointsBinInner,
    const std::vector<int>& rNumberKnotSpans, 
    const Vector& rKnotVectorU, 
    const Vector& rKnotVectorV,
    const Vector& rStartingPositionUV,
    std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable,
    ModelPart& rSurrogateModelPartInner
    ) 
{
    // Snake 2D works with a raycasting technique from each of the two directions

    const double knot_step_u = rKnotVectorU[1]-rKnotVectorU[0];
    const double knot_step_v = rKnotVectorV[1]-rKnotVectorV[0];
    
    IndexType id_surrogate_first_node; 
    if (rSurrogateModelPartInner.NumberOfNodes() == 0)
    {
        id_surrogate_first_node = rSurrogateModelPartInner.GetRootModelPart().NumberOfNodes() + 1;
        IndexType idSurrogateNode = id_surrogate_first_node;
        for (int j = 0; j < rNumberKnotSpans[1]; j++) {
            for (int i = 0; i < rNumberKnotSpans[0]; i++) {
                rSurrogateModelPartInner.CreateNewNode(idSurrogateNode, rKnotVectorU[i], rKnotVectorV[j], 0.0);
                idSurrogateNode++;
            }
        }
    } else 
    {
        id_surrogate_first_node = rSurrogateModelPartInner.GetRootModelPart().NumberOfNodes() - rNumberKnotSpans[1]*rNumberKnotSpans[0] + 1;
    }
    
    auto p_cond_prop = rSurrogateModelPartInner.pGetProperties(0);
    
    // Direction parallel to x
    IndexType id_surrogate_condition = rSurrogateModelPartInner.GetRootModelPart().NumberOfConditions() + 1;
    IndexType id_surrogate_first_condition = id_surrogate_condition;
    for (int j = 0; j < rNumberKnotSpans[1]; j++) {
        bool check_next_point = false;
        /*  
            Formula to connect i,j to the id of the model_part
            id = id_surrogate_first_node + [i + j*(n_knot_spans_uv[0]) + 1];
        */
        // move in the x direction
        for (int i = 0; i < rNumberKnotSpans[0]; i++) {
            if (check_next_point) {
                // Check i+1 point using isPointInsideSkinBoundary3D
                Point knot_span_center_point = Point(rStartingPositionUV[0] + (i + 0.5)*knot_step_u, rStartingPositionUV[1] + (j + 0.5)*knot_step_v, 0.0);
                bool is_exiting = false;
                if ( rKnotSpansAvailable[IdMatrix][j][i] == 1 ) {
                    // the knot span was already been checked very well
                    continue;
                }
                else if (IsPointInsideSkinBoundary(knot_span_center_point, rPointsBinInner, rSkinModelPartInner)) {
                    // STILL INSIDE --> do not save nothing and update rKnotSpansAvailable 
                    if ( rKnotSpansAvailable[IdMatrix][j][i] == -1) {
                        is_exiting = true;
                    }
                    else {
                        rKnotSpansAvailable[IdMatrix][j][i] = 1;
                    }
                }
                else {
                    is_exiting = true;
                }
                if (is_exiting) {
                    /* EXITING --> save last segment in direction x. i-th is the knot value. */
                    int node1_i = i; int node1_j = j;   
                    int node2_i = i; int node2_j = j+1; 

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*rNumberKnotSpans[0];
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*rNumberKnotSpans[0];

                    auto p_condition = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    p_condition->Set(BOUNDARY, false);

                    // surrogate_model_part_inner.AddCondition(p_condition);
                    id_surrogate_condition++;
                    check_next_point = false;
                }
                
            }
            else if (rKnotSpansAvailable[IdMatrix][j][i] == 1) {
                // ENTERING --> save first face in direction
                int node1_i = i; int node1_j = j;   
                int node2_i = i; int node2_j = j+1; 

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*rNumberKnotSpans[0]; 
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*rNumberKnotSpans[0];

                auto p_condition = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                p_condition->Set(BOUNDARY, true);
            }
        }
    }
    
    // Do the same for y direction, without isPointInsideSkinBoundary, since we have already done it
    // And it is not necessary do it again
    for (int i = 0; i < rNumberKnotSpans[0]; i++) {
        
        bool check_next_point = false;
        /*  
            Formula to connect i,j,k to the id of the model_part
            id = id_surrogate_first_node + [i + j*(n_knot_spans_uv[0]) + 1];
        */
        // move in the y direction
        for (int j = 0; j < rNumberKnotSpans[1]; j++) {
            if (check_next_point) {
                if (rKnotSpansAvailable[IdMatrix][j][i] != 1) {
                    /* EXITING --> save last face in direction x. i-th is the knot value. */
                    int node1_i = i;   int node1_j = j;
                    int node2_i = i+1; int node2_j = j;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*rNumberKnotSpans[0];
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*rNumberKnotSpans[0];
  
                    auto p_condition = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    p_condition->Set(BOUNDARY, false);

                    // surrogate_model_part_inner.AddCondition(p_cond);
                    id_surrogate_condition++;
                    check_next_point = false;
                } 
            }
            else if (rKnotSpansAvailable[IdMatrix][j][i] == 1) {
                // ENTERING --> save first face in direction
                int node1_i = i;   int node1_j = j;
                int node2_i = i+1; int node2_j = j;

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]);
                auto p_condition = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                // surrogate_model_part_inner.AddCondition(p_cond);
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                p_condition->Set(BOUNDARY, true);

            }
        }
    }

    // Create "fictituos element" to store starting and ending condition id for each surrogate boundary loop
    IndexType elem_id = rSurrogateModelPartInner.GetRootModelPart().NumberOfElements()+1;
    IndexType id_surrogate_last_condition = id_surrogate_condition-1;
    std::vector<ModelPart::IndexType> elem_nodes{id_surrogate_first_condition, id_surrogate_last_condition};
    rSurrogateModelPartInner.CreateNewElement("Element2D2N", elem_id, elem_nodes, p_cond_prop);
}



void SnakeSbmProcess::CreateSurrogateBuondaryFromSnakeOuter(
    const int IdMatrix, 
    const ModelPart& rSkinModelPartOuter,
    DynamicBins& rPointsBinOuter, 
    const std::vector<int>& rNumberKnotSpans, 
    const Vector& knot_vector_u, 
    const Vector& knot_vector_v, 
    const Vector& rStartingPositionUV,
    std::vector<std::vector<std::vector<int>>> & rKnotSpansAvailable,
    ModelPart& rSurrogateModelPartOuter)
{
    // CHECK ALL THE EXTERNAL KNOT SPANS

    // LEFT BOUNDARY
    double knot_step_u = knot_vector_u[1]-knot_vector_u[0];
    double knot_step_v = knot_vector_v[1]-knot_vector_v[0];

    for (int i = 0; i<2; i++) {
        for (int j = 0; j < (rNumberKnotSpans[0]); j++ ) {
            Point centroidKnotSpan = Point((j+0.5)*knot_step_u+rStartingPositionUV[0], (i+0.5)*knot_step_v+rStartingPositionUV[1], 0);
            if (IsPointInsideSkinBoundary(centroidKnotSpan, rPointsBinOuter, rSkinModelPartOuter) && rKnotSpansAvailable[IdMatrix][i][j] != -1) {
                rKnotSpansAvailable[IdMatrix][i][j] = 1;
                }
        }
    }
    // TOP BOUNDARY
    for (int j = int (rKnotSpansAvailable[IdMatrix][0].size()-1); j > int (rKnotSpansAvailable[IdMatrix][0].size()-3); j--) {
        for (int i = 0; i < (rNumberKnotSpans[1]); i++) {
            Point centroidKnotSpan = Point((j+0.5)*knot_step_u+rStartingPositionUV[0], (i+0.5)*knot_step_v+rStartingPositionUV[1], 0);
            if (IsPointInsideSkinBoundary(centroidKnotSpan, rPointsBinOuter, rSkinModelPartOuter) && rKnotSpansAvailable[IdMatrix][i][j] != -1) {
                rKnotSpansAvailable[IdMatrix][i][j] = 1;
                }
        }
    }
    // RIGHT BOUNDARY
    for (int i = int (rKnotSpansAvailable[IdMatrix].size()-1); i > int (rKnotSpansAvailable[IdMatrix].size()-3); i--) {
        for (int j = rNumberKnotSpans[0]-1; j > -1; j-- ) {
            Point centroidKnotSpan = Point((j+0.5)*knot_step_u+rStartingPositionUV[0], (i+0.5)*knot_step_v+rStartingPositionUV[1], 0);
            if (IsPointInsideSkinBoundary(centroidKnotSpan, rPointsBinOuter, rSkinModelPartOuter) && rKnotSpansAvailable[IdMatrix][i][j] != -1) {
                rKnotSpansAvailable[IdMatrix][i][j] = 1;
                }
        }
    }
    // BOTTOM BOUNDARY
    for (int j = 0; j<2; j++) {
        for (int i = rNumberKnotSpans[1]-1; i > -1 ; i--) {
            Point centroidKnotSpan = Point((j+0.5)*knot_step_u+rStartingPositionUV[0], (i+0.5)*knot_step_v+rStartingPositionUV[1], 0);
            if (IsPointInsideSkinBoundary(centroidKnotSpan, rPointsBinOuter, rSkinModelPartOuter) && rKnotSpansAvailable[IdMatrix][i][j] != -1) {
                rKnotSpansAvailable[IdMatrix][i][j] = 1;
                }
        }
    }
    
    // Snake 2D works with a raycasting technique from each of the two directions
    IndexType id_surrogate_first_node = rSurrogateModelPartOuter.GetRootModelPart().NumberOfNodes() + 1;
    IndexType idSurrogateNode = id_surrogate_first_node;
    for (int j = 0; j < rNumberKnotSpans[1]+1; j++) {
        for (int i = 0; i < rNumberKnotSpans[0]+1; i++) {
            rSurrogateModelPartOuter.CreateNewNode(idSurrogateNode, knot_vector_u[i], knot_vector_v[j], 0.0);
            idSurrogateNode++;
        }
    }
    
    // Direction parallel to x
    
    IndexType id_surrogate_condition = rSurrogateModelPartOuter.GetRootModelPart().NumberOfConditions() + 1;
    auto p_cond_prop = rSurrogateModelPartOuter.pGetProperties(0);
    
    for (int j = 0; j < rNumberKnotSpans[1]; j++) {
        
        bool check_next_point = false;
        /*  
            Formula to connect i,j,k to the id of the model_part
            id = id_surrogate_first_node + [i + j*(n_knot_spans_uv[0]) + 1];
        */
        // move in the x direction
        for (int i = 0; i < rNumberKnotSpans[0]; i++) {

            int node1_i; int node1_j;   
            int node2_i; int node2_j; 

            if (check_next_point) {
                // Check i+1 point using isPointInsideSkinBoundary
                Point knot_span_center_point = Point(rStartingPositionUV[0] + (i + 0.5)*knot_step_u, rStartingPositionUV[1] + (j + 0.5)*knot_step_v, 0.0);
                // FIXME:
                // auto p_center_point = Kratos::make_shared<>();
                bool is_exiting = false;
                node1_i = i; node1_j = j;   
                node2_i = i; node2_j = j+1; 
                if ( rKnotSpansAvailable[IdMatrix][j][i] == 1 ) {
                    // the knot span has already been checked very well
                }
                else if (IsPointInsideSkinBoundary(knot_span_center_point, rPointsBinOuter, rSkinModelPartOuter)) {
                    // STILL INSIDE --> do not save nothing and update knot_spans_available 
                    if ( rKnotSpansAvailable[IdMatrix][j][i] == -1) {
                        is_exiting = true;
                    }
                    else {
                        rKnotSpansAvailable[IdMatrix][j][i] = 1;
                    }
                }
                else {
                    is_exiting = true;
                }
                if (is_exiting) {
                    /* EXITING --> save last face in direction x. i-th is the knot value. */

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                        
                    auto p_condition = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    p_condition->Set(BOUNDARY, false);

                    id_surrogate_condition++;
                    check_next_point = false;    
                }
                
            }
            else if (rKnotSpansAvailable[IdMatrix][j][i] == 1) {
                // ENTERING --> save first face in direction
                int node1_i = i; int node1_j = j;   
                int node2_i = i; int node2_j = j+1; 

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                    
                auto p_condition = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                p_condition->Set(BOUNDARY, true);
            }

            if (rKnotSpansAvailable[IdMatrix][j][i] == 1 && i == rNumberKnotSpans[0]-1) 
            {
                // Check if we are at the end of the patch -> if yes close the surrogate boundary
                int node1_i = i+1; int node1_j = j;   
                int node2_i = i+1; int node2_j = j+1; 

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                    
                auto p_condition = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                p_condition->Set(BOUNDARY, false);
                id_surrogate_condition++;
                check_next_point = false;
            }
        }
    }

    // Do the same for y, without isPointInsideSkinBoundary, since we have already done it
    // And it is not necessary do it again
    for (int i = 0; i < rNumberKnotSpans[0]; i++) {
        
        bool check_next_point = false;
        /*  
            Formula to connect i,j,k to the id of the model_part
            i + j*(nKnotSpansUV[0]) + k*(nKnotSpansUV[1])*(nKnotSpansUV[0]);
        */
        // move in the y direction
        for (int j = 0; j < rNumberKnotSpans[1]; j++) {
            if (check_next_point) {
                int node1_i; int node1_j;
                int node2_i; int node2_j;
                if (rKnotSpansAvailable[IdMatrix][j][i] != 1) {
                    /* EXITING --> save last face in direction x. i-th is the knot value. */
                    node1_i = i;   node1_j = j;
                    node2_i = i+1; node2_j = j;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                        
                    auto p_condition = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    p_condition->Set(BOUNDARY, false);
                    id_surrogate_condition++;
                    check_next_point = false;
                } 
            }
            else if (rKnotSpansAvailable[IdMatrix][j][i] == 1) {
                // ENTERING --> save first face in direction
                int node1_i = i;   int node1_j = j;
                int node2_i = i+1; int node2_j = j;

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                auto p_condition = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                p_condition->Set(BOUNDARY, true);
            }

            if (rKnotSpansAvailable[IdMatrix][j][i] == 1 && j == rNumberKnotSpans[1]-1) 
            {
                // Check if we are at the end of the patch -> if yes close the surrogate boundary
                int node1_i = i;   int node1_j = j+1;   
                int node2_i = i+1; int node2_j = j+1; 

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                    
                auto p_condition = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                p_condition->Set(BOUNDARY, false);
                id_surrogate_condition++;
                check_next_point = false;
            }
        }
    }
}

bool SnakeSbmProcess::IsInside(
    const std::vector<std::vector<int>>& rKnotSpanUV,
    const std::vector<int>& rNumberKnotSpansUV) 
{
    return (rKnotSpanUV[0][0] < 0 || rKnotSpanUV[0][0] >= rNumberKnotSpansUV[0] ||
            rKnotSpanUV[1][0] < 0 || rKnotSpanUV[1][0] >= rNumberKnotSpansUV[1] ||
            rKnotSpanUV[0][1] < 0 || rKnotSpanUV[0][1] >= rNumberKnotSpansUV[0] ||
            rKnotSpanUV[1][1] < 0 || rKnotSpanUV[1][1] >= rNumberKnotSpansUV[1]); 
}

template <bool TIsInnerLoop>
void SnakeSbmProcess::KeepLargestZeroIsland(std::vector<std::vector<int>>& rGrid) 
{
    const int row_count = static_cast<int>(rGrid.size());
    if (row_count == 0) {
        return;
    }
    const int column_count = static_cast<int>(rGrid[0].size());
    if (column_count == 0) {
        return;
    }

    // Label map: -1 = unvisited/non-zero, 0..K = component id for zero components
    std::vector<std::vector<int>> label(row_count, std::vector<int>(column_count, -1));
    std::vector<int> component_size;

    static const int delta_row_8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    static const int delta_col_8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    // Use a snapshot to avoid chain reactions during this pass
    const auto original_grid = rGrid;

    if constexpr (TIsInnerLoop)
    {
        for (int row = 0; row < row_count; ++row) {
            for (int column = 0; column < column_count; ++column) {
                if (original_grid[row][column] != 0) {
                    continue;
                }

                bool has_zero_neighbor = false;
                for (int direction = 0; direction < 8; ++direction) {
                    const int neighbor_row = row + delta_row_8[direction];
                    const int neighbor_column = column + delta_col_8[direction];
                    if (0 <= neighbor_row && neighbor_row < row_count &&
                        0 <= neighbor_column && neighbor_column < column_count &&
                        original_grid[neighbor_row][neighbor_column] == 0) {
                        has_zero_neighbor = true;
                        break;
                    }
                }
                if (!has_zero_neighbor) {
                    rGrid[row][column] = 1; // isolated 0 → 1
                }
            }
        }
    }


    // 4-neighborhood
    const int delta_row_4[4] = {-1, 1, 0, 0};
    const int delta_col_4[4] = { 0, 0,-1, 1};

    int component_id = 0;
    int largest_component_id = -1;
    int largest_component_size = 0;

    for (int row = 0; row < row_count; ++row) {
        for (int column = 0; column < column_count; ++column) {
            if (rGrid[row][column] == 0 && label[row][column] == -1) {
                // BFS to label this zero-component
                std::queue<std::pair<int,int>> queue_zero_nodes;
                queue_zero_nodes.push({row, column});
                label[row][column] = component_id;
                int current_component_size = 0;

                while (!queue_zero_nodes.empty()) {
                    auto [current_row, current_column] = queue_zero_nodes.front();
                    queue_zero_nodes.pop();
                    ++current_component_size;

                    for (int direction = 0; direction < 4; ++direction) {
                        const int neighbor_row = current_row + delta_row_4[direction];
                        const int neighbor_column = current_column + delta_col_4[direction];
                        if (0 <= neighbor_row && neighbor_row < row_count &&
                            0 <= neighbor_column && neighbor_column < column_count &&
                            rGrid[neighbor_row][neighbor_column] == 0 &&
                            label[neighbor_row][neighbor_column] == -1) {
                            label[neighbor_row][neighbor_column] = component_id;
                            queue_zero_nodes.push({neighbor_row, neighbor_column});
                        }
                    }
                }

                component_size.push_back(current_component_size);
                // Track largest; if tie, keep first encountered
                if (current_component_size > largest_component_size) {
                    largest_component_size = current_component_size;
                    largest_component_id = component_id;
                }
                ++component_id;
            }
        }
    }

    if (largest_component_id == -1) {
        return; // no zeros at all
    }

    // Flip all zeros that are NOT in the largest component to 1
    if constexpr (TIsInnerLoop) {
        for (int row = 0; row < row_count; ++row) {
            for (int column = 0; column < column_count; ++column) {
                if (rGrid[row][column] == 0 && label[row][column] != largest_component_id) {
                    rGrid[row][column] = 1;
                }
            }
        }
    } 
    else { 
        for (int row = 0; row < row_count; ++row) {
            for (int column = 0; column < column_count; ++column) {
                if (rGrid[row][column] == 0 && label[row][column] == largest_component_id) {
                    rGrid[row][column] = 1;
                }
            }
        }
    }
}

}  // namespace Kratos.
