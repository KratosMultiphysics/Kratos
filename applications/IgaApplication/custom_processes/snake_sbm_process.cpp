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
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    mEchoLevel = mThisParameters["echo_level"].GetInt();
    mLambdaInner = mThisParameters["lambda_inner"].GetDouble();
    mLambdaOuter = mThisParameters["lambda_outer"].GetDouble();
    mNumberOfInnerLoops = mThisParameters["number_of_inner_loops"].GetInt();

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


void SnakeSbmProcess::CreateTheSnakeCoordinates()
{   
    // Vector know_w = mpIgaModelPart->GetValue(KNOT_VECTOR_W);
    if (mpIgaModelPart->GetValue(KNOT_VECTOR_W).size() == 0) {
        // 2D case
        CreateTheSnakeCoordinates2D();
    } else {
        // 3D case
        CreateTheSnakeCoordinates3D();
    }
}   

void SnakeSbmProcess::CreateTheSnakeCoordinates2D()
{   
    // Initilize the property of skin_model_part_in and out
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0) {
        if (!mpSkinModelPartInnerInitial->HasProperties(0)) mpSkinModelPartInnerInitial->CreateNewProperties(0);
        if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
        // template argument IsInnerLoop set true
        CreateTheSnakeCoordinates<true>(*mpSkinModelPartInnerInitial, mNumberOfInnerLoops, mLambdaInner, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart);
            
    }
    if (mpSkinModelPartOuterInitial->NumberOfNodes()>0) {
        if (!mpSkinModelPartOuterInitial->HasProperties(0)) mpSkinModelPartOuterInitial->CreateNewProperties(0);
        if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
        // template argument IsInnerLoop set false
        CreateTheSnakeCoordinates<false>(*mpSkinModelPartOuterInitial, 1, mLambdaOuter, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart);
    }
}   



template <bool TIsInnerLoop>
void SnakeSbmProcess::CreateTheSnakeCoordinates(
    const ModelPart& rSkinModelPartInitial,
    const std::size_t NumberOfLoops,
    const double Lambda,
    IndexType EchoLevel,
    ModelPart& rIgaModelPart,
    ModelPart& rSkinModelPart) 
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
    auto& surrogate_model_part = rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name);
    surrogate_model_part.GetParentModelPart().SetValue(KNOT_SPAN_SIZES, mesh_sizes_uv);

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
    
    surrogate_model_part.GetParentModelPart().SetValue(PARAMETER_SPACE_CORNERS, parameter_external_coordinates);

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
    }
    
    // Optimized Snake -> for inner loops
    int id_matrix_knot_spans_available = 0;
    IndexType id_first_node;
    bool new_inner_loop = true;
    
    if (EchoLevel >  0)
    {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::",  is_inner) << "Inner :: Starting SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Starting SnakeStep" << std::endl;
    }
    

    KRATOS_WARNING_IF("::[SnakeSbmProcess]::", rSkinModelPartInitial.NumberOfConditions() == 0) 
                    << "Reference Skin model part for SBM has no conditions." << std::endl;
            
    if (rSkinModelPartInitial.NumberOfConditions()> 0) {
        
        // CREATE FIRST NODE FOR SKIN SUB MODEL PART
        auto initial_condition = rSkinModelPartInitial.ConditionsBegin();
        const double x_true_boundary0 = initial_condition->GetGeometry()[0].X();
        const double y_true_boundary0 = initial_condition->GetGeometry()[0].Y();
        
        const int id_new_node = rSkinModelPart.GetRootModelPart().NumberOfNodes()+1; 
        r_skin_sub_model_part.CreateNewNode(id_new_node, x_true_boundary0, y_true_boundary0, 0.0);

        for (auto &i_cond : rSkinModelPartInitial.Conditions()) {  
            if (new_inner_loop) {
                id_first_node = i_cond.GetGeometry()[0].Id();
                new_inner_loop = false;
            }
            // Collect the coordinates of the points of the i_cond
            const auto& r_coords_true_boundary1 = i_cond.GetGeometry()[0].Coordinates();
            const auto& r_coords_true_boundary2 = i_cond.GetGeometry()[1].Coordinates();

            std::vector<std::vector<double>> xy_coord_i_cond(2);
            xy_coord_i_cond[0].resize(2); xy_coord_i_cond[1].resize(2); 
            
            xy_coord_i_cond[0][0] = r_coords_true_boundary1[0];
            xy_coord_i_cond[1][0] = r_coords_true_boundary1[1];
            xy_coord_i_cond[0][1] = r_coords_true_boundary2[0];
            xy_coord_i_cond[1][1] = r_coords_true_boundary2[1];
            
            // Collect the intersections of the skin boundary with the knot values
            std::vector<std::vector<int>> knot_span_uv(2);
            knot_span_uv[0].resize(2); knot_span_uv[1].resize(2);

            knot_span_uv[0][0] = (r_coords_true_boundary1[0]-starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_1st_point
            knot_span_uv[1][0] = (r_coords_true_boundary1[1]-starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_1st_point
            knot_span_uv[0][1] = (r_coords_true_boundary2[0]-starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_2nd_point
            knot_span_uv[1][1] = (r_coords_true_boundary2[1]-starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_2nd_point

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
            
            SnakeStep(id_matrix_knot_spans_available, knot_span_uv, xy_coord_i_cond, knot_step_uv, starting_pos_uv, 
                        r_skin_sub_model_part, knot_spans_available);
            
            if (i_cond.GetGeometry()[1].Id() == id_first_node) {
                id_matrix_knot_spans_available++;
                new_inner_loop = true;
            }
        }
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

    for (IndexType i = 0; i < NumberOfLoops; i++) {
        IndexType id_inner_loop = i;
        // Mark the knot_spans_available's for inner and outer loops
        MarkKnotSpansAvailable(id_inner_loop, points_bin, r_skin_sub_model_part, Lambda, 
                                n_knot_spans_uv, knot_step_uv, starting_pos_uv, knot_spans_available);  
    
        if (EchoLevel >  0) {
            KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Ending MarkKnotSpansAvailable" << std::endl;
            KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Ending MarkKnotSpansAvailable" << std::endl;
        }
        
        if (is_inner) {
            CreateSurrogateBuondaryFromSnakeInner(id_inner_loop, r_skin_sub_model_part, points_bin, n_knot_spans_uv, 
                                                    knot_vector_u, knot_vector_v, knot_spans_available, r_surrogate_sub_model_part);
            
            if (EchoLevel >  0)
                KRATOS_INFO("::[SnakeSbmProcess]::") << "Inner :: Snake process has finished" << std::endl;
        }
        else {
            CreateSurrogateBuondaryFromSnakeOuter (id_inner_loop, r_skin_sub_model_part, points_bin, n_knot_spans_uv, knot_vector_u,
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
    const std::vector<std::vector<double>>& rConditionCoord, 
    const Vector rKnotStepUV, 
    const Vector rStartingPosition,
    ModelPart& rSkinModelPart, 
    std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable)
{
    bool isSplitted = false;

    if (rKnotSpansUV[0][0] != rKnotSpansUV[0][1] || rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
        // Check if we are jumping some cut knot spans. If yes we split the true segment
        if (std::abs(rKnotSpansUV[1][0]-rKnotSpansUV[1][1]) > 1 || std::abs(rKnotSpansUV[0][0]-rKnotSpansUV[0][1]) > 1 || 
                (rKnotSpansUV[0][0] != rKnotSpansUV[0][1] && rKnotSpansUV[1][0] != rKnotSpansUV[1][1]) ) {
            isSplitted = true;

            // Split the segment and do it recursively
            double x_true_boundary_split = (rConditionCoord[0][0]+rConditionCoord[0][1]) / 2;
            double y_true_boundary_split = (rConditionCoord[1][0]+rConditionCoord[1][1]) / 2;
            int knot_span_u_point_split = (x_true_boundary_split-rStartingPosition[0]) / rKnotStepUV[0] ;
            int knot_span_v_point_split = (y_true_boundary_split-rStartingPosition[1]) / rKnotStepUV[1] ;

            if (knot_span_u_point_split == int (rKnotSpansAvailable[IdMatrix][0].size())) knot_span_u_point_split--;
            if (knot_span_v_point_split == int (rKnotSpansAvailable[IdMatrix].size())) knot_span_v_point_split--;

            // update xy_coord for the first split segment
            std::vector<std::vector<double>> xy_coord_i_cond_split(2);
            xy_coord_i_cond_split[0].resize(2); xy_coord_i_cond_split[1].resize(2); 
            xy_coord_i_cond_split[0][0] = rConditionCoord[0][0]; // x_true_boundary1
            xy_coord_i_cond_split[1][0] = rConditionCoord[1][0]; // y_true_boundary1
            xy_coord_i_cond_split[0][1] = x_true_boundary_split; // x_true_boundary_split
            xy_coord_i_cond_split[1][1] = y_true_boundary_split; // y_true_boundary_split
            // update knot_span_uv for the first split segment
            std::vector<std::vector<int>> knot_span_uv_split(2);
            knot_span_uv_split[0].resize(2); knot_span_uv_split[1].resize(2); 
            knot_span_uv_split[0][0] = rKnotSpansUV[0][0]; // knot_span_u_1st_point
            knot_span_uv_split[1][0] = rKnotSpansUV[1][0]; // knot_span_v_1st_point
            knot_span_uv_split[0][1] = knot_span_u_point_split; // knot_span_u_point_split
            knot_span_uv_split[1][1] = knot_span_v_point_split; // knot_span_v_point_split
            
            // __We do it recursively first split__
            SnakeStep(IdMatrix, knot_span_uv_split, xy_coord_i_cond_split, rKnotStepUV, rStartingPosition, 
                        rSkinModelPart, rKnotSpansAvailable);

            // update xy_coord for the second split segment
            xy_coord_i_cond_split[0][0] = x_true_boundary_split; // x_true_boundary_split
            xy_coord_i_cond_split[1][0] = y_true_boundary_split; // y_true_boundary_split
            xy_coord_i_cond_split[0][1] = rConditionCoord[0][1]; // x_true_boundary2
            xy_coord_i_cond_split[1][1] = rConditionCoord[1][1]; // y_true_boundary2
            // update knot_span_uv for the first split segment
            knot_span_uv_split[0][0] = knot_span_u_point_split; // knot_span_u_point_split
            knot_span_uv_split[1][0] = knot_span_v_point_split; // knot_span_v_point_split
            knot_span_uv_split[0][1] = rKnotSpansUV[0][1]; // knot_span_u_2nd_point
            knot_span_uv_split[1][1] = rKnotSpansUV[1][1]; // knot_span_v_2nd_point

            // __We do it recursively second split__
            SnakeStep(IdMatrix, knot_span_uv_split, xy_coord_i_cond_split, rKnotStepUV, rStartingPosition, 
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
    if (!isSplitted) {
        // Call the root model part for the Ids of the node
        auto idNode1 = (rSkinModelPart.GetRootModelPart().NodesEnd()-1)->Id();
        auto idNode2 = idNode1+1;
        // Create two nodes and two conditions for each skin condition
        rSkinModelPart.CreateNewNode(idNode2, (rConditionCoord[0][0]+rConditionCoord[0][1] ) / 2, (rConditionCoord[1][0]+rConditionCoord[1][1] ) / 2, 0.0);
        rSkinModelPart.CreateNewNode(idNode2+1, rConditionCoord[0][1], rConditionCoord[1][1], 0.0);
        auto p_cond_prop = rSkinModelPart.pGetProperties(0);
        auto p_cond1 = rSkinModelPart.CreateNewCondition("LineCondition2D2N", idNode1, {{idNode1, idNode2}}, p_cond_prop );
        auto p_cond2 = rSkinModelPart.CreateNewCondition("LineCondition2D2N", idNode2, {{idNode2, idNode2+1}}, p_cond_prop );
        rSkinModelPart.AddCondition(p_cond1);
        rSkinModelPart.AddCondition(p_cond2);
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
    auto nearest_condition_1 = rSkinModelPart.GetCondition(id_1);
    // Check if the condition is the first one and therefore the previous one does not exist
    IndexType id_2 = id_1 - 1;
    if (id_1 == rSkinModelPart.ConditionsBegin()->Id()) {
        int number_conditions = rSkinModelPart.NumberOfConditions();
        id_2 = id_1 + number_conditions - 1; 
    }
    auto nearest_condition_2 = rSkinModelPart.GetCondition(id_2);
    // The two candidates nodes
    const auto& r_coords_candidate_point_1 = nearest_condition_1.GetGeometry()[1].Coordinates();
    const auto& r_coords_candidate_point_2 = nearest_condition_2.GetGeometry()[0].Coordinates();
    
    array_1d<double,3> v_1;
    array_1d<double,3> v_2;

    if (MathUtils<double>::Norm(r_coords_candidate_point_1-rPoint1) > MathUtils<double>::Norm(r_coords_candidate_point_2-rPoint1)){
        // Need to invert the order to preserve the positivity of the area
        v_1 = r_coords_candidate_point_2 - rPoint1;
        v_2 = nearest_condition_1.GetGeometry()[0] - rPoint1;
    } else 
    {
        v_1 = nearest_condition_1.GetGeometry()[0] - rPoint1;
        v_2 = r_coords_candidate_point_1 - rPoint1;
    }

    array_1d<double,3> cross_product;
    MathUtils<double>::CrossProduct(cross_product, v_1, v_2);

    return cross_product[2] > 0;
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

                // Create 25 "fake" gauss_points to check if the majority are inside or outside
                const int num_fake_gauss_points = 5;
                int number_of_inside_gaussian_points = 0;
                for (IndexType i_GPx = 0; i_GPx < num_fake_gauss_points; i_GPx++){
                    double x_coord = j*rKnotStepUV[0] + rKnotStepUV[0]/(num_fake_gauss_points+1)*(i_GPx+1) + rStartingPosition[0];

                    // NOTE:: The v-knot spans are upside down in the matrix!!
                    for (IndexType i_GPy = 0; i_GPy < num_fake_gauss_points; i_GPy++) 
                    {
                        double y_coord = i*rKnotStepUV[1] + rKnotStepUV[1]/(num_fake_gauss_points+1)*(i_GPy+1) + rStartingPosition[1];
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
    const Vector& knot_vector_u, 
    const Vector& knot_vector_v,
    std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable,
    ModelPart& rSurrogateModelPartInner
    ) 
{
    // Snake 2D works with a raycasting technique from each of the two directions

    const double knot_step_u = knot_vector_u[1]-knot_vector_u[0];
    const double knot_step_v = knot_vector_v[1]-knot_vector_v[0];
    
    IndexType id_surrogate_first_node; 
    if (rSurrogateModelPartInner.NumberOfNodes() == 0)
    {
        id_surrogate_first_node = rSurrogateModelPartInner.GetRootModelPart().NumberOfNodes() + 1;
        IndexType idSurrogateNode = id_surrogate_first_node;
        for (int j = 0; j < rNumberKnotSpans[1]; j++) {
            for (int i = 0; i < rNumberKnotSpans[0]; i++) {
                rSurrogateModelPartInner.CreateNewNode(idSurrogateNode, knot_vector_u[i], knot_vector_v[j], 0.0);
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
                Point centerPoint = Point((i + 0.5)*knot_step_u, (j + 0.5)*knot_step_v, 0.0);
                bool is_exiting = false;
                if ( rKnotSpansAvailable[IdMatrix][j][i] == 1 ) {
                    // the knot span was already been checked very well
                }
                else if (IsPointInsideSkinBoundary(centerPoint, rPointsBinInner, rSkinModelPartInner)) {
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
                        
                    auto pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, false);

                    // surrogate_model_part_inner.AddCondition(pcond);
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
                    
                auto pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                pcond->Set(BOUNDARY, true);
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
                        
                    auto pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, false);

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
                auto pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                // surrogate_model_part_inner.AddCondition(p_cond);
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                pcond->Set(BOUNDARY, true);

            }
        }
    }

    // Create "fictituos element" to store starting and ending condition id for each surrogate boundary loop
    IndexType elem_id = rSurrogateModelPartInner.NumberOfElements()+1;
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
                Point centerPoint = Point((i + 0.5)*knot_step_u, (j + 0.5)*knot_step_v, 0.0);
                // FIXME:
                // auto p_center_point = Kratos::make_shared<>();
                bool is_exiting = false;
                node1_i = i; node1_j = j;   
                node2_i = i; node2_j = j+1; 
                if ( rKnotSpansAvailable[IdMatrix][j][i] == 1 ) {
                    // the knot span has already been checked very well
                }
                else if (IsPointInsideSkinBoundary(centerPoint, rPointsBinOuter, rSkinModelPartOuter)) {
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
                        
                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, false);

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
                    
                auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                pcond->Set(BOUNDARY, true);
            }

            if (rKnotSpansAvailable[IdMatrix][j][i] == 1 && i == rNumberKnotSpans[0]-1) 
            {
                // Check if we are at the end of the patch -> if yes close the surrogate boundary
                int node1_i = i+1; int node1_j = j;   
                int node2_i = i+1; int node2_j = j+1; 

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                    
                auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                pcond->Set(BOUNDARY, false);
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
                        
                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, false);
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
                auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_2, id_node_1}}, p_cond_prop );
                id_surrogate_condition++;
                check_next_point = true;

                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                pcond->Set(BOUNDARY, true);
            }

            if (rKnotSpansAvailable[IdMatrix][j][i] == 1 && j == rNumberKnotSpans[1]-1) 
            {
                // Check if we are at the end of the patch -> if yes close the surrogate boundary
                int node1_i = i;   int node1_j = j+1;   
                int node2_i = i+1; int node2_j = j+1; 

                IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1);
                IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1);
                    
                auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
                // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                pcond->Set(BOUNDARY, false);
                id_surrogate_condition++;
                check_next_point = false;
            }
        }
    }
}

bool SnakeSbmProcess::IsInside(
    const std::vector<std::vector<int>>& rKnotSpanUV,
    const std::vector<int>& NumberKnotSpansUV) 
{
    return (rKnotSpanUV[0][0] < 0 || rKnotSpanUV[0][0] >= NumberKnotSpansUV[0] ||
            rKnotSpanUV[1][0] < 0 || rKnotSpanUV[1][0] >= NumberKnotSpansUV[1] ||
            rKnotSpanUV[0][1] < 0 || rKnotSpanUV[0][1] >= NumberKnotSpansUV[0] ||
            rKnotSpanUV[1][1] < 0 || rKnotSpanUV[1][1] >= NumberKnotSpansUV[1]); 
}


void SnakeSbmProcess::CreateTheSnakeCoordinates3D()
{   
    // Initilize the property of skin_model_part_in and out
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0) {
        if (!mpSkinModelPartInnerInitial->HasProperties(0)) mpSkinModelPartInnerInitial->CreateNewProperties(0);
        if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
        // template argument IsInnerLoop set true
        CreateTheSnakeCoordinates3D<true>(*mpSkinModelPartInnerInitial, mNumberOfInnerLoops, mLambdaInner, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart);
            
    }
    if (mpSkinModelPartOuterInitial->NumberOfNodes()>0) {
        if (!mpSkinModelPartOuterInitial->HasProperties(0)) mpSkinModelPartOuterInitial->CreateNewProperties(0);
        if (!mpSkinModelPart->HasProperties(0)) mpSkinModelPart->CreateNewProperties(0);
        // template argument IsInnerLoop set false
        CreateTheSnakeCoordinates3D<false>(*mpSkinModelPartOuterInitial, 1, mLambdaOuter, mEchoLevel, *mpIgaModelPart, *mpSkinModelPart);
    }
}


template <bool TIsInnerLoop>
void SnakeSbmProcess::CreateTheSnakeCoordinates3D(
    const ModelPart& rSkinModelPartInitial,
    const std::size_t NumberOfLoops,
    const double Lambda,
    IndexType EchoLevel,
    ModelPart& rIgaModelPart,
    ModelPart& rSkinModelPart) 
{ 
    KRATOS_ERROR_IF(rIgaModelPart.GetValue(KNOT_VECTOR_U).size() == 0) << "::[SnakeSbmProcess]::" 
                << "The iga model part has KNOT_VECTOR_U of size 0" << std::endl;
    KRATOS_ERROR_IF(rIgaModelPart.GetValue(KNOT_VECTOR_V).size() == 0) << "::[SnakeSbmProcess]::" 
                << "The iga model part has KNOT_VECTOR_V of size 0" << std::endl;
    KRATOS_ERROR_IF(rIgaModelPart.GetValue(KNOT_VECTOR_W).size() == 0) << "::[SnakeSbmProcess]::" 
                << "The iga model part has KNOT_VECTOR_W of size 0" << std::endl;
    
    Vector knot_vector_u = rIgaModelPart.GetValue(KNOT_VECTOR_U);
    Vector knot_vector_v = rIgaModelPart.GetValue(KNOT_VECTOR_V);
    Vector knot_vector_w = rIgaModelPart.GetValue(KNOT_VECTOR_W);
    
    const bool is_inner = TIsInnerLoop;

    std::string surrogate_sub_model_part_name; 
    std::string skin_sub_model_part_name; 
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

    array_1d<double, 3> knot_step_uvw;
    KRATOS_WATCH(knot_vector_u)
    knot_step_uvw[0] = std::abs(knot_vector_u[std::ceil(knot_vector_u.size()/2) + 1] - knot_vector_u[std::ceil(knot_vector_u.size()/2)]);
    KRATOS_WATCH(knot_vector_v)
    knot_step_uvw[1] = std::abs(knot_vector_v[std::ceil(knot_vector_v.size()/2) + 1] - knot_vector_v[std::ceil(knot_vector_v.size()/2)]);
    KRATOS_WATCH(knot_vector_w)
    knot_step_uvw[2] = std::abs(knot_vector_w[std::ceil(knot_vector_w.size()/2) + 1] - knot_vector_w[std::ceil(knot_vector_w.size()/2)]);

    KRATOS_WATCH(knot_step_uvw)

    KRATOS_WATCH("1")

    // Set KNOT_SPAN_SIZES
    Vector mesh_sizes_uvw(3);
    mesh_sizes_uvw[0] = knot_step_uvw[0];
    mesh_sizes_uvw[1] = knot_step_uvw[1];
    mesh_sizes_uvw[2] = knot_step_uvw[2];

    auto& surrogate_model_part = rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name);
    surrogate_model_part.GetParentModelPart().SetValue(KNOT_SPAN_SIZES, mesh_sizes_uvw);

    // Starting positions
    array_1d<double, 3> starting_pos_uvw;
    starting_pos_uvw[0] = knot_vector_u[0];
    starting_pos_uvw[1] = knot_vector_v[0];
    starting_pos_uvw[2] = knot_vector_w[0];

    // Set PARAMETER_SPACE_CORNERS
    std::vector<Vector> parameter_external_coordinates(3);
    parameter_external_coordinates[0].resize(2); // U
    parameter_external_coordinates[1].resize(2); // V
    parameter_external_coordinates[2].resize(2); // W

    parameter_external_coordinates[0][0] = knot_vector_u[0];
    parameter_external_coordinates[1][0] = knot_vector_v[0];
    parameter_external_coordinates[2][0] = knot_vector_w[0];
    parameter_external_coordinates[0][1] = knot_vector_u[knot_vector_u.size()-1];
    parameter_external_coordinates[1][1] = knot_vector_v[knot_vector_v.size()-1];
    parameter_external_coordinates[2][1] = knot_vector_w[knot_vector_w.size()-1];

    surrogate_model_part.GetParentModelPart().SetValue(PARAMETER_SPACE_CORNERS, parameter_external_coordinates);

    KRATOS_WATCH("2")

    // Create the matrix of active/inactive knot spans, one for inner and one for outer loop
    std::vector<int> n_knot_spans_uvw(3);
    n_knot_spans_uvw[0] = knot_vector_u.size()-1; 
    n_knot_spans_uvw[1] = knot_vector_v.size()-1;
    n_knot_spans_uvw[2] = knot_vector_w.size()-1;

    std::vector<std::vector<std::vector<std::vector<int>>>> knot_spans_available;
    knot_spans_available.resize(NumberOfLoops);

    for (IndexType i = 0; i < NumberOfLoops; ++i) {
        knot_spans_available[i].resize(n_knot_spans_uvw[2]); // W
        for (int j = 0; j < n_knot_spans_uvw[2]; ++j) {
            knot_spans_available[i][j].resize(n_knot_spans_uvw[1]); // V
            for (int k = 0; k < n_knot_spans_uvw[1]; ++k) {
                knot_spans_available[i][j][k].resize(n_knot_spans_uvw[0], 0); // U
            }
        }
    }
    
    int id_matrix_knot_spans_available = 0;
    IndexType id_first_node;
    bool new_inner_loop = true;
    
    if (EchoLevel >  0)
    {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::",  is_inner) << "Inner :: Starting SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Starting SnakeStep" << std::endl;
    }
    
    KRATOS_WARNING_IF("::[SnakeSbmProcess]::", rSkinModelPartInitial.NumberOfConditions() == 0) 
                    << "Reference Skin model part for SBM has no conditions." << std::endl;
    
    KRATOS_WATCH("3")
    if (rSkinModelPartInitial.NumberOfConditions()> 0) {
        
        // Copy all the nodes of the initial_skin_model_part to the skin_model_part
        for (auto &i_node : rSkinModelPartInitial.Nodes()) {
            r_skin_sub_model_part.CreateNewNode(i_node.Id(), i_node.X(), i_node.Y(), i_node.Z());
        } 

        for (auto &i_cond : rSkinModelPartInitial.Conditions()) {  

            if (new_inner_loop) {
                id_first_node = i_cond.GetGeometry()[0].Id();
                new_inner_loop = false;
            }
            // Collect the coordinates of the points of the i_cond (3-node surface condition)
            const auto& r_coords_true_boundary1 = i_cond.GetGeometry()[0].Coordinates();
            const auto& r_coords_true_boundary2 = i_cond.GetGeometry()[1].Coordinates();
            const auto& r_coords_true_boundary3 = i_cond.GetGeometry()[2].Coordinates();

            // Initialize the container [x/y/z][point index]
            std::vector<std::vector<double>> xyz_coord_i_cond(3);
            xyz_coord_i_cond[0].resize(3);  // x
            xyz_coord_i_cond[1].resize(3);  // y
            xyz_coord_i_cond[2].resize(3);  // z

            // Fill values for node 1
            xyz_coord_i_cond[0][0] = r_coords_true_boundary1[0];  // x1
            xyz_coord_i_cond[1][0] = r_coords_true_boundary1[1];  // y1
            xyz_coord_i_cond[2][0] = r_coords_true_boundary1[2];  // z1

            // Fill values for node 2
            xyz_coord_i_cond[0][1] = r_coords_true_boundary2[0];  // x2
            xyz_coord_i_cond[1][1] = r_coords_true_boundary2[1];  // y2
            xyz_coord_i_cond[2][1] = r_coords_true_boundary2[2];  // z2

            // Fill values for node 3
            xyz_coord_i_cond[0][2] = r_coords_true_boundary3[0];  // x3
            xyz_coord_i_cond[1][2] = r_coords_true_boundary3[1];  // y3
            xyz_coord_i_cond[2][2] = r_coords_true_boundary3[2];  // z3

            // Collect the intersections of the skin boundary with the knot values in 3D
            std::vector<std::vector<int>> knot_span_uvw(3);
            knot_span_uvw[0].resize(3); // u indices
            knot_span_uvw[1].resize(3); // v indices
            knot_span_uvw[2].resize(3); // w indices

            // Point 1 (node 0)
            knot_span_uvw[0][0] = (r_coords_true_boundary1[0] - starting_pos_uvw[0]) / knot_step_uvw[0]; // u
            knot_span_uvw[1][0] = (r_coords_true_boundary1[1] - starting_pos_uvw[1]) / knot_step_uvw[1]; // v
            knot_span_uvw[2][0] = (r_coords_true_boundary1[2] - starting_pos_uvw[2]) / knot_step_uvw[2]; // w

            // Point 2 (node 1)
            knot_span_uvw[0][1] = (r_coords_true_boundary2[0] - starting_pos_uvw[0]) / knot_step_uvw[0];
            knot_span_uvw[1][1] = (r_coords_true_boundary2[1] - starting_pos_uvw[1]) / knot_step_uvw[1];
            knot_span_uvw[2][1] = (r_coords_true_boundary2[2] - starting_pos_uvw[2]) / knot_step_uvw[2];

            // Point 3 (node 2)
            knot_span_uvw[0][2] = (r_coords_true_boundary3[0] - starting_pos_uvw[0]) / knot_step_uvw[0];
            knot_span_uvw[1][2] = (r_coords_true_boundary3[1] - starting_pos_uvw[1]) / knot_step_uvw[1];
            knot_span_uvw[2][2] = (r_coords_true_boundary3[2] - starting_pos_uvw[2]) / knot_step_uvw[2];


            // In the inner case : check is the immersed object is inside the rectangular domain
            if (is_inner && IsInside3D(knot_span_uvw, n_knot_spans_uvw))
                KRATOS_ERROR << "[SnakeSbmProcess]:: The skin boundary provided is bigger than the background geometry in the parameter space." << std::endl;
            
            // In the outer case: additional check for knot_span_uvw computation on the domain border
            if (!is_inner)
            {
                for (int point = 0; point < 3; ++point) {
                    for (int d = 0; d < 3; ++d) {
                        if (knot_span_uvw[d][point] == n_knot_spans_uvw[d]) {
                            knot_span_uvw[d][point]--;
                        }
                    }
                }
            }
            array_1d<IndexType, 3> ordered_ids;
            ordered_ids[0] = i_cond.GetGeometry()[0].Id();
            ordered_ids[1] = i_cond.GetGeometry()[1].Id();
            ordered_ids[2] = i_cond.GetGeometry()[2].Id();

            // KRATOS_WATCH("\n \n \n")
            // KRATOS_WATCH("SnakeStep3D")

            SnakeStep3D(id_matrix_knot_spans_available, knot_span_uvw, xyz_coord_i_cond, knot_step_uvw, starting_pos_uvw, 
                        r_skin_sub_model_part, knot_spans_available, ordered_ids);
            
            // KRATOS_WATCH("SnakeStep3D -> end")
            // if (i_cond.GetGeometry()[1].Id() == id_first_node) {
            //     id_matrix_knot_spans_available++;
            //     new_inner_loop = true;
            // }
        }
    }

    KRATOS_WATCH("end snake step")
    KRATOS_WATCH(r_skin_sub_model_part)

    PointVector points;
    for (auto &i_cond : r_skin_sub_model_part.Conditions()) {
        points.push_back(Kratos::make_intrusive<PointType>(
            i_cond.Id(),
            i_cond.GetGeometry().Center().X(),
            i_cond.GetGeometry().Center().Y(),
            i_cond.GetGeometry().Center().Z()
            // i_cond.Id(),
            // i_cond.GetGeometry()[0].X(),
            // i_cond.GetGeometry()[0].Y(),
            // i_cond.GetGeometry()[0].Z()
        ));
    }
    DynamicBins points_bin(points.begin(), points.end());

    if (EchoLevel >  0) {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Ending SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Ending SnakeStep" << std::endl;
        
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Starting MarkKnotSpansAvailable" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Starting MarkKnotSpansAvailable" << std::endl;
    }
    // print r_skin_sub_model_part
    std::ofstream outputFile("skin_condition_nodes_coordinates.txt");
    outputFile << std::setprecision(16);
    for (const auto& cond : r_skin_sub_model_part.Conditions()) {
        const auto& geometry = cond.GetGeometry();
        outputFile << "Condition ID: " << cond.Id() << std::endl;
        for (std::size_t i = 0; i < geometry.PointsNumber(); ++i) {
            const auto& node = geometry[i];
            outputFile << "Node ID: " << node.Id() << "  "
                    << node.X() << " " << node.Y() << " " << node.Z() << std::endl;
        }
        outputFile << std::endl;
    }
    outputFile.close();

    // std::ofstream knotFile("knot_spans_available_before.txt", std::ios::app);
    // knotFile << std::setprecision(16);
    // // Dump in a readable format: loop through all 3D entries
    // for (std::size_t loop = 0; loop < knot_spans_available.size(); ++loop) {
    //     knotFile << "Loop ID: " << loop << "\n";
    //     for (std::size_t k = 0; k < knot_spans_available[loop].size(); ++k) {
    //         for (std::size_t j = 0; j < knot_spans_available[loop][k].size(); ++j) {
    //             for (std::size_t i = 0; i < knot_spans_available[loop][k][j].size(); ++i) {
    //                 knotFile << knot_spans_available[loop][k][j][i] << " ";
    //             }
    //             knotFile << "\n";
    //         }
    //         knotFile << "\n";
    //     }
    //     knotFile << "\n";
    // }
    // knotFile.close();

    for (IndexType i = 0; i < NumberOfLoops; i++) {
        IndexType id_inner_loop = i;
        // Mark the knot_spans_available's for inner and outer loops
        MarkKnotSpansAvailable3D(id_inner_loop, points_bin, r_skin_sub_model_part, Lambda,
                                n_knot_spans_uvw, knot_step_uvw, starting_pos_uvw, knot_spans_available);  
        
        // std::ofstream knotFile("knot_spans_available_after.txt", std::ios::app);
        // knotFile << std::setprecision(16);
        // // Dump in a readable format: loop through all 3D entries
        // for (std::size_t loop = 0; loop < knot_spans_available.size(); ++loop) {
        //     knotFile << "Loop ID: " << loop << "\n";
        //     for (std::size_t k = 0; k < knot_spans_available[loop].size(); ++k) {
        //         for (std::size_t j = 0; j < knot_spans_available[loop][k].size(); ++j) {
        //             for (std::size_t i = 0; i < knot_spans_available[loop][k][j].size(); ++i) {
        //                 knotFile << knot_spans_available[loop][k][j][i] << " ";
        //             }
        //             knotFile << "\n";
        //         }
        //         knotFile << "\n";
        //     }
        //     knotFile << "\n";
        // }
        // knotFile.close();

        // Make a copy of the original knot spans to safely modify
        auto cleaned_knot_spans = knot_spans_available;
        for (std::size_t loop = 0; loop < knot_spans_available.size(); ++loop) {
            auto& grid = knot_spans_available[loop]; // 3D grid for the current loop
            for (std::size_t k = 0; k < grid.size(); ++k) {
                for (std::size_t j = 0; j < grid[k].size(); ++j) {
                    for (std::size_t i = 0; i < grid[k][j].size(); ++i) {
                        // Process only cells with value 1
                        if (grid[k][j][i] == 1) {
                            bool has_neighbor = false;
                            // Check the 6 direct neighbors: ±x, ±y, ±z
                            const std::vector<std::tuple<int, int, int>> directions = {
                                {1, 0, 0}, {-1, 0, 0},   // x+ and x-
                                {0, 1, 0}, {0, -1, 0},   // y+ and y-
                                {0, 0, 1}, {0, 0, -1}    // z+ and z-
                            };
                            // Loop through all 6 directions to find at least one neighbor with value 1
                            for (const auto& [dx, dy, dz] : directions) {
                                int ni = static_cast<int>(i) + dx;
                                int nj = static_cast<int>(j) + dy;
                                int nk = static_cast<int>(k) + dz;

                                // Ensure neighbor indices are within bounds
                                if (nk >= 0 && nk < static_cast<int>(grid.size()) &&
                                    nj >= 0 && nj < static_cast<int>(grid[nk].size()) &&
                                    ni >= 0 && ni < static_cast<int>(grid[nk][nj].size())) {
                                    if (grid[nk][nj][ni] == 1) {
                                        has_neighbor = true;
                                        break; // No need to check further if a neighbor is found
                                    }
                                }
                            }
                            // If no neighbor is found, set this cell to 0 in the cleaned grid
                            if (!has_neighbor) {
                                cleaned_knot_spans[loop][k][j][i] = 0;
                            }
                        }
                    }
                }
            }
        }
        // Update the original grid with the cleaned one
        knot_spans_available = cleaned_knot_spans;



    
        if (EchoLevel >  0) {
            KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Ending MarkKnotSpansAvailable" << std::endl;
            KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Ending MarkKnotSpansAvailable" << std::endl;
        }
        
        if (is_inner) {
            CreateSurrogateBuondaryFromSnakeInner3D(id_inner_loop, r_skin_sub_model_part, points_bin, n_knot_spans_uvw, 
                                                    knot_vector_u, knot_vector_v, knot_vector_w, 
                                                    knot_spans_available, r_surrogate_sub_model_part);
            
            std::ofstream outputFile("surrogate_condition_nodes_coordinates.txt", std::ios::app);
            outputFile << std::setprecision(16);
            for (auto& cond : r_surrogate_sub_model_part.Conditions()) {
                auto geometry = cond.GetGeometry();
                outputFile << "Condition ID: " << cond.Id() << std::endl;
                for (std::size_t i = 0; i < geometry.PointsNumber(); ++i) {
                    const auto& node = geometry[i];
                    outputFile << "Node ID: " << node.Id() << "  "
                                << node.X() << " " << node.Y() << " " << node.Z() << std::endl;
                }
                outputFile << std::endl;
            }
            outputFile.close();
                                                                
            if (EchoLevel >  0)
                KRATOS_INFO("::[SnakeSbmProcess]::") << "Inner :: Snake process has finished" << std::endl;
        }
        else {
            CreateSurrogateBuondaryFromSnakeOuter3D(id_inner_loop, r_skin_sub_model_part, points_bin, n_knot_spans_uvw, 
                                                    knot_vector_u, knot_vector_v, knot_vector_w, starting_pos_uvw, 
                                                    knot_spans_available, r_surrogate_sub_model_part);
                                                    
            
            std::ofstream outputFile("surrogate_condition_nodes_coordinates.txt", std::ios::app);
            outputFile << std::setprecision(16);
            for (auto& cond : r_surrogate_sub_model_part.Conditions()) {
                auto geometry = cond.GetGeometry();
                outputFile << "Condition ID: " << cond.Id() << std::endl;
                for (std::size_t i = 0; i < geometry.PointsNumber(); ++i) {
                    const auto& node = geometry[i];
                    outputFile << "Node ID: " << node.Id() << "  "
                                << node.X() << " " << node.Y() << " " << node.Z() << std::endl;
                }
                outputFile << std::endl;
            }
            outputFile.close();                                       
            if (EchoLevel >  0)
                KRATOS_INFO("::[SnakeSbmProcess]::") << "Outer :: Snake process has finished" << std::endl;
        }
    }

    if (EchoLevel >  0) {
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", is_inner) << "Inner :: Loop finished" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmProcess]::", !is_inner) << "Outer :: Loop finished" << std::endl;
    }
}


bool SnakeSbmProcess::IsInside3D(
    const std::vector<std::vector<int>>& rKnotSpanUVW,
    const std::vector<int>& NumberKnotSpansUVW)
{
    for (int point = 0; point < 3; ++point) {
        if (rKnotSpanUVW[0][point] < 0 || rKnotSpanUVW[0][point] >= NumberKnotSpansUVW[0] ||  // u
            rKnotSpanUVW[1][point] < 0 || rKnotSpanUVW[1][point] >= NumberKnotSpansUVW[1] ||  // v
            rKnotSpanUVW[2][point] < 0 || rKnotSpanUVW[2][point] >= NumberKnotSpansUVW[2])    // w
        {
            return true;  // At least one node is outside
        }
    }

    return false;  // All nodes are inside
}


void SnakeSbmProcess::SnakeStep3D(
    const int IdMatrix, 
    const std::vector<std::vector<int>>& rKnotSpansUVW, 
    const std::vector<std::vector<double>>& rConditionCoord, 
    const Vector rKnotStepUVW, 
    const Vector rStartingPosition,
    ModelPart& rSkinModelPart, 
    std::vector<std::vector<std::vector<std::vector<int>>>>& rKnotSpansAvailable,
    array_1d<IndexType, 3>& ordered_ids)
{
    bool isSplitted = false;
    if (rKnotSpansUVW[0][0] != rKnotSpansUVW[0][1] || rKnotSpansUVW[0][0] != rKnotSpansUVW[0][2] || rKnotSpansUVW[0][1] != rKnotSpansUVW[0][2] || 
        rKnotSpansUVW[1][0] != rKnotSpansUVW[1][1] || rKnotSpansUVW[1][0] != rKnotSpansUVW[1][2] || rKnotSpansUVW[1][1] != rKnotSpansUVW[1][2] || 
        rKnotSpansUVW[2][0] != rKnotSpansUVW[2][1] || rKnotSpansUVW[2][0] != rKnotSpansUVW[2][2] || rKnotSpansUVW[2][1] != rKnotSpansUVW[2][2]) 
    {
        auto squared_distance = [&](int i, int j) {
            double dx = rConditionCoord[0][i] - rConditionCoord[0][j];
            double dy = rConditionCoord[1][i] - rConditionCoord[1][j];
            double dz = rConditionCoord[2][i] - rConditionCoord[2][j];
            return dx*dx + dy*dy + dz*dz;
        };
        if (std::abs(rKnotSpansUVW[0][0] - rKnotSpansUVW[0][1]) > 1 || std::abs(rKnotSpansUVW[0][0] - rKnotSpansUVW[0][2]) > 1 || std::abs(rKnotSpansUVW[0][1] - rKnotSpansUVW[0][2]) > 1 ||
            std::abs(rKnotSpansUVW[1][0] - rKnotSpansUVW[1][1]) > 1 || std::abs(rKnotSpansUVW[1][0] - rKnotSpansUVW[1][2]) > 1 || std::abs(rKnotSpansUVW[1][1] - rKnotSpansUVW[1][2]) > 1 ||
            std::abs(rKnotSpansUVW[2][0] - rKnotSpansUVW[2][1]) > 1 || std::abs(rKnotSpansUVW[2][0] - rKnotSpansUVW[2][2]) > 1 || std::abs(rKnotSpansUVW[2][1] - rKnotSpansUVW[2][2]) > 1 ||
            // std::abs(rConditionCoord[0][0] - rConditionCoord[0][1]) > rKnotStepUVW[0]/50 ||
            // std::abs(rConditionCoord[0][1] - rConditionCoord[0][2]) > rKnotStepUVW[0]/50 ||
            // std::abs(rConditionCoord[0][2] - rConditionCoord[0][0]) > rKnotStepUVW[0]/50 ||
            // std::abs(rConditionCoord[1][0] - rConditionCoord[1][1]) > rKnotStepUVW[1]/50 ||
            // std::abs(rConditionCoord[1][1] - rConditionCoord[1][2]) > rKnotStepUVW[1]/50 ||
            // std::abs(rConditionCoord[1][2] - rConditionCoord[1][0]) > rKnotStepUVW[1]/50 ||
            // std::abs(rConditionCoord[2][0] - rConditionCoord[2][1]) > rKnotStepUVW[2]/50 ||
            // std::abs(rConditionCoord[2][1] - rConditionCoord[2][2]) > rKnotStepUVW[2]/50 ||
            // std::abs(rConditionCoord[2][2] - rConditionCoord[2][0]) > rKnotStepUVW[2]/50
            squared_distance(0, 1) > rKnotStepUVW[0]*rKnotStepUVW[0]/5 ||
            squared_distance(0, 2) > rKnotStepUVW[0]*rKnotStepUVW[0]/5 ||
            squared_distance(1, 2) > rKnotStepUVW[0]*rKnotStepUVW[0]/5
        )
        {
            isSplitted = true;
            // KRATOS_INFO("::[SnakeSBMUtilities]::") << "SnakeStep :: Splitting a 3D condition" << std::endl;
            // KRATOS_WATCH(rConditionCoord)

            // Midpoints
            double x12 = (rConditionCoord[0][0] + rConditionCoord[0][1]) / 2.0;
            double y12 = (rConditionCoord[1][0] + rConditionCoord[1][1]) / 2.0;
            double z12 = (rConditionCoord[2][0] + rConditionCoord[2][1]) / 2.0;

            double x23 = (rConditionCoord[0][1] + rConditionCoord[0][2]) / 2.0;
            double y23 = (rConditionCoord[1][1] + rConditionCoord[1][2]) / 2.0;
            double z23 = (rConditionCoord[2][1] + rConditionCoord[2][2]) / 2.0;

            double x31 = (rConditionCoord[0][2] + rConditionCoord[0][0]) / 2.0;
            double y31 = (rConditionCoord[1][2] + rConditionCoord[1][0]) / 2.0;
            double z31 = (rConditionCoord[2][2] + rConditionCoord[2][0]) / 2.0;

            // Compute in which knot spans they lie
            int ku12 = (x12 - rStartingPosition[0]) / rKnotStepUVW[0];
            int kv12 = (y12 - rStartingPosition[1]) / rKnotStepUVW[1];
            int kw12 = (z12 - rStartingPosition[2]) / rKnotStepUVW[2];

            int ku23 = (x23 - rStartingPosition[0]) / rKnotStepUVW[0];
            int kv23 = (y23 - rStartingPosition[1]) / rKnotStepUVW[1];
            int kw23 = (z23 - rStartingPosition[2]) / rKnotStepUVW[2];

            int ku31 = (x31 - rStartingPosition[0]) / rKnotStepUVW[0];
            int kv31 = (y31 - rStartingPosition[1]) / rKnotStepUVW[1];
            int kw31 = (z31 - rStartingPosition[2]) / rKnotStepUVW[2];

            // CREATE three NEW NODES in the middle of each side
            auto idNode_12 = rSkinModelPart.NumberOfNodes() + 1;
            auto idNode_23 = idNode_12 + 1;
            auto idNode_31 = idNode_23 + 1;
            rSkinModelPart.CreateNewNode(idNode_12, x12, y12, z12);
            rSkinModelPart.CreateNewNode(idNode_23, x23, y23, z23);
            rSkinModelPart.CreateNewNode(idNode_31, x31, y31, z31);

            // We keep the same order so that the direction of the normal remains the same (inside/outiside keeps working)
            array_1d<IndexType, 3> ordered_ids_1;
            ordered_ids_1[0] = ordered_ids[0];
            ordered_ids_1[1] = idNode_12;
            ordered_ids_1[2] = idNode_31;
            array_1d<IndexType, 3> ordered_ids_2;
            ordered_ids_2[0] = ordered_ids[1];
            ordered_ids_2[1] = idNode_23;
            ordered_ids_2[2] = idNode_12;
            array_1d<IndexType, 3> ordered_ids_3;
            ordered_ids_3[0] = ordered_ids[2];
            ordered_ids_3[1] = idNode_31;
            ordered_ids_3[2] = idNode_23;
            array_1d<IndexType, 3> ordered_ids_4;
            ordered_ids_4[0] = idNode_12;
            ordered_ids_4[1] = idNode_23;
            ordered_ids_4[2] = idNode_31;

            // KRATOS_WATCH("TRIANGLES DIVISION")

            // const double tol = 1e-7;
            // if (std::abs(rConditionCoord[0][0] - rConditionCoord[0][1]) < tol &&
            //     std::abs(rConditionCoord[0][1] - rConditionCoord[0][2]) < tol &&
            //     std::abs(rConditionCoord[1][0] - rConditionCoord[1][1]) < tol &&
            //     std::abs(rConditionCoord[1][1] - rConditionCoord[1][2]) < tol &&
            //     std::abs(rConditionCoord[2][0] - rConditionCoord[2][1]) < tol &&
            //     std::abs(rConditionCoord[2][1] - rConditionCoord[2][2]) < tol)
            // {
            //     KRATOS_WARNING("SnakeStep3D") << "Degenerate triangle detected (within tolerance), skipping..." << std::endl;
            //     exit(0);
            // }

            // We do it recursively for all the new 4 tringular conditions
            // First subtriangle 
            // KRATOS_WATCH("first_subtriangle")
            SnakeStep3D(
                IdMatrix,
                {{rKnotSpansUVW[0][0], ku12, ku31},{rKnotSpansUVW[1][0], kv12, kv31},{rKnotSpansUVW[2][0], kw12, kw31}},
                {{rConditionCoord[0][0], x12, x31},{rConditionCoord[1][0], y12, y31},{rConditionCoord[2][0], z12, z31}},
                rKnotStepUVW,
                rStartingPosition,
                rSkinModelPart,
                rKnotSpansAvailable,
                ordered_ids_1
            );

            // Second subtriangle
            // KRATOS_WATCH("second_subtriangle")
            SnakeStep3D(
                IdMatrix,
                {{rKnotSpansUVW[0][1],  ku23, ku12},{rKnotSpansUVW[1][1], kv23, kv12},{rKnotSpansUVW[2][1], kw23, kw12}},
                {{rConditionCoord[0][1], x23, x12},{rConditionCoord[1][1], y23, y12},{rConditionCoord[2][1], z23, z12}},
                rKnotStepUVW,
                rStartingPosition,
                rSkinModelPart,
                rKnotSpansAvailable,
                ordered_ids_2
            );

            // Third subtriangle
            // KRATOS_WATCH("third_subtriangle")
            SnakeStep3D(
                IdMatrix,
                {{rKnotSpansUVW[0][2], ku31, ku23},{rKnotSpansUVW[1][2], kv31, kv23},{rKnotSpansUVW[2][2], kw31, kw23}},
                {{rConditionCoord[0][2], x31, x23},{rConditionCoord[1][2], y31, y23},{rConditionCoord[2][2], z31, z23}},
                rKnotStepUVW,
                rStartingPosition,
                rSkinModelPart,
                rKnotSpansAvailable,
                ordered_ids_3
            );

            // Central subtriangle --> Is it necessary?
            // KRATOS_WATCH("fourth_subtriangle")
            SnakeStep3D(
                IdMatrix,
                {{ku12, ku23, ku31},{kv12, kv23, kv31},{kw12, kw23, kw31}},
                {{x12,  x23,  x31},{y12,  y23,  y31},{z12, z23, z31}},
                rKnotStepUVW,
                rStartingPosition,
                rSkinModelPart,
                rKnotSpansAvailable,
                ordered_ids_4
            );

        }
        /* The condition it is NOT split
            Check if the true boundary crosses an u, v, or w knot value */
        // Points between node 1 & 2 
        // The condition is NOT split: Check if the true boundary crosses a u, v, or w knot value

        // Points between node 1 & 2
        else if (rKnotSpansUVW[0][0] != rKnotSpansUVW[0][1]) { // u knot crossed between node 1 & 2
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][1]] = 2;
        }
        else if (rKnotSpansUVW[1][0] != rKnotSpansUVW[1][1]) { // v knot crossed
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][1]][rKnotSpansUVW[0][0]] = 2;
        }
        else if (rKnotSpansUVW[2][0] != rKnotSpansUVW[2][1]) { // w knot crossed
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][1]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
        }
        // Points between node 1 & 3
        else if (rKnotSpansUVW[0][0] != rKnotSpansUVW[0][2]) {
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][2]] = 2;
        }
        else if (rKnotSpansUVW[1][0] != rKnotSpansUVW[1][2]) {
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][2]][rKnotSpansUVW[0][0]] = 2;
        }
        else if (rKnotSpansUVW[2][0] != rKnotSpansUVW[2][2]) {
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][0]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][2]][rKnotSpansUVW[1][0]][rKnotSpansUVW[0][0]] = 2;
        }
        // Points between node 3 & 2
        else if (rKnotSpansUVW[0][2] != rKnotSpansUVW[0][1]) {
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][2]][rKnotSpansUVW[1][2]][rKnotSpansUVW[0][2]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][2]][rKnotSpansUVW[1][2]][rKnotSpansUVW[0][1]] = 2;
        }
        else if (rKnotSpansUVW[1][2] != rKnotSpansUVW[1][1]) {
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][2]][rKnotSpansUVW[1][2]][rKnotSpansUVW[0][2]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][2]][rKnotSpansUVW[1][1]][rKnotSpansUVW[0][2]] = 2;
        }
        else if (rKnotSpansUVW[2][2] != rKnotSpansUVW[2][1]) {
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][2]][rKnotSpansUVW[1][2]][rKnotSpansUVW[0][2]] = 2;
            rKnotSpansAvailable[IdMatrix][rKnotSpansUVW[2][1]][rKnotSpansUVW[1][2]][rKnotSpansUVW[0][2]] = 2;
        }
        else {
            // TODO: delete this
            KRATOS_WATCH("Something went wrong");
            exit(0);
        }

    }
    if (!isSplitted) {
        // TO DO -> do you want to create a new skin_model_part with the split triangles?
        //           How do you manage the creation of the nodes? They might already exist in the model part
        //           Maybe is better to add nodes and condition to the initial_skin_model_part 

        //// Create 1 new conditions for each skin condition 
        Properties::Pointer p_cond_prop = rSkinModelPart.pGetProperties(0);

        // KRATOS_WATCH("SnakeStep3D -> create new condition")
        // KRATOS_WATCH("\n")

        IndexType last_condition_id;
        if (rSkinModelPart.NumberOfConditions() == 0) {
            last_condition_id = 0;
        } else {
            last_condition_id = (rSkinModelPart.GetRootModelPart().ConditionsEnd()-1)->Id();
        }

        Condition::Pointer p_cond1 = rSkinModelPart.CreateNewCondition("SurfaceCondition3D3N", last_condition_id+1, 
                                                                        {{ordered_ids[0], ordered_ids[1], ordered_ids[2]}}, p_cond_prop );
        
        // KRATOS_WATCH("SnakeStep3D -> create END condition")
    }
}


/**
* Marking process:
*   1) We set to 2 all the cut knot spans
*   2) We check the 26 neighbor knot spans and set them to 1|0  if inside|outside
*   3) We check all the cut knot spans and set them     to 1|-1 if inside|outside 
*/
void SnakeSbmProcess::MarkKnotSpansAvailable3D(
    const int IdMatrix,
    DynamicBins& rPointsBin, 
    const ModelPart& rSkinModelPart,
    const double Lambda, 
    const std::vector<int>& rNumberKnotSpans, 
    const array_1d<double, 3>& rKnotStepUVW,
    const Vector& rStartingPosition,
    std::vector<std::vector<std::vector<std::vector<int>>>> & rKnotSpansAvailable) 
{
    for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
        for (int i = 0; i < rNumberKnotSpans[1]; ++i) {
            for (int j = 0; j < rNumberKnotSpans[0]; ++j) {
                if (rKnotSpansAvailable[IdMatrix][k][i][j] == 2) {
    
                    // Loop through 26 neighbors (3x3x3 - 1)
                    for (int dk = -1; dk <= 1; ++dk) {
                        for (int di = -1; di <= 1; ++di) {
                            for (int dj = -1; dj <= 1; ++dj) {
                                if (dk == 0 && di == 0 && dj == 0)
                                    continue;
    
                                int kk = k + dk;
                                int ii = i + di;
                                int jj = j + dj;
                                
                                if (kk >= 0 && kk < rNumberKnotSpans[2] &&
                                    ii >= 0 && ii < rNumberKnotSpans[1] &&
                                    jj >= 0 && jj < rNumberKnotSpans[0] &&
                                    rKnotSpansAvailable[IdMatrix][kk][ii][jj] == 0) {
    
                                    double x_gp = (jj + 0.5) * rKnotStepUVW[0] + rStartingPosition[0];
                                    double y_gp = (ii + 0.5) * rKnotStepUVW[1] + rStartingPosition[1];
                                    double z_gp = (kk + 0.5) * rKnotStepUVW[2] + rStartingPosition[2];
                                    Point gauss_point(x_gp, y_gp, z_gp);
    
                                    if (IsPointInsideSkinBoundary3D(gauss_point, rPointsBin, rSkinModelPart)) {
                                        rKnotSpansAvailable[IdMatrix][kk][ii][jj] = 3; // originally was == 1
                                    }
                                }
                            }
                        }
                    }


                }
            }
        }
    }


    for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
        for (int i = 0; i < rNumberKnotSpans[1]; ++i) {
            for (int j = 0; j < rNumberKnotSpans[0]; ++j) {
                if (rKnotSpansAvailable[IdMatrix][k][i][j] > 1) { // originally was == 2
                    // Create fake Gauss points to check if the majority of them are inside
                    const int num_fake_gauss_points = 4;
                    int number_of_inside_gaussian_points = 0;

                    for (IndexType gp_x = 0; gp_x < num_fake_gauss_points; ++gp_x) {
                        double x_coord = j * rKnotStepUVW[0] + rKnotStepUVW[0] / (num_fake_gauss_points + 1) * (gp_x + 1) + rStartingPosition[0];

                        for (IndexType gp_y = 0; gp_y < num_fake_gauss_points; ++gp_y) {
                            double y_coord = i * rKnotStepUVW[1] + rKnotStepUVW[1] / (num_fake_gauss_points + 1) * (gp_y + 1) + rStartingPosition[1];

                            for (IndexType gp_z = 0; gp_z < num_fake_gauss_points; ++gp_z) {
                                double z_coord = k * rKnotStepUVW[2] + rKnotStepUVW[2] / (num_fake_gauss_points + 1) * (gp_z + 1) + rStartingPosition[2];

                                Point gauss_point(x_coord, y_coord, z_coord);

                                if (IsPointInsideSkinBoundary3D(gauss_point, rPointsBin, rSkinModelPart)) {
                                    number_of_inside_gaussian_points++;
                                }
                            }
                        }
                    }
    
                    // Classify as active or inactive
                    const double threshold = Lambda * std::pow(num_fake_gauss_points, 3);

                    if (number_of_inside_gaussian_points < threshold) {
                        rKnotSpansAvailable[IdMatrix][k][i][j] = -1;
                    } else {
                        rKnotSpansAvailable[IdMatrix][k][i][j] = 1;
                    }
                }
            }
        }
    }
}


bool SnakeSbmProcess::IsPointInsideSkinBoundary3D(
    const Point& rPoint1, 
    DynamicBins& rPointsBin, 
    const ModelPart& rSkinModelPart)
{    
    // Get the nearest point of the true boundary
    DynamicBinsPointerType p_point_to_search = DynamicBinsPointerType(new PointType(1, rPoint1.X(), rPoint1.Y(), rPoint1.Z()));
    DynamicBinsPointerType p_nearest_point = rPointsBin.SearchNearestPoint(*p_point_to_search);
    
    // Get the closest Condition
    IndexType id_1 = p_nearest_point->Id();
    auto nearest_condition = rSkinModelPart.GetCondition(id_1);
    // KRATOS_WATCH(nearest_condition)

    // Point0 -> pointToSearch
    Point point_1 = nearest_condition.GetGeometry()[0]; // FIRST POINT IN TRUE GEOM
    Point point_2 = nearest_condition.GetGeometry()[1]; // SECOND POINT IN TRUE GEOM
    Point point_3 = nearest_condition.GetGeometry()[2]; // THIRD POINT IN TRUE GEOM

    array_1d<double,3> v_1 = point_2 - point_1;
    array_1d<double,3> v_2 = point_3 - point_1;
    array_1d<double,3> normal;
    MathUtils<double>::CrossProduct(normal, v_1, v_2);
    
    normal = normal/norm_2(normal);
    
    array_1d<double,3> center_to_point_to_search = rPoint1 - nearest_condition.GetGeometry().Center() ;

    center_to_point_to_search = center_to_point_to_search/norm_2(center_to_point_to_search);

    // Scalar product between the nornal and the center_to_point_to_search
    bool is_inside = false;
    if (inner_prod(normal, center_to_point_to_search) > 0) {is_inside = true;}

    return is_inside; // should be "return is_inside;""

}


// bool SnakeSbmProcess::IsPointInsideSkinBoundary3D(
//     const Point& rPoint1,
//     DynamicBins& rPointsBin,
//     const ModelPart& rSkinModelPart)
// {
//     constexpr double search_radius = 1e-2; // puoi adattarlo alla scala del problema
//     constexpr std::size_t max_results = 10;

//     std::vector<DynamicBinsPointerType> results(max_results);
//     std::vector<double> list_of_distances(max_results);

//     using PointType = Point;
//     using PointerType = typename PointType::Pointer;
//     using SizeType = std::size_t;
//     using IndexType = std::size_t;

//     auto p_point_to_search = PointerType(new PointType(1, rPoint1.X(), rPoint1.Y(), rPoint1.Z()));
//     SizeType obtained_results = rPointsBin.SearchInRadius(*p_point_to_search, search_radius, results.begin(), list_of_distances.begin(), max_results);

//     KRATOS_ERROR_IF(obtained_results == 0) 
//         << "::[SnakeSbmProcess]:: No nearby points found in radius search for point: " << *p_point_to_search << std::endl;

//     double min_distance = 1e12;
//     bool found_coherent_result = false;
//     bool is_inside = false;

//     for (IndexType i = 0; i < obtained_results; ++i) {
//         auto& p_candidate = results[i];
//         if (!p_candidate) continue;

//         IndexType cond_id = p_candidate->Id();
//         auto it_cond = rSkinModelPart.Conditions().find(cond_id);
//         if (it_cond == rSkinModelPart.Conditions().end())
//             continue;

//         const auto& cond = *it_cond;
//         const auto& geom = cond.GetGeometry();
//         if (geom.size() != 3) continue;

//         // Compute normal of triangle
//         const Point& p1 = geom[0];
//         const Point& p2 = geom[1];
//         const Point& p3 = geom[2];

//         array_1d<double, 3> v1 = p2 - p1;
//         array_1d<double, 3> v2 = p3 - p1;
//         array_1d<double, 3> normal;
//         MathUtils<double>::CrossProduct(normal, v1, v2);
//         const double norm = norm_2(normal);
//         if (norm < 1e-12) continue;
//         normal /= norm;

//         array_1d<double, 3> center_to_point = rPoint1 - geom.Center();
//         double signed_distance = inner_prod(normal, center_to_point);

//         if (std::abs(signed_distance) < 1e-12) continue; // tangente, incerto

//         // Usa solo i triangoli che danno un risultato coerente (positivo o negativo)
//         if (list_of_distances[i] < min_distance) {
//             min_distance = list_of_distances[i];
//             is_inside = (signed_distance > 0.0); // se >0 → dentro
//             found_coherent_result = true;
//         }
//     }

//     KRATOS_ERROR_IF_NOT(found_coherent_result) 
//         << "::[SnakeSbmProcess]:: No valid orientation found for any nearby condition for point: " << *p_point_to_search << std::endl;

//     return is_inside;
// }


/**
* summary of knot_spans_available:
    " 1"  -> interior knot spans                                  
    "-1"  -> exterior knot spans well checked
    " 0"  -> exterior knot spans OR very interior knot spans (more 
                than one ks away from surrogate boundary)
*/
void SnakeSbmProcess::CreateSurrogateBuondaryFromSnakeInner3D(
    const int IdMatrix, 
    const ModelPart& rSkinModelPartInner, 
    DynamicBins& rPointsBinInner,
    const std::vector<int>& rNumberKnotSpans, 
    const Vector& rKnotVectorU, 
    const Vector& rKnotVectorV,
    const Vector& rKnotVectorW,
    std::vector<std::vector<std::vector<std::vector<int>>>>& rKnotSpansAvailable,
    ModelPart& rSurrogateModelPartInner
    ) 
{
    const double knot_step_u = rKnotVectorU[1] - rKnotVectorU[0];
    const double knot_step_v = rKnotVectorV[1] - rKnotVectorV[0];
    const double knot_step_w = rKnotVectorW[1] - rKnotVectorW[0];

    IndexType id_surrogate_first_node = rSurrogateModelPartInner.GetParentModelPart().NumberOfNodes() + 1;
    IndexType id_node = id_surrogate_first_node;

    // Create nodes for each knot span position
    for (int k = 0; k < rNumberKnotSpans[2]+1; ++k) {
        for (int j = 0; j < rNumberKnotSpans[1]+1; ++j) {
            for (int i = 0; i < rNumberKnotSpans[0]+1; ++i) {
                rSurrogateModelPartInner.CreateNewNode(id_node, rKnotVectorU[i], rKnotVectorV[j], rKnotVectorW[k]);
                ++id_node;
            }
        }
    }

    auto p_cond_prop = rSurrogateModelPartInner.pGetProperties(0);
    IndexType id_surrogate_condition = rSurrogateModelPartInner.NumberOfConditions() + 1;
    IndexType id_surrogate_first_condition = id_surrogate_condition;

    // Sweep along u-direction
    for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
        for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
            bool check_next_point = false;
            for (int i = 0; i < rNumberKnotSpans[0]; ++i) {
                if (check_next_point) {
                    Point center_point((i + 0.5) * knot_step_u, (j + 0.5) * knot_step_v, (k + 0.5) * knot_step_w);
                    bool is_exiting = false;
                    if (rKnotSpansAvailable[IdMatrix][k][j][i] == 1) {
                        // Already marked active
                    } else if (IsPointInsideSkinBoundary3D(center_point, rPointsBinInner, rSkinModelPartInner)) {
                        // STILL INSIDE --> do not save nothing and update knot_spans_available 
                        if (rKnotSpansAvailable[IdMatrix][k][j][i] == -1)
                            is_exiting = true;
                        else
                            rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
                    } else {
                        is_exiting = true;
                    }

                    if (is_exiting) {
                        int node1_i = i; int node1_j = j;   int node1_k = k;
                        int node2_i = i; int node2_j = j+1; int node2_k = k;
                        int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                        int node4_i = i; int node4_j = j;   int node4_k = k+1;
                        /*  
                            Formula to connect i,j,k to the id of the model_part
                            i + j*(rNumberKnotSpans[0]+1) + k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1) + 1;
                        */
                        IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                        IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                        IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                        IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                        auto pcond = rSurrogateModelPartInner.CreateNewCondition("SurfaceCondition3D4N", id_surrogate_condition++, {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                        pcond->Set(BOUNDARY, false);
                        check_next_point = false;
                    }
                } else if (rKnotSpansAvailable[IdMatrix][k][j][i] == 1) {
                    int node1_i = i; int node1_j = j;   int node1_k = k;
                    int node2_i = i; int node2_j = j+1; int node2_k = k;
                    int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                    int node4_i = i; int node4_j = j;   int node4_k = k+1;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartInner.CreateNewCondition("SurfaceCondition3D4N", id_surrogate_condition++, {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, true);
                    check_next_point = true;
                }
            }
        }
    }

    // Sweep along v-direction (no boundary check needed)
    for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
        for (int i = 0; i < rNumberKnotSpans[0]; ++i) {

            bool was_active = false;

            for (int j = 0; j < rNumberKnotSpans[1]; ++j) {

                bool is_active = (rKnotSpansAvailable[IdMatrix][k][j][i] == 1);

                if (!was_active && is_active) {
                    // ENTERING an active region → generate INLET face
                    int node1_i = i;   int node1_j = j;   int node1_k = k;
                    int node2_i = i+1; int node2_j = j;   int node2_k = k;
                    int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                    int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartInner.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, true);
                }

                else if (was_active && !is_active) {
                    // EXITING an active region → generate OUTLET face
                    int node1_i = i;   int node1_j = j;   int node1_k = k;
                    int node2_i = i+1; int node2_j = j;   int node2_k = k;
                    int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                    int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartInner.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, false);
                }

                // Update state for next step
                was_active = is_active;
            }
        }
    }


    // Sweep along w-direction (no boundary check needed)
    for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
        for (int i = 0; i < rNumberKnotSpans[0]; ++i) {

            bool was_active = false;

            for (int k = 0; k < rNumberKnotSpans[2]; ++k) {

                bool is_active = (rKnotSpansAvailable[IdMatrix][k][j][i] == 1);

                if (!was_active && is_active) {
                    // ENTERING an active region → generate INLET face (normal along -w)
                    int node1_i = i;   int node1_j = j;     int node1_k = k;
                    int node2_i = i+1; int node2_j = j;     int node2_k = k;
                    int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                    int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartInner.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, true);
                }

                else if (was_active && !is_active) {
                    // EXITING an active region → generate OUTLET face (normal along +w)
                    int node1_i = i;   int node1_j = j;     int node1_k = k;
                    int node2_i = i+1; int node2_j = j;     int node2_k = k;
                    int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                    int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartInner.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, false);
                }

                was_active = is_active;
            }
        }
    }


    // Create dummy element to keep track of first and last surrogate condition
    IndexType id_surrogate_last_condition = id_surrogate_condition - 1;
    IndexType id_surrogate_elem = rSurrogateModelPartInner.NumberOfElements() + 1;
    std::vector<ModelPart::IndexType> elem_nodes{id_surrogate_first_condition, id_surrogate_last_condition};
    rSurrogateModelPartInner.CreateNewElement("Element3D2N", id_surrogate_elem, elem_nodes, p_cond_prop);
}


/**
* summary of knot_spans_available:
    " 1"  -> interior knot spans                                  
    "-1"  -> exterior knot spans well checked
    " 0"  -> exterior knot spans OR very interior knot spans (more 
                than one ks away from surrogate boundary)
*/
void SnakeSbmProcess::CreateSurrogateBuondaryFromSnakeOuter3D(
    const int IdMatrix, 
    const ModelPart& rSkinModelPartOuter, 
    DynamicBins& rPointsBinOuter,
    const std::vector<int>& rNumberKnotSpans, 
    const Vector& rKnotVectorU, 
    const Vector& rKnotVectorV,
    const Vector& rKnotVectorW,
    const Vector& rStartingPositionUVW,
    std::vector<std::vector<std::vector<std::vector<int>>>>& rKnotSpansAvailable,
    ModelPart& rSurrogateModelPartOuter
    ) 
{
    // CHECK ALL THE EXTERNAL KNOT SPANS – includes an additional layer
    const double knot_step_u = rKnotVectorU[1] - rKnotVectorU[0];
    const double knot_step_v = rKnotVectorV[1] - rKnotVectorV[0];
    const double knot_step_w = rKnotVectorW[1] - rKnotVectorW[0];

    // // LEFT BOUNDARY (u)
    // for (int i = 0; i < 2; ++i) {
    //     for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
    //         for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
    //             Point centroid((i + 0.5) * knot_step_u + rStartingPositionUVW[0],
    //                            (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
    //                            (k + 0.5) * knot_step_w + rStartingPositionUVW[2]);
    //             if (IsPointInsideSkinBoundary3D(centroid, rPointsBinOuter, rSkinModelPartOuter) &&
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] != -1) {
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
    //             }
    //         }
    //     }
    // }

    // // RIGHT BOUNDARY (u)
    // for (int i = rNumberKnotSpans[0] - 1; i >= rNumberKnotSpans[0] - 2; --i) {
    //     for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
    //         for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
    //             Point centroid((i + 0.5) * knot_step_u + rStartingPositionUVW[0],
    //                         (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
    //                         (k + 0.5) * knot_step_w + rStartingPositionUVW[2]);
    //             if (IsPointInsideSkinBoundary3D(centroid, rPointsBinOuter, rSkinModelPartOuter) &&
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] != -1) {
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
    //             }
    //         }
    //     }
    // }

    // // FRONT BOUNDARY (v)
    // for (int j = 0; j < 2; ++j) {
    //     for (int i = 0; i < rNumberKnotSpans[0]; ++i) {
    //         for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
    //             Point centroid((i + 0.5) * knot_step_u + rStartingPositionUVW[0],
    //                         (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
    //                         (k + 0.5) * knot_step_w + rStartingPositionUVW[2]);
    //             if (IsPointInsideSkinBoundary3D(centroid, rPointsBinOuter, rSkinModelPartOuter) &&
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] != -1) {
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
    //             }
    //         }
    //     }
    // }

    // // BACK BOUNDARY (v)
    // for (int j = rNumberKnotSpans[1] - 1; j >= rNumberKnotSpans[1] - 2; --j) {
    //     for (int i = 0; i < rNumberKnotSpans[0]; ++i) {
    //         for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
    //             Point centroid((i + 0.5) * knot_step_u + rStartingPositionUVW[0],
    //                         (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
    //                         (k + 0.5) * knot_step_w + rStartingPositionUVW[2]);
    //             if (IsPointInsideSkinBoundary3D(centroid, rPointsBinOuter, rSkinModelPartOuter) &&
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] != -1) {
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
    //             }
    //         }
    //     }
    // }

    // // BOTTOM BOUNDARY (w)
    // for (int k = 0; k < 2; ++k) {
    //     for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
    //         for (int i = 0; i < rNumberKnotSpans[0]; ++i) {
    //             Point centroid((i + 0.5) * knot_step_u + rStartingPositionUVW[0],
    //                         (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
    //                         (k + 0.5) * knot_step_w + rStartingPositionUVW[2]);
    //             if (IsPointInsideSkinBoundary3D(centroid, rPointsBinOuter, rSkinModelPartOuter) &&
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] != -1) {
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
    //             }
    //         }
    //     }
    // }

    // // TOP BOUNDARY (w)
    // for (int k = rNumberKnotSpans[2] - 1; k >= rNumberKnotSpans[2] - 2; --k) {
    //     for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
    //         for (int i = 0; i < rNumberKnotSpans[0]; ++i) {
    //             Point centroid((i + 0.5) * knot_step_u + rStartingPositionUVW[0],
    //                         (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
    //                         (k + 0.5) * knot_step_w + rStartingPositionUVW[2]);
    //             if (IsPointInsideSkinBoundary3D(centroid, rPointsBinOuter, rSkinModelPartOuter) &&
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] != -1) {
    //                 rKnotSpansAvailable[IdMatrix][k][j][i] = 1;
    //             }
    //         }
    //     }
    // }


    IndexType id_surrogate_first_node = rSurrogateModelPartOuter.GetParentModelPart().NumberOfNodes() + 1;
    IndexType id_node = id_surrogate_first_node;

    // Create nodes for each knot span position
    for (int k = 0; k < rNumberKnotSpans[2]+1; ++k) {
        for (int j = 0; j < rNumberKnotSpans[1]+1; ++j) {
            for (int i = 0; i < rNumberKnotSpans[0]+1; ++i) {
                rSurrogateModelPartOuter.CreateNewNode(id_node, rKnotVectorU[i], rKnotVectorV[j], rKnotVectorW[k]);
                ++id_node;
            }
        }
    }

    auto p_cond_prop = rSurrogateModelPartOuter.pGetProperties(0);
    IndexType id_surrogate_condition = rSurrogateModelPartOuter.GetRootModelPart().NumberOfConditions() + 1;


    // Sweep along u-direction
    for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
        for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
            // KRATOS_WATCH("\n")
            // KRATOS_WATCH(rKnotSpansAvailable[IdMatrix][k][j])
            bool check_next_point = false;
            for (int i = 0; i < rNumberKnotSpans[0]; ++i) {
                // KRATOS_WATCH(rKnotSpansAvailable[IdMatrix][k][j][i])
                Point center_point(
                    (i + 0.5) * knot_step_u + rStartingPositionUVW[0],
                    (j + 0.5) * knot_step_v + rStartingPositionUVW[1],
                    (k + 0.5) * knot_step_w + rStartingPositionUVW[2]
                );

                if (check_next_point) {
                    // KRATOS_WATCH("check_next_point true")
                    bool is_exiting = false;
                    int& current_value = rKnotSpansAvailable[IdMatrix][k][j][i];

                    if (current_value == 1) {
                        // KRATOS_WATCH("current_value == 1")
                        // Already active, OK
                    } else {
                        if (current_value == -1) {
                            // KRATOS_WATCH("current_value == -1")
                            // Already -1, exit
                            is_exiting = true;
                        } else {
                            // current_value == 0
                            // KRATOS_WATCH("current_value == 0")
                            bool is_inside = IsPointInsideSkinBoundary3D(center_point, rPointsBinOuter, rSkinModelPartOuter);
                            // KRATOS_WATCH(is_inside)
                            if (is_inside) { // ??????
                                rKnotSpansAvailable[IdMatrix][k][j][i] = 1; // <-- IMPORTANT!! Inside I want to have all 1's
                                // KRATOS_WATCH(center_point)
                                // exit(0);
                            } else {
                                // exit(0);
                                is_exiting = true;
                            }
                        }
                    }

                    if (is_exiting) {
                        // KRATOS_WATCH(" save exit")
                        int node1_i = i; int node1_j = j;   int node1_k = k;
                        int node2_i = i; int node2_j = j+1; int node2_k = k;
                        int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                        int node4_i = i; int node4_j = j;   int node4_k = k+1;

                        IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                        IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                        IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                        IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                        auto pcond = rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", id_surrogate_condition++, {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                        pcond->Set(BOUNDARY, false);
                        check_next_point = false; // Exit surface
                    }
                } else if (rKnotSpansAvailable[IdMatrix][k][j][i] == 1) {
                    // KRATOS_WATCH(" save entry")
                    int node1_i = i; int node1_j = j;   int node1_k = k;
                    int node2_i = i; int node2_j = j+1; int node2_k = k;
                    int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                    int node4_i = i; int node4_j = j;   int node4_k = k+1;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", id_surrogate_condition++, {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, true);
                    check_next_point = true; // Enter surface
                }
            }
            // KRATOS_WATCH(rKnotSpansAvailable[IdMatrix][k][j])
        }
    }

    // Sweep along v-direction (no boundary check needed)
    for (int k = 0; k < rNumberKnotSpans[2]; ++k) {
        for (int i = 0; i < rNumberKnotSpans[0]; ++i) {

            bool was_active = false;

            for (int j = 0; j < rNumberKnotSpans[1]; ++j) {

                bool is_active = (rKnotSpansAvailable[IdMatrix][k][j][i] == 1);

                if (!was_active && is_active) {
                    // ENTERING an active region → generate INLET face
                    int node1_i = i;   int node1_j = j;   int node1_k = k;
                    int node2_i = i+1; int node2_j = j;   int node2_k = k;
                    int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                    int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, true);
                }

                else if (was_active && !is_active) {
                    // EXITING an active region → generate OUTLET face
                    int node1_i = i;   int node1_j = j;   int node1_k = k;
                    int node2_i = i+1; int node2_j = j;   int node2_k = k;
                    int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                    int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, false);
                }

                // Update state for next step
                was_active = is_active;
            }
        }
    }

    
    // Sweep along w-direction (no boundary check needed)
    for (int j = 0; j < rNumberKnotSpans[1]; ++j) {
        for (int i = 0; i < rNumberKnotSpans[0]; ++i) {

            bool was_active = false;

            for (int k = 0; k < rNumberKnotSpans[2]; ++k) {

                bool is_active = (rKnotSpansAvailable[IdMatrix][k][j][i] == 1);

                if (!was_active && is_active) {
                    // ENTERING an active region → generate INLET face (normal along -w)
                    int node1_i = i;   int node1_j = j;     int node1_k = k;
                    int node2_i = i+1; int node2_j = j;     int node2_k = k;
                    int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                    int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, true);
                }

                else if (was_active && !is_active) {
                    // EXITING an active region → generate OUTLET face (normal along +w)
                    int node1_i = i;   int node1_j = j;     int node1_k = k;
                    int node2_i = i+1; int node2_j = j;     int node2_k = k;
                    int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                    int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                    IndexType id_node_1 = id_surrogate_first_node + node1_i + node1_j*(rNumberKnotSpans[0]+1) + node1_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_2 = id_surrogate_first_node + node2_i + node2_j*(rNumberKnotSpans[0]+1) + node2_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_3 = id_surrogate_first_node + node3_i + node3_j*(rNumberKnotSpans[0]+1) + node3_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);
                    IndexType id_node_4 = id_surrogate_first_node + node4_i + node4_j*(rNumberKnotSpans[0]+1) + node4_k*(rNumberKnotSpans[1]+1)*(rNumberKnotSpans[0]+1);

                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition(
                        "SurfaceCondition3D4N", id_surrogate_condition++, 
                        {{id_node_1, id_node_2, id_node_3, id_node_4}}, p_cond_prop);
                    pcond->Set(BOUNDARY, false);
                }

                was_active = is_active;
            }
        }
    }
}
    




}  // namespace Kratos.
