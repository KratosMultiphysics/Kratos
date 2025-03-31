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
                        
                    auto pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );

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
                    
                auto pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
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
                        
                    auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
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
                auto pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", id_surrogate_condition, {{id_node_1, id_node_2}}, p_cond_prop );
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

}  // namespace Kratos.
