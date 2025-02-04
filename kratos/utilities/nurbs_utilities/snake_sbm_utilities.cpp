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
#include "snake_sbm_utilities.h"

namespace Kratos
{

    void SnakeSbmUtilities::CreateTheSnakeCoordinates(ModelPart& rIgaModelPart,
                                                    ModelPart& rSkinModelPartInnerInitial,
                                                    ModelPart& rSkinModelPartOuterInitial, 
                                                    ModelPart& rSkinModelPart, 
                                                    int rEchoLevel, 
                                                    Vector& knotVectorU, 
                                                    Vector& knotVectorV,
                                                    const Parameters mParameters) { 
        
        // Check
        KRATOS_ERROR_IF_NOT(mParameters.Has("sbm_parameters")) << "sbm_parameters has not been defined in the nurbs modeler" << std::endl;
        
        // Initilize the property of skin_model_part_in and out
        if (rSkinModelPartInnerInitial.Nodes().size()>0) {
        
            CreateTheSnakeCoordinates(rIgaModelPart, rSkinModelPartInnerInitial, rSkinModelPart, rEchoLevel, knotVectorU, knotVectorV, mParameters, true);
                
        }
        if (rSkinModelPartOuterInitial.Nodes().size()>0) {

            CreateTheSnakeCoordinates(rIgaModelPart, rSkinModelPartOuterInitial, rSkinModelPart, rEchoLevel, knotVectorU, knotVectorV, mParameters, false);
        }
    }   

    void SnakeSbmUtilities::CreateTheSnakeCoordinates(ModelPart& rIgaModelPart, 
                                                      ModelPart& rSkinModelPartInitial, 
                                                      ModelPart& rSkinModelPart,
                                                      int rEchoLevel, 
                                                      Vector& knotVectorU, 
                                                      Vector& knotVectorV,
                                                      const Parameters mParameters,
                                                      bool isInner) 
    { 
        
        std::string surrogate_sub_model_part_name; 
        std::string skin_sub_model_part_name; 
        // ModelPart skin_sub_model_part; 
        if (isInner) {
            surrogate_sub_model_part_name = "surrogate_inner";
            skin_sub_model_part_name = "inner";
        }
        else {
            surrogate_sub_model_part_name = "surrogate_outer";
            skin_sub_model_part_name = "outer";
        }
        
        ModelPart& skin_sub_model_part = rSkinModelPart.GetSubModelPart(skin_sub_model_part_name);
        ModelPart& surrogate_sub_model_part = rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name);

        if (!rSkinModelPartInitial.HasProperties(0)) rSkinModelPartInitial.CreateNewProperties(0);
        if (!rSkinModelPart.HasProperties(0)) rSkinModelPart.CreateNewProperties(0);
        Properties::Pointer p_cond_prop_in = rSkinModelPartInitial.pGetProperties(0);
        rSkinModelPartInitial.AddProperties(p_cond_prop_in);

        array_1d<double, 2> knot_step_uv(2);
        knot_step_uv[0] = std::abs(knotVectorU[int(knotVectorU.size()/2) +1]  - knotVectorU[int(knotVectorU.size()/2)] ) ;
        knot_step_uv[1] = std::abs(knotVectorV[int(knotVectorV.size()/2) +1]  - knotVectorV[int(knotVectorV.size()/2)] ) ;

        Vector meshSizes_uv(2);
        meshSizes_uv[0] = knot_step_uv[0]; 
        meshSizes_uv[1] = knot_step_uv[1];
        ModelPart& surrogate_model_part = rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name);
        surrogate_model_part.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uv);

        array_1d<double, 2> starting_pos_uv;
        starting_pos_uv[0] = knotVectorU[0];
        starting_pos_uv[1] = knotVectorV[0];

        Vector parameterExternalCoordinates(4);
        parameterExternalCoordinates[0] = knotVectorU[0];
        parameterExternalCoordinates[1] = knotVectorV[0];
        parameterExternalCoordinates[2] = knotVectorU[knotVectorU.size()-1];
        parameterExternalCoordinates[3] = knotVectorV[knotVectorV.size()-1];
        
        surrogate_model_part.GetProcessInfo().SetValue(LOAD_MESHES, parameterExternalCoordinates);

        // Create the matrix of active/inactive knot spans, one for inner and one for outer loop
        unsigned int numberOfLoops;
        if (isInner)
            numberOfLoops = mParameters["sbm_parameters"]["number_of_inner_loops"].GetInt();
        else 
            numberOfLoops = 1;
        
        std::vector<int> n_knot_spans_uv(2);
        n_knot_spans_uv[0] = knotVectorU.size()-1; 
        n_knot_spans_uv[1] = knotVectorV.size()-1;

        std::vector<std::vector<std::vector<int>>> knot_spans_available;
        knot_spans_available.reserve(numberOfLoops);

        for (IndexType i = 0; i < numberOfLoops; ++i) {
            std::vector<std::vector<int>> matrix; 
            matrix.reserve(n_knot_spans_uv[1]);
            for (int j = 0; j <= n_knot_spans_uv[1]-1; ++j) {
                std::vector<int> row(n_knot_spans_uv[0]); 
                matrix.push_back(row); 
            }
            knot_spans_available.push_back(matrix);
        }
        
        // Optimized Snake -> for inner loops
        int idMatrixKnotSpansAvailable = 0;
        IndexType idFirstNode;
        bool newInnerLoop = true;

        
        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && isInner) << "Inner :: Starting SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && !isInner) << "Outer :: Starting SnakeStep" << std::endl;
                
        if (rSkinModelPartInitial.Conditions().size() > 0) {
            
            // CREATE FIRST NODE FOR SKIN SUB MODEL PART
            auto initial_condition = rSkinModelPartInitial.GetCondition(1);
            double x_true_boundary0 = initial_condition.GetGeometry()[0].X();
            double y_true_boundary0 = initial_condition.GetGeometry()[0].Y();

            const int id_first_node = rSkinModelPart.GetRootModelPart().Nodes().size()+1;
            skin_sub_model_part.CreateNewNode(id_first_node, x_true_boundary0, y_true_boundary0, 0.0);

            for (auto &i_cond : rSkinModelPartInitial.Conditions()) {  
                if (newInnerLoop) {
                    idFirstNode = i_cond.GetGeometry()[0].Id();
                    newInnerLoop = false;
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

                // In the inner case, check is the immersed object is inside the rectangular domain
                if (isInner &&
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
                
                SnakeStep(skin_sub_model_part, knot_spans_available, idMatrixKnotSpansAvailable, 
                         knot_span_uv, xy_coord_i_cond, knot_step_uv, starting_pos_uv);
              
                if (i_cond.GetGeometry()[1].Id() == idFirstNode) {
                    idMatrixKnotSpansAvailable++;
                    newInnerLoop = true;
                }
            }
        }

        PointVector points;
        for (auto &i_cond : skin_sub_model_part.Conditions()) {
            points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry()[0].X(), i_cond.GetGeometry()[0].Y(), i_cond.GetGeometry()[0].Z())));
        }
        DynamicBins testBins(points.begin(), points.end());

        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && isInner) << "Inner :: Ending SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && !isInner) << "Outer :: Ending SnakeStep" << std::endl;
        
        // Read lambda parameters: 0.0 -> External,  0.5 -> Optimal
        double lambda;
        if (isInner) 
            if (mParameters["sbm_parameters"].Has("lambda_inner")) 
                lambda = mParameters["sbm_parameters"]["lambda_inner"].GetDouble();
            else 
                lambda = 0.5;
        else 
            if (mParameters["sbm_parameters"].Has("lambda_outer")) 
                lambda = mParameters["sbm_parameters"]["lambda_outer"].GetDouble();
            else 
                lambda = 0.5;
        
        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && isInner) << "Inner :: MarkKnotSpansAvailable" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && !isInner) << "Outer :: MarkKnotSpansAvailable" << std::endl;

        for (IndexType i = 0; i < numberOfLoops; i++) {
            IndexType idInnerLoop = i;
            // Mark the knot_spans_available's for inner and outer loops
            MarkKnotSpansAvailable(knot_spans_available, idInnerLoop, testBins, skin_sub_model_part, lambda, 
                                   n_knot_spans_uv, knot_step_uv, starting_pos_uv);  
            
            if (isInner) {
                CreateSurrogateBuondaryFromSnakeInner (knot_spans_available, idInnerLoop, surrogate_sub_model_part, 
                                    n_knot_spans_uv, knotVectorU, knotVectorV, starting_pos_uv );
                KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0) << "Inner :: Snake process has finished" << std::endl;
            }
            else {
                CreateSurrogateBuondaryFromSnakeOuter (testBins, skin_sub_model_part, knot_spans_available, idInnerLoop, surrogate_sub_model_part, 
                                    n_knot_spans_uv, knotVectorU, knotVectorV, starting_pos_uv);
                KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0) << "Outer :: Snake process has finished" << std::endl;
            }
        }

        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && isInner) << "Inner :: Loop finished" << std::endl;
        KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0 && !isInner) << "Outer :: Loop finished" << std::endl;
    }


    void SnakeSbmUtilities::SnakeStep(ModelPart& rSkinModelPart, 
                            std::vector<std::vector<std::vector<int>>> &knotSpansAvailable, 
                            int idMatrix, 
                            std::vector<std::vector<int>> knotSpansUV, 
                            std::vector<std::vector<double>> conditionCoord, 
                            Vector knotStepUV, 
                            Vector startingPos)
                            {
        bool isSplitted = false;

        if (knotSpansUV[0][0] != knotSpansUV[0][1] || knotSpansUV[1][0] != knotSpansUV[1][1]) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
            // Check if we are jumping some cut knot spans. If yes we split the true segment
            if (std::abs(knotSpansUV[1][0]-knotSpansUV[1][1]) > 1 || std::abs(knotSpansUV[0][0]-knotSpansUV[0][1]) > 1 || 
                    (knotSpansUV[0][0] != knotSpansUV[0][1] && knotSpansUV[1][0] != knotSpansUV[1][1]) ) {
                isSplitted = true;

                // Split the segment and do it recursively
                double x_true_boundary_split = (conditionCoord[0][0]+conditionCoord[0][1]) / 2;
                double y_true_boundary_split = (conditionCoord[1][0]+conditionCoord[1][1]) / 2;
                int knot_span_u_point_split = (x_true_boundary_split-startingPos[0]) / knotStepUV[0] ;
                int knot_span_v_point_split = (y_true_boundary_split-startingPos[1]) / knotStepUV[1] ;

                if (knot_span_u_point_split == int (knotSpansAvailable[idMatrix][0].size())) knot_span_u_point_split--;
                if (knot_span_v_point_split == int (knotSpansAvailable[idMatrix].size())) knot_span_v_point_split--;

                // update xy_coord for the first split segment
                std::vector<std::vector<double>> xy_coord_i_cond_split(2);
                xy_coord_i_cond_split[0].resize(2); xy_coord_i_cond_split[1].resize(2); 
                xy_coord_i_cond_split[0][0] = conditionCoord[0][0]; // x_true_boundary1
                xy_coord_i_cond_split[1][0] = conditionCoord[1][0]; // y_true_boundary1
                xy_coord_i_cond_split[0][1] = x_true_boundary_split; // x_true_boundary_split
                xy_coord_i_cond_split[1][1] = y_true_boundary_split; // y_true_boundary_split
                // update knot_span_uv for the first split segment
                std::vector<std::vector<int>> knot_span_uv_split(2);
                knot_span_uv_split[0].resize(2); knot_span_uv_split[1].resize(2); 
                knot_span_uv_split[0][0] = knotSpansUV[0][0]; // knot_span_u_1st_point
                knot_span_uv_split[1][0] = knotSpansUV[1][0]; // knot_span_v_1st_point
                knot_span_uv_split[0][1] = knot_span_u_point_split; // knot_span_u_point_split
                knot_span_uv_split[1][1] = knot_span_v_point_split; // knot_span_v_point_split
                
                // __We do it recursively first split__
                SnakeStep(rSkinModelPart, knotSpansAvailable, idMatrix, knot_span_uv_split, 
                         xy_coord_i_cond_split, knotStepUV, startingPos );

                // update xy_coord for the second split segment
                xy_coord_i_cond_split[0][0] = x_true_boundary_split; // x_true_boundary_split
                xy_coord_i_cond_split[1][0] = y_true_boundary_split; // y_true_boundary_split
                xy_coord_i_cond_split[0][1] = conditionCoord[0][1]; // x_true_boundary2
                xy_coord_i_cond_split[1][1] = conditionCoord[1][1]; // y_true_boundary2
                // update knot_span_uv for the first split segment
                knot_span_uv_split[0][0] = knot_span_u_point_split; // knot_span_u_point_split
                knot_span_uv_split[1][0] = knot_span_v_point_split; // knot_span_v_point_split
                knot_span_uv_split[0][1] = knotSpansUV[0][1]; // knot_span_u_2nd_point
                knot_span_uv_split[1][1] = knotSpansUV[1][1]; // knot_span_v_2nd_point

                // __We do it recursively second split__
                SnakeStep(rSkinModelPart, knotSpansAvailable, idMatrix, knot_span_uv_split, 
                         xy_coord_i_cond_split, knotStepUV, startingPos);
            }
            // Check if the true boundary crosses an u or a v knot value
            else if (knotSpansUV[0][0] != knotSpansUV[0][1]) { // u knot value is crossed
                // Find the "knotSpansAvailable" using the intersection
                knotSpansAvailable[idMatrix][knotSpansUV[1][0]][knotSpansUV[0][0]] = 2;
                knotSpansAvailable[idMatrix][knotSpansUV[1][0]][knotSpansUV[0][1]] = 2;
            }
            else if (knotSpansUV[1][0] != knotSpansUV[1][1]) { // v knot value is crossed
                // Find the "knotSpansAvailable" using the intersection (Snake_coordinate classic -> External Boundary)
                knotSpansAvailable[idMatrix][knotSpansUV[1][0]][knotSpansUV[0][0]] = 2;
                knotSpansAvailable[idMatrix][knotSpansUV[1][1]][knotSpansUV[0][0]] = 2;
            }
        }
        if (!isSplitted) {
            // Call the root model part for the Ids of the node
            auto idNode1 = rSkinModelPart.GetRootModelPart().Nodes().size();
            auto idNode2 = idNode1+1;
            // Create two nodes and two conditions for each skin condition
            rSkinModelPart.CreateNewNode(idNode2, (conditionCoord[0][0]+conditionCoord[0][1] ) / 2, (conditionCoord[1][0]+conditionCoord[1][1] ) / 2, 0.0);
            rSkinModelPart.CreateNewNode(idNode2+1, conditionCoord[0][1], conditionCoord[1][1], 0.0);
            Properties::Pointer p_cond_prop = rSkinModelPart.pGetProperties(0);
            Condition::Pointer p_cond1 = rSkinModelPart.CreateNewCondition("LineCondition2D2N", idNode1, {{idNode1, idNode2}}, p_cond_prop );
            Condition::Pointer p_cond2 = rSkinModelPart.CreateNewCondition("LineCondition2D2N", idNode2, {{idNode2, idNode2+1}}, p_cond_prop );
            rSkinModelPart.AddCondition(p_cond1);
            rSkinModelPart.AddCondition(p_cond2);
        }
    }


    bool SnakeSbmUtilities::IsPointInsideSkinBoundary(Point& point1, DynamicBins& testBins, ModelPart& rSkinModelPart)
    {
        // Get the nearest point of the true boundary
        PointerType pointToSearch = PointerType(new PointType(1000000, point1.X(), point1.Y(), 0.0));
        PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);
        
        // Get the closest Condition the initial_skin_model_part_in.Conditions
        IndexType id1 = nearestPoint->Id();
        auto nearestCondition1 = rSkinModelPart.GetCondition(id1);
        // Check if the condition is the first one and therefore the previous one does not exist
        IndexType id2 = id1 - 1;
        if (id1 == rSkinModelPart.ConditionsBegin()->Id()) {
            int nConditions = rSkinModelPart.Conditions().size();
            id2 = id1 + nConditions - 1; 
        }
        auto nearestCondition2 = rSkinModelPart.GetCondition(id2);
        // The two candidates nodes
        Point candidatePoint1 = nearestCondition1.GetGeometry()[1];
        Point candidatePoint2 = nearestCondition2.GetGeometry()[0];

        Point point2 = nearestCondition1.GetGeometry()[0]; // FIRST POINT IN TRUE GEOM
        Point point3 = candidatePoint1;// SECOND POINT IN TRUE GEOM
        if (MathUtils<double>::Norm(candidatePoint1-point1) > MathUtils<double>::Norm(candidatePoint2-point1)){
            // Need to invert the order to preserve the positivity of the area
            point3 = point2;
            point2 = candidatePoint2;
        }

        array_1d<double,3> v_1 = point2 - point1;
        array_1d<double,3> v_2 = point3 - point1;
        array_1d<double,3> crossProduct;
        MathUtils<double>::CrossProduct(crossProduct, v_1, v_2);

        bool isInside = false;
        if (crossProduct[2] > 0) {isInside = true;}

        return isInside;
    }



    void SnakeSbmUtilities::MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, int idMatrix,DynamicBins& testBin, 
                                                   ModelPart& rSkinModelPart, double lambda, std::vector<int>& nKnotSpans, 
                                                   array_1d<double, 2>& knotStepUV, const Vector& startingPos) {
        for (int i = 0; i < nKnotSpans[1]; i++) {
            for (int j = 0; j < nKnotSpans[0]; j++) {
                if (knotSpansAvailable[idMatrix][i][j] == 2) {
                    // Check the 8 neighbor knot spans -> Is there any completely inside? Note that we can just check 1 point.
                    
                    // right node
                    if (i != nKnotSpans[1]-1)
                        if (knotSpansAvailable[idMatrix][i+1][j] == 0) { 
                            Point gaussPoint = Point((j+0.5) * knotStepUV[0] + startingPos[0], (i+1+0.5) * knotStepUV[1] +startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i+1][j] = 1;}
                        }
                    // left node    
                    if (i != 0)
                        if (knotSpansAvailable[idMatrix][i-1][j] == 0) { 
                            Point gaussPoint = Point((j+0.5) * knotStepUV[0]+startingPos[0], (i-1+0.5) * knotStepUV[1] + startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i-1][j] = 1;}
                        }
                    // up node
                    if (j != nKnotSpans[0]-1)
                        if (knotSpansAvailable[idMatrix][i][j+1] == 0) { 
                            Point gaussPoint = Point((j+1+0.5) * knotStepUV[0]+startingPos[0], (i+0.5) * knotStepUV[1]+startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i][j+1] = 1;}
                        }
                    //down node
                    if (j != 0)
                        if (knotSpansAvailable[idMatrix][i][j-1] == 0) { 
                            Point gaussPoint = Point((j-1+0.5) * knotStepUV[0]+startingPos[0], (i+0.5) * knotStepUV[1]+startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i][j-1] = 1;}
                        } 

                    // corner right-down node
                    if (j != 0 && i != nKnotSpans[1]-1)
                        if (knotSpansAvailable[idMatrix][i+1][j-1] == 0) {
                            Point gaussPoint = Point((j-1+0.5) * knotStepUV[0]+startingPos[0], (i+1+0.5) * knotStepUV[1]+startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i+1][j-1] = 1;}
                        }
                    // corner left-down node
                    if (j != 0 && i != 0)
                        if (knotSpansAvailable[idMatrix][i-1][j-1] == 0) {
                            Point gaussPoint = Point((j-1+0.5) * knotStepUV[0]+startingPos[0], (i-1+0.5) * knotStepUV[1]+startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i-1][j-1] = 1;}
                        }
                    // corner right-up node
                    if (j != nKnotSpans[0]-1 && i != nKnotSpans[1]-1)
                        if (knotSpansAvailable[idMatrix][i+1][j+1] == 0) {
                            Point gaussPoint = Point((j+1+0.5) * knotStepUV[0]+startingPos[0], (i+1+0.5) * knotStepUV[1]+startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i+1][j+1] = 1;}
                        }
                    // corner left-up node
                    if (j != nKnotSpans[0]-1 && i != 0)
                        if (knotSpansAvailable[idMatrix][i-1][j+1] == 0) {
                            Point gaussPoint = Point((j+1+0.5) * knotStepUV[0]+startingPos[0], (i-1+0.5) * knotStepUV[1]+startingPos[1], 0);
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {knotSpansAvailable[idMatrix][i-1][j+1] = 1;}
                        }

                    // Create 25 "fake" GaussPoints to check if the majority are inside or outside
                    const int numFakeGaussPoints = 5;
                    int numberOfInsideGaussianPoints = 0;
                    for (IndexType i_GPx = 0; i_GPx < numFakeGaussPoints; i_GPx++){
                        double x_coord = j*knotStepUV[0] + knotStepUV[0]/(numFakeGaussPoints+1)*(i_GPx+1) + startingPos[0];

                        // NOTE:: The v-knot spans are upside down in the matrix!!
                        for (IndexType i_GPy = 0; i_GPy < numFakeGaussPoints; i_GPy++) 
                        {
                            double y_coord = i*knotStepUV[1] + knotStepUV[1]/(numFakeGaussPoints+1)*(i_GPy+1) + startingPos[1];
                            Point gaussPoint = Point(x_coord, y_coord, 0);  // GAUSSIAN POINT
                            if (IsPointInsideSkinBoundary(gaussPoint, testBin, rSkinModelPart)) {
                                // Sum over the number of numFakeGaussPoints per knot span
                                numberOfInsideGaussianPoints++;
                            }
                        }
                        
                    }
                
                    // Mark the knot span as available or not depending on the number of Gauss Points Inside/Outside
                    if (numberOfInsideGaussianPoints < lambda*numFakeGaussPoints*numFakeGaussPoints) {
                        knotSpansAvailable[idMatrix][i][j] = -1; // Cut knot spans that have been checked
                    }
                    else{
                        knotSpansAvailable[idMatrix][i][j] = 1; // The knot span is considered DEACTIVE
                    }
                }
            }
        }
    }


    void SnakeSbmUtilities::CreateSurrogateBuondaryFromSnakeInner (std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, int idMatrix, 
                            ModelPart& rSurrogateModelPartInner, std::vector<int>& nKnotSpans, 
                            Vector& knotVectorU, Vector&  knotVectorV,
                            const Vector& startingPos){
        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        for (int i = 0; i < (nKnotSpans[0]); i++) {
            for (int j = 0; j < (nKnotSpans[1]); j++ ) {
                if (knotSpansAvailable[idMatrix][j][i] == 1 ) {
                    // Check if one of the neighouring knot span is also == 1,
                    // otherwise the knot span is isolated I cannot be taken.
                    if ((knotSpansAvailable[idMatrix][j+1][i] == 1) || (knotSpansAvailable[idMatrix][j][i+1] == 1))
                    {
                        start_i = i;
                        start_j = j;
                        break;
                    }     
                }
            }
            if (knotSpansAvailable[idMatrix][start_j][start_i] == 1 ) { break; }
        }
        
        if (!rSurrogateModelPartInner.HasProperties(0)) {rSurrogateModelPartInner.CreateNewProperties(0);}

        Properties::Pointer p_cond_prop = rSurrogateModelPartInner.pGetProperties(0);
        Node snakeNode(1 , knotVectorU[start_i], knotVectorV[start_j], 0.0);
        
        const int id_first_node = rSurrogateModelPartInner.GetRootModelPart().Nodes().size() + 1;
        rSurrogateModelPartInner.CreateNewNode(id_first_node, snakeNode);
        IndexType idSnakeNode = id_first_node+1;

        // Follow the clockwise loop
        bool end = false;
        // We are going horizontally
        int direction = 0;
        // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int i = start_i; int j = start_j;
        int I = start_i; int J = start_j;
        int steps = 0;

        const int max_number_of_steps = 1e5;
        while (!end && steps < max_number_of_steps) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knotSpansAvailable[idMatrix][J][I] == 1) {
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; reset and move
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}
                if (knotSpansAvailable[idMatrix][J][I] == 1) {
                    // going straight
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    Node snakeNode(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                    rSurrogateModelPartInner.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                    idSnakeNode++;
                }
                else {
                    // Need to search to hte right; reset and move
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}
                    if (knotSpansAvailable[idMatrix][J][I] == 1) {
                        // We are moving to the right -> First move straight, store, the move to the right (i,j), store again
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        Node snakeNode1(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                        rSurrogateModelPartInner.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        Node snakeNode2(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                        rSurrogateModelPartInner.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        direction = direction + 1;
                        if (direction == 4) {direction = 0;}
                    }
                    else { // Special case of "isolated" knot span to be circumnavigated
                        is_special_case = 0;
                        // reset and move (I,J) backward
                        if (direction == 0) {I-- ; J-- ;} 
                        if (direction == 1) {J++ ; I-- ;}
                        if (direction == 2) {I++ ; J++ ;}
                        if (direction == 3) {J-- ; I++ ;}
                        // Check
                        if (knotSpansAvailable[idMatrix][J][I] == 1) {
                            // First straight, store, then move to the right, store again, then move to the right and store again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }

                            Node snakeNode(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                            rSurrogateModelPartInner.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;

                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }

                            Node snakeNode2(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                            rSurrogateModelPartInner.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;

                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }

                            Node snakeNode3(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                            rSurrogateModelPartInner.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;

                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else
                            KRATOS_ERROR << "SnakeSbmUtilities:: " <<  "error in the Snakes Coordinates" << std::endl;
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i && J == start_j && is_special_case == 1 ) {
                // End of the while loop
                end = true;
                KRATOS_INFO("number of steps in the snake step") << steps << std::endl;
            }
        }
        KRATOS_ERROR_IF(steps >= max_number_of_steps-1) << "SnakeSbmUtilities:: " <<  "reached maximum number of steps" << std::endl;

        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        IndexType initialId = id_first_node;
        if (rSurrogateModelPartInner.Elements().size()>0) {
            // Check if it is not the first inner loop
            initialId = rSurrogateModelPartInner.GetElement(rSurrogateModelPartInner.GetRootModelPart().Elements().size()).GetGeometry()[1].Id()+1;
        }
        std::vector<ModelPart::IndexType> elem_nodes{initialId, idSnakeNode-1};
        rSurrogateModelPartInner.CreateNewElement("Element2D2N", rSurrogateModelPartInner.GetRootModelPart().Elements().size()+1, elem_nodes, p_cond_prop);
        }


    void SnakeSbmUtilities::CreateSurrogateBuondaryFromSnakeOuter (DynamicBins& testBinOuter, ModelPart& rInitialSkinModelPartOuter, 
                            std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, int idMatrix, 
                            ModelPart& rSurrogateModelPartOuter, std::vector<int>& nKnotSpans, 
                            Vector& knotVectorU, Vector&  knotVectorV,
                            const Vector& startingPos){
                                             
        // CHECK ALL THE EXTERNAL KNOT SPANS

        // LEFT BOUNDARY
        double knot_step_u = knotVectorU[1]-knotVectorU[0];
        double knot_step_v = knotVectorV[1]-knotVectorV[0];

        for (int i = 0; i<2; i++) {
            for (int j = 0; j < (nKnotSpans[0]); j++ ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPos[0], (i+0.5)*knot_step_v+startingPos[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rInitialSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }
        // TOP BOUNDARY
        for (int j = int (knotSpansAvailable[idMatrix][0].size()-1); j > int (knotSpansAvailable[idMatrix][0].size()-3); j--) {
            for (int i = 0; i < (nKnotSpans[1]); i++) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPos[0], (i+0.5)*knot_step_v+startingPos[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rInitialSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }
        // RIGHT BOUNDARY
        for (int i = int (knotSpansAvailable[idMatrix].size()-1); i > int (knotSpansAvailable[idMatrix].size()-3); i--) {
            for (int j = nKnotSpans[0]-1; j > -1; j-- ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPos[0], (i+0.5)*knot_step_v+startingPos[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rInitialSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }
        // BOTTOM BOUNDARY
        for (int j = 0; j<2; j++) {
            for (int i = nKnotSpans[1]-1; i > -1 ; i--) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPos[0], (i+0.5)*knot_step_v+startingPos[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rInitialSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }

        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        for (int i = 0; i < (nKnotSpans[0]); i++) {
            for (int j = 0; j < (nKnotSpans[1]); j++ ) {
                if (knotSpansAvailable[idMatrix][j][i] == 1 ) {
                    // Check if one of the neighouring knot span is also == 1,
                    // otherwise the knot span is isolated I cannot be taken.
                    if ((knotSpansAvailable[idMatrix][j+1][i] == 1) || (knotSpansAvailable[idMatrix][j][i+1] == 1))
                    {
                        start_i = i;
                        start_j = j;
                        break;
                    }     
                }
            }
            if (knotSpansAvailable[idMatrix][start_j][start_i] == 1 ) { break; }
        }
        
        // EXTEND THE MATRIX
        //-------------------
        std::vector<std::vector<int>> knot_spans_available_extended(nKnotSpans[1]+2, std::vector<int>(nKnotSpans[0]+2));

        for (IndexType i = 0; i < knotSpansAvailable[idMatrix].size(); i++){
            for (IndexType j = 0; j < knotSpansAvailable[idMatrix][0].size(); j++) {
                knot_spans_available_extended[i+1][j+1] = knotSpansAvailable[idMatrix][i][j]; 
            }
        }  
        
        Properties::Pointer p_cond_prop = rSurrogateModelPartOuter.CreateNewProperties(1001);
        Node snakeNode(1 , knotVectorU[start_i], knotVectorV[start_j], 0.0);

        const int id_first_node = rSurrogateModelPartOuter.GetRootModelPart().Nodes().size() + 1;
        
        rSurrogateModelPartOuter.CreateNewNode(id_first_node, snakeNode);
        IndexType idSnakeNode = id_first_node+1;
        
        // Follow the clockwise loop
        bool end = false;
        // We are going orizontally
        int direction = 0 ;      // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int I = start_i+1;
        int J = start_j+1; 
        int i = start_i;
        int j = start_j;
        int steps = 0;
        const int max_number_of_steps = 1e5;
        while (!end && steps < max_number_of_steps) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knot_spans_available_extended[J][I] == 1) {
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; reset and move
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}

                if (knot_spans_available_extended[J][I] == 1) {

                    // Going straight! -> Don't store and move (i,j)
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    Node snakeNode(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                    rSurrogateModelPartOuter.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                    idSnakeNode++;
                }
                else {
                    // Need to search to the right; Reset and move
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}

                    if (knot_spans_available_extended[J][I] == 1) {

                        // Going to the right! -> first go forward, store, more right (i,j), store again
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        Node snakeNode1(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                        rSurrogateModelPartOuter.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        Node snakeNode2(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                        rSurrogateModelPartOuter.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        direction = direction + 1;
                        if (direction == 4) {direction = 0;}
                    }
                    else { // Special case of "isolated" knot span to be circumnavigated
                        is_special_case = 0;
                        // Reset and move (I,J) "backward"
                        if (direction == 0) {I-- ; J-- ;} 
                        if (direction == 1) {J++ ; I-- ;}
                        if (direction == 2) {I++ ; J++ ;}
                        if (direction == 3) {J-- ; I++ ;}
                        // Check
                        if (knot_spans_available_extended[J][I] == 1) {
                            // First go forward, store, then move to right, then store again, then move to the right and store again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }
                            Node snakeNode(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                            rSurrogateModelPartOuter.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }
                            Node snakeNode2(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                            rSurrogateModelPartOuter.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }
                            Node snakeNode3(idSnakeNode , knotVectorU[i], knotVectorV[j], 0.0);
                            rSurrogateModelPartOuter.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else
                            KRATOS_ERROR << "SnakeSbmUtilities:: " <<  "error in the Snakes Coordinates" << std::endl;
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i+1 && J == start_j+1 && is_special_case == 1 ) {
                // End of the while loop
                end = true;
                }
        }
        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        std::vector<ModelPart::IndexType> elem_nodes{1, idSnakeNode-1};
        const ModelPart& r_root = rSurrogateModelPartOuter.GetRootModelPart();
        IndexType elem_id = r_root.Elements().empty()
                            ? 1
                            : r_root.Elements().back().Id() + 1;
        rSurrogateModelPartOuter.CreateNewElement("Element2D2N", elem_id, elem_nodes, p_cond_prop);
    }
}  // namespace Kratos.
