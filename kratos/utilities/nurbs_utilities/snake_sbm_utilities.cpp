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
        DynamicBins testBinInner(points.begin(), points.end());

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
            MarkKnotSpansAvailable(knot_spans_available, idInnerLoop, testBinInner, skin_sub_model_part, lambda, 
                                   n_knot_spans_uv, knot_step_uv, starting_pos_uv);  
            
            if (isInner) {
                CreateSurrogateBuondaryFromSnakeInner(knot_spans_available, idInnerLoop, surrogate_sub_model_part, 
                                                      skin_sub_model_part, testBinInner, n_knot_spans_uv, 
                                                      knotVectorU, knotVectorV);

                KRATOS_INFO_IF("::[SnakeSbmUtilities]::", rEchoLevel > 0) << "Inner :: Snake process has finished" << std::endl;
            }
            else {
                CreateSurrogateBuondaryFromSnakeOuter (knot_spans_available, idInnerLoop, surrogate_sub_model_part,
                                                        skin_sub_model_part, testBinInner, n_knot_spans_uv, knotVectorU,
                                                        knotVectorV, starting_pos_uv);
                                                        
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


    /**
     * Marking process:
     *   1) We set to 2 all the cut knot spans
     *   2) We check the 8 neighbor knot spans and set them to 1|0  if inside|outside
     *   3) We check all the cut knot spans and set them    to 1|-1 if inside|outside 
     */
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

    /**
     * summary of knot_spans_available:
            " 1"  -> interior knot spans                                  
            "-1"  -> exterior knot spans well checked
            " 0"  -> exterior knot spans OR very interior knot spans (more 
                     than one ks away from surrogate boundary)
     */
    void SnakeSbmUtilities::CreateSurrogateBuondaryFromSnakeInner(std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, int idMatrix, 
                            ModelPart& rSurrogateModelPartInner, ModelPart& rSkinModelPartInner, DynamicBins& testBinInner,
                            std::vector<int>& nKnotSpans, Vector& knotVectorU, Vector&  knotVectorV) {
        
        KRATOS_INFO("::[SnakeSBMUtilities]::") << "Inner :: Check layers in 2D" << std::endl;

        // Snake 2D works with a raycasting technique from each of the two directions

        const double knot_step_u = knotVectorU[1]-knotVectorU[0];
        const double knot_step_v = knotVectorV[1]-knotVectorV[0];
        
        
        IndexType idSurrogateFirstNode; 
        if (rSurrogateModelPartInner.NumberOfNodes() == 0)
        {
            idSurrogateFirstNode = rSurrogateModelPartInner.GetRootModelPart().NumberOfNodes() + 1;
            IndexType idSurrogateNode = idSurrogateFirstNode;
            for (int j = 0; j < nKnotSpans[1]; j++) {
                for (int i = 0; i < nKnotSpans[0]; i++) {
                    rSurrogateModelPartInner.CreateNewNode(idSurrogateNode, knotVectorU[i], knotVectorV[j], 0.0);
                    idSurrogateNode++;
                }
            }
        } else 
        {
            idSurrogateFirstNode = rSurrogateModelPartInner.GetRootModelPart().NumberOfNodes() - nKnotSpans[1]*nKnotSpans[0] + 1;
        }
        
        Properties::Pointer p_cond_prop = rSurrogateModelPartInner.pGetProperties(0);
        
        // Direction parallel to x
        IndexType idSurrogateCondition = rSurrogateModelPartInner.GetRootModelPart().NumberOfConditions() + 1;
        for (int j = 0; j < nKnotSpans[1]; j++) {
            bool checkNextPoint = false;
            /*  
                Formula to connect i,j to the id of the model_part
                id = idSurrogateFirstNode + [i + j*(n_knot_spans_uv[0]) + 1];
            */
            // move in the x direction
            for (int i = 0; i < nKnotSpans[0]; i++) {
                if (checkNextPoint) {
                    // Check i+1 point using isPointInsideSkinBoundary3D
                    Point centerPoint = Point((i + 0.5)*knot_step_u, (j + 0.5)*knot_step_v, 0.0);
                    bool isExiting = false;
                    if ( knotSpansAvailable[idMatrix][j][i] == 1 ) {
                        // the knot span was already been checked very well
                    }
                    else if (IsPointInsideSkinBoundary(centerPoint, testBinInner, rSkinModelPartInner)) {
                        // STILL INSIDE --> do not save nothing and update knotSpansAvailable 
                        if ( knotSpansAvailable[idMatrix][j][i] == -1) {
                            isExiting = true;
                        }
                        else {
                            knotSpansAvailable[idMatrix][j][i] = 1;
                        }
                    }
                    else {
                        isExiting = true;
                    }
                    if (isExiting) {
                        /* EXITING --> save last segment in direction x. i-th is the knot value. */
                        int node1_i = i; int node1_j = j;   
                        int node2_i = i; int node2_j = j+1; 

                        IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*nKnotSpans[0];
                        IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*nKnotSpans[0];
                            
                        Condition::Pointer pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, false);

                        // surrogate_model_part_inner.AddCondition(pcond);
                        idSurrogateCondition++;
                        checkNextPoint = false;
                    }
                    
                }
                else if (knotSpansAvailable[idMatrix][j][i] == 1) {
                    // ENTERING --> save first face in direction
                    int node1_i = i; int node1_j = j;   
                    int node2_i = i; int node2_j = j+1; 

                    IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*nKnotSpans[0]; 
                    IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*nKnotSpans[0];
                        
                    Condition::Pointer pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                    idSurrogateCondition++;
                    checkNextPoint = true;

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, true);
                }
            }
        }
        
        // Do the same for y direction, without isPointInsideSkinBoundary, since we have already done it
        // And it is not necessary do it again
        for (int i = 0; i < nKnotSpans[0]; i++) {
            
            bool checkNextPoint = false;
            /*  
                Formula to connect i,j,k to the id of the model_part
                id = idSurrogateFirstNode + [i + j*(n_knot_spans_uv[0]) + 1];
            */
            // move in the y direction
            for (int j = 0; j < nKnotSpans[1]; j++) {
                if (checkNextPoint) {
                    if (knotSpansAvailable[idMatrix][j][i] != 1) {
                        /* EXITING --> save last face in direction x. i-th is the knot value. */
                        int node1_i = i;   int node1_j = j;
                        int node2_i = i+1; int node2_j = j;

                        IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*nKnotSpans[0];
                        IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*nKnotSpans[0];
                            
                        Condition::Pointer pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, false);

                        // surrogate_model_part_inner.AddCondition(p_cond);
                        idSurrogateCondition++;
                        checkNextPoint = false;
                    } 
                }
                else if (knotSpansAvailable[idMatrix][j][i] == 1) {
                    // ENTERING --> save first face in direction
                    int node1_i = i;   int node1_j = j;
                    int node2_i = i+1; int node2_j = j;

                    IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]);
                    IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]);
                    Condition::Pointer pcond = rSurrogateModelPartInner.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                    // surrogate_model_part_inner.AddCondition(p_cond);
                    idSurrogateCondition++;
                    checkNextPoint = true;

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, true);

                }
            }
        }

        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        IndexType elem_id = rSurrogateModelPartInner.NumberOfElements()+1;
        std::vector<ModelPart::IndexType> elem_nodes{elem_id, idSurrogateCondition};
        rSurrogateModelPartInner.CreateNewElement("Element2D2N", rSurrogateModelPartInner.Elements().size()+1, elem_nodes, p_cond_prop);
    }



    void SnakeSbmUtilities::CreateSurrogateBuondaryFromSnakeOuter(std::vector<std::vector<std::vector<int>>> & knotSpansAvailable,
                                int idMatrix, ModelPart& rSurrogateModelPartOuter, ModelPart& rSkinModelPartOuter,
                                DynamicBins& testBinOuter, std::vector<int>& nKnotSpans, 
                                Vector& knotVectorU, Vector& knotVectorV, const Vector& startingPositionUV)
    {

        KRATOS_INFO("::[SnakeSbmUtilities]::") << "Outer :: Check layers in 2D" << std::endl;

        // CHECK ALL THE EXTERNAL KNOT SPANS

        // LEFT BOUNDARY
        double knot_step_u = knotVectorU[1]-knotVectorU[0];
        double knot_step_v = knotVectorV[1]-knotVectorV[0];

        for (int i = 0; i<2; i++) {
            for (int j = 0; j < (nKnotSpans[0]); j++ ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPositionUV[0], (i+0.5)*knot_step_v+startingPositionUV[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }
        // TOP BOUNDARY
        for (int j = int (knotSpansAvailable[idMatrix][0].size()-1); j > int (knotSpansAvailable[idMatrix][0].size()-3); j--) {
            for (int i = 0; i < (nKnotSpans[1]); i++) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPositionUV[0], (i+0.5)*knot_step_v+startingPositionUV[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }
        // RIGHT BOUNDARY
        for (int i = int (knotSpansAvailable[idMatrix].size()-1); i > int (knotSpansAvailable[idMatrix].size()-3); i--) {
            for (int j = nKnotSpans[0]-1; j > -1; j-- ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPositionUV[0], (i+0.5)*knot_step_v+startingPositionUV[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }
        // BOTTOM BOUNDARY
        for (int j = 0; j<2; j++) {
            for (int i = nKnotSpans[1]-1; i > -1 ; i--) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+startingPositionUV[0], (i+0.5)*knot_step_v+startingPositionUV[1], 0);
                if (IsPointInsideSkinBoundary(centroidKnotSpan, testBinOuter, rSkinModelPartOuter) && knotSpansAvailable[idMatrix][i][j] != -1) {
                    knotSpansAvailable[idMatrix][i][j] = 1;
                    }
            }
        }

        KRATOS_INFO("::[SnakeSbmUtilities]::") << "Outer :: Starting Creation of Surrogate_Model_Part_Outer" << std::endl;
        
        // Snake 2D works with a raycasting technique from each of the two directions
        IndexType idSurrogateFirstNode = rSurrogateModelPartOuter.GetRootModelPart().NumberOfNodes() + 1;
        IndexType idSurrogateNode = idSurrogateFirstNode;
        for (int j = 0; j < nKnotSpans[1]+1; j++) {
            for (int i = 0; i < nKnotSpans[0]+1; i++) {
                rSurrogateModelPartOuter.CreateNewNode(idSurrogateNode, knotVectorU[i], knotVectorV[j], 0.0);
                idSurrogateNode++;
            }
        }
        
        // Direction parallel to x
       
        IndexType idSurrogateCondition = rSurrogateModelPartOuter.GetRootModelPart().NumberOfConditions() + 1;
        Properties::Pointer p_cond_prop = rSurrogateModelPartOuter.pGetProperties(0);
        
        for (int j = 0; j < nKnotSpans[1]; j++) {
            
            bool checkNextPoint = false;
            /*  
                Formula to connect i,j,k to the id of the model_part
                id = idSurrogateFirstNode + [i + j*(n_knot_spans_uv[0]) + 1];
            */
            // move in the x direction
            for (int i = 0; i < nKnotSpans[0]; i++) {

                int node1_i; int node1_j;   
                int node2_i; int node2_j; 

                if (checkNextPoint) {
                    // Check i+1 point using isPointInsideSkinBoundary
                    Point centerPoint = Point((i + 0.5)*knot_step_u, (j + 0.5)*knot_step_v, 0.0);
                    bool isExiting = false;
                    node1_i = i; node1_j = j;   
                    node2_i = i; node2_j = j+1; 
                    if ( knotSpansAvailable[idMatrix][j][i] == 1 ) {
                        // the knot span has already been checked very well
                    }
                    else if (IsPointInsideSkinBoundary(centerPoint, testBinOuter, rSkinModelPartOuter)) {
                        // STILL INSIDE --> do not save nothing and update knot_spans_available 
                        if ( knotSpansAvailable[idMatrix][j][i] == -1) {
                            isExiting = true;
                        }
                        else {
                            knotSpansAvailable[idMatrix][j][i] = 1;
                        }
                    }
                    else {
                        isExiting = true;
                    }
                    if (isExiting) {
                        /* EXITING --> save last face in direction x. i-th is the knot value. */

                        IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]+1);
                        IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]+1);
                            
                        Condition::Pointer pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, false);

                        idSurrogateCondition++;
                        checkNextPoint = false;    
                    }
                    
                }
                else if (knotSpansAvailable[idMatrix][j][i] == 1) {
                    // ENTERING --> save first face in direction
                    int node1_i = i; int node1_j = j;   
                    int node2_i = i; int node2_j = j+1; 

                    IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]+1);
                    IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]+1);
                        
                    Condition::Pointer pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                    idSurrogateCondition++;
                    checkNextPoint = true;

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, true);
                }

                if (knotSpansAvailable[idMatrix][j][i] == 1 && i == nKnotSpans[0]-1) 
                {
                    // Check if we are at the end of the patch -> if yes close the surrogate boundary
                    int node1_i = i+1; int node1_j = j;   
                    int node2_i = i+1; int node2_j = j+1; 

                    IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]+1);
                    IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]+1);
                        
                    Condition::Pointer pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, false);
                    idSurrogateCondition++;
                    checkNextPoint = false;
                }
            }
        }

        // Do the same for y, without isPointInsideSkinBoundary, since we have already done it
        // And it is not necessary do it again
        for (int i = 0; i < nKnotSpans[0]; i++) {
            
            bool checkNextPoint = false;
            /*  
                Formula to connect i,j,k to the id of the model_part
                i + j*(nKnotSpansUV[0]) + k*(nKnotSpansUV[1])*(nKnotSpansUV[0]);
            */
            // move in the y direction
            for (int j = 0; j < nKnotSpans[1]; j++) {
                if (checkNextPoint) {
                    int node1_i; int node1_j;
                    int node2_i; int node2_j;
                    if (knotSpansAvailable[idMatrix][j][i] != 1) {
                        /* EXITING --> save last face in direction x. i-th is the knot value. */
                        node1_i = i;   node1_j = j;
                        node2_i = i+1; node2_j = j;

                        IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]+1);
                        IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]+1);
                            
                        Condition::Pointer pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, false);
                        idSurrogateCondition++;
                        checkNextPoint = false;
                    } 
                }
                else if (knotSpansAvailable[idMatrix][j][i] == 1) {
                    // ENTERING --> save first face in direction
                    int node1_i = i;   int node1_j = j;
                    int node2_i = i+1; int node2_j = j;

                    IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]+1);
                    IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]+1);
                    Condition::Pointer pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                    idSurrogateCondition++;
                    checkNextPoint = true;

                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, true);
                }

                if (knotSpansAvailable[idMatrix][j][i] == 1 && j == nKnotSpans[1]-1) 
                {
                    // Check if we are at the end of the patch -> if yes close the surrogate boundary
                    int node1_i = i;   int node1_j = j+1;   
                    int node2_i = i+1; int node2_j = j+1; 

                    IndexType id_node_1 = idSurrogateFirstNode + node1_i + node1_j*(nKnotSpans[0]+1);
                    IndexType id_node_2 = idSurrogateFirstNode + node2_i + node2_j*(nKnotSpans[0]+1);
                        
                    Condition::Pointer pcond = rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", idSurrogateCondition, {{id_node_1, id_node_2}}, p_cond_prop );
                    // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                    pcond->Set(BOUNDARY, false);
                    idSurrogateCondition++;
                    checkNextPoint = false;
                }
            }
        }
    }

}  // namespace Kratos.
