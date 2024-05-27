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

    void SnakeSBMUtilities::CreateTheSnakeCoordinates(ModelPart& iga_model_part, ModelPart& skin_model_part_in, ModelPart& skin_model_part_out, 
                                                      ModelPart& initial_skin_model_part_in, ModelPart& initial_skin_model_part_out,  int rEchoLevel, 
                                                      Vector& knot_vector_u_complete, Vector& knot_vector_v_complete, double& knot_step_u, double& knot_step_v, 
                                                      const Parameters refinements_parameters, const Parameters mParameters, ModelPart& surrogate_model_part_inner, 
                                                      ModelPart& surrogate_model_part_outer) {

        
        // Initilize the property of skin_model_part_in and out
        if (!initial_skin_model_part_in.HasProperties(0)) {initial_skin_model_part_in.CreateNewProperties(0);}
        if (!initial_skin_model_part_out.HasProperties(0)) {initial_skin_model_part_out.CreateNewProperties(0);}
        Properties::Pointer p_cond_prop_in = initial_skin_model_part_in.pGetProperties(0);
        skin_model_part_in.AddProperties(p_cond_prop_in);
        Properties::Pointer p_cond_prop_out = initial_skin_model_part_out.pGetProperties(0);
        skin_model_part_out.AddProperties(p_cond_prop_out);
        

        int insert_nb_per_span_u = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        int insert_nb_per_span_v = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_v"].GetInt();

        int insert_nb_per_span_u_refined_out = insert_nb_per_span_u;
        int insert_nb_per_span_v_refined_out = insert_nb_per_span_v; 

        int insert_nb_per_span_u_refined_in = insert_nb_per_span_u;
        int insert_nb_per_span_v_refined_in = insert_nb_per_span_v; 

        // Additional refinement (LOCAL & TENSORIAL)
        if (refinements_parameters["refinements"][0]["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
            insert_nb_per_span_u_refined_in = (insert_nb_per_span_u+1) * refinements_parameters["refinements"][0]["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt() - 1 ;
            insert_nb_per_span_v_refined_in = (insert_nb_per_span_v +1) * refinements_parameters["refinements"][0]["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt() - 1 ;
        }
        
        // Build the knot_vector without the repetitive knot values

        
        //------------ INNER REFINED -------------------
        std::vector<double> knot_vector_u_refined_in(1);
        knot_vector_u_refined_in[0] = knot_vector_u_complete[0];
        std::vector<double> knot_vector_v_refined_in(1);
        knot_vector_v_refined_in[0] = knot_vector_v_complete[0];
        double knot_step_u_refined_in = (knot_vector_u_complete[knot_vector_u_complete.size()-1]-knot_vector_u_refined_in[0]) / (insert_nb_per_span_u_refined_in+1) ;
        double knot_step_v_refined_in = (knot_vector_v_complete[knot_vector_v_complete.size()-1]-knot_vector_v_refined_in[0]) / (insert_nb_per_span_v_refined_in+1) ;
        Vector meshSizes_uv(2);
        meshSizes_uv[0] = knot_step_u_refined_in; 
        meshSizes_uv[1] = knot_step_v_refined_in;
        surrogate_model_part_inner.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uv);
        for (IndexType j = 1; j < insert_nb_per_span_u_refined_in + 1; ++j) {
            knot_vector_u_refined_in.push_back(knot_vector_u_refined_in[0]+ knot_step_u_refined_in * j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_v_refined_in + 1; ++j) {
            knot_vector_v_refined_in.push_back(knot_vector_v_refined_in[0]+ knot_step_v_refined_in *j);
        }
        knot_vector_u_refined_in.push_back(knot_vector_u_complete[knot_vector_u_complete.size()-1]) ;
        knot_vector_v_refined_in.push_back(knot_vector_v_complete[knot_vector_v_complete.size()-1]) ;


        //------------ OUTER REFINED -------------------
        std::vector<double> knot_vector_u_refined_out(1);
        knot_vector_u_refined_out[0] = knot_vector_u_complete[0];
        std::vector<double> knot_vector_v_refined_out(1);
        knot_vector_v_refined_out[0] = knot_vector_v_complete[0];
        double knot_step_u_refined_out = (knot_vector_u_complete[knot_vector_u_complete.size()-1]-knot_vector_u_refined_out[0]) / (insert_nb_per_span_u_refined_out+1) ;
        double knot_step_v_refined_out = (knot_vector_v_complete[knot_vector_v_complete.size()-1]-knot_vector_v_refined_out[0]) / (insert_nb_per_span_v_refined_out+1) ;
        
        meshSizes_uv[0] = knot_step_u_refined_out; 
        meshSizes_uv[1] = knot_step_v_refined_out;
        surrogate_model_part_outer.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uv);
        for (IndexType j = 1; j < insert_nb_per_span_u_refined_out + 1; ++j) {
            knot_vector_u_refined_out.push_back(knot_vector_u_refined_out[0]+ knot_step_u_refined_out * j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_v_refined_out + 1; ++j) {
            knot_vector_v_refined_out.push_back(knot_vector_v_refined_out[0]+ knot_step_v_refined_out *j);
        }
        knot_vector_u_refined_out.push_back(knot_vector_u_complete[knot_vector_u_complete.size()-1]) ;
        knot_vector_v_refined_out.push_back(knot_vector_v_complete[knot_vector_v_complete.size()-1]) ;


        Vector parameterExternalCoordinates(4);
        parameterExternalCoordinates[0] = knot_vector_u_refined_out[0];
        parameterExternalCoordinates[1] = knot_vector_v_refined_out[0];
        parameterExternalCoordinates[2] = knot_vector_u_refined_out[knot_vector_u_refined_out.size()-1];
        parameterExternalCoordinates[3] = knot_vector_v_refined_out[knot_vector_v_refined_out.size()-1];
        surrogate_model_part_inner.GetProcessInfo().SetValue(LOAD_MESHES, parameterExternalCoordinates);
        surrogate_model_part_outer.GetProcessInfo().SetValue(LOAD_MESHES, parameterExternalCoordinates);


        //------------------------ INNER -----------------------------
        PointVector points_in;
        for (auto &i_cond : initial_skin_model_part_in.Conditions()) {
            points_in.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry()[0].X(), i_cond.GetGeometry()[0].Y(), i_cond.GetGeometry()[0].Z())));
        }
        DynamicBins testBins_in(points_in.begin(), points_in.end());

        // Create the matrix of active/inactive knot spans, one for inner and one for outer loop
        const int numberOfInnerLoops = mParameters["sbm_parameters"]["number_of_inner_loops"].GetInt();

        std::vector<std::vector<std::vector<int>>> knot_spans_available_inner;
        knot_spans_available_inner.reserve(numberOfInnerLoops);
        for (int i = 0; i < numberOfInnerLoops; ++i) {
            std::vector<std::vector<int>> matrix; 
            matrix.reserve(insert_nb_per_span_u_refined_in+1);
            for (int j = 0; j <= insert_nb_per_span_u_refined_in; ++j) {
                std::vector<int> row(insert_nb_per_span_v_refined_in+1); 
                matrix.push_back(row); 
            }
            knot_spans_available_inner.push_back(matrix); // Add the matrix to your 3D vector
        }
        
        // Optimized Snake -> for inner loops
        int idMatrixKnotSpansAvailable = 0;
        int idFirstNode;
        bool newInnerLoop = true;
        KRATOS_INFO("::[SnakeSBMUtilities]::") << "Inner :: Starting SnakeStep" << std::endl;

        if (initial_skin_model_part_in.Conditions().size() > 0) {
            
            // CREATE FIRST NODE FOR SKIN MODEL PART
            auto initial_condition_in = initial_skin_model_part_in.GetCondition(1);
            double x_true_boundary0_in = initial_condition_in.GetGeometry()[0].X();
            double y_true_boundary0_in = initial_condition_in.GetGeometry()[0].Y();
            skin_model_part_in.CreateNewNode(1, x_true_boundary0_in, y_true_boundary0_in, 0.0);


            for (auto &i_cond : initial_skin_model_part_in.Conditions()) {  
                if (newInnerLoop) {
                    idFirstNode = i_cond.GetGeometry()[0].Id();
                    newInnerLoop = false;
                }

                double x_true_boundary1 = i_cond.GetGeometry()[0].X();
                double y_true_boundary1 = i_cond.GetGeometry()[0].Y();
                double x_true_boundary2 = i_cond.GetGeometry()[1].X();
                double y_true_boundary2 = i_cond.GetGeometry()[1].Y();
                
                // Find the intersections of the skin boundary with the knot values
                int knot_span_u_1st_point = x_true_boundary1 / knot_step_u_refined_in ;
                int knot_span_u_2nd_point = x_true_boundary2 / knot_step_u_refined_in ;
                int knot_span_v_1st_point = y_true_boundary1 / knot_step_v_refined_in ;
                int knot_span_v_2nd_point = y_true_boundary2 / knot_step_v_refined_in ;

                SnakeStep(skin_model_part_in, knot_spans_available_inner, idMatrixKnotSpansAvailable, knot_span_u_1st_point, knot_span_u_2nd_point, knot_span_v_1st_point, knot_span_v_2nd_point,  x_true_boundary1, x_true_boundary2, y_true_boundary1, y_true_boundary2, knot_step_u_refined_in, knot_step_v_refined_in );
                
                if (i_cond.GetGeometry()[1].Id() == idFirstNode) {
                    idMatrixKnotSpansAvailable++;
                    newInnerLoop = true;
                }
            }
        }
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner :: Ending SnakeStep" << std::endl;
        // Read lambda parameters: 0.0 -> External,  0.5 -> Optimal
        const double lambda_inner = mParameters["sbm_parameters"]["lambda_inner"].GetDouble();
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner-> MarkKnotSpansAvailable" << std::endl;
        for (int i = 0; i < numberOfInnerLoops; i++) {
            int idInnerLoop = i;
            // Mark the knot_spans_available's for inner and outer loops
            MarkKnotSpansAvailable(knot_spans_available_inner, idInnerLoop, testBins_in, initial_skin_model_part_in, lambda_inner, insert_nb_per_span_u_refined_in, insert_nb_per_span_v_refined_in, knot_step_u_refined_in, knot_step_v_refined_in);
            CreateSurrogateBuondaryFromSnake_inner (knot_spans_available_inner,idInnerLoop, surrogate_model_part_inner, insert_nb_per_span_u_refined_in, insert_nb_per_span_v_refined_in, knot_vector_u_refined_in, knot_vector_v_refined_in);
            KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner-> Loop finished" << std::endl;
        }
        
        double box_refinement_u_min = 1e16;
        double box_refinement_u_max = -1e16;
        double box_refinement_v_min = 1e16;
        double box_refinement_v_max = -1e16;
        for(ModelPart::NodesContainerType::iterator i_node =  surrogate_model_part_inner.NodesBegin() ; i_node != surrogate_model_part_inner.NodesEnd() ; i_node++)
	    {
	      if( i_node-> X() < box_refinement_u_min ) {box_refinement_u_min = i_node-> X();}
          if( i_node-> X() > box_refinement_u_max ) {box_refinement_u_max = i_node-> X();}
          if( i_node-> Y() < box_refinement_v_min ) {box_refinement_v_min = i_node-> Y();}
          if( i_node-> Y() > box_refinement_v_max ) {box_refinement_v_max = i_node-> Y();}
	    }
        
        Vector box_refinement_u_v_min_max(4);
        box_refinement_u_v_min_max[0] = box_refinement_u_min; 
        box_refinement_u_v_min_max[1] = box_refinement_u_max;
        box_refinement_u_v_min_max[2] = box_refinement_v_min; 
        box_refinement_u_v_min_max[3] = box_refinement_v_max;
        
        iga_model_part.GetProcessInfo().SetValue(MARKER_MESHES, box_refinement_u_v_min_max);
        
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner-> Last loop has finished" << std::endl;

        //--------------------------------------------------------------------------------
        //---------------------------------- OUTER ---------------------------------------
        //--------------------------------------------------------------------------------
        PointVector points_out;
        for (auto &i_cond : initial_skin_model_part_out.Conditions()) {
            points_out.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry()[0].X(), i_cond.GetGeometry()[0].Y(), i_cond.GetGeometry()[0].Z())));
        }
        DynamicBins testBins_out(points_out.begin(), points_out.end());
        const int numberOfOuterLoops = 1;
        std::vector<std::vector<std::vector<int>>> knot_spans_available_outer;
        knot_spans_available_outer.reserve(numberOfOuterLoops);
        for (int i = 0; i < numberOfOuterLoops; ++i) {
            std::vector<std::vector<int>> matrix; 
            matrix.reserve(insert_nb_per_span_u_refined_out+1);
            for (int j = 0; j <= insert_nb_per_span_u_refined_out; ++j) {
                std::vector<int> row(insert_nb_per_span_v_refined_out+1); 
                matrix.push_back(row); 
            }
            knot_spans_available_outer.push_back(matrix); 
        }

        // Id conditions which cut the parameter space
        std::vector<int> idConditionOuterIntersection;

        // CREATE FIRST NODE FOR SKIN MODEL PART
        if (initial_skin_model_part_out.Conditions().size() > 0) {
            
            // CREATE FIRST NODE FOR SKIN MODEL PART
            auto initial_condition_out = initial_skin_model_part_out.GetCondition(1);
            double x_true_boundary0_out = initial_condition_out.GetGeometry()[0].X();
            double y_true_boundary0_out = initial_condition_out.GetGeometry()[0].Y();
            skin_model_part_out.CreateNewNode(1, x_true_boundary0_out, y_true_boundary0_out, 0.0);
            
            // Same for outer loops
            for (auto &i_cond : initial_skin_model_part_out.Conditions()) {
                double x_true_boundary1 = i_cond.GetGeometry()[0].X();
                double y_true_boundary1 = i_cond.GetGeometry()[0].Y();
                double x_true_boundary2 = i_cond.GetGeometry()[1].X();
                double y_true_boundary2 = i_cond.GetGeometry()[1].Y();
                
                // Check if point1 and point2 or the current condition is outside of the parameter space (ESSENTIAL)
                bool isPoint1Out = ((x_true_boundary1 >= knot_vector_u_refined_out[knot_vector_u_refined_out.size()-1]) || (x_true_boundary1 <= knot_vector_u_refined_out[0]) || 
                                    (y_true_boundary1 >= knot_vector_v_refined_out[knot_vector_v_refined_out.size()-1]) || (y_true_boundary1 <= knot_vector_v_refined_out[0]));
                
                bool isPoint2Out = ((x_true_boundary2 >= knot_vector_u_refined_out[knot_vector_u_refined_out.size()-1]) || (x_true_boundary2 <= knot_vector_u_refined_out[0]) || 
                                    (y_true_boundary2 >= knot_vector_v_refined_out[knot_vector_v_refined_out.size()-1]) || (y_true_boundary2 <= knot_vector_v_refined_out[0]));
                // Don't consider all the conditions which are outside the parameter space (ESSENTIAL)
                if (isPoint1Out || isPoint2Out) {
                    // Point 2 out Point 1 in
                    if  (isPoint2Out && !isPoint1Out) {
                        idConditionOuterIntersection.push_back(i_cond.Id());
                    }
                    // Point 1 out Point 2 in
                    else if (isPoint1Out && !isPoint2Out){
                        idConditionOuterIntersection.push_back(i_cond.Id());
                    }
                    continue;
                }
                // Find the intersections of the skin boundary with the knot values
                int knot_span_u_1st_point = x_true_boundary1 / knot_step_u_refined_out ;
                int knot_span_u_2nd_point = x_true_boundary2 / knot_step_u_refined_out ;
                int knot_span_v_1st_point = y_true_boundary1 / knot_step_v_refined_out ;
                int knot_span_v_2nd_point = y_true_boundary2 / knot_step_v_refined_out ;
                SnakeStep(skin_model_part_out, knot_spans_available_outer, 0,  knot_span_u_1st_point, knot_span_u_2nd_point, knot_span_v_1st_point, knot_span_v_2nd_point,  x_true_boundary1, x_true_boundary2, y_true_boundary1, y_true_boundary2, knot_step_u_refined_out, knot_step_v_refined_out);
            }      
            // Read lambda parameters: 0.0 -> External,  0.5 -> Optimal
            const double lambda_outer = mParameters["sbm_parameters"]["lambda_outer"].GetDouble();
            idMatrixKnotSpansAvailable = 0;
            if (initial_skin_model_part_out.Nodes().size() > 0) {
                MarkKnotSpansAvailable(knot_spans_available_outer, idMatrixKnotSpansAvailable, testBins_out, initial_skin_model_part_out, lambda_outer, insert_nb_per_span_u_refined_out, insert_nb_per_span_v_refined_out, knot_step_u_refined_out, knot_step_v_refined_out);
                CreateSurrogateBuondaryFromSnake_outer (testBins_out, initial_skin_model_part_out, knot_spans_available_outer, idMatrixKnotSpansAvailable, surrogate_model_part_outer, insert_nb_per_span_u_refined_out, insert_nb_per_span_v_refined_out, knot_vector_u_refined_out, knot_vector_v_refined_out);
            }
        }
    }   


    void SnakeSBMUtilities::SnakeStep(ModelPart& skin_model_part, std::vector<std::vector<std::vector<int>>> &knot_spans_available, int idMatrix, int knot_span_u_1st_point, int knot_span_u_2nd_point, int knot_span_v_1st_point,int knot_span_v_2nd_point,
                double& x_true_boundary1, double& x_true_boundary2, double& y_true_boundary1, double& y_true_boundary2, double& knot_step_u, double& knot_step_v) {
        bool isSplitted = false;
        if (knot_span_u_1st_point != knot_span_u_2nd_point || knot_span_v_1st_point != knot_span_v_2nd_point) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
            // Check if we are jumping some cut knot spans. If yes we split the true segment
            if (std::abs(knot_span_v_1st_point-knot_span_v_2nd_point) > 1 || std::abs(knot_span_u_1st_point-knot_span_u_2nd_point) > 1 || 
                    (knot_span_u_1st_point != knot_span_u_2nd_point && knot_span_v_1st_point != knot_span_v_2nd_point) ) {
                isSplitted = true;
                double x_true_boundary_split = (x_true_boundary1+x_true_boundary2 ) / 2;
                double y_true_boundary_split = (y_true_boundary1+y_true_boundary2 ) / 2;
                int knot_span_u_point_split = x_true_boundary_split / knot_step_u ;
                int knot_span_v_point_split = y_true_boundary_split / knot_step_v ;
                // We do it recursively
                SnakeStep(skin_model_part, knot_spans_available, idMatrix, knot_span_u_1st_point,   knot_span_u_point_split, knot_span_v_1st_point  , knot_span_v_point_split, x_true_boundary1     , x_true_boundary_split, y_true_boundary1     , y_true_boundary_split, knot_step_u, knot_step_v );
                SnakeStep(skin_model_part, knot_spans_available, idMatrix, knot_span_u_point_split, knot_span_u_2nd_point,   knot_span_v_point_split, knot_span_v_2nd_point,   x_true_boundary_split, x_true_boundary2,      y_true_boundary_split, y_true_boundary2     , knot_step_u, knot_step_v );
            }
            // Check if the true boundary crosses an u or a v knot value
            else if (knot_span_u_1st_point != knot_span_u_2nd_point) { // u knot value is crossed
                // Find the "knot_spans_available" using the intersection
                knot_spans_available[idMatrix][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_v_1st_point][knot_span_u_2nd_point] = 2;
            }
            else if (knot_span_v_1st_point != knot_span_v_2nd_point) { // v knot value is crossed
                // Find the "knot_spans_available" using the intersection (Snake_coordinate classic -> External Boundary)
                knot_spans_available[idMatrix][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_v_2nd_point][knot_span_u_1st_point] = 2;
            }
        }
        if (!isSplitted) {
            // int idNode1 = skin_model_part.Nodes().size();
            // int idNode2 = idNode1+1;
            // skin_model_part.CreateNewNode(idNode2, x_true_boundary2, y_true_boundary2, 0.0);
            // Properties::Pointer p_cond_prop = skin_model_part.pGetProperties(0);
            // IndexType idSnakeNode = skin_model_part.Nodes().size()+1;
            // Condition::Pointer p_cond = skin_model_part.CreateNewCondition("LineCondition2D2N", idNode1, {{idNode1, idNode2}}, p_cond_prop );
            // skin_model_part.AddCondition(p_cond);

            // Create at least two conditions for each skin condition
            auto idNode1 = skin_model_part.Nodes().size();
            auto idNode2 = idNode1+1;
            skin_model_part.CreateNewNode(idNode2, (x_true_boundary1+x_true_boundary2 ) / 2, (y_true_boundary1+y_true_boundary2 ) / 2, 0.0);
            skin_model_part.CreateNewNode(idNode2+1, x_true_boundary2, y_true_boundary2, 0.0);
            Properties::Pointer p_cond_prop = skin_model_part.pGetProperties(0);
            Condition::Pointer p_cond1 = skin_model_part.CreateNewCondition("LineCondition2D2N", idNode1, {{idNode1, idNode2}}, p_cond_prop );
            Condition::Pointer p_cond2 = skin_model_part.CreateNewCondition("LineCondition2D2N", idNode2, {{idNode2, idNode2+1}}, p_cond_prop );
            skin_model_part.AddCondition(p_cond1);
            skin_model_part.AddCondition(p_cond2);
        }
    }


    bool SnakeSBMUtilities::isPointInsideSkinBoundary(Point& point1, DynamicBins& testBins, ModelPart& skin_model_part)
    {
        // Get the nearest point of the true boundary
        PointerType pointToSearch = PointerType(new PointType(10000, point1.X(), point1.Y(), 0.0));
        PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);

        // SearchInRadius -> TOO SLOW
        // _________________________________________________________________________________________________________________________________
        // Vector meshSizes_uv = rModel->GetModelPart("surrogate_model_part").GetProcessInfo().GetValue(MARKER_MESHES);
        // double meshSize = meshSizes_uv[0];
        // if (meshSizes_uv[1] > meshSize) {meshSize = meshSizes_uv[1];}
        // const double radius = sqrt(2.5)*(meshSize)+1e-15; 
        // const int numberOfResults = 1e6; 
        // ModelPart::NodesContainerType::ContainerType Results(numberOfResults);
        // std::vector<double> list_of_distances(numberOfResults);
        // int obtainedResults = testBins.SearchInRadius(*pointToSearch, radius, Results.begin(), list_of_distances.begin(), numberOfResults);
        // double minimum_distance=1e10;
        // int nearestNodeId;
        // for (int i_distance = 0; i_distance < obtainedResults; i_distance++) {
        //     double new_distance = list_of_distances[i_distance];   
        //     if (new_distance < minimum_distance) { 
        //         minimum_distance = new_distance;
        //         nearestNodeId = i_distance;
        //         }
        // }
        // if (obtainedResults == 0) { KRATOS_WATCH("0 POINTS FOUND: EXIT");  exit(0);}
        // int id1 = Results[nearestNodeId]->Id();
        // _________________________________________________________________________________________________________________________________
        
        // Get the closest Condition the initial_skin_model_part_in.Conditions
        int id1 = nearestPoint->Id();
        auto nearestCondition1 = skin_model_part.GetCondition(id1);
        // Check if the condition is the first one and therefore the previous one does not exist
        int id2 = id1 - 1;
        if (id1 == 1) {
            int nConditions = skin_model_part.Conditions().size();
            id2 = nConditions; 
        }
        auto nearestCondition2 = skin_model_part.GetCondition(id2);
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



    void SnakeSBMUtilities::MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idMatrix,DynamicBins& testBin, ModelPart& skin_model_part, double lambda, int insert_nb_per_span_u, int insert_nb_per_span_v, double knot_step_u, double knot_step_v) {
        
        for (int i = 0; i < insert_nb_per_span_v+1; i++) {
            for (int j = 0; j < insert_nb_per_span_u+1; j++) {
                if (knot_spans_available[idMatrix][i][j] == 2) {
                    // Check the 8 neighbor knot spans -> Is there any completely inside? Note that we can just check 1 point.
                    if (i != 0 && i != insert_nb_per_span_v) {
                        if (knot_spans_available[idMatrix][i+1][j] == 0) { // right node
                            Point gaussPoint = Point((j+0.5)*knot_step_u, (i+1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i+1][j] = 1;}
                        }
                        if (knot_spans_available[idMatrix][i-1][j] == 0) { // left node
                            Point gaussPoint = Point((j+0.5)*knot_step_u, (i-1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i-1][j] = 1;}
                        }
                    }
                    if (j != 0 && j != insert_nb_per_span_u) {
                        if (knot_spans_available[idMatrix][i][j+1] == 0) { // up node
                            Point gaussPoint = Point((j+1+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i][j+1] = 1;}
                        }
                        if (knot_spans_available[idMatrix][i][j-1] == 0) { //down node
                            Point gaussPoint = Point((j-1+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i][j-1] = 1;}
                        } 
                    }
                    if ((i != 0 && i != insert_nb_per_span_v) && (j != 0 && j != insert_nb_per_span_u)){
                        if (knot_spans_available[idMatrix][i+1][j-1] == 0) { // corner right-down node
                            Point gaussPoint = Point((j-1+0.5)*knot_step_u, (i+1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i+1][j-1] = 1;}
                        }
                        if (knot_spans_available[idMatrix][i-1][j-1] == 0) { // corner left-down node
                            Point gaussPoint = Point((j-1+0.5)*knot_step_u, (i-1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i-1][j-1] = 1;}
                        }
                        if (knot_spans_available[idMatrix][i+1][j+1] == 0) { // corner right-up node
                            Point gaussPoint = Point((j+1+0.5)*knot_step_u, (i+1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i+1][j+1] = 1;}
                        }
                        if (knot_spans_available[idMatrix][i-1][j+1] == 0) { // corner left-up node
                            Point gaussPoint = Point((j+1+0.5)*knot_step_u, (i-1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i-1][j+1] = 1;}
                        }
                    }

                    // Create 25 "fake" GaussPoints to check if the majority are inside or outside
                    const int numFakeGaussPoints = 4;
                    int numberOfInsideGaussianPoints = 0;
                    for (int i_GPx = 0; i_GPx < numFakeGaussPoints; i_GPx++){
                        double x_coord = j*knot_step_u + knot_step_u/(numFakeGaussPoints+1)*(i_GPx+1);

                        // NOTE:: The v-knot spans are upside down in the matrix!!
                        for (int i_GPy = 0; i_GPy < numFakeGaussPoints; i_GPy++) 
                        {
                            double y_coord = i*knot_step_v + knot_step_v/(numFakeGaussPoints+1)*(i_GPy+1);
                            Point gaussPoint = Point(x_coord, y_coord, 0);  // GAUSSIAN POINT
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {
                                // Sum over the number of numFakeGaussPoints per knot span
                                numberOfInsideGaussianPoints++;
                            }
                        }
                        
                    }
                    // Mark the knot span as available or not depending on the number of Gauss Points Inside/Outside
                    if (numberOfInsideGaussianPoints < lambda*numFakeGaussPoints*numFakeGaussPoints) {
                        knot_spans_available[idMatrix][i][j] = -1; // Cut knot spans that have been checked
                    }
                    else{
                        knot_spans_available[idMatrix][i][j] = 1; // The knot span is considered DEACTIVE
                    }
                            
                }
            }
        }
    }


    void SnakeSBMUtilities::CreateSurrogateBuondaryFromSnake_inner (std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idMatrix, ModelPart& surrogate_model_part_inner, int insert_nb_per_span_u, int insert_nb_per_span_v, std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v){
        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        for (int i = 0; i < (insert_nb_per_span_u+1); i++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                if (knot_spans_available[idMatrix][j][i] == 1 ) {
                    // Check if one of the neighouring knot span is also == 1,
                    // otherwise the knot span is isolated I cannot be taken.
                    if ((knot_spans_available[idMatrix][j+1][i] == 1) || (knot_spans_available[idMatrix][j][i+1] == 1))
                    {
                        start_i = i;
                        start_j = j;
                        break;
                    }     
                }
            }
            if (knot_spans_available[idMatrix][start_j][start_i] == 1 ) { break; }
        }
        
        std::ofstream outputFile("txt_files/Snake_coordinates.txt", std::ios::app);
        outputFile << std::setprecision(16);
        outputFile << knot_vector_u[start_i] << " " << knot_vector_v[start_j] << std::endl;
        
        if (!surrogate_model_part_inner.HasProperties(0)) {surrogate_model_part_inner.CreateNewProperties(0);}

        Properties::Pointer p_cond_prop = surrogate_model_part_inner.pGetProperties(0);

        Node snakeNode(1 , knot_vector_u[start_i], knot_vector_v[start_j], 0.0);
        

        surrogate_model_part_inner.CreateNewNode(surrogate_model_part_inner.Nodes().size()+1, snakeNode);

        IndexType idSnakeNode = surrogate_model_part_inner.Nodes().size()+1;


        // Follow the clockwise loop
        int end = 0;
        // We are going orizontally
        int direction = 0 ;      // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int i = start_i;         // Indici del percorso lungo i knot values
        int j = start_j;         // Indici del percorso lungo i knot values
        int I = start_i;         // Indici del percorso lungo i knot spans
        int J = start_j;         // Indici del percorso lungo i knot spans
        int steps = 0;
        while (end == 0 && steps < 100000) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // KRATOS_WATCH(i)
            // KRATOS_WATCH(j)
            // KRATOS_WATCH(I)
            // KRATOS_WATCH(J)
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knot_spans_available[idMatrix][J][I] == 1) {
                // KRATOS_WATCH("trovato, sinistra") 
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; risetto e muovo
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}
                // KRATOS_WATCH(I)
                // KRATOS_WATCH(J)
                if (knot_spans_available[idMatrix][J][I] == 1) {
                    // KRATOS_WATCH("trovato, dritto")

                    // Stiamo andando a Dritti! -> Non scrivo nulla e muovo (i,j)
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                    Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                    surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                    idSnakeNode++;
                
                    // KRATOS_WATCH(knot_vector_u[i])
                    // KRATOS_WATCH(knot_vector_v[j])
                }
                else {
                    // Dobbiamo cercare a destra; resetto e muovo
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}

                    if (knot_spans_available[idMatrix][J][I] == 1) {
                        // KRATOS_WATCH("trovato, destra")

                        // Stiamo andando a DX! -> Prima passo dritto, poi stampo, poi passo a destra (i,j), poi scrivo di nuovo
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                        Node snakeNode1(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        // KRATOS_WATCH(knot_vector_u[i])
                        // KRATOS_WATCH(knot_vector_v[j])
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                        Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        // KRATOS_WATCH(knot_vector_u[i])
                        // KRATOS_WATCH(knot_vector_v[j])
                        direction = direction + 1;
                        if (direction == 4) {direction = 0;}
                    }
                    else { // Super special case of "isolated" knot span to be circumnavigated
                        // KRATOS_WATCH("Super special case of "isolated" knot span")
                        is_special_case = 0;
                        // Resetto e muovo (I,J) "indietro"
                        if (direction == 0) {I-- ; J-- ;} 
                        if (direction == 1) {J++ ; I-- ;}
                        if (direction == 2) {I++ ; J++ ;}
                        if (direction == 3) {J-- ; I++ ;}
                        // Check
                        if (knot_spans_available[idMatrix][J][I] == 1) {
                            // First passo dritto, poi print, then move to right, then print again, then move to the right and print again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                            Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                            Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            // surrogate_model_part_inner.AddCondition(pcond2);
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                            Node snakeNode3(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            // surrogate_model_part.AddCondition(pcond3);
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else{ KRATOS_WATCH("errore nello Snakes Coordinates"); 
                            // exit(0);
                            }  
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i && J == start_j && is_special_case == 1 ) {
                // End of the while loop
                end = 1;
                KRATOS_WATCH(steps)
                // Symbol on output file to end a inner loop
                outputFile << "end_inner_loop" << std::endl; //DELETE
                }
                

        }
        if (steps > 90000) { 
            KRATOS_ERROR << "errore nello SNAKE: STEPS>100000";
            exit(0);
        }
        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        IndexType initialId = 1;
        if (surrogate_model_part_inner.Elements().size()>0) {
            // Check if it is not the first inner loop
            initialId = surrogate_model_part_inner.GetElement(surrogate_model_part_inner.Elements().size()).GetGeometry()[1].Id()+1;
        }
        std::vector<ModelPart::IndexType> elem_nodes{initialId, idSnakeNode-1};
        surrogate_model_part_inner.CreateNewElement("Element2D2N", surrogate_model_part_inner.Elements().size()+1, elem_nodes, p_cond_prop);

        // KRATOS_WA<TCH(initialId)
        // KRATOS_WATCH(idSnakeNode-1)
        outputFile.close();
        KRATOS_WATCH("Inner Snake process has finished")
    }


    void SnakeSBMUtilities::CreateSurrogateBuondaryFromSnake_outer (DynamicBins& testBin_out, ModelPart& initial_skin_model_part_out, std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idMatrix, ModelPart& surrogate_model_part_outer, int insert_nb_per_span_u, int insert_nb_per_span_v, std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v){

        KRATOS_WATCH("comincia Outer")
        // CHECK ALL THE EXTERNAL KNOT SPANS -> Better also an additional layer
        // LEFT BOUNDARY
        double knot_step_u = knot_vector_u[1]-knot_vector_u[0];
        double knot_step_v = knot_vector_v[1]-knot_vector_v[0];
        // int i = 0;
        for (int i = 0; i<2; i++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }
        // TOP BOUNDARY
        // int j = knot_spans_available.size()-1;
        for (int j = knot_spans_available[idMatrix][0].size()-1; j > knot_spans_available[idMatrix][0].size()-3; j--) {
            for (int i = 0; i < (insert_nb_per_span_u+1); i++) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }
        // RIGHT BOUNDARY
        // i = knot_spans_available[0].size()-1;
        for (int i = knot_spans_available[idMatrix][0].size()-1; i > knot_spans_available[idMatrix][0].size()-3; i--) {
            for (int j = insert_nb_per_span_v; j > -1; j-- ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }
        // BOTTOM BOUNDARY
        // j = 0;
        for (int j = 0; j<2; j++) {
            for (int i = insert_nb_per_span_u; i > -1 ; i--) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }


        // MODIFIED 09/02/2024
        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        IndexType idSnakeNode = 1;
        for (int i = 0; i < (insert_nb_per_span_u+1); i++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                if (knot_spans_available[idMatrix][j][i] == 1 ) {
                    // Check if one of the neighouring knot span is also == 1,
                    // otherwise the knot span is isolated I cannot be taken.
                    if ((knot_spans_available[idMatrix][j+1][i] == 1) || (knot_spans_available[idMatrix][j][i+1] == 1))
                    {
                        start_i = i;
                        start_j = j;
                        break;
                    }     
                }
            }
            // KRATOS_WATCH(knot_spans_available[start_j][start_i])
            if (knot_spans_available[idMatrix][start_j][start_i] == 1 ) { break; }
        }
        
        // EXTEND THE MATRIX
        //-------------------
        std::vector<std::vector<int>> knot_spans_available_extended(insert_nb_per_span_v+1+2, std::vector<int>(insert_nb_per_span_u+1+2));

        for (int i = 0; i < knot_spans_available[idMatrix].size(); i++){
            for (int j = 0; j<knot_spans_available[idMatrix][0].size(); j++) {
                knot_spans_available_extended[i+1][j+1] = knot_spans_available[idMatrix][i][j]; 
            }
        }  
        //--------------------

        std::ofstream outputFile("txt_files/Snake_coordinates_outer.txt");
        outputFile << std::setprecision(16);
        outputFile << knot_vector_u[start_i] << " " << knot_vector_v[start_j] << std::endl;
        
        Properties::Pointer p_cond_prop = surrogate_model_part_outer.CreateNewProperties(0);
        Node snakeNode(1 , knot_vector_u[start_i], knot_vector_v[start_j], 0.0);
        
        // surrogate_model_part.AddNode(&snakeNode);

        surrogate_model_part_outer.CreateNewNode(1, snakeNode);
        idSnakeNode++;
        // KRATOS_WATCH(knot_vector_u[start_i])
        // KRATOS_WATCH(knot_vector_v[start_j])

        
        // Follow the clockwise loop
        int end = 0;
        // We are going orizontally
        int direction = 0 ;      // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int I = start_i+1;         // Indici del percorso lungo i knot spans
        int J = start_j+1;         // Indici del percorso lungo i knot spans
        int i = start_i;
        int j = start_j;
        int steps = 0;
        while (end == 0 && steps < 100000) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // KRATOS_WATCH(i)
            // KRATOS_WATCH(j)
            // KRATOS_WATCH(I)
            // KRATOS_WATCH(J)
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knot_spans_available_extended[J][I] == 1) {
                // KRATOS_WATCH("trovato, sinistra") 
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; risetto e muovo
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}
                // KRATOS_WATCH(I)
                // KRATOS_WATCH(J)
                if (knot_spans_available_extended[J][I] == 1) {
                    // KRATOS_WATCH("trovato, dritto")

                    // Stiamo andando a Dritti! -> Non scrivo nulla e muovo (i,j)
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                    Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                    surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                    idSnakeNode++;
                
                    // KRATOS_WATCH(knot_vector_u[i])
                    // KRATOS_WATCH(knot_vector_v[j])
                }
                else {
                    // Dobbiamo cercare a destra; resetto e muovo
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}

                    if (knot_spans_available_extended[J][I] == 1) {
                        // KRATOS_WATCH("trovato, destra")

                        // Stiamo andando a DX! -> Prima passo dritto, poi stampo, poi passo a destra (i,j), poi scrivo di nuovo
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                        Node snakeNode1(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        // KRATOS_WATCH(knot_vector_u[i])
                        // KRATOS_WATCH(knot_vector_v[j])
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                        Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        // KRATOS_WATCH(knot_vector_u[i])
                        // KRATOS_WATCH(knot_vector_v[j])
                        direction = direction + 1;
                        if (direction == 4) {direction = 0;}
                    }
                    else { // Super special case of "isolated" knot span to be circumnavigated
                        // KRATOS_WATCH("Super special case of "isolated" knot span")
                        is_special_case = 0;
                        // Resetto e muovo (I,J) "indietro"
                        if (direction == 0) {I-- ; J-- ;} 
                        if (direction == 1) {J++ ; I-- ;}
                        if (direction == 2) {I++ ; J++ ;}
                        if (direction == 3) {J-- ; I++ ;}
                        // Check
                        if (knot_spans_available_extended[J][I] == 1) {
                            // First passo dritto, poi print, then move to right, then print again, then move to the right and print again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                            Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                            Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                            Node snakeNode3(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else{ KRATOS_WATCH("errore nello Snakes Coordinates"); 
                            exit(0);
                            }  
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i+1 && J == start_j+1 && is_special_case == 1 ) {
                // End of the while loop
                end = 1;
                KRATOS_WATCH(steps)
                }
            // if ((I == end_i && J == end_j) && is_special_case == 1 ) {
            //     // End of the while loop
            //     end = 1;
            //     KRATOS_WATCH(steps)
            //     // Need to write the last condition in surrogate_model_part_outer
            //     if (J==0) {j--;}
            //     if (I==0) {i--;}
            //     if (J==knot_vector_v.size()-2) {j++;}
            //     if (I==knot_vector_u.size()-2) {i++;}
            //     outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl;
            //     Node snakeNode3(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
            //     surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode3);
            //     Condition::Pointer pcond3 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
            //     idSnakeNode++;
            //     }
                

        }
        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        std::vector<ModelPart::IndexType> elem_nodes{1, idSnakeNode-1};
        surrogate_model_part_outer.CreateNewElement("Element2D2N", 1, elem_nodes, p_cond_prop);
        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl;

        outputFile.close();
        KRATOS_WATCH("Outer Snake process has finished")
    }

    // ######################################################################################
    // 3 - DIMENSIONAL SNAKE PROCESS
    // #######################################################################################
    void SnakeSBMUtilities::CreateTheSnakeCoordinates3D(ModelPart& iga_model_part, ModelPart& skin_model_part_in, ModelPart& skin_model_part_out, 
                                                      ModelPart& initial_skin_model_part_in, ModelPart& initial_skin_model_part_out,  int rEchoLevel, 
                                                      Vector& knot_vector_u_complete, Vector& knot_vector_v_complete, Vector& knot_vector_w_complete, double& knot_step_u, double& knot_step_v, double& knot_step_w,
                                                      const Parameters refinements_parameters, const Parameters mParameters, ModelPart& surrogate_model_part_inner, 
                                                      ModelPart& surrogate_model_part_outer) {

        
        // Initilize the property of skin_model_part_in and out
        if (!initial_skin_model_part_in.HasProperties(0)) {initial_skin_model_part_in.CreateNewProperties(0);}
        if (!initial_skin_model_part_out.HasProperties(0)) {initial_skin_model_part_out.CreateNewProperties(0);}
        Properties::Pointer p_cond_prop_in = initial_skin_model_part_in.pGetProperties(0);
        skin_model_part_in.AddProperties(p_cond_prop_in);
        Properties::Pointer p_cond_prop_out = initial_skin_model_part_out.pGetProperties(0);
        skin_model_part_out.AddProperties(p_cond_prop_out);
        

        int insert_nb_per_span_u = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        int insert_nb_per_span_v = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_v"].GetInt();
        int insert_nb_per_span_w = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_w"].GetInt();

        int insert_nb_per_span_u_refined_out = insert_nb_per_span_u;
        int insert_nb_per_span_v_refined_out = insert_nb_per_span_v; 
        int insert_nb_per_span_w_refined_out = insert_nb_per_span_w; 

        int insert_nb_per_span_u_refined_in = insert_nb_per_span_u;
        int insert_nb_per_span_v_refined_in = insert_nb_per_span_v; 
        int insert_nb_per_span_w_refined_in = insert_nb_per_span_w; 

        // Additional refinement (LOCAL & TENSORIAL)
        if (refinements_parameters["refinements"][0]["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
            insert_nb_per_span_u_refined_in = (insert_nb_per_span_u + 1) * refinements_parameters["refinements"][0]["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt() - 1 ;
            insert_nb_per_span_v_refined_in = (insert_nb_per_span_v + 1) * refinements_parameters["refinements"][0]["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt() - 1 ;
            insert_nb_per_span_w_refined_in = (insert_nb_per_span_w + 1) * refinements_parameters["refinements"][0]["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt() - 1 ;
        }
        
        // Build the knot_vector without the repetitive knot values

        
        //------------ INNER REFINED -------------------
        std::vector<double> knot_vector_u_refined_in(1);    knot_vector_u_refined_in[0] = knot_vector_u_complete[0];
        std::vector<double> knot_vector_v_refined_in(1);    knot_vector_v_refined_in[0] = knot_vector_v_complete[0];
        std::vector<double> knot_vector_w_refined_in(1);    knot_vector_w_refined_in[0] = knot_vector_w_complete[0];
        
        double knot_step_u_refined_in = (knot_vector_u_complete[knot_vector_u_complete.size()-1]-knot_vector_u_refined_in[0]) / (insert_nb_per_span_u_refined_in+1) ;
        double knot_step_v_refined_in = (knot_vector_v_complete[knot_vector_v_complete.size()-1]-knot_vector_v_refined_in[0]) / (insert_nb_per_span_v_refined_in+1) ;
        double knot_step_w_refined_in = (knot_vector_w_complete[knot_vector_w_complete.size()-1]-knot_vector_w_refined_in[0]) / (insert_nb_per_span_w_refined_in+1) ;
        
        Vector meshSizes_uvw(3);
        meshSizes_uvw[0] = knot_step_u_refined_in; 
        meshSizes_uvw[1] = knot_step_v_refined_in;
        meshSizes_uvw[2] = knot_step_w_refined_in;
        surrogate_model_part_inner.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uvw);
        for (IndexType j = 1; j < insert_nb_per_span_u_refined_in + 1; ++j) {
            knot_vector_u_refined_in.push_back(knot_vector_u_refined_in[0]+ knot_step_u_refined_in * j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_v_refined_in + 1; ++j) {
            knot_vector_v_refined_in.push_back(knot_vector_v_refined_in[0]+ knot_step_v_refined_in *j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_w_refined_in + 1; ++j) {
            knot_vector_w_refined_in.push_back(knot_vector_w_refined_in[0]+ knot_step_w_refined_in *j);
        }
        knot_vector_u_refined_in.push_back(knot_vector_u_complete[knot_vector_u_complete.size()-1]) ;
        knot_vector_v_refined_in.push_back(knot_vector_v_complete[knot_vector_v_complete.size()-1]) ;
        knot_vector_w_refined_in.push_back(knot_vector_w_complete[knot_vector_w_complete.size()-1]) ;


        //------------ OUTER REFINED -------------------
        std::vector<double> knot_vector_u_refined_out(1);    knot_vector_u_refined_out[0] = knot_vector_u_complete[0];
        std::vector<double> knot_vector_v_refined_out(1);    knot_vector_v_refined_out[0] = knot_vector_v_complete[0];
        std::vector<double> knot_vector_w_refined_out(1);    knot_vector_w_refined_out[0] = knot_vector_w_complete[0];
        double knot_step_u_refined_out = (knot_vector_u_complete[knot_vector_u_complete.size()-1]-knot_vector_u_refined_out[0]) / (insert_nb_per_span_u_refined_out+1) ;
        double knot_step_v_refined_out = (knot_vector_v_complete[knot_vector_v_complete.size()-1]-knot_vector_v_refined_out[0]) / (insert_nb_per_span_v_refined_out+1) ;
        double knot_step_w_refined_out = (knot_vector_w_complete[knot_vector_w_complete.size()-1]-knot_vector_w_refined_out[0]) / (insert_nb_per_span_w_refined_out+1) ;
        
        meshSizes_uvw[0] = knot_step_u_refined_out; 
        meshSizes_uvw[1] = knot_step_v_refined_out;
        meshSizes_uvw[2] = knot_step_w_refined_out;
        surrogate_model_part_outer.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uvw);
        for (IndexType j = 1; j < insert_nb_per_span_u_refined_out + 1; ++j) {
            knot_vector_u_refined_out.push_back(knot_vector_u_refined_out[0]+ knot_step_u_refined_out * j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_v_refined_out + 1; ++j) {
            knot_vector_v_refined_out.push_back(knot_vector_v_refined_out[0]+ knot_step_v_refined_out *j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_w_refined_out + 1; ++j) {
            knot_vector_w_refined_out.push_back(knot_vector_w_refined_out[0]+ knot_step_w_refined_out *j);
        }
        knot_vector_u_refined_out.push_back(knot_vector_u_complete[knot_vector_u_complete.size()-1]) ;
        knot_vector_v_refined_out.push_back(knot_vector_v_complete[knot_vector_v_complete.size()-1]) ;
        knot_vector_w_refined_out.push_back(knot_vector_w_complete[knot_vector_w_complete.size()-1]) ;


        Vector parameterExternalCoordinates(6);
        parameterExternalCoordinates[0] = knot_vector_u_refined_out[0];
        parameterExternalCoordinates[1] = knot_vector_v_refined_out[0];
        parameterExternalCoordinates[2] = knot_vector_w_refined_out[0];
        parameterExternalCoordinates[3] = knot_vector_u_refined_out[knot_vector_u_refined_out.size()-1];
        parameterExternalCoordinates[4] = knot_vector_v_refined_out[knot_vector_v_refined_out.size()-1];
        parameterExternalCoordinates[5] = knot_vector_w_refined_out[knot_vector_w_refined_out.size()-1];
        surrogate_model_part_inner.GetProcessInfo().SetValue(LOAD_MESHES, parameterExternalCoordinates);
        surrogate_model_part_outer.GetProcessInfo().SetValue(LOAD_MESHES, parameterExternalCoordinates);


        //------------------------ INNER -----------------------------
        PointVector points_in;
        for (auto &i_cond : initial_skin_model_part_in.Conditions()) {
            // points_in.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry()[0].X(), i_cond.GetGeometry()[0].Y(), i_cond.GetGeometry()[0].Z())));
            points_in.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        }
        DynamicBins testBins_in(points_in.begin(), points_in.end());

        // Create the matrix of active/inactive knot spans, one for inner and one for outer loop
        const int numberOfInnerLoops = mParameters["sbm_parameters"]["number_of_inner_loops"].GetInt();

        std::vector<std::vector<std::vector<std::vector<int>>>> knot_spans_available_inner;
        knot_spans_available_inner.reserve(numberOfInnerLoops);
        for (int i = 0; i < numberOfInnerLoops; ++i) {
            std::vector<std::vector<std::vector<int>>> matrix3D; 
            matrix3D.reserve(insert_nb_per_span_u_refined_in+1);
            for (int j = 0; j <= insert_nb_per_span_u_refined_in; ++j) {
                std::vector<std::vector<int>> matrix;
                matrix.reserve(insert_nb_per_span_v_refined_in+1); 
                for (int k = 0; k <= insert_nb_per_span_v_refined_in; ++k) {
                    std::vector<int> row(insert_nb_per_span_w_refined_in+1);
                    matrix.push_back(row);
                }
                matrix3D.push_back(matrix); 
            }
            knot_spans_available_inner.push_back(matrix3D); 
        }
        
        // Optimized Snake -> for inner loops
        int idMatrixKnotSpansAvailable = 0;
        int idFirstNode;
        bool newInnerLoop = true;
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner :: Starting SnakeStep" << std::endl;

        if (initial_skin_model_part_in.Conditions().size() > 0) {
            
            // Copy all the nodes of the initial_skin_model_part to the skin_model_part
            for (auto &i_node : initial_skin_model_part_in.Nodes()) {
                skin_model_part_in.CreateNewNode(i_node.Id(), i_node.X(), i_node.Y(), i_node.Z());
            } 
            
            for (auto &i_cond : initial_skin_model_part_in.Conditions()) {  
                // if (newInnerLoop) {
                //     idFirstNode = i_cond.GetGeometry()[0].Id();
                //     newInnerLoop = false;
                // }
                double x_true_boundary1 = i_cond.GetGeometry()[0].X();
                double y_true_boundary1 = i_cond.GetGeometry()[0].Y();
                double z_true_boundary1 = i_cond.GetGeometry()[0].Z();
                double x_true_boundary2 = i_cond.GetGeometry()[1].X();
                double y_true_boundary2 = i_cond.GetGeometry()[1].Y();
                double z_true_boundary2 = i_cond.GetGeometry()[1].Z();
                double x_true_boundary3 = i_cond.GetGeometry()[2].X();
                double y_true_boundary3 = i_cond.GetGeometry()[2].Y();
                double z_true_boundary3 = i_cond.GetGeometry()[2].Z();
                
                // Find the intersections of the skin boundary with the knot values
                int knot_span_u_1st_point = x_true_boundary1 / knot_step_u_refined_in ;
                int knot_span_u_2nd_point = x_true_boundary2 / knot_step_u_refined_in ;
                int knot_span_u_3rd_point = x_true_boundary3 / knot_step_u_refined_in ;

                int knot_span_v_1st_point = y_true_boundary1 / knot_step_v_refined_in ;
                int knot_span_v_2nd_point = y_true_boundary2 / knot_step_v_refined_in ;
                int knot_span_v_3rd_point = y_true_boundary3 / knot_step_v_refined_in ;

                int knot_span_w_1st_point = z_true_boundary1 / knot_step_w_refined_in ;
                int knot_span_w_2nd_point = z_true_boundary2 / knot_step_w_refined_in ;
                int knot_span_w_3rd_point = z_true_boundary3 / knot_step_w_refined_in ;

                array_1d<IndexType, 3> ordered_ids;
                ordered_ids[0] = i_cond.GetGeometry()[0].Id();
                ordered_ids[1] = i_cond.GetGeometry()[1].Id();
                ordered_ids[2] = i_cond.GetGeometry()[2].Id();

                SnakeStep3D(skin_model_part_in, knot_spans_available_inner, idMatrixKnotSpansAvailable,
                            knot_span_u_1st_point, knot_span_u_2nd_point, knot_span_u_3rd_point,
                            knot_span_v_1st_point, knot_span_v_2nd_point, knot_span_v_3rd_point,
                            knot_span_w_1st_point, knot_span_w_2nd_point, knot_span_w_3rd_point,
                            x_true_boundary1, x_true_boundary2, x_true_boundary3, 
                            y_true_boundary1, y_true_boundary2, y_true_boundary3, 
                            z_true_boundary1, z_true_boundary2, z_true_boundary3, 
                            knot_step_u_refined_in, knot_step_v_refined_in, knot_step_w_refined_in,
                            ordered_ids);
                
                // if (i_cond.GetGeometry()[1].Id() == idFirstNode) {
                //     idMatrixKnotSpansAvailable++;
                //     newInnerLoop = true;
                // }
            }
        }

        // KRATOS_WATCH(knot_spans_available_inner)
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner :: Ending SnakeStep" << std::endl;
        // Read lambda parameters: 0.0 -> External,  0.5 -> Optimal
        const double lambda_inner = mParameters["sbm_parameters"]["lambda_inner"].GetDouble();
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner-> MarkKnotSpansAvailable" << std::endl;
        for (int i = 0; i < numberOfInnerLoops; i++) {
            int idInnerLoop = i;
            // Mark the knot_spans_available's for inner and outer loops
            MarkKnotSpansAvailable3D(knot_spans_available_inner, idInnerLoop, testBins_in, initial_skin_model_part_in, lambda_inner, 
                                    insert_nb_per_span_u_refined_in, insert_nb_per_span_v_refined_in, insert_nb_per_span_w_refined_in, 
                                    knot_step_u_refined_in, knot_step_v_refined_in, knot_step_w_refined_in);
            // KRATOS_WATCH(knot_spans_available_inner)

            CreateSurrogateBuondaryFromSnake_inner3D(knot_spans_available_inner,idInnerLoop, 
                                                    surrogate_model_part_inner, initial_skin_model_part_in, testBins_in,
                                                    insert_nb_per_span_u_refined_in, insert_nb_per_span_v_refined_in, insert_nb_per_span_w_refined_in,
                                                    knot_vector_u_refined_in, knot_vector_v_refined_in, knot_vector_w_refined_in,
                                                    knot_step_u_refined_in, knot_step_v_refined_in, knot_step_w_refined_in);
            KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner-> Loop finished" << std::endl;
        }
        

        // TO DO: Cannot do this anymore -> in the surrogate model part all possible nodes have been added

        double box_refinement_u_min = 1e16;
        double box_refinement_u_max = -1e16;
        double box_refinement_v_min = 1e16;
        double box_refinement_v_max = -1e16;
        double box_refinement_w_min = 1e16;
        double box_refinement_w_max = -1e16;

        for(ModelPart::ConditionsContainerType::iterator i_cond =  surrogate_model_part_inner.ConditionsBegin() ; i_cond != surrogate_model_part_inner.ConditionsEnd() ; i_cond++)
	    {
          if (i_cond->GetGeometry()[0].X() < box_refinement_u_min ){box_refinement_u_min = i_cond->GetGeometry()[0].X();}
          if (i_cond->GetGeometry()[2].X() > box_refinement_u_max ){box_refinement_u_max = i_cond->GetGeometry()[2].X();}

          if (i_cond->GetGeometry()[0].Y() < box_refinement_v_min ){box_refinement_v_min = i_cond->GetGeometry()[0].Y();}
          if (i_cond->GetGeometry()[2].Y() > box_refinement_v_max ){box_refinement_v_max = i_cond->GetGeometry()[2].Y();}

          if (i_cond->GetGeometry()[0].Z() < box_refinement_w_min ){box_refinement_w_min = i_cond->GetGeometry()[0].Z();}
          if (i_cond->GetGeometry()[2].Z() > box_refinement_w_max ){box_refinement_w_max = i_cond->GetGeometry()[2].Z();}
          
	    }
        Vector box_refinement_u_v_w_min_max(6);
        box_refinement_u_v_w_min_max[0] = box_refinement_u_min; 
        box_refinement_u_v_w_min_max[1] = box_refinement_u_max;
        box_refinement_u_v_w_min_max[2] = box_refinement_v_min; 
        box_refinement_u_v_w_min_max[3] = box_refinement_v_max;
        box_refinement_u_v_w_min_max[4] = box_refinement_w_min; 
        box_refinement_u_v_w_min_max[5] = box_refinement_w_max;
        iga_model_part.GetProcessInfo().SetValue(MARKER_MESHES, box_refinement_u_v_w_min_max);
        
    
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner-> Last loop has finished" << std::endl;

        // //--------------------------------------------------------------------------------
        // //---------------------------------- OUTER ---------------------------------------
        // //--------------------------------------------------------------------------------
        PointVector points_out;
        for (auto &i_cond : initial_skin_model_part_out.Conditions()) {
            points_out.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        }

        DynamicBins testBins_out(points_out.begin(), points_out.end());
        const int numberOfOuterLoops = 1;

        std::vector<std::vector<std::vector<std::vector<int>>>> knot_spans_available_outer;
        knot_spans_available_outer.reserve(numberOfOuterLoops);
        for (int i = 0; i < numberOfOuterLoops; ++i) {
            std::vector<std::vector<std::vector<int>>> matrix3D; 
            matrix3D.reserve(insert_nb_per_span_u_refined_out+1);
            for (int j = 0; j <= insert_nb_per_span_u_refined_out; ++j) {
                std::vector<std::vector<int>> matrix;
                matrix.reserve(insert_nb_per_span_v_refined_out+1); 
                for (int k = 0; k <= insert_nb_per_span_v_refined_out; ++k) {
                    std::vector<int> row(insert_nb_per_span_w_refined_out+1);
                    matrix.push_back(row);
                }
                matrix3D.push_back(matrix); 
            }
            knot_spans_available_outer.push_back(matrix3D); 
        }

        // Optimized Snake -> for outer loops
        idMatrixKnotSpansAvailable = 0;
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Outer :: Starting SnakeStep" << std::endl;

        if (initial_skin_model_part_out.Conditions().size() > 0) {
            
            // Copy all the nodes of the initial_skin_model_part to the skin_model_part
            for (auto &i_node : initial_skin_model_part_out.Nodes()) {
                skin_model_part_out.CreateNewNode(i_node.Id(), i_node.X(), i_node.Y(), i_node.Z());
            } 

            for (auto &i_cond : initial_skin_model_part_out.Conditions()) {  
                // if (newInnerLoop) {
                //     idFirstNode = i_cond.GetGeometry()[0].Id();
                //     newInnerLoop = false;
                // }
                double x_true_boundary1 = i_cond.GetGeometry()[0].X();
                double y_true_boundary1 = i_cond.GetGeometry()[0].Y();
                double z_true_boundary1 = i_cond.GetGeometry()[0].Z();
                double x_true_boundary2 = i_cond.GetGeometry()[1].X();
                double y_true_boundary2 = i_cond.GetGeometry()[1].Y();
                double z_true_boundary2 = i_cond.GetGeometry()[1].Z();
                double x_true_boundary3 = i_cond.GetGeometry()[2].X();
                double y_true_boundary3 = i_cond.GetGeometry()[2].Y();
                double z_true_boundary3 = i_cond.GetGeometry()[2].Z();
                
                // Find the intersections of the skin boundary with the knot values
                int knot_span_u_1st_point = x_true_boundary1 / knot_step_u_refined_out ;
                int knot_span_u_2nd_point = x_true_boundary2 / knot_step_u_refined_out ;
                int knot_span_u_3rd_point = x_true_boundary3 / knot_step_u_refined_out ;

                int knot_span_v_1st_point = y_true_boundary1 / knot_step_v_refined_out ;
                int knot_span_v_2nd_point = y_true_boundary2 / knot_step_v_refined_out ;
                int knot_span_v_3rd_point = y_true_boundary3 / knot_step_v_refined_out ;

                int knot_span_w_1st_point = z_true_boundary1 / knot_step_w_refined_out ;
                int knot_span_w_2nd_point = z_true_boundary2 / knot_step_w_refined_out ;
                int knot_span_w_3rd_point = z_true_boundary3 / knot_step_w_refined_out ;

                array_1d<IndexType, 3> ordered_ids;
                ordered_ids[0] = i_cond.GetGeometry()[0].Id();
                ordered_ids[1] = i_cond.GetGeometry()[1].Id();
                ordered_ids[2] = i_cond.GetGeometry()[2].Id();

                SnakeStep3D(skin_model_part_out, knot_spans_available_outer, idMatrixKnotSpansAvailable,
                            knot_span_u_1st_point, knot_span_u_2nd_point, knot_span_u_3rd_point,
                            knot_span_v_1st_point, knot_span_v_2nd_point, knot_span_v_3rd_point,
                            knot_span_w_1st_point, knot_span_w_2nd_point, knot_span_w_3rd_point,
                            x_true_boundary1, x_true_boundary2, x_true_boundary3, 
                            y_true_boundary1, y_true_boundary2, y_true_boundary3, 
                            z_true_boundary1, z_true_boundary2, z_true_boundary3, 
                            knot_step_u_refined_out, knot_step_v_refined_out, knot_step_w_refined_out, ordered_ids);
                
                // if (i_cond.GetGeometry()[1].Id() == idFirstNode) {
                //     idMatrixKnotSpansAvailable++;
                //     newInnerLoop = true;
                // }
            }

            // KRATOS_WATCH(knot_spans_available_outer)
            KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Outer :: Ending SnakeStep" << std::endl;
            // Read lambda parameters: 0.0 -> External,  0.5 -> Optimal
            const double lambda_outer = mParameters["sbm_parameters"]["lambda_outer"].GetDouble();
            KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Outer-> MarkKnotSpansAvailable" << std::endl;
            for (int i = 0; i < numberOfOuterLoops; i++) {
                int idOuterLoop = i;
                // Mark the knot_spans_available's for inner and outer loops
                MarkKnotSpansAvailable3D(knot_spans_available_outer, idOuterLoop, testBins_out, initial_skin_model_part_out, lambda_outer, 
                                        insert_nb_per_span_u_refined_out, insert_nb_per_span_v_refined_out, insert_nb_per_span_w_refined_out, 
                                        knot_step_u_refined_out, knot_step_v_refined_out, knot_step_w_refined_out);

                CreateSurrogateBuondaryFromSnake_outer3D(testBins_out, initial_skin_model_part_out, knot_spans_available_outer, idMatrixKnotSpansAvailable,
                                                        surrogate_model_part_outer,
                                                        insert_nb_per_span_u_refined_out, insert_nb_per_span_v_refined_out, insert_nb_per_span_w_refined_out,
                                                        knot_vector_u_refined_out, knot_vector_v_refined_out, knot_vector_w_refined_out);

                KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "OUTER-> Loop finished" << std::endl;
            }
        }
    }   






    
    void SnakeSBMUtilities::SnakeStep3D(ModelPart& skin_model_part, std::vector<std::vector<std::vector<std::vector<int>>>> &knot_spans_available, int idMatrix, 
                int knot_span_u_1st_point, int knot_span_u_2nd_point, int knot_span_u_3rd_point,
                int knot_span_v_1st_point, int knot_span_v_2nd_point, int knot_span_v_3rd_point,
                int knot_span_w_1st_point, int knot_span_w_2nd_point, int knot_span_w_3rd_point,
                double& x_true_boundary1, double& x_true_boundary2, double& x_true_boundary3,
                double& y_true_boundary1, double& y_true_boundary2, double& y_true_boundary3,
                double& z_true_boundary1, double& z_true_boundary2, double& z_true_boundary3,
                double& knot_step_u,      double& knot_step_v,      double& knot_step_w, 
                array_1d<IndexType, 3>& ordered_ids) {

        bool isSplitted = false;
        if (knot_span_u_1st_point != knot_span_u_2nd_point || knot_span_u_1st_point != knot_span_u_3rd_point || knot_span_u_2nd_point != knot_span_u_3rd_point || 
            knot_span_v_1st_point != knot_span_v_2nd_point || knot_span_v_1st_point != knot_span_v_3rd_point || knot_span_v_2nd_point != knot_span_v_3rd_point || 
            knot_span_w_1st_point != knot_span_w_2nd_point || knot_span_w_1st_point != knot_span_w_3rd_point || knot_span_w_2nd_point != knot_span_w_3rd_point ) 
            { 
            // Check if we are jumping some cut knot spans. If yes we split the true segment
            if (std::abs(knot_span_u_1st_point-knot_span_u_2nd_point) > 1 || std::abs(knot_span_u_1st_point-knot_span_u_3rd_point) > 1 || std::abs(knot_span_u_2nd_point-knot_span_u_3rd_point) > 1 || 
                std::abs(knot_span_v_1st_point-knot_span_v_2nd_point) > 1 || std::abs(knot_span_v_1st_point-knot_span_v_3rd_point) > 1 || std::abs(knot_span_v_2nd_point-knot_span_v_3rd_point) > 1 || 
                std::abs(knot_span_w_1st_point-knot_span_w_2nd_point) > 1 || std::abs(knot_span_w_1st_point-knot_span_w_3rd_point) > 1 || std::abs(knot_span_w_2nd_point-knot_span_w_3rd_point) > 1) {
                //  (knot_span_u_1st_point != knot_span_u_2nd_point && knot_span_v_1st_point != knot_span_v_2nd_point) 
                isSplitted = true;
                KRATOS_INFO("::[SnakeSBMUtilities]::") << "SnakeStep :: Spitting a 3D condition" << std::endl;
                // Compute in which knot spans they lie
                // Points between node 1 & 2
                double x_true_boundary_split12 = (x_true_boundary1+x_true_boundary2 ) / 2;
                double y_true_boundary_split12 = (y_true_boundary1+y_true_boundary2 ) / 2;
                double z_true_boundary_split12 = (z_true_boundary1+z_true_boundary2 ) / 2;
                int knot_span_u_point_split12 = x_true_boundary_split12 / knot_step_u ;
                int knot_span_v_point_split12 = y_true_boundary_split12 / knot_step_v ;
                int knot_span_w_point_split12 = z_true_boundary_split12 / knot_step_w ;
                // Points between node 2 & 3
                double x_true_boundary_split23 = (x_true_boundary2+x_true_boundary3 ) / 2;
                double y_true_boundary_split23 = (y_true_boundary2+y_true_boundary3 ) / 2;
                double z_true_boundary_split23 = (z_true_boundary2+z_true_boundary3 ) / 2;
                int knot_span_u_point_split23 = x_true_boundary_split23 / knot_step_u ;
                int knot_span_v_point_split23 = y_true_boundary_split23 / knot_step_v ;
                int knot_span_w_point_split23 = z_true_boundary_split23 / knot_step_w ;
                // Points between node 3 & 1
                double x_true_boundary_split31 = (x_true_boundary3+x_true_boundary1 ) / 2;
                double y_true_boundary_split31 = (y_true_boundary3+y_true_boundary1 ) / 2;
                double z_true_boundary_split31 = (z_true_boundary3+z_true_boundary1 ) / 2;
                int knot_span_u_point_split31 = x_true_boundary_split31 / knot_step_u ;
                int knot_span_v_point_split31 = y_true_boundary_split31 / knot_step_v ;
                int knot_span_w_point_split31 = z_true_boundary_split31 / knot_step_w ;

                // CREATE three NEW NODES in the middle of each side
                auto idNode_12 = skin_model_part.Nodes().size()+1;
                auto idNode_23 = idNode_12+1;
                auto idNode_31 = idNode_23+1;
                skin_model_part.CreateNewNode(idNode_12, (x_true_boundary1+x_true_boundary2 ) / 2, (y_true_boundary1+y_true_boundary2 ) / 2, (z_true_boundary1+z_true_boundary2 ) / 2);
                skin_model_part.CreateNewNode(idNode_23, (x_true_boundary3+x_true_boundary2 ) / 2, (y_true_boundary3+y_true_boundary2 ) / 2, (z_true_boundary3+z_true_boundary2 ) / 2);
                skin_model_part.CreateNewNode(idNode_31, (x_true_boundary1+x_true_boundary3 ) / 2, (y_true_boundary1+y_true_boundary3 ) / 2, (z_true_boundary1+z_true_boundary3 ) / 2);
                
                array_1d<IndexType, 3> ordered_ids_1;
                array_1d<IndexType, 3> ordered_ids_2;
                array_1d<IndexType, 3> ordered_ids_3;
                array_1d<IndexType, 3> ordered_ids_4;
                // We keep the same order so that the direction of the normal remains the same (inside/outiside keeps working)
                ordered_ids_1[0] = ordered_ids[0];     ordered_ids_1[1] = idNode_12;       ordered_ids_1[2] = idNode_31;
                ordered_ids_2[0] = ordered_ids[1];     ordered_ids_2[1] = idNode_23;       ordered_ids_2[2] = idNode_12;
                ordered_ids_3[0] = ordered_ids[2];     ordered_ids_3[1] = idNode_31;       ordered_ids_3[2] = idNode_23;
                ordered_ids_4[0] = idNode_12;          ordered_ids_4[1] = idNode_23;       ordered_ids_4[2] = idNode_31;
                

                // We do it recursively for all the new 4 tringular conditions
                SnakeStep3D(skin_model_part, knot_spans_available, idMatrix, knot_span_u_1st_point,   knot_span_u_point_split12, knot_span_u_point_split31,
                                                                           knot_span_v_1st_point,   knot_span_v_point_split12, knot_span_v_point_split31,
                                                                           knot_span_w_1st_point,   knot_span_w_point_split12, knot_span_w_point_split31,
                                                                           x_true_boundary1     ,   x_true_boundary_split12, x_true_boundary_split31, 
                                                                           y_true_boundary1     ,   y_true_boundary_split12, y_true_boundary_split31,
                                                                           z_true_boundary1     ,   z_true_boundary_split12, z_true_boundary_split31,
                                                                           knot_step_u, knot_step_v, knot_step_w, ordered_ids_1 );
                SnakeStep3D(skin_model_part, knot_spans_available, idMatrix, knot_span_u_2nd_point,   knot_span_u_point_split12, knot_span_u_point_split23,
                                                                           knot_span_v_2nd_point,   knot_span_v_point_split12, knot_span_v_point_split23,
                                                                           knot_span_w_2nd_point,   knot_span_w_point_split12, knot_span_w_point_split23,
                                                                           x_true_boundary2     ,   x_true_boundary_split12, x_true_boundary_split23, 
                                                                           y_true_boundary2     ,   y_true_boundary_split12, y_true_boundary_split23,
                                                                           z_true_boundary2     ,   z_true_boundary_split12, z_true_boundary_split23,
                                                                           knot_step_u, knot_step_v, knot_step_w, ordered_ids_2  );
                SnakeStep3D(skin_model_part, knot_spans_available, idMatrix, knot_span_u_3rd_point,   knot_span_u_point_split31, knot_span_u_point_split23,
                                                                           knot_span_v_3rd_point,   knot_span_v_point_split31, knot_span_v_point_split23,
                                                                           knot_span_w_3rd_point,   knot_span_w_point_split31, knot_span_w_point_split23,
                                                                           x_true_boundary3     ,   x_true_boundary_split31, x_true_boundary_split23, 
                                                                           y_true_boundary3     ,   y_true_boundary_split31, y_true_boundary_split23,
                                                                           z_true_boundary3     ,   z_true_boundary_split31, z_true_boundary_split23,
                                                                           knot_step_u, knot_step_v, knot_step_w, ordered_ids_3 );
                // Last tringle connecting all the split points
                // TO DO ... remove? 
                SnakeStep3D(skin_model_part, knot_spans_available, idMatrix, knot_span_u_point_split12, knot_span_u_point_split31, knot_span_u_point_split23,
                                                                           knot_span_v_point_split12, knot_span_v_point_split31, knot_span_v_point_split23,
                                                                           knot_span_w_point_split12, knot_span_w_point_split31, knot_span_w_point_split23,
                                                                           x_true_boundary_split12,   x_true_boundary_split31, x_true_boundary_split23, 
                                                                           y_true_boundary_split12,   y_true_boundary_split31, y_true_boundary_split23,
                                                                           z_true_boundary_split12,   z_true_boundary_split31, z_true_boundary_split23,
                                                                           knot_step_u, knot_step_v, knot_step_w, ordered_ids_4  );
            }
            /* The condition it is NOT split
               Check if the true boundary crosses an u, v, or w knot value */
            // Points between node 1 & 2 
            else if (knot_span_u_1st_point != knot_span_u_2nd_point) { // u knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_2nd_point] = 2;
            }
            else if (knot_span_v_1st_point != knot_span_v_2nd_point) { // v knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_2nd_point][knot_span_u_1st_point] = 2;
            }
            else if (knot_span_w_1st_point != knot_span_w_2nd_point) { // w knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_2nd_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
            }
            // Points between node 1 & 3
            else if (knot_span_u_1st_point != knot_span_u_3rd_point) { // u knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_3rd_point] = 2;
            }
            else if (knot_span_v_1st_point != knot_span_v_3rd_point) { // v knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_3rd_point][knot_span_u_1st_point] = 2;
            }
            else if (knot_span_w_1st_point != knot_span_w_3rd_point) { // w knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_1st_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_3rd_point][knot_span_v_1st_point][knot_span_u_1st_point] = 2;
            }
            // Points between node 3 & 2 // it never enters here?? Might be useless?
            else if (knot_span_u_3rd_point != knot_span_u_2nd_point) { // u knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_3rd_point][knot_span_v_3rd_point][knot_span_u_3rd_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_3rd_point][knot_span_v_3rd_point][knot_span_u_2nd_point] = 2;
            }
            else if (knot_span_v_3rd_point != knot_span_v_2nd_point) { // v knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_3rd_point][knot_span_v_3rd_point][knot_span_u_3rd_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_3rd_point][knot_span_v_2nd_point][knot_span_u_3rd_point] = 2;
            }
            else if (knot_span_w_3rd_point != knot_span_w_2nd_point) { // w knot value is crossed
                knot_spans_available[idMatrix][knot_span_w_3rd_point][knot_span_v_3rd_point][knot_span_u_3rd_point] = 2;
                knot_spans_available[idMatrix][knot_span_w_2nd_point][knot_span_v_3rd_point][knot_span_u_3rd_point] = 2;
            }
            else {KRATOS_WATCH("Something went wrong"); exit(0);}
        }
        if (!isSplitted) {
            // TO DO -> do you want to create a new skin_model_part with the split triangles?
            //           How do you manage the creation of the nodes? They might already exist in the model part
            //           Maybe is better to add nodes and condition to the initial_skin_model_part 

            //// Create 1 new conditions for each skin condition 
            Properties::Pointer p_cond_prop = skin_model_part.pGetProperties(0);
            auto lastIdCondition = skin_model_part.Conditions().size();

            Condition::Pointer p_cond1 = skin_model_part.CreateNewCondition("SurfaceCondition3D3N", lastIdCondition+1, 
                                                                            {{ordered_ids[0], ordered_ids[1], ordered_ids[2]}}, p_cond_prop );
        }
    }


    bool SnakeSBMUtilities::isPointInsideSkinBoundary3D(Point& point0, DynamicBins& testBins, ModelPart& skin_model_part)
    {
        // Get the nearest point of the true boundary
        PointerType pointToSearch = PointerType(new PointType(10000, point0.X(), point0.Y(), point0.Z()));
        PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);

        // Get the closest Condition the initial_skin_model_part_in.Conditions
        int id1 = nearestPoint->Id();
        auto nearestCondition = skin_model_part.GetCondition(id1);

        // Point0 -> pointToSearch
        Point point1 = nearestCondition.GetGeometry()[0]; // FIRST POINT IN TRUE GEOM
        Point point2 = nearestCondition.GetGeometry()[1]; // SECOND POINT IN TRUE GEOM
        Point point3 = nearestCondition.GetGeometry()[2]; // THIRD POINT IN TRUE GEOM

        array_1d<double,3> v_1 = point2 - point1;
        array_1d<double,3> v_2 = point3 - point1;
        array_1d<double,3> normal;
        MathUtils<double>::CrossProduct(normal, v_1, v_2);
        
        normal = normal/norm_2(normal);
        
        array_1d<double,3> CenterToPointToSearch = point0 - nearestCondition.GetGeometry().Center() ;

        CenterToPointToSearch = CenterToPointToSearch/norm_2(CenterToPointToSearch);
        // Scalar product between the nornal and the CenterToPointToSearch
        bool isInside = false;
        if (inner_prod(normal, CenterToPointToSearch) > 0) {isInside = true;}

        return isInside;
    }



    void SnakeSBMUtilities::MarkKnotSpansAvailable3D(std::vector<std::vector<std::vector<std::vector<int>>>> & knot_spans_available, int idMatrix,DynamicBins& testBin, ModelPart& skin_model_part, double lambda, 
                                                        int insert_nb_per_span_u, int insert_nb_per_span_v, int insert_nb_per_span_w, 
                                                        double knot_step_u, double knot_step_v, double knot_step_w) {
        
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                    if (knot_spans_available[idMatrix][k][j][i] == 2) {
                        // Check the 26(+1) neighbor knot spans -> Is there any completely inside? Note that we can just check 1 point.
                        std::vector<int> move_vector = {-1, 0, 1};
                        // Need all the combination in the three directions
                        for (int i_move : move_vector){
                            if (i != 0 && i != insert_nb_per_span_u) {
                                //---------
                                for (int j_move : move_vector){
                                    if (j != 0 && j != insert_nb_per_span_v) {
                                        //----------
                                        for (int k_move : move_vector){
                                            if (k != 0 && k != insert_nb_per_span_w) {
                                                //---------
                                                if (knot_spans_available[idMatrix][k+k_move][j+j_move][i + i_move] == 0) {
                                                    Point gaussPoint = Point((i + i_move + 0.5)*knot_step_u, (j+ j_move + 0.5)*knot_step_v, (k + k_move + 0.5)*knot_step_w);
                                                    if (isPointInsideSkinBoundary3D(gaussPoint, testBin, skin_model_part)) {
                                                            knot_spans_available[idMatrix][k+k_move][j+j_move][i+i_move] = 1;
                                                    }
                                                }
                                                //--------
                                            }
                                        }
                                        //----------
                                    }
                                }
                                //--------
                            }
                        }
                        
                        // Create "fake" GaussPoints to check if the majority are inside or outside
                        const int numFakeGaussPoints = 4;
                        int numberOfInsideGaussianPoints = 0;
                        for (int i_GPx = 0; i_GPx < numFakeGaussPoints; i_GPx++){
                            double x_coord = i*knot_step_u + knot_step_u/(numFakeGaussPoints+1)*(i_GPx+1);

                            // NOTE:: The v-knot spans are upside down in the matrix!!
                            for (int i_GPy = 0; i_GPy < numFakeGaussPoints; i_GPy++) 
                            {   
                                double y_coord = j*knot_step_v + knot_step_v/(numFakeGaussPoints+1)*(i_GPy+1);
                                for (int i_GPz = 0; i_GPz < numFakeGaussPoints; i_GPz++) 
                                {
                                    double z_coord = k*knot_step_w + knot_step_w/(numFakeGaussPoints+1)*(i_GPz+1);
                                    Point gaussPoint = Point(x_coord, y_coord, z_coord);  // GAUSSIAN POINT
                                    
                                    if (isPointInsideSkinBoundary3D(gaussPoint, testBin, skin_model_part)) {
                                        // Sum over the number of numFakeGaussPoints per knot span
                                        numberOfInsideGaussianPoints++;
                                    }
                                }
                            }
                            
                        }
                        // Mark the knot span as available or not depending on the number of Gauss Points Inside/Outside
                        if (numberOfInsideGaussianPoints < lambda*numFakeGaussPoints*numFakeGaussPoints*numFakeGaussPoints) {
                            knot_spans_available[idMatrix][k][j][i] = -1; // Cut knot spans that have been checked
                        }
                        else{
                            knot_spans_available[idMatrix][k][j][i] = 1; // The knot span is considered DEACTIVE
                        }
                                
                    }
                }
            }
        }
    }


    void SnakeSBMUtilities::CreateSurrogateBuondaryFromSnake_inner3D (std::vector<std::vector<std::vector<std::vector<int>>>> & knot_spans_available, int idMatrix, 
                         ModelPart& surrogate_model_part_inner, ModelPart& skin_model_part_inner, DynamicBins& testBins,
                         int insert_nb_per_span_u, int insert_nb_per_span_v,  int insert_nb_per_span_w,
                         std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v, std::vector<double>&  knot_vector_w,
                         double knot_step_u, double knot_step_v, double knot_step_w){
        
        std::ofstream outputFile("txt_files/Snake_coordinates.txt", std::ios::app);
        outputFile << std::setprecision(16);
        // Snake 3D works with a raycasting technique from each of the three directions
        int idSurrogateNode = 1;
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                    surrogate_model_part_inner.CreateNewNode(idSurrogateNode, knot_vector_u[i], knot_vector_v[j], knot_vector_w[k]);
                    idSurrogateNode++;
                }
            }
        }
        
        // Direction parallel to x
        bool knotSpansIsActive = true;
        
        uint idSurrogateCondition = 1;
        Properties::Pointer p_cond_prop = surrogate_model_part_inner.pGetProperties(0);
        
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                
                bool checkNextPoint = false;
                /*  
                    Formula to connect i,j,k to the id of the model_part
                    i + j*(insert_nb_per_span_u+1) + k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1) + 1;
                */
                // move in the x direction
                for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                    if (checkNextPoint) {
                        // Check i+1 point using isPointInsideSkinBoundary3D
                        Point centerPoint = Point((i + 0.5)*knot_step_u, (j + 0.5)*knot_step_v, (k + 0.5)*knot_step_w);
                        bool isExiting = false;
                        if ( knot_spans_available[idMatrix][k][j][i] == 1 ) {
                            // the knot span was already been checked very well
                        }
                        else if (isPointInsideSkinBoundary3D(centerPoint, testBins, skin_model_part_inner)) {
                            // STILL INSIDE --> do not save nothing and update knot_spans_available 
                            if ( knot_spans_available[idMatrix][k][j][i] == -1) {
                                isExiting = true;
                            }
                            else {
                                knot_spans_available[idMatrix][k][j][i] = 1;
                            }
                        }
                        else {
                            isExiting = true;
                        }
                        if (isExiting) {
                            /* EXITING --> save last face in direction x. i-th is the knot value. */
                            int node1_i = i; int node1_j = j;   int node1_k = k;
                            int node2_i = i; int node2_j = j+1; int node2_k = k;
                            int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                            int node4_i = i; int node4_j = j;   int node4_k = k+1;

                            IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                                
                            Condition::Pointer pcond = surrogate_model_part_inner.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                            
                            // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                            pcond->Set(BOUNDARY, false);

                            // surrogate_model_part_inner.AddCondition(pcond);
                            idSurrogateCondition++;
                            checkNextPoint = false;

                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k+1] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl;    
                        }
                        
                    }
                    else if (knot_spans_available[idMatrix][k][j][i] == 1) {
                        // ENTERING --> save first face in direction
                        int node1_i = i; int node1_j = j;   int node1_k = k;
                        int node2_i = i; int node2_j = j+1; int node2_k = k;
                        int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                        int node4_i = i; int node4_j = j;   int node4_k = k+1;

                        IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            
                        Condition::Pointer pcond = surrogate_model_part_inner.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                        // surrogate_model_part_inner.AddCondition(pcond);
                        idSurrogateCondition++;
                        checkNextPoint = true;

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, true);

                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k+1] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                    }
                }
            }
        }
        // Do the same for y and z direction, without isPointInsideSkinBoundary3D, since we have already done it
        // And it is not necessary do it again
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                
                bool checkNextPoint = false;
                /*  
                    Formula to connect i,j,k to the id of the model_part
                    i + j*(insert_nb_per_span_u+1) + k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                */
                // move in the x direction
                for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                    if (checkNextPoint) {
                        if (knot_spans_available[idMatrix][k][j][i] != 1) {
                            /* EXITING --> save last face in direction x. i-th is the knot value. */
                            int node1_i = i;   int node1_j = j;   int node1_k = k;
                            int node2_i = i+1; int node2_j = j;   int node2_k = k;
                            int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                            int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                            IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                                
                            Condition::Pointer pcond = surrogate_model_part_inner.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                            // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                            pcond->Set(BOUNDARY, false);

                            // surrogate_model_part_inner.AddCondition(p_cond);
                            idSurrogateCondition++;
                            checkNextPoint = false;

                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                        } 
                    }
                    else if (knot_spans_available[idMatrix][k][j][i] == 1) {
                        // ENTERING --> save first face in direction
                        int node1_i = i;   int node1_j = j;   int node1_k = k;
                        int node2_i = i+1; int node2_j = j;   int node2_k = k;
                        int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                        int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                        IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        Condition::Pointer pcond = surrogate_model_part_inner.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                        // surrogate_model_part_inner.AddCondition(p_cond);
                        idSurrogateCondition++;
                        checkNextPoint = true;

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, true);

                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                    }
                }
            }
        }
        // z direction
        for (int j = 0; j < insert_nb_per_span_v+1; j++) {
            for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                
                bool checkNextPoint = false;
                /*  
                    Formula to connect i,j,k to the id of the model_part
                    i + j*(insert_nb_per_span_u+1) + k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                */
                // move in the x direction
                for (int k = 0; k < insert_nb_per_span_w+1; k++) {
                    if (checkNextPoint) {
                        if (knot_spans_available[idMatrix][k][j][i] != 1) {
                            /* EXITING --> save last face in direction x. i-th is the knot value. */
                            int node1_i = i;   int node1_j = j;     int node1_k = k;
                            int node2_i = i+1; int node2_j = j;     int node2_k = k;
                            int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                            int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                            IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                                
                            Condition::Pointer pcond = surrogate_model_part_inner.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                            // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                            pcond->Set(BOUNDARY, false);
                            
                            // surrogate_model_part_inner.AddCondition(p_cond);
                            idSurrogateCondition++;
                            checkNextPoint = false;

                            outputFile << knot_vector_u[i]   << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i]   << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                        } 
                    }
                    else if (knot_spans_available[idMatrix][k][j][i] == 1) {
                        // ENTERING --> save first face in direction
                        int node1_i = i;   int node1_j = j;     int node1_k = k;
                        int node2_i = i+1; int node2_j = j;     int node2_k = k;
                        int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                        int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                        IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        Condition::Pointer pcond = surrogate_model_part_inner.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                        // surrogate_model_part_inner.AddCondition(p_cond);
                        idSurrogateCondition++;
                        checkNextPoint = true;

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, true);

                        outputFile << knot_vector_u[i]   << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i]   << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                    }
                }
            }
        }
        outputFile.close();
        KRATOS_WATCH("Inner Snake process has finished")


        // if (surrogate_model_part_inner.Elements().size()>0) {
        //     // Check if it is not the first inner loop
        //     initialId = surrogate_model_part_inner.GetElement(surrogate_model_part_inner.Elements().size()).GetGeometry()[1].Id()+1;
        // }
        // std::vector<ModelPart::IndexType> elem_nodes{initialId, idSurrogateCondition};
        // surrogate_model_part_inner.CreateNewElement("Element2D2N", surrogate_model_part_inner.Elements().size()+1, elem_nodes, p_cond_prop);

    }


    void SnakeSBMUtilities::CreateSurrogateBuondaryFromSnake_outer3D (DynamicBins& testBin_out, ModelPart& initial_skin_model_part_out, 
                                std::vector<std::vector<std::vector<std::vector<int>>>> & knot_spans_available, int idMatrix, ModelPart& surrogate_model_part_outer, 
                                int insert_nb_per_span_u, int insert_nb_per_span_v, int insert_nb_per_span_w, 
                                std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v, std::vector<double>&  knot_vector_w){

        KRATOS_INFO("::[SnakeSBMUtilities]::") << "Outer :: Check External layers in 3D" << std::endl;
        // CHECK ALL THE EXTERNAL KNOT SPANS -> Better also an additional layer
        
        double knot_step_u = knot_vector_u[1]-knot_vector_u[0];
        double knot_step_v = knot_vector_v[1]-knot_vector_v[0];
        double knot_step_w = knot_vector_w[1]-knot_vector_w[0];
        
        // LEFT BOUNDARY (u)
        for (int i = 0; i<2; i++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                for (int k = 0; k < (insert_nb_per_span_w+1); k++ ) {
                    Point centroidKnotSpan = Point((i+0.5)*knot_step_u, (j+0.5)*knot_step_v, (k+0.5)*knot_step_w);
                    if (isPointInsideSkinBoundary3D(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][k][j][i] != -1) {
                        knot_spans_available[idMatrix][k][j][i] = 1;
                    }
                }
            }
        }
        // RIGHT BOUNDARY (u)
        for (int i = knot_spans_available[idMatrix][0][0].size()-1; i > knot_spans_available[idMatrix][0][0].size()-3; i--) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ){
                for (int k = 0; k < (insert_nb_per_span_w+1); k++ ) {
                    Point centroidKnotSpan = Point((i+0.5)*knot_step_u, (j+0.5)*knot_step_v, (k+0.5)*knot_step_w);
                    if (isPointInsideSkinBoundary3D(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][k][j][i] != -1) {
                        knot_spans_available[idMatrix][k][j][i] = 1;
                    }
                }
            }
        }

        // FRONT BOUNDARY (v)
        for (int j = 0; j<2; j++) {
            for (int i = 0; i < (insert_nb_per_span_u+1); i++ ) {
                for (int k = 0; k < (insert_nb_per_span_w+1); k++ ) {
                    Point centroidKnotSpan = Point((i+0.5)*knot_step_u, (j+0.5)*knot_step_v, (k+0.5)*knot_step_w);
                    if (isPointInsideSkinBoundary3D(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][k][j][i] != -1) {
                        knot_spans_available[idMatrix][k][j][i] = 1;
                    }
                }
            }
        }
        // BACK BOUNDARY (v)
        for (int j = knot_spans_available[idMatrix][0].size()-1; j > knot_spans_available[idMatrix][0].size()-3; j--) {
            for (int i = 0; i < (insert_nb_per_span_u+1); i++ ) {
                for (int k = 0; k < (insert_nb_per_span_w+1); k++ ) {
                    Point centroidKnotSpan = Point((i+0.5)*knot_step_u, (j+0.5)*knot_step_v, (k+0.5)*knot_step_w);
                    if (isPointInsideSkinBoundary3D(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][k][j][i] != -1) {
                        knot_spans_available[idMatrix][k][j][i] = 1;
                    }
                }
            }
        }
        // BOTTOM BOUNDARY (w)
        for (int k = 0; k<2; k++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                for (int i = 0; i < (insert_nb_per_span_u+1); i++ ) {
                    Point centroidKnotSpan = Point((i+0.5)*knot_step_u, (j+0.5)*knot_step_v, (k+0.5)*knot_step_w);
                    if (isPointInsideSkinBoundary3D(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][k][j][i] != -1) {
                        knot_spans_available[idMatrix][k][j][i] = 1;
                    }
                }
            }
        }
        // TOP BOUNDARY (w)
        for (int k = knot_spans_available[idMatrix].size()-1; k > knot_spans_available[idMatrix].size()-3; k--) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                for (int i = 0; i < (insert_nb_per_span_u+1); i++ ) {
                    Point centroidKnotSpan = Point((i+0.5)*knot_step_u, (j+0.5)*knot_step_v, (k+0.5)*knot_step_w);
                    if (isPointInsideSkinBoundary3D(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][k][j][i] != -1) {
                        knot_spans_available[idMatrix][k][j][i] = 1;
                    }
                }
            }
        }

        KRATOS_INFO("::[SnakeSBMUtilities]::") << "Outer :: Starting Creation of Surrogate_Model_Part_Outer" << std::endl;
        
        // Create the surrogate boundary outer
        std::ofstream outputFile("txt_files/Snake_coordinates_outer.txt", std::ios::app);
        outputFile << std::setprecision(16);
        // Snake 3D works with a raycasting technique from each of the three directions
        int idSurrogateNode = 1;
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                    surrogate_model_part_outer.CreateNewNode(idSurrogateNode, knot_vector_u[i], knot_vector_v[j], knot_vector_w[k]);
                    idSurrogateNode++;
                }
            }
        }
        
        // Direction parallel to x
        bool knotSpansIsActive = true;
        
        uint idSurrogateCondition = 1;
        Properties::Pointer p_cond_prop = surrogate_model_part_outer.pGetProperties(0);
        
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                
                bool checkNextPoint = false;
                /*  
                    Formula to connect i,j,k to the id of the model_part
                    i + j*(insert_nb_per_span_u+1) + k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1) + 1;
                */
                // move in the x direction
                for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                    if (checkNextPoint) {
                        // Check i+1 point using isPointInsideSkinBoundary3D
                        Point centerPoint = Point((i + 0.5)*knot_step_u, (j + 0.5)*knot_step_v, (k + 0.5)*knot_step_w);
                        bool isExiting = false;
                        if ( knot_spans_available[idMatrix][k][j][i] == 1 ) {
                            // the knot span was already been checked very well
                        }
                        else if (isPointInsideSkinBoundary3D(centerPoint, testBin_out, initial_skin_model_part_out)) {
                            // STILL INSIDE --> do not save nothing and update knot_spans_available 
                            if ( knot_spans_available[idMatrix][k][j][i] == -1) {
                                isExiting = true;
                            }
                            else {
                                knot_spans_available[idMatrix][k][j][i] = 1;
                            }
                        }
                        else {
                            isExiting = true;
                        }
                        if (isExiting) {
                            /* EXITING --> save last face in direction x. i-th is the knot value. */
                            int node1_i = i; int node1_j = j;   int node1_k = k;
                            int node2_i = i; int node2_j = j+1; int node2_k = k;
                            int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                            int node4_i = i; int node4_j = j;   int node4_k = k+1;

                            IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                                
                            Condition::Pointer pcond = surrogate_model_part_outer.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                            
                            // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                            pcond->Set(BOUNDARY, false);

                            // surrogate_model_part_inner.AddCondition(pcond);
                            idSurrogateCondition++;
                            checkNextPoint = false;

                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k+1] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl;    
                        }
                        
                    }
                    else if (knot_spans_available[idMatrix][k][j][i] == 1) {
                        // ENTERING --> save first face in direction
                        int node1_i = i; int node1_j = j;   int node1_k = k;
                        int node2_i = i; int node2_j = j+1; int node2_k = k;
                        int node3_i = i; int node3_j = j+1; int node3_k = k+1; 
                        int node4_i = i; int node4_j = j;   int node4_k = k+1;

                        IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            
                        Condition::Pointer pcond = surrogate_model_part_outer.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                        // surrogate_model_part_inner.AddCondition(pcond);
                        idSurrogateCondition++;
                        checkNextPoint = true;

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, true);

                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k+1] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                    }
                }
            }
        }
        // Do the same for y and z direction, without isPointInsideSkinBoundary3D, since we have already done it
        // And it is not necessary do it again
        for (int k = 0; k < insert_nb_per_span_w+1; k++) {
            for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                
                bool checkNextPoint = false;
                /*  
                    Formula to connect i,j,k to the id of the model_part
                    i + j*(insert_nb_per_span_u+1) + k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                */
                // move in the x direction
                for (int j = 0; j < insert_nb_per_span_v+1; j++) {
                    if (checkNextPoint) {
                        if (knot_spans_available[idMatrix][k][j][i] != 1) {
                            /* EXITING --> save last face in direction x. i-th is the knot value. */
                            int node1_i = i;   int node1_j = j;   int node1_k = k;
                            int node2_i = i+1; int node2_j = j;   int node2_k = k;
                            int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                            int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                            IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                                
                            Condition::Pointer pcond = surrogate_model_part_outer.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                            // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                            pcond->Set(BOUNDARY, false);

                            // surrogate_model_part_inner.AddCondition(p_cond);
                            idSurrogateCondition++;
                            checkNextPoint = false;

                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                        } 
                    }
                    else if (knot_spans_available[idMatrix][k][j][i] == 1) {
                        // ENTERING --> save first face in direction
                        int node1_i = i;   int node1_j = j;   int node1_k = k;
                        int node2_i = i+1; int node2_j = j;   int node2_k = k;
                        int node3_i = i+1; int node3_j = j;   int node3_k = k+1; 
                        int node4_i = i;   int node4_j = j;   int node4_k = k+1;

                        IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        Condition::Pointer pcond = surrogate_model_part_outer.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                        // surrogate_model_part_inner.AddCondition(p_cond);
                        idSurrogateCondition++;
                        checkNextPoint = true;

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, true);

                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << " " << knot_vector_w[k+1] << std::endl; 
                    }
                }
            }
        }
        // z direction
        for (int j = 0; j < insert_nb_per_span_v+1; j++) {
            for (int i = 0; i < insert_nb_per_span_u+1; i++) {
                
                bool checkNextPoint = false;
                /*  
                    Formula to connect i,j,k to the id of the model_part
                    i + j*(insert_nb_per_span_u+1) + k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                */
                // move in the x direction
                for (int k = 0; k < insert_nb_per_span_w+1; k++) {
                    if (checkNextPoint) {
                        if (knot_spans_available[idMatrix][k][j][i] != 1) {
                            /* EXITING --> save last face in direction x. i-th is the knot value. */
                            int node1_i = i;   int node1_j = j;     int node1_k = k;
                            int node2_i = i+1; int node2_j = j;     int node2_k = k;
                            int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                            int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                            IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                            IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                                
                            Condition::Pointer pcond = surrogate_model_part_outer.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                            // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                            pcond->Set(BOUNDARY, false);
                            
                            // surrogate_model_part_inner.AddCondition(p_cond);
                            idSurrogateCondition++;
                            checkNextPoint = false;

                            outputFile << knot_vector_u[i]   << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                            outputFile << knot_vector_u[i]   << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                        } 
                    }
                    else if (knot_spans_available[idMatrix][k][j][i] == 1) {
                        // ENTERING --> save first face in direction
                        int node1_i = i;   int node1_j = j;     int node1_k = k;
                        int node2_i = i+1; int node2_j = j;     int node2_k = k;
                        int node3_i = i+1; int node3_j = j+1;   int node3_k = k; 
                        int node4_i = i;   int node4_j = j+1;   int node4_k = k;

                        IndexType id_node_1 = 1 + node1_i + node1_j*(insert_nb_per_span_u+1) + node1_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_2 = 1 + node2_i + node2_j*(insert_nb_per_span_u+1) + node2_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_3 = 1 + node3_i + node3_j*(insert_nb_per_span_u+1) + node3_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        IndexType id_node_4 = 1 + node4_i + node4_j*(insert_nb_per_span_u+1) + node4_k*(insert_nb_per_span_v+1)*(insert_nb_per_span_u+1);
                        Condition::Pointer pcond = surrogate_model_part_outer.CreateNewCondition("SurfaceCondition3D4N", idSurrogateCondition, {{id_node_1, id_node_2, id_node_3, id_node_4 }}, p_cond_prop );
                        // surrogate_model_part_inner.AddCondition(p_cond);
                        idSurrogateCondition++;
                        checkNextPoint = true;

                        // BOUNDARY true means that the condition (i.e. the sbm face) is entering looking from x,y,z positive
                        pcond->Set(BOUNDARY, true);

                        outputFile << knot_vector_u[i]   << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j]   << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i+1] << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                        outputFile << knot_vector_u[i]   << " " << knot_vector_v[j+1] << " " << knot_vector_w[k] << std::endl; 
                    }
                }
            }
        }
        outputFile.close();
        KRATOS_INFO("::[SnakeSBMUtilities]::") << "Outer :: Outer Snake process has finished" << std::endl;

    }

}  // namespace Kratos.
