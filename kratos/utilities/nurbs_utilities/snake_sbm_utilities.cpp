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

    void SnakeSBMUtilities::CreateTheSnakeCoordinates(ModelPart& iga_model_part,
                                                    ModelPart& skin_model_part_inner_initial,
                                                    ModelPart& skin_model_part_outer_initial, 
                                                    ModelPart& skin_model_part, 
                                                    int rEchoLevel, 
                                                    Vector& knot_vector_u, 
                                                    Vector& knot_vector_v,
                                                    const Parameters mParameters) { 
        
        // Check
        KRATOS_ERROR_IF_NOT(mParameters.Has("sbm_parameters")) << "sbm_parameters has not been defined in the nurbs modeler" << std::endl;
        
        // Initilize the property of skin_model_part_in and out
        if (skin_model_part_inner_initial.NumberOfNodes()>0 || skin_model_part_inner_initial.NumberOfGeometries()>0) {
        
            CreateTheSnakeCoordinates(iga_model_part, skin_model_part_inner_initial, skin_model_part, rEchoLevel, knot_vector_u, knot_vector_v, mParameters, true);
                
        }
        if (skin_model_part_outer_initial.NumberOfNodes()>0 || skin_model_part_outer_initial.NumberOfGeometries()>0) {

            CreateTheSnakeCoordinates(iga_model_part, skin_model_part_outer_initial, skin_model_part, rEchoLevel, knot_vector_u, knot_vector_v, mParameters, false);
        }
    }   

    void SnakeSBMUtilities::CreateTheSnakeCoordinates(ModelPart& iga_model_part, 
                                                      ModelPart& skin_model_part_initial, 
                                                      ModelPart& skin_model_part,
                                                      int rEchoLevel, 
                                                      Vector& knot_vector_u, 
                                                      Vector& knot_vector_v,
                                                      const Parameters mParameters,
                                                      bool is_inner) 
    { 
        
        std::string surrogate_sub_model_part_name; 
        std::string skin_sub_model_part_name; 
        // ModelPart skin_sub_model_part; 
        if (is_inner) {
            surrogate_sub_model_part_name = "surrogate_inner";
            skin_sub_model_part_name = "inner";
        }
        else {
            surrogate_sub_model_part_name = "surrogate_outer";
            skin_sub_model_part_name = "outer";
        }
        
        ModelPart& skin_sub_model_part = skin_model_part.GetSubModelPart(skin_sub_model_part_name);
        ModelPart& surrogate_sub_model_part = iga_model_part.GetSubModelPart(surrogate_sub_model_part_name);

        const bool is_initial_skin_nurbs = skin_model_part_initial.Geometries().size()>0;
        
        Properties::Pointer p_cond_prop_in = skin_model_part_initial.pGetProperties(0);
        skin_model_part_initial.AddProperties(p_cond_prop_in);

        Vector knot_step_uv(2);
        knot_step_uv[0] = std::abs(knot_vector_u[int(knot_vector_u.size()/2) +1]  - knot_vector_u[int(knot_vector_u.size()/2)] ) ;
        knot_step_uv[1] = std::abs(knot_vector_v[int(knot_vector_v.size()/2) +1]  - knot_vector_v[int(knot_vector_v.size()/2)] ) ;

        Vector meshSizes_uv(2);
        meshSizes_uv[0] = knot_step_uv[0]; 
        meshSizes_uv[1] = knot_step_uv[1];
        ModelPart& surrogate_model_part = iga_model_part.GetSubModelPart(surrogate_sub_model_part_name);
        surrogate_model_part.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uv);

        Vector starting_pos_uv(2);
        starting_pos_uv[0] = knot_vector_u[0];
        starting_pos_uv[1] = knot_vector_v[0];

        Vector parameterExternalCoordinates(4);
        parameterExternalCoordinates[0] = knot_vector_u[0];
        parameterExternalCoordinates[1] = knot_vector_v[0];
        parameterExternalCoordinates[2] = knot_vector_u[knot_vector_u.size()-1];
        parameterExternalCoordinates[3] = knot_vector_v[knot_vector_v.size()-1];
        
        surrogate_model_part.GetProcessInfo().SetValue(LOAD_MESHES, parameterExternalCoordinates);


        // Create the matrix of active/inactive knot spans, one for inner and one for outer loop
        unsigned int numberOfLoops;
        if (is_inner)
            numberOfLoops = mParameters["sbm_parameters"]["number_of_inner_loops"].GetInt();
        else 
            numberOfLoops = 1;
        
        std::vector<int> n_knot_spans_uv(2);
        n_knot_spans_uv[0] = knot_vector_u.size()-1; 
        n_knot_spans_uv[1] = knot_vector_v.size()-1;

        std::vector<std::vector<std::vector<int>>> knot_spans_available;
        knot_spans_available.reserve(numberOfLoops);

        for (int i = 0; i < numberOfLoops; ++i) {
            std::vector<std::vector<int>> matrix; 
            matrix.reserve(n_knot_spans_uv[1]);
            for (int j = 0; j <= n_knot_spans_uv[1]; ++j) {
                std::vector<int> row(n_knot_spans_uv[0]); 
                matrix.push_back(row); 
            }
            knot_spans_available.push_back(matrix);
        }
        
        // Optimized Snake -> for inner loops
        int idMatrixKnotSpansAvailable = 0;
        int idFirstNode;
        bool newInnerLoop = true;

        
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && is_inner)  << "Inner :: Starting SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && !is_inner) << "Outer :: Starting SnakeStep" << std::endl;

        if (skin_model_part_initial.Conditions().size() > 0) {
            
            // CREATE FIRST NODE FOR SKIN SUB MODEL PART
            auto initial_condition_in = skin_model_part_initial.GetCondition(1);
            double x_true_boundary0_in = initial_condition_in.GetGeometry()[0].X();
            double y_true_boundary0_in = initial_condition_in.GetGeometry()[0].Y();

            const int id_first_node = skin_model_part.GetRootModelPart().Nodes().size()+1;
            skin_sub_model_part.CreateNewNode(id_first_node, x_true_boundary0_in, y_true_boundary0_in, 0.0);

            for (auto &i_cond : skin_model_part_initial.Conditions()) {  
                if (newInnerLoop) {
                    idFirstNode = i_cond.GetGeometry()[0].Id();
                    newInnerLoop = false;
                }
                // Collect the coordinates of the points of the i_cond
                double x_true_boundary1 = i_cond.GetGeometry()[0].X();
                double y_true_boundary1 = i_cond.GetGeometry()[0].Y();
                double x_true_boundary2 = i_cond.GetGeometry()[1].X();
                double y_true_boundary2 = i_cond.GetGeometry()[1].Y();

                std::vector<std::vector<double>> xy_coord_i_cond(2);
                xy_coord_i_cond[0].resize(2); xy_coord_i_cond[1].resize(2); 
                
                xy_coord_i_cond[0][0] = i_cond.GetGeometry()[0].X(); // x_true_boundary1
                xy_coord_i_cond[1][0] = i_cond.GetGeometry()[0].Y(); // y_true_boundary1
                xy_coord_i_cond[0][1] = i_cond.GetGeometry()[1].X(); // x_true_boundary2
                xy_coord_i_cond[1][1] = i_cond.GetGeometry()[1].Y(); // y_true_boundary2
                
                // Collect the intersections of the skin boundary with the knot values
                std::vector<std::vector<int>> knot_span_uv(2);
                knot_span_uv[0].resize(2); knot_span_uv[1].resize(2);

                knot_span_uv[0][0] = (x_true_boundary1-starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_1st_point
                knot_span_uv[1][0] = (y_true_boundary1-starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_1st_point
                knot_span_uv[0][1] = (x_true_boundary2-starting_pos_uv[0]) / knot_step_uv[0]; // knot_span_u_2nd_point
                knot_span_uv[1][1] = (y_true_boundary2-starting_pos_uv[1]) / knot_step_uv[1]; // knot_span_v_2nd_point

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

                SnakeStep(skin_sub_model_part, knot_spans_available, idMatrixKnotSpansAvailable, 
                         knot_span_uv, xy_coord_i_cond, knot_step_uv, starting_pos_uv);
                          
                if (i_cond.GetGeometry()[1].Id() == idFirstNode) {
                    idMatrixKnotSpansAvailable++;
                    newInnerLoop = true;
                }
            }
        }
        else if (is_initial_skin_nurbs) {
            const int n_initial_points_for_side = 100000;
            int first_node_id = skin_model_part.GetRootModelPart().NumberOfNodes()+1;
            const SizeType n_boundary_curves = skin_model_part_initial.NumberOfGeometries();
            bool newInnerLoop = true;
            for (IndexType i_boundary_curve = 0; i_boundary_curve < n_boundary_curves; i_boundary_curve++) 
            {
                NurbsCurveGeometryPointerType p_curve = std::dynamic_pointer_cast<Kratos::NurbsCurveGeometry<2, Kratos::PointerVector<Kratos::Node>>>(skin_model_part_initial.pGetGeometry(i_boundary_curve));
                if (!p_curve) 
                    KRATOS_ERROR << "NURBS curve " << i_boundary_curve << " not defined in the initial Model Part. Check the importNurbsSbmModeler." << std::endl;

                // first point
                CoordinatesArrayType first_point_coords(3);
                Vector first_point_local_coord = ZeroVector(3);
                p_curve->GlobalCoordinates(first_point_coords, first_point_local_coord);
                // check the first point of the curve
                if (newInnerLoop) 
                {
                    Node::Pointer node = new Node(first_node_id, first_point_coords[0], first_point_coords[1], first_point_coords[2]);

                    // Create two nodes and two conditions for each skin condition
                    std::string layer_name = p_curve->GetValue(IDENTIFIER);

                    ModelPart& skin_layer_sub_model_part = skin_model_part.HasSubModelPart(layer_name) ? 
                                                           skin_model_part.GetSubModelPart(layer_name) : skin_model_part.CreateSubModelPart(layer_name);
            
                    skin_sub_model_part.AddNode(node);
                    newInnerLoop = false;
                } else 
                {
                    const int last_node_id = skin_model_part.GetRootModelPart().NumberOfNodes();
                    KRATOS_ERROR_IF(norm_2(first_point_coords - skin_sub_model_part.GetNode(last_node_id)) > 1e-15)
                                    << "ImportNurbsSbmModeler: error in the json NURBS file. The boundary curves" 
                                    << " are not correctly ordered." << std::endl;;
                }
                // add the specified number of points
                Vector second_point_local_coord = ZeroVector(3);
                    CoordinatesArrayType second_point_coords(3);
                for (int i = 1; i < n_initial_points_for_side; i++)
                {
                    second_point_local_coord[0] = (double) i/(n_initial_points_for_side-1);
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
                    if (knot_span_uv[0][0] == n_knot_spans_uv[0]) knot_span_uv[0][0]--; if (knot_span_uv[1][0] == n_knot_spans_uv[1]) knot_span_uv[1][0]--;
                    if (knot_span_uv[0][1] == n_knot_spans_uv[0]) knot_span_uv[0][1]--; if (knot_span_uv[1][1] == n_knot_spans_uv[1]) knot_span_uv[1][1]--;

                    std::vector<double> local_coords{first_point_local_coord[0], second_point_local_coord[0]};

                    // KRATOS_WATCH(knot_span_uv)
                    // KRATOS_WATCH(knot_spans_available.size())
                    // KRATOS_WATCH(first_point_coords)
                    // KRATOS_WATCH(second_point_coords)

                    SnakeStepNurbs(skin_sub_model_part, knot_spans_available, idMatrixKnotSpansAvailable, 
                                   knot_span_uv, xy_coord_i_cond, knot_step_uv, starting_pos_uv, local_coords, p_curve);
                    
                    first_point_local_coord = second_point_local_coord;
                    first_point_coords = second_point_coords;
                    // ********************************************************

                    // // compute normal and store value to the point
                    // std::vector<CoordinatesArrayType> global_space_derivatives;
                    // SizeType derivative_order = 1;
                    // p_curve->GlobalSpaceDerivatives(global_space_derivatives, local_coord, derivative_order);
                    // CoordinatesArrayType tangent_vector = global_space_derivatives[1];

                    // double tangent_magnitude = norm_2(tangent_vector);
                    // tangent_vector /= tangent_magnitude;
                    // Vector normal_vector = ZeroVector(3);
                    // normal_vector[0] = tangent_vector[1];
                    // normal_vector[1] = -tangent_vector[0];

                    // // node->SetValue(NORMAL, node)
                    // skin_model_part_initial.AddNode(node);

                    //-------------------------------------------------
                }
                // check the last point of the curve
                if (norm_2(second_point_coords - skin_sub_model_part.GetNode(first_node_id)) < 1e-15)
                {
                    first_node_id = skin_sub_model_part.GetRootModelPart().NumberOfNodes()+1;
                    newInnerLoop = true;
                    idMatrixKnotSpansAvailable++;
                }
            }
            // exit(0);
        }
        else {
            KRATOS_ERROR << "::[SnakeSbmUtilities]:: The initial skin model part is empty.";
        }
        PointVector points;
        for (auto &i_cond : skin_sub_model_part.Conditions()) {
            points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry()[0].X(), i_cond.GetGeometry()[0].Y(), i_cond.GetGeometry()[0].Z())));
        }
        DynamicBins testBins(points.begin(), points.end());

        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && is_inner) << "Inner :: Ending SnakeStep" << std::endl;
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && !is_inner) << "Outer :: Ending SnakeStep" << std::endl;
        
        // Read lambda parameters: 0.0 -> External,  0.5 -> Optimal
        double lambda;
        if (is_inner) 
            if (mParameters["sbm_parameters"].Has("lambda_inner")) 
                lambda = mParameters["sbm_parameters"]["lambda_inner"].GetDouble();
            else 
                lambda = 0.5;
        else 
            if (mParameters["sbm_parameters"].Has("lambda_outer")) 
                lambda = mParameters["sbm_parameters"]["lambda_outer"].GetDouble();
            else 
                lambda = 0.5;

        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && is_inner) << "Inner :: MarkKnotSpansAvailable" << std::endl;
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && !is_inner) << "Outer :: MarkKnotSpansAvailable" << std::endl;

        for (int i = 0; i < numberOfLoops; i++) {
            int idInnerLoop = i;
            // Mark the knot_spans_available's for inner and outer loops
            MarkKnotSpansAvailable(knot_spans_available, idInnerLoop, testBins, skin_sub_model_part, lambda, 
                                   n_knot_spans_uv, knot_step_uv, starting_pos_uv);  
                  
            if (is_inner) {
                CreateSurrogateBuondaryFromSnake_inner (knot_spans_available, idInnerLoop, surrogate_sub_model_part, 
                                    n_knot_spans_uv, knot_vector_u, knot_vector_v, starting_pos_uv );
                KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Inner :: Snake process has finished" << std::endl;
            }
            else {
                CreateSurrogateBuondaryFromSnake_outer (testBins, skin_sub_model_part, knot_spans_available, idInnerLoop, surrogate_sub_model_part, 
                                    n_knot_spans_uv, knot_vector_u, knot_vector_v, starting_pos_uv);
                KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0) << "Outer :: Snake process has finished" << std::endl;
            }
        }
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && is_inner) << "Inner :: Loop finished" << std::endl;
        KRATOS_INFO_IF("::[SnakeSBMUtilities]::", rEchoLevel > 0 && !is_inner) << "Outer :: Loop finished" << std::endl;
    }


    void SnakeSBMUtilities::SnakeStep(ModelPart& skin_model_part, 
                            std::vector<std::vector<std::vector<int>>> &knot_spans_available, 
                            int idMatrix, 
                            std::vector<std::vector<int>> knot_spans_uv, 
                            std::vector<std::vector<double>> xy_coord_i_cond, 
                            Vector knot_step_uv, 
                            Vector starting_pos_uv)
                            {
        bool isSplitted = false;

        if (knot_spans_uv[0][0] != knot_spans_uv[0][1] || knot_spans_uv[1][0] != knot_spans_uv[1][1]) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
            // Check if we are jumping some cut knot spans. If yes we split the true segment
            if (std::abs(knot_spans_uv[1][0]-knot_spans_uv[1][1]) > 1 || std::abs(knot_spans_uv[0][0]-knot_spans_uv[0][1]) > 1 || 
                    (knot_spans_uv[0][0] != knot_spans_uv[0][1] && knot_spans_uv[1][0] != knot_spans_uv[1][1]) ) {
                isSplitted = true;

                // Split the segment and do it recursively
                double x_true_boundary_split = (xy_coord_i_cond[0][0]+xy_coord_i_cond[0][1]) / 2;
                double y_true_boundary_split = (xy_coord_i_cond[1][0]+xy_coord_i_cond[1][1]) / 2;
                int knot_span_u_point_split = (x_true_boundary_split-starting_pos_uv[0]) / knot_step_uv[0] ;
                int knot_span_v_point_split = (y_true_boundary_split-starting_pos_uv[1]) / knot_step_uv[1] ;

                // update xy_coord for the first split segment
                std::vector<std::vector<double>> xy_coord_i_cond_split(2);
                xy_coord_i_cond_split[0].resize(2); xy_coord_i_cond_split[1].resize(2); 
                xy_coord_i_cond_split[0][0] = xy_coord_i_cond[0][0]; // x_true_boundary1
                xy_coord_i_cond_split[1][0] = xy_coord_i_cond[1][0]; // y_true_boundary1
                xy_coord_i_cond_split[0][1] = x_true_boundary_split; // x_true_boundary_split
                xy_coord_i_cond_split[1][1] = y_true_boundary_split; // y_true_boundary_split
                // update knot_span_uv for the first split segment
                std::vector<std::vector<int>> knot_span_uv_split(2);
                knot_span_uv_split[0].resize(2); knot_span_uv_split[1].resize(2); 
                knot_span_uv_split[0][0] = knot_spans_uv[0][0]; // knot_span_u_1st_point
                knot_span_uv_split[1][0] = knot_spans_uv[1][0]; // knot_span_v_1st_point
                knot_span_uv_split[0][1] = knot_span_u_point_split; // knot_span_u_point_split
                knot_span_uv_split[1][1] = knot_span_v_point_split; // knot_span_v_point_split
                
                // __We do it recursively first split__
                SnakeStep(skin_model_part, knot_spans_available, idMatrix, knot_span_uv_split, 
                         xy_coord_i_cond_split, knot_step_uv, starting_pos_uv );

                // update xy_coord for the second split segment
                xy_coord_i_cond_split[0][0] = x_true_boundary_split; // x_true_boundary_split
                xy_coord_i_cond_split[1][0] = y_true_boundary_split; // y_true_boundary_split
                xy_coord_i_cond_split[0][1] = xy_coord_i_cond[0][1]; // x_true_boundary2
                xy_coord_i_cond_split[1][1] = xy_coord_i_cond[1][1]; // y_true_boundary2
                // update knot_span_uv for the first split segment
                knot_span_uv_split[0][0] = knot_span_u_point_split; // knot_span_u_point_split
                knot_span_uv_split[1][0] = knot_span_v_point_split; // knot_span_v_point_split
                knot_span_uv_split[0][1] = knot_spans_uv[0][1]; // knot_span_u_2nd_point
                knot_span_uv_split[1][1] = knot_spans_uv[1][1]; // knot_span_v_2nd_point

                // __We do it recursively second split__
                SnakeStep(skin_model_part, knot_spans_available, idMatrix, knot_span_uv_split, 
                         xy_coord_i_cond_split, knot_step_uv, starting_pos_uv);
            }
            // Check if the true boundary crosses an u or a v knot value
            else if (knot_spans_uv[0][0] != knot_spans_uv[0][1]) { // u knot value is crossed
                // Find the "knot_spans_available" using the intersection
                knot_spans_available[idMatrix][knot_spans_uv[1][0]][knot_spans_uv[0][0]] = 2;
                knot_spans_available[idMatrix][knot_spans_uv[1][0]][knot_spans_uv[0][1]] = 2;

            }
            else if (knot_spans_uv[1][0] != knot_spans_uv[1][1]) { // v knot value is crossed
                // Find the "knot_spans_available" using the intersection (Snake_coordinate classic -> External Boundary)
                knot_spans_available[idMatrix][knot_spans_uv[1][0]][knot_spans_uv[0][0]] = 2;
                knot_spans_available[idMatrix][knot_spans_uv[1][1]][knot_spans_uv[0][0]] = 2;
            }
        }
        if (!isSplitted) {
            // Call the root model part for the Ids of the node
            auto idNode1 = skin_model_part.GetRootModelPart().Nodes().size();
            auto idNode2 = idNode1+1;
            // Create two nodes and two conditions for each skin condition
            skin_model_part.CreateNewNode(idNode2, (xy_coord_i_cond[0][0]+xy_coord_i_cond[0][1] ) / 2, (xy_coord_i_cond[1][0]+xy_coord_i_cond[1][1] ) / 2, 0.0);
            skin_model_part.CreateNewNode(idNode2+1, xy_coord_i_cond[0][1], xy_coord_i_cond[1][1], 0.0);
            Properties::Pointer p_cond_prop = skin_model_part.pGetProperties(0);
            Condition::Pointer p_cond1 = skin_model_part.CreateNewCondition("LineCondition2D2N", idNode1, {{idNode1, idNode2}}, p_cond_prop );
            Condition::Pointer p_cond2 = skin_model_part.CreateNewCondition("LineCondition2D2N", idNode2, {{idNode2, idNode2+1}}, p_cond_prop );
            skin_model_part.AddCondition(p_cond1);
            skin_model_part.AddCondition(p_cond2);
        }
    }

    void SnakeSBMUtilities::SnakeStepNurbs(ModelPart& skin_model_part, 
                            std::vector<std::vector<std::vector<int>>> &knot_spans_available, 
                            int idMatrix, 
                            std::vector<std::vector<int>> knot_spans_uv, 
                            std::vector<std::vector<double>> xy_coord, 
                            Vector knot_step_uv, 
                            const Vector starting_pos_uv,
                            const std::vector<double> local_coords,
                            const NurbsCurveGeometryPointerType p_curve)
                            {
        bool isSplitted = false;

        if (knot_spans_uv[0][0] != knot_spans_uv[0][1] || knot_spans_uv[1][0] != knot_spans_uv[1][1]) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
            // Check if we are jumping some cut knot spans. If yes we split the true segment
            if (std::abs(knot_spans_uv[1][0]-knot_spans_uv[1][1]) > 1 || std::abs(knot_spans_uv[0][0]-knot_spans_uv[0][1]) > 1 || 
                    (knot_spans_uv[0][0] != knot_spans_uv[0][1] && knot_spans_uv[1][0] != knot_spans_uv[1][1]) ) {
                isSplitted = true;
                
                // Split the segment and do it recursively
                Vector local_coords_split = ZeroVector(3);
                local_coords_split[0] = (local_coords[0] + local_coords[1]) / 2;
                CoordinatesArrayType xy_true_boundary_split;
                p_curve->GlobalCoordinates(xy_true_boundary_split, local_coords_split);

                int knot_span_u_point_split = (xy_true_boundary_split[0]-starting_pos_uv[0]) / knot_step_uv[0] ;
                int knot_span_v_point_split = (xy_true_boundary_split[1]-starting_pos_uv[1]) / knot_step_uv[1] ;

                // update xy_coord for the first split segment
                std::vector<std::vector<double>> xy_coord_i_cond_split(2);
                xy_coord_i_cond_split[0].resize(2); xy_coord_i_cond_split[1].resize(2); 
                xy_coord_i_cond_split[0][0] = xy_coord[0][0]; // x_true_boundary1
                xy_coord_i_cond_split[1][0] = xy_coord[1][0]; // y_true_boundary1
                xy_coord_i_cond_split[0][1] = xy_true_boundary_split[0]; // x_true_boundary_split
                xy_coord_i_cond_split[1][1] = xy_true_boundary_split[1]; // y_true_boundary_split
                // update knot_span_uv for the first split segment
                std::vector<std::vector<int>> knot_span_uv_split(2);
                knot_span_uv_split[0].resize(2); knot_span_uv_split[1].resize(2); 
                knot_span_uv_split[0][0] = knot_spans_uv[0][0]; // knot_span_u_1st_point
                knot_span_uv_split[1][0] = knot_spans_uv[1][0]; // knot_span_v_1st_point
                knot_span_uv_split[0][1] = knot_span_u_point_split; // knot_span_u_point_split
                knot_span_uv_split[1][1] = knot_span_v_point_split; // knot_span_v_point_split

                std::vector<double> local_coords_split_segment1{local_coords[0], local_coords_split[0]};
                
                // KRATOS_WATCH(knot_spans_uv)
                // KRATOS_WATCH(local_coords)
                // KRATOS_WATCH(local_coords_split)
                // KRATOS_WATCH(xy_coord)
                // exit(0);
                // __We do it recursively first split__
                SnakeStepNurbs(skin_model_part, knot_spans_available, idMatrix, knot_span_uv_split, 
                         xy_coord_i_cond_split, knot_step_uv, starting_pos_uv, local_coords_split_segment1, p_curve);

                // update xy_coord for the second split segment
                xy_coord_i_cond_split[0][0] = xy_true_boundary_split[0]; // x_true_boundary_split
                xy_coord_i_cond_split[1][0] = xy_true_boundary_split[1]; // y_true_boundary_split
                xy_coord_i_cond_split[0][1] = xy_coord[0][1]; // x_true_boundary2
                xy_coord_i_cond_split[1][1] = xy_coord[1][1]; // y_true_boundary2
                // update knot_span_uv for the first split segment
                knot_span_uv_split[0][0] = knot_span_u_point_split; // knot_span_u_point_split
                knot_span_uv_split[1][0] = knot_span_v_point_split; // knot_span_v_point_split
                knot_span_uv_split[0][1] = knot_spans_uv[0][1]; // knot_span_u_2nd_point
                knot_span_uv_split[1][1] = knot_spans_uv[1][1]; // knot_span_v_2nd_point

                std::vector<double> local_coords_split_segment2{local_coords_split[0], local_coords[1]};
                // __We do it recursively second split__
                SnakeStepNurbs(skin_model_part, knot_spans_available, idMatrix, knot_span_uv_split, 
                         xy_coord_i_cond_split, knot_step_uv, starting_pos_uv, local_coords_split_segment2, p_curve);



            }
            // Check if the true boundary crosses an u or a v knot value
            else if (knot_spans_uv[0][0] != knot_spans_uv[0][1]) { // u knot value is crossed
                // Find the "knot_spans_available" using the intersection
                knot_spans_available[idMatrix][knot_spans_uv[1][0]][knot_spans_uv[0][0]] = 2;
                knot_spans_available[idMatrix][knot_spans_uv[1][0]][knot_spans_uv[0][1]] = 2;

            }
            else if (knot_spans_uv[1][0] != knot_spans_uv[1][1]) { // v knot value is crossed
                // Find the "knot_spans_available" using the intersection (Snake_coordinate classic -> External Boundary)
                knot_spans_available[idMatrix][knot_spans_uv[1][0]][knot_spans_uv[0][0]] = 2;
                knot_spans_available[idMatrix][knot_spans_uv[1][1]][knot_spans_uv[0][0]] = 2;
            }
        }
        if (!isSplitted) {
            // Call the root model part for the Ids of the node
            auto idNode1 = skin_model_part.GetRootModelPart().Nodes().size();
            auto idNode2 = idNode1+1;
            // Create two nodes and two conditions for each skin condition
            auto node = new Node(idNode2, xy_coord[0][1], xy_coord[1][1], 0.0);

            std::string layer_name = p_curve->GetValue(IDENTIFIER);

            ModelPart& skin_layer_sub_model_part = skin_model_part.HasSubModelPart(layer_name) ? 
                                                skin_model_part.GetSubModelPart(layer_name) : skin_model_part.CreateSubModelPart(layer_name);
            
            // compute normal and tangent informations at the local coord of the point 
            skin_layer_sub_model_part.AddNode(node);

            Properties::Pointer p_cond_prop = skin_model_part.pGetProperties(0);
            Condition::Pointer p_cond = skin_model_part.CreateNewCondition("LineCondition2D2N", idNode1, {{idNode1, idNode2}}, p_cond_prop );
            skin_model_part.AddCondition(p_cond);
        }
    }


    bool SnakeSBMUtilities::isPointInsideSkinBoundary(Point& point1, DynamicBins& testBins, ModelPart& skin_model_part)
    {
        // Get the nearest point of the true boundary
        PointerType pointToSearch = PointerType(new PointType(1000000, point1.X(), point1.Y(), 0.0));
        PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);
        
        // Get the closest Condition the initial_skin_model_part_in.Conditions
        int id1 = nearestPoint->Id();
        auto nearestCondition1 = skin_model_part.GetCondition(id1);
        // Check if the condition is the first one and therefore the previous one does not exist
        int id2 = id1 - 1;
        if (id1 == skin_model_part.ConditionsBegin()->Id()) {
            int nConditions = skin_model_part.Conditions().size();
            id2 = id1 + nConditions - 1; 
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



    void SnakeSBMUtilities::MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idMatrix,DynamicBins& testBin, 
                                                   ModelPart& skin_model_part, double lambda, std::vector<int>& n_knot_spans_uv, 
                                                   Vector& knot_step_uv, const Vector& starting_pos_uv) {
        for (int i = 0; i < n_knot_spans_uv[1]; i++) {
            for (int j = 0; j < n_knot_spans_uv[0]; j++) {
                if (knot_spans_available[idMatrix][i][j] == 2) {
                    // Check the 8 neighbor knot spans -> Is there any completely inside? Note that we can just check 1 point.
                    
                    // right node
                    if (i != n_knot_spans_uv[1]-1)
                        if (knot_spans_available[idMatrix][i+1][j] == 0) { 
                            Point gaussPoint = Point((j+0.5) * knot_step_uv[0] + starting_pos_uv[0], (i+1+0.5) * knot_step_uv[1] +starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i+1][j] = 1;}
                        }
                    // left node    
                    if (i != 0)
                        if (knot_spans_available[idMatrix][i-1][j] == 0) { 
                            Point gaussPoint = Point((j+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i-1+0.5) * knot_step_uv[1] + starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i-1][j] = 1;}
                        }
                    // up node
                    if (j != n_knot_spans_uv[0]-1)
                        if (knot_spans_available[idMatrix][i][j+1] == 0) { 
                            Point gaussPoint = Point((j+1+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i+0.5) * knot_step_uv[1]+starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i][j+1] = 1;}
                        }
                    //down node
                    if (j != 0)
                        if (knot_spans_available[idMatrix][i][j-1] == 0) { 
                            Point gaussPoint = Point((j-1+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i+0.5) * knot_step_uv[1]+starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i][j-1] = 1;}
                        } 

                    // corner right-down node
                    if (j != 0 && i != n_knot_spans_uv[1]-1)
                        if (knot_spans_available[idMatrix][i+1][j-1] == 0) {
                            Point gaussPoint = Point((j-1+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i+1+0.5) * knot_step_uv[1]+starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i+1][j-1] = 1;}
                        }
                    // corner left-down node
                    if (j != 0 && i != 0)
                        if (knot_spans_available[idMatrix][i-1][j-1] == 0) {
                            Point gaussPoint = Point((j-1+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i-1+0.5) * knot_step_uv[1]+starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i-1][j-1] = 1;}
                        }
                    // corner right-up node
                    if (j != n_knot_spans_uv[0]-1 && i != n_knot_spans_uv[1]-1)
                        if (knot_spans_available[idMatrix][i+1][j+1] == 0) {
                            Point gaussPoint = Point((j+1+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i+1+0.5) * knot_step_uv[1]+starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i+1][j+1] = 1;}
                        }
                    // corner left-up node
                    if (j != n_knot_spans_uv[0]-1 && i != 0)
                        if (knot_spans_available[idMatrix][i-1][j+1] == 0) {
                            Point gaussPoint = Point((j+1+0.5) * knot_step_uv[0]+starting_pos_uv[0], (i-1+0.5) * knot_step_uv[1]+starting_pos_uv[1], 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBin, skin_model_part)) {knot_spans_available[idMatrix][i-1][j+1] = 1;}
                        }

                    // Create 25 "fake" GaussPoints to check if the majority are inside or outside
                    const int numFakeGaussPoints = 8;
                    int numberOfInsideGaussianPoints = 0;
                    for (int i_GPx = 0; i_GPx < numFakeGaussPoints; i_GPx++){
                        double x_coord = j*knot_step_uv[0] + knot_step_uv[0]/(numFakeGaussPoints+1)*(i_GPx+1) + starting_pos_uv[0];

                        // NOTE:: The v-knot spans are upside down in the matrix!!
                        for (int i_GPy = 0; i_GPy < numFakeGaussPoints; i_GPy++) 
                        {
                            double y_coord = i*knot_step_uv[1] + knot_step_uv[1]/(numFakeGaussPoints+1)*(i_GPy+1) + starting_pos_uv[1];
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

                    // exit(0);
                            
                }
            }
        }
    }


    void SnakeSBMUtilities::CreateSurrogateBuondaryFromSnake_inner (std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idMatrix, 
                            ModelPart& surrogate_model_part_inner, std::vector<int>& n_knot_spans_uv, 
                            Vector& knot_vector_u, Vector&  knot_vector_v,
                            const Vector& starting_pos_uv){
        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        for (int i = 0; i < (n_knot_spans_uv[0]); i++) {
            for (int j = 0; j < (n_knot_spans_uv[1]); j++ ) {
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
        
        if (!surrogate_model_part_inner.HasProperties(0)) {surrogate_model_part_inner.CreateNewProperties(1000);}

        Properties::Pointer p_cond_prop = surrogate_model_part_inner.pGetProperties(0);
        Node snakeNode(1 , knot_vector_u[start_i], knot_vector_v[start_j], 0.0);
        
        const int id_first_node = surrogate_model_part_inner.GetRootModelPart().Nodes().size() + 1;
        surrogate_model_part_inner.CreateNewNode(id_first_node, snakeNode);
        IndexType idSnakeNode = id_first_node+1;

        // Follow the clockwise loop
        int end = 0;
        // We are going orizontally
        int direction = 0;
        // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int i = start_i; int j = start_j;
        int I = start_i; int J = start_j;
        int steps = 0;

        const int max_number_of_steps = 1e5;
        while (end == 0 && steps < max_number_of_steps) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knot_spans_available[idMatrix][J][I] == 1) {
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; risetto e muovo
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}
                if (knot_spans_available[idMatrix][J][I] == 1) {
                    // going straight
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                    surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                    idSnakeNode++;
                }
                else {
                    // Need to search to hte right; rese e muovo
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}
                    if (knot_spans_available[idMatrix][J][I] == 1) {
                        // We are moving to the right -> First move straight, store, the move to the right (i,j), store again
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        Node snakeNode1(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
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
                        if (knot_spans_available[idMatrix][J][I] == 1) {
                            // First straight, store, then move to the right, store again, then move to the right and store again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }
                            Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }

                            Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;

                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }

                            Node snakeNode3(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_inner.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = surrogate_model_part_inner.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;

                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else
                            KRATOS_ERROR << "SnakeSBMUtilities:: " <<  "error in the Snakes Coordinates" << std::endl;
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i && J == start_j && is_special_case == 1 ) {
                // End of the while loop
                end = 1;
                KRATOS_INFO("number of steps in the snake step") << steps << std::endl;
            }
        }
        KRATOS_ERROR_IF(steps >= max_number_of_steps-1) << "SnakeSBMUtilities:: " <<  "reached maximum number of steps" << std::endl;

        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        IndexType initialId = id_first_node;
        if (surrogate_model_part_inner.Elements().size()>0) {
            // Check if it is not the first inner loop
            initialId = surrogate_model_part_inner.GetElement(surrogate_model_part_inner.GetRootModelPart().Elements().size()).GetGeometry()[1].Id()+1;
        }
        std::vector<ModelPart::IndexType> elem_nodes{initialId, idSnakeNode-1};
        surrogate_model_part_inner.CreateNewElement("Element2D2N", surrogate_model_part_inner.GetRootModelPart().Elements().size()+1, elem_nodes, p_cond_prop);
    }


    void SnakeSBMUtilities::CreateSurrogateBuondaryFromSnake_outer (DynamicBins& testBin_out, ModelPart& initial_skin_model_part_out, 
                            std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idMatrix, 
                            ModelPart& surrogate_model_part_outer, std::vector<int>& n_knot_spans_uv, 
                            Vector& knot_vector_u, Vector&  knot_vector_v,
                            const Vector& starting_pos_uv){
                                             
        // CHECK ALL THE EXTERNAL KNOT SPANS

        // LEFT BOUNDARY
        double knot_step_u = knot_vector_u[1]-knot_vector_u[0];
        double knot_step_v = knot_vector_v[1]-knot_vector_v[0];

        for (int i = 0; i<2; i++) {
            for (int j = 0; j < (n_knot_spans_uv[1]); j++ ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+starting_pos_uv[0], (i+0.5)*knot_step_v+starting_pos_uv[1], 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }
        // TOP BOUNDARY
        for (int j = knot_spans_available[idMatrix][0].size()-1; j > knot_spans_available[idMatrix][0].size()-3; j--) {
            for (int i = 0; i < (n_knot_spans_uv[0]); i++) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+starting_pos_uv[0], (i+0.5)*knot_step_v+starting_pos_uv[1], 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }
        // RIGHT BOUNDARY
        for (int i = knot_spans_available[idMatrix][0].size()-1; i > knot_spans_available[idMatrix][0].size()-3; i--) {
            for (int j = n_knot_spans_uv[1]-1; j > -1; j-- ) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+starting_pos_uv[0], (i+0.5)*knot_step_v+starting_pos_uv[1], 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }
        // BOTTOM BOUNDARY
        for (int j = 0; j<2; j++) {
            for (int i = n_knot_spans_uv[0]-1; i > -1 ; i--) {
                Point centroidKnotSpan = Point((j+0.5)*knot_step_u+starting_pos_uv[0], (i+0.5)*knot_step_v+starting_pos_uv[1], 0);
                if (isPointInsideSkinBoundary(centroidKnotSpan, testBin_out, initial_skin_model_part_out) && knot_spans_available[idMatrix][i][j] != -1) {
                    knot_spans_available[idMatrix][i][j] = 1;
                    }
            }
        }

        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        for (int i = 0; i < (n_knot_spans_uv[0]); i++) {
            for (int j = 0; j < (n_knot_spans_uv[1]); j++ ) {
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
        
        // EXTEND THE MATRIX
        //-------------------
        std::vector<std::vector<int>> knot_spans_available_extended(n_knot_spans_uv[1]+2, std::vector<int>(n_knot_spans_uv[0]+2));

        for (int i = 0; i < knot_spans_available[idMatrix].size(); i++){
            for (int j = 0; j < knot_spans_available[idMatrix][0].size(); j++) {
                knot_spans_available_extended[i+1][j+1] = knot_spans_available[idMatrix][i][j]; 
            }
        }  
        
        Properties::Pointer p_cond_prop = surrogate_model_part_outer.CreateNewProperties(1001);
        Node snakeNode(1 , knot_vector_u[start_i], knot_vector_v[start_j], 0.0);

        const int id_first_node = surrogate_model_part_outer.GetRootModelPart().Nodes().size() + 1;
        
        surrogate_model_part_outer.CreateNewNode(id_first_node, snakeNode);
        IndexType idSnakeNode = id_first_node+1;
        
        // Follow the clockwise loop
        int end = 0;
        // We are going orizontally
        int direction = 0 ;      // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int I = start_i+1;
        int J = start_j+1; 
        int i = start_i;
        int j = start_j;
        int steps = 0;
        const int max_number_of_steps = 1e5;
        while (end == 0 && steps < max_number_of_steps) {
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

                    Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                    surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
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
                        Node snakeNode1(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
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
                            Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }
                            Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }
                            Node snakeNode3(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part_outer.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = surrogate_model_part_outer.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else
                            KRATOS_ERROR << "SnakeSBMUtilities:: " <<  "error in the Snakes Coordinates" << std::endl;
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i+1 && J == start_j+1 && is_special_case == 1 ) {
                // End of the while loop
                end = 1;
                }
        }
        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        std::vector<ModelPart::IndexType> elem_nodes{1, idSnakeNode-1};
        int elem_id = surrogate_model_part_outer.GetRootModelPart().Elements().size()+1;
        surrogate_model_part_outer.CreateNewElement("Element2D2N", elem_id, elem_nodes, p_cond_prop);
    }
}  // namespace Kratos.
