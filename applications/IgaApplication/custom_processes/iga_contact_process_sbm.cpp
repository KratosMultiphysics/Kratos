//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

// Project includes
#include "iga_contact_process_sbm.h"
#include "geometries/nurbs_volume_geometry.h"

#include "utilities/variable_utils.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"
#include "utilities/nurbs_utilities/projection_nurbs_contact_utilities.h"

namespace Kratos
{

    IgaContactProcessSbm::IgaContactProcessSbm(
        Model& rModel, Parameters ThisParameters) : 
        Process(), 
        mpModel(&rModel), 
        mParameters(ThisParameters)
    {
        KRATOS_WATCH("inittttttttttttttttttttttttttttttttttttttttttttt")
        // TODO: SBM->Body Fitted case or viceversa
        // the class only does the case SBM->SBM at the moment
        
        mEchoLevel = mParameters["echo_level"].GetInt();

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("analysis_model_part_name")) << "::[IgaContactProcessSbm]::" 
                            << " Missing \"analysis_model_part_name\" parameter. "<< std::endl; 
        
        KRATOS_ERROR_IF_NOT(ThisParameters.Has("contact_sub_model_part_name")) << "::[IgaContactProcessSbm]::" 
                            << " Missing \"contact_sub_model_part_name\" parameter. "<< std::endl; 

        //-----------------------------------------------------------------------------------------------
        // Obtain SLAVE interface b_reps
        std::string slave_model_part_name = mParameters["contact_parameters"]["slave_model_part"]["sub_model_part_name"].GetString();
        const std::string slave_layer_name = mParameters["contact_parameters"]["slave_model_part"]["layer_name"].GetString();
        slave_model_part_name += "." + slave_layer_name;


        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_model_part_name)) << "ERROR: SLAVE MODEL PART " 
                                                << slave_model_part_name << "NOT CREATED" << std::endl; 

        mrSlaveModelPart = &(mpModel->GetModelPart(slave_model_part_name));


        const IndexType slave_property_id = mParameters["contact_parameters"]["slave_model_part"]["property_id"].GetInt();
        mpPropSlave = mrSlaveModelPart->pGetProperties(slave_property_id);

        // Obtain the slave skin model part
        std::string slave_skin_model_part_name = mParameters["contact_parameters"]["slave_model_part"]["initial_skin_model_part_name"].GetString();

        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_skin_model_part_name)) << "ERROR: SLAVE SKIN MODEL PART " 
                                                << slave_skin_model_part_name << "NOT CREATED" << std::endl; 

        mrSlaveSkinModelPart = &(mpModel->GetModelPart(slave_skin_model_part_name));

        //--------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------
    
        // Obtain MASTER interface b_reps
        std::string master_model_part_name = mParameters["contact_parameters"]["master_model_part"]["sub_model_part_name"].GetString();
        const std::string master_layer_name = mParameters["contact_parameters"]["master_model_part"]["layer_name"].GetString();
        master_model_part_name += "." + master_layer_name;


        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_model_part_name)) << "ERROR: MASTER MODEL PART " 
                                                << master_model_part_name << "NOT CREATED" << std::endl; 

        mrMasterModelPart = &(mpModel->GetModelPart(master_model_part_name));

        const IndexType master_property_id = mParameters["contact_parameters"]["master_model_part"]["property_id"].GetInt();
        mpPropMaster = mrMasterModelPart->pGetProperties(master_property_id);

        // Obtain the slave skin model part
        std::string master_skin_model_part_name = mParameters["contact_parameters"]["master_model_part"]["initial_skin_model_part_name"].GetString();

        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_skin_model_part_name)) << "ERROR: MASTER SKIN MODEL PART " 
                                                << master_skin_model_part_name << "NOT CREATED" << std::endl; 

        mrMasterSkinModelPart = &(mpModel->GetModelPart(master_skin_model_part_name));

        //------------------------------------------------------------------------
        std::string analysis_model_part_name = mParameters["analysis_model_part_name"].GetString();
        std::string contact_sub_model_part_name = mParameters["contact_sub_model_part_name"].GetString();

        std::string contact_model_part_name = analysis_model_part_name + ".ContactInterface." + contact_sub_model_part_name;

        mrContactModelPart = &(mpModel->CreateModelPart(contact_model_part_name));

        //********************************************************************** */
        const double epsilon = 1e-12;
        // master
        Vector master_knot_step_uv = mrMasterModelPart->GetParentModelPart().GetValue(MARKER_MESHES);
        Vector master_parameter_space_extremes = mrMasterModelPart->GetParentModelPart().GetValue(LOAD_MESHES);
        double master_domain_length_u = master_parameter_space_extremes[2] - master_parameter_space_extremes[0];
        double master_domain_length_v = master_parameter_space_extremes[3] - master_parameter_space_extremes[1];

        const SizeType master_number_of_knot_spans_u = floor((master_domain_length_u + epsilon)/master_knot_step_uv[0]);
        const SizeType master_number_of_knot_spans_v = floor((master_domain_length_v + epsilon)/master_knot_step_uv[1]);
        mSparseBrepMatrixMaster.resize(master_number_of_knot_spans_u+1, master_number_of_knot_spans_v+1, false);

        for (auto& r_geometry_master : mrMasterModelPart->Geometries()) {

            int master_brep_id = r_geometry_master.Id();
            auto p_master_geometry = mrMasterModelPart->pGetGeometry(master_brep_id);

            // compute the center of each surrogate brep
            Vector brep_domain_interval;
            p_master_geometry->DomainInterval(brep_domain_interval);
            CoordinatesArrayType brep_center = ZeroVector(3);
            CoordinatesArrayType brep_center_local_coords = ZeroVector(3);
            brep_center_local_coords[0] = (brep_domain_interval[0]+brep_domain_interval[1])/2;
            p_master_geometry->GlobalCoordinates(brep_center, brep_center_local_coords);

            IndexType knot_span_surrounded_u_id = floor((brep_center[0]- master_parameter_space_extremes[0]+epsilon)/master_knot_step_uv[0]);
            IndexType knot_span_surrounded_v_id = floor((brep_center[1]- master_parameter_space_extremes[1]+epsilon)/master_knot_step_uv[1]);
            mSparseBrepMatrixMaster(knot_span_surrounded_u_id, knot_span_surrounded_v_id) = master_brep_id;
        }

        //---------------------
        // slave
        Vector slave_knot_step_uv = mrSlaveModelPart->GetParentModelPart().GetValue(MARKER_MESHES);
        Vector slave_parameter_space_extremes = mrSlaveModelPart->GetParentModelPart().GetValue(LOAD_MESHES);
        double slave_domain_length_u = slave_parameter_space_extremes[2] - slave_parameter_space_extremes[0];
        double slave_domain_length_v = slave_parameter_space_extremes[3] - slave_parameter_space_extremes[1];

        const SizeType slave_number_of_knot_spans_u = floor((slave_domain_length_u + epsilon)/slave_knot_step_uv[0]);
        const SizeType slave_number_of_knot_spans_v = floor((slave_domain_length_v + epsilon)/slave_knot_step_uv[1]);
        mSparseBrepMatrixSlave.resize(slave_number_of_knot_spans_u+1, slave_number_of_knot_spans_v+1,false);

        for (auto& r_geometry_slave : mrSlaveModelPart->Geometries()) {

            int slave_brep_id = r_geometry_slave.Id();
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);

            // compute the center of each surrogate brep
            Vector brep_domain_interval;
            p_slave_geometry->DomainInterval(brep_domain_interval);
            CoordinatesArrayType brep_center = ZeroVector(3);
            CoordinatesArrayType brep_center_local_coords = ZeroVector(3);
            brep_center_local_coords[0] = (brep_domain_interval[0]+brep_domain_interval[1])/2;
            p_slave_geometry->GlobalCoordinates(brep_center, brep_center_local_coords);
            IndexType knot_span_surrounded_u_id = floor((brep_center[0]- slave_parameter_space_extremes[0]+epsilon)/slave_knot_step_uv[0]);
            IndexType knot_span_surrounded_v_id = floor((brep_center[1]- slave_parameter_space_extremes[1]+epsilon)/slave_knot_step_uv[1]);

            mSparseBrepMatrixSlave(knot_span_surrounded_u_id, knot_span_surrounded_v_id) = slave_brep_id;
        }

        // KRATOS_WATCH(master_knot_step_uv)
        // KRATOS_WATCH(master_number_of_knot_spans_u)
        // KRATOS_WATCH(master_number_of_knot_spans_v)
        // KRATOS_WATCH(slave_knot_step_uv)
        // KRATOS_WATCH(slave_number_of_knot_spans_u)
        // KRATOS_WATCH(slave_number_of_knot_spans_v)
        // KRATOS_WATCH(SparseBrepMatrixMaster)
        // KRATOS_WATCH(SparseBrepMatrixSlave)
        // exit(0);

   }

    void IgaContactProcessSbm::Execute(){

        std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::trunc);

        const std::string master_layer_name = mParameters["contact_parameters"]["master_model_part"]["layer_name"].GetString();
        const std::string slave_layer_name = mParameters["contact_parameters"]["slave_model_part"]["layer_name"].GetString();

        ConditionsArrayType& r_conditions_array = mrContactModelPart->GetParentModelPart().Conditions();
        KRATOS_TRACE_IF("Empty model part", r_conditions_array.size() == 0) << "YOUR CONTACT MODEL PART IS EMPTY" << std::endl;
        block_for_each(r_conditions_array, [&](Condition& rCond) {
                rCond.Set(TO_ERASE, true);
        });

        NodesArrayType& r_nodes_array = mrContactModelPart->GetParentModelPart().Nodes();
        block_for_each(r_nodes_array, [&](Node& rNode) {
                rNode.Set(TO_ERASE, true);
        });

        mrContactModelPart->GetParentModelPart().RemoveConditionsFromAllLevels(TO_ERASE);
        mrContactModelPart->GetParentModelPart().RemoveNodesFromAllLevels(TO_ERASE);

        SizeType id = 1;
        if (mrContactModelPart->GetRootModelPart().Conditions().size() > 0)
            id = mrContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;
        
        KRATOS_ERROR_IF_NOT(mParameters.Has("name"))
            << "\"name\" need to be specified." << std::endl;
        std::string name = mParameters["name"].GetString();

        //************************** */
        // PREPARE TESTBINS WITH INTEGRATION POINTS OF THE SLAVE SIDE
        //------------------------------------------------------------
        PointVector points;
        Vector meshSizes_uv = mrSlaveModelPart->GetParentModelPart().GetValue(MARKER_MESHES);
        IndexType count_slave_brep = 0;
        for (auto& r_geometry_slave : mrSlaveModelPart->Geometries()) {

            int slave_brep_id = r_geometry_slave.Id();
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);
            p_slave_geometry->SetValue(IDENTIFIER, "non_active");

            // add the center of each surrogate brep
            Vector brep_domain_interval;
            p_slave_geometry->DomainInterval(brep_domain_interval);
            CoordinatesArrayType brep_center = ZeroVector(3);
            CoordinatesArrayType brep_center_local_coords = ZeroVector(3);
            brep_center_local_coords[0] = (brep_domain_interval[0]+brep_domain_interval[1])/2;
            p_slave_geometry->GlobalCoordinates(brep_center, brep_center_local_coords);

            points.push_back(PointTypePointer(new PointType(p_slave_geometry->Id(), brep_center)));

            count_slave_brep++;
        }

        
        
        DynamicBins testBins(points.begin(), points.end());
        const int numberOfResults = 1e2; 
        ModelPart::NodesContainerType::ContainerType Results(numberOfResults);
        std::vector<double> list_of_distances(numberOfResults);
        double meshSize = meshSizes_uv[0];
        if (meshSizes_uv[1] > meshSize) {meshSize = meshSizes_uv[1];}
        double radius = 3*sqrt(3)*(meshSize); 

        // PointerType pointToSearch = PointerType(new PointType(1, 0.0195833, 0, 0));

        Vector master_knot_step_uv = mrMasterModelPart->GetParentModelPart().GetValue(MARKER_MESHES);

        // int obtainedResults = testBins.SearchInRadius(*pointToSearch, radius, Results.begin(), list_of_distances.begin(), numberOfResults);

        //********************************************************
        // ********************************************************** */
        //FIXME:
        bool is_surrogate_correct = false;
        for (auto& r_geometry_master : mrMasterModelPart->Geometries()) {
            int master_brep_id = r_geometry_master.Id();
            auto p_master_geometry = mrMasterModelPart->pGetGeometry(master_brep_id);
            auto master_brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_master_geometry);

            KRATOS_ERROR_IF(!master_brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << master_brep_id << "is not a Brep." << std::endl;

            GeometriesArrayType gp_list_master;
            master_brep_geometry->GetQuadraturePointGeometries(gp_list_master);

            // initialize the resulting geometries array
            GeometriesArrayType rResultGeometries(gp_list_master.size());

            // compute the projection for each closest point to each gauss point
            IndexType count_gp_master = 0;
            bool is_converged_at_least_once_in_brep = false;

            for (IndexType i_gp_in_brep = 0; i_gp_in_brep < gp_list_master.size(); i_gp_in_brep++) {
                auto& gp_in_brep = gp_list_master(i_gp_in_brep);
                // collect the skin projection
                Node& skin_node = gp_in_brep->GetValue(NEIGHBOUR_NODES)[0];

                //FIXME:
                if (gp_in_brep->Center()[0] < 1e-1 && gp_in_brep->Center()[1]< 8.11)
                {   
                    if (skin_node.X() < 1e-1 && skin_node.Y()< 8.1)
                        is_surrogate_correct = true;
                    
                    // KRATOS_WATCH(gp_in_brep->Center())
                    // KRATOS_WATCH(skin_node)
                }

                // project the skin_node to the skin boundary on the slave side
                double best_curve_distance = 1e16;
                CoordinatesArrayType best_curve_projection;
                CoordinatesArrayType best_curve_projection_local_coord;
                bool is_converged_at_least_once = false;
                IndexType best_slave_curve_id;

                CoordinatesArrayType skin_node_deformed_coordinates;
                GetDeformedPosition(skin_node, *mrMasterModelPart, mSparseBrepMatrixMaster, skin_node_deformed_coordinates);

                for (auto &i_slave_curve : mrSlaveSkinModelPart->Geometries())
                { 
                    if (i_slave_curve.GetValue(IDENTIFIER) != slave_layer_name) continue;
                    CoordinatesArrayType local_coord = ZeroVector(3);
                    local_coord[0] = 0; //initial guess
                    CoordinatesArrayType projected_point;
                    double distance;

                    // bool is_converged = ProjectionNurbsContactUtilities<PointType, PointerVector<NodeType>>::NewtonRaphsonCurve(
                    //                     local_coord,
                    //                     skin_node,
                    //                     projected_point,
                    //                     i_slave_curve,
                    //                     distance,
                    //                     20,
                    //                     1e-9);


                    int slave_nurbs_curve_id = i_slave_curve.Id();
                    auto p_slave_nurbs_curve_geometry = mrSlaveSkinModelPart->pGetGeometry(slave_nurbs_curve_id);
                    auto slave_nurbs_curve_geometry = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_slave_nurbs_curve_geometry);

                    KRATOS_ERROR_IF(!slave_nurbs_curve_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << slave_nurbs_curve_id 
                                            << " is not a NurbsCurveGeometryType." << std::endl;

                    bool is_converged = NewtonRaphsonCurveOnDeformed(
                                        local_coord,
                                        skin_node_deformed_coordinates,
                                        projected_point,
                                        *slave_nurbs_curve_geometry,
                                        distance,
                                        20,
                                        10,
                                        1e-9);

                    if (is_converged || distance < 1e-2)
                    {
                        is_converged_at_least_once = true;
                    
                        if (distance < best_curve_distance) {
                            best_curve_distance = distance;
                            best_curve_projection = projected_point;
                            best_curve_projection_local_coord = local_coord; 
                            best_slave_curve_id = i_slave_curve.Id();

                            std::vector<array_1d<double, 3>> ppp;
                            //FIXME:
                            slave_nurbs_curve_geometry->GlobalSpaceDerivatives(ppp, best_curve_projection_local_coord, 1);
                            best_curve_projection = ppp[0];
                        }
                    } 
                }
                if (!is_converged_at_least_once) continue;
                else is_converged_at_least_once_in_brep = true;
                // associate the best skin_node with its projection
                Node& best_skin_node = gp_in_brep->GetValue(NEIGHBOUR_NODES)[0];

                IndexType idNewNode = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size()+1;
                auto new_slave_skin_node = new Node(idNewNode, best_curve_projection);


                // std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::app);
                //         outputFile <<  best_skin_node.X() << " " << best_skin_node.Y() << " "  << new_slave_skin_node->X() << " " << new_slave_skin_node->Y() <<"\n";
                //         outputFile.close();

                std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::app);
                        outputFile <<  best_skin_node.X() << " " << best_skin_node.Y() << " "  << gp_in_brep->Center().X() << " " << gp_in_brep->Center().Y() <<"\n";
                        outputFile.close();

                //-------------------------------------------------------------
                //  START SEARCH
                //-------------------------------------------------------------
                NodePointerVector neighbour_vector;
                neighbour_vector.push_back(new_slave_skin_node); 
                best_skin_node.SetValue(NEIGHBOUR_NODES, neighbour_vector);

                // compute normals and useful info
                std::vector<CoordinatesArrayType> global_space_derivatives;
                SizeType derivative_order = 1;
                mrSlaveSkinModelPart->pGetGeometry(best_slave_curve_id)->GlobalSpaceDerivatives(global_space_derivatives, best_curve_projection_local_coord, derivative_order);
                CoordinatesArrayType tangent_vector = global_space_derivatives[1];
                double tangent_magnitude = norm_2(tangent_vector);
                tangent_vector /= tangent_magnitude;
                Vector normal_vector = ZeroVector(3);
                normal_vector[0] = tangent_vector[1];
                normal_vector[1] = -tangent_vector[0];

                new_slave_skin_node->SetValue(NORMAL, normal_vector);
                new_slave_skin_node->SetValue(LOCAL_TANGENT, tangent_vector);


                // find the closest slave gauss point to the projected skin node on the slave boundary
                PointerType pointToSearch = PointerType(new PointType(1, new_slave_skin_node->Coordinates()));

                int obtainedResults = testBins.SearchInRadius(*pointToSearch, radius, Results.begin(), list_of_distances.begin(), numberOfResults);

                double minimum_distance=1e10;
                int nearestNodeId;
                for (int i_distance = 0; i_distance < obtainedResults; i_distance++) {
                    double new_distance = list_of_distances[i_distance];   
                    if (new_distance < minimum_distance) { 
                        minimum_distance = new_distance;
                        nearestNodeId = i_distance;
                        }
                }
                
                KRATOS_ERROR_IF(obtainedResults == 0) << "::[IgaContactProcessSbm]:: Zero points found in search!!" 
                    << *new_slave_skin_node <<  best_skin_node << std::endl;

                IndexType best_brep_surrogate_slave_id = Results[nearestNodeId]->Id();
                auto best_slave_brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(mrSlaveModelPart->pGetGeometry(best_brep_surrogate_slave_id));
                KRATOS_ERROR_IF(!best_slave_brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << best_brep_surrogate_slave_id << "is not a Brep." << std::endl;
                // find the closest integration point in the brep
                GeometriesArrayType gp_list_slave;
                best_slave_brep_geometry->GetQuadraturePointGeometries(gp_list_slave);
                double best_distance_skin_surrogate_slave = 1e16;
                IndexType best_surrogate_node_array_position;
                IndexType count_gp_slave = 0;
                for (auto& gp_in_brep : gp_list_slave) {
                    double curr_distance = norm_2(gp_in_brep.Center() - *new_slave_skin_node);

                    if (curr_distance < best_distance_skin_surrogate_slave)
                    {
                        best_distance_skin_surrogate_slave = curr_distance;
                        best_surrogate_node_array_position = count_gp_slave;
                    }
                    count_gp_slave++;
                }
                new_slave_skin_node->SetValue(NEIGHBOUR_GEOMETRY, gp_list_slave(best_surrogate_node_array_position));

                rResultGeometries(count_gp_master) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCouplingGeometry2D(
                                                        gp_in_brep, gp_list_slave(best_surrogate_node_array_position));

                
                gp_in_brep->SetValue(IDENTIFIER, "active");
                gp_list_slave(best_surrogate_node_array_position)->SetValue(IDENTIFIER, "active");


                rResultGeometries(count_gp_master)->SetValue(MARKER_MESHES, master_knot_step_uv);
                
                count_gp_master++;
            }
            rResultGeometries.resize(count_gp_master);
            if (is_converged_at_least_once_in_brep)
            {
                this->CreateConditions(
                                rResultGeometries.ptr_begin(), rResultGeometries.ptr_end(),
                                *mrContactModelPart, name, id, mpPropMaster, mpPropSlave);
            }
        }

        // KRATOS_ERROR_IF_NOT(is_surrogate_correct) << "SURROGATE NOT CORRECT FOR CONTACT SBM PROBLEM" << std::endl;

        // Set non active slave breps to SBM load condition to zero
        std::string default_condition_name = "SBMLoadSolid2DCondition";
        SizeType rIdCounter = 1;
        if (mrContactModelPart->GetRootModelPart().Conditions().size() > 0)
            rIdCounter = mrContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;

        for (auto& r_geometry_slave : mrSlaveModelPart->Geometries()) {

            int slave_brep_id = r_geometry_slave.Id();
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);

            auto slave_brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);

            KRATOS_ERROR_IF(!slave_brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << slave_brep_id << "is not a Brep." << std::endl;

            GeometriesArrayType gp_list_slave;
            slave_brep_geometry->GetQuadraturePointGeometries(gp_list_slave);

            // CREATE SBM LOAD CONDITIONS 
            const Condition& rReferenceCondition = KratosComponents<Condition>::Get(default_condition_name);
            ModelPart::ConditionsContainerType new_condition_list;

            KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
                << "Creating conditions of type " << default_condition_name
                << " in " << mrContactModelPart->GetParentModelPart().Name() << "-SubModelPart." << std::endl;

            IndexType count_cond = 0;
            for (auto it = gp_list_slave.ptr_begin(); it != gp_list_slave.ptr_end(); ++it) {
                
                if ((*it)->GetValue(IDENTIFIER) == "active") continue;

                new_condition_list.push_back(
                    rReferenceCondition.Create(rIdCounter, (*it), mpPropSlave));
                
                new_condition_list.GetContainer()[count_cond]->SetValue(MARKER_MESHES, meshSizes_uv);
                                
                for (SizeType i = 0; i < (*it)->size(); ++i) {
                    // These are the control points associated with the basis functions involved in the condition we are creating
                    // rModelPart.AddNode((*it)->pGetPoint(i));
                    mrContactModelPart->GetParentModelPart().Nodes().push_back((*it)->pGetPoint(i));
                }
                rIdCounter++;
                count_cond++;
            }
            mrContactModelPart->GetParentModelPart().AddConditions(new_condition_list.begin(), new_condition_list.end());
        }
            //-----------------------------------------------------------------------------------------------------
            // master inactive breps
        if (mrContactModelPart->GetRootModelPart().Conditions().size() > 0)
            rIdCounter = mrContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;
        for (auto& r_geometry_master : mrMasterModelPart->Geometries()) {

            int master_brep_id = r_geometry_master.Id();
            auto p_master_geometry = mrMasterModelPart->pGetGeometry(master_brep_id);

            auto master_brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_master_geometry);

            KRATOS_ERROR_IF(!master_brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << master_brep_id << "is not a Brep." << std::endl;

            GeometriesArrayType gp_list_master;
            master_brep_geometry->GetQuadraturePointGeometries(gp_list_master);

            // CREATE SBM LOAD CONDITIONS 
            const Condition& rReferenceCondition = KratosComponents<Condition>::Get(default_condition_name);
            ModelPart::ConditionsContainerType new_condition_list;

            KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
                << "Creating conditions of type " << default_condition_name
                << " in " << mrContactModelPart->GetParentModelPart().Name() << "-SubModelPart." << std::endl;

            IndexType count_cond = 0;
            for (auto it = gp_list_master.ptr_begin(); it != gp_list_master.ptr_end(); ++it) {
                
                if ((*it)->GetValue(IDENTIFIER) == "active") continue;

                new_condition_list.push_back(rReferenceCondition.Create(rIdCounter, (*it), mpPropMaster));
                
                new_condition_list.GetContainer()[count_cond]->SetValue(MARKER_MESHES, meshSizes_uv);
                new_condition_list.GetContainer()[count_cond]->SetValue(IDENTIFIER, "outer");
                                
                for (SizeType i = 0; i < (*it)->size(); ++i) {
                    // These are the control points associated with the basis functions involved in the condition we are creating
                    // rModelPart.AddNode((*it)->pGetPoint(i));
                    mrContactModelPart->GetParentModelPart().Nodes().push_back((*it)->pGetPoint(i));
                }
                rIdCounter++;
                count_cond++;
            }

            mrContactModelPart->GetParentModelPart().AddConditions(new_condition_list.begin(), new_condition_list.end());
        }

        EntitiesUtilities::InitializeEntities<Condition>(mrContactModelPart->GetParentModelPart());

        KRATOS_ERROR_IF(mrContactModelPart->NumberOfConditions() == 0) << "YOUR CONTACT MODEL PART IS EMPTY" << std::endl;
    }
    



     void IgaContactProcessSbm::CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pPropMaster,
        PropertiesPointerType pPropSlave) const
    {
        // const SupportContact2DCondition rReferenceCondition = SupportContact2DCondition();
        // const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;

        KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
            << "Creating conditions of type " << rConditionName
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_condition_list.push_back(
                Kratos::make_intrusive<SbmContact2DCondition>(
                rIdCounter, (*it), pPropMaster, pPropSlave));
            
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                // rModelPart.AddNode((*it)->pGetPoint(i));
                rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            }
            rIdCounter++;
        }
        
        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());

        EntitiesUtilities::InitializeEntities<Condition>(rModelPart);

    }

    void IgaContactProcessSbm::GetDeformedPosition(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& SparseBrepMatrix,
        CoordinatesArrayType& rPointDeformedCoordinates)
    {
        const double epsilon = 1e-12;
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(MARKER_MESHES);
        Vector parameter_space_extremes = rModelPart.GetParentModelPart().GetValue(LOAD_MESHES);

        // compute the knot span where the point lies
        int knot_span_surrounded_u_id = floor((rPointGlobalCoordinates[0] - parameter_space_extremes[0]+epsilon)/knot_step_uv[0]);
        int knot_span_surrounded_v_id = floor((rPointGlobalCoordinates[1] - parameter_space_extremes[1]+epsilon)/knot_step_uv[1]);

        // check where the surrogate brep is in the boundary knot span

        IndexType brep_id = 0;
        bool index_found = false;
        for (int i = 0; i < 2; i++)
        {
            if (knot_span_surrounded_u_id+i < 0) continue;
            for (int j = 0; j < 2; j++)
            {
                if (knot_span_surrounded_v_id+j < 0) continue;
                if (SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j) != 0)
                {
                    brep_id = SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j);
                    index_found = true;
                    break;
                }
            }
            if (index_found) break;
        }
        int first_i = 0; int first_j = 0;
        while (!index_found)
        {
            first_i--; first_j--;
            for (int i = first_i; i < 2; i++)
            {
                if (knot_span_surrounded_u_id+i < 0) continue;
                for (int j = first_j; j < 2; j++)
                {
                    if (knot_span_surrounded_v_id+j < 0) continue;
                    if (SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j) != 0)
                    {
                        brep_id = SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j);
                        index_found = true;
                        break;
                    }
                }
                if (index_found) break;
            }
        }
        KRATOS_ERROR_IF(brep_id == 0) << "::[IgaContactProcessSbm]:: no boundary surrogate brep around true point: " <<
                                         rPointGlobalCoordinates << " found in knot span: (" << 
                                         knot_span_surrounded_u_id << ", " << knot_span_surrounded_v_id << ")" << std::endl;
        // retrieve the brep
        auto p_geometry = rModelPart.pGetGeometry(brep_id);
        auto brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_geometry);

        KRATOS_ERROR_IF(!brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << brep_id << "is not a Brep." << std::endl;

        GeometriesArrayType gp_list;
        brep_geometry->GetQuadraturePointGeometries(gp_list);

        //----------------------------
        // get displacement of the reference integration point on the surrogate
        auto& surrogate_point_geometry = gp_list(1);
        SizeType number_of_control_points = surrogate_point_geometry->size();
        const SizeType mat_size = number_of_control_points * 2;
        Vector control_point_displacement_coefficients = ZeroVector(mat_size);
        GetValuesVector(*surrogate_point_geometry, control_point_displacement_coefficients);


        const Matrix& H = surrogate_point_geometry->ShapeFunctionsValues();
        std::vector<Matrix> n_shape_function_derivatives; 

        const GeometryType::ShapeFunctionsGradientsType& DN_De = surrogate_point_geometry->ShapeFunctionsLocalGradients();
        int basis_functions_order = std::sqrt(DN_De[0].size1()) - 1;
        for (int n = 1; n <= basis_functions_order; n++) {
            n_shape_function_derivatives.push_back(surrogate_point_geometry->ShapeFunctionDerivatives(n, 0));
        }

        // compute displacement extension matrix
        Vector distance_vector = rPointGlobalCoordinates - surrogate_point_geometry->Center();
        Matrix H_sum = ZeroMatrix(1, number_of_control_points);
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            double H_taylor_term = 0.0; // Reset for each node
            for (int n = 1; n <= basis_functions_order; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = n_shape_function_derivatives[n-1];
                for (int k = 0; k <= n; k++) {
                    int n_k = n - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term += ComputeTaylorTerm(derivative, distance_vector[0], n_k, distance_vector[1], k);
                }
            }
            H_sum(0,i) = H_taylor_term + H(0, i);;
        }
        //----------------------------------
        // extend the value to the true boundary
        // reset matrix for matrix to vector product
        const SizeType dim = 2;
        Matrix H_projection_matrix = ZeroMatrix(dim, dim*number_of_control_points);
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            for (IndexType i_dim = 0; i_dim < dim; i_dim++) {
                H_projection_matrix(i_dim, dim*i+i_dim) = H_sum(0, i);
            }
        }

        Vector point_displacement = prod(H_projection_matrix,control_point_displacement_coefficients);

        rPointDeformedCoordinates[0] = rPointGlobalCoordinates[0] + point_displacement[0];
        rPointDeformedCoordinates[1] = rPointGlobalCoordinates[1] + point_displacement[1];
    }



    void IgaContactProcessSbm::GetDeformedGradient(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& SparseBrepMatrix,
        Matrix& rPointGradientDeformation,
        Matrix& rPointHessianDeformation)
    {
        const double epsilon = 1e-12;
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(MARKER_MESHES);
        Vector parameter_space_extremes = rModelPart.GetParentModelPart().GetValue(LOAD_MESHES);

        // compute the knot span where the point lies
        int knot_span_surrounded_u_id = floor((rPointGlobalCoordinates[0] - parameter_space_extremes[0]+epsilon)/knot_step_uv[0]);
        int knot_span_surrounded_v_id = floor((rPointGlobalCoordinates[1] - parameter_space_extremes[1]+epsilon)/knot_step_uv[1]);

        // check where the surrogate brep is in the boundary knot span
        IndexType brep_id = 0;
        bool index_found = false;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j) != 0)
                {
                    brep_id = SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j);
                    index_found = true;
                    break;
                }
            }
            if (index_found) break;
        }
        int first_i = 0; int first_j = 0;
        while (!index_found)
        {
            first_i--; first_j--;
            for (int i = first_i; i < 2; i++)
            {
                if (knot_span_surrounded_u_id+i < 0) continue;
                for (int j = first_j; j < 2; j++)
                {
                    if (knot_span_surrounded_v_id+j < 0) continue;
                    if (SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j) != 0)
                    {
                        brep_id = SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j);
                        index_found = true;
                        break;
                    }
                }
                if (index_found) break;
            }
        }
        
        KRATOS_ERROR_IF(brep_id == 0) << "::[IgaContactProcessSbm]:: no boundary surrogate brep around true point: " <<
                                         rPointGlobalCoordinates << " found in knot span: (" << 
                                         knot_span_surrounded_u_id << ", " << knot_span_surrounded_v_id << ")" << std::endl;
        // retrieve the brep
        auto p_geometry = rModelPart.pGetGeometry(brep_id);
        auto brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_geometry);

        KRATOS_ERROR_IF(!brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << brep_id << "is not a Brep." << std::endl;

        GeometriesArrayType gp_list;
        brep_geometry->GetQuadraturePointGeometries(gp_list);

        //----------------------------
        // get displacement of the reference integration point on the surrogate
        auto& surrogate_point_geometry = gp_list(1);
        SizeType number_of_control_points = surrogate_point_geometry->size();
        const SizeType mat_size = number_of_control_points * 2;
        Vector control_point_displacement_coefficients = ZeroVector(mat_size);
        GetValuesVector(*surrogate_point_geometry, control_point_displacement_coefficients);


        const Matrix& H = surrogate_point_geometry->ShapeFunctionsValues();
        std::vector<Matrix> n_shape_function_derivatives; 

        const GeometryType::ShapeFunctionsGradientsType& DN_De = surrogate_point_geometry->ShapeFunctionsLocalGradients();

        int basis_functions_order = std::sqrt(DN_De[0].size1()) - 1;
        for (int n = 1; n <= basis_functions_order; n++) {
            n_shape_function_derivatives.push_back(surrogate_point_geometry->ShapeFunctionDerivatives(n, 0));
        }

        // compute displacement gradient extension matrix
        Vector distance_vector = rPointGlobalCoordinates - surrogate_point_geometry->Center();
        Matrix H_grad = ZeroMatrix(number_of_control_points, 2);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            double H_taylor_term_X = 0.0; // Reset for each node
            double H_taylor_term_Y = 0.0; 
            for (int n = 2; n <= basis_functions_order; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = n_shape_function_derivatives[n-1];
                for (int k = 0; k <= n-1; k++) {
                    int n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_X += ComputeTaylorTerm(derivative, distance_vector[0], n_k, distance_vector[1], k);
                }
                for (int k = 0; k <= n-1; k++) {
                    int n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, distance_vector[0], n_k, distance_vector[1], k);
                }
            }
            H_grad(i,0) = H_taylor_term_X + n_shape_function_derivatives[0](i,0);
            H_grad(i,1) = H_taylor_term_Y + n_shape_function_derivatives[0](i,1);
        }         
        //----------------------------------
        // extend the value to the true boundary
        // reset matrix for matrix to vector product
        const SizeType dim = 2;
        Matrix H_grad_u_projection_matrix = ZeroMatrix(2, dim*number_of_control_points);
        Matrix H_grad_v_projection_matrix = ZeroMatrix(2, dim*number_of_control_points);
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            H_grad_u_projection_matrix(0, dim*i) = H_grad(i,0);
            H_grad_u_projection_matrix(1, dim*i) = H_grad(i,1);

            H_grad_v_projection_matrix(0, dim*i+1) = H_grad(i,0);
            H_grad_v_projection_matrix(1, dim*i+1) = H_grad(i,1);
        }

        Vector point_gradient_deformation_u = prod(H_grad_u_projection_matrix,control_point_displacement_coefficients);
        Vector point_gradient_deformation_v = prod(H_grad_v_projection_matrix,control_point_displacement_coefficients);

        if (rPointGradientDeformation.size1() != dim || rPointGradientDeformation.size2() != dim)
        {
            rPointGradientDeformation.resize(dim, dim);
        }
        noalias(row(rPointGradientDeformation, 0)) = point_gradient_deformation_u;
        noalias(row(rPointGradientDeformation, 1)) = point_gradient_deformation_v;
        //-------------------------------------------------------------------------
        // COMPUTE HESSIAN EXTENSION MATRIX
        Matrix H_hess = ZeroMatrix(number_of_control_points, 3);
        if (basis_functions_order >= 2)
        {
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            double H_taylor_term_X = 0.0; // Reset for each node
            double H_taylor_term_XY = 0.0; 
            double H_taylor_term_Y = 0.0; 
            for (int n = 3; n <= basis_functions_order; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = n_shape_function_derivatives[n-1];
                for (int k = 0; k <= n-2; k++) {
                    int n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_X += ComputeTaylorTerm(derivative, distance_vector[0], n_k, distance_vector[1], k);
                }
                for (int k = 0; k <= n-2; k++) {
                    int n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_XY += ComputeTaylorTerm(derivative, distance_vector[0], n_k, distance_vector[1], k);
                }
                for (int k = 0; k <= n-2; k++) {
                    int n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+2); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, distance_vector[0], n_k, distance_vector[1], k);
                }
            }
            H_hess(i,0) = H_taylor_term_X + n_shape_function_derivatives[1](i,0);
            H_hess(i,1) = H_taylor_term_XY + n_shape_function_derivatives[1](i,1);
            H_hess(i,2) = H_taylor_term_Y + n_shape_function_derivatives[1](i,2);
        }     
        }  

        // extend the value to the true boundary

        Matrix H_hess_u_projection_matrix = ZeroMatrix(3, dim*number_of_control_points);
        Matrix H_hess_v_projection_matrix = ZeroMatrix(3, dim*number_of_control_points);
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            H_hess_u_projection_matrix(0, dim*i) = H_hess(i,0);
            H_hess_u_projection_matrix(1, dim*i) = H_hess(i,1);
            H_hess_u_projection_matrix(2, dim*i) = H_hess(i,2);
            
            H_hess_v_projection_matrix(0, dim*i+1) = H_hess(i,0);
            H_hess_v_projection_matrix(1, dim*i+1) = H_hess(i,1);
            H_hess_v_projection_matrix(2, dim*i+1) = H_hess(i,2);
        }

        Vector point_hessian_deformation_u = prod(H_hess_u_projection_matrix,control_point_displacement_coefficients);
        Vector point_hessian_deformation_v = prod(H_hess_v_projection_matrix,control_point_displacement_coefficients);
        if (rPointHessianDeformation.size1() != dim || rPointHessianDeformation.size2() != 3)
        {
            rPointHessianDeformation.resize(dim, 3);
        }

        noalias(row(rPointHessianDeformation, 0)) = point_hessian_deformation_u;
        noalias(row(rPointHessianDeformation, 1)) = point_hessian_deformation_v;
    }

    void IgaContactProcessSbm::GetValuesVector(
        const Geometry<PointType>& rGeometry,
        Vector& rValues)
    {
        const SizeType number_of_control_points = rGeometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = rGeometry(i)->GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }



    bool IgaContactProcessSbm::NewtonRaphsonCurveOnDeformed(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinates, // already deformed
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const NurbsCurveGeometryType& rPairedGeometry, 
        double& distance,
        const int rNumberOfInitialGuesses,
        const int MaxIterations,
        const double Accuracy)
    {
        // Intialize variables
        double residual = 0.0, delta_t = 0.0;
        distance = 0.0;

        std::vector<array_1d<double, 3>> curve_derivatives(3, ZeroVector(3));
        array_1d<double, 3> distance_vector = ZeroVector(3);
        
        Vector gradient_derivatives_updated = ZeroVector(3);
        Vector hessian_derivatives_updated = ZeroVector(3);

        CoordinatesArrayType projected_point_deformed_global_coordinates = ZeroVector(3);
        Matrix projected_point_gradient_deformation = ZeroMatrix(2, 2);
        Matrix projected_point_hessian_deformation = ZeroMatrix(2, 3);

        CoordinatesArrayType current_point_global_coordinates = ZeroVector(3);

        Vector curve_interval(2);
        auto interval = rPairedGeometry.DomainInterval();
        curve_interval[0] = interval.GetT0();
        curve_interval[1] = interval.GetT1();
        bool is_converged_at_least_once = false;
        double best_distance = 1e12;

        for (IndexType i_guess = 0; i_guess < rNumberOfInitialGuesses; i_guess++)
        {   
            residual = Accuracy + 1;
            CoordinatesArrayType t = ZeroVector(3);
            t[0] = curve_interval[0] +  (curve_interval[1]- curve_interval[0])/(rNumberOfInitialGuesses-1) * float(i_guess);

            // Loop over all Newton-Raphson iterations
            for (int i = 0; i < MaxIterations; ++i)
            {
                // Compute the position, the base and the acceleration vector in the curve space
                rPairedGeometry.GlobalSpaceDerivatives(
                    curve_derivatives,
                    t,
                    2);
                    
                current_point_global_coordinates = curve_derivatives[0]; // undeformed

                // deformed position

                if (current_point_global_coordinates[0] < 0)
                {
                    KRATOS_WATCH(current_point_global_coordinates)
                    exit(0);
                }
                GetDeformedPosition(current_point_global_coordinates, *mrSlaveModelPart, mSparseBrepMatrixSlave, projected_point_deformed_global_coordinates);

                // // NEW
                GetDeformedGradient(current_point_global_coordinates, *mrSlaveModelPart, mSparseBrepMatrixSlave, projected_point_gradient_deformation, projected_point_hessian_deformation);

                for (int i_dim = 0; i_dim < 2; i_dim++) {
                    gradient_derivatives_updated[i_dim] = curve_derivatives[1][i_dim] +
                                                          projected_point_gradient_deformation(i_dim,0)*curve_derivatives[1][0] +
                                                          projected_point_gradient_deformation(i_dim,1)*curve_derivatives[1][1];
                }

                for (int i_dim = 0; i_dim < 2; i_dim++) {
                    hessian_derivatives_updated[i_dim] = curve_derivatives[2][i_dim] +
                                                          projected_point_hessian_deformation(i_dim,0)*curve_derivatives[1][0]*curve_derivatives[1][0] +
                                                          projected_point_hessian_deformation(i_dim,2)*curve_derivatives[1][1]*curve_derivatives[1][1] +
                                                          projected_point_hessian_deformation(i_dim,0)*curve_derivatives[2][0] +
                                                          projected_point_hessian_deformation(i_dim,2)*curve_derivatives[2][2] +
                                                          2*projected_point_hessian_deformation(i_dim,1)*curve_derivatives[1][0]*curve_derivatives[1][1];
                }


                // Compute the distance vector between the point and its
                // projection on the curve
                distance_vector = projected_point_deformed_global_coordinates - rPointGlobalCoordinates;
                distance = norm_2(distance_vector);

                if (distance < Accuracy) // Acc
                {
                    rProjectedPointLocalCoordinates = t;
                    rProjectedPointGlobalCoordinates = current_point_global_coordinates;
                    return true;
                }

                // Compute the residual
                residual = inner_prod(distance_vector, gradient_derivatives_updated);
                if (std::abs(residual) < Accuracy) // Acc
                {
                    if (std::isnan(residual)) break;
                    else 
                    {
                        is_converged_at_least_once = true;

                        if (distance < best_distance) {
                            best_distance = distance;
                            rProjectedPointGlobalCoordinates = current_point_global_coordinates;
                            rProjectedPointLocalCoordinates = t;
                        }
                        break;
                    }
                }

                // Compute the increment
                double denom = inner_prod(hessian_derivatives_updated, distance_vector) + pow(norm_2(gradient_derivatives_updated), 2);
                if (std::abs(denom) < 1e-10) denom = 1e-10;  // Avoid division by zero
                delta_t = residual / denom;
                
                // Increment the parametric coordinate
                t[0] -= delta_t;

                // Check if the increment is too small and if yes return true
                if (norm_2(delta_t *gradient_derivatives_updated) < Accuracy) // Acc
                {
                    double to_check = norm_2(delta_t * gradient_derivatives_updated);
                    if (std::isnan(to_check)) break;
                    else 
                    {
                        if (distance < best_distance) {
                            best_distance = distance;
                            rProjectedPointLocalCoordinates = t;

                            rPairedGeometry.GlobalSpaceDerivatives(
                                                                    curve_derivatives,
                                                                    t,
                                                                    2);
                            rProjectedPointGlobalCoordinates = curve_derivatives[0];
                        }
                        break;
                    }
                }
                    

                // Check if the parameter gets out of its interval of definition and if so clamp it
                // back to the boundaries
                // CoordinatesArrayType closest_t = ZeroVector(3);
                // int check = rPairedGeometry.ClosestPointLocalToLocalSpace(
                //     t, t,2);

                int check = 0;
                double curve_interval_min; double curve_interval_max; 
                // KRATOS_WATCH(t)
                if (curve_interval[0] < curve_interval[1]) {
                    curve_interval_min = curve_interval[0];
                    curve_interval_max = curve_interval[1];
                }
                else {
                    curve_interval_min = curve_interval[1];
                    curve_interval_max = curve_interval[0];
                }
                //----
                if (t[0] > curve_interval_min && t[0] < curve_interval_max) check = 1;
                else if (std::abs(t[0] - curve_interval_min) < 1e-7) {
                    t[0] = curve_interval_min; 
                    check = 2;
                }
                else if (std::abs(t[0] - curve_interval_max) < 1e-7) {
                    t[0] = curve_interval_max; 
                    check = 2;
                }

                if (check == 0) {
                    // if (projection_reset_to_boundary) { return false; }
                    // else { projection_reset_to_boundary = true; }

                    break;
                }
            }
        }

        if (is_converged_at_least_once) 
        {
            distance = best_distance;
            return true;
        }
        // Return false if the Newton-Raphson iterations did not converge
        
        return false;
    }


    // Function to compute a single term in the Taylor expansion
    double IgaContactProcessSbm::ComputeTaylorTerm(double derivative, double dx, int n_k, double dy, int k)
    {
        return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (Factorial(k) * Factorial(n_k));    
    }

    unsigned long long IgaContactProcessSbm::Factorial(int n) 
    {
        if (n == 0) return 1;
        unsigned long long result = 1;
        for (int i = 2; i <= n; ++i) result *= i;
        return result;
    }

} // End namespace Kratos
