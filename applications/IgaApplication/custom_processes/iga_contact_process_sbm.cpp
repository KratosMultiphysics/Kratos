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
        Vector master_knot_step_uv = mrMasterModelPart->GetParentModelPart().GetValue(MARKER_MESHES)/2; // need half the size of the edge
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
        Vector slave_knot_step_uv = mrSlaveModelPart->GetParentModelPart().GetValue(MARKER_MESHES)/2; // need half the size of the edge
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
   }

    void IgaContactProcessSbm::Execute(){

        std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::trunc);

        std::string master_layer_name = mParameters["contact_parameters"]["master_model_part"]["layer_name"].GetString();
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
        

        Vector master_knot_step_uv = mrMasterModelPart->GetParentModelPart().GetValue(MARKER_MESHES);

        // FIXME: new mortar implementation

        // Map all brep vertices from slave to master with backward projection
        
        // 1) retrieve list of slave vertices to project back
        const double epsilon = 1e-12;
        Vector slave_knot_step_uv = mrSlaveModelPart->GetParentModelPart().GetValue(MARKER_MESHES)/2; // need half the size of the edge
        Vector slave_parameter_space_extremes = mrSlaveModelPart->GetParentModelPart().GetValue(LOAD_MESHES);
        double slave_domain_length_u = slave_parameter_space_extremes[2] - slave_parameter_space_extremes[0];
        double slave_domain_length_v = slave_parameter_space_extremes[3] - slave_parameter_space_extremes[1];

        const SizeType slave_number_of_knot_spans_u = floor((slave_domain_length_u + epsilon)/slave_knot_step_uv[0]);
        const SizeType slave_number_of_knot_spans_v = floor((slave_domain_length_v + epsilon)/slave_knot_step_uv[1]);
        UblasSpace<int, CompressedMatrix, Vector>::MatrixType is_slave_vertex_projected(slave_number_of_knot_spans_u+1, slave_number_of_knot_spans_v+1,false);

        // Lista di Point
        std::vector<CoordinatesArrayType> brep_slave_vertices;
        std::vector<array_1d<IndexType, 2>> brep_id_of_vertices;
        IndexType count_vertex = 0;
        for (auto& r_geometry_slave : mrSlaveModelPart->Geometries()) {

            int slave_brep_id = r_geometry_slave.Id();
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);
            p_slave_geometry->SetValue(ACTIVATION_LEVEL, 0.0);

            // compute the vertex of each surrogate brep
            Vector brep_domain_interval;
            p_slave_geometry->DomainInterval(brep_domain_interval);
            CoordinatesArrayType brep_vertex_1 = ZeroVector(3); CoordinatesArrayType brep_vertex_1_local_coords = ZeroVector(3);
            CoordinatesArrayType brep_vertex_2 = ZeroVector(3); CoordinatesArrayType brep_vertex_2_local_coords = ZeroVector(3);
            brep_vertex_1_local_coords[0] = brep_domain_interval[0];
            brep_vertex_2_local_coords[0] = brep_domain_interval[1];

            p_slave_geometry->GlobalCoordinates(brep_vertex_1, brep_vertex_1_local_coords);
            p_slave_geometry->GlobalCoordinates(brep_vertex_2, brep_vertex_2_local_coords);

            // check if vertex already projected
            IndexType knot_span_u_1 = floor((brep_vertex_1[0]- slave_parameter_space_extremes[0]+epsilon)/slave_knot_step_uv[0]);
            IndexType knot_span_v_1 = floor((brep_vertex_1[1]- slave_parameter_space_extremes[1]+epsilon)/slave_knot_step_uv[1]);

            IndexType knot_span_u_2 = floor((brep_vertex_2[0]- slave_parameter_space_extremes[0]+epsilon)/slave_knot_step_uv[0]);
            IndexType knot_span_v_2 = floor((brep_vertex_2[1]- slave_parameter_space_extremes[1]+epsilon)/slave_knot_step_uv[1]);

            if (is_slave_vertex_projected(knot_span_u_1, knot_span_v_1) == 0) {
                is_slave_vertex_projected(knot_span_u_1, knot_span_v_1) = count_vertex+1;
                //add brep_vertex_1 to point list
                brep_slave_vertices.push_back(brep_vertex_1);

                array_1d<IndexType, 2> vertex_brep_ids;
                vertex_brep_ids[0] = slave_brep_id;
                vertex_brep_ids[1] = 0;
                brep_id_of_vertices.push_back(vertex_brep_ids);
                count_vertex++;
            }
            else {
                int vertex_id = is_slave_vertex_projected(knot_span_u_1, knot_span_v_1)-1;
                brep_id_of_vertices[vertex_id][1] = slave_brep_id;
            }
            if (is_slave_vertex_projected(knot_span_u_2, knot_span_v_2) == 0) {
                is_slave_vertex_projected(knot_span_u_2, knot_span_v_2) = count_vertex+1;
                //add brep_vertex_2 to point list
                brep_slave_vertices.push_back(brep_vertex_2);

                array_1d<IndexType, 2> vertex_brep_ids;
                vertex_brep_ids[0] = slave_brep_id;
                vertex_brep_ids[1] = 0;
                brep_id_of_vertices.push_back(vertex_brep_ids);
                count_vertex++;
            }
            else {
                int vertex_id = is_slave_vertex_projected(knot_span_u_2, knot_span_v_2)-1;
                
                brep_id_of_vertices[vertex_id][1] = slave_brep_id;
            }
        }

        // 2) project the brep slave vertices to the slave true boundary

        std::vector<CoordinatesArrayType> projected_slave_vertices_on_true;
        std::vector<array_1d<double, 3>> first_derivatives_projected_slave_vertices_on_true;

        for (const auto& vertex : brep_slave_vertices) {
            std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));
            CoordinatesArrayType best_projected_vertex;
            
            std::string projection_layer_name = "";
            ProjectToSkinBoundary(mrSlaveSkinModelPart, projection_layer_name, vertex, best_projected_vertex, curve_derivatives, 10); // FIXME:

            projected_slave_vertices_on_true.push_back(best_projected_vertex);
            first_derivatives_projected_slave_vertices_on_true.push_back(curve_derivatives[1]);
        }

        // KRATOS_WATCH(brep_slave_vertices)
        // KRATOS_WATCH(projected_slave_vertices_on_true)

        // 3) project the slave_true_points to the master true boundary via ray tracing
        count_vertex = 0;
        std::vector<CoordinatesArrayType> master_projections_on_true_of_slave_vertices;
        std::vector<Vector> normal_projections_on_master_of_slave_vertices;
        for (auto& r_geometry : mrSlaveModelPart->Geometries()) {
            r_geometry.SetValue(ACTIVATION_LEVEL, 0);
            r_geometry.SetValue(IDENTIFIER, "");
        }
        for (const auto& skin_vertex_slave : projected_slave_vertices_on_true) {
            
            // get the deformed position of the slave vertex and the normal
            CoordinatesArrayType skin_vertex_slave_deformed;
            Matrix skin_vertex_gradient_deformation = ZeroMatrix(2, 2);
            Matrix skin_vertex_hessian_deformation = ZeroMatrix(2, 3);

            Vector slave_tangent_on_vertex_deformed = ZeroVector(3);

            GetDeformedPosition(skin_vertex_slave, *mrSlaveModelPart, mSparseBrepMatrixSlave, skin_vertex_slave_deformed);

            GetDeformedGradient(skin_vertex_slave, *mrSlaveModelPart, mSparseBrepMatrixSlave, skin_vertex_gradient_deformation, skin_vertex_hessian_deformation);

            for (int i_dim = 0; i_dim < 2; i_dim++) {
                slave_tangent_on_vertex_deformed[i_dim] = first_derivatives_projected_slave_vertices_on_true[count_vertex][i_dim] +
                                                        skin_vertex_gradient_deformation(i_dim,0)*first_derivatives_projected_slave_vertices_on_true[count_vertex][0] +
                                                        skin_vertex_gradient_deformation(i_dim,1)*first_derivatives_projected_slave_vertices_on_true[count_vertex][1];
            }

            slave_tangent_on_vertex_deformed /= norm_2(slave_tangent_on_vertex_deformed);
            CoordinatesArrayType slave_normal_on_vertex_deformed = ZeroVector(3);
            slave_normal_on_vertex_deformed[0] = -slave_tangent_on_vertex_deformed[1];
            slave_normal_on_vertex_deformed[1] = slave_tangent_on_vertex_deformed[0];

            // CHATGPT CODE
            CoordinatesArrayType best_projected_point_master = ZeroVector(3);
            Vector normal_vector = ZeroVector(3);
            double best_distance_to_master = 1e12;

            bool has_converged_at_least_once = false;

            for (auto& i_master_curve : mrMasterSkinModelPart->Geometries()) {
                if (i_master_curve.GetValue(IDENTIFIER) != master_layer_name) continue;

                int master_id = i_master_curve.Id();
                auto p_master_geom = mrMasterSkinModelPart->pGetGeometry(master_id);
                auto master_curve = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_master_geom);
                KRATOS_ERROR_IF(!master_curve) << ":::[IgaContactProcessSbm]::: Curve id " << master_id << " is not a NurbsCurveGeometryType." << std::endl;

                CoordinatesArrayType projected_point;
                CoordinatesArrayType projected_point_local_coords = ZeroVector(3);
                double projection_distance;
                bool converged = ProjectPointViaRayTracingToMasterCurveDeformed(
                    skin_vertex_slave_deformed,
                    slave_normal_on_vertex_deformed,
                    *master_curve,
                    *mrMasterModelPart,
                    mSparseBrepMatrixMaster,
                    projected_point,
                    projected_point_local_coords,
                    projection_distance);

                if (converged && projection_distance < best_distance_to_master && projection_distance < 2e0) {
                    has_converged_at_least_once = true;
                    best_distance_to_master = projection_distance;
                    best_projected_point_master = projected_point;


                    // compute normals 
                    std::vector<CoordinatesArrayType> global_space_derivatives;
                    SizeType derivative_order = 2;
                    master_curve->GlobalSpaceDerivatives(global_space_derivatives, projected_point_local_coords, derivative_order);
                    CoordinatesArrayType tangent_vector = global_space_derivatives[1];
                    double tangent_magnitude = norm_2(tangent_vector);
                    tangent_vector /= tangent_magnitude;
                    
                    normal_vector[0] = tangent_vector[1];
                    normal_vector[1] = -tangent_vector[0];
                }
            }

            if (has_converged_at_least_once)
            {
                master_projections_on_true_of_slave_vertices.push_back(best_projected_point_master);
                normal_projections_on_master_of_slave_vertices.push_back(normal_vector);

                for (int k = 0; k < 2; k++)
                {
                    IndexType slave_brep_id = brep_id_of_vertices[count_vertex][k];
                    if (slave_brep_id == 0) continue;
                    auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);
                    int activation_level = p_slave_geometry->GetValue(ACTIVATION_LEVEL) + 1.0;
                    p_slave_geometry->SetValue(ACTIVATION_LEVEL, activation_level);

                    // retrieve vertices another time
                    if (p_slave_geometry->GetValue(IDENTIFIER) == "last_vertex_activated" || p_slave_geometry->GetValue(IDENTIFIER) == "first_vertex_activated") {
                        p_slave_geometry->SetValue(IDENTIFIER, "completely_activated");
                    }
                    else 
                    {
                        Vector brep_domain_interval;
                        p_slave_geometry->DomainInterval(brep_domain_interval);
                        CoordinatesArrayType brep_vertex_1 = ZeroVector(3); CoordinatesArrayType brep_vertex_1_local_coords = ZeroVector(3);
                        CoordinatesArrayType brep_vertex_2 = ZeroVector(3); CoordinatesArrayType brep_vertex_2_local_coords = ZeroVector(3);
                        brep_vertex_1_local_coords[0] = brep_domain_interval[0];
                        brep_vertex_2_local_coords[0] = brep_domain_interval[1];

                        p_slave_geometry->GlobalCoordinates(brep_vertex_1, brep_vertex_1_local_coords);
                        p_slave_geometry->GlobalCoordinates(brep_vertex_2, brep_vertex_2_local_coords);
                        
                        if (norm_2(brep_slave_vertices[count_vertex] - brep_vertex_1) < epsilon)
                        {
                            p_slave_geometry->SetValue(IDENTIFIER, "first_vertex_activated");
                        }   
                        else if (norm_2(brep_slave_vertices[count_vertex] - brep_vertex_2) < epsilon)
                        {
                            p_slave_geometry->SetValue(IDENTIFIER, "last_vertex_activated");
                        }
                        else
                        {
                            KRATOS_WATCH(brep_slave_vertices[count_vertex])
                            KRATOS_WATCH(brep_vertex_1)
                            KRATOS_WATCH(brep_vertex_2)
                            KRATOS_ERROR << "to better check " << std::endl;
                        }
                    }
                }
            }
            // Usa best_projection_point_on_master come risultato della proiezione
            
            count_vertex++;
        }

        // KRATOS_WATCH(master_projections_on_true_of_slave_vertices)
        // KRATOS_WATCH(normal_projections_on_master_of_slave_vertices)
        
        // 4) project the master_true_points to the master surrogate boundary via ray tracing
        count_vertex = 0;
        std::vector<IndexType> brep_id_backward_projections;
        std::vector<double> brep_local_coord_backward_projections;
        for (const auto& skin_vertex_master : master_projections_on_true_of_slave_vertices) {
            
            Vector normal_vector = normal_projections_on_master_of_slave_vertices[count_vertex];

            IndexType brep_id = -1;
            double projection_brep_local_variable = 0.0;

            ProjectBackToSurrogateBoundary(*mrMasterModelPart, skin_vertex_master, mSparseBrepMatrixMaster, normal_vector,
                                            brep_id, projection_brep_local_variable);
            
            // KRATOS_ERROR_IF(best_brep_id == -1) << "::[IgaContactProcessSbm]:: No brep found for the master skin node " 
            //                                 << skin_vertex_master << std::endl;

            if (brep_id == -1)
            {
                continue;
            }

            brep_local_coord_backward_projections.push_back(projection_brep_local_variable);
            brep_id_backward_projections.push_back(brep_id);

            auto p_geometry = mrMasterModelPart->pGetGeometry(brep_id);
            CoordinatesArrayType point_on_master = ZeroVector(3); CoordinatesArrayType local_coords = ZeroVector(3);
            local_coords[0] = projection_brep_local_variable;

            p_geometry->GlobalCoordinates(point_on_master, local_coords);
            count_vertex++;
        }
        
        // 5) create the contact model part (quadrature points from master + projections on slave)
        SizeType count_contact_gp_master = 0;
        SizeType count_neumann_gp_master = 0;

        double integral = 0;
        GeometriesArrayType result_geometries_contact;
        GeometriesArrayType result_geometries_neumann;

        std::vector<IndexType> brep_id_forward_projections;
        std::vector<double> brep_local_coord_forward_projections;
        
        for (auto& r_geometry_master : mrMasterModelPart->Geometries()) {
            int master_brep_id = r_geometry_master.Id();
            auto p_master_geometry = mrMasterModelPart->pGetGeometry(master_brep_id);
            auto master_brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_master_geometry);

            KRATOS_ERROR_IF(!master_brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << master_brep_id << "is not a Brep." << std::endl;

            std::vector<double> spans;
            master_brep_geometry->SpansLocalSpace(spans);

            for (IndexType i = 0; i < brep_local_coord_backward_projections.size(); i++) 
            {
                if (brep_id_backward_projections[i] == master_brep_id) {
                    spans.push_back(brep_local_coord_backward_projections[i]);
                }
            }

            std::sort(spans.begin(), spans.end());

            // remove duplicates with tolerance
            auto last = std::unique(spans.begin(), spans.end(), [](double a, double b) {
                return std::abs(a - b) < 1e-4;
            });
            spans.erase(last, spans.end());
            
            ///////////////////////////////////////////
            GeometriesArrayType geometries;
            SizeType shape_function_derivatives_order = 3;
            if (mParameters.Has("shape_function_derivatives_order")) {
                shape_function_derivatives_order = mParameters["shape_function_derivatives_order"].GetInt();
            }
            else {
                KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
                    << "shape_function_derivatives_order is not provided and thus being considered as 3. " << std::endl;
            }

            std::string quadrature_method = mParameters.Has("quadrature_method")
                ? mParameters["integration_rule"].GetString()
                : "GAUSS";
            IntegrationInfo integration_info = master_brep_geometry->GetDefaultIntegrationInfo();

            if (mParameters.Has("number_of_integration_points_per_span")) {
                for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                    integration_info.SetNumberOfIntegrationPointsPerSpan(i, mParameters["number_of_integration_points_per_span"].GetInt());
                }
            }
            
            for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                if (quadrature_method == "GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
                }
                else if (quadrature_method == "EXTENDED_GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
                }
                else if (quadrature_method == "GRID") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
                }
                else {
                    KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                        << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
                }
            }


            IntegrationPointsArrayType integration_points;
            IntegrationPointUtilities::CreateIntegrationPoints1D(integration_points, spans, integration_info);

            master_brep_geometry->CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_points, integration_info, false);
            
            // FIXME: delete
            // KRATOS_WATCH(geometries(0)->Center())
            // KRATOS_WATCH(geometries(0)->ShapeFunctionsLocalGradients())
            // GeometriesArrayType gp_list_master;
            // master_brep_geometry->GetQuadraturePointGeometries(gp_list_master);
            // KRATOS_WATCH(gp_list_master(0)->Center())
            // KRATOS_WATCH(gp_list_master(0)->ShapeFunctionsLocalGradients())

            // exit(0);
            // initialize the resulting geometries array
            result_geometries_contact.resize(count_contact_gp_master + geometries.size());
            result_geometries_neumann.resize(count_neumann_gp_master + geometries.size());

            bool project_vertex_to_slave = false;

            // for each master quadrature point compute the projection on the slave side
            for (IndexType i_gp_in_brep = 0; i_gp_in_brep < geometries.size(); i_gp_in_brep++) {
                auto& gp_in_brep = geometries(i_gp_in_brep);

                // KRATOS_WATCH(gp_in_brep->Center())
                
                // compute projection on skin master
                std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));
                CoordinatesArrayType gp_in_skin_master = ZeroVector(3);
                
                bool is_projected = ProjectToSkinBoundary(mrMasterSkinModelPart, master_layer_name, gp_in_brep->Center(), gp_in_skin_master, curve_derivatives, 10);
                
                if (!is_projected) {
                    // save it as Neumann quadrature point and continue
                    bool is_projected_to_skin_master = false;
                    std::string projection_layer_name = "";
                    is_projected_to_skin_master = ProjectToSkinBoundary(mrMasterSkinModelPart, projection_layer_name, gp_in_brep->Center(), gp_in_skin_master, curve_derivatives, 10);
                    

                    IndexType id_new_node_master = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size()+1;
                    auto new_master_skin_node = new Node(id_new_node_master, gp_in_skin_master);

                    NodePointerVector empty_vector;
                    empty_vector.push_back(new_master_skin_node); // Just it_node-plane neighbours
                    gp_in_brep->SetValue(NEIGHBOUR_NODES, empty_vector);

                    // compute normals and useful info
                    std::vector<CoordinatesArrayType> global_space_derivatives_master;
                    CoordinatesArrayType tangent_vector_master = curve_derivatives[1];
                    double tangent_magnitude_master = norm_2(tangent_vector_master);
                    tangent_vector_master /= tangent_magnitude_master;
                    Vector normal_vector_master = ZeroVector(3);
                    normal_vector_master[0] = tangent_vector_master[1];
                    normal_vector_master[1] = -tangent_vector_master[0];

                    new_master_skin_node->SetValue(NORMAL, normal_vector_master);
                    new_master_skin_node->SetValue(LOCAL_TANGENT, tangent_vector_master);


                    result_geometries_neumann(count_neumann_gp_master) = gp_in_brep;
                    count_neumann_gp_master++;
                    
                    continue;
                };
                IndexType id_new_node_master = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size()+1;
                auto new_master_skin_node = new Node(id_new_node_master, gp_in_skin_master);

                NodePointerVector empty_vector;
                empty_vector.push_back(new_master_skin_node); // Just it_node-plane neighbours
                gp_in_brep->SetValue(NEIGHBOUR_NODES, empty_vector);


                // compute normals and useful info
                std::vector<CoordinatesArrayType> global_space_derivatives_master;
                CoordinatesArrayType tangent_vector_master = curve_derivatives[1];
                double tangent_magnitude_master = norm_2(tangent_vector_master);
                tangent_vector_master /= tangent_magnitude_master;
                Vector normal_vector_master = ZeroVector(3);
                normal_vector_master[0] = tangent_vector_master[1];
                normal_vector_master[1] = -tangent_vector_master[0];

                new_master_skin_node->SetValue(NORMAL, normal_vector_master);
                new_master_skin_node->SetValue(LOCAL_TANGENT, tangent_vector_master);

                // compute the curvature
                CoordinatesArrayType curve_first_derivative_vector_master = curve_derivatives[1];
                CoordinatesArrayType curve_second_derivative_vector_master = curve_derivatives[2];

                double curvature_master = 
                                    norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector_master, curve_second_derivative_vector_master)) 
                                    / pow(norm_2(curve_first_derivative_vector_master), 3);
                new_master_skin_node->SetValue(CURVATURE, curvature_master);

                // project the skin_node to the skin boundary on the slave side
                double best_curve_distance = 1e16;
                CoordinatesArrayType best_curve_projection;
                CoordinatesArrayType best_curve_projection_local_coord;
                bool is_converged_at_least_once = false;
                IndexType best_slave_curve_id;

                CoordinatesArrayType skin_node_deformed_coordinates;
                GetDeformedPosition(gp_in_skin_master, *mrMasterModelPart, mSparseBrepMatrixMaster, skin_node_deformed_coordinates);

                for (auto &i_slave_curve : mrSlaveSkinModelPart->Geometries())
                { 
                    if (i_slave_curve.GetValue(IDENTIFIER) != slave_layer_name) continue;
                    CoordinatesArrayType local_coord = ZeroVector(3);
                    local_coord[0] = 0; //initial guess
                    CoordinatesArrayType projected_point;
                    double distance;

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

                    if ((is_converged || distance < 1e-2) && distance < 2e0)
                    {
                        is_converged_at_least_once = true;
                    
                        if (distance < best_curve_distance) {
                            best_curve_distance = distance;
                            best_curve_projection = projected_point;
                            best_curve_projection_local_coord = local_coord; 
                            best_slave_curve_id = slave_nurbs_curve_id;

                            std::vector<array_1d<double, 3>> ppp;
                            //FIXME:
                            slave_nurbs_curve_geometry->GlobalSpaceDerivatives(ppp, best_curve_projection_local_coord, 1);
                            best_curve_projection = ppp[0];
                        }
                    } 
                }
                if (!is_converged_at_least_once) 
                {
                    // save it as Neumann quadrature point and continue
                    result_geometries_neumann(count_neumann_gp_master) = gp_in_brep;
                    count_neumann_gp_master++;
                    continue;
                };
               
                // associate the best skin_node with its projection
                IndexType idNewNode = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size()+1;
                auto new_slave_skin_node = new Node(idNewNode, best_curve_projection);

                NodePointerVector neighbour_vector;
                neighbour_vector.push_back(new_slave_skin_node); 
                new_master_skin_node->SetValue(NEIGHBOUR_NODES, neighbour_vector);

                // compute normals and useful info
                std::vector<CoordinatesArrayType> global_space_derivatives;
                SizeType derivative_order = 2;
                mrSlaveSkinModelPart->pGetGeometry(best_slave_curve_id)->GlobalSpaceDerivatives(global_space_derivatives, best_curve_projection_local_coord, derivative_order);
                CoordinatesArrayType tangent_vector = global_space_derivatives[1];
                double tangent_magnitude = norm_2(tangent_vector);
                tangent_vector /= tangent_magnitude;
                Vector normal_vector = ZeroVector(3);
                normal_vector[0] = tangent_vector[1];
                normal_vector[1] = -tangent_vector[0];

                new_slave_skin_node->SetValue(NORMAL, normal_vector);
                new_slave_skin_node->SetValue(LOCAL_TANGENT, tangent_vector);

                // compute the curvature
                CoordinatesArrayType curve_first_derivative_vector = global_space_derivatives[1];
                CoordinatesArrayType curve_second_derivative_vector = global_space_derivatives[2];

                double curvature = norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector, curve_second_derivative_vector)) / pow(norm_2(curve_first_derivative_vector), 3);
                new_slave_skin_node->SetValue(CURVATURE, curvature);

                //-------------------------------------------------------------
                //  project skin slave to surrogate slave
                //-------------------------------------------------------------
                // find the closest slave gauss point to the projected skin node on the slave boundary
                // perform ray tracing to the slave geometry starting from the new_slave_skin_node with its normal

                // find the closest breps (want to try the ray tracing only on the potential breps)
                CoordinatesArrayType slave_skin_node_coords = new_slave_skin_node->Coordinates();
                
                IndexType best_brep_id = -1;
                double brep_intersection_local_variable = 0;
                ProjectBackToSurrogateBoundary(*mrSlaveModelPart, slave_skin_node_coords, mSparseBrepMatrixSlave, normal_vector,
                                            best_brep_id, brep_intersection_local_variable);

                // FIXME:
                if (best_brep_id == -1) continue;


                KRATOS_ERROR_IF(best_brep_id == -1) << "::[IgaContactProcessSbm]:: No brep found for the slave skin node " 
                                            << *new_slave_skin_node << std::endl;gp_in_brep->SetValue(IDENTIFIER, "active");

                auto p_brep = mrSlaveModelPart->pGetGeometry(best_brep_id);
                auto brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep);

                if (brep_geometry->GetValue(ACTIVATION_LEVEL) == 1)
                {
                    project_vertex_to_slave = true;
                    // continue;
                }
                if (brep_geometry->GetValue(ACTIVATION_LEVEL) == 0)
                {
                    continue;
                }

                std::ofstream outputFile("txt_files/Contact_Projection_Coordinates.txt", std::ios::app);
                        outputFile <<  new_master_skin_node->X() << " " << new_master_skin_node->Y() << " "  << new_slave_skin_node->X() << " " << new_slave_skin_node->Y() <<"\n";
                        outputFile.close();


                IntegrationPoint<1> integration_point_slave(brep_intersection_local_variable);
                IndexType number_of_shape_functions_derivatives = 5;
                IntegrationPointsArrayType surrogate_integration_points_list; // only contains the new gauss point

                surrogate_integration_points_list.push_back(integration_point_slave);
                IntegrationInfo integration_info_slave = brep_geometry->GetDefaultIntegrationInfo();

                GeometriesArrayType quadrature_point_list; // only contains the new gauss point

                brep_geometry->pGetCurveOnSurface()->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives, 
                                                            surrogate_integration_points_list, integration_info_slave);

                
                new_slave_skin_node->SetValue(NEIGHBOUR_GEOMETRY, quadrature_point_list(0));

                result_geometries_contact(count_contact_gp_master) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCouplingGeometry2D(
                                                    gp_in_brep, quadrature_point_list(0));
                
                integral += gp_in_brep->IntegrationPoints()[0].Weight();
                                
                result_geometries_contact(count_contact_gp_master)->SetValue(MARKER_MESHES, master_knot_step_uv);
                
                count_contact_gp_master++;
            }


            // FIXME: projection to slave of master vertices
            if (project_vertex_to_slave)
            {
                Vector brep_domain_interval;
                master_brep_geometry->DomainInterval(brep_domain_interval);
                std::vector<CoordinatesArrayType> brep_vertex(2); brep_vertex[0] = ZeroVector(3); brep_vertex[1] = ZeroVector(3);
                CoordinatesArrayType brep_vertex_1_local_coords = ZeroVector(3); CoordinatesArrayType brep_vertex_2_local_coords = ZeroVector(3);
                brep_vertex_1_local_coords[0] = brep_domain_interval[0];
                brep_vertex_2_local_coords[0] = brep_domain_interval[1];

                master_brep_geometry->GlobalCoordinates(brep_vertex[0], brep_vertex_1_local_coords);
                master_brep_geometry->GlobalCoordinates(brep_vertex[1], brep_vertex_2_local_coords);

                // KRATOS_WATCH(brep_vertex)

                ///////////////////////////////////////////////////////////////////////////////////////////////
                // need to know which one of the vertices is connvected to the fully contact element
                for (IndexType i_vertex = 0; i_vertex < 2; i_vertex++)
                {
                     // compute projection on skin master
                    std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));
                    CoordinatesArrayType gp_in_skin_master = ZeroVector(3);
                    
                    bool is_projected = ProjectToSkinBoundary(mrMasterSkinModelPart, master_layer_name, brep_vertex[i_vertex], gp_in_skin_master, curve_derivatives, 10);
                    
                    if (!is_projected) {
                        continue;
                    };

                    // project the skin_node to the skin boundary on the slave side
                    double best_curve_distance = 1e16;
                    CoordinatesArrayType best_curve_projection;
                    CoordinatesArrayType best_curve_projection_local_coord;
                    bool is_converged_at_least_once = false;
                    IndexType best_slave_curve_id;

                    CoordinatesArrayType skin_node_deformed_coordinates;
                    GetDeformedPosition(gp_in_skin_master, *mrMasterModelPart, mSparseBrepMatrixMaster, skin_node_deformed_coordinates);

                    for (auto &i_slave_curve : mrSlaveSkinModelPart->Geometries())
                    { 
                        if (i_slave_curve.GetValue(IDENTIFIER) != slave_layer_name) continue;
                        CoordinatesArrayType local_coord = ZeroVector(3);
                        local_coord[0] = 0; //initial guess
                        CoordinatesArrayType projected_point;
                        double distance;

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

                        if ((is_converged || distance < 1e-2))// && distance < 2e0)
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
                    // compute normals and useful info
                    std::vector<CoordinatesArrayType> global_space_derivatives;
                    SizeType derivative_order = 2;
                    mrSlaveSkinModelPart->pGetGeometry(best_slave_curve_id)->GlobalSpaceDerivatives(global_space_derivatives, best_curve_projection_local_coord, derivative_order);
                    CoordinatesArrayType tangent_vector = global_space_derivatives[1];
                    double tangent_magnitude = norm_2(tangent_vector);
                    tangent_vector /= tangent_magnitude;
                    Vector normal_vector = ZeroVector(3);
                    normal_vector[0] = tangent_vector[1];
                    normal_vector[1] = -tangent_vector[0];
                    //-------------------------------------------------------------
                    //  project skin slave to surrogate slave
                    //-------------------------------------------------------------
                    // std::ofstream outputFile2("txt_files/Cut_Elements_Coordinates.txt", std::ios::app);
                    //     outputFile2 <<  best_curve_projection[0] << " " << best_curve_projection[1] << " "  << 0 << " " << 0 <<"\n";
                    //     outputFile2.close();

                    IndexType best_brep_id = -1;
                    double brep_intersection_local_variable = 0;
                    ProjectBackToSurrogateBoundary(*mrSlaveModelPart, best_curve_projection, mSparseBrepMatrixSlave, normal_vector,
                                                best_brep_id, brep_intersection_local_variable);
        
                    KRATOS_ERROR_IF(best_brep_id == -1) << "::[IgaContactProcessSbm]:: No brep found for the slave skin node "
                                                << best_curve_projection << std::endl;
                                                
                    
                    brep_id_forward_projections.push_back(best_brep_id);
                    brep_local_coord_forward_projections.push_back(brep_intersection_local_variable);
                    //FIXME: actually doing it double the time for repeated vertices

                    // KRATOS_WATCH(best_brep_id)
                    // KRATOS_WATCH(best_curve_projection)
                }

            }
        }

        // create contact conditions
        result_geometries_contact.resize(count_contact_gp_master);
        this->CreateConditions(
                        result_geometries_contact.ptr_begin(), result_geometries_contact.ptr_end(),
                        *mrContactModelPart, name, id, mpPropMaster, mpPropSlave);
        
        // create neumann conditions
        std::string default_condition_name = "SBMLoadSolid2DCondition";
        result_geometries_neumann.resize(count_neumann_gp_master);
        Vector meshSizes_uv = mrSlaveModelPart->GetParentModelPart().GetValue(MARKER_MESHES);

        if (mrContactModelPart->GetRootModelPart().Conditions().size() > 0)
            id = mrContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;
        
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(default_condition_name);
        ModelPart::ConditionsContainerType new_condition_list;

        KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
            << "Creating conditions of type " << default_condition_name
            << " in " << mrContactModelPart->GetParentModelPart().Name() << "-SubModelPart." << std::endl;

        IndexType count_cond = 0;
        for (auto it = result_geometries_neumann.ptr_begin(); it != result_geometries_neumann.ptr_end(); ++it) {
            
            auto neigh_nodes = (*it)->GetValue(NEIGHBOUR_NODES);

            std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
            outputFile <<  neigh_nodes(0)->X() << " " << neigh_nodes(0)->Y() << " "  << (*it)->Center().X() << " " << (*it)->Center().Y() <<"\n";
            outputFile.close();

            new_condition_list.push_back(rReferenceCondition.Create(id, (*it), mpPropMaster));
            
            new_condition_list.GetContainer()[count_cond]->SetValue(MARKER_MESHES, meshSizes_uv);
            new_condition_list.GetContainer()[count_cond]->SetValue(IDENTIFIER, "outer");
                            
            id++;
            count_cond++;
        }
        mrContactModelPart->GetParentModelPart().AddConditions(new_condition_list.begin(), new_condition_list.end());

        // // Set non active slave breps to SBM load condition to zero
        // std::string default_condition_name = "SBMLoadSolid2DCondition";
        SizeType rIdCounter = 1;
        if (mrContactModelPart->GetRootModelPart().Conditions().size() > 0)
            rIdCounter = mrContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;

        for (auto& r_geometry_slave : mrSlaveModelPart->Geometries()) {
            
            if (r_geometry_slave.GetValue(ACTIVATION_LEVEL) == 2) continue;
            int slave_brep_id = r_geometry_slave.Id();
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);

            auto slave_brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);

            KRATOS_ERROR_IF(!slave_brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << slave_brep_id << "is not a Brep." << std::endl;

            if (r_geometry_slave.GetValue(ACTIVATION_LEVEL) == 1) 
            {
                std::vector<double> spans;
                slave_brep_geometry->SpansLocalSpace(spans);

                for (IndexType i = 0; i < brep_local_coord_forward_projections.size(); i++) 
                {
                    if (brep_id_forward_projections[i] == slave_brep_id) {
                        spans.push_back(brep_local_coord_forward_projections[i]);
                    }
                }

                std::sort(spans.begin(), spans.end());

                // remove duplicates with tolerance
                auto last = std::unique(spans.begin(), spans.end(), [](double a, double b) {
                    return std::abs(a - b) < 1e-4;
                });
                spans.erase(last, spans.end());

                // remove duplicates with tolerance
                std::vector<double> spans_final(2);

                if (slave_brep_geometry->GetValue(IDENTIFIER) == "first_vertex_activated")
                {
                    spans_final[0] = spans[spans.size()-1];
                    spans_final[1] = spans[spans.size()-2];
                }
                else if (slave_brep_geometry->GetValue(IDENTIFIER) == "last_vertex_activated")
                {
                    spans_final[0] = spans[0];
                    spans_final[1] = spans[1];
                }
                else
                {
                    KRATOS_ERROR << "ERROR IN THE CREATION OF THE MORTAR SPACE FOR THE SLAVE BOUNDARY " << std::endl;
                }
                
                ///////////////////////////////////////////
                GeometriesArrayType geometries;
                SizeType shape_function_derivatives_order = 3;
                if (mParameters.Has("shape_function_derivatives_order")) {
                    shape_function_derivatives_order = mParameters["shape_function_derivatives_order"].GetInt();
                }
                else {
                    KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
                        << "shape_function_derivatives_order is not provided and thus being considered as 3. " << std::endl;
                }

                std::string quadrature_method = mParameters.Has("quadrature_method")
                    ? mParameters["integration_rule"].GetString()
                    : "GAUSS";
                IntegrationInfo integration_info = slave_brep_geometry->GetDefaultIntegrationInfo();

                if (mParameters.Has("number_of_integration_points_per_span")) {
                    for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                        integration_info.SetNumberOfIntegrationPointsPerSpan(i, mParameters["number_of_integration_points_per_span"].GetInt());
                    }
                }
                
                for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                    if (quadrature_method == "GAUSS") {
                        integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
                    }
                    else if (quadrature_method == "EXTENDED_GAUSS") {
                        integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
                    }
                    else if (quadrature_method == "GRID") {
                        integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
                    }
                    else {
                        KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                            << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
                    }
                }


                IntegrationPointsArrayType integration_points;
                IntegrationPointUtilities::CreateIntegrationPoints1D(integration_points, spans_final, integration_info);

                slave_brep_geometry->CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_points, integration_info, false);
            
                for (auto it = geometries.ptr_begin(); it != geometries.ptr_end(); ++it) {
                    
                    std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));
                    CoordinatesArrayType best_projection = ZeroVector(3);

                    std::string projection_layer_name = "";

                    ProjectToSkinBoundary(mrSlaveSkinModelPart, projection_layer_name, (*it)->Center(), best_projection, curve_derivatives, 10);

                    IndexType id_new_node_slave = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size()+1;
                    auto new_slave_skin_node = new Node(id_new_node_slave, best_projection);
                    // mrSlaveSkinModelPart->GetSubModelPart(projection_layer_name).AddNode(new_slave_skin_node);

                    //FIXME: To do automatically in some way
                    Vector force_vector = ZeroVector(3); 
                    // force_vector[1] = -0.3;
                    new_slave_skin_node->SetValue(FORCE, force_vector);

                    NodePointerVector empty_vector;
                    empty_vector.push_back(new_slave_skin_node); // Just it_node-plane neighbours
                    (*it)->SetValue(NEIGHBOUR_NODES, empty_vector);

                    std::ofstream outputFile2("txt_files/Cut_Elements_Coordinates.txt", std::ios::app);
                        outputFile2 <<  new_slave_skin_node->X() << " " << new_slave_skin_node->Y() << " "  << (*it)->Center().X() << " " << (*it)->Center().Y() <<"\n";
                        outputFile2.close();

                    // compute normals and useful info
                    std::vector<CoordinatesArrayType> global_space_derivatives_slave;
                    CoordinatesArrayType tangent_vector_slave = curve_derivatives[1];
                    double tangent_magnitude_slave = norm_2(tangent_vector_slave);
                    tangent_vector_slave /= tangent_magnitude_slave;
                    Vector normal_vector_slave = ZeroVector(3);
                    normal_vector_slave[0] = tangent_vector_slave[1];
                    normal_vector_slave[1] = -tangent_vector_slave[0];

                    new_slave_skin_node->SetValue(NORMAL, normal_vector_slave);
                    new_slave_skin_node->SetValue(LOCAL_TANGENT, tangent_vector_slave);

                    std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
                    outputFile <<  new_slave_skin_node->X() << " " << new_slave_skin_node->Y() << " "  << (*it)->Center().X() << " " << (*it)->Center().Y() <<"\n";
                    outputFile.close();

                    new_condition_list.push_back(
                        rReferenceCondition.Create(rIdCounter, (*it), mpPropSlave));
                    
                    new_condition_list.GetContainer()[count_cond]->SetValue(MARKER_MESHES, meshSizes_uv);
                                    
                    // for (SizeType i = 0; i < (*it)->size(); ++i) {
                    //     // These are the control points associated with the basis functions involved in the condition we are creating
                    //     // rModelPart.AddNode((*it)->pGetPoint(i));
                    //     mrContactModelPart->GetParentModelPart().Nodes().push_back((*it)->pGetPoint(i));
                    // }
                    rIdCounter++;
                    count_cond++;
                }

                mrContactModelPart->GetParentModelPart().AddConditions(new_condition_list.begin(), new_condition_list.end());
            
            }
            else
            {
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
                    
                    std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));
                    CoordinatesArrayType best_projection = ZeroVector(3);

                    std::string projection_layer_name = "";
                    ProjectToSkinBoundary(mrSlaveSkinModelPart, projection_layer_name, (*it)->Center(), best_projection, curve_derivatives, 10);
                        

                    IndexType id_new_node_slave = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size()+1;
                    auto new_slave_skin_node = new Node(id_new_node_slave, best_projection);
                    // mrSlaveSkinModelPart->GetSubModelPart(projection_layer_name).AddNode(new_slave_skin_node);

                    //FIXME: To do automatically in some way
                    Vector force_vector = ZeroVector(3); 
                    // force_vector[1] = -0.3;
                    new_slave_skin_node->SetValue(FORCE, force_vector);

                    NodePointerVector empty_vector;
                    empty_vector.push_back(new_slave_skin_node); // Just it_node-plane neighbours
                    (*it)->SetValue(NEIGHBOUR_NODES, empty_vector);

                    // compute normals and useful info
                    std::vector<CoordinatesArrayType> global_space_derivatives_slave;
                    CoordinatesArrayType tangent_vector_slave = curve_derivatives[1];
                    double tangent_magnitude_slave = norm_2(tangent_vector_slave);
                    tangent_vector_slave /= tangent_magnitude_slave;
                    Vector normal_vector_slave = ZeroVector(3);
                    normal_vector_slave[0] = tangent_vector_slave[1];
                    normal_vector_slave[1] = -tangent_vector_slave[0];

                    new_slave_skin_node->SetValue(NORMAL, normal_vector_slave);
                    new_slave_skin_node->SetValue(LOCAL_TANGENT, tangent_vector_slave);

                    std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
                    outputFile <<  new_slave_skin_node->X() << " " << new_slave_skin_node->Y() << " "  << (*it)->Center().X() << " " << (*it)->Center().Y() <<"\n";
                    outputFile.close();

                    new_condition_list.push_back(
                        rReferenceCondition.Create(rIdCounter, (*it), mpPropSlave));
                    
                    new_condition_list.GetContainer()[count_cond]->SetValue(MARKER_MESHES, meshSizes_uv);
                                    
                    // for (SizeType i = 0; i < (*it)->size(); ++i) {
                    //     // These are the control points associated with the basis functions involved in the condition we are creating
                    //     // rModelPart.AddNode((*it)->pGetPoint(i));
                    //     mrContactModelPart->GetParentModelPart().Nodes().push_back((*it)->pGetPoint(i));
                    // }
                    rIdCounter++;
                    count_cond++;
                }

                mrContactModelPart->GetParentModelPart().AddConditions(new_condition_list.begin(), new_condition_list.end());
            }
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
            
            // for (SizeType i = 0; i < (*it)->size(); ++i) {
            //     // rModelPart.AddNode((*it)->pGetPoint(i));
            //     rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            // }
            rIdCounter++;
        }
        
        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());

        // EntitiesUtilities::InitializeEntities<Condition>(rModelPart);

    }

    void IgaContactProcessSbm::GetDeformedPosition(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& SparseBrepMatrix,
        CoordinatesArrayType& rPointDeformedCoordinates)
    {
        const double epsilon = 1e-12;
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(MARKER_MESHES)/2;
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
        int max_i = 2; int max_j = 2;
        while (!index_found)
        {
            first_i--; first_j--;
            max_i++; max_j++;
            for (int i = first_i; i < max_i; i++)
            {
                if (knot_span_surrounded_u_id+i < 0 || knot_span_surrounded_u_id+i > SparseBrepMatrix.size1()-1) continue;
                for (int j = first_j; j < max_j; j++)
                {
                    if (knot_span_surrounded_v_id+j < 0 || knot_span_surrounded_v_id+j > SparseBrepMatrix.size2()-1)  continue;
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
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(MARKER_MESHES)/2;
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
        int max_i = 2; int max_j = 2;
        while (!index_found)
        {
            first_i--; first_j--;
            max_i++; max_j++;
            for (int i = first_i; i < max_i; i++)
            {
                if (knot_span_surrounded_u_id+i < 0 || knot_span_surrounded_u_id+i > SparseBrepMatrix.size1()-1) continue;
                for (int j = first_j; j < max_j; j++)
                {
                    if (knot_span_surrounded_v_id+j < 0 || knot_span_surrounded_v_id+j > SparseBrepMatrix.size2()-1)  continue;
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

        double curve_interval_min; double curve_interval_max; 
        if (curve_interval[0] < curve_interval[1]) {
            curve_interval_min = curve_interval[0];
            curve_interval_max = curve_interval[1];
        }
        else {
            curve_interval_min = curve_interval[1];
            curve_interval_max = curve_interval[0];
        }

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

                 // Check if the parameter gets out of its interval of definition and if so clamp it
                // back to the boundaries
                // CoordinatesArrayType closest_t = ZeroVector(3);
                // int check = rPairedGeometry.ClosestPointLocalToLocalSpace(
                //     t, t,2);

                int check = 0;
                
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
                    break;
                }

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



    std::vector<IndexType> IgaContactProcessSbm::FindClosestBrepId(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& SparseBrepMatrix)
    {
        const double epsilon = 1e-12;
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(MARKER_MESHES)/2;
        Vector parameter_space_extremes = rModelPart.GetParentModelPart().GetValue(LOAD_MESHES);

        // compute the knot span where the point lies
        int knot_span_surrounded_u_id = floor((rPointGlobalCoordinates[0] - parameter_space_extremes[0]+epsilon)/knot_step_uv[0]);
        int knot_span_surrounded_v_id = floor((rPointGlobalCoordinates[1] - parameter_space_extremes[1]+epsilon)/knot_step_uv[1]);

        // check where the surrogate breps are in the boundary knot span
        std::vector<IndexType> list_closest_brep_id;

        bool index_found = false;
        int number_of_close_span_to_check = 6;
        for (int i = -number_of_close_span_to_check; i < number_of_close_span_to_check; i++)
        {
            if (knot_span_surrounded_u_id+i < 0) continue;
            for (int j = -number_of_close_span_to_check; j < number_of_close_span_to_check; j++)
            {
                if (knot_span_surrounded_v_id+j < 0) continue;
                if (SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j) != 0)
                {
                    list_closest_brep_id.push_back(SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j));
                    index_found = true;
                }
            }
        }
        int first_i = -number_of_close_span_to_check; int first_j = -number_of_close_span_to_check;
        int max_i = number_of_close_span_to_check; int max_j = number_of_close_span_to_check;
        while (!index_found)
        {
            first_i--; first_j--;
            max_i++; max_j++;
            for (int i = first_i; i < max_i; i++)
            {
                if (knot_span_surrounded_u_id+i < 0 || knot_span_surrounded_u_id+i > SparseBrepMatrix.size1()-1) continue;
                for (int j = first_j; j < max_j; j++)
                {
                    if (knot_span_surrounded_v_id+j < 0 || knot_span_surrounded_v_id+j > SparseBrepMatrix.size2()-1) continue;
                    if (SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j) != 0)
                    {
                        list_closest_brep_id.push_back(SparseBrepMatrix(knot_span_surrounded_u_id+i, knot_span_surrounded_v_id+j));
                        index_found = true;
                        break;
                    }
                }
                if (index_found) break;
            }
        }
        KRATOS_ERROR_IF(list_closest_brep_id.size() == 0) << "::[IgaContactProcessSbm]:: no boundary surrogate brep around true point: " <<
                                         rPointGlobalCoordinates << " found in knot span: (" << 
                                         knot_span_surrounded_u_id << ", " << knot_span_surrounded_v_id << ")" << std::endl;

        return list_closest_brep_id;
    }


    bool IgaContactProcessSbm::ProjectPointViaRayTracingToMasterCurveDeformed(
        const CoordinatesArrayType& point_slave_deformed,
        const Vector& normal_slave_deformed,
        const NurbsCurveGeometryType& master_curve,
        const ModelPart& master_model_part,
        const SparseMatrixType& sparse_matrix_master,
        CoordinatesArrayType& rBestProjectedPoint,
        CoordinatesArrayType& rBestProjectedPointLocalCoords,
        double& rBestDistance,
        double tolerance,
        int max_iter,
        int n_initial_guesses)
    {
        CoordinatesArrayType local_coords = ZeroVector(3);
        CoordinatesArrayType global_coords = ZeroVector(3);
        CoordinatesArrayType global_coords_updated = ZeroVector(3);
        Vector tangent_deformed = ZeroVector(3);

        auto interval = master_curve.DomainInterval();
        double t_min = interval.GetT0();
        double t_max = interval.GetT1();

        int larger = std::abs(normal_slave_deformed[0]) > std::abs(normal_slave_deformed[1]) ? 0 : 1;
        int smaller = 1 - larger;

        bool is_converged = false;
        rBestDistance = 1e12;

        for (int i_guess = 0; i_guess < n_initial_guesses; ++i_guess) {
            double res = tolerance + 1.0;
            int iter = 0;

            double t = t_min + (t_max - t_min) * double(i_guess) / double(n_initial_guesses - 1);
            local_coords[0] = t;

            master_curve.GlobalCoordinates(global_coords, local_coords);

            // Deformazione punto master
            CoordinatesArrayType global_coords_deformed;
            GetDeformedPosition(global_coords, master_model_part, sparse_matrix_master, global_coords_deformed);

            // Derivate deformate (prima derivata curva + gradiente spostamento)
            Matrix gradient_deformation = ZeroMatrix(2,2);
            Matrix hessian_deformation = ZeroMatrix(2,3);
            GetDeformedGradient(global_coords, master_model_part, sparse_matrix_master, gradient_deformation, hessian_deformation);

            std::vector<CoordinatesArrayType> global_derivatives;
            master_curve.GlobalSpaceDerivatives(global_derivatives, local_coords, 1);

            for (int i_dim = 0; i_dim < 2; ++i_dim) {
                tangent_deformed[i_dim] = global_derivatives[1][i_dim]
                                        + gradient_deformation(i_dim, 0)*global_derivatives[1][0]
                                        + gradient_deformation(i_dim, 1)*global_derivatives[1][1];
            }

            double s_line = (global_coords_deformed[larger] - point_slave_deformed[larger]) / normal_slave_deformed[larger];
            double smaller_line = s_line * normal_slave_deformed[smaller] + point_slave_deformed[smaller];
            res = smaller_line - global_coords_deformed[smaller];

            while (std::abs(res) > tolerance && iter < max_iter) {
                double df = normal_slave_deformed[smaller]/normal_slave_deformed[larger] * tangent_deformed[larger]
                        - tangent_deformed[smaller];
                if (std::abs(df) < 1e-10) break;

                double delta_t = res / df;
                t -= delta_t;

                if (t < t_min || t > t_max) break;

                local_coords[0] = t;
                master_curve.GlobalCoordinates(global_coords, local_coords);
                GetDeformedPosition(global_coords, master_model_part, sparse_matrix_master, global_coords_deformed);

                GetDeformedGradient(global_coords, master_model_part, sparse_matrix_master, gradient_deformation, hessian_deformation);
                master_curve.GlobalSpaceDerivatives(global_derivatives, local_coords, 1);

                for (int i_dim = 0; i_dim < 2; ++i_dim) {
                    tangent_deformed[i_dim] = global_derivatives[1][i_dim]
                                            + gradient_deformation(i_dim, 0)*global_derivatives[1][0]
                                            + gradient_deformation(i_dim, 1)*global_derivatives[1][1];
                }

                s_line = (global_coords_deformed[larger] - point_slave_deformed[larger]) / normal_slave_deformed[larger];
                smaller_line = s_line * normal_slave_deformed[smaller] + point_slave_deformed[smaller];
                res = smaller_line - global_coords_deformed[smaller];
                iter++;
            }

            if (std::abs(res) < tolerance) {
                double dist = norm_2(global_coords_deformed - point_slave_deformed);
                if (dist < rBestDistance) {
                    rBestDistance = dist;
                    rBestProjectedPoint = global_coords_deformed;
                    rBestProjectedPointLocalCoords = local_coords;
                    is_converged = true;
                }
            }
        }

        return is_converged;
    }





    bool IgaContactProcessSbm::ProjectToSkinBoundary(
        const ModelPart* pSkinModelPart,
        std::string &rLayerName,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rProjectedPoint,
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
            if  (rLayerName != "" && i_curve.GetValue(IDENTIFIER) != rLayerName) continue;

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

                bool is_projected = nurbs_curve_geometry->ProjectionPointGlobalToLocalSpace(rPoint, projected_point_local);

                if (!is_projected) continue;

                nurbs_curve_geometry->GlobalCoordinates(projected_point, projected_point_local);

                double curr_distance = norm_2(rPoint - projected_point);

                if (curr_distance < best_distance) {
                    best_distance = curr_distance;
                    rProjectedPoint = projected_point;
                    nurbs_curve_geometry->GlobalSpaceDerivatives(best_curve_derivatives, projected_point_local, 2);
                    is_projected_at_least_once = true;
                    best_layer_name = i_curve.GetValue(IDENTIFIER);
                }
            }
        }

        rCurveDerivatives = best_curve_derivatives;

        if (rLayerName == "") rLayerName = best_layer_name;

        if (!is_projected_at_least_once)
        {
            KRATOS_WARNING("::[IgaContactProcessSbm]:: no projection found on the skin boundary with layer ") << rLayerName 
                        << " for the point: " << rPoint << std::endl;
        }

        return is_projected_at_least_once;
    }


    // bool IgaContactProcessSbm::SkinRayTracingProjection( 
    //     const ModelPart &rOriginSurrogateModelPart,
    //     const ModelPart &rDestinationSkinModelPart,
    //     const CoordinatesArrayType& rOriginPoint,
    //     const CoordinatesArrayType& rDestinationPoint,
    //     const SparseMatrixType& rOriginSparseBrepMatrix,
    //     const SparseMatrixType& rDestinationSparseBrepMatrix)
    // {
    //     // get the deformed position of the slave vertex and the normal
    //     CoordinatesArrayType skin_vertex_slave_deformed;
    //     Matrix skin_vertex_gradient_deformation = ZeroMatrix(2, 2);
    //     Matrix skin_vertex_hessian_deformation = ZeroMatrix(2, 3);

    //     Vector slave_tangent_on_vertex_deformed = ZeroVector(3);

    //     GetDeformedPosition(rOriginPoint, rOriginSurrogateModelPart, rOriginSparseBrepMatrix, skin_vertex_slave_deformed);

    //     GetDeformedGradient(rOriginPoint, rOriginSurrogateModelPart, rOriginSparseBrepMatrix, skin_vertex_gradient_deformation, skin_vertex_hessian_deformation);

    //     for (int i_dim = 0; i_dim < 2; i_dim++) {
    //         slave_tangent_on_vertex_deformed[i_dim] = first_derivatives_projected_slave_vertices_on_true[count_vertex][i_dim] +
    //                                                 skin_vertex_gradient_deformation(i_dim,0)*first_derivatives_projected_slave_vertices_on_true[count_vertex][0] +
    //                                                 skin_vertex_gradient_deformation(i_dim,1)*first_derivatives_projected_slave_vertices_on_true[count_vertex][1];
    //     }

    //     slave_tangent_on_vertex_deformed /= norm_2(slave_tangent_on_vertex_deformed);
    //     CoordinatesArrayType slave_normal_on_vertex_deformed = ZeroVector(3);
    //     slave_normal_on_vertex_deformed[0] = -slave_tangent_on_vertex_deformed[1];
    //     slave_normal_on_vertex_deformed[1] = slave_tangent_on_vertex_deformed[0];

    //     // CHATGPT CODE
    //     CoordinatesArrayType best_projected_point_master = ZeroVector(3);
    //     Vector normal_vector = ZeroVector(3);
    //     double best_distance_to_master = 1e12;

    //     bool has_converged_at_least_once = false;

    //     for (auto& i_master_curve : mrMasterSkinModelPart->Geometries()) {
    //         if (i_master_curve.GetValue(IDENTIFIER) != master_layer_name) continue;

    //         int master_id = i_master_curve.Id();
    //         auto p_master_geom = mrMasterSkinModelPart->pGetGeometry(master_id);
    //         auto master_curve = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_master_geom);
    //         KRATOS_ERROR_IF(!master_curve) << ":::[IgaContactProcessSbm]::: Curve id " << master_id << " is not a NurbsCurveGeometryType." << std::endl;

    //         CoordinatesArrayType projected_point;
    //         CoordinatesArrayType projected_point_local_coords = ZeroVector(3);
    //         double projection_distance;
    //         bool converged = ProjectPointViaRayTracingToMasterCurveDeformed(
    //             skin_vertex_slave_deformed,
    //             slave_normal_on_vertex_deformed,
    //             *master_curve,
    //             *mrMasterModelPart,
    //             mSparseBrepMatrixMaster,
    //             projected_point,
    //             projected_point_local_coords,
    //             projection_distance);

    //         if (converged && projection_distance < best_distance_to_master && projection_distance < 2e0) {
    //             has_converged_at_least_once = true;
    //             best_distance_to_master = projection_distance;
    //             best_projected_point_master = projected_point;


    //             // compute normals 
    //             std::vector<CoordinatesArrayType> global_space_derivatives;
    //             SizeType derivative_order = 2;
    //             master_curve->GlobalSpaceDerivatives(global_space_derivatives, projected_point_local_coords, derivative_order);
    //             CoordinatesArrayType tangent_vector = global_space_derivatives[1];
    //             double tangent_magnitude = norm_2(tangent_vector);
    //             tangent_vector /= tangent_magnitude;
                
    //             normal_vector[0] = tangent_vector[1];
    //             normal_vector[1] = -tangent_vector[0];
    //         }
    //     }

    //     return has_converged_at_least_once;
    // }


    bool IgaContactProcessSbm::ProjectBackToSurrogateBoundary(
        const ModelPart &rSurrogateModelPart,
        const CoordinatesArrayType& rSkinPoint,
        const SparseMatrixType& rSparseBrepMatrix,
        const Vector & rSkinNormal,
        IndexType& rBrepId,
        double& rProjectionBrepLocalVariable
    )
    {
        std::vector closest_brep_id = FindClosestBrepId(rSkinPoint, rSurrogateModelPart, rSparseBrepMatrix);
        rBrepId = -1;
        CoordinatesArrayType best_surrogate_projection = ZeroVector(3);
        rProjectionBrepLocalVariable = 0;
        for (IndexType brep_id: closest_brep_id)
        {
            // KRATOS_WATCH("*******************************************************")
            auto p_brep = rSurrogateModelPart.pGetGeometry(brep_id);
            auto brep_geometry = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep);

            KRATOS_ERROR_IF(!brep_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << p_brep 
                                    << " is not a BrepCurveOnSurfaceGeometry." << std::endl;

            CoordinatesArrayType first_vertex = ZeroVector(3); CoordinatesArrayType first_vertex_local_coords = ZeroVector(3);
            CoordinatesArrayType second_vertex = ZeroVector(3); CoordinatesArrayType second_vertex_local_coords = ZeroVector(3);
            first_vertex_local_coords[0] = brep_geometry->DomainInterval().GetT0(); 
            second_vertex_local_coords[0] = brep_geometry->DomainInterval().GetT1();

            brep_geometry->GlobalCoordinates(first_vertex, first_vertex_local_coords);
            brep_geometry->GlobalCoordinates(second_vertex, second_vertex_local_coords);

            rProjectionBrepLocalVariable = 0;

            const double eps = 1e-10;
            // find brep direction (x parallel or y parallel)
            if (std::abs(first_vertex[0] - second_vertex[0]) < eps)
            {

                if (std::abs(rSkinNormal[0]) < eps) continue;
                // y aligned (x constant)
                double t_ray_tracing_intersection = (first_vertex[0] - rSkinPoint[0])/ rSkinNormal[0];

                //compute y of the intersection
                double y_intersection = rSkinPoint[1] + t_ray_tracing_intersection * rSkinNormal[1];

                // check if the intersection is inside the brep
                if (y_intersection > std::max(first_vertex[1], second_vertex[1]) + eps || y_intersection < std::min(first_vertex[1], second_vertex[1]) - eps)
                {
                    // the intersection is outside the brep
                    continue;
                }

                if (first_vertex[1] > second_vertex[1])
                    rProjectionBrepLocalVariable = first_vertex[1] - (y_intersection - second_vertex[1]);
                else
                    rProjectionBrepLocalVariable = y_intersection;

            }
            else if (std::abs(first_vertex[1] - second_vertex[1]) < eps)
            {
                if (std::abs(rSkinNormal[1]) < eps) continue;
                // x aligned (y constant)
                double t_ray_tracing_intersection = (first_vertex[1] - rSkinPoint[1])/ rSkinNormal[1];

                //compute x of the intersection
                double x_intersection = rSkinPoint[0] + t_ray_tracing_intersection * rSkinNormal[0];
                // check if the intersection is inside the brep
                
                if (x_intersection > std::max(first_vertex[0], second_vertex[0]) + eps || x_intersection < std::min(first_vertex[0], second_vertex[0]) - eps)
                {
                    // the intersection is outside the brep
                    continue;
                }

                if (first_vertex[0] > second_vertex[0])
                    rProjectionBrepLocalVariable = first_vertex[0] - (x_intersection - second_vertex[0]);
                else
                    rProjectionBrepLocalVariable = x_intersection;
            }
            else
            {
                KRATOS_ERROR << ":::[IgaContactProcessSbm]::: the brep with id " << brep_id 
                                    << " is not a aligned neither with x or y." << std::endl;
            }
            rBrepId = brep_id;
            break;
        }

        if (rBrepId == -1) return false;
        else return true;

        
    }


} // End namespace Kratos
