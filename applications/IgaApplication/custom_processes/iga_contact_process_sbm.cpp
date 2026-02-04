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
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "utilities/variable_utils.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"
#include "utilities/nurbs_utilities/projection_nurbs_contact_utilities.h"
#include "integration/integration_point_utilities.h"
#include "utilities/math_utils.h"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Kratos
{

    IgaContactProcessSbm::IgaContactProcessSbm(
        Model& rModel, Parameters ThisParameters) : 
        Process(), 
        mpModel(&rModel), 
        mParameters(ThisParameters)
    {
        // TODO: SBM->Body Fitted case or viceversa
        // the class only does the case SBM->SBM at the moment
        
        mEchoLevel = mParameters["echo_level"].GetInt();

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("analysis_model_part_name")) << "::[IgaContactProcessSbm]::" 
                            << " Missing \"analysis_model_part_name\" parameter. "<< std::endl; 
        
        KRATOS_ERROR_IF_NOT(ThisParameters.Has("contact_sub_model_part_name")) << "::[IgaContactProcessSbm]::" 
                            << " Missing \"contact_sub_model_part_name\" parameter. "<< std::endl; 

        if (mParameters.Has("integrate_on_true_boundary"))
            mIntegrateOnTrueBoundary = mParameters["integrate_on_true_boundary"].GetBool();

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
        Vector master_knot_step_uv = mrMasterModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES)/2; // need half the size of the edge
        auto master_parameter_space_extremes = mrMasterModelPart->GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);
        double master_domain_length_u = master_parameter_space_extremes[0][1] - master_parameter_space_extremes[0][0];
        double master_domain_length_v = master_parameter_space_extremes[1][1] - master_parameter_space_extremes[1][0];

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

            IndexType knot_span_surrounded_u_id = floor((brep_center[0]- master_parameter_space_extremes[0][0]+epsilon)/master_knot_step_uv[0]);
            IndexType knot_span_surrounded_v_id = floor((brep_center[1]- master_parameter_space_extremes[1][0]+epsilon)/master_knot_step_uv[1]);
            mSparseBrepMatrixMaster(knot_span_surrounded_u_id, knot_span_surrounded_v_id) = master_brep_id;
        }

        //---------------------
        // slave
        Vector slave_knot_step_uv = mrSlaveModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES)/2; // need half the size of the edge
        auto slave_parameter_space_extremes = mrSlaveModelPart->GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);
        double slave_domain_length_u = slave_parameter_space_extremes[0][1] - slave_parameter_space_extremes[0][0];
        double slave_domain_length_v = slave_parameter_space_extremes[1][1] - slave_parameter_space_extremes[1][0];

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
            IndexType knot_span_surrounded_u_id = floor((brep_center[0]- slave_parameter_space_extremes[0][0]+epsilon)/slave_knot_step_uv[0]);
            IndexType knot_span_surrounded_v_id = floor((brep_center[1]- slave_parameter_space_extremes[1][0]+epsilon)/slave_knot_step_uv[1]);

            mSparseBrepMatrixSlave(knot_span_surrounded_u_id, knot_span_surrounded_v_id) = slave_brep_id;
        }



        // create submodelparts of contact to store the neumann conditions on the master and slave side
        //TODO: create them as submoelparts of the contactModelPart, as well as for the InterfaceModelPart
        ModelPart& r_master_neumann_sub_model_part = mrContactModelPart->CreateSubModelPart("master");

        ModelPart& r_slave_neumann_sub_model_part = mrContactModelPart->CreateSubModelPart("slave");

        ModelPart& r_contact_sub_model_part = mrContactModelPart->CreateSubModelPart("contact");
   }

    
    void IgaContactProcessSbm::Execute(){
        using Clock = std::chrono::steady_clock;
        double time_project_to_skin = 0.0;
        double time_newton = 0.0;
        double time_back_projection = 0.0;
        double time_create_quadrature = 0.0;
        double time_create_quadrature_sbm = 0.0;
        double time_create_coupling = 0.0;
        
        std::string master_layer_name = mParameters["contact_parameters"]["master_model_part"]["layer_name"].GetString();
        const std::string slave_layer_name = mParameters["contact_parameters"]["slave_model_part"]["layer_name"].GetString();

        ConditionsArrayType& r_conditions_array = mrContactModelPart->GetParentModelPart().Conditions();
        KRATOS_TRACE_IF("Empty model part", r_conditions_array.size() == 0) << "YOUR CONTACT MODEL PART IS EMPTY" << std::endl;

        // if (r_conditions_array.size() > 0) return;
        block_for_each(r_conditions_array, [&](Condition& rCond) {
                rCond.Set(TO_ERASE, true);
        });

        NodesArrayType& r_nodes_array = mrContactModelPart->GetParentModelPart().Nodes();
        block_for_each(r_nodes_array, [&](Node& rNode) {
                rNode.Set(TO_ERASE, true);
        });

        mrContactModelPart->GetParentModelPart().RemoveConditionsFromAllLevels(TO_ERASE);
        mrContactModelPart->GetParentModelPart().RemoveNodesFromAllLevels(TO_ERASE);

        SizeType next_condition_id = 1;
        if (mrContactModelPart->GetRootModelPart().Conditions().size() > 0)
            next_condition_id = mrContactModelPart->GetRootModelPart().Conditions().back().Id() + 1;
        
        KRATOS_ERROR_IF_NOT(mParameters.Has("name"))
            << "\"name\" need to be specified." << std::endl;
        std::string name = mParameters["name"].GetString();

        const Vector master_knot_step_uv = mrMasterModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

        ModelPart& r_contact_root = mrContactModelPart->GetRootModelPart();

        IndexType next_node_id = 1;
        next_node_id = r_contact_root.Nodes().back().Id() + 1;

        for (auto& r_slave_geometry : mrSlaveModelPart->Geometries()) {
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(r_slave_geometry.Id());
            auto p_slave_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);
            p_slave_brep->SetValue(ACTIVATION_LEVEL, 0);
        }

    std::vector<MasterSegmentData> master_segments;
    master_segments.reserve(mrMasterModelPart->NumberOfGeometries());

    std::vector<NodeType::Pointer> master_nodes;
    master_nodes.reserve(2 * mrMasterModelPart->NumberOfGeometries());

    const double coordinate_tolerance = master_knot_step_uv[0]/1E9;

    for (auto& r_geometry_master : mrMasterModelPart->Geometries()) {
        auto p_master_geometry = mrMasterModelPart->pGetGeometry(r_geometry_master.Id());
        auto p_brep_curve = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_master_geometry);
        KRATOS_ERROR_IF(!p_brep_curve) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << r_geometry_master.Id() << " is not a Brep." << std::endl;

        Vector domain_interval;
        p_brep_curve->DomainInterval(domain_interval);

        CoordinatesArrayType local_begin = ZeroVector(3);
        CoordinatesArrayType local_end = ZeroVector(3);
        local_begin[0] = domain_interval[0];
        local_end[0] = domain_interval[1];

        CoordinatesArrayType global_begin = ZeroVector(3);
        CoordinatesArrayType global_end = ZeroVector(3);
        p_brep_curve->GlobalCoordinates(global_begin, local_begin);
        p_brep_curve->GlobalCoordinates(global_end, local_end);

        NodeType::Pointer p_node_begin = CreateOrRetrieveContactNode(global_begin, master_nodes, next_node_id, coordinate_tolerance);
        NodeType::Pointer p_node_end = CreateOrRetrieveContactNode(global_end, master_nodes, next_node_id, coordinate_tolerance);

        MasterSegmentData segment_data;
        segment_data.p_brep = p_brep_curve;
        segment_data.local_begin = domain_interval[0];
        segment_data.local_end = domain_interval[1];
        segment_data.p_node_begin = p_node_begin;
        segment_data.p_node_end = p_node_end;

        master_segments.push_back(segment_data);
    }

    const double projection_distance_limit = master_knot_step_uv[0]*1;
    const double projection_distance_fallback = master_knot_step_uv[0]/2;

    for (NodeType::Pointer p_node : master_nodes) {
        Vector projection_coordinates = ZeroVector(3);

        std::vector<array_1d<double, 3>> master_curve_derivatives(3, ZeroVector(3));
        CoordinatesArrayType master_skin_projection = ZeroVector(3);
        std::string projection_layer_name = master_layer_name;
        const auto t_proj_begin = Clock::now();
        bool is_projected_master = ProjectToSkinBoundary(
            mrMasterSkinModelPart,
            projection_layer_name,
            p_node->Coordinates(),
            master_skin_projection,
            master_curve_derivatives,
            10);
        time_project_to_skin += std::chrono::duration<double>(Clock::now() - t_proj_begin).count();

        if (!is_projected_master) {
            projection_layer_name.clear();
            const auto t_proj_begin_retry = Clock::now();
            is_projected_master = ProjectToSkinBoundary(
                mrMasterSkinModelPart,
                projection_layer_name,
                p_node->Coordinates(),
                master_skin_projection,
                master_curve_derivatives,
                10);
            time_project_to_skin += std::chrono::duration<double>(Clock::now() - t_proj_begin_retry).count();
        }

        if (!is_projected_master) {
            KRATOS_ERROR << "::[IgaModelerSbm]:: master surrogate point " << p_node->Coordinates() << "could not be projected to the skin master"
                         << "the best projection was: " << master_skin_projection << std::endl;
        }

        CoordinatesArrayType master_skin_deformed = ZeroVector(3);
        GetDeformedPosition(master_skin_projection, *mrMasterModelPart, mSparseBrepMatrixMaster, master_skin_deformed);

        bool best_projection_found = false;
        double best_projection_distance = std::numeric_limits<double>::max();
        CoordinatesArrayType best_slave_projection = ZeroVector(3);
        CoordinatesArrayType best_slave_local_coords = ZeroVector(3);
        IndexType best_slave_curve_id = 0;
        bool best_is_converged = false;
        const bool is_target_node =
            (std::abs(p_node->X()) < 1.0e-6) &&
            (std::abs(p_node->Y() - 20.96) < 1.0e-6) &&
            (std::abs(p_node->Z()) < 1.0e-6);

        for (auto& r_slave_curve_geometry : mrSlaveSkinModelPart->Geometries()) {
            if (r_slave_curve_geometry.GetValue(IDENTIFIER) != slave_layer_name) {
                continue;
            }

            auto p_slave_curve_geometry = mrSlaveSkinModelPart->pGetGeometry(r_slave_curve_geometry.Id());
            auto p_slave_nurbs_curve = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_slave_curve_geometry);
            KRATOS_ERROR_IF(!p_slave_nurbs_curve) <<  ":::[IgaContactProcessSbm]::: the geometry with id "
                << r_slave_curve_geometry.Id() << " is not a NurbsCurveGeometryType." << std::endl;

            CoordinatesArrayType local_coords = ZeroVector(3);
            CoordinatesArrayType projected_point = ZeroVector(3);
            double distance = std::numeric_limits<double>::max();

            const auto t_newton_begin = Clock::now();
            bool is_converged = NewtonRaphsonCurveOnDeformed(
                local_coords,
                master_skin_deformed,
                projected_point,
                *p_slave_nurbs_curve,
                distance,
                20,
                10,
                1e-9);
            time_newton += std::chrono::duration<double>(Clock::now() - t_newton_begin).count();

            if ((is_converged || distance < projection_distance_fallback) && distance < best_projection_distance) {
                best_projection_found = true;
                best_projection_distance = distance;
                best_slave_local_coords = local_coords;
                best_slave_curve_id = r_slave_curve_geometry.Id();
                best_is_converged = is_converged;

                std::vector<array_1d<double, 3>> derivatives(2, ZeroVector(3));
                p_slave_nurbs_curve->GlobalSpaceDerivatives(derivatives, local_coords, 1);
                best_slave_projection = derivatives[0];
            }
        }

        if (!best_is_converged) {
            if (p_node->X() < 0.05) {
                projection_coordinates[0] = p_node->X();
                projection_coordinates[1] = p_node->Y();
                projection_coordinates[2] = p_node->Z();
                p_node->SetValue(ACTIVATION_LEVEL, 1.0);
                p_node->SetValue(PROJECTION_NODE_COORDINATES, projection_coordinates);
            } else {
                p_node->SetValue(ACTIVATION_LEVEL, 0.0);
            }
            continue;
        }

        const bool is_active = best_projection_found &&
            ((best_projection_distance < projection_distance_limit) || (best_projection_distance < projection_distance_fallback));

        // if (is_target_node) {
        //     KRATOS_INFO("IgaContactProcessSbm")
        //         << "Target master node coords " << p_node->Coordinates()
        //         << " best_projection_found " << best_projection_found
        //         << " best_is_converged " << best_is_converged
        //         << " best_projection_distance " << best_projection_distance
        //         << " limit " << projection_distance_limit
        //         << " fallback " << projection_distance_fallback
        //         << " best_slave_curve_id " << best_slave_curve_id
        //         << " best_slave_local_coords " << best_slave_local_coords
        //         << " best_slave_projection " << best_slave_projection
        //         << std::endl;
        // }

        if (!is_active) {
            p_node->SetValue(ACTIVATION_LEVEL, 0.0);
            // KRATOS_INFO("IgaContactProcessSbm") << "Master node inactive at coords " << p_node->Coordinates()
            //                                     << " activation " << p_node->GetValue(ACTIVATION_LEVEL) << std::endl;
            continue;
        }

        std::vector<array_1d<double, 3>> slave_global_derivatives(2, ZeroVector(3));
        mrSlaveSkinModelPart->pGetGeometry(best_slave_curve_id)->GlobalSpaceDerivatives(slave_global_derivatives, best_slave_local_coords, 1);
        CoordinatesArrayType slave_tangent = slave_global_derivatives[1];
        const double slave_tangent_norm = norm_2(slave_tangent);
        if (slave_tangent_norm > 0.0) {
            slave_tangent /= slave_tangent_norm;
        }
        Vector slave_normal = ZeroVector(3);
        slave_normal[0] = slave_tangent[1];
        slave_normal[1] = -slave_tangent[0];

        IndexType slave_brep_id = std::numeric_limits<IndexType>::max();
        double slave_brep_local_parameter = 0.0;
        const auto t_back_begin = Clock::now();
        ProjectBackToSurrogateBoundary(
            *mrSlaveModelPart,
            best_slave_projection,
            mSparseBrepMatrixSlave,
            slave_normal,
            slave_brep_id,
            slave_brep_local_parameter);
        time_back_projection += std::chrono::duration<double>(Clock::now() - t_back_begin).count();

        if (slave_brep_id == std::numeric_limits<IndexType>::max()) {
            p_node->SetValue(ACTIVATION_LEVEL, 0.0);
            KRATOS_INFO("IgaContactProcessSbm") << "Master node inactive (no slave brep) at coords " << p_node->Coordinates()
                                                << " activation " << p_node->GetValue(ACTIVATION_LEVEL) << std::endl;
            continue;
        }

        auto p_slave_brep_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);
        auto p_slave_brep_curve = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_brep_geometry);
        KRATOS_ERROR_IF(!p_slave_brep_curve) << ":::[IgaContactProcessSbm]::: the geometry with id " << slave_brep_id << " is not a Brep." << std::endl;

        CoordinatesArrayType surrogate_local_coords = ZeroVector(3);
        surrogate_local_coords[0] = slave_brep_local_parameter;
        CoordinatesArrayType surrogate_point = ZeroVector(3);
        p_slave_brep_curve->GlobalCoordinates(surrogate_point, surrogate_local_coords);

        for (IndexType i = 0; i < 3; ++i) {
            projection_coordinates[i] = surrogate_point[i];
        }

        p_node->SetValue(ACTIVATION_LEVEL, 1.0);
        p_node->SetValue(PROJECTION_NODE_COORDINATES, projection_coordinates); // FIXME: rename variable once new variable is available
        // KRATOS_INFO("IgaContactProcessSbm") << "Master node active at coords " << p_node->Coordinates()
        //                                     << " activation " << p_node->GetValue(ACTIVATION_LEVEL) << std::endl;
    }

    for (auto& r_segment : master_segments) {
        const double begin_activation = r_segment.p_node_begin->GetValue(ACTIVATION_LEVEL);
        const double end_activation = r_segment.p_node_end->GetValue(ACTIVATION_LEVEL);

        const bool is_begin_active = begin_activation > 0.5;
        const bool is_end_active = end_activation > 0.5;

        if (is_begin_active == is_end_active) {
            continue;
        }

        NodeType::Pointer p_active_node = is_begin_active ? r_segment.p_node_begin : r_segment.p_node_end;
        const Vector& r_projection_vector = p_active_node->GetValue(PROJECTION_NODE_COORDINATES);

        if (r_projection_vector.size() == 0) {
            continue;
        }

        const CoordinatesArrayType projection_point = ConvertVectorToCoordinates(r_projection_vector);

        for (auto& r_slave_geometry : mrSlaveModelPart->Geometries()) {
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(r_slave_geometry.Id());
            auto p_slave_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);
            if (!p_slave_brep) {
                continue;
            }

            Vector slave_domain_interval;
            p_slave_brep->DomainInterval(slave_domain_interval);

            CoordinatesArrayType slave_local_begin = ZeroVector(3);
            CoordinatesArrayType slave_local_end = ZeroVector(3);
            slave_local_begin[0] = slave_domain_interval[0];
            slave_local_end[0] = slave_domain_interval[1];

            CoordinatesArrayType slave_vertex_begin = ZeroVector(3);
            CoordinatesArrayType slave_vertex_end = ZeroVector(3);
            p_slave_brep->GlobalCoordinates(slave_vertex_begin, slave_local_begin);
            p_slave_brep->GlobalCoordinates(slave_vertex_end, slave_local_end);

            if (!IsPointOnAxisAlignedSegment(projection_point, slave_vertex_begin, slave_vertex_end, coordinate_tolerance)) {
                continue;
            }

            Vector brep_projection_vector = ZeroVector(3);
            for (IndexType i = 0; i < std::min<IndexType>(brep_projection_vector.size(), r_projection_vector.size()); ++i) {
                brep_projection_vector[i] = r_projection_vector[i];
            }

            p_slave_brep->SetValue(PROJECTION_NODE_COORDINATES, brep_projection_vector); // FIXME: rename variable once new variable is available
            p_slave_brep->SetValue(ACTIVATION_LEVEL, 2.0);
            break;
        }
    }

    std::vector<CoordinatesArrayType> slave_vertices;
    CollectUniqueSlaveVertices(slave_vertices, coordinate_tolerance);
    SplitMasterSegmentsWithSlaveVertices(
        slave_vertices,
        master_segments,
        master_nodes,
        next_node_id,
        coordinate_tolerance,
        master_layer_name,
        slave_layer_name);

    std::vector<GeometryType::Pointer> contact_geometries;
    std::vector<GeometryType::Pointer> neumann_geometries_master;
    std::vector<GeometryType::Pointer> neumann_geometries_slave;

    for (const auto& r_segment : master_segments) {
        const bool begin_active = r_segment.p_node_begin->GetValue(ACTIVATION_LEVEL) == 1;
        const bool end_active = r_segment.p_node_end->GetValue(ACTIVATION_LEVEL) == 1;

        NurbsInterval segment_interval(r_segment.local_begin, r_segment.local_end);
        auto p_curve_on_surface = r_segment.p_brep->pGetCurveOnSurface();
        auto p_segment_brep = Kratos::make_shared<BrepCurveOnSurfaceType>(p_curve_on_surface, segment_interval, r_segment.p_brep->HasSameCurveDirection());

        GeometriesArrayType quadrature_geometries;
        const std::vector<double> custom_spans{r_segment.local_begin, r_segment.local_end};
        const auto t_quad_begin = Clock::now();
        CreateQuadratureGeometries(*p_segment_brep, quadrature_geometries, &custom_spans);
        time_create_quadrature += std::chrono::duration<double>(Clock::now() - t_quad_begin).count();

        if (!begin_active || !end_active) {
            for (IndexType i = 0; i < quadrature_geometries.size(); ++i) {
                auto p_gp = quadrature_geometries(i);
                // const auto& r_gp_center = p_gp->Center();
                // const bool is_target_gp = (std::abs(r_gp_center[0] + 0.305) < 1.0e-2) &&
                //                           (std::abs(r_gp_center[1] - 36.205) < 1.0e-2);
                // if (is_target_gp) {
                //     KRATOS_WATCH("TARGET_GP_BEGIN_END_INACTIVE")
                //     KRATOS_WATCH(r_gp_center)
                //     KRATOS_WATCH(begin_active)
                //     KRATOS_WATCH(end_active)
                // }

                std::vector<array_1d<double, 3>> master_curve_derivatives(3, ZeroVector(3));
                CoordinatesArrayType master_skin_point = ZeroVector(3);
                std::string projection_layer_name = master_layer_name;
                const auto t_proj_begin = Clock::now();
                bool is_projected_master = ProjectToSkinBoundary(
                    mrMasterSkinModelPart,
                    projection_layer_name,
                    p_gp->Center(),
                    master_skin_point,
                    master_curve_derivatives,
                    10);
                time_project_to_skin += std::chrono::duration<double>(Clock::now() - t_proj_begin).count();

                if (!is_projected_master) {
                    projection_layer_name.clear();
                    const auto t_proj_begin_retry = Clock::now();
                    ProjectToSkinBoundary(
                        mrMasterSkinModelPart,
                        projection_layer_name,
                        p_gp->Center(),
                        master_skin_point,
                        master_curve_derivatives,
                        10);
                    time_project_to_skin += std::chrono::duration<double>(Clock::now() - t_proj_begin_retry).count();
                }
                if (!is_projected_master) {
                    KRATOS_WARNING("") << "::[IgaModelerSbm]:: master surrogate point " << p_gp->Center() << "could not be projected to the skin master. \n"
                                << "The best projection was: " << master_skin_point << std::endl;
                    // if (is_target_gp) {
                    //     KRATOS_WATCH("TARGET_GP_MASTER_PROJECTION_FAILED_INACTIVE_BRANCH")
                    //     KRATOS_WATCH(master_skin_point)
                    // }
                }

                IndexType new_node_id = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size() + 1;
                auto p_master_skin_node = new Node(new_node_id, master_skin_point);
                NodePointerVector neighbour_nodes;

                CoordinatesArrayType tangent_vector = master_curve_derivatives[1];
                const double tangent_norm = norm_2(tangent_vector);
                if (tangent_norm > 0.0) {
                    tangent_vector /= tangent_norm;
                }
                Vector normal_vector = ZeroVector(3);
                normal_vector[0] = tangent_vector[1];
                normal_vector[1] = -tangent_vector[0];
                p_master_skin_node->SetValue(NORMAL, normal_vector);
                

                CoordinatesArrayType curve_first_derivative_vector = master_curve_derivatives[1];
                CoordinatesArrayType curve_second_derivative_vector = master_curve_derivatives[2];
                double curvature = norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector, curve_second_derivative_vector)) / pow(norm_2(curve_first_derivative_vector), 3);
                p_master_skin_node->SetValue(CURVATURE, curvature);
                
                p_master_skin_node->SetValue(LOCAL_TANGENT, tangent_vector);
                
                neighbour_nodes.push_back(p_master_skin_node);
                p_gp->SetValue(NEIGHBOUR_NODES, neighbour_nodes);

                neumann_geometries_master.push_back(p_gp);
                // if (is_target_gp) {
                //     KRATOS_WATCH("TARGET_GP_ASSIGNED_NEUMANN_INACTIVE_BRANCH")
                // }
            }
            continue;
        }

        for (IndexType i = 0; i < quadrature_geometries.size(); ++i) {
            auto p_gp = quadrature_geometries(i);
            // const auto& r_gp_center = p_gp->Center();
            // const bool is_target_gp = (std::abs(r_gp_center[0] + 0.305) < 1.0e-2) &&
            //                           (std::abs(r_gp_center[1] - 36.205) < 1.0e-2);
            // if (is_target_gp) {
            //     KRATOS_WATCH("TARGET_GP_ACTIVE_BRANCH_START")
            //     KRATOS_WATCH(r_gp_center)
            // }

            std::vector<array_1d<double, 3>> master_curve_derivatives(3, ZeroVector(3));
            CoordinatesArrayType master_skin_point = ZeroVector(3);
            std::string projection_layer_name = master_layer_name;
            const auto t_proj_begin = Clock::now();
            bool is_projected_master = ProjectToSkinBoundary(
                mrMasterSkinModelPart,
                projection_layer_name,
                p_gp->Center(),
                master_skin_point,
                master_curve_derivatives,
                10);
            time_project_to_skin += std::chrono::duration<double>(Clock::now() - t_proj_begin).count();

            // if (is_target_gp) {
            //     KRATOS_WATCH("TARGET_GP_BEGIN_END_INACTIVE")
            //     KRATOS_WATCH(r_gp_center)
            //     KRATOS_WATCH(master_skin_point)
            //     // exit(0);
            // }

            if (!is_projected_master) {
                projection_layer_name.clear();
                const auto t_proj_begin_retry = Clock::now();
                is_projected_master = ProjectToSkinBoundary(
                    mrMasterSkinModelPart,
                    projection_layer_name,
                    p_gp->Center(),
                    master_skin_point,
                    master_curve_derivatives,
                    10);
                time_project_to_skin += std::chrono::duration<double>(Clock::now() - t_proj_begin_retry).count();
            }

            if (!is_projected_master && norm_2(master_skin_point - p_gp->Center()) > projection_distance_fallback) {
                neumann_geometries_master.push_back(p_gp);
                // if (is_target_gp) {
                //     KRATOS_WATCH(projection_distance_fallback)
                //     KRATOS_WATCH(norm_2(master_skin_point - p_gp->Center()))
                //     KRATOS_WATCH("TARGET_GP_ASSIGNED_NEUMANN_MASTER_PROJECTION_FAILED")
                // }
                continue;
            }

            IndexType master_skin_node_id = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size() + 1;
            auto p_master_skin_node = new Node(master_skin_node_id, master_skin_point);

            CoordinatesArrayType master_tangent = master_curve_derivatives[1];
            const double master_tangent_norm = norm_2(master_tangent);
            if (master_tangent_norm > 0.0) {
                master_tangent /= master_tangent_norm;
            }
            Vector master_normal = ZeroVector(3);
            master_normal[0] = master_tangent[1];
            master_normal[1] = -master_tangent[0];
            p_master_skin_node->SetValue(NORMAL, master_normal);
            p_master_skin_node->SetValue(LOCAL_TANGENT, master_tangent);

            CoordinatesArrayType curve_first_derivative_vector = master_curve_derivatives[1];
            CoordinatesArrayType curve_second_derivative_vector = master_curve_derivatives[2];
            double curvature = norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector, curve_second_derivative_vector)) / pow(norm_2(curve_first_derivative_vector), 3);
            p_master_skin_node->SetValue(CURVATURE, curvature);

            NodePointerVector master_neighbours;
            master_neighbours.push_back(p_master_skin_node);
            p_gp->SetValue(NEIGHBOUR_NODES, master_neighbours);

            CoordinatesArrayType master_skin_deformed = ZeroVector(3);
            GetDeformedPosition(master_skin_point, *mrMasterModelPart, mSparseBrepMatrixMaster, master_skin_deformed);

            bool best_projection_found = false;
            double best_projection_distance = std::numeric_limits<double>::max();
            CoordinatesArrayType best_slave_projection = ZeroVector(3);
            CoordinatesArrayType best_slave_local_coords = ZeroVector(3);
            IndexType best_slave_curve_id = 0;

            for (auto& r_slave_curve_geometry : mrSlaveSkinModelPart->Geometries()) {
                if (r_slave_curve_geometry.GetValue(IDENTIFIER) != slave_layer_name) {
                    continue;
                }

                auto p_slave_curve_geometry = mrSlaveSkinModelPart->pGetGeometry(r_slave_curve_geometry.Id());
                auto p_slave_nurbs_curve = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_slave_curve_geometry);
                KRATOS_ERROR_IF(!p_slave_nurbs_curve) <<  ":::[IgaContactProcessSbm]::: the geometry with id "
                    << r_slave_curve_geometry.Id() << " is not a NurbsCurveGeometryType." << std::endl;

                CoordinatesArrayType local_coords = ZeroVector(3);
                CoordinatesArrayType projected_point = ZeroVector(3);
                double distance = std::numeric_limits<double>::max();

                const auto t_newton_begin = Clock::now();
                bool is_converged = NewtonRaphsonCurveOnDeformed(
                    local_coords,
                    master_skin_deformed,
                    projected_point,
                    *p_slave_nurbs_curve,
                    distance,
                    20,
                    10,
                    1e-9
                );
                time_newton += std::chrono::duration<double>(Clock::now() - t_newton_begin).count();
                // if (is_target_gp) {
                //     KRATOS_WATCH("TARGET_GP_SLAVE_CURVE_TRY")
                //     KRATOS_WATCH(r_slave_curve_geometry.Id())
                //     KRATOS_WATCH(master_skin_point)
                //     KRATOS_WATCH(master_skin_deformed)
                //     KRATOS_WATCH(projected_point)
                //     KRATOS_WATCH(slave_layer_name)
                //     KRATOS_WATCH(is_converged)
                //     KRATOS_WATCH(distance)
                // }

                if ((is_converged || distance < projection_distance_fallback) && distance < best_projection_distance) {
                    best_projection_found = true;
                    best_projection_distance = distance;
                    best_slave_local_coords = local_coords;
                    best_slave_curve_id = r_slave_curve_geometry.Id();

                    std::vector<array_1d<double, 3>> derivatives(2, ZeroVector(3));
                    p_slave_nurbs_curve->GlobalSpaceDerivatives(derivatives, local_coords, 1);
                    best_slave_projection = derivatives[0];
                }
            }

            if (!(best_projection_found &&
                (best_projection_distance < projection_distance_limit || best_projection_distance < projection_distance_fallback))) {
                neumann_geometries_master.push_back(p_gp);
                // if (is_target_gp) {
                //     KRATOS_WATCH("TARGET_GP_ASSIGNED_NEUMANN_SLAVE_PROJECTION_FAILED_OR_DISTANCE")
                //     KRATOS_WATCH(best_projection_found)
                //     KRATOS_WATCH(best_projection_distance)
                //     KRATOS_WATCH(projection_distance_limit)
                //     KRATOS_WATCH(projection_distance_fallback)
                // }
                continue;
            }

            IndexType slave_skin_node_id = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size() + 1;
            auto p_slave_skin_node = new Node(slave_skin_node_id, best_slave_projection);
            NodePointerVector slave_neighbours;
            slave_neighbours.push_back(p_slave_skin_node);
            p_master_skin_node->SetValue(NEIGHBOUR_NODES, slave_neighbours);

            std::vector<array_1d<double, 3>> slave_global_derivatives(2, ZeroVector(3));
            mrSlaveSkinModelPart->pGetGeometry(best_slave_curve_id)->GlobalSpaceDerivatives(slave_global_derivatives, best_slave_local_coords, 2);
            CoordinatesArrayType slave_tangent = slave_global_derivatives[1];
            const double slave_tangent_norm = norm_2(slave_tangent);
            if (slave_tangent_norm > 0.0) {
                slave_tangent /= slave_tangent_norm;
            }
            Vector slave_normal = ZeroVector(3);
            slave_normal[0] = slave_tangent[1];
            slave_normal[1] = -slave_tangent[0];
            p_slave_skin_node->SetValue(NORMAL, slave_normal);
            p_slave_skin_node->SetValue(LOCAL_TANGENT, slave_tangent);

            curve_first_derivative_vector = slave_global_derivatives[1];
            curve_second_derivative_vector = slave_global_derivatives[2];
            curvature = norm_2(MathUtils<double>::CrossProduct(curve_first_derivative_vector, curve_second_derivative_vector)) / pow(norm_2(curve_first_derivative_vector), 3);
            p_slave_skin_node->SetValue(CURVATURE, curvature);

            IndexType slave_brep_id = std::numeric_limits<IndexType>::max();
            double slave_brep_local_parameter = 0.0;
            const auto t_back_begin = Clock::now();
            ProjectBackToSurrogateBoundary(
                *mrSlaveModelPart,
                best_slave_projection,
                mSparseBrepMatrixSlave,
                slave_normal,
                slave_brep_id,
                slave_brep_local_parameter);
            time_back_projection += std::chrono::duration<double>(Clock::now() - t_back_begin).count();

            if (slave_brep_id == std::numeric_limits<IndexType>::max()) {
                neumann_geometries_master.push_back(p_gp);
                // if (is_target_gp) {
                //     KRATOS_WATCH("TARGET_GP_ASSIGNED_NEUMANN_NO_SLAVE_BREP")
                // }
                continue;
            }

            auto p_slave_brep_geometry = mrSlaveModelPart->pGetGeometry(slave_brep_id);
            auto p_slave_brep_curve = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_brep_geometry);
            KRATOS_ERROR_IF(!p_slave_brep_curve) <<  ":::[IgaContactProcessSbm]::: the geometry with id "
                << slave_brep_id << " is not a Brep." << std::endl;

            IntegrationPoint<1> surrogate_integration_point(slave_brep_local_parameter);
            IntegrationPointsArrayType surrogate_integration_points_list;
            surrogate_integration_points_list.push_back(surrogate_integration_point);
            GeometriesArrayType surrogate_quadrature_list;

            const int number_of_shape_function_derivatives = 5; //FIXME:
            IntegrationInfo surrogate_integration_info = p_slave_brep_curve->GetDefaultIntegrationInfo();
            const auto t_quad_sbm_begin = Clock::now();
            p_slave_brep_curve->pGetCurveOnSurface()->CreateQuadraturePointGeometriesSBM(
                surrogate_quadrature_list,
                number_of_shape_function_derivatives,
                surrogate_integration_points_list,
                surrogate_integration_info);
            time_create_quadrature_sbm += std::chrono::duration<double>(Clock::now() - t_quad_sbm_begin).count();

            // FIXME: is it really necessary?
            std::vector<Geometry<Node>::Pointer> neighbour_geometries;
            neighbour_geometries.push_back(surrogate_quadrature_list(0));
            p_slave_skin_node->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);

            const auto t_coupling_begin = Clock::now();
            auto p_coupling_geometry = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCouplingGeometry2D(
                p_gp,
                surrogate_quadrature_list(0));
            time_create_coupling += std::chrono::duration<double>(Clock::now() - t_coupling_begin).count();
            p_coupling_geometry->SetValue(KNOT_SPAN_SIZES, master_knot_step_uv);
                contact_geometries.push_back(p_coupling_geometry);
                // if (is_target_gp) {
                //     KRATOS_WATCH("TARGET_GP_ASSIGNED_CONTACT")
                // }

                auto p_slave_brep_for_contact = mrSlaveModelPart->pGetGeometry(slave_brep_id);
                if (p_slave_brep_for_contact->GetValue(ACTIVATION_LEVEL) != 2.0) {
                    p_slave_brep_for_contact->SetValue(ACTIVATION_LEVEL, 1.0);
                }
        }
    }

    std::vector<CoordinatesArrayType> level_one_vertices;
    const double level_one_tolerance = 1.0e-6;
    for (auto& r_slave_geometry : mrSlaveModelPart->Geometries()) {
        auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(r_slave_geometry.Id());
        auto p_slave_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);

        if (std::abs(p_slave_brep->GetValue(ACTIVATION_LEVEL) - 1.0) < 1.0e-12) {
            auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(r_slave_geometry.Id());
            auto p_slave_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);
            if (!p_slave_brep) {
                continue;
            }

            Vector slave_domain_interval;
            p_slave_brep->DomainInterval(slave_domain_interval);

            CoordinatesArrayType local_begin = ZeroVector(3);
            CoordinatesArrayType local_end = ZeroVector(3);
            local_begin[0] = slave_domain_interval[0];
            local_end[0] = slave_domain_interval[1];

            CoordinatesArrayType vertex_begin = ZeroVector(3);
            CoordinatesArrayType vertex_end = ZeroVector(3);
            p_slave_brep->GlobalCoordinates(vertex_begin, local_begin);
            p_slave_brep->GlobalCoordinates(vertex_end, local_end);

            AppendUniqueVertex(level_one_vertices, vertex_begin, level_one_tolerance);
            AppendUniqueVertex(level_one_vertices, vertex_end, level_one_tolerance);
        }
    }

    for (auto& r_slave_geometry : mrSlaveModelPart->Geometries()) {
        auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(r_slave_geometry.Id());
        auto p_slave_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);
        KRATOS_ERROR_IF(!p_slave_brep) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << r_slave_geometry.Id() << " is not a Brep." << std::endl;

        if (p_slave_brep->GetValue(ACTIVATION_LEVEL) == 0.0) {
            GeometriesArrayType surrogate_quadrature;
            // CreateQuadratureGeometries(*p_slave_brep, surrogate_quadrature); //FIXME:
            p_slave_brep->GetQuadraturePointGeometries(surrogate_quadrature);
            for (IndexType i = 0; i < surrogate_quadrature.size(); ++i) {
                // SetSurrogateNeighbourNodes(surrogate_quadrature(i), p_slave_brep, slave_layer_name);
                neumann_geometries_slave.push_back(surrogate_quadrature(i)); //FIXME:
            }
        } else if (p_slave_brep->GetValue(ACTIVATION_LEVEL) == 1.0) {
            // nothing to do, already included in the contact part
        } else if (p_slave_brep->GetValue(ACTIVATION_LEVEL) == 2.0) {
            const Vector& r_projection_vector = p_slave_brep->GetValue(PROJECTION_NODE_COORDINATES);
            if (r_projection_vector.size() == 0) {
                KRATOS_ERROR << "::[IgaContactProcess]:: projection vector not defined for slave brep" << std::endl;
            }


            const CoordinatesArrayType projection_point = ConvertVectorToCoordinates(r_projection_vector);

            Vector slave_domain_interval;
            p_slave_brep->DomainInterval(slave_domain_interval);

            CoordinatesArrayType local_begin = ZeroVector(3);
            CoordinatesArrayType local_end = ZeroVector(3);
            local_begin[0] = slave_domain_interval[0];
            local_end[0] = slave_domain_interval[1];

            CoordinatesArrayType vertex_begin = ZeroVector(3);
            CoordinatesArrayType vertex_end = ZeroVector(3);
            p_slave_brep->GlobalCoordinates(vertex_begin, local_begin);
            p_slave_brep->GlobalCoordinates(vertex_end, local_end);

            const bool begin_is_active = ContainsVertex(level_one_vertices, vertex_begin, level_one_tolerance);
            const bool end_is_active = ContainsVertex(level_one_vertices, vertex_end, level_one_tolerance);

            CoordinatesArrayType target_vertex = begin_is_active && !end_is_active ? vertex_end : vertex_begin;
            double target_parameter = begin_is_active && !end_is_active ? slave_domain_interval[1] : slave_domain_interval[0];
            if (begin_is_active && end_is_active) {
                KRATOS_WARNING("::[IgaContactProcessSbm]:: both vertices active in cut slave brep") << std::endl;
                continue;
            }

            CoordinatesArrayType projection_local = ZeroVector(3);

            int is_projected_on_local_space = p_slave_brep->ProjectionPointGlobalToLocalSpace(projection_point, projection_local);

            const double parameter_projection = projection_local[0];
            const double parameter_target = target_parameter;

            if (std::abs(parameter_projection - parameter_target) < 1.0e-8) {
                continue;
            }

            const double trimmed_lower = std::min(parameter_projection, parameter_target);
            const double trimmed_upper = std::max(parameter_projection, parameter_target);

            NurbsInterval trimmed_interval(trimmed_lower, trimmed_upper);

            auto p_trimmed_brep = Kratos::make_shared<BrepCurveOnSurfaceType>(
                p_slave_brep->pGetCurveOnSurface(),
                trimmed_interval,
                p_slave_brep->HasSameCurveDirection());

            GeometriesArrayType surrogate_quadrature;
            const std::vector<double> custom_spans{trimmed_lower, trimmed_upper};
            CreateQuadratureGeometries(*p_trimmed_brep, surrogate_quadrature, &custom_spans);
            for (IndexType i = 0; i < surrogate_quadrature.size(); ++i) {
                SetSurrogateNeighbourNodes(surrogate_quadrature(i), p_trimmed_brep, slave_layer_name);
                neumann_geometries_slave.push_back(surrogate_quadrature(i));
            }
        }
    }

    GeometriesArrayType contact_geometries_array(contact_geometries.size());
    for (IndexType i = 0; i < contact_geometries.size(); ++i) {
        contact_geometries_array(i) = contact_geometries[i];
    }
    
    GeometriesArrayType neumann_geometries_array_master(neumann_geometries_master.size());
    for (IndexType i = 0; i < neumann_geometries_master.size(); ++i) {
        neumann_geometries_array_master(i) = neumann_geometries_master[i];
    }

    GeometriesArrayType neumann_geometries_array_slave(neumann_geometries_slave.size());
    for (IndexType i = 0; i < neumann_geometries_slave.size(); ++i) {
        neumann_geometries_array_slave(i) = neumann_geometries_slave[i];
    }
    this->CreateConditions(
        contact_geometries_array.ptr_begin(),
        contact_geometries_array.ptr_end(),
        mrContactModelPart->GetSubModelPart("contact"),
        name,
        next_condition_id,
        mpPropMaster,
        mpPropSlave);

    const std::string default_condition_name = "SbmLoadSolidCondition";
    const Condition& r_reference_condition = KratosComponents<Condition>::Get(default_condition_name);
    Vector mesh_sizes_uv = mrSlaveModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    ModelPart::ConditionsContainerType new_condition_list_master;
    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << default_condition_name
        << " in " << mrContactModelPart->GetParentModelPart().Name() << ".master-SubModelPart." << std::endl;

    IndexType count_cond = 0;
    for (auto it = neumann_geometries_array_master.ptr_begin(); it != neumann_geometries_array_master.ptr_end(); ++it) {
        new_condition_list_master.push_back(r_reference_condition.Create(next_condition_id++, (*it), mpPropMaster));
        new_condition_list_master.GetContainer()[count_cond]->SetValue(KNOT_SPAN_SIZES, mesh_sizes_uv);
        new_condition_list_master.GetContainer()[count_cond]->SetValue(IDENTIFIER, "outer");
        ++count_cond;
    }

    (mrContactModelPart->GetSubModelPart("master")).AddConditions(new_condition_list_master.begin(), new_condition_list_master.end());
    
    ModelPart::ConditionsContainerType new_condition_list_slave;
    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << default_condition_name
        << " in " << mrContactModelPart->GetParentModelPart().Name() << ".slave-SubModelPart." << std::endl;


    count_cond = 0;
    for (auto it = neumann_geometries_array_slave.ptr_begin(); it != neumann_geometries_array_slave.ptr_end(); ++it) {
        new_condition_list_slave.push_back(r_reference_condition.Create(next_condition_id++, (*it), mpPropSlave));
        new_condition_list_slave.GetContainer()[count_cond]->SetValue(KNOT_SPAN_SIZES, mesh_sizes_uv);
        new_condition_list_slave.GetContainer()[count_cond]->SetValue(IDENTIFIER, "outer");
        ++count_cond;
    }

    (mrContactModelPart->GetSubModelPart("slave")).AddConditions(new_condition_list_slave.begin(), new_condition_list_slave.end());
    
    EntitiesUtilities::InitializeEntities<Condition>(mrContactModelPart->GetParentModelPart());

    KRATOS_ERROR_IF(mrContactModelPart->NumberOfConditions() == 0) << "YOUR CONTACT MODEL PART IS EMPTY" << std::endl;


    // Obtain the slave skin model part //TODO: read from input
    std::string slave_skin_model_part_name = "skin_Body1.outer.ContactSide";
    std::string master_skin_model_part_name = "skin_Body2.outer.ContactSide";

    // std::string master_skin_model_part_name = "skin_Body1.outer.Bottom";
    // std::string slave_skin_model_part_name = "skin_Body2.outer.Top_1";

    // std::string master_skin_model_part_name = "skin_Body1.outer.bottom";
    // std::string slave_skin_model_part_name = "skin_Body2.outer.top";

    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_skin_model_part_name)) << "ERROR: SLAVE SKIN MODEL PART " 
                                            << slave_skin_model_part_name << "NOT CREATED" << std::endl; 

    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_skin_model_part_name)) << "ERROR: MASTER SKIN MODEL PART " 
                                            << master_skin_model_part_name << "NOT CREATED" << std::endl; 

    auto& r_slave_skin_sub_model_part = (mpModel->GetModelPart(slave_skin_model_part_name));
    auto& r_master_skin_sub_model_part = (mpModel->GetModelPart(master_skin_model_part_name));

    KRATOS_INFO("IgaContactProcessSbmTiming") << "ProjectToSkinBoundary: " << time_project_to_skin
        << "s, NewtonRaphsonCurveOnDeformed: " << time_newton
        << "s, ProjectBackToSurrogateBoundary: " << time_back_projection
        << "s, CreateQuadratureGeometries: " << time_create_quadrature
        << "s, CreateQuadraturePointGeometriesSBM: " << time_create_quadrature_sbm
        << "s, CreateQuadraturePointCouplingGeometry2D: " << time_create_coupling
        << "s" << std::endl;

    KRATOS_INFO("[IgaContactProcess]:: finished projections \n");

    if (mIntegrateOnTrueBoundary)
        PrepareIntegrationOnTrueBoundary(r_master_skin_sub_model_part, r_slave_skin_sub_model_part);
    
    KRATOS_INFO("[IgaContactProcess]:: finished preparation for integration on true \n");
    
}






    Node::Pointer IgaContactProcessSbm::FindExistingNode(
        const std::vector<NodeType::Pointer>& rNodes,
        const CoordinatesArrayType& rCoordinates,
        double tolerance) const
    {
        const double tolerance_squared = tolerance * tolerance;
        for (NodeType::Pointer p_candidate : rNodes) {
            double distance_squared = 0.0;
            for (IndexType i = 0; i < 3; ++i) {
                const double diff = p_candidate->Coordinates()[i] - rCoordinates[i];
                distance_squared += diff * diff;
            }
            if (distance_squared <= tolerance_squared) {
                return p_candidate;
            }
        }
        return nullptr;
    }

    Node::Pointer IgaContactProcessSbm::CreateOrRetrieveContactNode(
        const CoordinatesArrayType& rCoordinates,
        std::vector<NodeType::Pointer>& rNodes,
        IndexType& rNextNodeId,
        double tolerance)
    {
        if (NodeType::Pointer p_existing = FindExistingNode(rNodes, rCoordinates, tolerance); p_existing != nullptr) {
            return p_existing;
        }

        auto p_new_node = mrContactModelPart->CreateNewNode(rNextNodeId++, rCoordinates[0], rCoordinates[1], rCoordinates[2]);
        p_new_node->SetValue(ACTIVATION_LEVEL, 0.0);
        Vector projection_coordinates = ZeroVector(3);
        p_new_node->SetValue(PROJECTION_NODE_COORDINATES, projection_coordinates);

        rNodes.push_back(p_new_node);
        return p_new_node;
    }

IgaContactProcessSbm::CoordinatesArrayType IgaContactProcessSbm::ConvertVectorToCoordinates(const Vector& rVector)
{
    CoordinatesArrayType result = ZeroVector(3);
    const IndexType size = std::min<IndexType>(rVector.size(), 3);
    for (IndexType i = 0; i < size; ++i) {
        result[i] = rVector[i];
    }
    return result;
}

bool IgaContactProcessSbm::IsPointOnAxisAlignedSegment(
    const CoordinatesArrayType& rPoint,
    const CoordinatesArrayType& rSegmentBegin,
    const CoordinatesArrayType& rSegmentEnd,
    double tolerance)
{
    int varying_axis = -1;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(rSegmentBegin[i] - rSegmentEnd[i]) > tolerance) {
            if (varying_axis != -1) {
                return false;
            }
            varying_axis = i;
        }
    }

    if (varying_axis == -1) {
        for (int i = 0; i < 3; ++i) {
            if (std::abs(rPoint[i] - rSegmentBegin[i]) > tolerance) {
                return false;
            }
        }
        return true;
    }

    for (int i = 0; i < 3; ++i) {
        if (i == varying_axis) {
            continue;
        }
        if (std::abs(rPoint[i] - rSegmentBegin[i]) > tolerance) {
            return false;
        }
    }

    const double min_value = std::min(rSegmentBegin[varying_axis], rSegmentEnd[varying_axis]) - tolerance;
    const double max_value = std::max(rSegmentBegin[varying_axis], rSegmentEnd[varying_axis]) + tolerance;
    return (rPoint[varying_axis] >= min_value) && (rPoint[varying_axis] <= max_value);
}

void IgaContactProcessSbm::CreateQuadratureGeometries(
    BrepCurveOnSurfaceType& rBrep,
    GeometriesArrayType& rGeometries,
    const std::vector<double>* pCustomSpans)
{
    SizeType shape_function_derivatives_order = 3;
    if (mParameters.Has("shape_function_derivatives_order")) {
        shape_function_derivatives_order = mParameters["shape_function_derivatives_order"].GetInt();
    } else {
        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
            << "shape_function_derivatives_order is not provided and thus being considered as 3. " << std::endl;
    }

    std::string quadrature_method = mParameters.Has("quadrature_method")
        ? mParameters["integration_rule"].GetString()
        : "GAUSS";

    IntegrationInfo integration_info = rBrep.GetDefaultIntegrationInfo();
    if (mParameters.Has("number_of_integration_points_per_span")) {
        for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
            integration_info.SetNumberOfIntegrationPointsPerSpan(i, mParameters["number_of_integration_points_per_span"].GetInt());
        }
    }

    for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
        if (quadrature_method == "GAUSS") {
            integration_info.SetQuadratureMethod(i, IntegrationInfo::QuadratureMethod::GAUSS);
        } else if (quadrature_method == "GRID") {
            integration_info.SetQuadratureMethod(i, IntegrationInfo::QuadratureMethod::GRID);
        } else {
            KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
        }
    }

    std::vector<double> spans;
    if (pCustomSpans) {
        spans = *pCustomSpans;
    } else {
        rBrep.SpansLocalSpace(spans);
        std::sort(spans.begin(), spans.end());
        spans.erase(std::unique(spans.begin(), spans.end(), [](double a, double b) {
            return std::abs(a - b) < 1e-4;
        }), spans.end());
    }

    IntegrationPointsArrayType integration_points;
    IntegrationPointUtilities::CreateIntegrationPoints1D(integration_points, spans, integration_info);

    rBrep.CreateQuadraturePointGeometries(rGeometries, shape_function_derivatives_order, integration_points, integration_info, false);
}

void IgaContactProcessSbm::CollectUniqueSlaveVertices(
    std::vector<CoordinatesArrayType>& rVertices,
    double tolerance) const
{
    rVertices.clear();
    rVertices.reserve(mrSlaveModelPart->NumberOfGeometries() * 2);

    for (auto& r_slave_geometry : mrSlaveModelPart->Geometries()) {
        auto p_slave_geometry = mrSlaveModelPart->pGetGeometry(r_slave_geometry.Id());
        auto p_slave_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_slave_geometry);
        if (!p_slave_brep) {
            continue;
        }

        Vector slave_domain_interval;
        p_slave_brep->DomainInterval(slave_domain_interval);

        CoordinatesArrayType local_begin = ZeroVector(3);
        CoordinatesArrayType local_end = ZeroVector(3);
        local_begin[0] = slave_domain_interval[0];
        local_end[0] = slave_domain_interval[1];

        CoordinatesArrayType vertex_begin = ZeroVector(3);
        CoordinatesArrayType vertex_end = ZeroVector(3);
        p_slave_brep->GlobalCoordinates(vertex_begin, local_begin);
        p_slave_brep->GlobalCoordinates(vertex_end, local_end);

        AppendUniqueVertex(rVertices, vertex_begin, tolerance);
        AppendUniqueVertex(rVertices, vertex_end, tolerance);
    }
}

void IgaContactProcessSbm::SplitMasterSegmentsWithSlaveVertices(
    const std::vector<CoordinatesArrayType>& rSlaveVertices,
    std::vector<MasterSegmentData>& rMasterSegments,
    std::vector<NodeType::Pointer>& rMasterNodes,
    IndexType& rNextNodeId,
    double coordinate_tolerance,
    const std::string& master_layer_name,
    const std::string& slave_layer_name)
{
    if (rSlaveVertices.empty()) {
        return;
    }

    const Vector& master_knot_span_sizes = mrMasterModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    const double projection_distance_limit = master_knot_span_sizes[0]*1;
    const double projection_distance_fallback = master_knot_span_sizes[0]/2;

    for (const auto& r_slave_vertex : rSlaveVertices) {
        CoordinatesArrayType slave_skin_projection = ZeroVector(3);
        std::vector<array_1d<double, 3>> slave_curve_derivatives(3, ZeroVector(3));
        std::string projection_layer_name = slave_layer_name;

        bool is_projected_slave = ProjectToSkinBoundary(
            mrSlaveSkinModelPart,
            projection_layer_name,
            r_slave_vertex,
            slave_skin_projection,
            slave_curve_derivatives,
            10);

        if (!is_projected_slave) {
            projection_layer_name.clear();
            is_projected_slave = ProjectToSkinBoundary(
                mrSlaveSkinModelPart,
                projection_layer_name,
                r_slave_vertex,
                slave_skin_projection,
                slave_curve_derivatives,
                10);
        }

        if (!is_projected_slave) {
            continue;
        }

        CoordinatesArrayType slave_deformed = ZeroVector(3);
        GetDeformedPosition(slave_skin_projection, *mrSlaveModelPart, mSparseBrepMatrixSlave, slave_deformed);

        Matrix slave_gradient = ZeroMatrix(2, 2);
        Matrix slave_hessian = ZeroMatrix(2, 3);
        GetDeformedGradient(slave_skin_projection, *mrSlaveModelPart, mSparseBrepMatrixSlave, slave_gradient, slave_hessian);

        CoordinatesArrayType slave_tangent_deformed = ZeroVector(3);
        const CoordinatesArrayType& slave_tangent = slave_curve_derivatives[1];
        for (int i_dim = 0; i_dim < 2; ++i_dim) {
            slave_tangent_deformed[i_dim] = slave_tangent[i_dim]
                + slave_gradient(i_dim, 0) * slave_tangent[0]
                + slave_gradient(i_dim, 1) * slave_tangent[1];
        }

        const double tangent_norm = norm_2(slave_tangent_deformed);
        if (tangent_norm < std::numeric_limits<double>::epsilon()) {
            continue;
        }
        slave_tangent_deformed /= tangent_norm;

        CoordinatesArrayType slave_normal_deformed = ZeroVector(3);
        slave_normal_deformed[0] = -slave_tangent_deformed[1];
        slave_normal_deformed[1] = slave_tangent_deformed[0];

        bool best_projection_found = false;
        double best_projection_distance = std::numeric_limits<double>::max();
        CoordinatesArrayType best_master_local_coords = ZeroVector(3);
        NurbsCurveGeometryType::Pointer p_best_master_curve = nullptr;

        for (auto& r_master_curve_geometry : mrMasterSkinModelPart->Geometries()) {
            if (r_master_curve_geometry.GetValue(IDENTIFIER) != master_layer_name) {
                continue;
            }

            auto p_master_curve_geometry = mrMasterSkinModelPart->pGetGeometry(r_master_curve_geometry.Id());
            auto p_master_curve = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_master_curve_geometry);
            KRATOS_ERROR_IF(!p_master_curve) <<  ":::[IgaContactProcessSbm]::: the geometry with id "
                << r_master_curve_geometry.Id() << " is not a NurbsCurveGeometryType." << std::endl;

            CoordinatesArrayType master_projection_deformed = ZeroVector(3);
            CoordinatesArrayType master_local_coords = ZeroVector(3);
            double projection_distance = std::numeric_limits<double>::max();

            bool converged = ProjectPointViaRayTracingToMasterCurveDeformed(
                slave_deformed,
                slave_normal_deformed,
                *p_master_curve,
                *mrMasterModelPart,
                mSparseBrepMatrixMaster,
                master_projection_deformed,
                master_local_coords,
                projection_distance);

            if (!converged) {
                continue;
            }

            if (projection_distance < best_projection_distance) {
                best_projection_found = true;
                best_projection_distance = projection_distance;
                best_master_local_coords = master_local_coords;
                p_best_master_curve = p_master_curve;
            }
        }

        if (!best_projection_found || p_best_master_curve == nullptr ||
            (best_projection_distance >= projection_distance_limit &&
             best_projection_distance >= projection_distance_fallback)) {
            continue;
        }

        CoordinatesArrayType master_skin_projection = ZeroVector(3);
        p_best_master_curve->GlobalCoordinates(master_skin_projection, best_master_local_coords);

        std::vector<array_1d<double, 3>> master_global_derivatives(2, ZeroVector(3));
        p_best_master_curve->GlobalSpaceDerivatives(master_global_derivatives, best_master_local_coords, 1);
        CoordinatesArrayType master_tangent = master_global_derivatives[1];
        const double master_tangent_norm = norm_2(master_tangent);
        if (master_tangent_norm < std::numeric_limits<double>::epsilon()) {
            continue;
        }
        master_tangent /= master_tangent_norm;

        Vector master_normal = ZeroVector(3);
        master_normal[0] = master_tangent[1];
        master_normal[1] = -master_tangent[0];

        IndexType master_brep_id = std::numeric_limits<IndexType>::max();
        double master_brep_local_parameter = 0.0;
        bool is_projected_back = ProjectBackToSurrogateBoundary(
            *mrMasterModelPart,
            master_skin_projection,
            mSparseBrepMatrixMaster,
            master_normal,
            master_brep_id,
            master_brep_local_parameter);

        if (!is_projected_back || master_brep_id == std::numeric_limits<IndexType>::max()) {
            continue;
        }

        auto p_master_brep_geometry = mrMasterModelPart->pGetGeometry(master_brep_id);
        auto p_master_brep = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_master_brep_geometry);
        KRATOS_ERROR_IF(!p_master_brep) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << master_brep_id << " is not a Brep." << std::endl;

        CoordinatesArrayType split_local_coords = ZeroVector(3);
        split_local_coords[0] = master_brep_local_parameter;
        CoordinatesArrayType split_global_coords = ZeroVector(3);
        p_master_brep->GlobalCoordinates(split_global_coords, split_local_coords);

        NodeType::Pointer p_split_node = CreateOrRetrieveContactNode(split_global_coords, rMasterNodes, rNextNodeId, coordinate_tolerance);

        SplitMasterSegment(master_brep_id, master_brep_local_parameter, p_split_node, rMasterSegments, coordinate_tolerance);
    }
}

void IgaContactProcessSbm::SplitMasterSegment(
    IndexType brep_id,
    double split_parameter,
    NodeType::Pointer p_split_node,
    std::vector<MasterSegmentData>& rMasterSegments,
    double tolerance)
{
    for (SizeType i = 0; i < rMasterSegments.size(); ++i) {
        auto& r_segment = rMasterSegments[i];
        if (!r_segment.p_brep || r_segment.p_brep->Id() != brep_id) {
            continue;
        }

        const double param_min = std::min(r_segment.local_begin, r_segment.local_end);
        const double param_max = std::max(r_segment.local_begin, r_segment.local_end);
        if (split_parameter <= param_min + tolerance || split_parameter >= param_max - tolerance) {
            continue;
        }

        const double orientation = (r_segment.local_end >= r_segment.local_begin) ? 1.0 : -1.0;
        if ((split_parameter - r_segment.local_begin) * orientation <= tolerance ||
            (r_segment.local_end - split_parameter) * orientation <= tolerance) {
            continue;
        }

        // FIXME: NEW TO CHECK
        if (r_segment.p_node_begin->GetValue(ACTIVATION_LEVEL)== 0 || r_segment.p_node_end->GetValue(ACTIVATION_LEVEL)== 0) return;
        p_split_node ->SetValue(ACTIVATION_LEVEL, 1);

        MasterSegmentData new_segment = r_segment;
        NodeType::Pointer p_original_end = r_segment.p_node_end;
        const double original_end_parameter = r_segment.local_end;

        r_segment.local_end = split_parameter;
        r_segment.p_node_end = p_split_node;

        new_segment.local_begin = split_parameter;
        new_segment.local_end = original_end_parameter;
        new_segment.p_node_begin = p_split_node;
        new_segment.p_node_end = p_original_end;

        rMasterSegments.insert(rMasterSegments.begin() + i + 1, new_segment);
        break;
    }
}

void IgaContactProcessSbm::AppendUniqueVertex(
    std::vector<CoordinatesArrayType>& rVertices,
    const CoordinatesArrayType& rVertex,
    double tolerance)
{
    if (!ContainsVertex(rVertices, rVertex, tolerance)) {
        rVertices.push_back(rVertex);
    }
}

bool IgaContactProcessSbm::ContainsVertex(
    const std::vector<CoordinatesArrayType>& rVertices,
    const CoordinatesArrayType& rVertex,
    double tolerance)
{
    for (const auto& r_existing : rVertices) {
        if (norm_2(r_existing - rVertex) < tolerance) {
            return true;
        }
    }
    return false;
}

void IgaContactProcessSbm::SetSurrogateNeighbourNodes(
    GeometryPointerType p_geometry,
    BrepCurveOnSurfaceType::Pointer p_brep_geometry,
    const std::string& rSlaveLayerName)
{
    std::vector<array_1d<double, 3>> curve_derivatives(3, ZeroVector(3));
    CoordinatesArrayType projected_point = ZeroVector(3);
    CoordinatesArrayType geometry_center = p_geometry->Center();

    std::string projection_layer_name = rSlaveLayerName;
    bool is_projected = ProjectToSkinBoundary(
        mrSlaveSkinModelPart,
        projection_layer_name,
        geometry_center,
        projected_point,
        curve_derivatives,
        10);

    if (!is_projected) {
        projection_layer_name.clear();
        is_projected = ProjectToSkinBoundary(
            mrSlaveSkinModelPart,
            projection_layer_name,
            geometry_center,
            projected_point,
            curve_derivatives,
            10);
    }

    if (!is_projected) {
        return;
    }

    IndexType new_node_id = mrSlaveSkinModelPart->GetRootModelPart().Nodes().size() + 1;
    auto p_slave_skin_node = new Node(new_node_id, projected_point);

    CoordinatesArrayType tangent_vector = curve_derivatives[1];
    const double tangent_norm = norm_2(tangent_vector);
    if (tangent_norm > std::numeric_limits<double>::epsilon()) {
        tangent_vector /= tangent_norm;
    }

    Vector normal_vector = ZeroVector(3);
    normal_vector[0] = tangent_vector[1];
    normal_vector[1] = -tangent_vector[0];

    double curvature = 0.0;
    const CoordinatesArrayType& first_derivative = curve_derivatives[1];
    const CoordinatesArrayType& second_derivative = curve_derivatives[2];
    const double first_norm = norm_2(first_derivative);
    if (first_norm > std::numeric_limits<double>::epsilon()) {
        curvature = norm_2(MathUtils<double>::CrossProduct(first_derivative, second_derivative)) / std::pow(first_norm, 3);
    }

    p_slave_skin_node->SetValue(NORMAL, normal_vector);
    p_slave_skin_node->SetValue(LOCAL_TANGENT, tangent_vector);
    p_slave_skin_node->SetValue(CURVATURE, curvature);

    NodePointerVector neighbour_nodes;
    neighbour_nodes.push_back(p_slave_skin_node);
    p_geometry->SetValue(NEIGHBOUR_NODES, neighbour_nodes);

    (void)p_brep_geometry;
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
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES)/2;
        auto parameter_space_extremes = rModelPart.GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);

        const double epsilon = std::max(knot_step_uv[0], knot_step_uv[1])/1e12;

        // compute the knot span where the point lies
        int knot_span_surrounded_u_id = floor((rPointGlobalCoordinates[0] - parameter_space_extremes[0][0]+epsilon)/knot_step_uv[0]);
        int knot_span_surrounded_v_id = floor((rPointGlobalCoordinates[1] - parameter_space_extremes[1][0]+epsilon)/knot_step_uv[1]);

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
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES)/2;
        auto parameter_space_extremes = rModelPart.GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);

        const double epsilon = std::max(knot_step_uv[0], knot_step_uv[1])/1e12;

        // compute the knot span where the point lies
        int knot_span_surrounded_u_id = floor((rPointGlobalCoordinates[0] - parameter_space_extremes[0][0]+epsilon)/knot_step_uv[0]);
        int knot_span_surrounded_v_id = floor((rPointGlobalCoordinates[1] - parameter_space_extremes[1][0]+epsilon)/knot_step_uv[1]);

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
        Vector knot_step_uv = rModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES)/2;

        const double epsilon = std::max(knot_step_uv[0], knot_step_uv[1])/1e12;

        auto parameter_space_extremes = rModelPart.GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);

        // compute the knot span where the point lies
        int knot_span_surrounded_u_id = floor((rPointGlobalCoordinates[0] - parameter_space_extremes[0][0]+epsilon)/knot_step_uv[0]);
        int knot_span_surrounded_v_id = floor((rPointGlobalCoordinates[1] - parameter_space_extremes[1][0]+epsilon)/knot_step_uv[1]);

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
        auto FindBestProjection = [&](const std::string& rLayerFilter,
                                      double& rBestDistance,
                                      std::vector<array_1d<double, 3>>& rBestDerivatives,
                                      std::string& rBestLayerName,
                                      bool& rFound) {
            for (auto &i_curve : pSkinModelPart->Geometries())
            {   
                if  (rLayerFilter != "" && i_curve.GetValue(IDENTIFIER) != rLayerFilter) continue;

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

                    projected_point_local[0] = t0 + (t1 - t0) * double(i_guess) / (nInitialGuesses - 1);

                    bool is_projected = nurbs_curve_geometry->ProjectionPointGlobalToLocalSpace(rPoint, projected_point_local, 1e-13);

                    if (!is_projected) continue;

                    nurbs_curve_geometry->GlobalCoordinates(projected_point, projected_point_local);

                    double curr_distance = norm_2(rPoint - projected_point);

                    if (curr_distance < rBestDistance) {
                        rBestDistance = curr_distance;
                        rProjectedPoint = projected_point;
                        nurbs_curve_geometry->GlobalSpaceDerivatives(rBestDerivatives, projected_point_local, 2);
                        rFound = true;
                        rBestLayerName = i_curve.GetValue(IDENTIFIER);
                    }
                }
            }
        };

        bool is_projected_at_least_once = false;
        double best_distance = 1e12;
        std::vector<array_1d<double, 3>> best_curve_derivatives(3, ZeroVector(3));
        std::string best_layer_name = "";

        FindBestProjection(rLayerName, best_distance, best_curve_derivatives, best_layer_name, is_projected_at_least_once);

        rCurveDerivatives = best_curve_derivatives;

        if (rLayerName == "") rLayerName = best_layer_name;

        /*
        // DynamicBins fallback (disabled)
        {
            Vector knot_span_sizes;
            if (pSkinModelPart->Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = pSkinModelPart->GetValue(KNOT_SPAN_SIZES);
            } else if (pSkinModelPart->GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = pSkinModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
            } else if (pSkinModelPart->GetRootModelPart().Has(KNOT_SPAN_SIZES)) {
                knot_span_sizes = pSkinModelPart->GetRootModelPart().GetValue(KNOT_SPAN_SIZES);
            }

            double mesh_size = 0.0;
            if (knot_span_sizes.size() >= 2) {
                mesh_size = std::max(knot_span_sizes[0], knot_span_sizes[1]);
            }
            if (mesh_size <= 0.0) {
                mesh_size = 1.0;
            }
            const double radius = std::sqrt(2.0) * mesh_size;

            PointVector points;
            points.reserve(pSkinModelPart->NumberOfConditions());
            for (auto& r_cond : pSkinModelPart->Conditions()) {
                if (rLayerName != "" && r_cond.GetValue(LAYER_NAME) != rLayerName) {
                    continue;
                }
                const auto& c = r_cond.GetGeometry().Center();
                points.push_back(Kratos::make_intrusive<PointType>(r_cond.Id(), c.X(), c.Y(), c.Z()));
            }

            if (!points.empty()) {
                DynamicBins bins(points.begin(), points.end());
                const int number_of_results = static_cast<int>(points.size());
                PointVector results(number_of_results);
                std::vector<double> list_of_distances(number_of_results);

                PointerType p_point_to_search(new PointType(10000, rPoint[0], rPoint[1], rPoint[2]));
                const int obtained_results = bins.SearchInRadius(*p_point_to_search, radius,
                    results.begin(), list_of_distances.begin(), number_of_results);

                if (obtained_results > 0) {
                    double minimum_distance = list_of_distances[0];
                    int nearest_index = 0;
                    for (int i_result = 1; i_result < obtained_results; ++i_result) {
                        if (list_of_distances[i_result] < minimum_distance) {
                            minimum_distance = list_of_distances[i_result];
                            nearest_index = i_result;
                        }
                    }

                    if (results[nearest_index] && (minimum_distance < best_distance || !is_projected_at_least_once)) {
                        best_distance = minimum_distance;
                        rProjectedPoint[0] = (*results[nearest_index])[0];
                        rProjectedPoint[1] = (*results[nearest_index])[1];
                        rProjectedPoint[2] = (*results[nearest_index])[2];
                        rCurveDerivatives = std::vector<array_1d<double, 3>>(3, ZeroVector(3));
                        is_projected_at_least_once = true;

                        if (rLayerName == "") {
                            const IndexType condition_id = results[nearest_index]->Id();
                            if (pSkinModelPart->HasCondition(condition_id)) {
                                rLayerName = pSkinModelPart->GetCondition(condition_id).GetValue(LAYER_NAME);
                            }
                        }
                    }
                }
            }

            if (!is_projected_at_least_once) {
                KRATOS_WARNING("::[IgaContactProcessSbm]:: no projection found on the skin boundary with layer ") << rLayerName 
                            << " for the point: " << rPoint << std::endl;
            }
        }
        */

        // Fallback: check only curve endpoints (t0, t1) without ProjectionPointGlobalToLocalSpace
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

            CoordinatesArrayType projected_point_local = ZeroVector(3);
            CoordinatesArrayType projected_point = ZeroVector(3);

            projected_point_local[0] = t0;
            nurbs_curve_geometry->GlobalCoordinates(projected_point, projected_point_local);
            double curr_distance = norm_2(rPoint - projected_point);
            if (curr_distance < best_distance) {
                best_distance = curr_distance;
                rProjectedPoint = projected_point;
                is_projected_at_least_once = true;
                nurbs_curve_geometry->GlobalSpaceDerivatives(rCurveDerivatives, projected_point_local, 2);
            }

            projected_point_local[0] = t1;
            nurbs_curve_geometry->GlobalCoordinates(projected_point, projected_point_local);
            curr_distance = norm_2(rPoint - projected_point);
            if (curr_distance < best_distance) {
                is_projected_at_least_once = true;
                best_distance = curr_distance;
                rProjectedPoint = projected_point;
                nurbs_curve_geometry->GlobalSpaceDerivatives(rCurveDerivatives, projected_point_local, 2);
            }
        }

        // Keep warning behavior tied only to initial FindBestProjection result
        if (!is_projected_at_least_once) {
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


void IgaContactProcessSbm::PrepareIntegrationOnTrueBoundary(
    ModelPart& rMasterSkinModelPart,
    ModelPart& rSlaveSkinModelPart) const
{
    auto& r_contact_sub_model_part = mrContactModelPart->GetSubModelPart("contact");
    auto& r_master_contact_sub_model_part = mrContactModelPart->GetSubModelPart("master");
    auto& r_slave_contact_sub_model_part = mrContactModelPart->GetSubModelPart("slave");

    for (auto& r_condition : r_contact_sub_model_part.Conditions()) {
        r_condition.SetValue(INTEGRATION_POINTS_MASTER, std::vector<Vector>());
        r_condition.SetValue(INTEGRATION_POINTS_SLAVE, std::vector<Vector>());
        r_condition.SetValue(INTEGRATION_WEIGHTS_MASTER, std::vector<double>());
        r_condition.SetValue(INTEGRATION_WEIGHTS_SLAVE, std::vector<double>());
    }

    PointVector points_master;
    PointVector points_slave;
    points_master.reserve(r_contact_sub_model_part.NumberOfConditions() + r_master_contact_sub_model_part.NumberOfConditions());
    points_slave.reserve(r_contact_sub_model_part.NumberOfConditions() + r_slave_contact_sub_model_part.NumberOfConditions());

    for (auto& r_condition : r_contact_sub_model_part.Conditions()) {
        const auto master_center = r_condition.GetGeometry().GetGeometryPart(0).Center();
        points_master.push_back(PointTypePointer(new PointType(r_condition.Id(), master_center.X(), master_center.Y(), master_center.Z())));

        const auto slave_center = r_condition.GetGeometry().GetGeometryPart(1).Center();
        points_slave.push_back(PointTypePointer(new PointType(r_condition.Id(), slave_center.X(), slave_center.Y(), slave_center.Z())));
    }

    for (auto& r_condition : r_master_contact_sub_model_part.Conditions()) {
        const auto center = r_condition.GetGeometry().Center();
        points_master.push_back(PointTypePointer(new PointType(r_condition.Id(), center.X(), center.Y(), center.Z())));
    }

    for (auto& r_condition : r_slave_contact_sub_model_part.Conditions()) {
        const auto center = r_condition.GetGeometry().Center();
        points_slave.push_back(PointTypePointer(new PointType(r_condition.Id(), center.X(), center.Y(), center.Z())));
    }

    if (points_master.empty() && points_slave.empty()) {
        KRATOS_WARNING("IgaContactProcessSbm") << "No contact conditions found to prepare integration on true boundary." << std::endl;
        return;
    }

    int order = 5;
    if (mParameters.Has("precision_order_on_integration")) {
        order = mParameters["precision_order_on_integration"].GetInt();
    }
    const SizeType num_points = order + 1;

    Vector knot_span_sizes = ZeroVector(3);
    if (mrContactModelPart->Has(KNOT_SPAN_SIZES)) {
        knot_span_sizes = mrContactModelPart->GetValue(KNOT_SPAN_SIZES);
    } else if (mrContactModelPart->GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
        knot_span_sizes = mrContactModelPart->GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    } else if (mrContactModelPart->GetRootModelPart().Has(KNOT_SPAN_SIZES)) {
        knot_span_sizes = mrContactModelPart->GetRootModelPart().GetValue(KNOT_SPAN_SIZES);
    }

    double mesh_size = 0.0;
    if (knot_span_sizes.size() >= 2) {
        mesh_size = std::max(knot_span_sizes[0], knot_span_sizes[1]);
    }
    if (mesh_size <= 0.0) {
        mesh_size = 1.0;
    }
    const double radius = std::sqrt(2.0) * mesh_size;

    const int number_of_results = 100000;
    ModelPart::NodesContainerType::ContainerType results(number_of_results);
    std::vector<double> list_of_distances(number_of_results);

    const auto& integration_point_list_u = IntegrationPointUtilities::s_gauss_legendre[num_points - 1];
    const double distance_tolerance = std::numeric_limits<double>::epsilon();

    auto process_skin = [&](ModelPart& rSkinModelPart,
                            DynamicBins& rBins,
                            const Variable<std::vector<Vector>>& rPointsVariable,
                            const Variable<std::vector<double>>& rWeightsVariable) {

        if (rSkinModelPart.NumberOfConditions() == 0) {
            return;
        }

        for (auto& r_skin_condition : rSkinModelPart.Conditions()) {
            const auto& r_point0 = r_skin_condition.GetGeometry()[0].Coordinates();
            const auto& r_point1 = r_skin_condition.GetGeometry()[1].Coordinates();

            array_1d<double, 3> distance_vector = r_point1 - r_point0;
            const double segment_length = norm_2(distance_vector);

            if (segment_length <= distance_tolerance) {
                continue;
            }

            for (SizeType point_index = 0; point_index < num_points; ++point_index) {
                array_1d<double, 3> current_point = r_point0;
                for (IndexType i = 0; i < 3; ++i) {
                    current_point[i] += distance_vector[i] * integration_point_list_u[point_index][0];
                }
                const double current_weight = integration_point_list_u[point_index][1] * segment_length;

                PointerType p_point_to_search(new PointType(10000, current_point[0], current_point[1], current_point[2]));
                const int obtained_results = rBins.SearchInRadius(*p_point_to_search, radius, results.begin(), list_of_distances.begin(), number_of_results);

                if (obtained_results == 0) {
                    continue;
                }

                double minimum_distance = list_of_distances[0];
                int nearest_index = 0;
                for (int i_result = 1; i_result < obtained_results; ++i_result) {
                    if (list_of_distances[i_result] < minimum_distance) {
                        minimum_distance = list_of_distances[i_result];
                        nearest_index = i_result;
                    }
                }
                                const IndexType condition_id = results[nearest_index]->Id();
                if (!r_contact_sub_model_part.HasCondition(condition_id)) {
                    continue;
                }

                auto& r_target_condition = r_contact_sub_model_part.GetCondition(condition_id);

                auto integration_points = r_target_condition.GetValue(rPointsVariable);
                Vector integration_point_vector = ZeroVector(3);
                for (IndexType i = 0; i < 3; ++i) {
                    integration_point_vector[i] = current_point[i];
                }
                integration_points.push_back(integration_point_vector);
                r_target_condition.SetValue(rPointsVariable, integration_points);

                auto integration_weights = r_target_condition.GetValue(rWeightsVariable);
                integration_weights.push_back(current_weight);
                r_target_condition.SetValue(rWeightsVariable, integration_weights);
            }
        }
    };

    if (!points_master.empty()) {
        DynamicBins test_bins_master(points_master.begin(), points_master.end());
        process_skin(rMasterSkinModelPart, test_bins_master, INTEGRATION_POINTS_MASTER, INTEGRATION_WEIGHTS_MASTER);
    }

    if (!points_slave.empty()) {
        DynamicBins test_bins_slave(points_slave.begin(), points_slave.end());
        process_skin(rSlaveSkinModelPart, test_bins_slave, INTEGRATION_POINTS_SLAVE, INTEGRATION_WEIGHTS_SLAVE);
    }
}


} // End namespace Kratos
