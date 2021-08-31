//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "structure_mpm_modeler.h"
#include "geometries/bounding_box.h"
#include "boost/geometry/geometry.hpp"
#include "integration/integration_point_utilities.h"


namespace Kratos
{
    ///@name Stages
    ///@{


    void StructureMpmModeler::SetupGeometryModel()
    {
        CheckParameters();
        
        const bool is_create_segmented_fem_quads = false;

        Model* p_model_mpm = (mIsOriginMpm) ? mpModelOrigin : mpModelDest;
        Model* p_model_fem = (mIsOriginMpm) ? mpModelDest : mpModelOrigin;

        ModelPart& coupling_model_part = (mpModelOrigin->HasModelPart("coupling"))
            ? mpModelOrigin->GetModelPart("coupling")
            : mpModelOrigin->CreateModelPart("coupling");

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;
        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }
        else
        {
            KRATOS_ERROR << "Automatic interface creation is not implemented yet" << std::endl;
        }

        ModelPart& background_grid_model_part = (p_model_mpm->HasModelPart("Background_Grid"))
            ? p_model_mpm->GetModelPart("Background_Grid")
            : p_model_mpm->CreateModelPart("Background_Grid");

        background_grid_model_part.GetProcessInfo().SetValue(IS_COSIM_COUPLED, true);

        KRATOS_ERROR_IF(background_grid_model_part.NumberOfElements() == 0) << "Background_Grid model part has zero elements!\n";

        // create coupling conditions on interface depending on the dimension
        ModelPart& r_fem_interface = (mIsOriginMpm)
            ? p_model_fem->GetModelPart(destination_interface_sub_model_part_name)
            : p_model_fem->GetModelPart(origin_interface_sub_model_part_name);
        std::vector<GeometryPointerType> interface_geoms;
        std::vector<GeometryPointerType> segmented_fem_quad_points;
        std::vector<GeometryPointerType> quads_structure;
        
        const IndexType gauss_order = mParameters["gauss_integration_order"].GetInt();
        GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod(gauss_order - 1);

        const IndexType dim = 2;
        if (dim == 2)
        {
            CreateInterfaceLineCouplingConditions(r_fem_interface, interface_geoms,
                background_grid_model_part);
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }


        if (is_create_segmented_fem_quads)
        {
            CreateExactlySegmentedStructureQuadraturePointGeometries(r_fem_interface,
                interface_geoms, background_grid_model_part, quads_structure);
        }
        else
        {
            CreateStructureQuadraturePointGeometries<
                std::vector<GeometryPointerType>, std::vector<GeometryPointerType>>(
                    interface_geoms, quads_structure, integration_method);
        }


        std::vector<GeometryPointerType> quads_mpm(quads_structure.size());
        CreateMpmQuadraturePointGeometries<2, std::vector<GeometryPointerType>>(
            quads_structure, quads_mpm, background_grid_model_part);

        // Transfer everything into the coupling modelpart
        ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");

        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");

        // Temporarily store coupling nodes in modelparts. Then we can setNodes instead of addNodes to avoid duplicates!
        ModelPart& mpm_coupling_nodes = (p_model_mpm->HasModelPart("coupling_nodes"))
            ? p_model_mpm->GetModelPart("coupling_nodes")
            : p_model_mpm->CreateModelPart("coupling_nodes");

        const double tolerance = mParameters["minimum_shape_function_value"].GetDouble();
        for (IndexType i = 0; i < quads_mpm.size(); ++i) {
            for (IndexType j = 0; j < quads_mpm[i]->size(); ++j) {
                if (quads_mpm[i]->ShapeFunctionValue(0,j) > tolerance)
                {
                    mpm_coupling_nodes.AddNode(quads_mpm[i]->pGetPoint(j));
                }
            }
        }
        ModelPart& fem_coupling_nodes = (p_model_fem->HasModelPart("coupling_nodes"))
            ? p_model_fem->GetModelPart("coupling_nodes")
            : p_model_fem->CreateModelPart("coupling_nodes");
        for (IndexType i = 0; i < quads_structure.size(); ++i) {
            for (IndexType j = 0; j < quads_structure[i]->size(); ++j) {
                fem_coupling_nodes.AddNode(quads_structure[i]->pGetPoint(j));
            }
        }

        if (mIsOriginMpm)
        {
            coupling_interface_origin.SetNodes(mpm_coupling_nodes.pNodes());
            coupling_interface_destination.SetNodes(fem_coupling_nodes.pNodes());
        }
        else
        {
            coupling_interface_origin.SetNodes(fem_coupling_nodes.pNodes());
            coupling_interface_destination.SetNodes(mpm_coupling_nodes.pNodes());
        }

        // We fix the interface nodes so they can receive the prescribed displacements from FEM.
        // not needed for penalty?!?!
        // FixMPMDestInterfaceNodes(mpm_coupling_nodes);

        std::vector<GeometryPointerType>& p_quads_origin = (mIsOriginMpm) ? quads_mpm : quads_structure;
        std::vector<GeometryPointerType>& p_quads_dest = (mIsOriginMpm) ? quads_structure : quads_mpm;

        // Determine next condition number
        IndexType condition_id = (coupling_model_part.GetRootModelPart().NumberOfConditions() == 0)
            ? 1 : (coupling_model_part.GetRootModelPart().ConditionsEnd() - 1)->Id() + 1;
        for (IndexType i = 0; i < p_quads_origin.size(); ++i) {
            coupling_model_part.AddCondition(Kratos::make_intrusive<Condition>(
                condition_id + i, Kratos::make_shared<CouplingGeometry<Node<3>>>(p_quads_origin[i], p_quads_dest[i])));
        }
        KRATOS_WATCH("2222222222222HEeeeeeeeeeeeeeeeeeeeeeeeeeeeeeLLO")
        
        ModelPart& material_model_model_part = (p_model_mpm->HasModelPart("MPM_Material"))
            ? p_model_mpm->GetModelPart("MPM_Material")
            : p_model_mpm->CreateModelPart("MPM_Material");

        auto mpm_material_model_part_name = "Coupling_Interface";
        ModelPart& mpm_support_mp = material_model_model_part.HasSubModelPart(mpm_material_model_part_name)
        ? material_model_model_part.GetSubModelPart(mpm_material_model_part_name)
        : material_model_model_part.CreateSubModelPart(mpm_material_model_part_name);
        // TO DO Check id properties
        Properties::Pointer p_properties = coupling_model_part.GetRootModelPart().CreateNewProperties(10000);
        const Condition& new_condition = KratosComponents<Condition>::Get("MPMParticlePenaltyCouplingInterfaceCondition2D3N");
        KRATOS_WATCH("HEeeeeeeeeeeeeeeeeeeeeeeeeeeeeeLLO")
        for (IndexType i = 0; i < p_quads_dest.size(); ++i) {

            Condition::Pointer p_condition = new_condition.Create(condition_id + i+p_quads_origin.size() , p_quads_dest[i], p_properties);
            ProcessInfo process_info = ProcessInfo();
            std::vector<array_1d<double, 3>> mpc_xg = { ZeroVector(3) };
            std::vector<array_1d<double, 3>> mpc_imposed_displacement = { ZeroVector(3) };

            p_quads_dest[i]->GlobalCoordinates(mpc_xg[0],0);
            p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, { mpc_imposed_displacement }, process_info);

            std::vector<double> mpc_penalty_factor(1);
            
            mpc_penalty_factor[0]=1e5;
            p_condition->SetValuesOnIntegrationPoints(PENALTY_FACTOR, mpc_penalty_factor , process_info);
            p_condition->SetValuesOnIntegrationPoints(MPC_COORD, mpc_xg , process_info);
            p_condition->Set(INTERFACE);
            // p_condition->Initialize(process_info);
            // p_condition->InitializeSolutionStep(process_info);
            
            // Setting particle condition's initial condition
            // p_condition->SetValuesOnIntegrationPoints(MPC_COORD, mpc_xg , process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_AREA,  mpc_area  , process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_NORMAL, { mpc_normal }, process_info);

            
            // p_condition->SetValuesOnIntegrationPoints(MPC_DISPLACEMENT, { mpc_displacement }, process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, { mpc_imposed_displacement }, process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_VELOCITY, { mpc_velocity }, process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_VELOCITY, { mpc_imposed_velocity }, process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_ACCELERATION, { mpc_acceleration }, process_info);
            // p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_ACCELERATION, { mpc_imposed_acceleration }, process_info);
                                // Mark as boundary condition
            p_condition->Set(BOUNDARY, false);
            
            mpm_support_mp.AddCondition(p_condition);
            
        }
        
        // KRATOS_WATCH(*p_model_mpm)
    }

    void StructureMpmModeler::UpdateGeometryModel()
    {
        KRATOS_WATCH("NEXT STEP")
        // The FEM coupling quad points naturally deform with the FEM domain and track the interface correctly.
        // The MPM coupling quad points need to be updated with the FEM quad point locations.

        Model* p_model_mpm = (mIsOriginMpm) ? mpModelOrigin : mpModelDest;
        const IndexType mpm_index = (mIsOriginMpm) ? 0 : 1;

        KRATOS_ERROR_IF_NOT(mpModelOrigin->HasModelPart("coupling"))
            << "The origin model does not have the 'coupling' model part. This should have been created during StructureMpmModeler::SetupGeometryModel\n";
        ModelPart& coupling_model_part = mpModelOrigin->GetModelPart("coupling");

        KRATOS_ERROR_IF(coupling_model_part.NumberOfConditions() == 0)
            << "The MPM model model part 'coupling_model_part' has no conditions, which should have been created in the StructureMpmModeler::SetupGeometryModel\n";

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;
        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }
        else
        {
            KRATOS_ERROR << "Automatic interface creation is not implemented yet" << std::endl;
        }

        ModelPart& mpm_background_grid_model_part = p_model_mpm->GetModelPart("Background_Grid");

        //UpdateMpmQuadraturePointGeometries<3,
        //    typename ModelPart::ConditionsContainerType>(
        //        coupling_model_part.Conditions(), mpm_background_grid_model_part);
        UpdateMpmQuadraturePointGeometriesWithFilter<3,
            typename ModelPart::ConditionsContainerType>(
                coupling_model_part.Conditions(), mpm_background_grid_model_part);

        ModelPart& mpm_coupling_nodes = (p_model_mpm->HasModelPart("coupling_nodes"))
            ? p_model_mpm->GetModelPart("coupling_nodes")
            : p_model_mpm->CreateModelPart("coupling_nodes");

        KRATOS_ERROR_IF(mpm_coupling_nodes.NumberOfNodes() == 0)
             << "The MPM model has no model part 'coupling_nodes', which should have been created in the structure_mpm_modeler setupGeometry";

        // Unfix interface nodes of previous timestep (just mpm, since the fem nodes stay the same),
        // ReleaseMPMDestInterfaceNodes(mpm_coupling_nodes);
        

        // // Set all mpm interface nodal forces to be zero
        // if (mIsOriginMpm && mParameters["is_gauss_seidel"].GetBool()) {
        //     for (auto interface_node : mpm_coupling_nodes.NodesArray())
        //     {
        //         array_1d<double, 3 >& point_load = (interface_node)->FastGetSolutionStepValue(POINT_LOAD);
        //         point_load.clear();
        //     }
        // }

        // Remove all old interface nodes from coupling nodes modelpart
        mpm_coupling_nodes.NodesArray().clear();

        ModelPart& coupling_interface_mpm = (mIsOriginMpm)
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.GetSubModelPart("interface_destination");
        coupling_interface_mpm.NodesArray().clear();

        // Add in new interface nodes
        const double tolerance = mParameters["minimum_shape_function_value"].GetDouble();
        for (auto cond_it : coupling_model_part.ConditionsArray())
        {
            auto quads_mpm = cond_it->GetGeometry().pGetGeometryPart(mpm_index);
            for (IndexType j = 0; j < quads_mpm->PointsNumber(); ++j) {
                if (quads_mpm->ShapeFunctionValue(0, j) > tolerance)
                {
                    mpm_coupling_nodes.AddNode(quads_mpm->pGetPoint(j));
                }
            }
        }
        // Set coupling interface nodes in coupling model part
        coupling_interface_mpm.SetNodes(mpm_coupling_nodes.pNodes());

        // // We fix the interface nodes so they can receive the prescribed displacements from FEM.
        // FixMPMDestInterfaceNodes(mpm_coupling_nodes);
    }

    void StructureMpmModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in StructureMpmModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in StructureMpmModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("is_interface_sub_model_parts_specified"))
            << "Missing \"is_interface_sub_model_parts_specified\" in StructureMpmModeler Parameters." << std::endl;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            KRATOS_ERROR_IF_NOT(mParameters.Has("origin_interface_sub_model_part_name"))
                << "Missing \"origin_interface_sub_model_part_name\" in StructureMpmModeler Parameters." << std::endl;

            KRATOS_ERROR_IF_NOT(mParameters.Has("destination_interface_sub_model_part_name"))
                << "Missing \"destination_interface_sub_model_part_name\" in StructureMpmModeler Parameters." << std::endl;
        }

        mParameters.AddMissingParameters(this->GetModelerDefaultSettings());

        const int gauss_order = mParameters["gauss_integration_order"].GetInt();
        KRATOS_ERROR_IF(gauss_order < 1 || gauss_order > 5)
            << "Invalid gauss_integration_order specified. It must be within 1 and 5 inclusive." << std::endl;

    }

	void StructureMpmModeler::CreateInterfaceLineCouplingConditions(
        ModelPart& rInterfaceModelPart,
        std::vector<GeometryPointerType>& rGeometries,
        ModelPart& rBackgroundGrid)
	{
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_geometries_model_part = rInterfaceModelPart.GetSubModelPart("coupling_conditions");
        const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();

        IndexType interface_node_id;
        IndexType trial_interface_node_id;
        IndexType trial_geom_node_id;

        for (size_t node_index = 0; node_index < rInterfaceModelPart.NumberOfNodes() - 1; ++node_index)
        {
            interface_node_id = (rInterfaceModelPart.NodesBegin() + node_index)->Id();
            std::vector< GeometryPointerType> p_geom_vec;
            for (auto& ele_it : root_mp.Elements())
            {
                auto p_geom = ele_it.pGetGeometry();
                for (size_t i = 0; i < p_geom->size(); i++)
                {
                    if ((*p_geom)[i].Id() == interface_node_id)
                    {
                        p_geom_vec.push_back(p_geom);
                    }
                }
            }

            if (p_geom_vec.size() == 0) KRATOS_ERROR << "Interface node not found in modelpart geom\n";

            // Loop over all geometries that have nodes on the interface
            for (size_t interface_geom_index = 0; interface_geom_index < p_geom_vec.size(); interface_geom_index++)
            {
                GeometryType& r_interface_geom = *(p_geom_vec[interface_geom_index]);

                // Loop over remaining interface nodes, see if any of them are nodes in the interface geom
                for (size_t geom_node_index = 0; geom_node_index < r_interface_geom.size(); geom_node_index++)
                {
                    trial_geom_node_id = r_interface_geom[geom_node_index].Id();

                    for (size_t trial_index = node_index + 1; trial_index < rInterfaceModelPart.NumberOfNodes(); ++trial_index)
                    {
                        trial_interface_node_id = (rInterfaceModelPart.NodesBegin() + trial_index)->Id();
                        if (trial_geom_node_id == trial_interface_node_id)
                        {
                            // Another interface node was found in the same geom, make line condition between them
                            Line2D2<NodeType>::Pointer p_line = Kratos::make_shared<Line2D2<NodeType>>(
                                rInterfaceModelPart.pGetNode(interface_node_id), rInterfaceModelPart.pGetNode(trial_interface_node_id));
                            coupling_geometries_model_part.AddGeometry(p_line);
                            rGeometries.push_back(p_line);
                        }
                    }
                }
            }
        }


	}

    void StructureMpmModeler::CreateExactlySegmentedStructureQuadraturePointGeometries(
        ModelPart& rInterfaceModelPart,
        std::vector<GeometryPointerType>& rFEMInterfaceGeometries,
        ModelPart& rBackgroundGrid,
        std::vector<GeometryPointerType>& rSegmentedQuadraturePoints)
    {

        double interface_length_physical = 0;
        double interface_length_boost_intersections = 0;
        double interface_length_quad_points = 0;

        // Get bounding box of FEM interface and make reduced set of grid boost geoms
        BoundingBox<Node<3>> fem_interface_bounding_box = BoundingBox<Node<3>>(rInterfaceModelPart.NodesBegin(), rInterfaceModelPart.NodesEnd());
        array_1d<double, 3>& max_coord = fem_interface_bounding_box.GetMaxPoint().Coordinates();
        array_1d<double, 3>& min_coord = fem_interface_bounding_box.GetMinPoint().Coordinates();
        double safety_buffer = 1.0 * norm_2(max_coord - min_coord);
        for (size_t i = 0; i < 3; i++) {
            max_coord[i] += safety_buffer;
            min_coord[i] -= safety_buffer;
        }

        typedef boost::geometry::model::d2::point_xy<double> BoostPoint;
        typedef boost::geometry::model::linestring<BoostPoint> BoostLinestring;
        typedef boost::geometry::model::multi_linestring<BoostLinestring> BoostMultiLinestring;
        typedef boost::geometry::model::polygon<BoostPoint> BoostPolygon;

        // Make boost linestring of total fem interface
        BoostLinestring fem_interface_line;
        for (size_t i = 0; i < rFEMInterfaceGeometries.size(); ++i)
        {
            auto& r_line = *(rFEMInterfaceGeometries[i].get());
            for (size_t j = 0; j < r_line.PointsNumber(); j++) fem_interface_line.push_back(BoostPoint(r_line.GetPoint(j).X(), r_line.GetPoint(j).Y()));
        }
        interface_length_physical = boost::geometry::length(fem_interface_line);

        // Make boost polygon representation of background grid around FEM interface
        std::vector< BoostPolygon> boost_interface_grid_geoms;
        block_for_each(rBackgroundGrid.Elements(), [&](ModelPart::ElementType& r_element)
            {
                auto& r_geom = r_element.GetGeometry();
                for (auto& r_node : r_geom.Points())
                {
                    array_1d<double, 3>& node_pos = r_node.Coordinates();
                    if (node_pos[0] >= min_coord[0] && node_pos[1] >= min_coord[1]) {
                        if (node_pos[0] <= max_coord[0] && node_pos[1] <= max_coord[1]) {
                            // make boost polygon
                            BoostPolygon boost_poly;
                            std::vector< BoostPoint> boost_polygon_points(r_geom.PointsNumber() + 1);
                            for (size_t i = 0; i < r_geom.PointsNumber(); ++i) boost_polygon_points[i] = BoostPoint(r_geom.GetPoint(i).X(), r_geom.GetPoint(i).Y());
                            boost_polygon_points[r_geom.PointsNumber()] = boost_polygon_points[0];
                            boost_poly.outer().assign(boost_polygon_points.begin(), boost_polygon_points.end());
                            if (boost::geometry::covered_by(fem_interface_line, boost_poly) || boost::geometry::crosses(fem_interface_line, boost_poly))
                            {
                                #pragma omp critical
                                boost_interface_grid_geoms.push_back(boost_poly);

                                break;
                            }

                        }
                    }
                }
            }
        );


        // Now we loop over all fem interface lines and split them over the mpm grid
        IndexPartition<>(rFEMInterfaceGeometries.size()).for_each([&](SizeType i)
            {
                BoostLinestring boost_line;
                std::vector< BoostLinestring> boost_individual_intersected_interface_lines;
                auto& r_line = *(rFEMInterfaceGeometries[i].get());

                for (size_t j = 0; j < r_line.PointsNumber(); j++) boost_line.push_back(BoostPoint(r_line.GetPoint(j).X(), r_line.GetPoint(j).Y()));

                for (size_t boost_grid_index = 0; boost_grid_index < boost_interface_grid_geoms.size(); ++boost_grid_index)
                {
                    bool add_geom = false;
                    std::vector<BoostPoint> intersection_points;
                    BoostLinestring intersection;
                    boost::geometry::intersection(boost_line, boost_interface_grid_geoms[boost_grid_index], intersection_points);
                    if (intersection_points.size() > 0)
                    {
                        if (intersection_points.size() == 1)
                        {
                            // A line end might be within the polygon
                            for (size_t k = 0; k < 2; k++) {
                                if (boost::geometry::covered_by(boost_line[k], boost_interface_grid_geoms[boost_grid_index]))
                                {
                                    add_geom = true;
                                    intersection.push_back(intersection_points[0]);
                                    intersection.push_back(boost_line[k]);
                                    break;
                                }
                            }
                        }
                        else if (intersection_points.size() == 2)
                        {
                            // A line at least spans the polygon
                            add_geom = true;
                            intersection = BoostLinestring(intersection_points.begin(), intersection_points.end());
                        }
                    }
                    else
                    {
                        // A line might be wholly inside the polygon
                        if (boost::geometry::covered_by(boost_line, boost_interface_grid_geoms[boost_grid_index]))
                        {
                            add_geom = true;
                            intersection = boost_line;
                        }
                    }


                    if (add_geom)
                    {
                        bool is_duplicate = false;
                        for (auto& check_line : boost_individual_intersected_interface_lines) {
                            if (boost::geometry::equals(intersection, check_line)) {
                                is_duplicate = true;
                                break;
                            }
                        }
                        if (!is_duplicate)
                        {
                            boost_individual_intersected_interface_lines.push_back(intersection);
                            #pragma omp critical
                            interface_length_boost_intersections += boost::geometry::length(intersection);
                        }
                    }

                }
                KRATOS_ERROR_IF(boost_individual_intersected_interface_lines.size() == 0) << "No intersected grid geoms found!\n";

                // Create quadrature points
                CoordinatesArrayType local_parameter_1 = ZeroVector(3);
                CoordinatesArrayType local_parameter_2 = ZeroVector(3);

                for (size_t j = 0; j < boost_individual_intersected_interface_lines.size(); ++j)
                {
                    BoostLinestring& r_intersected_line = boost_individual_intersected_interface_lines[j];
                    CoordinatesArrayType intersection_low = ZeroVector(3);
                    intersection_low[0] += r_intersected_line[0].get<0>();
                    intersection_low[1] += r_intersected_line[0].get<1>();

                    CoordinatesArrayType intersection_high = ZeroVector(3);
                    intersection_high[0] += r_intersected_line[1].get<0>();
                    intersection_high[1] += r_intersected_line[1].get<1>();
                    r_line.PointLocalCoordinates(local_parameter_1, intersection_low);
                    r_line.PointLocalCoordinates(local_parameter_2, intersection_high);

                    const SizeType IntegrationPointsPerSpan = 2; // TODO this should depend on the basis order
                    IntegrationPointsArrayType integration_points(IntegrationPointsPerSpan);
                    typename IntegrationPointsArrayType::iterator integration_point_iterator = integration_points.begin();

                    IntegrationPointUtilities::IntegrationPoints1D(
                        integration_point_iterator, IntegrationPointsPerSpan,
                        local_parameter_1[0], local_parameter_2[0]);

                    GeometriesArrayType quadrature_point_geometries_fem(IntegrationPointsPerSpan);

                    // CreateQuadraturePointsUtility::Create re-written below to give geometry pointer type
                    auto default_method = r_line.GetDefaultIntegrationMethod();
                    const size_t NumberOfShapeFunctionDerivatives = 1;
                    Vector N;
                    Matrix DN_De;
                    for (IndexType k = 0; k < integration_points.size(); ++k)
                    {
                        r_line.ShapeFunctionsValues(N, integration_points[k]);

                        Matrix N_matrix = ZeroMatrix(1, N.size());
                        if (NumberOfShapeFunctionDerivatives >= 0) {
                            for (IndexType m = 0; m < N.size(); ++m) N_matrix(0, m) = N[m];
                        }

                        /// Get Shape Function Derivatives DN_De, ...
                        if (NumberOfShapeFunctionDerivatives > 0) {
                            r_line.ShapeFunctionsLocalGradients(DN_De, integration_points[k]);
                        }

                        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                            default_method, integration_points[k],
                            N_matrix, DN_De);

                        GeometryPointerType quad_geom = CreateQuadraturePointsUtility<Node<3>>::CreateQuadraturePoint(
                            r_line.WorkingSpaceDimension(), r_line.LocalSpaceDimension(),
                            data_container, r_line);

                        #pragma omp critical
                        {
                            rSegmentedQuadraturePoints.push_back(quad_geom);
                            interface_length_quad_points += quad_geom->IntegrationPoints()[0].Weight() * quad_geom->DeterminantOfJacobian(0);
                        }
                    }
                }
            }
        );


        KRATOS_CHECK_NEAR(interface_length_physical, interface_length_boost_intersections, 1e-9);
        KRATOS_CHECK_NEAR(interface_length_physical, interface_length_quad_points, 1e-9);
    }

}
