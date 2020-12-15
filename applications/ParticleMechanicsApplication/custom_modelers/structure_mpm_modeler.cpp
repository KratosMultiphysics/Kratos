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


namespace Kratos
{
    ///@name Stages
    ///@{

    void StructureMpmModeler::SetupGeometryModel()
    {
        CheckParameters();

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

        KRATOS_ERROR_IF(background_grid_model_part.NumberOfElements() == 0) << "Background_Grid model part has zero elements!\n";

        // create coupling conditions on interface depending on the dimension
        ModelPart& r_fem_interface = (mIsOriginMpm)
            ? p_model_fem->GetModelPart(destination_interface_sub_model_part_name)
            : p_model_fem->GetModelPart(origin_interface_sub_model_part_name);
        std::vector<GeometryPointerType> interface_geoms;
        const IndexType dim = 2;
        if (dim == 2)
        {
            CreateInterfaceLineCouplingConditions(r_fem_interface, interface_geoms);
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }
        const IndexType gauss_order = mParameters["gauss_integration_order"].GetInt();
        GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod(gauss_order-1);

        std::vector<GeometryPointerType> quads_structure;
        CreateStructureQuadraturePointGeometries<
            std::vector<GeometryPointerType>, std::vector<GeometryPointerType>>(
                interface_geoms, quads_structure, integration_method);

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
        FixMPMDestInterfaceNodes(mpm_coupling_nodes);

        std::vector<GeometryPointerType>& p_quads_origin = (mIsOriginMpm) ? quads_mpm : quads_structure;
        std::vector<GeometryPointerType>& p_quads_dest = (mIsOriginMpm) ? quads_structure : quads_mpm;

        // Determine next condition number
        IndexType condition_id = (coupling_model_part.GetRootModelPart().NumberOfConditions() == 0)
            ? 1 : (coupling_model_part.GetRootModelPart().ConditionsEnd() - 1)->Id() + 1;
        for (IndexType i = 0; i < p_quads_origin.size(); ++i) {
            coupling_model_part.AddCondition(Kratos::make_intrusive<Condition>(
                condition_id + i, Kratos::make_shared<CouplingGeometry<Node<3>>>(p_quads_origin[i], p_quads_dest[i])));
        }
    }

    void StructureMpmModeler::UpdateGeometryModel()
    {
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

        UpdateMpmQuadraturePointGeometries<3,
            typename ModelPart::ConditionsContainerType>(
                coupling_model_part.Conditions(), mpm_background_grid_model_part);

        ModelPart& mpm_coupling_nodes = (p_model_mpm->HasModelPart("coupling_nodes"))
            ? p_model_mpm->GetModelPart("coupling_nodes")
            : p_model_mpm->CreateModelPart("coupling_nodes");

        KRATOS_ERROR_IF(mpm_coupling_nodes.NumberOfNodes() == 0)
             << "The MPM model has no model part 'coupling_nodes', which should have been created in the structure_mpm_modeler setupGeometry";

        // Unfix interface nodes of previous timestep (just mpm, since the fem nodes stay the same),
        ReleaseMPMDestInterfaceNodes(mpm_coupling_nodes);

        // Set all mpm interface nodal forces to be zero
        if (mIsOriginMpm && mParameters["is_gauss_seidel"].GetBool()) {
            for (auto interface_node : mpm_coupling_nodes.NodesArray())
            {
                array_1d<double, 3 >& point_load = (interface_node)->FastGetSolutionStepValue(POINT_LOAD);
                point_load.clear();
            }
        }

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

        // We fix the interface nodes so they can receive the prescribed displacements from FEM.
        FixMPMDestInterfaceNodes(mpm_coupling_nodes);
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
        std::vector<GeometryPointerType>& rGeometries)
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
}
