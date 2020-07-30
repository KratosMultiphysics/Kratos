//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "fill_cfd_modelpart_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void FillCfdModelpartUtilities::FillPartsFluid(ModelPart& OriginModelPart, const std::string& ElementModelName) {
        // mrModelPart = CfdModelPart

        if (!mrModelPart.HasSubModelPart("Parts_Fluid"))
            mrModelPart.CreateSubModelPart("Parts_Fluid");

        auto& r_cfd_nodes_array = mrModelPart.Nodes();
        auto& r_cfd_elems_array = mrModelPart.Elements();

        auto& r_elem_model_part = OriginModelPart.GetSubModelPart(ElementModelName);
        auto& r_nodes_array = r_elem_model_part.Nodes();
        auto& r_elements_array = r_elem_model_part.Elements();

        Properties::Pointer p_prop = OriginModelPart.pGetProperties(0);
        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> elem_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            const auto node = r_nodes_array.begin() + i_node;
            auto existing_node_it = r_cfd_nodes_array.find(node->Id());
            if (existing_node_it == mrModelPart.NodesEnd()) {
                // node does not exist
                mrModelPart.CreateNewNode(node->Id(),
                                          node->Coordinates()[0],
                                          node->Coordinates()[1],
                                          node->Coordinates()[2]);
            }
            // [CHECK IT] in any case we add the node into the sub model part
            node_ids.push_back(node->Id());
        }

        for (int i_elem = 0; i_elem < static_cast<int>(r_elements_array.size()); ++i_elem) {
            const auto elem = r_elements_array.begin() + i_elem;
            auto existing_element_it = r_cfd_elems_array.find(elem->Id());
            if (existing_element_it == mrModelPart.ElementsEnd()) {
                // element does not exist
                auto& r_geom = elem->GetGeometry();  // nodes
                mrModelPart.CreateNewElement("Element3D4N",
                                             elem->Id(),
                                             {r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id()},
                                             p_prop);
            }
            // [CHECK IT] in any case we add the element into the sub model part
            elem_ids.push_back(elem->Id());
        }

        auto& r_fluid_model_part = mrModelPart.GetSubModelPart("Parts_Fluid");
        r_fluid_model_part.AddNodes(node_ids);
        r_fluid_model_part.AddElements(elem_ids);

        // for Debugging purpose only
        KRATOS_INFO("*** [PARTS_FLUID] In total ") << node_ids.size() << " nodes are added!" << std::endl;
        KRATOS_INFO("*** [PARTS_FLUID] In total ") << elem_ids.size() << " elements are added!" << std::endl;
    }


    void FillCfdModelpartUtilities::FillInlet(ModelPart& OriginModelPart, const std::string& ConditionModelName) {
        // mrModelPart = CfdModelPart

        if (!mrModelPart.HasSubModelPart("Inlet"))
            mrModelPart.CreateSubModelPart("Inlet");

        auto& r_cfd_nodes_array = mrModelPart.Nodes();
        auto& r_cfd_conds_array = mrModelPart.Conditions();

        auto& r_cond_model_part = OriginModelPart.GetSubModelPart(ConditionModelName);
        auto& r_nodes_array = r_cond_model_part.Nodes();
        auto& r_conditions_array = r_cond_model_part.Conditions();

        Properties::Pointer p_prop = OriginModelPart.pGetProperties(0);
        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> cond_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            const auto node = r_nodes_array.begin() + i_node;
            auto existing_node_it = r_cfd_nodes_array.find(node->Id());
            if (existing_node_it == mrModelPart.NodesEnd()) {
				// node does not exist
                auto n = mrModelPart.CreateNewNode( node->Id(),
                                                    node->Coordinates()[0],
                                                    node->Coordinates()[1],
                                                    node->Coordinates()[2]);
            }
			// [CHECK IT] in any case we add the node into the sub model part
			node_ids.push_back(node->Id());
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_conditions_array.size()); ++i_cond) {
            const auto cond = r_conditions_array.begin() + i_cond;
            auto existing_condition_it = r_cfd_conds_array.find(cond->Id());
            if (existing_condition_it == mrModelPart.ConditionsEnd()) {
				// condition does not exist
                auto& r_geom = cond->GetGeometry();  // nodes
                auto c = mrModelPart.CreateNewCondition("SurfaceCondition3D3N",
                                                        cond->Id(),
                                                        {r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id()},
                                                        p_prop);
            }
            // [CHECK IT] in any case we add the condition into the sub model part
            cond_ids.push_back(cond->Id());
        }

        auto& r_fluid_model_part = mrModelPart.GetSubModelPart("Inlet");
        r_fluid_model_part.AddNodes(node_ids);
        r_fluid_model_part.AddConditions(cond_ids);

        // for Debugging purpose only
        KRATOS_INFO("*** [INLET] In total ") << node_ids.size() << " nodes are added!" << std::endl;
        KRATOS_INFO("*** [INLET] In total ") << cond_ids.size() << " conditions are added!" << std::endl;
    }


    void FillCfdModelpartUtilities::FillOutlet(ModelPart& OriginModelPart, const std::string& ConditionModelName) {
        // mrModelPart = CfdModelPart

        if (!mrModelPart.HasSubModelPart("Outlet"))
            mrModelPart.CreateSubModelPart("Outlet");

        auto& r_cfd_nodes_array = mrModelPart.Nodes();
        auto& r_cfd_conds_array = mrModelPart.Conditions();

        auto& r_cond_model_part = OriginModelPart.GetSubModelPart(ConditionModelName);
        auto& r_nodes_array = r_cond_model_part.Nodes();
        auto& r_conditions_array = r_cond_model_part.Conditions();

        Properties::Pointer p_prop = OriginModelPart.pGetProperties(0);
        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> cond_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            const auto node = r_nodes_array.begin() + i_node;
            auto existing_node_it = r_cfd_nodes_array.find(node->Id());
            if (existing_node_it == mrModelPart.NodesEnd()) {
				// node does not exist
                auto n = mrModelPart.CreateNewNode( node->Id(),
                                                    node->Coordinates()[0],
                                                    node->Coordinates()[1],
                                                    node->Coordinates()[2]);
            }
			// [CHECK IT] in any case we add the node into the sub model part
			node_ids.push_back(node->Id());
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_conditions_array.size()); ++i_cond) {
            const auto cond = r_conditions_array.begin() + i_cond;
            auto existing_condition_it = r_cfd_conds_array.find(cond->Id());
            if (existing_condition_it == mrModelPart.ConditionsEnd()) {
				// condition does not exist
                auto& r_geom = cond->GetGeometry();  // nodes
                auto c = mrModelPart.CreateNewCondition("SurfaceCondition3D3N",
                                                        cond->Id(),
                                                        {r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id()},
                                                        p_prop);
            }
            // [CHECK IT] in any case we add the condition into the sub model part
            cond_ids.push_back(cond->Id());
        }

        auto& r_fluid_model_part = mrModelPart.GetSubModelPart("Outlet");
        r_fluid_model_part.AddNodes(node_ids);
        r_fluid_model_part.AddConditions(cond_ids);

        // for Debugging purpose only
        KRATOS_INFO("*** [OUTLET] In total ") << node_ids.size() << " nodes are added!" << std::endl;
        KRATOS_INFO("*** [OUTLET] In total ") << cond_ids.size() << " conditions are added!" << std::endl;
    }


    void FillCfdModelpartUtilities::FillSlip(ModelPart& OriginModelPart, const std::string& ConditionModelName) {
        // mrModelPart = CfdModelPart

        if (!mrModelPart.HasSubModelPart("Slip"))
            mrModelPart.CreateSubModelPart("Slip");

        auto& r_cfd_nodes_array = mrModelPart.Nodes();
        auto& r_cfd_conds_array = mrModelPart.Conditions();

        auto& r_cond_model_part = OriginModelPart.GetSubModelPart(ConditionModelName);
        auto& r_nodes_array = r_cond_model_part.Nodes();
        auto& r_conditions_array = r_cond_model_part.Conditions();

        Properties::Pointer p_prop = OriginModelPart.pGetProperties(0);
        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> cond_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            const auto node = r_nodes_array.begin() + i_node;
            auto existing_node_it = r_cfd_nodes_array.find(node->Id());
            if (existing_node_it == mrModelPart.NodesEnd()) {
				// node does not exist
                auto n = mrModelPart.CreateNewNode( node->Id(),
                                                    node->Coordinates()[0],
                                                    node->Coordinates()[1],
                                                    node->Coordinates()[2]);
            }
			// [CHECK IT] in any case we add the node into the sub model part
			node_ids.push_back(node->Id());
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_conditions_array.size()); ++i_cond) {
            const auto cond = r_conditions_array.begin() + i_cond;
            auto existing_condition_it = r_cfd_conds_array.find(cond->Id());
            if (existing_condition_it == mrModelPart.ConditionsEnd()) {
				// condition does not exist
                auto& r_geom = cond->GetGeometry();  // nodes
                auto c = mrModelPart.CreateNewCondition("SurfaceCondition3D3N",
                                                        cond->Id(),
                                                        {r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id()},
                                                        p_prop);
            }
            // [CHECK IT] in any case we add the condition into the sub model part
            cond_ids.push_back(cond->Id());
        }

        auto& r_fluid_model_part = mrModelPart.GetSubModelPart("Slip");
        r_fluid_model_part.AddNodes(node_ids);
        r_fluid_model_part.AddConditions(cond_ids);

        // for Debugging purpose only
        KRATOS_INFO("*** [OUTLET] In total ") << node_ids.size() << " nodes are added!" << std::endl;
        KRATOS_INFO("*** [OUTLET] In total ") << cond_ids.size() << " conditions are added!" << std::endl;
    }


    void FillCfdModelpartUtilities::FillNoslip(ModelPart& OriginModelPart, const std::string& ConditionModelName) {
        // mrModelPart = CfdModelPart

        if (!mrModelPart.HasSubModelPart("NoSlip"))
            mrModelPart.CreateSubModelPart("NoSlip");

        auto& r_cfd_nodes_array = mrModelPart.Nodes();
        auto& r_cfd_conds_array = mrModelPart.Conditions();

        auto& r_cond_model_part = OriginModelPart.GetSubModelPart(ConditionModelName);
        auto& r_nodes_array = r_cond_model_part.Nodes();
        auto& r_conditions_array = r_cond_model_part.Conditions();

        Properties::Pointer p_prop = OriginModelPart.pGetProperties(0);
        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> cond_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            const auto node = r_nodes_array.begin() + i_node;
            auto existing_node_it = r_cfd_nodes_array.find(node->Id());
            if (existing_node_it == mrModelPart.NodesEnd()) {
				// node does not exist
                auto n = mrModelPart.CreateNewNode( node->Id(),
                                                    node->Coordinates()[0],
                                                    node->Coordinates()[1],
                                                    node->Coordinates()[2]);
            }
			// [CHECK IT] in any case we add the node into the sub model part
			node_ids.push_back(node->Id());
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_conditions_array.size()); ++i_cond) {
            const auto cond = r_conditions_array.begin() + i_cond;
            auto existing_condition_it = r_cfd_conds_array.find(cond->Id());
            if (existing_condition_it == mrModelPart.ConditionsEnd()) {
				// condition does not exist
                auto& r_geom = cond->GetGeometry();  // nodes
                auto c = mrModelPart.CreateNewCondition("SurfaceCondition3D3N",
                                                        cond->Id(),
                                                        {r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id()},
                                                        p_prop);
            }
            // [CHECK IT] in any case we add the condition into the sub model part
            cond_ids.push_back(cond->Id());
        }

        auto& r_fluid_model_part = mrModelPart.GetSubModelPart("NoSlip");
        r_fluid_model_part.AddNodes(node_ids);
        r_fluid_model_part.AddConditions(cond_ids);

        // for Debugging purpose only
        KRATOS_INFO("*** [OUTLET] In total ") << node_ids.size() << " nodes are added!" << std::endl;
        KRATOS_INFO("*** [OUTLET] In total ") << cond_ids.size() << " conditions are added!" << std::endl;
    }



    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const FillCfdModelpartUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
