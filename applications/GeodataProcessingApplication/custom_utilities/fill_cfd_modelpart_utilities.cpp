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

    void FillCfdModelpartUtilities::FillModelPart(ModelPart& NewModelPart) {
        // mrModelPart = ModelPart
        //      TopModelPart
        //      BottomModelPart
        //      LateralSector_... (from 1 to n)
        //      SKIN_ISOSURFACE
        //      Parts_Fluid
        //      NO LateralModelPart
        // NewModelPart = Modelpart for the CFD analysis. It must be empty

        /*
        * TODO: check in NewModelPart. MUST BE EMPTY
        */

       KRATOS_INFO("FillModelPart") << "START" << std::endl;

        auto& r_nodes_array = mrModelPart.Nodes();
        auto& r_elems_array = mrModelPart.Elements();

        Properties::Pointer p_prop = NewModelPart.pGetProperties(0);

        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> elem_ids;

        // Parts_Fluid (Nodes and Elements)
        ModelPart& r_fluid_model_part = NewModelPart.CreateSubModelPart("Parts_Fluid");
        KRATOS_INFO("FillModelPart") << "Parts_Fluid CREATED!" << std::endl;

        // we create new nodes in "Parts_Fluid"
        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            const auto p_node = r_nodes_array.begin() + i_node;
            r_fluid_model_part.CreateNewNode(   p_node->Id(),
                                                p_node->Coordinates()[0],
                                                p_node->Coordinates()[1],
                                                p_node->Coordinates()[2]);
        }
        KRATOS_INFO("FillModelPart") << "Parts_Fluid Nodes CREATED!" << std::endl;

        for (int i_elem = 0; i_elem < static_cast<int>(r_elems_array.size()); ++i_elem) {
            const auto p_elem = r_elems_array.begin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();   // Nodes

            r_fluid_model_part.CreateNewElement("Element3D4N",
                                                p_elem->Id(),
                                                {{r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id()}},
                                                p_prop);
        }
        KRATOS_INFO("FillModelPart") << "Parts_Fluid Elements CREATED!" << std::endl;

        // TopModelPart, BottomModelPart, LateralSector_, SKIN_ISOSURFACE (Nodes and Conditions)
        for (auto& submodelpart : mrModelPart.SubModelParts()) {

            if (submodelpart.NumberOfElements() > 0)     // only Parts_Fluid has elements. we skip it
                continue;
            
            std::string submodelpart_name = submodelpart.Name();    // name of current sub model part
            if (submodelpart_name == "LateralModelPart") {
                KRATOS_INFO("FillModelPart") << "LateralModelPart SKIPPED!" << std::endl;
                continue;
            }

            ModelPart& r_current_sub_model_part = NewModelPart.CreateSubModelPart(submodelpart_name);   // new sub model part with the name of submodelpart_name
            KRATOS_INFO("FillModelPart") << submodelpart_name << " CREATED!" << std::endl;

            // Nodes
            std::vector<ModelPart::IndexType> node_ids;

            for (int i_node = 0; i_node < static_cast<int>(submodelpart.NumberOfNodes()); ++i_node) {
                const auto p_node = submodelpart.NodesBegin() + i_node;
                node_ids.push_back(p_node->Id());
            }

            r_current_sub_model_part.AddNodes(node_ids);    // Nodes are added in current sub model part
            KRATOS_INFO("FillModelPart") << submodelpart_name << " Nodes ADDED!" << std::endl;

            // Conditions
            for (int i_cond = 0; i_cond < static_cast<int>(submodelpart.NumberOfConditions()); ++i_cond) {
                const auto p_cond = submodelpart.ConditionsBegin() + i_cond;
                auto& r_geom = p_cond->GetGeometry();   // Nodes

                r_current_sub_model_part.CreateNewCondition("SurfaceCondition3D3N",
                                                            p_cond->Id(),
                                                            {{r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id()}},
                                                            p_prop);
            }
            KRATOS_INFO("FillModelPart") << submodelpart_name << " Conditions CREATED!" << std::endl;
        }

        KRATOS_INFO("FillModelPart") << "END" << std::endl;

    }

/*
    void FillCfdModelpartUtilities::FillPartsFluid(const std::string& ElementModelName) {
        // mrModelPart = ModelPart

        KRATOS_INFO("FillPartsFluid") << "ElementModelName: " << ElementModelName << std::endl;

        if (!mrModelPart.HasSubModelPart("Parts_Fluid")) {
            KRATOS_INFO("FillPartsFluid") << "No Parts_Fluid Model Part is present!" << std::endl;
            return;
        }

        auto& r_elem_model_part = mrModelPart.GetSubModelPart(ElementModelName);
        auto& r_fluid_model_part = mrModelPart.GetSubModelPart("Parts_Fluid");

        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> elem_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_elem_model_part.NumberOfNodes()); ++i_node) {
            auto p_node = r_elem_model_part.NodesBegin() + i_node;
            node_ids.push_back(p_node->Id());
        }

        for (int i_elem = 0; i_elem < static_cast<int>(r_elem_model_part.NumberOfElements()); ++i_elem) {
            auto p_elem = r_elem_model_part.ElementsBegin() + i_elem;
            elem_ids.push_back(p_elem->Id());
        }

        // nodes and elements are added in "Parts_Fluid" sub model part
        r_fluid_model_part.AddNodes(node_ids);
        r_fluid_model_part.AddElements(elem_ids);

        // delete sub model part
        mrModelPart.RemoveSubModelPart(ElementModelName);

        // for Debugging purpose only
        KRATOS_INFO("\t[Parts_Fluid] In total ") << node_ids.size() << " nodes are moved from " << ElementModelName << " to \"Parts_Fluid\" sub model part!" << std::endl;
        KRATOS_INFO("\t[Parts_Fluid] In total ") << elem_ids.size() << " conditions are moved from " << ElementModelName << " to \"Parts_Fluid\" sub model part!" << std::endl;
    }
*/
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


    void FillCfdModelpartUtilities::FillSlip(const std::string& ConditionModelName) {
        // mrModelPart = ModelPart

        KRATOS_INFO("FillSlip") << "ConditionModelName: " << ConditionModelName << std::endl;

        if (!mrModelPart.HasSubModelPart("Slip")) {
            KRATOS_INFO("FillSlip") << "No Slip Model Part is present!" << std::endl;
            return;
        }

        auto& r_cond_model_part = mrModelPart.GetSubModelPart(ConditionModelName);
        auto& r_slip_model_part = mrModelPart.GetSubModelPart("Slip");

        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> cond_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_cond_model_part.NumberOfNodes()); ++i_node) {
            auto p_node = r_cond_model_part.NodesBegin() + i_node;
            node_ids.push_back(p_node->Id());
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_cond_model_part.NumberOfConditions()); ++i_cond) {
            auto p_cond = r_cond_model_part.ConditionsBegin() + i_cond;
            cond_ids.push_back(p_cond->Id());
        }

        // nodes and conditions are added in "Slip" sub model part
        r_slip_model_part.AddNodes(node_ids);
        r_slip_model_part.AddConditions(cond_ids);

        // delete sub model part
        mrModelPart.RemoveSubModelPart(ConditionModelName);

        // for Debugging purpose only
        KRATOS_INFO("\t[Slip] In total ") << node_ids.size() << " nodes are moved from " << ConditionModelName << " to \"Slip\" sub model part!" << std::endl;
        KRATOS_INFO("\t[Slip] In total ") << cond_ids.size() << " conditions are moved from " << ConditionModelName << " to \"Slip\" sub model part!" << std::endl;
    }


    void FillCfdModelpartUtilities::FillNoslip(const std::string& ConditionModelName) {
        // mrModelPart = ModelPart

        KRATOS_INFO("FillNoslip") << "ConditionModelName: " << ConditionModelName << std::endl;

        if (!mrModelPart.HasSubModelPart("NoSlip")) {
            KRATOS_INFO("FillNoslip") << "No NoSlip Model Part is present!" << std::endl;
            return;
        }

        auto& r_cond_model_part = mrModelPart.GetSubModelPart(ConditionModelName);
        auto& r_noslip_model_part = mrModelPart.GetSubModelPart("NoSlip");

        std::vector<ModelPart::IndexType> node_ids;
        std::vector<ModelPart::IndexType> cond_ids;

        for (int i_node = 0; i_node < static_cast<int>(r_cond_model_part.NumberOfNodes()); ++i_node) {
            auto p_node = r_cond_model_part.NodesBegin() + i_node;
            node_ids.push_back(p_node->Id());
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_cond_model_part.NumberOfConditions()); ++i_cond) {
            auto p_cond = r_cond_model_part.ConditionsBegin() + i_cond;
            cond_ids.push_back(p_cond->Id());
        }

        // nodes and conditions are added in "NoSlip" sub model part
        r_noslip_model_part.AddNodes(node_ids);
        r_noslip_model_part.AddConditions(cond_ids);

        // delete sub model part
        mrModelPart.RemoveSubModelPart(ConditionModelName);

        // for Debugging purpose only
        KRATOS_INFO("\t[NoSlip] In total ") << node_ids.size() << " nodes are moved from " << ConditionModelName << " to \"NoSlip\" sub model part!" << std::endl;
        KRATOS_INFO("\t[NoSlip] In total ") << cond_ids.size() << " conditions are moved from " << ConditionModelName << " to \"NoSlip\" sub model part!" << std::endl;
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
