//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                               Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/extend_pressure_condition_process.h"
#include "processes/find_elements_neighbours_process.h"

namespace Kratos {

template <SizeType TDim>
ExtendPressureConditionProcess<TDim>::ExtendPressureConditionProcess(
    ModelPart& r_model_part)
    : mrModelPart(r_model_part)
{
}

template <>
void ExtendPressureConditionProcess<2>::Execute()
{
    // We search the neighbours for the 3 node generation of line loads
    auto find_neigh = FindElementalNeighboursProcess(mrModelPart, 2, 5);
    find_neigh.Execute();

    // Remove previous line loads
    this->RemovePreviousLineLoads();
    this->CreateNewConditions();
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::RemovePreviousLineLoads()
{
    // We remove only the line loads of all the SubModels
    std::vector<std::string> submodel_parts_names = mrModelPart.GetSubModelPartNames();
    std::vector<std::string> pressure_sub_models;
    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        if (submodel_parts_names[i].substr(0, 11) == "Normal_Load") {
            // Remove the line loads
            auto& r_sub_model = mrModelPart.GetSubModelPart(submodel_parts_names[i]);
            for (auto it_cond = r_sub_model.ConditionsBegin(); it_cond != r_sub_model.ConditionsEnd(); it_cond++) {
                it_cond->Set(TO_ERASE, true);
            }
        }
    }
    mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<2>::CreateNewConditions()
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    int maximum_condition_id;
    this->GetMaximumConditionIdOnSubmodelPart(maximum_condition_id);

    // Loop over the elements (all active, the inactive have been removed in GeneratingDEM)
    for (auto it_elem = mrModelPart.Elements().ptr_begin();  it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {

        // We count how many nodes are wet
        auto& r_geometry = (*it_elem)->GetGeometry();
        int wet_nodes_counter = 0, non_wet_local_id_node = 10, pressure_id;

        for (IndexType local_id = 0; local_id < r_geometry.PointsNumber(); ++local_id) {
            if (r_geometry[local_id].GetValue(PRESSURE_ID) != 0) {
                wet_nodes_counter++;
                pressure_id = r_geometry[local_id].GetValue(PRESSURE_ID);
            } else {
                non_wet_local_id_node = local_id;
            }
        }

        if (wet_nodes_counter == 2) {
            this->GenerateLineLoads2Nodes(non_wet_local_id_node, pressure_id, maximum_condition_id, it_elem);
            r_process_info[ITER] = 1;
        } else if (wet_nodes_counter == 3) {
            this->GenerateLineLoads3Nodes(pressure_id, maximum_condition_id, it_elem);
            r_process_info[ITER] = 1;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ExtendPressureConditionProcess<2>::GenerateLineLoads2Nodes(
    const int NonWetLocalIdNode,
    const int PressureId,
    int& rMaximumConditionId,
    ModelPart::ElementsContainerType::ptr_iterator itElem
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
    ModelPart::PropertiesType::Pointer p_properties = r_sub_model_part.pGetProperties(1);
    auto& r_geom = (*itElem)->GetGeometry();

    const IndexType id_1 = NonWetLocalIdNode == 0 ? 0 : NonWetLocalIdNode == 1 ? 1 : 2;
    const IndexType id_2 = NonWetLocalIdNode == 0 ? 1 : NonWetLocalIdNode == 1 ? 2 : 0;
    const IndexType id_3 = NonWetLocalIdNode == 0 ? 2 : NonWetLocalIdNode == 1 ? 0 : 1;

    std::vector<IndexType> condition_nodes_id(2);
    condition_nodes_id[0] = r_geom[id_3].Id();
    condition_nodes_id[1] = r_geom[id_2].Id();
	rMaximumConditionId++;

    // We create the Line Load Condition
    const auto& r_line_condition = r_sub_model_part.CreateNewCondition(
					                    "LineLoadCondition2D2N",
					                    rMaximumConditionId,
					                    condition_nodes_id,
					                    p_properties, 0);

    // Adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_condition);

    // Adding the nodes to the SubModelPart
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_3].Id()));
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_2].Id()));
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ExtendPressureConditionProcess<2>::GenerateLineLoads3Nodes(
    const int PressureId,
    int& rMaximumConditionId,
    ModelPart::ElementsContainerType::ptr_iterator itElem
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
    ModelPart::PropertiesType::Pointer p_properties = r_sub_model_part.pGetProperties(1);

    // We get the neighbour elements
    WeakPointerVector<Element>& r_elem_neigb = (*itElem)->GetValue(NEIGHBOUR_ELEMENTS);

    IndexType alone_edge_local_id = 10;
    for (IndexType i = 0; i < r_elem_neigb.size(); ++i) {
        if ((*itElem)->Id() == r_elem_neigb[i].Id()) {
            alone_edge_local_id = i;
        }
    }
    KRATOS_ERROR_IF(alone_edge_local_id == 10) << "Unexpected error in extrapolating the pressure load..." << std::endl;
    
    const IndexType id_1 = alone_edge_local_id == 0 ? 0 : alone_edge_local_id == 1 ? 1 : 2;
    const IndexType id_2 = alone_edge_local_id == 0 ? 1 : alone_edge_local_id == 1 ? 2 : 0;
    const IndexType id_3 = alone_edge_local_id == 0 ? 2 : alone_edge_local_id == 1 ? 0 : 1;

    std::vector<IndexType> condition_nodes_id(2);
    auto& r_geom = (*itElem)->GetGeometry();
    condition_nodes_id[0] = r_geom[id_2].Id();
    condition_nodes_id[1] = r_geom[id_1].Id();
    rMaximumConditionId++;
    const auto& r_line_cond1 = r_sub_model_part.CreateNewCondition(
                                    "LineLoadCondition2D2N",
                                    rMaximumConditionId,
                                    condition_nodes_id,
                                    p_properties, 0);    

    condition_nodes_id[0] = r_geom[id_1].Id();
    condition_nodes_id[1] = r_geom[id_3].Id();
    rMaximumConditionId++;
    const auto& r_line_cond2 = r_sub_model_part.CreateNewCondition(
                                    "LineLoadCondition2D2N",
                                    rMaximumConditionId,
                                    condition_nodes_id,
                                    p_properties, 0);

    // Adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond1);
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond2);

    // Adding the nodes to the SubModelPart
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_3].Id()));
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_2].Id()));
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_1].Id()));
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<3>::Execute()
{
    // todo implementation in 3D
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<2>::GetMaximumConditionIdOnSubmodelPart(
      int& rMaximumConditionId
)
{
    rMaximumConditionId = 0;
    for (auto it_cond = mrModelPart.ConditionsBegin();
         it_cond != mrModelPart.ConditionsEnd(); it_cond++) {
        if (((*it_cond)).Id() > rMaximumConditionId) rMaximumConditionId = ((*it_cond)).Id();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::CalculateNumberOfElementsOnNodes()
{
    // Reset the Flag
    for (auto it_node = mrModelPart.Nodes().ptr_begin(); it_node != mrModelPart.Nodes().ptr_end(); ++it_node) {
        int& number_of_elems = (*it_node)->GetValue(NUMBER_OF_ACTIVE_ELEMENTS);
        number_of_elems = 0;
    }
    // Add the active elements
    for (auto itElem = mrModelPart.Elements().ptr_begin(); itElem != mrModelPart.Elements().ptr_end(); ++itElem) {
        bool condition_is_active = true;
        if ((*itElem)->IsDefined(ACTIVE)) {
            condition_is_active = (*itElem)->Is(ACTIVE);
        }
        if (condition_is_active) {
            auto& r_geom = (*itElem)->GetGeometry();
            for (IndexType i = 0; i <  r_geom.size(); ++i) {
                int& number_of_elems = r_geom[i].GetValue(NUMBER_OF_ACTIVE_ELEMENTS);
                number_of_elems++;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class ExtendPressureConditionProcess<2>;
template class ExtendPressureConditionProcess<3>;

}  // namespace Kratos
