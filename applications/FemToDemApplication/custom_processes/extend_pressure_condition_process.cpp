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
void ExtendPressureConditionProcess<2>::CreateAndAddPressureConditions2Nodes(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    const unsigned int LocalId,
    const int PressureId,
	int& rMaximumConditionId,
    std::vector<IndexType>& ToEraseConditionsId
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);

    std::vector<IndexType> condition_nodes_id(2);
    ModelPart::PropertiesType::Pointer p_properties = r_sub_model_part.pGetProperties(1);

    // We add the node to the submodel
    auto& r_geom = (*itElem)->GetGeometry();
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[LocalId].Id()));

    // Set the vectors
    // mNodeIdContainer.push_back(r_geom[LocalId].Id());
    // mNodePressureIdContainer.push_back(PressureId);
    // test
    mrModelPart.pGetNode(r_geom[LocalId].Id())->SetValue(PRESSURE_ID, PressureId);

    KRATOS_WATCH(r_geom[LocalId].Id())

    const IndexType id_1 = LocalId == 0 ? 0 : LocalId == 1 ? 1 : 2;
    const IndexType id_2 = LocalId == 0 ? 1 : LocalId == 1 ? 2 : 0;
    const IndexType id_3 = LocalId == 0 ? 2 : LocalId == 1 ? 0 : 1;

    condition_nodes_id[0] = r_geom[id_1].Id();
    condition_nodes_id[1] = r_geom[id_2].Id();
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

    // adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond1);
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond2);

    // We remove the condition regarding the erased edge...
    for (auto it = mrModelPart.Conditions().ptr_begin(); it != mrModelPart.Conditions().ptr_end(); ++it) {
        // Nodes of the condition
        if ((*it)->GetGeometry().size() > 1) {
            const IndexType Id1 = (*it)->GetGeometry()[0].Id();
            const IndexType Id2 = (*it)->GetGeometry()[1].Id();
            if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_3].Id()) ||
                (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_3].Id())) {
                ToEraseConditionsId.push_back((*it)->Id());
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ExtendPressureConditionProcess<2>::CreateAndAddPressureConditions3Nodes(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    const int PressureId,
	int& rMaximumConditionId,
    std::vector<IndexType>& ToEraseConditionsId
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);

    std::vector<IndexType> condition_nodes_id(2);
    ModelPart::PropertiesType::Pointer p_properties = r_sub_model_part.pGetProperties(1);
    auto& r_geom = (*itElem)->GetGeometry();

    IndexType local_id;
    int aux_counter = 0;
    std::vector<IndexType> inactive_nodes_id;
    std::vector<int> inactive_nodes_local_id;

    const auto& r_process_info = mrModelPart.GetProcessInfo();
	const int counter_of_affected_nodes = r_process_info[ITER];
    if (counter_of_affected_nodes != 1) this->CalculateNumberOfElementsOnNodes();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        if (r_geom[i].GetValue(NUMBER_OF_ACTIVE_ELEMENTS) == 1) {
            local_id = i;
            aux_counter++;
        } else if (r_geom[i].GetValue(NUMBER_OF_ACTIVE_ELEMENTS) == 0) {
            inactive_nodes_id.push_back(r_geom[i].Id());
            inactive_nodes_local_id.push_back(i);
        }
    }

    if (aux_counter == 1) { // common case
        const IndexType id_1 = local_id == 0 ? 0 : local_id == 1 ? 1 : 2;
        const IndexType id_2 = local_id == 0 ? 1 : local_id == 1 ? 2 : 0;
        const IndexType id_3 = local_id == 0 ? 2 : local_id == 1 ? 0 : 1;

        condition_nodes_id[0] = r_geom[id_2].Id();
        condition_nodes_id[1] = r_geom[id_3].Id();
        rMaximumConditionId++;
        const auto& r_line_cond = r_sub_model_part.CreateNewCondition(
                                        "LineLoadCondition2D2N",
                                        rMaximumConditionId,
                                        condition_nodes_id,
                                        p_properties, 0);
        // adding the conditions to the computing model part
        mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond);

        // We remove the condition regarding the erased edges...
        for (auto it = mrModelPart.Conditions().ptr_begin(); it != mrModelPart.Conditions().ptr_end(); ++it) {
            auto& r_it_geometry = (*it)->GetGeometry();
            if (r_it_geometry.size() > 1) { // avoid nodal forces
                const IndexType Id1 = r_it_geometry[0].Id();
                const IndexType Id2 = r_it_geometry[1].Id();

                // Remove the old conditions
                if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_1].Id()) ||
                    (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_1].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                } else if ((Id1 == r_geom[id_1].Id() && Id2 == r_geom[id_3].Id()) ||
                           (Id2 == r_geom[id_1].Id() && Id1 == r_geom[id_3].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                }
            }
        }
    // One inactive node
    } else if (inactive_nodes_id.size() == 1) {
        const IndexType id_1 = inactive_nodes_local_id[0] == 0 ? 0 : inactive_nodes_local_id[0] == 1 ? 1 : 2;
        const IndexType id_2 = inactive_nodes_local_id[0] == 0 ? 1 : inactive_nodes_local_id[0] == 1 ? 2 : 0;
        const IndexType id_3 = inactive_nodes_local_id[0] == 0 ? 2 : inactive_nodes_local_id[0] == 1 ? 0 : 1;

        condition_nodes_id[0] = r_geom[id_2].Id();
        condition_nodes_id[1] = r_geom[id_3].Id();
        rMaximumConditionId++;
        const auto& r_line_cond = r_sub_model_part.CreateNewCondition(
                                           "LineLoadCondition2D2N",
                                           rMaximumConditionId,
                                           condition_nodes_id,
                                           p_properties, 0);

        // adding the conditions to the computing model part
        mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond);

        // We remove the condition regarding the erased edges...
        for (auto it = mrModelPart.Conditions().ptr_begin(); it != mrModelPart.Conditions().ptr_end(); ++it) {
            auto& r_it_geom = (*it)->GetGeometry();
            if (r_it_geom.size() > 1) {
                const IndexType Id1 = r_it_geom[0].Id();
                const IndexType Id2 = r_it_geom[1].Id();

                if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_1].Id()) ||
                    (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_1].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                } else if ((Id1 == r_geom[id_1].Id() && Id2 == r_geom[id_3].Id()) ||
                        (Id2 == r_geom[id_1].Id() && Id1 == r_geom[id_3].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                }
            }
        }
    } else if (inactive_nodes_id.size() == 3) { // elem and nodes are removed afterwards
        // We remove the condition regarding the erased edges...
        for (auto it = mrModelPart.Conditions().ptr_begin(); it != mrModelPart.Conditions().ptr_end(); ++it) {
            auto& r_it_geom = (*it)->GetGeometry();
            if (r_it_geom.size() > 1) {
                const IndexType Id1 = r_it_geom[0].Id();
                const IndexType Id2 = r_it_geom[1].Id();

                const IndexType id_1 = inactive_nodes_local_id[0] == 0 ? 0 : inactive_nodes_local_id[0] == 1 ? 1 : 2;
                const IndexType id_2 = inactive_nodes_local_id[0] == 0 ? 1 : inactive_nodes_local_id[0] == 1 ? 2 : 0;
                const IndexType id_3 = inactive_nodes_local_id[0] == 0 ? 2 : inactive_nodes_local_id[0] == 1 ? 0 : 1;

                if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_1].Id()) ||
                    (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_1].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                } else if ((Id1 == r_geom[id_1].Id() && Id2 == r_geom[id_3].Id()) ||
                           (Id2 == r_geom[id_1].Id() && Id1 == r_geom[id_3].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                } else if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_3].Id()) ||
                           (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_3].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                }
            }
        }
    } else if (aux_counter == 0 && inactive_nodes_id.size() == 0) {

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

        condition_nodes_id[0] = r_geom[id_1].Id();
        condition_nodes_id[1] = r_geom[id_2].Id();
        rMaximumConditionId++;
        const auto& r_line_cond1 = r_sub_model_part.CreateNewCondition(
                                        "LineLoadCondition2D2N",
                                        rMaximumConditionId,
                                        condition_nodes_id,
                                        p_properties, 0);
        // adding the conditions to the computing model part
        mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond1);

        condition_nodes_id[0] = r_geom[id_3].Id();
        condition_nodes_id[1] = r_geom[id_1].Id();
        rMaximumConditionId++;
        const auto& r_line_cond2 = r_sub_model_part.CreateNewCondition(
                                        "LineLoadCondition2D2N",
                                        rMaximumConditionId,
                                        condition_nodes_id,
                                        p_properties, 0);
        // adding the conditions to the computing model part
        mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond2);

        // We remove the condition regarding the erased edge...
        for (auto it = mrModelPart.Conditions().ptr_begin(); it != mrModelPart.Conditions().ptr_end(); ++it) {
            auto& r_it_geometry = (*it)->GetGeometry();
            if (r_it_geometry.size() > 1) { // avoid nodal forces
                const IndexType Id1 = r_it_geometry[0].Id();
                const IndexType Id2 = r_it_geometry[1].Id();

                // Remove the old conditions
                if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_3].Id()) ||
                    (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_3].Id())) {
                    ToEraseConditionsId.push_back((*it)->Id());
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ExtendPressureConditionProcess<2>::CreateAndAddPressureConditions1Node(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    const int PressureId,
	int& rMaximumConditionId,
    std::vector<IndexType>& ToEraseConditionsId
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);

    std::vector<IndexType> condition_nodes_id(2);
    auto p_properties = r_sub_model_part.pGetProperties(1);

    IndexType local_id;
    auto& r_geom = (*itElem)->GetGeometry();
    // Let's identify the node with the pressure load (wet)
	for (IndexType i = 0; i < r_geom.PointsNumber(); ++i) {
        if (r_geom[i].GetValue(PRESSURE_ID) != 0) {
            local_id = i;
            break;
        }
    }

    const IndexType id_1 = local_id == 0 ? 0 : local_id == 1 ? 1 : 2;
    const IndexType id_2 = local_id == 0 ? 1 : local_id == 1 ? 2 : 0;
    const IndexType id_3 = local_id == 0 ? 2 : local_id == 1 ? 0 : 1;

    // We add the nodes to the submodel
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_2].Id()));
    r_sub_model_part.AddNode(mrModelPart.pGetNode(r_geom[id_3].Id()));

    // Set the 2 new nodes to wet nodes
    // mNodeIdContainer.push_back(r_geom[id_2].Id());
    // mNodeIdContainer.push_back(r_geom[id_3].Id());
    // mNodePressureIdContainer.push_back(PressureId);
    // mNodePressureIdContainer.push_back(PressureId);
    mrModelPart.pGetNode(r_geom[id_2].Id())->SetValue(PRESSURE_ID, PressureId);
    mrModelPart.pGetNode(r_geom[id_3].Id())->SetValue(PRESSURE_ID, PressureId);

    // We create the new pressure conditions
    condition_nodes_id[0] = r_geom[id_1].Id();
    condition_nodes_id[1] = r_geom[id_2].Id();
	rMaximumConditionId++;
    const auto& r_line_cond1 = r_sub_model_part.CreateNewCondition(
					                    "LineLoadCondition2D2N",
					                    rMaximumConditionId,
					                    condition_nodes_id,
					                    p_properties, 0);
    condition_nodes_id[0] = r_geom[id_2].Id();
    condition_nodes_id[1] = r_geom[id_3].Id();
	rMaximumConditionId++;
    const auto& r_line_cond2 = r_sub_model_part.CreateNewCondition(
					                    "LineLoadCondition2D2N",
					                    rMaximumConditionId,
					                    condition_nodes_id,
					                    p_properties, 0);
    condition_nodes_id[0] = r_geom[id_3].Id();
    condition_nodes_id[1] = r_geom[id_1].Id();
	rMaximumConditionId++;
    const auto& r_line_cond3 = r_sub_model_part.CreateNewCondition(
					                    "LineLoadCondition2D2N",
					                    rMaximumConditionId,
					                    condition_nodes_id,
					                    p_properties, 0);

    // adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond1);
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond2);
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(r_line_cond3);
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

template <>
void ExtendPressureConditionProcess<2>::Execute()
{
    // We search the neighbours for some purposes
    auto find_neigh = FindElementalNeighboursProcess(mrModelPart, 2, 5);
    find_neigh.Execute();

	int maximum_condition_id, counter_of_affected_nodes = 0;
    std::vector<IndexType> ToEraseConditionsId;
    this->GetMaximumConditionIdOnSubmodelPart(maximum_condition_id);

    // Loop over the elements in order to extrapolate the pressure load on its nodes if necessary
    for (auto it_elem = mrModelPart.Elements().ptr_begin();  it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
        bool condition_is_active = true;
        if ((*it_elem)->IsDefined(ACTIVE)) {
            condition_is_active = (*it_elem)->Is(ACTIVE);
        }
        // it_elem's going to be removed
        if (condition_is_active == false && (*it_elem)->GetValue(SMOOTHING) == false) {
            KRATOS_WATCH((*it_elem)->Id())
            unsigned int local_id, counter = 0, pressure_id;
            auto& r_geometry = (*it_elem)->GetGeometry();
            // Loop over nodes in order to check if there's pressure on nodes
            for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                if (r_geometry.GetPoint(i).GetValue(PRESSURE_ID) != 0) {
                    pressure_id = r_geometry.GetPoint(i).GetValue(PRESSURE_ID);
                    counter++;
                } else {
                    local_id = i;
                }
            }
            if (counter == 2) {
                std::cout << "22222222222222" << std::endl;
                this->CreateAndAddPressureConditions2Nodes(it_elem, local_id, pressure_id, maximum_condition_id, ToEraseConditionsId);
                counter_of_affected_nodes++;
                // We use this flag to enter once on each element
                (*it_elem)->SetValue(SMOOTHING, true);
            } else if (counter == 1) {
                std::cout << "11111111111111111111" << std::endl;
                this->CreateAndAddPressureConditions1Node(it_elem, pressure_id, maximum_condition_id, ToEraseConditionsId);
                counter_of_affected_nodes++;
                // We use this flag to enter once on each element
                (*it_elem)->SetValue(SMOOTHING, true);
            } else if (counter == 3) {
                this->CreateAndAddPressureConditions3Nodes(it_elem, pressure_id, maximum_condition_id, ToEraseConditionsId);
                counter_of_affected_nodes++;
                // We use this flag to enter once on each element
                (*it_elem)->SetValue(SMOOTHING, true);
            }
        }
    }
    auto& r_process_info = mrModelPart.GetProcessInfo();
    r_process_info[ITER] = counter_of_affected_nodes;

    // for (IndexType i = 0; i < mNodeIdContainer.size(); ++i) {
    //     mrModelPart.GetNode(mNodeIdContainer[i]).SetValue(PRESSURE_ID, mNodePressureIdContainer[i]);
    // }
    // mNodeIdContainer.clear();
    // mNodePressureIdContainer.clear();

    for (IndexType i = 0; i < ToEraseConditionsId.size(); ++i) {
        mrModelPart.RemoveConditionFromAllLevels(ToEraseConditionsId[i]);
    }
    ToEraseConditionsId.clear();
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
template <SizeType TDim>
bool ExtendPressureConditionProcess<TDim>::CheckIfHasConditionId(const IndexType Id)
{
    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); it_cond++) {
        if ((*it_cond).Id() == Id) {
            return true;
        }
    }
    return false;
}
/***********************************************************************************/
/***********************************************************************************/

template class ExtendPressureConditionProcess<2>;
template class ExtendPressureConditionProcess<3>;

}  // namespace Kratos
