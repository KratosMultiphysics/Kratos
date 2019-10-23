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
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];
    mPressureName = (dimension == 2) ? "Normal_Load" : "Pressure_Load";
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::GetMaximumConditionIdOnSubmodelPart(
    int& rMaximumConditionId
)
{
    rMaximumConditionId = 0;
    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); it_cond++) {
        rMaximumConditionId = (((*it_cond)).Id() > rMaximumConditionId) ? ((*it_cond)).Id() : rMaximumConditionId;
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<3>::CreatePressureLoads(
    const int Id1,
    const int Id2,
    const int Id3,
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    ModelPart& rSubModelPart,
    ModelPart::PropertiesType::Pointer pProperties,
    int& rMaximumConditionId
    )
{
    auto& r_geom = (*itElem)->GetGeometry();
    std::vector<IndexType> condition_nodes_id(3);
    condition_nodes_id[0] = r_geom[Id1].Id();
    condition_nodes_id[1] = r_geom[Id2].Id();
    condition_nodes_id[2] = r_geom[Id3].Id();
    rMaximumConditionId++;

    // Adding the nodes to the SubModelPart
    rSubModelPart.AddNodes(condition_nodes_id);

    // We create the Line Load Condition
    const auto p_pressure_condition = rSubModelPart.CreateNewCondition(
                                        "SurfaceLoadCondition3D3N",
                                        rMaximumConditionId,
                                        condition_nodes_id,
                                        pProperties, 0);

    // Adding the conditions to the computing model part
    mrModelPart.GetSubModelPart("computing_domain").AddCondition(p_pressure_condition); 
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<2>::CreateLineLoads(
    const int Id1,
    const int Id2,
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    ModelPart& rSubModelPart,
    ModelPart::PropertiesType::Pointer pProperties,
    int& rMaximumConditionId
    )
{
    std::vector<IndexType> condition_nodes_id(2);
    auto& r_geom = (*itElem)->GetGeometry();
    condition_nodes_id[0] = r_geom[Id1].Id();
    condition_nodes_id[1] = r_geom[Id2].Id();
    rMaximumConditionId++;
    rSubModelPart.AddNodes(condition_nodes_id);

    const auto p_line_cond = rSubModelPart.CreateNewCondition(
                                    "LineLoadCondition2D2N",
                                    rMaximumConditionId,
                                    condition_nodes_id,
                                    pProperties, 0);

    mrModelPart.GetSubModelPart("computing_domain").AddCondition(p_line_cond);
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
    sub_model_name = mPressureName + "-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
    ModelPart::PropertiesType::Pointer p_properties = mpPropertiesVector[PressureId - 1];
    auto& r_geom = (*itElem)->GetGeometry();

    // We check some things...
    auto& r_elem_neigb = (*itElem)->GetValue(NEIGHBOUR_ELEMENTS);
    if (r_elem_neigb[NonWetLocalIdNode].Id() == (*itElem)->Id()) {
        const IndexType id_1 = (NonWetLocalIdNode == 0) ? 0 : (NonWetLocalIdNode == 1) ? 1 : 2;
        const IndexType id_2 = (NonWetLocalIdNode == 0) ? 1 : (NonWetLocalIdNode == 1) ? 2 : 0;
        const IndexType id_3 = (NonWetLocalIdNode == 0) ? 2 : (NonWetLocalIdNode == 1) ? 0 : 1;
        this->CreateLineLoads(id_2, id_3, itElem, r_sub_model_part, p_properties, rMaximumConditionId);
    }
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
    sub_model_name = mPressureName + "-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
    ModelPart::PropertiesType::Pointer p_properties = mpPropertiesVector[PressureId - 1];

    // We get the neighbour elements
    GlobalPointersVector<Element>& r_elem_neigb = (*itElem)->GetValue(NEIGHBOUR_ELEMENTS);

    IndexType alone_edge_local_id = 10;
    int number_of_free_edges = 0, non_free_edge;
    for (IndexType i = 0; i < r_elem_neigb.size(); ++i) {
        if ((*itElem)->Id() == r_elem_neigb[i].Id()) {
            alone_edge_local_id = i;
            number_of_free_edges++;
        } else {
            non_free_edge = i;
        }
    }
    if (number_of_free_edges == 2) {
        const IndexType id_1 = (non_free_edge == 0) ? 0 : (non_free_edge == 1) ? 1 : 2;
        const IndexType id_2 = (non_free_edge == 0) ? 1 : (non_free_edge == 1) ? 2 : 0;
        const IndexType id_3 = (non_free_edge == 0) ? 2 : (non_free_edge == 1) ? 0 : 1;
        this->CreateLineLoads(id_1, id_2, itElem, r_sub_model_part, p_properties, rMaximumConditionId);
        this->CreateLineLoads(id_3, id_1, itElem, r_sub_model_part, p_properties, rMaximumConditionId);
    } else if (number_of_free_edges == 1) {
        const IndexType id_1 = (alone_edge_local_id == 0) ? 0 : (alone_edge_local_id == 1) ? 1 : 2;
        const IndexType id_2 = (alone_edge_local_id == 0) ? 1 : (alone_edge_local_id == 1) ? 2 : 0;
        const IndexType id_3 = (alone_edge_local_id == 0) ? 2 : (alone_edge_local_id == 1) ? 0 : 1;
        this->CreateLineLoads(id_3, id_2, itElem, r_sub_model_part, p_properties, rMaximumConditionId);
    }
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ExtendPressureConditionProcess<3>::GeneratePressureLoads4WetNodes(
    const int PressureId,
    int& rMaximumConditionId,
    ModelPart::ElementsContainerType::ptr_iterator itElem
    )
{
    std::string sub_model_name;
    sub_model_name = mPressureName + "-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
    ModelPart::PropertiesType::Pointer p_properties = mpPropertiesVector[PressureId - 1];
    const int id = (*itElem)->Id();

    // We only create pressure loads when the surface is skin
    auto& r_elem_neigb = (*itElem)->GetValue(NEIGHBOUR_ELEMENTS);

    // Loop over the faces
    for (int i = 0; i < r_elem_neigb.size(); i++) {
        if (r_elem_neigb[i].Id() == id) { // it is skin face
            // The associated node to this skin face is excluded
            const int excluded_local_id_node = i;

            const IndexType id_1 = (excluded_local_id_node == 0) ? 3 : (excluded_local_id_node == 1) ? 0 : (excluded_local_id_node == 2) ? 3 : 0;
            const IndexType id_2 = (excluded_local_id_node == 0) ? 2 : (excluded_local_id_node == 1) ? 2 : (excluded_local_id_node == 2) ? 1 : 1;
            const IndexType id_3 = (excluded_local_id_node == 0) ? 1 : (excluded_local_id_node == 1) ? 3 : (excluded_local_id_node == 2) ? 0 : 2;
            this->CreatePressureLoads(id_1, id_2, id_3, itElem, r_sub_model_part, p_properties, rMaximumConditionId);
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/
template<>
void ExtendPressureConditionProcess<3>::GeneratePressureLoads3WetNodes(
    const int NonWetLocalIdNode,
    const int PressureId,
    int& rMaximumConditionId,
    ModelPart::ElementsContainerType::ptr_iterator itElem
    )
{
    std::string sub_model_name;
    sub_model_name = mPressureName + "-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mrModelPart.GetSubModelPart(sub_model_name);
    ModelPart::PropertiesType::Pointer p_properties = mpPropertiesVector[PressureId - 1];
    auto& r_geom = (*itElem)->GetGeometry();

    const IndexType id_1 = (NonWetLocalIdNode == 0) ? 0 : (NonWetLocalIdNode == 1) ? 1 : (NonWetLocalIdNode == 2) ? 2 : 3;
    const IndexType id_2 = (NonWetLocalIdNode == 0) ? 3 : (NonWetLocalIdNode == 1) ? 0 : (NonWetLocalIdNode == 2) ? 3 : 0;
    const IndexType id_3 = (NonWetLocalIdNode == 0) ? 2 : (NonWetLocalIdNode == 1) ? 2 : (NonWetLocalIdNode == 2) ? 1 : 1;
    const IndexType id_4 = (NonWetLocalIdNode == 0) ? 1 : (NonWetLocalIdNode == 1) ? 3 : (NonWetLocalIdNode == 2) ? 0 : 2;

    // We only create pressure loads when the surface is skin
    auto& r_elem_neigb = (*itElem)->GetValue(NEIGHBOUR_ELEMENTS);
    if (r_elem_neigb[NonWetLocalIdNode].Id() == (*itElem)->Id()) {     
        this->CreatePressureLoads(id_2, id_3, id_4, itElem, r_sub_model_part, p_properties, rMaximumConditionId);  
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<2>::CreateNewConditions()
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    int maximum_condition_id;
    this->GetMaximumConditionIdOnSubmodelPart(maximum_condition_id);
    r_process_info[INTERNAL_PRESSURE_ITERATION] = 0;

    // Loop over the elements (all active, the inactive have been removed in GeneratingDEM)
    for (auto it_elem = mrModelPart.Elements().ptr_begin(); it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
        if (!(*it_elem)->GetValue(SMOOTHING)) {
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
                r_process_info[INTERNAL_PRESSURE_ITERATION] = 10;
                (*it_elem)->SetValue(SMOOTHING, true);
            } else if (wet_nodes_counter == 3) {
                this->GetPressureId(it_elem, pressure_id);
                this->GenerateLineLoads3Nodes(pressure_id, maximum_condition_id, it_elem);
                r_process_info[INTERNAL_PRESSURE_ITERATION] = 10;
                (*it_elem)->SetValue(SMOOTHING, true);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<3>::CreateNewConditions()
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    int maximum_condition_id;
    this->GetMaximumConditionIdOnSubmodelPart(maximum_condition_id);
    r_process_info[INTERNAL_PRESSURE_ITERATION] = 0;

    // Loop over the elements (all active, the inactive have been removed in GeneratingDEM)
    for (auto it_elem = mrModelPart.Elements().ptr_begin(); it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
        if (!(*it_elem)->GetValue(SMOOTHING)) {
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
            if (wet_nodes_counter == 3) {
                this->GeneratePressureLoads3WetNodes(non_wet_local_id_node, pressure_id, maximum_condition_id, it_elem);
                r_process_info[INTERNAL_PRESSURE_ITERATION] = 10;
                (*it_elem)->SetValue(SMOOTHING, true);
            } else if (wet_nodes_counter == 4) {
                this->GeneratePressureLoads4WetNodes(pressure_id, maximum_condition_id, it_elem);
                r_process_info[INTERNAL_PRESSURE_ITERATION] = 10;
                (*it_elem)->SetValue(SMOOTHING, true);
            }
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::Execute()
{
    // We search the neighbours for the generation of line loads
    auto find_neigh = FindElementalNeighboursProcess(mrModelPart, TDim, 5);
    find_neigh.Execute();
    auto& r_process_info = mrModelPart.GetProcessInfo();

    // Remove previous line loads-> Only the 1st iteration
    if (r_process_info[INTERNAL_PRESSURE_ITERATION] == 1) {
        // We fill the properties vectors to be reassigned afterwards
        this->SavePreviousProperties();
        this->RemovePreviousLineLoads();
        this->ResetFlagOnElements();
    }
    // Generate the new ones
    this->CreateNewConditions();
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::SavePreviousProperties()
{
    const std::vector<std::string> submodel_parts_names = mrModelPart.GetSubModelPartNames();
    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        if (submodel_parts_names[i].substr(0, 8) == mPressureName.substr(0, 8)) {
            auto& r_sub_model_part = mrModelPart.GetSubModelPart(submodel_parts_names[i]);
            ModelPart::ConditionIterator it_cond = r_sub_model_part.ConditionsBegin();
            ModelPart::PropertiesType::Pointer p_properties = it_cond->pGetProperties();
            mpPropertiesVector.push_back(p_properties);
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::RemovePreviousLineLoads()
{
    // We remove only the line loads of all the SubModels
    const std::vector<std::string> submodel_parts_names = mrModelPart.GetSubModelPartNames();
    for (IndexType i = 0; i < submodel_parts_names.size(); ++i) {
        if (submodel_parts_names[i].substr(0, 8) == mPressureName.substr(0, 8)) {
            // Remove the line loads
            auto& r_sub_model = mrModelPart.GetSubModelPart(submodel_parts_names[i]);

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_sub_model.Conditions().size()); i++) {
                auto it_cond = r_sub_model.ConditionsBegin() + i;
                it_cond->Set(TO_ERASE, true);
            }
        }
    }
    mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::ResetFlagOnElements()
{
    auto it_elem_begin = mrModelPart.ElementsBegin();
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = it_elem_begin + i;
        it_elem->SetValue(SMOOTHING, false);
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void ExtendPressureConditionProcess<TDim>::GetPressureId(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    int& rPressureId
    )
{
    auto& r_geometry = (*itElem)->GetGeometry();
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        if (r_geometry[i].GetValue(PRESSURE_ID) != 0) {
            rPressureId = r_geometry[i].GetValue(PRESSURE_ID);
            break;
        }
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
