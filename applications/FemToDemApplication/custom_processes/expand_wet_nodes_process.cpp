//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//  Kratos default license:
//  kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/expand_wet_nodes_process.h"

namespace Kratos {

ExpandWetNodesProcess::ExpandWetNodesProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void ExpandWetNodesProcess::Execute() 
{
    int extrapolated_elements = 1;
    int pressure_id;
    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];
    while (extrapolated_elements > 0) {
        extrapolated_elements = 0;

        const auto it_elem_begin = mrModelPart.ElementsBegin();
        //#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
            auto it_elem = it_elem_begin + i;
            
            bool element_done = it_elem->GetValue(PRESSURE_EXPANDED);
            bool condition_is_active = true;
            if (it_elem->IsDefined(ACTIVE)) {
                condition_is_active = it_elem->Is(ACTIVE);
            }
            int number_of_wet_nodes;
            const bool has_wet_nodes = this->ElementHasWetNodes(it_elem, pressure_id, number_of_wet_nodes);
            if (number_of_wet_nodes > dimension - 1 && !condition_is_active && !element_done) {
                this->ExpandWetNodes(it_elem, pressure_id);
                extrapolated_elements++;
                it_elem->SetValue(PRESSURE_EXPANDED, true);
            }
        }
    }
    
    #pragma omp parallel for
    for(int i = 0; i<static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        it_elem->SetValue(PRESSURE_EXPANDED, false);
    }

    if (dimension == 2) 
        this->ExpandWetNodesIfTheyAreSkin();
}

/***********************************************************************************/
/***********************************************************************************/

bool ExpandWetNodesProcess::ElementHasWetNodes(
    ElementIterator itElem,
    int& rPressureId,
    int& rNumberOfWetNodes
    )
{
    rNumberOfWetNodes = 0;
    bool auxiliar = false;
    auto& r_geometry = itElem->GetGeometry();
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        const int pressure_id = r_geometry[i].GetValue(PRESSURE_ID);
        if (pressure_id != 0) {
            rNumberOfWetNodes++;
            rPressureId = pressure_id;
            auxiliar = true;
        }
    }
    return auxiliar;
}

/***********************************************************************************/
/***********************************************************************************/

void ExpandWetNodesProcess::ExpandWetNodes(
    ElementIterator itElem,
    const int PressureId
    )
{
    auto& r_geometry = itElem->GetGeometry();
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        r_geometry[i].SetValue(PRESSURE_ID, PressureId);
    }

    // Indicator to reconstruct the Pressure afterwards
    auto& r_process_info = mrModelPart.GetProcessInfo();
    r_process_info[RECONSTRUCT_PRESSURE_LOAD] = 1;
}

/***********************************************************************************/
/***********************************************************************************/

void ExpandWetNodesProcess::ExpandWetNodesIfTheyAreSkin()
{
    Parameters skin_process_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part"              : "SkinModelPart",
        "name_auxiliar_condition"               : "Condition",
        "list_model_parts_to_assign_conditions" : [],
        "echo_level"                            : 0
    })");

    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];
    
    // Evaluating according dimension
    if (dimension == 2) {
        auto skin_process = SkinDetectionProcess<2>(mrModelPart, skin_process_parameters);
        skin_process.Execute();
    } else {
        auto skin_process = SkinDetectionProcess<3>(mrModelPart, skin_process_parameters);
        skin_process.Execute();
    }
    
    auto& r_sub_model_part = mrModelPart.GetSubModelPart("SkinModelPart");

    auto it_node_begin = r_sub_model_part.NodesBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_sub_model_part.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(IS_SKIN, true);
    }

    int expanded_elements = 1;
    while (expanded_elements > 0) {
        expanded_elements = 0;

        //#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
            auto it_elem = mrModelPart.ElementsBegin() + i;
            auto& r_geometry = it_elem->GetGeometry();

            bool element_done = it_elem->GetValue(PRESSURE_EXPANDED);
            int number_of_wet_nodes, node_pressure_id;
            // Indicator to reconstruct the Pressure afterwards
            auto& r_process_info = mrModelPart.GetProcessInfo();

            if (this->ElementHasWetNodes(it_elem, node_pressure_id, number_of_wet_nodes) && !element_done) {
                // Loop over the nodes
                for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                   auto& r_node = r_geometry[i];
                   const int reference_pressure_id = node_pressure_id;
                   node_pressure_id = r_node.GetValue(PRESSURE_ID);

                    if (node_pressure_id == 0 && r_node.GetValue(IS_SKIN)) {
                       r_node.SetValue(PRESSURE_ID, reference_pressure_id);
                       expanded_elements++;
                       it_elem->SetValue(PRESSURE_EXPANDED, true);
                       r_process_info[RECONSTRUCT_PRESSURE_LOAD] = 1;
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i<static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        it_elem->SetValue(PRESSURE_EXPANDED, false);
    }
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos
