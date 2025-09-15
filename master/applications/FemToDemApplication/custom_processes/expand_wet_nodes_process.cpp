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
    auto& r_process_info = mrModelPart.GetProcessInfo();
    std::size_t dimension = r_process_info[DOMAIN_SIZE];
    mPressureName = (dimension == 2) ? "Normal_Load" : "Pressure_Load";
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
            
            const bool element_done = it_elem->GetValue(PRESSURE_EXPANDED);
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

    const auto elem_begin = mrModelPart.ElementsBegin();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = elem_begin + i;
        it_elem->SetValue(PRESSURE_EXPANDED, false);
    }

    if (dimension == 2) {
        this->ExpandWetNodesIfTheyAreSkin();
        // this->ExpandWetNodesWithLatestPressureId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool ExpandWetNodesProcess::ElementHasWetNodes(
    ElementIterator itElem,
    int& rPressureId,
    int& rNumberOfWetNodes
    )
{   // Assigns the maximum Pressure Id 
    rNumberOfWetNodes = 0;
    rPressureId = 0;
    bool auxiliar = false;
    auto& r_geometry = itElem->GetGeometry();
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        const int pressure_id = r_geometry[i].GetValue(PRESSURE_ID);
        if (pressure_id != 0) {
            rNumberOfWetNodes++;
            rPressureId = (pressure_id > rPressureId) ? pressure_id : rPressureId;
            auxiliar = true;
        }
    }
    return auxiliar;
}

/***********************************************************************************/
/***********************************************************************************/

bool ExpandWetNodesProcess::ElementHasWetNodes2(
    ElementIterator itElem,
    int& rPressureId,
    int& rNumberOfWetNodes
    )
{   // Assigns the dominant Pressure Id 
    rNumberOfWetNodes = 0;
    rPressureId = 0;
    bool auxiliar = false;
    auto& r_geometry = itElem->GetGeometry();
    const auto number_points = r_geometry.PointsNumber();
    Vector pressures_ids(number_points);
    for (IndexType i = 0; i < number_points; ++i) {
        const int pressure_id = r_geometry[i].GetValue(PRESSURE_ID);

        if (pressure_id != 0) {
            pressures_ids(i) = pressure_id;
            rNumberOfWetNodes++;
            rPressureId = (pressure_id > rPressureId) ? pressure_id : rPressureId;
            auxiliar = true;
        }
    }

    if (rNumberOfWetNodes == number_points) {
        if (pressures_ids(0) == pressures_ids(1)) {
            rPressureId = pressures_ids(0);
        } else if (pressures_ids(0) == pressures_ids(2)) {
            rPressureId = pressures_ids(0);
        } else if (pressures_ids(1) == pressures_ids(2)) {
            rPressureId = pressures_ids(1);
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
    double initial_vol_pressure, vol_pressure;
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        if (r_geometry[i].GetValue(PRESSURE_INITIAL_VOLUME) != 0.0) {
            initial_vol_pressure = r_geometry[i].GetValue(PRESSURE_INITIAL_VOLUME);
            vol_pressure = r_geometry[i].GetValue(PRESSURE_VOLUME);
            break;
        }
    }

    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        r_geometry[i].SetValue(PRESSURE_ID, PressureId);
        r_geometry[i].SetValue(PRESSURE_INITIAL_VOLUME, initial_vol_pressure);
        r_geometry[i].SetValue(PRESSURE_VOLUME, vol_pressure);
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
        auto it_elem_begin = mrModelPart.ElementsBegin();

        //#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
            auto it_elem = it_elem_begin + i;
            auto& r_geometry = it_elem->GetGeometry();

            const bool element_done = it_elem->GetValue(PRESSURE_EXPANDED);
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

void ExpandWetNodesProcess::ExpandWetNodesWithLatestPressureId()
{
    int expanded_elements = 1;
    auto& r_process_info = mrModelPart.GetProcessInfo();
    const std::size_t dimension = r_process_info[DOMAIN_SIZE];
    bool volumes_updated = false;

    while (expanded_elements > 0) {

        expanded_elements = 0;
        const auto it_elem_begin = mrModelPart.ElementsBegin();

        //#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
            auto it_elem = it_elem_begin + i;
            auto& r_geometry = it_elem->GetGeometry();

            bool condition_is_active = true;
            if (it_elem->IsDefined(ACTIVE)) {
                condition_is_active = it_elem->Is(ACTIVE);
            }

            const bool element_done = it_elem->GetValue(PRESSURE_EXPANDED);
            int number_of_wet_nodes, node_pressure_id, pressure_id;
            // Indicator to reconstruct the Pressure afterwards
            auto& r_process_info = mrModelPart.GetProcessInfo();

            int non_pressure_node;
            for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                auto& r_node = r_geometry[i];
                node_pressure_id = r_node.GetValue(PRESSURE_ID);
                if (node_pressure_id == 0)
                    non_pressure_node = i;
            }


            if (this->ElementHasWetNodes(it_elem, pressure_id, number_of_wet_nodes) && !element_done && condition_is_active) {
                // Loop over the nodes
                for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                   auto& r_node = r_geometry[i];
                   const int reference_pressure_id = pressure_id;
                   node_pressure_id = r_node.GetValue(PRESSURE_ID);

                    if (node_pressure_id != 0 && 
                        node_pressure_id < reference_pressure_id && 
                        number_of_wet_nodes < r_geometry.PointsNumber()) {

                        std::string old_blast_sub_model_name, new_blast_sub_model_name;
                        old_blast_sub_model_name = mPressureName + "-auto-" + std::to_string(node_pressure_id);
                        new_blast_sub_model_name = mPressureName + "-auto-" + std::to_string(reference_pressure_id);

                        auto& r_new_blast_sub_model_part = mrModelPart.GetSubModelPart(new_blast_sub_model_name);
                        auto& r_old_blast_sub_model_part = mrModelPart.GetSubModelPart(old_blast_sub_model_name);
                        const auto it_old_blast_node_begin = r_old_blast_sub_model_part.NodesBegin();
                        const auto it_new_blast_node_begin = r_new_blast_sub_model_part.NodesBegin();

                        const double old_pressure_volume = it_old_blast_node_begin->GetValue(PRESSURE_VOLUME);
                        const double new_pressure_volume = it_new_blast_node_begin->GetValue(PRESSURE_VOLUME);

                        const double new_pressure_volume_initial = it_new_blast_node_begin->GetValue(PRESSURE_INITIAL_VOLUME);
                        
                        if (!volumes_updated) {
                            const auto it_node_begin = mrModelPart.NodesBegin();
                            #pragma omp parallel for // We update the volume as the sum of the 2 (merging volumes)
                            for (int j = 0; j < static_cast<int>(mrModelPart.Nodes().size()); j++) {
                                auto it_node = it_node_begin + j;
                                if (it_node->GetValue(PRESSURE_ID) == node_pressure_id ||
                                    it_node->GetValue(PRESSURE_ID) == reference_pressure_id) {
                                        // it_node->SetValue(PRESSURE_VOLUME, old_pressure_volume + new_pressure_volume);
                                        // it_node->SetValue(PRESSURE_INITIAL_VOLUME, new_pressure_volume_initial);
                                }
                            }
                        }


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
