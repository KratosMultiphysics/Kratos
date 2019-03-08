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
    while (extrapolated_elements > 0) {
        extrapolated_elements = 0;
        for (auto it_elem = mrModelPart.Elements().ptr_begin(); it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
            
            bool element_done = (*it_elem)->GetValue(PRESSURE_EXPANDED);
            bool condition_is_active = true;
            if ((*it_elem)->IsDefined(ACTIVE)) {
                condition_is_active = (*it_elem)->Is(ACTIVE);
            }
            int number_of_wet_nodes;
            bool has_wet_nodes = this->ElementHasWetNodes(it_elem, pressure_id, number_of_wet_nodes);
            if (number_of_wet_nodes > 1 && condition_is_active == false && element_done == false) {
                this->ExpandWetNodes(it_elem, pressure_id);
                extrapolated_elements++;
                (*it_elem)->SetValue(PRESSURE_EXPANDED, true);
            }
        }
    }
    for (auto it_elem = mrModelPart.Elements().ptr_begin(); it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
        (*it_elem)->SetValue(PRESSURE_EXPANDED, false);
    }

	this->ExpandWetNodesIfTheyAreSkin();



    for (auto it_node = mrModelPart.Nodes().ptr_begin(); it_node != mrModelPart.Nodes().ptr_end(); ++it_node) {
        if ((*it_node)->GetValue(PRESSURE_ID) != 0) {
            KRATOS_WATCH((*it_node)->Id())
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool ExpandWetNodesProcess::ElementHasWetNodes(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    int& rPressureId,
    int& rNumberOfWetNodes
    )
{
    rNumberOfWetNodes = 0;
    bool auxiliar = false;
    auto& r_geometry = (*itElem)->GetGeometry();
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
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    const int PressureId
    )
{
    auto& r_geometry = (*itElem)->GetGeometry();
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

    auto skin_process = SkinDetectionProcess<2>(mrModelPart, skin_process_parameters); // TODO what's up in 3d?
    skin_process.Execute();
	auto& r_sub_model_part = mrModelPart.GetSubModelPart("SkinModelPart");




    // for (auto it_node = r_sub_model_part.Nodes().ptr_begin(); it_node != r_sub_model_part.Nodes().ptr_end(); ++it_node) {
    //     KRATOS_WATCH((*it_node)->Id())
    // }

    // for (auto it_node = mrModelPart.Nodes().ptr_begin(); it_node != mrModelPart.Nodes().ptr_end(); ++it_node) {
    //     if ((*it_node)->GetValue(PRESSURE_ID) != 0) {
    //         KRATOS_WATCH((*it_node)->Id())
    //     }
    // }





    // Let's assign the flag IS_SKIN to the nodes
    for (auto it_node = r_sub_model_part.Nodes().ptr_begin(); it_node != r_sub_model_part.Nodes().ptr_end(); ++it_node) {
        (*it_node)->SetValue(IS_SKIN, true);
    }

    int expanded_elements = 1;
    while (expanded_elements > 0) {
        expanded_elements = 0;
        for (auto it_elem = mrModelPart.Elements().ptr_begin(); it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
            auto& r_geometry = (*it_elem)->GetGeometry();

            // if ((*it_elem)->Id() == 158) {
            //     int number_of_wet_nodes, pressure_id;
            //     bool aux = this->ElementHasWetNodes(it_elem, pressure_id, number_of_wet_nodes);
            //     KRATOS_WATCH(aux)
            //     for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
            //         auto& r_node = r_geometry[i];
            //         pressure_id = r_node.GetValue(PRESSURE_ID);
            //         KRATOS_WATCH(pressure_id)
            //         KRATOS_WATCH(r_node.GetValue(IS_SKIN))
            //         if (pressure_id == 0 && r_node.GetValue(IS_SKIN)) {
            //             // KRATOS_WATCH(pressure_id)
            //             // KRATOS_WATCH(r_node.GetValue(IS_SKIN))
            //         }
            //     }
            // }

            bool element_done = (*it_elem)->GetValue(PRESSURE_EXPANDED);
            int number_of_wet_nodes, pressure_id;
            // Indicator to reconstruct the Pressure afterwards
            auto& r_process_info = mrModelPart.GetProcessInfo();

            if (this->ElementHasWetNodes(it_elem, pressure_id, number_of_wet_nodes) && element_done == false) {
                // Loop over the nodes
                for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                    auto& r_node = r_geometry[i];
                    pressure_id = r_node.GetValue(PRESSURE_ID);
                    if (pressure_id == 0 && r_node.GetValue(IS_SKIN)) {
                        r_node.SetValue(PRESSURE_ID, pressure_id);
                        expanded_elements++;
                        (*it_elem)->SetValue(PRESSURE_EXPANDED, true);
                        r_process_info[RECONSTRUCT_PRESSURE_LOAD] = 1;


                        // if ((*it_elem)->Id() == 158) {
                        //     KRATOS_WATCH(pressure_id)
                        //     KRATOS_WATCH(r_node.GetValue(IS_SKIN))
                        // }



                    }
                }
            }
        }
    }
    for (auto it_elem = mrModelPart.Elements().ptr_begin(); it_elem != mrModelPart.Elements().ptr_end(); ++it_elem) {
        (*it_elem)->SetValue(PRESSURE_EXPANDED, false);
    }
}



} // namespace Kratos