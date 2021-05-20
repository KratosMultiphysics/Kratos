//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Marc Nunez
//


#include "define_embedded_wake_process_3d.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"


namespace Kratos
{
// Constructor for DefineEmbeddedWake3DProcess Process
DefineEmbeddedWake3DProcess::DefineEmbeddedWake3DProcess(ModelPart& rModelPart,
                    ModelPart& rWakeModelPart
                ):
    Process(),
    mrModelPart(rModelPart),
    mrWakeModelPart(rWakeModelPart)
{}

void DefineEmbeddedWake3DProcess::Execute()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]!=3) << "DOMAIN_SIZE is not 3. DefineEmbeddedWake3DProcess is only implemented for 3D cases!" << std::endl;

    ComputeDistanceToWake();
    MarkWakeElements();
    ComputeTrailingEdgeNode();
    // MarkKuttaWakeElements();

    KRATOS_CATCH("");
}


void DefineEmbeddedWake3DProcess::ComputeDistanceToWake(){

    CalculateDiscontinuousDistanceToSkinProcess<3> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();
}

void DefineEmbeddedWake3DProcess::MarkWakeElements(){

    ModelPart& deactivated_model_part = mrModelPart.CreateSubModelPart("deactivated_model_part");
    std::vector<std::size_t> deactivated_ids;

    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i;

        BoundedVector<double, 4> nodal_distances_to_wake = it_elem->GetValue(ELEMENTAL_DISTANCES);
        // double epsilon = 1e-12;
        // for (std::size_t j = 0; j < nodal_distances_to_wake.size(); j++)
        // {
        //     if (std::abs(nodal_distances_to_wake[j]) < epsilon) {
        //         nodal_distances_to_wake[j] = epsilon;
        //     }

        // }
        it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);

        // Selecting the cut (wake) elements
        const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<3,4>(nodal_distances_to_wake);
        // const bool is_wake_element = this->Is(TO_SPLIT) ;

        BoundedVector<double,4> geometry_distances;
        for(unsigned int i_node = 0; i_node< 4; i_node++){
            geometry_distances[i_node] = it_elem->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<3,4>(geometry_distances);

        if (it_elem->Id() == 1194){
            KRATOS_WATCH(is_wake_element)
            KRATOS_WATCH(is_embedded)
            KRATOS_WATCH(nodal_distances_to_wake)
        }
        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element) {
            if (is_embedded){
                deactivated_ids.push_back(it_elem->Id());
                it_elem->Set(ACTIVE, false);
            }
            else{
                it_elem->SetValue(WAKE, true);
                auto r_geometry = it_elem->GetGeometry();
                for (unsigned int i_node = 0; i_node < it_elem->GetGeometry().size(); i_node++) {
                    r_geometry[i_node].SetLock();
                    r_geometry[i_node].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i_node));
                    r_geometry[i_node].UnSetLock();
                }
            }
        }
    }
    deactivated_model_part.AddElements(deactivated_ids);
}

void DefineEmbeddedWake3DProcess::ComputeTrailingEdgeNode(){

    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");

    // Find furthest deactivated element to the wake origin
    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;
        for(unsigned int i_node = 0; i_node < 4; i_node++){
            it_elem-> GetGeometry()[i_node].GetValue(TRAILING_EDGE) = true;
        }
    }

    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;
        for(unsigned int i_node = 0; i_node < 4; i_node++){
            const GlobalPointersVector<Element>& r_node_elem_neighbour = it_elem -> GetGeometry()[i_node].GetValue(NEIGHBOUR_ELEMENTS);
            for (std::size_t j = 0; j < r_node_elem_neighbour.size(); j++) {
                auto p_neighbour_element = r_node_elem_neighbour(j);
                const int wake = p_neighbour_element -> GetValue(WAKE);
                const bool active = p_neighbour_element -> Is(ACTIVE);
                int counter = 0;
                for (unsigned int sub_node = 0; sub_node < 4; sub_node++) {
                    if (p_neighbour_element -> GetGeometry()[sub_node].GetValue(TRAILING_EDGE) == 1) {
                        counter++;
                    }
                }

                if ((wake == 1) && active && (counter > 2)) {
                    p_neighbour_element -> Set(STRUCTURE);
                    auto wake_elemental_distances =   p_neighbour_element -> GetValue(WAKE_ELEMENTAL_DISTANCES);

                    // for (unsigned int iterator = 0; iterator<wake_elemental_distances.size(); iterator++){
                    //     if (wake_elemental_distances[iterator] > 0.0){
                    //         p_neighbour_element-> GetGeometry()[iterator].GetValue(TRAILING_EDGE) = 0;
                    //     }
                    // }

                }
                // if ((wake == 1) && active && (counter < 3)) {
                //     auto wake_elemental_distances =   p_neighbour_element -> GetValue(WAKE_ELEMENTAL_DISTANCES);
                //     int npos = 0;
                //     int non_te_nodes = 0;
                //     int nneg =0;
                //     for (unsigned int iterator = 0; iterator<wake_elemental_distances.size(); iterator++){
                //         if (p_neighbour_element->GetGeometry()[i_node].GetValue(TRAILING_EDGE) == 1) {
                //             non_te_nodes++;
                //         }
                //         if (wake_elemental_distances[iterator] < 0.0){
                //             nneg++;
                //         }else{
                //             npos++;
                //         }
                //     }
                //     if (nneg > non_te_nodes - 1 ){
                //         p_neighbour_element -> SetValue(WAKE, false);
                //         p_neighbour_element -> SetValue(KUTTA, true);
                //     }
                // }
            }
        }
    }

    mrModelPart.RemoveSubModelPart("deactivated_model_part");
}

void DefineEmbeddedWake3DProcess::MarkKuttaWakeElements(){

    // // Find elements that touch the furthes deactivated element and that are part of the wake.
    // for (std::size_t i = 0; i < mKuttaWakeElementCandidates.size(); i++)
    // {
    //     auto& r_geometry = mKuttaWakeElementCandidates[i].GetGeometry();
    //     if (mKuttaWakeElementCandidates[i].GetValue(WAKE) && mKuttaWakeElementCandidates[i].Is(ACTIVE)) {
    //         for (std::size_t i_node= 0; i_node < r_geometry.size(); i_node++) {
    //             if(r_geometry[i_node].GetValue(TRAILING_EDGE)){
    //                 mKuttaWakeElementCandidates[i].Set(STRUCTURE);
    //                 break;
    //             }
    //         }
    //     }

    // }
}
}// Namespace Kratos
