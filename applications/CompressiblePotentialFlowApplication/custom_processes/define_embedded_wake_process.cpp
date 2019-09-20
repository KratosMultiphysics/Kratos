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


#include "define_embedded_wake_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"


namespace Kratos
{
// Constructor for DefineEmbeddedWakeProcess Process
DefineEmbeddedWakeProcess::DefineEmbeddedWakeProcess(ModelPart& rModelPart,
                    ModelPart& rWakeModelPart
                ):
    Process(),
    mrModelPart(rModelPart),
    mrWakeModelPart(rWakeModelPart)
{}

void DefineEmbeddedWakeProcess::Execute()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]>2) << "DOMAIN_SIZE is greater than 2. DefineEmbeddedWakeProcess is only implemented for 2D cases!" << std::endl;

    ComputeDistanceToWake();
    MarkWakeElements();
    ComputeTrailingEdgeNode();
    MarkKuttaWakeElements();

    KRATOS_CATCH("");
}


void DefineEmbeddedWakeProcess::ComputeDistanceToWake(){

    CalculateDiscontinuousDistanceToSkinProcess<2> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();
}

void DefineEmbeddedWakeProcess::MarkWakeElements(){

    ModelPart& deactivated_model_part = mrModelPart.CreateSubModelPart("deactivated_model_part");
    std::vector<std::size_t> deactivated_ids;
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i;

        BoundedVector<double, 3> nodal_distances_to_wake = it_elem->GetValue(ELEMENTAL_DISTANCES);
        it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);

        // Selecting the cut (wake) elements
        const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(nodal_distances_to_wake);

        BoundedVector<double,3> geometry_distances;
        for(unsigned int i_node = 0; i_node<3; i_node++){
            geometry_distances[i_node] = it_elem->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);

        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element) {
            if (is_embedded){
                deactivated_ids.push_back(it_elem->Id());
                it_elem->Set(ACTIVE, false);
            }
            else{
                it_elem->SetValue(WAKE, true);
                auto r_geometry = it_elem->GetGeometry();
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].UnSetLock();
                }
            }
        }
    }
    deactivated_model_part.AddElements(deactivated_ids);
}

void DefineEmbeddedWakeProcess::ComputeTrailingEdgeNode(){

    double max_distance = 0.0;
    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");
    Element::Pointer p_max_elem;

    auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];

    // Find furthest deactivated element to the wake origin
    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;

        BoundedVector<double,2> distance_vector;
        distance_vector[0] = wake_origin[0] - it_elem->GetGeometry().Center().X();
        distance_vector[1] = wake_origin[1] - it_elem->GetGeometry().Center().Y();
        double norm = norm_2(distance_vector);
        if(norm>max_distance){
            max_distance = norm;
            p_max_elem = mrModelPart.pGetElement(it_elem->Id());
        }
    }

    // Mark nodes of the furthest deactivated element and store its neighbour elements
    for (unsigned int i_node= 0; i_node < p_max_elem->GetGeometry().size(); i_node++) {
        p_max_elem->GetGeometry()[i_node].SetValue(TRAILING_EDGE,true);
        const GlobalPointersVector<Element>& r_node_elem_candidates = p_max_elem -> GetGeometry()[i_node].GetValue(NEIGHBOUR_ELEMENTS);
        for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
            mKuttaWakeElementCandidates.push_back(r_node_elem_candidates(j));
        }
    }

    mrModelPart.RemoveSubModelPart("deactivated_model_part");
}

void DefineEmbeddedWakeProcess::MarkKuttaWakeElements(){

    // Find elements that touch the furthes deactivated element and that are part of the wake.
    for (std::size_t i = 0; i < mKuttaWakeElementCandidates.size(); i++)
    {
        auto& r_geometry = mKuttaWakeElementCandidates[i].GetGeometry();
        if (mKuttaWakeElementCandidates[i].GetValue(WAKE) && mKuttaWakeElementCandidates[i].Is(ACTIVE)) {
            for (std::size_t i_node= 0; i_node < r_geometry.size(); i_node++) {
                if(r_geometry[i_node].GetValue(TRAILING_EDGE)){
                    mKuttaWakeElementCandidates[i].Set(STRUCTURE);
                    break;
                }
            }
        }

    }
}
}// Namespace Kratos
