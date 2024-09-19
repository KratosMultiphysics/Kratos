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

    ExecuteInitialize();
    ComputeDistanceToWake();
    MarkWakeElements();
    ComputeTrailingEdgeNode();

    KRATOS_CATCH("");
}

void DefineEmbeddedWakeProcess::ExecuteInitialize() {

    KRATOS_TRY;

    block_for_each(mrModelPart.Elements(), [&](Element& rElem) {
        rElem.SetValue(WAKE, false);
    });
    block_for_each(mrModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        rNode.SetValue(WAKE_DISTANCE, 0.0);
        rNode.SetValue(WAKE, false);
        rNode.SetValue(KUTTA, false);
    });

    const auto free_stream_velocity = mrModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
    KRATOS_ERROR_IF(free_stream_velocity.size() != 3)
        << "The free stream velocity should be a vector with 3 components!"
        << std::endl;
    const double norm = norm_2(free_stream_velocity);
    KRATOS_ERROR_IF(norm < std::numeric_limits<double>::epsilon())
        << "The norm of the free stream velocity should be different than 0."
        << std::endl;
    // The wake direction is the free stream direction
    const auto wake_direction = free_stream_velocity / norm;
    array_1d<double, 3> wake_normal;

    wake_normal[0] = -wake_direction[1];
    wake_normal[1] = wake_direction[0];
    wake_normal[2] = 0.0;
    mrModelPart.GetRootModelPart().GetProcessInfo()[WAKE_NORMAL] = wake_normal;

    KRATOS_CATCH("");
}

void DefineEmbeddedWakeProcess::ComputeDistanceToWake(){

    KRATOS_TRY;

    CalculateDiscontinuousDistanceToSkinProcess<2> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();

    KRATOS_CATCH("");
}

void DefineEmbeddedWakeProcess::MarkWakeElements(){

    KRATOS_TRY;

    block_for_each(mrModelPart.Elements(), [&](Element& rElem)
    {
        auto& r_geometry = rElem.GetGeometry();
        BoundedVector<double, 3> nodal_distances_to_wake = rElem.GetValue(ELEMENTAL_DISTANCES);
        rElem.SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);

        // Selecting the cut (wake) elements
        const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(nodal_distances_to_wake);

        BoundedVector<double,3> geometry_distances;
        for(unsigned int i_node = 0; i_node< r_geometry.size(); i_node++){
            geometry_distances[i_node] = r_geometry[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);

        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element && rElem.Is(ACTIVE)) {
            rElem.SetValue(WAKE, true);
            if (is_embedded){
                rElem.Set(STRUCTURE, true);
                for (unsigned int i = 0; i < r_geometry.size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(KUTTA, true);
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].UnSetLock();
                }
            }
            else{
                for (unsigned int i = 0; i < r_geometry.size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].SetValue(WAKE, true);
                    r_geometry[i].UnSetLock();
                }
            }
        }
    });

    KRATOS_CATCH("");

}

ModelPart::NodeType::Pointer DefineEmbeddedWakeProcess::pGetTrailingEdgeNode(){

    for (auto& r_node : mrModelPart.Nodes()) {
        bool is_positive = r_node.GetValue(WAKE_DISTANCE) <= 0.0;
        bool is_wake = r_node.GetValue(WAKE);
        bool is_kutta = r_node.GetValue(KUTTA);
        if (is_positive && is_wake && is_kutta) {
            r_node.SetValue(TRAILING_EDGE, true);
            return &r_node;
        }
    }
    KRATOS_ERROR << "No trailing edge node was found!" << std::endl;
    return nullptr;

}

void DefineEmbeddedWakeProcess::ComputeTrailingEdgeNode(){

    KRATOS_TRY;

    auto p_node = pGetTrailingEdgeNode();

    std::vector<std::size_t> trailing_edge_node_list;
    trailing_edge_node_list.push_back(p_node->Id());

    if (mrModelPart.HasSubModelPart("trailing_edge_sub_model_part")){
        mrModelPart.RemoveSubModelPart("trailing_edge_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("trailing_edge_sub_model_part");

    std::sort(trailing_edge_node_list.begin(),
              trailing_edge_node_list.end());
    mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").AddNodes(trailing_edge_node_list);

    KRATOS_CATCH("");
}
}
