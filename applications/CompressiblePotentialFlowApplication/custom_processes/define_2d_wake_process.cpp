//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

// System includes

// External includes

// Project includes
#include "define_2d_wake_process.h"

namespace Kratos {

/***********************************************************************************/
/***********************************************************************************/

void Define2DWakeProcess::ExecuteInitialize()
{
    std::cout << "Entering ExecuteInitialize" << std::endl;
    SetWakeDirectionAndNormal();
    // Save the trailing edge for further computations
    SaveTrailingEdgeNode();
    // Check which elements are cut and mark them as wake
    MarkWakeElements();
}

void Define2DWakeProcess::SetWakeDirectionAndNormal()
{
    auto free_stream_velocity = mrModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
    KRATOS_ERROR_IF(free_stream_velocity.size() != 3)
        << "The free stream velocity should be a vector with 3 components!" << std::endl;

    double norm = sqrt(inner_prod(free_stream_velocity, free_stream_velocity));
    mWakeDirection = free_stream_velocity/norm;
    mWakeNormal(0) = -mWakeDirection(1);
    mWakeNormal(1) = mWakeDirection(0);
    mWakeNormal(2) = 0.0;
}

BoundedVector<double, 3> Define2DWakeProcess::GetWakeDirection()
{
    return mWakeDirection;
}

void Define2DWakeProcess::SaveTrailingEdgeNode()
{
    std::cout << "Entering SaveTrailingEdgeNode" << std::endl;
    double max_x_coordinate = -1e30;
    NodeIteratorType trailing_edge_node;
    for (auto it_node = mrModelPart.NodesBegin();
         it_node != mrModelPart.NodesEnd(); ++it_node) {
        if (it_node->X() > max_x_coordinate) {
            max_x_coordinate = it_node->X();
            trailing_edge_node = it_node;
        }
    }
    trailing_edge_node->SetValue(TRAILING_EDGE, true);
    SetTrailingEdgeNode(trailing_edge_node);
}

void Define2DWakeProcess::SetTrailingEdgeNode(NodeIteratorType TrailingEdgeNode)
{
    mTrailingEdgeNode = TrailingEdgeNode;
}

ModelPart::NodeIterator Define2DWakeProcess::GetTrailingEdgeNode()
{
    return mTrailingEdgeNode;
}

// This function checks which elements are cut by the wake
// and marks them as wake elements
void Define2DWakeProcess::MarkWakeElements()
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    bool potentially_wake = false;
    bool is_wake_element = false;
    BoundedVector<double, 3> distances_to_wake = ZeroVector(3);

    //#pragma omp parallel for
    for (int i=0; i< static_cast<int>(root_model_part.Elements().size()); i++){
        ModelPart::ElementIterator it_elem = root_model_part.ElementsBegin() + i;

        // Elements downstream the trailing edge can be wake elements
        potentially_wake = CheckIfPotentiallyWakeElement(it_elem);

        if(potentially_wake){
            // Compute the nodal distances of the element to the wake
            distances_to_wake = ComputeDistancesToWake(it_elem);

            // Selecting the cut (wake) elements
            is_wake_element = CheckIfWakeElement(distances_to_wake);

            if(is_wake_element){
                it_elem->SetValue(WAKE, true);
                it_elem->SetValue(ELEMENTAL_DISTANCES, distances_to_wake);
            }
        }
    }
}

// This function selects the elements downstream the
// trailing edge as potentially wake elements
bool Define2DWakeProcess::CheckIfPotentiallyWakeElement(ElementIteratorType& rElement)
{
    auto trailing_edge_node = GetTrailingEdgeNode();

    // Compute the distance from the element's center to
    // the trailing edge
    BoundedVector<double, 3> distance_to_te = ZeroVector(3);
    distance_to_te(0) =
        rElement->GetGeometry().Center().X() - trailing_edge_node->X();
    distance_to_te(1) =
        rElement->GetGeometry().Center().Y() - trailing_edge_node->Y();

    // Compute the projection of the distance in the wake direction
    auto wake_direction = GetWakeDirection();
    double projection_on_wake = inner_prod(distance_to_te, wake_direction);

    return projection_on_wake > 0.0;
}

BoundedVector<double, 3> Define2DWakeProcess::ComputeDistancesToWake(ElementIteratorType& rElement)
{
    BoundedVector<double, 3> nodal_distances_to_wake = ZeroVector(3);
    BoundedVector<double, 3> distance_to_te = ZeroVector(3);
    double distance_to_wake = 1.0;
    for (unsigned int i = 0; i < rElement->GetGeometry().size(); i++) {
        distance_to_te(0) =
            rElement->GetGeometry().pGetPoint(i)->X() - mTrailingEdgeNode->X();
        distance_to_te(1) =
            rElement->GetGeometry().pGetPoint(i)->Y() - mTrailingEdgeNode->Y();

        distance_to_wake = inner_prod(distance_to_te, mWakeNormal);

        // Nodes laying on the wake have a positive distance
        if (std::abs(distance_to_wake) < mTolerance) {
            distance_to_wake = mTolerance;
        }
        nodal_distances_to_wake[i] = distance_to_wake;
    }
    return nodal_distances_to_wake;
}

bool Define2DWakeProcess::CheckIfWakeElement(BoundedVector<double, 3>& rDistancesToWake)
{
    unsigned int number_of_nodes_with_positive_distance = 0;
    unsigned int number_of_nodes_with_negative_distance = 0;

    for (unsigned int i = 0; i < rDistancesToWake.size(); i++) {
        if (rDistancesToWake(i) < 0.0) {
            number_of_nodes_with_negative_distance += 1;
        }
        else {
            number_of_nodes_with_positive_distance += 1;
        }
    }

    return number_of_nodes_with_negative_distance > 0 &&
           number_of_nodes_with_positive_distance > 0;
}

// void Define2DWakeProcess::SetTrailingEdgeNode(NodeType node)
// {
//     mTrailingEdgeNode = node;
// }

}  // namespace Kratos.



