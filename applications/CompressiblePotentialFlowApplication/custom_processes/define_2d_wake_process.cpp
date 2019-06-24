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

// Project includes
#include "define_2d_wake_process.h"
#include "utilities/variable_utils.h"

namespace Kratos {

// Constructor for Define2DWakeProcess Process
Define2DWakeProcess::Define2DWakeProcess(ModelPart& rBodyModelPart, const double Tolerance)
    : Process(), mrBodyModelPart(rBodyModelPart), mTolerance(Tolerance)
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    root_model_part.CreateSubModelPart("trailing_edge_sub_model_part");
}

void Define2DWakeProcess::ExecuteInitialize()
{
    // Initialize element and node variables for the process
    //InitializeVariables();
    // Set the wake direction and normal for further computations
    SetWakeDirectionAndNormal();
    // Save the trailing edge for further computations
    SaveTrailingEdgeNode();
    // Check which elements are cut and mark them as wake
    MarkWakeElements();
    // Mark the elements touching the trailing edge from below as kutta
    MarkKuttaElements();
    // Mark the trailing edge element that is further downstream as wake
    MartkWakeTrailingEdgeElement();
}

const void Define2DWakeProcess::InitializeVariables()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    // Initialize element variables
    VariableUtils().SetNonHistoricalVariable(WAKE, false, root_model_part.Elements());
    VariableUtils().SetNonHistoricalVariable(KUTTA, false, root_model_part.Elements());
    VariableUtils().SetNonHistoricalVariable(TRAILING_EDGE, false, root_model_part.Elements());
    VariableUtils().SetFlag(STRUCTURE, false, root_model_part.Elements());
    VariableUtils().SetNonHistoricalVariable(ELEMENTAL_DISTANCES, ZeroVector(3), root_model_part.Elements());

    // Initialize node variables
    VariableUtils().SetNonHistoricalVariable(TRAILING_EDGE, false, root_model_part.Nodes());
    VariableUtils().SetVariable(DISTANCE, 1.0, root_model_part.Nodes());
}

// This function sets the wake direction and normal for further computations
const void Define2DWakeProcess::SetWakeDirectionAndNormal()
{
    // Reading the free_stream_velocity from the properties
    const auto free_stream_velocity = mrBodyModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
    KRATOS_ERROR_IF(free_stream_velocity.size() != 3)
        << "The free stream velocity should be a vector with 3 components!"
        << std::endl;

    // Computing the norm of the free_stream_velocity vector
    const double norm = sqrt(inner_prod(free_stream_velocity, free_stream_velocity));

    // The wake direction is the free stream direction
    mWakeDirection = free_stream_velocity / norm;

    // Rotating 90° to obtain the wake normal
    mWakeNormal(0) = -mWakeDirection(1);
    mWakeNormal(1) = mWakeDirection(0);
    mWakeNormal(2) = 0.0;
}

// This function finds and saves the trailing edge for further computations
const void Define2DWakeProcess::SaveTrailingEdgeNode()
{
    double max_x_coordinate = -1e30;
    NodeIteratorType trailing_edge_node;
    for (auto it_node = mrBodyModelPart.NodesBegin();
         it_node != mrBodyModelPart.NodesEnd(); ++it_node) {
        if (it_node->X() > max_x_coordinate) {
            max_x_coordinate = it_node->X();
            trailing_edge_node = it_node;
        }
    }
    trailing_edge_node->SetValue(TRAILING_EDGE, true);
    mTrailingEdgeNode = trailing_edge_node;
}

// This function checks which elements are cut by the wake and marks them as
// wake elements
const void Define2DWakeProcess::MarkWakeElements()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    bool potentially_wake = false;
    bool is_wake_element = false;
    BoundedVector<double, 3> nodal_distances_to_wake = ZeroVector(3);

    #pragma omp parallel for private(potentially_wake, is_wake_element, nodal_distances_to_wake)
    for (int i = 0; i < static_cast<int>(root_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = root_model_part.ElementsBegin() + i;

        // Check if the element is touching the trailing edge
        CheckIfTrailingEdgeElement(it_elem);

        // Elements downstream the trailing edge can be wake elements
        potentially_wake = CheckIfPotentiallyWakeElement(it_elem);

        if (potentially_wake) {
            // Compute the nodal distances of the element to the wake
            nodal_distances_to_wake = ComputeNodalDistancesToWake(it_elem);

            // Selecting the cut (wake) elements
            is_wake_element = CheckIfWakeElement(nodal_distances_to_wake);

            // Mark wake element and save their nodal distances to the wake
            if (is_wake_element) {
                #pragma omp critical
                {
                    it_elem->SetValue(WAKE, true);
                    it_elem->SetValue(ELEMENTAL_DISTANCES, nodal_distances_to_wake);
                    for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                        it_elem->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE) = nodal_distances_to_wake(i);
                    }
                }
            }
        }
    }
    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeElements();
}

// This function checks if the element is touching the trailing edge
const void Define2DWakeProcess::CheckIfTrailingEdgeElement(ElementIteratorType& rElement)
{
    // Loop over element nodes
    for (unsigned int i = 0; i < rElement->GetGeometry().size(); i++) {
        // Elements touching the trailing edge are trailing edge elements
        if (rElement->GetGeometry()[i].Id() == mTrailingEdgeNode->Id()) {
            MarkTrailingEdgeElement(rElement);
        }
    }
}

// This function marks the elements touching the trailing edge as trailing edge
const void Define2DWakeProcess::MarkTrailingEdgeElement(ElementIteratorType& rElement)
{
    #pragma omp critical
    {
    rElement->SetValue(TRAILING_EDGE, true);
    mTrailingEdgeElementsOrderedIds.push_back(rElement->Id());
    }
}

// This function selects the elements downstream the trailing edge as
// potentially wake elements
const bool Define2DWakeProcess::CheckIfPotentiallyWakeElement(ElementIteratorType& rElement)
{
    // Compute the distance from the trailing edge to the element's center
    BoundedVector<double, 3> distance_to_element_center = ZeroVector(3);
    distance_to_element_center =
        ComputeDistanceFromTrailingEdgeToPoint(rElement->GetGeometry().Center());

    // Compute the projection of the distance in the wake direction
    double projection_on_wake = inner_prod(distance_to_element_center, mWakeDirection);

    return projection_on_wake > 0.0;
}

// This function computes the distance of the element nodes to the wake
const BoundedVector<double, 3> Define2DWakeProcess::ComputeNodalDistancesToWake(ElementIteratorType& rElement)
{
    BoundedVector<double, 3> nodal_distances_to_wake = ZeroVector(3);
    BoundedVector<double, 3> distance_from_te_to_node = ZeroVector(3);
    double distance_to_wake = 1.0;
    for (unsigned int i = 0; i < rElement->GetGeometry().size(); i++) {
        // Compute the distance from the trailing edge to the node
        distance_from_te_to_node =
            ComputeDistanceFromTrailingEdgeToPoint(rElement->GetGeometry()[i]);

        // Compute the projection of the distance vector in the wake normal
        // direction
        distance_to_wake = inner_prod(distance_from_te_to_node, mWakeNormal);

        // Nodes laying on the wake have a positive distance
        if (std::abs(distance_to_wake) < mTolerance) {
            distance_to_wake = mTolerance;
        }
        nodal_distances_to_wake[i] = distance_to_wake;
    }
    return nodal_distances_to_wake;
}

// This function checks whether the element is cut by the wake
const bool Define2DWakeProcess::CheckIfWakeElement(BoundedVector<double, 3>& rNodalDistancesToWake)
{
    // Initialize counters
    unsigned int number_of_nodes_with_positive_distance = 0;
    unsigned int number_of_nodes_with_negative_distance = 0;

    // Count how many element nodes are above and below the wake
    for (unsigned int i = 0; i < rNodalDistancesToWake.size(); i++) {
        if (rNodalDistancesToWake(i) < 0.0) {
            number_of_nodes_with_negative_distance += 1;
        }
        else {
            number_of_nodes_with_positive_distance += 1;
        }
    }

    // Elements with nodes above and below the wake are wake elements
    return number_of_nodes_with_negative_distance > 0 &&
           number_of_nodes_with_positive_distance > 0;
}

// This function adds the trailing edge elements in the
// trailing_edge_sub_model_part
const void Define2DWakeProcess::AddTrailingEdgeElements()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    std::sort(mTrailingEdgeElementsOrderedIds.begin(),
              mTrailingEdgeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("trailing_edge_sub_model_part").AddElements(mTrailingEdgeElementsOrderedIds);
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
const void Define2DWakeProcess::MarkKuttaElements()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_sub_model_part");

    BoundedVector<double, 3> distance_to_element_center = ZeroVector(3);
    double distance_to_wake = 1.0;

    for (int i = 0;
         i < static_cast<int>(trailing_edge_sub_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem =
            trailing_edge_sub_model_part.ElementsBegin() + i;

        // Compute the distance from the element's center to the trailing edge
        distance_to_element_center =
            ComputeDistanceFromTrailingEdgeToPoint(it_elem->GetGeometry().Center());

        // Compute the projection of the distance vector in the wake normal
        // direction
        distance_to_wake = inner_prod(distance_to_element_center, mWakeNormal);

        if (distance_to_wake < 0.0) {
            it_elem->SetValue(KUTTA, true);
        }
    }
}

// This function finds the trailing edge element that is further downstream and
// marks it as wake trailing edge element. The rest of trailing edge elements
// are unassigned from the wake.
const void Define2DWakeProcess::MartkWakeTrailingEdgeElement()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_sub_model_part");

    for (int i = 0;
         i < static_cast<int>(trailing_edge_sub_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem =
            trailing_edge_sub_model_part.ElementsBegin() + i;
        if(it_elem->GetValue(WAKE)){
            // Trailing edge wake element
            if(CheckIfTrailingEdgeElementIsCutByWake(it_elem)){
                it_elem->Set(STRUCTURE);
                it_elem->SetValue(KUTTA, false);
            }
            //Rest of elements touching the trailing edge but not part of the wake
            else{
                it_elem->SetValue(WAKE, false);
            }
        }
    }
}

// This function checks if the element is cut by the wake
const bool Define2DWakeProcess::CheckIfTrailingEdgeElementIsCutByWake(ElementIteratorType& rElement)
{
    unsigned int number_of_nodes_with_negative_distance = 0;
    // REMINDER: In 3D the elemental_distances may not be match with the nodal
    // distances if CalculateDistanceToSkinProcess is used.
    auto nodal_distances_to_wake = rElement->GetValue(ELEMENTAL_DISTANCES);

    // Count how many element nodes are above and below the wake
    for (unsigned int i = 0; i < nodal_distances_to_wake.size(); i++) {
        if (nodal_distances_to_wake(i) < 0.0) {
            number_of_nodes_with_negative_distance += 1;
        }
    }

    return number_of_nodes_with_negative_distance == 1;
}

// This function computes the distance from the trailing edge to an input point
const BoundedVector<double, 3> Define2DWakeProcess::ComputeDistanceFromTrailingEdgeToPoint(Point InputPoint)
{
    BoundedVector<double, 3> distance_to_point = ZeroVector(3);

    distance_to_point(0) = InputPoint.X() - mTrailingEdgeNode->X();
    distance_to_point(1) = InputPoint.Y() - mTrailingEdgeNode->Y();

    return distance_to_point;
}

} // namespace Kratos.
