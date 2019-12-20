//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

// Project includes
#include "define_2d_wake_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {

// Constructor for Define2DWakeProcess Process
Define2DWakeProcess::Define2DWakeProcess(ModelPart& rBodyModelPart, const double Tolerance)
    : Process(), mrBodyModelPart(rBodyModelPart), mTolerance(Tolerance)
{}

void Define2DWakeProcess::ExecuteInitialize()
{
    // Initialize submodelparts
    InitializeTrailingEdgeSubModelpart();
    InitializeWakeSubModelpart();
    // Set the wake direction and normal for further computations
    SetWakeDirectionAndNormal();
    // Save the trailing edge for further computations
    SaveTrailingEdgeNode();
    // Check which elements are cut and mark them as wake
    MarkWakeElements();
    // Mark the elements touching the trailing edge from below as kutta
    MarkKuttaElements();
    // Mark the trailing edge element that is further downstream as wake
    MarkWakeTrailingEdgeElement();
}

// This function initializes the variables and removes all of the elements of
// the trailing edge submodelpart
void Define2DWakeProcess::InitializeTrailingEdgeSubModelpart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("trailing_edge_sub_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // trailing_edge_sub_model_part
        ModelPart& trailing_edge_sub_model_part =
            root_model_part.GetSubModelPart("trailing_edge_sub_model_part");

        for (auto& r_element : trailing_edge_sub_model_part.Elements()){
            r_element.SetValue(TRAILING_EDGE, false);
            r_element.SetValue(KUTTA, false);
            r_element.Reset(STRUCTURE);
            r_element.Set(TO_ERASE, true);
        }
        trailing_edge_sub_model_part.RemoveElements(TO_ERASE);
    }
    else{
        // Creating the trailing_edge_sub_model_part
        root_model_part.CreateSubModelPart("trailing_edge_sub_model_part");
    }
}

// This function initializes the variables and removes all of the elements of
// the wake submodelpart
void Define2DWakeProcess::InitializeWakeSubModelpart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("wake_sub_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // wake_sub_model_part
        ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_sub_model_part");

        for (auto& r_element : wake_sub_model_part.Elements()){
            r_element.SetValue(WAKE, false);
            r_element.SetValue(WAKE_ELEMENTAL_DISTANCES, ZeroVector(3));
            r_element.Set(TO_ERASE, true);
        }
        wake_sub_model_part.RemoveElements(TO_ERASE);
    }
    else{
        // Creating the wake_sub_model_part
        root_model_part.CreateSubModelPart("wake_sub_model_part");
    }
}

// This function sets the wake direction and normal for further computations
void Define2DWakeProcess::SetWakeDirectionAndNormal()
{
    // Reading the free_stream_velocity from the properties
    const auto free_stream_velocity = mrBodyModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
    KRATOS_ERROR_IF(free_stream_velocity.size() != 3)
        << "The free stream velocity should be a vector with 3 components!"
        << std::endl;

    // Computing the norm of the free_stream_velocity vector
    const double norm = std::sqrt(inner_prod(free_stream_velocity, free_stream_velocity));

    const double eps = std::numeric_limits<double>::epsilon();
    KRATOS_ERROR_IF(norm < eps)
        << "The norm of the free stream velocity should be different than 0."
        << std::endl;

    // The wake direction is the free stream direction
    mWakeDirection = free_stream_velocity / norm;

    // Rotating 90° to obtain the wake normal
    mWakeNormal(0) = -mWakeDirection(1);
    mWakeNormal(1) = mWakeDirection(0);
    mWakeNormal(2) = 0.0;
}

// This function finds and saves the trailing edge for further computations
void Define2DWakeProcess::SaveTrailingEdgeNode()
{
    KRATOS_ERROR_IF(mrBodyModelPart.NumberOfNodes() == 0) << "There are no nodes in the body_model_part!"<< std::endl;

    double max_x_coordinate = std::numeric_limits<double>::lowest();
    auto p_trailing_edge_node = &*mrBodyModelPart.NodesBegin();

    for (auto& r_node : mrBodyModelPart.Nodes()) {
        if (r_node.X() > max_x_coordinate) {
            max_x_coordinate = r_node.X();
            p_trailing_edge_node = &r_node;
        }
    }

    p_trailing_edge_node->SetValue(TRAILING_EDGE, true);
    mpTrailingEdgeNode = p_trailing_edge_node;
}

// This function checks which elements are cut by the wake and marks them as
// wake elements
void Define2DWakeProcess::MarkWakeElements()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    std::vector<std::size_t> wake_elements_ordered_ids;

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(root_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = root_model_part.ElementsBegin() + i;

        // Check if the element is touching the trailing edge
        CheckIfTrailingEdgeElement(*it_elem);

        // Elements downstream the trailing edge can be wake elements
        bool potentially_wake = CheckIfPotentiallyWakeElement(*it_elem);

        if (potentially_wake) {
            // Compute the nodal distances of the element to the wake
            BoundedVector<double, 3> nodal_distances_to_wake = ComputeNodalDistancesToWake(*it_elem);

            // Selecting the cut (wake) elements
            const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(nodal_distances_to_wake);;

            // Mark wake element and save their nodal distances to the wake
            if (is_wake_element) {
                it_elem->SetValue(WAKE, true);
                it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);
                #pragma omp critical
                {
                    wake_elements_ordered_ids.push_back(it_elem->Id());
                }
                auto r_geometry = it_elem->GetGeometry();
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].UnSetLock();
                }
            }
        }
    }
    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeAndWakeElements(wake_elements_ordered_ids);
}

// This function checks if the element is touching the trailing edge
void Define2DWakeProcess::CheckIfTrailingEdgeElement(Element& rElement)
{
    // Loop over element nodes
    for (unsigned int i = 0; i < rElement.GetGeometry().size(); i++) {
        // Elements touching the trailing edge are trailing edge elements
        if (rElement.GetGeometry()[i].Id() == mpTrailingEdgeNode->Id()) {
            rElement.SetValue(TRAILING_EDGE, true);
            #pragma omp critical
            {
                mTrailingEdgeElementsOrderedIds.push_back(rElement.Id());
            }
        }
    }
}

// This function selects the elements downstream the trailing edge as
// potentially wake elements
bool Define2DWakeProcess::CheckIfPotentiallyWakeElement(const Element& rElement) const
{
    // Compute the distance from the trailing edge to the element's center
    const auto distance_to_element_center =
        ComputeDistanceFromTrailingEdgeToPoint(rElement.GetGeometry().Center());

    // Compute the projection of the distance in the wake direction
    return inner_prod(distance_to_element_center, mWakeDirection) > 0.0 ? true : false;
}

// This function computes the distance of the element nodes to the wake
const BoundedVector<double, 3> Define2DWakeProcess::ComputeNodalDistancesToWake(const Element& rElement) const
{
    BoundedVector<double, 3> nodal_distances_to_wake = ZeroVector(3);

    for (unsigned int i = 0; i < rElement.GetGeometry().size(); i++) {
        // Compute the distance from the trailing edge to the node
        const auto distance_from_te_to_node =
            ComputeDistanceFromTrailingEdgeToPoint(rElement.GetGeometry()[i]);

        // Compute the projection of the distance vector in the wake normal
        // direction
        double distance_to_wake = inner_prod(distance_from_te_to_node, mWakeNormal);

        // Nodes laying on the wake have a positive distance
        if (std::abs(distance_to_wake) < mTolerance) {
            distance_to_wake = mTolerance;
        }
        nodal_distances_to_wake[i] = distance_to_wake;
    }
    return nodal_distances_to_wake;
}

// This function adds the trailing edge elements in the
// trailing_edge_sub_model_part
void Define2DWakeProcess::AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds)
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    std::sort(rWakeElementsOrderedIds.begin(),
              rWakeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("wake_sub_model_part").AddElements(rWakeElementsOrderedIds);

    std::sort(mTrailingEdgeElementsOrderedIds.begin(),
              mTrailingEdgeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("trailing_edge_sub_model_part").AddElements(mTrailingEdgeElementsOrderedIds);
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
void Define2DWakeProcess::MarkKuttaElements() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_sub_model_part");

    for (auto& r_element : trailing_edge_sub_model_part.Elements()) {
        // Compute the distance from the element's center to the trailing edge
        const auto distance_to_element_center =
            ComputeDistanceFromTrailingEdgeToPoint(r_element.GetGeometry().Center());

        // Compute the projection of the distance vector in the wake normal
        // direction
        const double distance_to_wake = inner_prod(distance_to_element_center, mWakeNormal);

        if (distance_to_wake < 0.0) {
            r_element.SetValue(KUTTA, true);
        }
    }
}

// This function finds the trailing edge element that is further downstream and
// marks it as wake trailing edge element. The rest of trailing edge elements
// are unassigned from the wake.
void Define2DWakeProcess::MarkWakeTrailingEdgeElement() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_sub_model_part");

    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_sub_model_part");

    for (auto& r_element : trailing_edge_sub_model_part.Elements()) {
        if(r_element.GetValue(WAKE)){
            // Trailing edge wake element
            if(CheckIfTrailingEdgeElementIsCutByWake(r_element)){
                r_element.Set(STRUCTURE);
                r_element.SetValue(KUTTA, false);
            }
            //Rest of elements touching the trailing edge but not part of the wake
            else{
                r_element.SetValue(WAKE, false);
                wake_sub_model_part.RemoveElement(r_element.Id());
            }
        }
    }
}

// This function checks if the element is cut by the wake
bool Define2DWakeProcess::CheckIfTrailingEdgeElementIsCutByWake(const Element& rElement) const
{
    unsigned int number_of_nodes_with_negative_distance = 0;
    // REMINDER: In 3D the elemental_distances may not be match with the nodal
    // distances if CalculateDistanceToSkinProcess is used.
    const auto nodal_distances_to_wake = rElement.GetValue(WAKE_ELEMENTAL_DISTANCES);

    // Count how many element nodes are above and below the wake
    for (unsigned int i = 0; i < nodal_distances_to_wake.size(); i++) {
        if (nodal_distances_to_wake(i) < 0.0) {
            number_of_nodes_with_negative_distance += 1;
        }
    }

    return number_of_nodes_with_negative_distance == 1;
}

// This function computes the distance from the trailing edge to an input point
const BoundedVector<double, 3> Define2DWakeProcess::ComputeDistanceFromTrailingEdgeToPoint(const Point& rInputPoint) const
{
    BoundedVector<double, 3> distance_to_point = ZeroVector(3);

    distance_to_point(0) = rInputPoint.X() - mpTrailingEdgeNode->X();
    distance_to_point(1) = rInputPoint.Y() - mpTrailingEdgeNode->Y();

    return distance_to_point;
}

} // namespace Kratos.
