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
#include "define_3d_wake_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "processes/calculate_distance_to_skin_process.h"

namespace Kratos {

// Constructor for Define3DWakeProcess Process
Define3DWakeProcess::Define3DWakeProcess(ModelPart& rTrailingEdgeModelPart,
                                         ModelPart& rBodyModelPart,
                                         ModelPart& rStlWakeModelPart,
                                         const double Tolerance,
                                         const Vector& rWakeNormal)
    : Process(),
      mrTrailingEdgeModelPart(rTrailingEdgeModelPart),
      mrBodyModelPart(rBodyModelPart),
      mrStlWakeModelPart(rStlWakeModelPart),
      mTolerance(Tolerance),
      mWakeNormal(rWakeNormal)
{
    KRATOS_ERROR_IF(mWakeNormal.size() != 3)
        << "The mWakeNormal should be a vector with 3 components!"
        << std::endl;
}

void Define3DWakeProcess::ExecuteInitialize()
{
    InitializeTrailingEdgeSubModelpart();

    InitializeWakeSubModelpart();

    SetWakeAndSpanDirections();

    MarkTrailingEdgeNodes();

    ComputeLowerSurfaceNormals();

    MarkWakeElements();

    MarkKuttaElements();

    AddWakeNodes();

    //Print();
}

void Define3DWakeProcess::Print() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    std::ofstream outfile;
    outfile.open("normal_elements_id.txt");
    std::ofstream outfile_wake;
    outfile_wake.open("wake_elements_id.txt");
    unsigned int normal_elements_counter = 0;
    unsigned int wake_elements_counter = 0;
    for (auto& r_element : trailing_edge_sub_model_part.Elements()){
        if(!r_element.GetValue(WAKE)){
            normal_elements_counter += 1;
            std::ofstream outfile;
            outfile.open("normal_elements_id.txt", std::ios_base::app);
            outfile << r_element.Id();
            outfile << "\n";
        }
        else{
            wake_elements_counter += 1;
            std::ofstream outfile_wake;
            outfile_wake.open("wake_elements_id.txt", std::ios_base::app);
            outfile_wake << r_element.Id();
            outfile_wake << "\n";
        }
    }
    KRATOS_WATCH(normal_elements_counter)
    KRATOS_WATCH(wake_elements_counter)
}

void Define3DWakeProcess::AddWakeNodes() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_elements_model_part");

    std::vector<std::size_t> wake_nodes_ordered_ids;
    for (auto& r_element : wake_sub_model_part.Elements()){
        for (unsigned int i = 0; i < r_element.GetGeometry().size(); i++){
            r_element.GetGeometry()[i].SetValue(WAKE, true);
            wake_nodes_ordered_ids.push_back(r_element.GetGeometry()[i].Id());
        }
    }

    std::sort(wake_nodes_ordered_ids.begin(),
              wake_nodes_ordered_ids.end());
    wake_sub_model_part.AddNodes(wake_nodes_ordered_ids);

}

// This function initializes the variables and removes all of the elements of
// the trailing edge submodelpart
void Define3DWakeProcess::InitializeTrailingEdgeSubModelpart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("trailing_edge_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // trailing_edge_sub_model_part
        ModelPart& trailing_edge_sub_model_part =
            root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

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
        root_model_part.CreateSubModelPart("trailing_edge_elements_model_part");
    }
}

// This function initializes the variables and removes all of the elements of
// the wake submodelpart
void Define3DWakeProcess::InitializeWakeSubModelpart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("wake_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // wake_sub_model_part
        ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_elements_model_part");

        for (auto& r_element : wake_sub_model_part.Elements()){
            r_element.SetValue(WAKE, false);
            r_element.SetValue(WAKE_ELEMENTAL_DISTANCES, ZeroVector(3));
            r_element.Set(TO_ERASE, true);
        }
        wake_sub_model_part.RemoveElements(TO_ERASE);
    }
    else{
        // Creating the wake_sub_model_part
        root_model_part.CreateSubModelPart("wake_elements_model_part");
    }
}

void Define3DWakeProcess::SetWakeAndSpanDirections()
{
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
    MathUtils<double>::CrossProduct(mSpanDirection, mWakeNormal, mWakeDirection);
}

void Define3DWakeProcess::MarkTrailingEdgeNodes()
{
    KRATOS_ERROR_IF(mrTrailingEdgeModelPart.NumberOfNodes() == 0) << "There are no nodes in the mrTrailingEdgeModelPart!"<< std::endl;

    double max_span_position = std::numeric_limits<double>::lowest();
    double min_span_position = std::numeric_limits<double>::max();

    auto p_right_wing_tip_node = &*mrTrailingEdgeModelPart.NodesBegin();
    auto p_left_wing_tip_node = &*mrTrailingEdgeModelPart.NodesBegin();

    for (auto& r_node : mrTrailingEdgeModelPart.Nodes()) {
        r_node.SetValue(TRAILING_EDGE, true);
        r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS, 0.0);
        r_node.SetValue(TE_ELEMENT_COUNTER, 0);
        const auto& r_coordinates = r_node.Coordinates();
        const double distance_projection = inner_prod(r_coordinates, mSpanDirection);

        if(distance_projection > max_span_position){
            p_right_wing_tip_node = &r_node;
            max_span_position = distance_projection;
        }
        else if(distance_projection < min_span_position){
            p_left_wing_tip_node = &r_node;
            min_span_position = distance_projection;
        }
    }

    mpRightWingTipNode = p_right_wing_tip_node;
    mpLeftWingTipNode = p_left_wing_tip_node;
}

void Define3DWakeProcess::ComputeLowerSurfaceNormals() const
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrBodyModelPart.Conditions().size()); i++) {
        ModelPart::ConditionIterator it_cond = mrBodyModelPart.ConditionsBegin() + i;

        auto r_geometry = it_cond->GetGeometry();
        auto surface_normal = r_geometry.UnitNormal(0);
        const double projection = inner_prod(surface_normal, mWakeNormal);

        if(projection > 0.0){
            for (unsigned int i = 0; i < it_cond->GetGeometry().size(); i++) {
                r_geometry[i].SetLock();
                r_geometry[i].SetValue(NORMAL, surface_normal);
                r_geometry[i].SetValue(LOWER_SURFACE, true);
                r_geometry[i].UnSetLock();
            }
        }
        else{
            for (unsigned int i = 0; i < it_cond->GetGeometry().size(); i++) {
                r_geometry[i].SetLock();
                r_geometry[i].SetValue(UPPER_SURFACE, true);
                r_geometry[i].UnSetLock();
            }
        }
    }

}

// This function checks which elements are cut by the wake and marks them as
// wake elements
void Define3DWakeProcess::MarkWakeElements()
{
    KRATOS_TRY;
    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements..." << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    std::vector<std::size_t> wake_elements_ordered_ids;

    CalculateDistanceToSkinProcess<3> distance_calculator(root_model_part, mrStlWakeModelPart);
    distance_calculator.Execute();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(root_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = root_model_part.ElementsBegin() + i;

        // Check if the element is touching the trailing edge
        CheckIfTrailingEdgeElement(*it_elem);

        if(it_elem->Is(TO_SPLIT)){
            it_elem->SetValue(WAKE, true);
            #pragma omp critical
            {
                wake_elements_ordered_ids.push_back(it_elem->Id());
            }
            auto wake_elemental_distances = it_elem->GetValue(ELEMENTAL_DISTANCES);
            auto r_geometry = it_elem->GetGeometry();
            for(unsigned int j = 0; j < wake_elemental_distances.size(); j++){
                if(std::abs(wake_elemental_distances[j] < mTolerance)){
                    if(wake_elemental_distances[j] < 0.0){
                        wake_elemental_distances[j] = - mTolerance;
                    }
                    else{
                        wake_elemental_distances[j] = mTolerance;
                    }
                }
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(WAKE_DISTANCE, wake_elemental_distances[j]);
                r_geometry[j].UnSetLock();
            }
            it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
    }
    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeAndWakeElements(wake_elements_ordered_ids);
    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements finished..." << std::endl;
    KRATOS_CATCH("");
}

// This function checks if the element is touching the trailing edge
void Define3DWakeProcess::CheckIfTrailingEdgeElement(Element& rElement)
{
    // Loop over element nodes
    for (unsigned int i = 0; i < rElement.GetGeometry().size(); i++) {
        // Elements touching the trailing edge are trailing edge elements
        const auto& r_node = rElement.GetGeometry()[i];
        if (r_node.GetValue(TRAILING_EDGE)) {
            rElement.SetValue(TRAILING_EDGE, true);
            #pragma omp critical
            {
                mTrailingEdgeElementsOrderedIds.push_back(rElement.Id());
            }
        }
    }
}

// This function adds the trailing edge elements in the
// trailing_edge_sub_model_part
void Define3DWakeProcess::AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds)
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    std::sort(rWakeElementsOrderedIds.begin(),
              rWakeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("wake_elements_model_part").AddElements(rWakeElementsOrderedIds);

    std::sort(mTrailingEdgeElementsOrderedIds.begin(),
              mTrailingEdgeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("trailing_edge_elements_model_part").AddElements(mTrailingEdgeElementsOrderedIds);
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
void Define3DWakeProcess::MarkKuttaElements()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");
    ModelPart& wake_sub_model_part =
        root_model_part.GetSubModelPart("wake_elements_model_part");

    for (int i = 0; i < static_cast<int>(trailing_edge_sub_model_part.Elements().size()); i++) {
        ModelPart::ElementsContainerType::ptr_iterator it_p_elem = trailing_edge_sub_model_part.Elements().ptr_begin() + i;
        auto p_elem = *it_p_elem;
        NodeType::Pointer p_trailing_edge_node;
        unsigned int number_of_te_nodes = 0;
        for (unsigned int j = 0; j < p_elem->GetGeometry().size(); j++){
            const auto& r_node = p_elem->GetGeometry()[j];
            if (r_node.GetValue(TRAILING_EDGE)){
                p_trailing_edge_node = p_elem->GetGeometry()(j);
                number_of_te_nodes += 1;
            }
        }

        KRATOS_ERROR_IF(number_of_te_nodes < 0.5) << "Number of trailing edge nodes must be larger than 0 " << p_elem->Id()
        << " number_of_te_nodes = " << number_of_te_nodes << std::endl;

        const unsigned int number_of_non_te_nodes = 4 - number_of_te_nodes;

        Vector nodal_distances_to_te = ZeroVector(number_of_non_te_nodes);
        ComputeNodalDistancesToWakeAndLowerSurface(*p_elem, p_trailing_edge_node, nodal_distances_to_te);

        unsigned int number_of_nodes_with_negative_distance = 0;
        unsigned int number_of_nodes_with_positive_distance = 0;

        for(unsigned int j = 0; j < nodal_distances_to_te.size(); j++){
            if(nodal_distances_to_te[j] < 0.0){
                number_of_nodes_with_negative_distance += 1;
            }
            else{
                number_of_nodes_with_positive_distance +=1;
            }
        }

        if(number_of_nodes_with_negative_distance > number_of_te_nodes - 1){
            p_elem->Set(STRUCTURE);
            p_elem->SetValue(WAKE, true);
            wake_sub_model_part.AddElement(*it_p_elem);
            BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
            unsigned int counter = 0;
            for(unsigned int j = 0; j < p_elem->GetGeometry().size(); j++){
                auto& r_node = p_elem->GetGeometry()[j];
                if(r_node.GetValue(TRAILING_EDGE)){
                    r_node.SetValue(WAKE_DISTANCE, mTolerance);
                    wake_elemental_distances[j] = mTolerance;
                    auto& r_number_of_neighbour_elements = r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                    r_number_of_neighbour_elements += 1;
                }
                else{
                    r_node.SetValue(WAKE_DISTANCE, nodal_distances_to_te[counter]);
                    wake_elemental_distances[j] = nodal_distances_to_te[counter];
                    counter += 1;
                }
            }
            p_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
        else if(number_of_nodes_with_positive_distance > 0 && number_of_nodes_with_negative_distance > 0 && p_elem->GetValue(WAKE)){
            p_elem->Set(STRUCTURE);
            p_elem->SetValue(WAKE, true);
            wake_sub_model_part.AddElement(*it_p_elem);
            BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
            unsigned int counter = 0;
            for(unsigned int j = 0; j < p_elem->GetGeometry().size(); j++){
                auto& r_node = p_elem->GetGeometry()[j];
                if(r_node.GetValue(TRAILING_EDGE)){
                    r_node.SetValue(WAKE_DISTANCE, mTolerance);
                    wake_elemental_distances[j] = mTolerance;
                    auto& r_number_of_neighbour_elements = r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                    r_number_of_neighbour_elements += 1;
                }
                else{
                    r_node.SetValue(WAKE_DISTANCE, nodal_distances_to_te[counter]);
                    wake_elemental_distances[j] = nodal_distances_to_te[counter];
                    counter += 1;
                }
            }
            p_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
        else{
            p_elem->SetValue(WAKE, false);
            wake_sub_model_part.RemoveElement(*it_p_elem);
        }
    }
}

void Define3DWakeProcess::ComputeNodalDistancesToWakeAndLowerSurface(Element& rElement, NodeType::Pointer pTrailingEdgeNode, Vector& rNodalDistancesToTe)
{
    unsigned int counter = 0;
    for (unsigned int i = 0; i < rElement.GetGeometry().size(); i++){
        const auto& r_node = rElement.GetGeometry()[i];
        if (!r_node.GetValue(TRAILING_EDGE)){
            // Compute the distance vector from the trailing edge to the node
            const array_1d<double,3> distance_vector = r_node.Coordinates() - pTrailingEdgeNode->Coordinates();

            // Compute the distance in the free stream direction
            const double free_stream_direction_distance = inner_prod(distance_vector, mWakeDirection);

            double distance;
            if(r_node.GetValue(LOWER_SURFACE)){
                distance = - mTolerance;
            }
            else if(r_node.GetValue(UPPER_SURFACE)){
                distance = mTolerance;
            }
            else if(free_stream_direction_distance < 0.0){
                distance = inner_prod(distance_vector, pTrailingEdgeNode->GetValue(NORMAL));
            }
            else{
                distance = inner_prod(distance_vector, mWakeNormal);
            }

            KRATOS_ERROR_IF(counter+1 > rNodalDistancesToTe.size()) << "Counter larger than size " << rElement.Id()
                << " counter = " << counter << " rNodalDistancesToTe.size() = " << rNodalDistancesToTe.size() << std::endl;

            if(std::abs(distance) < mTolerance){
                distance = - mTolerance;
            }

            rNodalDistancesToTe[counter] = distance;
            counter += 1;
        }
    }
}

// This function sets the wake direction and normal for further computations
void Define3DWakeProcess::SetWakeDirectionAndNormal()
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
    mWakeNormalOld(0) = -mWakeDirection(1);
    mWakeNormalOld(1) = mWakeDirection(0);
    mWakeNormalOld(2) = 0.0;
}

// This function finds and saves the trailing edge for further computations
void Define3DWakeProcess::SaveTrailingEdgeNode()
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

// This function selects the elements downstream the trailing edge as
// potentially wake elements
const bool Define3DWakeProcess::CheckIfPotentiallyWakeElement(const Element& rElement) const
{
    // Compute the distance from the trailing edge to the element's center
    const auto distance_to_element_center =
        ComputeDistanceFromTrailingEdgeToPoint(rElement.GetGeometry().Center());

    // Compute the projection of the distance in the wake direction
    return inner_prod(distance_to_element_center, mWakeDirection) > 0.0 ? true : false;
}

// This function computes the distance of the element nodes to the wake
const BoundedVector<double, 3> Define3DWakeProcess::ComputeNodalDistancesToWake(const Element& rElement) const
{
    BoundedVector<double, 3> nodal_distances_to_wake = ZeroVector(3);

    for (unsigned int i = 0; i < rElement.GetGeometry().size(); i++) {
        // Compute the distance from the trailing edge to the node
        const auto distance_from_te_to_node =
            ComputeDistanceFromTrailingEdgeToPoint(rElement.GetGeometry()[i]);

        // Compute the projection of the distance vector in the wake normal
        // direction
        double distance_to_wake = inner_prod(distance_from_te_to_node, mWakeNormalOld);

        // Nodes laying on the wake have a positive distance
        if (std::abs(distance_to_wake) < mTolerance) {
            distance_to_wake = mTolerance;
        }
        nodal_distances_to_wake[i] = distance_to_wake;
    }
    return nodal_distances_to_wake;
}

// This function finds the trailing edge element that is further downstream and
// marks it as wake trailing edge element. The rest of trailing edge elements
// are unassigned from the wake.
void Define3DWakeProcess::MarkWakeTrailingEdgeElement() const
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
const bool Define3DWakeProcess::CheckIfTrailingEdgeElementIsCutByWake(const Element& rElement) const
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
const BoundedVector<double, 3> Define3DWakeProcess::ComputeDistanceFromTrailingEdgeToPoint(const Point& rInputPoint) const
{
    BoundedVector<double, 3> distance_to_point = ZeroVector(3);

    distance_to_point(0) = rInputPoint.X() - mpTrailingEdgeNode->X();
    distance_to_point(1) = rInputPoint.Y() - mpTrailingEdgeNode->Y();

    return distance_to_point;
}

} // namespace Kratos.
