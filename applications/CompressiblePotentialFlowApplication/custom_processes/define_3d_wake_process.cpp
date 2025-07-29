//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Marc Nunez
//

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "compressible_potential_flow_application_variables.h"
#include "define_3d_wake_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "utilities/builtin_timer.h"

namespace Kratos {

// Constructor for Define3DWakeProcess Process
Define3DWakeProcess::Define3DWakeProcess(ModelPart& rTrailingEdgeModelPart,
                                         ModelPart& rBodyModelPart,
                                         ModelPart& rStlWakeModelPart,
                                         Parameters ThisParameters)
    : Process(),
      mrTrailingEdgeModelPart(rTrailingEdgeModelPart),
      mrBodyModelPart(rBodyModelPart),
      mrStlWakeModelPart(rStlWakeModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "tolerance"                            : 1e-9,
        "wake_normal"                          : [0.0,0.0,1.0],
        "wake_direction"                       : [1.0,0.0,0.0],
        "switch_wake_normal"                   : false,
        "count_elements_number"                : false,
        "write_elements_ids_to_file"           : false,
        "shed_wake_from_trailing_edge"         : false,
        "shedded_wake_distance"                : 12.5,
        "shedded_wake_element_size"            : 0.2,
        "decrease_wake_width_at_the_wing_tips" : false,
        "echo_level": 1
    })" );
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mTolerance = ThisParameters["tolerance"].GetDouble();
    mWakeNormal = ThisParameters["wake_normal"].GetVector();
    mWakeDirection = ThisParameters["wake_direction"].GetVector();
    mSwitchWakeDirection = ThisParameters["switch_wake_normal"].GetBool();
    mCountElementsNumber = ThisParameters["count_elements_number"].GetBool();
    mWriteElementsIdsToFile = ThisParameters["write_elements_ids_to_file"].GetBool();
    mShedWakeFromTrailingEdge = ThisParameters["shed_wake_from_trailing_edge"].GetBool();
    mSheddedWakeDistance = ThisParameters["shedded_wake_distance"].GetDouble();
    mSheddedWakeElementSize = ThisParameters["shedded_wake_element_size"].GetDouble();
    mDecreaseWakeWidthAtTheWingTips = ThisParameters["decrease_wake_width_at_the_wing_tips"].GetBool();
    mEchoLevel = ThisParameters["echo_level"].GetInt();

    KRATOS_ERROR_IF(mWakeNormal.size() != 3)
        << "The mWakeNormal should be a vector with 3 components!"
        << std::endl;
}

void Define3DWakeProcess::ExecuteInitialize()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    block_for_each(root_model_part.Nodes(), [&](Node& r_nodes)
    {
        r_nodes.SetValue(UPPER_SURFACE, false);
        r_nodes.SetValue(LOWER_SURFACE, false);
        r_nodes.SetValue(TRAILING_EDGE, false);
        r_nodes.SetValue(WAKE_DISTANCE, 0.0);
    });
    auto& r_elements = root_model_part.Elements();
    VariableUtils().SetNonHistoricalVariable(WAKE, 0, r_elements);

    InitializeTrailingEdgeSubModelpart();

    InitializeWakeSubModelpart();

    // Compute span direction as the cross product: mWakeNormal x mWakeDirection
    MathUtils<double>::CrossProduct(mSpanDirection, mWakeNormal, mWakeDirection);

    MarkTrailingEdgeNodesAndFindWingtipNodes();

    ComputeWingLowerSurfaceNormals();

    ComputeAndSaveLocalWakeNormal();

    if(mShedWakeFromTrailingEdge){
        KRATOS_INFO("Define3DWakeProcess") << " Shedding wake from the trailing edge" << std::endl;
        ShedWakeSurfaceFromTheTrailingEdge();
    }

    MarkWakeElements();

    RecomputeNodalDistancesToWakeOrWingLowerSurface();

    MarkKuttaElements();

    SaveLocalWakeNormalInElements();

    AddWakeNodesToWakeModelPart();

    if(mCountElementsNumber){
        CountElementsNumber();
    }

    if(mWriteElementsIdsToFile){
        WriteElementIdsToFile();
    }
}

// This function initializes the variables and removes all of the elements and
// nodes of the trailing edge submodelpart
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
        auto& r_nodes = trailing_edge_sub_model_part.Nodes();
        VariableUtils().SetFlag(TO_ERASE, false, r_nodes);
        trailing_edge_sub_model_part.RemoveElements(TO_ERASE);
        trailing_edge_sub_model_part.RemoveNodes(TO_ERASE);
    }
    else{
        // Creating the trailing_edge_sub_model_part
        root_model_part.CreateSubModelPart("trailing_edge_elements_model_part");
    }
}

// This function initializes the variables and removes all of the elements and
// nodes of the wake submodelpart
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
        auto& r_nodes = wake_sub_model_part.Nodes();
        VariableUtils().SetFlag(TO_ERASE, false, r_nodes);
        wake_sub_model_part.RemoveElements(TO_ERASE);
        wake_sub_model_part.RemoveNodes(TO_ERASE);
    }
    else{
        // Creating the wake_sub_model_part
        root_model_part.CreateSubModelPart("wake_elements_model_part");
    }
}

// This function marks the trailing edge nodes and finds the wingtip nodes
void Define3DWakeProcess::MarkTrailingEdgeNodesAndFindWingtipNodes()
{
    double max_span_position = std::numeric_limits<double>::lowest();
    double min_span_position = std::numeric_limits<double>::max();

    auto p_right_wing_tip_node = &*mrTrailingEdgeModelPart.NodesBegin();
    auto p_left_wing_tip_node = &*mrTrailingEdgeModelPart.NodesBegin();

    for (auto& r_node : mrTrailingEdgeModelPart.Nodes()) {
        r_node.SetValue(TRAILING_EDGE, true);
        const auto& r_coordinates = r_node.Coordinates();
        const double distance_projection = inner_prod(r_coordinates, mSpanDirection);

        // Wingtip nodes have maximum and minimum span positions
        if(distance_projection > max_span_position){
            p_right_wing_tip_node = &r_node;
            max_span_position = distance_projection;
        }
        if(distance_projection < min_span_position){
            p_left_wing_tip_node = &r_node;
            min_span_position = distance_projection;
        }
    }

    // Marking wingtip nodes
    p_right_wing_tip_node->SetValue(WING_TIP, true);
    p_left_wing_tip_node->SetValue(WING_TIP, true);
}

// This function computes the wing lower surface normals and marks the upper
// and lower surfaces. The wing lower surface normals are used later in
// RecomputeComputeNodalDistancesToWakeOrWingLowerSurface inside the
// MarkKuttaElements function to check whether nodes are above or below the wake
// TODO: Think a better way of doing this.
void Define3DWakeProcess::ComputeWingLowerSurfaceNormals() const
{
    // TO DISCUSS: Is it worth it to run these loops in parallel?
    // So far it has not been noted much difference in terms of speed.
    // Also, when ran in parallel, the result can be random.
    // Mark upper surface
    for (auto& r_cond : mrBodyModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        const auto& surface_normal = r_geometry.UnitNormal(0);
        const double projection = inner_prod(surface_normal, mWakeNormal);

        if(!(projection > 0.0)){
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(UPPER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }

    // Mark lower surface
    for (auto& r_cond : mrBodyModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        const auto& surface_normal = r_geometry.UnitNormal(0);
        const double projection = inner_prod(surface_normal, mWakeNormal);

        if(projection > 0.0){
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(NORMAL, surface_normal);
                r_geometry[j].SetValue(LOWER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }
}


// This function computes the local wake normal at each trailing edge node
// by avaraging the local wake normals of the surrounding conditions
void Define3DWakeProcess::ComputeAndSaveLocalWakeNormal() const
{
    // If mrTrailingEdgeModelPart doesn't have conditions use the global wake normal
    if(mrTrailingEdgeModelPart.NumberOfConditions() < 1){
        // Use global wake normal
        auto& r_nodes = mrTrailingEdgeModelPart.Nodes();
        VariableUtils().SetNonHistoricalVariable(WAKE_NORMAL, mWakeNormal, r_nodes);
    }
    else{
        for (auto& r_cond : mrTrailingEdgeModelPart.Conditions()){
            auto& r_geometry = r_cond.GetGeometry();
            const auto& coordinates1 = r_geometry[0].Coordinates();
            const auto& coordinates2 = r_geometry[1].Coordinates();

            const auto& local_span_direction = coordinates2 - coordinates1;
            array_1d<double, 3> local_wake_normal = ZeroVector(3);
            MathUtils<double>::CrossProduct(local_wake_normal, mWakeDirection, local_span_direction);

            const double projection = inner_prod(local_wake_normal, mWakeNormal);

            if(projection < 0.0){
                local_wake_normal *= -1.0;
            }

            for (unsigned int i = 0; i < r_geometry.size(); i++){
                r_geometry[i].GetValue(WAKE_NORMAL) += local_wake_normal;
            }
        }

        for (auto& r_node: mrTrailingEdgeModelPart.Nodes()){
            auto& local_wake_normal = r_node.GetValue(WAKE_NORMAL);
            const double norm = MathUtils<double>::Norm3(local_wake_normal);
            local_wake_normal /= norm;
        }
    }
}

// This function creates the wake surface automatically by shedding it from the
// trailing edge in the direction of the free stream velocity (mWakeDirection).
// The user can decide how much distance is to be shedded and the element size
// of the wake surface in the wake direction. Note that the element size in span
// direction is predetermined by the size of the conditions constituting the
// trailing edge.
void Define3DWakeProcess::ShedWakeSurfaceFromTheTrailingEdge() const
{
    const Properties::Pointer pElemProp = mrStlWakeModelPart.pGetProperties(0);
    const double number_of_elements = mSheddedWakeDistance / mSheddedWakeElementSize;
    const unsigned int number_of_elements_in_wake_direction = int(number_of_elements);

    IndexType node_index = 0;
    IndexType element_index = 0;

    array_1d<double,3> coordinates1;
    array_1d<double,3> coordinates2;
    array_1d<double,3> coordinates3;
    array_1d<double,3> coordinates4;

    for (auto& r_cond : mrTrailingEdgeModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        coordinates1 = r_geometry[0].Coordinates();
        coordinates2 = r_geometry[1].Coordinates();

        if(mDecreaseWakeWidthAtTheWingTips){
            if(r_geometry[0].GetValue(WING_TIP) != 0){
                DecreaseWakeWidthAtTheWingTips(coordinates1, coordinates2);
            }
            else if(r_geometry[1].GetValue(WING_TIP) != 0){
                DecreaseWakeWidthAtTheWingTips(coordinates2, coordinates1);
            }
        }

        coordinates3 = coordinates1 + mSheddedWakeElementSize * mWakeDirection;
        coordinates4 = coordinates2 + mSheddedWakeElementSize * mWakeDirection;

        CreateWakeSurfaceNodesAndElements(node_index, coordinates1, coordinates2, coordinates3,
                                   coordinates4, element_index, pElemProp);

        for (unsigned int j = 0; j < number_of_elements_in_wake_direction; j++){
            coordinates1 = coordinates3;
            coordinates2 = coordinates4;
            coordinates3 = coordinates1 + mSheddedWakeElementSize * mWakeDirection;
            coordinates4 = coordinates2 + mSheddedWakeElementSize * mWakeDirection;

            CreateWakeSurfaceNodesAndElements(node_index, coordinates1, coordinates2, coordinates3,
                                   coordinates4, element_index, pElemProp);
        }
    }
}


// This function decreases the wake width, avoiding cutting elements outside of the wake.
void Define3DWakeProcess::DecreaseWakeWidthAtTheWingTips(array_1d<double, 3>& rPoint1,
                                                         const array_1d<double, 3>& rPoint2) const
{
    const auto& te_edge = rPoint1 - rPoint2;
    const double projection = inner_prod(te_edge, mSpanDirection);

    if (projection > 0.0)
    {
        rPoint1 -= 1e-6 * te_edge;
    }
    else
    {
        rPoint1 += 1e-6 * te_edge;
    }
}

void Define3DWakeProcess::CreateWakeSurfaceNodesAndElements(
    IndexType& rNode_index,
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4,
    IndexType& rElement_index,
    const Properties::Pointer pElemProp) const
{
    const std::array<ModelPart::IndexType, 4> nodes_ids = CreateWakeSurfaceNodes(
        rNode_index, rCoordinates1, rCoordinates2, rCoordinates3, rCoordinates4);

    const double normal_projection = ComputeFaceNormalProjectionToWakeNormal(
        rCoordinates1, rCoordinates2, rCoordinates3, rCoordinates4);

    CreateWakeSurfaceElements(normal_projection, rElement_index, nodes_ids, pElemProp);
}

std::array<ModelPart::IndexType, 4> Define3DWakeProcess::CreateWakeSurfaceNodes(
    IndexType& rNode_index,
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4) const
{
    const auto& p_node1 = mrStlWakeModelPart.CreateNewNode(
        ++rNode_index, rCoordinates1[0], rCoordinates1[1], rCoordinates1[2]);
    const auto& p_node2 = mrStlWakeModelPart.CreateNewNode(
        ++rNode_index, rCoordinates2[0], rCoordinates2[1], rCoordinates2[2]);
    const auto& p_node3 = mrStlWakeModelPart.CreateNewNode(
        ++rNode_index, rCoordinates3[0], rCoordinates3[1], rCoordinates3[2]);
    const auto& p_node4 = mrStlWakeModelPart.CreateNewNode(
        ++rNode_index, rCoordinates4[0], rCoordinates4[1], rCoordinates4[2]);

    return {p_node1->Id(), p_node2->Id(), p_node3->Id(), p_node4->Id()};
}

double Define3DWakeProcess::ComputeFaceNormalProjectionToWakeNormal(
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4) const
{
    const auto& side1 = rCoordinates2 - rCoordinates1;
    const auto& side2 = rCoordinates3 - rCoordinates1;
    array_1d<double, 3> face_normal = ZeroVector(3);
    MathUtils<double>::CrossProduct(face_normal, side1, side2);
    return inner_prod(face_normal, mWakeNormal);
}

void Define3DWakeProcess::CreateWakeSurfaceElements(const double normal_projection,
                                                    IndexType& rElement_index,
                                                    const std::array<ModelPart::IndexType, 4>& rNodes_ids,
                                                    const Properties::Pointer pElemProp) const
{
    if(normal_projection > 0.0){
        const std::vector<ModelPart::IndexType> elem_nodes_1{rNodes_ids[0], rNodes_ids[1], rNodes_ids[2]};
        const std::vector<ModelPart::IndexType> elem_nodes_2{rNodes_ids[1], rNodes_ids[3], rNodes_ids[2]};
        mrStlWakeModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_1, pElemProp);
        mrStlWakeModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_2, pElemProp);
    }
    else{
        const std::vector<ModelPart::IndexType> elem_nodes_1{rNodes_ids[0], rNodes_ids[2], rNodes_ids[1]};
        const std::vector<ModelPart::IndexType> elem_nodes_2{rNodes_ids[1], rNodes_ids[2], rNodes_ids[3]};
        mrStlWakeModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_1, pElemProp);
        mrStlWakeModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_2, pElemProp);
    }
}

// This function checks which elements are cut by the wake and marks them as
// wake elements
void Define3DWakeProcess::MarkWakeElements() const
{
    KRATOS_TRY;
    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements..." << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    BuiltinTimer timer;

    CalculateDiscontinuousDistanceToSkinProcess<3> distance_calculator(root_model_part, mrStlWakeModelPart);
    distance_calculator.Execute();

    KRATOS_INFO_IF("MarkWakeElements", mEchoLevel > 0)
        << "distance_calculator took " << timer.ElapsedSeconds() << " [sec]" << std::endl;

    // This variable allows to inverse the distances computed with the distance
    // process, it is useful if the user makes a mistake and defines the wake
    // surface with the normal pointing downwards.
    double wake_normal_switching_factor = 1.0;
    if(mSwitchWakeDirection){
        KRATOS_INFO("MarkWakeElements") << " Switching wake element distances!" << std::endl;
        wake_normal_switching_factor = -1.0;
    }

    // Variable to store element ids
    moodycamel::ConcurrentQueue<std::size_t> wake_elements_ordered_ids_concurrent_queue;
    moodycamel::ConcurrentQueue<std::size_t> trailing_edge_elements_ordered_ids_concurrent_queue;

    block_for_each(root_model_part.Elements(), [&](Element& rElement)
    {
        // Check if the element is touching the trailing edge
        auto& r_geometry = rElement.GetGeometry();
        CheckIfTrailingEdgeElement(rElement, r_geometry, trailing_edge_elements_ordered_ids_concurrent_queue);

        // Mark wake elements, save their ids, save the elemental distances in
        // the element and in the nodes, and save wake nodes ids
        if (rElement.Is(TO_SPLIT))
        {
            // Mark wake elements
            rElement.SetValue(WAKE, true);

            wake_elements_ordered_ids_concurrent_queue.enqueue(rElement.Id());

            // Save elemental distances in the element
            array_1d<double, 4> wake_elemental_distances = ZeroVector(4);
            wake_elemental_distances = wake_normal_switching_factor *
                                       rElement.GetValue(ELEMENTAL_DISTANCES);
            rElement.SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);

            // Save elemental distances in the nodes
            for (unsigned int j = 0; j < wake_elemental_distances.size(); j++)
            {
                if (std::abs(wake_elemental_distances[j]) < mTolerance)
                {
                    if (wake_elemental_distances[j] < 0.0)
                    {
                        wake_elemental_distances[j] = -mTolerance;
                    }
                    else
                    {
                        wake_elemental_distances[j] = mTolerance;
                    }
                }
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(WAKE_DISTANCE, wake_elemental_distances[j]);
                r_geometry[j].UnSetLock();
            }
        }
    });

    // Variable to store element ids
    std::vector<std::size_t> wake_elements_ordered_ids;
    std::vector<std::size_t> trailing_edge_elements_ordered_ids;

    bool found_value;
    std::size_t id;

    while ( (found_value = wake_elements_ordered_ids_concurrent_queue.try_dequeue(id)) ) {
        wake_elements_ordered_ids.push_back(id);
    }

    while ( (found_value = trailing_edge_elements_ordered_ids_concurrent_queue.try_dequeue(id)) ) {
        trailing_edge_elements_ordered_ids.push_back(id);
    }

    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeAndWakeElements(wake_elements_ordered_ids, trailing_edge_elements_ordered_ids);

    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements finished..." << std::endl;
    KRATOS_CATCH("");
}

// This function checks if the element is touching the trailing edge
void Define3DWakeProcess::CheckIfTrailingEdgeElement(
    Element& rElement,
    const Geometry<NodeType>& rGeometry,
    moodycamel::ConcurrentQueue<std::size_t>& rTrailingEdgeElementsOrderedIds) const
{
    // Loop over element nodes
    for (unsigned int i = 0; i < rGeometry.size(); i++)
    {
        // Elements touching the trailing edge are trailing edge elements
        if (rGeometry[i].GetValue(TRAILING_EDGE))
        {
            rElement.SetValue(TRAILING_EDGE, true);
            rTrailingEdgeElementsOrderedIds.enqueue(rElement.Id());
        }
    }
}

// This function adds the trailing edge elements in the
// trailing_edge_sub_model_part
void Define3DWakeProcess::AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds,
                                                         std::vector<std::size_t>& rTrailingEdgeElementsOrderedIds) const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    std::sort(rWakeElementsOrderedIds.begin(),
              rWakeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("wake_elements_model_part").AddElements(rWakeElementsOrderedIds);

    std::sort(rTrailingEdgeElementsOrderedIds.begin(),
              rTrailingEdgeElementsOrderedIds.end());
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");
    trailing_edge_sub_model_part.AddElements(rTrailingEdgeElementsOrderedIds);

    std::vector<std::size_t> trailing_edge_nodes_ordered_ids;
    trailing_edge_nodes_ordered_ids.reserve(trailing_edge_sub_model_part.NumberOfElements());

    // Add also the nodes of the elements touching the trailing edge
    for (auto& r_elem : trailing_edge_sub_model_part.Elements())
    {
        for (auto& r_node : r_elem.GetGeometry()){
            trailing_edge_nodes_ordered_ids.push_back(r_node.Id());
        }
    }

    std::sort(trailing_edge_nodes_ordered_ids.begin(),
              trailing_edge_nodes_ordered_ids.end());
    trailing_edge_sub_model_part.AddNodes(trailing_edge_nodes_ordered_ids);
}

// This function recomputes the wake distances from the nodes belonging to the
// elements. These distances are used later to decide which elements are KUTTA,
// which are WAKE (STRUCTURE), and which are NORMAL
void Define3DWakeProcess::RecomputeNodalDistancesToWakeOrWingLowerSurface() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    block_for_each(trailing_edge_sub_model_part.Nodes(), [&](Node& r_node)
    {
        // Trailing edge nodes are assigned a positive distance
        if(r_node.GetValue(TRAILING_EDGE)){
            r_node.SetValue(WAKE_DISTANCE, mTolerance);
        }
        // Nodes in the lower surface are assigned a negative distance
        else if(r_node.GetValue(LOWER_SURFACE)){
            r_node.SetValue(WAKE_DISTANCE, -mTolerance);
        }
        // Nodes in the upper surface are assigned a positive distance
        else if(r_node.GetValue(UPPER_SURFACE)){
            r_node.SetValue(WAKE_DISTANCE, mTolerance);
        }
        // For the rest of the nodes the distance is recomputed:
        else{
            // TODO: improve the distance calculation.
            // Find closest trailing edge node
            NodeType::Pointer p_closest_te_node = *mrTrailingEdgeModelPart.NodesBegin().base();
            FindClosestTrailingEdgeNode(p_closest_te_node, r_node);
            RecomputeDistance(p_closest_te_node, r_node);
        }
    });
}

// This function finds the closest trailing edge node to the given point
void Define3DWakeProcess::FindClosestTrailingEdgeNode(NodeType::Pointer& pClosest_te_node,
                                                      const array_1d<double, 3>& rPoint) const
{
    double min_distance_to_te = std::numeric_limits<double>::max();
    for (auto& r_te_node : mrTrailingEdgeModelPart.Nodes())
    {
        const auto& distance_vector = rPoint - r_te_node;
        const double distance_to_te = inner_prod(distance_vector, distance_vector);
        if (distance_to_te < min_distance_to_te)
        {
            min_distance_to_te = distance_to_te;
            pClosest_te_node = &r_te_node;
        }
    }
}

// This function recomputes the distances from the given node to the wake or to
// the wing lower surface, depending whether the node is behind or in front of
// the trailing edge according to the wake direction (free stream velocity
// direction).
void Define3DWakeProcess::RecomputeDistance(NodeType::Pointer& pClosest_te_node,
                                            NodeType& rNode) const
{
    // Compute the distance vector from the closest trailing edge node to the
    // current node
    const auto& distance_vector = rNode - *pClosest_te_node;

    // Compute the distance in the free stream direction
    const double free_stream_direction_distance = inner_prod(distance_vector, mWakeDirection);

    //  For the nodes in front of the selected trailing edge node, the
    //  distance is computed according to the lower surface normal
    if(free_stream_direction_distance < 0.0){
        double distance = inner_prod(distance_vector, pClosest_te_node->GetValue(NORMAL));
        // Nodes close to the wing are given a negative distance
        if(std::abs(distance) < mTolerance){
            distance = - mTolerance;
        }
        rNode.SetValue(WAKE_DISTANCE, distance);
    }
    // For the nodes behind of the selected trailing edge node, the
    // distance is computed according to the wake normal
    else{
        const auto& local_wake_normal = pClosest_te_node->GetValue(WAKE_NORMAL);
        double distance = inner_prod(distance_vector, local_wake_normal);
        // Nodes slightly below and above the wake are given a positive distance (wake down)
        if(std::abs(distance) < mTolerance){
            distance = mTolerance;
        }
        rNode.SetValue(WAKE_DISTANCE, distance);
    }
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
void Define3DWakeProcess::MarkKuttaElements() const
{
    KRATOS_INFO("MarkKuttaElements") << "...Selecting kutta elements..." << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    std::vector<std::size_t> wake_elements_ordered_ids;

    // TO DISCUSS: Is it worth it to run this loop in parallel?
    // So far not much different in terms of speed has been noted.
    block_for_each(trailing_edge_sub_model_part.Elements(), [&](Element& rElement)
    {
        auto& r_geometry = rElement.GetGeometry();
        // Counting number of trailing edge nodes in the element
        unsigned int number_of_te_nodes = CountNumberOfTrailindEdgeNodesInElement(r_geometry);

        KRATOS_ERROR_IF(number_of_te_nodes < 1)
            << "Number of trailing edge nodes must be 1 or larger. Element Id: "
            << rElement.Id() << " number_of_te_nodes = " << number_of_te_nodes
            << std::endl;

        // Checking only distances from nodes that are not trailing edge
        unsigned int number_of_nodes_with_negative_distance = 0;
        unsigned int number_of_nodes_with_positive_distance = 0;

        CountNumberOfPositiveAndNegativeDistances(
            r_geometry, number_of_nodes_with_negative_distance,
            number_of_nodes_with_positive_distance);

        SelectElementType(rElement, r_geometry, number_of_te_nodes,
                          number_of_nodes_with_negative_distance,
                          number_of_nodes_with_positive_distance);
    });

    // Remove elements that were wake and now are either kutta or normal.
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");
    wake_sub_model_part.RemoveElements(TO_ERASE);
    KRATOS_INFO("MarkKuttaElements") << "...Selecting kutta elements finished..." << std::endl;
}

// This function returns the number of trailing edge nodes
unsigned int Define3DWakeProcess::CountNumberOfTrailindEdgeNodesInElement(
    const Geometry<NodeType>& rGeometry) const
{
    unsigned int number_of_te_nodes = 0;

    for (unsigned int j = 0; j < rGeometry.size(); j++)
    {
        const auto& r_node = rGeometry[j];
        if (r_node.GetValue(TRAILING_EDGE))
        {
            number_of_te_nodes += 1;
        }
    }

    return number_of_te_nodes;
}

// This function returns the number of non trailing edge nodes with positive and
// negative distance
void Define3DWakeProcess::CountNumberOfPositiveAndNegativeDistances(
    const Geometry<NodeType>& rGeometry,
    unsigned int& number_of_nodes_with_negative_distance,
    unsigned int& number_of_nodes_with_positive_distance) const
{
    for (unsigned int j = 0; j < rGeometry.size(); j++){
        const auto& r_node = rGeometry[j];
        if (!r_node.GetValue(TRAILING_EDGE)){
            const auto& distance = r_node.GetValue(WAKE_DISTANCE);
            if(distance < 0.0){
                number_of_nodes_with_negative_distance += 1;
            }
            else{
                number_of_nodes_with_positive_distance +=1;
            }
        }
    }
}

// This function selects the element type (wake, kutta, or normal) according to
// the number of non trailing edge nodes with positive and negative distance
void Define3DWakeProcess::SelectElementType(Element& rElement,
                                            const Geometry<NodeType>& rGeometry,
                                            const unsigned int number_of_te_nodes,
                                            const unsigned int number_of_nodes_with_negative_distance,
                                            const unsigned int number_of_nodes_with_positive_distance) const
{
    // Wake structure elements (cut)
    if (number_of_nodes_with_positive_distance > 0 &&
        number_of_nodes_with_negative_distance > 0 && rElement.GetValue(WAKE))
    {
        rElement.Set(STRUCTURE);
        BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
        for (unsigned int j = 0; j < rGeometry.size(); j++)
        {
            const auto& r_node = rGeometry[j];
            const auto& distance = r_node.GetValue(WAKE_DISTANCE);
            wake_elemental_distances[j] = distance;
        }
        rElement.SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
    }

    // Kutta elements (below). Kutta elements have all non_te_nodes with
    // negative distance:
    // 1 te node  -> 3 non te nodes with negative distance > 2 = 3 - 1
    // 2 te nodes -> 2 non te nodes with negative distance > 1 = 3 - 2
    else if (number_of_nodes_with_negative_distance > 3 - number_of_te_nodes)
    {
        rElement.SetValue(KUTTA, true);
        rElement.SetValue(WAKE, false);
        rElement.Set(TO_ERASE, true);
    }

    // Normal elements (above). Normal elements have all nodes with positive
    // distance.
    else
    {
        rElement.SetValue(WAKE, false);
        rElement.Set(TO_ERASE, true);
    }
}

// This function saves the local wake normal in the element to be used later
// inside the element to apply the wake conditions
void Define3DWakeProcess::SaveLocalWakeNormalInElements() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");

    for (auto& r_elem : wake_sub_model_part.Elements()){
        NodeType::Pointer p_closest_te_node = *mrTrailingEdgeModelPart.NodesBegin().base();
        FindClosestTrailingEdgeNode(p_closest_te_node, r_elem.GetGeometry().Center());

        r_elem.GetValue(WAKE_NORMAL) = p_closest_te_node->GetValue(WAKE_NORMAL);
    }
}

void Define3DWakeProcess::AddWakeNodesToWakeModelPart() const
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

// This function counts the number of elements of each type. Useful for
// debugging purposes.
void Define3DWakeProcess::CountElementsNumber() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    // Initialize counters
    unsigned int normal_elements_counter = 0;
    unsigned int wake_elements_counter = 0;
    unsigned int kutta_elements_counter = 0;
    unsigned int structure_elements_counter = 0;
    for (auto& r_element : trailing_edge_sub_model_part.Elements()){
        if(!r_element.GetValue(WAKE)){
            if(r_element.GetValue(KUTTA)){
                kutta_elements_counter += 1;
            }
            else{
                normal_elements_counter += 1;
            }
        }
        else{
            wake_elements_counter += 1;
            if(r_element.Is(STRUCTURE)){
                structure_elements_counter += 1;
            }
        }
    }

    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");
    const unsigned int all_wake_elements_counter = wake_sub_model_part.NumberOfElements();
    KRATOS_WATCH(normal_elements_counter)
    KRATOS_WATCH(kutta_elements_counter)
    KRATOS_WATCH(wake_elements_counter)
    KRATOS_WATCH(structure_elements_counter)
    KRATOS_WATCH(all_wake_elements_counter)
}

// This function prints the elements ids into a file. Useful for debugging
// purposes.
void Define3DWakeProcess::WriteElementIdsToFile() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    std::ofstream outfile;
    outfile.open("normal_elements_id.txt");
    std::ofstream outfile_wake;
    outfile_wake.open("wake_elements_id.txt");
    std::ofstream outfile_structure;
    outfile_structure.open("structure_elements_id.txt");
    std::ofstream outfile_kutta;
    outfile_kutta.open("kutta_elements_id.txt");
        for (auto& r_element : trailing_edge_sub_model_part.Elements()){
        if(!r_element.GetValue(WAKE)){
            if(r_element.GetValue(KUTTA)){
                outfile_kutta << r_element.Id();
                outfile_kutta << "\n";
            }
            else{
                outfile << r_element.Id();
                outfile << "\n";
            }
        }
        else{
            outfile_wake << r_element.Id();
            outfile_wake << "\n";
            if(r_element.Is(STRUCTURE)){
                outfile_structure << r_element.Id();
                outfile_structure << "\n";
            }
        }
    }
    outfile_kutta.close();
    outfile.close();
    outfile_structure.close();
    outfile_wake.close();

    // Loop over all wake elements
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");

    std::ofstream outfile_all_wake;
    outfile_all_wake.open("all_wake_elements_id.txt");
    for (auto& r_element : wake_sub_model_part.Elements()){
        outfile_all_wake << r_element.Id();
        outfile_all_wake << "\n";
    }
    outfile_all_wake.close();
}
} // namespace Kratos.
