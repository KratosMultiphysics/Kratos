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
#include "custom_processes/move_model_part_process.h"

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

bool DefineEmbeddedWakeProcess::CheckIfValid(NodeType* p_trailing_edge_candidate) {
    auto neighbour_elements_list = p_trailing_edge_candidate->GetValue(NEIGHBOUR_ELEMENTS);
    for (auto neighbour_element : neighbour_elements_list) {
        if (neighbour_element.IsNot(ACTIVE)) {
            return true;
        }
        auto r_neighbour_geometry = neighbour_element.GetGeometry();
        for (unsigned int i_node = 0; i_node < r_neighbour_geometry.size(); i_node++) {
            auto neighbour_distance_value = r_neighbour_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
            if (neighbour_distance_value < 0.0 && r_neighbour_geometry[i_node].X() < p_trailing_edge_candidate->X()) {
                if (CheckIfValid(&r_neighbour_geometry[i_node])) {
                    return true;
                }
            }
        }
    }

    return false;
}

Node<3>* DefineEmbeddedWakeProcess::FindNode(double restarted_search_maximum) {
    double max_x_coordinate = -1e30;
    NodeType* p_trailing_edge_node;
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        if (it_elem->Is(TO_SPLIT) && it_elem->Is(ACTIVE)) {
            auto& r_geometry = it_elem->GetGeometry();
            for (unsigned int i_node = 0; i_node < r_geometry.size(); i_node++) {
                auto x_value = r_geometry[i_node].X();
                auto distance_value = r_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
                if (x_value>max_x_coordinate && x_value < restarted_search_maximum && distance_value < 0.0) {
                    max_x_coordinate = x_value;
                    p_trailing_edge_node = &r_geometry[i_node];
                }

            }
        }
    }
    return p_trailing_edge_node;
}

void DefineEmbeddedWakeProcess::Execute()
{
    KRATOS_TRY;
    double restarted_search_maximum = 1e30;
    NodeType* trailing_edge_candidate = FindNode(restarted_search_maximum);
    bool is_valid = CheckIfValid(trailing_edge_candidate);

    while (!is_valid){
        KRATOS_INFO("FAILED_NODE") << trailing_edge_candidate->Id() << std::endl;
        trailing_edge_candidate->FastGetSolutionStepValue(GEOMETRY_DISTANCE) = 1e-9;
        trailing_edge_candidate = FindNode(trailing_edge_candidate->X());
        is_valid = CheckIfValid(trailing_edge_candidate);
        if (!is_valid) {
            KRATOS_INFO("FAILED_NODE") << trailing_edge_candidate->Id() << std::endl;
            trailing_edge_candidate->FastGetSolutionStepValue(GEOMETRY_DISTANCE) = 1e-9;
        }
    }


    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->SetValue(WING_TIP, false);
        it_node->SetValue(TRAILING_EDGE, false);
        it_node->SetValue(AIRFOIL, false);
        it_node->SetValue(KUTTA, false);
        it_node->SetValue(LOWER_SURFACE, false);
        it_node->SetValue(UPPER_SURFACE, false);
    }

    double max_value = -1e30;
    for(unsigned int i = 0; i<3; i++){
        mWakeOrigin[i]=0.0;
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        if (it_elem->Is(TO_SPLIT)) {
            it_elem->Set(BOUNDARY, true);
        }
        if (it_elem->IsNot(ACTIVE)) {
            auto element_center = it_elem ->GetGeometry().Center();
            if (element_center[0] > max_value) {
                max_value = element_center[0];
                mWakeOrigin = element_center;
            }

        }
        auto geometry_elemental_distances = it_elem->GetValue(ELEMENTAL_DISTANCES);
        it_elem->SetValue(GEOMETRY_ELEMENTAL_DISTANCES, geometry_elemental_distances);
    }

    // Parameters moving_parameters;
    // // // mWakeOrigin[0] += -0.001;
    // moving_parameters.AddEmptyValue("origin");
    // moving_parameters["origin"].SetVector(mWakeOrigin);
    // moving_parameters.AddEmptyValue("rotation_angle");
    // double angle=-mrModelPart.GetProcessInfo()[ROTATION_ANGLE]*Globals::Pi/180;
    // mrModelPart.GetProcessInfo()[WAKE_ORIGIN] = mWakeOrigin;
    // moving_parameters["rotation_angle"].SetDouble(angle);
    // // moving_parameters["rotation_angle"].SetDouble(0.0);
    // MoveModelPartProcess move_process(mrWakeModelPart, moving_parameters);
    // move_process.Execute();

    KRATOS_ERROR_IF(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]>2) << "DOMAIN_SIZE is greater than 2. DefineEmbeddedWakeProcess is only implemented for 2D cases!" << std::endl;

    ExecuteInitialize();
    ComputeDistanceToWake();
    MarkWakeElements();
    ComputeTrailingEdgeNode();
    // MarkKuttaWakeElements();



    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        it_elem->Set(TO_SPLIT, false);
        if (it_elem->Is(BOUNDARY)) {
            it_elem->Set(TO_SPLIT, true);
        }
        auto geometry_elemental_distances = it_elem->GetValue(GEOMETRY_ELEMENTAL_DISTANCES);
        it_elem->SetValue(ELEMENTAL_DISTANCES, geometry_elemental_distances);
    }

    RedefineWake();

    std::unordered_set<std::size_t> touching_elements;

    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");
    for (auto& r_elem : deactivated_model_part.Elements()) {
        if (r_elem.Is(STRUCTURE) && touching_elements.find(r_elem.Id()) == touching_elements.end()) {
            // KRATOS_WATCH("***********************************")

            std::unordered_set<std::size_t> visited_elements;
            bool touches_wake = TouchesWake(r_elem, visited_elements);
            // KRATOS_WATCH("***********************************")
            // KRATOS_WATCH(touches_wake)
            // KRATOS_WATCH("***********************************")
            if (!touches_wake) {
                r_elem.Set(MARKER, true);
                r_elem.SetValue(WAKE, false);
                r_elem.Set(STRUCTURE, false);
                // for (auto& r_node : r_elem.GetGeometry()) {
                //     r_node.SetValue(WING_TIP, false);
                // }
            } else {
                for (auto elem_id : visited_elements) {
                    touching_elements.insert(elem_id);
                }
                for (auto& r_node : r_elem.GetGeometry()) {
                    r_node.SetValue(WAKE, true);
                }
            }
        }
    }

    for (auto& r_elem : deactivated_model_part.Elements()) {
        if (r_elem.Is(MARKER)) {
            for (auto& r_node : r_elem.GetGeometry()) {
                // if(r_node.GetValue(UPPER_SURFACE) && r_node.GetValue(LOWER_SURFACE) && !r_node.GetValue(WAKE)){
                if(r_node.GetValue(WAKE_DISTANCE) < 0.0) {
                // if(r_node.GetValue(UPPER_SURFACE) && r_node.GetValue(LOWER_SURFACE)){
                    r_node.SetValue(TRAILING_EDGE, true);
                    r_elem.SetValue(KUTTA, true);
                }
            }
        }
    }

    // for (auto& r_elem : mrModelPart.Elements()) {
    //     std::size_t upper_nodes = 0;
    //     bool touches_te = false;
    //     for (auto& r_node : r_elem.GetGeometry()) {
    //         if(r_node.GetValue(UPPER_SURFACE)){
    //             upper_nodes++;
    //         }
    //         if(r_node.GetValue(TRAILING_EDGE))
    //             touches_te = true;
    //     }
    //     bool is_upper = upper_nodes == 3;
    //     if (is_upper && touches_te) {
    //         r_elem.SetValue(KUTTA, true);
    //     }
    //     // if (touches_te && r_elem.Is(ACTIVE)) {
    //     //     for (auto& r_node : r_elem.GetGeometry()) {
    //     //         r_node.SetValue(WING_TIP, true);
    //     //     }
    //     // }
    // }

    Element::Pointer airfoil_pointer = mrModelPart.pGetElement(1);
    Element::Pointer kutta_wake_pointer = mrModelPart.pGetElement(1);
    for (auto& r_elem : mrModelPart.Elements()) {
        std::size_t counter = 0;
        for (auto& r_node : r_elem.GetGeometry()) {
            bool is_airfoil = r_node.GetValue(AIRFOIL);
            bool is_kutta_wake = r_node.GetValue(KUTTA);
            if (is_airfoil && is_kutta_wake)
                counter++;
        }
        if (counter==2 &&  r_elem.IsNot(ACTIVE)) {
            airfoil_pointer = &r_elem;
        }
        else if (counter==2 && r_elem.Is(STRUCTURE)) {
            kutta_wake_pointer = &r_elem;
        } else if (counter>0) {
            for (auto& r_node : r_elem.GetGeometry()) {
                r_node.SetValue(WING_TIP, true);
            }
        }

    }
    std::size_t x_nodes = 0;
    //CHECKING airfoil element
    double airfoil_center_x = airfoil_pointer->GetGeometry().Center().X();
    for (auto& r_node : airfoil_pointer->GetGeometry()) {
        if (r_node.X() > airfoil_center_x)
            x_nodes++;
    }
    // if (x_nodes == 1) {
    if (true) {
        double max_x = -1e10;
        Node<3>::Pointer p_node = nullptr;
        for (auto& r_node : kutta_wake_pointer->GetGeometry()) {
            if (r_node.GetValue(AIRFOIL) && r_node.GetValue(KUTTA) && r_node.X() > max_x) {
                max_x = r_node.X();
                p_node = &r_node;
            }
        }

        for (auto& r_neigh_elem : p_node->GetValue(NEIGHBOUR_ELEMENTS)) {
            for (auto& r_node : r_neigh_elem.GetGeometry()) {
                r_node.SetValue(WING_TIP, true);
            }
        }

    } else {
        KRATOS_WATCH("*************")
        KRATOS_WATCH("DOING NOTHING")
        KRATOS_WATCH("*************")
    }


    // std::vector<Element::Pointer> elem_vector;
    // for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
    //     ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i;
    //     auto& r_geometry = it_elem->GetGeometry();


    //     BoundedVector<double,3> geometry_distances;
    //     for(unsigned int i_node = 0; i_node<3; i_node++){
    //         geometry_distances[i_node] = r_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
    //     }
    //     const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);
    //     bool is_touching_marker = false;
    //     // if (is_embedded){
    //     if (true){
    //         for (auto& r_node : r_geometry) {
    //             if (r_node.GetValue(WING_TIP) == 1) {
    //                 // std::cout << "Found embedded elem " << it_elem->Id() << " touching marker" << std::endl;
    //                 is_touching_marker = true;
    //                 // break;
    //             }
    //         }
    //     }

    //     if (is_touching_marker) {
    //         elem_vector.push_back(&(*it_elem));
    //     }
    // }
    // for (auto p_elem : elem_vector) {
    //     for (auto& r_node : p_elem->GetGeometry()) {
    //         r_node.SetValue(WING_TIP, true);
    //     }
    // }



    mrModelPart.RemoveSubModelPart("deactivated_model_part");
    mrModelPart.RemoveSubModelPart("upper_surface_sub_model_part");
    mrModelPart.RemoveSubModelPart("lower_surface_sub_model_part");


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

    ModelPart& deactivated_model_part = mrModelPart.CreateSubModelPart("deactivated_model_part");
    std::vector<std::size_t> deactivated_elements_id_list;
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


        if (is_embedded){
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), Vector(geometry_distances));
            // Computing Normal
            std::vector<array_1d<double,3>> cut_normal;
            pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
            double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
            array_1d<double,3> unit_normal = cut_normal[0]/norm_normal;
            it_elem->SetValue(VELOCITY_LOWER,unit_normal);
        }
        if (it_elem->IsNot(ACTIVE)) {
            for(unsigned int i_node = 0; i_node<3; i_node++){
                it_elem->GetGeometry()[i_node].SetValue(AIRFOIL, true);
            }
        }
        // if (is_wake_element && it_elem->IsNot(ACTIVE)) {
        //     auto& r_geometry = it_elem->GetGeometry();

        //     for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
        //         r_geometry[i].SetLock();
        //         r_geometry[i].SetValue(WING_TIP, true);
        //         r_geometry[i].UnSetLock();
        //     }
        // }
        // Mark wake element and save their nodal distances to the wake

        if (is_wake_element && it_elem->Is(ACTIVE)) {
            it_elem->SetValue(WAKE, true);

            if (is_embedded){
                #pragma omp critical
                {
                    deactivated_elements_id_list.push_back(it_elem->Id());
                }
                // it_elem->Set(ACTIVE, false);
                // it_elem->SetValue(WAKE, true);
                it_elem->Set(STRUCTURE, true);
                auto& r_geometry = it_elem->GetGeometry();
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WING_TIP, true);
                    r_geometry[i].SetValue(KUTTA, true);
                    r_geometry[i].Set(MARKER, true);
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].UnSetLock();
                }
            }
            else{
                // it_elem->SetValue(WAKE, true);
                auto& r_geometry = it_elem->GetGeometry();
                for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                    r_geometry[i].SetLock();
                    r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                    r_geometry[i].SetValue(WAKE, true);
                    r_geometry[i].UnSetLock();
                }
            }
        }
    }
    deactivated_model_part.AddElements(deactivated_elements_id_list);
}



void DefineEmbeddedWakeProcess::ComputeTrailingEdgeNode(){

    // double max_distance = 0.0;
    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");
    Element::Pointer p_max_elem;

    auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];

    // Find furthest deactivated element to the wake origin
    // Find deactivated element with Dim nodes shared with a wake element
    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;
        std::size_t wake_nodes_counter = 0;
        auto& r_geometry = it_elem->GetGeometry();
        for (unsigned int i_node= 0; i_node < it_elem->GetGeometry().size(); i_node++) {
            if (r_geometry[i_node].GetValue(WAKE)) {
                wake_nodes_counter++;
            }
        }
        // std::cout << it_elem->Id() << " " << wake_nodes_counter << std::endl;
        if (wake_nodes_counter == 2) {
            p_max_elem = mrModelPart.pGetElement(it_elem->Id());
        }
    }
    auto angle_in_deg = -1*mrModelPart.GetProcessInfo()[ROTATION_ANGLE];
    BoundedVector<double, 3> wake_direction;
    wake_direction[0] = cos(angle_in_deg*Globals::Pi/180);
    wake_direction[1] = sin(angle_in_deg*Globals::Pi/180);
    wake_direction[2] = 0.0;

    // Mark nodes of the furthest deactivated element and store its neighbour elements
    std::size_t max_node_i = -1;
    double max_x = -1e10;
    auto& r_max_element_geometry=p_max_elem->GetGeometry();


    BoundedVector<double, 3> base_vector;
    base_vector[0] = 1000.0;
    base_vector[1] = 0.0;


    for (unsigned int i_node= 0; i_node < r_max_element_geometry.size(); i_node++) {
        // r_max_element_geometry[i_node].SetValue(AIRFOIL,true);
        // BoundedVector<double, 3> nodal_position_vector = r_max_element_geometry[i_node].Coordinates() - wake_origin;
        BoundedVector<double, 3> nodal_position_vector = base_vector+r_max_element_geometry[i_node].Coordinates();
        double projection = std::abs(inner_prod(nodal_position_vector, wake_direction));
        // KRATOS_WATCH(projection);
        // KRATOS_WATCH(r_max_element_geometry[i_node].Id());
        // if (r_max_element_geometry[i_node].X()>max_x) {
        // if (projection>max_x) {
        if (projection>max_x && r_max_element_geometry[i_node].GetValue(WAKE_DISTANCE)  > 0.0) {
        // if (projection>max_x) {
            max_x = projection;
            max_node_i = i_node;
        }

        const GlobalPointersVector<Element>& r_node_elem_candidates = r_max_element_geometry[i_node].GetValue(NEIGHBOUR_ELEMENTS);
        for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
            mKuttaWakeElementCandidates.push_back(r_node_elem_candidates(j));
        }
    }

    // for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
    //     ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;
    //     // it_elem->SetValue(KUTTA, true);
    //     // it_elem->Set(MARKER, true);
    //     it_elem->Set(STRUCTURE, true);
    //     // it_elem->SetValue(WAKE, true);
    //     auto& r_geometry = it_elem->GetGeometry();
    //     for (unsigned int i_node= 0; i_node < it_elem->GetGeometry().size(); i_node++) {
    //         r_geometry[i_node].Set(MARKER, true);
    //         // std::cout << r_geometry[i_node].Id() << std::endl;
    //         // if (r_geometry[i_node].GetValue(WAKE_DISTANCE) > 0.0) {
    //         //     // r_geometry[i_node].SetValue(TRAILING_EDGE, true);
    //         // }
    //     }
    // }



    KRATOS_ERROR_IF(max_node_i < 0) << "No trailing edge nodes were found" << std::endl;
    // KRATOS_WATCH(r_max_element_geometry[max_node_i].Id());

    // r_max_element_geometry[max_node_i].SetValue(AIRFOIL, true);
    // r_max_element_geometry[max_node_i].SetValue(TRAILING_EDGE, true);
    mrTrailingEdgeNode = r_max_element_geometry[max_node_i];

    std::vector<std::size_t> trailing_edge_node_list;
    trailing_edge_node_list.push_back(r_max_element_geometry[max_node_i].Id());

    if (mrModelPart.HasSubModelPart("trailing_edge_sub_model_part")){
        mrModelPart.RemoveSubModelPart("trailing_edge_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("trailing_edge_sub_model_part");

    std::sort(trailing_edge_node_list.begin(),
              trailing_edge_node_list.end());
    mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").AddNodes(trailing_edge_node_list);
}

bool DefineEmbeddedWakeProcess::TouchesWake(Element& rElem, std::unordered_set<std::size_t> visited_elements) {

    // KRATOS_WATCH(rElem.Id())
    auto& r_geometry = rElem.GetGeometry();
    visited_elements.insert(rElem.Id());
    std::unordered_set<std::size_t> failed_elements;
    for (auto& r_node : r_geometry) {
        for (auto& r_elem_neigh : r_node.GetValue(NEIGHBOUR_ELEMENTS)) {
            if (failed_elements.find(r_elem_neigh.Id()) == failed_elements.end()) {
                bool is_struct = r_elem_neigh.Is(STRUCTURE);
                bool is_wake = r_elem_neigh.GetValue(WAKE);
                std::size_t matching_nodes = 0;
                for (auto& r_node_neigh : r_elem_neigh.GetGeometry()) {
                    for (auto& r_node_main : rElem.GetGeometry()) {
                        if (r_node_main.Id() == r_node_neigh.Id()) {
                            matching_nodes++;
                        }
                    }

                }
                bool is_matching = matching_nodes == 2;
                bool is_visited = visited_elements.find(r_elem_neigh.Id()) != visited_elements.end();
                // std::cout<<rElem.Id() << " NEIGH id, struct,wake "<< r_elem_neigh.Id() << " " << is_struct << " " << is_wake << " "<< is_visited<< " "<< is_matching <<  std::endl;
                // std::cout<<"visited elements: ";
                // for (auto& value : visited_elements) {
                //     std::cout << value << " ";
                // }
                // std::cout << std::endl;
                if (is_matching && is_wake && !is_struct) {
                    return true;
                }
                if (is_matching && is_wake && is_struct && !is_visited) {
                    if (TouchesWake(r_elem_neigh, visited_elements)) {
                        return true;
                    } else {
                        failed_elements.insert(r_elem_neigh.Id());
                    }
                }
            }
        }
    }
    return false;

}

void DefineEmbeddedWakeProcess::RedefineWake(){


    if (mrModelPart.HasSubModelPart("lower_surface_sub_model_part")){
        mrModelPart.RemoveSubModelPart("lower_surface_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("lower_surface_sub_model_part");

    if (mrModelPart.HasSubModelPart("upper_surface_sub_model_part")){
        mrModelPart.RemoveSubModelPart("upper_surface_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("upper_surface_sub_model_part");

    std::vector<std::size_t> lower_surface_elements_list;
    std::vector<std::size_t> upper_surface_elements_list;

    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        BoundedVector<double, 3> geometry_elemental_distances = it_elem->GetValue(GEOMETRY_ELEMENTAL_DISTANCES);
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_elemental_distances);

        auto angle_in_deg = -1*mrModelPart.GetProcessInfo()[ROTATION_ANGLE];
        auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];
        BoundedVector<double, 3> wake_direction;
        // wake_direction[0] = cos(angle_in_deg*Globals::Pi/180);
        // wake_direction[1] = sin(angle_in_deg*Globals::Pi/180);
        // wake_direction[2] = 0.0;
        wake_direction[0] = sin(angle_in_deg*Globals::Pi/180);
        wake_direction[1] = -cos(angle_in_deg*Globals::Pi/180);
        wake_direction[2] = 0.0;
        auto& r_geometry = it_elem->GetGeometry();
        if (is_embedded && it_elem->Is(ACTIVE)) {

            BoundedVector<double, 3> unit_normal = it_elem->GetValue(VELOCITY_LOWER);
            // it_elem->SetValue(VELOCITY_LOWER, unit_normal);
            double projection = inner_prod(unit_normal, wake_direction);

            if (projection > 0.0) {
                upper_surface_elements_list.push_back(it_elem->Id());
                for (std::size_t i_node = 0; i_node<it_elem->GetGeometry().size();i_node++) {
                    r_geometry[i_node].SetValue(UPPER_SURFACE, true);
                }
            } else {
                lower_surface_elements_list.push_back(it_elem->Id());
                for (std::size_t i_node = 0; i_node<it_elem->GetGeometry().size();i_node++) {
                    r_geometry[i_node].SetValue(LOWER_SURFACE, true);
                }
            }
        }
    }

    mrModelPart.GetSubModelPart("lower_surface_sub_model_part").AddElements(lower_surface_elements_list);
    mrModelPart.GetSubModelPart("upper_surface_sub_model_part").AddElements(upper_surface_elements_list);

    auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];
    auto reference_chord = mrModelPart.GetProcessInfo()[REFERENCE_CHORD];

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        bool is_upper = it_node->GetValue(UPPER_SURFACE);
        bool is_lower = it_node->GetValue(LOWER_SURFACE);
        // bool is_right_from_te = it_node->X() > wake_origin[0]-reference_chord/2;
        bool is_right_from_te = it_node->X() > 0.0;
        if (is_lower && is_upper && is_right_from_te && it_node->IsNot(MARKER)) {
            it_node->SetValue(TRAILING_EDGE, true);
        }

        // Also copying distance field to temperature:
        double distance_value = it_node -> FastGetSolutionStepValue(GEOMETRY_DISTANCE);
        it_node -> SetValue(TEMPERATURE, distance_value);

    }

    std::string sub_model_part_to_mark_as_kutta_name;
    if (mrTrailingEdgeNode.GetValue(WAKE_DISTANCE) > 0.0) {
        sub_model_part_to_mark_as_kutta_name = "lower_surface_sub_model_part";
    } else {
        sub_model_part_to_mark_as_kutta_name = "upper_surface_sub_model_part";
    }

    for (int i = 0; i < static_cast<int>(mrModelPart.GetSubModelPart(sub_model_part_to_mark_as_kutta_name).Elements().size()); i++) {
        auto it_elem = mrModelPart.GetSubModelPart(sub_model_part_to_mark_as_kutta_name).ElementsBegin() + i;
        bool is_wake = it_elem -> GetValue(WAKE);
        bool is_active = it_elem -> Is(ACTIVE);
        auto& r_geometry = it_elem->GetGeometry();
        if (!is_wake && is_active) {
            for (std::size_t i_node = 0; i_node<it_elem->GetGeometry().size();i_node++) {
                double distance_value = r_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);

                if (r_geometry[i_node].GetValue(TRAILING_EDGE) && distance_value<0.0) {
                    // it_elem->SetValue(KUTTA, true);
                    r_geometry[i_node].SetValue(TEMPERATURE, distance_value);
                }
            }
        }
    }
}

ModifiedShapeFunctions::Pointer DefineEmbeddedWakeProcess::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

}// Namespace Kratos


// ModifiedShapeFunctions::Pointer p_modified_geometry = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), Vector(geometry_elemental_distances));
// Matrix geometry_sh_func_values;
// ModifiedShapeFunctions::ShapeFunctionsGradientsType geometry_sh_func_gradients;
// Vector geometry_weights;
// p_modified_geometry -> ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
//     geometry_sh_func_values,
//     geometry_sh_func_gradients,
//     geometry_weights,
//     GeometryData::GI_GAUSS_1);

// Vector interface_coordinates = ZeroVector(3);
// for (std::size_t i_node = 0; i_node<3; i_node++) {
//     auto nodal_coordinates = it_elem->GetGeometry()[i_node].Coordinates();
//     for (std::size_t i_dim = 0; i_dim<3; i_dim++) {
//         interface_coordinates[i_dim] += geometry_sh_func_values(0,i_node)*nodal_coordinates[i_dim];
//     }

// }
// auto interface_vector = interface_coordinates-wake_origin;
