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
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{
// Constructor for DefineEmbeddedWakeProcess3D Process
DefineEmbeddedWakeProcess3D::DefineEmbeddedWakeProcess3D(ModelPart& rModelPart,
                    ModelPart& rWakeModelPart
                ):
    Process(),
    mrModelPart(rModelPart),
    mrWakeModelPart(rWakeModelPart)
{}

void DefineEmbeddedWakeProcess3D::Execute()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.GetProcessInfo()[DOMAIN_SIZE] != 3) << "DOMAIN_SIZE is not 3. DefineEmbeddedWakeProcess3D is only implemented for 3D cases!" << std::endl;

    ExecuteInitialize();
    ComputeDistanceToWake();
    MarkWakeElements();

    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = mrModelPart.ElementsBegin() + i;
        BoundedVector<double, 4> geometry_distances;
        auto& r_geometry = it_elem->GetGeometry();
        for(unsigned int i_node = 0; i_node< r_geometry.size(); i_node++){
            geometry_distances[i_node] = r_geometry[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<3,4>(geometry_distances);

        auto angle_in_deg = -1*mrModelPart.GetProcessInfo()[ROTATION_ANGLE];
        auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];
        BoundedVector<double, 3> wake_normal;
        // wake_direction[0] = cos(angle_in_deg*Globals::Pi/180);
        // wake_direction[1] = sin(angle_in_deg*Globals::Pi/180);
        // wake_direction[2] = 0.0;
        wake_normal[0] = sin(angle_in_deg*Globals::Pi/180);
        wake_normal[1] = 0.0;
        wake_normal[2] = -cos(angle_in_deg*Globals::Pi/180);
        if (is_embedded && it_elem->Is(ACTIVE)) {

            // INNER UNIT NORMAL (points inside volume)
            BoundedVector<double, 3> unit_normal = it_elem->GetValue(VELOCITY_LOWER);
            // it_elem->SetValue(VELOCITY_LOWER, unit_normal);
            double projection = inner_prod(unit_normal, wake_normal);

            it_elem->SetValue(TEMPERATURE, std::abs(unit_normal[1])); //multiply projection to x-axis and z-axis

            double side_limit = 0.001;
            if (projection > side_limit) {
                // upper_surface_elements_list.push_back(it_elem->Id());
                for (std::size_t i_node = 0; i_node<it_elem->GetGeometry().size();i_node++) {
                    r_geometry[i_node].SetValue(UPPER_SURFACE, true);
                }
            }
            if (projection < -side_limit) {
                // lower_surface_elements_list.push_back(it_elem->Id());
                for (std::size_t i_node = 0; i_node<it_elem->GetGeometry().size();i_node++) {
                    r_geometry[i_node].SetValue(LOWER_SURFACE, true);
                }
            }
        }
    }

    for (auto& r_elem : mrModelPart.Elements()) {
        bool is_upper = false;
        bool is_lower = false;
        bool is_wake = false;
        auto& r_geometry = r_elem.GetGeometry();
        for (auto& r_node : r_geometry) {
            if (r_node.GetValue(UPPER_SURFACE) == 1) {
                is_upper = true;
            }
            if (r_node.GetValue(LOWER_SURFACE) == 1) {
                is_lower = true;
            }
            if (r_node.GetValue(WAKE) == 1) {
                is_wake = true;
            }
        }

        bool is_struct = r_elem.Is(STRUCTURE);
        double y_projection = r_elem.GetValue(TEMPERATURE);
        if (is_struct){
            if (y_projection < 0.99999) {
                if ((is_upper && is_lower) or is_wake) {
                    for (auto& r_node : r_geometry) {
                        r_node.SetValue(KUTTA, true);
                    }
                } else {
                    r_elem.Set(STRUCTURE, false);
                    // r_elem.Set(MARKER, true);
                    r_elem.SetValue(WAKE, false);
                    auto& wake_elemental_distances = r_elem.GetValue(WAKE_ELEMENTAL_DISTANCES);
                    if (is_upper) {
                        r_elem.SetValue(KUTTA, true);
                        for (IndexType i=0; i<r_geometry.size(); i++) {
                            if (wake_elemental_distances[i] < 0.0) {
                                r_geometry[i].SetValue(TRAILING_EDGE, true);
                            }
                        }
                    }
                    if (is_lower) {
                        r_elem.SetValue(KUTTA, true);
                        for (IndexType i=0; i<r_geometry.size(); i++) {
                            if (wake_elemental_distances[i] > 0.0) {
                                r_geometry[i].SetValue(TRAILING_EDGE, true);
                            }
                        }
                    }
                }
            } else {
                r_elem.Set(STRUCTURE, false);
                r_elem.SetValue(WAKE, false);
            }
        }
    }


    // KRATOS_INFO("DefineEmbeddedWakeProcess3D")
    //     << "Deactivating unconnected structure elements" << std::endl;
    // std::unordered_set<std::size_t> touching_elements;

    // // ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");
    // for (auto& r_elem : mrModelPart.Elements()) {
    //     if (r_elem.Is(STRUCTURE) &&
    //         touching_elements.find(r_elem.Id()) == touching_elements.end()) {
    //         // KRATOS_WATCH("***********************************")
    //         // KRATOS_WATCH(r_elem.Id())
    //         std::unordered_set<std::size_t> visited_elements;
    //         bool touches_wake = TouchesWake(r_elem, visited_elements);
    //         // KRATOS_WATCH("***********************************")
    //         // KRATOS_WATCH(touches_wake)
    //         // KRATOS_WATCH("***********************************")
    //         if (!touches_wake) {
    //             r_elem.Set(MARKER, true);
    //             r_elem.SetValue(WAKE, false);
    //             r_elem.Set(STRUCTURE, false);
    //             // for (auto& r_node : r_elem.GetGeometry()) {
    //             //     r_node.SetValue(WING_TIP, false);
    //             // }
    //         }
    //         else {
    //             for (auto elem_id : visited_elements) {
    //                 touching_elements.insert(elem_id);
    //             }
    //             for (auto& r_node : r_elem.GetGeometry()) {
    //                 r_node.SetValue(WAKE, true);
    //             }
    //         }
    //     }
    // }
    // KRATOS_INFO("DefineEmbeddedWakeProcess3D")
    //     << "Finished deactivating unconnected structure elements" << std::endl;

    // for (auto& r_elem : mrModelPart.Elements()) {
    //     bool is_x = r_elem.GetGeometry().Center().X() > 0.0;
    //     // if (r_elem.Is(MARKER)) {
    //     if (is_x) {
    //         for (auto& r_node : r_elem.GetGeometry()) {
    //             // if (r_node.GetValue(UPPER_SURFACE) && r_node.GetValue(LOWER_SURFACE) && !r_node.GetValue(WAKE)) {
    //             if (r_node.GetValue(UPPER_SURFACE) && r_node.GetValue(LOWER_SURFACE)) {
    //                 r_node.SetValue(TRAILING_EDGE, true);
    //             }
    //             // if (r_node.GetValue(WAKE_DISTANCE) > 0.0) {
    //             //     r_node.SetValue(TRAILING_EDGE, true);
    //             // }
    //         }
    //     }
    // }

    // for (auto& r_elem : mrModelPart.Elements()) {
    //     std::size_t upper_nodes = 0;
    //     std::size_t lower_nodes = 0;
    //     bool touches_te = false;
    //     for (auto& r_node : r_elem.GetGeometry()) {
    //         if (r_node.GetValue(UPPER_SURFACE)) {
    //             upper_nodes++;
    //         }
    //         if (r_node.GetValue(LOWER_SURFACE)) {
    //             lower_nodes++;
    //         }
    //         if (r_node.GetValue(TRAILING_EDGE))
    //             touches_te = true;
    //     }
    //     bool is_lower =lower_nodes >= upper_nodes;
    //     // if (is_upper && touches_te) {
    //     if (!is_lower && touches_te) {
    //         r_elem.SetValue(KUTTA, true);
    //     }
    //     // if (touches_te && r_elem.Is(ACTIVE)) {
    //     //     for (auto& r_node : r_elem.GetGeometry()) {
    //     //         r_node.SetValue(WING_TIP, true);
    //     //     }
    //     // }
    // }

    // for (auto& r_elem : mrModelPart.Elements()) {
    //     bool is_found1 = false;
    //     for (auto& r_node : r_elem.GetGeometry()) {
    //         if (r_node.Id() == 102034)
    //             is_found1 = true;
    //     }
    //     bool is_found2 = false;
    //     for (auto& r_node : r_elem.GetGeometry()) {
    //         if (r_node.Id() == 134913)
    //             is_found2 = true;

    //     }
    //     bool is_found3 = false;
    //     for (auto& r_node : r_elem.GetGeometry()) {
    //         if (r_node.Id() == 193859)
    //             is_found3 = true;
    //     }
    //     if (is_found1 && is_found2 && is_found3) {
    //         KRATOS_WATCH("FOUND")
    //         KRATOS_WATCH(r_elem.Id())
    //     }
    // }

    ComputeTrailingEdgeNode();

    // for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
    //     ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i;
    //     auto& r_geometry = it_elem->GetGeometry();


    //     BoundedVector<double,4> geometry_distances;
    //     for(unsigned int i_node = 0; i_node<r_geometry.size(); i_node++){
    //         geometry_distances[i_node] = r_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
    //     }
    //     const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<3,4>(geometry_distances);
    //     bool is_touching_marker = false;
    //     if (is_embedded){
    //         for (auto& r_node : r_geometry) {
    //             if (r_node.GetValue(KUTTA) == 1) {
    //                 // std::cout << "Found embedded elem " << it_elem->Id() << " touching marker" << std::endl;
    //                 is_touching_marker = true;
    //                 break;
    //             }
    //         }
    //     }

    //     if (is_touching_marker) {
    //         for (auto& r_node : r_geometry) {
    //             r_node.SetValue(KUTTA, true);
    //         }
    //     }
    // }


    std::ofstream outfile_wake;
    outfile_wake.open("wake_elements_id.txt");
    std::ofstream outfile_marker;
    outfile_marker.open("marker_elements_id.txt");
    std::ofstream outfile_structure;
    outfile_structure.open("structure_elements_id.txt");
    std::ofstream outfile_kutta;
    outfile_kutta.open("kutta_elements_id.txt");
    std::ofstream outfile_upperlower;
    outfile_upperlower.open("upperlower_elements_id.txt");
    for (auto& r_element : mrModelPart.Elements()){
        if(r_element.Is(MARKER)){
            outfile_marker << r_element.Id();
            outfile_marker << " ";
        }
        if(!r_element.GetValue(WAKE)){
            if(r_element.GetValue(KUTTA)){
                outfile_kutta << r_element.Id();
                outfile_kutta << " ";
            }
        }
        else{
            outfile_wake << r_element.Id();
            outfile_wake << " ";
            if(r_element.Is(STRUCTURE)){
                outfile_structure << r_element.Id();
                outfile_structure << " ";
            }
        }
        std::size_t upper_nodes = 0;
        std::size_t lower_nodes = 0;
        for (auto& r_node : r_element.GetGeometry()) {
            if (r_node.GetValue(UPPER_SURFACE)) {
                upper_nodes++;
            }
            if (r_node.GetValue(LOWER_SURFACE)) {
                lower_nodes++;
            }
        }
        if (upper_nodes==4 && lower_nodes==4 && r_element.GetValue(WAKE)) {
                outfile_upperlower << r_element.Id();
                outfile_upperlower << " ";
        }
    }
    outfile_kutta.close();
    outfile_marker.close();
    outfile_structure.close();
    outfile_wake.close();
    outfile_upperlower.close();
    KRATOS_CATCH("");
}

void DefineEmbeddedWakeProcess3D::ExecuteInitialize() {

    KRATOS_TRY;

    block_for_each(mrModelPart.Elements(), [&](Element& rElem) {
        rElem.SetValue(WAKE, false);
    });
    block_for_each(mrModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        rNode.SetValue(WAKE_DISTANCE, 0.0);
        rNode.SetValue(WAKE, false);
        rNode.SetValue(TRAILING_EDGE, false);
        rNode.SetValue(KUTTA, false);
    });

    KRATOS_CATCH("");
}

void DefineEmbeddedWakeProcess3D::ComputeDistanceToWake(){

    KRATOS_TRY;

    CalculateDiscontinuousDistanceToSkinProcess<3> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();

    KRATOS_CATCH("");
}

void DefineEmbeddedWakeProcess3D::MarkWakeElements(){

    KRATOS_TRY;

    block_for_each(mrModelPart.Elements(), [&](Element& rElem)
    {
        auto& r_geometry = rElem.GetGeometry();
        BoundedVector<double, 4> nodal_distances_to_wake = rElem.GetValue(ELEMENTAL_DISTANCES);
        rElem.SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);

        // Selecting the cut (wake) elements
        const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<3,4>(nodal_distances_to_wake);

        BoundedVector<double, 4> geometry_distances;
        for(unsigned int i_node = 0; i_node< r_geometry.size(); i_node++){
            geometry_distances[i_node] = r_geometry[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<3,4>(geometry_distances);

        if (is_embedded){
            Tetrahedra3D4ModifiedShapeFunctions tetrehedra_shape_func(rElem.pGetGeometry(), Vector(geometry_distances));
            // Computing Normal
            std::vector<array_1d<double,3>> cut_normal;
            tetrehedra_shape_func.ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
            double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
            array_1d<double,3> unit_normal = cut_normal[0]/norm_normal;
            rElem.SetValue(VELOCITY_LOWER,unit_normal);
        }

        if (rElem.Id()==2899496) {
            KRATOS_WATCH(nodal_distances_to_wake)
            KRATOS_WATCH(rElem.Is(TO_SPLIT))

        }
        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element && rElem.Is(ACTIVE)) {
            rElem.SetValue(WAKE, true);
            if (is_embedded){
                rElem.Set(STRUCTURE, true);
                for (unsigned int i = 0; i < r_geometry.size(); i++) {
                    r_geometry[i].SetLock();
                    // r_geometry[i].SetValue(KUTTA, true);
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


bool DefineEmbeddedWakeProcess3D::TouchesWake(Element& rElem,
                                              std::unordered_set<std::size_t> visited_elements)
{
    // KRATOS_WATCH(rElem.Id())
    auto& r_geometry = rElem.GetGeometry();
    visited_elements.insert(rElem.Id());
    std::unordered_set<std::size_t> failed_elements;
    double global_max_node_x = -1e10;
    for (auto& r_node : rElem.GetGeometry()) {
        global_max_node_x = std::max(global_max_node_x, r_node.X());
    }
    for (auto& r_node : r_geometry) {
        std::vector<Element> vector_wake_elem;
        std::vector<Element> sorted_vector_wake_elem;
        for (auto& r_elem_neigh : r_node.GetValue(NEIGHBOUR_ELEMENTS)) {
            std::size_t matching_nodes = 0;

            for (auto& r_node_neigh : r_elem_neigh.GetGeometry()) {
                for (auto& r_node_main : rElem.GetGeometry()) {
                    if (r_node_main.Id() == r_node_neigh.Id()) {
                        matching_nodes++;
                    }
                }
            }
            bool is_matching = matching_nodes == 3;
            bool is_wake = r_elem_neigh.GetValue(WAKE);
            if (is_matching && is_wake) {
                vector_wake_elem.push_back(r_elem_neigh);
            }
        }
        while (vector_wake_elem.size()>0) {
            IndexType max_i = 0;
            double max_x_node = -1e10;
            for (IndexType i =0 ; i<vector_wake_elem.size(); i++) {
                for (auto& r_node : vector_wake_elem[i].GetGeometry()) {
                    if (r_node.X()>max_x_node) {
                        max_x_node = r_node.X();
                        max_i = i;
                    }
                }
            }
            sorted_vector_wake_elem.push_back(vector_wake_elem[max_i]);
            vector_wake_elem.erase(vector_wake_elem.begin()+ max_i);
            // KRATOS_WATCH(vector_wake_elem.size())
        }

        for (auto& r_elem_neigh : sorted_vector_wake_elem) {
            if (failed_elements.find(r_elem_neigh.Id()) == failed_elements.end() &&
                r_elem_neigh.GetValue(WAKE)) {
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
                bool is_matching = matching_nodes ==  3;
                bool is_visited = visited_elements.find(r_elem_neigh.Id()) !=
                                  visited_elements.end();
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
                    }
                    else {
                        failed_elements.insert(r_elem_neigh.Id());
                    }
                }
            }
        }
    }
    return false;

}

// std::vector<IndexType> DefineEmbeddedWakeProcess3D::pGetTrailingEdgeNode(){

//     for (auto& r_node : mrModelPart.Nodes()) {
//         bool is_positive = r_node.GetValue(WAKE_DISTANCE) > 0.0;
//         bool is_wake = r_node.GetValue(WAKE);
//         bool is_kutta = r_node.GetValue(KUTTA);
//         if (is_positive && is_wake && is_kutta) {
//             r_node.SetValue(TRAILING_EDGE, true);
//             return &r_node;
//         }
//     }
//     KRATOS_ERROR << "No trailing edge node was found!" << std::endl;
//     return nullptr;

// }

void DefineEmbeddedWakeProcess3D::ComputeTrailingEdgeNode(){

    KRATOS_TRY;

    // auto p_node = pGetTrailingEdgeNode();

    std::vector<std::size_t> trailing_edge_node_list;

    // for (auto& r_node : mrModelPart.Nodes()) {
    //     bool is_positive = r_node.GetValue(WAKE_DISTANCE) > 0.0;
    //     bool is_wake = r_node.GetValue(WAKE);
    //     bool is_kutta = r_node.GetValue(KUTTA);
    //     if (is_positive && is_wake && is_kutta) {
    //         r_node.SetValue(TRAILING_EDGE, true);
    //         // return &r_node;
    //         trailing_edge_node_list.push_back(r_node.Id());
    //     }
    // }
    for (auto& r_node : mrModelPart.Nodes()) {
        bool is_te = r_node.GetValue(TRAILING_EDGE);
        if (is_te) {
            trailing_edge_node_list.push_back(r_node.Id());
        }
    }
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
