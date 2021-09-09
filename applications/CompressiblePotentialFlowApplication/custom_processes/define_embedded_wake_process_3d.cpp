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
        BoundedVector<double, 3> wake_direction;
        // wake_direction[0] = cos(angle_in_deg*Globals::Pi/180);
        // wake_direction[1] = sin(angle_in_deg*Globals::Pi/180);
        // wake_direction[2] = 0.0;
        wake_direction[0] = sin(angle_in_deg*Globals::Pi/180);
        wake_direction[1] = 0.0;
        wake_direction[2] = -cos(angle_in_deg*Globals::Pi/180);
        if (is_embedded && it_elem->Is(ACTIVE)) {

            BoundedVector<double, 3> unit_normal = it_elem->GetValue(VELOCITY_LOWER);
            // it_elem->SetValue(VELOCITY_LOWER, unit_normal);
            double projection = inner_prod(unit_normal, wake_direction);

            it_elem->SetValue(TEMPERATURE, unit_normal[0]*unit_normal[2]); //multiply projection to x-axis and z-axis

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
        if (is_struct){
            if ((is_upper && is_lower) or is_wake) {
                for (auto& r_node : r_geometry) {
                    r_node.SetValue(KUTTA, true);
                }
            } else {
                r_elem.Set(STRUCTURE, false);
                // r_elem.Set(MARKER, true);
                r_elem.SetValue(WAKE, false);
                if (is_lower) {
                    r_elem.SetValue(WAKE, false);
                    r_elem.Set(STRUCTURE, false);
                    r_elem.SetValue(KUTTA, true);
                    for (auto& r_node : r_geometry) {
                        if (r_node.GetValue(WAKE_DISTANCE) > 0.0) {
                            r_node.SetValue(TRAILING_EDGE, true);
                        }
                    }
                }
            }
        } else {
            // bool is_x = r_elem.GetGeometry().Center().X() > 0.0;
            if (is_upper && is_lower) {
                bool is_neighbour = false;
                for (auto& r_node : r_geometry) {
                    for (auto& r_elem_neigh : r_node.GetValue(NEIGHBOUR_ELEMENTS)) {
                        if (r_elem_neigh.Is(STRUCTURE)) {
                            is_neighbour = true;
                            break;
                        }
                    }
                }
                if (is_neighbour) {
                    r_elem.Set(MARKER, true);
                    // r_elem.SetValue(WAKE, false);
                    // r_elem.Set(STRUCTURE, false);
                    // r_elem.SetValue(KUTTA, true);
                    // for (auto& r_node : r_geometry) {
                    //     if (r_node.GetValue(WAKE_DISTANCE) > 0.0) {
                    //         r_node.SetValue(TRAILING_EDGE, true);
                    //     }
                    // }
                }

            }
        }
    }

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
    outfile_wake.open("marker_elements_id.txt");
    std::ofstream outfile_structure;
    outfile_structure.open("structure_elements_id.txt");
    std::ofstream outfile_kutta;
    outfile_kutta.open("kutta_elements_id.txt");
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
    }
    outfile_kutta.close();
    outfile_marker.close();
    outfile_structure.close();
    outfile_wake.close();

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
            std::vector<Vector> cut_normal;
            tetrehedra_shape_func.ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
            double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
            array_1d<double,3> unit_normal = cut_normal[0]/norm_normal;
            rElem.SetValue(VELOCITY_LOWER,unit_normal);
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
