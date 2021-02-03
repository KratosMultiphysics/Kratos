//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "nicola.h"


namespace Kratos
{
    /* Public functions *******************************************************/

    // void CleaningUtilities::CleanIsolatedNodes(){

    //     const int initial_num = mrModelPart.Nodes().size();
    //     auto& r_nodes_array = mrModelPart.Nodes();
    //     const auto& r_elem_array = mrModelPart.Elements();

    //     // marking all nodes as "superfluous"
    //     #pragma omp parallel for
    //     for( int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node ){
    //         auto node = r_nodes_array.begin() + i_node;
    //         node->Set(TO_ERASE, true);
    //     }

    //     // saving the nodes that belong to an element
    //     #pragma omp parallel for
    //     for( int i_elem = 0; i_elem < static_cast<int>(r_elem_array.size()); ++i_elem ){
    //         const auto elem = r_elem_array.begin() + i_elem;
    //         auto& r_geom = elem->GetGeometry();

    //         for (unsigned int i = 0; i < r_geom.size(); ++i){
    //             r_geom[i].Set(TO_ERASE, false);
    //         }
    //     }

    //     mrModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    //     const int final_num = mrModelPart.Nodes().size();
    //     KRATOS_INFO("CleaningUtilities") << "In total " << (initial_num - final_num) <<" superfluous nodes were cleared" << std::endl;

    // }

    void Nicola::test_nicola(){
        int a = 1;
        int b = 2;
        int c = a+b;
        KRATOS_INFO("Nicola") << "test_nicola!" << std::endl;
    }


    // ModelPart& Nicola::CheckIfInternal( ModelPart& VolumeModelPart, ModelPart& GeometryModelPart ){
    ModelPart& Nicola::CheckIfInternal( ModelPart& VolumeModelPart ){

        auto& r_nodes_array = VolumeModelPart.Nodes();
        
        return r_nodes_array;

        // // double array_z[VolumeModelPart.NumberOfNodes()] = {};       // length of array = VolumeModelPart.NumberOfNodes()
        // std::vector<std::size_t> coord_z;

        // // loop to compute the maximum value of z
        // for (auto& node : VolumeModelPart.Nodes()){
        //     coord_z.push_back(node->Coordinates()[2]);
        // }
        // const max_z = *std::max_element(coord_z.begin(), coord_z.end());    // max value of z

        // for (auto& node : VolumeModelPart.Nodes()){
        //     const array_1d<double,3> p_coords = node.Coordinates();
        //     const array_1d<double,3> t_coords = {p_coords[0],
        //                                          p_coords[1],
        //                                          p_coords[2]+max_z};
        // }
    }
}
